#ifndef KALMAN_FILTER_H_
#define KALMAN_FILTER_H_

#include "Eigen/Dense"
#include "Eigen/LU"

/**
 * Generic Kalman Filter class.
 *
 * The filter maintains the system state and the variance as its state. The
 * updates to the state are managed by Sensor objects. Multiple Sensor objects
 * can reference the same KalmanFilter.
 *
 * The template class is parameterised by the dimension of the state vector.
 * Each Sensor template class is also parameterised by the dimension of the
 * measurement vector. This makes it possible to use fixed-size matrix math
 * for improved compile-time correctness checking and (sometimes) better
 * performance.
 */
template <int StateSize>
struct KalmanFilter {
  typedef Eigen::Matrix<double, StateSize, 1> StateVector;
  typedef Eigen::Matrix<double, StateSize, StateSize> StateMatrix;

  /**
   * A generic sensor for a Kalman Filter. The Sensor is responsible for
   * updating the state of the filter based on its measurements.
   */
  template <int MeasurementSize>
  struct Sensor {
    typedef Eigen::Matrix<double, MeasurementSize, 1> MeasurementVector;
    typedef Eigen::Matrix<double, MeasurementSize, MeasurementSize>
      MeasurementMatrix;
    typedef Eigen::Matrix<double, MeasurementSize, StateSize>
      MeasurementStateMatrix;
    typedef Eigen::Matrix<double, StateSize, MeasurementSize>
      StateMeasurementMatrix;

    /**
     * Constructor: provide a reference to the KalmanFilter that this Sensor
     * updates.
     */
    Sensor(KalmanFilter<StateSize> &filter) : filter_(filter) { }

    /**
     * Update a linear Kalman Filter with a new measurement.
     *
     * @param z measurement
     * @param H measurement matrix
     * @param R measurement covariance matrix
     */
    void Update(
      const MeasurementVector &z,
      const MeasurementStateMatrix &H,
      const MeasurementMatrix &R) {
      UpdateEKF(z, H * filter_.state(), H, R);
    }

    /**
     * Update an Extended Kalman Filter with a new measurement.
     *
     * @param z measurement
     * @param h measurement predicted (possibly non-linearly) from the state
     * @param H measurement matrix (possibly linearized)
     * @param R measurement covariance matrix
     */
    void UpdateEKF(
      const MeasurementVector &z,
      const MeasurementVector &h,
      const MeasurementStateMatrix &H,
      const MeasurementMatrix &R)
    {
      const StateMatrix &P = filter_.covariance();
      MeasurementVector y = z - h;
      MeasurementMatrix S = H * P * H.transpose() + R;

      // Avoid explicit matrix inversion. Starting from
      // K = P H^T S^-1
      // postmultiply by S to get
      // K S = P H^T
      // which is equivalent to
      // S^T K^T = H P^T
      // which is in standard form (Ax=B)
      Eigen::FullPivLU<MeasurementMatrix> lu(S.transpose());
      StateMeasurementMatrix K = lu.solve(H * P.transpose()).transpose();

      filter_.Update(K * y, K * H);
    }

  protected:
    KalmanFilter<StateSize> &filter_;
  };

  /**
   * Constructor: the state and covariance matrices are uninitialized (may
   * contain values from uninitialized memory).
   */
  KalmanFilter() : KalmanFilter(StateVector(), StateMatrix()) {}

  /**
   * Constructor: set the initial state and covariance matrix.
   */
  KalmanFilter(const StateVector &x, const StateMatrix &P)
    : I_(StateMatrix::Identity()), x_(x), P_(P)
  { }

  const StateVector &state() const { return x_; }

  const StateMatrix &covariance() const { return P_; }

  /**
   * Predict the next state
   *
   * @param F state transition matrix
   * @param Q process covariance matrix
   */
  void Predict(const StateMatrix &F, const StateMatrix &Q) {
    x_ = F * x_;
    P_ = F * P_ * F.transpose() + Q;
  }

  /**
   * Initialize (or reinitialize) the state of the filter.
   *
   * @param x state
   * @param P covariance matrix
   */
  void Initialize(const StateVector &x, const StateMatrix &P) {
    x_ = x;
    P_ = P;
  }

  /**
   * Update the state and covariance estimates. The caller is responsible for
   * computing the Kalman gain and using it to scale the residual and the
   * measurement matrix. This allows this method to be independent of the
   * dimension of the measurement vector.
   *
   * @param Ky residual, premultiplied by the Kalman gain, K
   * @param KH measurement matrix, premultiplied by the Kalman gain, K
   */
  void Update(const StateVector &Ky, const StateMatrix &KH) {
    x_ = x_ + Ky;
    P_ = (I_ - KH) * P_;
  }

private:
  // Identity matrix (for Update)
  StateMatrix I_;

  // State vector
  StateVector x_;

  // State covariance matrix
  StateMatrix P_;
};

#endif /* KALMAN_FILTER_H_ */

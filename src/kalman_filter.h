#ifndef KALMAN_FILTER_H_
#define KALMAN_FILTER_H_

#include "Eigen/Dense"
#include "Eigen/LU"

template <int StateSize>
struct KalmanFilter {
  typedef Eigen::Matrix<double, StateSize, 1> StateVector;
  typedef Eigen::Matrix<double, StateSize, StateSize> StateMatrix;

  template <int MeasurementSize>
  struct Sensor {
    typedef Eigen::Matrix<double, MeasurementSize, 1> MeasurementVector;
    typedef Eigen::Matrix<double, MeasurementSize, MeasurementSize>
      MeasurementMatrix;
    typedef Eigen::Matrix<double, MeasurementSize, StateSize>
      MeasurementStateMatrix;
    typedef Eigen::Matrix<double, StateSize, MeasurementSize>
      StateMeasurementMatrix;

    Sensor(KalmanFilter<StateSize> &filter) : filter_(filter) { }

    void Update(
      const MeasurementVector &z,
      const MeasurementStateMatrix &H,
      const MeasurementMatrix &R) {
      UpdateEKF(z, H * filter_.state(), H, R);
    }

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

  KalmanFilter() : KalmanFilter(StateVector(), StateMatrix()) {}

  KalmanFilter(const StateVector &x, const StateMatrix &P)
    : I_(StateMatrix::Identity()), x_(x), P_(P)
  { }

  const StateVector &state() const { return x_; }

  const StateMatrix &covariance() const { return P_; }

  void Predict(const StateMatrix &F, const StateMatrix &Q) {
    x_ = F * x_;
    P_ = F * P_ * F.transpose() + Q;
  }

  void Initialize(const StateVector &x, const StateMatrix &P) {
    x_ = x;
    P_ = P;
  }

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

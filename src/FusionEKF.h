#ifndef FusionEKF_H_
#define FusionEKF_H_

#include "measurement_package.h"
#include "Eigen/Dense"
#include <vector>
#include <string>
#include <fstream>
#include "kalman_filter.h"
#include "tools.h"

/**
 * Sensor-fusing Extended Kalman Filter for pedestrian tracking with laser
 * (LIDAR) and radar data.
 */
class FusionEKF {
public:
  /**
  * Constructor.
  */
  FusionEKF();

  /**
  * Destructor.
  */
  virtual ~FusionEKF();

  /**
  * Run the whole flow of the Kalman Filter from here.
  */
  void ProcessMeasurement(const MeasurementPackage &measurement_pack);

  // State accessors:
  const Eigen::Vector4d &state() const { return ekf_.state(); }
  double px() const { return state()(0); }
  double py() const { return state()(1); }
  double vx() const { return state()(2); }
  double vy() const { return state()(3); }

private:
  typedef KalmanFilter<4> Filter;

  /**
   * The radar sensor: handles the nonlinear measurement function and Jacobian
   * calculation for updating the filter. If the first measurement is a radar
   * measurement, this class also handles initializing the filter using that
   * first radar measurement.
   *
   * The radar returns a three-dimensional measurement vector: range (rho),
   * angle (phi) and radial speed (rho_dot).
   */
  struct Radar : public Filter::Sensor<3>
  {
    explicit Radar(Filter &filter);

    void Initialize(const MeasurementVector &z);

    void Update(const MeasurementVector &z);

  private:
    // Measurement covariance matrix for radar
    MeasurementMatrix R_;

    // Approximate measurement matrix for radar
    MeasurementStateMatrix Hj_;
  };

  /**
   * The laser sensor: much simpler than the Radar sensor, because the
   * measurement function is linear. If the first measurement is a laser
   * measurement, this class also handles initializing the filter using that
   * first laser measurement.
   *
   * The laser returns a two-dimensional measurement vector: x and y position.
   */
  struct Laser : public Filter::Sensor<2>
  {
    explicit Laser(Filter &filter);

    void Initialize(const MeasurementVector &z);

    void Update(const MeasurementVector &z);

  private:
    // Measurement covariance matrix for laser
    MeasurementMatrix R_;

    // Measurement matrix for laser
    MeasurementStateMatrix H_;
  };

  // check whether the tracking toolbox was initiallized or not (first measurement)
  bool is_initialized_;

  // previous timestamp
  long previous_timestamp_;

  // We'll have a single Filter object, and it will be updated by the two
  // Sensor objects, Radar and Laser, which represent our two sensors to be
  // fused.
  Filter ekf_;
  Radar radar_;
  Laser laser_;

  // state transistion matrix
  Eigen::Matrix4d F_;

  // process covariance matrix
  Eigen::Matrix4d Q_;

  void Initialize(const MeasurementPackage &measurement_pack);
  void Predict(const MeasurementPackage &measurement_pack);
  void Update(const MeasurementPackage &measurement_pack);
};

#endif /* FusionEKF_H_ */

#ifndef FusionEKF_H_
#define FusionEKF_H_

#include "measurement_package.h"
#include "Eigen/Dense"
#include <vector>
#include <string>
#include <fstream>
#include "kalman_filter.h"
#include "tools.h"

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

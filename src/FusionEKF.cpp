#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
        0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
        0, 0.0009, 0,
        0, 0, 0.09;

  /**
    * Finish initializing the FusionEKF.
    * Set the process and measurement noises
  */
  ekf_.F_ = MatrixXd::Identity(4, 4); // transition (dt values updated later)

  H_laser_ <<
    1, 0, 0, 0,
    0, 1, 0, 0;
}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::Initialize(const MeasurementPackage &measurement_pack) {
  /**
    * Initialize the state ekf_.x_ with the first measurement.
    * Create the covariance matrix.
  */
  double rho, phi, x, y;
  switch (measurement_pack.sensor_type_) {
    case MeasurementPackage::RADAR:
      // Convert radar from polar to cartesian coordinates.
      rho = measurement_pack.raw_measurements_(0);
      phi = measurement_pack.raw_measurements_(1);
      x = rho * cos(phi);
      y = rho * sin(phi);
      break;
    case MeasurementPackage::LASER:
      x = measurement_pack.raw_measurements_(0);
      y = measurement_pack.raw_measurements_(1);
      break;
    default:
      cerr << "bad measurement pack sensor type " <<
        measurement_pack.sensor_type_ << endl;
      exit(EXIT_FAILURE);
  }

  ekf_.x_ = VectorXd(4);
  ekf_.x_ << x, y, 0, 0;

  ekf_.P_ = MatrixXd(4, 4);
	ekf_.P_ <<
    1, 0,    0,    0,
	  0, 1,    0,    0,
	  0, 0, 1000,    0,
	  0, 0,    0, 1000;
}

void FusionEKF::Predict(const MeasurementPackage &measurement_pack) {
  /**
   * Update the state transition matrix F according to the new elapsed time.
   * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */
  // Get elapsed time in seconds (from microseconds).
  double dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1e6;

  //1. Modify the F matrix so that the time is integrated
  cout << "dt=" << dt << endl;
  ekf_.F_(0, 2) = dt;
  ekf_.F_(1, 3) = dt;

	//2. Set the process covariance matrix Q
	double dt2 = dt * dt;
	double dt3 = dt2 * dt / 2.0f;
	double dt4 = dt2 * dt2 / 4.0f;
  double noise_ax = 9;
  double noise_ay = 9;
	ekf_.Q_ = MatrixXd(4, 4);
	ekf_.Q_ <<
	  dt4 * noise_ax,              0, dt3 * noise_ax,              0,
	               0, dt4 * noise_ay,              0, dt3 * noise_ay,
	  dt3 * noise_ax,              0, dt2 * noise_ax,              0,
	               0, dt3 * noise_ay,              0, dt2 * noise_ay;

  ekf_.Predict();
}

void FusionEKF::Update(const MeasurementPackage &measurement_pack) {
  /**
   * Use the sensor type to perform the update step.
   * Update the state and covariance matrices.
   */
  switch(measurement_pack.sensor_type_) {
    case MeasurementPackage::RADAR:
      UpdateWithRadar(measurement_pack);
    break;
    case MeasurementPackage::LASER:
      UpdateWithLaser(measurement_pack);
      break;
    default:
      cerr << "bad measurement pack sensor type " <<
        measurement_pack.sensor_type_ << endl;
      exit(EXIT_FAILURE);
  }
}

void FusionEKF::UpdateWithRadar(const MeasurementPackage &measurement_pack) {
  ekf_.H_ = tools.CalculateJacobian(ekf_.x_);
  ekf_.R_ = R_radar_;

  double px = ekf_.x_(0);
  double py = ekf_.x_(1);
  double vx = ekf_.x_(2);
  double vy = ekf_.x_(3);
  double phi = atan2(py, px);
  double rho = sqrt(px*px + py*py);
  if (abs(rho) < 1e-9) {
    return;
  }
  double rho_dot = (px * vx + py * vy) / rho;
  VectorXd h(3);
  h << rho, phi, rho_dot;

  ekf_.UpdateEKF(measurement_pack.raw_measurements_, h);
}

void FusionEKF::UpdateWithLaser(const MeasurementPackage &measurement_pack) {
  ekf_.H_ = H_laser_;
  ekf_.R_ = R_laser_;
  ekf_.Update(measurement_pack.raw_measurements_);
}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
  if (is_initialized_) {
    Predict(measurement_pack);
    Update(measurement_pack);
  } else {
    cout << "EKF: " << endl;
    Initialize(measurement_pack);
    is_initialized_ = true;
  }

  previous_timestamp_ = measurement_pack.timestamp_;

  cout << "x_ = " << endl << ekf_.x_ << endl;
  cout << "P_ = " << endl << ekf_.P_ << endl;
}

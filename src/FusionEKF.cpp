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

double norm(double x, double y) {
  return sqrt(x*x + y*y);
}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {


  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    /**
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
    // first measurement
    cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);

    float x, y;
    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      float rho = measurement_pack.raw_measurements_(0);
      float phi = measurement_pack.raw_measurements_(1);
      x = rho * cos(phi);
      y = rho * sin(phi);
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      x = measurement_pack.raw_measurements_(0);
      y = measurement_pack.raw_measurements_(1);
    } else {
      cerr << "bad measurement pack sensor type " <<
        measurement_pack.sensor_type_ << endl;
      exit(EXIT_FAILURE);
    }
    ekf_.x_ << x, y, 0, 0;

    ekf_.P_ = MatrixXd(4, 4);
  	ekf_.P_ <<
      1, 0,    0,    0,
		  0, 1,    0,    0,
		  0, 0, 1000,    0,
		  0, 0,    0, 1000;

    // done initializing, no need to predict or update
    is_initialized_ = true;
    previous_timestamp_ = measurement_pack.timestamp_;
    cout << "x_ = " << endl << ekf_.x_ << endl;
    cout << "P_ = " << endl << ekf_.P_ << endl;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */
  //1. Modify the F matrix so that the time is integrated
  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1e6;
  cout << "dt=" << dt << endl;
  ekf_.F_(0, 2) = dt;
  ekf_.F_(1, 3) = dt;

	//2. Set the process covariance matrix Q
	float dt2 = dt * dt;
	float dt3 = dt2 * dt / 2.0f;
	float dt4 = dt2 * dt2 / 4.0f;
  float noise_ax = 9;
  float noise_ay = 9;
	ekf_.Q_ = MatrixXd(4, 4);
	ekf_.Q_ <<
	  dt4 * noise_ax,              0, dt3 * noise_ax,              0,
	               0, dt4 * noise_ay,              0, dt3 * noise_ay,
	  dt3 * noise_ax,              0, dt2 * noise_ax,              0,
	               0, dt3 * noise_ay,              0, dt2 * noise_ay;

  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
    ekf_.H_ = tools.CalculateJacobian(ekf_.x_);
    ekf_.R_ = R_radar_;
    double phi = atan2(ekf_.x_(1), ekf_.x_(0));
    double rho = norm(ekf_.x_(0), ekf_.x_(1));
    double rho_dot = (ekf_.x_(0) * ekf_.x_(2) + ekf_.x_(1) * ekf_.x_(3)) / rho;
    VectorXd h(3);
    h << rho, phi, rho_dot;
    ekf_.UpdateEKF(measurement_pack.raw_measurements_, h);
  } else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
    // Laser updates
    ekf_.H_ = H_laser_;
    ekf_.R_ = R_laser_;
    ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  cout << "x_ = " << endl << ekf_.x_ << endl;
  cout << "P_ = " << endl << ekf_.P_ << endl;
  previous_timestamp_ = measurement_pack.timestamp_;
}

#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

// Tolerance for avoiding division by zero.
const double EPSILON = 1e-6;

// Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
const double NOISE_AX = 9;
const double NOISE_AY = 9;

FusionEKF::Radar::Radar(FusionEKF::Filter &filter) :
  FusionEKF::Filter::Sensor<3>(filter)
{
  // Measurement covariance matrix:
  R_ <<
    0.09, 0, 0,
    0, 0.0009, 0,
    0, 0, 0.09;
}

void FusionEKF::Radar::Update(const MeasurementVector &z) {
  const Filter::StateVector &x = filter_.state();

  double px = x(0);
  double py = x(1);
  double vx = x(2);
  double vy = x(3);
  double phi = atan2(py, px);
  double rho = sqrt(px*px + py*py);
  if (abs(rho) < EPSILON) {
    return;
  }
  double rho_dot = (px * vx + py * vy) / rho;
  Eigen::Vector3d h;
  h << rho, phi, rho_dot;
  MeasurementStateMatrix H = tools.CalculateJacobian(x);
  FusionEKF::Filter::Sensor<3>::UpdateEKF(z, h, H, R_);
}

FusionEKF::Laser::Laser(FusionEKF::Filter &filter)
  : FusionEKF::Filter::Sensor<2>(filter)
{
  // Measurement covariance matrix:
  R_ <<
    0.0225, 0,
    0, 0.0225;

  // Measurement matrix:
  H_ <<
    1, 0, 0, 0,
    0, 1, 0, 0;
}

void FusionEKF::Laser::Update(const MeasurementVector &z) {
  FusionEKF::Filter::Sensor<2>::Update(z, H_, R_);
}

/*
 * Constructor.
 */
FusionEKF::FusionEKF() :
  is_initialized_(false),
  previous_timestamp_(0),
  radar_(ekf_),
  laser_(ekf_),
  F_(Eigen::Matrix4d::Identity())  // transition matrix (dt terms added later)
{ }

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::Initialize(const MeasurementPackage &measurement_pack) {
  /**
    * Initialize the state ekf_.x_ with the first measurement.
    * Create the covariance matrix.
  */
  double rho, phi, px, py;
  switch (measurement_pack.sensor_type_) {
    case MeasurementPackage::RADAR:
      // Convert radar from polar to cartesian coordinates.
      rho = measurement_pack.raw_measurements_(0);
      phi = measurement_pack.raw_measurements_(1);
      px = rho * cos(phi);
      py = rho * sin(phi);
      break;
    case MeasurementPackage::LASER:
      px = measurement_pack.raw_measurements_(0);
      py = measurement_pack.raw_measurements_(1);
      break;
    default:
      cerr << "bad measurement pack sensor type " <<
        measurement_pack.sensor_type_ << endl;
      exit(EXIT_FAILURE);
  }

  Filter::StateVector x;
  x << px, py, 0, 0;

  Filter::StateMatrix P;
	P <<
    1, 0,    0,    0,
	  0, 1,    0,    0,
	  0, 0, 1000,    0,
	  0, 0,    0, 1000;

  ekf_.Initialize(x, P);
}

void FusionEKF::Predict(const MeasurementPackage &measurement_pack) {
  /**
   * Update the state transition matrix F according to the new elapsed time.
   */
  // Get elapsed time in seconds (from microseconds).
  double dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1e6;

  //1. Modify the F matrix so that the time is integrated
  F_(0, 2) = dt;
  F_(1, 3) = dt;

	//2. Set the process covariance matrix Q
	double dt2 = dt * dt;
	double dt3 = dt2 * dt / 2.0f;
	double dt4 = dt2 * dt2 / 4.0f;
	Q_ <<
	  dt4 * NOISE_AX,              0, dt3 * NOISE_AX,              0,
	               0, dt4 * NOISE_AY,              0, dt3 * NOISE_AY,
	  dt3 * NOISE_AX,              0, dt2 * NOISE_AX,              0,
	               0, dt3 * NOISE_AY,              0, dt2 * NOISE_AY;

  ekf_.Predict(F_, Q_);
}

void FusionEKF::Update(const MeasurementPackage &measurement_pack) {
  /**
   * Use the sensor type to perform the update step.
   * Update the state and covariance matrices.
   */
  switch(measurement_pack.sensor_type_) {
    case MeasurementPackage::RADAR:
      radar_.Update(measurement_pack.raw_measurements_);
    break;
    case MeasurementPackage::LASER:
      laser_.Update(measurement_pack.raw_measurements_);
      break;
    default:
      cerr << "bad measurement pack sensor type " <<
        measurement_pack.sensor_type_ << endl;
      exit(EXIT_FAILURE);
  }
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

  cout << "x_ = " << endl << ekf_.state() << endl;
  cout << "P_ = " << endl << ekf_.covariance() << endl;
}

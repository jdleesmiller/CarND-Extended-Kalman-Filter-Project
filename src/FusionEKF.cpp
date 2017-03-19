#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;

// Tolerance for avoiding division by zero.
const double EPSILON = 1e-6;

// We don't infer anything about speed when initializing, so the variance is
// large.
const double INITIAL_VELOCITY_VARIANCE = 25;

// Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
const double NOISE_AX = 9;
const double NOISE_AY = 9;

//
// Radar Sensor
//

FusionEKF::Radar::Radar(FusionEKF::Filter &filter) : Filter::Sensor<3>(filter) {
  // Measurement covariance matrix:
  R_ <<
    0.09, 0, 0,
    0, 0.0009, 0,
    0, 0, 0.09;

  // Approximate measurement matrix: will be reinitialized later, but zero it
  // for safety.
  Hj_ = MeasurementStateMatrix::Zero();
}

void FusionEKF::Radar::Initialize(const MeasurementVector &z) {
  double rho = z(0);
  double phi = z(1);
  double var_rho = R_(0, 0);
  double var_phi = R_(1, 1);

  double px = rho * cos(phi);
  double py = rho * sin(phi);

  // In addition to the position, we can also obtain an estimate of the variance
  // due to measurement error in the radar sensor. Let's take the variance of px
  // as an example; the variance for py is computed similarly.
  //
  // We need to propagate error through a product, rho * cos(phi), and a cosine.
  // Fortunately, Wikipedia tells us exactly how to do this:
  // https://en.wikipedia.org/wiki/Propagation_of_uncertainty#Example_formulas
  // (But note that all of these rules are approximate; they hold only when
  // the measurement error is small relative to the value being measured.)
  //
  // The radar R_ matrix is diagonal, which tells us that measurements of rho
  // and phi are uncorrelated. The formula for propagating error through the
  // product AB for normal random variables A and B therefore simplifies to:
  //
  // Var(AB) = Var(A) B^2 + Var(B) A^2
  //
  // The rule for cos(C) for normal random variable C is:
  //
  // Var(cos(C)) = Var(C) sin^2(C)
  //
  // Combining the above, we obtain
  //
  // Var(rho cos(phi))
  //   = Var(rho) cos^2(phi) + Var(cos(phi)) rho^2
  //   =                 ... + Var(phi) sin^2(phi) rho^2
  //   =                 ... + Var(phi) (py)^2
  //
  // Note: it turns out is possible to obtain an exact expression for
  // Var(cos(phi)): http://nbviewer.jupyter.org/gist/dougalsutherland/8513749
  // so we could try that if this approximation causes problems.
  double cos2phi = cos(phi) * cos(phi);
  double sin2phi = sin(phi) * sin(phi);
  double var_px = var_rho * cos2phi + var_phi * py * py;
  double var_py = var_rho * sin2phi + var_phi * px * px;

  Filter::StateVector x;
  x << px, py, 0, 0;

  Filter::StateMatrix P;
  double var_v = INITIAL_VELOCITY_VARIANCE;
  P <<
        var_px,      0,     0,     0,
             0, var_py,     0,     0,
             0,      0, var_v,     0,
             0,      0,     0, var_v;

  filter_.Initialize(x, P);
}

void FusionEKF::Radar::Update(const MeasurementVector &z) {
  const Filter::StateVector &x = filter_.state();

  double px = x(0);
  double py = x(1);
  double vx = x(2);
  double vy = x(3);

  // Project the current state into measurement space using the nonlinear 'h'.
  // (See course notes.)
  double phi = atan2(py, px);
  double rho = sqrt(px*px + py*py);
  // Check division by zero: ignore the measurement.
  if (abs(rho) < EPSILON) {
    return;
  }
  double rho_dot = (px * vx + py * vy) / rho;
  Eigen::Vector3d h;
  h << rho, phi, rho_dot;

  // Calculate the linearized projection Hj for calculating the Kalman gain and
  // updating the covariance. (See course notes.)
  double hp = px * px + py * py;
  double dp = sqrt(hp);
  double p32 = dp * hp;
  double c = vx * py - vy * px;
  // Check division by zero: ignore the measurement.
  if (abs(hp) < EPSILON) {
    return;
  }
  Hj_ <<
           px / dp,       py / dp,       0,       0,
          -py / hp,       px / hp,       0,       0,
      py * c / p32, -px * c / p32, px / dp, py / dp;

  FusionEKF::Filter::Sensor<3>::UpdateEKF(z, h, Hj_, R_);
}

//
// Laser Sensor
//

FusionEKF::Laser::Laser(FusionEKF::Filter &filter) : Filter::Sensor<2>(filter) {
  // Measurement covariance matrix:
  R_ <<
    0.0225, 0,
    0, 0.0225;

  // Measurement matrix:
  H_ <<
    1, 0, 0, 0,
    0, 1, 0, 0;
}

void FusionEKF::Laser::Initialize(const MeasurementVector &z) {
  Filter::StateVector x;
  x << z(0), z(1), 0, 0;

  Filter::StateMatrix P;
  double var_v = INITIAL_VELOCITY_VARIANCE;
  P <<
    R_(0),     0,     0,     0,
        0, R_(1),     0,     0,
        0,     0, var_v,     0,
        0,     0,     0, var_v;

  filter_.Initialize(x, P);
}

void FusionEKF::Laser::Update(const MeasurementVector &z) {
  Filter::Sensor<2>::Update(z, H_, R_);
}

//
// Main Fusion Filter
//

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
  switch (measurement_pack.sensor_type_) {
    case MeasurementPackage::RADAR:
      radar_.Initialize(measurement_pack.raw_measurements_);
      break;
    case MeasurementPackage::LASER:
      laser_.Initialize(measurement_pack.raw_measurements_);
      break;
    default:
      cerr << "bad measurement pack sensor type " <<
        measurement_pack.sensor_type_ << endl;
      exit(EXIT_FAILURE);
  }
}

void FusionEKF::Predict(const MeasurementPackage &measurement_pack) {
  // Get elapsed time in seconds (from microseconds).
  double dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1e6;

  // 1. Modify the F matrix so that the time is integrated
  F_(0, 2) = dt;
  F_(1, 3) = dt;

  // 2. Set the process covariance matrix Q. (See course notes.)
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

#include "kalman_filter.h"

#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  // MINE
  x_ = F_ * x_;
  P_ = F_ * P_ * F_.transpose() + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  // MINE
  UpdateEKF(z, H_ * x_);
}

void KalmanFilter::UpdateEKF(const VectorXd &z, const VectorXd &h) {
  // MINE
  std::cout << "z=" << z << std::endl;
  std::cout << "h=" << h << std::endl;
  VectorXd y = z - h;
  MatrixXd S = H_ * P_ * H_.transpose() + R_;
  MatrixXd K = P_ * H_.transpose() * S.inverse();
  MatrixXd I = MatrixXd::Identity(x_.size(), x_.size());

  x_ = x_ + K * y;
  P_ = (I - K * H_) * P_;
}

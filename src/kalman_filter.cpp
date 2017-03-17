#include "kalman_filter.h"
#include "Eigen/LU"

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
  VectorXd y = z - h;
  MatrixXd S = H_ * P_ * H_.transpose() + R_;
  MatrixXd I = MatrixXd::Identity(x_.size(), x_.size());

  // Avoid explicit matrix inversion. Starting from
  // K = P H^T S^-1
  // postmultiply by S to get
  // K S = P H^T
  // which is equivalent to
  // S^T K^T = H P^T
  // which is in standard form (Ax=B)
  Eigen::FullPivLU<MatrixXd> lu(S.transpose());
  MatrixXd K = lu.solve(H_ * P_.transpose()).transpose();

  x_ = x_ + K * y;
  P_ = (I - K * H_) * P_;
}

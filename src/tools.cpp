#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  // MINE
  VectorXd rmse(4);
  rmse << 0,0,0,0;

  // check the validity of the following inputs:
  //  * the estimation vector size should not be zero
  //  * the estimation vector size should equal ground truth vector size
  if (estimations.size() != ground_truth.size() || estimations.size() == 0) {
  	std::cerr << "Invalid estimation or ground_truth data" << std::endl;
  	return rmse;
  }

  //accumulate squared residuals
  for (size_t i=0; i < estimations.size(); ++i) {
  	VectorXd residual = estimations[i] - ground_truth[i];

  	//coefficient-wise multiplication
  	residual = residual.array()*residual.array();
  	rmse += residual;
  }

  //calculate the mean
  rmse = rmse / estimations.size();

  //calculate the squared root
  rmse = rmse.array().sqrt();

  //return the result
  return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
  // MINE
  MatrixXd Hj(3,4);

  //recover state parameters
	float px = x_state(0);
	float py = x_state(1);
	float vx = x_state(2);
	float vy = x_state(3);

  float hp = px * px + py * py;
  float dp = sqrt(hp);
  float p32 = dp * hp;
  float c = vx * py - vy * px;

  //check division by zero
  if (hp == 0) {
    return Hj;
  }

  //compute the Jacobian matrix
  Hj <<
           px / dp,       py / dp,       0,       0,
          -py / hp,       px / hp,       0,       0,
      py * c / p32, -px * c / p32, px / dp, py / dp;

  return Hj;
}

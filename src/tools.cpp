#include "tools.h"
#include <iostream>
#include <math.h>

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  VectorXd rmse(4);
  rmse << 0, 0, 0, 0;

  // check the validity of the following inputs
  if (estimations.size() != ground_truth.size() || estimations.size() == 0) {
    std::cout << "Invalid estimation or ground_truth data" << std::endl;
    return rmse;
  }

  // accumulate squared residuals
  for (unsigned int i = 0; i < estimations.size(); ++i) {

    VectorXd residual = estimations[i] - ground_truth[i];

    // coefficient-wise multiplication
    residual = residual.array() * residual.array();
    rmse += residual;
  }

  // calculate the mean
  rmse = rmse / estimations.size();
  std::cout << rmse << std::endl;

  // calculate the squared root
  rmse = rmse.array().sqrt();

  return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {

  // recover state parameters
  float px = x_state(0);
  float py = x_state(1);
  float vx = x_state(2);
  float vy = x_state(3);
  float px2py2 = px * px + py * py;
  if (px2py2 < 0.0001) {  // check division by zero
    std::cout << "CalculateJakobian () - Error - Division by Zero" << std::endl;
    return MatrixXd::Zero(3, 4);
  }

  MatrixXd Hj(3, 4);
  float sqrt_px2py2 = sqrt(px2py2);
  float px2py2_pow32 = px2py2 * sqrt_px2py2;
  Hj << px/sqrt_px2py2,                py/sqrt_px2py2,                0.,             0.,
        -py/px2py2,                    px/px2py2,                     0.,             0.,
        py*(vx*py-vy*px)/px2py2_pow32, px*(vy*px-vx*py)/px2py2_pow32, px/sqrt_px2py2, py/sqrt_px2py2;
  return Hj;
}

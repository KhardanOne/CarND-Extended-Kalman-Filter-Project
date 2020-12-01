#include "kalman_filter.h"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

/* 
 * Please note that the Eigen library does not initialize 
 *   VectorXd or MatrixXd objects with zeros upon creation.
 */

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(const VectorXd &x_in, const MatrixXd &P_in, const MatrixXd &F_in,
                        /*const MatrixXd &H_in,*/ const MatrixXd &R_in, const MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  // H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
  I_.setIdentity(P_.rows(), P_.cols());
}

void KalmanFilter::Predict() {
  x_ = F_ * x_;
  P_ = F_ * P_ * F_.transpose() + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
   * DONE: update the state by using Kalman Filter equations
   */
  VectorXd y = z - H_ * x_;
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd K = P_ * Ht * S.inverse();
  x_ += K * y;
  P_ = (I_ - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
   * TODO: update the state by using Extended Kalman Filter equations
   */
  VectorXd y = z - CartesianToPolar(x_);
  y[1] = atan2(sin(y[1]), cos(y[1]));  // crop phi between -pi and +pi
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd K = P_ * Ht * S.inverse();
  x_ += K * y;
  P_ = (I_ - K * H_) * P_;
}

Eigen::VectorXd KalmanFilter::CartesianToPolar(const Eigen::VectorXd &in) {
  double ppx = in(0);
  double ppy = in(1);
  double vpx = in(2);
  double vpy = in(3);
  double rho = sqrt(ppx * ppx + ppy * ppy);
  if (rho < 0.0001) {
    std::cout << "rho too low, setting to 0.0001 ====================================================================================================" << std::endl;
    rho = 0.0001;
  }
  double phi = atan2(ppy, ppx);
  double rhodot = (ppx * vpx + ppy * vpy) / rho;
  VectorXd res(3);
  res << rho, phi, rhodot;
  return res;
}

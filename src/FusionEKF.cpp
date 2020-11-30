#include "FusionEKF.h"
#include <iostream>
#include "Eigen/Dense"
#include "tools.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::cout;
using std::endl;
using std::vector;

/**
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
   * TODO: Finish initializing the FusionEKF.
   * TODO: Set the process and measurement noises
   */
  H_laser_ << 1, 0, 0, 0,
              0, 1, 0, 0;
  double noise_ax = 9.;
  double noise_ay = 9.;
  double ax2 = noise_ax * noise_ax;
  double ay2 = noise_ay * noise_ay;
  Q_base_ = MatrixXd(4, 4);
  Q_base_ <<  ax2/4, 0, ax2/2, 0,
              0, ay2/4, 0, ay2/2,
              ax2/2, 0, ax2, 0,
              0, ay2/2, 0, ay2;
}

/**
 * Destructor.
 */
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
  /**
   * Initialization
   */
  if (!is_initialized_) {
    /**
     * TODO: Initialize the state ekf_.x_ with the first measurement.
     * TODO: Create the covariance matrix.
     * You'll need to convert radar from polar to cartesian coordinates.
     */

    // first measurement
    cout << "EKF: " << endl;
    MatrixXd P(4, 4);
    P <<  1, 0, 0, 0,
          0, 1, 0, 0,
          0, 0, 1000, 0,
          0, 0, 0, 1000;
    MatrixXd F = MatrixXd::Identity(4, 4);

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      // TODO: Convert radar from polar to cartesian coordinates and initialize state.
      double rho = measurement_pack.raw_measurements_(int(DataIdxR::rho_measured));
      double phi = measurement_pack.raw_measurements_(int(DataIdxR::phi_measured));
      VectorXd x(4);
      x(0) = rho*cos(phi);
      x(1) = -rho*sin(phi);
      x(2) = 0.;
      x(3) = 0.;
      ekf_.Init(x, P, F, R_radar_, Q_base_);
      previous_timestamp_ = measurement_pack.raw_measurements_(int(DataIdxR::timestamp));
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      // TODO: Initialize state.
      VectorXd x(4);
      x(0) = measurement_pack.raw_measurements_(int(DataIdxL::x_measured));
      x(1) = measurement_pack.raw_measurements_(int(DataIdxL::y_measured));
      x(2) = 0.;
      x(3) = 0.;
      ekf_.Init(x, P, F, R_radar_, Q_base_);
      previous_timestamp_ = measurement_pack.raw_measurements_(int(DataIdxL::timestamp));
    }

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /**
   * Prediction
   */

  /**
   * DONE: Update the state transition matrix F according to the new elapsed time.
   * Time is measured in seconds.
   */ 
  double dt = measurement_pack.timestamp_ - previous_timestamp_;
  double dt2 = dt * dt;
  double dt3 = dt * dt2;
  double dt4 = dt * dt3;
  previous_timestamp_ = measurement_pack.timestamp_;
  ekf_.F_(0, 2) = dt;
  ekf_.F_(1, 3) = dt;

  /**
   * DONE: Update the process noise covariance matrix.
   * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */
  ekf_.Q_(0, 0) = Q_base_(0, 0) * dt4;
  ekf_.Q_(0, 2) = Q_base_(0, 2) * dt3;
  ekf_.Q_(1, 1) = Q_base_(1, 1) * dt4;
  ekf_.Q_(1, 3) = Q_base_(1, 3) * dt3;
  ekf_.Q_(2, 0) = Q_base_(2, 0) * dt3;
  ekf_.Q_(2, 2) = Q_base_(2, 2) * dt2;
  ekf_.Q_(3, 1) = Q_base_(3, 1) * dt3;
  ekf_.Q_(3, 3) = Q_base_(3, 3) * dt2;

  ekf_.Predict();

  /**
   * Update
   */

  /**
   * TODO:
   * - Use the sensor type to perform the update step.
   * - Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // TODO: Radar updates
    double rho = measurement_pack.raw_measurements_(int(DataIdxR::rho_measured));
    double phi = measurement_pack.raw_measurements_(int(DataIdxR::phi_measured));
    double rhodot = measurement_pack.raw_measurements_(int(DataIdxR::rhodot_measured));
    VectorXd z(3);
    z << rho, phi, rhodot;
    ekf_.UpdateEKF(z);
  } else {
    // TODO: Laser updates
    double x = measurement_pack.raw_measurements_(int(DataIdxL::x_measured));
    double y = measurement_pack.raw_measurements_(int(DataIdxL::y_measured));
    VectorXd z(4);
    z << x, y, 0, 0;
    ekf_.Update(z);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}

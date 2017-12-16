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

  // Measurement covariance matrix R - laser
  R_laser_ << 0.0225, 0,
              0, 0.0225;

  // Measurement covariance matrix R - radar
  R_radar_ << 0.09, 0, 0,
              0, 0.0009, 0,
              0, 0, 0.09;

  // Measurement function H - radar
  H_laser_ << 1, 0, 0, 0,
              0, 1, 0, 0;
  
  // State covariance matrix P 
  ekf_.P_ = MatrixXd(4, 4);
  ekf_.P_ <<  1, 0, 0,    0,
              0, 1, 0,    0,
              0, 0, 1000, 0,
              0, 0, 0,    1000;
  
  // Transition matrix F, 
  // the values at (1, 3) and (2, 4) will correspond to dt after first update
  ekf_.F_ = MatrixXd(4, 4);
  ekf_.F_ <<  1, 0, 1, 0,
              0, 1, 0, 1,
              0, 0, 1, 0,
              0, 0, 0, 1;

  // Set noise components from task as per project spec
  noise_ax = 9;
  noise_ay = 9;
}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    // State vector X => x, y, vx, vy
    ekf_.x_ = VectorXd(4);

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      float rho = measurement_pack.raw_measurements_[0];
      float phi = measurement_pack.raw_measurements_[1];
      float drho_dt = measurement_pack.raw_measurements_[2];

      ekf_.x_ <<  cos(phi) * rho, 
                  sin(phi) * rho, 
                  drho_dt * cos(phi), 
                  drho_dt * sin(phi);
    }
    else  if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      ekf_.x_ <<  measurement_pack.raw_measurements_[0],
                  measurement_pack.raw_measurements_[1],
                  0,
                  0;
    }
    if (ekf_.x_(0) < 0.0001) {
      ekf_.x_(0) = 0.0001;
    }

    previous_timestamp_ = measurement_pack.timestamp_;

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  // Elapsed time between measurements dt in seconds
  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
  cout << dt << endl;
  previous_timestamp_ = measurement_pack.timestamp_;

  // dt at different powers
  float dt_2 = dt * dt;
  float dt_3 = dt_2 * dt;
  float dt_4 = dt_3 * dt;
  float dt_4_4 = dt_4 / 4;
  float dt_3_2 = dt_3 / 2;


  // Update state transition matrix F with dt
  ekf_.F_(0, 2) = dt;
  ekf_.F_(1, 3) = dt;

  // Compute process covariance matrix Q
  ekf_.Q_ = MatrixXd(4, 4);
  ekf_.Q_ <<  dt_4_4 * noise_ax,    0,                  dt_3_2 * noise_ax,   0,
              0,                    dt_4_4 * noise_ay,  0,                   dt_3_2 * noise_ay,
              dt_3_2 * noise_ax,    0,                  dt_2 * noise_ax,     0,
              0,                    dt_3_2 * noise_ay,  0,                   dt_2 * noise_ay;
  
  // Run prediction
  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  // Laser
  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    ekf_.H_ = tools.CalculateJacobian(ekf_.x_);
    ekf_.R_ = R_radar_;
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);
  }
  // Radar
  else {
    ekf_.H_ = H_laser_;
    ekf_.R_ = R_laser_;
    ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}

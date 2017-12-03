#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;
#define SMALL_VALUE 0.0001

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
}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

Eigen::VectorXd FusionEKF::PolarToCartesian(const long rho, const long phi, const long drho_dt)
{
  Eigen::VectorXd x = Eigen::VectorXd(4);
  x <<  rho * cos(phi),
        rho * sin(phi),
        drho_dt * cos(phi),
        drho_dt * sin(phi);
  return x;
}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    // State vector X => x, y, vx, vy
    ekf_.x_ = VectorXd(4);

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      long rho = measurement_pack.raw_measurements_[0];
      long phi = measurement_pack.raw_measurements_[1];
      long drho_dt = measurement_pack.raw_measurements_[2];
      ekf_.x_ = PolarToCartesian(rho, phi, drho_dt);
    }
    else  {
      // LIDAR
      ekf_.x_ <<  measurement_pack.raw_measurements_[0],
                  measurement_pack.raw_measurements_[1],
                  0,
                  0;
    }
    if (fabs(ekf_.x_(0)) < SMALL_VALUE && fabs(ekf_.x_(1)) < SMALL_VALUE)
    {
      ekf_.x_(0) = SMALL_VALUE;
      ekf_.x_(1) = SMALL_VALUE;
    }
    // State covariance matrix P 
    ekf_.P_ = MatrixXd(4, 4);
    ekf_.P_ <<  1, 0, 0,    0,
                0, 1, 0,    0,
                0, 0, 1000, 0,
                0, 0, 0,    1000;
    
    // Transition matrix F, the values at (1, 3) and (2, 4) will correspond to dt after first update
    ekf_.F_ = MatrixXd(4, 4);
    ekf_.F_ <<  1, 0, 1, 0,
                0, 1, 0, 1,
                0, 0, 1, 0,
                0, 0, 0, 1;

    previous_timestamp_ = measurement_pack.timestamp_;

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  // Elapsed time between measurements dt in seconds
  long long dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000;
  previous_timestamp_ = measurement_pack.timestamp_;

  // dt at different powers
  long long dt_2 = dt * dt;
  long long dt_3 = dt_2 * dt;
  long long dt_4 = dt_3 * dt;
  long long dt_4_4 = dt_4 / 4;
  long long dt_3_2 = dt_3 / 2;

  // Set noise components from task, values specified by project spec
  float noise_ax = 9.0;
  float noise_ay = 9.0;

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

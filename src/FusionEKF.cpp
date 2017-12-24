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
  ekf_ = KalmanFilter();

  // Measurement covariance matrix R - laser
  R_laser_ = MatrixXd(2, 2);
  R_laser_ << 0.0225, 0,
              0, 0.0225;

  // Measurement covariance matrix R - radar
  R_radar_ = MatrixXd(3, 3);
  R_radar_ << 0.09, 0, 0,
              0, 0.0009, 0,
              0, 0, 0.09;

  // Measurement function H - radar
  H_laser_ = MatrixXd(2, 4);
  H_laser_ << 1, 0, 0, 0,
              0, 1, 0, 0;
  
  // State covariance matrix P 
  ekf_.P_ = MatrixXd(4, 4);
  ekf_.P_ <<  1, 0, 0,    0,
              0, 1, 0,    0,
              0, 0, 1000, 0,
              0, 0, 0,    1000;
  
  // Transition matrix F, 
  // the values at (0, 2) and (1, 3) will correspond to dt after first update
  ekf_.F_ = MatrixXd(4, 4);
  ekf_.F_ <<  1, 0, 1, 0,
              0, 1, 0, 1,
              0, 0, 1, 0,
              0, 0, 0, 1;

  // Process matrix Q
  ekf_.Q_ = MatrixXd(4, 4);
  ekf_.Q_ <<  0, 0, 0, 0,
              0, 0, 0, 0,
              0, 0, 0, 0,
              0, 0, 0, 0;

  // State vector X => x, y, vx, vy
  ekf_.x_ = VectorXd(4);
  ekf_.x_ << 1, 1, 1, 1;

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
    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      double rho = measurement_pack.raw_measurements_(0);
      double phi = measurement_pack.raw_measurements_(1);
      double drho_dt = measurement_pack.raw_measurements_(2);

      ekf_.x_ <<  cos(phi) * rho, 
                  sin(phi) * rho, 
                  0,
                  0;
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      ekf_.x_ <<  measurement_pack.raw_measurements_(0),
                  measurement_pack.raw_measurements_(1),
                  0,
                  0;
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
  previous_timestamp_ = measurement_pack.timestamp_;

  // Update state transition matrix F with dt
  ekf_.F_(0, 2) = dt;
  ekf_.F_(1, 3) = dt;
  
  float dt_2 = pow(dt, 2);
  float dt_3 = pow(dt, 3);
  float dt_4 = pow(dt, 4);
  float dt_4_4 = dt_4 / 4;
  float dt_3_2 = dt_3 / 2;
  
  // Compute process covariance matrix Q
  ekf_.Q_ << dt_4_4 * noise_ax, 0,                    dt_3_2 * noise_ax,  0,
             0,                   dt_4_4 * noise_ay,  0,                   dt_3_2 * noise_ay,
             dt_3_2 * noise_ax, 0,                    dt_2 * noise_ax,    0,
	           0,                   dt_3_2 * noise_ay,  0,                   dt_2 * noise_ay;
  
  // Run prediction
  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  // Radar
  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    try {
      ekf_.H_ = tools.CalculateJacobian(ekf_.x_);
      ekf_.R_ = R_radar_;
      ekf_.UpdateEKF(measurement_pack.raw_measurements_);
    } catch (std::invalid_argument e) {
      cout << "Error during radar update: " << e.what() << '\n';
    }
  }
  // Laser
  else {
    ekf_.H_ = H_laser_;
    ekf_.R_ = R_laser_;
    ekf_.Update(measurement_pack.raw_measurements_);
  }
}

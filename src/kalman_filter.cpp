#include "kalman_filter.h"
#include <iostream>
#include <math.h>

using namespace std;
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
  x_ = F_ * x_;
  MatrixXd F_transpose = F_.transpose();
  P_ = F_ * P_ * F_transpose + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  VectorXd y = z - H_ * x_;
  MatrixXd I = MatrixXd::Identity(x_.size(), x_.size());
  MatrixXd H_transpose = H_.transpose();
  MatrixXd S = H_ * P_ * H_transpose + R_;
  MatrixXd K = P_ * H_transpose * S.inverse();
  x_ += K * y;
  // update
  P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  float px, py, vx, vy, rho, phi, drho_dt;
  px = x_[0];
  py = x_[1];
  vx = x_[2];
  vy = x_[3];
  rho = sqrt(px * px + py*py);
  phi = atan2(py, px);
  if (rho < 0.0001)
  {
    rho = 0.0001;
  } 
  drho_dt = (px * vx + py * vy) / rho;
  VectorXd prediction(3);
  prediction << rho, phi, drho_dt;
  VectorXd y = z - prediction;
  while (y(1) > M_PI) {
    y(1) -= 2 * M_PI;
  }
  while (y(1) < -M_PI) {
    y(1) += 2 * M_PI;
  }
  MatrixXd I = MatrixXd::Identity(x_.size(), x_.size());
  MatrixXd H_transpose = H_.transpose();
  MatrixXd S = H_ * P_ * H_transpose + R_;
  MatrixXd K = P_ * H_transpose * S.inverse();
  x_ += K * y;
  // update
  P_ = (I - K * H_) * P_;
}
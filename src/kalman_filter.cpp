#include "kalman_filter.h"

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
  MatrixXd H_transpose = H_.transpose();
  MatrixXd S = H_ * P_ * H_transpose + R_;
  MatrixXd K = P_ * H_transpose * S.inverse();
  MatrixXd I = MatrixXd::Identity(x_.size(), x_.size());
  x_ += K * y;

  // Update
  P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  float rho = sqrt(x_(0) * x_(0) + x_(1) * x_(1));
  float phi = atan2(x_(1), x_(0));
  float drho_dt;
  if (fabs(rho) < 0.0000001)
  {
    drho_dt = 0;
  } else
  {
    drho_dt = (x_(0) * x_(2) + x_(1) * x_(3)) / rho;
  }
  VectorXd prediction(3);
  prediction << rho, phi, drho_dt;
  VectorXd y = z - prediction;
  MatrixXd H_transpose = H_.transpose();
  MatrixXd S = H_ * P_ * H_transpose + R_;
  MatrixXd K = P_ * H_transpose * S.inverse();
  MatrixXd I = MatrixXd::Identity(x_.size(), x_.size());
  x_ += K * y;

  // Update
  P_ = (I - K * H_) * P_;
}
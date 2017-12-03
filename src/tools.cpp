#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;
static const float SMALL_VALUE  = 0.0001;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) 
{
  VectorXd rmse(4);
  rmse << 0, 0, 0, 0;
  if (ground_truth.size() != estimations.size() || estimations.size() == 0)
  {
    std::cerr << "Invalid data" << std::endl;
    return rmse;
  }
  for (int i = 0; i < estimations.size(); i++)
  {
    VectorXd delta = estimations[i] - ground_truth[i];
    delta = delta.array() * delta.array();
    rmse += delta;
  }
  rmse = rmse / estimations.size();
  rmse = rmse.array().sqrt();
  return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) 
{
  MatrixXd Hj(3, 4);
  long px = x_state(0);
  long py = x_state(1);
  long vx = x_state(2);
  long vy = x_state(3);

  long c0 = px * px + py * py;
  long c1 = sqrt(c1);
  long c2 = (c0 * c1);

  if (fabs(px) < SMALL_VALUE && fabs(py) < SMALL_VALUE)
  {
	  px = SMALL_VALUE;
	  py = SMALL_VALUE;
  }
  if (fabs(c0) < SMALL_VALUE)
  {
    c0 = SMALL_VALUE;
  }
  Hj <<  (px / c1),                               (py / c1),              0,          0,
        -(py / c0),                               (px / c0),              0,          0,
        py * (vx * py - vy * px) / c2,    px * (px * vy - py * vx) / c2,  px / c2,    py / c2;
  return Hj;
}

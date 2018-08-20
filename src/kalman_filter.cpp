#include "kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

#define PI 3.14159265

// Please note that the Eigen library does not initialize 
// VectorXd or MatrixXd objects with zeros upon creation.

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
  // No external control function so consider control vector = 0;
  // In prediction step, we are using linear model, so we still use regular KF equation
  // rather than EKF (used if non-linear model is used)
  // x = F*x + u ; u = 0
  x_ = F_ * x_ ; 
  MatrixXd Ft = F_.transpose();
  P_ = F_ * P_ * Ft + Q_; // calculate Process  Covariance matrix
}

void KalmanFilter::Update(const VectorXd &z) {
  // Calculate Error function for laser
  VectorXd y = z - H_ * x_; // error calculation
  
  // Update State measurement
  // Calculate Update Matrices
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;

  //new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
}



void KalmanFilter::UpdateEKF(const VectorXd &z) {
  // Calculate Error function for Radar (y = z = hbar(x))
  VectorXd h = VectorXd(3); // h(x_)

  double rho = sqrt(x_(0)*x_(0) + x_(1)*x_(1));
  double phi = atan2(x_(1), x_(0));
  double rhd = (x_(0)*x_(2) + x_(1)*x_(3)) / rho;
  
  h << rho, phi, rhd;
  
  VectorXd y = z - h;
  
  // Normalize the phi value to be in the valid range [-PI, PI]
  if (y(1) > PI)  { y(1) = 2.0*PI - y(1); }
  else if (y(1) < (-1.0 * PI)) {y(1) = 2.0*PI + y(1); }
  
  // Update State measurement
  // Calculate Update Matrices
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;

  //new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
}


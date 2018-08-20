#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

#define EPS 0.0001

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
  Hj_      = MatrixXd(3, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
              0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
              0, 0.0009, 0,
              0, 0, 0.09;

  /**
  TODO:
    * Finish initializing the FusionEKF.
    * Set the process and measurement noises
  */
  ekf_.x_ = VectorXd(4);
  ekf_.P_ = MatrixXd(4, 4);
  ekf_.F_ = MatrixXd(4, 4);
  //state covariance matrix P
  ekf_.P_ << 1, 0, 0, 0,
		     0, 1, 0, 0,
		     0, 0, 1000, 0,
		     0, 0, 0, 1000;


  //measurement matrix - laser
  H_laser_ << 1, 0, 0, 0,
			  0, 1, 0, 0;

  //the initial transition matrix F_
  ekf_.F_ << 1, 0, 1, 0,
	         0, 1, 0, 1,
		     0, 0, 1, 0,
		     0, 0, 0, 1;

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
    /**
    TODO:
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
    // first measurement
    cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
	  // Read in the Radar measurements rho = Radial distance, 
	  // phi = Angular Direction, rho_dot = Radial velocity
	  double rho = measurement_pack.raw_measurements_[0]; 
	  double phi = measurement_pack.raw_measurements_[1]; 
	  double rhd = measurement_pack.raw_measurements_[2]; 
	  // Coordinates convertion from polar to cartesian
	  double x = rho * cos(phi); 
	  double y = rho * sin(phi);
	  if (x < EPS) { x = EPS; }
	  if (y < EPS) { y = EPS; }
	  double vx = rhd * cos(phi);
	  double vy = rhd * sin(phi);
	  ekf_.x_ << x, y, vx , vy;  
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
	  //set the state with the initial location and zero velocity
	  ekf_.x_ << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0;
    }
	
	
    // done initializing, no need to predict or update
	previous_timestamp_ = measurement_pack.timestamp_;
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
   TODO:
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */
   //compute the time elapsed between the current and previous measurements
	double dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;	//dt - expressed in seconds
	previous_timestamp_ = measurement_pack.timestamp_;
	
	//1. Modify the F matrix so that the time is integrated
	ekf_.F_ << 1, 0, dt, 0,
			   0, 1, 0, dt,
			   0, 0, 1, 0,
			   0, 0, 0, 1;
 
	//2. Set the process covariance matrix Q
	double noise_ax = 9;
	double noise_ay = 9;
	double dt4 = std::pow(dt,4)/4;
	double dt3 = std::pow(dt,3)/2;
	double dt2 = std::pow(dt,2)/1;
	ekf_.Q_ = MatrixXd(4, 4);
	ekf_.Q_ << dt4*noise_ax, 0, dt3*noise_ax, 0,
			   0, dt4*noise_ay, 0, dt3*noise_ay,
			   dt3*noise_ax, 0, dt2*noise_ax, 0,
			   0, dt3*noise_ay, 0, dt2*noise_ay;
	
    ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/
   // Call the Kalman Filter update() function with the most recent raw measurements_
    VectorXd z;
    z = measurement_pack.raw_measurements_;
  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
	// For measurement update function use Jacobian matrix
	ekf_.H_ = tools.CalculateJacobian(ekf_.x_);
	ekf_.R_ = R_radar_;
	ekf_.UpdateEKF(z);
  } else {
    // Laser measurement update
	ekf_.H_ = H_laser_;
	ekf_.R_ = R_laser_;
	ekf_.Update(z);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}

#include "ukf.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

Tools tools;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {

  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);
  x_ << 0, 0, 0, 0, 0;

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 30;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 30;

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;

  // initial covariance matrix
  P_ = MatrixXd(5, 5);
  P_ << 1, 0, 0, 0, 0,
        0, 1, 0, 0, 0,
        0, 0, 1, 0, 0,
        0, 0, 0, 1, 0,
        0, 0, 0, 0, 1; // TODO std...

    // NEU

  //set state dimension
  n_x = 5;

  //set augmented dimension
  n_aug = 7;

  //set measurement dimension
  n_z_radar = 3;
  n_z_laser = 2;

  //define spreading parameter
  lambda = 3 - n_aug;

  /**
  Complete the initialization. See ukf.h for other member properties.
  Hint: one or more values initialized above might be wildly off...
  */
  is_initialized_ = false;
  previous_timestamp_ = 0;

//  long long time_us_;

  // set weights
  weights = VectorXd(2 * n_aug + 1);
  weights.setConstant( 0.5 / (lambda + n_aug) );
  weights[0] =  lambda / (lambda + n_aug);

  // Sensor matrices
  R_laser_ = MatrixXd(2, 2);
  R_laser_ << std_laspx_ * std_laspx_, 0,
              0, std_laspy_ * std_laspy_;

  R_radar_ = MatrixXd(3, 3);
  R_radar_ << std_radr_ * std_radr_, 0, 0,
              0, std_radphi_ * std_radphi_, 0,
              0, 0, std_radrd_ * std_radrd_;
}

UKF::~UKF() {}


/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
// void UKF::ProcessMeasurement(MeasurementPackage meas_package) { }
void UKF::ProcessMeasurement(MeasurementPackage measurement_pack) {

  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    /**
      * Initialize the state x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
    // first measurement
    x_ << 1, 1, 1, 1, 1;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      // RADAR: meas_package.raw_measurements_ << ro, phi, ro_dot;
        float rho     = measurement_pack.raw_measurements_[0];
        float phi     = measurement_pack.raw_measurements_[1];
        float rho_dot = measurement_pack.raw_measurements_[2];

        float x = rho * cos(phi);  
        float y = rho * sin(phi);

        // state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
        x_ << x, y, rho_dot, phi, 0;
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      set the state with the initial location and zero velocity
      */
     float rho = measurement_pack.raw_measurements_[0]; 
     float phi = measurement_pack.raw_measurements_[1]; 

      // state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
      x_ << rho, phi, 0, 0, 0;
    }

    // remember previous timestamp
    previous_timestamp_ = measurement_pack.timestamp_;

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction (identical for laser and radar)
   ****************************************************************************/

  /**
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */

    //compute the time elapsed between the current and previous measurements
    //dt - expressed in seconds
    float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0; 
    previous_timestamp_ = measurement_pack.timestamp_;

    //predict
    Prediction(dt);

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */
  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
    // measurement update
    UpdateRadar(measurement_pack);

    /*VectorXd z_pred = VectorXd(4);
    Hj_ = tools.CalculateJacobian(ekf_.x_);
    z_pred = tools.fromCartesianToPolar(ekf_.x_);
    SetMeasurementMatrix(Hj_);
    SetMeasurementNoise(R_radar_);
    UpdateLidar(measurement_pack.raw_measurements_, z_pred);*/
  } else {
    // Laser updates
    UpdateLidar(measurement_pack);
    
    // measurement update
    /*SetMeasurementMatrix(H_laser_);
    SetMeasurementNoise(R_laser_);*/
  }

  // print the output
  cout << "x_ = " << x_ << endl;
  cout << "P_ = " << P_ << endl;
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */
  /*****************************************************************************
   * predict the state
   * Sigma points prediction
   ****************************************************************************/
  // generate sigma points
  MatrixXd Xsig; // Augmented sigma points
//  MatrixXd Xsig_pred; // predicted sigma points

  AugmentedSigmaPoints(&Xsig, x_, P_);
  SigmaPointPrediction(&Xsig_pred_, Xsig);
  PredictMeanAndCovariance(&x_, &P_, Xsig_pred_);
  // used to be:
  //    x_ = F_ * x_;
  //    P_ = F_ * P_ * F_.transpose() + Q_;
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage measurement_pack) {
  /**
  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
  float px = measurement_pack.raw_measurements_[0]; 
  float py = measurement_pack.raw_measurements_[1]; 

  VectorXd z = VectorXd(n_z_laser);
  z << py, py; 

  VectorXd z_pred;
  MatrixXd S;
  MatrixXd Zsig;

  PredictLaserMeasurement(&z_pred, &S, Xsig_pred_);

  //new estimate
  UpdateState(&x_, &P_, Xsig_pred_, x_, P_, Zsig, z_pred, S, z);

  // TODO: NIS
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage measurement_pack) {
  /**
  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
  float rho  = measurement_pack.raw_measurements_[0]; 
  float phi  = measurement_pack.raw_measurements_[1]; 
  float rhod = measurement_pack.raw_measurements_[2]; 

  VectorXd z = VectorXd(n_z_radar);
  z << rho, phi, rhod;

  VectorXd z_pred;
  MatrixXd S;
  MatrixXd Zsig;
  PredictRadarMeasurement(&z_pred, &S, Xsig_pred_);

  //new estimate
  UpdateState(&x_, &P_, Xsig_pred_, x_, P_, Zsig, z_pred, S, z);
}

/**
 * TODO
 * @param TODO
 */
void UKF::AugmentedSigmaPoints(MatrixXd* Xsig_out, const VectorXd& x, const MatrixXd& P) {

 /* //set state dimension
  int n_x = 5; // TODO class

  //set augmented dimension
  int n_aug = 7;

  //Process noise standard deviation longitudinal acceleration in m/s^2
  double std_a = 0.2;

  //Process noise standard deviation yaw acceleration in rad/s^2
  double std_yawdd = 0.2;

  //define spreading parameter
  double lambda = 3 - n_aug;
*/

  //x: set example state
  //P: create example covariance matrix

  //create augmented mean vector
  VectorXd x_aug = VectorXd(n_aug);

  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(n_aug, n_aug);

  //create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug, 2 * n_aug + 1);
 
  //create augmented mean state
  x_aug.fill(0.0);
  x_aug.head(5) = x;
  
  //create augmented covariance matrix
  P_aug.fill(0.0);
  P_aug.topLeftCorner(5, 5) = P;
  P_aug(5,5) = std_a_*std_a_;
  P_aug(6,6) = std_yawdd_*std_yawdd_;
  // std::cout << "P_aug = " << std::endl << P_aug << std::endl;  
  
  //create square root matrix
  MatrixXd A_aug = P_aug.llt().matrixL();
  
  //create augmented sigma points
  
  //set first column of sigma point matrix
  Xsig_aug.col(0)  = x_aug;

  //set remaining sigma points
  for (int i = 0; i < n_aug; i++)
  {
    Xsig_aug.col(i+1)       = x_aug + sqrt(lambda+n_aug) * A_aug.col(i);
    Xsig_aug.col(i+1+n_aug) = x_aug - sqrt(lambda+n_aug) * A_aug.col(i);
  }  
  
  //print result
  std::cout << "Xsig_aug = " << std::endl << Xsig_aug << std::endl;

  //write result
  *Xsig_out = Xsig_aug;
}

void UKF::SigmaPointPrediction(MatrixXd* Xsig_out, const MatrixXd& Xsig_aug) {

 /* //set state dimension
  int n_x = 5;

  //set augmented dimension
  int n_aug = 7;
*/
  // MatrixXd Xsig_aug = MatrixXd(n_aug, 2 * n_aug + 1);
  //create example sigma point matrix

  //create matrix with predicted sigma points as columns
  MatrixXd Xsig_pred = MatrixXd(n_x, 2 * n_aug + 1);

  double delta_t = 0.1; //time diff in sec

  for(int i = 0; i < 2 * n_aug + 1; ++i) {
    //extract values for better readability
    double p_x      = Xsig_aug(0,i);
    double p_y      = Xsig_aug(1,i);
    double v        = Xsig_aug(2,i);
    double yaw      = Xsig_aug(3,i);
    double yawd     = Xsig_aug(4,i);
    double nu_a     = Xsig_aug(5,i);
    double nu_yawdd = Xsig_aug(6,i);

    //predicted state values
    double px_p, py_p, v_p, yaw_p, yawd_p;

    //avoid division by zero
    if (fabs(yawd) < 0.001) {
        px_p = p_x + v*delta_t*cos(yaw);
        py_p = p_y + v*delta_t*sin(yaw);
    }
    else {
        px_p = p_x + v/yawd * ( sin (yaw + yawd*delta_t) - sin(yaw));
        py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) );
    }

    v_p = v;
    yaw_p = yaw + yawd*delta_t;
    yawd_p = yawd;

    //add noise
    px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
    py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
    v_p = v_p + nu_a*delta_t;

    yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
    yawd_p = yawd_p + nu_yawdd*delta_t;

    //write predicted sigma point into right column
    Xsig_pred.col(i) << px_p, py_p, v_p, yaw_p, yawd_p;
  }
  
  //print result
  std::cout << "Xsig_pred = " << std::endl << Xsig_pred << std::endl;

  //write result
  *Xsig_out = Xsig_pred;
}


void UKF::PredictMeanAndCovariance(VectorXd* x_out, MatrixXd* P_out, const MatrixXd& Xsig_pred) {

  /*//set state dimension
  int n_x = 5;

  //set augmented dimension
  int n_aug = 7;

  //define spreading parameter
  double lambda = 3 - n_aug;
*/
  // MatrixXd Xsig_pred = MatrixXd(n_x, 2 * n_aug + 1);
  //create example matrix with predicted sigma points

  //create vector for weights
  VectorXd weights = VectorXd(2*n_aug+1);
  
  //create vector for predicted state
  VectorXd x = VectorXd(n_x);

  //create covariance matrix for prediction
  MatrixXd P = MatrixXd(n_x, n_x);

  //set weights
  weights.setConstant( 0.5 / (lambda + n_aug) );
  weights[0] =  lambda / (lambda + n_aug);
  std::cout << "weights = " << std::endl << weights << std::endl;  
  
  //predict state mean
  x.setZero();
  for(int i = 0; i < 2*n_aug + 1; ++i) {
    x += weights[i] * Xsig_pred.col(i);
  }
  
  //predict state covariance matrix
  P.setZero();
  for(int i = 0; i < 2*n_aug + 1; ++i) {
      VectorXd x_diff = Xsig_pred.col(i) - x;

      //angle normalization
      while (x_diff(3) >  M_PI) x_diff(3) -= 2.0*M_PI;
      while (x_diff(3) < -M_PI) x_diff(3) += 2.0*M_PI;
      
      P += weights[i] * x_diff * x_diff.transpose();
  }

  //print result
  std::cout << "Predicted state" << std::endl;
  std::cout << x << std::endl;
  std::cout << "Predicted covariance matrix" << std::endl;
  std::cout << P << std::endl;

  //write result
  *x_out = x;
  *P_out = P;
}

void UKF::PredictLaserMeasurement(VectorXd* z_out, MatrixXd* S_out, const MatrixXd& Xsig_pred) {
}

void UKF::PredictRadarMeasurement(VectorXd* z_out, MatrixXd* S_out, const MatrixXd& Xsig_pred) {
// TODO return Zsig

 /* //set state dimension
  int n_x = 5;

  //set augmented dimension
  int n_aug = 7;

  //set measurement dimension, radar can measure r, phi, and r_dot
  int n_z = 3;

  //define spreading parameter
  double lambda = 3 - n_aug;

  //set vector for weights
  VectorXd weights = VectorXd(2*n_aug+1);
  double weight_0 = lambda/(lambda+n_aug);
  weights(0) = weight_0;
  for (int i=1; i<2*n_aug+1; i++) {  
    double weight = 0.5/(n_aug+lambda);
    weights(i) = weight;
  } // TODO init
*/

  //radar measurement noise standard deviation radius in m
  double std_radr = 0.3;

  //radar measurement noise standard deviation angle in rad
  double std_radphi = 0.0175;

  //radar measurement noise standard deviation radius change in m/s
  double std_radrd = 0.1;

  //create example matrix with predicted sigma points
  // MatrixXd Xsig_pred = MatrixXd(n_x, 2 * n_aug + 1);

  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z_radar, 2 * n_aug + 1);

  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z_radar);
  
  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z_radar, n_z_radar);
 
  //transform sigma points into measurement space
  for(int i = 0; i < 2 * n_aug + 1; i++) {
      double p_x = Xsig_pred(0,i);
      double p_y = Xsig_pred(1,i);
      double v   = Xsig_pred(2,i);
      double yaw = Xsig_pred(3,i);
      
      double rho = sqrt(p_x*p_x + p_y*p_y);
      double phi = atan2(p_y,p_x);
      double rhod = (p_x * cos(yaw) * v + p_y * sin(yaw) * v)/rho;
      
      Zsig.col(i) << rho, phi, rhod;
  } 
  
  //calculate mean predicted measurement
  z_pred.setZero();                           
  for(int i = 0; i < 2*n_aug + 1; i++) {
    z_pred += weights[i] * Zsig.col(i);
  }
 
  //calculate measurement covariance matrix S
  S.setZero();
  for(int i = 0; i < 2*n_aug + 1; i++) {
    VectorXd z_diff = Zsig.col(i) - z_pred;
    
    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    S += weights[i] * z_diff * z_diff.transpose();
  }
  
  MatrixXd R = MatrixXd(n_z_radar, n_z_radar);
  R.setZero();
  R(0,0) = std_radr*std_radr;
  R(1,1) = std_radphi*std_radphi;
  R(2,2) = std_radrd*std_radrd;
  
  S = S + R;

  //print result
  std::cout << "z_pred: " << std::endl << z_pred << std::endl;
  std::cout << "S: " << std::endl << S << std::endl;

  //write result
  *z_out = z_pred;
  *S_out = S;
}

void UKF::UpdateState(VectorXd* x_out, MatrixXd* P_out, const MatrixXd& Xsig_pred, const VectorXd& x, const MatrixXd& P,
    const MatrixXd& Zsig, const VectorXd& z_pred, const MatrixXd& S, const VectorXd& z) {

  /*//set state dimension
  int n_x = 5;

  //set augmented dimension
  int n_aug = 7;

  //set measurement dimension, radar can measure r, phi, and r_dot
  int n_z = 3;

  //define spreading parameter
  double lambda = 3 - n_aug;

  //set vector for weights
  VectorXd weights = VectorXd(2*n_aug+1);
  double weight_0 = lambda/(lambda+n_aug);
  weights(0) = weight_0;
  for (int i=1; i<2*n_aug+1; i++) {  
    double weight = 0.5/(n_aug+lambda);
    weights(i) = weight;
  }*/

  int n_z = n_z_radar;
  // TODO

  //create example matrix with predicted sigma points in state space
  // MatrixXd Xsig_pred = MatrixXd(n_x, 2 * n_aug + 1);

  //create example vector for predicted state mean
  VectorXd x_new = VectorXd(n_x);

  //create example matrix for predicted state covariance
  MatrixXd P_new = MatrixXd(n_x,n_x);

  //create example matrix with sigma points in measurement space
  // MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug + 1);

  //create example vector for mean predicted measurement
  // VectorXd z_pred = VectorXd(n_z);

  //create example matrix for predicted measurement covariance
  // MatrixXd S = MatrixXd(n_z,n_z);

  //create example vector for incoming radar measurement
  // VectorXd z = VectorXd(n_z);

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x, n_z);

  //calculate cross correlation matrix
  Tc.setZero();
  for(int i = 0; i < 2*n_aug + 1; i++) {
      VectorXd x_diff = Xsig_pred.col(i) - x;
      VectorXd z_diff = Zsig.col(i) - z_pred;
      
    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;
      
    Tc += weights[i] * x_diff * z_diff.transpose();
  }
  
  //calculate Kalman gain K;
  MatrixXd K = Tc * S.inverse();
  
  VectorXd z_diff = z - z_pred;
  while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
  
  //update state mean and covariance matrix
  x_new = x + K * z_diff;
  P_new = P - K * S * K.transpose();

  //print result
  std::cout << "Updated state x: " << std::endl << x << std::endl;
  std::cout << "Updated state covariance P: " << std::endl << P << std::endl;

  //write result
  *x_out = x_new;
  *P_out = P_new;
}

/*  //measurement covariance matrix - laser
  R_laser_ = MatrixXd(2, 2);
  R_laser_ << 0.0225, 0,
              0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ = MatrixXd(3, 3);
  R_radar_ << 0.09, 0, 0,
              0, 0.0009, 0,
              0, 0, 0.09;

  //measurement matrix
  H_laser_ = MatrixXd(2, 4);
  H_laser_ << 1, 0, 0, 0,
              0, 1, 0, 0;

  //the initial transition matrix F_
  F_ = MatrixXd(4, 4);
  F_ << 1, 0, 1, 0,
        0, 1, 0, 1,
        0, 0, 1, 0,
        0, 0, 0, 1;

  // TODO: Noise matrix
  Q_ = MatrixXd(4, 4); */

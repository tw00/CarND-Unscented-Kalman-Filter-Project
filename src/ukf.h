#ifndef UKF_H
#define UKF_H

#include "measurement_package.h"
#include "Eigen/Dense"
#include <vector>
#include <string>
#include <fstream>
#include "tools.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

class UKF {
public:

  ///* set to true to print results after each step
  bool print_result_;

  ///* initially set to false, set to true in first call of ProcessMeasurement
  bool is_initialized_;

  ///* if this is false, laser measurements will be ignored (except for init)
  bool use_laser_;

  ///* if this is false, radar measurements will be ignored (except for init)
  bool use_radar_;

  ///* state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
  VectorXd x_;

  ///* state covariance matrix
  MatrixXd P_;

  ///* augmented sigma points matrix
  MatrixXd Xsig_aug_;

  ///* predicted sigma points matrix
  MatrixXd Xsig_pred_;

  // Sensor Noise 
  MatrixXd R_laser_;
  MatrixXd R_radar_;

  // Mean predicted measurement for the Radar data
  VectorXd z_pred; 

  ///* weights...
  VectorXd weights;

  ///* time when the state is true, in us
  //long long time_us_;
  long long previous_timestamp_;

  ///* Process noise standard deviation longitudinal acceleration in m/s^2
  double std_a_;

  ///* Process noise standard deviation yaw acceleration in rad/s^2
  double std_yawdd_;

  ///* Laser measurement noise standard deviation position1 in m
  double std_laspx_;

  ///* Laser measurement noise standard deviation position2 in m
  double std_laspy_;

  ///* Radar measurement noise standard deviation radius in m
  double std_radr_;

  ///* Radar measurement noise standard deviation angle in rad
  double std_radphi_;

  ///* Radar measurement noise standard deviation radius change in m/s
  double std_radrd_ ;

  ///* Weights of sigma points
  VectorXd weights_;

  ///* State dimension
  int n_x_;

  ///* Augmented state dimension
  int n_aug_;

  ///* Sigma point spreading parameter
  double lambda_;

  ///* the current NIS for radar
  double NIS_radar_;

  ///* the current NIS for laser
  double NIS_laser_;

  //set state dimension
  int n_x;

  //set augmented dimension
  int n_aug;

  //set measurement dimension, radar can measure r, phi, and r_dot
  int n_z_radar;

  //set measurement dimension, laser can measure x, y
  int n_z_laser;

  //define spreading parameter
  double lambda;

  /**
   * Constructor
   */
  UKF();

  /**
   * Destructor
   */
  virtual ~UKF();

  /**
   * ProcessMeasurement
   * @param meas_package The latest measurement data of either radar or laser
   */
  void ProcessMeasurement(MeasurementPackage measurement_package);

  /**
   * Prediction Predicts sigma points, the state, and the state covariance
   * matrix
   * @param delta_t Time between k and k+1 in s
   */
  void Prediction(double delta_t);

  /**
   * Updates the state and the state covariance matrix using a laser measurement
   * @param meas_package The measurement at k+1
   */
  void UpdateLaser(MeasurementPackage measurement_package);

  /**
   * Updates the state and the state covariance matrix using a radar measurement
   * @param meas_package The measurement at k+1
   */
  void UpdateRadar(MeasurementPackage measurement_package);

  void AugmentedSigmaPoints(MatrixXd* Xsig_out, const VectorXd& x, const MatrixXd& P);
  void SigmaPointPrediction(MatrixXd* Xsig_out, const MatrixXd& Xsig_aug, double delta_t);
  void PredictMeanAndCovariance(VectorXd* x_out, MatrixXd* P_out, const MatrixXd& Xsig_pred);
  void PredictRadarMeasurement(VectorXd* z_out, MatrixXd* S_out, MatrixXd* Zsig_out, const MatrixXd& Xsig_pred);
  void PredictLaserMeasurement(VectorXd* z_out, MatrixXd* S_out, MatrixXd* Zsig_out, const MatrixXd& Xsig_pred);
  void UpdateState(VectorXd* x_out, MatrixXd* P_out, const MatrixXd& Xsig_pred, const VectorXd& x, const MatrixXd& P,
    const MatrixXd& Zsig, const VectorXd& z_pred, const MatrixXd& S, const VectorXd& z, MeasurementPackage::SensorType sensor_type);
};

#endif /* UKF_H */

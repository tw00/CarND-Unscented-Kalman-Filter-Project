#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

double CalculateNIS(const VectorXd &bla){
    return 0.0;
}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth){
  /**
    * Calculate the RMSE here.
  */

	VectorXd rmse(4);
	rmse << 0,0,0,0;

	// check the validity of the following inputs:
	//  * the estimation vector size should not be zero
	//  * the estimation vector size should equal ground truth vector size
	if(estimations.size() != ground_truth.size()
			|| estimations.size() == 0){
        std::cout << "Invalid estimation or ground_truth data" << std::endl;
		return rmse;
	}

	//accumulate squared residuals
	for(unsigned int i=0; i < estimations.size(); ++i){

		VectorXd residual = estimations[i] - ground_truth[i];

		//coefficient-wise multiplication
		residual = residual.array()*residual.array();
		rmse += residual;
	}

	//calculate the mean
	rmse = rmse/estimations.size();

	//calculate the squared root
	rmse = rmse.array().sqrt();

	//return the result
	return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd &x_state) {
  /**
    * Calculate a Jacobian here.
  */

	MatrixXd Hj(3,4);
	//recover state parameters
	float px = x_state(0);
	float py = x_state(1);
	float vx = x_state(2);
	float vy = x_state(3);

	//check division by zero
	float L = px*px + py*py;
	float Ls = sqrt(L);
	float L3 = sqrt(L*L*L);
	if (fabs(L) < 0.00001) {
        std::cout << "Error - Division by zero (Tools::CalculateJacobian)" << std::endl;
        return Hj;
	}
	
	//compute the Jacobian matrix
	Hj << px/Ls, py/Ls, 0, 0,
	      -py/L, px/L, 0, 0,
	      py*(vx*py-vy*px)/L3, px*(vy*px-vx*py)/L3, px/Ls, py/Ls;

	return Hj;
}


VectorXd Tools::fromCartesianToPolar(const VectorXd& x_state) {
  // TODO: COPY AND PASTE

  VectorXd z_pred(3);

  const double distance = sqrt(pow(x_state[0], 2) + pow(x_state[1], 2));

  if(distance < 1e-4)
  {
      // Set angle and distance change rate to zero for very small movement
      z_pred << distance, 0.0, 0.0;
  }
  else
  {
      z_pred << distance,
                atan2(x_state[1], x_state[0]),
                (x_state[0] * x_state[2] + x_state[1] * x_state[3]) / distance;
  }

  return z_pred;
}

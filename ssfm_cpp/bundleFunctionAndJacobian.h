#pragma once 

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SVD>
#include <vector>
#include <string>
#include <unordered_set>



#include <iostream>
#include <unordered_map>
#include <assert.h>
#include <math.h>

#include <tuple>

using namespace Eigen;

using namespace std;


MatrixXd functionForRotationAndTransition(const MatrixXd& parameters,const MatrixXd& variables);

MatrixXd jacobianForPoint(const MatrixXd& parameters,const MatrixXd& variables);

MatrixXd jacobianForRotationAndTransition(const MatrixXd& parameters,const MatrixXd& variables);


MatrixXd functionForRotationAndTransitionUnitLength(const MatrixXd& parameters,const MatrixXd& variables);

MatrixXd jacobianForPointUnitLength(const MatrixXd& parameters,const MatrixXd& variables);

MatrixXd jacobianForRotationAndTransitionUnitLength(const MatrixXd& parameters,const MatrixXd& variables);






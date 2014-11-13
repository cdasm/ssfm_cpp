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


MatrixXd functionForRotationAndTransition(const MatrixXd& ,const MatrixXd&, const MatrixXd&);

MatrixXd jacobianForPoint(const MatrixXd& ,const MatrixXd&, const MatrixXd&);

MatrixXd jacobianForRotationAndTransition(const MatrixXd& ,const MatrixXd& , const MatrixXd&);


MatrixXd functionForRotationAndTransitionUnitLength(const MatrixXd& ,const MatrixXd& , const MatrixXd&);

MatrixXd jacobianForPointUnitLength(const MatrixXd& ,const MatrixXd& , const MatrixXd&);

MatrixXd jacobianForRotationAndTransitionUnitLength(const MatrixXd& ,const MatrixXd&,  const MatrixXd&);






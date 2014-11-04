#include <Eigen/Dense>

#include <Eigen/Sparse>
#include "basicSettings.inl"


using namespace Eigen;


#pragma once

typedef MatrixXd (*funcType)(const MatrixXd& variables,const MatrixXd& parameters);

MatrixXd levenbergM_simple(MatrixXd& dataset,const MatrixXd& obj_vals,funcType func,funcType jfunc,MatrixXd& initParameters,int maxiter_times=170);
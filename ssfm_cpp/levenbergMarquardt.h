#include <Eigen/Dense>

#include <Eigen/Sparse>
#include "basicSettings.inl"
#include <vector>

using namespace std;

using namespace Eigen;

using namespace lM;

#pragma once

typedef MatrixXd (*funcType)(const MatrixXd& variables,const MatrixXd& parameters);

typedef MatrixXd (*funcType2)(const vector<MatrixXd>& variables);

MatrixXd levenbergM_simple(MatrixXd& dataset,const MatrixXd& obj_vals,funcType func,funcType jfunc,MatrixXd& initParameters,int maxiter_times=170);

MatrixXd levenbergM_advanced(MatrixXd& dataset,MatrixXd& assistantPara,const vector<vector<int> >& funcDataMap,const vector<vector<int> >& jfuncDataMap,const MatrixXd& obj_vals,vector<funcType2>& funcs,vector<funcType2>& jfuncs,MatrixXd& initParameters,int maxiter_times=170);


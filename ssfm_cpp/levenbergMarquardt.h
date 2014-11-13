#include <Eigen/Dense>

#include <Eigen/Sparse>
#include "basicSettings.inl"
#include <vector>

using namespace std;

using namespace Eigen;

using namespace lM;

#pragma once

typedef MatrixXd (*funcTypeadhoc)(const MatrixXd&,const MatrixXd&,const MatrixXd& );

typedef MatrixXd (*funcType2)(const vector<MatrixXd>& variables);

MatrixXd levenbergM_adhoc(const vector<MatrixXd>& dataset,const MatrixXd& obj_vals,funcTypeadhoc func,funcTypeadhoc jfunc,MatrixXd& initParameters,int maxiter_times=170);

//MatrixXd levenbergM_advanced(MatrixXd& dataset,MatrixXd& assistantPara,const vector<vector<int> >& funcDataMap,const vector<vector<int> >& jfuncDataMap,const MatrixXd& obj_vals,vector<funcType2>& funcs,vector<funcType2>& jfuncs,MatrixXd& initParameters,int maxiter_times=260);

enum dataORParameter{_data,_parameter};

//MatrixXd levenbergM_advanced(const vector<MatrixXd>& dataset,const vector<pair<int,vector<pair<dataORParameter,pair<int,int>>>>>& funcDataMap,const vector<pair<int,vector<pair<dataORParameter,pair<int,int>>>>>& jfuncDataMap,const MatrixXd& obj_vals,vector<funcType2>& funcs,vector<funcType2>& jfuncs,MatrixXd& initParameters,int dataNumber,int observationNumber,int maxiter_times=260);
MatrixXd levenbergM_advanced(const vector<MatrixXd>& dataset,const vector<pair<pair<int,int>,vector<pair<dataORParameter,pair<int,int>>>>>& funcDataMap,const vector<pair<pair<int,pair<int,int>>,vector<pair<dataORParameter,pair<int,int>>>>>& jfuncDataMap,const MatrixXd& obj_vals,vector<funcType2>& funcs,vector<funcType2>& jfuncs,MatrixXd& initParameters,int dataNumber,int observationNumber,int maxiter_times=260);
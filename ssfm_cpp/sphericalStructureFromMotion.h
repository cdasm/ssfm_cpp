#pragma once 

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SVD>
#include <vector>
#include <string>
#include <unordered_set>
#include "basicSettings.inl"
#include "../FileIO/FileInOut.h"
#include "levenbergMarquardt.h"

#include <iostream>
#include <unordered_map>
#include <assert.h>
#include <math.h>
#include <tuple>

using namespace Eigen;

using namespace sSfm;
using namespace std;



//this is the function which converts from 2D image coordinates to 3D spherical coordinate with reference to image size.
template<class T>
MatrixXd imageCordinate2Phere(vector<T> coordinates,vector<int> imageSize)
{
	double a=(double)coordinates[0];
	double b=(double)coordinates[1];
	double c=(double)imageSize[0];
	double d=(double)imageSize[1];
	MatrixXd y(1,3);

	y(0,0) = cos((PI*a*2.0)/c)*sin((PI*b)/d);
	y(0,1) = -sin((PI*a*2.0)/c)*sin((PI*b)/d);
	y(0,2) = cos((PI*b)/d);

	return y;
}



auto collectFromPairwiseMatching(const string& featureLst,const string& matchLst,vector<vector<vector<int>>>& feas,vector<map<int,int> >& sth, vector<unordered_set<int> >& contain,vector<vector<int> >& featureIsPoint )-> vector<pair<vector<int> ,  vector<int> > >;

//these are the types for how the cameras are defined during bundle adjustment
//_static means the projections on the camera are considered while the position of the camera will not change
//_unitLength means the norm of this cameras transition is 1, and the camera is 5-degree free
//_ordinary means the transition and rotation will be estimated for the camera, and the camera is 6-degree free
enum cameraType{_static,_unitLength,_ordinary};



vector<pair<MatrixXd,MatrixXd> > transitionAndRotationFromEssential(const MatrixXd& essential);


pair<MatrixXd,vector<double> > bestPoints(const MatrixXd& spnts1,const vector<int>& ind1,const MatrixXd& spnts2,const vector<int>& ind2,const vector<MatrixXd>& transitions,const vector<MatrixXd>& rotations);




//auto geometricReconstructionFrom2Frames(const MatrixXd& sphericalPoints1,const vector<int>& ind1,const MatrixXd& sphericalPoints2,const vector<int>& ind2,vector<double>& transition,vector<double >& rotation)->pair<MatrixXd,vector<double> >;
auto geometricReconstructionFrom2Frames(const MatrixXd& sphericalPoints1,const vector<int>& ind1,const MatrixXd& sphericalPoints2,const vector<int>& ind2,MatrixXd& transition,MatrixXd& rotation)->pair<MatrixXd,vector<double> >;


auto threeDimensionReconstruction(const string& featureFileName,const string& matchFileName)->tuple<MatrixXd,MatrixXd,MatrixXd>;



MatrixXd bestPoint(const MatrixXd& p, const MatrixXd& u);



MatrixXd functionForRotationAndTransition(const MatrixXd& parameters,const MatrixXd& variables);

MatrixXd jacobianForPoint(const MatrixXd& parameters,const MatrixXd& variables);

MatrixXd jacobianForRotationAndTransition(const MatrixXd& parameters,const MatrixXd& variables);

MatrixXd functionForRotationAndTransition(const vector<MatrixXd>& input);

MatrixXd jacobianForRotationAndTransition(const vector<MatrixXd>& input);

MatrixXd jacobianForPoint(const vector<MatrixXd>& input);


MatrixXd functionForRotationAndTransitionUnitLength(const MatrixXd& parameters,const MatrixXd& variables);

MatrixXd jacobianForPointUnitLength(const MatrixXd& parameters,const MatrixXd& variables);

MatrixXd jacobianForRotationAndTransitionUnitLength(const MatrixXd& parameters,const MatrixXd& variables);

MatrixXd functionForRotationAndTransitionUnitLength(const vector<MatrixXd>& input);

MatrixXd jacobianForRotationAndTransitionUnitLength(const vector<MatrixXd>& input);

MatrixXd jacobianForPointUnitLength(const vector<MatrixXd>& input);


MatrixXd estimateCameraParameter(const MatrixXd& projPoints,const MatrixXd& points);

MatrixXd estimateCameraParameter(const MatrixXd& projPoints,const vector<int>& ind1,const MatrixXd& points,const vector<int>& ind2);

/*
template<class T>
auto geometricReconstructionFrom2Frames(const vector<vector<T> >& pnts1,const vector<vector<T> >& pnts2,vector<double>& transition,vector<double >& rotation,const vector<int>& imageSize)->pair<MatrixXd,vector<double> >
{
	assert(pnts1.size()==pnts2.size());

	MatrixXd sphericalPoints1(pnts1.size(),3);
	MatrixXd sphericalPoints2(pnts2.size(),3);

	for (int i = 0; i < pnts1.size(); i++)
	{
		sphericalPoints1.row(i)=imageCordinate2Phere(pnts1[i],imageSize);
		sphericalPoints2.row(i)=imageCordinate2Phere(pnts2[i],imageSize);
	}
	return geometricReconstructionFrom2Frames(sphericalPoints1,sphericalPoints2,transition,rotation);
	
}*/

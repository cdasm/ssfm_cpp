#pragma once 

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <vector>
#include <string>
#include <unordered_set>
#include "basicSettings.inl"
#include "../FileIO/FileInOut.h"

#include <iostream>
#include <assert.h>
using namespace Eigen;

using namespace sSfm;
using namespace std;



//this is the function which converts from 2D image coordinates to 3D spherical coordinate with reference to image size.
template<class T>
vector<double> imageCordinate2Phere(vector<T> coordinates,vector<T> imageSize)
{
	double a=(double)coordinates[0];
	double b=(double)coordinates[1];
	double c=(double)imageSize[0];
	double d=(double)imageSize[1];
	vector<double> y(3);

	y[0] = cos((PI*a*2.0)/c)*sin((PI*b)/d);
	y[1] = -sin((PI*a*2.0)/c)*sin((PI*b)/d);
	y[2] = cos((PI*b)/d);

	return y;
}



auto collectFromPairwiseMatching(const string& featureLst,const string& matchLst,vector<vector<vector<int>>>& feas,vector<map<int,int> >& sth, vector<unordered_set<int> >& contain,vector<vector<int> >& featureIsPoint )-> pair<vector<pair<int,int> >,  vector<vector<int> > >;

//these are the types for how the cameras are defined during bundle adjustment
//_static means the projections on the camera are considered while the position of the camera will not change
//_unitLength means the norm of this cameras transition is 1, and the camera is 5-degree free
//_ordinary means the transition and rotation will be estimated for the camera, and the camera is 6-degree free
enum cameraType{_static,_unitLength,_ordinary};

template<class T>
auto geometricReconstructionFrom2Frames(const vector<vector<T> >& pnts1,const vector<vector<T> >& pnts2,vector<double>& transition,vector<vector<double> > rotation)->vector<vector<double> >
{
	assert(pnts1.size()==pnts2.size());

	vector<vector<double> > sphericalPoints1(pnts1.size());
	vector<vector<double> > sphericalPoints2(pnts2.size());

	for (int i = 0; i < pnts1.size(); i++)
	{
		sphericalPoints1[i]=
	}

}
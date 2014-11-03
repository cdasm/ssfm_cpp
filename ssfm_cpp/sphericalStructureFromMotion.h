#pragma once 

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SVD>
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



vector<pair<MatrixXd,MatrixXd> > transitionAndRotationFromEssential(const MatrixXd& essential);

template<class T>
auto geometricReconstructionFrom2Frames(const vector<vector<T> >& pnts1,const vector<vector<T> >& pnts2,vector<double>& transition,vector<vector<double> >& rotation,const vector<T>& imageSize)->pair<vector<vector<double> >,vector<double> >
{
	assert(pnts1.size()==pnts2.size());

	vector<vector<double> > sphericalPoints1(pnts1.size());
	vector<vector<double> > sphericalPoints2(pnts2.size());

	for (int i = 0; i < pnts1.size(); i++)
	{
		sphericalPoints1[i]=imageCordinate2Phere(pnts1[i],imageSize);
		sphericalPoints2[i]=imageCordinate2Phere(pnts2[i],imageSize);
	}

	vector<vector<double> > points(pnts1.size(),vector<double>(3,0.0));
	vector<double> errors(pnts1.size(),-1);


	MatrixXd observation=MatrixXd::Zero(points.size(),9);
	
	for (int i = 0; i < points.size(); i++)
	{
		//a,&b,&c,&d,&e,&f;
		double &a=sphericalPoints1[i][0];
		double &b=sphericalPoints1[i][1];
		double &c=sphericalPoints1[i][2];
		double &d=sphericalPoints2[i][0];
		double &e=sphericalPoints2[i][1];
		double &f=sphericalPoints2[i][2];

		observation(i,0) = a*d;
		observation(i,1) = b*d;
		observation(i,2) = c*d;
		observation(i,3) = a*e;
		observation(i,4) = b*e;
		observation(i,5) = c*e;
		observation(i,6) = a*f;
		observation(i,7) = b*f;
		observation(i,8) = c*f;


	}

	JacobiSVD<decltype(observation)> svd(observation, Eigen::ComputeFullU |
                                        Eigen::ComputeFullV);

	//std::cout << "singular values" << std::endl
    //      << svd.singularValues() << std::endl;
	//std::cout << "matrix U" << std::endl << svd.matrixU() << std::endl;
	std::cout << "matrix V" << std::endl << svd.matrixV() << std::endl;

	vector<MatrixXd> essentialMatrixes(2,MatrixXd::Zero(3,3));

	auto matrixV=svd.matrixV();

	for (int i = 0; i < 2; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			for (int k = 0; k < 3; k++)
			{
				essentialMatrixes[i](j,k)=matrixV (j*3+k,7+i);
			}
		}
	}

	vector<pair<MatrixXd,MatrixXd> > transtionAndRotations;
	for (int i = 0; i < 2; i++)
	{
		auto ttar=transitionAndRotationFromEssential(essentialMatrixes[i]);
		transtionAndRotations.insert(transtionAndRotations.end(),ttar.begin(),ttar.end());
	}
	

	//cout<<observation;
	return make_pair(points,errors);

}
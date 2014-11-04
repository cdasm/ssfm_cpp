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
#include <assert.h>
#include <math.h>

using namespace Eigen;

using namespace sSfm;
using namespace std;



//this is the function which converts from 2D image coordinates to 3D spherical coordinate with reference to image size.
template<class T>
MatrixXd imageCordinate2Phere(vector<T> coordinates,vector<T> imageSize)
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



auto collectFromPairwiseMatching(const string& featureLst,const string& matchLst,vector<vector<vector<int>>>& feas,vector<map<int,int> >& sth, vector<unordered_set<int> >& contain,vector<vector<int> >& featureIsPoint )-> pair<vector<pair<int,int> >,  vector<vector<int> > >;

//these are the types for how the cameras are defined during bundle adjustment
//_static means the projections on the camera are considered while the position of the camera will not change
//_unitLength means the norm of this cameras transition is 1, and the camera is 5-degree free
//_ordinary means the transition and rotation will be estimated for the camera, and the camera is 6-degree free
enum cameraType{_static,_unitLength,_ordinary};



vector<pair<MatrixXd,MatrixXd> > transitionAndRotationFromEssential(const MatrixXd& essential);


pair<MatrixXd,vector<bool> > bestPoints(const MatrixXd& spnts1,const MatrixXd& spnts2,const vector<MatrixXd>& transitions,const vector<MatrixXd>& rotations);


template<class T>
auto geometricReconstructionFrom2Frames(const vector<vector<T> >& pnts1,const vector<vector<T> >& pnts2,vector<double>& transition,vector<double >& rotation,const vector<T>& imageSize)->pair<vector<vector<double> >,vector<double> >
{
	assert(pnts1.size()==pnts2.size());

	MatrixXd sphericalPoints1(pnts1.size(),3);
	MatrixXd sphericalPoints2(pnts2.size(),3);

	for (int i = 0; i < pnts1.size(); i++)
	{
		sphericalPoints1.row(i)=imageCordinate2Phere(pnts1[i],imageSize);
		sphericalPoints2.row(i)=imageCordinate2Phere(pnts2[i],imageSize);
	}

	vector<vector<double> > points(pnts1.size(),vector<double>(3,0.0));
	vector<double> errors(pnts1.size(),-1);


	MatrixXd observation(points.size(),9);

	auto convert=[](const MatrixXd& from)->MatrixXd
	{
		MatrixXd A0(1,9);
		double a,b,c,d,e,f,g,h,i;
		a=from(0,0) ;
		b=from(0,1) ;
		c=from(0,2) ;
		d=from(1,0) ;
		e=from(1,1) ;
		f=from(1,2) ;
		g=from(2,0) ;
		h=from(2,1) ;
		i=from(2,2) ;
		A0(0,0) = a;
		A0(0,1) = d;
		A0(0,2) = g;
		A0(0,3) = b;
		A0(0,4) = e;
		A0(0,5) = h;
		A0(0,6) = c;
		A0(0,7) = f;
		A0(0,8) = i;
		return A0;
	};
	
	for (int i = 0; i < points.size(); i++)
	{
		//a,&b,&c,&d,&e,&f;
		//double &a=sphericalPoints1[i][0];
		//double &b=sphericalPoints1[i][1];
		//double &c=sphericalPoints1[i][2];
		//double &d=sphericalPoints2[i][0];
		//double &e=sphericalPoints2[i][1];
		//double &f=sphericalPoints2[i][2];/

		MatrixXd tob=sphericalPoints1.row(i).transpose()*sphericalPoints2.row(i);

		observation.row(i)=convert(tob);
		/*
		observation(i,0) = a*d;
		observation(i,1) = b*d;
		observation(i,2) = c*d;
		observation(i,3) = a*e;
		observation(i,4) = b*e;
		observation(i,5) = c*e;
		observation(i,6) = a*f;
		observation(i,7) = b*f;
		observation(i,8) = c*f;
		*/

	}

	JacobiSVD<decltype(observation)> svd(observation, Eigen::ComputeFullU |
                                        Eigen::ComputeFullV);
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
	
	vector<pair<MatrixXd,vector<bool> > > bestPointCandidates(transtionAndRotations.size());

	for (int i = 0; i < transtionAndRotations.size(); i++)
	{
		vector<MatrixXd> trans(2,MatrixXd::Zero(1,3));
		trans[1]=transtionAndRotations[i].first;
		vector<MatrixXd> rots(2,MatrixXd::Identity(3,3));
		rots[1]=transtionAndRotations[i].second;
		bestPointCandidates[i]=bestPoints(sphericalPoints1,sphericalPoints2,trans,rots);
	}
	//cout<<observation;
	return make_pair(points,errors);

}
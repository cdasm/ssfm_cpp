#include <iostream>
#include <Windows.h>

#include <vector>
#include <map>
#include <string>
#include <direct.h>
#include <Eigen/Dense>
#include <algorithm>
#include "sphericalStructureFromMotion.h"
#include "utinity.h"
using namespace Eigen;
using namespace std;

int main_()
{
  

	auto test=imageCordinate2Phere(vector<double> (2,10),vector<int> (2,28));
	

	MatrixXd m = MatrixXd::Random(3,3);
	cout<<m.cols()<<endl;
	cout<<m.rows()<<endl;

	cout<<m;
	m = (m + MatrixXd::Constant(3,3,1.2)) ;
	cout << "m =" << endl << m << endl;
	VectorXd v(3);
	v << 1, 2, 3;
	
	cout << "m * v =" << endl << m * v << endl;

	cout<<"v*v'=\n"<<v*v.transpose()<<endl;
	cout<<"v'*v=\n"<<v.transpose()*v<<endl;

	//getchar();
	return 0;
}


void testTwoFrame()
{
	MatrixXd sp1(30,3);
	MatrixXd sp2(30,3);

	sp1<<0.6248,0.4874,0.6100
	,0.3779,0.6295,0.6789
	,0.6384,0.4412,0.6307
	,0.6208,0.5933,0.5125
	,0.5969,0.5685,0.5662
	,0.6550,0.4276,0.6230
	,0.6395,0.5994,0.4814
	,0.6450,0.6019,0.4708
	,0.4480,0.5155,0.7304
	,0.5790,0.5367,0.6137
	,0.4922,0.6328,0.5978
	,0.4394,0.4795,0.7596
	,0.6674,0.6109,0.4258
	,0.5876,0.4783,0.6526
	,0.5846,0.5758,0.5716
	,0.5991,0.6539,0.4620
	,0.4611,0.6695,0.5824
	,0.7100,0.4919,0.5040
	,0.6750,0.5977,0.4326
	,0.6289,0.4401,0.6409
	,0.5891,0.5082,0.6282
	,0.6026,0.6735,0.4280
	,0.4861,0.7611,0.4294
	,0.6335,0.6817,0.3660
	,0.6275,0.6159,0.4764
	,0.6254,0.6495,0.4325
	,0.6282,0.5811,0.5174
	,0.4045,0.6298,0.6631
	,0.5369,0.7175,0.4437
	,0.5034,0.6270,0.5946;

	sp2<< 0.6952,0.3884,0.6049
	,0.3870,0.5682,0.7262
	,0.7195,0.3059,0.6235
	,0.6908,0.5309,0.4908
	,0.6604,0.4984,0.5616
	,0.7337,0.2989,0.6101
	,0.7081,0.5414,0.4533
	,0.7172,0.5439,0.4356
	,0.4721,0.4006,0.7853
	,0.6560,0.4288,0.6211
	,0.5342,0.5797,0.6152
	,0.4625,0.3543,0.8128
	,0.7429,0.5562,0.3724
	,0.6534,0.3676,0.6618
	,0.6471,0.5068,0.5696
	,0.6694,0.6157,0.4158
	,0.4944,0.6286,0.6004
	,0.7997,0.3902,0.4563
	,0.7492,0.5395,0.3843
	,0.6985,0.3252,0.6374
	,0.6529,0.4157,0.6332
	,0.6641,0.6413,0.3844
	,0.5121,0.7618,0.3968
	,0.6938,0.6517,0.3065
	,0.6961,0.5625,0.4461
	,0.6912,0.6083,0.3901
	,0.6998,0.5145,0.4955
	,0.4249,0.5706,0.7028
	,0.5811,0.6987,0.4173
	,0.5485,0.5722,0.6097;

	vector<double> tt(3,0);
	vector<double> rr(3,0);

	geometricReconstructionFrom2Frames(sp1,sp2,tt,rr);
}

void testJacobian()
{
	MatrixXd para1(1,6);
	para1<<1, 2, 3, 1, 1, 1;
	MatrixXd var1(1,6);
	var1<<3, 2, 1, 2, 1, 2;
	cout<<functionForRotationAndTransition(para1,var1)<<endl;
	cout<<functionForRotationAndTransition2(para1,var1)<<endl;
	cout<<jacobianForRotationAndTransition(para1,var1)<<endl;
	cout<<jacobianForRotationAndTransition2(para1,var1)<<endl;
	cout<<jacobianForPoint(para1,var1)<<endl;
	cout<<jacobianForPoint2(para1,var1)<<endl;
}

int main()
{
	testJacobian();
	 Eigen::MatrixXf m(4,4);
	  m <<  1, 2, 3, 4,
			5, 6, 7, 8,
			9,10,11,12,
		   13,14,15,16;
	  cout << "Block in the middle" << endl;
	  cout << m.block<2,2>(1,1) << endl << endl;
	  for (int i = 1; i <= 4; ++i)
	  {
		cout << "Block of size " << i << "x" << i << endl;
		cout << m.block(0,0,i,i) << endl << endl;
	  }

	MatrixXd test=MatrixXd::Random(3,4);

	cout<<test<<endl;

	MatrixXd vv=MatrixXd::Zero(1,3);
	test.block(0,1,1,3)=vv;
	cout<<test;
	//testTwoFrame();
	MatrixXd p(2,3);
	MatrixXd u(2,3);

	p<<0,0,0,1,1,1;
	u<<2,4,5,3,2,1;

	auto pnt=bestPoint(p,u);

	cout<<pnt<<endl;

	cout<<MatrixXd::Identity(3,3)<<endl;
	cout<<acos(-0.5)<<endl;
	cout<<acos(0)<<endl;
	cout<<acos(1)<<endl;
	

	main_();
	
	_chdir("D:\\agood");

	vector<vector<vector<int> > > features;
	vector<map<int,int>> correspondences;

	vector<unordered_set<int> > contain;

	vector<vector<int> > featureIsPoint;

	auto trajectories=collectFromPairwiseMatching("orb.lst","match.lst",features,correspondences,contain,featureIsPoint);



	vector<vector<double> > transitions(features.size(),vector<double>(3,0.0));



	vector<vector<double > > rotations(features.size(),vector<double>(3,0));

	vector<vector<double> > reconstructedPoints(trajectories.first.size(),vector<double>(3,0.0));
	
	vector<bool> alreadyReconstructed(trajectories.first.size(),false);

	vector<vector<int> > points1(correspondences[0].size()),points2(correspondences[0].size());

	vector<int> index(correspondences[0].size());

	int count=0;
	for ( auto&k:correspondences[0])
	{
		points1[count]=features[0][k.first];
		points2[count]=features[1][k.second];
		index[count]=featureIsPoint[0][k.first];
		++count;
	}

	vector<int> imageSize(2);
	imageSize[0]=512;imageSize[1]=256;

	geometricReconstructionFrom2Frames(points1,points2,transitions[1],rotations[1],imageSize);
	//set_intersection(contain[0].begin(),contain[0].end(),contain[1].begin(),contain[1].end(),back_inserter(pointOnFirst2Frames));
	
	//auto pointOnFirst2Frames=set_intersect(contain[0],contain[1]);



	return 0;

}
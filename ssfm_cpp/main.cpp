#include <iostream>
#include <Windows.h>

#include <vector>
#include <map>
#include <string>
#include <direct.h>
#include <Eigen/Dense>
#include <algorithm>
#include <fstream>
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

	MatrixXd tt(1,3);
	MatrixXd rr(1,3);

	vector<int> ind(30);
	for (int i = 0; i < 30; i++)
	{
		ind[i]=i;
	}


	geometricReconstructionFrom2Frames(sp1,ind,sp2,ind,tt,rr);
}

void testJacobian()
{
	MatrixXd para1(1,6);
	para1<<1, 2, 3, 1, 1, 1;
	MatrixXd var1(1,6);
	var1<<3, 2, 1, 2, 1, 2;
	cout<<functionForRotationAndTransition(para1,var1)<<endl;
	cout<<functionForRotationAndTransitionUnitLength(para1,var1)<<endl;
	cout<<jacobianForRotationAndTransition(para1,var1)<<endl;
	cout<<jacobianForRotationAndTransitionUnitLength(para1,var1)<<endl;
	cout<<jacobianForPoint(para1,var1)<<endl;
	cout<<jacobianForPointUnitLength(para1,var1)<<endl;
}


void testLM()
{
	MatrixXd pnts(30,3),prjp(30,3);
	pnts<<18.1472,17.0605,17.5127
		,19.0579,10.3183,12.5510
		,11.2699,12.7692,15.0596
		,19.1338,10.4617,16.9908
		,16.3236,10.9713,18.9090
		,10.9754,18.2346,19.5929
		,12.7850,16.9483,15.4722
		,15.4688,13.1710,11.3862
		,19.5751,19.5022,11.4929
		,19.6489,10.3445,12.5751
		,11.5761,14.3874,18.4072
		,19.7059,13.8156,12.5428
		,19.5717,17.6552,18.1428
		,14.8538,17.9520,12.4352
		,18.0028,11.8687,19.2926
		,11.4189,14.8976,13.4998
		,14.2176,14.4559,11.9660
		,19.1574,16.4631,12.5108
		,17.9221,17.0936,16.1604
		,19.5949,17.5469,14.7329
		,16.5574,12.7603,13.5166
		,10.3571,16.7970,18.3083
		,18.4913,16.5510,15.8526
		,19.3399,11.6261,15.4972
		,16.7874,11.1900,19.1719
		,17.5774,14.9836,12.8584
		,17.4313,19.5974,17.5720
		,13.9223,13.4039,17.5373
		,16.5548,15.8527,13.8045
		,11.7119,12.2381,15.6782;

	prjp<<0.5415, 0.6478, 0.5358
		, 0.7348, 0.4979, 0.4607
		, 0.4396, 0.6619, 0.6072
		, 0.6579, 0.4665, 0.5913
		, 0.5540, 0.4931, 0.6708
		, 0.3022, 0.7195, 0.6252
		, 0.4227, 0.7365, 0.5280
		, 0.6245, 0.6554, 0.4248
		, 0.6032, 0.7278, 0.3262
		, 0.7440, 0.4903, 0.4540
		, 0.3779, 0.6470, 0.6623
		, 0.6890, 0.5933, 0.4163
		, 0.5575, 0.6378, 0.5314
		, 0.5053, 0.7639, 0.4014
		, 0.5778, 0.4997, 0.6454
		, 0.4341, 0.7405, 0.5131
		, 0.5561, 0.7052, 0.4399
		, 0.6343, 0.6665, 0.3917
		, 0.5528, 0.6650, 0.5021
		, 0.6002, 0.6657, 0.4433
		, 0.6258, 0.6035, 0.4941
		, 0.3095, 0.7163, 0.6254
		, 0.5774, 0.6489, 0.4955
		, 0.6716, 0.5126, 0.5350
		, 0.5586, 0.4925, 0.6674
		, 0.6210, 0.6534, 0.4330
		, 0.4909, 0.7058, 0.5107
		, 0.4789, 0.6083, 0.6330
		, 0.5691, 0.6805, 0.4615
		, 0.4534, 0.6298, 0.6306;

	cout<<estimateCameraParameter(prjp,pnts);
}

void testInter()

{
	testLM();
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
}


void testInverseWithSolve()
{
	MatrixXd A=MatrixXd::Random(3,3);
	VectorXd b=VectorXd::Random(3);

//	auto an1=A.colPivHouseholderQr().solve(b);
	VectorXd c=A.row(0);
	cout<<A<<endl;
	cout<<c<<endl;

	MatrixXd d(1,3);

	d=c.transpose();
	cout<<d<<endl;

	VectorXd x=A.transpose().colPivHouseholderQr().solve(b);

//	dp=d*J*H_lm.inverse();
//	VectorXd tosdj=d*J;
//	VectorXd solution=H_lm..transpose().colPivHouseholderQr().solve(tosdj);
//	dp=solution.transpose();

	cout << x <<"\n"<< endl;
	cout<<b.transpose()*A.inverse()<<endl;
}

int main()
{
	cout<<1e10/1e8<<endl;
//	testInverseWithSolve();

	
	
	_chdir("D:\\agood");
//set_intersection(contain[0].begin(),contain[0].end(),contain[1].begin(),contain[1].end(),back_inserter(pointOnFirst2Frames));
	
	//auto pointOnFirst2Frames=set_intersect(contain[0],contain[1]);

	auto pointsCameras=threeDimensionReconstruction("orb.lst","match.lst");
	ofstream f1;
	f1.open("points.txt");
	f1<<get<2>(pointsCameras);
	f1.close();
	f1.open("transitions.txt");
	f1<<get<0>(pointsCameras);
	f1.close();
	f1.open("rotations.txt");
	f1<<get<1>(pointsCameras);
	f1.close();

	

	return 0;

}
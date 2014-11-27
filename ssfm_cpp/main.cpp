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


//	geometricReconstructionFrom2Frames(sp1,ind,sp2,ind,tt,rr);
}
/*
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
}*/


void testLM()
{
	MatrixXd pnts(30,3),prjp(30,3);
	pnts<<15.6882,10.8444,13.3772
	,14.6939,13.9978,19.0005
	,10.1190,12.5987,13.6925
	,13.3712,18.0007,11.1120
	,11.6218,14.3141,17.8025
	,17.9428,19.1065,13.8974
	,13.1122,11.8185,12.4169
	,15.2853,12.6380,14.0391
	,11.6565,11.4554,10.9645
	,16.0198,11.3607,11.3197
	,12.6297,18.6929,19.4205
	,16.5408,15.7970,19.5613
	,16.8921,15.4986,15.7521
	,17.4815,11.4495,10.5978
	,14.5054,18.5303,12.3478
	,10.8382,16.2206,13.5316
	,12.2898,13.5095,18.2119
	,19.1334,15.1325,10.1540
	,11.5238,14.0181,10.4302
	,18.2582,10.7597,11.6899
	,15.3834,12.3992,16.4912
	,19.9613,11.2332,17.3172
	,10.7818,11.8391,16.4775
	,14.4268,12.3995,14.5092
	,11.0665,14.1727,15.4701
	,19.6190,10.4965,12.9632
	,10.0463,19.0272,17.4469
	,17.7491,19.4479,11.8896
	,18.1730,14.9086,16.8678
	,18.6869,14.8925,11.8351
	;

	prjp<< 0.7245, 0.3448, 0.5969
		, 0.5344, 0.4054, 0.7417
		, 0.4945, 0.5245, 0.6931
		, 0.5733, 0.6831, 0.4525
		, 0.4487, 0.4774, 0.7555
		, 0.6429, 0.5832, 0.4966
		, 0.6590, 0.4539, 0.5997
		, 0.6724, 0.4203, 0.6093
		, 0.6523, 0.4997, 0.5698
		, 0.7720, 0.3895, 0.5023
		, 0.4164, 0.5683, 0.7097
		, 0.5597, 0.4344, 0.7057
		, 0.6399, 0.4729, 0.6058
		, 0.8136, 0.3737, 0.4455
		, 0.5822, 0.6536, 0.4835
		, 0.4718, 0.6431, 0.6032
		, 0.4738, 0.4310, 0.7680
		, 0.7853, 0.4872, 0.3820
		, 0.6010, 0.6221, 0.5017
		, 0.8153, 0.3185, 0.4836
		, 0.6236, 0.3730, 0.6870
		, 0.7239, 0.2620, 0.6382
		, 0.4680, 0.4164, 0.7795
		, 0.6405, 0.4181, 0.6441
		, 0.4741, 0.5286, 0.7041
		, 0.8139, 0.2755, 0.5115
		, 0.3471, 0.6466, 0.6793
		, 0.6590, 0.6188, 0.4275
		, 0.6578, 0.4199, 0.6253
		, 0.7563, 0.4706, 0.4544;


	vector<int> ind(30);
	for (int i = 0; i < 30; i++)
	{
		ind[i]=i;
	}

	cout<<estimateCameraParameter(prjp,ind,pnts,ind);
}

void testInter()

{
	testLM();
//	testJacobian();
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

//	auto pnt= reconstructPoint (p,u);

//	cout<<pnt<<endl;

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

void testBestPoint()
{
	MatrixXd a(2,3),b(2,3);
	a<< 0  ,       0 ,        0,
    0.3793,    0.2698 ,   2.6540;

	b<<  0.5960 ,   0.5603  ,  0.5752,
    0.6211  ,  0.5869 ,   0.5194;
//	reconstructPoint(a,b);
}

int main()
{

	
	_chdir("D:\\agood");

	auto pointsCameras=threeDimensionReconstruction("torb.lst","tmatch.lst",512,256);
	ofstream f1;
	f1.open("points.txt");
	f1<<(pointsCameras.second);
	f1.close();
	//p1 p2 p3
	f1.open("cameras.txt");
	f1<<(pointsCameras.first);

	//r1 r2 r3 t1 t2 t3
	f1.close();


	

	return 0;

}
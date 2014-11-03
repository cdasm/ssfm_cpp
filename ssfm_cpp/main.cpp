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
  

	auto test=imageCordinate2Phere(vector<double> (2,10),vector<double> (2,28));
	

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

	getchar();
	return 0;
}


int main()
{
	main_();
	
	_chdir("D:\\agood");

	vector<vector<vector<int> > > features;
	vector<map<int,int>> correspondences;

	vector<unordered_set<int> > contain;

	vector<vector<int> > featureIsPoint;

	auto trajectories=collectFromPairwiseMatching("orb.lst","match.lst",features,correspondences,contain,featureIsPoint);



	vector<vector<double> > transitions(features.size(),vector<double>(3,0.0));

	vector<vector<double> > eyeI(3,vector<double>(3,0));
	eyeI[0][0]=1;eyeI[1][1]=1;eyeI[2][2]=1;

	vector<vector<vector<double> > > rotations(features.size(),eyeI);

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
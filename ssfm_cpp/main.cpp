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
	m = (m + MatrixXd::Constant(3,3,1.2)) * 50;
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
	//set_intersection(contain[0].begin(),contain[0].end(),contain[1].begin(),contain[1].end(),back_inserter(pointOnFirst2Frames));
	
	//auto pointOnFirst2Frames=set_intersect(contain[0],contain[1]);



	return 0;

}
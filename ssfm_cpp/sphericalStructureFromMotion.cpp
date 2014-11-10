#include "sphericalStructureFromMotion.h"


template<class Pnt>
auto incrementalTrajectoryDetect(const vector<vector<Pnt> >& features, vector<map<int,int>>& sth)->   vector< pair< vector<int> ,  vector<int>>>
{



	vector<pair<vector<int> ,vector<int> > > result;

	vector<vector<bool> > markers(features.size());	

	for (size_t i = 0; i < features.size(); i++)
	{
		markers[i].resize(features[i].size(),false);
	}

	for (int i = 0; i < features.size()-1; i++)
	{
		if(i%30==0)
			cout<<"processing frame " <<i<<endl;
		for (int j = 0; j < markers[i].size(); j++)
		{
			if(!markers[i][j])
			{
				vector<int> stend;
				vector<int> traj;
		
				stend.push_back(i);
				traj.push_back(j);
				markers[i][j]=true;
				
				int cur_frame=i;
				int cur_index=j;
		
				while( (cur_frame<features.size()-1) && (sth[cur_frame].count(cur_index)))
				{
					int indd=sth[cur_frame][cur_index];
					stend.push_back(cur_frame+1);		
					traj.push_back(indd);

					markers[cur_frame+1][indd]=true;
				
					cur_index=indd;
					++cur_frame;			
				}
				if(stend.size()>1)
				{

					result.push_back(make_pair(stend,traj));
				}
			}
		}
	}
	cout<<"feature tracing finished"<<endl;
	return result;
}



auto collectFromPairwiseMatching(const string& featureLst,const string& matchLst,vector<vector<vector<int>>>& feas,vector<map<int,int> >& sth,vector<unordered_set<int> >& contain,vector<vector<int> >& featureIsPoint)->  vector<pair<vector<int> ,  vector<int> > >
{

	auto feaNames=fileIOclass::InVectorString(featureLst);
	auto mathces=fileIOclass::InVectorString(matchLst);

	vector<vector<vector<int> > > correpss(mathces.size());
	feas.resize(feaNames.size());
	featureIsPoint.resize(feaNames.size());

	for (int i = 0; i < feaNames.size(); i++)
	{
		feas[i]=fileIOclass::InVectorSInt(feaNames[i],2);
		featureIsPoint[i].resize(feas[i].size(),-1);
	}

	for (int i = 0; i < correpss.size(); i++)
	{
		correpss[i]=fileIOclass::InVectorSInt(mathces[i],2);
	}


	sth.resize(correpss.size());

	for(int i=0;i<correpss.size();++i)
	{
		for(int j=0;j<correpss[i].size();++j)
			sth[i][correpss[i][j][0]]=correpss[i][j][1];
	}

	auto trajs=incrementalTrajectoryDetect(feas,sth);
	
	contain.resize(feas.size());

	for (int i = 0; i < trajs.size(); i++)
	{
		for (int j =0;j < trajs[i].first.size(); j++)
		{
			contain[ trajs[i].first[j] ].insert(i);
			featureIsPoint[ trajs[i].first[j] ] [ trajs[i].second[j] ]=i;
		}
	}

	return trajs;
}


MatrixXd rotationThomason(const MatrixXd& para)
{
	double a,b,c;
	a=para(0,0);
	b=para(0,1);
	c=para(0,2);
	MatrixXd A0=MatrixXd::Zero(3,3);
	double t2 ,t3 ,t4 ,t5 ,t6 ,t7 ,t8 ,t9 ,t10 ;
	t2 = a*a;
	t3 = b*b;
	t4 = c*c;
	t5 = t2+t3+t4+4.0;
	t6 = 1.0/t5;
	t7 = c*2.0;
	t8 = b*2.0;
	t9 = a*c;
	t10 = a*2.0;
	A0(0,0) = t6*(t2*2.0+8.0)-1.0;
	A0(0,1) = t6*(t7-a*b)*-2.0;
	A0(0,2) = t6*(t8+t9)*2.0;
	A0(1,0) = t6*(t7+a*b)*2.0;
	A0(1,1) = t6*(t3*2.0+8.0)-1.0;
	A0(1,2) = t6*(t10-b*c)*-2.0;
	A0(2,0) = t6*(t8-t9)*-2.0;
	A0(2,1) = t6*(t10+b*c)*2.0;
	A0(2,2) = t6*(t4*2.0+8.0)-1.0;


	return A0;
}

MatrixXd rotationThomasonPara(const MatrixXd& mtr)
{

	MatrixXd m=MatrixXd::Zero(4,4);
	m.block(1,1,3,3)=mtr;
	double a,b,c,d,e,x,y,z;
	MatrixXd result=MatrixXd::Zero(1,3);

	a=m(1,1)+m(2,2)+m(3,3);
	b=16/(a+1);

	c=(m(1,3)-m(3,1))*b;

	y=c/8;

	d=(m(2,1)-m(1,2))*b;
	z=d/8;

	e=(m(3,2)-m(2,3))*b;

	x=e/8;

	result(0,0)=x;
	result(0,1)=y;
	result(0,2)=z;
	return result;
}

MatrixXd transitionFrom2Para(const MatrixXd& inp)
{
	double a,b;
	a=inp(0,0) ;
	b=inp(0,1) ;

	MatrixXd A0(1,3);
	double t2 ;
	t2 = cos(a);
	A0(0,0) = t2*cos(b);
	A0(0,1) = t2*sin(b);
	A0(0,2) = sin(a);
	return A0;
}

MatrixXd transition2Para(const MatrixXd& inp)
{
	double t1,t2,t3;
	t1=inp(0,0);
	t2=inp(0,1);
	t3=inp(0,2);
	double a,b;

	b=atan2(t2,t1);
	a=atan2(t3,sqrt(t1*t1+t2*t2));
	MatrixXd A0(1,2);
	A0(0,0)=a;
	A0(0,1)=b;
	return A0;
}

double distanceBetweenPointLine(const MatrixXd&x,const MatrixXd&p,const MatrixXd& u)
{
	double x1,x2,x3,p1,p2,p3,u1,u2,u3;
	x1=x(0,0);
	x2=x(0,1);
	x3=x(0,2);
	p1=p(0,0);
	p2=p(0,1);
	p3=p(0,2);
	u1=u(0,0);
	u2=u(0,1);
	u3=u(0,2);
	double t3 ,t4 ,t5 ,t6 ,t7 ,t8 ,t9 ,t10 ,t11 ,t12 ,t13 ,t14 ,t2 ,t15 ,t16 ,t0 ;
	t3 = p1-x1;
	t4 = t3*u1;
	t5 = p2-x2;
	t6 = t5*u2;
	t7 = p3-x3;
	t8 = t7*u3;
	t9 = t4+t6+t8;
	t10 = u1*u1;
	t11 = u2*u2;
	t12 = u3*u3;
	t13 = t10+t11+t12;
	t14 = 1.0/t13;
	t2 = -p1+x1+t9*t14*u1;
	t15 = -p2+x2+t9*t14*u2;
	t16 = -p3+x3+t9*t14*u3;
	t0 = t2*t2+t15*t15+t16*t16;
	return t0;
}
MatrixXd transitionFromCrossMatrix(const MatrixXd& crossM)
{
	MatrixXd result=MatrixXd::Zero(1,3);

	result(0,0)=crossM(2,1);
	result(0,1)=crossM(0,2);
	result(0,2)=crossM(1,0);
	return result;
}

vector<pair<MatrixXd,MatrixXd> > transitionAndRotationFromEssential(const MatrixXd& essential)
{
	JacobiSVD<MatrixXd> svd(essential,Eigen::ComputeFullU |
                                        Eigen::ComputeFullV);


	MatrixXd V=	svd.matrixU();
	MatrixXd U=	svd.matrixV();


	MatrixXd W=MatrixXd::Zero(3,3);

	W<<0, -1, 0, 1, 0, 0, 0, 0, 1;

	MatrixXd Z=MatrixXd::Zero(3,3);

	Z<<0, 1, 0, -1, 0, 0, 0, 0, 0;

	MatrixXd Tx=U*Z*U.transpose();
	MatrixXd R1=U*W*V.transpose();
	if(R1.determinant()<0)
		R1*=-1;

	MatrixXd R2=U*W.transpose()*V.transpose();
	if(R2.determinant()<0)
		R2*=-1;

//	getchar();
	vector<pair<MatrixXd,MatrixXd> > result;

	MatrixXd Ty=transitionFromCrossMatrix(Tx);

	result.push_back(make_pair(Ty,R1));
	result.push_back(make_pair(Ty*-1,R1));
	result.push_back(make_pair(Ty,R2));
	result.push_back(make_pair(Ty*-1,R2));

	return result;

}


MatrixXd bestPointCoefficient(double a,double b,double c,double d,double e,double f)
{
	MatrixXd A0=MatrixXd::Zero(3,4);

	double t2 ,t3 ,t4 ,t5 ,t6 ,t7 ,t8 ,t9 ,t10 ,t11 ,t12 ,t13 ;
	t2 = e*e;
	t3 = f*f;
	t4 = d*d;
	t5 = t2+t3+t4;
	t6 = 1.0/t5;
	t7 = t3*2.0;
	t8 = t4*2.0;
	t9 = t2*2.0;
	t10 = a*d;
	t11 = b*e;
	t12 = c*f;
	t13 = t10+t11+t12;
	A0(0,0) = t6*(t7+t9);
	A0(0,1) = d*e*t6*-2.0;
	A0(0,2) = d*f*t6*-2.0;
	A0(0,3) = -t6*(a*(t2+t3)*2.0-d*(b*e*2.0+c*f*2.0));
	A0(1,0) = d*e*t6*-2.0;
	A0(1,1) = t6*(t7+t8);
	A0(1,2) = e*f*t6*-2.0;
	A0(1,3) = b*-2.0+e*t6*t13*2.0;
	A0(2,0) = d*f*t6*-2.0;
	A0(2,1) = e*f*t6*-2.0;
	A0(2,2) = t6*(t8+t9);
	A0(2,3) = c*-2.0+f*t6*t13*2.0;

	return A0;
}

MatrixXd bestPoint(const MatrixXd& p, const MatrixXd& u)
{
	MatrixXd coefficients= MatrixXd::Zero(3,4);

	assert(p.cols()==3);
	assert(u.cols()==3);
	assert(p.rows()==u.rows());
	assert(p.rows()>1);

	for (int i = 0; i < p.rows(); i++)
	{
		double a,b,c,d,e,f;
		a=p(i,0);
		b=p(i,1);
		c=p(i,2);

		d=u(i,0);
		e=u(i,1);
		f=u(i,2);

		MatrixXd curcoe=bestPointCoefficient(a,b,c,d,e,f);
		coefficients+=curcoe;
	}


	double a,b,c,d,e,f,g,h,i,j,k,l;

	a	=	coefficients(0,0);
	b	=	coefficients(0,1);
	c	=	coefficients(0,2);
	d	=	coefficients(0,3);
	e	=	coefficients(1,0);
	f	=	coefficients(1,1);
	g	=	coefficients(1,2);
	h	=	coefficients(1,3);
	i	=	coefficients(2,0);
	j	=	coefficients(2,1);
	k	=	coefficients(2,2);
	l	=	coefficients(2,3);

	double x,y,z;
	double r=1.0/(a*f*k - a*g*j - b*e*k + b*g*i + c*e*j - c*f*i);
	x=-(b*g*l - b*h*k - c*f*l + c*h*j + d*f*k - d*g*j)*r;
	y=(a*g*l - a*h*k - c*e*l + c*h*i + d*e*k - d*g*i)*r;
	z=-(a*f*l - a*h*j - b*e*l + b*h*i + d*e*j - d*f*i)*r;

	MatrixXd result=MatrixXd(1,3);
	result(0,0)=x;
	result(0,1)=y;
	result(0,2)=z;
//	return result;

}

double length(const MatrixXd& a)
{
	return sqrt((a*a.transpose())(0,0));
}

double angleBetween(const MatrixXd& a,const MatrixXd& b)
{
	return acos((a*b.transpose())(0,0)/length(a)/length(b));
}



bool pointBeforeCamera(const MatrixXd&x,const MatrixXd&p,const MatrixXd& u)
{
	MatrixXd tx=x-p;

	auto core=[](double a,double b)->bool
	{
		return (a*b)>0 && abs(a)>abs(b);
	};
	return core(tx(0,0),u(0,0)) && core(tx(0,1),u(0,1)) && core(tx(0,2),u(0,2));
}

pair<MatrixXd,vector<double> > bestPoints(const MatrixXd& spnts1,const vector<int>& ind1,const MatrixXd& spnts2,const vector<int>& ind2,const vector<MatrixXd>& transitions,const vector<MatrixXd>& rotations)
{
	assert(ind1.size()==ind2.size());

	MatrixXd points(ind1.size(),3);
	vector<double> error(ind1.size(),-1);

	MatrixXd curP(2,3);
	curP.row(0)=transitions[0];
	curP.row(1)=transitions[1];
	for (int i = 0; i < ind1.size(); i++)
	{
		
		MatrixXd curU(2,3);
		curU.row(0)=(rotations[0]*spnts1.row(ind1[i]).transpose()).transpose();		
		curU.row(1)=(rotations[1]*spnts2.row(ind2[i]).transpose()).transpose();
		points.row(i)=bestPoint(curP,curU);

	/*	if(angleBetween(points.row(i)-curP.row(0),curU.row(0))<constrain_on_goodPoint && angleBetween(points.row(i)-curP.row(1),curU.row(1))<constrain_on_goodPoint )
		{
	
			goodlabel[i]=true;
		}*/

		if(pointBeforeCamera(points.row(i),curP.row(0),curU.row(0)) && pointBeforeCamera(points.row(i),curP.row(1),curU.row(1)))
		{
	
			//goodlabel[i]=true;
			error[i]=distanceBetweenPointLine(points.row(i),curP.row(0),curU.row(0))+distanceBetweenPointLine(points.row(i),curP.row(1),curU.row(1));
		}
	}

	return make_pair(points,error);
}




auto geometricReconstructionFrom2Frames(const MatrixXd& sphericalPoints1,const vector<int>& ind1,const MatrixXd& sphericalPoints2,const vector<int>& ind2,MatrixXd& transition,MatrixXd& rotation)->pair<MatrixXd,vector<double> >
//auto geometricReconstructionFrom2Frames(const MatrixXd& sphericalPoints1,const MatrixXd& sphericalPoints2,vector<double>& transition,vector<double >& rotation)->pair<MatrixXd ,vector<double> >
{
//	vector<vector<double> > points(sphericalPoints1.rows(),vector<double>(3,0.0));
//	vector<double> errors(sphericalPoints1.rows(),-1);

	assert(ind1.size()==ind2.size());
	MatrixXd observation(ind1.size(),9);

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
	
	for (int i = 0; i < observation.rows(); i++)
	{

		MatrixXd tob=sphericalPoints1.row(ind1[i]).transpose()*sphericalPoints2.row(ind2[i]);

		observation.row(i)=convert(tob);
	}

	//cout<<observation<<endl<<endl;

	JacobiSVD<decltype(observation)> svd(observation, Eigen::ComputeFullU |
                                        Eigen::ComputeFullV);
	//std::cout << "matrix V" << endl << svd.matrixV() <<"\n"<< endl;

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
	
	vector<pair<MatrixXd,vector<double> > > bestPointCandidates(transtionAndRotations.size());

	auto countDouble=[](const vector<double>& toc)->int
	{
		int sum=0;
		for(auto s:toc)
			if(s>0)
				++sum;

		return sum;
	};

	int bestIndex;
	int bestCount;
	vector<int> goodPointCount(transtionAndRotations.size(),0);

	//cout<<"transition:\n"<<endl;
	for (int i = 0; i < transtionAndRotations.size(); i++)
	{
		//cout<<transtionAndRotations[i].first<<endl;
		vector<MatrixXd> trans(2,MatrixXd::Zero(1,3));
		trans[1]=transtionAndRotations[i].first;
		vector<MatrixXd> rots(2,MatrixXd::Identity(3,3));
		rots[1]=transtionAndRotations[i].second;
		bestPointCandidates[i]=bestPoints(sphericalPoints1,ind1,sphericalPoints2,ind2,trans,rots);
		goodPointCount[i]=countDouble(bestPointCandidates[i].second);

		if(i==0)
		{
			bestIndex=i;
			bestCount=goodPointCount[i];
		}
		else
		{
			if(goodPointCount[i]>bestCount)
			{
				bestIndex=i;
				bestCount=goodPointCount[i];
			}
		}

	}
	transition=transtionAndRotations[bestIndex].first;
	rotation=rotationThomasonPara (transtionAndRotations[bestIndex].second.inverse());

//	bestIndex=3;
 //	cout<<"the best index:"<<bestIndex<<endl;
//	cout<<"the best transition:\n"<<transtionAndRotations[bestIndex].first <<endl;
//	cout<<"the best rotation:\n"<<transtionAndRotations[bestIndex].second<<endl;

	


	//cout<<observation;
	return bestPointCandidates[bestIndex];

}




auto threeDimensionReconstruction(const string& featureFileName,const string& matchFileName)->tuple<MatrixXd,MatrixXd,MatrixXd>
{
	
	vector<vector<vector<int> > > features;
	vector<map<int,int>> correspondences;
	vector<unordered_set<int> > contain;
	vector<vector<int> > featureIsPoint;

	auto trajectories=collectFromPairwiseMatching(featureFileName,matchFileName,features,correspondences,contain,featureIsPoint);


	MatrixXd transitions(features.size(),3);

	MatrixXd rotations(features.size(),3);

	vector<bool> alreadyEstimated(features.size(),false);

	transitions.row(0)=MatrixXd::Zero(1,3);
	rotations.row(0)=MatrixXd::Zero(1,3);
	alreadyEstimated[0]=true;

	MatrixXd reconstructedPoints(trajectories.size(),3);

	vector<bool> alreadyReconstructed(trajectories.size(),false);

	vector<int> imageSize(2);
	imageSize[0]=512;imageSize[1]=256;

	auto convertFeature=[&](const vector<vector<int>> &fea)->MatrixXd
	{
		MatrixXd sfea(fea.size(),3);
		for (int i = 0; i < fea.size(); i++)
		{
			sfea.row(i)=imageCordinate2Phere(fea[i],imageSize);
		}
		return sfea;
	};
	vector<MatrixXd> sphericalFeatures(features.size());
	
	for (int i = 0; i < sphericalFeatures.size(); i++)
	{
		sphericalFeatures[i]=convertFeature(features[i]);
	}

//	vector<vector<int> > points1(correspondences[0].size()),points2(correspondences[0].size());

	vector<int> index(correspondences[0].size());

	vector<int> ind1,ind2;

	int count=0;
	for ( auto&k:correspondences[0])
	{
		ind1.push_back(k.first);
		ind2.push_back(k.second);
		index[count]=featureIsPoint[0][k.first];
		++count;
	}

	MatrixXd curTransition,curRotation;
	auto pointsError2Frame=geometricReconstructionFrom2Frames(sphericalFeatures[0],ind1,sphericalFeatures[1],ind2,curTransition,curRotation);
	transitions.row(1)=curTransition;
	rotations.row(1)=curRotation;
	alreadyEstimated[1]=true;


	auto setReconstructedPoints=[&reconstructedPoints,&alreadyReconstructed](const pair<MatrixXd,vector<double> >& mpointsError2Frame,const vector<int>& mindex)
	{
		assert(mpointsError2Frame.second.size()==mindex.size());


		for (int i = 0; i < mpointsError2Frame.second.size(); i++)
		{
			if (mpointsError2Frame.second[i]>0 && mpointsError2Frame.second[i]<constrain_on_goodPoint)
			{
				reconstructedPoints.row(mindex[i])=mpointsError2Frame.first.row(i);
				alreadyReconstructed[mindex[i]]=true;

			}
			else
			{
				alreadyReconstructed[mindex[i]]=false;
			}

		}
	};

	setReconstructedPoints(pointsError2Frame,index);

	auto reconstructPoints=[&](int lcameraIndex)
	{
		assert(lcameraIndex>=1);
		assert(alreadyEstimated[lcameraIndex-1] && alreadyEstimated[lcameraIndex]);
		int lCount=0;
		vector<int> lIndex,lIndex1,lIndex2;
		for (auto& s:correspondences[lcameraIndex-1])
		{
//			if(!alreadyReconstructed[featureIsPoint[lcameraIndex-1][s.first]])
			{
				lIndex.push_back(featureIsPoint[lcameraIndex-1][s.first]);
				lIndex1.push_back(s.first);
				lIndex2.push_back(s.second);
			}
			++lCount;
		}
		vector<MatrixXd> ltransitions(2);
		vector<MatrixXd> lrotations(2);
		ltransitions[0]=transitions.row(lcameraIndex-1);
		ltransitions[1]=transitions.row(lcameraIndex);

		lrotations[0]=rotationThomason(rotations.row(lcameraIndex-1)).inverse();
		lrotations[1]=rotationThomason(rotations.row(lcameraIndex)).inverse();
		auto lpointsError2Frame=bestPoints(sphericalFeatures[lcameraIndex-1],lIndex1,sphericalFeatures[lcameraIndex],lIndex2,ltransitions,lrotations);
		setReconstructedPoints(lpointsError2Frame,lIndex);

	};
	vector<funcType2> functions(2);
	functions[0]=&functionForRotationAndTransition;
	functions[1]=&functionForRotationAndTransitionUnitLength;
	vector<funcType2> jacabianFunctions(4);
	jacabianFunctions[0]=&jacobianForRotationAndTransition;
	jacabianFunctions[1]=&jacobianForRotationAndTransitionUnitLength;
	jacabianFunctions[2]=&jacobianForPoint;
	jacabianFunctions[3]=&jacobianForPointUnitLength;

	auto bundleAdjustment=[&](unordered_map<int,cameraType>& lbundlePara)
	{
		unordered_set<int> pntIndx;

		unordered_map<int,pair<int,int>> cameraParaLookUp;
		unordered_map<int,pair<int,int>> pointParaLookUp;
		unordered_map<int,int> staticCameraParaLookUp;

		vector<vector<int> > funcDataMap;
		vector<vector<int> > jfuncDataMap;
		int totalProjCount=0;
		int totalParameterCount=0;
		vector<pair<int,int> > allProjs;

		for (auto& s:lbundlePara)
		{
			assert(alreadyEstimated[s.first]);
			if (s.second!=cameraType::_static)
			{
				for (auto&w: contain[s.first])
				{
					if(alreadyReconstructed[w] && !pntIndx.count(w))
					{
						int prjCount=0;

						unordered_set<int> curPrj;
						for (int i = 0; i < trajectories[w].first.size(); i++)
						{
							if(lbundlePara.count(trajectories[w].first[i]))
							{
								++prjCount;
								curPrj.insert(i);
							}
						}

						if(prjCount>1)
						{
							if(!pointParaLookUp.count(w))
							{
								pointParaLookUp[w]=make_pair(totalParameterCount,3);
								totalParameterCount+=3;
							}
							for (auto& x:curPrj)
							{
								int curPrjCamera=trajectories[w].first[x];

								



								switch (lbundlePara[curPrjCamera])
								{
								case cameraType::_static:
									{
										if(!staticCameraParaLookUp.count(curPrjCamera))
										{
											int _t=staticCameraParaLookUp.size();
											staticCameraParaLookUp[curPrjCamera]=_t;
										}
										
										vector<int> static_funcMap(7);
										vector<int> static_jfuncMap(8);
										static_funcMap[0]=0;
										static_funcMap[1]=totalProjCount*3;
										static_funcMap[2]=totalProjCount;
										static_funcMap[3]=1;
										static_funcMap[4]=staticCameraParaLookUp[curPrjCamera];
										static_funcMap[5]=pointParaLookUp[w].first;
										static_funcMap[6]=pointParaLookUp[w].second;
										funcDataMap.push_back(static_funcMap);
										
										static_jfuncMap[0]=2;
										static_jfuncMap[1]=totalProjCount*3;
										static_jfuncMap[2]=pointParaLookUp[w].first;
										static_jfuncMap[3]=totalProjCount;
										static_jfuncMap[4]=1;
										static_jfuncMap[5]=staticCameraParaLookUp[curPrjCamera];
										static_jfuncMap[6]=pointParaLookUp[w].first;
										static_jfuncMap[7]=pointParaLookUp[w].second;
										jfuncDataMap.push_back(static_jfuncMap);
									}
									break;
								case cameraType::_unitLength:
									{
										if(!cameraParaLookUp.count(curPrjCamera))
										{
											cameraParaLookUp[curPrjCamera]=make_pair(totalParameterCount,5);
											totalParameterCount+=5;
										}
										vector<int> funcMap(9);	
										vector<int> jfuncMap(10);
										funcMap[0]=1;
										funcMap[1]=totalProjCount*3;
										funcMap[2]=totalProjCount;
										funcMap[3]=-1;
										funcMap[4]=-1;
										funcMap[5]=cameraParaLookUp[curPrjCamera].first;
										funcMap[6]=cameraParaLookUp[curPrjCamera].second;
										funcMap[7]=pointParaLookUp[w].first;
										funcMap[8]=pointParaLookUp[w].second;
										funcDataMap.push_back(funcMap);
										
										jfuncMap[0]=1;
										jfuncMap[1]=totalProjCount*3;
										jfuncMap[2]=cameraParaLookUp[curPrjCamera].first;
										jfuncMap[3]=totalProjCount;
										jfuncMap[4]=-1;
										jfuncMap[5]=-1;
										jfuncMap[6]=cameraParaLookUp[curPrjCamera].first;
										jfuncMap[7]=cameraParaLookUp[curPrjCamera].second;
										jfuncMap[8]=pointParaLookUp[w].first;
										jfuncMap[9]=pointParaLookUp[w].second;
										jfuncDataMap.push_back(jfuncMap);

										jfuncMap[0]=3;
									
										jfuncMap[2]=pointParaLookUp[w].first;
									
										jfuncDataMap.push_back(jfuncMap);

									}
									break;
								case cameraType::_ordinary:
									{
										if(!cameraParaLookUp.count(curPrjCamera))
										{
											cameraParaLookUp[curPrjCamera]=make_pair(totalParameterCount,6);
											totalParameterCount+=6;
										}
										vector<int> funcMap(9);	
										vector<int> jfuncMap(10);
										funcMap[0]=0;
										funcMap[1]=totalProjCount*3;
										funcMap[2]=totalProjCount;
										funcMap[3]=-1;
										funcMap[4]=-1;
										funcMap[5]=cameraParaLookUp[curPrjCamera].first;
										funcMap[6]=cameraParaLookUp[curPrjCamera].second;
										funcMap[7]=pointParaLookUp[w].first;
										funcMap[8]=pointParaLookUp[w].second;
										funcDataMap.push_back(funcMap);
										
										jfuncMap[0]=0;
										jfuncMap[1]=totalProjCount*3;
										jfuncMap[2]=cameraParaLookUp[curPrjCamera].first;
										jfuncMap[3]=totalProjCount;
										jfuncMap[4]=-1;
										jfuncMap[5]=-1;
										jfuncMap[6]=cameraParaLookUp[curPrjCamera].first;
										jfuncMap[7]=cameraParaLookUp[curPrjCamera].second;
										jfuncMap[8]=pointParaLookUp[w].first;
										jfuncMap[9]=pointParaLookUp[w].second;
										jfuncDataMap.push_back(jfuncMap);

										jfuncMap[0]=2;
									
										jfuncMap[2]=pointParaLookUp[w].first;
									
										jfuncDataMap.push_back(jfuncMap);
									}
									break;
								default:
									break;
								}
								++totalProjCount;
								allProjs.push_back(make_pair( trajectories[w].first[x],trajectories[w].second[x]));
							}
							pntIndx.insert(w);
						}
					}
				}
			}
		}


		MatrixXd assistantPara(staticCameraParaLookUp.size(),6);

		MatrixXd dataSet(totalProjCount,3);

		MatrixXd obj_vals=MatrixXd::Zero(1,totalProjCount*3);
		MatrixXd intit_parameters(1,totalParameterCount);

		for ( auto& sc: staticCameraParaLookUp)
		{
			assistantPara.block(sc.second,0,1,3)=rotations.row(sc.first);
			assistantPara.block(sc.second,3,1,3)=transitions.row(sc.first);
		}
		for (int i=0;i<allProjs.size();++i)
		{
			dataSet.row(i)=sphericalFeatures[allProjs[i].first].row(allProjs[i].second);
		}
		for (auto& cc:cameraParaLookUp)
		{
			
			MatrixXd paraM;
			if(lbundlePara[cc.first]==cameraType::_unitLength)
			{
				paraM.resize(1,5);
				paraM.block(0,3,1,2)=transition2Para(transitions.row(cc.first));
			}

			if(lbundlePara[cc.first]==cameraType::_ordinary)
			{
				paraM.resize(1,6);
				paraM.block(0,3,1,3)=transitions.row(cc.first);
			}
			paraM.block(0,0,1,3)=rotations.row(cc.first);
			intit_parameters.block(0,cc.second.first,1,cc.second.second)=paraM;
		}

		for(auto& p:pointParaLookUp)
		{

			intit_parameters.block(0,p.second.first,1,p.second.second)=reconstructedPoints.row(p.first);
		}
		auto est_parameters=levenbergM_advanced(dataSet,assistantPara,funcDataMap,jfuncDataMap,obj_vals,functions,jacabianFunctions,intit_parameters);


		for (auto& cc:cameraParaLookUp)
		{
			
			MatrixXd paraM=est_parameters.block(0,cc.second.first,1,cc.second.second);
			if(lbundlePara[cc.first]==cameraType::_unitLength)
			{
				
				transitions.row(cc.first)=transitionFrom2Para(paraM.block(0,3,1,2));
			}

			if(lbundlePara[cc.first]==cameraType::_ordinary)
			{
				
				transitions.row(cc.first)=paraM.block(0,3,1,3);
			}
			rotations.row(cc.first)=paraM.block(0,0,1,3);
			//=paraM;
		}

		for(auto& p:pointParaLookUp)
		{

			reconstructedPoints.row(p.first)=est_parameters.block(0,p.second.first,1,p.second.second);
		}


	};

	unordered_map<int,cameraType> bundlePara1;

	bundlePara1[0]=cameraType::_static;
	bundlePara1[1]=cameraType::_unitLength;

	cout<<"bundle adjustment for cameras and points before camera number "<<"1"<<endl;
	bundleAdjustment(bundlePara1);
	reconstructPoints(1);

	for (int cameraIndex = 2; cameraIndex < features.size(); cameraIndex++)
	{
		vector<int> curInd1,curInd2;
		for (int i = 0; i < featureIsPoint[cameraIndex].size(); i++)
		{
			if (featureIsPoint[cameraIndex][i]>=0 && alreadyReconstructed[featureIsPoint[cameraIndex][i]])
			{
				curInd1.push_back(i);
				curInd2.push_back(featureIsPoint[cameraIndex][i]);
			}
		}
		cout<<curInd1.size()<<" projections matching "<<curInd2.size()<<" points for camera position estimation"<<"\n\n"<<endl;

		if(cameraIndex==12)
		{
			cout<<"break point";
		}
		auto cameraPara=estimateCameraParameter(sphericalFeatures[cameraIndex],curInd1,reconstructedPoints,curInd2);

		cout<<"estimated camera parameters of camera number "<<cameraIndex<<": "<<cameraPara<<endl;

		rotations.row(cameraIndex)=cameraPara.block(0,0,1,3);
		transitions.row(cameraIndex)=cameraPara.block(0,3,1,3);
		alreadyEstimated[cameraIndex]=true;
		reconstructPoints(cameraIndex);




		unordered_map<int,cameraType> bundlePara;

		bundlePara[0]=cameraType::_static;
		bundlePara[1]=cameraType::_unitLength;
		for (int i = 2; i <=cameraIndex; i++)
		{
			bundlePara[i]=cameraType::_ordinary;
		}
		cout<<"bundle adjustment for cameras and points before camera number "<<cameraIndex<<endl;
		bundleAdjustment(bundlePara);
		reconstructPoints(cameraIndex);
	}


	return make_tuple(transitions,rotations,reconstructedPoints);
	
}

MatrixXd functionForRotationAndTransition(const MatrixXd& parameters,const MatrixXd& variables)
{
	double a,b,c,d,e,f,t1,t2,t3,r1,r2,r3;//% a b c coordinates of projpoints, d e f coordinates of points r1 r2 r3, rotation parameters t1 t2 t3 transition parameters
	a=parameters(0,0) ;
	b=parameters(0,1) ;
	c=parameters(0,2) ;
	d=parameters(0,3) ;
	e=parameters(0,4) ;
	f=parameters(0,5) ;

	r1=variables(0,0) ;
	r2=variables(0,1) ;
	r3=variables(0,2) ;
	t1=variables(0,3) ;
	t2=variables(0,4) ;
	t3=variables(0,5) ;
	MatrixXd A0(1,3);
	double t5 ,t6 ,t7 ,t8 ,t9 ,t10 ,t11 ,t12 ,t13 ,t14 ,t15 ,t16 ,t17 ;
	t5 = r1*r1;
	t6 = r2*r2;
	t7 = r3*r3;
	t8 = 1.0/b;
	t9 = d*d;
	t10 = e*e;
	t11 = f*f;
	t12 = t9+t10+t11;
	t13 = 1.0/sqrt(t12);
	t14 = t5+t6+t7+4.0;
	t15 = 1.0/t14;
	t16 = 1.0/a;
	t17 = 1.0/c;
	A0(0,0) = -t8*t13*t15*t16*(a*e*4.0-b*d*4.0-a*t2*4.0+b*t1*4.0+a*d*r3*4.0-a*f*r1*4.0+b*e*r3*4.0-b*f*r2*4.0-a*e*t5-b*d*t5+a*e*t6+b*d*t6-a*e*t7+b*d*t7+a*r1*t3*4.0-a*r3*t1*4.0+b*r2*t3*4.0-b*r3*t2*4.0+a*t2*t5-a*t2*t6+a*t2*t7+b*t1*t5-b*t1*t6-b*t1*t7+a*d*r1*r2*2.0-b*e*r1*r2*2.0+a*f*r2*r3*2.0-b*f*r1*r3*2.0-a*r1*r2*t1*2.0-a*r2*r3*t3*2.0+b*r1*r2*t2*2.0+b*r1*r3*t3*2.0);
	A0(0,1) = -t8*t13*t15*t17*(b*f*4.0-c*e*4.0-b*t3*4.0+c*t2*4.0-b*d*r2*4.0+b*e*r1*4.0-c*d*r3*4.0+c*f*r1*4.0-b*f*t5+c*e*t5-b*f*t6-c*e*t6+b*f*t7+c*e*t7-b*r1*t2*4.0+b*r2*t1*4.0-c*r1*t3*4.0+c*r3*t1*4.0+b*t3*t5+b*t3*t6-b*t3*t7-c*t2*t5+c*t2*t6-c*t2*t7+b*d*r1*r3*2.0-c*d*r1*r2*2.0+b*e*r2*r3*2.0-c*f*r2*r3*2.0-b*r1*r3*t1*2.0-b*r2*r3*t2*2.0+c*r1*r2*t1*2.0+c*r2*r3*t3*2.0);
	A0(0,2) = t13*t15*t16*t17*(a*f*4.0-c*d*4.0-a*t3*4.0+c*t1*4.0-a*d*r2*4.0+a*e*r1*4.0+c*e*r3*4.0-c*f*r2*4.0-a*f*t5-c*d*t5-a*f*t6+c*d*t6+a*f*t7+c*d*t7-a*r1*t2*4.0+a*r2*t1*4.0+c*r2*t3*4.0-c*r3*t2*4.0+a*t3*t5+a*t3*t6-a*t3*t7+c*t1*t5-c*t1*t6-c*t1*t7+a*d*r1*r3*2.0+a*e*r2*r3*2.0-c*e*r1*r2*2.0-c*f*r1*r3*2.0-a*r1*r3*t1*2.0-a*r2*r3*t2*2.0+c*r1*r2*t2*2.0+c*r1*r3*t3*2.0);

	return A0;
}


MatrixXd functionForRotationAndTransition_(const MatrixXd& parameters,const MatrixXd& variables)
{
	double a,b,c,d,e,f,t1,t2,t3,r1,r2,r3;//% a b c coordinates of projpoints, d e f coordinates of points r1 r2 r3, rotation parameters t1 t2 t3 transition parameters
	a=parameters(0,0) ;
	b=parameters(0,1) ;
	c=parameters(0,2) ;
	d=parameters(0,3) ;
	e=parameters(0,4) ;
	f=parameters(0,5) ;

	r1=variables(0,0) ;
	r2=variables(0,1) ;
	r3=variables(0,2) ;
	t1=variables(0,3) ;
	t2=variables(0,4) ;
	t3=variables(0,5) ;

	double t5 ,t6 ,t7 ,t8 ,t9 ,t10 ,t11 ,t12 ;
	t5 = r1*r1;
	t6 = r2*r2;
	t7 = r3*r3;
	t8 = 1.0/b;
	t9 = t5+t6+t7+4.0;
	t10 = 1.0/t9;
	t11 = 1.0/a;
	t12 = 1.0/c;
	MatrixXd A0(1,3);
	A0(0,0) = -t8*t10*t11*(a*e*4.0-b*d*4.0-a*t2*4.0+b*t1*4.0+a*d*r3*4.0-a*f*r1*4.0+b*e*r3*4.0-b*f*r2*4.0-a*e*t5-b*d*t5+a*e*t6+b*d*t6-a*e*t7+b*d*t7+a*r1*t3*4.0-a*r3*t1*4.0+b*r2*t3*4.0-b*r3*t2*4.0+a*t2*t5-a*t2*t6+a*t2*t7+b*t1*t5-b*t1*t6-b*t1*t7+a*d*r1*r2*2.0-b*e*r1*r2*2.0+a*f*r2*r3*2.0-b*f*r1*r3*2.0-a*r1*r2*t1*2.0-a*r2*r3*t3*2.0+b*r1*r2*t2*2.0+b*r1*r3*t3*2.0);
	A0(0,1) = -t8*t10*t12*(b*f*4.0-c*e*4.0-b*t3*4.0+c*t2*4.0-b*d*r2*4.0+b*e*r1*4.0-c*d*r3*4.0+c*f*r1*4.0-b*f*t5+c*e*t5-b*f*t6-c*e*t6+b*f*t7+c*e*t7-b*r1*t2*4.0+b*r2*t1*4.0-c*r1*t3*4.0+c*r3*t1*4.0+b*t3*t5+b*t3*t6-b*t3*t7-c*t2*t5+c*t2*t6-c*t2*t7+b*d*r1*r3*2.0-c*d*r1*r2*2.0+b*e*r2*r3*2.0-c*f*r2*r3*2.0-b*r1*r3*t1*2.0-b*r2*r3*t2*2.0+c*r1*r2*t1*2.0+c*r2*r3*t3*2.0);
	A0(0,2) = t10*t11*t12*(a*f*4.0-c*d*4.0-a*t3*4.0+c*t1*4.0-a*d*r2*4.0+a*e*r1*4.0+c*e*r3*4.0-c*f*r2*4.0-a*f*t5-c*d*t5-a*f*t6+c*d*t6+a*f*t7+c*d*t7-a*r1*t2*4.0+a*r2*t1*4.0+c*r2*t3*4.0-c*r3*t2*4.0+a*t3*t5+a*t3*t6-a*t3*t7+c*t1*t5-c*t1*t6-c*t1*t7+a*d*r1*r3*2.0+a*e*r2*r3*2.0-c*e*r1*r2*2.0-c*f*r1*r3*2.0-a*r1*r3*t1*2.0-a*r2*r3*t2*2.0+c*r1*r2*t2*2.0+c*r1*r3*t3*2.0);
	return A0;
}



MatrixXd functionForRotationAndTransitionUnitLength_(const MatrixXd& parameters,const MatrixXd& variables2)
{
	MatrixXd variable(1,6);
	variable.block(0,0,1,3)=variables2.block(0,0,1,3);
	variable.block(0,3,1,3)=transitionFrom2Para(variables2.block(0,3,1,2));

	return functionForRotationAndTransition(parameters,variable);
}


MatrixXd functionForRotationAndTransitionUnitLength(const MatrixXd& parameters,const MatrixXd& variables2)
{
	MatrixXd variable(1,6);
	variable.block(0,0,1,3)=variables2.block(0,0,1,3);
	variable.block(0,3,1,3)=transitionFrom2Para(variables2.block(0,3,1,2));

	return functionForRotationAndTransition(parameters,variable);
}

MatrixXd jacobianForPointUnitLength(const MatrixXd& parameters,const MatrixXd& variables2)
{
	MatrixXd variable(1,6);
	variable.block(0,0,1,3)=variables2.block(0,0,1,3);
	variable.block(0,3,1,3)=transitionFrom2Para(variables2.block(0,3,1,2));

	return jacobianForPoint(parameters,variable);
}

MatrixXd jacobianForPoint(const MatrixXd& parameters,const MatrixXd& variables)
{

	double a,b,c,d,e,f,t1,t2,t3,r1,r2,r3;
	a=parameters(0,0) ;
	b=parameters(0,1) ;
	c=parameters(0,2) ;
	d=parameters(0,3) ;
	e=parameters(0,4) ;
	f=parameters(0,5) ;

	r1=variables(0,0) ;
	r2=variables(0,1) ;
	r3=variables(0,2) ;
	t1=variables(0,3) ;
	t2=variables(0,4) ;
	t3=variables(0,5) ;
	MatrixXd A0(3,3);
	double t5 ,t6 ,t7 ,t8 ,t9 ,t10 ,t11 ,t12 ,t13 ,t14 ,t15 ,t16 ,t17 ,t18 ,t19 ,t20 ,t21 ,t22 ,t23 ,t24 ,t25 ,t26 ,t27 ,t28 ,t29 ,t30 ,t31 ,t32 ,t33 ,t35 ,t36 ,t37 ,t38 ,t39 ,t40 ,t41 ,t42 ,t43 ,t44 ,t45 ,t46 ,t47 ,t48 ,t49 ,t50 ,t34 ,t51 ,t52 ,t53 ,t54 ,t55 ,t56 ,t57 ,t58 ,t59 ,t60 ,t61 ,t62 ,t63 ,t64 ,t65 ,t66 ,t67 ,t68 ,t69 ,t73 ,t74 ,t75 ,t76 ,t77 ,t78 ,t79 ,t80 ,t81 ,t82 ,t83 ,t84 ,t85 ,t86 ,t87 ,t88 ,t70 ,t71 ,t72 ,t89 ,t90 ,t91 ,t92 ,t93 ,t94 ,t95 ,t96 ,t97 ,t98 ,t99 ,t100 ,t101 ,t102 ,t103 ,t104 ,t105 ,t106 ,t107 ,t108 ,t109 ,t113 ,t114 ,t115 ,t116 ,t117 ,t118 ,t119 ,t120 ,t121 ,t122 ,t123 ,t124 ,t125 ,t126 ,t127 ,t128 ,t110 ,t111 ,t112 ;
	t5 = r1*r1;
	t6 = r2*r2;
	t7 = r3*r3;
	t8 = 1.0/a;
	t9 = 1.0/b;
	t10 = d*d;
	t11 = e*e;
	t12 = f*f;
	t13 = t10+t11+t12;
	t14 = t5+t6+t7+4.0;
	t15 = 1.0/t14;
	t16 = 1.0/sqrt(t13);
	t17 = 1.0/pow(t13,3.0/2.0);
	t18 = a*e*4.0;
	t19 = b*t1*4.0;
	t20 = a*t2*t5;
	t21 = a*t2*t7;
	t22 = b*t1*t5;
	t23 = a*d*r3*4.0;
	t24 = b*e*r3*4.0;
	t25 = a*r1*t3*4.0;
	t26 = b*r2*t3*4.0;
	t27 = a*e*t6;
	t28 = b*d*t6;
	t29 = b*d*t7;
	t30 = a*d*r1*r2*2.0;
	t31 = a*f*r2*r3*2.0;
	t32 = b*r1*r2*t2*2.0;
	t33 = b*r1*r3*t3*2.0;
	t35 = b*d*4.0;
	t36 = a*t2*4.0;
	t37 = a*t2*t6;
	t38 = b*t1*t6;
	t39 = b*t1*t7;
	t40 = a*f*r1*4.0;
	t41 = b*f*r2*4.0;
	t42 = a*r3*t1*4.0;
	t43 = b*r3*t2*4.0;
	t44 = a*e*t5;
	t45 = b*d*t5;
	t46 = a*e*t7;
	t47 = b*e*r1*r2*2.0;
	t48 = b*f*r1*r3*2.0;
	t49 = a*r1*r2*t1*2.0;
	t50 = a*r2*r3*t3*2.0;
	t34 = t18+t19+t20+t21+t22+t23+t24+t25+t26+t27+t28+t29+t30+t31+t32+t33-t35-t36-t37-t38-t39-t40-t41-t42-t43-t44-t45-t46-t47-t48-t49-t50;
	t51 = b*r2*4.0;
	t52 = b*r1*r3*2.0;
	t53 = 1.0/c;
	t54 = b*f*4.0;
	t55 = c*t2*4.0;
	t56 = b*t3*t5;
	t57 = b*t3*t6;
	t58 = c*t2*t6;
	t59 = b*e*r1*4.0;
	t60 = c*f*r1*4.0;
	t61 = b*r2*t1*4.0;
	t62 = c*r3*t1*4.0;
	t63 = c*e*t5;
	t64 = b*f*t7;
	t65 = c*e*t7;
	t66 = b*d*r1*r3*2.0;
	t67 = b*e*r2*r3*2.0;
	t68 = c*r1*r2*t1*2.0;
	t69 = c*r2*r3*t3*2.0;
	t73 = c*e*4.0;
	t74 = b*t3*4.0;
	t75 = b*t3*t7;
	t76 = c*t2*t5;
	t77 = c*t2*t7;
	t78 = b*d*r2*4.0;
	t79 = c*d*r3*4.0;
	t80 = b*r1*t2*4.0;
	t81 = c*r1*t3*4.0;
	t82 = b*f*t5;
	t83 = b*f*t6;
	t84 = c*e*t6;
	t85 = c*d*r1*r2*2.0;
	t86 = c*f*r2*r3*2.0;
	t87 = b*r1*r3*t1*2.0;
	t88 = b*r2*r3*t2*2.0;
	t70 = t54+t55+t56+t57+t58+t59+t60+t61+t62+t63+t64+t65+t66+t67+t68+t69-t73-t74-t75-t76-t77-t78-t79-t80-t81-t82-t83-t84-t85-t86-t87-t88;
	t71 = b*t6;
	t72 = b*t7;
	t89 = c*t5;
	t90 = c*t7;
	t91 = a*r1*4.0;
	t92 = c*r3*4.0;
	t93 = c*r1*r2*2.0;
	t94 = a*f*4.0;
	t95 = c*t1*4.0;
	t96 = a*t3*t5;
	t97 = a*t3*t6;
	t98 = c*t1*t5;
	t99 = a*e*r1*4.0;
	t100 = c*e*r3*4.0;
	t101 = a*r2*t1*4.0;
	t102 = c*r2*t3*4.0;
	t103 = c*d*t6;
	t104 = a*f*t7;
	t105 = c*d*t7;
	t106 = a*d*r1*r3*2.0;
	t107 = a*e*r2*r3*2.0;
	t108 = c*r1*r2*t2*2.0;
	t109 = c*r1*r3*t3*2.0;
	t113 = c*d*4.0;
	t114 = a*t3*4.0;
	t115 = a*t3*t7;
	t116 = c*t1*t6;
	t117 = c*t1*t7;
	t118 = a*d*r2*4.0;
	t119 = c*f*r2*4.0;
	t120 = a*r1*t2*4.0;
	t121 = c*r3*t2*4.0;
	t122 = a*f*t5;
	t123 = c*d*t5;
	t124 = a*f*t6;
	t125 = c*e*r1*r2*2.0;
	t126 = c*f*r1*r3*2.0;
	t127 = a*r1*r3*t1*2.0;
	t128 = a*r2*r3*t2*2.0;
	t110 = t94+t95+t96+t97+t98+t99+t100+t101+t102+t103+t104+t105+t106+t107+t108+t109-t113-t114-t115-t116-t117-t118-t119-t120-t121-t122-t123-t124-t125-t126-t127-t128;
	t111 = a*4.0;
	t112 = a*t6;
	A0(0,0) = -t8*t9*t15*t16*(b*-4.0+t71+t72+a*r3*4.0-b*t5+a*r1*r2*2.0)+d*t8*t9*t15*t17*t34;
	A0(0,1) = -t8*t9*t15*t16*(t111+t112+b*r3*4.0-a*t5-a*t7-b*r1*r2*2.0)+e*t8*t9*t15*t17*t34;
	A0(0,2) = t8*t9*t15*t16*(t51+t52+t91-a*r2*r3*2.0)+f*t8*t9*t15*t17*t34;
	A0(1,0) = t9*t15*t16*t53*(t51-t52+t92+t93)+d*t9*t15*t17*t53*t70;
	A0(1,1) = -t9*t15*t16*t53*(c*-4.0+t89+t90+b*r1*4.0-c*t6+b*r2*r3*2.0)+e*t9*t15*t17*t53*t70;
	A0(1,2) = -t9*t15*t16*t53*(b*4.0-t71+t72+c*r1*4.0-b*t5-c*r2*r3*2.0)+f*t9*t15*t17*t53*t70;
	A0(2,0) = -t8*t15*t16*t53*(c*4.0+t89-t90+a*r2*4.0-c*t6-a*r1*r3*2.0)-d*t8*t15*t17*t53*t110;
	A0(2,1) = t8*t15*t16*t53*(t91+t92-t93+a*r2*r3*2.0)-e*t8*t15*t17*t53*t110;
	A0(2,2) = -t8*t15*t16*t53*(-t111+t112+c*r2*4.0+a*t5-a*t7+c*r1*r3*2.0)-f*t8*t15*t17*t53*t110;

	return A0;
}
MatrixXd jacobianForPoint_(const MatrixXd& parameters,const MatrixXd& variables)
{

	double a,b,c,d,e,f,t1,t2,t3,r1,r2,r3;
	a=parameters(0,0) ;
	b=parameters(0,1) ;
	c=parameters(0,2) ;
	d=parameters(0,3) ;
	e=parameters(0,4) ;
	f=parameters(0,5) ;

	r1=variables(0,0) ;
	r2=variables(0,1) ;
	r3=variables(0,2) ;
	t1=variables(0,3) ;
	t2=variables(0,4) ;
	t3=variables(0,5) ;
	MatrixXd jsym(3,3);

	double mt4=1.0/(r1*r1+r2*r2+r3*r3+4.0);
	double mt1=mt4/a/b;
	double mt2=mt4/a/c;
	double mt3=mt4/b/c;
	

	jsym(0,0) = -(b*-4.0+a*r3*4.0-b*(r1*r1)+b*(r2*r2)+b*(r3*r3)+a*r1*r2*2.0)*mt1;
	jsym(0,1) = -(a*4.0+b*r3*4.0-a*(r1*r1)+a*(r2*r2)-a*(r3*r3)-b*r1*r2*2.0)*mt1;
	jsym(0,2) = (a*r1*4.0+b*r2*4.0-a*r2*r3*2.0+b*r1*r3*2.0)*mt1;
	jsym(1,0) = (b*r2*4.0+c*r3*4.0-b*r1*r3*2.0+c*r1*r2*2.0)*mt3;
	jsym(1,1) = -(c*-4.0+b*r1*4.0+c*(r1*r1)-c*(r2*r2)+c*(r3*r3)+b*r2*r3*2.0)*mt3;
	jsym(1,2) = -(b*4.0+c*r1*4.0-b*(r1*r1)-b*(r2*r2)+b*(r3*r3)-c*r2*r3*2.0)*mt3;
	jsym(2,0) = -(c*4.0+a*r2*4.0+c*(r1*r1)-c*(r2*r2)-c*(r3*r3)-a*r1*r3*2.0)*mt2;
	jsym(2,1) = (a*r1*4.0+c*r3*4.0+a*r2*r3*2.0-c*r1*r2*2.0)*mt2;
	jsym(2,2) = -(a*-4.0+c*r2*4.0+a*(r1*r1)+a*(r2*r2)-a*(r3*r3)+c*r1*r3*2.0)*mt2;
	return jsym;
}
MatrixXd jacobianForRotationAndTransitionUnitLength(const MatrixXd& parameters,const MatrixXd& variables)
{

	double a,b,c,d,e,f,t1,t2,r1,r2,r3;
	a=parameters(0,0) ;
	b=parameters(0,1) ;
	c=parameters(0,2) ;
	d=parameters(0,3) ;
	e=parameters(0,4) ;
	f=parameters(0,5) ;

	r1=variables(0,0) ;
	r2=variables(0,1) ;
	r3=variables(0,2) ;
	t1=variables(0,3) ;
	t2=variables(0,4) ;
	MatrixXd A0(3,5);
	double t4 ,t5 ,t6 ,t7 ,t8 ,t9 ,t10 ,t11 ,t12 ,t13 ,t14 ,t28 ,t15 ,t16 ,t17 ,t18 ,t26 ,t19 ,t20 ,t21 ,t22 ,t23 ,t24 ,t25 ,t27 ,t29 ,t30 ,t31 ,t32 ,t33 ,t34 ,t35 ,t36 ,t37 ,t38 ,t39 ,t40 ,t41 ,t45 ,t42 ,t43 ,t44 ,t46 ,t47 ,t48 ,t49 ,t50 ,t51 ,t52 ,t53 ,t54 ,t55 ,t56 ,t57 ,t60 ,t58 ,t59 ,t61 ,t62 ,t63 ,t64 ,t65 ,t66 ,t67 ,t68 ,t69 ,t70 ,t71 ,t72 ,t73 ,t74 ,t75 ,t76 ,t77 ,t78 ,t79 ,t80 ,t81 ,t82 ,t83 ,t84 ,t85 ,t86 ,t87 ,t88 ,t89 ,t90 ,t91 ;
	t4 = r1*r1;
	t5 = t4*(1.0/4.0);
	t6 = r2*r2;
	t7 = t6*(1.0/4.0);
	t8 = r3*r3;
	t9 = t8*(1.0/4.0);
	t10 = t5+t7+t9+1.0;
	t11 = 1.0/(t10*t10);
	t12 = 1.0/t10;
	t13 = cos(t1);
	t14 = cos(t2);
	t28 = t13*t14;
	t15 = b-t28;
	t16 = r1*t12*(1.0/2.0);
	t17 = t5+t7+t9-1.0;
	t18 = sin(t2);
	t26 = t13*t18;
	t19 = b-t26;
	t20 = r1*r3*t11*(1.0/2.0);
	t21 = r2*t4*t11*(1.0/4.0);
	t22 = sin(t1);
	t23 = b-t22;
	t24 = 1.0/a;
	t25 = r1*r2*r3*t11*(1.0/4.0);
	t27 = r1*t6*t11*(1.0/4.0);
	t29 = r2*t12*(1.0/2.0);
	t30 = 1.0/b;
	t31 = r2*r3*t11*(1.0/2.0);
	t32 = r1*r2*t11*(1.0/2.0);
	t33 = d*d;
	t34 = e*e;
	t35 = f*f;
	t36 = t33+t34+t35;
	t37 = 1.0/sqrt(t36);
	t38 = r3*t12*(1.0/2.0);
	t39 = r3*t4*t11*(1.0/4.0);
	t40 = t4+t6+t8+4.0;
	t41 = 1.0/t40;
	t45 = t4*t11*(1.0/2.0);
	t42 = t12+t25-t45;
	t43 = t23*t42;
	t44 = r1*t11*t17*(1.0/2.0);
	t46 = r1*t8*t11*(1.0/4.0);
	t47 = 1.0/c;
	t48 = r3*t6*t11*(1.0/4.0);
	t49 = t6*t11*(1.0/2.0);
	t50 = r2*t11*t17*(1.0/2.0);
	t51 = -t16+t27+t31;
	t52 = t29+t50-r2*t6*t11*(1.0/4.0);
	t53 = t19*t52;
	t54 = r2*t8*t11*(1.0/4.0);
	t55 = t8*t11*(1.0/2.0);
	t56 = -t12+t25+t55;
	t57 = t15*t56;
	t60 = r3*t11*t17*(1.0/2.0);
	t58 = t38+t48-t60;
	t59 = t19*t58;
	t61 = t16+t44-r1*t4*t11*(1.0/4.0);
	t62 = t15*t61;
	t63 = t20-t21+t29;
	t64 = t19*t63;
	t65 = t32+t38-t39;
	t66 = -t12+t25+t45;
	t67 = t19*t66;
	t68 = t16-t44+t46;
	t69 = t23*t68;
	t70 = t67+t69-t15*t65;
	t71 = t47*t70;
	t72 = t32-t38+t48;
	t73 = t19*t72;
	t74 = t12+t25-t49;
	t75 = t15*t74;
	t76 = t29-t50+t54;
	t77 = t23*t76;
	t78 = t73+t75+t77;
	t79 = t47*t78;
	t80 = -t12+t25+t49;
	t81 = t23*t80;
	t82 = t16-t27+t31;
	t83 = -t16+t31+t46;
	t84 = t23*t83;
	t85 = t16+t31-t46;
	t86 = t15*t85;
	t87 = t20-t29+t54;
	t88 = t38+t60-r3*t8*t11*(1.0/4.0);
	t89 = t23*t88;
	t90 = t86+t89-t19*t87;
	t91 = t47*t90;
	A0(0,0) = t37*(t24*(t62+t64-t23*(t32+t39-r3*t12*(1.0/2.0)))+t30*(t43+t15*(t20+t21-r2*t12*(1.0/2.0))+t19*(t16+t27-r1*t11*t17*(1.0/2.0))));
	A0(0,1) = -t37*(t30*(t53-t15*t51+t23*(t32+t38-r3*t6*t11*(1.0/4.0)))+t24*(t81-t19*t82+t15*(t21+t29-r2*t11*t17*(1.0/2.0))));
	A0(0,2) = -t37*(t24*(t84+t19*(t12+t25-t8*t11*(1.0/2.0))+t15*(t38+t39-r3*t11*t17*(1.0/2.0)))-t30*(t57+t59-t23*(t20+t29-r2*t8*t11*(1.0/4.0))));
	A0(0,3) = -t24*t30*t37*t41*(a*r1*t13*4.0+b*r2*t13*4.0+a*t18*t22*4.0-b*t14*t22*4.0-a*r2*r3*t13*2.0+b*r1*r3*t13*2.0+a*r3*t14*t22*4.0+b*r3*t18*t22*4.0-a*t4*t18*t22+a*t6*t18*t22-a*t8*t18*t22-b*t4*t14*t22+b*t6*t14*t22+b*t8*t14*t22+a*r1*r2*t14*t22*2.0-b*r1*r2*t18*t22*2.0);
	A0(0,4) = -t13*t24*t30*t37*t41*(a*t14*-4.0-b*t18*4.0+a*r3*t18*4.0-b*r3*t14*4.0+a*t4*t14-a*t6*t14+a*t8*t14-b*t4*t18+b*t6*t18+b*t8*t18+a*r1*r2*t18*2.0+b*r1*r2*t14*2.0);
	A0(1,0) = t37*(t71-t30*(t43+t15*(t20+t21-t29)+t19*(t16+t27-t44)));
	A0(1,1) = t37*(t79+t30*(t53-t15*t51+t23*(t32+t38-t48)));
	A0(1,2) = -t37*(t91+t30*(t57+t59-t23*(t20+t29-t54)));
	A0(1,3) = t30*t37*t41*t47*(b*t13*4.0+c*r1*t13*4.0-b*t4*t13-b*t6*t13+b*t8*t13+c*t18*t22*4.0-c*r2*r3*t13*2.0+b*r2*t14*t22*4.0-b*r1*t18*t22*4.0+c*r3*t14*t22*4.0-c*t4*t18*t22+c*t6*t18*t22-c*t8*t18*t22-b*r1*r3*t14*t22*2.0-b*r2*r3*t18*t22*2.0+c*r1*r2*t14*t22*2.0);
	A0(1,4) = t13*t30*t37*t41*t47*(c*t14*-4.0+b*r1*t14*4.0+b*r2*t18*4.0+c*r3*t18*4.0+c*t4*t14-c*t6*t14+c*t8*t14+b*r2*r3*t14*2.0-b*r1*r3*t18*2.0+c*r1*r2*t18*2.0);
	A0(2,0) = -t37*(t71+t24*(t62+t64-t23*(t32-t38+t39)));
	A0(2,1) = -t37*(t79-t24*(t81-t19*t82+t15*(t21+t29-t50)));
	A0(2,2) = t37*(t91+t24*(t84+t19*(t12+t25-t55)+t15*(t38+t39-t60)));
	A0(2,3) = t24*t37*t41*t47*(a*t13*-4.0+c*r2*t13*4.0+a*t4*t13+a*t6*t13-a*t8*t13-c*t14*t22*4.0+c*r1*r3*t13*2.0-a*r2*t14*t22*4.0+a*r1*t18*t22*4.0+c*r3*t18*t22*4.0-c*t4*t14*t22+c*t6*t14*t22+c*t8*t14*t22+a*r1*r3*t14*t22*2.0+a*r2*r3*t18*t22*2.0-c*r1*r2*t18*t22*2.0);
	A0(2,4) = -t13*t24*t37*t41*t47*(c*t18*4.0+a*r1*t14*4.0+a*r2*t18*4.0+c*r3*t14*4.0+c*t4*t18-c*t6*t18-c*t8*t18+a*r2*r3*t14*2.0-a*r1*r3*t18*2.0-c*r1*r2*t14*2.0);


	return A0;
}

MatrixXd jacobianForRotationAndTransitionUnitLength_(const MatrixXd& parameters,const MatrixXd& variables)
{

	double a,b,c,d,e,f,t1,t2,r1,r2,r3;
	a=parameters(0,0) ;
	b=parameters(0,1) ;
	c=parameters(0,2) ;
	d=parameters(0,3) ;
	e=parameters(0,4) ;
	f=parameters(0,5) ;

	r1=variables(0,0) ;
	r2=variables(0,1) ;
	r3=variables(0,2) ;
	t1=variables(0,3) ;
	t2=variables(0,4) ;
	//t3=variables(0,5) ;
	MatrixXd A0(3,5);
	double t4 ,t5 ,t6 ,t7 ,t8 ,t9 ,t10 ,t11 ,t12 ,t13 ,t14 ,t15 ,t16 ,t17 ,t18 ,t19 ,t20 ,t21 ,t22 ,t23 ,t24 ,t25 ,t26 ,t27 ,t28 ,t29 ,t30 ,t31 ,t32 ,t33 ,t34 ,t35 ,t36 ,t37 ,t38 ,t39 ,t40 ,t41 ,t42 ,t43 ,t44 ,t45 ,t46 ,t47 ,t48 ,t49 ,t50 ,t51 ,t52 ,t53 ,t54 ,t55 ,t56 ,t57 ,t58 ,t59 ,t60 ,t61 ,t62 ,t63 ,t64 ,t65 ,t66 ,t67 ,t68 ,t69 ,t70 ,t71 ,t72 ,t73 ,t74 ,t75 ,t76 ,t77 ,t78 ,t79 ,t80 ,t81 ,t82 ,t114 ,t121 ,t83 ,t84 ,t85 ,t86 ,t87 ,t88 ,t89 ,t90 ,t91 ,t92 ,t93 ,t94 ,t95 ,t96 ,t97 ,t98 ,t99 ,t100 ,t101 ,t102 ,t103 ,t104 ,t105 ,t106 ,t107 ,t108 ,t109 ,t110 ,t111 ,t112 ,t113 ,t115 ,t116 ,t117 ,t118 ,t119 ,t120 ;
	t4 = r1*r1;
	t5 = sin(t1);
	t6 = r2*r2;
	t7 = r3*r3;
	t8 = cos(t1);
	t9 = cos(t2);
	t10 = sin(t2);
	t11 = t4+t6+t7+4.0;
	t12 = 1.0/(t11*t11);
	t13 = 1.0/a;
	t14 = f*1.6E1;
	t15 = t5*1.6E1;
	t16 = t4*t5*4.0;
	t17 = f*t6*4.0;
	t18 = f*t7*4.0;
	t19 = f*r1*r2*r3*4.0;
	t20 = 1.0/b;
	t21 = f*r3*8.0;
	t22 = f*r3*t7*2.0;
	t23 = f*r3*t6*2.0;
	t24 = r1*r2*t5*8.0;
	t25 = r3*t4*t5*2.0;
	t26 = 1.0/t11;
	t27 = t5*t6*4.0;
	t28 = e*r1*1.6E1;
	t29 = f*t4*4.0;
	t30 = d*r2*t4*2.0;
	t31 = e*r1*t6*4.0;
	t32 = r2*t8*t9*8.0;
	t33 = d*r1*r3*8.0;
	t34 = r2*t6*t8*t9*2.0;
	t35 = r2*t7*t8*t9*2.0;
	t36 = r1*r2*r3*t5*4.0;
	t37 = e*1.6E1;
	t38 = e*t4*4.0;
	t39 = e*t6*4.0;
	t40 = t7*t8*t10*4.0;
	t41 = e*r1*r2*r3*4.0;
	t42 = d*r1*8.0;
	t43 = d*r1*t4*2.0;
	t44 = d*r1*t7*2.0;
	t45 = e*r2*t4*4.0;
	t46 = e*r2*t7*4.0;
	t47 = f*r3*t4*2.0;
	t48 = r3*t5*t6*2.0;
	t49 = f*r1*r2*8.0;
	t50 = r1*t6*t8*t9*2.0;
	t51 = r2*r3*t8*t9*8.0;
	t52 = 1.0/c;
	t53 = d*1.6E1;
	t54 = d*t4*4.0;
	t55 = d*t6*4.0;
	t56 = t7*t8*t9*4.0;
	t57 = r1*r2*r3*t8*t9*4.0;
	t58 = e*r2*8.0;
	t59 = e*r2*t6*2.0;
	t60 = e*r2*t7*2.0;
	t61 = e*r1*r3*8.0;
	t62 = r2*t4*t8*t10*2.0;
	t63 = f*r2*8.0;
	t64 = d*t7*4.0;
	t65 = f*r2*t6*2.0;
	t66 = f*r2*t4*2.0;
	t67 = r2*t5*t7*2.0;
	t68 = r3*t8*t10*1.6E1;
	t69 = f*r1*r3*8.0;
	t70 = t6*t8*t9*4.0;
	t71 = d*r1*r2*r3*4.0;
	t72 = r3*t6*t8*t10*4.0;
	t73 = d*r3*8.0;
	t74 = d*r3*t7*2.0;
	t75 = e*t7*4.0;
	t76 = r1*t5*1.6E1;
	t77 = d*r3*t6*2.0;
	t78 = r1*t5*t7*4.0;
	t79 = d*r1*r2*8.0;
	t80 = t4*t8*t10*4.0;
	t81 = r3*t4*t8*t9*2.0;
	t82 = r1*r2*r3*t8*t10*4.0;
	t114 = t8*t10*1.6E1;
	t121 = t6*t8*t10*4.0;
	t83 = t37-t38+t39-t40-t41+t73+t74+t75+t76+t77+t78+t79+t80+t81+t82-t114-t121-f*r1*1.6E1-d*r3*t4*2.0-f*r1*t7*4.0-r3*t8*t9*8.0-r1*r2*t8*t9*8.0-r3*t6*t8*t9*2.0-r3*t7*t8*t9*2.0;
	t84 = d*r1*t6*4.0;
	t85 = d*r1*t7*4.0;
	t86 = e*r2*t4*2.0;
	t87 = r1*r3*t8*t10*8.0;
	t88 = r2*t7*t8*t10*2.0;
	t89 = e*r1*8.0;
	t90 = e*r1*t4*2.0;
	t91 = e*r1*t7*2.0;
	t92 = r2*t8*t9*1.6E1;
	t93 = e*r2*r3*8.0;
	t94 = r2*t4*t8*t9*4.0;
	t95 = r1*t6*t8*t10*2.0;
	t96 = t14-t15-t16-t17+t18-t19+t27+t29+t36+t89+t90+t91+t92+t93+t94+t95-d*r2*1.6E1-t5*t7*4.0-d*r2*t4*4.0-e*r1*t6*2.0-r1*t8*t10*8.0-r2*r3*t8*t10*8.0-r1*t4*t8*t10*2.0-r1*t7*t8*t10*2.0;
	t97 = t12*t13*t96;
	t98 = f*r2*1.6E1;
	t99 = e*r3*t6*2.0;
	t100 = f*r2*t7*4.0;
	t101 = r3*t8*t10*8.0;
	t102 = e*r1*r2*8.0;
	t103 = r3*t7*t8*t10*2.0;
	t104 = r3*t4*t8*t10*2.0;
	t105 = t53+t54-t55-t56-t57+t64+t70+t71+t98+t99+t100+t101+t102+t103+t104-e*r3*8.0-r2*t5*1.6E1-t8*t9*1.6E1-e*r3*t4*2.0-e*r3*t7*2.0-r2*t5*t7*4.0-t4*t8*t9*4.0-r1*r2*t8*t10*8.0-r3*t6*t8*t10*2.0;
	t106 = t12*t52*t105;
	t107 = d*r1*t6*2.0;
	t108 = f*r3*t4*4.0;
	t109 = f*r3*t6*4.0;
	t110 = d*r2*r3*8.0;
	t111 = r1*t7*t8*t9*2.0;
	t112 = t42+t43-t44-t50-t51+t58+t59-t60-t61-t62+t86+t87+t88+t107+t108+t109+t110+t111-r3*t4*t5*4.0-r3*t5*t6*4.0-r1*t8*t9*8.0-r2*t8*t10*8.0-r1*t4*t8*t9*2.0-r2*t6*t8*t10*2.0;
	t113 = r1*t4*t5*2.0;
	t115 = d*r3*1.6E1;
	t116 = r1*t5*8.0;
	t117 = d*r3*t4*4.0;
	t118 = f*r1*t7*2.0;
	t119 = r1*t5*t6*2.0;
	t120 = f*r2*r3*8.0;
	A0(0,0) = t12*t13*(t21+t22+t23+t24+t25+t58+t59+t60+t61+t62+t84+t85-r3*t5*8.0-f*r1*r2*8.0-e*r2*t4*2.0-f*r3*t4*2.0-r3*t5*t6*2.0-r3*t5*t7*2.0-r2*t8*t10*8.0-r1*r3*t8*t10*8.0-r1*t6*t8*t9*4.0-r1*t7*t8*t9*4.0-r2*t6*t8*t10*2.0-r2*t7*t8*t10*2.0)+t12*t20*(t14-t15+t16+t17+t18+t19+t28+t30+t31+t32+t33+t34+t35-d*r2*8.0-f*t4*4.0-t5*t6*4.0-t5*t7*4.0-d*r2*t6*2.0-d*r2*t7*2.0-r1*t8*t10*1.6E1-r1*r2*r3*t5*4.0-r1*r3*t8*t9*8.0-r2*t4*t8*t9*2.0-r1*t6*t8*t10*4.0);
	A0(0,1) = t97-t12*t20*(t21+t22-t23-t24-t25+t42+t43+t44+t45+t46+t47+t48+t49+t50+t51-r3*t5*8.0-d*r2*r3*8.0-d*r1*t6*2.0-r3*t5*t7*2.0-r1*t8*t9*8.0-r1*t4*t8*t9*2.0-r2*t4*t8*t10*4.0-r1*t7*t8*t9*2.0-r2*t7*t8*t10*4.0);
	A0(0,2) = -t12*t20*(t53+t54+t55+t56+t57+t63+t65+t66+t67+t68+t69+t72-e*r3*1.6E1-d*t7*4.0-r2*t5*8.0-t8*t9*1.6E1-e*r3*t6*4.0-f*r2*t7*2.0-r1*r3*t5*8.0-r2*t4*t5*2.0-r2*t5*t6*2.0-t4*t8*t9*4.0-t6*t8*t9*4.0-d*r1*r2*r3*4.0)-t12*t13*(t37+t38+t39+t40+t41+t113+t115+t116+t117+t118+t119+t120-f*r1*8.0-e*t7*4.0-t8*t10*1.6E1-f*r1*t4*2.0-f*r1*t6*2.0-r2*r3*t5*8.0-r1*t5*t7*2.0-r3*t8*t9*1.6E1-t4*t8*t10*4.0-t6*t8*t10*4.0-r3*t4*t8*t9*4.0-r1*r2*r3*t8*t10*4.0);
	A0(0,3) = -t13*t20*t26*(a*r1*t8*4.0+b*r2*t8*4.0+a*t5*t10*4.0-b*t5*t9*4.0-a*r2*r3*t8*2.0+b*r1*r3*t8*2.0+a*r3*t5*t9*4.0+b*r3*t5*t10*4.0-a*t4*t5*t10+a*t5*t6*t10-a*t5*t7*t10-b*t4*t5*t9+b*t5*t6*t9+b*t5*t7*t9+a*r1*r2*t5*t9*2.0-b*r1*r2*t5*t10*2.0);
	A0(0,4) = -t8*t13*t20*t26*(a*t9*-4.0-b*t10*4.0+a*r3*t10*4.0-b*r3*t9*4.0+a*t4*t9-a*t6*t9+a*t7*t9-b*t4*t10+b*t6*t10+b*t7*t10+a*r1*r2*t10*2.0+b*r1*r2*t9*2.0);
	A0(1,0) = -t12*t20*(t14-t15+t16+t17+t18+t19-t27+t28-t29+t30+t31+t32+t33+t34+t35-t36-d*r2*8.0-t5*t7*4.0-d*r2*t6*2.0-d*r2*t7*2.0-r1*t8*t10*1.6E1-r1*r3*t8*t9*8.0-r2*t4*t8*t9*2.0-r1*t6*t8*t10*4.0)-t12*t52*t83;
	A0(1,1) = t106+t12*t20*(t21+t22-t23-t24-t25+t42+t43+t44+t45+t46+t47+t48+t49+t50+t51-r3*t5*8.0-d*r2*r3*8.0-d*r1*t6*2.0-r3*t5*t7*2.0-r1*t8*t9*8.0-r1*t4*t8*t9*2.0-r2*t4*t8*t10*4.0-r1*t7*t8*t9*2.0-r2*t7*t8*t10*4.0);
	A0(1,2) = t12*t20*(t53+t54+t55+t56+t57+t63-t64+t65+t66+t67+t68+t69-t70-t71+t72-e*r3*1.6E1-r2*t5*8.0-t8*t9*1.6E1-e*r3*t6*4.0-f*r2*t7*2.0-r1*r3*t5*8.0-r2*t4*t5*2.0-r2*t5*t6*2.0-t4*t8*t9*4.0)-t12*t52*t112;
	A0(1,3) = t20*t26*t52*(b*t8*4.0+c*r1*t8*4.0-b*t4*t8-b*t6*t8+b*t7*t8+c*t5*t10*4.0-c*r2*r3*t8*2.0-b*r1*t5*t10*4.0+b*r2*t5*t9*4.0+c*r3*t5*t9*4.0-c*t4*t5*t10+c*t5*t6*t10-c*t5*t7*t10-b*r1*r3*t5*t9*2.0-b*r2*r3*t5*t10*2.0+c*r1*r2*t5*t9*2.0);
	A0(1,4) = t8*t20*t26*t52*(c*t9*-4.0+b*r1*t9*4.0+b*r2*t10*4.0+c*r3*t10*4.0+c*t4*t9-c*t6*t9+c*t7*t9-b*r1*r3*t10*2.0+b*r2*r3*t9*2.0+c*r1*r2*t10*2.0);
	A0(2,0) = -t12*t13*(t21+t22+t23+t24+t25-t47-t48-t49+t58+t59+t60+t61+t62+t84+t85-t86-t87-t88-r3*t5*8.0-r3*t5*t7*2.0-r2*t8*t10*8.0-r1*t6*t8*t9*4.0-r1*t7*t8*t9*4.0-r2*t6*t8*t10*2.0)+t12*t52*t83;
	A0(2,1) = -t97-t106;
	A0(2,2) = t12*t13*(t37+t38+t39+t40+t41-t75-t80-t82+t113-t114+t115+t116+t117+t118+t119+t120-t121-f*r1*8.0-f*r1*t4*2.0-f*r1*t6*2.0-r2*r3*t5*8.0-r1*t5*t7*2.0-r3*t8*t9*1.6E1-r3*t4*t8*t9*4.0)+t12*t52*t112;
	A0(2,3) = t13*t26*t52*(a*t8*-4.0+c*r2*t8*4.0+a*t4*t8+a*t6*t8-a*t7*t8-c*t5*t9*4.0+c*r1*r3*t8*2.0+a*r1*t5*t10*4.0-a*r2*t5*t9*4.0+c*r3*t5*t10*4.0-c*t4*t5*t9+c*t5*t6*t9+c*t5*t7*t9+a*r1*r3*t5*t9*2.0+a*r2*r3*t5*t10*2.0-c*r1*r2*t5*t10*2.0);
	A0(2,4) = -t8*t13*t26*t52*(c*t10*4.0+a*r1*t9*4.0+a*r2*t10*4.0+c*r3*t9*4.0+c*t4*t10-c*t6*t10-c*t7*t10-a*r1*r3*t10*2.0+a*r2*r3*t9*2.0-c*r1*r2*t9*2.0);
	return A0;
}


MatrixXd jacobianForRotationAndTransition(const MatrixXd& parameters,const MatrixXd& variables)
{

	double a,b,c,d,e,f,t1,t2,t3,r1,r2,r3;
	a=parameters(0,0) ;
	b=parameters(0,1) ;
	c=parameters(0,2) ;
	d=parameters(0,3) ;
	e=parameters(0,4) ;
	f=parameters(0,5) ;

	r1=variables(0,0) ;
	r2=variables(0,1) ;
	r3=variables(0,2) ;
	t1=variables(0,3) ;
	t2=variables(0,4) ;
	t3=variables(0,5) ;
	

	
	MatrixXd A0(3,6);
	double t5 ,t6 ,t7 ,t8 ,t9 ,t10 ,t11 ,t12 ,t13 ,t14 ,t15 ,t16 ,t17 ,t18 ,t19 ,t20 ,t21 ,t22 ,t23 ,t24 ,t25 ,t26 ,t27 ,t28 ,t29 ,t30 ,t31 ,t32 ,t33 ,t35 ,t36 ,t37 ,t38 ,t39 ,t40 ,t41 ,t42 ,t43 ,t44 ,t45 ,t46 ,t47 ,t48 ,t49 ,t50 ,t34 ,t51 ,t52 ,t53 ,t54 ,t55 ,t56 ,t57 ,t58 ,t59 ,t60 ,t61 ,t62 ,t63 ,t64 ,t65 ,t66 ,t67 ,t68 ,t69 ,t70 ,t75 ,t76 ,t77 ,t78 ,t79 ,t80 ,t81 ,t82 ,t83 ,t84 ,t85 ,t86 ,t87 ,t88 ,t89 ,t90 ,t71 ,t72 ,t73 ,t74 ,t91 ,t92 ,t93 ,t94 ,t95 ,t96 ,t97 ,t98 ,t99 ,t100 ,t101 ,t102 ,t103 ,t104 ,t105 ,t106 ,t107 ,t108 ,t109 ,t110 ,t111 ,t112 ,t113 ,t114 ,t115 ,t116 ,t117 ,t118 ,t119 ,t124 ,t125 ,t126 ,t127 ,t128 ,t129 ,t130 ,t131 ,t132 ,t133 ,t134 ,t135 ,t136 ,t137 ,t138 ,t139 ,t120 ,t121 ,t122 ,t123 ,t140 ,t141 ,t142 ,t143 ,t144 ,t145 ,t146 ;
	t5 = 1.0/a;
	t6 = 1.0/b;
	t7 = d*d;
	t8 = e*e;
	t9 = f*f;
	t10 = t7+t8+t9;
	t11 = 1.0/sqrt(t10);
	t12 = r1*r1;
	t13 = r2*r2;
	t14 = r3*r3;
	t15 = t12+t13+t14+4.0;
	t16 = 1.0/t15;
	t17 = 1.0/(t15*t15);
	t18 = a*e*4.0;
	t19 = b*t1*4.0;
	t20 = a*t2*t12;
	t21 = a*t2*t14;
	t22 = b*t1*t12;
	t23 = a*d*r3*4.0;
	t24 = b*e*r3*4.0;
	t25 = a*r1*t3*4.0;
	t26 = b*r2*t3*4.0;
	t27 = a*e*t13;
	t28 = b*d*t13;
	t29 = b*d*t14;
	t30 = a*d*r1*r2*2.0;
	t31 = a*f*r2*r3*2.0;
	t32 = b*r1*r2*t2*2.0;
	t33 = b*r1*r3*t3*2.0;
	t35 = b*d*4.0;
	t36 = a*t2*4.0;
	t37 = a*t2*t13;
	t38 = b*t1*t13;
	t39 = b*t1*t14;
	t40 = a*f*r1*4.0;
	t41 = b*f*r2*4.0;
	t42 = a*r3*t1*4.0;
	t43 = b*r3*t2*4.0;
	t44 = a*e*t12;
	t45 = b*d*t12;
	t46 = a*e*t14;
	t47 = b*e*r1*r2*2.0;
	t48 = b*f*r1*r3*2.0;
	t49 = a*r1*r2*t1*2.0;
	t50 = a*r2*r3*t3*2.0;
	t34 = t18+t19+t20+t21+t22+t23+t24+t25+t26+t27+t28+t29+t30+t31+t32+t33-t35-t36-t37-t38-t39-t40-t41-t42-t43-t44-t45-t46-t47-t48-t49-t50;
	t51 = b*e*4.0;
	t52 = b*d*r3*2.0;
	t53 = b*r1*t3*2.0;
	t54 = 1.0/c;
	t55 = b*f*4.0;
	t56 = c*t2*4.0;
	t57 = b*t3*t12;
	t58 = b*t3*t13;
	t59 = c*t2*t13;
	t60 = b*e*r1*4.0;
	t61 = c*f*r1*4.0;
	t62 = b*r2*t1*4.0;
	t63 = c*r3*t1*4.0;
	t64 = c*e*t12;
	t65 = b*f*t14;
	t66 = c*e*t14;
	t67 = b*d*r1*r3*2.0;
	t68 = b*e*r2*r3*2.0;
	t69 = c*r1*r2*t1*2.0;
	t70 = c*r2*r3*t3*2.0;
	t75 = c*e*4.0;
	t76 = b*t3*4.0;
	t77 = b*t3*t14;
	t78 = c*t2*t12;
	t79 = c*t2*t14;
	t80 = b*d*r2*4.0;
	t81 = c*d*r3*4.0;
	t82 = b*r1*t2*4.0;
	t83 = c*r1*t3*4.0;
	t84 = b*f*t12;
	t85 = b*f*t13;
	t86 = c*e*t13;
	t87 = c*d*r1*r2*2.0;
	t88 = c*f*r2*r3*2.0;
	t89 = b*r1*r3*t1*2.0;
	t90 = b*r2*r3*t2*2.0;
	t71 = t55+t56+t57+t58+t59+t60+t61+t62+t63+t64+t65+t66+t67+t68+t69+t70-t75-t76-t77-t78-t79-t80-t81-t82-t83-t84-t85-t86-t87-t88-t89-t90;
	t72 = b*d*r1*2.0;
	t73 = b*e*r2*2.0;
	t74 = b*f*r3*2.0;
	t91 = b*r2*4.0;
	t92 = b*r1*r3*2.0;
	t93 = b*t13;
	t94 = b*t14;
	t95 = c*r1*t1*2.0;
	t96 = c*r2*t2*2.0;
	t97 = c*r3*t3*2.0;
	t98 = a*f*4.0;
	t99 = c*t1*4.0;
	t100 = a*d*4.0;
	t101 = c*f*4.0;
	t102 = a*f*r2*2.0;
	t103 = c*e*r1*2.0;
	t104 = a*r3*t2*2.0;
	t105 = c*r2*t1*2.0;
	t106 = a*t3*t12;
	t107 = a*t3*t13;
	t108 = c*t1*t12;
	t109 = a*e*r1*4.0;
	t110 = c*e*r3*4.0;
	t111 = a*r2*t1*4.0;
	t112 = c*r2*t3*4.0;
	t113 = c*d*t13;
	t114 = a*f*t14;
	t115 = c*d*t14;
	t116 = a*d*r1*r3*2.0;
	t117 = a*e*r2*r3*2.0;
	t118 = c*r1*r2*t2*2.0;
	t119 = c*r1*r3*t3*2.0;
	t124 = c*d*4.0;
	t125 = a*t3*4.0;
	t126 = a*t3*t14;
	t127 = c*t1*t13;
	t128 = c*t1*t14;
	t129 = a*d*r2*4.0;
	t130 = c*f*r2*4.0;
	t131 = a*r1*t2*4.0;
	t132 = c*r3*t2*4.0;
	t133 = a*f*t12;
	t134 = c*d*t12;
	t135 = a*f*t13;
	t136 = c*e*r1*r2*2.0;
	t137 = c*f*r1*r3*2.0;
	t138 = a*r1*r3*t1*2.0;
	t139 = a*r2*r3*t2*2.0;
	t120 = t98+t99+t106+t107+t108+t109+t110+t111+t112+t113+t114+t115+t116+t117+t118+t119-t124-t125-t126-t127-t128-t129-t130-t131-t132-t133-t134-t135-t136-t137-t138-t139;
	t121 = a*r1*t1*2.0;
	t122 = a*r2*t2*2.0;
	t123 = a*r3*t3*2.0;
	t140 = c*t12;
	t141 = c*t14;
	t142 = a*r1*4.0;
	t143 = c*r3*4.0;
	t144 = c*r1*r2*2.0;
	t145 = a*4.0;
	t146 = a*t13;
	A0(0,0) = t5*t6*t11*t16*(t72+t73+t74+t98-a*t3*4.0-a*d*r2*2.0+a*e*r1*2.0-a*r1*t2*2.0+a*r2*t1*2.0-b*r1*t1*2.0-b*r2*t2*2.0-b*r3*t3*2.0)+r1*t5*t6*t11*t17*t34*2.0;
	A0(0,1) = t5*t6*t11*t16*(t55+t121+t122+t123-b*t3*4.0-a*d*r1*2.0-a*e*r2*2.0-b*d*r2*2.0+b*e*r1*2.0-a*f*r3*2.0-b*r1*t2*2.0+b*r2*t1*2.0)+r2*t5*t6*t11*t17*t34*2.0;
	A0(0,2) = -t5*t6*t11*t16*(t51+t52+t53+t100+t102+t104-a*t1*4.0-b*t2*4.0-a*e*r3*2.0-b*f*r1*2.0-a*r2*t3*2.0-b*r3*t1*2.0)+r3*t5*t6*t11*t17*t34*2.0;
	A0(0,3) = t5*t6*t11*t16*(b*-4.0+t93+t94+a*r3*4.0-b*t12+a*r1*r2*2.0);
	A0(0,4) = t5*t6*t11*t16*(t145+t146+b*r3*4.0-a*t12-a*t14-b*r1*r2*2.0);
	A0(0,5) = -t5*t6*t11*t16*(t91+t92+t142-a*r2*r3*2.0);
	A0(1,0) = -t6*t11*t16*t54*(t51+t52+t53+t101+t103+t105-b*t2*4.0-c*t3*4.0-c*d*r2*2.0-b*f*r1*2.0-b*r3*t1*2.0-c*r1*t2*2.0)+r1*t6*t11*t17*t54*t71*2.0;
	A0(1,1) = -t6*t11*t16*t54*(t19-t35+t95+t96+t97-c*d*r1*2.0+b*e*r3*2.0-b*f*r2*2.0-c*e*r2*2.0-c*f*r3*2.0+b*r2*t3*2.0-b*r3*t2*2.0)+r2*t6*t11*t17*t54*t71*2.0;
	A0(1,2) = -t6*t11*t16*t54*(t72+t73+t74+t99-c*d*4.0+c*e*r3*2.0-c*f*r2*2.0-b*r1*t1*2.0-b*r2*t2*2.0-b*r3*t3*2.0+c*r2*t3*2.0-c*r3*t2*2.0)+r3*t6*t11*t17*t54*t71*2.0;
	A0(1,3) = -t6*t11*t16*t54*(t91-t92+t143+t144);
	A0(1,4) = t6*t11*t16*t54*(c*-4.0+t140+t141+b*r1*4.0-c*t13+b*r2*r3*2.0);
	A0(1,5) = t6*t11*t16*t54*(b*4.0-t93+t94+c*r1*4.0-b*t12-c*r2*r3*2.0);
	A0(2,0) = t5*t11*t16*t54*(t18-t36+t95+t96+t97+a*d*r3*2.0-a*f*r1*2.0-c*d*r1*2.0-c*e*r2*2.0-c*f*r3*2.0+a*r1*t3*2.0-a*r3*t1*2.0)-r1*t5*t11*t17*t54*t120*2.0;
	A0(2,1) = -t5*t11*t16*t54*(t100+t101+t102+t103+t104+t105-a*t1*4.0-c*t3*4.0-a*e*r3*2.0-c*d*r2*2.0-a*r2*t3*2.0-c*r1*t2*2.0)-r2*t5*t11*t17*t54*t120*2.0;
	A0(2,2) = -t5*t11*t16*t54*(t56-t75+t121+t122+t123-a*d*r1*2.0-a*e*r2*2.0-a*f*r3*2.0-c*d*r3*2.0+c*f*r1*2.0-c*r1*t3*2.0+c*r3*t1*2.0)-r3*t5*t11*t17*t54*t120*2.0;
	A0(2,3) = t5*t11*t16*t54*(c*4.0+t140-t141+a*r2*4.0-c*t13-a*r1*r3*2.0);
	A0(2,4) = -t5*t11*t16*t54*(t142+t143-t144+a*r2*r3*2.0);
	A0(2,5) = t5*t11*t16*t54*(-t145+t146+c*r2*4.0+a*t12-a*t14+c*r1*r3*2.0);


	return A0;
}
MatrixXd jacobianForRotationAndTransition_(const MatrixXd& parameters,const MatrixXd& variables)
{

	double a,b,c,d,e,f,t1,t2,t3,r1,r2,r3;
	a=parameters(0,0) ;
	b=parameters(0,1) ;
	c=parameters(0,2) ;
	d=parameters(0,3) ;
	e=parameters(0,4) ;
	f=parameters(0,5) ;

	r1=variables(0,0) ;
	r2=variables(0,1) ;
	r3=variables(0,2) ;
	t1=variables(0,3) ;
	t2=variables(0,4) ;
	t3=variables(0,5) ;
	double t5 ,t6 ,t7 ,t8 ,t9 ,t10 ,t11 ,t12 ,t13 ,t14 ,t15 ,t16 ,t17 ,t18 ,t19 ,t20 ,t21 ,t22 ,t23 ,t24 ,t25 ,t26 ,t27 ,t28 ,t30 ,t31 ,t32 ,t33 ,t34 ,t35 ,t36 ,t37 ,t38 ,t39 ,t40 ,t41 ,t42 ,t43 ,t44 ,t45 ,t29 ,t46 ,t47 ,t48 ,t49 ,t50 ,t51 ,t52 ,t53 ,t54 ,t55 ,t56 ,t57 ,t58 ,t59 ,t60 ,t61 ,t62 ,t63 ,t64 ,t65 ,t70 ,t71 ,t72 ,t73 ,t74 ,t75 ,t76 ,t77 ,t78 ,t79 ,t80 ,t81 ,t82 ,t83 ,t84 ,t85 ,t66 ,t67 ,t68 ,t69 ,t86 ,t87 ,t88 ,t89 ,t90 ,t91 ,t92 ,t93 ,t94 ,t95 ,t96 ,t97 ,t98 ,t99 ,t100 ,t101 ,t102 ,t103 ,t104 ,t105 ,t106 ,t107 ,t108 ,t109 ,t110 ,t111 ,t112 ,t113 ,t114 ,t119 ,t120 ,t121 ,t122 ,t123 ,t124 ,t125 ,t126 ,t127 ,t128 ,t129 ,t130 ,t131 ,t132 ,t133 ,t134 ,t115 ,t116 ,t117 ,t118 ,t135 ,t136 ,t137 ,t138 ,t139 ,t140 ,t141 ;
	t5 = 1.0/a;
	t6 = 1.0/b;
	t7 = r1*r1;
	t8 = r2*r2;
	t9 = r3*r3;
	t10 = t7+t8+t9+4.0;
	t11 = 1.0/t10;
	t12 = 1.0/(t10*t10);
	t13 = a*e*4.0;
	t14 = b*t1*4.0;
	t15 = a*t2*t7;
	t16 = a*t2*t9;
	t17 = b*t1*t7;
	t18 = a*d*r3*4.0;
	t19 = b*e*r3*4.0;
	t20 = a*r1*t3*4.0;
	t21 = b*r2*t3*4.0;
	t22 = a*e*t8;
	t23 = b*d*t8;
	t24 = b*d*t9;
	t25 = a*d*r1*r2*2.0;
	t26 = a*f*r2*r3*2.0;
	t27 = b*r1*r2*t2*2.0;
	t28 = b*r1*r3*t3*2.0;
	t30 = b*d*4.0;
	t31 = a*t2*4.0;
	t32 = a*t2*t8;
	t33 = b*t1*t8;
	t34 = b*t1*t9;
	t35 = a*f*r1*4.0;
	t36 = b*f*r2*4.0;
	t37 = a*r3*t1*4.0;
	t38 = b*r3*t2*4.0;
	t39 = a*e*t7;
	t40 = b*d*t7;
	t41 = a*e*t9;
	t42 = b*e*r1*r2*2.0;
	t43 = b*f*r1*r3*2.0;
	t44 = a*r1*r2*t1*2.0;
	t45 = a*r2*r3*t3*2.0;
	t29 = t13+t14+t15+t16+t17+t18+t19+t20+t21+t22+t23+t24+t25+t26+t27+t28-t30-t31-t32-t33-t34-t35-t36-t37-t38-t39-t40-t41-t42-t43-t44-t45;
	t46 = b*e*4.0;
	t47 = b*d*r3*2.0;
	t48 = b*r1*t3*2.0;
	t49 = 1.0/c;
	t50 = b*f*4.0;
	t51 = c*t2*4.0;
	t52 = b*t3*t7;
	t53 = b*t3*t8;
	t54 = c*t2*t8;
	t55 = b*e*r1*4.0;
	t56 = c*f*r1*4.0;
	t57 = b*r2*t1*4.0;
	t58 = c*r3*t1*4.0;
	t59 = c*e*t7;
	t60 = b*f*t9;
	t61 = c*e*t9;
	t62 = b*d*r1*r3*2.0;
	t63 = b*e*r2*r3*2.0;
	t64 = c*r1*r2*t1*2.0;
	t65 = c*r2*r3*t3*2.0;
	t70 = c*e*4.0;
	t71 = b*t3*4.0;
	t72 = b*t3*t9;
	t73 = c*t2*t7;
	t74 = c*t2*t9;
	t75 = b*d*r2*4.0;
	t76 = c*d*r3*4.0;
	t77 = b*r1*t2*4.0;
	t78 = c*r1*t3*4.0;
	t79 = b*f*t7;
	t80 = b*f*t8;
	t81 = c*e*t8;
	t82 = c*d*r1*r2*2.0;
	t83 = c*f*r2*r3*2.0;
	t84 = b*r1*r3*t1*2.0;
	t85 = b*r2*r3*t2*2.0;
	t66 = t50+t51+t52+t53+t54+t55+t56+t57+t58+t59+t60+t61+t62+t63+t64+t65-t70-t71-t72-t73-t74-t75-t76-t77-t78-t79-t80-t81-t82-t83-t84-t85;
	t67 = b*d*r1*2.0;
	t68 = b*e*r2*2.0;
	t69 = b*f*r3*2.0;
	t86 = b*r2*4.0;
	t87 = b*r1*r3*2.0;
	t88 = b*t8;
	t89 = b*t9;
	t90 = c*r1*t1*2.0;
	t91 = c*r2*t2*2.0;
	t92 = c*r3*t3*2.0;
	t93 = a*f*4.0;
	t94 = c*t1*4.0;
	t95 = a*d*4.0;
	t96 = c*f*4.0;
	t97 = a*f*r2*2.0;
	t98 = c*e*r1*2.0;
	t99 = a*r3*t2*2.0;
	t100 = c*r2*t1*2.0;
	t101 = a*t3*t7;
	t102 = a*t3*t8;
	t103 = c*t1*t7;
	t104 = a*e*r1*4.0;
	t105 = c*e*r3*4.0;
	t106 = a*r2*t1*4.0;
	t107 = c*r2*t3*4.0;
	t108 = c*d*t8;
	t109 = a*f*t9;
	t110 = c*d*t9;
	t111 = a*d*r1*r3*2.0;
	t112 = a*e*r2*r3*2.0;
	t113 = c*r1*r2*t2*2.0;
	t114 = c*r1*r3*t3*2.0;
	t119 = c*d*4.0;
	t120 = a*t3*4.0;
	t121 = a*t3*t9;
	t122 = c*t1*t8;
	t123 = c*t1*t9;
	t124 = a*d*r2*4.0;
	t125 = c*f*r2*4.0;
	t126 = a*r1*t2*4.0;
	t127 = c*r3*t2*4.0;
	t128 = a*f*t7;
	t129 = c*d*t7;
	t130 = a*f*t8;
	t131 = c*e*r1*r2*2.0;
	t132 = c*f*r1*r3*2.0;
	t133 = a*r1*r3*t1*2.0;
	t134 = a*r2*r3*t2*2.0;
	t115 = t93+t94+t101+t102+t103+t104+t105+t106+t107+t108+t109+t110+t111+t112+t113+t114-t119-t120-t121-t122-t123-t124-t125-t126-t127-t128-t129-t130-t131-t132-t133-t134;
	t116 = a*r1*t1*2.0;
	t117 = a*r2*t2*2.0;
	t118 = a*r3*t3*2.0;
	t135 = c*t7;
	t136 = c*t9;
	t137 = a*r1*4.0;
	t138 = c*r3*4.0;
	t139 = c*r1*r2*2.0;
	t140 = a*4.0;
	t141 = a*t8;
	MatrixXd A0(3,6);
	A0(0,0) = t5*t6*t11*(t67+t68+t69+t93-a*t3*4.0-a*d*r2*2.0+a*e*r1*2.0-a*r1*t2*2.0+a*r2*t1*2.0-b*r1*t1*2.0-b*r2*t2*2.0-b*r3*t3*2.0)+r1*t5*t6*t12*t29*2.0;
	A0(0,1) = t5*t6*t11*(t50+t116+t117+t118-b*t3*4.0-a*d*r1*2.0-a*e*r2*2.0-b*d*r2*2.0+b*e*r1*2.0-a*f*r3*2.0-b*r1*t2*2.0+b*r2*t1*2.0)+r2*t5*t6*t12*t29*2.0;
	A0(0,2) = -t5*t6*t11*(t46+t47+t48+t95+t97+t99-a*t1*4.0-b*t2*4.0-a*e*r3*2.0-b*f*r1*2.0-a*r2*t3*2.0-b*r3*t1*2.0)+r3*t5*t6*t12*t29*2.0;
	A0(0,3) = t5*t6*t11*(b*-4.0+t88+t89+a*r3*4.0-b*t7+a*r1*r2*2.0);
	A0(0,4) = t5*t6*t11*(t140+t141+b*r3*4.0-a*t7-a*t9-b*r1*r2*2.0);
	A0(0,5) = -t5*t6*t11*(t86+t87+t137-a*r2*r3*2.0);
	A0(1,0) = -t6*t11*t49*(t46+t47+t48+t96+t98+t100-b*t2*4.0-c*t3*4.0-c*d*r2*2.0-b*f*r1*2.0-b*r3*t1*2.0-c*r1*t2*2.0)+r1*t6*t12*t49*t66*2.0;
	A0(1,1) = -t6*t11*t49*(t14-t30+t90+t91+t92-c*d*r1*2.0+b*e*r3*2.0-b*f*r2*2.0-c*e*r2*2.0-c*f*r3*2.0+b*r2*t3*2.0-b*r3*t2*2.0)+r2*t6*t12*t49*t66*2.0;
	A0(1,2) = -t6*t11*t49*(t67+t68+t69+t94-c*d*4.0+c*e*r3*2.0-c*f*r2*2.0-b*r1*t1*2.0-b*r2*t2*2.0-b*r3*t3*2.0+c*r2*t3*2.0-c*r3*t2*2.0)+r3*t6*t12*t49*t66*2.0;
	A0(1,3) = -t6*t11*t49*(t86-t87+t138+t139);
	A0(1,4) = t6*t11*t49*(c*-4.0+t135+t136+b*r1*4.0-c*t8+b*r2*r3*2.0);
	A0(1,5) = t6*t11*t49*(b*4.0-t88+t89+c*r1*4.0-b*t7-c*r2*r3*2.0);
	A0(2,0) = t5*t11*t49*(t13-t31+t90+t91+t92+a*d*r3*2.0-a*f*r1*2.0-c*d*r1*2.0-c*e*r2*2.0-c*f*r3*2.0+a*r1*t3*2.0-a*r3*t1*2.0)-r1*t5*t12*t49*t115*2.0;
	A0(2,1) = -t5*t11*t49*(t95+t96+t97+t98+t99+t100-a*t1*4.0-c*t3*4.0-a*e*r3*2.0-c*d*r2*2.0-a*r2*t3*2.0-c*r1*t2*2.0)-r2*t5*t12*t49*t115*2.0;
	A0(2,2) = -t5*t11*t49*(t51-t70+t116+t117+t118-a*d*r1*2.0-a*e*r2*2.0-a*f*r3*2.0-c*d*r3*2.0+c*f*r1*2.0-c*r1*t3*2.0+c*r3*t1*2.0)-r3*t5*t12*t49*t115*2.0;
	A0(2,3) = t5*t11*t49*(c*4.0+t135-t136+a*r2*4.0-c*t8-a*r1*r3*2.0);
	A0(2,4) = -t5*t11*t49*(t137+t138-t139+a*r2*r3*2.0);
	A0(2,5) = t5*t11*t49*(-t140+t141+c*r2*4.0+a*t7-a*t9+c*r1*r3*2.0);
	return A0;

}




MatrixXd functionForRotationAndTransition(const vector<MatrixXd>& input)
{
	MatrixXd parameters(1,6);
	parameters.block(0,0,1,3)=input[0];
	parameters.block(0,3,1,3)=input[2];
	MatrixXd variables=input[1];
	return functionForRotationAndTransition(parameters,variables);
}

MatrixXd jacobianForRotationAndTransition(const vector<MatrixXd>& input)
{
	MatrixXd parameters(1,6);
	parameters.block(0,0,1,3)=input[0];
	parameters.block(0,3,1,3)=input[2];
	MatrixXd variables=input[1];
	return jacobianForRotationAndTransition(parameters,variables);

}

MatrixXd jacobianForPoint(const vector<MatrixXd>& input)
{
	MatrixXd parameters(1,6);
	parameters.block(0,0,1,3)=input[0];
	parameters.block(0,3,1,3)=input[2];
	MatrixXd variables=input[1];
	return jacobianForPoint(parameters,variables);

}



MatrixXd functionForRotationAndTransitionUnitLength(const vector<MatrixXd>& input)
{
	MatrixXd parameters(1,6);
	parameters.block(0,0,1,3)=input[0];
	parameters.block(0,3,1,3)=input[2];
	MatrixXd variables=input[1];
	return functionForRotationAndTransitionUnitLength(parameters,variables);
}

MatrixXd jacobianForRotationAndTransitionUnitLength(const vector<MatrixXd>& input)
{
	MatrixXd parameters(1,6);
	parameters.block(0,0,1,3)=input[0];
	parameters.block(0,3,1,3)=input[2];
	MatrixXd variables=input[1];
	return jacobianForRotationAndTransitionUnitLength(parameters,variables);

}

MatrixXd jacobianForPointUnitLength(const vector<MatrixXd>& input)
{
	MatrixXd parameters(1,6);
	parameters.block(0,0,1,3)=input[0];
	parameters.block(0,3,1,3)=input[2];
	MatrixXd variables=input[1];
	return jacobianForPointUnitLength(parameters,variables);

}
/*
MatrixXd estimateCameraParameter(const MatrixXd& projPoints,const MatrixXd& points)
{
	assert(projPoints.rows()==points.rows());
	MatrixXd dataset(projPoints.rows(),6);

	dataset.block(0,0,projPoints.rows(),3)=projPoints;
	dataset.block(0,3,projPoints.rows(),3)=points;
	MatrixXd obj_vals=MatrixXd::Zero(1,projPoints.rows()*3);
	MatrixXd para=MatrixXd::Zero(1,6);
	return	levenbergM_simple(dataset,obj_vals,functionForRotationAndTransition,jacobianForRotationAndTransition,para);
}*/


MatrixXd estimateCameraParameter(const MatrixXd& projPoints,const vector<int>& ind1,const MatrixXd& points,const vector<int>& ind2)
{
	assert(ind1.size()==ind2.size());

	MatrixXd dataset(ind1.size(),6);

	for (int i = 0; i < ind1.size(); i++)
	{
		dataset.block(i,0,1,3)=projPoints.row(ind1[i]);
		dataset.block(i,3,1,3)=points.row(ind2[i]);
	}

	
	MatrixXd obj_vals=MatrixXd::Zero(1,ind1.size()*3);
	MatrixXd para=MatrixXd::Zero(1,6);
	return	levenbergM_simple(dataset,obj_vals,functionForRotationAndTransition,jacobianForRotationAndTransition,para);
}
#include "sphericalStructureFromMotion.h"


template<class Pnt>
   vector< unordered_map<int,int> > incrementalTrajectoryDetect(const vector<vector<Pnt> >& features, vector<map<int,int>>& twoFrameCorrespondences)
{



	vector<unordered_map<int,int> > result;

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
		
				unordered_map<int,int> pntme;
			
				pntme[i]=j;
				markers[i][j]=true;
				
				int cur_frame=i;
				int cur_index=j;
		
				while( (cur_frame<features.size()-1) && (twoFrameCorrespondences[cur_frame].count(cur_index)))
				{
					int indd=twoFrameCorrespondences[cur_frame][cur_index];
					pntme[cur_frame+1]=indd;
					markers[cur_frame+1][indd]=true;
				
					cur_index=indd;
					++cur_frame;			
				}
				if(pntme.size()>1)
				{

					result.push_back(pntme);
				}
			}
		}
	}
	
	
	cout<<"feature tracing finished"<<endl;
	return result;
}



auto collectFromPairwiseMatching(const string& featureLst,const string& matchLst,vector<vector<vector<int>>>& feas,vector<map<int,int> >& sth,vector<unordered_set<int> >& contain,vector<vector<int> >& featureIsPoint)->  vector<unordered_map<int,int> >
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

	   vector< unordered_map<int,int> > trajs=incrementalTrajectoryDetect(feas,sth);
	
	contain.resize(feas.size());

	for (int i = 0; i < trajs.size(); i++)
	{

		for (auto &s:trajs[i])
		{
			contain[s.first].insert(i);
			featureIsPoint[s.first][s.second]=i;
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



MatrixXd pointProjectError(const MatrixXd& pU,const MatrixXd& pnt)
{
	double p1,p2,p3,u1,u2,u3;
	double x,y,z;

	p1=pU(0,0);
	p2=pU(0,1);
	p3=pU(0,2);
	
	u1=pU(0,3);
	u2=pU(0,4);
	u3=pU(0,5);

	x=pnt(0,0);
	y=pnt(0,1);
	z=pnt(0,2);

	MatrixXd A0(1,3);
	double t2 ,t3 ,t4 ,t5 ,t6 ,t7 ,t8 ,t9 ;
	t2 = p1-x;
	t3 = p2-y;
	t4 = p3-z;
	t5 = t2*t2;
	t6 = t3*t3;
	t7 = t4*t4;
	t8 = t5+t6+t7;
	t9 = 1.0/sqrt(t8);
	A0(0,0) = -u1-t2*t9;
	A0(0,1) = -u2-t3*t9;
	A0(0,2) = -u3-t4*t9;

	return A0;
}



MatrixXd pointProjectErrorJac0bian(const MatrixXd& pU,const MatrixXd& pnt)
{
	double p1,p2,p3,u1,u2,u3;
	double x,y,z;

	p1=pU(0,0);
	p2=pU(0,1);
	p3=pU(0,2);
	
	u1=pU(0,3);
	u2=pU(0,4);
	u3=pU(0,5);

	x=pnt(0,0);
	y=pnt(0,1);
	z=pnt(0,2);

	MatrixXd A0(3,3);
	double t2 ,t3 ,t4 ,t5 ,t6 ,t7 ,t8 ,t9 ,t10 ,t17 ,t11 ,t12 ,t13 ,t18 ,t14 ,t15 ,t19 ,t16 ;
	t2 = p1-x;
	t3 = p2-y;
	t4 = p3-z;
	t5 = t2*t2;
	t6 = t3*t3;
	t7 = t4*t4;
	t8 = t5+t6+t7;
	t9 = 1.0/pow(t8,3.0/2.0);
	t10 = p1*2.0;
	t17 = x*2.0;
	t11 = t10-t17;
	t12 = 1.0/sqrt(t8);
	t13 = p2*2.0;
	t18 = y*2.0;
	t14 = t13-t18;
	t15 = p3*2.0;
	t19 = z*2.0;
	t16 = t15-t19;
	A0(0,0) = t12-t2*t9*t11*(1.0/2.0);
	A0(0,1) = t2*t9*t14*(-1.0/2.0);
	A0(0,2) = t2*t9*t16*(-1.0/2.0);
	A0(1,0) = t3*t9*t11*(-1.0/2.0);
	A0(1,1) = t12-t3*t9*t14*(1.0/2.0);
	A0(1,2) = t3*t9*t16*(-1.0/2.0);
	A0(2,0) = t4*t9*t11*(-1.0/2.0);
	A0(2,1) = t4*t9*t14*(-1.0/2.0);
	A0(2,2) = t12-t4*t9*t16*(1.0/2.0);

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

	MatrixXd dataset(p.rows(),6);
	
	dataset.block(0,0,p.rows(),p.cols())=p;
	dataset.block(0,p.cols(),u.rows(),u.cols())=u;

	MatrixXd obj_vals=MatrixXd(1,3*p.rows());

	//cout<<"\n"<<dataset<<endl;
	auto good=levenbergM_simple(dataset,obj_vals,pointProjectError,pointProjectErrorJac0bian,result);
	//cout<<"from "<<result<<" become "<<good<<endl;
	//getchar();
	return good;

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

		MatrixXd pU(curP.rows(),curP.cols()+curU.cols());

		pU.block(0,0,curP.rows(),curP.cols())=curP;
		pU.block(0,curP.cols(),curU.rows(),curU.cols())=curU;
		if(pointBeforeCamera(points.row(i),curP.row(0),curU.row(0)) && pointBeforeCamera(points.row(i),curP.row(1),curU.row(1)))
		{
	
			//goodlabel[i]=true;
			//error[i]=distanceBetweenPointLine(points.row(i),curP.row(0),curU.row(0))+distanceBetweenPointLine(points.row(i),curP.row(1),curU.row(1));
			error[i]=0;
			for(int pui=0;pui<pU.rows();++pui)
			{
				MatrixXd terror=pointProjectError( pU.row(pui),points.row(i));
				error[i]+= (terror*terror.transpose())(0,0);
			}
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

    vector< unordered_map<int,int> > trajectories=collectFromPairwiseMatching(featureFileName,matchFileName,features,correspondences,contain,featureIsPoint);

	vector<unordered_map<int,bool> > projectionsValid(trajectories.size());

	for ( int i=0;i<trajectories.size();++i)
	{
		for(auto&w:trajectories[i])
			projectionsValid[i][w.first]=false;

	}

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
			if (mpointsError2Frame.second[i]>0 && mpointsError2Frame.second[i]<constrain_on_point_error &&length(mpointsError2Frame.first.row(i))<constrain_on_point_length )
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
			if(!alreadyReconstructed[featureIsPoint[lcameraIndex-1][s.first]])
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
						for (auto &mi: trajectories[w])
						{
							if(lbundlePara.count(mi.first))
							{
								++prjCount;
								curPrj.insert(mi.first);
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
							//	int curPrjCamera=trajectories[w].first[x];
								int curPrjCamera=x;
								



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
								allProjs.push_back(make_pair( x,trajectories[w][x]));
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
		//getchar();
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
#include "sphericalStructureFromMotion.h"





//typedef vector<int> Pnt;

template<class Pnt>
auto incrementalTrajectoryDetect(const vector<vector<Pnt> >& features, vector<map<int,int>>& sth)-> pair<vector<pair<int,int> >,  vector<vector<int> > >
{

	vector<pair<int,int> > sttend;
	vector<vector<int> > trajs;
	vector<vector<bool> > markers(features.size());	

	for (size_t i = 0; i < features.size(); i++)
	{
		markers[i].resize(features[i].size(),false);
	}

	for (int i = 0; i < features.size()-1; i++)
	{
		cout<<"processing frame " <<i<<endl;
		for (int j = 0; j < markers[i].size(); j++)
		{
			if(!markers[i][j])
			{
				pair<int,int> stend;
				vector<int> traj;
		
				stend.first=i;
				traj.push_back(j);
				markers[i][j]=true;
				
				int cur_frame=i;
				int cur_index=j;
		
				while( (cur_frame<features.size()-1) && (sth[cur_frame].count(cur_index)))
				{
					int indd=sth[cur_frame][cur_index];
				
					traj.push_back(indd);

					markers[cur_frame+1][indd]=true;
				
					cur_index=indd;
					++cur_frame;			
				}
				stend.second=cur_frame;
				if(stend.second>stend.first)
				{
					trajs.push_back(traj);
					sttend.push_back(stend);}
				}
		}
	}
	cout<<"feature tracing finished"<<endl;
	return make_pair(sttend,trajs);
}



auto collectFromPairwiseMatching(const string& featureLst,const string& matchLst,vector<vector<vector<int>>>& feas,vector<map<int,int> >& sth,vector<unordered_set<int> >& contain,vector<vector<int> >& featureIsPoint)-> pair<vector<pair<int,int> >,  vector<vector<int> > >
{

	auto feaNames=fileIOclass::InVectorString(featureLst);
	auto mathces=fileIOclass::InVectorString(matchLst);

	vector<vector<vector<int> > > correpss(mathces.size());
	//vector<vector<vector<int> > > 
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


	//vector<map<int,int> > 
	sth.resize(correpss.size());

	for(int i=0;i<correpss.size();++i)
	{
		for(int j=0;j<correpss[i].size();++j)
			sth[i][correpss[i][j][0]]=correpss[i][j][1];
	}

	auto trajs=incrementalTrajectoryDetect(feas,sth);
	
	contain.resize(feas.size());

	for (int i = 0; i < trajs.first.size(); i++)
	{
		for (int j = trajs.first[i].first; j <= trajs.first[i].second; j++)
		{
			contain[j].insert(i);
			featureIsPoint[j][trajs.second[i][j-trajs.first[i].first]]=i;
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

	//cout<<essential<<endl;

	MatrixXd V=	svd.matrixU();
	MatrixXd U=	svd.matrixV();

	cout<<"V\n"<<endl;

	cout<<V<<endl;

	cout<<"U\n"<<endl;

	cout<<U<<endl;
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

	cout<<Tx<<endl;
	cout<<R1<<endl;
	cout<<R2<<endl;
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
	return result;

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

pair<MatrixXd,vector<bool> > bestPoints(const MatrixXd& spnts1,const MatrixXd& spnts2,const vector<MatrixXd>& transitions,const vector<MatrixXd>& rotations)
{
	assert(spnts1.rows()==spnts2.rows());

	MatrixXd points(spnts1.rows(),3);
	vector<bool> goodlabel(spnts1.rows(),false);

	MatrixXd curP(2,3);
	curP.row(0)=transitions[0];
	curP.row(1)=transitions[1];
	for (int i = 0; i < spnts1.rows(); i++)
	{
		
		MatrixXd curU(2,3);
		curU.row(0)=(rotations[0]*spnts1.row(i).transpose()).transpose();		
		curU.row(1)=(rotations[1]*spnts2.row(i).transpose()).transpose();
		points.row(i)=bestPoint(curP,curU);

	/*	if(angleBetween(points.row(i)-curP.row(0),curU.row(0))<constrain_on_goodPoint && angleBetween(points.row(i)-curP.row(1),curU.row(1))<constrain_on_goodPoint )
		{
	
			goodlabel[i]=true;
		}*/

		if(pointBeforeCamera(points.row(i),curP.row(0),curU.row(0)) && pointBeforeCamera(points.row(i),curP.row(1),curU.row(1)))
		{
	
			goodlabel[i]=true;
		}
	}

	return make_pair(points,goodlabel);
}





pair<vector<vector<double> >,vector<double> > geometricReconstructionFrom2Frames(const MatrixXd& sphericalPoints1,const MatrixXd& sphericalPoints2,vector<double>& transition,vector<double >& rotation)
{
	vector<vector<double> > points(sphericalPoints1.rows(),vector<double>(3,0.0));
	vector<double> errors(sphericalPoints1.rows(),-1);


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

		MatrixXd tob=sphericalPoints1.row(i).transpose()*sphericalPoints2.row(i);

		observation.row(i)=convert(tob);
	}

	cout<<observation<<endl<<endl;

	JacobiSVD<decltype(observation)> svd(observation, Eigen::ComputeFullU |
                                        Eigen::ComputeFullV);
	std::cout << "matrix V" << endl << svd.matrixV() <<"\n"<< endl;

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

	auto countBool=[](const vector<bool>& toc)->int
	{
		int sum=0;
		for(auto s:toc)
			if(s)
				++sum;

		return sum;
	};

	int bestIndex;
	int bestCount;
	vector<int> goodPointCount(transtionAndRotations.size(),0);

	cout<<"transition:\n"<<endl;
	for (int i = 0; i < transtionAndRotations.size(); i++)
	{
		cout<<transtionAndRotations[i].first<<endl;
		vector<MatrixXd> trans(2,MatrixXd::Zero(1,3));
		trans[1]=transtionAndRotations[i].first;
		vector<MatrixXd> rots(2,MatrixXd::Identity(3,3));
		rots[1]=transtionAndRotations[i].second;
		bestPointCandidates[i]=bestPoints(sphericalPoints1,sphericalPoints2,trans,rots);
		goodPointCount[i]=countBool(bestPointCandidates[i].second);

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

//	bestIndex=3;
	cout<<"the best index:"<<bestIndex<<endl;
	cout<<"the best transition:\n"<<transtionAndRotations[bestIndex].first <<endl;
	cout<<"the best rotation:\n"<<transtionAndRotations[bestIndex].second<<endl;

	bestIndex=6;
	cout<<"the best index:"<<bestIndex<<endl;
	cout<<"the best transition:\n"<<transtionAndRotations[bestIndex].first <<endl;
	cout<<"the best rotation:\n"<<transtionAndRotations[bestIndex].second<<endl;


	//cout<<observation;
	return make_pair(points,errors);

}
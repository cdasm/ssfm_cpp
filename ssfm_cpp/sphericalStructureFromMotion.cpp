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

pair<MatrixXd,vector<double> > bestPoints(const MatrixXd& spnts1,const MatrixXd& spnts2,const vector<MatrixXd>& transitions,const vector<MatrixXd>& rotations)
{
	assert(spnts1.rows()==spnts2.rows());

	MatrixXd points(spnts1.rows(),3);
	vector<double> error(spnts1.rows(),-1);

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
	
			//goodlabel[i]=true;
			error[i]=distanceBetweenPointLine(points.row(i),curP.row(0),curU.row(0))+distanceBetweenPointLine(points.row(i),curP.row(1),curU.row(1));
		}
	}

	return make_pair(points,error);
}





auto geometricReconstructionFrom2Frames(const MatrixXd& sphericalPoints1,const MatrixXd& sphericalPoints2,vector<double>& transition,vector<double >& rotation)->pair<MatrixXd ,vector<double> >
{
//	vector<vector<double> > points(sphericalPoints1.rows(),vector<double>(3,0.0));
//	vector<double> errors(sphericalPoints1.rows(),-1);


	MatrixXd observation(sphericalPoints1.rows(),9);

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
	
	for (int i = 0; i < sphericalPoints1.rows(); i++)
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

	cout<<"transition:\n"<<endl;
	for (int i = 0; i < transtionAndRotations.size(); i++)
	{
		cout<<transtionAndRotations[i].first<<endl;
		vector<MatrixXd> trans(2,MatrixXd::Zero(1,3));
		trans[1]=transtionAndRotations[i].first;
		vector<MatrixXd> rots(2,MatrixXd::Identity(3,3));
		rots[1]=transtionAndRotations[i].second;
		bestPointCandidates[i]=bestPoints(sphericalPoints1,sphericalPoints2,trans,rots);
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

//	bestIndex=3;
	cout<<"the best index:"<<bestIndex<<endl;
	cout<<"the best transition:\n"<<transtionAndRotations[bestIndex].first <<endl;
	cout<<"the best rotation:\n"<<transtionAndRotations[bestIndex].second<<endl;

	


	//cout<<observation;
	return bestPointCandidates[bestIndex];

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



MatrixXd functionForRotationAndTransition2(const MatrixXd& parameters,const MatrixXd& variables2)
{
	MatrixXd variable(1,6);
	variable.block(0,0,1,3)=variables2.block(0,0,1,3);
	variable.block(0,3,1,3)=transitionFrom2Para(variables2.block(0,3,1,2));

	return functionForRotationAndTransition(parameters,variable);
}

MatrixXd jacobianForPoint2(const MatrixXd& parameters,const MatrixXd& variables2)
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

MatrixXd jacobianForRotationAndTransition2(const MatrixXd& parameters,const MatrixXd& variables)
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



MatrixXd functionForRotationAndTransition2(const vector<MatrixXd>& input)
{
	MatrixXd parameters(1,6);
	parameters.block(0,0,1,3)=input[0];
	parameters.block(0,3,1,3)=input[2];
	MatrixXd variables=input[1];
	return functionForRotationAndTransition2(parameters,variables);
}

MatrixXd jacobianForRotationAndTransition2(const vector<MatrixXd>& input)
{
	MatrixXd parameters(1,6);
	parameters.block(0,0,1,3)=input[0];
	parameters.block(0,3,1,3)=input[2];
	MatrixXd variables=input[1];
	return jacobianForRotationAndTransition2(parameters,variables);

}

MatrixXd jacobianForPoint2(const vector<MatrixXd>& input)
{
	MatrixXd parameters(1,6);
	parameters.block(0,0,1,3)=input[0];
	parameters.block(0,3,1,3)=input[2];
	MatrixXd variables=input[1];
	return jacobianForPoint2(parameters,variables);

}
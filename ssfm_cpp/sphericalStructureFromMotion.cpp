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

MatrixXd pointProjectError_(const MatrixXd& pU,const MatrixXd& pnt)
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
	double t2 ,t3 ,t4 ,t5 ,t6 ,t7 ,t8 ,t9 ,t10 ,t11 ,t12 ,t13 ;
	t2 = p1-x;
	t3 = p2-y;
	t4 = p3-z;
	t5 = 1.0/u2;
	t6 = t3*t5;
	t7 = t2*t2;
	t8 = t3*t3;
	t9 = t4*t4;
	t10 = t7+t8+t9;
	t11 = 1.0/sqrt(t10);
	t12 = 1.0/u1;
	t13 = 1.0/u3;
	A0(0,0) = t11*(t6-t2*t12);
	A0(0,1) = -t11*(t6-t4*t13);
	A0(0,2) = t11*(t2*t12-t4*t13);

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
	t2 = 1.0/u2;
	t3 = p2-y;
	t4 = t2*t3;
	t5 = 1.0/u1;
	t6 = p1-x;
	t7 = 1.0/u3;
	t8 = p3-z;
	t9 = t7*t8;
	A0(0,0) = t4-t5*t6;
	A0(0,1) = -t4+t9;
	A0(0,2) = -t9+t5*t6;

	return A0;
}


MatrixXd pointProjectErrorJac0bian_(const MatrixXd& pU,const MatrixXd& pnt)
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
	double t2 ,t3 ,t4 ;
	t2 = 1.0/u2;
	t3 = 1.0/u1;
	t4 = 1.0/u3;
	A0(0,0) = t3;
	A0(0,1) = -t2;
	A0(1,1) = t2;
	A0(1,2) = -t4;
	A0(2,0) = -t3;
	A0(2,2) = t4;
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
	double t2 ,t3 ,t4 ,t5 ,t6 ,t7 ,t8 ,t9 ,t10 ,t11 ,t12 ,t15 ,t13 ,t14 ,t16 ,t26 ,t17 ,t18 ,t27 ,t19 ,t20 ,t24 ,t21 ,t22 ,t29 ,t23 ,t25 ,t28 ;
	t2 = p1-x;
	t3 = p2-y;
	t4 = p3-z;
	t5 = 1.0/u1;
	t6 = t2*t2;
	t7 = t3*t3;
	t8 = t4*t4;
	t9 = t6+t7+t8;
	t10 = 1.0/u2;
	t11 = 1.0/sqrt(t9);
	t12 = t2*t5;
	t15 = t3*t10;
	t13 = t12-t15;
	t14 = 1.0/pow(t9,3.0/2.0);
	t16 = p1*2.0;
	t26 = x*2.0;
	t17 = t16-t26;
	t18 = p2*2.0;
	t27 = y*2.0;
	t19 = t18-t27;
	t20 = 1.0/u3;
	t24 = t4*t20;
	t21 = t15-t24;
	t22 = p3*2.0;
	t29 = z*2.0;
	t23 = t22-t29;
	t25 = t5*t11;
	t28 = t12-t24;
	A0(0,0) = t25-t13*t14*t17*(1.0/2.0);
	A0(0,1) = -t10*t11-t13*t14*t19*(1.0/2.0);
	A0(0,2) = t13*t14*t23*(-1.0/2.0);
	A0(1,0) = t14*t17*t21*(-1.0/2.0);
	A0(1,1) = t10*t11-t14*t19*t21*(1.0/2.0);
	A0(1,2) = -t11*t20-t14*t21*t23*(1.0/2.0);
	A0(2,0) = -t25+t14*t17*t28*(1.0/2.0);
	A0(2,1) = t14*t19*t28*(1.0/2.0);
	A0(2,2) = t11*t20+t14*t23*t28*(1.0/2.0);

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

	auto good=levenbergM_simple(dataset,obj_vals,pointProjectError,pointProjectErrorJac0bian,result);
	cout<<"from "<<result<<" become "<<good<<endl;
	getchar();
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

			MatrixXd terror=pointProjectError( pU,points.row(i));
			error[i]= (terror*terror.transpose())(0,0);
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
	getchar();
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

MatrixXd functionForRotationAndTransition__(const MatrixXd& parameters,const MatrixXd& variables)
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
	double t5 ,t6 ,t7 ,t8 ,t9 ,t10 ,t11 ,t12 ,t13 ,t14 ,t15 ,t16 ,t17 ,t18 ,t19 ,t20 ,t21 ,t22 ,t31 ,t23 ,t24 ,t32 ,t25 ,t26 ,t27 ,t28 ,t29 ,t30 ,t33 ,t34 ,t35 ,t36 ,t37 ,t38 ,t39 ,t40 ,t41 ,t42 ,t43 ,t44 ,t45 ,t46 ,t47 ,t48 ;
	t5 = r1*r1;
	t6 = t5*(1.0/4.0);
	t7 = r2*r2;
	t8 = t7*(1.0/4.0);
	t9 = r3*r3;
	t10 = t9*(1.0/4.0);
	t11 = t6+t8+t10+1.0;
	t12 = 1.0/t11;
	t13 = r3*t12;
	t14 = d-t1;
	t15 = f-t3;
	t16 = t6+t8+t10-1.0;
	t17 = e-t2;
	t18 = 1.0/b;
	t19 = r1*r2*t12*(1.0/2.0);
	t20 = t13+t19;
	t21 = t14*t20;
	t22 = r1*t12;
	t31 = r2*r3*t12*(1.0/2.0);
	t23 = t22-t31;
	t24 = t7*t12*(1.0/2.0);
	t32 = t12*t16;
	t25 = t24-t32;
	t26 = t17*t25;
	t27 = t21+t26-t15*t23;
	t28 = t18*t27;
	t29 = r2*t12;
	t30 = r1*r3*t12*(1.0/2.0);
	t33 = t14*t14;
	t34 = t17*t17;
	t35 = t15*t15;
	t36 = t33+t34+t35;
	t37 = 1.0/sqrt(t36);
	t38 = 1.0/a;
	t39 = t29+t30;
	t40 = t15*t39;
	t41 = t5*t12*(1.0/2.0);
	t42 = 1.0/c;
	t43 = t29-t30;
	t44 = t14*t43;
	t45 = t22+t31;
	t46 = t32-t9*t12*(1.0/2.0);
	t47 = t15*t46;
	t48 = t42*(t44+t47-t17*t45);
	A0(0,0) = -t37*(t28-t38*(t40-t17*(t13-r1*r2*t12*(1.0/2.0))+t14*(t41-t12*t16)));
	A0(0,1) = t37*(t28+t48);
	A0(0,2) = -t37*(t48-t38*(-t40+t17*(t13-t19)+t14*(t32-t41)));
	return A0;
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



MatrixXd functionForRotationAndTransitionUnitLength(const MatrixXd& parameters,const MatrixXd& variables2)
{
	MatrixXd variable(1,6);
	variable.block(0,0,1,3)=variables2.block(0,0,1,3);
	variable.block(0,3,1,3)=transitionFrom2Para(variables2.block(0,3,1,2));

	return functionForRotationAndTransition(parameters,variable);
}


MatrixXd functionForRotationAndTransitionUnitLength__(const MatrixXd& parameters,const MatrixXd& variables2)
{
	MatrixXd variable(1,6);
	variable.block(0,0,1,3)=variables2.block(0,0,1,3);
	variable.block(0,3,1,3)=transitionFrom2Para(variables2.block(0,3,1,2));

	return functionForRotationAndTransition(parameters,variable);
}

MatrixXd jacobianForPointUnitLength(const MatrixXd& parameters,const MatrixXd& variables)
{
	double a,b,c,d,e,f,t1,t200,r1,r2,r3;
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
	t200=variables(0,4) ;
	MatrixXd A0(3,3);
	double t2 ,t3 ,t4 ,t5 ,t6 ,t7 ,t8 ,t9 ,t10 ,t11 ,t12 ,t13 ,t14 ,t15 ,t16 ,t17 ,t18 ,t19 ,t20 ;
	t2 = r1*r1;
	t3 = r2*r2;
	t4 = r3*r3;
	t5 = 1.0/a;
	t6 = 1.0/b;
	t7 = t2+t3+t4+4.0;
	t8 = 1.0/t7;
	t9 = b*r2*2.0;
	t10 = b*r1*r3;
	t11 = 1.0/c;
	t12 = b*t3;
	t13 = b*t4;
	t14 = c*t2;
	t15 = c*t4;
	t16 = a*r1*2.0;
	t17 = c*r3*2.0;
	t18 = c*r1*r2;
	t19 = a*4.0;
	t20 = a*t3;
	A0(0,0) = -t5*t6*t8*(b*-4.0+t12+t13+a*r3*4.0-b*t2+a*r1*r2*2.0);
	A0(0,1) = -t5*t6*t8*(t19+t20+b*r3*4.0-a*t2-a*t4-b*r1*r2*2.0);
	A0(0,2) = t5*t6*t8*(t9+t10+t16-a*r2*r3)*2.0;
	A0(1,0) = t6*t8*t11*(t9-t10+t17+t18)*2.0;
	A0(1,1) = -t6*t8*t11*(c*-4.0+t14+t15+b*r1*4.0-c*t3+b*r2*r3*2.0);
	A0(1,2) = -t6*t8*t11*(b*4.0-t12+t13+c*r1*4.0-b*t2-c*r2*r3*2.0);
	A0(2,0) = -t5*t8*t11*(c*4.0+t14-t15+a*r2*4.0-c*t3-a*r1*r3*2.0);
	A0(2,1) = t5*t8*t11*(t16+t17-t18+a*r2*r3)*2.0;
	A0(2,2) = -t5*t8*t11*(-t19+t20+a*t2+c*r2*4.0-a*t4+c*r1*r3*2.0);
	return A0;
}

MatrixXd jacobianForPointUnitLength__(const MatrixXd& parameters,const MatrixXd& variables)
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
	MatrixXd A0(3,3);
	double t4 ,t5 ,t6 ,t7 ,t8 ,t9 ,t10 ,t11 ,t12 ,t18 ,t13 ,t20 ,t14 ,t15 ,t16 ,t17 ,t19 ,t21 ,t22 ,t23 ,t24 ,t25 ,t26 ,t27 ,t28 ,t29 ,t30 ,t31 ,t45 ,t46 ,t47 ,t48 ,t49 ,t50 ,t51 ,t52 ,t32 ,t33 ,t34 ,t35 ,t36 ,t37 ,t38 ,t39 ,t40 ,t54 ,t55 ,t56 ,t57 ,t58 ,t59 ,t60 ,t61 ,t41 ,t53 ,t62 ,t42 ,t43 ,t44 ,t63 ,t94 ,t64 ,t65 ,t66 ,t67 ,t68 ,t69 ,t70 ,t71 ,t72 ,t73 ,t74 ,t75 ,t76 ,t81 ,t82 ,t83 ,t84 ,t85 ,t86 ,t87 ,t77 ,t88 ,t78 ,t79 ,t98 ,t80 ,t89 ,t90 ,t91 ,t92 ,t93 ,t95 ,t96 ,t97 ,t99 ,t100 ,t101 ,t102 ,t103 ;
	t4 = cos(t1);
	t5 = cos(t2);
	t6 = r1*r1;
	t7 = r2*r2;
	t8 = r3*r3;
	t9 = sin(t1);
	t10 = sin(t2);
	t11 = t6+t7+t8+4.0;
	t12 = 1.0/t11;
	t18 = t4*t5;
	t13 = d-t18;
	t20 = t4*t10;
	t14 = e-t20;
	t15 = f-t9;
	t16 = 1.0/a;
	t17 = 1.0/b;
	t19 = t13*t13;
	t21 = t14*t14;
	t22 = t15*t15;
	t23 = t19+t21+t22;
	t24 = d*4.0;
	t25 = f*r2*4.0;
	t26 = d*t6;
	t27 = r3*t4*t10*4.0;
	t28 = e*r1*r2*2.0;
	t29 = f*r1*r3*2.0;
	t30 = t4*t5*t7;
	t31 = t4*t5*t8;
	t45 = t4*t5*4.0;
	t46 = e*r3*4.0;
	t47 = d*t7;
	t48 = d*t8;
	t49 = r2*t9*4.0;
	t50 = r1*r3*t9*2.0;
	t51 = t4*t5*t6;
	t52 = r1*r2*t4*t10*2.0;
	t32 = t24+t25+t26+t27+t28+t29+t30+t31-t45-t46-t47-t48-t49-t50-t51-t52;
	t33 = e*4.0;
	t34 = d*r3*4.0;
	t35 = e*t7;
	t36 = r1*t9*4.0;
	t37 = d*r1*r2*2.0;
	t38 = f*r2*r3*2.0;
	t39 = t4*t6*t10;
	t40 = t4*t8*t10;
	t54 = t4*t10*4.0;
	t55 = f*r1*4.0;
	t56 = e*t6;
	t57 = e*t8;
	t58 = r2*r3*t9*2.0;
	t59 = r3*t4*t5*4.0;
	t60 = t4*t7*t10;
	t61 = r1*r2*t4*t5*2.0;
	t41 = t33+t34+t35+t36+t37+t38+t39+t40-t54-t55-t56-t57-t58-t59-t60-t61;
	t53 = t12*t16*t32;
	t62 = t12*t17*t41;
	t42 = t53-t62;
	t43 = 1.0/pow(t23,3.0/2.0);
	t44 = 1.0/sqrt(t23);
	t63 = d*2.0;
	t94 = t4*t5*2.0;
	t64 = t63-t94;
	t65 = 1.0/c;
	t66 = b*r2*2.0;
	t67 = b*r1*r3;
	t68 = f*4.0;
	t69 = t9*4.0;
	t70 = t6*t9;
	t71 = t7*t9;
	t72 = e*r1*4.0;
	t73 = f*t8;
	t74 = r2*t4*t5*4.0;
	t75 = d*r1*r3*2.0;
	t76 = e*r2*r3*2.0;
	t81 = t8*t9;
	t82 = d*r2*4.0;
	t83 = f*t6;
	t84 = f*t7;
	t85 = r1*t4*t10*4.0;
	t86 = r2*r3*t4*t10*2.0;
	t87 = r1*r3*t4*t5*2.0;
	t77 = t68-t69+t70+t71+t72+t73+t74+t75+t76-t81-t82-t83-t84-t85-t86-t87;
	t88 = t12*t65*t77;
	t78 = t62-t88;
	t79 = e*2.0;
	t98 = t4*t10*2.0;
	t80 = t79-t98;
	t89 = f*2.0;
	t90 = t9*2.0;
	t91 = t89-t90;
	t92 = b*t7;
	t93 = b*t8;
	t95 = c*t6;
	t96 = c*t8;
	t97 = t53-t88;
	t99 = a*r1*2.0;
	t100 = c*r3*2.0;
	t101 = c*r1*r2;
	t102 = a*4.0;
	t103 = a*t7;
	A0(0,0) = t42*t43*t64*(-1.0/2.0)-t12*t16*t17*t44*(b*-4.0+t92+t93+a*r3*4.0-b*t6+a*r1*r2*2.0);
	A0(0,1) = t42*t43*t80*(-1.0/2.0)-t12*t16*t17*t44*(t102+t103+b*r3*4.0-a*t6-a*t8-b*r1*r2*2.0);
	A0(0,2) = t42*t43*t91*(-1.0/2.0)+t12*t16*t17*t44*(t66+t67+t99-a*r2*r3)*2.0;
	A0(1,0) = t43*t64*t78*(-1.0/2.0)+t12*t17*t44*t65*(t66-t67+t100+t101)*2.0;
	A0(1,1) = t43*t78*t80*(-1.0/2.0)-t12*t17*t44*t65*(c*-4.0+t95+t96+b*r1*4.0-c*t7+b*r2*r3*2.0);
	A0(1,2) = t43*t78*t91*(-1.0/2.0)-t12*t17*t44*t65*(b*4.0-t92+t93+c*r1*4.0-b*t6-c*r2*r3*2.0);
	A0(2,0) = t43*t64*t97*(1.0/2.0)-t12*t16*t44*t65*(c*4.0+t95-t96+a*r2*4.0-c*t7-a*r1*r3*2.0);
	A0(2,1) = t43*t80*t97*(1.0/2.0)+t12*t16*t44*t65*(t99+t100-t101+a*r2*r3)*2.0;
	A0(2,2) = t43*t91*t97*(1.0/2.0)-t12*t16*t44*t65*(-t102+t103+c*r2*4.0+a*t6-a*t8+c*r1*r3*2.0);

	return A0;
}

MatrixXd jacobianForPoint__(const MatrixXd& parameters,const MatrixXd& variables)
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
	double t5 ,t6 ,t7 ,t8 ,t9 ,t10 ,t11 ,t12 ,t13 ,t14 ,t15 ,t16 ,t17 ,t18 ,t19 ,t20 ,t24 ,t21 ,t22 ,t23 ,t25 ,t26 ,t27 ,t28 ,t29 ,t42 ,t30 ,t31 ,t32 ,t33 ,t34 ,t35 ,t36 ,t48 ,t37 ,t38 ,t39 ,t47 ,t40 ,t41 ,t43 ,t49 ,t50 ,t44 ,t45 ,t46 ,t51 ,t52 ,t53 ,t54 ,t55 ,t56 ,t57 ,t64 ,t58 ,t59 ,t60 ,t61 ,t62 ,t63 ,t65 ,t66 ,t67 ,t68 ,t69 ,t70 ,t71 ,t75 ,t72 ,t73 ,t74 ;
	t5 = r1*r1;
	t6 = t5*(1.0/4.0);
	t7 = r2*r2;
	t8 = t7*(1.0/4.0);
	t9 = r3*r3;
	t10 = t9*(1.0/4.0);
	t11 = t6+t8+t10+1.0;
	t12 = 1.0/t11;
	t13 = d-t1;
	t14 = e-t2;
	t15 = f-t3;
	t16 = 1.0/a;
	t17 = r3*t12;
	t18 = r1*r2*t12*(1.0/2.0);
	t19 = t5*t12*(1.0/2.0);
	t20 = t6+t8+t10-1.0;
	t24 = t12*t20;
	t21 = t19-t24;
	t22 = 1.0/b;
	t23 = t17+t18;
	t25 = t13*t13;
	t26 = t14*t14;
	t27 = t15*t15;
	t28 = t25+t26+t27;
	t29 = t17-t18;
	t42 = t7*t12*(1.0/2.0);
	t30 = t24-t42;
	t31 = 1.0/sqrt(t28);
	t32 = r2*t12;
	t33 = r1*r3*t12*(1.0/2.0);
	t34 = t32+t33;
	t35 = t15*t34;
	t36 = t13*t21;
	t48 = t14*t29;
	t37 = t35+t36-t48;
	t38 = t16*t37;
	t39 = r1*t12;
	t47 = r2*r3*t12*(1.0/2.0);
	t40 = t39-t47;
	t41 = t15*t40;
	t43 = t14*t30;
	t49 = t13*t23;
	t50 = t41+t43-t49;
	t44 = t22*t50;
	t45 = t38+t44;
	t46 = 1.0/pow(t28,3.0/2.0);
	t51 = 1.0/c;
	t52 = t32-t33;
	t53 = d*2.0;
	t54 = t1*2.0;
	t55 = t53-t54;
	t56 = t39+t47;
	t57 = t13*t52;
	t64 = t9*t12*(1.0/2.0);
	t58 = t24-t64;
	t59 = t15*t58;
	t60 = e*2.0;
	t61 = t2*2.0;
	t62 = t60-t61;
	t63 = t22*t40;
	t65 = f*2.0;
	t66 = t3*2.0;
	t67 = t65-t66;
	t68 = t51*t52;
	t69 = t16*t21;
	t70 = t16*t29;
	t71 = t51*t56;
	t75 = t14*t56;
	t72 = t51*(t57+t59-t75);
	t73 = t38+t72;
	t74 = t16*t34;
	A0(0,0) = t31*(t69-t22*t23)-t45*t46*t55*(1.0/2.0);
	A0(0,1) = -t31*(t70-t22*t30)-t45*t46*t62*(1.0/2.0);
	A0(0,2) = t31*(t63+t74)-t45*t46*t67*(1.0/2.0);
	A0(1,0) = t31*(t68+t22*t23)+t46*t55*(t44-t51*(t57+t59-t14*t56))*(1.0/2.0);
	A0(1,1) = -t31*(t71+t22*t30)+t46*t62*(t44-t51*(t57+t59-t14*t56))*(1.0/2.0);
	A0(1,2) = -t31*(t63-t51*t58)+t46*t67*(t44-t51*(t57+t59-t14*t56))*(1.0/2.0);
	A0(2,0) = -t31*(t68+t69)+t46*t55*t73*(1.0/2.0);
	A0(2,1) = t31*(t70+t71)+t46*t62*t73*(1.0/2.0);
	A0(2,2) = -t31*(t74+t51*t58)+t46*t67*t73*(1.0/2.0);

	return A0;
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
MatrixXd jacobianForRotationAndTransitionUnitLength__(const MatrixXd& parameters,const MatrixXd& variables)
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
	double t4 ,t5 ,t6 ,t7 ,t8 ,t9 ,t10 ,t11 ,t12 ,t29 ,t13 ,t31 ,t14 ,t15 ,t16 ,t17 ,t18 ,t19 ,t20 ,t21 ,t22 ,t23 ,t24 ,t25 ,t26 ,t27 ,t28 ,t30 ,t32 ,t33 ,t34 ,t35 ,t36 ,t37 ,t38 ,t39 ,t40 ,t41 ,t42 ,t43 ,t44 ,t182 ,t183 ,t184 ,t185 ,t186 ,t187 ,t188 ,t189 ,t45 ,t46 ,t47 ,t48 ,t49 ,t50 ,t51 ,t52 ,t53 ,t109 ,t110 ,t111 ,t112 ,t113 ,t114 ,t115 ,t116 ,t54 ,t117 ,t190 ,t55 ,t56 ,t57 ,t58 ,t59 ,t60 ,t61 ,t62 ,t63 ,t64 ,t65 ,t66 ,t67 ,t68 ,t69 ,t70 ,t71 ,t72 ,t73 ,t74 ,t75 ,t76 ,t77 ,t78 ,t79 ,t80 ,t81 ,t91 ,t93 ,t94 ,t96 ,t97 ,t146 ,t147 ,t82 ,t83 ,t84 ,t85 ,t86 ,t87 ,t88 ,t89 ,t90 ,t92 ,t95 ,t98 ,t99 ,t100 ,t101 ,t102 ,t103 ,t104 ,t105 ,t106 ,t107 ,t108 ,t118 ,t119 ,t199 ,t120 ,t121 ,t122 ,t123 ,t124 ,t125 ,t126 ,t127 ,t128 ,t129 ,t191 ,t192 ,t193 ,t194 ,t195 ,t196 ,t197 ,t130 ,t198 ,t131 ,t132 ,t201 ,t133 ,t134 ,t135 ,t136 ,t137 ,t138 ,t139 ,t140 ,t141 ,t142 ,t143 ,t174 ,t181 ,t144 ,t145 ,t148 ,t149 ,t150 ,t151 ,t152 ,t153 ,t154 ,t155 ,t156 ,t157 ,t158 ,t159 ,t160 ,t161 ,t162 ,t163 ,t164 ,t165 ,t166 ,t167 ,t168 ,t169 ,t170 ,t171 ,t172 ,t173 ,t175 ,t176 ,t177 ,t178 ,t179 ,t180 ,t200 ;
	t4 = r1*r1;
	t5 = sin(t1);
	t6 = r2*r2;
	t7 = r3*r3;
	t8 = cos(t1);
	t9 = cos(t2);
	t10 = sin(t2);
	t11 = t4+t6+t7+4.0;
	t12 = 1.0/(t11*t11);
	t29 = t8*t9;
	t13 = d-t29;
	t31 = t8*t10;
	t14 = e-t31;
	t15 = f-t5;
	t16 = 1.0/a;
	t17 = f*1.6E1;
	t18 = t5*1.6E1;
	t19 = t4*t5*4.0;
	t20 = f*t6*4.0;
	t21 = f*t7*4.0;
	t22 = f*r1*r2*r3*4.0;
	t23 = 1.0/b;
	t24 = f*r3*8.0;
	t25 = f*r3*t7*2.0;
	t26 = f*r3*t6*2.0;
	t27 = r1*r2*t5*8.0;
	t28 = r3*t4*t5*2.0;
	t30 = t13*t13;
	t32 = t14*t14;
	t33 = t15*t15;
	t34 = t30+t32+t33;
	t35 = 1.0/sqrt(t34);
	t36 = 1.0/t11;
	t37 = d*4.0;
	t38 = f*r2*4.0;
	t39 = d*t4;
	t40 = r3*t8*t10*4.0;
	t41 = e*r1*r2*2.0;
	t42 = f*r1*r3*2.0;
	t43 = t6*t8*t9;
	t44 = t7*t8*t9;
	t182 = t8*t9*4.0;
	t183 = e*r3*4.0;
	t184 = d*t6;
	t185 = d*t7;
	t186 = r2*t5*4.0;
	t187 = r1*r3*t5*2.0;
	t188 = t4*t8*t9;
	t189 = r1*r2*t8*t10*2.0;
	t45 = t37+t38+t39+t40+t41+t42+t43+t44-t182-t183-t184-t185-t186-t187-t188-t189;
	t46 = e*4.0;
	t47 = d*r3*4.0;
	t48 = e*t6;
	t49 = r1*t5*4.0;
	t50 = d*r1*r2*2.0;
	t51 = f*r2*r3*2.0;
	t52 = t4*t8*t10;
	t53 = t7*t8*t10;
	t109 = t8*t10*4.0;
	t110 = f*r1*4.0;
	t111 = e*t4;
	t112 = e*t7;
	t113 = r2*r3*t5*2.0;
	t114 = r3*t8*t9*4.0;
	t115 = t6*t8*t10;
	t116 = r1*r2*t8*t9*2.0;
	t54 = t46+t47+t48+t49+t50+t51+t52+t53-t109-t110-t111-t112-t113-t114-t115-t116;
	t117 = t23*t36*t54;
	t190 = t16*t36*t45;
	t55 = -t117+t190;
	t56 = 1.0/pow(t34,3.0/2.0);
	t57 = t5*t6*4.0;
	t58 = e*r1*1.6E1;
	t59 = f*t4*4.0;
	t60 = d*r2*t4*2.0;
	t61 = e*r1*t6*4.0;
	t62 = r2*t8*t9*8.0;
	t63 = d*r1*r3*8.0;
	t64 = r2*t6*t8*t9*2.0;
	t65 = r2*t7*t8*t9*2.0;
	t66 = r1*r2*r3*t5*4.0;
	t67 = e*1.6E1;
	t68 = e*t4*4.0;
	t69 = e*t6*4.0;
	t70 = t7*t8*t10*4.0;
	t71 = e*r1*r2*r3*4.0;
	t72 = d*r1*8.0;
	t73 = d*r1*t4*2.0;
	t74 = d*r1*t7*2.0;
	t75 = e*r2*t4*4.0;
	t76 = e*r2*t7*4.0;
	t77 = f*r3*t4*2.0;
	t78 = r3*t5*t6*2.0;
	t79 = f*r1*r2*8.0;
	t80 = r1*t6*t8*t9*2.0;
	t81 = r2*r3*t8*t9*8.0;
	t91 = d*r1*t6*2.0;
	t93 = r1*t8*t9*8.0;
	t94 = d*r2*r3*8.0;
	t96 = r1*t4*t8*t9*2.0;
	t97 = r1*t7*t8*t9*2.0;
	t146 = r3*t5*t7*2.0;
	t147 = r3*t5*8.0;
	t82 = t24+t25-t26-t27-t28+t72+t73+t74+t75+t76+t77+t78+t79+t80+t81-t91-t93-t94-t96-t97-t146-t147-r2*t4*t8*t10*4.0-r2*t7*t8*t10*4.0;
	t83 = 1.0/c;
	t84 = d*1.6E1;
	t85 = d*t4*4.0;
	t86 = d*t6*4.0;
	t87 = t7*t8*t9*4.0;
	t88 = r1*r2*r3*t8*t9*4.0;
	t89 = e*r2*8.0;
	t90 = e*r2*t6*2.0;
	t92 = e*r2*t7*2.0;
	t95 = e*r1*r3*8.0;
	t98 = r2*t4*t8*t10*2.0;
	t99 = f*r2*8.0;
	t100 = d*t7*4.0;
	t101 = f*r2*t6*2.0;
	t102 = f*r2*t4*2.0;
	t103 = r2*t5*t7*2.0;
	t104 = r3*t8*t10*1.6E1;
	t105 = f*r1*r3*8.0;
	t106 = t6*t8*t9*4.0;
	t107 = d*r1*r2*r3*4.0;
	t108 = r3*t6*t8*t10*4.0;
	t118 = t5*t10*t14*2.0;
	t119 = t5*t9*t13*2.0;
	t199 = t8*t15*2.0;
	t120 = t118+t119-t199;
	t121 = f*4.0;
	t122 = t5*4.0;
	t123 = t4*t5;
	t124 = t5*t6;
	t125 = e*r1*4.0;
	t126 = f*t7;
	t127 = r2*t8*t9*4.0;
	t128 = d*r1*r3*2.0;
	t129 = e*r2*r3*2.0;
	t191 = t5*t7;
	t192 = d*r2*4.0;
	t193 = f*t4;
	t194 = f*t6;
	t195 = r1*t8*t10*4.0;
	t196 = r2*r3*t8*t10*2.0;
	t197 = r1*r3*t8*t9*2.0;
	t130 = t121-t122+t123+t124+t125+t126+t127+t128+t129-t191-t192-t193-t194-t195-t196-t197;
	t198 = t36*t83*t130;
	t131 = t117-t198;
	t132 = d*t10;
	t201 = e*t9;
	t133 = t132-t201;
	t134 = d*r3*8.0;
	t135 = d*r3*t7*2.0;
	t136 = e*t7*4.0;
	t137 = r1*t5*1.6E1;
	t138 = d*r3*t6*2.0;
	t139 = r1*t5*t7*4.0;
	t140 = d*r1*r2*8.0;
	t141 = t4*t8*t10*4.0;
	t142 = r3*t4*t8*t9*2.0;
	t143 = r1*r2*r3*t8*t10*4.0;
	t174 = t8*t10*1.6E1;
	t181 = t6*t8*t10*4.0;
	t144 = t67-t68+t69-t70-t71+t134+t135+t136+t137+t138+t139+t140+t141+t142+t143-t174-t181-f*r1*1.6E1-d*r3*t4*2.0-f*r1*t7*4.0-r3*t8*t9*8.0-r1*r2*t8*t9*8.0-r3*t6*t8*t9*2.0-r3*t7*t8*t9*2.0;
	t145 = t12*t83*t144;
	t148 = d*r1*t6*4.0;
	t149 = d*r1*t7*4.0;
	t150 = e*r2*t4*2.0;
	t151 = r1*r3*t8*t10*8.0;
	t152 = r2*t7*t8*t10*2.0;
	t153 = e*r1*8.0;
	t154 = e*r1*t4*2.0;
	t155 = e*r1*t7*2.0;
	t156 = r2*t8*t9*1.6E1;
	t157 = e*r2*r3*8.0;
	t158 = r2*t4*t8*t9*4.0;
	t159 = r1*t6*t8*t10*2.0;
	t160 = t17-t18-t19-t20+t21-t22+t57+t59+t66+t153+t154+t155+t156+t157+t158+t159-d*r2*1.6E1-t5*t7*4.0-d*r2*t4*4.0-e*r1*t6*2.0-r1*t8*t10*8.0-r2*r3*t8*t10*8.0-r1*t4*t8*t10*2.0-r1*t7*t8*t10*2.0;
	t161 = f*r2*1.6E1;
	t162 = e*r3*t6*2.0;
	t163 = f*r2*t7*4.0;
	t164 = r3*t8*t10*8.0;
	t165 = e*r1*r2*8.0;
	t166 = r3*t7*t8*t10*2.0;
	t167 = r3*t4*t8*t10*2.0;
	t168 = t84+t85-t86-t87-t88+t100+t106+t107+t161+t162+t163+t164+t165+t166+t167-e*r3*8.0-r2*t5*1.6E1-t8*t9*1.6E1-e*r3*t4*2.0-e*r3*t7*2.0-r2*t5*t7*4.0-t4*t8*t9*4.0-r1*r2*t8*t10*8.0-r3*t6*t8*t10*2.0;
	t169 = t12*t83*t168;
	t170 = f*r3*t4*4.0;
	t171 = f*r3*t6*4.0;
	t172 = t72+t73-t74-t80-t81+t89+t90+t91-t92-t93+t94-t95-t96+t97-t98+t150+t151+t152+t170+t171-r3*t4*t5*4.0-r3*t5*t6*4.0-r2*t8*t10*8.0-r2*t6*t8*t10*2.0;
	t173 = r1*t4*t5*2.0;
	t175 = d*r3*1.6E1;
	t176 = r1*t5*8.0;
	t177 = d*r3*t4*4.0;
	t178 = f*r1*t7*2.0;
	t179 = r1*t5*t6*2.0;
	t180 = f*r2*r3*8.0;
	t200 = t190-t198;
	A0(0,0) = t35*(t12*t16*(t24+t25+t26+t27+t28+t89+t90+t92+t95+t98+t148+t149-r3*t5*8.0-f*r1*r2*8.0-e*r2*t4*2.0-f*r3*t4*2.0-r3*t5*t6*2.0-r3*t5*t7*2.0-r2*t8*t10*8.0-r1*r3*t8*t10*8.0-r1*t6*t8*t9*4.0-r1*t7*t8*t9*4.0-r2*t6*t8*t10*2.0-r2*t7*t8*t10*2.0)+t12*t23*(t17-t18+t19+t20+t21+t22+t58+t60+t61+t62+t63+t64+t65-d*r2*8.0-f*t4*4.0-t5*t6*4.0-t5*t7*4.0-d*r2*t6*2.0-d*r2*t7*2.0-r1*t8*t10*1.6E1-r1*r2*r3*t5*4.0-r1*r3*t8*t9*8.0-r2*t4*t8*t9*2.0-r1*t6*t8*t10*4.0));
	A0(0,1) = -t35*(t12*t23*t82-t12*t16*t160);
	A0(0,2) = -t35*(t12*t23*(t84+t85+t86+t87+t88+t99+t101+t102+t103+t104+t105+t108-e*r3*1.6E1-d*t7*4.0-r2*t5*8.0-t8*t9*1.6E1-e*r3*t6*4.0-f*r2*t7*2.0-r1*r3*t5*8.0-r2*t4*t5*2.0-r2*t5*t6*2.0-t4*t8*t9*4.0-t6*t8*t9*4.0-d*r1*r2*r3*4.0)+t12*t16*(t67+t68+t69+t70+t71+t173+t175+t176+t177+t178+t179+t180-f*r1*8.0-e*t7*4.0-t8*t10*1.6E1-f*r1*t4*2.0-f*r1*t6*2.0-r2*r3*t5*8.0-r1*t5*t7*2.0-r3*t8*t9*1.6E1-t4*t8*t10*4.0-t6*t8*t10*4.0-r3*t4*t8*t9*4.0-r1*r2*r3*t8*t10*4.0));
	A0(0,3) = t55*t56*t120*(-1.0/2.0)-t16*t23*t35*t36*(a*r1*t8*4.0+b*r2*t8*4.0+a*t5*t10*4.0-b*t5*t9*4.0-a*r2*r3*t8*2.0+b*r1*r3*t8*2.0+a*r3*t5*t9*4.0+b*r3*t5*t10*4.0-a*t4*t5*t10+a*t5*t6*t10-a*t5*t7*t10-b*t4*t5*t9+b*t5*t6*t9+b*t5*t7*t9+a*r1*r2*t5*t9*2.0-b*r1*r2*t5*t10*2.0);
	A0(0,4) = -t8*t55*t56*t133-t8*t16*t23*t35*t36*(a*t9*-4.0-b*t10*4.0+a*r3*t10*4.0-b*r3*t9*4.0+a*t4*t9-a*t6*t9+a*t7*t9-b*t4*t10+b*t6*t10+b*t7*t10+a*r1*r2*t10*2.0+b*r1*r2*t9*2.0);
	A0(1,0) = -t35*(t145+t12*t23*(t17-t18+t19+t20+t21+t22-t57+t58-t59+t60+t61+t62+t63+t64+t65-t66-d*r2*8.0-t5*t7*4.0-d*r2*t6*2.0-d*r2*t7*2.0-r1*t8*t10*1.6E1-r1*r3*t8*t9*8.0-r2*t4*t8*t9*2.0-r1*t6*t8*t10*4.0));
	A0(1,1) = t35*(t169+t12*t23*t82);
	A0(1,2) = t35*(t12*t23*(t84+t85+t86+t87+t88+t99-t100+t101+t102+t103+t104+t105-t106-t107+t108-e*r3*1.6E1-r2*t5*8.0-t8*t9*1.6E1-e*r3*t6*4.0-f*r2*t7*2.0-r1*r3*t5*8.0-r2*t4*t5*2.0-r2*t5*t6*2.0-t4*t8*t9*4.0)-t12*t83*t172);
	A0(1,3) = t56*t120*t131*(-1.0/2.0)+t23*t35*t36*t83*(b*t8*4.0+c*r1*t8*4.0-b*t4*t8-b*t6*t8+b*t7*t8+c*t5*t10*4.0-c*r2*r3*t8*2.0-b*r1*t5*t10*4.0+b*r2*t5*t9*4.0+c*r3*t5*t9*4.0-c*t4*t5*t10+c*t5*t6*t10-c*t5*t7*t10-b*r1*r3*t5*t9*2.0-b*r2*r3*t5*t10*2.0+c*r1*r2*t5*t9*2.0);
	A0(1,4) = -t8*t56*t131*t133+t8*t23*t35*t36*t83*(c*t9*-4.0+b*r1*t9*4.0+b*r2*t10*4.0+c*r3*t10*4.0+c*t4*t9-c*t6*t9+c*t7*t9-b*r1*r3*t10*2.0+b*r2*r3*t9*2.0+c*r1*r2*t10*2.0);
	A0(2,0) = t35*(t145-t12*t16*(t24+t25+t26+t27+t28-t77-t78-t79+t89+t90+t92+t95+t98-t146-t147+t148+t149-t150-t151-t152-r2*t8*t10*8.0-r1*t6*t8*t9*4.0-r1*t7*t8*t9*4.0-r2*t6*t8*t10*2.0));
	A0(2,1) = -t35*(t169+t12*t16*t160);
	A0(2,2) = t35*(t12*t16*(t67+t68+t69+t70+t71-t136-t141-t143+t173-t174+t175+t176+t177+t178+t179+t180-t181-f*r1*8.0-f*r1*t4*2.0-f*r1*t6*2.0-r2*r3*t5*8.0-r1*t5*t7*2.0-r3*t8*t9*1.6E1-r3*t4*t8*t9*4.0)+t12*t83*t172);
	A0(2,3) = t56*t120*t200*(1.0/2.0)+t16*t35*t36*t83*(a*t8*-4.0+c*r2*t8*4.0+a*t4*t8+a*t6*t8-a*t7*t8-c*t5*t9*4.0+c*r1*r3*t8*2.0+a*r1*t5*t10*4.0-a*r2*t5*t9*4.0+c*r3*t5*t10*4.0-c*t4*t5*t9+c*t5*t6*t9+c*t5*t7*t9+a*r1*r3*t5*t9*2.0+a*r2*r3*t5*t10*2.0-c*r1*r2*t5*t10*2.0);
	A0(2,4) = t8*t56*t200*(t132-t201)-t8*t16*t35*t36*t83*(c*t10*4.0+a*r1*t9*4.0+a*r2*t10*4.0+c*r3*t9*4.0+c*t4*t10-c*t6*t10-c*t7*t10-a*r1*r3*t10*2.0+a*r2*r3*t9*2.0-c*r1*r2*t9*2.0);



	return A0;
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
	//t3=variables(0,5) ;
	MatrixXd A0(3,6);
	
	return A0;
}


MatrixXd jacobianForRotationAndTransition__(const MatrixXd& parameters,const MatrixXd& variables)
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
	double t5 ,t6 ,t7 ,t8 ,t9 ,t10 ,t11 ,t12 ,t13 ,t14 ,t15 ,t16 ,t17 ,t18 ,t19 ,t20 ,t21 ,t22 ,t23 ,t24 ,t25 ,t26 ,t27 ,t28 ,t31 ,t29 ,t30 ,t32 ,t44 ,t33 ,t34 ,t35 ,t36 ,t37 ,t38 ,t50 ,t39 ,t40 ,t41 ,t49 ,t42 ,t43 ,t45 ,t51 ,t63 ,t46 ,t47 ,t48 ,t52 ,t53 ,t54 ,t55 ,t56 ,t57 ,t58 ,t59 ,t60 ,t61 ,t62 ,t64 ,t65 ,t66 ,t67 ,t68 ,t69 ,t78 ,t70 ,t71 ,t72 ,t73 ,t74 ,t79 ,t75 ,t83 ,t76 ,t77 ,t80 ,t81 ,t82 ,t84 ,t85 ,t86 ,t87 ,t88 ,t89 ,t90 ,t91 ,t92 ,t93 ,t94 ,t95 ,t96 ,t97 ,t98 ,t99 ,t100 ,t101 ,t102 ,t103 ,t104 ,t105 ,t106 ,t107 ,t108 ,t109 ;
	t5 = d-t1;
	t6 = e-t2;
	t7 = f-t3;
	t8 = r1*r1;
	t9 = r2*r2;
	t10 = r3*r3;
	t11 = 1.0/a;
	t12 = 1.0/b;
	t13 = t5*t5;
	t14 = t6*t6;
	t15 = t7*t7;
	t16 = t13+t14+t15;
	t17 = 1.0/sqrt(t16);
	t18 = t8+t9+t10+4.0;
	t19 = 1.0/(t18*t18);
	t20 = t8*(1.0/4.0);
	t21 = t9*(1.0/4.0);
	t22 = t10*(1.0/4.0);
	t23 = t20+t21+t22+1.0;
	t24 = 1.0/t23;
	t25 = r3*t24;
	t26 = r1*r2*t24*(1.0/2.0);
	t27 = t8*t24*(1.0/2.0);
	t28 = t20+t21+t22-1.0;
	t31 = t24*t28;
	t29 = t27-t31;
	t30 = t25+t26;
	t32 = t25-t26;
	t44 = t9*t24*(1.0/2.0);
	t33 = t31-t44;
	t34 = r2*t24;
	t35 = r1*r3*t24*(1.0/2.0);
	t36 = t34+t35;
	t37 = t7*t36;
	t38 = t5*t29;
	t50 = t6*t32;
	t39 = t37+t38-t50;
	t40 = t11*t39;
	t41 = r1*t24;
	t49 = r2*r3*t24*(1.0/2.0);
	t42 = t41-t49;
	t43 = t7*t42;
	t45 = t6*t33;
	t51 = t5*t30;
	t63 = t43+t45-t51;
	t46 = t12*t63;
	t47 = t40+t46;
	t48 = 1.0/pow(t16,3.0/2.0);
	t52 = b*e*1.6E1;
	t53 = b*t2*t10*4.0;
	t54 = b*e*t8*4.0;
	t55 = b*e*t9*4.0;
	t56 = b*e*r1*r2*r3*4.0;
	t57 = 1.0/c;
	t58 = b*e*r2*8.0;
	t59 = b*e*r2*t9*2.0;
	t60 = b*e*r1*r3*8.0;
	t61 = b*e*r2*t10*2.0;
	t62 = b*r2*t2*t8*2.0;
	t64 = t34-t35;
	t65 = d*2.0;
	t66 = t1*2.0;
	t67 = t65-t66;
	t68 = t41+t49;
	t69 = t5*t64;
	t78 = t10*t24*(1.0/2.0);
	t70 = t31-t78;
	t71 = t7*t70;
	t72 = e*2.0;
	t73 = t2*2.0;
	t74 = t72-t73;
	t79 = t6*t68;
	t75 = t69+t71-t79;
	t83 = t57*t75;
	t76 = t46-t83;
	t77 = t12*t42;
	t80 = f*2.0;
	t81 = t3*2.0;
	t82 = t80-t81;
	t84 = c*f*r3*8.0;
	t85 = c*f*r3*t10*2.0;
	t86 = c*f*r1*r2*8.0;
	t87 = c*f*r3*t8*2.0;
	t88 = c*r3*t3*t9*2.0;
	t89 = a*d*1.6E1;
	t90 = c*f*1.6E1;
	t91 = a*t1*t10*4.0;
	t92 = c*t3*t8*4.0;
	t93 = a*d*t8*4.0;
	t94 = a*d*t9*4.0;
	t95 = c*f*t9*4.0;
	t96 = c*f*t10*4.0;
	t97 = c*f*r1*r2*r3*4.0;
	t98 = a*r1*r2*r3*t1*4.0;
	t99 = a*r1*t1*t8*2.0;
	t100 = a*r1*t1*8.0;
	t101 = a*d*r2*r3*8.0;
	t102 = a*d*r1*t9*2.0;
	t103 = a*r1*t1*t10*2.0;
	t104 = t57*t64;
	t105 = t11*t29;
	t106 = t11*t32;
	t107 = t57*t68;
	t108 = t40+t83;
	t109 = t11*t36;
	A0(0,0) = t11*t12*t17*t19*(t58+t59+t60+t61+t62+a*f*1.6E1-a*t3*1.6E1-a*d*r2*8.0+a*e*r1*1.6E1+b*f*r3*8.0-a*f*t8*4.0+a*f*t9*4.0+a*f*t10*4.0-a*r1*t2*1.6E1+a*r2*t1*8.0-b*r2*t2*8.0-b*r3*t3*8.0+a*t3*t8*4.0-a*t3*t9*4.0-a*t3*t10*4.0+a*d*r1*r3*8.0+a*d*r2*t8*2.0-a*d*r2*t9*2.0-a*d*r2*t10*2.0-b*f*r1*r2*8.0+a*e*r1*t9*4.0+b*d*r1*t9*4.0+b*d*r1*t10*4.0-b*e*r2*t8*2.0-b*f*r3*t8*2.0+b*f*r3*t9*2.0+b*f*r3*t10*2.0-a*r1*r3*t1*8.0+b*r1*r2*t3*8.0-b*r1*r3*t2*8.0-a*r2*t1*t8*2.0-a*r1*t2*t9*4.0+a*r2*t1*t9*2.0+a*r2*t1*t10*2.0-b*r1*t1*t9*4.0-b*r1*t1*t10*4.0-b*r2*t2*t9*2.0-b*r2*t2*t10*2.0+b*r3*t3*t8*2.0-b*r3*t3*t9*2.0-b*r3*t3*t10*2.0+a*f*r1*r2*r3*4.0-a*r1*r2*r3*t3*4.0);
	A0(0,1) = t11*t12*t17*t19*(t99+t100+t101+t102+t103+b*f*1.6E1-b*t3*1.6E1-a*d*r1*8.0-b*d*r2*1.6E1+b*e*r1*8.0-a*f*r3*8.0+b*f*t8*4.0-b*f*t9*4.0+b*f*t10*4.0+a*r3*t3*8.0-b*r1*t2*8.0+b*r2*t1*1.6E1-b*t3*t8*4.0+b*t3*t9*4.0-b*t3*t10*4.0-a*f*r1*r2*8.0+b*e*r2*r3*8.0-a*d*r1*t8*2.0-a*d*r1*t10*2.0-a*e*r2*t8*4.0-b*d*r2*t8*4.0-a*e*r2*t10*4.0+b*e*r1*t8*2.0-b*e*r1*t9*2.0-a*f*r3*t8*2.0+b*e*r1*t10*2.0+a*f*r3*t9*2.0-a*f*r3*t10*2.0+a*r1*r2*t3*8.0-a*r2*r3*t1*8.0-b*r2*r3*t2*8.0-a*r1*t1*t9*2.0+a*r2*t2*t8*4.0+a*r2*t2*t10*4.0+a*r3*t3*t8*2.0-a*r3*t3*t9*2.0+a*r3*t3*t10*2.0-b*r1*t2*t8*2.0+b*r2*t1*t8*4.0+b*r1*t2*t9*2.0-b*r1*t2*t10*2.0-b*f*r1*r2*r3*4.0+b*r1*r2*r3*t3*4.0);
	A0(0,2) = -t11*t12*t17*t19*(t52+t53+t54+t55+t56+t89+t91+t93+t94+t98-a*t1*1.6E1-b*t2*1.6E1-a*e*r3*1.6E1+b*d*r3*1.6E1+a*f*r2*8.0-a*d*t10*4.0-b*f*r1*8.0-b*e*t10*4.0-a*r2*t3*8.0+a*r3*t2*1.6E1+b*r1*t3*8.0-b*r3*t1*1.6E1-a*t1*t8*4.0-a*t1*t9*4.0-b*t2*t8*4.0-b*t2*t9*4.0+a*f*r1*r3*8.0+b*f*r2*r3*8.0+b*d*r3*t8*4.0-a*e*r3*t9*4.0+a*f*r2*t8*2.0+a*f*r2*t9*2.0-a*f*r2*t10*2.0-b*f*r1*t8*2.0-b*f*r1*t9*2.0+b*f*r1*t10*2.0-a*r1*r3*t3*8.0-b*r2*r3*t3*8.0-a*r2*t3*t8*2.0-a*r2*t3*t9*2.0+a*r3*t2*t9*4.0+a*r2*t3*t10*2.0+b*r1*t3*t8*2.0-b*r3*t1*t8*4.0+b*r1*t3*t9*2.0-b*r1*t3*t10*2.0-a*d*r1*r2*r3*4.0-b*r1*r2*r3*t2*4.0);
	A0(0,3) = -t17*(t105-t12*t30)+t47*t48*t67*(1.0/2.0);
	A0(0,4) = t17*(t106-t12*t33)+t47*t48*t74*(1.0/2.0);
	A0(0,5) = -t17*(t77+t109)+t47*t48*t82*(1.0/2.0);
	A0(1,0) = -t12*t17*t19*t57*(t52-t53-t54+t55-t56+t90+t92+t95+t96+t97-b*t2*1.6E1-c*t3*1.6E1+b*d*r3*8.0-c*d*r2*8.0-b*f*r1*1.6E1+c*e*r1*1.6E1+b*e*t10*4.0-c*f*t8*4.0+b*r1*t3*1.6E1-b*r3*t1*8.0-c*r1*t2*1.6E1+c*r2*t1*8.0+b*t2*t8*4.0-b*t2*t9*4.0-c*t3*t9*4.0-c*t3*t10*4.0+b*d*r1*r2*8.0+c*d*r1*r3*8.0-b*d*r3*t8*2.0+b*d*r3*t9*2.0+b*d*r3*t10*2.0+c*d*r2*t8*2.0-c*d*r2*t9*2.0-c*d*r2*t10*2.0+c*e*r1*t9*4.0-b*f*r1*t10*4.0-b*r1*r2*t1*8.0-c*r1*r3*t1*8.0+b*r3*t1*t8*2.0-b*r3*t1*t9*2.0+b*r1*t3*t10*4.0-b*r3*t1*t10*2.0-c*r2*t1*t8*2.0-c*r1*t2*t9*4.0+c*r2*t1*t9*2.0+c*r2*t1*t10*2.0+b*r1*r2*r3*t2*4.0-c*r1*r2*r3*t3*4.0);
	A0(1,1) = t12*t17*t19*t57*(t84+t85+t86+t87+t88+b*d*1.6E1-b*t1*1.6E1+c*d*r1*8.0-b*e*r3*8.0+b*f*r2*1.6E1+b*d*t8*4.0-b*d*t9*4.0+b*d*t10*4.0-b*r2*t3*1.6E1+b*r3*t2*8.0-c*r1*t1*8.0-c*r3*t3*8.0-b*t1*t8*4.0+b*t1*t9*4.0-b*t1*t10*4.0+b*e*r1*r2*8.0-c*d*r2*r3*8.0+c*d*r1*t8*2.0-c*d*r1*t9*2.0-b*e*r3*t8*2.0+c*d*r1*t10*2.0+b*e*r3*t9*2.0-b*e*r3*t10*2.0+c*e*r2*t8*4.0+b*f*r2*t10*4.0+c*e*r2*t10*4.0-c*f*r3*t9*2.0-b*r1*r2*t2*8.0-c*r1*r2*t3*8.0+c*r2*r3*t1*8.0+b*r3*t2*t8*2.0-b*r3*t2*t9*2.0-b*r2*t3*t10*4.0+b*r3*t2*t10*2.0-c*r1*t1*t8*2.0+c*r1*t1*t9*2.0-c*r1*t1*t10*2.0-c*r2*t2*t8*4.0-c*r2*t2*t10*4.0-c*r3*t3*t8*2.0-c*r3*t3*t10*2.0+b*d*r1*r2*r3*4.0-b*r1*r2*r3*t1*4.0);
	A0(1,2) = -t12*t17*t19*t57*(t58+t59-t60-t61-t62-c*d*1.6E1+c*t1*1.6E1+b*d*r1*8.0+c*e*r3*1.6E1-c*f*r2*8.0-c*d*t8*4.0-c*d*t9*4.0+c*d*t10*4.0-b*r1*t1*8.0-b*r2*t2*8.0+c*r2*t3*8.0-c*r3*t2*1.6E1+c*t1*t8*4.0+c*t1*t9*4.0-c*t1*t10*4.0+b*d*r2*r3*8.0+b*d*r1*t8*2.0+b*d*r1*t9*2.0-b*d*r1*t10*2.0-c*f*r1*r3*8.0+b*e*r2*t8*2.0+b*f*r3*t8*4.0+b*f*r3*t9*4.0+c*e*r3*t9*4.0-c*f*r2*t8*2.0-c*f*r2*t9*2.0+c*f*r2*t10*2.0+b*r1*r3*t2*8.0-b*r2*r3*t1*8.0+c*r1*r3*t3*8.0-b*r1*t1*t8*2.0-b*r1*t1*t9*2.0+b*r1*t1*t10*2.0-b*r2*t2*t9*2.0+b*r2*t2*t10*2.0-b*r3*t3*t8*4.0-b*r3*t3*t9*4.0+c*r2*t3*t8*2.0+c*r2*t3*t9*2.0-c*r3*t2*t9*4.0-c*r2*t3*t10*2.0+c*d*r1*r2*r3*4.0-c*r1*r2*r3*t1*4.0);
	A0(1,3) = -t17*(t104+t12*t30)-t48*t67*t76*(1.0/2.0);
	A0(1,4) = t17*(t107+t12*(t31-t44))-t48*t74*t76*(1.0/2.0);
	A0(1,5) = t17*(t77-t57*t70)-t48*t76*t82*(1.0/2.0);
	A0(2,0) = -t11*t17*t19*t57*(t84+t85-t86-t87-t88-a*e*1.6E1+a*t2*1.6E1-a*d*r3*8.0+a*f*r1*1.6E1+c*e*r2*8.0+a*e*t8*4.0-a*e*t9*4.0-a*e*t10*4.0-a*r1*t3*1.6E1+a*r3*t1*8.0-c*r2*t2*8.0-c*r3*t3*8.0-a*t2*t8*4.0+a*t2*t9*4.0+a*t2*t10*4.0-a*d*r1*r2*8.0+a*d*r3*t8*2.0-a*d*r3*t9*2.0-a*d*r3*t10*2.0+c*e*r1*r3*8.0+c*d*r1*t9*4.0+a*f*r1*t10*4.0+c*d*r1*t10*4.0-c*e*r2*t8*2.0+c*e*r2*t9*2.0+c*e*r2*t10*2.0+c*f*r3*t9*2.0+a*r1*r2*t1*8.0+c*r1*r2*t3*8.0-c*r1*r3*t2*8.0-a*r3*t1*t8*2.0+a*r3*t1*t9*2.0-a*r1*t3*t10*4.0+a*r3*t1*t10*2.0-c*r1*t1*t9*4.0-c*r1*t1*t10*4.0+c*r2*t2*t8*2.0-c*r2*t2*t9*2.0-c*r2*t2*t10*2.0+c*r3*t3*t8*2.0-c*r3*t3*t10*2.0+a*e*r1*r2*r3*4.0-a*r1*r2*r3*t2*4.0);
	A0(2,1) = -t11*t17*t19*t57*(t89+t90-t91-t92+t93-t94-t95+t96-t97-t98-a*t1*1.6E1-c*t3*1.6E1-a*e*r3*8.0+a*f*r2*1.6E1-c*d*r2*1.6E1+a*d*t10*4.0+c*e*r1*8.0+c*f*t8*4.0-a*r2*t3*1.6E1+a*r3*t2*8.0-c*r1*t2*8.0+c*r2*t1*1.6E1-a*t1*t8*4.0+a*t1*t9*4.0+c*t3*t9*4.0-c*t3*t10*4.0+a*e*r1*r2*8.0+c*e*r2*r3*8.0-a*e*r3*t8*2.0+a*e*r3*t9*2.0-a*e*r3*t10*2.0-c*d*r2*t8*4.0+a*f*r2*t10*4.0+c*e*r1*t8*2.0-c*e*r1*t9*2.0+c*e*r1*t10*2.0-a*r1*r2*t2*8.0-c*r2*r3*t2*8.0+a*r3*t2*t8*2.0-a*r3*t2*t9*2.0-a*r2*t3*t10*4.0+a*r3*t2*t10*2.0-c*r1*t2*t8*2.0+c*r2*t1*t8*4.0+c*r1*t2*t9*2.0-c*r1*t2*t10*2.0+a*d*r1*r2*r3*4.0+c*r1*r2*r3*t3*4.0);
	A0(2,2) = -t11*t17*t19*t57*(t99+t100-t101-t102-t103-c*e*1.6E1+c*t2*1.6E1-a*d*r1*8.0-a*e*r2*8.0-c*d*r3*1.6E1+c*f*r1*8.0-c*e*t8*4.0-c*e*t9*4.0+c*e*t10*4.0+a*r2*t2*8.0-c*r1*t3*8.0+c*r3*t1*1.6E1+c*t2*t8*4.0+c*t2*t9*4.0-c*t2*t10*4.0+a*e*r1*r3*8.0-a*d*r1*t8*2.0+a*d*r1*t10*2.0-a*e*r2*t8*2.0-a*e*r2*t9*2.0+a*e*r2*t10*2.0-c*f*r2*r3*8.0-a*f*r3*t8*4.0-c*d*r3*t8*4.0-a*f*r3*t9*4.0+c*f*r1*t8*2.0+c*f*r1*t9*2.0-c*f*r1*t10*2.0-a*r1*r3*t2*8.0+a*r2*r3*t1*8.0+c*r2*r3*t3*8.0+a*r1*t1*t9*2.0+a*r2*t2*t8*2.0+a*r2*t2*t9*2.0-a*r2*t2*t10*2.0+a*r3*t3*t8*4.0+a*r3*t3*t9*4.0-c*r1*t3*t8*2.0+c*r3*t1*t8*4.0-c*r1*t3*t9*2.0+c*r1*t3*t10*2.0-c*e*r1*r2*r3*4.0+c*r1*r2*r3*t2*4.0);
	A0(2,3) = t17*(t104+t105)-t48*t67*t108*(1.0/2.0);
	A0(2,4) = -t17*(t106+t107)-t48*t74*t108*(1.0/2.0);
	A0(2,5) = t17*(t109+t57*(t31-t78))-t48*t82*t108*(1.0/2.0);



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
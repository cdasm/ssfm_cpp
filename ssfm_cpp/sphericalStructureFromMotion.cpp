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

MatrixXd convertTransitionRotationToCameraPosition(const MatrixXd& t,const MatrixXd& r)
{
	MatrixXd result(1,6);
	result.block(0,0,1,3)=r;
	result.block(0,3,1,3)=t;
	return result;
}

vector<MatrixXd> transitionAndRotationFromEssential(const MatrixXd& essential)
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

	vector<MatrixXd > result;

	MatrixXd Ty=transitionFromCrossMatrix(Tx);

	R1=rotationThomasonPara(R1.inverse());
	R2=rotationThomasonPara(R2.inverse());

	

	result.push_back(convertTransitionRotationToCameraPosition(Ty,R1));
	result.push_back(convertTransitionRotationToCameraPosition(Ty*-1,R1));
	result.push_back(convertTransitionRotationToCameraPosition(Ty,R2));
	result.push_back(convertTransitionRotationToCameraPosition(Ty*-1,R2));

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





MatrixXd projectionErrorForPointReconstruction(const MatrixXd& var1,const MatrixXd& var2,const MatrixXd& var3)
{
	return  projectionError (var1,var3,var2);
}


MatrixXd jacobianForPointReconstruction(const MatrixXd& var1,const MatrixXd& var2,const MatrixXd& var3)
{
	return jacobianForPoint (var1,var3,var2);
}

MatrixXd coeffcientLinear(const MatrixXd& var0,const MatrixXd& var1)
{
	double b0,b1,b2;
	b0=var0(0,0);
	b1=var0(0,1);
	b2=var0(0,2);
	double c0,c1,c2;
	c0=var1(0,0);
	c1=var1(0,1);
	c2=var1(0,2);
	MatrixXd A0(3,4);
	double t2 ,t3 ,t4 ,t5 ,t6 ,t7 ,t8 ,t9 ,t10 ,t11 ;
	t2 = c1*c1;
	t3 = c2*c2;
	t4 = t3*2.0;
	t5 = c0*c0;
	t6 = t5*2.0;
	t7 = t2*2.0;
	t8 = b0*c0;
	t9 = b1*c1;
	t10 = b2*c2;
	t11 = t8+t9+t10;
	A0(0,0) = t4+t7;
	A0(0,1) = c0*c1*-2.0;
	A0(0,2) = c0*c2*-2.0;
	A0(0,3) = b0*(t2+t3)*-2.0+c0*(b1*c1*2.0+b2*c2*2.0);
	A0(1,0) = c0*c1*-2.0;
	A0(1,1) = t4+t6;
	A0(1,2) = c1*c2*-2.0;
	A0(1,3) = b1*-2.0+c1*t11*2.0;
	A0(2,0) = c0*c2*-2.0;
	A0(2,1) = c1*c2*-2.0;
	A0(2,2) = t6+t7;
	A0(2,3) = b2*-2.0+c2*t11*2.0;
	return A0;
}

MatrixXd reconstructPointLinear(const MatrixXd& projections, const MatrixXd& cameras)
{
	assert(projections.rows()>1);
	assert(projections.rows()==cameras.rows());

	auto calculateSingleCoefficientLinear=[](const MatrixXd& prj,const MatrixXd& cam)->MatrixXd
	{
		MatrixXd p=cam.block(0,3,1,3);
		MatrixXd u=(rotationThomason(cam.block(0,0,1,3)).transpose()*prj.transpose()).transpose();
		return coeffcientLinear(p,u);
	};

	MatrixXd coefficients(3,4);
	coefficients*=0.0;

	for (int i = 0; i < projections.rows(); i++)
	{
		coefficients+=calculateSingleCoefficientLinear(projections.row(i),cameras.row(i));
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
/*
MatrixXd reconstructPoint(const MatrixXd& projections,const MatrixXd& cameras)
{	
	MatrixXd obj_vals=MatrixXd::Zero(1,projections.rows()*observationDim);
	MatrixXd initPara=MatrixXd::Random(1,pointDim);
	vector<MatrixXd> data(2);
	data[0]=projections;
	data[1]=cameras;
	return levenbergM_adhoc(data,obj_vals,&projectionErrorForPointReconstruction,&jacobianForPointReconstruction,initPara);
}*/

double length(const MatrixXd& a)
{
	return sqrt((a*a.transpose())(0,0));
}

double projectionErrorValue(const MatrixXd& projection,const MatrixXd& point,const MatrixXd& camera)
{
	auto terror=projectionError(projection,point,camera);

	return sqrt((terror*terror.transpose())(0,0));
}



bool pointBeforeCamera_xpu(const MatrixXd&x,const MatrixXd&p,const MatrixXd& u)
{
	MatrixXd tx=x-p;

	auto core=[](double a,double b)->bool
	{
		return (a*b)>0 && abs(a)>abs(b);
	};
	return core(tx(0,0),u(0,0)) && core(tx(0,1),u(0,1)) && core(tx(0,2),u(0,2));
}

bool pointBeforeCamera(const MatrixXd& point,const MatrixXd& projection,const MatrixXd& cameraPos)
{
	MatrixXd transition=cameraPos.block(0,3,1,3);
	MatrixXd rotation=rotationThomason(cameraPos.block(0,0,1,3));

	MatrixXd coordinate=(rotation*(point-transition).transpose()).transpose();

	return pointBeforeCamera_xpu(coordinate,MatrixXd::Zero(1,3),projection) && length(coordinate)>1.0;
}

bool reconstructPoint(const MatrixXd& projections,const MatrixXd& cameras,vector<bool>& flags,MatrixXd& pnt)
{
	assert(projections.rows()>1);
	assert(projections.rows()==cameras.rows());
	assert(projections.rows()==flags.size());
	flags.resize(flags.size(),false);




	int goodcamera=0;
	int usedcamera;
	MatrixXd tprojections=projections;
	MatrixXd tcameras=cameras;
	vector<int> match(flags.size());
	for(int i=0;i<match.size();++i) match[i]=i;

	do
	{
		pnt=reconstructPointLinear(tprojections,tcameras);

		usedcamera=tprojections.rows();
		goodcamera=0;	
		
		vector<int> tmatch;
		for (int i = 0; i < tprojections.rows(); i++)
		{
			if(pointBeforeCamera(pnt,tprojections.row(i),tcameras.row(i)) && projectionErrorValue(pnt,tprojections.row(i),tcameras.row(i))<constrain_on_point_error )
			{
				tmatch.push_back(match[i]);
				++goodcamera;
			}
		}

		match.swap(tmatch);

		tprojections=MatrixXd(goodcamera,projections.cols());
		tcameras=MatrixXd(goodcamera,cameras.cols());

		for (int i = 0; i < match.size(); i++)
		{
			tprojections.row(i)=projections.row(match[i]);
			tcameras.row(i)=cameras.row(match[i]);
		}


	}while(goodcamera<usedcamera && goodcamera>1);

	for (int i = 0; i < match.size(); i++)
	{
		flags[match[i]]=true;
	}

	if(goodcamera<2)
		return false;
	
	return true;
	


}

/*



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

	return good;

}*/



double angleBetween(const MatrixXd& a,const MatrixXd& b)
{
	return acos((a*b.transpose())(0,0)/length(a)/length(b));
}






pair<MatrixXd,vector<double> > reconstructPointsFor2Cameras(const MatrixXd& spnts1,const vector<int>& ind1,const MatrixXd& spnts2,const vector<int>& ind2,const MatrixXd&  cameraPositions)
{
	assert(ind1.size()==ind2.size());

	assert(cameraPositions.rows()==2);

	MatrixXd points(ind1.size(),3);
	vector<double> error(ind1.size(),-1);


	for (int i = 0; i < ind1.size(); i++)
	{
		
		cout<<i<<"\t";
		MatrixXd projs(2,3);
		projs.row(0)=spnts1.row(ind1[i]);
		projs.row(1)=spnts2.row(ind2[i]);
		points.row(i)=reconstructPointLinear(projs,cameraPositions);
		if(pointBeforeCamera(points.row(i),projs.row(0),cameraPositions.row(0) ) && pointBeforeCamera(points.row(i),projs.row(1),cameraPositions.row(1) ))
		{
			
			error[i]=0;
			for(int pui=0;pui<2;++pui)
			{
				MatrixXd terror=projectionError(projs.row(pui),points.row(i),cameraPositions.row(pui));
				error[i]+= (terror*terror.transpose())(0,0);
			}
			error[i]/=2.0;
		}
	}

	return make_pair(points,error);
}




auto geometricReconstructionFrom2Frames(const MatrixXd& sphericalPoints1,const vector<int>& ind1,const MatrixXd& sphericalPoints2,const vector<int>& ind2,MatrixXd& cameraPosition)->pair<MatrixXd,vector<double> >
{

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

	vector<MatrixXd > transtionAndRotations;
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
		MatrixXd cameras(2,transtionAndRotations[0].cols());
		cameras.row(0)*=0;
		cameras.row(1)=transtionAndRotations[i];

		bestPointCandidates[i]=reconstructPointsFor2Cameras (sphericalPoints1,ind1,sphericalPoints2,ind2,cameras);
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
	cameraPosition=transtionAndRotations[bestIndex];


	return bestPointCandidates[bestIndex];

}


auto threeDimensionReconstruction(const string& featureFileName,const string& matchFileName,int imageWidth,int imageHeight)->pair<MatrixXd,MatrixXd>
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


	

	MatrixXd points(trajectories.size(),3);

	vector<bool> alreadyReconstructed(trajectories.size(),false);

	vector<int> imageSize(2);
	imageSize[0]=imageWidth;imageSize[1]=imageHeight;

	auto convertFeature=[&](const vector<vector<int>> &fea)->MatrixXd
	{
		MatrixXd sfea(fea.size(),3);
		for (int i = 0; i < fea.size(); i++)
		{
			sfea.row(i)=imageCordinate2Phere(fea[i],imageSize);
		}
		return sfea;
	};
	vector<MatrixXd> sphericalFeaturesAndCameras(features.size()+1);
	
	for (int i = 0; i < sphericalFeaturesAndCameras.size()-1; i++)
	{
		sphericalFeaturesAndCameras[i]=convertFeature(features[i]);
	}


	MatrixXd& cameraPosition=sphericalFeaturesAndCameras[sphericalFeaturesAndCameras.size()-1];

	cameraPosition=MatrixXd(features.size(),6);

	vector<bool> alreadyEstimated(features.size(),false);

	cameraPosition.row(0)*=0;

	alreadyEstimated[0]=true;

	

	
	auto buildFromNothingForTwoFrames=[&correspondences,&sphericalFeaturesAndCameras,&featureIsPoint,&cameraPosition,&alreadyEstimated,&alreadyReconstructed,&points](int cameraIndex)
	{
		assert(cameraIndex==0);
		assert(cameraIndex<correspondences.size());
		vector<int> index(correspondences[cameraIndex].size());

		vector<int> ind1,ind2;

		int count=0;
		for ( auto&k:correspondences[cameraIndex])
		{
			ind1.push_back(k.first);
			ind2.push_back(k.second);
			index[count]=featureIsPoint[cameraIndex][k.first];
			++count;
		}

		MatrixXd curCamerapos;
		auto pointsError2Frame=geometricReconstructionFrom2Frames(sphericalFeaturesAndCameras[0],ind1,sphericalFeaturesAndCameras[1],ind2,curCamerapos)	;


		cameraPosition.row(cameraIndex+1)=curCamerapos;

		alreadyEstimated[cameraIndex+1]=true;

		
		for (int i = 0; i < pointsError2Frame.second.size(); i++)
		{
			if (pointsError2Frame.second[i]>0 && pointsError2Frame.second[i]<constrain_on_point_error  )
			{
				points.row(index[i])=pointsError2Frame.first.row(i);
				alreadyReconstructed[index[i]]=true;

			}
			else
			{
				alreadyReconstructed[index[i]]=false;
			}

		}

	};
	
	buildFromNothingForTwoFrames(0);
	alreadyEstimated[1]=true;

	

	
	//to be rewrite
	auto reconstructPoints=[&](int lcameraIndex)
	{
		assert(lcameraIndex<features.size());

		for (auto& s:contain[lcameraIndex])
		{
			vector<int> curCameras;

			for (auto& w:trajectories[s])
			{
				if(alreadyEstimated[w.first])
				{
					curCameras.push_back(w.first);
				}
			}

			if(curCameras.size()>1)
			{
				MatrixXd curPrj(curCameras.size(),pointDim),curCams(curCameras.size(),cameraDim);
				vector<bool> curflags(curCameras.size());

				for (int i = 0; i < curCameras.size(); i++)
				{
					curPrj.row(i)=sphericalFeaturesAndCameras[curCameras[i] ].row( trajectories[s][curCameras[i]] );
					curCams.row(i)=cameraPosition.row(curCameras[i]);
				}
				MatrixXd pnt;
				if(reconstructPoint(curPrj,curCams,curflags,pnt))
				{
					alreadyReconstructed[s]=true;
					points.row(s)=pnt;
					for (int i = 0; i < curflags.size(); i++)
					{
						if(curflags[i])
						{
							projectionsValid[s][curCameras[i]]=true;
						}
					}
				}
			}
		}
	};

	vector<funcType2> functions(2);
	functions[0]=&projectionError;
	functions[1]=&projectionErrorUnitLength;
	vector<funcType2> jacabianFunctions(4);
	jacabianFunctions[0]=&jacobianForCamera;
	jacabianFunctions[1]=&jacobianForCameraUnitLength;
	jacabianFunctions[2]=&jacobianForPoint;
	jacabianFunctions[3]=&jacobianForPointUnitLength;

	auto insersectSpecial=[](unordered_map<int,bool>& a,unordered_map<int,cameraType>& b)->unordered_set<int>
	{
		unordered_set<int> c;
		for(auto& d:b)
		{
			if(a.count(d.first) && d.second)
				c.insert(d.first);
		}

		return c;
	};
	auto cameraPosUnitLenthConvert=[](const MatrixXd& v,int flag)
	{
		MatrixXd r;
		if(flag<0)
		{
			r=MatrixXd(1,5);
			r.block(0,3,1,2)=transition2Para(v.block(0,3,1,3));
		}
		else
		{
			r=MatrixXd(1,6);
			r.block(0,3,1,3)=transitionFrom2Para(v.block(0,3,1,2));
		
		}
		r.block(0,0,1,3)=v.block(0,0,1,3);
		return r;
	};


	auto bundleAdjustment=[&](unordered_map<int,cameraType>& lbundlePara)
	{
		unordered_set<int> pntIndx;

		unordered_map<MatrixXd*, unordered_map<int,pair<int,int> > > parameterLookUp;

		vector<pair<int,int> > projections;

		unordered_set<int> cameraForBundle;
		unordered_set<int> pointForBundle;

		for (auto& s:lbundlePara)
		{
			assert(alreadyEstimated[s.first]);
			if (s.second!=cameraType::_static)
			{
				for (auto&w: contain[s.first])
				{
					if( !pntIndx.count(w))
					{
						
						pntIndx.insert(w);
						unordered_set<int> curPrj;
						
						curPrj=insersectSpecial(projectionsValid[w],lbundlePara);

						if(curPrj.size()>1)
						{
							pointForBundle.insert(w);
							for(auto x:curPrj)
							{
								projections.push_back(make_pair(w,x));
								cameraForBundle.insert(x);
							}
						}
						
					}
				}
			}
		}

		int parameterNumber=0;
		for (auto& camera:cameraForBundle)
		{
			if(lbundlePara[camera]==cameraType::_unitLength)
			{
				parameterLookUp[&cameraPosition][camera]=make_pair(parameterNumber,cameraDimUnitLength);
				parameterNumber+=cameraDimUnitLength;
			}
			else if(lbundlePara[camera]==cameraType::_ordinary)
			{
				parameterLookUp[&cameraPosition][camera]=make_pair(parameterNumber,cameraDim);
				parameterNumber+=cameraDim;
			}
		}

		for(auto& pnt:pointForBundle)
		{
			parameterLookUp[&points][pnt]=make_pair(parameterNumber,pointDim);
			parameterNumber+=3;
		}

		/*
		vector<funcType2> functions(2);
		functions[0]=&functionForRotationAndTransition;
		functions[1]=&functionForRotationAndTransitionUnitLength;
		vector<funcType2> jacabianFunctions(4);
		jacabianFunctions[0]=&jacobianForRotationAndTransition;
		jacabianFunctions[1]=&jacobianForRotationAndTransitionUnitLength;
		jacabianFunctions[2]=&jacobianForPoint;
		jacabianFunctions[3]=&jacobianForPointUnitLength;
		*/
		unordered_map<cameraType,int> funcIndex;
		funcIndex[cameraType::_static]=0;
		funcIndex[cameraType::_ordinary]=0;		
		funcIndex[cameraType::_unitLength]=1;
		unordered_map<cameraType,int> jfuncForPointIndex;
		jfuncForPointIndex[cameraType::_ordinary]=2;		
		jfuncForPointIndex[cameraType::_unitLength]=3;
		unordered_map<cameraType,int> jfuncForCameraIndex;
		jfuncForCameraIndex[cameraType::_ordinary]=0;		
		jfuncForCameraIndex[cameraType::_unitLength]=1;


		vector<pair<pair<int,int>,vector<pair<dataORParameter,pair<int,int>>>>> funcDataMap;
		vector<pair<pair<int,pair<int,int>>,vector<pair<dataORParameter,pair<int,int>>>>> jfuncDataMap;

		for (int i = 0; i < projections.size(); i++)
		{
			auto& proj=projections[i];
			auto _camType=lbundlePara[proj.second];
			pair<int,int> funcPosmap;		
			funcPosmap.first=funcIndex[_camType];
			funcPosmap.second=i*observationDim;
			vector<pair<dataORParameter,pair<int,int>>> datamap(3);
			pair<int,pair<int,int>> jfuncPosmap;

			

			datamap[0]=make_pair(dataORParameter::_data,make_pair(proj.second,trajectories[proj.first][proj.second]));
			datamap[1]=make_pair(dataORParameter::_parameter,parameterLookUp[&points][proj.first]);

			jfuncPosmap.first=jfuncForPointIndex[_camType];
			jfuncPosmap.second=make_pair(i*observationDim,parameterLookUp[&points][proj.first].first);

			if (_camType==cameraType::_static)
			{								
					datamap[2]=make_pair(dataORParameter::_data,make_pair(sphericalFeaturesAndCameras.size()-1,proj.second));
					funcDataMap.push_back(make_pair(funcPosmap,datamap));
					
					jfuncDataMap.push_back(make_pair(jfuncPosmap,datamap));
			}
			else
			{							
					datamap[2]=make_pair(dataORParameter::_parameter,parameterLookUp[&cameraPosition][proj.second]);
					funcDataMap.push_back(make_pair(funcPosmap,datamap));

					jfuncDataMap.push_back(make_pair(jfuncPosmap,datamap));
					
					jfuncPosmap.first=jfuncForCameraIndex[_camType];
					jfuncPosmap.second=make_pair(i*observationDim,parameterLookUp[&cameraPosition][proj.second].first);
					jfuncDataMap.push_back(make_pair(jfuncPosmap,datamap));
																	
			}

		}



		MatrixXd obj_vals=MatrixXd::Zero(1,projections.size()*observationDim);
		MatrixXd intit_parameters(1,parameterNumber);
		

		for(auto& paratop:parameterLookUp)
		{
			for(auto& para:paratop.second)
			{
				if(paratop.first==(&cameraPosition) && lbundlePara[ para.first]==cameraType::_unitLength)
				{
					intit_parameters.block(0,para.second.first,1,para.second.second)=cameraPosUnitLenthConvert((* paratop.first).row(para.first),1);
				}
				else
				{
					intit_parameters.block(0,para.second.first,1,para.second.second)=(* paratop.first).row(para.first);
				}
			}
		}

		auto est_parameters=levenbergM_advanced(sphericalFeaturesAndCameras,funcDataMap,jfuncDataMap,obj_vals,functions,jacabianFunctions,intit_parameters,projections.size(),observationDim);

		for(auto& paratop:parameterLookUp)
		{
			for(auto& para:paratop.second)
			{
				if(paratop.first==(&cameraPosition) && lbundlePara[ para.first]==cameraType::_unitLength)
				{
					//intit_parameters.block(0,para.second.first,1,para.second.second)=cameraPosUnitLenthConvert((* para.first.first).row(para.first.second),1);
					(* paratop.first).row(para.first)=cameraPosUnitLenthConvert(est_parameters.block(0,para.second.first,1,para.second.second),-1);
				}
				else
				{
					//intit_parameters.block(0,para.second.first,1,para.second.second)=(* para.first.first).row(para.first.second);
					(* paratop.first).row(para.first)=est_parameters.block(0,para.second.first,1,para.second.second);
				}
			}
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
		auto cameraPara=estimateCameraParameter(sphericalFeaturesAndCameras[cameraIndex],curInd1,points,curInd2);

		cout<<"estimated camera parameters of camera number "<<cameraIndex<<": "<<cameraPara<<endl;
		//getchar();
		cameraPosition.row(cameraIndex)=cameraPara;

		alreadyEstimated[cameraIndex]=true;
		reconstructPoints(cameraIndex);




		unordered_map<int,cameraType> bundlePara;

		if(cameraIndex<15)
		{
			bundlePara[0]=cameraType::_static;
			bundlePara[1]=cameraType::_unitLength;
			for (int i = 2; i <=cameraIndex; i++)
			{
				bundlePara[i]=cameraType::_ordinary;
			}
		}
		else
		{
			for (int i = cameraIndex-10; i < cameraIndex-5; i++)
			{
				bundlePara[i]=cameraType::_static;
			}
			for (int i = cameraIndex-5; i <= cameraIndex; i++)
			{
				bundlePara[i]=cameraType::_ordinary;
			}
		}
		cout<<"bundle adjustment for cameras and points before camera number "<<cameraIndex<<endl;
		bundleAdjustment(bundlePara);
		reconstructPoints(cameraIndex);
	}


	return make_pair(cameraPosition,points);
	
}








MatrixXd projectionError(const vector<MatrixXd>& input)
{
//	MatrixXd parameters(1,6);
//	parameters.block(0,0,1,3)=input[0];
//	parameters.block(0,3,1,3)=input[2];
//	MatrixXd variables=input[1];
	assert(input.size()==3);
	return projectionError(input[0],input[1],input[2]);
}

MatrixXd jacobianForCamera(const vector<MatrixXd>& input)
{
//	MatrixXd parameters(1,6);
//	parameters.block(0,0,1,3)=input[0];
//	parameters.block(0,3,1,3)=input[2];
//	MatrixXd variables=input[1];
	assert(input.size()==3);
	return jacobianForCamera(input[0],input[1],input[2]);

}

MatrixXd jacobianForPoint(const vector<MatrixXd>& input)
{
//	MatrixXd parameters(1,6);
//	parameters.block(0,0,1,3)=input[0];
//	parameters.block(0,3,1,3)=input[2];
//	MatrixXd variables=input[1];
	assert(input.size()==3);
	return jacobianForPoint(input[0],input[1],input[2]);

}



MatrixXd projectionErrorUnitLength(const vector<MatrixXd>& input)
{
//	MatrixXd parameters(1,6);
//	parameters.block(0,0,1,3)=input[0];
//	parameters.block(0,3,1,3)=input[2];
//	MatrixXd variables=input[1];
	assert(input.size()==3);
	return projectionErrorUnitLength(input[0],input[1],input[2]);
}

MatrixXd jacobianForCameraUnitLength(const vector<MatrixXd>& input)
{
//	MatrixXd parameters(1,6);
//	parameters.block(0,0,1,3)=input[0];
//	parameters.block(0,3,1,3)=input[2];
//	MatrixXd variables=input[1];
	assert(input.size()==3);
	return jacobianForCameraUnitLength(input[0],input[1],input[2]);

}

MatrixXd jacobianForPointUnitLength(const vector<MatrixXd>& input)
{
//	MatrixXd parameters(1,6);
//	parameters.block(0,0,1,3)=input[0];
//	parameters.block(0,3,1,3)=input[2];
//	MatrixXd variables=input[1];
	assert(input.size()==3);
	return jacobianForPointUnitLength(input[0],input[1],input[2]);

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

	vector<MatrixXd> ddata(2);

	ddata[0]=MatrixXd(ind1.size(),projPoints.cols());
	ddata[1]=MatrixXd(ind1.size(),points.cols());
//	MatrixXd dataset(ind1.size(),6);

	for (int i = 0; i < ind1.size(); i++)
	{
		ddata[0].block(i,0,1,3)=projPoints.row(ind1[i]);
		ddata[1].block(i,3,1,3)=points.row(ind2[i]);
	}

	
	MatrixXd obj_vals=MatrixXd::Zero(1,ind1.size()*3);
	MatrixXd para=MatrixXd::Zero(1,6);

	return	levenbergM_adhoc(ddata,obj_vals,projectionError,jacobianForCamera,para);
}
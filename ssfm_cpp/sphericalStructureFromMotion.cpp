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

	for (size_t i = 0; i < features.size()-1; i++)
	{
		cout<<"processing frame " <<i<<endl;
		for (size_t j = 0; j < markers[i].size(); j++)
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
	MatrixXd rot=MatrixXd::Zero(3,3);
	double r=(a*a+b*b+c*c+4.0);

	rot(0,0) = ((a*a)*2.0+8.0)/r-1.0;
	rot(0,1) = (c*-4.0+a*b*2.0)/r;
	rot(0,2) = (b*4.0+a*c*2.0)/r;
	rot(1,0) = (c*4.0+a*b*2.0)/r;
	rot(1,1) = ((b*b)*2.0+8.0)/r-1.0;
	rot(1,2) = (a*-4.0+b*c*2.0)/r;
	rot(2,0) = (b*-4.0+a*c*2.0)/r;
	rot(2,1) = (a*4.0+b*c*2.0)/r;
	rot(2,2) = ((c*c)*2.0+8.0)/r-1.0;

	return rot;
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

	MatrixXd U=	svd.matrixU();
	MatrixXd V=	svd.matrixV();

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

	//cout<<Tx<<endl;
	//cout<<R1<<endl;
	//cout<<R2<<endl;
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
	MatrixXd re=MatrixXd::Zero(3,4);

	double r=d*d+e*e+f*f;
	re(0,0) = ((d*d)*-2.0)/(r)+2.0;
	re(0,1) = (d*e*-2.0)/(r);
	re(0,2) = (d*f*-2.0)/(r);
	re(0,3) = (d*(b*e*2.0+c*f*2.0)-a*(e*e+f*f)*2.0)/(r);
	re(1,0) = (d*e*-2.0)/(r);
	re(1,1) = ((e*e)*-2.0)/(r)+2.0;
	re(1,2) = (e*f*-2.0)/(r);
	re(1,3) = b*-2.0+(e*(a*d+b*e+c*f)*2.0)/(r);
	re(2,0) = (d*f*-2.0)/(r);
	re(2,1) = (e*f*-2.0)/(r);
	re(2,2) = ((f*f)*-2.0)/(r)+2.0;
	re(2,3) = c*-2.0+(f*(a*d+b*e+c*f)*2.0)/(r);
	return re;
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

	a	=	coefficients(0,0)	;
	b	=	coefficients(0,1)	;
	c	=	coefficients(0,2)	;
	d	=	coefficients(0,3)	;
	e	=	coefficients(1,0)	;
	f	=	coefficients(1,1)	;
	g	=	coefficients(1,2)	;
	h	=	coefficients(1,3)	;
	i	=	coefficients(2,0)	;
	j	=	coefficients(2,1)	;
	k	=	coefficients(2,2)	;
	l	=	coefficients(2,3)	;

	double x,y,z;
	double r=(a*f*k - a*g*j - b*e*k + b*g*i + c*e*j - c*f*i);
	x=-(b*g*l - b*h*k - c*f*l + c*h*j + d*f*k - d*g*j)/r;
	y=(a*g*l - a*h*k - c*e*l + c*h*i + d*e*k - d*g*i)/r;
	z=-(a*f*l - a*h*j - b*e*l + b*h*i + d*e*j - d*f*i)/r;

	MatrixXd result=MatrixXd(1,3);
	result(0,0)=x;
	result(0,1)=y;
	result(0,2)=z;
	return result;

}

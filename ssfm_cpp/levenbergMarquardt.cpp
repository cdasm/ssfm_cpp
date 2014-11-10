#include "levenbergMarquardt.h"


#include <iostream>

MatrixXd levenbergM_simple(MatrixXd& dataset,const MatrixXd& obj_vals,funcType func,funcType jfunc,MatrixXd& initParameters,int maxiter_times)
{
	int dataNumber=dataset.rows();
	int parameterNumber=initParameters.cols();
	int observationNumber=obj_vals.cols()/dataNumber;

	double lambda=0.01;

	bool updateJ=true;
	auto para_est=initParameters;
	auto para_lm=initParameters;
	MatrixXd dis_init(1,dataNumber*observationNumber);

	MatrixXd J;//(dataNumber,parameterNumber);

	MatrixXd d=obj_vals;
	MatrixXd dp(1,parameterNumber);
	MatrixXd H(parameterNumber,parameterNumber);
	MatrixXd H_lm(parameterNumber,parameterNumber);
	double e=-1;

	for (int it = 0; it < maxiter_times; it++)
	{
		if (updateJ)
		{
			J=MatrixXd::Zero(dataNumber*observationNumber,parameterNumber);
			for (int i = 0; i < dataNumber; i++)
			{
				J.block(i*observationNumber,0,observationNumber,parameterNumber)=jfunc(dataset.row(i),para_est);
			}

			for (int i = 0; i < dataNumber; i++)
			{
				dis_init.block(0,i*observationNumber,1,observationNumber)=func(dataset.row(i),para_est);
			}
			d=obj_vals-dis_init;
			H=J.transpose()*J;
			if(it==0)
				e=(d*d.transpose())(0,0);
		}
		H_lm=H+MatrixXd::Identity(parameterNumber,parameterNumber)*lambda;
		dp=d*J*H_lm.inverse();
		para_lm=para_est+dp;
		for (int i = 0; i < dataNumber; i++)
		{
			dis_init.block(0,i*observationNumber,1,observationNumber)=func(dataset.row(i),para_lm);
		}
		d=obj_vals-dis_init;
		double e_lm=(d*d.transpose())(0,0);

		if(e_lm<e)
		{
			lambda/=10;
			para_est=para_lm;

			if(e-e_lm<constrain_on_delta_error)
				break;
			e=e_lm;
			updateJ=true;

			
		}
		else
		{
			updateJ=false;
			lambda*=10;
			if(lambda>lambda_limit)
				break;
		}
	}


	return para_est;
}


vector<Triplet<double> > tosetSparseMatrix(const MatrixXd& dm,int startrow,int startcolumn)
{
	vector<Triplet<double> > toset; 
	for (int i = 0; i < dm.rows(); i++)
	{
		for (int j = 0; j < dm.cols(); j++)
		{
			toset.push_back(Triplet<double>(startrow+i,startcolumn+j,dm(i,j)));
		}
	}
	return toset;
//	sm.setFromTriplets(toset.begin(),toset.end());
}


MatrixXd levenbergM_advanced(MatrixXd& dataset,MatrixXd& assistantPara,const vector<vector<int> >& funcDataMap,const vector<vector<int> >& jfuncDataMap,const MatrixXd& obj_vals,vector<funcType2>& funcs,vector<funcType2>& jfuncs,MatrixXd& initParameters,int maxiter_times)
{
	int dataNumber=dataset.rows();
	int parameterNumber=initParameters.cols();
	int observationNumber=obj_vals.cols()/dataNumber;

	double lambda=0.01;

	bool updateJ=true;
	auto para_est=initParameters;
	auto para_lm=initParameters;
	MatrixXd dis_init(1,dataNumber*observationNumber);

	SparseMatrix<double> J;//(dataNumber*observationNumber,parameterNumber);//(dataNumber,parameterNumber);

	MatrixXd d=obj_vals;
	MatrixXd dp(1,parameterNumber);
	MatrixXd H(parameterNumber,parameterNumber);
	MatrixXd H_lm(parameterNumber,parameterNumber);
	double e=-1;

	auto update_dis_init =[&](const MatrixXd& paraHere)
	{
		dis_init*=0;
		for (int i = 0; i < funcDataMap.size(); i++)
		{
			vector<MatrixXd> curparas;
			int funcind=funcDataMap[i][0];
			int sc;
			sc=funcDataMap[i][1];
			//sc=funcDataMap[i][2];

			curparas.push_back(dataset.row(funcDataMap[i][2]));

			if(funcDataMap[i][3]>0)
				curparas.push_back(assistantPara.row(funcDataMap[i][4]));

				
			for (int j = 5; j < funcDataMap[i].size(); j+=2)
			{
				curparas.push_back(paraHere.block(0,funcDataMap[i][j],1,funcDataMap[i][j+1]));					
			}
			MatrixXd upd=funcs[funcind](curparas);
			//cout<<i<<endl;
			dis_init.block(0,sc,1,upd.cols())=upd;
		}
	};

	for (int it = 0; it < maxiter_times; it++)
	{
		cout<<"iteration Number "<<it<<endl;
		
		if (updateJ)
		{
			J=SparseMatrix<double>(dataNumber*observationNumber,parameterNumber);
			vector<Triplet<double> > tosetJ;
			//for (int i = 0; i < dataNumber; i++)
			//{
			//	J.block(i*observationNumber,0,observationNumber,parameterNumber)=jfunc(dataset.row(i),para_est);
			//}

			for (int i = 0; i < jfuncDataMap.size(); i++)
			{
				vector<MatrixXd> curparas;
				int funcind=jfuncDataMap[i][0];
				int sr, sc;
				sr=jfuncDataMap[i][1];
				sc=jfuncDataMap[i][2];

				curparas.push_back(dataset.row(jfuncDataMap[i][3]));

				if(jfuncDataMap[i][4]>0)
					curparas.push_back(assistantPara.row(jfuncDataMap[i][5]));

				
				for (int j = 6; j < jfuncDataMap[i].size(); j+=2)
				{
					curparas.push_back(para_est.block(0,jfuncDataMap[i][j],1,jfuncDataMap[i][j+1]));					
				}
				MatrixXd upd=jfuncs[funcind](curparas);
				//cout<<upd<<endl;
//				J.block(sr,sc,upd.rows(),upd.cols())=upd;
			//	setSparseMatrix(J,upd,sr,sc);
				auto cursetJ=tosetSparseMatrix(upd,sr,sc);
				tosetJ.insert(tosetJ.end(),cursetJ.begin(),cursetJ.end());
			}
			J.setFromTriplets(tosetJ.begin(),tosetJ.end());
			update_dis_init(para_est);

			cout<<dis_init.transpose()<<endl;
			getchar();

			d=obj_vals-dis_init;
	//		J.transpose();
			SparseMatrix<double> sH=J.transpose()*J;
			H= MatrixXd( sH);
			if(it==0)
				e=(d*d.transpose())(0,0);
		}
		cout<<"current error is:"<<e<<endl;
		getchar();
		H_lm=H+MatrixXd::Identity(parameterNumber,parameterNumber)*lambda;
		dp=d*J*H_lm.inverse();
		//VectorXd tosdj=d*J;
		//VectorXd solution=H_lm.transpose().colPivHouseholderQr().solve(tosdj);
		//dp=solution.transpose();

		para_lm=para_est+dp;
		update_dis_init(para_lm);

		d=obj_vals-dis_init;
		double e_lm=(d*d.transpose())(0,0);

		if(e_lm<e || e<0)
		{
		
			lambda/=10;
			para_est=para_lm;
		if(e-e_lm<constrain_on_delta_error)
				break;
			e=e_lm;
			updateJ=true;
		}
		else
		{
			updateJ=false;
			lambda*=10;
			if(lambda>lambda_limit)
				break;
		}
	}


	return para_est;
}
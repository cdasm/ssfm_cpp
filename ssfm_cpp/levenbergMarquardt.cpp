#include "levenbergMarquardt.h"




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
		}
	}


	return para_est;
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

	MatrixXd J;//(dataNumber,parameterNumber);

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
			J.block(0,sc,1,sc+upd.rows())=upd;
		}
	};

	for (int it = 0; it < maxiter_times; it++)
	{
		if (updateJ)
		{
			J*=0;
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

				curparas.push_back(dataset.row(funcDataMap[i][3]));

				if(jfuncDataMap[i][4]>0)
					curparas.push_back(assistantPara.row(jfuncDataMap[i][5]));

				
				for (int j = 6; j < jfuncDataMap[i].size(); j+=2)
				{
					curparas.push_back(para_est.block(0,jfuncDataMap[i][j],1,jfuncDataMap[i][j+1]));					
				}
				MatrixXd upd=jfuncs[funcind](curparas);
				J.block(sr,sc,sr+upd.rows(),sc+upd.cols())=upd;
			}
			update_dis_init(para_est);


			d=obj_vals-dis_init;
			H=J.transpose()*J;
			if(it==0)
				e=(d*d.transpose())(0,0);
		}
		H_lm=H+MatrixXd::Identity(parameterNumber,parameterNumber)*lambda;
		dp=d*J*H_lm.inverse();
		para_lm=para_est+dp;
		update_dis_init(para_lm);

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
		}
	}


	return para_est;
}
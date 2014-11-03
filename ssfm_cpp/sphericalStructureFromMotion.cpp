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


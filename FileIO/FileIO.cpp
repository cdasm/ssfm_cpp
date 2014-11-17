#include "FileInOut.h"

#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS
#endif
#pragma warning(disable: 4996) /* Disable deprecation */
vector<int>  fileIOclass::InVectorInt(string s)
{
	vector<int> result;
	FILE* fp;
	fp=fopen(s.c_str(),"r");
	int size;
	fscanf(fp,"%d\n",&size);
	result.resize(size,0);
	for (int i=0;i<size;i++)
	{

		fscanf(fp,"%d\n",&result[i]);
	}

	fclose(fp);

	return result;
};
void  fileIOclass::OutVectorInt(string s,vector<int> v)
{
	FILE* fp;
	fp=fopen(s.c_str(),"w");
	int size(v.size());
	fprintf(fp,"%d\n",size);

	for (int i=0;i<size;i++)
	{

		fprintf(fp,"%d\n",v[i]);
	}

	fclose(fp);
};



vector<double>  fileIOclass::InVectorDouble(string s)
{
	vector<double> result;
	FILE* fp;
	fp=fopen(s.c_str(),"r");
	int size;
	fscanf(fp,"%d\n",&size);
	result.resize(size,0);
	for (int i=0;i<size;i++)
	{

		fscanf(fp,"%lf\n",&result[i]);
	}

	fclose(fp);

	return result;
};
void  fileIOclass::OutVectorDouble(string s,vector<double> v)
{
	FILE* fp;
	fp=fopen(s.c_str(),"w");
	int size(v.size());
	fprintf(fp,"%d\n",size);

	for (int i=0;i<size;i++)
	{

		fprintf(fp,"%lf\n",v[i]);
	}

	fclose(fp);
};


vector<string>  fileIOclass::InVectorString(string s)
{
	vector<string> result;
	FILE* fp;
	fp=fopen(s.c_str(),"r");
	int size;
	fscanf(fp,"%d\n",&size);
	result.resize(size,"");
	char tem[100];

	for (int i=0;i<size;i++)
	{

		fscanf(fp,"%s\n",&tem);
		string s(tem);
		result[i]=s;
	}

	fclose(fp);

	return result;
};
void  fileIOclass::OutVectorString(string s,vector<string> v)
{
	FILE* fp;
	fp=fopen(s.c_str(),"w");
	int size(v.size());
	fprintf(fp,"%d\n",size);

	for (int i=0;i<size;i++)
	{

		fprintf(fp,"%s\n",v[i].c_str());
	}

	fclose(fp);
};


vector<vector<int> > fileIOclass::InVectorSInt(string s)
{
	vector<vector<int> > result;
	FILE* fp;
	fp=fopen(s.c_str(),"r");
	int size,s2;

	fscanf(fp,"%d %d\n",&size,&s2);
	result.resize(size,vector<int>(s2,0));
	for (int i=0;i<size;i++)
	{
		for (int j=0;j<s2;j++)
		{
			fscanf(fp,"%d ",&result[i][j]);
		}
		fscanf(fp,"\n");
	}

	fclose(fp);

	return result;
};
void fileIOclass::OutVectorSInt(string s,vector<vector<int> > v,bool secondDim)
{
	FILE* fp;
	fp=fopen(s.c_str(),"w");
	int size,s2(-1);
	size=v.size();
	if (size>0)
	{
		s2=v[0].size();
	}
	if(secondDim)
		fprintf(fp,"%d %d\n",size,s2);
	else
		fprintf(fp,"%d\n",size);
	
	for (int i=0;i<size;i++)
	{
		for (int j=0;j<s2;j++)
		{
			fprintf(fp,"%d ",v[i][j]);
		}
		fprintf(fp,"\n");
	}

	fclose(fp);
};

vector<vector<double> > fileIOclass::InVectorSDouble(string s)
{
	vector<vector<double> > result;
	FILE* fp;
	fp=fopen(s.c_str(),"r");
	int size,s2;

	fscanf(fp,"%d %d\n",&size,&s2);
	result.resize(size,vector<double>(s2,0));
	for (int i=0;i<size;i++)
	{
		for (int j=0;j<s2;j++)
		{
			fscanf(fp,"%lf ",&result[i][j]);
		}
		fscanf(fp,"\n");
	}

	fclose(fp);

	return result;
};
void fileIOclass::OutVectorSDouble(string s,vector<vector<double> > v,bool secondDim)
{
	FILE* fp;
	fp=fopen(s.c_str(),"w");
	int size,s2(-1);
	size=v.size();
	if (size>0)
	{
		s2=v[0].size();
	}
	if(secondDim)
		fprintf(fp,"%d %d\n",size,s2);
	else
		fprintf(fp,"%d\n",size);

	for (int i=0;i<size;i++)
	{
		for (int j=0;j<s2;j++)
		{
			fprintf(fp,"%lf ",v[i][j]);
		}
		fprintf(fp,"\n");
	}

	fclose(fp);
};


vector<vector<int> > fileIOclass::InVectorSInt(string s, int dim)
{
	vector<vector<int> > result;
	FILE* fp;
	fp=fopen(s.c_str(),"r");
	int size,s2;

	fscanf(fp,"%d\n",&size);
	s2=dim;
	result.resize(size,vector<int>(s2,0));
	for (int i=0;i<size;i++)
	{
		for (int j=0;j<s2;j++)
		{
			fscanf(fp,"%d ",&result[i][j]);
		}
		fscanf(fp,"\n");
	}

	fclose(fp);

	return result;
};

vector<vector<double> > fileIOclass::InVectorSDouble(string s,int dim)
{
	vector<vector<double> > result;
	FILE* fp;
	fp=fopen(s.c_str(),"r");
	int size,s2;

	fscanf(fp,"%d\n",&size);
	s2=dim;
	result.resize(size,vector<double>(s2,0));
	for (int i=0;i<size;i++)
	{
		for (int j=0;j<s2;j++)
		{
			fscanf(fp,"%lf ",&result[i][j]);
		}
		fscanf(fp,"\n");
	}

	fclose(fp);

	return result;
};
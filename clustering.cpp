// clustering.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
//#include "LloydsClustering.h"
#include "EMClustering.h"
#include "FineClustering.h"
#include "Problems.h"
#include <vector>
#include <iostream>
#include <conio.h>
using namespace std;

void PrintCluster( vector<vector<EMPoint<double> > > &points)
{
	int x = 0;
	for(int i = 0; i < points.size(); i++)
	{
		cout<<"Cluster "<<i<<" having points : "<<points[i].size()<<endl;
		x += points[i].size();
		double Xm = 0;
		double Ym = 0;
		int size = points[i].size();
		for(int j = 0; j < size; j++)
		{
			cout<<"X:"<<points[i][j].X()<<" Y:"<<points[i][j].Y()<<endl;
			Xm += points[i][j].X();
			Ym += points[i][j].Y();
		}
		cout<<endl;
		cout<<"Xm "<<Xm/size<<" and Ym "<< Ym/size<<endl;
		cout<<endl;
	}
	cout<<"Total Points = "<<x<<endl;
}


int _tmain(int argc, _TCHAR* argv[])
{
	Problems problem;
	//One - 50 wala
	//Two - 10 wala
	vector<LPoint> points = problem.GetProblem(Two);
	int totCluster = 10;
	for(int i = 0; i < points.size(); i++)
	{
		cout<<i<<"- X:"<<points[i].X<<" Y:"<<points[i].Y<<endl;
	}
	cout<<endl;
	vector<EMPoint<double> > ip(points.size());
	for(int i = 0; i < points.size(); i++)
	{
		ip[i].X() = (double)points[i].X;
		ip[i].Y() = (double)points[i].Y;
	}
	char ans = 'e';
	cout<<"Total Clusters desired."<<endl;
	cin>>totCluster;
	cout<<"Which Clustering method you want? Type -"<<endl;
	cout<<"- f for fine clustering"<<endl;
	cout<<"- e for em clustering"<<endl;
	cin>>ans;
	if(ans == 'e')
	{
		EMClustering EMCluster(ip, totCluster, 100);
		cout<<"Simple Version or Filtered Version? Type"<<endl;
		ans = 's';
		cout<<"a) Simple Version - 's'"<<endl;
		cout<<"b) Filtered Version - 'f'"<<endl;
		cin>>ans;
		if(ans == 'f')
		{
			vector<vector<EMPoint<double> > > EMclusteredpts = EMCluster.Cluster();
			PrintCluster(EMclusteredpts);
		}
		else
		{
			vector<vector<EMPoint<double> > > EMclusteredpts = EMCluster.Cluster(Simple);
			PrintCluster(EMclusteredpts);
		}
	}
	else
	{
		FineClustering FCluster;
		vector<vector<EMPoint<double> > > Fclusteredpts = FCluster.Cluster(ip, totCluster);
		PrintCluster(Fclusteredpts);
	}
	_getch();
	return 0;
}


#ifndef _LLOYDSCLUSTERING___
#define _LLOYDSCLUSTERING___

#pragma once
#include <vector>
#include <map>
#include <iostream>
#include "Definitions.h"

using namespace std;


class KMean
{
	// _changeInDistortion is the percentage distortion beyond which the algorithm stops.
	const double _boundaryDistortion;
	int _totalSets;
	const int _totalNumStages;

	vector<EMPoint<double> > _centres;
	vector<EMPoint<double> > _points;
	vector<int> _pointCentres;
	vector<double> _distortion;
	double _totalDistortion;
	multimap<int, int> _sets;
	typedef multimap<int, int>::iterator SetItr;
	typedef pair<SetItr, SetItr> SetPair;
	typedef pair<int, int> InsertPair;
private:
public:
	KMean(int TotalNumberofClusters, vector<EMPoint<double> > points, vector<EMPoint<double> >* centres);
	virtual ~KMean(void);
	vector<vector<EMPoint<double> > > Cluster();

private:
	void CalculateCentres();
	void CalculateInitialCentres();
	void CalculateInitialSets();
	double Distortion(int CentreIndex);
	void CalculateInitialDistortion();
	void CalculateDistortion()
	{
		CalculateInitialDistortion();
	}
	void ReOrganizeClusters();
	double CalculateDistance(EMPoint<double>  point1, EMPoint<double>  point2)
	{
		double X = fabs(point1.X() - point2.X());
		double Y = fabs(point1.Y() - point2.Y());
		return fabs(sqrt(X*X + Y*Y));
	}
	void PrintSets()
	{
		int x = 0;
		for(int i = 0; i < _totalSets; i++)
		{
			SetPair setpair = _sets.equal_range(i);
			SetItr itr;
			cout<<"Set - "<<i<<" Centre = "<<_centres[i].X()<<" "<<_centres[i].Y()<<endl;
			for(itr = setpair.first; itr != setpair.second; itr++)
			{
				int pindex = (*itr).second;
				int X = _points[pindex].X();
				int Y = _points[pindex].Y();
				EMPoint<double>  pt;
				pt.X() = X; pt.Y() = Y;
				cout<<x<<" X "<<X<<" Y "<<Y<<" Distance "<<CalculateDistance(_centres[i], pt)<<endl;
				x++;
			}
		}
		cout<<"Total PTs : "<<x<<endl;
	}
	void PrintCentre()
	{
		for(int i = 0; i < _totalSets; i++)
		{
			cout<<"Centre "<<i<<" "<<_centres[i].X()<<" "<<_centres[i].Y()<<endl;
		}
	}
};



#endif
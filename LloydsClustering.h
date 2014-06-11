#ifndef _LLOYDSCLUSTERING___
#define _LLOYDSCLUSTERING___

#pragma once
#include <vector>
#include <map>
#include <iostream>
#include "Definitions.h"

using namespace std;


class CLloydsClustering
{
	// _changeInDistortion is the percentage distortion beyond which the algorithm stops.
	const double _boundaryDistortion;
	int _totalSets;
	const int _totalNumStages;

	vector<LPoint> _centres;
	vector<LPoint> _points;
	vector<int> _pointCentres;
	vector<double> _distortion;
	int _totalDistortion;
	multimap<int, int> _sets;
	typedef multimap<int, int>::iterator SetItr;
	typedef pair<SetItr, SetItr> SetPair;
	typedef pair<int, int> InsertPair;
private:
public:
	CLloydsClustering(int TotalNumberofClusters, vector<LPoint> points);
	virtual ~CLloydsClustering(void);
	vector<vector<LPoint> > Cluster();

private:
	void CalculateCentres();
	void CalculateInitialCentres();
	void CalculateInitialSets();
	void CalculateInitialSetsNew();
	double Distortion(int CentreIndex);
	void CalculateInitialDistortion();
	void CalculateDistortion()
	{
		CalculateInitialDistortion();
	}
	void ReOrganizeClusters();
	double CalculateDistance(LPoint point1, LPoint point2)
	{
		int X = (point1.X - point2.X);
		int Y = (point1.Y - point2.Y);
		return sqrt((double)X*X + Y*Y);
	}
	void PrintSets()
	{
		int x = 0;
		for(int i = 0; i < _totalSets; i++)
		{
			SetPair setpair = _sets.equal_range(i);
			SetItr itr;
			cout<<"Set - "<<i<<" Centre = "<<_centres[i].X<<" "<<_centres[i].Y<<endl;
			for(itr = setpair.first; itr != setpair.second; itr++)
			{
				int pindex = (*itr).second;
				int X = _points[pindex].X;
				int Y = _points[pindex].Y;
				LPoint pt;
				pt.X = X; pt.Y = Y;
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
			cout<<"Centre "<<i<<" "<<_centres[i].X<<" "<<_centres[i].Y<<endl;
		}
	}
	void PrintDistortion()
	{
		cout<<"Distortion = "<<_totalDistortion<<endl;
	}
	void PrintBB()
	{
		double MinX = DBL_MAX;
		double MaxX = DBL_MIN;
		double MinY = DBL_MAX;
		double MaxY = DBL_MIN;
		for(int i = 0; i < _points.size(); i++)
		{
			if(MinX > _points[i].X)
				MinX = _points[i].X;
			if(MaxX < _points[i].X)
				MaxX = _points[i].X;
			if(MinY > _points[i].Y)
				MinY = _points[i].Y;
			if(MaxY < _points[i].Y)
				MaxY = _points[i].Y;
		}
		cout<<"Xmin "<<MinX<<" Xmax "<<MaxX<<" Ymin "<<MinY<<" Ymax "<<MaxY<<endl;
	}
};



#endif
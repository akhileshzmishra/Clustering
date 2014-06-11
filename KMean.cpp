#include "StdAfx.h"
#include "KMean.h"
#include <cstdlib>			// standard C++ includes
#include <math.h>	
#include <random>
#include <algorithm>

#define DONE 1
#define NOTDONE 0
#define BADINDEX -1


static double RandomInt();
static double RandomGauss();
static double UniformRand(double lo, double hi);
//--------------------------------------------------------------------------------------------------------------
KMean::KMean(int TotalNumberofClusters, vector<EMPoint<double> > points, vector<EMPoint<double> >* centres):
_boundaryDistortion(0.0001),
_totalSets(TotalNumberofClusters),
_totalNumStages(TotalNumberofClusters*TotalNumberofClusters), 
_centres(_totalSets), 
_points(points), 
_pointCentres(points.size()), 
_distortion(_totalSets),
_totalDistortion(0.0)
{
	if(centres == 0)
	{
		CalculateInitialCentres();
	}
	else
	{
		for(int i = 0; i < _totalSets; i++)
		{
			_centres[i] = (*centres)[i];
		}
	}
	CalculateInitialSets();
	CalculateInitialDistortion();
}
KMean::~KMean(void)
{
}
vector<vector<EMPoint<double> > > KMean::Cluster()
{
	double currentDistortionChange = 100.0; 
	int currStage = 0;
	PrintCentre();
	double prevDistortion = _totalDistortion;
	while(currentDistortionChange > _boundaryDistortion && currStage < _totalNumStages)
	{
		PrintSets();
		CalculateCentres();
		PrintCentre();
		ReOrganizeClusters();
		CalculateDistortion();
		currentDistortionChange = (prevDistortion - _totalDistortion)/prevDistortion;
		prevDistortion = _totalDistortion;
		currStage++;
	}
	PrintSets();
	vector<EMPoint<double> > Set;
	vector<vector<EMPoint<double> > > Clusters(_totalSets, Set);
	for(int i = 0; i < _totalSets; i++)
	{
		SetPair setpair = _sets.equal_range(i);
		SetItr itr;
		for(itr = setpair.first; itr != setpair.second; itr++)
		{
			int pindex = (*itr).second;
			EMPoint<double>  point = _points[pindex];
			Clusters[i].push_back(point);
		}
	}
	return Clusters;

}
void KMean::CalculateCentres()
{
	for(int i = 0; i < _totalSets; i++)
	{
		SetPair setpair = _sets.equal_range(i);
		SetItr itr;
		double X = 0.0;
		double Y = 0.0;
		int num = 0;
		for(itr = setpair.first; itr != setpair.second; itr++)
		{
			int pindex = (*itr).second;
			X += _points[pindex].X();
			Y += _points[pindex].Y();
			num++;
		}
		if(num != 0)
		{
			X /= num;
			Y /= num;
			_centres[i].X() = X;
			_centres[i].Y() = Y;
		}
		else
		{
		}
	}
}
void KMean::ReOrganizeClusters()
{
	for(int i = 0; i < _points.size(); i++)
	{
		double MinDist = DBL_MAX;
		int MinIndex = -1;
		for(int j = 0; j < _totalSets; j++)
		{
			double dist = CalculateDistance(_centres[j], _points[i]);
			if(dist < MinDist)
			{
				MinDist = dist;
				MinIndex = j;
			}
		}
		if(MinIndex != -1)
		{
			if(_pointCentres[i] != MinIndex)
			{
				SetPair setpair = _sets.equal_range(MinIndex);
				SetItr itr;
				for(itr = setpair.first; itr != setpair.second; itr++)
				{
					int pindex = (*itr).second;
					if(pindex != i)
						continue;
					_sets.erase(itr);
					_sets.insert(InsertPair(MinIndex, i));
					_pointCentres[i] = MinIndex;
					break;
				}
			}
		}
	}
}
void KMean::CalculateInitialCentres()
{
	double MinX = DBL_MAX;
	double MaxX = DBL_MIN;
	double MinY = DBL_MAX;
	double MaxY = DBL_MIN;
	for(int i = 0; i < _points.size(); i++)
	{
		if(MinX > _points[i].X())
			MinX = _points[i].X();
		if(MaxX < _points[i].X())
			MaxX = _points[i].X();
		if(MinY > _points[i].Y())
			MinY = _points[i].Y();
		if(MaxY < _points[i].Y())
			MaxY = _points[i].Y();
	}
	for(int i = 0; i < _totalSets; i++)
	{
		int X = UniformRand(MinX, MaxX) + RandomGauss();
		_centres[i].X() = X % (int)(MaxX - MinX);
		int Y = UniformRand(MinY, MaxY) + RandomGauss();
		_centres[i].Y() = Y % (int)(MaxY - MinY);
	}
}
void KMean::CalculateInitialSets()
{
	for(int i = 0; i < _points.size(); i++)
	{
		double MinDist = DBL_MAX;
		int MinIndex = -1;
		for(int j = 0; j < _totalSets; j++)
		{
			double dist = CalculateDistance(_centres[j], _points[i]);
			if(dist < MinDist)
			{
				MinDist = dist;
				MinIndex = j;
			}
		}
		if(MinIndex != -1)
		{
			_sets.insert(InsertPair(MinIndex, i));
			_pointCentres[i] = MinIndex;
		}
	}
}
void KMean::CalculateInitialDistortion()
{
	_totalDistortion = 0.0;
	for(int i = 0 ; i < _totalSets; i++)
	{
		_distortion[i] = Distortion(i);
		_totalDistortion += _distortion[i];
	}
}
double KMean::Distortion(int CentreIndex)
{
	SetPair setpair = _sets.equal_range(CentreIndex);
	SetItr itr;
	double Distortion = 0.0;
	for(itr = setpair.first; itr != setpair.second; itr++)
	{
		int PointIndex = (*itr).second;
		Distortion += CalculateDistance(_centres[CentreIndex], _points[PointIndex]);
	}
	return Distortion;
}

















static double MyRandom()
{
	return rand();
}
static double RandomInt()
{
    int j;
	static int kmIdum = 0;
    static double y, maxran, v[98];
    static int iff = 0;
    if (kmIdum < 0 || iff == 0) 
	{
#ifdef WIN32				
		maxran = RAND_MAX;
#else
		unsigned i, k;
		i = 2;
		do 
		{
			k = i;
			i <<= 1;
		} 
		while (i);
		maxran = (double) k;
#endif
 		iff = 1;  
		for (j = 1; j <= 97; j++)	// exercise the system routine
			MyRandom();			// (value intentionally ignored)

		for (j = 1; j <= 97; j++)	// Then save 97 values and a 98th
			v[j] = MyRandom();
		y = MyRandom();
     }
    j = 1 + (int) (97.0 * (y / maxran));
    y = v[j];
    v[j] = MyRandom();			
    return(y / maxran);
}
static double RandomGauss()
{
    static int iset=0;
    static double gset;
    if (iset == 0) 
	{			
		double v1, v2;
		double r = 2.0;
		while (r >= 1.0) 
		{
			v1 = UniformRand(-1, 1);
			v2 = UniformRand(-1, 1);
			r = v1 * v1 + v2 * v2;
		}
		double fac = sqrt(-2.0 * log(r) / r);
		gset = v1 * fac;
		iset = 1;		    	// set flag
		return v2 * fac;
    }
    else 
	{				
		iset = 0;			
		return gset;			
    }
}
static double UniformRand(double lo, double hi)
{
	 return RandomInt()*(hi-lo) + lo;
}

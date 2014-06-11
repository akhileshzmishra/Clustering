#include "StdAfx.h"
#include "LloydsClustering.h"
#include <cstdlib>			// standard C++ includes
#include <math.h>	
#include <random>
#include <algorithm>

#define DONE 1
#define NOTDONE 0
#define BADINDEX -1

//Declarations
class DistanceSet
{
private:
	int _pointIndex;
	double _distance;
public:
	DistanceSet()
	{
		_pointIndex = BADINDEX;
		_distance = 0.0;
	}
	int& Index()
	{
		return _pointIndex;
	}
	double& Distance()
	{
		return _distance;
	}
	bool IsSetGood()
	{
		return ((_pointIndex >= 0) && (_distance >= 0));
	}
};
//for algorithms predicate(X, Y) X > Y returns true. necessary for vectors
bool DistancePredicate(DistanceSet set1, DistanceSet set2)
{
	return (set1.Distance() > set2.Distance());
}
bool PointsSortPredicate(LPoint p1, LPoint p2)
{
	return (p1.X > p2.X);
}

//Standard deviation and random functions taken from gnu- kmclustering algorithm series
static double RandomInt();
static double RandomGauss();
static double UniformRand(double lo, double hi);
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

//--------------------------------------------------------------------------------------------------------------
CLloydsClustering::CLloydsClustering(int TotalNumberofClusters, vector<LPoint> points):
_boundaryDistortion(1.0), 
_totalSets(TotalNumberofClusters),
_totalNumStages(TotalNumberofClusters*2), 
_centres(_totalSets), 
_points(points), 
_pointCentres(points.size()), 
_distortion(_totalSets),
_totalDistortion(0.0)
{
	//CalculateInitialCentres();
	CalculateInitialSetsNew();
	CalculateInitialDistortion();
}
CLloydsClustering::~CLloydsClustering(void)
{
}
vector<vector<LPoint> > CLloydsClustering::Cluster()
{
	double currentDistortionChange = 100.0; 
	int currStage = 0;
	PrintBB();
	PrintCentre();
	while(currentDistortionChange > _boundaryDistortion && currStage < _totalNumStages)
	{
		PrintSets();
		CalculateCentres();
		PrintCentre();
		ReOrganizeClusters();
		double prevDistortion = _totalDistortion;
		_totalDistortion = 0.0;
		CalculateDistortion();
		cout<<"Current Distortion = "<<_totalDistortion<<" previous Distortion = "<<prevDistortion<<endl;
		currentDistortionChange = ((prevDistortion - _totalDistortion)*100.0)/prevDistortion;
		cout<<"Current Distortion Change = "<<currentDistortionChange<<endl;
		currStage++;
	}
	PrintSets();
	vector<LPoint> Set;
	vector<vector<LPoint> > Clusters(_totalSets, Set);
	for(int i = 0; i < _totalSets; i++)
	{
		SetPair setpair = _sets.equal_range(i);
		SetItr itr;
		for(itr = setpair.first; itr != setpair.second; itr++)
		{
			int pindex = (*itr).second;
			LPoint point = _points[pindex];
			Clusters[i].push_back(point);
		}
	}
	return Clusters;

}
void CLloydsClustering::CalculateCentres()
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
			X += _points[pindex].X;
			Y += _points[pindex].Y;
			num++;
		}
		if(num != 0)
		{
			X /= num;
			Y /= num;
			_centres[i].X = X;
			_centres[i].Y = Y;
		}
		else
		{
		}
	}
}

void CLloydsClustering::CalculateInitialSetsNew()
{
	//In this we would first find the first point from the x axis;
	sort(_points.begin(), _points.end(), PointsSortPredicate);
	vector<DistanceSet>distances(_points.size());
	//vector of all the distances between any two points which is sorted by the distance it makes
	vector<vector<distanceSet> >AllDistances(_points.size(), distances);
	vector<int> DistanceInfo(_points.size(), NOTDONE);
	vector<vector<int> >AllDistanceInfo(_points.size(), DistanceInfo);
	for(int i = 0; i < _points.size(); i++)
	{
		for(int j = 0; j < _points.size(); j++)
		{
			if(i == j)
				continue;
			AllDistances[i][j].Distance() = CalculateDistance(_points[i], _points[j]);
			AllDistances[i][j].Index() = j;
			sort(AllDistances[i].begin(), AllDistance[i].end(), DistancePredicate);
		}
	}
	int numberofCentresFound = 0;
	vector<LPoint> CentresFound;
	int pointsProcessed = 0;
	//To tell if a point has already been counted into a set or not
	vector<int> PointsInfo(_points.size(), NOTDONE);
	int currentIndex = 0;
	while(pointsProcessed < _points.size())
	{
		int currentSubIndex = 1;
		int numpts = 1;
		LPoint& currentPt = _points[currentIndex];
		int Xm = currentPt.X;
		int Ym = currentPt.Y;
		LPoint& itrPt = _points[AllDistances[currentIndex][currentSubIndex].Index()];
		numpts++;
		int newXm = (Xm + itrPt.X)/numpts;
		int newYm = (Ym + itrPt.Y)/numpts;
		
		Xm 
	}
	
}
void CLloydsClustering::ReOrganizeClusters()
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
						itr++;
					_sets.erase(itr);
					break;
				}
				_sets.insert(InsertPair(MinIndex, i));

			}
		}
	}
}
void CLloydsClustering::CalculateInitialCentres()
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
	for(int i = 0; i < _totalSets; i++)
	{
		int X = UniformRand(MinX, MaxX) + RandomGauss();
		_centres[i].X = X % (int)(MaxX - MinX);
		int Y = UniformRand(MinY, MaxY) + RandomGauss();
		_centres[i].Y = Y % (int)(MaxY - MinY);
		cout<<"Centre "<< i<<" = "<<X<<" , "<<Y<<endl;
	}
}
void CLloydsClustering::CalculateInitialSets()
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
void CLloydsClustering::CalculateInitialDistortion()
{
	_totalDistortion = 0.0;
	for(int i = 0 ; i < _totalSets; i++)
	{
		_distortion[i] = Distortion(i);
		_totalDistortion += _distortion[i];
	}
}
double CLloydsClustering::Distortion(int CentreIndex)
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
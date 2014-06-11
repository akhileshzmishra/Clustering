#include "StdAfx.h"
#include "EMClustering.h"
#include <cmath>
#include <cstdlib>
#include <random>
#include "KMean.h"

const double IPI = 3.14152;

EMClustering::EMClustering(vector<EMPoint<double> > points, int numofClusters, int numitr):
_numClusters(numofClusters),
_mean(numofClusters),
_coVariance(numofClusters),
_pi(numofClusters),
_likelihood(points.size(), _pi),
_inputPoints(points),
_hiddenVars(points.size(), _pi),
_numIterations(numitr)
{
	_gauBVConst = 1.0/(2* IPI);
	CalculateBB();
	Init();
}
void EMClustering::PrintComps()
{
	cout<<"--- ----- --------- -------- ----------- ------ ----"<<endl;
	for(int i = 0; i < _numClusters; i++)
	{
		cout<<"******************************************"<<endl;
		cout<<"Mean "<<_mean[i].X()<<", "<<_mean[i].Y()<<" Pi "<<_pi[i]<<endl;
		//cout<<"Hidden Vars or weights and likelihood "<<endl;
		//for(int j = 0; j < _inputPoints.size(); j++)
		//{
		//	cout<<_hiddenVars[j][i]<<"/"<<_likelihood[j][i]<<" , ";
		//}
		//cout<<endl;
	}
	cout<<"******************************************"<<endl;
}

EMClustering::~EMClustering(void)
{
}

double EMClustering::CalculateGaussian(EMPoint<double>& point, int CentreIndex)
{
	EMPoint<double> pt(point.X() - _mean[CentreIndex].X(), point.Y() - _mean[CentreIndex].Y());
	CovMat inverse = _coVariance[CentreIndex].Inverse();
	//double inverseDet = fabs(inverse.Determinant());
	double Det = fabs(_coVariance[CentreIndex].Determinant());
	double xsx = pt.X() * (pt.X()*inverse.CoVariance[0][0] + pt.Y()*inverse.CoVariance[1][0]) +
                 pt.Y() * (pt.X()*inverse.CoVariance[0][1] + pt.Y()*inverse.CoVariance[1][1]);
	double exponent = exp(-0.5*xsx);
	//return _gauBVConst*inverseDet*exponent;
	double retval = _gauBVConst*exponent/sqrt(Det);
	return retval;
}

void EMClustering::Estimate()
{
	// update the likelihood
    for (int ptIdx = 0; ptIdx < _inputPoints.size(); ptIdx++)
	{
        for (int clustCentreIdx = 0; clustCentreIdx < _numClusters; clustCentreIdx++)
		{
            _likelihood[ptIdx][clustCentreIdx] = CalculateGaussian(_inputPoints[ptIdx],  clustCentreIdx);
		}
	}
	 // update the hidden variables - These are weights assigned to the points which we do not know
    for (int ptIdx = 0; ptIdx < _inputPoints.size(); ptIdx++)
	{
        double sum = 0.0;
        for (int clustCentreIdx = 0; clustCentreIdx < _numClusters; clustCentreIdx++)
		{
            _hiddenVars[ptIdx][clustCentreIdx] = (_likelihood[ptIdx][clustCentreIdx] * _pi[clustCentreIdx]);
            sum += (_hiddenVars[ptIdx][clustCentreIdx]);
        }
        if (DBLCompare::IsZero(sum, precisionDBL))  // uniform distribution
		{
            for (int clustCentreIdx = 0; clustCentreIdx < _numClusters; clustCentreIdx++)
			{
                _hiddenVars[ptIdx][clustCentreIdx] = 1.0 / (double)_numClusters;
			}
		}
        else                    // normalization
		{
            for (int clustCentreIdx = 0; clustCentreIdx < _numClusters; clustCentreIdx++)
			{
                _hiddenVars[ptIdx][clustCentreIdx] /= sum;
			}
		}
    }

}
void EMClustering:: Maximize()
{
	// Here we would average out the weights of the points of different clusters.
	for (int clustCentreIdx = 0; clustCentreIdx < _numClusters; clustCentreIdx++)	
	{
        double sum = 0.0;
		EMPoint<double> currpt;
        for (int ptIdx = 0; ptIdx < _inputPoints.size(); ptIdx++)
		{
            //_pi[clustCentreIdx] += _hiddenVars[ptIdx][clustCentreIdx];
            sum += _hiddenVars[ptIdx][clustCentreIdx];
			currpt.X() += _hiddenVars[ptIdx][clustCentreIdx]*_inputPoints[ptIdx].X();
			currpt.Y() += _hiddenVars[ptIdx][clustCentreIdx]*_inputPoints[ptIdx].Y();
        }
		if(DBLCompare::Greater(sum, ZERO))
		{
			//now the _pi stores average hidden variable or avg weight for each point in the cluster.
			_pi[clustCentreIdx] = sum /(double)_inputPoints.size();
			//We have to calculate the weighted mean
			_mean[clustCentreIdx].X() = currpt.X() / sum;
			_mean[clustCentreIdx].Y() = currpt.Y() / sum;
			// update the covariance matrix
			double ROxsq = 0.0;
			double ROxy = 0.0;
			double ROysq = 0.0;
			int size = 1;//_inputPoints.size();
            for (int i = 0; i < _inputPoints.size(); i++)
			{
				EMPoint<double> pt(_inputPoints[i].X() - _mean[clustCentreIdx].X(),
					_inputPoints[i].Y() - _mean[clustCentreIdx].Y());
                ROxsq += _hiddenVars[i][clustCentreIdx] * pt.X() * pt.X();//orthox square
                ROxy  += _hiddenVars[i][clustCentreIdx] * pt.X() * pt.Y();//orthox X orthoy
                ROysq += _hiddenVars[i][clustCentreIdx] * pt.Y() * pt.Y();//orthoy square
			}
			_coVariance[clustCentreIdx].CoVariance[0][0] = ROxsq / (sum*size);
            _coVariance[clustCentreIdx].CoVariance[0][1] = ROxy / (sum*(size));
			_coVariance[clustCentreIdx].CoVariance[1][0] = _coVariance[clustCentreIdx].CoVariance[0][1];
            _coVariance[clustCentreIdx].CoVariance[1][1] = ROysq / (sum*size);
		}
		else
		{
			Init(clustCentreIdx);
		}
	}
}
void EMClustering::CalculateBB()
{
	double MinX = DBL_MAX;
	double MaxX = DBL_MIN;
	double MinY = DBL_MAX;
	double MaxY = DBL_MIN;
	for(int i = 0; i < _inputPoints.size(); i++)
	{
		if(MinX > _inputPoints[i].X())
			MinX = _inputPoints[i].X();
		if(MaxX < _inputPoints[i].X())
			MaxX = _inputPoints[i].X();
		if(MinY > _inputPoints[i].Y())
			MinY = _inputPoints[i].Y();
		if(MaxY < _inputPoints[i].Y())
			MaxY = _inputPoints[i].Y();
	}
	_boundingBox.Set(MaxX, MinX, MaxY, MinY);
}
void EMClustering::Init()
{
	_widthRatio = double(RAND_MAX)/_boundingBox.Width();
	_heightRatio = double(RAND_MAX)/_boundingBox.Height();
	_covRatio = double(RAND_MAX)/10000.0;
	for (int clustCentreIdx = 0; clustCentreIdx < _numClusters; clustCentreIdx++)	
	{
		Init(clustCentreIdx);
	}
}
void EMClustering::Init(int cIndex)
{
	_pi[cIndex] = 1.0/(double)_numClusters;
	int Width = (double)_boundingBox.Width();
	if(Width == 0)
		Width = 1;
	int Height = (double)_boundingBox.Height();
	if(Height == 0)
		Height = 1;
	int CovRat = _covRatio;
	if(CovRat == 0)
		CovRat = 1;
    _mean[cIndex].X() = fabs((double)((int)rand() % Width));
    _mean[cIndex].Y() = fabs((double)((int)rand() % Height));
    _coVariance[cIndex].CoVariance[0][0] = 8000;fabs((double)((int)rand() / CovRat));
    _coVariance[cIndex].CoVariance[1][1] = 8000;fabs((double)((int)rand() / CovRat));
    _coVariance[cIndex].CoVariance[0][1] = 0.0;
	_coVariance[cIndex].CoVariance[1][0] = 0.0;
}
vector<vector<EMPoint<double> > >EMClustering::Cluster(TypeofEmCluster type)
{
	for(int i = 0; i < _numIterations; i++)
	{
		Estimate();
		Maximize();
		PrintComps();
	}
	if(type == Filtered)
	{
		KMean SecondStageCluster(_numClusters, _inputPoints, &_mean);
		return SecondStageCluster.Cluster();
	}
	return CalculateSets();
}
vector<vector<EMPoint<double> > > EMClustering::CalculateSets()
{
	vector<vector<EMPoint<double> > >retval(_numClusters);
	for(int i = 0; i < _inputPoints.size(); i++)
	{
		double MinDist = DBL_MAX;
		int MinIndex = -1;
		for(int j = 0; j < _numClusters; j++)
		{
			double dist = CalculateDistance(_mean[j], _inputPoints[i]);
			if(DBLCompare::Less(dist, MinDist))
			{
				MinDist = dist;
				MinIndex = j;
			}
		}
		if(MinIndex != -1)
		{
			retval[MinIndex].push_back(_inputPoints[i]);
		}
	}
	return retval;
}

vector<vector<EMPoint<double> > > EMClustering::CalculateSetsNew()
{
	vector<vector<EMPoint<double> > >retval(_numClusters);
	for(int i = 0; i < _inputPoints.size(); i++)
	{
		double minProb = 0.0;
		int MinIndex = -1;
		for(int j = 0; j < _numClusters; j++)
		{
			double newProb = _hiddenVars[i][j];
			if(DBLCompare::Less(minProb, newProb))
			{
				minProb = newProb;
				MinIndex = j;
			}
		}
		if(MinIndex != -1)
		{
			retval[MinIndex].push_back(_inputPoints[i]);
		}
	}
	return retval;
}

#pragma once

#ifndef _EM_CLUSTERING___
#define _EM_CLUSTERING___

#include "stdafx.h"
#include "Definitions.h"
#include <vector>
#include <math.h>
#include <iostream>
using namespace std;


#define MAX_ITR 100

enum TypeofEmCluster
{
	Simple,
	Filtered
};
const double ZERO = 1e-32;
const int precisionDBL = 8;
//Hidden variables are actually weights assigned to the points.
//The mean we would take is weighted mean.
class EMClustering
{
	const int _numClusters;
	vector<EMPoint<double> > _mean;
	struct CovMat
	{
		double CoVariance[2][2];
		double Determinant()
		{
			double retval = (CoVariance[0][0]* CoVariance[1][1] - 
						  CoVariance[0][1]*CoVariance[1][0]);

			if(DBLCompare::IsZero(retval, precisionDBL))
			{
				CovMat pcov;
				pcov.CoVariance[0][0] = CoVariance[0][0] + 1e-3;
				pcov.CoVariance[0][1] = CoVariance[0][1];
				pcov.CoVariance[1][0] = CoVariance[1][0];
				pcov.CoVariance[1][1] = CoVariance[1][1] + 1e-3;
				return pcov.Determinant();
			}
			else if(DBLCompare::Less(retval, ZERO))
			{
				cout<<"the determinant of covariant is negative"<<endl;
			}
			return retval;
		}
		// compute the inverse matrix of a 3x3 covariance matrix
		CovMat Inverse()
		{
			CovMat retval;
			double determinant = Determinant();
			retval.CoVariance[0][0] = CoVariance[1][1]/determinant;
			retval.CoVariance[1][1] = CoVariance[0][0]/determinant;
			retval.CoVariance[1][0] = -CoVariance[1][0]/determinant;
			retval.CoVariance[0][1] = -CoVariance[0][1]/determinant;
			return retval;
		}
	};
	vector<CovMat> _coVariance;
	vector<double> _pi;
	vector<vector<double> > _likelihood;
	vector<EMPoint<double> > _inputPoints;
	vector<vector<double> >_hiddenVars;

	// constant variables
	double _gauBVConst;

	// Clustering control
	const int _numIterations;
	// Cluster parameters from point
	BB _boundingBox;
	double _widthRatio;
	double _heightRatio;
	double _covRatio;
public:
	EMClustering(vector<EMPoint<double> > points, int numofClusters, int numitr = MAX_ITR);
	~EMClustering(void);
	vector<vector<EMPoint<double> > >Cluster(TypeofEmCluster type = Filtered);
private:
	double CalculateGaussian(EMPoint<double>& point, int CentreIndex);
	void Estimate();
	void Maximize();
	void CalculateBB();
	void Init();
	void Init(int cIndex);
	vector<vector<EMPoint<double> > > CalculateSets();
	vector<vector<EMPoint<double> > > CalculateSetsNew();
	double CalculateDistance(EMPoint<double> point1, EMPoint<double> point2)
	{
		double X = (double)(point1.X() - point2.X());
		double Y = (double)(point1.Y() - point2.Y());
		return sqrt(X*X + Y*Y);
	}

	//Print Functions
	void PrintComps();
};

#endif
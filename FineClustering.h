#pragma once
#ifndef _FINECLUSTERING___
#define _FINECLUSTERING___

#include "Definitions.h"
#include <vector>
#include "EMClustering.h"
#include "KMean.h"
using namespace std;

class FineClustering
{
public:
	FineClustering(void);
	~FineClustering(void);
	vector<vector<EMPoint<double> > > Cluster(vector<EMPoint<double> >& point, int clusters);
};

#endif
#include "StdAfx.h"
#include "FineClustering.h"


FineClustering::FineClustering(void)
{
}


FineClustering::~FineClustering(void)
{
}
vector<vector<EMPoint<double> > > FineClustering::Cluster(vector<EMPoint<double> >& point, int clusters)
{
	vector<vector<EMPoint<double> > > retval;
	EMClustering emcluster(point, clusters, 50);
	vector<vector<EMPoint<double> > > val = emcluster.Cluster(Simple);
	for(int i = 0; i < val.size(); i++)
	{
		for(int j = i + 1; j < val.size(); j++)
		{
			vector<EMPoint<double> > pair(val[i].size() + val[j].size());
			vector<EMPoint<double> > centrepair(2);
			int k = 0;
			for(k = 0; k < val[i].size(); k++)
			{
				pair[k] = val[i][k];
				centrepair[0].Translate(pair[k]);
			}
			k = val[i].size();
			for(int k1 = 0; k1 < val[j].size(); k1++)
			{
				pair[k] = val[j][k1];
				centrepair[1].Translate(pair[k]);
				k++;
			}
			centrepair[0].Reduce(val[i].size());
			centrepair[1].Reduce(val[j].size());
			EMClustering tempcluster(pair, 2, 3);
			//KMean tempcluster(2, pair, &centrepair);
			vector<vector<EMPoint<double> > > res = tempcluster.Cluster();
			val[i].clear();
			val[i] = res[0];
			val[j].clear();
			val[j] = res[1];
		}
	}
	return val;
}
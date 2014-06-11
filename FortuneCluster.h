#pragma once
#include <iostream>
#include <vector>
using namespace std;
struct FPoint
{
	int X;
	int Y;
};
class CFortuneCluster
{
public:
	CFortuneCluster(void);
	virtual ~CFortuneCluster(void);
	void Cluster(vector<FPoint>& points);
};


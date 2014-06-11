#pragma once
#ifndef __POINTPROBLEMS___
#define __POINTPROBLEMS___
#include <vector>
#include "Definitions.h"
using namespace std;

enum ProblemNum
{
	One,
	Two,

	Last
};

class Problems
{
public:
	Problems(void);
	~Problems(void);

private:
	vector<LPoint> CreateProblem1();//50wala
	vector<LPoint> CreateProblem2();//10 wala

public:
	vector<LPoint> GetProblem(ProblemNum num);

};

#endif
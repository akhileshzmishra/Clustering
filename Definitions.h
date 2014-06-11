
#ifndef _CLUSTER_DEFINITIONS___
#define _CLUSTER_DEFINITIONS___

#include "stdafx.h"
#include <cmath>
#include <iostream>
using namespace std;
template<class T>
class EMPoint
{
	T _X;
	T _Y;
public:
	EMPoint():_X(0), _Y(0)
	{
	}
	EMPoint(T X, T Y): _X(X), _Y(Y)
	{
	}
	T& X()
	{
		return _X;
	}
	T& Y()
	{
		return _Y;
	}
	T VectorVal()
	{
		return (_X*_X + _Y*_Y);
	}
	void Translate(T x, T y)
	{
		_X += x;
		_Y += y
	}
	void Translate(EMPoint<T> anotherPt)
	{
		_X += anotherPt.X();
		_Y += anotherPt.Y();
	}
	void Reduce(int amount)
	{
		if(amount <= 0)
			return;

		_X /= amount;
		_Y /= amount;
	}
	
};

class BB
{
	double _MinX;// = DBL_MAX;
	double _MaxX;// = DBL_MIN;
	double _MinY;// = DBL_MAX;
	double _MaxY;// = DBL_MIN;
public:
	BB()
	{
		_MinX = DBL_MAX;
		_MaxX = DBL_MIN;
		_MinY = DBL_MAX;
		_MaxY = DBL_MIN;
	}
	BB(double MaxX, double MinX, double MaxY, double MinY)
	{
		_MinX = MinX;
		_MaxX = MaxX;
		_MinY = MinY;
		_MaxY = MaxY;
	}
	void Set(double MaxX, double MinX, double MaxY, double MinY)
	{
		_MinX = MinX;
		_MaxX = MaxX;
		_MinY = MinY;
		_MaxY = MaxY;
	}
	double Width()
	{
		return fabs(_MaxX - _MinX);
	}
	double Height()
	{
		return fabs(_MaxY - _MinY);
	}
};
struct LPoint
{
	int X;
	int Y;
};

namespace DBLCompare
{
	int 	Compare(double d1, double d2);
	bool	Equals(double d1, double d2);
	bool	Greater(double d1, double d2);
	bool	GreaterOrEqual(double d1, double d2);
	bool	Less(double d1, double d2);
	bool	LessOrEqual(double d1, double d2);
	bool    IsZero(double d, int precision);
};


#endif
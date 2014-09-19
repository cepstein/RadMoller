#ifndef __GENINTEGRATOR_H__
#define __GENINTEGRATOR_H__
#include "TMath.h"
#include "TRandom3.h"
#include "RandGen.h"
#include <assert.h>

class GenIntegrator {

public:
	GenIntegrator(){}
	virtual ~GenIntegrator();

	void setBins(int);

	double 	*getRand2D();
	void 	setArray(int, int, double);

	double setRowVals(double,double);
	double setColVals(double,double);
    void setRandom(RandGen*);

private:

	int nBins;
	double *dimXY;

	RandGen *random;

	double rLow;
	double rHigh;
	double cLow;
	double cHigh;

	double rRange;
	double cRange;

	double getVal(double, double *);

	double getVal2D(double , double , double *);

	double tInt(double *);

	double rVal(double);

	double cVal(double);

	double getRand(double *);

	double cdfInvert(double, double *);

	double* collapseVec(double *);
		
	double* interpolVec(double , int , double *);


};

#endif

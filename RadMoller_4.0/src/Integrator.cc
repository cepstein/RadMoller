#include "GenIntegrator.h"
#include "Riostream.h"

void GenIntegrator::setBins(int bins){
	nBins = bins;
	dimXY = new double[nBins*nBins];
}

GenIntegrator::~GenIntegrator(){
	delete [] dimXY;
}

// TRandom3 *randGen2 = new TRandom3(0);
void GenIntegrator::setRandom(RandGen *rGen){
	random = rGen;
}

double GenIntegrator::setRowVals(double LowR, double HighR){
	rLow = LowR;
	rHigh = HighR;
	rRange = rHigh-rLow;
}

double GenIntegrator::setColVals(double LowC, double HighC){
	cLow = LowC;
	cHigh = HighC;
	cRange = cHigh-cLow;
}

void GenIntegrator::setArray(int r, int c, double val){
	assert(r<nBins&&r>-1);
	assert(c<nBins&&c>-1);
	dimXY[r*nBins+c] = val;
}


double GenIntegrator::getVal(double coord, double *vec){
	int lBin = (int) floor(coord);
	double diff = coord - lBin;
	double val = vec[lBin]+(vec[lBin+1]-vec[lBin])*diff;
	return val;
}

double GenIntegrator::getVal2D(double rCoord, double cCoord, double *vec){//DOES THIS WORK??????
	int lRow = (int) floor(rCoord);
	double diff = rCoord - lRow;
	double *intVec = interpolVec(diff, lRow, vec);
	double val = getVal(cCoord, intVec);
	delete[] intVec;
	return val;
}	

double GenIntegrator::tInt(double *vec){
	double intgrl = 0.;
	for(int i=1; i<nBins; i++){
		intgrl += 0.5*(vec[i-1]+vec[i]);
	}
	return intgrl;
}


double GenIntegrator::rVal(double bin){
	return rLow + (bin*rRange)/(nBins-1.);
}

double GenIntegrator::cVal(double bin){
	return cLow + (bin*cRange)/(nBins-1.);
}

double GenIntegrator::getRand(double *vec){
	double iVal = random->randomOne();
	double coord = cdfInvert(iVal, vec);
	return coord;
}

double GenIntegrator::cdfInvert(double iVal, double *vec){
	double lVal = 0;
	double uVal = nBins-1.;
	double val;

	double total = tInt(vec);

	double *CDF = new double [nBins];
	CDF[0] = 0;
	for(int i=1; i<nBins; i++){
		CDF[i]=0.5*(vec[i]+vec[i-1])/total+CDF[i-1];
	}
	if(0){
		// cout<<"Old method"<<endl;
		while((uVal-lVal)/(uVal+lVal)>1.e-3){
			val = 0.5*(lVal+uVal);
			int ind = floor(val);
			double diff = val - ind;
			double evalCDF = CDF[ind]+diff*0.5*(vec[ind]+getVal(val,vec))/total;

			if(evalCDF<iVal){lVal=val;}
			else if(evalCDF>iVal){uVal=val;}
		}
	}
	if(1){
		int indL;
		// cout<<"Starting"<<endl;
		for (int i = 0; i<nBins;i++){
			if (CDF[i]>iVal){
				indL = i-1;
				break;
			}
		}
		// cout<<"Found bin"<<endl;
		if (vec[indL]==vec[indL+1]){
			val = (iVal-CDF[indL])/(vec[indL]/total)+indL;
			// cout<<"Same"<<endl;
		}
		else if (abs((vec[indL]-vec[indL+1])/(vec[indL]+vec[indL+1]))<1.e-10){
			val = (iVal-CDF[indL])/(vec[indL]/total)+indL;
			// cout<<"Diff ish"<<vec[indL]<<" and "<<vec[indL+1]<<endl;

		}
		else{
			val = (vec[indL]/total-sqrt(pow(vec[indL]/total,2.)-2.*(iVal-CDF[indL])*(vec[indL]/total-vec[indL+1]/total)))/(vec[indL]/total-vec[indL+1]/total)+indL;
			// cout<<"Diff "<<vec[indL]<<" and "<<vec[indL+1]<<endl;
		}
	}
	delete [] CDF;
	return val;
}	
double* GenIntegrator::collapseVec(double *vec){
	double *collapsedVec = new double[nBins];
	for(int i=0;i<nBins;i++){
		double val = tInt(&vec[nBins*i]);
		collapsedVec[i] = val;
	}
	return collapsedVec;
}	
	
double* GenIntegrator::interpolVec(double diff, int lRow, double *vec){ //DOES THIS WORK????????
	double *iVec = new double[nBins];
	for(int i=0;i<nBins;i++){
		iVec[i]=vec[lRow*nBins+i]+diff*(vec[(lRow+1)*nBins+i]-vec[lRow*nBins+i]);
	}
	return iVec;
}	

double* GenIntegrator::getRand2D(){// DOES THIS WORK???????

	double *collapsedVec = collapseVec(dimXY);
	double rCoord = getRand(collapsedVec);
	double totInt = tInt(collapsedVec);

	int lRow = (int) floor(rCoord);
	double diff = rCoord - lRow;
	double *intVec = interpolVec(diff, lRow, dimXY);
	double cCoord = getRand(intVec);
	double val = getVal2D(rCoord,cCoord,dimXY);

	double weight = pow((val*(nBins-1.)*(nBins-1.))/(totInt*rRange*cRange),-1);

	delete [] collapsedVec;
	delete [] intVec;

	double *coords = new double[3];
	coords[0]=rVal(rCoord);
	coords[1]=cVal(cCoord);
	coords[2]=weight;
	return coords;
}

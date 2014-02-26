#include "Riostream.h"
#include "TFile.h"
#include "TFoam.h"
#include "TNtuple.h"
#include "TLorentzVector.h"
#include "TF1.h"
#include "TF2.h"
#include "TH1.h"
#include "TH2D.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TRandom3.h"
#include "TSystem.h"
#include "TStyle.h"
#include <complex>
#include <stdio.h>


long nEve;
const double me = 0.510998910;
const double Ebeam = Tbeam+me;
const double sw2 = 0.2216189;
const double cw2 = 1.-sw2;
const double alpha = 1./137.035999074;

const double pi = 4.0*atan(1.0);
const double twopi = 2.*pi;

const double Pbeam = sqrt(pow(Ebeam,2.)-pow(me,2.));
const double betacm = Pbeam/(Ebeam+me);
const double gammacm = 1./sqrt(1.-pow(betacm,2.));

const double Ecm = gammacm*Ebeam - gammacm*betacm*Pbeam + gammacm*me;
const double Pcm = sqrt(pow(Ecm/2.,2)-pow(me,2.));//momentum of either beam in CM frame
const double Ecmp = Ecm/2.; //Ecm per particle


const double Pcmp = sqrt(pow(Ecmp,2)-pow(me,2.));//momentum of either beam in CM frame
const double r0 = alpha*(1.97e-11)/me;


const double ec= sqrt(4.*pi*alpha); //electron charge
const double mz = 91188.; // Z boson mass
const double mz2 = pow(mz,2.);
const double mw2 = mz2*cw2;

const double gv = -0.25+sw2; //vector coupling
const double ga = 0.25; //axial coupling

const double se = 4.*Ecmp*Ecmp; //Elastic Mandelstam S ("s" was unavailable)
const double dE = dE_frac*sqrt(se);

const double EkMax = (Ecm*Ecm-4.*me*me)/(2.*Ecm);

//Mandelstam T 
double te(double x)
    {
    return -4.*Ecmp*Ecmp*pow(sin(x/2.),2.);
    }

//Mandelstam U
double ue(double x)
    {
    return -4.*Ecmp*Ecmp*pow(cos(x/2.),2.);
    }

double sqr(double x)
    {
    return x*x;
    }

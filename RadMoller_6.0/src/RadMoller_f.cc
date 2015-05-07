//Function implementations for RadMoller_Gen class
#include "RadMoller.h"
#include "Riostream.h"
#include "TF1.h"
#include "Math/WrappedTF1.h"
#include "Math/BrentRootFinder.h"
#include "Math/Functor.h"
#include "Math/RootFinder.h"
#include "TSystem.h"
#include "TComplex.h"
#include "TFoamIntegrand.h"
#include "TFoam.h"
void RadMoller_Gen::setRandom(RandGen *rGen){
    random = rGen;
}


void RadMoller_Gen::SetMoller(){
    mb_flag = 1;
}

void RadMoller_Gen::SetBhabha(){
    mb_flag = 0;
}

void RadMoller_Gen::SetpRes(int res){
    pRes = res;
}

void RadMoller_Gen::SetpCells(int bins){
    pCells = bins;
}

void RadMoller_Gen::SetRadFrac(double rdf){
    radFrac = rdf;
}

// void RadMoller_Gen::SetCM(){
//     CM_flag = 1;
// }
// void RadMoller_Gen::SetLab()){
//     CM_flag = 0;
// }

// double RadMoller_Gen::CMTheta(double theta){
//     return TMath::ATan()

// }

double RadMoller_Gen::symWeight(double tq2, double pq2){
    double symFactor;
    while(pq2<0){
        pq2 = pq2 + twopi;
    }

    if((tq2 < tqrCut1) && (tq2 > tqrCut0) && (pq2 < phiq1) && (pq2 > phiq0)){
        symFactor = 1.;
    }
    else{
        symFactor = 2.;
    }
    // cout<<symFactor<<endl;
    return symFactor;
}


void RadMoller_Gen::SetTCuts(double tqrc0, double tqrc1, double pq0, double pq1){
    if(mb_flag == 1){
        cout<<"Moller Events Requested!"<<endl;
        cout<<"You Requested: "<<endl;
        // cout<<"Min. Photon Theta = "<<tkc0<<endl;
        // cout<<"Max. Photon Theta = "<<tkc1<<endl;
        cout<<"Min Electron_1 Theta = "<<tqrc0<<endl;
        cout<<"Max Electron_1 Theta = "<<tqrc1<<endl;
        cout<<"Min Electron_1 Phi = "<<pq0<<endl;
        cout<<"Max Electron_1 Phi = "<<pq1<<endl;
        
        tkCut0 = 0.;//(tkc0);
        tqrCut0 = tqrc0;
        tkCut1 = TMath::Pi();//(tkc1);
        tqrCut1 = tqrc1;
        phiq0 = pq0;
        phiq1 = pq1;



    }
    else{
        cout<<"Bhabha Events Requested!"<<endl;
        cout<<"You Requested: "<<endl;
        // cout<<"Min. Photon Theta = "<<tkc0<<endl;
        // cout<<"Max. Photon Theta = "<<tkc1<<endl;
        cout<<"Min Positron Theta = "<<tqrc0<<endl;
        cout<<"Max Positron Theta = "<<tqrc1<<endl;
        cout<<"Min Positron Phi = "<<pq0<<endl;
        cout<<"Max Positron Phi = "<<pq1<<endl;

        tkCut0 = 0.;//tkc0;
        tqrCut0 = tqrc0;
        tkCut1 = TMath::Pi();//tkc1;
        tqrCut1 = tqrc1;
        phiq0 = pq0;
        phiq1 = pq1;

    }


}

void RadMoller_Gen::SetECut(double dEf){
    dE_frac = dEf;
}

void RadMoller_Gen::SetLumi(double lum){
    Lumi = lum;
}

void RadMoller_Gen::SetTBeam(double tb){
    Tbeam = tb;
}

//Mandelstam T 
double RadMoller_Gen::te(double x)
    {
    return -4.*Ecmp*Ecmp*pow(sin(x/2.),2.);
    }

//Mandelstam U
double RadMoller_Gen::ue(double x)
    {
    return -4.*Ecmp*Ecmp*pow(cos(x/2.),2.);
    }

double RadMoller_Gen::sqr(double x)
    {
    return x*x;
    }

double RadMoller_Gen::eps1u(TVector3 *q1, TLorentzVector *k){ //me-> units of MeV    
    double Cosa = TMath::Cos(k->Angle(*q1));
    double E_k = k->E();
    return me*(-((Cosa*E_k*sqrt(-pow((2*Ecmp)/me - E_k/me,2) + 
            (pow(Cosa,2)*pow(E_k,2))/pow(me,2) + 
            (4*pow(Ecmp,2)*pow(Ecmp/me - E_k/me,2))/pow(me,2)))/me) + 
     (2*Ecmp*(Ecmp/me - E_k/me)*((2*Ecmp)/me - E_k/me))/me)/
   (pow((2*Ecmp)/me - E_k/me,2) - (pow(Cosa,2)*pow(E_k,2))/pow(me,2));
}

double RadMoller_Gen::eps1l(TVector3 *q1, TLorentzVector *k){ //me-> units of MeV    
    double Cosa = TMath::Cos(k->Angle(*q1));
    double E_k = k->E();

    return me*((Cosa*E_k*sqrt(-pow((2*Ecmp)/me - E_k/me,2) + 
          (pow(Cosa,2)*pow(E_k,2))/pow(me,2) + 
          (4*pow(Ecmp,2)*pow(Ecmp/me - E_k/me,2))/pow(me,2)))/me + 
     (2*Ecmp*(Ecmp/me - E_k/me)*((2*Ecmp)/me - E_k/me))/me)/
   (pow((2*Ecmp)/me - E_k/me,2) - (pow(Cosa,2)*pow(E_k,2))/pow(me,2));
}

double RadMoller_Gen::q1u(TVector3 *q1, TLorentzVector *k){ //me-> units of MeV
    double Cosa = TMath::Cos(k->Angle(*q1));
    double E_k = k->E();

    return me*(sqrt(-pow((2*Ecmp)/me - E_k/me,2) + (pow(Cosa,2)*pow(E_k,2))/pow(me,2) + 
        (4*pow(Ecmp,2)*pow(Ecmp/me - E_k/me,2))/pow(me,2))*
      ((2*Ecmp)/me - E_k/me) - (2*Cosa*Ecmp*E_k*(Ecmp/me - E_k/me))/pow(me,2))/
   (pow((2*Ecmp)/me - E_k/me,2) - (pow(Cosa,2)*pow(E_k,2))/pow(me,2));
}

double RadMoller_Gen::q1l(TVector3 *q1, TLorentzVector *k){ //me-> units of MeV
    double Cosa = TMath::Cos(k->Angle(*q1));
    double E_k = k->E();
    return me*(sqrt(-pow((2*Ecmp)/me - E_k/me,2) + (pow(Cosa,2)*pow(E_k,2))/pow(me,2) + 
        (4*pow(Ecmp,2)*pow(Ecmp/me - E_k/me,2))/pow(me,2))*
      ((-2*Ecmp)/me + E_k/me) - (2*Cosa*Ecmp*E_k*(Ecmp/me - E_k/me))/pow(me,2))/
   (pow((2*Ecmp)/me - E_k/me,2) - (pow(Cosa,2)*pow(E_k,2))/pow(me,2));
}

//Squared Tree-Level Matrix element - does not take m->0 limit
double RadMoller_Gen::M2(double x)
    {
    return 64.0*pow(pi,2.)*pow(alpha,2.)*(pow(me,4.)/pow(te(x),2)*
        ((pow(se,2)+pow(ue(x),2))/(2.*pow(me,4))+4.*te(x)/pow(me,2)-4.0)+pow(me,4)/
        pow(ue(x),2)*((pow(se,2)+pow(te(x),2))/(2.0*pow(me,4))+4.0*ue(x)/pow(me,2)-4.)+
        pow(me,4)/(ue(x)*te(x))*(se/pow(me,2)-2.0)*(se/pow(me,2)-6.0));
    }

//Squared Tree-Level Bhabha Matrix element
double RadMoller_Gen::M2b(double x)
    {
    return 64.0*pow(pi,2.)*pow(alpha,2.)*(pow(me,4.)/pow(te(x),2)*
        ((pow(ue(x),2)+pow(se,2))/(2.*pow(me,4))+4.*te(x)/pow(me,2)-4.0)+pow(me,4)/
        pow(se,2)*((pow(ue(x),2)+pow(te(x),2))/(2.0*pow(me,4))+4.0*se/pow(me,2)-4.)+
        pow(me,4)/(se*te(x))*(ue(x)/pow(me,2)-2.0)*(ue(x)/pow(me,2)-6.0));
    }

Double_t RadMoller_Gen::SoftPhoton_Moller_Integrand(Double_t *x, Double_t *par){
    double var = x[0];

   //  return (sqrt(SD)*(-2 + TD)*log((sqrt(SD) + sqrt(-4 + SD + 4*TD*(1 - var)*var))/
   //      (sqrt(SD) - sqrt(-4 + SD + 4*TD*(1 - var)*var))))/
   //  ((1 - TD*(1 - var)*var)*sqrt(-4 + SD + 4*TD*(1 - var)*var)) + 
   // (sqrt(SD)*(-2 + UD)*log((sqrt(SD) + sqrt(-4 + SD + 4*UD*(1 - var)*var))/
   //      (sqrt(SD) - sqrt(-4 + SD + 4*UD*(1 - var)*var))))/
   //  ((1 - UD*(1 - var)*var)*sqrt(-4 + SD + 4*UD*(1 - var)*var));
        return ((sqrt(SD)*(-2 + TD)*log((sqrt(SD) + sqrt(-4 + SD + 4*TD*(1 - var)*var))/
        (sqrt(SD) - sqrt(-4 + SD + 4*TD*(1 - var)*var))))/
    ((1 - TD*(1 - var)*var)*sqrt(-4 + SD + 4*TD*(1 - var)*var)) + 
   (sqrt(SD)*(-2 + UD)*log((sqrt(SD) + sqrt(-4 + SD + 4*UD*(1 - var)*var))/
        (sqrt(SD) - sqrt(-4 + SD + 4*UD*(1 - var)*var))))/
    ((1 - UD*(1 - var)*var)*sqrt(-4 + SD + 4*UD*(1 - var)*var)));

}

double RadMoller_Gen::SoftPhoton_Moller_Integral(){
// gSystem->Load("libMathMore");
   // cout<<(sqrt(SD)*(-2 + TD)*log((sqrt(SD) + sqrt(-4 + SD + 4*TD*(1 - 0.2)*0.2))/
   //      (sqrt(SD) - sqrt(-4 + SD + 4*TD*(1 - 0.2)*0.2))))/
   //  ((1 - TD*(1 - 0.2)*0.2)*sqrt(-4 + SD + 4*TD*(1 - 0.2)*0.2)) + 
   // (sqrt(SD)*(-2 + UD)*log((sqrt(SD) + sqrt(-4 + SD + 4*UD*(1 - 0.2)*0.2))/
   //      (sqrt(SD) - sqrt(-4 + SD + 4*UD*(1 - 0.2)*0.2))))/
   //  ((1 - UD*(1 - 0.2)*0.2)*sqrt(-4 + SD + 4*UD*(1 - 0.2)*0.2))<<endl;
   TF1 f("Integrand", this,&RadMoller_Gen::SoftPhoton_Moller_Integrand,0,0.5,0,"RadMoller_Gen","SoftPhoton_Moller_Integrand");
   ROOT::Math::WrappedTF1 wf1(f);
   // cout<<f(0.2)<<endl;
   // Create the Integrator
   ROOT::Math::GSLIntegrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVE);
   // Set parameters of the integration
   ig.SetFunction(wf1);
   ig.SetRelTolerance(0.0001);
   return ig.Integral(0, 0.5);
}

double RadMoller_Gen::SoftPhoton_Moller(double x, double dE){
    SD = se/(me*me);
    UD = ue(x)/(me*me);
    TD = te(x)/(me*me);
    double III = SoftPhoton_Moller_Integral();
    // return (2*alpha*(III + (2*sqrt(SD)*log((sqrt(-4 + SD) + sqrt(SD))/2.))/sqrt(-4 + SD) + 
    //    4*log(me/(2.*dE))*(0.5 + ((-2 + SD)*log((sqrt(-4 + SD) + sqrt(SD))/2.))/
    //        (sqrt(-4 + SD)*sqrt(SD)) + 
    //       ((-2 + TD)*log((sqrt(4 - TD) + sqrt(-TD))/2.))/
    //        (sqrt(4 - TD)*sqrt(-TD)) + 
    //       ((-2 + UD)*log((sqrt(4 - UD) + sqrt(-UD))/2.))/
    //        (sqrt(4 - UD)*sqrt(-UD))) + 
    //    ((-4 + 2*SD)*(pow(pi,2)/6. + 
    //         ((4 - SD)*pow(log((sqrt(-4 + SD) + sqrt(SD))/2.),2))/(-4 + SD) + 
    //         log((sqrt(-4 + SD) + sqrt(SD))/2.)*log(-4 + SD) - 
    //         TMath::DiLog((-2 + SD - sqrt(-4*SD + pow(SD,2)))/2.)))/
    //     sqrt(-4*SD + pow(SD,2))))/pi;
    return ((2*alpha*(III + (2*sqrt(SD)*log((sqrt(-4 + SD) + sqrt(SD))/2.))/sqrt(-4 + SD) + 
       4*log(me/(2.*dE))*(0.5 + ((-2 + SD)*log((sqrt(-4 + SD) + sqrt(SD))/2.))/
           (sqrt(-4 + SD)*sqrt(SD)) + 
          ((-2 + TD)*log((sqrt(4 - TD) + sqrt(-TD))/2.))/
           (sqrt(4 - TD)*sqrt(-TD)) + 
          ((-2 + UD)*log((sqrt(4 - UD) + sqrt(-UD))/2.))/
           (sqrt(4 - UD)*sqrt(-UD))) + 
       ((-4 + 2*SD)*(pow(pi,2)/6. + 
            ((4 - SD)*pow(log((sqrt(-4 + SD) + sqrt(SD))/2.),2))/(-4 + SD) + 
            log((sqrt(-4 + SD) + sqrt(SD))/2.)*log(-4 + SD) - 
            TMath::DiLog((-2 + SD - sqrt(-4*SD + pow(SD,2)))/2.)))/
        sqrt(-4*SD + pow(SD,2))))/pi);

}

Double_t RadMoller_Gen::SoftPhoton_Bhabha_Integrand(Double_t *x, Double_t *par){
    double var = x[0];

    return ((TComplex::Sqrt(SD)*(-2 + TD)*TComplex::Log((TComplex::Sqrt(SD) + TComplex::Sqrt(-4 + SD + 4*TD*(1 - var)*var))/
        (TComplex::Sqrt(SD) - TComplex::Sqrt(-4 + SD + 4*TD*(1 - var)*var))))/
    ((1 - TD*(1 - var)*var)*TComplex::Sqrt(-4 + SD + 4*TD*(1 - var)*var)) + 
   (TComplex::Sqrt(SD)*(2 - UD)*TComplex::Log((TComplex::Sqrt(SD) + TComplex::Sqrt(-4 + SD + 4*UD*(1 - var)*var))/
        (TComplex::Sqrt(SD) - TComplex::Sqrt(-4 + SD + 4*UD*(1 - var)*var))))/
    ((1 - UD*(1 - var)*var)*TComplex::Sqrt(-4 + SD + 4*UD*(1 - var)*var))).Rho();
}

double RadMoller_Gen::SoftPhoton_Bhabha_Integral(){
// gSystem->Load("libMathMore");
   // cout<<SoftPhoton_Bhabha_Integrand(0.2)<<endl;
   TF1 f("Integrand", this,&RadMoller_Gen::SoftPhoton_Bhabha_Integrand,0,1,0,"RadMoller_Gen","SoftPhoton_Bhabha_Integrand");
   ROOT::Math::WrappedTF1 wf1(f);
   // Create the Integrator
   ROOT::Math::GSLIntegrator ig(ROOT::Math::IntegrationOneDim::kADAPTIVE);
   // Set parameters of the integration
   ig.SetFunction(wf1);
   ig.SetRelTolerance(0.0001);
   return ig.Integral(0, 0.5);
}

double RadMoller_Gen::SoftPhoton_Bhabha(double x, double dE){
    SD = se/(me*me);
    UD = ue(x)/(me*me);
    TD = te(x)/(me*me);
    double III = SoftPhoton_Bhabha_Integral();
    return ((2*alpha*(III + (2*sqrt(SD)*TComplex::Log((TComplex::Sqrt(-4 + SD) + TComplex::Sqrt(SD))/2.))/TComplex::Sqrt(-4 + SD) + 
       4*log(me/(2.*dE))*(0.5 + ((2 - SD)*TComplex::Log((TComplex::Sqrt(-4 + SD) + TComplex::Sqrt(SD))/2.))/
           (TComplex::Sqrt(-4 + SD)*TComplex::Sqrt(SD)) + 
          ((-2 + TD)*TComplex::Log((TComplex::Sqrt(4 - TD) + TComplex::Sqrt(-TD))/2.))/
           (TComplex::Sqrt(4 - TD)*TComplex::Sqrt(-TD)) + 
          ((2 - UD)*TComplex::Log((TComplex::Sqrt(4 - UD) + TComplex::Sqrt(-UD))/2.))/(TComplex::Sqrt(4 - UD)*TComplex::Sqrt(-UD))
          ) + ((-4 + 2*SD)*(-pow(pi,2)/6. + 
            TMath::Power(TComplex::Log((TComplex::Sqrt(-4 + SD) + TComplex::Sqrt(SD))/2.),2) - 
            TComplex::Log((TComplex::Sqrt(-4 + SD) + TComplex::Sqrt(SD))/2.)*TComplex::Log(-4 + SD) + 
            TMath::DiLog((-2 + SD - TComplex::Sqrt(-4*SD + pow(SD,2)))/2.)))/
        TComplex::Sqrt(-4*SD + pow(SD,2))))/pi).Rho();
}

//Construct the Tree-Level Cross Section (CMS)
double RadMoller_Gen::tree_cs(double x)
    {
    return 0.5*pow(8.0*pi,-2)*M2(x)*pow(Ecm,-2);
    }

double RadMoller_Gen::tree_cs_b(double x)
    {
    return 0.5*pow(8.0*pi,-2)*M2b(x)*pow(Ecm,-2);
    }

//Construct the Soft-Brehmsstralung-Corrected Cross-Section (CMS)
double RadMoller_Gen::mCSfunc(double x, double dE)
    {
        return Lumi*hbarc2*tree_cs(x)*(1+SoftPhoton_Moller(x,dE));//(1.+soft_cs_tsai(x,dE));
    }

double RadMoller_Gen::bCSfunc(double x, double dE)
    {
        //Extra factor of two because e+,e- are distinguishable 
        return 2.0*Lumi*hbarc2*tree_cs_b(x)*(1.+SoftPhoton_Bhabha(x,dE));
    }


//Construct the e-e- Bremsstrahlung Cross Section 
double RadMoller_Gen::bremCS(double E_k, 
  TLorentzVector *q1, TLorentzVector *k)
    {
        double bCS;
        double Cosa = cos(k->Angle(q1->Vect()));
        bCS = (E_k*me)/
        (2.*Ecmp*sqrt(-pow((2.*Ecmp)/me - E_k/me,2.) + 
       (pow(Cosa,2.)*pow(E_k,2.))/pow(me,2.) + 
       (4.*pow(Ecmp,2.)*pow(Ecmp/me - E_k/me,2.))/pow(me,2.))*Pcmp)*pow(q1->P(),2.);
        
        return 0.5*bCS*Lumi*Mh2(p1,p2,q1,k)*pow(2*pi,-5)*pow(me,-3)*hbarc2/32;
    }

//Construct the e-e+ Bremsstrahlung Cross Section 
double RadMoller_Gen::bremCSb(double E_k, 
  TLorentzVector *q1, TLorentzVector *k)
    {
        double bCS;
        double Cosa = cos(k->Angle(q1->Vect()));
        bCS = (E_k*me)/
        (2.*Ecmp*sqrt(-pow((2.*Ecmp)/me - E_k/me,2.) + 
       (pow(Cosa,2.)*pow(E_k,2.))/pow(me,2.) + 
       (4.*pow(Ecmp,2.)*pow(Ecmp/me - E_k/me,2.))/pow(me,2.))*Pcmp)*pow(q1->P(),2.);
        
        return bCS*Lumi*Mh2b(p1,p2,q1,k)*pow(2*pi,-5)*pow(me,-3)*hbarc2/32;
    }


class RadMoller_Gen::TFDISTRAD: public TFoamIntegrand{
public:
  double tqr0;
  RadMoller_Gen* RG;
  TFDISTRAD(RadMoller_Gen *RG0, double thetaq){
    tqr0=thetaq;
    RG = RG0;
  };
  Double_t Density(Int_t nDim, Double_t *Xarg){
    double Ek0 = Xarg[0]*(RG->k0*(1.-1.e-3)-RG->dE)+RG->dE;
    double tk0 = Xarg[1]*RG->pi;
    double phik0 = Xarg[2]*RG->twopi;
    double val = RG->bremInt(Ek0, tk0, phik0, tqr0,0.);
    return val;
  }
};

class RadMoller_Gen::GenRandTR: public TRandom{
public:
  RandGen *RR;
  GenRandTR(Int_t i=0){};
  void SetRandom(RandGen *R0){RR = R0;}
  Double_t Uniform(Double_t x1=1){
      return x1*RR->randomOne();
  }
  Double_t Uniform(Double_t x1, Double_t x2){
    return x1+(x2-x1)*RR->randomOne();
  }
  Double_t Rndm(){
    return RR->randomOne();
  }
  void RndmArray(Int_t n, Double_t *array){
    for(int i=0;i<n;i++){
      array[i] = RR->randomOne();
    }
  }
  void RndmArray(Int_t n, Float_t *array){
    for(int i=0;i<n;i++){
      array[i] = RR->randomOne();
    }
  }
};

double RadMoller_Gen::bremInt(double Ek0, double tk0, double phik0, double tqr0,double pqr0){
    TLorentzVector *kcm0 = new TLorentzVector;
    TLorentzVector *q1cm0 = new TLorentzVector;
    TVector3 *q1r0 = new TVector3(0,0,0);

    q1r0->SetXYZ(sin(tqr0)*cos(pqr0),sin(tqr0)*sin(pqr0),cos(tqr0));
    kcm0->SetPxPyPzE(Ek0*sin(tk0)*cos(phik0),Ek0*sin(tk0)*sin(phik0),Ek0*cos(tk0),Ek0);
    q1cm0->SetPxPyPzE(q1u(q1r0,kcm0)*q1r0->X(),q1u(q1r0,kcm0)*q1r0->Y(),
        q1u(q1r0,kcm0)*q1r0->Z(),eps1u(q1r0,kcm0));

    double val;
    if(mb_flag==1){
        val = bremCS(Ek0,q1cm0,kcm0);
    }
    else if(mb_flag==0){
        val = bremCSb(Ek0,q1cm0,kcm0);
    }
    // cout<<"val: "<<val<<endl;
    if(val<0){

    cout<<"ERROR VAL < 0: "<<val<<endl;
    // cout<<"Ek: "<<Ek0<<endl;
    // cout<<"kcm: "<<kcm0->X()<<", "<<kcm0->Y()<<", "<<kcm0->Z()<<", "<<kcm0->E()<<endl;
    // cout<<"q1cm: "<<q1cm0->X()<<", "<<q1cm0->Y()<<", "<<q1cm0->Z()<<", "<<q1cm0->E()<<endl;
    }
    if(std::isnan(val)){
        cout<<"ERROR NAN in FILL"<<endl;
        cout<<"Ek: "<<Ek0<<endl;
        cout<<"..vs k0: "<<k0<<endl;
        cout<<"kcm: "<<kcm0->X()<<", "<<kcm0->Y()<<", "<<kcm0->Z()<<", "<<kcm0->E()<<endl;
        cout<<"q1cm: "<<q1cm0->X()<<", "<<q1cm0->Y()<<", "<<q1cm0->Z()<<", "<<q1cm0->E()<<endl;

    }
    return val*sin(tk0);

}



void RadMoller_Gen::InitGenerator_RadMoller(){
    CM_flag = 1;
    me = 0.510998910;
    Ebeam = Tbeam+me;
    alpha = 1./137.035999074;
    hbarc2 = pow(1.97326e-11,2.);
    pi = 4.0*atan(1.0);
    twopi = 2.*pi;
    Pbeam = sqrt(pow(Ebeam,2.)-pow(me,2.));
    // betacm = Pbeam/(Ebeam+me);
    // gammacm = 1./sqrt(1.-pow(betacm,2.));
    // Ecm = gammacm*Ebeam - gammacm*betacm*Pbeam + gammacm*me;
    cm = new TLorentzVector(0.,0.,Pbeam,Ebeam+me);
    p1 = new TLorentzVector(0,0,Pbeam,Ebeam);
    p2 = new TLorentzVector(0,0,0,me);
    p1->Boost(-cm->BoostVector());
    p2->Boost(-cm->BoostVector());
    Ecm = p1->E()+p2->E();

    Pcm = sqrt(pow(Ecm/2.,2)-pow(me,2.));//momentum of either 
    Ecmp = Ecm/2.; //Ecm per particle
    Pcmp = sqrt(pow(Ecmp,2)-pow(me,2.));//momentum of either b
    ec = sqrt(4.*pi*alpha); //electron charge
    se = Ecm*Ecm; //Elastic Mandelstam S ("s" was unavail`
    dE = dE_frac*Ecm;
    EkMax = (Ecm*Ecm-4.*me*me)/(2.*Ecm);
    k0 = me*(2*Ecmp*(-1 + Ecmp/me))/((-1 + (2*Ecmp)/me)*me); //units of MeV

    cout<<"EkMax: "<<EkMax<<endl;
    cout<<"Ecm: "<<Ecm<<endl;


    photonInt = new TFoam*[pRes];
    PseRan = new GenRandTR();
    PseRan->SetRandom(random);

    RHO = new TFDISTRAD*[pRes];

    MCvect =new Double_t[3];

    for(int i = 0;i<pRes;i++){//for every tqr
            double tqr0 = tqrCut0 + (tqrCut1-tqrCut0)*double(i)/double(pRes-1);
            double pqr0 = 0.;
            RHO[i] = new TFDISTRAD(this,tqr0);
            photonInt[i] = new TFoam("PhotonIntegrator");
            photonInt[i]->SetkDim(3);         // No. of dimensions, obligatory!
            photonInt[i]->SetnCells(pCells);     // Optionally No. of cells, default=2000
            photonInt[i]->SetRho(RHO[i]);  // Set 2-dim distribution, included below
            photonInt[i]->SetPseRan(PseRan);  // Set random number generator
            photonInt[i]->Initialize();       // Initialize simulator, may take time...
    }

}

double RadMoller_Gen::tqrFunc(double tqr){
    // return 4/(pi*(1 - pow(-1 + (2*tqr)/pi,2))*
    //  (log(pi - tqrCut0) - log(tqrCut0) - log(pi - tqrCut1) + log(tqrCut1)));
    return 4/(3.*pi*(1 - pow(-1 + (2*tqr)/(3.*pi),2))*
     (log(3*pi - tqrCut0) - log(tqrCut0) - log(3*pi - tqrCut1) + log(tqrCut1)));
    // return 1./((1 - tqr)*(log(1 - tqrCut0) - log(1 - tqrCut1)));
}

double RadMoller_Gen::tqrFunc_Moller(double tqr){
    return 4/(pi*(1 - pow(-1 + (2*tqr)/pi,2))*
     (log(pi - tqrCut0) - log(tqrCut0) - log(pi - tqrCut1) + log(tqrCut1)));

    // return 1./((1 - pow(tqr,2))*(-TMath::ATanH(tqrCut0) + TMath::ATanH(tqrCut1)));
}


double RadMoller_Gen::tqriCDF(double rand){
   //  return (pi*pow(pi - tqrCut0,rand)*tqrCut0*pow(tqrCut1,rand))/
   // (pi*pow(tqrCut0,rand)*pow(pi - tqrCut1,rand) - 
   //   pow(tqrCut0,1 + rand)*pow(pi - tqrCut1,rand) + 
   //   pow(pi - tqrCut0,rand)*tqrCut0*pow(tqrCut1,rand));
    return (3*pi*pow(3*pi - tqrCut0,rand)*tqrCut0*pow(tqrCut1,rand))/
   (3*pi*pow(tqrCut0,rand)*pow(3*pi - tqrCut1,rand) - 
     pow(tqrCut0,1 + rand)*pow(3*pi - tqrCut1,rand) + 
     pow(3*pi - tqrCut0,rand)*tqrCut0*pow(tqrCut1,rand));
    return (pow(-pow(tqrCut0,3) + pow(tqrCut1,3),1/3.)*
     pow(-pow(tqrCut0,3) + rand*pow(tqrCut0,3) - rand*pow(tqrCut1,3),
      1/3.))/pow(pow(tqrCut0,3) - pow(tqrCut1,3),1/3.);
    // return 1. - exp(-(rand*log(1 - tqrCut0)) + rand*log(1 - tqrCut1))*(1 - tqrCut0);
}

double RadMoller_Gen::tkFunc(double cos_t, double beta){
    return (1 - pow(cos_t,2))/(pow(1/beta - cos_t,2)*(-4 + (2*log((1 + beta)/(1 - beta)))/beta));
}

double RadMoller_Gen::tkCDF(double cos_t, double beta){
    return -(-beta + (-1 + pow(beta,2))/(-1 - beta) + 2*log(1 + beta))/
    (2.*(2*beta - log((1 + beta)/(1 - beta)))) + 
   (beta*cos_t + (-1 + pow(beta,2))/(-1 + beta*cos_t) + 2*log(1 - beta*cos_t))/
    (2.*(2*beta - log((1 + beta)/(1 - beta))));
}

double RadMoller_Gen::tkInvert(double y, double beta){

    double min = -1;
    double max = 1;
    double integrator;
    for(int i = 0;i<25;i++){
        integrator = (min+max)/2.;
        if (tkCDF(integrator,beta)>y){max = integrator;}
        else if (tkCDF(integrator,beta)<y){min = integrator;}
    }
    // cout<<integrator<<endl;
    return integrator;
}

double RadMoller_Gen::tqriCDF_Moller(double rand){
    return (pi*pow(pi - tqrCut0,rand)*tqrCut0*pow(tqrCut1,rand))/
   (pi*pow(tqrCut0,rand)*pow(pi - tqrCut1,rand) - 
     pow(tqrCut0,1 + rand)*pow(pi - tqrCut1,rand) + 
     pow(pi - tqrCut0,rand)*tqrCut0*pow(tqrCut1,rand));
   //  return (-1 + exp(2*rand*(TMath::ATanH(tqrCut0) - TMath::ATanH(tqrCut1))) - tqrCut0 - 
   //   exp(2*rand*(TMath::ATanH(tqrCut0) - TMath::ATanH(tqrCut1)))*tqrCut0)/
   // (-1 - exp(2*rand*(TMath::ATanH(tqrCut0) - TMath::ATanH(tqrCut1))) - tqrCut0 + 
   //   exp(2*rand*(TMath::ATanH(tqrCut0) - TMath::ATanH(tqrCut1)))*tqrCut0);

}


double RadMoller_Gen::ekFunc(double ek){
    return 2./(exp((2*ek)/EkMax)*(-exp(-2) + exp((-2*dE)/EkMax))*EkMax);
}

double RadMoller_Gen::ekiCDF(double rand){
    return (EkMax*log(exp(2 + (2*dE)/EkMax)/
       (exp(2) - exp(2)*rand + exp((2*dE)/EkMax)*rand)))/2.;
}

void RadMoller_Gen::Generate_Event(){
    

    tqr = tqrCut0+(tqrCut1-tqrCut0)*random->randomOne();
    tqrWeight = sin(tqr)*(tqrCut1-tqrCut0);

    pqr = phiq0+(phiq1-phiq0)*random->randomOne();
    pqrWeight = phiq1-phiq0;

    pickProc = random->randomOne();

    delete q1cm;
    delete q2cm;
    delete q1;
    delete q2;

    if (pickProc<radFrac){//Bremsstrahlung

        Ek   = dE+(EkMax-dE)*random->randomOne();//ekiCDF(random->randomOne());//
        ekWeight   = EkMax-dE;//1./ekFunc(Ek);//
        // ekWeight = 1.;
        delete qcm;
        delete q1r;
        delete k;
        delete kcm;


        q1r = new TVector3(sin(tqr)*cos(pqr),sin(tqr)*sin(pqr),cos(tqr));
        // kcm = new TLorentzVector(Ek*sin(tk)*cos(phik),Ek*sin(tk)*sin(phik),Ek*cos(tk),Ek);

        if (Ek>k0){

            double aRand = random->randomOne();
            delete newAxis;
            newAxis = new TVector3(-sin(tqr)*cos(pqr),-sin(tqr)*sin(pqr),-cos(tqr));
            cosaMax = -((sqrt(pow((2*Ecmp)/me - Ek/me,2) - 
                (4*pow(Ecmp,2)*pow(Ecmp/me - Ek/me,2))/pow(me,2))*me)/Ek);

            if(1){//Uniform
                        tk  = TMath::ACos(1.-(1.+cosaMax)*random->randomOne());
                        phik = twopi*random->randomOne();

                        kcm = new TLorentzVector(Ek*sin(tk)*cos(phik),Ek*sin(tk)*sin(phik),Ek*cos(tk),Ek);
                        kcm->RotateUz(*newAxis);

                        tkWeight  = (1.+cosaMax);
                        phikWeight = twopi;

                        if (aRand>0.5){
                            q1cm = new TLorentzVector(q1u(q1r,kcm)*q1r->X(),q1u(q1r,kcm)*q1r->Y(),
                                q1u(q1r,kcm)*q1r->Z(),eps1u(q1r,kcm));   
                        }
                        else if (aRand <0.5){
                            q1cm = new TLorentzVector(q1l(q1r,kcm)*q1r->X(),q1l(q1r,kcm)*q1r->Y(),
                                q1l(q1r,kcm)*q1r->Z(),eps1l(q1r,kcm));
                        }
            }



            aFlag = 2.; //weight up by 2.
            q2cm = new TLorentzVector(*p1+*p2-*q1cm-*kcm);
            if(mb_flag==1){//Moller
                weight = bremCS(Ek,q1cm,kcm)*symWeight(q2cm->Theta(),q2cm->Phi())/radFrac\
                *ekWeight*tkWeight*phikWeight*tqrWeight*pqrWeight*aFlag;
                // cout<<"ek "<<ekWeight<<" tkw "<<tkWeight<<" pkw "<<phikWeight<<" tqrw "<<tqrWeight<<" pqrw "<<pqrWeight<<endl;
            }

            if(mb_flag==0){//Bhabha
                weight = bremCSb(Ek,q1cm,kcm)/radFrac\
                *ekWeight*tkWeight*phikWeight*tqrWeight*pqrWeight*aFlag;
            }


        }

        else if (Ek<k0){
            kcm = new TLorentzVector;
            q1cm = new TLorentzVector;
            double tBinD = (tqr-tqrCut0)/(tqrCut1-tqrCut0)*(double(pRes)-1.);
            int tBin = (int) floor(tBinD+0.5);


            photonInt[tBin]->MakeEvent();            // generate MC event
            photonInt[tBin]->GetMCvect(MCvect);

            double EkF = MCvect[0]*(k0*(1.-1.e-3)-dE)+dE;
            tk = MCvect[1]*pi;
            phik = MCvect[2]*twopi+pqr;


            kcm->SetPxPyPzE(EkF*sin(tk)*cos(phik),EkF*sin(tk)*sin(phik),EkF*cos(tk),EkF);
            q1cm->SetPxPyPzE(q1u(q1r,kcm)*q1r->X(),q1u(q1r,kcm)*q1r->Y(),
                q1u(q1r,kcm)*q1r->Z(),eps1u(q1r,kcm));

            // delete intPhoton;
            // photonInt[tBin]->GetMCwt(MCwt);
            // cout<<"WEIGHT: "<<MCwt<<endl;
            photonInt[tBin]->GetIntegMC(MCResult,MCError);
            Double_t *XX = new Double_t[3];
            XX[0] = MCvect[0];
            XX[1] = MCvect[1];
            XX[2] = MCvect[2];
            aFlag = 1.;  //no reweighting
            q2cm = new TLorentzVector(*p1+*p2-*q1cm-*kcm);
            if(mb_flag==1){//Moller
              weight = MCResult*(bremCS(EkF,q1cm,kcm)/RHO[tBin]->Density(3,XX))*\
              symWeight(q2cm->Theta(),q2cm->Phi())/radFrac*tqrWeight*pqrWeight*ekWeight*pi*twopi*sin(tk)*aFlag;
            }//sin(tk)*

            if(mb_flag==0){//Bhabha
              weight = MCResult*(bremCSb(EkF,q1cm,kcm)/RHO[tBin]->Density(3,XX))*\
              1./radFrac*tqrWeight*pqrWeight*ekWeight*pi*sin(tk)*twopi*aFlag;
            }

         }


        q1 = new TLorentzVector(*q1cm);
        q2 = new TLorentzVector(*q2cm);
        k = new TLorentzVector(*kcm);

        q1->Boost(cm->BoostVector());
        q2->Boost(cm->BoostVector());
        k->Boost(cm->BoostVector());

        elFlag = 1;
    }

    if (pickProc >radFrac){//Elastic Kinematics

        q1cm = new TLorentzVector(Pcmp*sin(tqr)*cos(pqr),\
            Pcmp*sin(tqr)*sin(pqr),Pcmp*cos(tqr),Ecmp);
        q2cm = new TLorentzVector(-Pcmp*sin(tqr)*cos(pqr),\
            -Pcmp*sin(tqr)*sin(pqr),-Pcmp*cos(tqr),Ecmp);

        q1 = new TLorentzVector(*q1cm);
        q2 = new TLorentzVector(*q2cm);
        if(mb_flag==1){//Moller
        weight = mCSfunc(tqr,dE)*symWeight(q2cm->Theta(),q2cm->Phi())*tqrWeight*pqrWeight/(1.-radFrac);//
        }

        if(mb_flag==0){//Bhabha
        weight = bCSfunc(tqr,dE)*tqrWeight*pqrWeight/(1.-radFrac);
        }

        q1->Boost(cm->BoostVector());
        q2->Boost(cm->BoostVector());
        elFlag = 0;
    }



}

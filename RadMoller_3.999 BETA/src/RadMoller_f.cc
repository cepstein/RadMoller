//Function implementations for RadMoller_Gen class
#include "RadMoller.h"
#include "Riostream.h"
#include "TF1.h"
#include "Math/WrappedTF1.h"
#include "Math/BrentRootFinder.h"
#include "Math/Functor.h"
#include "Math/RootFinder.h"
#include "GenIntegrator.h"

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

void RadMoller_Gen::SetpBins(int bins){
    pBins = bins;
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

//Tsai's Soft Corrections - Moller
double RadMoller_Gen::soft_cs_tsai(double x, double dE){
    return (pow(te(x),2)*pow(ue(x),2)*(-((pow(me,2)*
            (pow(te(x),4) + pow(ue(x),4) + pow(se,2)*pow(te(x) + ue(x),2)))/
          (pow(te(x),2)*pow(ue(x),2))) - 
       (alpha*pow(me,2)*(2*pow(se,2)*te(x)*(se + te(x))*log(-1 - se/te(x)) - 
            (se + te(x))*(2*se + te(x))*(2*pow(se,2) + 3*se*te(x) + 2*pow(te(x),2))*
             pow(log(-1 - se/te(x)),2) + 
            2*te(x)*pow(se + te(x),3)*log(-(se/te(x))) + 
            (se + te(x))*(2*se + te(x))*(2*pow(se,2) + se*te(x) + pow(te(x),2))*
             pow(log(-(se/te(x))),2) + 
            te(x)*log(-(te(x)/(se + te(x))))*
             (2*pow(se,2)*(se + te(x)) + 
               (pow(se,3) + se*pow(te(x),2) - 2*pow(te(x),3))*
                log(-(te(x)/(se + te(x))))) + 
            2*pow(te(x),3)*(se + te(x))*log(1/(1 + te(x)/se)) + 
            (-2*pow(se,3)*te(x) + pow(se,2)*pow(te(x),2) + pow(te(x),4))*
             pow(log(1/(1 + te(x)/se)),2)))/
        (2.*pi*pow(te(x),2)*pow(se + te(x),2)) + 
       (pow(me,2)*(pow(se,2) + (pow(se,2)*te(x))/ue(x) + pow(ue(x),2))*
          (1 + (alpha*(-46 + 33*log(-(te(x)/pow(me,2))) - 
                 36*log(Ecmp/dE)*(-1 + log((te(x)*ue(x))/(pow(me,2)*se)))))/(9.*pi)
            ))/pow(te(x),2) + (pow(me,2)*
          (pow(se,2) + pow(te(x),2) + (pow(se,2)*ue(x))/te(x))*
          (1 + (alpha*(-46 + 33*log(-(ue(x)/pow(me,2))) - 
                 36*log(Ecmp/dE)*(-1 + log((te(x)*ue(x))/(pow(me,2)*se)))))/(9.*pi)
            ))/pow(ue(x),2)))/
   (pow(me,2)*(pow(te(x),4) + pow(ue(x),4) + pow(se,2)*pow(te(x) + ue(x),2)));
}

//Soft corrections to Bhabha scattering (A.B. Arbuzov, E.S. Scherbakova) & Glover, Tausk, van der Bij
double RadMoller_Gen::soft_bhabha(double x, double dE){
    return (alpha*(-2*(6 + pow(pi,2)) - 6*TMath::DiLog(-(te(x)/se)) + 6*TMath::DiLog(1 + te(x)/se) + 
       9*log(se/pow(me,2)) + 12*log(dE/Ecmp)*
        (-1 + log(se/pow(me,2)) + log(-(te(x)/se)) + log(1/(1 + te(x)/se))) + 
       6*log(se/pow(me,2))*(log(-1 - se/te(x)) + log(-(te(x)/(se*(1 + te(x)/se))))) + 
       (3*((pow(pi,2)*(4 + (8*te(x))/se + (27*pow(te(x),2))/pow(se,2) + 
                 (26*pow(te(x),3))/pow(se,3) + (16*pow(te(x),4))/pow(se,4)))/
             12. + ((6 + (8*te(x))/se + (9*pow(te(x),2))/pow(se,2) + 
                 (3*pow(te(x),3))/pow(se,3))*log(-(te(x)/se)))/2. - 
            (te(x)*(3 + te(x)/se - (3*pow(te(x),2))/pow(se,2) - 
                 (4*pow(te(x),3))/pow(se,3))*pow(log(-(te(x)/se)),2))/(4.*se) + 
            (te(x)*(1 + pow(te(x),2)/pow(se,2))*log(1 + te(x)/se))/(2.*se) + 
            ((4 + (8*te(x))/se + (7*pow(te(x),2))/pow(se,2) + 
                 (2*pow(te(x),3))/pow(se,3))*log(-(te(x)/se))*log(1 + te(x)/se))/2.\
             + ((-2 - (5*te(x))/se - (7*pow(te(x),2))/pow(se,2) - 
                 (5*pow(te(x),3))/pow(se,3) - (2*pow(te(x),4))/pow(se,4))*
               pow(log(1 + te(x)/se),2))/2.))/
        pow(1 + te(x)/se + pow(te(x),2)/pow(se,2),2)))/(3.*pi);
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
        return Lumi*hbarc2*tree_cs(x)*(1.+soft_cs_tsai(x,dE));
    }

double RadMoller_Gen::bCSfunc(double x, double dE)
    {
        //Extra factor of two because e+,e- are distinguishable 
        return 2.0*Lumi*hbarc2*tree_cs_b(x)*(1.+soft_bhabha(x,dE));
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


void RadMoller_Gen::InitGenerator_RadMoller(){
    CM_flag = 0;
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
    se = Ecm*Ecm; //Elastic Mandelstam S ("s" was unavail
    dE = dE_frac*Ecm;
    EkMax = (Ecm*Ecm-4.*me*me)/(2.*Ecm);
    k0 = me*(2*Ecmp*(-1 + Ecmp/me))/((-1 + (2*Ecmp)/me)*me); //units of MeV

    cout<<"EkMax: "<<EkMax<<endl;
    cout<<"Ecm: "<<Ecm<<endl;

    // pRes = 25;
    // pBins = 75;
    // cout<<"pRes: "<<pRes<<"pBins: "<<pBins<<endl;
    photonArray = new GenIntegrator*[pRes*pRes];
    TLorentzVector *kcm0 = new TLorentzVector;
    TLorentzVector *q1cm0 = new TLorentzVector;
    TVector3 *q1r0 = new TVector3(0,0,0);
    cout<<"Integrating Photons..."<<endl;

    for(int i = 0; i<pRes; i++){ //for every Ek, i
        double Ek0 = dE+(k0*(1.-1.e-3)-dE)*double(i)/double(pRes-1);
        for(int j = 0;j<pRes;j++){//for every tqr, i
            double tqr0 = tqrCut0 + (tqrCut1-tqrCut0)*double(j)/double(pRes-1);
            double pqr0 = 0.; //phi always 0 (use symmetry to simplify)
            photonArray[i*pRes+j] = new GenIntegrator;
            // cout<<"1!"<<endl;
            photonArray[i*pRes+j]->setBins(pBins);

            photonArray[i*pRes+j]->setRandom(random);

            // cout<<"2!"<<endl;
            photonArray[i*pRes+j]->setRowVals(0,pi);
            // cout<<"3!"<<endl;
            photonArray[i*pRes+j]->setColVals(0,2.*pi);
            // cout<<"4!"<<endl;
            q1r0->SetXYZ(sin(tqr0)*cos(pqr0),sin(tqr0)*sin(pqr0),cos(tqr0));
            // cout<<"Set"<<endl;
                    for(int k=0;k<pBins;k++){
                        for(int l=0;l<pBins;l++){
                            double tk0 = double(k)/double(pBins-1)*pi;
                            double phik0 = double(l)/double(pBins-1)*2.*pi;
                            // cout<<"fill"<<endl;
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
                            photonArray[i*pRes+j]->setArray(k,l,val);
                            // cout<<"5!"<<endl;
                        }
                    }


        }
    cout<<"[]"<<flush;

    }
    cout<<endl;
    cout<<"Finished photon integration."<<endl;
    delete kcm0;
    delete q1r0;
    delete q1cm0;
}

// double RadMoller_Gen::tkFunc(double tk){
//     // return 4/(pi*(1 - pow(-1 + (2*tk)/pi,2))*
//     //  (log(pi - tkCut0) - log(tkCut0) - log(pi - tkCut1) + log(tkCut1)));
//     return 1./((1 - pow(tk,2))*(-TMath::ATanH(tkCut0) + TMath::ATanH(tkCut1)));

// }

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

// double RadMoller_Gen::tkiCDF(double rand){
//     return (-1 + exp(2*rand*(TMath::ATanH(tkCut0) - TMath::ATanH(tkCut1))) - tkCut0 - 
//      exp(2*rand*(TMath::ATanH(tkCut0) - TMath::ATanH(tkCut1)))*tkCut0)/
//    (-1 - exp(2*rand*(TMath::ATanH(tkCut0) - TMath::ATanH(tkCut1))) - tkCut0 + 
//      exp(2*rand*(TMath::ATanH(tkCut0) - TMath::ATanH(tkCut1)))*tkCut0);
// }

// double RadMoller_Gen::tkFunc(double tk){
//     return 4/(pi*(1 - pow(-1 + (2*tk)/pi,2))*
//      (log(pi - tkCut0) - log(tkCut0) - log(pi - tkCut1) + log(tkCut1)));
// }

// double RadMoller_Gen::tkiCDF(double rand){
//     return (pi*pow(pi - tkCut0,rand)*tkCut0*pow(tkCut1,rand))/
//    (pi*pow(tkCut0,rand)*pow(pi - tkCut1,rand) - 
//      pow(tkCut0,1 + rand)*pow(pi - tkCut1,rand) + 
//      pow(pi - tkCut0,rand)*tkCut0*pow(tkCut1,rand));
// }

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
    
    // if(mb_flag==1){//Moller
    //     tqr  = (tqriCDF_Moller(random->randomOne()));//tqrCut0+(tqrCut1-tqrCut0)*random->randomOne();//
    //     tqrWeight  = sin(tqr)/tqrFunc_Moller(tqr);//tqrCut1-tqrCut0;//1./tqrFunc(tqr);

    // }
    // else{//Bhabha
    //     tqr  = (tqriCDF(random->randomOne()));//tqrCut0+(tqrCut1-tqrCut0)*random->randomOne();//
    //     tqrWeight  = sin(tqr)/tqrFunc(tqr);//tqrCut1-tqrCut0;//1./tqrFunc(tqr);

    //     // tqr  = TMath::ACos(tqriCDF(random->randomOne()));//tqrCut0+(tqrCut1-tqrCut0)*random->randomOne();//
    //     // tqrWeight  = 1./tqrFunc(TMath::Cos(tqr));//tqrCut1-tqrCut0;//1./tqrFunc(tqr);
    // // tqr = TMath::ACos(tqrCut0+(tqrCut1-tqrCut0)*random->randomOne());
    // // tqrWeight = tqrCut1-tqrCut0;
    // }

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
        ekWeight   =EkMax-dE;//1./ekFunc(Ek);//
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
                        // cout<<"cosa "<<cos(kcm->Angle(*q1r))<<" max "<<cosaMax<<endl;
                        // tk  = (pi-TMath::ACos(cosaMax))*random->randomOne();

                        // phik = twopi*random->randomOne();

                        // kcm = new TLorentzVector(Ek*sin(tk)*cos(phik),Ek*sin(tk)*sin(phik),Ek*cos(tk),Ek);
                        // kcm->RotateUz(*newAxis);

                        // tkWeight  = sin(tk)*(pi-TMath::ACos(cosaMax));//(1.+cosaMax)*sin(tk);
                        // phikWeight = twopi;

                        if (aRand>0.5){
                            q1cm = new TLorentzVector(q1u(q1r,kcm)*q1r->X(),q1u(q1r,kcm)*q1r->Y(),
                                q1u(q1r,kcm)*q1r->Z(),eps1u(q1r,kcm));   
                        }
                        else if (aRand <0.5){
                            q1cm = new TLorentzVector(q1l(q1r,kcm)*q1r->X(),q1l(q1r,kcm)*q1r->Y(),
                                q1l(q1r,kcm)*q1r->Z(),eps1l(q1r,kcm));
                        }
            }

            if(0){//Integrator
                // cout<<"HI"<<endl;
                double tkM = pi-TMath::ACos(cosaMax);

                kcm = new TLorentzVector;
                q1cm = new TLorentzVector;

                GenIntegrator *intPhoton = new GenIntegrator;
                intPhoton->setBins(pBins);
                intPhoton->setRandom(random);
                intPhoton->setRowVals(0,(1.-1.e-4)*tkM);
                intPhoton->setColVals(0,2.*pi);

                for(int i = 0;i<pBins;i++){
                    for(int j=0;j<pBins;j++){

                        double tk0 = double(i)/double(pBins-1)*(1.-1.e-4)*tkM;
                        double phik0 = double(j)/double(pBins-1)*2.*pi;

                        kcm->SetPxPyPzE(Ek*sin(tk0)*cos(phik0),Ek*sin(tk0)*sin(phik0),Ek*cos(tk0),Ek);
                        kcm->RotateUz(*newAxis);
                        if(aRand>0.5){
                            q1cm->SetPxPyPzE(q1u(q1r,kcm)*q1r->X(),q1u(q1r,kcm)*q1r->Y(),
                                q1u(q1r,kcm)*q1r->Z(),eps1u(q1r,kcm));
                        }
                        if(aRand<0.5){
                            q1cm->SetPxPyPzE(q1l(q1r,kcm)*q1r->X(),q1l(q1r,kcm)*q1r->Y(),
                                q1l(q1r,kcm)*q1r->Z(),eps1l(q1r,kcm));
                        }

                        double val = bremCS(Ek,q1cm,kcm);
                        intPhoton->setArray(i,j,val);
                    }
                }
                photonCoords = intPhoton->getRand2D();

                tk = photonCoords[0];
                phik = photonCoords[1];
                tkWeight = photonCoords[2]*sin(tk);
                phikWeight = 1.;

                kcm->SetPxPyPzE(Ek*sin(tk)*cos(phik),Ek*sin(tk)*sin(phik),Ek*cos(tk),Ek);
                kcm->RotateUz(*newAxis);

                if (aRand>0.5){
                    q1cm->SetPxPyPzE(q1u(q1r,kcm)*q1r->X(),q1u(q1r,kcm)*q1r->Y(),
                        q1u(q1r,kcm)*q1r->Z(),eps1u(q1r,kcm));   
                }
                else if (aRand <0.5){
                    q1cm->SetPxPyPzE(q1l(q1r,kcm)*q1r->X(),q1l(q1r,kcm)*q1r->Y(),
                        q1l(q1r,kcm)*q1r->Z(),eps1l(q1r,kcm));
                }

                delete intPhoton;
            }

            aFlag = 2.; //weight up by 2.


        }

        else if (Ek<k0){
            // phik = twopi*random->randomOne();
            // phikWeight = twopi;
            // tk   = tkCut0+(tkCut1-tkCut0)*random->randomOne();//tkiCDF(random->randomOne());//
            // tkWeight = sin(tk)*(tkCut1-tkCut0);
            // // // newAxis = new TVector3(sin(tqr)*cos(pqr),sin(tqr)*sin(pqr),cos(tqr));

                

            // phik = twopi*random->randomOne();
            // phikWeight= twopi;
                                                    

            kcm = new TLorentzVector;
            q1cm = new TLorentzVector;

            // int bins = 5;
            // GenIntegrator *intPhoton = new GenIntegrator;
            // intPhoton->setBins(bins);

            // for(int i =0;i<bins;i++){
            //     for(int j=0;j<bins;j++){
            //         intPhoton->setRowVals(0,pi);
            //         intPhoton->setColVals(0,2.*pi);
            //         double tk0 = double(i)/double(bins-1)*pi;
            //         double phik0 = double(j)/double(bins-1)*2.*pi;
            //         // cout<<"tk: "<<tk<<endl;
            //         // cout<<"phik: "<<phik<<endl;

            //         kcm->SetPxPyPzE(Ek*sin(tk0)*cos(phik0),Ek*sin(tk0)*sin(phik0),Ek*cos(tk0),Ek);
            //         q1cm->SetPxPyPzE(q1u(q1r,kcm)*q1r->X(),q1u(q1r,kcm)*q1r->Y(),
            //             q1u(q1r,kcm)*q1r->Z(),eps1u(q1r,kcm));
            //         double val = bremCS(Ek,q1cm,kcm);
            //         // cout<<"val: "<<val<<endl;
            //         intPhoton->setArray(i,j,val);
            //     }
            // }
            double eBinD = (Ek-dE)/(k0-dE)*(double(pRes)-1.);
            double tBinD = (tqr-tqrCut0)/(tqrCut1-tqrCut0)*(double(pRes)-1.);
            // cout<<"eBin_d: "<<eBinD<<endl;
            int eBin = (int) floor(eBinD+0.5);
            int tBin = (int) floor(tBinD+0.5);

            photonCoords = photonArray[eBin*pRes+tBin]->getRand2D();
            // cout<<"coords 0: "<<photonCoords[0]<<endl;
            // cout<<"coords 1: "<<photonCoords[1]<<endl;

            tk = photonCoords[0];
            phik = photonCoords[1]+pqr;
            tkWeight = photonCoords[2]*sin(tk);
            phikWeight = 1.;


            kcm->SetPxPyPzE(Ek*sin(tk)*cos(phik),Ek*sin(tk)*sin(phik),Ek*cos(tk),Ek);
            q1cm->SetPxPyPzE(q1u(q1r,kcm)*q1r->X(),q1u(q1r,kcm)*q1r->Y(),
                q1u(q1r,kcm)*q1r->Z(),eps1u(q1r,kcm));

            // delete intPhoton;

            aFlag = 1.;  //no reweighting

        }

        q2cm = new TLorentzVector(*p1+*p2-*q1cm-*kcm);

        q1 = new TLorentzVector(*q1cm);
        q2 = new TLorentzVector(*q2cm);
        k = new TLorentzVector(*kcm);

        q1->Boost(cm->BoostVector());
        q2->Boost(cm->BoostVector());
        k->Boost(cm->BoostVector());

        if(mb_flag==1){//Moller
            weight = bremCS(Ek,q1cm,kcm)*symWeight(q2cm->Theta(),q2cm->Phi())/radFrac\
            *ekWeight*tkWeight*phikWeight*tqrWeight*pqrWeight*aFlag;
            // cout<<"ek "<<ekWeight<<" tkw "<<tkWeight<<" pkw "<<phikWeight<<" tqrw "<<tqrWeight<<" pqrw "<<pqrWeight<<endl;
        }

        if(mb_flag==0){//Bhabha
            weight = bremCSb(Ek,q1cm,kcm)/radFrac\
            *ekWeight*tkWeight*phikWeight*tqrWeight*pqrWeight*aFlag;
        }

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

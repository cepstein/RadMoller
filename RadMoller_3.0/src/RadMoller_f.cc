//Function implementations for RadMoller_Gen class
#include "RadMoller.h"
#include "Riostream.h"


void RadMoller_Gen::SetMoller(){
    mb_flag = 1;
}

void RadMoller_Gen::SetBhabha(){
    mb_flag = 0;
}

void RadMoller_Gen::SetRadFrac(double rdf){
    radFrac = rdf;
}


double RadMoller_Gen::symWeight(double tq2, double pq2){
    double symFactor;
    if((tq2 < tqrCut1) && (tq2 > tqrCut0) && (pq2 < phiq1) && (pq2 > phiq0)){
        symFactor = 1.;
    }
    else{
        symFactor = 2.;
    }
    return symFactor;
}


void RadMoller_Gen::SetTCuts(double tqrc0, double tqrc1, double pq0, double pq1){
    if(mb_flag == 1){
        cout<<"Moller Events Requested!"<<endl;
    }
    else{
        cout<<"Bhabha Events Requested!"<<endl;
    }

        cout<<"You Requested: "<<endl;
        cout<<"Min Electron_1 Theta = "<<tqrc0<<endl;
        cout<<"Max Electron_1 Theta = "<<tqrc1<<endl;
        cout<<"Min Electron_1 Phi = "<<pq0<<endl;
        cout<<"Max Electron_1 Phi = "<<pq1<<endl;
        
        tkCut0 = 0.;
        tkCut1 = TMath::Pi();
        tqrCut0 = tqrc0;
        tqrCut1 = tqrc1;
        phiq0 = pq0;
        phiq1 = pq1;

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
    return me*(-((Cosa*Ek*sqrt(-pow((2*Ecmp)/me - Ek/me,2) + 
            (pow(Cosa,2)*pow(Ek,2))/pow(me,2) + 
            (4*pow(Ecmp,2)*pow(Ecmp/me - Ek/me,2))/pow(me,2)))/me) + 
     (2*Ecmp*(Ecmp/me - Ek/me)*((2*Ecmp)/me - Ek/me))/me)/
   (pow((2*Ecmp)/me - Ek/me,2) - (pow(Cosa,2)*pow(Ek,2))/pow(me,2));
}

double RadMoller_Gen::eps1l(TVector3 *q1, TLorentzVector *k){ //me-> units of MeV    
    double Cosa = TMath::Cos(k->Angle(*q1));
    return me*((Cosa*Ek*sqrt(-pow((2*Ecmp)/me - Ek/me,2) + 
          (pow(Cosa,2)*pow(Ek,2))/pow(me,2) + 
          (4*pow(Ecmp,2)*pow(Ecmp/me - Ek/me,2))/pow(me,2)))/me + 
     (2*Ecmp*(Ecmp/me - Ek/me)*((2*Ecmp)/me - Ek/me))/me)/
   (pow((2*Ecmp)/me - Ek/me,2) - (pow(Cosa,2)*pow(Ek,2))/pow(me,2));
}

double RadMoller_Gen::q1u(TVector3 *q1, TLorentzVector *k){ //me-> units of MeV
    double Cosa = TMath::Cos(k->Angle(*q1));
    return me*(sqrt(-pow((2*Ecmp)/me - Ek/me,2) + (pow(Cosa,2)*pow(Ek,2))/pow(me,2) + 
        (4*pow(Ecmp,2)*pow(Ecmp/me - Ek/me,2))/pow(me,2))*
      ((2*Ecmp)/me - Ek/me) - (2*Cosa*Ecmp*Ek*(Ecmp/me - Ek/me))/pow(me,2))/
   (pow((2*Ecmp)/me - Ek/me,2) - (pow(Cosa,2)*pow(Ek,2))/pow(me,2));
}

double RadMoller_Gen::q1l(TVector3 *q1, TLorentzVector *k){ //me-> units of MeV
    double Cosa = TMath::Cos(k->Angle(*q1));
    return me*(sqrt(-pow((2*Ecmp)/me - Ek/me,2) + (pow(Cosa,2)*pow(Ek,2))/pow(me,2) + 
        (4*pow(Ecmp,2)*pow(Ecmp/me - Ek/me,2))/pow(me,2))*
      ((-2*Ecmp)/me + Ek/me) - (2*Cosa*Ecmp*Ek*(Ecmp/me - Ek/me))/pow(me,2))/
   (pow((2*Ecmp)/me - Ek/me,2) - (pow(Cosa,2)*pow(Ek,2))/pow(me,2));
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
    return (alpha*(-46*pow(se,2)*pow(te(x) + ue(x),2) - 46*(pow(te(x),4) + pow(ue(x),4)) + 
       33*(pow(ue(x),4) + pow(se,2)*ue(x)*(te(x) + ue(x)))*log(-(te(x)/pow(me,2))) + 
       33*te(x)*(pow(te(x),3) + pow(se,2)*(te(x) + ue(x)))*log(-(ue(x)/pow(me,2))) - 
       36*(pow(te(x),4) + pow(ue(x),4) + pow(se,2)*pow(te(x) + ue(x),2))*log(Ecmp/dE)*
        (-1 + log((te(x)*ue(x))/(pow(me,2)*se)))))/
   (9.*pi*(pow(te(x),4) + pow(ue(x),4) + pow(se,2)*pow(te(x) + ue(x),2)));
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
        //Exponentiating: CS_soft = CS_tree x exp(delta), with delta = soft cs correction / tree cs
        //Instead of 1 + delta, it is exp(delta), to account for multiple soft photons
        //Includes factor of 2*pi*sin(theta) to make it d_sigma/d_theta
        return Lumi*hbarc2*tree_cs(x)*exp(soft_cs_tsai(x,dE));
    }

double RadMoller_Gen::bCSfunc(double x, double dE)
    {
        //Extra factor of two because e+,e- are distinguishable 
        return 2.0*Lumi*hbarc2*tree_cs_b(x)*exp(soft_bhabha(x,dE));
    }


//Construct the e-e- Bremsstrahlung Cross Section 
double RadMoller_Gen::bremCS(double Ek, 
  TLorentzVector *q1, TLorentzVector *k)
    {
        double bCS;
        double Cosa = cos(k->Angle(q1->Vect()));
        bCS = (Ek*me)/
        (2.*Ecmp*sqrt(-pow((2.*Ecmp)/me - Ek/me,2.) + 
       (pow(Cosa,2.)*pow(Ek,2.))/pow(me,2.) + 
       (4.*pow(Ecmp,2.)*pow(Ecmp/me - Ek/me,2.))/pow(me,2.))*Pcmp)*pow(q1->P(),2.);
        
        return 0.5*bCS*Lumi*Mh2(p1,p2,q1,k)*pow(2*pi,-5)*pow(me,-3)*hbarc2/32;
    }

//Construct the e-e+ Bremsstrahlung Cross Section 
double RadMoller_Gen::bremCSb(double Ek, 
  TLorentzVector *q1, TLorentzVector *k)
    {
        double bCS;
        double Cosa = cos(k->Angle(q1->Vect()));
        bCS = (Ek*me)/
        (2.*Ecmp*sqrt(-pow((2.*Ecmp)/me - Ek/me,2.) + 
       (pow(Cosa,2.)*pow(Ek,2.))/pow(me,2.) + 
       (4.*pow(Ecmp,2.)*pow(Ecmp/me - Ek/me,2.))/pow(me,2.))*Pcmp)*pow(q1->P(),2.);
        
        return bCS*Lumi*Mh2b(p1,p2,q1,k)*pow(2*pi,-5)*pow(me,-3)*hbarc2/32;
    }


void RadMoller_Gen::InitGenerator_RadMoller(){
    me = 0.510998910;
    Ebeam = Tbeam+me;
    alpha = 1./137.035999074;
    hbarc2 = pow(1.97326e-11,2.);
    pi = 4.0*atan(1.0);
    twopi = 2.*pi;
    Pbeam = sqrt(pow(Ebeam,2.)-pow(me,2.));
    betacm = Pbeam/(Ebeam+me);
    gammacm = 1./sqrt(1.-pow(betacm,2.));
    Ecm = gammacm*Ebeam - gammacm*betacm*Pbeam + gammacm*me;
    Pcm = sqrt(pow(Ecm/2.,2)-pow(me,2.));//momentum of either 
    Ecmp = Ecm/2.; //Ecm per particle
    Pcmp = sqrt(pow(Ecmp,2)-pow(me,2.));//momentum of either b
    ec = sqrt(4.*pi*alpha); //electron charge
    se = 4.*Ecmp*Ecmp; //Elastic Mandelstam S ("s" was unavail
    dE = dE_frac*sqrt(se);
    EkMax = (Ecm*Ecm-4.*me*me)/(2.*Ecm);
    k0 = me*(2*Ecmp*(-1 + Ecmp/me))/((-1 + (2*Ecmp)/me)*me); //units of MeV

    cm = new TLorentzVector(0.,0.,Pbeam,Ebeam+me);
    p1 = new TLorentzVector(0,0,Pbeam,Ebeam);
    p2 = new TLorentzVector(0,0,0,me);
    p1->Boost(-cm->BoostVector());
    p2->Boost(-cm->BoostVector());

}


double RadMoller_Gen::tqrFunc(double tqr){
    return 4/(3.*pi*(1 - pow(-1 + (2*tqr)/(3.*pi),2))*
     (log(3*pi - tqrCut0) - log(tqrCut0) - log(3*pi - tqrCut1) + log(tqrCut1)));
}

double RadMoller_Gen::tqrFunc_Moller(double tqr){
    return 4/(pi*(1 - pow(-1 + (2*tqr)/pi,2))*
     (log(pi - tqrCut0) - log(tqrCut0) - log(pi - tqrCut1) + log(tqrCut1)));
}

double RadMoller_Gen::tkFunc(double tk){
    return 4/(pi*(1 - pow(-1 + (2*tk)/pi,2))*
     (log(pi - tkCut0) - log(tkCut0) - log(pi - tkCut1) + log(tkCut1)));
}

double RadMoller_Gen::tkiCDF(double rand){
    return (pi*pow(pi - tkCut0,rand)*tkCut0*pow(tkCut1,rand))/
   (pi*pow(tkCut0,rand)*pow(pi - tkCut1,rand) - 
     pow(tkCut0,1 + rand)*pow(pi - tkCut1,rand) + 
     pow(pi - tkCut0,rand)*tkCut0*pow(tkCut1,rand));
}

double RadMoller_Gen::tqriCDF(double rand){
    return (3*pi*pow(3*pi - tqrCut0,rand)*tqrCut0*pow(tqrCut1,rand))/
   (3*pi*pow(tqrCut0,rand)*pow(3*pi - tqrCut1,rand) - 
     pow(tqrCut0,1 + rand)*pow(3*pi - tqrCut1,rand) + 
     pow(3*pi - tqrCut0,rand)*tqrCut0*pow(tqrCut1,rand));
}

double RadMoller_Gen::tqriCDF_Moller(double rand){
    return (pi*pow(pi - tqrCut0,rand)*tqrCut0*pow(tqrCut1,rand))/
   (pi*pow(tqrCut0,rand)*pow(pi - tqrCut1,rand) - 
     pow(tqrCut0,1 + rand)*pow(pi - tqrCut1,rand) + 
     pow(pi - tqrCut0,rand)*tqrCut0*pow(tqrCut1,rand));
}

double RadMoller_Gen::ekFunc(double ek){
    return 2./(exp((2*ek)/EkMax)*(-exp(-2) + exp((-2*dE)/EkMax))*EkMax);
}

double RadMoller_Gen::ekiCDF(double rand){
    return (EkMax*log(exp(2 + (2*dE)/EkMax)/
       (exp(2) - exp(2)*rand + exp((2*dE)/EkMax)*rand)))/2.;
}

void RadMoller_Gen::Generate_Event(){
    
    if(mb_flag==1){//Moller
        tqr  = tqriCDF_Moller(randomGen());
        tqrWeight  = sin(tqr)/tqrFunc_Moller(tqr);
    }
    else{//Bhabha
        tqr  = tqriCDF(randomGen());
        tqrWeight  = sin(tqr)/tqrFunc(tqr);
    }

    pqr = phiq0+(phiq1-phiq0)*randomGen();
    pqrWeight = phiq1-phiq0;

    pickProc = randomGen();


    if (pickProc<radFrac){//Bremsstrahlung

        Ek = dE+(EkMax-dE)*randomGen();//ekiCDF(randomGen());//
        ekWeight = EkMax-dE;//1./ekFunc(Ek);//

        delete qcm;
        delete q1r;
        delete q1cm;
        delete q2cm;
        delete q1;
        delete q2;
        delete k;
        delete kcm;

        double aRand = randomGen();
        // cout<<"test point"<<endl;
        q1r = new TVector3(sin(tqr)*cos(pqr),sin(tqr)*sin(pqr),cos(tqr));

        if (Ek>k0){
            cosaMax = -((sqrt(pow((2*Ecmp)/me - Ek/me,2) - 
                (4*pow(Ecmp,2)*pow(Ecmp/me - Ek/me,2))/pow(me,2))*me)/Ek);
            
            delete newAxis;
            newAxis = new TVector3(-sin(tqr)*cos(pqr),-sin(tqr)*sin(pqr),-cos(tqr));

            tk  = (pi-TMath::ACos(cosaMax))*randomGen();

            phik = twopi*randomGen();

            kcm = new TLorentzVector(Ek*sin(tk)*cos(phik),Ek*sin(tk)*sin(phik),Ek*cos(tk),Ek);
            kcm->RotateUz(*newAxis);

            tkWeight  = (1.+cosaMax)*sin(tk);
            phikWeight = twopi;

            if (aRand>0.5){
                q1cm = new TLorentzVector(q1u(q1r,kcm)*q1r->X(),q1u(q1r,kcm)*q1r->Y(),
                    q1u(q1r,kcm)*q1r->Z(),eps1u(q1r,kcm));   
            }
            else if (aRand <0.5){
                q1cm = new TLorentzVector(q1l(q1r,kcm)*q1r->X(),q1l(q1r,kcm)*q1r->Y(),
                    q1l(q1r,kcm)*q1r->Z(),eps1l(q1r,kcm));
            }
            aFlag = 2.; //weight up by 2.

        }

        else if (Ek<k0){
            phik = twopi*randomGen();
            phikWeight = twopi;

            tk   = tkCut0+(tkCut1-tkCut0)*randomGen();
            tkWeight = sin(tk)*(tkCut1-tkCut0);

            kcm = new TLorentzVector(Ek*sin(tk)*cos(phik),Ek*sin(tk)*sin(phik),Ek*cos(tk),Ek);

            q1cm = new TLorentzVector(q1u(q1r,kcm)*q1r->X(),q1u(q1r,kcm)*q1r->Y(),
                q1u(q1r,kcm)*q1r->Z(),eps1u(q1r,kcm));
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
        }

        if(mb_flag==0){//Bhabha
            weight = bremCSb(Ek,q1cm,kcm)/radFrac\
            *ekWeight*tkWeight*phikWeight*tqrWeight*pqrWeight*aFlag;
        }

        elFlag = 1;
    }

    if (pickProc >radFrac){//Elastic Kinematics
        delete q1cm;
        delete q2cm;
        delete q1;
        delete q2;

        if(mb_flag==1){//Moller
        weight = mCSfunc(tqr,dE)*symWeight(q2cm->Theta(),q2cm->Phi())*tqrWeight*pqrWeight/(1.-radFrac);
        }

        if(mb_flag==0){//Bhabha
        weight = bCSfunc(tqr,dE)*tqrWeight*pqrWeight/(1.-radFrac);
        }

        q1cm = new TLorentzVector(Pcmp*sin(tqr)*cos(pqr),\
            Pcmp*sin(tqr)*sin(pqr),Pcmp*cos(tqr),Ecmp);
        q2cm = new TLorentzVector(-Pcmp*sin(tqr)*cos(pqr),\
            -Pcmp*sin(tqr)*sin(pqr),-Pcmp*cos(tqr),Ecmp);

        q1 = new TLorentzVector(*q1cm);
        q2 = new TLorentzVector(*q2cm);

        q1->Boost(cm->BoostVector());
        q2->Boost(cm->BoostVector());
        elFlag = 0;
    }



}

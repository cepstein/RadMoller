//Function implementations for RadMoller_Gen class

double RadMoller_Gen::GetHistRandom(TH1D *hist){
    //Adapted from ROOT GetRandom to use class-implemented random number generator.
    //Return a random number distributed according the histogram bin contents.
    delete[] fIntegral;
    int nbinsx = hist->GetNbinsX();
    double integral = 0;
    fIntegral = new double [nbinsx+2];
    fIntegral = hist->GetIntegral();
    double fEntries = hist->GetEntries();

    double r1 = randomGen();
    int ibin = TMath::BinarySearch(nbinsx,fIntegral,r1);
    double x = hist->GetBinLowEdge(ibin+1);
    if (r1 > fIntegral[ibin]) x +=
      hist->GetBinWidth(ibin+1)*(r1-fIntegral[ibin])/(fIntegral[ibin+1] - fIntegral[ibin]);
    return x;

}

void RadMoller_Gen::SetMoller(){
    mb_flag = 1;
}

void RadMoller_Gen::SetBhabha(){
    mb_flag = 0;
}

void RadMoller_Gen::SetOutputFlags(int rf, int tf){
    root_flag = rf;
    txt_flag = tf;
}

void RadMoller_Gen::SetRadFrac(double rdf){
    radFrac = rdf;
}

void RadMoller_Gen::SetTCuts(double tkc, double tqrc, double xec){
    tkCut = tkc;
    tqrCut = tqrc;
    xeCut = xec;
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

//Squared Tree-Level Matrix element - does not take m->0 limit
double RadMoller_Gen::M2(Double_t x)
    {
    return 64.0*pow(pi,2.)*pow(alpha,2.)*(pow(me,4.)/pow(te(x),2)*
        ((pow(se,2)+pow(ue(x),2))/(2.*pow(me,4))+4.*te(x)/pow(me,2)-4.0)+pow(me,4)/
        pow(ue(x),2)*((pow(se,2)+pow(te(x),2))/(2.0*pow(me,4))+4.0*ue(x)/pow(me,2)-4.)+
        pow(me,4)/(ue(x)*te(x))*(se/pow(me,2)-2.0)*(se/pow(me,2)-6.0));
    }

//Squared Tree-Level Bhabha Matrix element
double RadMoller_Gen::M2b(Double_t x)
    {
    return 64.0*pow(pi,2.)*pow(alpha,2.)*(pow(me,4.)/pow(te(x),2)*
        ((pow(ue(x),2)+pow(se,2))/(2.*pow(me,4))+4.*te(x)/pow(me,2)-4.0)+pow(me,4)/
        pow(se,2)*((pow(ue(x),2)+pow(te(x),2))/(2.0*pow(me,4))+4.0*se/pow(me,2)-4.)+
        pow(me,4)/(se*te(x))*(ue(x)/pow(me,2)-2.0)*(ue(x)/pow(me,2)-6.0));
    }

//Tsai's Soft Corrections - Moller - expressed as Delta(dsigma)
double RadMoller_Gen::soft_cs_tsai(double x, double dE){
    return -(alpha*pow(me,2)*pow(alpha/me,2)*(46*
         (pow(te(x),4) + pow(ue(x),4) + pow(se,2)*pow(te(x) + ue(x),2)) - 
        33*ue(x)*(pow(ue(x),3) + pow(se,2)*(te(x) + ue(x)))*log(-(te(x)/pow(me,2))) - 
        33*te(x)*(pow(te(x),3) + pow(se,2)*(te(x) + ue(x)))*log(-(ue(x)/pow(me,2))) + 
        36*(pow(te(x),4) + pow(ue(x),4) + pow(se,2)*pow(te(x) + ue(x),2))*
         log(Ecmp/dE)*(-1 + log((te(x)*ue(x))/(pow(me,2)*se)))))/
   (72.*pow(Ecmp,2)*pi*pow(te(x),2)*pow(ue(x),2));
}

//Soft corrections to Bhabha scattering (Glover, Tausk, Bij) - expresed as dsigma/dsigma0
double RadMoller_Gen::soft_bhabha(double x, double dE){
    return (pow(alpha,2)*((4.5 + 12*log(dE/Ecmp) + 8*pow(log(dE/Ecmp),2))*
        pow(log(se/pow(me,2)),2) + 
       log(se/pow(me,2))*(-11.625 - (5*pow(pi,2))/2. + 
          pow(log(dE/Ecmp),2)*(-16 + 16*log(-(te(x)/(se*(1 + te(x)/se))))) + 
          (3*((pow(pi,2)*(4 + (8*te(x))/se + (27*pow(te(x),2))/pow(se,2) + 
                    (26*pow(te(x),3))/pow(se,3) + (16*pow(te(x),4))/pow(se,4))
                  )/12. + ((6 + (8*te(x))/se + (9*pow(te(x),2))/pow(se,2) + 
                    (3*pow(te(x),3))/pow(se,3))*log(-(te(x)/se)))/2. - 
               (te(x)*(3 + te(x)/se - (3*pow(te(x),2))/pow(se,2) - 
                    (4*pow(te(x),3))/pow(se,3))*pow(log(-(te(x)/se)),2))/(4.*se)
                 + (te(x)*(1 + pow(te(x),2)/pow(se,2))*log(1 + te(x)/se))/(2.*se) + 
               ((4 + (8*te(x))/se + (7*pow(te(x),2))/pow(se,2) + 
                    (2*pow(te(x),3))/pow(se,3))*log(-(te(x)/se))*log(1 + te(x)/se))/
                2. + ((-2 - (5*te(x))/se - (7*pow(te(x),2))/pow(se,2) - 
                    (5*pow(te(x),3))/pow(se,3) - (2*pow(te(x),4))/pow(se,4))*
                  pow(log(1 + te(x)/se),2))/2.))/
           pow(1 + te(x)/se + pow(te(x),2)/pow(se,2),2) - 
          6*TMath::DiLog(-(te(x)/se)) + 6*TMath::DiLog(1 + te(x)/se) + 
          log(dE/Ecmp)*(-28 - (8*pow(pi,2))/3. + 
             12*log(-(te(x)/(se*(1 + te(x)/se)))) + 
             (4*((pow(pi,2)*(4 + (8*te(x))/se + (27*pow(te(x),2))/pow(se,2) + 
                       (26*pow(te(x),3))/pow(se,3) + 
                       (16*pow(te(x),4))/pow(se,4)))/12. + 
                  ((6 + (8*te(x))/se + (9*pow(te(x),2))/pow(se,2) + 
                       (3*pow(te(x),3))/pow(se,3))*log(-(te(x)/se)))/2. - 
                  (te(x)*(3 + te(x)/se - (3*pow(te(x),2))/pow(se,2) - 
                       (4*pow(te(x),3))/pow(se,3))*pow(log(-(te(x)/se)),2))/
                   (4.*se) + (te(x)*(1 + pow(te(x),2)/pow(se,2))*log(1 + te(x)/se))/
                   (2.*se) + ((4 + (8*te(x))/se + (7*pow(te(x),2))/pow(se,2) + 
                       (2*pow(te(x),3))/pow(se,3))*log(-(te(x)/se))*log(1 + te(x)/se)
                     )/2. + ((-2 - (5*te(x))/se - (7*pow(te(x),2))/pow(se,2) - 
                       (5*pow(te(x),3))/pow(se,3) - (2*pow(te(x),4))/pow(se,4)
                       )*pow(log(1 + te(x)/se),2))/2.))/
              pow(1 + te(x)/se + pow(te(x),2)/pow(se,2),2) - 
             8*TMath::DiLog(-(te(x)/se)) + 8*TMath::DiLog(1 + te(x)/se)) + 6*1.20206)))/
   pow(pi,2);
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
        return 2.0*pi*sin(x)*Lumi*pow(1.97e-11,2)*tree_cs(x)*exp(soft_cs_tsai(x,dE)/tree_cs(x));
    }

double RadMoller_Gen::bCSfunc(double x, double dE)
    {
        //Extra factor of two because e+,e- are distinguishable 
        return 2.0*2.0*pi*sin(x)*Lumi*pow(1.97e-11,2)*tree_cs_b(x)*exp(soft_bhabha(x,dE)-1);
    }


//Construct the e-e- Bremsstrahlung Cross Section 
double RadMoller_Gen::dSigmahdEkdTkdTqr(double Ek , double tk, double tqr, double pqr, 
  TLorentzVector *q1, TLorentzVector *q2, TLorentzVector *k)
    {
    double bremCS = pow(1.97e-11,2.)*Lumi*sin(tk)*sin(tqr)*Ek*Mh2(p1,p2,q1,q2,k)/(64.*se*pow(2.*pi,4.));
    return bremCS;
    }

//Construct the e-e+ Bremsstrahlung Cross Section 
double RadMoller_Gen::dSigmahdEkdTkdTqr_b(double Ek , double tk, double tqr, double pqr, 
  TLorentzVector *q1, TLorentzVector *q2, TLorentzVector *k)
    {
    //Extra factor of two because e+,e- are distinguishable 
    double bremCSb = 2.0*pow(1.97e-11,2.)*Lumi*sin(tk)*sin(tqr)*Ek*Mh2b(p1,p2,q1,q2,k)/(64.*se*pow(2.*pi,4.));
    return bremCSb;
    }


void RadMoller_Gen::InitGenerator_RadMoller(){
    me = 0.510998910;
    Ebeam = Tbeam+me;
    alpha = 1./137.035999074;
    pi = 4.0*atan(1.0);
    twopi = 2.*pi;
    Pbeam = sqrt(pow(Ebeam,2.)-pow(me,2.));
    betacm = Pbeam/(Ebeam+me);
    gammacm = 1./sqrt(1.-pow(betacm,2.));
    Ecm = gammacm*Ebeam - gammacm*betacm*Pbeam + gammacm*me;
    Pcm = sqrt(pow(Ecm/2.,2)-pow(me,2.));//momentum of either 
    Ecmp = Ecm/2.; //Ecm per particle
    Pcmp = sqrt(pow(Ecmp,2)-pow(me,2.));//momentum of either b
    ec= sqrt(4.*pi*alpha); //electron charge
    se = 4.*Ecmp*Ecmp; //Elastic Mandelstam S ("s" was unavail
    dE = dE_frac*sqrt(se);
    EkMax = (Ecm*Ecm-4.*me*me)/(2.*Ecm);

    ekDist = new TH1D("ekDist","Ek Distribution",1000,dE,EkMax);
    tkDist = new TH1D("tkDist","Theta K Distribution",1000,tkCut,pi-tkCut);
    tqrDist = new TH1D("tqrDist","Theta K Distribution",1000,tqrCut,pi-tqrCut);

    for (int i = 0; i<ekDist->GetNbinsX(); i++)
        {
        ekDist->SetBinContent(i,1./(i+100.));
        }
    ekDist->Scale(1./ekDist->Integral("width"));

    for (int i = 0; i<tkDist->GetNbinsX(); i++)
        {
        tkDist->SetBinContent(i,(sqrt(pow(pi/2.,2.)-pow((i/1000.)*pi-pi/2.,2.))));
        }
    tkDist->Scale(1./tkDist->Integral("width"));

    for (int i = 0; i<tqrDist->GetNbinsX(); i++)
        {
        tqrDist->SetBinContent(i,1./(sqrt(pow(pi/2.,2.)-pow((i/1000.)*pi-pi/2.,2.))));
        }
    tqrDist->Scale(1./tqrDist->Integral("width"));


    cm = new TLorentzVector(0.,0.,Pbeam,Ebeam+me);
    p1 = new TLorentzVector(0,0,Pbeam,Ebeam);
    p2 = new TLorentzVector(0,0,0,me);


}

void RadMoller_Gen::Generate_Event(){
    pqr = twopi*randomGen();
    Ek   = GetHistRandom(ekDist);//dE+(EkMax-dE)*randomGen();
    tk   = tkCut+(pi-2.*tkCut)*randomGen();//GetHistRandom(tkDist);
    tqr  = GetHistRandom(tqrDist);//tqrCut+(pi-2.*tqrCut)*randomGen();
    ekWeight   =1./ekDist->GetBinContent(ekDist->FindBin(Ek));//EkMax-dE;
    tkWeight   = pi-2.*tkCut;//1./tkDist->GetBinContent(tkDist->FindBin(tk));
    tqrWeight  = 1./tqrDist->GetBinContent(tqrDist->FindBin(tqr));//pi-2.*tqrCut;
    pqrWeight = twopi;
    phik = twopi*randomGen();
    xe = xeCut+(pi-xeCut)*randomGen(); //Elastic Angle

    pickProc = randomGen();


    if (pickProc<radFrac){//Bremsstrahlung

        delete qcm;
        delete qr;
        delete q1r;
        delete q2r;
        delete q1cm;
        delete q2cm;
        delete q1;
        delete q2;
        delete k;


        qcm =  new TLorentzVector(-Ek*sin(tk)*cos(phik),-Ek*sin(tk)*sin(phik),-Ek*cos(tk),Ecm-Ek);
        qr = new TLorentzVector(*qcm);
        qr->Boost(-qcm->BoostVector());

        pqcm = sqrt(pow(qr->E()/2.,2.)-me*me);

        q1r = new TLorentzVector(pqcm*sin(tqr)*cos(pqr),\
            pqcm*sin(tqr)*sin(pqr),pqcm*cos(tqr),qr->E()/2.);
        q2r = new TLorentzVector(-pqcm*sin(tqr)*cos(pqr),\
            -pqcm*sin(tqr)*sin(pqr),-pqcm*cos(tqr),qr->E()/2.);

        q1cm = new TLorentzVector(*q1r);
        q2cm = new TLorentzVector(*q2r);
        q1cm->Boost(qcm->BoostVector());
        q2cm->Boost(qcm->BoostVector());


        q1 = new TLorentzVector(*q1cm);
        q2 = new TLorentzVector(*q2cm);
        q1->Boost(cm->BoostVector());
        q2->Boost(cm->BoostVector());
        k = new TLorentzVector(*cm-*q1-*q2);

        if(mb_flag==1){//Moller
            weight = dSigmahdEkdTkdTqr(Ek,tk,tqr,pqr,q1,q2,k)/radFrac\
            *ekWeight*tkWeight*tqrWeight*pqrWeight;
        }
        if(mb_flag==0){//Bhabha
            weight = dSigmahdEkdTkdTqr_b(Ek,tk,tqr,pqr,q1,q2,k)/radFrac\
            *ekWeight*tkWeight*tqrWeight*pqrWeight;
        }

        elFlag = 1;
    }

    if (pickProc >radFrac){//Elastic Kinematics
        delete q1cm;
        delete q2cm;
        delete q1;
        delete q2;

        if(mb_flag==1){//Moller
        weight = mCSfunc(xe,dE)*(pi-2.*xeCut)/(1.-radFrac);
        }

        if(mb_flag==0){//Bhabha
        weight = bCSfunc(xe,dE)*(pi-2.*xeCut)/(1.-radFrac);
        }

        q1cm = new TLorentzVector(Pcmp*sin(xe)*cos(phik),\
            Pcmp*sin(xe)*sin(phik),Pcmp*cos(xe),Ecmp);
        q2cm = new TLorentzVector(-Pcmp*sin(xe)*cos(phik),\
            -Pcmp*sin(xe)*sin(phik),-Pcmp*cos(xe),Ecmp);

        q1 = new TLorentzVector(*q1cm);
        q2 = new TLorentzVector(*q2cm);

        q1->Boost(cm->BoostVector());
        q2->Boost(cm->BoostVector());
        elFlag = 0;
    }



}

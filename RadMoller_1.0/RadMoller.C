//================================================================================
//================================================================================
//
//      Radiative Moller Generator
//      Charles Epstein, MIT, Fall 2013
//      Soft Radiative Corrections based on Denner & Pozzorini, 1999 
//            or Tsai, 1960 (User Choice)
//      Hard Brehmsstralung calculated using FeynArts & FormCalc
//            for unpolarized ee>eey scattering
//
//================================================================================
//================================================================================

#include "vars.h"
#include "./include/RadMoller.h"
#include "./include/MSqBrem.h"
#include "./include/DenPozSoft.h"
#include "./include/Tsai.h"

//Squared Tree-Level Matrix element - does not take m->0 limit
double M2(Double_t x)
    {
    return 64.0*pow(pi,2.)*pow(alpha,2.)*(pow(me,4.)/pow(te(x),2)*
        ((pow(se,2)+pow(ue(x),2))/(2.*pow(me,4))+4.*te(x)/pow(me,2)-4.0)+pow(me,4)/
        pow(ue(x),2)*((pow(se,2)+pow(te(x),2))/(2.0*pow(me,4))+4.0*ue(x)/pow(me,2)-4.)+
        pow(me,4)/(ue(x)*te(x))*(se/pow(me,2)-2.0)*(se/pow(me,2)-6.0));
    }


//Construct the Tree-Level Cross Section (CMS)
double tree_cs(double x)
    {
    return 0.5*pow(8.0*pi,-2)*M2(x)*pow(Ecm,-2);
    }

//Construct the Soft-Brehmsstralung-Corrected Cross-Section (CMS)
double mCSfunc(double x, double dE)
    {
    double mCS;
    if (soft_flag == 0){
        mCS = 2.0*pi*sin(x)*Lumi*pow(1.97e-11,2)*(tree_cs(x)+corr_soft(x,dE));
      }
    else if (soft_flag == 1){
            //Exponentiating: CS_soft = CS_tree x exp(delta), with delta = soft correction / CS_tree
            //Instead of 1 + delta, it is exp(delta), to account for multiple soft photons
            //Includes factor of 2*pi*sin(theta) to make it d_sigma/d_theta
        mCS = 2.0*pi*sin(x)*Lumi*pow(1.97e-11,2)*tree_cs(x)*exp(corr_soft_tsai(x,dE));
    }
    else {
      cout<<"Error: soft_flag not properly set."<<endl;
      mCS = 0;
    }
    return mCS;
    }


//Construct the e-e Bremsstrahlung Cross Section 
double dSigmahdEkdTkdTqr(double Ek , double tk, double tqr, double pqr, 
  TLorentzVector *q1, TLorentzVector *q2, TLorentzVector *k)
    {
    TLorentzVector *p1 = new TLorentzVector(0,0,Pbeam,Ebeam);
    TLorentzVector *p2 = new TLorentzVector(0,0,0,me);
    double bremCS = twopi*pow(1.97e-11,2.)*Lumi*sin(tk)*sin(tqr)*Ek*Mh2(p1,p2,q1,q2,k)/(64.*se*pow(2.*pi,4.));
    return bremCS;
    }


int main(int argc, char* argv[])
    {
        if (argc<2){
            cout<<"Usage: ./RadMoller [nEve]"<<endl;
            return 1;
        }
        nEve = atoi(argv[1]); 

    //=========================================================
    //=========================================================
    TFile *f = new TFile("histos.root","recreate");
    TFile *fout = new TFile("radEve.root","recreate");
    FILE *eFile;
    eFile = fopen("radEve.txt","w");

    if(txt_flag){
        fprintf(eFile,"Weight\tP1x\tP1y\tP1z\tP2z\tP2y\tP2z\tKx\tKy\tKz\n");
    }

    TH2D  *hst_xy = new TH2D("hst_xy" ,  "Scattered Electron Energy vs. Angle", 90,0,90,501,0.5,Ebeam);
    TH2D *ptpz = new TH2D("ptpz", "Electron Pt vs Pz",104,0,Ebeam,50,0,Ebeam/10);
    TH2D *ptpzg = new TH2D("ptpzg", "Photon Pt vs Pz",104,0,Ebeam,50,0,Ebeam/10);
    TH2D *aa = new TH2D("aa","e1 Angle vs e2 Angle",180,0,90,180,0,90);
    TH2D *pp = new TH2D("pp","e1 Momentum vs e2 Momentum",199,0.5,Ebeam,199,0.5,Ebeam);
    TH2D *eg = new TH2D("eg","e1 Angle vs Photon Angle",360,0,180,180,0,90);
    TH1D *kk = new TH1D("kk","Ek distribution",50,0,Ebeam);
    TH1D *CSint = new TH1D("CSint","Photon ",100,0.004,0.008);
    TH1D *w_int = new TH1D("w_int","Electron Angular Dist",89,0.5,89.5);
    TH1D *ph_angles = new TH1D("ph_angles","Photon Angular Dist",89,0.5,89.5);
    TH2D *ph_erg = new TH2D("ph_erg","Photon Energy Angular Dist",89,0.5,89.5,101,0,Ebeam);
        
    TNtuple *radEve_t = new TNtuple("radEve_t", "Radiative Events", "p1x:p1y:p1z:p2x:p2y:p2z:k1x:k1y:k1z:weight");
    TNtuple *softEve_t = new TNtuple("softEve_t", "Soft Corrected Events", "p1x:p1y:p1z:p2x:p2y:p2z:weight");

        
    hst_xy->Sumw2();
    ptpz->Sumw2();
    ptpzg->Sumw2();
    aa->Sumw2();
    pp->Sumw2();
    eg->Sumw2();
    kk->Sumw2();
    CSint->Sumw2();
    w_int->Sumw2();
    ph_angles->Sumw2();
    ph_erg->Sumw2();
    

    hst_xy->SetXTitle("Scattered Electron Angle");
    hst_xy->SetYTitle("Scattered Electron Energy");

    ptpz->SetXTitle("Pz");
    ptpz->SetYTitle("Pt");
    ptpzg->SetXTitle("Pz");
    ptpzg->SetYTitle("Pt");

    aa->SetXTitle("e1 Angle");
    aa->SetYTitle("e2 Angle");

    pp->SetXTitle("e1 Momentum");
    pp->SetYTitle("e2 Momentum");

    eg->SetYTitle("e1 Angle");
    eg->SetXTitle("Photon Angle");

    kk->SetXTitle("Photon Energy [MeV]");


    TH1D *ekDist = new TH1D("ekDist","Ek Distribution",1000,dE,EkMax);
    for (int i = 0; i<ekDist->GetNbinsX(); i++)
        {
        ekDist->SetBinContent(i,1./(i+100.));
        }
    ekDist->Scale(1./ekDist->Integral("width"));


    TH1D *tkDist = new TH1D("tkDist","Theta K Distribution",1000,tkCut,pi-tkCut);
    TH1D *tqrDist = new TH1D("tqrDist","Theta K Distribution",1000,tqrCut,pi-tqrCut);

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


    TRandom *randGen   = new TRandom3();  // Create random number generator
    randGen->SetSeed(0);
    //=========================================================
    //======================================================================
    cout<<"Generating Events: ";
    int nshow = nEve/10;
    for(long loop=0; loop<nEve; loop++)
        {
        if (loop == nshow){
            nshow += nEve/10;
            cout<<"[][][][]";
        }
        double pqr = twopi*randGen->Uniform();
        double weight;

        double Ek   = ekDist->GetRandom();//dE+(EkMax-dE)*randGen->Uniform();//Ek=//Must be less than 5.016  //
        double tk   = tkCut+(pi-2.*tkCut)*randGen->Uniform();//tkDist->GetRandom();//tkCut+(pi-tkCut)*PseRan->Uniform();
        double tqr  = tqrDist->GetRandom();//tqrCut+(pi-2.*tqrCut)*randGen->Uniform();//tqrCut+(pi-tqrCut)*PseRan->Uniform();

        double ekWeight   = 1./ekDist->GetBinContent(ekDist->FindBin(Ek));//EkMax-dE;//// replaces range of Ek
        double tkWeight   = pi-2.*tkCut;//1./tkDist->GetBinContent(tkDist->FindBin(tk));//replaces range of tk
        double tqrWeight  = 1./tqrDist->GetBinContent(tqrDist->FindBin(tqr));//pi-2.*tqrCut;////replaces range of tk

        double pqrWeight = twopi;

        double phik = twopi*randGen->Uniform();
        double xe = xeCut+(pi-xeCut)*randGen->Uniform(); //Elastic Angle
        // double tmin = 0.;//0.001;
        // double tmax = pi;//pi/2.-tmin;
        // double Ecut = 0.;//0.001;
        double pickProc = randGen->Uniform();


        TLorentzVector *cm = new TLorentzVector(0.,0.,Pbeam,Ebeam+me);

        
        if (pickProc<0.5)
            {
            TLorentzVector *qcm = new TLorentzVector(-Ek*sin(tk)*cos(phik),-Ek*sin(tk)*sin(phik),-Ek*cos(tk),Ecm-Ek);
            TLorentzVector *qr = new TLorentzVector(*qcm);
            qr->Boost(-qcm->BoostVector());

            double pqcm = sqrt(pow(qr->E()/2.,2.)-me*me);

            TLorentzVector *q1r = new TLorentzVector(pqcm*sin(tqr)*cos(pqr),\
                pqcm*sin(tqr)*sin(pqr),pqcm*cos(tqr),qr->E()/2.);
            TLorentzVector *q2r = new TLorentzVector(-pqcm*sin(tqr)*cos(pqr),\
                -pqcm*sin(tqr)*sin(pqr),-pqcm*cos(tqr),qr->E()/2.);

            TLorentzVector *q1cm = new TLorentzVector(*q1r);
            TLorentzVector *q2cm = new TLorentzVector(*q2r);
            q1cm->Boost(qcm->BoostVector());
            q2cm->Boost(qcm->BoostVector());


            TLorentzVector *q1 = new TLorentzVector(*q1cm);
            TLorentzVector *q2 = new TLorentzVector(*q2cm);
            q1->Boost(cm->BoostVector());
            q2->Boost(cm->BoostVector());
            TLorentzVector *k = new TLorentzVector(*cm-*q1-*q2);

            weight = 2.*dSigmahdEkdTkdTqr(Ek,tk,tqr,pqr,q1,q2,k)/nEve\
            *ekWeight*tkWeight*tqrWeight*pqrWeight;
                
            if (weight<0){cout<<"Error! Negative Weight: "<<weight<<endl;}

            if (root_flag){
                    radEve_t->Fill(q1->Px(),q1->Py(),q1->Pz(),q2->Px(),q2->Py(),
                        q2->Pz(),k->Px(),k->Py(),k->Pz(),weight);
            }
            if (txt_flag){
                   fprintf(eFile,"%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n",
                    weight,q1->Px(),q1->Py(),q1->Pz(),q2->Px(),q2->Py(),q2->Pz(),k->Px(),k->Py(),k->Pz());
            }

            hst_xy->Fill(q1->Theta()*180./pi,q1->E(),weight);
            //hst_xy->Fill(q2->Theta()*180./pi,q2->E(),weight);
            aa->Fill(q1->Theta()*180./pi,q2->Theta()*180./pi,weight);
            pp->Fill(q1->P(),q2->P(),weight);
            ptpz->Fill(q1->Pz(),q1->Pt(),weight);
            //ptpz->Fill(q2->Pz(),q2->Pt(),weight);
            ptpzg->Fill(k->Pz(),k->Pt(),weight);
            //cout<<"q1 E "<<q1->E()<<endl;
            eg->Fill(k->Theta()*180./pi,q1->Theta()*180./pi,weight);
            CSint->Fill(q1->Theta(),weight);
            w_int->Fill(q1->Theta()*180./pi,weight);
            //if (k->E()>1.0){
            kk->Fill(k->E(),weight);
            ph_angles->Fill(k->Theta()*180./pi,weight);
            ph_erg->Fill(k->Theta()*180./pi,k->E(),weight);
                       

            }
            
        if (pickProc >0.5){

            weight = 2.*mCSfunc(xe,dE)*(pi-2.*xeCut)/nEve;
                
            TLorentzVector *q1cm = new TLorentzVector(Pcmp*sin(xe)*cos(phik),\
                Pcmp*sin(xe)*sin(phik),Pcmp*cos(xe),Ecmp);
            TLorentzVector *q2cm = new TLorentzVector(-Pcmp*sin(xe)*cos(phik),\
                -Pcmp*sin(xe)*sin(phik),-Pcmp*cos(xe),Ecmp);

            TLorentzVector *q1 = new TLorentzVector(*q1cm);
            TLorentzVector *q2 = new TLorentzVector(*q2cm);

            q1->Boost(cm->BoostVector());
            q2->Boost(cm->BoostVector());


            //Save To Files
            if (root_flag){
                softEve_t->Fill(q1->Px(),q1->Py(),q1->Pz(),q2->Px(),q2->Py(),q2->Pz(),weight);
            }
            if (txt_flag){
                fprintf(eFile,"%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n",
                    weight,q1->Px(),q1->Py(),q1->Pz(),q2->Px(),q2->Py(),q2->Pz());
            }

            //Plot in Histograms
            hst_xy->Fill(q1->Theta()*180./pi,q1->E(),weight);
            ptpz->Fill(q1->Pz(),q1->Pt(),weight);
            //hst_xy->Fill(q2->Theta()*180./pi,q2->E(),weight);
            //ptpz->Fill(q2->Pz(),q2->Pt(),weight);
            aa->Fill(q1->Theta()*180./pi,q2->Theta()*180./pi,weight);
            pp->Fill(q1->P(),q2->P(),weight);
            CSint->Fill(q1->Theta(),weight);
            w_int->Fill(q1->Theta()*180./pi,weight);

            }


}
    cout<<endl;
    cout<<"Average (CrossSection)*(Luminosity) of Generated Events "<<w_int->GetBinContent(2)/nEve\
    <<" +/- "<<w_int->GetBinError(2)/nEve<<endl;
    
    TH1D *ph_proj;
    for (int i = 0;i<w_int->GetNbinsX();i++){
        cout<<w_int->GetBinCenter(i)<<" degree electron rate (Hz):\t"<<\
        w_int->GetBinContent(i)/nEve<<"\t+/-\t"<<w_int->GetBinError(i)/nEve<<endl;
    }

    for (int i = 0;i<ph_erg->GetNbinsX();i++){
        ph_proj = ph_erg->ProjectionY("ph_proj",i,i);
        cout<<ph_angles->GetBinCenter(i)<<" degree mean photon energy (MeV):\t"<<\
        ph_proj->GetMean()<<"\tStdDev: "<<ph_proj->GetRMS()<<endl;
        cout<<ph_angles->GetBinCenter(i)<<" degree photon rate (Hz):\t"<<\
        ph_angles->GetBinContent(i)/nEve<<"\t+/-\t"<<ph_angles->GetBinError(i)/nEve<<endl;
    }
    f->cd();
    hst_xy->Write();
    ptpz->Write();
    ptpzg->Write();
    pp->Write();
    ph_angles->Write();
    ph_erg->Write();
    kk->Write();
    w_int->Write();
    f->Close();

    if(root_flag){
        fout->cd();
        radEve_t->Write();
        softEve_t->Write();
        radEve_t->Print();
        softEve_t->Print();
        fout->Close();
    }
 
    
    if(txt_flag){
        fclose(eFile);
    }
    cout<<"Finished Generating Events!  View Histograms by executing 'root -l ../plotHistos.C'"<<endl;

    return EXIT_SUCCESS;
    } 


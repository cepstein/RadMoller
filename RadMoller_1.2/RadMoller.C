//================================================================================
//================================================================================
//
//      Radiative Moller Generator
//      Charles Epstein, MIT, Spring 2014
//      Soft Radiative Corrections based on Denner & Pozzorini, 1999 
//            or Tsai, 1960 (User Choice)
//      Hard Brehmsstralung calculated using FeynArts & FormCalc
//            for unpolarized ee>eey scattering
//
//================================================================================
//================================================================================

#include "RadMoller.h"
#include "MSqBrem.h"
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

#include "RadMoller_f.h"



TRandom3 *randGen = new TRandom3(0);
double RadMoller_Gen::randomGen(){
    return randGen->Uniform();
}


int main(int argc, char* argv[])
    {
    if (argc<2){
        cout<<"Usage: ./RadMoller [nEve]"<<endl;
        return 1;
    }
    long nEve = atoi(argv[1]); 

    //Constants that must be set:
    //Output Flags
    int root_flag = 1; //1 = Output TNtuple of event parameters
    int txt_flag = 1;  //1 = Output txt file of event parameters
    double radFrac = 0.75; //Fraction of events that are radiative (purely statistical)

    //Phase-space cuts on outgoing particles
    double tkCut = 0.001;
    double tqrCut = 0.001;
    double xeCut = 0.001;

    double dE_frac = 1.e-3;//CMS photon energy cutoff as fraction of S

    //const double Lumi = 1.e30; //cm^2 s^-1 - gives CS in microbarns
    double Lumi = 1.e30;

    //Beam Kinetic Energy
    double Tbeam = 100.; 
    
    double pi = 4.0*atan(1.0);

    //=========================================================
    //  Initialize Generator
    RadMoller_Gen* rMollerGen = new RadMoller_Gen;
    rMollerGen->SetOutputFlags(root_flag,txt_flag);
    rMollerGen->SetRadFrac(radFrac);
    rMollerGen->SetTCuts(tkCut,tqrCut,xeCut);
    rMollerGen->SetECut(dE_frac);
    rMollerGen->SetLumi(Lumi);
    rMollerGen->SetTBeam(Tbeam);
    rMollerGen->InitGenerator_RadMoller();


    //
    //======================================================================

    //=========================================================
    //=========================================================
    TFile *f = new TFile("histos.root","recreate");
    TFile *fout = new TFile("radEve.root","recreate");
    FILE *eFile;
    FILE *oFile;
    eFile = fopen("radEve.txt","w");
    oFile = fopen("summary.txt","w");
    
    fprintf(oFile,"Generator Summary:\n");
    fprintf(oFile,"Number of Events: %.0f\n",nEve);
    fprintf(oFile,"Beam Energy: %.1f\tDelta E Fraction: %.3f\tLuminosity: %.2f\tRadiative Fraction: %.3f\n",Tbeam,dE_frac,Lumi,radFrac);
    fprintf(oFile,"Photon Theta Cut: %.3f\t Rad Electron Theta Cut: %.3f\tElastic Eletron Theta Cut: %.3f\n",tkCut,tqrCut,xeCut);
    
    if(txt_flag){
        fprintf(eFile,"Weight\tP1x\tP1y\tP1z\tP2z\tP2y\tP2z\tKx\tKy\tKz\n");
    }

    TH2D  *hst_xy = new TH2D("hst_xy" ,  "Scattered Electron Energy vs. Angle", 90,0,90,501,0.5,Tbeam);
    TH2D *ptpz = new TH2D("ptpz", "Electron Pt vs Pz",104,0,Tbeam,50,0,Tbeam/10);
    TH2D *ptpzg = new TH2D("ptpzg", "Photon Pt vs Pz",104,0,Tbeam,50,0,Tbeam/10);
    TH2D *aa = new TH2D("aa","e1 Angle vs e2 Angle",180,0,90,180,0,90);
    TH2D *pp = new TH2D("pp","e1 Momentum vs e2 Momentum",199,0.5,Tbeam,199,0.5,Tbeam);
    TH2D *eg = new TH2D("eg","e1 Angle vs Photon Angle",360,0,180,180,0,90);
    TH1D *kk = new TH1D("kk","Ek distribution",50,0,Tbeam);
    TH1D *CSint = new TH1D("CSint","Photon ",100,0.004,0.008);
    TH1D *w_int = new TH1D("w_int","Electron Angular Dist",89,0.5,89.5);
    TH1D *ph_angles = new TH1D("ph_angles","Photon Angular Dist",89,0.5,89.5);
    TH2D *ph_erg = new TH2D("ph_erg","Photon Energy Angular Dist",89,0.5,89.5,101,0,Tbeam);
        
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
    hst_xy->SetYTitle("Scattered Electron Energy [MeV]");

    ptpz->SetXTitle("Pz [MeV/c]");
    ptpz->SetYTitle("Pt [MeV/c]");
    ptpzg->SetXTitle("Pz [MeV/c]");
    ptpzg->SetYTitle("Pt [MeV/c]");

    aa->SetXTitle("e1 Angle");
    aa->SetYTitle("e2 Angle");

    pp->SetXTitle("e1 Momentum [MeV/c]");
    pp->SetYTitle("e2 Momentum [MeV/c]");

    eg->SetYTitle("e1 Angle");
    eg->SetXTitle("Photon Angle");

    kk->SetXTitle("Photon Energy [MeV]");

    ph_angles->SetXTitle("Angle");

    ph_erg->SetXTitle("Angle");
    ph_erg->SetYTitle("Energy [MeV]");

    w_int->SetXTitle("Angle");
    //=====================================================
    cout<<"Generating Events: ";
    int nshow = nEve/10;
    for(long loop=0; loop<nEve; loop++)
        {
        if (loop == nshow){
            nshow += nEve/10;
            cout<<"[][][][]";
        }
        rMollerGen->Generate_Event();
        int elFlag = rMollerGen->GetElFlag();
        double weight = rMollerGen->GetWeight();

        if(elFlag){    
            if (weight<0){cout<<"Error! Negative Weight: "<<weight<<endl;}
            TLorentzVector *q1 = rMollerGen->Getq1(); 
            TLorentzVector *q2 = rMollerGen->Getq2(); 
            TLorentzVector *k = rMollerGen->Getk(); 

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

        else{
            //Save To Files
            TLorentzVector *q1 = rMollerGen->Getq1(); 
            TLorentzVector *q2 = rMollerGen->Getq2(); 

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
    fprintf(oFile,"Average (CrossSection)*(Luminosity) of Generated Events %.2f +/- %.2f \n",w_int->GetBinContent(2)/nEve\
    ,w_int->GetBinError(2)/nEve);
    printf("Average (CrossSection)*(Luminosity) of Generated Events %.2f +/- %.2f \n",w_int->GetBinContent(2)/nEve\
    ,w_int->GetBinError(2)/nEve);

    
    TH1D *ph_proj;
    for (int i = 0;i<w_int->GetNbinsX();i++){
        fprintf(oFile,"%.0f degree electron rate (Hz):\t %.2f\t+/-\t%.2f\n",\
            w_int->GetBinCenter(i),w_int->GetBinContent(i)/nEve,w_int->GetBinError(i)/nEve);
    }

    for (int i = 0;i<ph_erg->GetNbinsX();i++){
        ph_proj = ph_erg->ProjectionY("ph_proj",i,i);
        fprintf(oFile,"%.0f degree mean photon energy (MeV):\t%.2f\tStdDev: %.2f\n",
            ph_angles->GetBinCenter(i),ph_proj->GetMean(),ph_proj->GetRMS());
        fprintf(oFile,"%.0f degree photon rate (Hz):\t %.2f\t+/-\t%.2f\n",ph_angles->GetBinCenter(i),\
        ph_angles->GetBinContent(i)/nEve,ph_angles->GetBinError(i)/nEve);
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
    fclose(oFile);
    cout<<"Finished Generating Events!  View Histograms by executing 'root -l ../plotHistos.C'"<<endl;

    return EXIT_SUCCESS;
    } 


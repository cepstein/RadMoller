//================================================================================
//================================================================================
//
//      Radiative Moller/Bhabha Generator
//      Charles Epstein, MIT, Spring 2014
//      Soft Radiative Corrections based on Tsai, 1960 (Moller); 
//            (A.B. Arbuzov, E.S. Scherbakova 2006) (Bhabha);
//      Hard Brehmsstralung calculated using FeynArts & FormCalc
//            for unpolarized e-e->e-e-y and e+e->e+e-y scattering
//
//================================================================================
//================================================================================

#include "RadMoller.h"
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



TRandom3 *randGen = new TRandom3(0);
double RadMoller_Gen::randomGen(){
    return randGen->Uniform();
}

double version = 3.0;

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
    double pi = 4.0*atan(1.0);

    //CMS Phase-space cuts on outgoing particles
    double tqrCut0 = 1*pi/180.;
    double tqrCut1 = 179*pi/180.;

    double phiq0 = 0.;
    double phiq1 = 2.*pi;

    double dE_frac = 1.e-3;//CMS photon energy cutoff as fraction of S

    //const double Lumi = 1.e30; //cm^2 s^-1 - gives CS in microbarns
    double Lumi = 6.e35;

    //Beam Kinetic Energy
    double Tbeam = 100.; 


    //=========================================================
    //=========================================================
    //=========================================================
    //  Initialize Generator
    RadMoller_Gen* rMollerGen = new RadMoller_Gen;
    rMollerGen->SetMoller();
    rMollerGen->SetRadFrac(radFrac);
    rMollerGen->SetTCuts(tqrCut0,tqrCut1,phiq0,phiq1);
    rMollerGen->SetECut(dE_frac);
    rMollerGen->SetLumi(Lumi);
    rMollerGen->SetTBeam(Tbeam);
    rMollerGen->InitGenerator_RadMoller();
    //
    //=========================================================
    //=========================================================
    //=========================================================


    TFile *f = new TFile("histos.root","recreate");
    TFile *fout = new TFile("radEve.root","recreate");
    FILE *eFile;
    FILE *oFile;
    eFile = fopen("radEve.txt","w");
    oFile = fopen("summary.txt","w");
    fprintf(oFile,"Radiative Moller/Bhabha Generator Version %.2f\n",version);
    fprintf(oFile,"Generator Summary:\n");
    fprintf(oFile,"Number of Events: %.1f\n",nEve);
    fprintf(oFile,"Beam Energy: %.1f\tDelta E Fraction: %.3f\tLuminosity: %.2f\tRadiative Fraction: %.3f\n",Tbeam,dE_frac,Lumi,radFrac);
    fprintf(oFile,"Lepton 1 Theta Cuts: %.3f\t %.3f\t Lepton 1 Phi Cuts: %.3f\t %.3f\n",tqrCut0,tqrCut1,phiq0,phiq1);
    if(txt_flag){
        fprintf(eFile,"Weight\tP1x\tP1y\tP1z\tP2z\tP2y\tP2z\tKx\tKy\tKz\n");
    }

    TH2D *hst_xy = new TH2D("hst_xy" ,  "Scattered Lepton Energy vs. Angle", 180,0,90,501,0.5,Tbeam);
    TH2D *ptpz = new TH2D("ptpz", "Lepton Pt vs Pz",200,0,Tbeam+2,100,0,Tbeam/10);
    TH2D *ptpzg = new TH2D("ptpzg", "Photon Pt vs Pz",200,0,Tbeam+2,100,0,Tbeam/10);
    TH2D *aa = new TH2D("aa","l1 Angle vs l2 Angle",180,0,90,180,0,90);
    TH2D *pp = new TH2D("pp","l1 Momentum vs l2 Momentum",199,0.5,Tbeam,199,0.5,Tbeam);
    TH2D *eg = new TH2D("eg","l1 Angle vs Photon Angle",360,0,180,180,0,90);
    TH1D *kk = new TH1D("kk","Ek distribution",50,0,Tbeam);
    TH1D *CSint = new TH1D("CSint","Photon ",100,0.004,0.008);

    TH1D *weights = new TH1D("weights","stuff",10000,0,2500000);
    TH1D *w_int = new TH1D("w_int","Cross-Section for Lepton 1",100,0,90);

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
    
    hst_xy->SetXTitle("Scattered Lepton Angle");
    hst_xy->SetYTitle("Scattered Lepton Energy [MeV]");

    ptpz->SetXTitle("Pz [MeV/c]");
    ptpz->SetYTitle("Pt [MeV/c]");
    ptpzg->SetXTitle("Pz [MeV/c]");
    ptpzg->SetYTitle("Pt [MeV/c]");

    aa->SetXTitle("l1 Angle");
    aa->SetYTitle("l2 Angle");

    pp->SetXTitle("l1 Momentum [MeV/c]");
    pp->SetYTitle("l2 Momentum [MeV/c]");

    eg->SetYTitle("l1 Angle");
    eg->SetXTitle("Photon Angle");

    kk->SetXTitle("Photon Energy [MeV]");

    ph_angles->SetXTitle("Angle");

    ph_erg->SetXTitle("Angle");
    ph_erg->SetYTitle("Energy [MeV]");

    w_int->SetXTitle("Theta [deg]");
    //=====================================================
    printf("Generating Events: \n");
    int nshow = nEve/10;
    for(long loop=0; loop<nEve; loop++)
        {

        if (loop == nshow){
            nshow += nEve/10;
            cout<<"[][][]"<<flush;
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

            hst_xy->Fill(q1->Theta()*180./pi,q1->E(),weight/nEve);
            hst_xy->Fill(q2->Theta()*180./pi,q2->E(),weight/nEve);
            aa->Fill(q1->Theta()*180./pi,q2->Theta()*180./pi,weight/nEve);
            pp->Fill(q1->P(),q2->P(),weight/nEve);
            ptpz->Fill(q1->Pz(),q1->Pt(),weight/nEve);
            ptpz->Fill(q2->Pz(),q2->Pt(),weight/nEve);
            ptpzg->Fill(k->Pz(),k->Pt(),weight/nEve);

            w_int->Fill(q1->Theta()*180./pi,weight/nEve);

            eg->Fill(k->Theta()*180./pi,q1->Theta()*180./pi,weight/nEve);
            CSint->Fill(q1->Theta(),weight/nEve);

            kk->Fill(k->E(),weight/nEve);
            ph_angles->Fill(k->Theta()*180./pi,weight/nEve);
            ph_erg->Fill(k->Theta()*180./pi,k->E(),weight/nEve);
            
            }

        else{
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
            hst_xy->Fill(q1->Theta()*180./pi,q1->E(),weight/nEve);
            ptpz->Fill(q1->Pz(),q1->Pt(),weight/nEve);
            hst_xy->Fill(q2->Theta()*180./pi,q2->E(),weight/nEve);
            ptpz->Fill(q2->Pz(),q2->Pt(),weight/nEve);
            aa->Fill(q1->Theta()*180./pi,q2->Theta()*180./pi,weight/nEve);
            pp->Fill(q1->P(),q2->P(),weight/nEve);
            CSint->Fill(q1->Theta(),weight/nEve);
            w_int->Fill(q1->Theta()*180./pi,weight/nEve);

            }
            
        }


    cout<<endl;

    f->cd();
    hst_xy->Write();
    ptpz->Write();
    ptpzg->Write();
    pp->Write();
    ph_angles->Write();
    ph_erg->Write();
    kk->Write();
    w_int->Write();
    weights->Write();
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

//================================================================================
//================================================================================
//================================================================================
//
//      Radiative Moller/Bhabha Generator
//      Charles Epstein, MIT, Spring 2014
//      Soft Radiative Corrections based on Tsai, 1960 (Moller); 
//            (A.B. Arbuzov, E.S. Scherbakova 2006) (Bhabha 1st Order);
//            Glover,Tausk, & Bij, 2001 (Bhabha 2nd order)
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
#include <TObjArray.h>
#include <TObjString.h>
#include <getopt.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <errno.h>
#include <assert.h>
#include "RandT3.h"

void argUsage(char  *argv0 , const char *message=NULL);
void argUsageLong(const char *message=NULL);
void setRandomSeed(char *cSeed);

// TRandom3 *randGen = new TRandom3(0);
// double RadMoller_Gen::randomGen(){
//     return randGen->Uniform();
// }

double version = 3.0;


//---------------------------------------------------------------------
//---------------------------------------------------------------------
//---------------------------------------------------------------------

int main(int argc,char** argv){
  char  *userName    = getenv("USER");

  long    nEve=3;
  int    root_flag = false;
  int    txt_flag = false;
  int   Moller = false;
  int   Bhabha = false;
  char  *thetaCut  = (char *)"lepTH_0.001_3.14059";
  char  *phiCut  = (char *)"lepPHI_0_6.2831853";
  double  Tbeam  = 100.;
  double  radFrac  = 0.75;
  int pRes = 25;
  int pCells = 1000;
  double  Lumi     = 2.4e36;//1.e36 cm^2 s^-1 - gives CS in picobarns
  double dE_frac   = 1.e-3;//CMS photon energy cutoff as fraction of S

  //============================================================
  {//========================== START argv parsing ============== for 100+ lines
      // working variables
  int    optInd   = 0;
  char   optChar  = 0;  

  static struct option optLong[] = {
    { "root_flag"  , 0 , &root_flag    , true },
    { "txt_flag"  , 0 , &txt_flag    , true },
    { "Moller"    , 0 , 0    , 'M' },
    { "Bhabha"    , 0 , 0    , 'B' },
    { "help"       , 0 , 0    , 'h' },
    { "longHelp"   , 0 , 0    , 'H' },
    { "nEve"     , 1 , 0    , 'n' }, // has one param & is equivalent to -n
    { "thetaCut"  , 1 , 0    , 'T' },
    { "phiCut"  , 1 , 0    , 'P' },
    { "Tbeam"  , 1 , 0    , 'E' },
    { "Lumi"  , 1 , 0    , 'L' },
    { "dE_frac"  , 1 , 0    , 'd' },
    { "radFrac"  , 1 , 0    , 'f' },
    { "pRes"  , 1 , 0    , 's' },
    { "pCells"  , 1 , 0    , 'b' },

    { 0, 0, 0, 0}
  };


  char  * argv0 = argv[0];
  
  opterr = 0;
  while((optChar = getopt_long(argc,argv,"rtT:P:n:s:b:E:f:L:d:MBhH",optLong,&optInd)) != EOF) {
    switch(optChar) {
    case  0  :   break;
    case 'r' : root_flag  = true;    break; 
    case 't' : txt_flag  = true;    break; 
    case 'T' : thetaCut  = optarg;  break; 
    case 'P' : phiCut  = optarg;  break; 
    case 'n' : nEve      = atoi(optarg);     break; 
    case 's' : pRes      = atoi(optarg);     break; 
    case 'b' : pCells      = atoi(optarg);     break; 
    case 'E' : Tbeam  = atof(optarg);  break; 
    case 'f' : radFrac  = atof(optarg);  break; 
    case 'L' : Lumi  = atof(optarg);  break; 
    case 'd' : dE_frac  = atof(optarg);  break; 
    case 'M' : Moller = true; break;
    case 'B' : Bhabha = true; break;
    case 'h' : argUsage(argv0); return(0);     break;
    case 'H' : argUsage(argv0);  argUsageLong(argv0);  return(0);       break;
    case '?':   
      if (isprint (optopt))
    fprintf (stderr, "\nUnknown option `-%c'\n\n", optopt);
      else
    fprintf (stderr,"\nUnknown option character `\\x%x'.\n\n",optopt);
    default  : argUsage(argv0,"unknown option");    break;
    };
  }
  
 /* Print any remaining command line arguments (not options).   */
  if (optind < argc)    {
      printf ("\n WARN, non-options ARGV-elements: ");
      while (optind < argc) 
    printf ("%s ", argv[optind++]);
      putchar ('\n');
      return -3;
  }

  printf("\n**** Final paramater choice made by user=%s  *** \n",userName);
  printf("Executing  %s  nEve=%d  \n",argv0,nEve);
  printf("physics:  thetaCut=%s  phiCut=%s  \n", thetaCut,phiCut);
  printf("pRes=%d, pCells=%d \n",pRes,pCells);
  printf("output: TNtuple=%d  txt=%d    \n",root_flag, txt_flag);

}

    double tqrCut0;
    double tqrCut1;

    double phiq0;
    double phiq1;

    TString *parList = new TString(thetaCut);
    TObjArray * objA=parList->Tokenize("_");
    TIter*  iter = new TIter(objA);  
    TObject* obj = 0;
    int k=0;
    while( (obj=iter->Next())) {    
      k++;
      TString ss= (( TObjString *)obj)->GetString();
      switch (k) {
      case 1: assert(ss=="lepTH"); break;            
      case 2: tqrCut0=atof(ss.Data()); break;
      case 3: tqrCut1=atof(ss.Data()); break;
      default:
        assert(345==33987); // unexpected
      }
    }// end of while
    assert(tqrCut1>0);
    assert(tqrCut1>tqrCut0);

    delete parList;
    delete iter;
    parList = new TString(phiCut);
    objA=parList->Tokenize("_");
    iter = new TIter(objA);  
    obj = 0;
    k=0;
    while( (obj=iter->Next())) {    
      k++;
      TString ss= (( TObjString *)obj)->GetString();
      switch (k) {
      case 1: assert(ss=="lepPHI"); break;            
      case 2: phiq0=atof(ss.Data()); break;
      case 3: phiq1=atof(ss.Data()); break;
      default:
        assert(345==33987); // unexpected
      }
    }// end of while
    assert(phiq1>0);
    assert(phiq1>phiq0);

    //Constants that must be set:
    double pi = 4.0*atan(1.0);
    
    
    double tkCut0 = 0;
    double tkCut1 = pi-tkCut0;
    

    //=========================================================
    //  Initialize Generator
    RadMoller_Gen* rMollerGen = new RadMoller_Gen;

    if(Moller==true){rMollerGen->SetMoller();}
    else if(Bhabha==true){rMollerGen->SetBhabha();}
    else{
        cout<<"Improper process selection. Pick '-M' or '-B'"<<endl;
        return 0;
        }

    RandT3 *randomGenerator = new RandT3;

    rMollerGen->setRandom(randomGenerator);
    rMollerGen->SetRadFrac(radFrac);
    rMollerGen->SetTCuts(tqrCut0,tqrCut1,phiq0,phiq1);    
    rMollerGen->SetECut(dE_frac);
    rMollerGen->SetLumi(Lumi);
    rMollerGen->SetTBeam(Tbeam);
    rMollerGen->SetpRes(pRes);
    rMollerGen->SetpCells(pCells);
    rMollerGen->InitGenerator_RadMoller();
    


    //
    //======================================================================
    //=========================================================
    //=========================================================
    TFile *f;
    TFile *fout;
    FILE *eFile;
    FILE *oFile;

    // if(Moller==true){
    //     f = new TFile("histos_Moller.root","recreate");
    //     fout = new TFile("radEve_Moller.root","recreate");
    //     eFile = fopen("radEve_Moller.txt","w");
    //     oFile = fopen("summary_Moller.txt","w");
    // }
    // else if(Bhabha==true){
        f = new TFile("histos.root","recreate");
        fout = new TFile("radEve.root","recreate");
        eFile = fopen("radEve.txt","w");
        oFile = fopen("summary.txt","w");
    // }

    fprintf(oFile,"Radiative Moller/Bhabha Generator Version %.2f\n",version);
    fprintf(oFile,"Generator Summary:\n");
    fprintf(oFile,"Number of Events: %.1f\n",nEve);
    fprintf(oFile,"Beam Energy: %.1f\tDelta E Fraction: %.3f\tLuminosity: %.2f\tRadiative Fraction: %.3f\n",Tbeam,dE_frac,Lumi,radFrac);
    fprintf(oFile,"Photon Theta Cuts: %.3f\t %.3f\t  Rad Lepton Theta Cuts: %.3f\t %.3f\t Elastic Lepton Theta Cut: %.3f\n",tkCut0,tkCut1,tqrCut0,tqrCut1,1);
    if(txt_flag){
        fprintf(eFile,"Weight\tP1x\tP1y\tP1z\tP2z\tP2y\tP2z\tKx\tKy\tKz\n");
    }

    TH2D *hst_xy = new TH2D("hst_xy" ,  "Scattered Lepton Energy vs. Angle", 180,0,90,501,0.5,Tbeam);
    TH2D *ptpz = new TH2D("ptpz", "Lepton Pt vs Pz",200,0,Tbeam,200,0,23);
    TH2D *ptpzg = new TH2D("ptpzg", "Photon Pt vs Pz",200,0,Tbeam,200,0,23);
    TH2D *aa = new TH2D("aa","l1 Angle vs l2 Angle",180,0,90,180,0,90);
    TH2D *pp = new TH2D("pp","l1 Momentum vs l2 Momentum",199,0.5,Tbeam,199,0.5,Tbeam);
    TH2D *eg = new TH2D("eg","l1 Angle vs Photon Angle",360,0,180,180,0,90);
    TH1D *kk = new TH1D("kk","Ek distribution",50,0,Tbeam);
    TH1D *CSint = new TH1D("CSint","Photon ",100,0.004,0.008);

    TH1D *conv = new TH1D("conv","Convergence",100,1.e-3*45.445,0.23);

    TH1D *weights = new TH1D("weights","stuff",1,1,1);
    TH1D *eweights = new TH1D("wweights","stuff",1,1,1);

    TH1D *w_int = new TH1D("w_int","CS vs dE at 20 deg",100,1.e-3*45.445,22.6675);
    TH1D *w_int1 = new TH1D("w_int1","CS vs dE at 40 deg",100,1.e-3*45.445,22.6675);
    TH1D *w_int2 = new TH1D("w_int2","CS vs dE at 70 deg",100,1.e-3*45.445,22.6675);

    TH1D *ph_angles = new TH1D("ph_angles","Photon Angular Dist",89,0.5,89.5);
    TH2D *ph_erg = new TH2D("ph_erg","Photon Energy Angular Dist",89,0.5,89.5,101,0,Tbeam);
        
    TNtuple *radEve_t = new TNtuple("radEve_t", "Radiative Events", "p1x:p1y:p1z:p2x:p2y:p2z:k1x:k1y:k1z:weight");
    TNtuple *softEve_t = new TNtuple("softEve_t", "Soft Corrected Events", "p1x:p1y:p1z:p2x:p2y:p2z:weight");
    
    conv->Sumw2();
    hst_xy->Sumw2();
    ptpz->Sumw2();
    ptpzg->Sumw2();
    aa->Sumw2();
    pp->Sumw2();
    eg->Sumw2();
    kk->Sumw2();
    CSint->Sumw2();
    w_int->Sumw2();
    w_int1->Sumw2();
    w_int2->Sumw2();
    weights->Sumw2();
    eweights->Sumw2();

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

    w_int->SetXTitle("dEr");
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

void argUsage(char  *argv0 , const char *message){
  if (message) fprintf(stderr,"%s: %s\n",argv0,message);
  fprintf(stderr,"usage: %s  [OPTIONS]\n",argv0);
  fprintf(stderr," -M | --Moller                     : generator:  Moller\n");
  fprintf(stderr," -B | --Bhabha                         : generator:  Bhabha\n");
  fprintf(stderr," -T | --thetaCut  <lepTH_[lower]_[upper]>  : theta cuts, use --longHelp for details  \n");
  fprintf(stderr," -P | --phiCut  <lepPHI_[lower]_[upper]>  : phi cuts, use --longHelp for details  \n");
  fprintf(stderr," -E | --Tbeam <100>                      : Beam kinetic Energy [MeV] \n");
  fprintf(stderr," -f | --radFrac <0.75>                      : Fraction of events that are radiative \n");
  fprintf(stderr," -L | --Lumi <1e30>                      : Luminosity in cm^-2s^-1 \n");
  fprintf(stderr," -d | --dE_frac <1e-4>                      : Delta-E as a fraction of sqrt(s) \n");
  fprintf(stderr," -n | --nEve <10>                      : number of generated events \n");
  fprintf(stderr," -r | --root_flag                      : output TNtuple? \n");
  fprintf(stderr," -t | --txt_flag                       : output txt file? \n");
  fprintf(stderr," -s | --pRes                       : number of photon integrators \n");
  fprintf(stderr," -b | --pCells                       : number of cells per TFoam \n");

  fprintf(stderr," -h | --help         : this short help\n");
  fprintf(stderr," -H | --longHelp     : detailed explanation of switches\n");
  if(message) exit(-1);
  return;
}

//----------------------------------------
void argUsageLong(const char *message){
  fprintf(stderr,"\n -------------------------------\nDetailed explanation of switches %s\n",message); 


  fprintf(stderr,"\n*)------ Generator Selection --------\n"
      "enabled by switch  --Moller or  --Bhabha \n"
      "Expected content of --thetaCut  lepTH_[ang1]_[ang2] \n" 
      "   the range of theta for the electron (positron) is [ang1, ang2], float values in radians\n" 
      " \n"
      " Example: to generate mollers for theta=[0.01,3.14] and phi=[0.0,6.28] use flags  \n"
      "           ...  -M  -T lepTH_0.01_3.14 -P lepPHI_0.0_6.28 .... \n"
      " or Bhabha   ...  -B  -T lepTH_0.01_3.14 -P lepPHI_0.0_6.28  .... \n"
      " Note - default phi range is [0,2pi]\n"
      "\n");



}

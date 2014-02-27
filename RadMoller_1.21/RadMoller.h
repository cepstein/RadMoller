#include "TLorentzVector.h"
#include "TH1.h"
#include "TMath.h"

class RadMoller_Gen {
    public:
        void InitGenerator_RadMoller();
        void SetOutputFlags(int,int);
        void SetRadFrac(double);
        void SetTCuts(double,double,double);
        void SetECut(double);
        void SetLumi(double);
        void SetTBeam(double);
        void Generate_Event();
        int GetElFlag(){return elFlag;}
        double GetWeight(){return weight;}
        double GetHistRandom(TH1D*);
        TLorentzVector* Getq1(){return q1;}
        TLorentzVector* Getq2(){return q2;}
        TLorentzVector* Getk(){return k;}

    private:
        double randomGen();
        //Output Flags
        int root_flag; 
        int txt_flag; 

        int elFlag;

        double radFrac;

        double tkCut;
        double tqrCut;
        double xeCut;

        double dE_frac;

        double Lumi;
        double Tbeam; 
        double M2(double);
        double tree_cs(double);
        double mCSfunc(double,double);
        double dSigmahdEkdTkdTqr(double,double,double,double,\
            TLorentzVector*,TLorentzVector*,TLorentzVector*);
        double corr_soft_tsai(double,double);
        double Mh2(TLorentzVector*, TLorentzVector*, TLorentzVector*, TLorentzVector*, TLorentzVector*);

        TH1D *ekDist;
        TH1D *tkDist;
        TH1D *tqrDist; 

        double pqr;
        double weight;

        double Ek;
        double tk;
        double tqr;
        double* fIntegral;
        double ekWeight;
        double tkWeight;
        double tqrWeight;

        double pqrWeight;

        double phik;
        double xe; 
        double pickProc;


        double pqcm;
        double kp1;
        double kp2;
        double kp3;
        double p1p2;
        double p1p3;


        double me;
        double Ebeam;
        double alpha;
        double pi;
        double twopi;
        double Pbeam;
        double betacm;
        double gammacm;
        double Ecm;
        double Pcm;//momentum of either beam in CM frame
        double Ecmp; //Ecm per particle
        double Pcmp;//momentum of either beam in CM frame
        double r0;
        double ec; //electron charge
        double se; //Elastic Mandelstam S ("s" was unavailable)
        double dE;
        double EkMax;


        TLorentzVector *cm;
        TLorentzVector *p1;
        TLorentzVector *p2;
        TLorentzVector *qcm;
        TLorentzVector *qr;

        TLorentzVector *q1r;
        TLorentzVector *q2r;

        TLorentzVector *q1cm;
        TLorentzVector *q2cm;

        TLorentzVector *q1;
        TLorentzVector *q2;
        TLorentzVector *k;

        double te(double);
        double ue(double);
        double sqr(double);




};

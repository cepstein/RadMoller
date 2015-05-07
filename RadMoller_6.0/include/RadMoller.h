#ifndef __RADMOLLER_H__
#define __RADMOLLER_H__

#include "TLorentzVector.h"
#include "TH1.h"
#include "TMath.h"
// #include "GenIntegrator.h"
#include "RandGen.h"
#include "TF1.h"
#include "Math/WrappedTF1.h"
#include "Math/GSLIntegrator.h"
#include "TFoam.h"
#include "TRandom.h"
class RadMoller_Gen {
    public:
        RadMoller_Gen(){}

        virtual ~RadMoller_Gen(){}
        
        void InitGenerator_RadMoller();
        void SetpRes(int);
        void SetpCells(int);

        void SetRadFrac(double);
        void SetTCuts(double,double,double,double);
        void SetECut(double);
        void SetLumi(double);
        void SetTBeam(double);
        void SetCM();
        void SetLab();
        void SetMoller();
        void SetBhabha();
        void Generate_Event();
        int GetElFlag(){return elFlag;}
        double GetWeight(){return weight;}
        // double GetHistRandom(TH1D*);
        TLorentzVector* Getq1(){
            if(CM_flag==1){
                return q1cm;}
            else{
                return q1;}
        }
        TLorentzVector* Getq2(){            
            if(CM_flag==1){
                return q2cm;}
            else{
                return q2;}
        }
        TLorentzVector* Getk(){
            if(CM_flag==1){
                return kcm;}
            else{
                return k;}
        }
        double mCSfunc(double,double);
        double dEr(){return Ek;}
        void setRandom(RandGen*);

    private:
        class TFDISTRAD;
        friend class TFDISTRAD;
        class GenRandTR;
        friend class GenRandTR;

        double bremInt(double,double,double,double,double);
        // GenIntegrator **photonArray;


        TFoam **photonInt;
        TFDISTRAD **RHO;

        RandGen *random;
        GenRandTR *PseRan;
        Double_t *MCvect;
        Double_t MCResult;
        Double_t MCError;
        // Double_t MCwt;

        double preWeight;

        double *photonCoords;
        int mb_flag;
        int CM_flag;
        int pRes;
        int pCells;

        int elFlag;
        double symWeight(double, double);
        double radFrac;
        double hbarc2;
        double tkCut0;
        double tqrCut0;
        double tkCut1;
        double tqrCut1;

        double phiq0;
        double phiq1;

        double xeWeight;

        double dE_frac;
        double k0;
        double cosaMax;

        double Lumi;
        double Tbeam; 
        double M2(double);
        double M2b(double);
        double tree_cs(double);
        double tree_cs_b(double);
        double bCSfunc(double,double);
        double bremCS(double,\
            TLorentzVector*,TLorentzVector*);
        double bremCSb(double,\
            TLorentzVector*,TLorentzVector*);

        // double soft_cs_tsai(double,double);
        // double soft_bhabha(double,double);
        double Mh2(TLorentzVector*, TLorentzVector*, TLorentzVector*, TLorentzVector*);
        double Mh2b(TLorentzVector*, TLorentzVector*, TLorentzVector*, TLorentzVector*);

        double eps1u(TVector3*, TLorentzVector*);
        double eps1l(TVector3*, TLorentzVector*);
        double q1u(TVector3*, TLorentzVector*);
        double q1l(TVector3*, TLorentzVector*);

        double ekFunc(double);
        double tkFunc(double,double);
        double tqrFunc(double);
        double tqrFunc_Moller(double);

        double ekiCDF(double);
        double tkCDF(double,double);
        double tqriCDF(double);
        double tqriCDF_Moller(double);
        double tkInvert(double,double);

        Double_t SoftPhoton_Moller_Integrand(Double_t*,Double_t*);
        double SoftPhoton_Moller_Integral();
        double SoftPhoton_Moller(double,double);

        Double_t SoftPhoton_Bhabha_Integrand(Double_t*,Double_t*);
        double SoftPhoton_Bhabha_Integral();
        double SoftPhoton_Bhabha(double,double);

        double pqr;
        double weight;
        double eFlag;
        double Ek;
        double tk;
        double tqr;
        double ekWeight;
        double tkWeight;
        double tqrWeight;
        double phikWeight;
        double pqrWeight;

        double phik;
        double xe; 
        double pickProc;

        double SD;
        double UD;
        double TD;

        double pqcm;
        double kp1;
        double kp2;
        double kp3;
        double p1p2;
        double p1p3;

        double dE;
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
        double ec; //electron charge
        double se; //Elastic Mandelstam S ("s" was unavailable)
        double EkMax;

        double aFlag; //ambiguity flag

        TLorentzVector *cm;
        TLorentzVector *p1;
        TLorentzVector *p2;
        TLorentzVector *qcm;
        TLorentzVector *qr;

        TVector3 *q1r;
        TVector3 *newAxis;

        TLorentzVector *q2r;

        TLorentzVector *q1cm;
        TLorentzVector *q2cm;

        TLorentzVector *q1;
        TLorentzVector *q2;
        TLorentzVector *k;
        TLorentzVector *kcm;

        double te(double);
        double ue(double);
        double sqr(double);




};
#endif

Radiative Moller/Bhabha Generator
Charles Epstein, MIT
Version 3.3, August 2014


===	Overview	===

This generator outputs two types of events: 

Soft-corrected elastic events-

Accounts for soft radiative corrections, ie σ'=(1+δ)σ where soft (undetectable) photons are emitted.  

The correction factor for Moller scattering consists of the corrections of Tsai (1960).  For Bhabha scattering, the soft corrections of (A.B. Arbuzov, E.S. Scherbakova) are taken into account.  

Bremsstrahlung events-

Events consisting of ee -> eey.  The squared matrix elements for radiative Moller and Bhabha scattering were calculated using the Mathematica plugins FeynArts and FormCalc.  The CM phase-space parametrization of Haug & Nakel (The Elementary Process of Bremsstrahlung, 2003) was used to translate these to a cross-section.  A phase-space factor of two for Bhabha scattering relative to Moller scattering is included both at radiative and tree-levels.

The user may specify what fraction of generated events are radiative events in order to improve statistics for one type of process.  However, all generated events must be examined in order to determine lepton rates, i.e., one cannot simply look at the bremsstrahlung events because these do not account for photons of energy less than dE.

===	Compiling and Running	===

You will need installations of ROOT and CMAKE.  In the directory RadMoller_###, run the commands:

mkdir build
cd build
cmake ..
make

=== Command-Line Usage ===

Operating help:

  usage: ./RadMoller  [OPTIONS]
   -M | --Moller                     : generator:  Moller
   -B | --Bhabha                         : generator:  Bhabha
   -T | --thetaCut  <lepTH_[lower]_[upper]>  : theta cuts, use --longHelp for details  
   -P | --phiCut  <lepPHI_[lower]_[upper]>  : phi cuts, use --longHelp for details  
   -E | --Tbeam <100>                      : Beam kinetic Energy [MeV] 
   -f | --radFrac <0.75>                      : Fraction of events that are radiative 
   -L | --Lumi <1e30>                      : Luminosity in cm^-2s^-1 
   -d | --dE_frac <1e-4>                      : Delta-E as a fraction of sqrt(s) 
   -n | --nEve <10>                      : number of generated events 
   -r | --root_flag                      : output TNtuple? 
   -t | --txt_flag                       : output txt file? 

   -h | --help         : this short help
   -H | --longHelp     : detailed explanation of switches

	------ Generator Selection --------
      enabled by switch  --Moller or  --Bhabha 
      Expected content of --thetaCut  lepTH_[ang1]_[ang2]  
         the range of theta for the electron (positron) is [ang1, ang2], float values in radians 
       
       Example: to generate mollers for theta=[0.01,3.14] and phi=[0.0,6.28] use flags  
                 ...  -M  -T lepTH_0.01_3.14 -P lepPHI_0.0_6.28 .... 
       or Bhabha   ...  -B  -T lepTH_0.01_3.14 -P lepPHI_0.0_6.28  .... 
       Note - default phi range is [0,2pi]

Example:

./RadMoller -M -n 10000 -E 2000 -f 1 -T lepTH_0.01_3.14 -P lepPHI_0.0_6.28 -d 1e-3 -L 1e36 -r -t

Or, using some defaults:

./RadMoller -M -n 10000 -E 100

(etc)

===	Further Usage Info ===

User-specified parameters are:

int root_flag = 1; //1 = Output TNtuple of event parameters
int txt_flag = 1;  //1 = Output txt file of event parameters

double radFrac = 0.5; //Fraction of events that are radiative 

—

Phase-space cuts on outgoing particles, in CM system coordinates:

tqrCut0: Lower theta cut for emitted electron 1 (or positron, for Bhabha)
tqrCut1: (Upper)
phiq0: Lower phi cut for emitted electron 1 (or positron, for Bhabha)
phiq1: (Upper)

—

double dE_frac = 1.e-4;

This sets the crossover between what is considered a hard and a soft photon.  It is expressed in units of sqrt(s), ie, the center of mass energy.

—
double Lumi = 6.e35;

This is the luminosity in cm^2s^-1 in order to convert the output from cross-section to rates. 

—

double Tbeam = 100.; 

This is the incident beam kinetic energy, in MeV.


===	Output File Summary	===

The generator has up to three outputs.

(1) Histogram file: histos.root.  By running 'root -l ../plotHistos.C' from the build directory, one can see various kinematic plots.

(2) Txt file: Containing 10 tab-delimited columns: weight, p1 (x,y,z), p2 (x,y,z), k (x,y,z).

(3) Root TNtuple: Containing a set of bremsstrahlung events, and a set of soft-corrected events (stored separately because the number of parameters are different - but all must be used to obtain physical results).  

===	Weighting Scheme	===

The weights follow the convention of (cross-section)*(luminosity): ie, they are not normalized to the number of events in the set.  With this scheme, the integrated output rates must be divided by the number of events in the set in order to obtain physical rates.  This is the same scheme as the DL MadGraph event sets but multiplied by the luminosity (which could be set by the user to be unity).

The generator will output all events for which "Lepton 1" is in the specified region.  For Bhabha events, this means the positron.  For Moller events, the generator will give all events that have _an_ electron in the region.  Events are ensured to not be overcounted by the internal function symWeight().  This checks whether the second electron falls into the region over which the first electron is being integrated - if it is outside the region, a factor of two is applied to account for the event in which the "primary" electron is outside and the secondary is inside.  
 
=== Implementing the Generator in Outside Applications ===

First, import "RadMoller.h" and add the three src files.  Provide an implementation of the random number generator.  For example: 

TRandom3 *randGen = new TRandom3(0);
double RadMoller_Gen::randomGen(){
    return randGen->Uniform();
}


Instantiate an instance of the generator by:

RadMoller_Gen* rMollerGen = new RadMoller_Gen;

Choose Moller Scattering:
rMollerGen->SetMoller();
or Bhabha scattering:
rMollerGen->SetBhabha();


Set the flags and parameters (defined above) by:
rMollerGen->SetRadFrac(radFrac);
rMollerGen->SetTCuts(tqrCut0,tqrCut1,phiq0,phiq1);
rMollerGen->SetECut(dE_frac);
rMollerGen->SetLumi(Lumi);
rMollerGen->SetTBeam(Tbeam);

Then initialize the generator:
rMollerGen->InitGenerator_RadMoller();

Generate an event by calling: 
rMollerGen->Generate_Event();

And access the output by:


rMollerGen->GetElFlag(); 
---> this is a flag of whether the event was elastic (2 e-) or radiative (2e-, 1y). (1 = radiative)
rMollerGen->GetWeight(); 
---> access the weight of the event

rMollerGen->Getq1(); 
--->Electron (positron) outgoing momentum

rMollerGen->Getq2();
--->Electron outgoing momentum

rMollerGen->Getk();
--->Photon momentum



===	Development		===

There are almost certainly bugs in the code that may or may not have physical consequences.  Please send physics concerns / discovered bugs to cepstein@mit.edu.
 

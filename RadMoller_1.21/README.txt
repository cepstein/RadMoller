Radiative Moller Generator
Charles Epstein, MIT
Version 1.21, February 2014


===	Overview	===

This generator outputs two types of events: 

Soft-corrected elastic events-

Accounts for soft radiative corrections, ie σ'=(1+δ)σ —-> σ'=exp(δ)σ where soft (undetectable) photons are emitted.  The exponentiation	accounts for the effects of multiple soft photons.

The correction factor consists of the corrections of Tsai (1960).  Previously included was the option for the corrections of Denner & Pozzorini (1999), which are temporarily removed from the release.  The latter is more relevant at high energies and includes weak interaction effects.

Bremsstrahlung events-

Events consisting of ee -> eey.  The squared matrix element was	calculated using the Mathematica plugins FeynArts and FormCalc.  The convenient phase-space parametrization of Petriello, 2003 was used to translate this to a cross-section.

The user may specify what fraction of generated events are radiative events in order to improve statistics for one type of process.  However, all generated events must be examined in order to determine electron rates, i.e., one cannot simply look at the bremsstrahlung events because these do not account for photons of energy less than dE.

===	Compiling and Running	===

You will need installations of ROOT and CMAKE.  In the directory RadMoller_1.0, run the commands:

mkdir build
cd build
cmake ..
make
./RadMoller [nEve]

where [nEve] is the desired number of events.  It is recommended to keep this in sets of less than ~(a few)*10^6 in order to prevent output file sizes from becoming unreasonable.  Generation of such a number of events typically takes less than a minute.


===	Usage		===

User-specified parameters are set in int main().  They are:

int root_flag = 1; //1 = Output TNtuple of event parameters
int txt_flag = 1;  //1 = Output txt file of event parameters

double radFrac = 0.5; //Fraction of events that are radiative 

—

Phase-space cuts on outgoing particles:

double tkCut = 0.001;
double tqrCut = 0.001;
double xeCut = 0.001;

These should not be used to place cuts on lab-frame phase-space, as they are defined for various non-lab-frame coordinate systems.  The first two refer to the bremsstrahlung events: tkCut is the cut in the CM polar photon angle, tqrCut is the cut in the polar angle of the electrons in their center of mass frame.  For elastic events, xeCut is the cut in the elastic electron angle in the CM system. They are non-zero to avoid singularities.

—

double dE_frac = 1.e-3;

This sets the crossover between what is considered a hard and a soft photon.  It is expressed in units of s, ie, Ecm^2.

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
 

===	Development		===

There are almost certainly bugs in the code that may or may not have physical consequences.  Please send physics concerns / discovered bugs to cepstein@mit.edu.


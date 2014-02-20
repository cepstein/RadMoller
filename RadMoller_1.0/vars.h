//Constants that must be set:

//Output Flags
const int soft_flag = 1; //0 = Denner & Pozzorini, 1 = Tsai
const int root_flag = 1; //1 = Output TNtuple of event parameters
const int txt_flag = 1;  //1 = Output txt file of event parameters

//Phase-space cuts on outgoing particles
const double tkCut = 0.001;
const double tqrCut = 0.001;
const double xeCut = 0.001;

const double dE_frac = 1.e-3;//fraction of S

//const double Lumi = 1.e30; //cm^2 s^-1 - gives CS in microbarns
const double Lumi = 6.e35;

//Beam Kinetic Energy
const double Tbeam = 100.; 

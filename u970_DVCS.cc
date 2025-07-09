#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <cmath>
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TLorentzRotation.h"
#include "Phast.h" 
#include "PaSetup.h"  
#include "PaAlgo.h"
#include "PaEvent.h" 
#include "PaHodoHelper.h"
#include "G3part.h" 
#include <utility>

#include "ecal_time_cuts.h"
#include "Tools_Camera.hh"
#include "PAMBH_interface_HepGEN.h"

//#include "/Users/gursimran/cern/phastPackages/xcheck_newTiS/tis_range.cc"
//#include "/Users/gursimran/cern/phastPackages/xcheck_newTiS/tis_range.h"
#include "Toololols_KinFitterInterface.hh"
#include "UConn_Tools.h" 

#include "/afs/cern.ch/user/g/gkainth/phastPackages/xcheck_newTiS/tis_range.cc"
#include "/afs/cern.ch/user/g/gkainth/phastPackages/xcheck_newTiS/tis_range.h"
//#include "/afs/cern.ch/user/g/gkainth/phastPackages/git_COMPASS/UConn_Tools.h"
//#include "/afs/cern.ch/user/g/gkainth/phastPackages/flux_files/flux_Johannes/functions/TiS_range.cc"

// ************************************************************************** //
// UserEvent for preselecting exclusive photon events	(DVCS)			            //
// In the selection all possible combinations of:					                    //
// Vertices (incoming and outgoing muons), exclusive photon and recoil proton //
// ************************************************************************** //

extern "C" float prob_(float&, int&); 

// Event statistic counter flags  
EventFlags eventFlags(EventFlags::DVCS); // Create an instance of flags and counters 

// Global user selection flags 
bool verbose_mode = false; // Create an instance of verbose_mode
// this is a good example of what can be an environmental variable 
bool leptoMC      = false; // Set to true for LEPTO MC to remove events with exclusive event topology 

//*****************************************************************************
void UserEvent970(PaEvent & e) { // begin event loop

    // Define constants
    static PaHodoHelper* HodoHelper = nullptr; 
    static TiSRange*     tis_range  = nullptr; 
    static PaCamera*     cam_inst   = nullptr;

    //const double M_pi    = G3partMass[8];  // Pion mass 
    //const double M_gamma = G3partMass[1];  // Photon mass
    const double M_mu = G3partMass[5];  // Muon mass
    const double M_p  = G3partMass[14]; // Proton mass 

    // Declare all objects globally
    // Add histograms as well if needed (ex. static TH1F* h97_Zprim  = NULL;)
    static TTree* tree(NULL);     // tree for sotring real data passed through the full event selection process 
    static TTree* tree_MC(NULL);  // tree for storing MC data passed through the full event selection process 
    static TTree* tree_gen(NULL); // tree for storing HEPGEN data prior to event selection (aceptance study)

    //
    // Variables to be stored into analysis Ntuple
    //
    // (if you want to extend the list:
    // - add variable here,
    // - add in Ntuple definition below 
    // and do not forget to fill the variable in the code
    //
    //*******************************************
    // Shared variables for real and MC data    
    static unsigned long long int Evt; // event number - unique evt number (based on run, spill and evt num in spill) 
    static int    Run;          // run number
    static int    LastRun = -1; // store the previous run number (used to reintialize Hodohelper, tis_range if there are multiple runs)
    static int    EvtInSpill;   // event number in spill 
    static int    Spill;        // spill number
    static double TimeInSpill;  // time in spill  
    static float  Chi2;         // Chi2 of the reconstructed vertex 
    static int    Nprim;        // Number of tracks in primary vertex (-1 in fot found)
    static int    Q_beam;       // Charge of the beam 
    static int    trig_mask;

    static TVector3 pVtx_vec;     // position vector for the primary vertex (X, Y, Z)
    static TVector3 posRingA_vec; // hit position of the proton in CAMERA ring A
    static TVector3 posRingB_vec; // hit position of the proton in CAMERA ring B

    static TLorentzVector inMu_TL;               // energy-momentum four vector of the beam muon
    static TLorentzVector outMu_TL;              // energy-momentum four vector of the scattered muon
    static TLorentzVector targ_TL(0, 0, 0, M_p); // energy momentum four-vector of the target proton 
    static TLorentzVector gamma_TL;              // energy-momentum four vector of the real photon 
    static TLorentzVector cluster_TL;            // vector used to store the position and energy information of the cluster
    static TLorentzVector p_camera_TL;           // energy-momentum four vector of the proton measured by camera
    static TLorentzVector gammaLow_TL;           // energy-momentum four vector of the low energy photon (used to find pi0)
    static TLorentzVector clusterLow_TL;         // vector used to store the position and energy information of the low-energy cluster

    static int low_calo;    // calorimeter ID for low energy cluster used for pi0 reconstruction 
    static int excl_calo;   // calorimeter ID for exclusive photon cluster  

    static double y;   // fractional energy loss of the incoming lepton 
    static double nu;  // energy of the virtual photon 
    static double Q2;  // four-momentum transfer sqaured (Q is the four-momentum transferred between the incoming and outgoing lepton)
    static double W2;  // effective mass of the final state hadron system squared 
    static double xbj; // measure of the elasticity of the scattering process 
    static double t;   // four-momentum transfer squared of the nucleon 

    // Exclusivity variables 
    static double delta_phi; 
    static double delta_pt;
    static double delta_Z;
    static double M2x;

    // Kinematic fit;
    static TVector3 pVtxFit_vec; 
    static TVector3 posRingAFit_vec; 
    static TVector3 posRingBFit_vec; 
    static TVector3 clusterFit_vec; // x, y and E of the cluster post kinematic fit 

    static TLorentzVector inMuFit_TL; 
    static TLorentzVector outMuFit_TL; 
    static TLorentzVector targetFit_TL;  
    static TLorentzVector gammaFit_TL; 
    static TLorentzVector protonFit_TL;  

    static double inMu_sigmaX; 
    static double inMu_sigmaY;
    static double inMu_sigmaPx;
    static double inMu_sigmaPy;
    static double inMu_sigmaPz;

    static double outMu_sigmaX;
    static double outMu_sigmaY;
    static double outMu_sigmaPx;
    static double outMu_sigmaPy;
    static double outMu_sigmaPz;

    static double gamma_sigmaX;
    static double gamma_sigmaY;
    static double gamma_sigmaE;

    static double proton_sigmaP; 
    static double proton_sigmaTheta;
    static double proton_sigmaPhi;

    static double ringA_sigmaR;
    static double ringA_sigmaPhi;
    static double ringA_sigmaZ;

    static double ringB_sigmaR;
    static double ringB_sigmaPhi;
    static double ringB_sigmaZ;

    static double chi2_fit; // chi2 of the fit  
    static int    ndf_fit;  // ndf of the fit 

    static double y_fit;
    static double nu_fit;
    static double Q2_fit; 
    static double t_fit;   

    //*******************************************
    // HEPGEN BH data (prior to event selection) -> used for acceptance study 
    static TLorentzVector inMu_gen_TL ; 
    static TLorentzVector outMu_gen_TL; 
    static TLorentzVector gamma_gen_TL; 
    static TLorentzVector proton_gen_TL;
    static TLorentzVector q_gen;         // four momentum of the virtual photon

    static double y_gen; 
    static double Q2_gen; 
    static double xbj_gen;
    static double nu_gen; 
    static double W2_gen; 
    static double t_gen;  
    static double phi_gg_gen; // azimuthal angle between the virtual and real photon production planes 

    static double weight_all; 
    static double weight_DVCS; 
    static double weight_BH; 
    static double weight_Interference; 
    static double phase_fac; 
    static double weight_PAMBH; 

    //*******************************************
    // Event selection flags (not the same as statistic counter flags)
    bool trig_MT   = false;
    bool trig_LT   = false;
    bool trig_OT   = false;
    bool trig_LAST = false;
    bool trig_flag = false; 
    bool TiS_flag  = true;
    bool fit_conv  = false; 
    bool save_evt  = false; 

    //*****************************************************************************
    static bool first(true);
    if (first) { // histograms and Ntuples booking block
        Phast::Ref().HistFileDir("UserEvent970");

        //
        // Ntuple definition 
        //
        //*******************************************
        tree = new TTree("USR970","User 970 DVCS NTuple - Real Data"); // name (has to be unique) and title of the Ntuple
        // Basic event information
        tree->Branch("Run", &Run, "Run/I");
        tree->Branch("Evt", &Evt, "Evt/l");
        tree->Branch("Chi2", &Chi2, "Chi2/F");
        tree->Branch("Spill", &Spill, "Spill/I");
        tree->Branch("Q_beam", &Q_beam, "Q_beam/I");
        // Trigger information
        tree->Branch("trig_MT", &trig_MT, "trig_MT/O");
        tree->Branch("trig_LT", &trig_LT, "trig_LT/O");
        tree->Branch("trig_OT", &trig_OT, "trig_OT/O");
        tree->Branch("trig_LAST", &trig_LAST, "trig_LAST/O");
        // Particle/vertex vectors 
        tree->Branch("inMu_TL", &inMu_TL);
        tree->Branch("outMu_TL", &outMu_TL);
        tree->Branch("gamma_TL", &gamma_TL);
        tree->Branch("cluster_TL", &cluster_TL);
        tree->Branch("p_camera_TL", &p_camera_TL);
        tree->Branch("gammaLow_TL", &gammaLow_TL);
        tree->Branch("clusterLow_TL", &clusterLow_TL);
        tree->Branch("pVtx_vec", &pVtx_vec); 
        tree->Branch("posRingA_vec", &posRingA_vec);
        tree->Branch("posRingB_vec", &posRingB_vec);
        // Kinematic variables 
        tree->Branch("y", &y, "y/D");
        tree->Branch("t", &t, "t/D");
        tree->Branch("nu", &nu, "nu/D");
        tree->Branch("Q2", &Q2, "Q2/D");
        tree->Branch("W2", &W2, "W2/D");
        tree->Branch("xbj", &xbj, "xbj/D");
        // Invariant mass of visible pi0
        tree->Branch("low_calo", &low_calo, "low_calo/I");
        tree->Branch("excl_calo", &excl_calo, "excl_calo/I");
        // Exclusivity variables      
        tree->Branch("M2x", &M2x, "M2x/D");
        tree->Branch("delta_Z", &delta_Z, "delta_Z/D");
        tree->Branch("delta_pt", &delta_pt, "delta_pt/D");
        tree->Branch("delta_phi", &delta_phi, "delta_phi/D");
        // Kinematic fit vectors 
        tree->Branch("inMuFit_TL", &inMuFit_TL);
        tree->Branch("outMuFit_TL", &outMuFit_TL);
        tree->Branch("gammaFit_TL", &gammaFit_TL);
        tree->Branch("protonFit_TL", &protonFit_TL);
        tree->Branch("pVtxFit_vec", &pVtxFit_vec);
        tree->Branch("clusterFit_vec", &clusterFit_vec);
        tree->Branch("posRingAFit_vec", &posRingAFit_vec);
        tree->Branch("posRingBFit_vec", &posRingBFit_vec);
        // Kinematic fit 
        tree->Branch("y_fit", &y_fit, "y_fit/D");
        tree->Branch("t_fit", &t_fit, "t_fit/D");
        tree->Branch("nu_fit", &nu_fit, "nu_fit/D");
        tree->Branch("Q2_fit", &Q2_fit, "Q2_fit/D");
        tree->Branch("ndf_fit", &ndf_fit, "ndf_fit/I");
        tree->Branch("chi2_fit", &chi2_fit, "chi2_fit/D");
        tree->Branch("fit_conv", &fit_conv, "fit_conv/O");
        tree->Branch("inMu_sigmaX", &inMu_sigmaX, "inMu_sigmaX/D");
        tree->Branch("inMu_sigmaY", &inMu_sigmaY, "inMu_sigmaY/D");
        tree->Branch("inMu_sigmaPx", &inMu_sigmaPx, "inMu_sigmaPx/D");
        tree->Branch("inMu_sigmaPy", &inMu_sigmaPy, "inMu_sigmaPy/D");
        tree->Branch("inMu_sigmaPz", &inMu_sigmaPz, "inMu_sigmaPz/D");
        tree->Branch("outMu_sigmaX", &outMu_sigmaX, "outMu_sigmaX/D");
        tree->Branch("outMu_sigmaY", &outMu_sigmaY, "outMu_sigmaY/D");
        tree->Branch("outMu_sigmaPx", &outMu_sigmaPx, "outMu_sigmaPx/D");
        tree->Branch("outMu_sigmaPy", &outMu_sigmaPy, "outMu_sigmaPy/D");
        tree->Branch("outMu_sigmaPz", &outMu_sigmaPz, "outMu_sigmaPz/D");
        tree->Branch("gamma_sigmaX", &gamma_sigmaX, "gamma_sigmaX/D");
        tree->Branch("gamma_sigmaY", &gamma_sigmaY, "gamma_sigmaY/D");
        tree->Branch("gamma_sigmaE", &gamma_sigmaE, "gamma_sigmaE/D");
        tree->Branch("proton_sigmaP", &proton_sigmaP, "proton_sigmaP/D");
        tree->Branch("proton_sigmaPhi", &proton_sigmaPhi, "proton_sigmaPhi/D");
        tree->Branch("proton_sigmaTheta", &proton_sigmaTheta, "proton_sigmaTheta/D");
        tree->Branch("ringA_sigmaR", &ringA_sigmaR, "ringA_sigmaR/D");
        tree->Branch("ringA_sigmaZ", &ringA_sigmaZ, "ringA_sigmaZ/D");
        tree->Branch("ringA_sigmaPhi", &ringA_sigmaPhi, "ringA_sigmaPhi/D");
        tree->Branch("ringB_sigmaR", &ringB_sigmaR, "ringB_sigmaR/D");
        tree->Branch("ringB_sigmaZ", &ringB_sigmaZ, "ringB_sigmaZ/D");
        tree->Branch("ringB_sigmaPhi", &ringB_sigmaPhi, "ringB_sigmaPhi/D");

        //*******************************************
        tree_MC = new TTree("USR970_MC","User 970 DVCS NTuple - MC Data"); // name (has to be unique) and title of the Ntuple
        // Basic event information
        tree_MC->Branch("Run", &Run, "Run/I");
        tree_MC->Branch("Evt", &Evt, "Evt/l");
        tree_MC->Branch("Chi2", &Chi2, "Chi2/F");
        tree_MC->Branch("Spill", &Spill, "Spill/I");
        tree_MC->Branch("Q_beam", &Q_beam, "Q_beam/I");
        // Trigger information
        tree_MC->Branch("trig_MT", &trig_MT, "trig_MT/O");
        tree_MC->Branch("trig_LT", &trig_LT, "trig_LT/O");
        tree_MC->Branch("trig_OT", &trig_OT, "trig_OT/O");
        tree_MC->Branch("trig_LAST", &trig_LAST, "trig_LAST/O");
        // Particle/vertex vectors 
        tree_MC->Branch("inMu_TL", &inMu_TL);
        tree_MC->Branch("outMu_TL", &outMu_TL);
        tree_MC->Branch("gamma_TL", &gamma_TL);
        tree_MC->Branch("cluster_TL", &cluster_TL);
        tree_MC->Branch("p_camera_TL", &p_camera_TL);
        tree_MC->Branch("gammaLow_TL", &gammaLow_TL);
        tree_MC->Branch("clusterLow_TL", &clusterLow_TL);
        tree_MC->Branch("pVtx_vec", &pVtx_vec);
        tree_MC->Branch("posRingA_vec", &posRingA_vec);
        tree_MC->Branch("posRingB_vec", &posRingB_vec);
        // Kinematic variables 
        tree_MC->Branch("y", &y, "y/D");
        tree_MC->Branch("t", &t, "t/D");
        tree_MC->Branch("nu", &nu, "nu/D");
        tree_MC->Branch("Q2", &Q2, "Q2/D");
        tree_MC->Branch("W2", &W2, "W2/D");
        tree_MC->Branch("xbj", &xbj, "xbj/D");
        // Invariant mass of visible pi0
        tree_MC->Branch("low_calo", &low_calo, "low_calo/I");
        tree_MC->Branch("excl_calo", &excl_calo, "excl_calo/I");
        // Exclusivity variables      
        tree_MC->Branch("M2x", &M2x, "M2x/D");
        tree_MC->Branch("delta_Z", &delta_Z, "delta_Z/D");
        tree_MC->Branch("delta_pt", &delta_pt, "delta_pt/D");
        tree_MC->Branch("delta_phi", &delta_phi, "delta_phi/D");
        // Kinematic fit vectors 
        tree_MC->Branch("inMuFit_TL", &inMuFit_TL);
        tree_MC->Branch("outMuFit_TL", &outMuFit_TL);
        tree_MC->Branch("gammaFit_TL", &gammaFit_TL);
        tree_MC->Branch("protonFit_TL", &protonFit_TL);
        tree_MC->Branch("pVtxFit_vec", &pVtxFit_vec);
        tree_MC->Branch("clusterFit_vec", &clusterFit_vec);
        tree_MC->Branch("posRingAFit_vec", &posRingAFit_vec);
        tree_MC->Branch("posRingBFit_vec", &posRingBFit_vec);
        // Kinematic fit
        tree_MC->Branch("y_fit", &y_fit, "y_fit/D");
        tree_MC->Branch("t_fit", &t_fit, "t_fit/D");
        tree_MC->Branch("nu_fit", &nu_fit, "nu_fit/D");
        tree_MC->Branch("Q2_fit", &Q2_fit, "Q2_fit/D");
        tree_MC->Branch("ndf_fit", &ndf_fit, "ndf_fit/I");
        tree_MC->Branch("chi2_fit", &chi2_fit, "chi2_fit/D");
        tree_MC->Branch("fit_conv", &fit_conv, "fit_conv/O");
        tree_MC->Branch("inMu_sigmaX", &inMu_sigmaX, "inMu_sigmaX/D");
        tree_MC->Branch("inMu_sigmaY", &inMu_sigmaY, "inMu_sigmaY/D");
        tree_MC->Branch("inMu_sigmaPx", &inMu_sigmaPx, "inMu_sigmaPx/D");
        tree_MC->Branch("inMu_sigmaPy", &inMu_sigmaPy, "inMu_sigmaPy/D");
        tree_MC->Branch("inMu_sigmaPz", &inMu_sigmaPz, "inMu_sigmaPz/D");
        tree_MC->Branch("outMu_sigmaX", &outMu_sigmaX, "outMu_sigmaX/D");
        tree_MC->Branch("outMu_sigmaY", &outMu_sigmaY, "outMu_sigmaY/D");
        tree_MC->Branch("outMu_sigmaPx", &outMu_sigmaPx, "outMu_sigmaPx/D");
        tree_MC->Branch("outMu_sigmaPy", &outMu_sigmaPy, "outMu_sigmaPy/D");
        tree_MC->Branch("outMu_sigmaPz", &outMu_sigmaPz, "outMu_sigmaPz/D");
        tree_MC->Branch("gamma_sigmaX", &gamma_sigmaX, "gamma_sigmaX/D");
        tree_MC->Branch("gamma_sigmaY", &gamma_sigmaY, "gamma_sigmaY/D");
        tree_MC->Branch("gamma_sigmaE", &gamma_sigmaE, "gamma_sigmaE/D");
        tree_MC->Branch("proton_sigmaP", &proton_sigmaP, "proton_sigmaP/D");
        tree_MC->Branch("proton_sigmaPhi", &proton_sigmaPhi, "proton_sigmaPhi/D");
        tree_MC->Branch("proton_sigmaTheta", &proton_sigmaTheta, "proton_sigmaTheta/D");
        tree_MC->Branch("ringA_sigmaR", &ringA_sigmaR, "ringA_sigmaR/D");
        tree_MC->Branch("ringA_sigmaZ", &ringA_sigmaZ, "ringA_sigmaZ/D");
        tree_MC->Branch("ringA_sigmaPhi", &ringA_sigmaPhi, "ringA_sigmaPhi/D");
        tree_MC->Branch("ringB_sigmaR", &ringB_sigmaR, "ringB_sigmaR/D");
        tree_MC->Branch("ringB_sigmaZ", &ringB_sigmaZ, "ringB_sigmaZ/D");
        tree_MC->Branch("ringB_sigmaPhi", &ringB_sigmaPhi, "ringB_sigmaPhi/D");

        //*******************************************
        tree_gen = new TTree("USR970_GEN","User 970 HEPGEN NTuple");
        // Basic event information
        tree_gen->Branch("Run", &Run, "Run/I");
        tree_gen->Branch("Evt", &Evt, "Evt/l");
        tree_gen->Branch("Spill", &Spill, "Spill/I");
        tree_gen->Branch("Q_beam", &Q_beam, "Q_beam/I");
        // Particle vectors 
        tree_gen->Branch("inMu_gen_TL", &inMu_gen_TL);
        tree_gen->Branch("outMu_gen_TL", &outMu_gen_TL);
        tree_gen->Branch("gamma_gen_TL", &gamma_gen_TL);
        tree_gen->Branch("proton_gen_TL", &proton_gen_TL);
        // Kinematic variables 
        tree_gen->Branch("y_gen", &y_gen, "y_gen/D");
        tree_gen->Branch("t_gen", &t_gen, "t_gen/D");
        tree_gen->Branch("Q2_gen", &Q2_gen, "Q2_gen/D");
        tree_gen->Branch("nu_gen", &nu_gen, "nu_gen/D");
        tree_gen->Branch("W2_gen", &W2_gen, "W2_gen/D");
        tree_gen->Branch("xbj_gen", &xbj_gen, "xbj_gen/D");
        tree_gen->Branch("phi_gg_gen", &phi_gg_gen, "phi_gg_gen/D");
        // Weights 
        tree_gen->Branch("phase_fac", &phase_fac, "phase_fac/D");
        tree_gen->Branch("weight_BH", &weight_BH, "weight_BH/D");
        tree_gen->Branch("weight_all", &weight_all, "weight_all/D");
        tree_gen->Branch("weight_DVCS", &weight_DVCS, "weight_DVCS/D");
        tree_gen->Branch("weight_PAMBH", &weight_PAMBH, "weight_PAMBH/D");
        tree_gen->Branch("weight_Interference", &weight_Interference, "weight_Interference/D");
        
        first = false;
    } // end of histogram booking

    //*****************************************************************************
    // Assign names to trigger bits
    enum trigger {   
      Tiger = 1<<0,
      MT = 1<<1,
      LT = 1<<2,
      OT = 1<<3,
      CT = 1<<4,
      IV = 1<<5,
      HaloT = 1<<6,
      BT = 1<<7,
      Tiger_only = 1<<8,
      LAST = 1<<9,
      TRand = 1<<10,
      NRand = 1<<11
    };

    trig_mask = e.TrigMask();
    trig_mask = trig_mask&2047;

    std::string trigCheck = ""; // empty string to store trigger information for debugging 
    if (trig_mask & MT) {
        trigCheck += "MT ";
        trig_flag = true;
        trig_MT = true; 
    }
    if (trig_mask & LT) {
        trigCheck += "LT ";
        trig_flag = true;
        trig_LT = true;
    }
    if (trig_mask & OT) {
        trigCheck += "OT ";
        trig_flag = true;
        trig_OT = true;
    }
    if (trig_mask & LAST) {
        trigCheck += "LAST ";
        trig_flag = true; 
        trig_LAST = true;
    }

    //*******************************************
    // Initialize variables, extra flags and check time in spill (time in spill cut is applied later not here)   
    // Check flags are only used to check statistics not for final event selection 
    eventFlags.createFlag("singleTrack_flag", "No. of events where primary vertex only has one outgoing track");
    eventFlags.createFlag("Q2_DIS_flag", "No. of events where Q2 > 0.5 (DIS)");
    //eventFlags.createFlag("y_DIS_flag", "No. of events where 0.01 < y < 0.99");

    eventFlags.createFlag("clAll_flag", "No. of events that have any clusters");
    eventFlags.createFlag("clNeutral_flag", "No. of events where clusters are not associated with charged tracks");
    eventFlags.createFlag("clTime_flag", "No. of events where cluster timing is within requirements");
    eventFlags.createFlag("lowECl_flag", "No. of events with low energy photons");
    eventFlags.createFlag("nCls_flag", "No. of events where high energy clusters are found in ECal 0, 1 or 2 only");
    eventFlags.createFlag("singleCl_flag", "No. of events where there is only a single high energy cluster in the ECals");

    eventFlags.createFlag("protonAll_flag", "No. of events with proton candidates");
    eventFlags.createFlag("proton_flag", "No. of events where proton candidates have 0.1 < beta < 1");
    eventFlags.createFlag("delta_pt_flag", "No. of events where |delta_pt| < 0.3 GeV/c");
    eventFlags.createFlag("delta_phi_flag", "No. of events where |delta_phi| < 0.4 rad");
    eventFlags.createFlag("delta_Z_flag", "No. of events where |delta_Z| < 16 cm");
    eventFlags.createFlag("M2x_flag", "No. of events where |(M_x)^2| < 0.3 (GeV/c^2)^2");

    eventFlags.createFlag("Q2Fit_flag", "No. of events where 1 < Q2_fit < 10");
    eventFlags.createFlag("yFit_flag", "No. of events where 0.05 < y_fit < 0.95"); //final DVCS cut is at 0.9 not 0.95 
    eventFlags.createFlag("tFit_flag", "No. of events where -0.08 < t_fit < -0.64");
    eventFlags.createFlag("nuFit_flag", "No. of events where 10 < nu_fit < 144");
    eventFlags.createFlag("kinFitAll_flag", "No. of events surviving all kinematic cuts (Q2, y, t, nu)");
    
    eventFlags.createFlag("nExclCombo_flag", "No. of events where at least one vertex, photon and proton combination satisfies 4/5 exclusivity cuts");
    eventFlags.createFlag("nExclComboPi0_flag", "No. of events with pi0 candidates which satisfy at least 4/5 exclsuvity cuts");
     
    eventFlags.resetFlags(); // Reset all event statistic counter flags to false

    Run         = e.RunNum();
    Evt         = e.UniqueEvNum();
    EvtInSpill  = e.EvInSpill(); 
    Spill       = e.SpillNum(); 
    TimeInSpill = e.TimeInSpill();  

    if (Run != LastRun) { // Reinitialize HodoHelper and tis_range only if the run number changes 
        HodoHelper = & PaHodoHelper::Init("/afs/cern.ch/user/g/gkainth/phast/dat/trigger_config/2016", true);  
        //tis_range  = new TiSRange("/Users/gursimran/cern/phastPackages/flux_files/flux_Johannes/2016/flux_files");
        tis_range  = new TiSRange("/afs/cern.ch/user/g/gkainth/phastPackages/flux_files/flux_Johannes/2016/flux_files");
        //Set_TiSrange("/afs/cern.ch/user/g/gkainth/phastPackages/flux_files/flux_Johannes/2016/flux_files", Run, Run);
        LastRun = Run;  // Update LastRun to the current run number
    }

    if (!e.IsMC()) {
      TiS_flag = tis_range->CheckWindow(Run, Spill, TimeInSpill); // check the time in spill for real data
    }

    cam_inst = & PaCamera::GetInstance();
    cam_inst->NewEvent(e, e.IsMC()); // only have it like this for the x-check for the full analysis use cam_inst->NewEvent(e); 
    //cam_inst->NewEvent(e);

    //*******************************************
    // Get MC PAM weights for HEPGEN and save to tree_gen 
    if (e.IsMC()) { // Begin loop over MC data (important only for HEPGEN BH but save it for all anyways) 
      NLUDATA ld; 
      if (e.MCgen(ld)) {
        weight_all  = ld.uservar[2]; 
        weight_DVCS = ld.uservar[15];  
        weight_BH   = ld.uservar[16];  
        phase_fac   = ld.uservar[9]; 
        weight_Interference = weight_all - weight_DVCS - weight_BH;
      }

      for (int iv = 0; iv < e.NMCvertex(); iv++) {
        const PaMCvertex &v = e.vMCvertex(iv);
        if (!v.IsPrimary()) continue; // must be a primary vertex 
        if (v.NMCtrack() != 4) continue; // must have exactly four outgoing tracks 

        const PaMCtrack &t_beam   = e.vMCtrack(v.iMCtrack(0));
			  const PaMCtrack &t_omu    = e.vMCtrack(v.iMCtrack(1));
			  const PaMCtrack &t_gamma  = e.vMCtrack(v.iMCtrack(2));
			  const PaMCtrack &t_proton = e.vMCtrack(v.iMCtrack(3));

        const PaTPar &par_beam   = t_beam.ParInVtx();
        const PaTPar &par_outMu  = t_omu.ParInVtx();
        const PaTPar &par_gamma  = t_gamma.ParInVtx();
        const PaTPar &par_proton = t_proton.ParInVtx();

        inMu_gen_TL   = par_beam.LzVec(M_mu); 
        outMu_gen_TL  = par_outMu.LzVec(M_mu); 
        gamma_gen_TL  = par_gamma.LzVec(0); 
        proton_gen_TL = par_proton.LzVec(M_p); 
        q_gen         = inMu_gen_TL - outMu_gen_TL;

        if (par_beam.Q() > 0) {
          Q_beam = 1;
        } else {Q_beam = -1;} 

        y_gen   = (inMu_gen_TL.E() - outMu_gen_TL.E()) / inMu_gen_TL.E();
        Q2_gen  = PaAlgo::Q2 (inMu_gen_TL, outMu_gen_TL); 
        xbj_gen = PaAlgo::xbj (inMu_gen_TL, outMu_gen_TL);
        nu_gen  = inMu_gen_TL.E() - outMu_gen_TL.E(); 
        W2_gen  = PaAlgo::W2 (inMu_gen_TL, outMu_gen_TL); 
        t_gen   = (targ_TL - proton_gen_TL) * (targ_TL - proton_gen_TL);  

        double E_gen  = inMu_gen_TL.E();

        phi_gg_gen = phiRV(inMu_gen_TL, outMu_gen_TL, proton_gen_TL, gamma_gen_TL, true);
        weight_PAMBH = Weight_Pam_BH(xbj_gen, Q2_gen, phi_gg_gen, t_gen, E_gen, phase_fac); 

        tree_gen->Fill(); 
      }

    } // End loop over MC data 

    // Check for exclusive event topology in the case of LEPTO
    enum class LujetCheckStatus {NOT_MC, I};
    if (leptoMC && !exclLepto(e)) {
        return;
    }

    //*******************************************  
    eventFlags.setFlagByName("allEvts_flag", true); 
    
    // Loop over reconstructed vertices in the event - REAL and MC RECONSTRUCTED DATA 
		for (int iv = 0; iv < e.NVertex(); iv++) { // begin loop over vertices
			//******************************************* 
			// Store info about primary vertex (if found) 
			const PaVertex &v = e.vVertex(iv); 
			if (!v.IsPrimary()) continue; // skip any vertices that are not primary 
      eventFlags.setFlagByName("pVtx_flag", true);
      pVtx_vec.SetXYZ(v.Pos(0), v.Pos(1), v.Pos(2));
      Chi2  = v.Chi2(); 
      Nprim = v.NOutParticles(); // number of tracks in vertex

			//*******************************************
    	// Store info about incoming muon beam (inMu)
      TVector3 SM2center(0, 0, 1825);
      const TVector3 SM2field = PaSetup::Ref().MagField(SM2center);
      Q_beam = SM2field(1) < 0 ? 1 : -1;

      static PaParticle beam; 
      static PaTrack beam_track; 
      static PaTPar par_beam;

			static BeamFluxParams beamParams; // Create an instance of BeamFluxParams
      // Will make cut on TiS later, just setting the event counter flag here 
	 		bool inMu_flag = beamFluxCheck(e, v, iv, Run, beamParams, beam, beam_track, par_beam, eventFlags);
      if (!inMu_flag) continue; // ignore all events that dont satisfy the beam flux requirements

      // TiS flag is used only to check DVCS statistics
      if (TiS_flag) {  
        eventFlags.setFlagByName("timeInSpill_flag", true); 
      }  
      
      //*******************************************
      // Store info about scattered muon (outMu)
      static PaParticle outMu; 
      static PaTrack outMu_track; 
      static PaTPar par_outMu; 

      static OutMuParams outMuParams; // Create an instance of OutMuParams 

      // outMu_flag will be true as long as any muon is found ( and all other conditions satisfied), even if it doesnt pass the full Hodoscope check 
      bool outMu_flag = outMuCheck(e, v, iv, Run, beam, HodoHelper, trig_flag, TiS_flag, outMuParams, outMu, outMu_track, par_outMu, eventFlags); 
      if (!outMu_flag) continue; // ignore events that dont satisfy requirements for scattered muon

      //*******************************************
      // Kinematic variables ... (1/2)
      inMu_TL  = par_beam.LzVec(M_mu); 
      outMu_TL = par_outMu.LzVec(M_mu);  

      Q2  = PaAlgo::Q2 (inMu_TL, outMu_TL); //Q2 = -(inMu_TL - outMu_TL).M2();
      y   = (inMu_TL.E() - outMu_TL.E()) / inMu_TL.E();
      nu  = (inMu_TL.E() - outMu_TL.E());
      W2  = PaAlgo::W2 (inMu_TL, outMu_TL);
      xbj = PaAlgo::xbj (inMu_TL, outMu_TL); //xbj = Q2/(2*q*targ_TL);

      //*******************************************
      // Exclusive selection starts here  ...  
      // Only one outgoing particle (scattered proton and photon are detected using ECals and CAMERA so will not be found here)
      if (Nprim != 1) continue; 
      eventFlags.setFlagByName("singleTrack_flag", true);

      // Current kinematic cuts will be tightened after the kinematically constrained fit is applied 
      if (Q2 < 0.5) continue; // inclusive Q2 cut
      eventFlags.setFlagByName("Q2_DIS_flag", true);   

      //if (y < 0.01 || y > 0.99) continue; // inclusive y cut
      //eventFlags.setFlagByName("y_DIS_flag", true);

      //*******************************************
      // Store info about real photon - check ECals for single neutral cluster (can be from any exclusive event (DVCS or BH))
      // ! CURRENTLY LOOKING FOR ANY NEUTRAL CLUSTERS, NOT SPECIFICALLY FOR ONE !
      int ecal0id = PaSetup::Ref().iCalorim("EC00P1__");
      int ecal1id = PaSetup::Ref().iCalorim("EC01P1__");
      int ecal2id = PaSetup::Ref().iCalorim("EC02P1__");

      std::vector<int> cl_id; // Store cluster index for valid exclusive clusters
      std::vector<int> pi0_cl_id; // Store cluster index for low-energy clusters 

      float EC0_thr = 4; // ECAL thresholds for exclusive event selection
      float EC1_thr = 5;
      float EC2_thr = 10;

/*       float EC0Low_thr = 0.5; // ECAL thresholds for visible pi0 identification (anything below this is noise)
      float EC1Low_thr = 0.63;  */

      for (int iclus = 0; iclus < e.NCaloClus(); iclus++) { // Begin loop over clusters 
          eventFlags.setFlagByName("clAll_flag", true);
          const PaCaloClus & cl = e.vCaloClus(iclus);
          int icalo = cl.iCalorim();

          if (cl.iTrack() != -1) continue; // Require neutral
          eventFlags.setFlagByName("clNeutral_flag", true);

          if (!EcalTimeCut(beam_track, cl, Run) && !e.IsMC()) continue;
          eventFlags.setFlagByName("clTime_flag", true);

          // Remove any clusters which do not have the required ECal ID 
          if ((icalo != ecal0id) && (icalo != ecal1id) && (icalo != ecal2id)) continue;

          float E = cl.E();
          if ((icalo == ecal0id && E >= EC0_thr) ||
              (icalo == ecal1id && E >= EC1_thr) ||
              (icalo == ecal2id && E >= EC2_thr)) {
            cl_id.push_back(iclus);
          } else {
            pi0_cl_id.push_back(iclus);
          }

      } // End loop over clusters 

      if (cl_id.size() != 0) {
        eventFlags.setFlagByName("nCls_flag", true);
      }

      if (pi0_cl_id.size() != 0) {
        eventFlags.setFlagByName("lowECl_flag", true);
      }

      // ! NOW CHECK THAT THERE IS ONLY ONE NEUTRAL CLUSTER !
      if (cl_id.size() != 1) continue; // skip event if there is more than one exclusive cluster
      eventFlags.setFlagByName("singleCl_flag", true);

      // Calculate TLorentz vectors for single high energy photon
      buildClusterVecs(e, v, cl_id[0], gamma_TL, cluster_TL);
      excl_calo = e.vCaloClus(cl_id[0]).iCalorim(); // get the ECal ID for the exclusive cluster

      //*******************************************
      // Store info about scattered proton candidates (check CAMERA)
      // ! CURRENTLY LOOKING FOR ANY PROTONS, NOT SPECIFICALLY FOR ONE !
      vector <CameraProton> proton_candidates = cam_inst->GetGoodCandidates(v); // vector holding all proton candidates

      TVector3 R_vtx;
      R_vtx.SetXYZ(v.Pos(0),v.Pos(1),v.Pos(2));

      int nCombs = 0;
      TLorentzVector p_miss_TL = targ_TL + inMu_TL - outMu_TL - gamma_TL;
      double pt_miss = p_miss_TL.Pt(); 
      double phi_miss = p_miss_TL.Phi();

      for (auto proton: proton_candidates) { // Begin loop over proton candidates
        eventFlags.setFlagByName("protonAll_flag", true); 
        bool DeltaPtPassed  = false; 
        bool DeltaPhiPassed = false; 
        bool DeltaZPassed   = false; 
        bool DeltaM2xPassed = false; 
        bool tPassed        = false; 
        bool betaPassed     = false; 

        double beta = proton.beta;
        if (beta >= 0.1 && beta <= 1) {betaPassed = true;}
        if (!betaPassed) continue;
        eventFlags.setFlagByName("proton_flag", true);

        p_camera_TL = proton.p4; 
        if (p_camera_TL.Mag() == 0) continue; // ignore events where there is no TL vector

        t = (targ_TL - p_camera_TL) * (targ_TL - p_camera_TL); 
        //double E_miss  = nu - gamma_TL.E() + t/(2 * M_p); // Missing energy (assuming proton)

        //*******************************************
        // Check that all combintations of the vertex, photon and proton satisfy exclusivity conditions 
        // 4 exclusivity variables: detla_phi, delta_pt, delta_Z, M2x
        double pt_camera = p_camera_TL.Pt(); 
        delta_pt = pt_camera - pt_miss; // transverse momentum
        double phi_camera = p_camera_TL.Phi();
        delta_phi = phi_camera - phi_miss; // azimuthal angle

        posRingA_vec = proton.Ahit.vec;
			  posRingB_vec = proton.Bhit.vec;
        double Z_inter = cam_inst->GetZA(R_vtx, posRingA_vec, posRingB_vec, phi_camera);
        delta_Z = posRingA_vec.Z() - Z_inter; // z position of the hits in the inner CAMERA ring

        M2x = (p_miss_TL - p_camera_TL) * (p_miss_TL - p_camera_TL);

        int protonPassCount = 0;
        if ((DeltaPtPassed  = std::fabs(delta_pt) <= 0.3)) protonPassCount++;
        if ((DeltaZPassed   = std::fabs(delta_Z) <= 16)) protonPassCount++;
        if ((DeltaM2xPassed = std::fabs(M2x) <= 0.3)) protonPassCount++;
        if ((tPassed = (t < -0.08 && t > -0.64))) protonPassCount++;

        if (delta_phi >= -0.4 && delta_phi <= 0.4) {
            DeltaPhiPassed = true;
            protonPassCount++;
        } else if ((delta_phi + 2 * TMath::Pi()) >= -0.4 && (delta_phi + 2 * TMath::Pi()) <= 0.4) {
            delta_phi += 2 * TMath::Pi();
            DeltaPhiPassed = true;
            protonPassCount++;
        } else if ((delta_phi - 2 * TMath::Pi()) >= -0.4 && (delta_phi - 2 * TMath::Pi()) <= 0.4) {
            delta_phi -= 2 * TMath::Pi();
            DeltaPhiPassed = true;
            protonPassCount++;
        }

        if (TiS_flag && (eventFlags.getFlag("passHodo_flag"))) {

          if (Q2 >= 1 && Q2 <= 10) {
            eventFlags.setFlagByName("Q2_DVCS_flag", true);

            if (y >= 0.05 && y <= 0.9) {
              eventFlags.setFlagByName("y_DVCS_flag", true);

              if (DeltaPtPassed) {
                eventFlags.setFlagByName("delta_pt_flag", true); 

                if (DeltaPhiPassed) {
                  eventFlags.setFlagByName("delta_phi_flag", true);

                  if (DeltaZPassed) {
                    eventFlags.setFlagByName("delta_Z_flag", true);

                    if (DeltaM2xPassed) {
                      eventFlags.setFlagByName("M2x_flag", true);
                    }

                  }

                }

              }

            }

          }

        }

        //*******************************************
        // Perform the kinematic fit for current combination and save results 
        static Fitter* FitInterface = &(Fitter::GetInstance());
        FitInterface->Init(R_vtx, beam_track, outMu_track, p_camera_TL, posRingA_vec, posRingB_vec); 
        FitInterface->Add_Photon(e.vCaloClus(cl_id[0]));
        FitInterface->SetupFit();
        FitInterface->DoFit(0, 1000);

        pVtxFit_vec     = *(FitInterface->GetVertex()->getCurr3Vec());
        posRingAFit_vec = *(FitInterface->GetHitA()->getCurr3Vec());
        posRingBFit_vec = *(FitInterface->GetHitB()->getCurr3Vec()); 
        clusterFit_vec  = *(FitInterface->GetOutPhotons()[0]->getCurr3Vec());

        inMuFit_TL   = *(FitInterface->GetMuonIn()->getCurr4Vec());
        outMuFit_TL  = *(FitInterface->GetMuonOut()->getCurr4Vec());
        protonFit_TL = *(FitInterface->GetProtonOut()->getCurr4Vec());
        targetFit_TL = *(FitInterface->GetProtonTarget()->getCurr4Vec());
        gammaFit_TL  = *(FitInterface->GetOutPhotons()[0]->getCurr4Vec());

        inMu_sigmaX  = TMath::Sqrt((*(FitInterface->GetMuonIn()->getCovMatrixDeltaY()))(0, 0));
        inMu_sigmaY  = TMath::Sqrt((*(FitInterface->GetMuonIn()->getCovMatrixDeltaY()))(1, 1));
        inMu_sigmaPx = TMath::Sqrt((*(FitInterface->GetMuonIn()->getCovMatrixDeltaY()))(2, 2));
        inMu_sigmaPy = TMath::Sqrt((*(FitInterface->GetMuonIn()->getCovMatrixDeltaY()))(3, 3));
        inMu_sigmaPz = TMath::Sqrt((*(FitInterface->GetMuonIn()->getCovMatrixDeltaY()))(4, 4));

        outMu_sigmaX  = TMath::Sqrt((*(FitInterface->GetMuonOut()->getCovMatrixDeltaY()))(0, 0));
        outMu_sigmaY  = TMath::Sqrt((*(FitInterface->GetMuonOut()->getCovMatrixDeltaY()))(1, 1));
        outMu_sigmaPx = TMath::Sqrt((*(FitInterface->GetMuonOut()->getCovMatrixDeltaY()))(2, 2));
        outMu_sigmaPy = TMath::Sqrt((*(FitInterface->GetMuonOut()->getCovMatrixDeltaY()))(3, 3));
        outMu_sigmaPz = TMath::Sqrt((*(FitInterface->GetMuonOut()->getCovMatrixDeltaY()))(4, 4));

        gamma_sigmaX = TMath::Sqrt((*(FitInterface->GetOutPhotons()[0]->getCovMatrixDeltaY()))(0, 0));
        gamma_sigmaY = TMath::Sqrt((*(FitInterface->GetOutPhotons()[0]->getCovMatrixDeltaY()))(1, 1));
        gamma_sigmaE = TMath::Sqrt((*(FitInterface->GetOutPhotons()[0]->getCovMatrixDeltaY()))(2, 2));

        proton_sigmaP     = TMath::Sqrt((*(FitInterface->GetProtonOut()->getCovMatrixDeltaY()))(0, 0));
        proton_sigmaTheta = TMath::Sqrt((*(FitInterface->GetProtonOut()->getCovMatrixDeltaY()))(1, 1));
        proton_sigmaPhi   = TMath::Sqrt((*(FitInterface->GetProtonOut()->getCovMatrixDeltaY()))(2, 2));

        ringA_sigmaR   = TMath::Sqrt((*(FitInterface->GetHitA()->getCovMatrixDeltaY()))(0, 0));
        ringA_sigmaPhi = TMath::Sqrt((*(FitInterface->GetHitA()->getCovMatrixDeltaY()))(1, 1));
        ringA_sigmaZ   = TMath::Sqrt((*(FitInterface->GetHitA()->getCovMatrixDeltaY()))(2, 2));

        ringB_sigmaR   = TMath::Sqrt((*(FitInterface->GetHitA()->getCovMatrixDeltaY()))(0, 0));
        ringB_sigmaPhi = TMath::Sqrt((*(FitInterface->GetHitA()->getCovMatrixDeltaY()))(1, 1));
        ringB_sigmaZ   = TMath::Sqrt((*(FitInterface->GetHitA()->getCovMatrixDeltaY()))(2, 2));
         
        fit_conv = FitInterface->GetFitOutput(chi2_fit, ndf_fit);  

        Q2_fit = PaAlgo::Q2 (inMuFit_TL, outMuFit_TL); 
        y_fit  = (inMuFit_TL.E() - outMuFit_TL.E()) / inMuFit_TL.E(); 
        nu_fit = (inMuFit_TL.E() - outMuFit_TL.E()); 
        t_fit  = (targetFit_TL - protonFit_TL) * (targetFit_TL - protonFit_TL);  

        bool Q2_cut = false; 
        bool y_cut  = false; 
        bool t_cut  = false; 
        bool nu_cut = false; 

        if ((Q2_fit > 1 && Q2_fit < 10) || std::isnan(Q2_fit)) {
          Q2_cut = true;
          eventFlags.setFlagByName("Q2Fit_flag", true);
        }
        if ((y_fit > 0.05 && y_fit < 0.95) || std::isnan(y_fit)) {
          y_cut = true; 
          eventFlags.setFlagByName("yFit_flag", true);
        }
        if ((t_fit < -0.08 && t_fit > -0.64) || std::isnan(t_fit)) {
          t_cut = true;
          eventFlags.setFlagByName("tFit_flag", true);
          std::cout << std::endl << "DEBUG:: " << Evt << "," << t_fit << std::endl;  
        }
/*         if (t_fit >= -0.08 || t_fit <= -0.64) {  
          t_cut = false; 
        } else {
          t_cut = true; 
          eventFlags.setFlagByName("tFit_flag", true);
        } */
        if ((nu_fit > 10 && nu_fit < 144) || std::isnan(nu_fit)) { 
          nu_cut = true; 
          eventFlags.setFlagByName("nuFit_flag", true);
        }
        if (Q2_cut && y_cut && t_cut && nu_cut) {
          eventFlags.setFlagByName("kinFitAll_flag", true);
        }

        nCombs++; // increment the counter 

        if (protonPassCount >= 4) { // Begin loop over candidats to save 
          eventFlags.setFlagByName("nExclCombo_flag", true);
          save_evt = true;
          if (!pi0_cl_id.empty()) { 
            for (auto iLow = std::size_t{0}; iLow < pi0_cl_id.size(); ++iLow) { // Begin loop over low-energy clusters 
              eventFlags.setFlagByName("nExclComboPi0_flag", true);
              const auto& cl_LowE = e.vCaloClus(pi0_cl_id[iLow]);
              low_calo = cl_LowE.iCalorim();

              // Build low-energy photon TLorentzVectors using function
              buildClusterVecs(e, v, pi0_cl_id[iLow], gammaLow_TL, clusterLow_TL);

              if (e.IsMC()) tree_MC->Fill();
              else tree->Fill();
            }
          } else {
            low_calo = -999;

            if (e.IsMC()) tree_MC->Fill();
            else tree->Fill();
          } // End loop over low-energy clusters

        } // End loop over candidats to save

      } // End loop over proton candidates

      //*******************************************
      // Debug statements ...
      printDebug("     ");
      printDebug("*** Run: " + std::to_string(Run) + ", spill: " + std::to_string(Spill) + ", event: " + std::to_string(EvtInSpill) + " ***");
      printDebug("    Vertex: (" + std::to_string(pVtx_vec.X()) + ", " + std::to_string(pVtx_vec.Y()) + ", " + std::to_string(pVtx_vec.Z()) + ")");
      printDebug("    mu: P: " + std::to_string(inMu_TL.P()) + " GeV/c, Charge: " + std::to_string(beam.Q()));
      printDebug("    mu': P: " + std::to_string(outMu_TL.P()) + " GeV/c, Charge: " + std::to_string(outMu.Q()));
      printDebug("    Kinematics: Q2: " + std::to_string(Q2) +  " GeV2, y: " + std::to_string(y) + ", W2: " + std::to_string(W2) + " GeV2, x: " + std::to_string(xbj));

      //*******************************************
		} // End loop over vertices 

  // Increment all counters whose flags are "true"
  eventFlags.incrementCounters(); 
  //Save the event
  if (save_evt) { 
    e.TagToSave();
  }   

} // End event loop 

void UserJobEnd970() {
  eventFlags.printFlags(); // Print to output stream 
}

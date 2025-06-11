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
bool verbose_mode      = false; // Create an instance of verbose_mode
bool plotECalEnergy    = false; // If this is true, DVCS threshold cut is turned off for ECals
bool kinfitPerformance = true;  // If this is true, additional data it stored in the tree for real data
bool leptoMC           = false; // Set to true for LEPTO MC to remove events with exclusive event topology 

//*****************************************************************************
void UserEvent970(PaEvent & e) { // begin event loop

    // Define constants
    static PaHodoHelper* HodoHelper = NULL; 
    static TiSRange*     tis_range  = NULL; 
    static PaCamera*     cam_inst   = NULL;

    const double M_pi    = G3partMass[8];  // Pion mass 
    const double M_gamma = G3partMass[1];  // Photon mass
    const double M_mu    = G3partMass[5];  // Muon mass
    const double M_p     = G3partMass[14]; // Proton mass 

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
    static int    Year;         // year data was taken 
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
    static TLorentzVector p_camera_TL;           // energy-momentum four vector of the proton measured by camera

    static int    pi0_calo;    // calorimeter ID for low energy clusters used in pi0 reconstruction -> same calorimeter as DVCS photon
    static double M_pi0;       // pi0 invariant mass 
    static double E_gammaLow;  // low-energy photon energy  

    static double y;   // fractional energy loss of the incoming lepton 
    static double nu;  // energy of the virtual photon 
    static double Q2;  // four-momentum transfer sqaured (Q is the four-momentum transferred between the incoming and outgoing lepton)
    static double W2;  // effective mass of the final state hadron system squared 
    static double xbj; // measure of the elasticity of the scattering process 
    static double t;   // four-momentum transfer squared of the nucleon 

    // Without using CAMERA information, the exclusive (DVCS) process can be selected through the 
    // detection of only the incident and outgoing muons and the photon, and using a cut on 
    // the missing energy or on the missing mass of a particle that is assumed to be a proton
    static double E_miss; 
    static double M2_miss; 

    // Exclusivity variables 
    static double delta_phi; 
    static double delta_pt;
    static double delta_Z;
    static double M2x;

    // Kinematic fit;
    static TVector3 pVtxFit_vec; 
    static TVector3 posRingAFit_vec; 
    static TVector3 posRingBFit_vec; 

    static TLorentzVector inMuFit_TL; 
    static TLorentzVector outMuFit_TL; 
    static TLorentzVector targetFit_TL;  
    static TLorentzVector gammaFit_TL; 
    static TLorentzVector protonFit_TL;  

    static TMatrixD Cov_inMu;   // covariance matrix for the beam muon  
    static TMatrixD Cov_outMu;  // covariance matrix for the scattered muon
    static TMatrixD Cov_gamma;  // covariance matrix for the photon
    static TMatrixD Cov_proton; // covariance matrix for the proton measured by camera 
    static TMatrixD Cov_ringA;  // covariance matrix for CAMERA ring A hit 
    static TMatrixD Cov_ringB;  // covariance matrix for CAMERA ring B hit

    static double chi2_fit; // chi2 of the fit  
    static int    ndf_fit;  // ndf of the fit 

    //*******************************************
    // HEPGEN data (prior to event selection) -> used for acceptance study 
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
    bool trig_flag = false; 
    bool TiS_flag  = false;
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
        tree->Branch("Evt", &Evt, "Evt/I");
        tree->Branch("Chi2", &Chi2, "Chi2/F");
        tree->Branch("Spill", &Spill, "Spill/I");
        tree->Branch("Q_beam", &Q_beam, "Q_beam/I");
        // Particle/vertex vectors 
        tree->Branch("inMu_TL", &inMu_TL);
        tree->Branch("outMu_TL", &outMu_TL);
        tree->Branch("gamma_TL", &gamma_TL);
        tree->Branch("pVtx_vec", &pVtx_vec);
        tree->Branch("p_camera_TL", &p_camera_TL);
        tree->Branch("posRingA_vec", &posRingA_vec);
        tree->Branch("posRingB_vec", &posRingB_vec);
        // Kinematic variables 
        tree->Branch("y", &y, "y/D");
        tree->Branch("t", &t, "t/D");
        tree->Branch("nu", &nu, "nu/D");
        tree->Branch("Q2", &Q2, "Q2/D");
        tree->Branch("W2", &W2, "W2/D");
        tree->Branch("xbj", &xbj, "xbj/D");
        // Missing energy/mass (assuming missing proton)
        tree->Branch("E_miss", &E_miss, "E_miss/D");
        tree->Branch("M2_miss", &M2_miss, "M2_miss/D");
        // Invariant mass of visible pi0
        tree->Branch("M_pi0", &M_pi0, "M_pi0/D");
        tree->Branch("pi0_calo", &pi0_calo, "pi0_calo/I");
        tree->Branch("E_gammaLow", &E_gammaLow, "E_gammaLow/D");  
        // Exclusivity variables      
        tree->Branch("M2x", &M2x, "M2x/D");
        tree->Branch("delta_Z", &delta_Z, "delta_Z/D");
        tree->Branch("delta_pt", &delta_pt, "delta_pt/D");
        tree->Branch("delta_phi", &delta_phi, "delta_phi/D");
        // Kinematic fit vectors 
        tree->Branch("inMuFit_TL", &inMuFit_TL);
        tree->Branch("outMuFit_TL", &outMuFit_TL);
        tree->Branch("gammaFit_TL", &gammaFit_TL);
        tree->Branch("pVtxFit_vec", &pVtxFit_vec);
        tree->Branch("protonFit_TL", &protonFit_TL);
        tree->Branch("posRingAFit_vec", &posRingAFit_vec);
        tree->Branch("posRingBFit_vec", &posRingBFit_vec);
        // Kinematic fit covariance matrices
        tree->Branch("Cov_inMu", &Cov_inMu);
        tree->Branch("Cov_outMu", &Cov_outMu);
        tree->Branch("Cov_gamma", &Cov_gamma);
        tree->Branch("Cov_ringA", &Cov_ringA);
        tree->Branch("Cov_ringB", &Cov_ringB);
        tree->Branch("Cov_proton", &Cov_proton);
        tree->Branch("ndf_fit", &ndf_fit, "ndf_fit/I");
        tree->Branch("chi2_fit", &chi2_fit, "chi2_fit/D");
        tree->Branch("fit_conv", &fit_conv, "fit_conv/O");

        //*******************************************
        tree_MC = new TTree("USR970_MC","User 970 DVCS NTuple - MC Data"); // name (has to be unique) and title of the Ntuple
        // Basic event information
        tree_MC->Branch("Run", &Run, "Run/I");
        tree_MC->Branch("Evt", &Evt, "Evt/I");
        tree_MC->Branch("Chi2", &Chi2, "Chi2/F");
        tree_MC->Branch("Spill", &Spill, "Spill/I");
        tree_MC->Branch("Q_beam", &Q_beam, "Q_beam/I");
        // Particle/vertex vectors 
        tree_MC->Branch("inMu_TL", &inMu_TL);
        tree_MC->Branch("outMu_TL", &outMu_TL);
        tree_MC->Branch("gamma_TL", &gamma_TL);
        tree_MC->Branch("pVtx_vec", &pVtx_vec);
        tree_MC->Branch("p_camera_TL", &p_camera_TL);
        tree_MC->Branch("posRingA_vec", &posRingA_vec);
        tree_MC->Branch("posRingB_vec", &posRingB_vec);
        // Kinematic variables 
        tree_MC->Branch("y", &y, "y/D");
        tree_MC->Branch("t", &t, "t/D");
        tree_MC->Branch("nu", &nu, "nu/D");
        tree_MC->Branch("Q2", &Q2, "Q2/D");
        tree_MC->Branch("W2", &W2, "W2/D");
        tree_MC->Branch("xbj", &xbj, "xbj/D");
        // Missing energy/mass (assuming missing proton)
        tree_MC->Branch("E_miss", &E_miss, "E_miss/D");
        tree_MC->Branch("M2_miss", &M2_miss, "M2_miss/D");
        // Invariant mass of visible pi0
        tree_MC->Branch("M_pi0", &M_pi0, "M_pi0/D");
        tree_MC->Branch("pi0_calo", &pi0_calo, "pi0_calo/I");
        tree_MC->Branch("E_gammaLow", &E_gammaLow, "E_gammaLow/D");  
        // Exclusivity variables      
        tree_MC->Branch("M2x", &M2x, "M2x/D");
        tree_MC->Branch("delta_Z", &delta_Z, "delta_Z/D");
        tree_MC->Branch("delta_pt", &delta_pt, "delta_pt/D");
        tree_MC->Branch("delta_phi", &delta_phi, "delta_phi/D");
        // Kinematic fit vectors 
        tree_MC->Branch("inMuFit_TL", &inMuFit_TL);
        tree_MC->Branch("outMuFit_TL", &outMuFit_TL);
        tree_MC->Branch("gammaFit_TL", &gammaFit_TL);
        tree_MC->Branch("pVtxFit_vec", &pVtxFit_vec);
        tree_MC->Branch("protonFit_TL", &protonFit_TL);
        tree_MC->Branch("posRingAFit_vec", &posRingAFit_vec);
        tree_MC->Branch("posRingBFit_vec", &posRingBFit_vec);
        // Kinematic fit covariance matrices
        tree_MC->Branch("Cov_inMu", &Cov_inMu);
        tree_MC->Branch("Cov_outMu", &Cov_outMu);
        tree_MC->Branch("Cov_gamma", &Cov_gamma);
        tree_MC->Branch("Cov_ringA", &Cov_ringA);
        tree_MC->Branch("Cov_ringB", &Cov_ringB);
        tree_MC->Branch("Cov_proton", &Cov_proton);
        tree_MC->Branch("ndf_fit", &ndf_fit, "ndf_fit/I");
        tree_MC->Branch("chi2_fit", &chi2_fit, "chi2_fit/D");
        tree_MC->Branch("fit_conv", &fit_conv, "fit_conv/O");

        //*******************************************
        tree_gen = new TTree("USR970_GEN","User 970 HEPGEN NTuple");
        // Basic event information
        tree_gen->Branch("Run", &Run, "Run/I");
        tree_gen->Branch("Evt", &Evt, "Evt/I");
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
    }
    if (trig_mask & LT) {
        trigCheck += "LT ";
        trig_flag = true;
    }
    if (trig_mask & OT) {
        trigCheck += "OT ";
        trig_flag = true;
    }
    if (trig_mask & LAST) {
        trigCheck += "LAST ";
        trig_flag = true;
    }

    //*******************************************
    // Initialize variables, extra flags and check time in spill (time in spill cut is applied later not here)   
    // Check flags are only used to check statistics not for final event selection 
    
    eventFlags.createFlag("Q2_DVCS_flag", "No. of events where 1 < Q2 < 10 (for stats only)");
    eventFlags.createFlag("y_DVCS_flag", "No. of events where 0.05 < y < 0.9 (for stats only)");
    eventFlags.createFlag("singleTrack_flag", "No. of events where primary vertex only has one outgoing track");
    eventFlags.createFlag("DVCS_Check", "STATS* No. of events which satisfy beam, scattered muon and DVCS kin cuts");
    eventFlags.createFlag("Q2_DIS_flag", "No. of events where Q2 > 0.5 (DIS)");
    //eventFlags.createFlag("y_DIS_flag", "No. of events where 0.01 < y < 0.99");

    eventFlags.createFlag("clAll_flag", "No. of events that have any clusters");
    eventFlags.createFlag("clNeutral_flag", "No. of events where clusters are not associated with charged tracks");
    eventFlags.createFlag("clTime_flag", "No. of events where cluster timing is within requirements");
    eventFlags.createFlag("nCls_flag", "No. of events where clusters are found in ECal 0, 1 or 2 only");
    //eventFlags.createFlag("lowECl_flag", "No. of events with low energy photons");
    eventFlags.createFlag("singleCl_flag", "No. of events where there is only a single high energy cluster in the ECals");

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
    
    eventFlags.createFlag("nCombo_flag", "No. of events where a vertex, photon and proton combination satisfies exclusivity conditions");
     
    eventFlags.resetFlags(); // Reset all event statistic counter flags to false

    Run         = e.RunNum();
    Evt         = e.UniqueEvNum();
    Year        = e.Year();
    EvtInSpill  = e.EvInSpill(); 
    Spill       = e.SpillNum(); 
    TimeInSpill = e.TimeInSpill();  

    if (Run != LastRun) { // Reinitialize HodoHelper and tis_range only if the run number changes 
        HodoHelper = & PaHodoHelper::Init("", true);  
        //tis_range  = new TiSRange("/Users/gursimran/cern/phastPackages/flux_files/flux_Johannes/2016/flux_files");
        tis_range  = new TiSRange("/afs/cern.ch/user/g/gkainth/phastPackages/flux_files/flux_Johannes/2016/flux_files");
        //Set_TiSrange("/afs/cern.ch/user/g/gkainth/phastPackages/flux_files/flux_Johannes/2016/flux_files", Run, Run);
        LastRun = Run;  // Update LastRun to the current run number
    }
 
    TiS_flag = tis_range->CheckWindow(Run, Spill, TimeInSpill); // check the time in spill 
    cam_inst= & PaCamera::GetInstance();
    cam_inst->NewEvent(e);

    //*******************************************
    // Get MC PAM weights for HEPGEN and save to tree_gen 
    if (e.IsMC()) { // Begin loop over MC data 
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

        const PaMCtrack & t_beam   = e.vMCtrack(v.iMCtrack(0));
			  const PaMCtrack & t_omu    = e.vMCtrack(v.iMCtrack(1));
			  const PaMCtrack & t_gamma  = e.vMCtrack(v.iMCtrack(2));
			  const PaMCtrack & t_proton = e.vMCtrack(v.iMCtrack(3));

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
    bool exclEvt = exclLepto(e, leptoMC); 
    if (exclEvt) return;  

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
      if (!eventFlags.getFlag("pVtx_flag")) continue; // primary vertices only 
      //if (!eventFlags.getFlag("inMuTrack_flag")) continue; // beam has track with parameters 

      // Checks for Flux statistics - not final event seletion 
      if (eventFlags.getFlag("zFirst_flag") && eventFlags.getFlag("momRange_flag") && eventFlags.getFlag("momErr_flag") && 
          eventFlags.getFlag("BMS_flag") && eventFlags.getFlag("FI_flag") && eventFlags.getFlag("SI_flag") && 
          eventFlags.getFlag("crossCells_flag") && eventFlags.getFlag("meantime_flag")) {
          eventFlags.setFlagByName("flux_Check", true);  
      }
      if (eventFlags.getFlag("flux_Check") && TiS_flag) {
        eventFlags.setFlagByName("fluxTiS_Check", true); 
      }

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
/*       if (eventFlags.getFlag("outMuFound_flag") != 1) continue; // scattered muon is found
      if (eventFlags.getFlag("vtxInTarget_flag") != 1) continue; // vertex must be in the target 
      if (eventFlags.getFlag("zFirstLast_flag") != 1) continue;  */

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

      if (eventFlags.getFlag("passHodo_flag") == 1) {
        if ( Q2 > 1 && Q2 < 10) {
          eventFlags.setFlagByName("Q2_DVCS_flag", true); // just used to check statistics - not an actual cut 
          if ( y > 0.05 && y < 0.9) {
            eventFlags.setFlagByName("y_DVCS_flag", true); // just used to check statistics - not an actual cut
          }
        }
      }

      //*******************************************
      // Exclusive selection starts here  ...  
      // Only one outgoing particle (scattered proton and photon are detected using ECals and CAMERA so will not be found here)
      if (Nprim != 1) continue; 
      if (eventFlags.getFlag("y_DVCS_flag") == 1) {
        eventFlags.setFlagByName("singleTrack_flag", true);
      }
      
      if ((eventFlags.getFlag("passHodoStat_flag") == 1) && (eventFlags.getFlag("Q2_DVCS_flag") == 1) && (eventFlags.getFlag("y_DVCS_flag") == 1)) { 
        eventFlags.setFlagByName("DVCS_Check", true);
      }

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

      float EC0Low_thr = 0.5; // ECAL thresholds for visible pi0 identification (anything below this is noise)
      float EC1Low_thr = 0.63; 

      int clusterCount = 0;
      for (int iclus = 0; iclus < e.NCaloClus(); iclus++) { // Begin loop over clusters 
          const PaCaloClus & cl = e.vCaloClus(iclus);
          int icalo = cl.iCalorim();

          eventFlags.setFlagByName("clAll_flag", true);
          if (cl.iTrack() != -1) continue; // Require neutral
          eventFlags.setFlagByName("clNeutral_flag", true);
          if (!EcalTimeCut(beam_track, cl, Run) && !e.IsMC()) continue;
          eventFlags.setFlagByName("clTime_flag", true);

          // Remove any clusters which do not have the required ECal ID 
          if ((icalo != ecal0id) && (icalo != ecal1id) && (icalo != ecal2id)) continue; 

          float E = cl.E();
          /*         if (plotECalEnergy) { // Fill the reconstructed energy distributions for each ECal
          if (icalo == ecal0id) h97_gamma_E_EC0->Fill(E);
          if (icalo == ecal1id) h97_gamma_E_EC1->Fill(E);
          if (icalo == ecal2id) h97_gamma_E_EC2->Fill(E);
        } */

          if (icalo == ecal0id) {
              if (E >= EC0_thr) {
                  cl_id.push_back(iclus);
                  clusterCount++;
              }
              else if (E >= EC0Low_thr) {
                  pi0_cl_id.push_back(iclus);
              }
          }
          else if (icalo == ecal1id) {
              if (E >= EC1_thr) {
                  cl_id.push_back(iclus);
                  clusterCount++;
              }
              else if (E >= EC1Low_thr) {
                  pi0_cl_id.push_back(iclus);
              }
          }
          else if (icalo == ecal2id) {
              if (E >= EC2_thr) {
                  cl_id.push_back(iclus);
                  clusterCount++;
              }
          }
      } // End loop over clusters 

      if (clusterCount == 0) continue; // skip event if there are no exclusive clusters
      eventFlags.setFlagByName("nCls_flag", true);

      if (pi0_cl_id.size() != 0) {
        eventFlags.setFlagByName("lowECl_flag", true);
      }

      // ! NOW CHECK THAT THERE IS ONLY ONE NEUTRAL CLUSTER !
      if (clusterCount != 1) continue; // skip event if there is more than one exclusive cluster
      eventFlags.setFlagByName("singleCl_flag", true);

     // Calculate TL vector for single high energy photon
      if (clusterCount != 0) { // Ensure there is at least one cluster
        const auto& cluster = e.vCaloClus(cl_id[0]); // Directly access the first (and only) cluster
        double phdz  = cluster.Z() - v.Z();
        double phdy  = cluster.Y() - v.Y();
        double phdx  = cluster.X() - v.X();
        double phL   = TMath::Sqrt(phdx * phdx + phdy * phdy + phdz * phdz);
        double cl_E  = cluster.E();
        gamma_TL.SetPxPyPzE(cl_E * phdx / phL, cl_E * phdy / phL, cl_E * phdz / phL, cl_E);
      }
 
      //*******************************************
      // Kinematic variables ... (2/2)
      E_miss  = nu - gamma_TL.E() + t/(2 * M_p); // Missing energy (assuming proton)

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
        double beta = proton.beta;
        if (beta < 0.1 || beta > 1) continue;

        p_camera_TL = proton.p4; 
        if (p_camera_TL.Mag() == 0) continue; // ignore events where there is no TL vector
        eventFlags.setFlagByName("proton_flag", true);

        //*******************************************
        // Check that all combintations of the vertex, photon and proton satisfy exclusivity conditions 
        // 4 exclusivity variables: detla_phi, delta_pt, delta_Z, M2x
        double pt_camera = p_camera_TL.Pt(); 
        delta_pt = pt_camera - pt_miss; // transverse momentum

        double phi_camera = p_camera_TL.Phi();
        delta_phi = phi_camera - phi_miss; // azimuthal angle

        TVector3 posRingA = proton.Ahit.vec;
			  TVector3 posRingB = proton.Bhit.vec;
        double Z_inter = cam_inst->GetZA(R_vtx, posRingA, posRingB, phi_camera);
        delta_Z = posRingA.Z() - Z_inter; // z position of the hits in the inner CAMERA ring

        M2x = (p_miss_TL - p_camera_TL) * (p_miss_TL - p_camera_TL);

        if (TiS_flag && (eventFlags.getFlag("passHodo_flag"))) {

          if (Q2 >= 1 && Q2 <= 10) {
            eventFlags.setFlagByName("Q2_DVCS_flag", true);

            if (y >= 0.05 && y <= 0.9) {
              eventFlags.setFlagByName("y_DVCS_flag", true);

              if (std::fabs(delta_pt) <= 0.3) {
                eventFlags.setFlagByName("delta_pt_flag", true);
                double delta_phi_norm = TVector2::Phi_mpi_pi(delta_phi);

                if (std::fabs(delta_phi_norm) <= 0.4) {
                  eventFlags.setFlagByName("delta_phi_flag", true);

                  if (std::fabs(delta_Z) <= 16) {
                    eventFlags.setFlagByName("delta_Z_flag", true);

                    if (std::fabs(M2x) <= 0.3) {
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
        FitInterface->Init(R_vtx, beam_track, outMu_track, p_camera_TL, posRingA, posRingB); 
        FitInterface->Add_Photon(e.vCaloClus(cl_id[0]));
        FitInterface->SetupFit();
        FitInterface->DoFit(0, 1000);

        pVtxFit_vec     = *(FitInterface->GetVertex()->getCurr3Vec());
        posRingAFit_vec = *(FitInterface->GetHitA()->getCurr3Vec());
        posRingBFit_vec = *(FitInterface->GetHitB()->getCurr3Vec());

        inMuFit_TL   = *(FitInterface->GetMuonIn()->getCurr4Vec());
        outMuFit_TL  = *(FitInterface->GetMuonOut()->getCurr4Vec());
        protonFit_TL = *(FitInterface->GetProtonOut()->getCurr4Vec());
        targetFit_TL = *(FitInterface->GetProtonTarget()->getCurr4Vec());
        gammaFit_TL  = *(FitInterface->GetOutPhotons()[0]->getCurr4Vec());

        Cov_inMu   = *(FitInterface->GetMuonIn()->getCovMatrixDeltaY());
        Cov_outMu  = *(FitInterface->GetMuonOut()->getCovMatrixDeltaY()); 
        Cov_gamma  = *(FitInterface->GetOutPhotons()[0]->getCovMatrixDeltaY());
        Cov_proton = *(FitInterface->GetProtonOut()->getCovMatrixDeltaY());
        Cov_ringA  = *(FitInterface->GetHitA()->getCovMatrixDeltaY());
        Cov_ringB  = *(FitInterface->GetHitB()->getCovMatrixDeltaY());
         
        fit_conv = FitInterface->GetFitOutput(chi2_fit, ndf_fit);  

        double Q2_fit = PaAlgo::Q2 (inMuFit_TL, outMuFit_TL); 
        double y_fit  = (inMuFit_TL.E() - outMuFit_TL.E()) / inMuFit_TL.E(); 
        double nu_fit = (inMuFit_TL.E() - outMuFit_TL.E()); 
        double t_fit  = (targetFit_TL - protonFit_TL) * (targetFit_TL - protonFit_TL);  
        t      = (targ_TL - p_camera_TL) * (targ_TL - p_camera_TL); 

        // ! NO CUT! Check for statistics only :D 
        if ((Q2_fit > 1 && Q2_fit < 10) || std::isnan(Q2_fit)) {
          eventFlags.setFlagByName("Q2Fit_flag", true);
        }
        if ((y_fit > 0.05 && y_fit < 0.95) || std::isnan(y_fit)) {
          eventFlags.setFlagByName("yFit_flag", true);
        }
        if ((fabs(t_fit) > 0.08 && fabs(t_fit) < 0.64) || std::isnan(t_fit)) {
          eventFlags.setFlagByName("tFit_flag", true);
        }
        if ((nu_fit > 10 && nu_fit < 144) || std::isnan(nu_fit)) { 
          eventFlags.setFlagByName("nuFit_flag", true);
        }

        if ((eventFlags.getFlag("Q2Fit_flag") == 1) && (eventFlags.getFlag("yFit_flag") == 1) && 
            (eventFlags.getFlag("tFit_flag") == 1) && (eventFlags.getFlag("nuFit_flag") == 1)) {
          eventFlags.setFlagByName("kinFitAll_flag", true);
        }

        nCombs++; // increment the counter 
        eventFlags.setFlagByName("nCombo_flag", true); 

      } // End loop over proton candidates

      if (nCombs > 0) { // only save events where there is at least one proton, vertex and exlusive photon combination 
        save_evt = true; 
      } 

/*       if (!((Q2 < 10 && Q2 > 1) && (y > 0.05 && y < 0.9) && (nu > 10 && nu < 144) && (t > -0.64 && t < -0.08))) continue; 
      for (auto proton: protons_excl) {
        TLorentzVector p_camera_TL = proton.p4;
        M2_miss = (inMu_TL + targ_TL - outMu_TL - gamma_TL - p_camera_TL)*(inMu_TL + targ_TL - outMu_TL - gamma_TL - p_camera_TL); 
      }    */

      // Find visible pi0 contamination in sample
/*       for (int iLow = 0; iLow < pi0_cl_id.size(); ++iLow) { // Begin loop over photon pairs 
        const auto& cl_LowE = e.vCaloClus(pi0_cl_id[iLow]);
        int DVCS_calo = e.vCaloClus(cl_id[0]).iCalorim(); 
        pi0_calo = cl_LowE.iCalorim();

        // Build low-energy photon TLorentzVector
        double phdz = cl_LowE.Z() - v.Z();
        double phdy = cl_LowE.Y() - v.Y();
        double phdx = cl_LowE.X() - v.X();
        double phL  = TMath::Sqrt(phdx * phdx + phdy * phdy + phdz * phdz);
        double cl_E = cl_LowE.E();

        TLorentzVector gammaLow_TL;
        gammaLow_TL.SetPxPyPzE(cl_E * phdx / phL, cl_E * phdy / phL, cl_E * phdz / phL, cl_E);

        // Combine to form a pi0 candidate
        TLorentzVector pi0Cand_TL = gamma_TL + gammaLow_TL;
        M_pi0 = pi0Cand_TL.M(); 
        E_gammaLow = gammaLow_TL.E(); 

      } // End loop over photon pairs  */

      //*******************************************
      // Fill histograms for ECal Energy after all exclusivity cuts
/*       if (plotECalEnergy) {
        if (clusterCount != 0) { // Ensure there is at least one cluster
          const auto& cluster = e.vCaloClus(cl_id.front()); // Directly access the first (and only) cluster
          if (cl_id.front() == ecal0id) {h97_gamma_E_EC0_ex->Fill(cluster.E());}
          if (cl_id.front() == ecal1id) {h97_gamma_E_EC1_ex->Fill(cluster.E());}
          if (cl_id.front() == ecal2id) {h97_gamma_E_EC2_ex->Fill(cluster.E());} 
        }
      } */

      //*******************************************
      // Debug statements ...
      printDebug("     ");
      printDebug("*** Run: " + std::to_string(Run) + ", spill: " + std::to_string(Spill) + ", event: " + std::to_string(EvtInSpill) + " ***");
      printDebug("    Vertex: (" + std::to_string(pVtx_vec.X()) + ", " + std::to_string(pVtx_vec.Y()) + ", " + std::to_string(pVtx_vec.Z()) + ")");
      printDebug("    mu: P: " + std::to_string(inMu_TL.P()) + " GeV/c, Charge: " + std::to_string(beam.Q()));
      printDebug("    mu': P: " + std::to_string(outMu_TL.P()) + " GeV/c, Charge: " + std::to_string(outMu.Q()));
      printDebug("    Kinematics: Q2: " + std::to_string(Q2) +  " GeV2, y: " + std::to_string(y) + ", W2: " + std::to_string(W2) + " GeV2, x: " + std::to_string(xbj));

      //*******************************************
      // Fill histrograms here if needed and save data to the tree(s)
      if (e.IsMC()) {
        tree_MC->Fill();
      } else {tree->Fill();}

		} // End loop over vertices 

  // Increment all counters whose flags are "true"
  eventFlags.incrementCounters(); 
  //Save the event
  if (save_evt) {e.TagToSave();}   

} // End event loop 

void UserJobEnd970() {
  if (plotECalEnergy) {
    std::cout << "! WARNING: ECAL ENERGY PLOT IS ON. TURN OFF TO IMPROVE EFFICIENCY. !"<< std::endl;
    std::cout << std::endl;
  }

  eventFlags.printFlags(); // Print to output stream 

}




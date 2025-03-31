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
#include "G3part.h" 
#include "PaHodoHelper.h"

#include "ecal_time_cuts.h"
#include "Tools_Camera.hh"

//#include "/Users/gursimran/cern/phastPackages/xcheck_newTiS/tis_range.cc"
//#include "/Users/gursimran/cern/phastPackages/xcheck_newTiS/tis_range.h"
//#include "UConn_Tools.h"

#include "/afs/cern.ch/user/g/gkainth/phastPackages/xcheck_newTiS/tis_range.cc"
#include "/afs/cern.ch/user/g/gkainth/phastPackages/xcheck_newTiS/tis_range.h"
#include "/afs/cern.ch/user/g/gkainth/phastPackages/git_COMPASS/UConn_Tools.h"
//#include "/afs/cern.ch/user/g/gkainth/phastPackages/flux_files/flux_Johannes/functions/TiS_range.cc"

// ************************************************************************** //
// UserEvent for preselecting exclusive photon events	(DVCS)			            //
// In the selection all possible combinations of:					                    //
// Vertices (incoming and outgoing muons), exclusive photon and recoil proton //
// ************************************************************************** //

extern "C" float prob_(float&, int&);

// Global flag for verbose mode 
bool verbose_mode = false; // Create an instance of verbose_mode

// Event statistic counter flags  
EventFlags eventFlags; // Create an instance of flags and counters 

// User selection flags 
bool plotECalEnergy = true; // if this is true DVCS threshold cut is turned off 

//*****************************************************************************
void UserEvent970(PaEvent & e) { // begin event loop

    // Define constants
    static PaHodoHelper* HodoHelper = NULL; 
    static TiSRange*     tis_range  = NULL; 
    static PaCamera*     cam_inst   = NULL;

    const double M_pi    = G3partMass[8]; // Pion mass 
    const double M_gamma = G3partMass[1]; // Photon mass
    const double M_mu    = G3partMass[5]; // Muon mass
    const double M_p     = G3partMass[14]; // Proton mass 

    // Declare all objects globally
    static TH1F* h97_Zprim  = NULL; 
    static TH1F* h97_Yprim  = NULL; 
    static TH1F* h97_Xprim  = NULL; 
    static TH2F* h97_XYprim = NULL;

    static TH1F* h97_inMu_p  = NULL;
    static TH1F* h97_inMu_py = NULL;
    static TH1F* h97_inMu_px = NULL; 

    static TH1F* h97_outMu_p  = NULL;
    static TH1F* h97_outMu_py = NULL;
    static TH1F* h97_outMu_px = NULL; 

    static TH1F* h97_y      = NULL;
    static TH1F* h97_nu     = NULL;
    static TH1F* h97_Q2     = NULL;
    static TH1F* h97_W2     = NULL;
    static TH1F* h97_xbj    = NULL;
    static TH2F* h97_Q2xbj  = NULL;
    static TH2F* h97_Chi2Q2 = NULL; 
    static TH1F* h97_t      = NULL;

    static TH1F* h97_gamma_E_EC0 = NULL;
    static TH1F* h97_gamma_E_EC1 = NULL;
    static TH1F* h97_gamma_E_EC2 = NULL;
    static TH1F* h97_E_miss      = NULL;
    static TH1F* h97_M2_miss     = NULL;

    static TH1F* h97_delta_phi = NULL; 
    static TH1F* h97_delta_pt  = NULL;
    static TH1F* h97_delta_Z   = NULL;
    static TH1F* h97_M2x       = NULL;

    static TH1F* h97_gamma_E_EC0_ex = NULL;
    static TH1F* h97_gamma_E_EC1_ex = NULL;
    static TH1F* h97_gamma_E_EC2_ex = NULL;

    static TTree* tree(NULL);

    //
    // Variables to be stored into analysis Ntuple
    //
    // (if you want to extend the list:
    // - add variable here,
    // - add in Ntuple definition below 
    // and do not forget to fill the variable in the code
    //
    static unsigned long long int Evt; // event number - unique evt number (based on run, spill and evt num in spill) 
    static int    Run;          // run number
    static int    LastRun = -1; // store the previous run number (used to reintialize Hodohelper, tis_range if there are multiple runs)
    static int    Year;         // year data was taken 
    static int    EvtInSpill;   // event number in spill 
    static int    Spill;        // spill number
    static double TimeInSpill;  // time in spill 
    static float  Zprim;        // Z coordinate of primary vertex (10000 if not found)
    static float  Yprim;        // Y coordinate of primary vertex (10000 if not found)
    static float  Xprim;        // X coordinate of primary vertex (10000 if not found)
    static float  Chi2;         // Chi2 of the reconstructed vertex 
    static int    Nprim;        // Number of tracks in primary vertex (-1 in fot found)
    static int    trig_mask;

    static double inMu_pz; // Z component of the beam muon momentum 
    static double inMu_py; // Y component of the beam muon momentum 
    static double inMu_px; // X component of the beam muon momentum
    static double inMu_p;  // Magnitude of the beam muon momentum 
    static double inMu_E;  // Magnitude of the beam muon energy
    
    static double outMu_pz; // Z component of the scattered muon momentum 
    static double outMu_py; // Y component of the scattered muon momentum
    static double outMu_px; // X component of the scattered muon momentum 
    static double outMu_p;  // Magnitude of the scattered muon momentum 
    static double outMu_E;  // Magnitude of the scattered muon energy

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

    // Event selection flags 
    bool trig_flag  = false; 
    bool TiS_flag   = false;
    bool flux_flag  = false;
    bool outMu_flag = false;  

    //*****************************************************************************
    static bool first(true);
    if (first) { // histograms and Ntuples booking block
        Phast::Ref().HistFileDir("UserEvent970");
    
        // 1D and 2D
        h97_Zprim  = new TH1F("h97_Zprim", "Primary Vertex Z (cm); Z [cm]; Events", 100, -250, 250);
        h97_Yprim  = new TH1F("h97_Yprim", "Primary Vertex Y (cm); Y [cm]; Events", 100, -3, 3);
        h97_Xprim  = new TH1F("h97_Xprim", "Primary Vertex X (cm); X [cm]; Events", 100, -3, 3);
        h97_XYprim = new TH2F("h97_XYprim", "Primary Vertex XY (cm); X [cm]; Y [cm]", 100, -3, 3, 100, -3, 3);

        h97_inMu_p  = new TH1F("h97_inMu_p", "P Incoming Muon (GeV/c); P [GeV/c]; Events", 100, 0, 200);
        h97_inMu_py = new TH1F("h97_inMu_py", "Py Incoming Muon (GeV/c); Py [GeV/c]; Events", 100, 0, 10);
        h97_inMu_px = new TH1F("h97_inMu_px", "Px Incoming Muon (GeV/c); Px [GeV/c]; Events", 100, 0, 10);

        h97_outMu_p  = new TH1F("h97_outMu_p", "P Scattered Muon (GeV/c); P [GeV/c]; Events", 100, 0, 200);
        h97_outMu_py = new TH1F("h97_outMu_py", "Py Scattered Muon (GeV/c); Py [GeV/c]; Events", 100, 0, 10);
        h97_outMu_px = new TH1F("h97_outMu_px", "Px Scattered Muon (GeV/c); Px [GeV/c]; Events", 100, 0, 10);

        h97_y      = new TH1F("h97_y", "Fractional Energy Loss of Incoming Muon (y); y; Events", 100, 0, 1);
        h97_nu     = new TH1F("h97_nu", "Energy of the virtual photon (#nu); #nu [GeV]; Events", 100, 0, 180);
        h97_Q2     = new TH1F("h97_Q2", "Four-momentum Transfer Squared (Lepton) (Q^{2}); Q^{2} [GeV^{2}]; Events", 100, 0, 10);
        h97_W2     = new TH1F("h97_W2", "Effective Mass of final state hadrons Squared (W^{2}); W^{2} [GeV^{2}]; Events", 100, 0, 350);
        h97_Q2xbj  = new TH2F("h97_Q2xbj", "Kinematic Coverage of Dataset; x_{bj}; Q^{2} [GeV^{2}]", 150, 0, 1, 150, 0, 100);
        h97_Chi2Q2 = new TH2F("h97_Chi2Q2", "Chi2 of Reconstruced Vertex vs. Q2; Chi^{2}; Q^{2} [GeV^{2}]", 100, 0, 10, 100, 0, 10);
        h97_t      = new TH1F("h97_t", "Four-momentum Transfer Squared (Nucleon) (t); t [GeV^{2}]; Events", 100, 0, 200);

        const int nBins = 100;      // Number of bins
        double xMin = 1e-3;         // Minimum x value (avoid 0 because log(0) is undefined)
        double xMax = 2.0;          // Maximum x value
        double binEdges[nBins + 1]; // Bin edges array
        for (int i = 0; i <= nBins; ++i) {
        binEdges[i] = xMin * pow(xMax / xMin, double(i) / nBins);
        }
        h97_xbj = new TH1F("h97_xbj", "Elasticity of the Scattering Process (x_{bj}); x_{bj}; Events", nBins, binEdges);

        h97_gamma_E_EC0 = new TH1F("h97_E_EC0", "Photon Energy ECal 0 - Before Excl. Cuts; E_{#gamma} [GeV]; Counts", 100, 0, 50);
        h97_gamma_E_EC1 = new TH1F("h97_E_EC1", "Photon Energy ECal 1 - Before Excl. Cuts; E_{#gamma} [GeV]; Counts", 100, 0, 100);
        h97_gamma_E_EC2 = new TH1F("h97_E_EC2", "Photon Energy ECal 2 - Before Excl. Cuts; E_{#gamma} [GeV]; Counts", 100, 0, 200);
        h97_E_miss      = new TH1F("h97_E_miss", "Missing Energy (All ECals); E_{#gamma} [GeV]; Counts", 100, -20, 20);
        h97_M2_miss     = new TH1F("h97_M2_miss", "Missing Mass Squared (All ECals); M^{2}_{#gamma} [GeV^{2}/c^{4}]; Counts", 100, 0, 300);

        h97_delta_phi = new TH1F("h97_delta_phi", "#Delta#phi = #phi^{cam} - #phi^{miss}; #Delta#phi [rad]; Counts", 100, -0.5, 0.5);
        h97_delta_pt = new TH1F("h97_delta_pt", "#DeltaP_{t} = P_{t}^{cam} - P_{t}^{miss}; #DeltaP_{t} [GeV/c]; Counts", 100, -0.4, 0.4);
        h97_delta_Z = new TH1F("h97_delta_Z", "#DeltaZ_{A} = Z_{A}^{cam} - Z_{A}^{miss}; #DeltaZ_{A} [cm]; Counts", 100, -20, 20);
        h97_M2x      = new TH1F("h97_M2x", "M^{2}_{undet} = (k + p - k'- q'- p')^{2}; M^{2}_{undet} [rad]; Counts", 100, -0.5, 0.5);

        h97_gamma_E_EC0_ex = new TH1F("h97_E_EC0_ex", "Photon Energy ECal 0 - After Excl. Cuts; E_{#gamma} [GeV]; Counts", 100, 0, 50);
        h97_gamma_E_EC1_ex = new TH1F("h97_E_EC1_ex", "Photon Energy ECal 1 - After Excl. Cuts; E_{#gamma} [GeV]; Counts", 100, 0, 100);
        h97_gamma_E_EC2_ex = new TH1F("h97_E_EC2_ex", "Photon Energy ECal 2 - After Excl. Cuts; E_{#gamma} [GeV]; Counts", 100, 0, 200);

        //
        // Ntuple definition 
        //
        tree = new TTree("USR97","User 97 DVCS NTuple"); // name (has to be unique) and title of the Ntuple
        
        tree->Branch("Run",     &Run,     "Run/I");
        tree->Branch("Evt",     &Evt,     "Evt/I");
        tree->Branch("Zprim",   &Zprim,   "Zprim/F");
        tree->Branch("Yprim",   &Yprim,   "Yprim/F");
        tree->Branch("Xprim",   &Xprim,   "Xprim/F");
        tree->Branch("Chi2",    &Chi2,    "Chi2/F");

        tree->Branch("inMu_pz", &inMu_pz, "inMu_pz/D");
        tree->Branch("inMu_py", &inMu_py, "inMu_py/D");
        tree->Branch("inMu_px", &inMu_px, "inMu_px/D");
        tree->Branch("inMu_E",  &inMu_E,  "inMu_E/D");
 
        tree->Branch("outMu_pz", &outMu_pz, "outMu_pz/D");
        tree->Branch("outMu_py", &outMu_py, "outMu_py/D");
        tree->Branch("outMu_px", &outMu_px, "outMu_px/D");
        tree->Branch("outMu_E",  &outMu_E,  "outMu_E/D");

        tree->Branch("y",   &y,   "y/D");
        tree->Branch("nu",  &nu,  "nu/D");
        tree->Branch("Q2",  &Q2,  "Q2/D");
        tree->Branch("W2",  &W2,  "W2/D");
        tree->Branch("xbj", &xbj, "xbj/D");
        tree->Branch("t",   &t,   "t/D");

        tree->Branch("E_miss",  &E_miss,  "E_miss/D");
        tree->Branch("M2_miss", &M2_miss, "M2_miss/D");

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
    //eventFlags.createFlag("Q2_DIS_flag", "No. of events where Q2 > 0.5");
    //eventFlags.createFlag("y_DIS_flag", "No. of events where 0.01 < y < 0.99");
    eventFlags.createFlag("Q2_flag", "No. of events where 1 < Q2 < 10");
    eventFlags.createFlag("y_flag", "No. of events where 0.05 < y < 0.9");

    eventFlags.createFlag("singleTrack_flag", "No. of events where primary vertex only has one outgoing track");
    eventFlags.createFlag("nCls_flag", "No. of events where multiple clusters are found in ECal 0, 1 or 2 only");
    eventFlags.createFlag("singleCl_flag", "No. of events where there is only a single cluster in the ECals");
    eventFlags.createFlag("proton_flag", "No. of events where proton candidates have 0.1 < beta < 0.95");
    eventFlags.createFlag("delta_pt_flag", "No. of events where |delta_pt| < 0.3 GeV/c");
    eventFlags.createFlag("delta_phi_flag", "No. of events where |delta_phi| < 0.4 rad");
    eventFlags.createFlag("delta_Z_flag", "No. of events where |delta_Z| < 16 cm");
    eventFlags.createFlag("M2x_flag", "No. of events where |(M_x)^2| < 0.3 (GeV/c^2)^2");
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
    // Debug statements ... (1/2)
    //printDebug("     ");
    //printDebug("*** Run: " + std::to_string(Run) + ", spill: " + std::to_string(Spill) + ", event: " + std::to_string(EvtInSpill) + " ***");
    
    //*******************************************  
    eventFlags.setFlagByName("allEvts_flag", true);
    
    // Loop over reconstructed vertices in the event 
		for (int iv = 0; iv < e.NVertex(); iv++) { // begin loop over vertices
			//******************************************* 
			// Store info about primary vertex (if found) 
			const PaVertex & v = e.vVertex(iv);
			if (!v.IsPrimary()) continue; // skip any vertices that are not primary 
      eventFlags.setFlagByName("pVtx_flag", true);
			Zprim = v.Pos(2);
			Yprim = v.Pos(1);
			Xprim = v.Pos(0);  
      Chi2  = v.Chi2(); 
      Nprim = v.NOutParticles(); // number of tracks in vertex

			//*******************************************
    	// Store info about incoming muon beam (inMu)
      static PaParticle beam; 
      static PaTrack beam_track; 
      static PaTPar Par_beam;

			static BeamFluxParams beamParams; // Create an instance of BeamFluxParams
			flux_flag = beamFluxCheck(e, v, iv, Run, TiS_flag, beamParams, beam, beam_track, Par_beam, eventFlags);
			if (!flux_flag) continue;  

      //*******************************************
      // Store info about scattered muon (outMu)
      static PaParticle outMu; 
      static PaTrack outMu_track; 
      static PaTPar Par_outMu;

      static OutMuParams outMuParams; // Create an instance of OutMuParams 
      outMu_flag = outMuCheck(e, v, iv, Run, beam, HodoHelper, trig_flag, outMuParams, 
                            outMu, outMu_track, Par_outMu, eventFlags);
      if (!outMu_flag) continue; 

      //*******************************************
      // Kinematic variables ... (1/2)
      TLorentzVector inMu_TL  = Par_beam.LzVec(M_mu); 
      TLorentzVector outMu_TL = Par_outMu.LzVec(M_mu);  
      TLorentzVector q        = (inMu_TL - outMu_TL); // four momentum of the virtual photon
      TLorentzVector targ_TL;

      targ_TL.SetPxPyPzE(0,0,0,M_p);
      Q2  = PaAlgo::Q2 (inMu_TL, outMu_TL); //Q2 = -(inMu_TL - outMu_TL).M2();
      y   = (inMu_TL.E() - outMu_TL.E()) / inMu_TL.E();
      nu  = (inMu_TL.E() - outMu_TL.E());
      W2  = PaAlgo::W2 (inMu_TL, outMu_TL);
      xbj = PaAlgo::xbj (inMu_TL, outMu_TL); //xbj = Q2/(2*q*targ_TL);

      // Current kinematic cuts may be tightened after the kinematically constrained fit is applied 
      bool test_flag = false; 
      //if (Q2 < 0.5) continue; // inclusive Q2 cut
      if (Q2 >= 0.5) {test_flag = true;} 
      if (!test_flag) continue; 
      eventFlags.setFlagByName("Q2_DIS_flag", true);
      //if (Q2 < 1 || Q2 > 10) continue; // exclusive Q2 cut 
      //eventFlags.setFlagByName("Q2_flag", true);

      //if (y < 0.01 || y > 0.99) continue; // inclusive y cut 
      //eventFlags.setFlagByName("y_DIS_flag", true);
      //if (y < 0.05 || y > 0.9) continue; // exclusive y cut 
      //eventFlags.setFlagByName("y_flag", true);

      double inMu_p = beam.ParInVtx(iv).Mom();
      double outMu_p = outMu.ParInVtx(iv).Mom();
       
      //*******************************************
      // Exclusive selection starts here  ...  
      // Only one outgoing particle (scattered proton and photon are detected using ECals and CAMERA so will not be found here)
      if (Nprim != 1) continue; 
      //std::cout << std::endl <<"DEBUG :: " << Evt << " " << Q2 << " " << y << std::endl;
      eventFlags.setFlagByName("singleTrack_flag", true);  

      //*******************************************
      // Store info about real photon - check ECals for single neutral cluster
      // ! CURRENTLY LOOKING FOR ANY NEUTRAL CLUSTERS, NOT SPECIFICALLY FOR ONE !
      //printDebug("     ", true);
      //printDebug("*** Run: " + std::to_string(Run) + ", spill: " + std::to_string(Spill) + ", event: " + std::to_string(EvtInSpill) + " ***", true);

      int ecal0id = PaSetup::Ref().iCalorim("EC00P1__");
      int ecal1id = PaSetup::Ref().iCalorim("EC01P1__");
      int ecal2id = PaSetup::Ref().iCalorim("EC02P1__");

      std::vector<PaCaloClus> Clus_gamma; // Store valid clusters dynamically
      std::vector<int> cl_id; // Store calorimeter IDs for valid clusters

      float EC0_thr = 4; // ECAL thresholds for DVCS event selection
      float EC1_thr = 5;
      float EC2_thr = 10;

      int clusterCount = 0;
      for (int iclus = 0; iclus < e.NCaloClus(); iclus++) { // Begin loop over clusters
        const PaCaloClus & cl = e.vCaloClus(iclus);
        int icalo = cl.iCalorim();

        if (cl.iTrack() != -1) continue; // cluster needs to have no associated charged track
        if (!EcalTimeCut(beam_track, cl, Run) && !e.IsMC()) continue; // check timing of the clusters
 
        if (plotECalEnergy) { // Fill the reconstructed energy distributions for each ECal
          if (icalo == ecal0id) h97_gamma_E_EC0->Fill(cl.E());
          if (icalo == ecal1id) h97_gamma_E_EC1->Fill(cl.E());
          if (icalo == ecal2id) h97_gamma_E_EC2->Fill(cl.E());
        }

        if ((icalo == ecal0id) && (cl.E() >= EC0_thr)) { // EC0
          Clus_gamma.push_back(cl); // Store the valid cluster
          cl_id.push_back(icalo);   // Store the ECal ID 
          //printDebug("    Evt: " + std::to_string(Evt) + ", ECal ID:" + std::to_string(icalo) + ", Energy: " + std::to_string(cl.E()), true);
          clusterCount++; 
        }

        else if ((icalo == ecal1id) && (cl.E() >= EC1_thr)) { // EC1
          Clus_gamma.push_back(cl); // Store the valid cluster
          cl_id.push_back(icalo);   // Store the ECal ID 
          //printDebug("    Evt: " + std::to_string(Evt) + ", ECal ID:" + std::to_string(icalo) + ", Energy: " + std::to_string(cl.E()), true);
          clusterCount++; 
        }

        else if ((icalo == ecal2id) && (cl.E() >= EC2_thr)) { //EC2
          Clus_gamma.push_back(cl); // Store the valid cluster
          cl_id.push_back(icalo);   // Store the ECal ID 
          //printDebug("    Evt: " + std::to_string(Evt) + ", ECal ID:" + std::to_string(icalo) + ", Energy: " + std::to_string(cl.E()), true);
          clusterCount++; 
        }

      } // End loop over clusters

      if (clusterCount == 0) continue; // skip event if there are no clusters
      eventFlags.setFlagByName("nCls_flag", true);

      // ! NOW CHECK THAT THERE IS ONLY ONE NEUTRAL CLUSTER !
      if (clusterCount != 1) continue; // skip event if there is more than one cluster
      eventFlags.setFlagByName("singleCl_flag", true);

      // Calculate TL vector for single photon
      TLorentzVector gamma_TL;

      if (!Clus_gamma.empty()) { // Ensure there is at least one cluster
        const auto& cluster = Clus_gamma.front(); // Directly access the first (and only) cluster
        int ecalid = cl_id.front(); // Retrieve cluster ID

        double phdz  = cluster.Z() - v.Z();
        double phdy  = cluster.Y() - v.Y();
        double phdx  = cluster.X() - v.X();
        double phL   = TMath::Sqrt(phdx * phdx + phdy * phdy + phdz * phdz);
        double cl_E  = cluster.E();

        gamma_TL.SetPxPyPzE(cl_E * phdx / phL, cl_E * phdy / phL, cl_E * phdz / phL, cl_E);
      }

      //*******************************************
      // Kinematic variables ... (2/2)
      t       = (q - gamma_TL).M2();
      E_miss  = nu - gamma_TL.E() + t/(2 * M_p); // Missing energy (assuming proton)
      M2_miss = (2 * M_p * (nu - gamma_TL.E())) + (M_p * M_p) + t; // Missing mass squared (assuming proton) 

      //*******************************************
      // Store info about scattered proton candidates (check CAMERA)
      // ! CURRENTLY LOOKING FOR ANY PROTONS, NOT SPECIFICALLY FOR ONE !
      vector <CameraProton> proton_candidates = cam_inst->GetGoodCandidates(v); // vector holding all proton candidates
      vector <CameraProton> protons; // vector for storing candidates that satisfy 0.1 < beta < 0.95 

      for (auto proton: proton_candidates) {
        double beta = proton.beta; 
        if (beta >= 0.1 || beta <= 0.95) {
          protons.emplace_back(proton);
        } 
      }

      if (protons.empty()) continue; // skip events with no good protons 
      eventFlags.setFlagByName("proton_flag", true); 

      //*******************************************
      // Check that all combintations of the vertex, photon and proton satisfy exclusivity conditions 
      // 4 exclusivity variables: detla_phi, delta_pt, delta_Z, M2x
      int nCombs = 0;
      vector <CameraProton> protons_excl; // vector for storing candidates which pass the exclusivity cuts 

      TLorentzVector p_miss_TL = targ_TL + inMu_TL - outMu_TL - gamma_TL;
      double phi_miss = p_miss_TL.Phi();
      double pt_miss  = p_miss_TL.Vect().Pt(); 

      TVector3 R_vtx;
      R_vtx.SetXYZ(v.Pos(0),v.Pos(1),v.Pos(2));

      // zero missing mass
      double M2x = (inMu_TL + targ_TL - outMu_TL - gamma_TL - p_miss_TL) * (inMu_TL + targ_TL - outMu_TL - gamma_TL - p_miss_TL); 

      for (auto proton: protons) { 
        TLorentzVector p_camera_TL = proton.p4; 
        double phi_camera = p_camera_TL.Phi();
        double delta_phi  = phi_camera - phi_miss; // azimuthal angle

        double pt_camera = p_camera_TL.Vect().Pt(); 
        double delta_pt  = std::abs(pt_camera) - std::abs(pt_miss); // transverse momentum

        TVector3 posRingA = proton.Ahit.vec;
			  TVector3 posRingB = proton.Bhit.vec;
        double Z_miss     = cam_inst->GetZA(R_vtx, posRingA, posRingB, phi_miss);
        double delta_Z    = posRingA.Z() - Z_miss; // z position of the hits in the inner CAMERA ring

        h97_delta_phi->Fill(delta_phi);
        h97_delta_pt->Fill(delta_pt);
        h97_delta_Z->Fill(delta_Z);
        h97_M2x->Fill(M2x);

        // Apply cuts on the exclusivity variables 
        if (std::fabs(delta_pt) > 0.3) continue; 
        eventFlags.setFlagByName("delta_pt_flag", true);

        if (std::fabs(delta_phi) > 0.4) continue; 
        eventFlags.setFlagByName("delta_phi_flag", true);

        if (std::fabs(delta_Z) > 16) continue; 
        eventFlags.setFlagByName("delta_Z_flag", true);

        if (std::fabs(M2x) > 0.3) continue; 
        eventFlags.setFlagByName("M2x_flag", true);
        nCombs++; 
      }

      if (protons_excl.empty()) continue;
      eventFlags.setFlagByName("nCombo_flag", true);

/*       printDebug("    ", true);
      printDebug("*** Run: " + std::to_string(Run) + ", spill: " + std::to_string(Spill) + ", event: " + std::to_string(EvtInSpill) + " ***", true);
      printDebug("    Vertex: (" + std::to_string(Xprim) + ", " + std::to_string(Yprim) + ", " + std::to_string(Zprim) + ")", true);
      printDebug("    No. protons: " + std::to_string(protons_excl.size()) + ", No. photons: " + std::to_string(Clus_gamma.size()), true); */

      //*******************************************
      // Fill histograms for ECal Energy after all exclusivity cuts
      if (plotECalEnergy) {
        if (!Clus_gamma.empty()) { // Ensure there is at least one cluster
          const auto& cluster = Clus_gamma.front(); // Directly access the first (and only) cluster
          int ecalid = cl_id.front(); // Retrieve cluster ID
          if (ecalid == ecal0id) {h97_gamma_E_EC0_ex->Fill(cluster.E());}
          if (ecalid == ecal1id) {h97_gamma_E_EC1_ex->Fill(cluster.E());}
          if (ecalid == ecal2id) {h97_gamma_E_EC2_ex->Fill(cluster.E());}
        }
      }

      //*******************************************
      // Debug statements ... (2/2)
      printDebug("     ");
      printDebug("*** Run: " + std::to_string(Run) + ", spill: " + std::to_string(Spill) + ", event: " + std::to_string(EvtInSpill) + " ***");
      printDebug("    Vertex: (" + std::to_string(Xprim) + ", " + std::to_string(Yprim) + ", " + std::to_string(Zprim) + ")");
      printDebug("    mu: P: " + std::to_string(inMu_p) + " GeV/c, Charge: " + std::to_string(beam.Q()));
      printDebug("    mu': P: " + std::to_string(outMu_p) + " GeV/c, Charge: " + std::to_string(outMu.Q()));
      printDebug("    Kinematics: Q2: " + std::to_string(Q2) +  " GeV2, y: " + std::to_string(y) + ", W2: " + std::to_string(W2) + " GeV2, x: " + std::to_string(xbj));

      //*******************************************
      // Fill histrograms 
      inMu_py = beam_track.vTPar(0).Py();
      inMu_px = beam_track.vTPar(0).Px();

      outMu_py = outMu_track.vTPar(0).Py();
      outMu_px = outMu_track.vTPar(0).Px();

      h97_Zprim->Fill(Zprim);
      h97_Yprim->Fill(Yprim);
      h97_Xprim->Fill(Xprim);
      h97_XYprim->Fill(Xprim,Yprim);

      h97_inMu_p->Fill(inMu_p);
      h97_inMu_py->Fill(inMu_py);
      h97_inMu_px->Fill(inMu_px);

      h97_outMu_p->Fill(outMu_p);
      h97_outMu_py->Fill(outMu_py);
      h97_outMu_px->Fill(outMu_px);

      h97_y->Fill(y);
      h97_nu->Fill(nu);
      h97_Q2->Fill(Q2);
      h97_W2->Fill(W2);
      h97_xbj->Fill(xbj); 
      h97_Q2xbj->Fill(xbj,Q2);
      h97_Chi2Q2->Fill(Chi2,Q2);
      h97_t->Fill(t);

      h97_E_miss->Fill(E_miss);
      h97_M2_miss->Fill(M2_miss);     

		} // end loop over vertices 

  // Increment all counters whose flags are "true"
  eventFlags.incrementCounters(); 

  //Save the event
  //if (saveEvent_flag) {e.TagToSave();}
  e.TagToSave();  

} // end event loop 

void UserJobEnd970() {
  if (plotECalEnergy) {
    std::cout << "! WARNING: ECAL ENERGY PLOT IS ON. TURN OFF TO IMPROVE EFFICIENCY. !"<< std::endl;
    std::cout << std::endl;
  }

  eventFlags.printFlags(); // Print to output stream

}




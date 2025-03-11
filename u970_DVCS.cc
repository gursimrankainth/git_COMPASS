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
//#include "Tools_Camera.hh"

#include "ecal_time_cuts.h"
#include "/Users/gursimran/cern/phastPackages/xcheck_newTiS/tis_range.cc"
#include "/Users/gursimran/cern/phastPackages/xcheck_newTiS/tis_range.h"
#include "/Users/gursimran/cern/phastPackages/xcheck_newTiS/xcheck_cuts_flags.cc"
#include "/Users/gursimran/cern/phastPackages/xcheck_newTiS/xcheck_cuts_flags.h"
//#include "/afs/cern.ch/user/g/gkainth/phastPackages/flux_files/flux_Johannes/functions/TiS_range.cc"

#include "UConn_Tools.h"

// ************************************************************************** //
// UserEvent for preselecting exclusive photon events	(DVCS)			      //
// In the selection all possible combinations of:					          //
// Vertices (incoming and outgoing muons), exclusive photon and recoil proton //
// ************************************************************************** //

extern "C" float prob_(float&, int&);

// Global flag for verbose mode 
bool verbose_mode = true; // Create an instance of verbose_mode

//*****************************************************************************
void UserEvent970(PaEvent & e) { // begin event loop

    // Define constants
    static PaHodoHelper* HodoHelper = NULL;
    static TiSRange*     tis_range  = NULL; 

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

    static TH1F* h97_y     = NULL;
    static TH1F* h97_Q2    = NULL;
    static TH1F* h97_W2    = NULL;
    static TH1F* h97_xbj   = NULL;
    static TH2F* h97_Q2xbj = NULL;

    static TH1F* h97_gamma_E_EC0 = NULL;
    static TH1F* h97_gamma_E_EC1 = NULL;
    static TH1F* h97_gamma_E_EC2 = NULL;
    static TH1F* h97_E_miss      = NULL;
    static TH1F* h97_M2_miss     = NULL;

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
		bool TiS_flag  = false;  
    bool trig_flag = false;

    struct EventFlags {
        //DIS cuts
        bool flux_flag  = false; 
        bool outMu_flag = false;
        bool Q2_flag    = false;
        bool y_flag     = false;
        bool DIS_flag   = false; // true if all DIS cuts have been applied 
        //Exclusive cuts 
        bool singleTrack_flag = false;
        bool checkCls_flag    = false; // true if all found clusters have a track, and pass the timing and ECal energy cuts 
        bool singleCl_flag    = false; // true after a single neutral cluster has been found in the ECals
        bool exclEvent_flag   = false; // true if all exclusive cuts have ben applied 
        bool saveEvent_flag   = false;
    };

    EventFlags eventFlags;   // Create an instance of event flags

    // Register XCheck flags to track cut statistics  
    XCHECK_REGISTER_FLAG(Event_AllEvents, "Total no. of events processed by PHAST user script");
    XCHECK_REGISTER_FLAG(Event_PVtx, "No. of events with a primary vertex");
		XCHECK_REGISTER_FLAG(Event_BeamTrack, "No. of events where beam has a track with parameters");
		XCHECK_REGISTER_FLAG(Event_ZFirst, "No. of events where beam was first measured before the target");
		XCHECK_REGISTER_FLAG(Event_BeamMom, "No. of events where beam momentum falls within acceptable range");
		XCHECK_REGISTER_FLAG(Event_BeamMomErr, "No. of events where beam momentum error falls within acceptable range");
		XCHECK_REGISTER_FLAG(Event_BMS, "No. of events where beam is detected by BMS");
		XCHECK_REGISTER_FLAG(Event_FI, "No. of events where beam is detected by SCIFI");
		XCHECK_REGISTER_FLAG(Event_SI, "No. of events where beam is detected by SI");
		XCHECK_REGISTER_FLAG(Event_CrossCells, "No. of events where the beam crosses full target length");
		XCHECK_REGISTER_FLAG(Event_Meantime, "No. of events where beam track meantime is within acceptable flux requirements");
		XCHECK_REGISTER_FLAG(Event_TiS, "No. of events where time in spill is within acceptable flux requirements");
    XCHECK_REGISTER_FLAG(Event_Flux, "No. of events where all flux requirements are satisfied");
    XCHECK_REGISTER_FLAG(Event_InTarget, "No. of events where scattered muon vertex is in target");
    XCHECK_REGISTER_FLAG(Event_Trigger, "No. of events with MT, LT, OT or LAST physics triggers");
    XCHECK_REGISTER_FLAG(Event_Hodo, "No. of events where scattered muon passes Hodoscope check");
    XCHECK_REGISTER_FLAG(Event_Charge, "No. of events where scattered muon has the same charge as the beam");
    XCHECK_REGISTER_FLAG(Event_ZFirstLast, "No. of events where first and last scattered muon z coord. are measured before and after SM1");
    XCHECK_REGISTER_FLAG(Event_outMu, "No. of events where scattered muon passes all checks");
    XCHECK_REGISTER_FLAG(Event_Q2, "No. of events where 1 < Q2 < 10");
    XCHECK_REGISTER_FLAG(Event_Y, "No. of events where 0.05 < y < 0.9");
    XCHECK_REGISTER_FLAG(Event_SingleTrack, "No. of events where primary vertex only has one outgoing track");
    XCHECK_REGISTER_FLAG(Event_MultTracks, "No. of events where primary vertex has multiple outgoing tracks");
    XCHECK_REGISTER_FLAG(Event_SingleCl, "No. of events where there is only a single neutral clsuter in the ECals");

    XCHECK_REGISTER_FLAG(Event_DIS, "No. of events that pass all DIS cuts"); // TBD 

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

        h97_y     = new TH1F("h97_y", "Fractional Energy Loss of Incoming Muon (y); y; Events", 100, 0, 1);
        h97_Q2    = new TH1F("h97_Q2", "Four-momentum Transfer Squared (Q^{2}); Q^{2} [GeV^{2}]; Events", 100, 0, 10);
        h97_W2    = new TH1F("h97_W2", "Effective Mass of final state hadrons Squared (W^{2}); W^{2} [GeV^{2}]; Events", 100, 0, 350);
        h97_Q2xbj = new TH2F("h97_Q2xbj", "Kinematic Coverage of Dataset; x_{bj}; Q^{2} [GeV^{2}]", 100, 0, 1, 100, 0, 100);

        const int nBins = 100;      // Number of bins
        double xMin = 1e-3;         // Minimum x value (avoid 0 because log(0) is undefined)
        double xMax = 2.0;          // Maximum x value
        double binEdges[nBins + 1]; // Bin edges array
        for (int i = 0; i <= nBins; ++i) {
        binEdges[i] = xMin * pow(xMax / xMin, double(i) / nBins);
        }
        h97_xbj = new TH1F("h97_xbj", "Elasticity of the Scattering Process (x_{bj}); x_{bj}; Events", nBins, binEdges);

        h97_gamma_E_EC0 = new TH1F("h97_E_EC0", "Photon Energy ECal 0 - Low Energy Cut; E_{#gamma} [GeV]; Counts", 100, 0, 50);
        h97_gamma_E_EC1 = new TH1F("h97_E_EC1", "Photon Energy ECal 1 - Low Energy Cut; E_{#gamma} [GeV]; Counts", 100, 0, 100);
        h97_gamma_E_EC2 = new TH1F("h97_E_EC2", "Photon Energy ECal 2 - Low Energy Cut; E_{#gamma} [GeV]; Counts", 100, 0, 200);
        h97_E_miss      = new TH1F("h97_E_miss", "Missing Energy (All ECals); E_{#gamma} [GeV]; Counts", 100, -20, 20);
        h97_M2_miss     = new TH1F("h97_M2_miss", "Missing Mass Squared (All ECals); M^{2}_{#gamma} [GeV^{2}/c^{4}]; Counts", 100, 0, 300);

        //
        // Ntuple definition 
        //
        tree = new TTree("USR97","User 97 DVCS NTuple"); // name (has to be unique) and title of the Ntuple
        
        tree->Branch("Run",     &Run,     "Run/I");
        tree->Branch("Evt",     &Evt,     "Evt/I");
        tree->Branch("Zprim",   &Zprim,   "Zprim/F");
        tree->Branch("Yprim",   &Yprim,   "Yprim/F");
        tree->Branch("Xprim",   &Xprim,   "Xprim/F");

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

        first=false;
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
    // Initialize variables and check time in spill (time in spill cut is applied later not here)
    Run         = e.RunNum();
    Evt         = e.UniqueEvNum();
    Year        = e.Year();
    EvtInSpill  = e.EvInSpill(); 
    Spill       = e.SpillNum(); 
    TimeInSpill = e.TimeInSpill();  

    if (Run != LastRun) { // Reinitialize HodoHelper and tis_range only if the run number changes 
        HodoHelper = & PaHodoHelper::Init("", true);  
        tis_range  = new TiSRange("/Users/gursimran/cern/phastPackages/flux_files/flux_Johannes/2016/flux_files");
        //Set_TiSrange("/afs/cern.ch/user/g/gkainth/phastPackages/flux_files/flux_Johannes/2016/flux_files", Run, Run);
        LastRun = Run;  // Update LastRun to the current run number
    }
 
    TiS_flag = tis_range->CheckWindow(Run, Spill, TimeInSpill); // check the time in spill 
    
    printDebug("     ");
    printDebug("*** Run: " + std::to_string(Run) + ", spill: " + std::to_string(Spill) + ", event: " + std::to_string(EvtInSpill) + " ***");
    
    //*******************************************  
    XCHECK_COUNT_FLAG(Event_AllEvents, "Total no. of events processed by PHAST user script"); 
    
    // Loop over reconstructed vertices in the event 
		for (int iv = 0; iv < e.NVertex(); iv++) { // begin loop over vertices
			//******************************************* 
			// Store info about primary vertex (if found) 
			const PaVertex & v = e.vVertex(iv);
			if (!v.IsPrimary()) continue;
			XCHECK_COUNT_FLAG(Event_PVtx, "No. of events with a primary vertex");
			Zprim = v.Pos(2);
			Yprim = v.Pos(1);
			Xprim = v.Pos(0); 
      Nprim = v.NOutParticles(); // number of tracks in vertex

			//*******************************************
    	// Store info about incoming muon beam (inMu)
      static PaParticle beam; 
      static PaTrack beam_track; 
      static PaTPar Par_beam;

			static BeamFluxParams beamParams; // Create an instance of BeamFluxParams
			eventFlags.flux_flag = beamFluxCheck(e, v, iv, Run, TiS_flag, beamParams, beam, beam_track, Par_beam);
			if (!eventFlags.flux_flag) continue; 
      XCHECK_COUNT_FLAG(Event_Flux, "No. of events where all flux requirements are satisfied");

      //*******************************************
      // Store info about scattered muon (outMu)
      static PaParticle outMu; 
      static PaTrack outMu_track; 
      static PaTPar Par_outMu;

      static OutMuParams outMuParams; // Create an instance of OutMuParams 
      eventFlags.outMu_flag = outMuCheck(e, v, iv, Run, beam, HodoHelper, trig_flag, outMuParams, outMu, outMu_track, Par_outMu);
      if (!eventFlags.outMu_flag) continue; 
      XCHECK_COUNT_FLAG(Event_outMu, "No. of events where scattered muon passes all checks");

      //*******************************************
      // Kinematic variables ... (1/2)
      TLorentzVector inMu_TL  = Par_beam.LzVec(M_mu); 
      TLorentzVector outMu_TL = Par_outMu.LzVec(M_mu); 
      //TLorentzVector targ_TL(0,0,0,M_p);
      TLorentzVector q = (inMu_TL - outMu_TL); // four momentum of the virtual photon

      Q2  = PaAlgo::Q2 (inMu_TL, outMu_TL); //Q2  = -(inMu_TL - outMu_TL).M2();
      y   = (inMu_TL.E() - outMu_TL.E()) / inMu_TL.E();
      nu  = (inMu_TL.E() - outMu_TL.E());
      W2  = PaAlgo::W2 (inMu_TL, outMu_TL);
      xbj = PaAlgo::xbj (inMu_TL, outMu_TL); //xbj = Q2/(2*q*targ_TL);

      // Current kinematic cuts may be tightened after the kinematically constrained fit is applied 
      if (Q2 < 1 || Q2 > 10) continue;
      eventFlags.Q2_flag = (Q2 > 1 && Q2 < 10);
      XCHECK_COUNT_FLAG(Event_Q2, "No. of events where 1 < Q2 < 10");

      if (y < 0.05 || y > 0.9) continue; 
      eventFlags.y_flag = (y > 0.05 && y < 0.9);
      XCHECK_COUNT_FLAG(Event_Y, "No. of events where 0.05 < y < 0.9");

      double inMu_mom = beam_track.vTPar(0).Mom();
      double outMu_mom = outMu_track.vTPar(0).Mom();
      //*******************************************
      // Exclusive selection starts here  ... 
      // Only one outgoing particle (scattered proton and photon are detected using ECals and CAMERA so will not be found here)
      if (Nprim == 1) {
        eventFlags.singleTrack_flag = true; 
        XCHECK_COUNT_FLAG(Event_SingleTrack, "No. of events where primary vertex only has one outgoing track");
      } 
      else {
        XCHECK_COUNT_FLAG(Event_MultTracks, "No. of events where primary vertex has multiple outgoing tracks");
        printDebug("    NPrim:" + std::to_string(Nprim) + ", mu: P: " + std::to_string(inMu_mom) + " GeV/c");
        printDebug("    NPrim:" + std::to_string(Nprim) + ", mu': P: " + std::to_string(outMu_mom) + " GeV/c");
      }

      //*******************************************

		} // end loop over vertices 

	//if (saveEvent_flag) {e.TagToSave();}
  e.TagToSave();  

} // end event loop 

void UserJobEnd970() {
  XCHECK_JOB_END(); // Print to output stream
}




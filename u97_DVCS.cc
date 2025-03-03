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
#include "PaHodoHelper.h"
#include "PaEvent.h" 
#include "G3part.h" 

#include "/Users/gursimran/cern/flux_files/flux_Johannes/functions/TiS_range.cc"
//#include "/afs/cern.ch/user/g/gkainth/flux_files/flux_Johannes/functions/TiS_range.cc"
#include "ecal_time_cuts.h"

/// ! THERE IS AN ISSUE WITH HODOHELPER SO YOU CAN HAVE ONE INPUT FILE AT A TIME!

// ************************************************************************** //
// UserEvent for preselecting exclusive photon events	(DVCS)			            //
// In the selection all possible combinations of:					                    //
// Vertices (incoming and outgoing muons), exclusive photon and recoil proton	//
// ************************************************************************** //

extern "C" float prob_(float&, int&);

// Global flag for verbose mode 
bool verbose = false; // Set to true for verbose output, false to suppress
void printDebug(const std::string &message, bool forcePrint = false) {
    if (verbose || forcePrint) {
        std::cout << message << std::endl;
    }
}

//*****************************************************************************
void UserEvent97(PaEvent & e) {

  // Define constants
  static PaHodoHelper* HodoHelper = NULL; 

  const double M_pi    = G3partMass[8]; // Pion mass 
  const double M_gamma = G3partMass[1]; // Photon mass
  const double M_mu    = G3partMass[5]; // Muon mass
  const double M_p     = G3partMass[14]; // Proton mass 

  // Declare all objects globally
  static TH1F *h97_cutStats = NULL;

  static TH1F* h97_Zprim  = NULL; 
  static TH1F* h97_Yprim  = NULL; 
  static TH1F* h97_Xprim  = NULL; 
  static TH2F* h97_XYprim = NULL;

  static TH1F* h97_inMu_p  = NULL;
  static TH1F* h97_inMu_py = NULL;
  static TH1F* h97_inMu_px = NULL; 
  //static TH1F* h97_inMu_theta = NULL;
  //static TH1F* h97_inMu_phi   = NULL;

  static TH1F* h97_outMu_p  = NULL;
  static TH1F* h97_outMu_py = NULL;
  static TH1F* h97_outMu_px = NULL; 
  //static TH1F* h97_outMu_theta = NULL;
  //static TH1F* h97_outMu_phi   = NULL;

  static TH1F* h97_y     = NULL;
  static TH1F* h97_Q2    = NULL;
  static TH1F* h97_W2    = NULL;
  static TH1F* h97_xbj   = NULL;
  static TH2F* h97_Q2xbj = NULL;

  static TH1F* h97_gamma_E_EC0 = NULL;
  static TH1F* h97_gamma_E_EC1 = NULL;
  static TH1F* h97_gamma_E_EC2 = NULL;
  static TH1F* h97_E_miss      = NULL;

  static TTree* tree(NULL);

  //
  // Variables to be stored into analysis Ntuple
  //
  // (if you want to extend the list:
  // - add variable here,
  // - add in Ntuple definition below 
  // and do not forget to fill the variable in the code
  //
  static int    Run;         // run number
  static int    LastRun = -1;     // store the previous run number (used to reintialize Hodohelper if there are multiple input files)
  static int    Evt;         // event number - unique evt number (based on run, spill and evt num in spill) 
  static int    Year;        // year data was taken 
  static int    EvtInSpill;  // event number in spill 
  static int    Spill;       // spill number
  static double TimeInSpill; // time in spill 
  static float  Zprim;       // Z coordinate of primary vertex (10000 if not found)
  static float  Yprim;       // Y coordinate of primary vertex (10000 if not found)
  static float  Xprim;       // X coordinate of primary vertex (10000 if not found)
  static int    Nprim;       // Number of tracks in primary vertex (-1 in fot found)
  static float  Ch2prim;     // Chi2 of primary vertex
  static int    trig_mask;

  static double inMu_pz; // Z component of the beam muon momentum 
  static double inMu_py; // Y component of the beam muon momentum 
  static double inMu_px; // X component of the beam muon momentum
  static double inMu_p;  // Magnitude of the beam muon momentum 
  static double inMu_E;  // Magnitude of the beam muon energy
  //static double inMu_theta; // Theta angle of beam track 
  //static double inMu_phi; // Phi angle of beam track 
  
  static int    Hodo_i_omu;
  static double outMu_pz; // Z component of the scattered muon momentum 
  static double outMu_py; // Y component of the scattered muon momentum
  static double outMu_px; // X component of the scattered muon momentum 
  static double outMu_p;  // Magnitude of the scattered muon momentum 
  static double outMu_E;  // Magnitude of the scattered muon energy
  //static double outMu_theta; // Theta angle of scattered muon track 
  //static double outMu_phi; // Phi angle of scattered muon track 

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
  struct EventFlags {
    //DIS cuts
    bool flux_flag = false;
    bool hodo_flag = false;
    bool outMuTrack_flag = false;
    bool vertex_flag = false;
    bool charge_flag = false;
    bool zFirstLast_flag = false;
    bool trig_flag = false;
    bool Q2_flag = false;
    bool y_flag = false;
    bool DIS_flag = false; // true if all DIS cuts have been applied 
    //Exclusive cuts 
    bool singleTrack_flag = false;
    bool checkCls_flag    = false; // true if all found clusters have a track, and pass the timing and ECal energy cuts 
    bool singleCl_flag    = false; // true after a single neutral cluster has been found in the ECals
    bool exclEvent_flag   = false; // true if all exclusive cuts have ben applied 
    bool saveEvent_flag   = false;
  };

  struct FluxFlags {
    bool allVtx_flag = false;
    bool pVtx_flag = false;
    bool inMu_flag = false;
    bool inMuTrack_flag = false;
    bool inMuPar_flag = false;
    bool passTarget_flag = false;
    bool zFirst_flag = false;
    bool BMS_flag = false;
    bool FI_flag = false;
    bool SI_flag = false;
    bool momRange_flag = false;
    bool momErr_flag = false;
    bool meantime_flag = false;
    bool TiS_flag = false;
  };

  EventFlags   eventFlags;   // Create an instance of event flags
  FluxFlags    fluxFlags;    // Create an instance of event flags

  // Selection statics (how many events remain after each cut)
  static std::map<std::string, std::pair<int, bool>> cutCounts = {
    {"Cut 00", {0, false}},  
    {"Cut 01", {0, false}},   
    {"Cut 02", {0, false}}, 
    {"Cut 03", {0, false}},   
    {"Cut 04", {0, false}},   
    {"Cut 05", {0, false}},     
    {"Cut 06", {0, false}},  
    {"Cut 07", {0, false}}, 
    {"Cut 08", {0, false}}, 
    {"Cut 09", {0, false}}, 
    {"Cut 10", {0, false}}, 
    {"Cut 11", {0, false}},  
    {"Cut 12", {0, false}},  
    {"Cut 13", {0, false}},
    {"Cut 14", {0, false}},
    {"Cut 15", {0, false}},
    {"Cut 16", {0, false}},
    {"Cut 17", {0, false}},
    {"Cut 18", {0, false}},
    {"Cut 19", {0, false}},
    {"Cut 20", {0, false}},
    {"Cut 21", {0, false}},
    {"Cut 22", {0, false}},
    {"Cut 23", {0, false}},
    {"Cut 24", {0, false}},   
    {"Cut 25", {0, false}},
  };

  //*****************************************************************************

  static bool first(true);
  if (first) { // histograms and Ntuples booking block
    Phast::Ref().HistFileDir("UserEvent97");

    if (HodoHelper == NULL) {
      std::cout << " *** HodoHelper is NULL. Initializing ... ***\n" << std::flush;
      HodoHelper = &PaHodoHelper::Init("", true);
      //HodoHelper = &PaHodoHelper::Init("/Users/gursimran/cern/phast.8.032/dat/trigger_config/2016",true);
      //HodoHelper = &PaHodoHelper::Init("/afs/cern.ch/user/g/gkainth/phast/dat/trigger_config/2016",true);
    }

    // 1D and 2D
    int nCuts = cutCounts.size();
    h97_cutStats = new TH1F("h97_cutStats", "Event Cuts Breakdown; Cut; Number of Events Removed", nCuts, 0, nCuts);

    h97_Zprim  = new TH1F("h97_Zprim", "Primary Vertex Z (cm); Z [cm]; Events", 100, -250, 250);
    h97_Yprim  = new TH1F("h97_Yprim", "Primary Vertex Y (cm); Y [cm]; Events", 100, -3, 3);
    h97_Xprim  = new TH1F("h97_Xprim", "Primary Vertex X (cm); X [cm]; Events", 100, -3, 3);
    h97_XYprim = new TH2F("h97_XYprim", "Primary Vertex XY (cm); X [cm]; Y [cm]", 100, -3, 3, 100, -3, 3);

    h97_inMu_p  = new TH1F("h97_inMu_p", "P Incoming Muon (GeV/c); P [GeV/c]; Events", 100, 0, 200);
    h97_inMu_py = new TH1F("h97_inMu_py", "Py Incoming Muon (GeV/c); Py [GeV/c]; Events", 100, 0, 10);
    h97_inMu_px = new TH1F("h97_inMu_px", "Px Incoming Muon (GeV/c); Px [GeV/c]; Events", 100, 0, 10);
    //h97_inMu_theta = new TH1F("h97_inMu_theta", "Theta Incoming Muon (rad); Theta [rad]; Events", 100, 0, 0.1);
    //h97_inMu_phi   = new TH1F("h97_inMu_phi", "Phi Incoming Muon (rad); Phi [rad]; Events", 100, -M_PI, M_PI);

    h97_outMu_p  = new TH1F("h97_outMu_p", "P Scattered Muon (GeV/c); P [GeV/c]; Events", 100, 0, 200);
    h97_outMu_py = new TH1F("h97_outMu_py", "Py Scattered Muon (GeV/c); Py [GeV/c]; Events", 100, 0, 10);
    h97_outMu_px = new TH1F("h97_outMu_px", "Px Scattered Muon (GeV/c); Px [GeV/c]; Events", 100, 0, 10);
    //h97_outMu_theta = new TH1F("h97_outMu_theta", "Theta Scattered Muon (rad); Theta [rad]; Events", 100, 0, 0.1);
    //h97_outMu_phi   = new TH1F("h97_outMu_phi", "Phi Scattered Muon (rad); Phi [rad]; Events", 100, -M_PI, M_PI);

    h97_y     = new TH1F("h97_y", "Fractional Energy Loss of Incoming Muon (y); y; Events", 100, 0, 1);
    h97_Q2    = new TH1F("h97_Q2", "Four-momentum Transfer Squared (Q2); Q2 [GeV2]; Events", 100, 0, 10);
    h97_W2    = new TH1F("h97_W2", "Effective Mass of final state hadrons Squared (W2); W2 [GeV2]; Events; W2 [GeV2]", 100, 0, 350);
    h97_Q2xbj = new TH2F("h97_Q2xbj", "Kinematic Coverage of Dataset; x_bj; Q2 [GeV2]", 100, 0, 1, 100, 0, 100);

    const int nBins = 100;       // Number of bins
    double xMin = 1e-3;          // Minimum x value (avoid 0 because log(0) is undefined)
    double xMax = 2.0;           // Maximum x value
    double binEdges[nBins + 1];  // Bin edges array
    for (int i = 0; i <= nBins; ++i) {
      binEdges[i] = xMin * pow(xMax / xMin, double(i) / nBins);
    }
    h97_xbj = new TH1F("h97_xbj", "Elasticity of the Scattering Process (x_bj); x_bj; Events", nBins, binEdges);

    h97_gamma_E_EC0 = new TH1F("h97_E_EC0", "Photon Energy ECal 0 - Low Energy Cut; E_{#gamma} [GeV]; Counts", 100, 0, 50);
    h97_gamma_E_EC1 = new TH1F("h97_E_EC1", "Photon Energy ECal 1 - Low Energy Cut; E_{#gamma} [GeV]; Counts", 100, 0, 100);
    h97_gamma_E_EC2 = new TH1F("h97_E_EC2", "Photon Energy ECal 2 - Low Energy Cut; E_{#gamma} [GeV]; Counts", 100, 0, 200);
    h97_E_miss      = new TH1F("h97_E_miss", "Missing Energy (All ECals); E_{#gamma} [GeV]; Counts", 100, -20, 20);

    //
    // Ntuple definition 
    //
    tree = new TTree("USR97","User 97 DVCS NTuple"); // name (has to be unique) and title of the Ntuple
    
    tree->Branch("Run",     &Run,     "Run/I");
    tree->Branch("Evt",     &Evt,     "Evt/I");
    tree->Branch("Zprim",   &Zprim,   "Zprim/F");
    tree->Branch("Yprim",   &Yprim,   "Yprim/F");
    tree->Branch("Xprim",   &Xprim,   "Xprim/F");
    tree->Branch("Nprim",   &Nprim,   "Nprim/I");
    tree->Branch("Ch2prim", &Ch2prim, "Chi2prim/F");

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

    tree->Branch("E_miss", &E_miss, "E_miss/D");

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
    eventFlags.trig_flag = true;
  }
  if (trig_mask & LT) {
    trigCheck += "LT ";
    eventFlags.trig_flag = true;
  }
  if (trig_mask & OT) {
    trigCheck += "OT ";
    eventFlags.trig_flag = true;
  }
  //if (trig_mask & LAST) {
  //  trigCheck += "LAST ";
  //  eventFlags.trig_flag = true;
  //}

  //*******************************************
  Run        = e.RunNum();
  Evt        = e.UniqueEvNum();
  Year       = e.Year();
  EvtInSpill = e.EvInSpill(); 
  Spill      = e.SpillNum(); 

  // Time in spill check 
  TimeInSpill  = e.TimeInSpill(); // check the time in spill 
  // Set_TiSrange("PATH_TO_FLUX_FILES", runMin, runMax);
  Set_TiSrange("/Users/gursimran/cern/flux_files/flux_Johannes/2016/flux_files", Run, Run);
  //Set_TiSrange("/afs/cern.ch/user/g/gkainth/flux_files/flux_Johannes/2016/flux_files", Run, Run);
  fluxFlags.TiS_flag = Check_TiS_window(Run,Spill,TimeInSpill);

  //*******************************************  
  // Reset the "already counted" flags for each event 
  for (auto& cut : cutCounts) {
    cut.second.second = false; 
  }
  cutCounts["Cut 00"].first++;

  // Loop over reconstructed vertices in the event 
  for (int iv = 0; iv < e.NVertex(); iv++) { // begin loop over vertices
  fluxFlags.allVtx_flag = true;
  if (fluxFlags.allVtx_flag && !cutCounts["Cut 01"].second) {
    cutCounts["Cut 01"].first++;
    cutCounts["Cut 01"].second = true;
  } 

    //******************************************* 
    // Store info about primary vertex (if found) 
    const PaVertex & v = e.vVertex(iv);
    if (!v.IsPrimary()) continue;
    fluxFlags.pVtx_flag = v.IsPrimary();
    if (fluxFlags.pVtx_flag && !cutCounts["Cut 02"].second) {
      cutCounts["Cut 02"].first++;
      cutCounts["Cut 02"].second = true;
    }
    Zprim   = v.Pos(2);
    Yprim   = v.Pos(1);
    Xprim   = v.Pos(0);
    Ch2prim = v.Chi2();
    Nprim   = v.NOutParticles(); // number of tracks in vertex 

    //*******************************************
    // Store info about incoming muon beam (inMu)
    int muBeam = v.InParticle(); 
    if (muBeam == -1) continue; // there is incoming particle associated with the primary vertex
    fluxFlags.inMu_flag = (muBeam != -1);  
    if (fluxFlags.inMu_flag && !cutCounts["Cut 03"].second) {
      cutCounts["Cut 03"].first++;
      cutCounts["Cut 03"].second = true;
    }

    const PaParticle & beam = e.vParticle(muBeam);
		int it_beam = beam.iTrack();
		if (it_beam == -1) continue; // the incoming particle has a track associated with it
    fluxFlags.inMuTrack_flag = (it_beam != -1);
    if (fluxFlags.inMuTrack_flag && !cutCounts["Cut 04"].second) {
      cutCounts["Cut 04"].first++;
      cutCounts["Cut 04"].second = true;
    }

    const PaTrack & beam_track = e.vTrack(muBeam);
		if (beam_track.NTPar() == 0) continue; // the track has parameters
    fluxFlags.inMuPar_flag = (beam_track.NTPar() != 0); 
    if (fluxFlags.inMuPar_flag && !cutCounts["Cut 05"].second) {
      cutCounts["Cut 05"].first++;
      cutCounts["Cut 05"].second = true;
    }

    if (beam_track.ZFirst() >= -78.5) continue; // incoming muon was first measured before the target
    fluxFlags.zFirst_flag = (beam_track.ZFirst() < -78.5);
    if (fluxFlags.zFirst_flag && !cutCounts["Cut 06"].second) {
      cutCounts["Cut 06"].first++;
      cutCounts["Cut 06"].second = true;
    }

    double inMu_mom = beam_track.vTPar(0).Mom();
    if (inMu_mom < 140.0 || inMu_mom > 180.0) continue; // momentum falls within acceptable range
    fluxFlags.momRange_flag = (inMu_mom >= 140.0 && inMu_mom <= 180.0);
    if (fluxFlags.momRange_flag && !cutCounts["Cut 07"].second) {
      cutCounts["Cut 07"].first++;
      cutCounts["Cut 07"].second = true;
    }
    
    double inMu_momErr = sqrt(beam_track.vTPar(0)(5,5))/(beam_track.vTPar(0)(5)*beam_track.vTPar(0)(5));
    if (inMu_momErr > 0.025*inMu_mom) continue; // momentum error falls within acceptable range  
    fluxFlags.momErr_flag = (inMu_momErr <= 0.025*inMu_mom);
    if (fluxFlags.momErr_flag && !cutCounts["Cut 08"].second) {
      cutCounts["Cut 08"].first++;
      cutCounts["Cut 08"].second = true;
    }

    // incoming muon is detected by detectors along the beamline 
    int nhits_BMS = beam_track.NHitsFoundInDetect("BM"); 
    int nhits_FI  = beam_track.NHitsFoundInDetect("FI"); 
		int nhits_SI  = beam_track.NHitsFoundInDetect("SI"); 

    if ((nhits_BMS < 3)) continue; 
    fluxFlags.BMS_flag = (nhits_BMS >= 3); 
    if (fluxFlags.BMS_flag && !cutCounts["Cut 09"].second) {
      cutCounts["Cut 09"].first++;
      cutCounts["Cut 09"].second = true;
    }

    if ((nhits_FI < 2)) continue; 
    fluxFlags.FI_flag = (nhits_FI >= 2);
    if (fluxFlags.FI_flag && !cutCounts["Cut 10"].second) {
      cutCounts["Cut 10"].first++;
      cutCounts["Cut 10"].second = true;
    }

		if ((nhits_SI < 3)) continue;
    fluxFlags.SI_flag = (nhits_SI >= 3);
    if (fluxFlags.SI_flag && !cutCounts["Cut 11"].second) {
      cutCounts["Cut 11"].first++;
      cutCounts["Cut 11"].second = true;
    }

    // PaAlgo::CrossCells(t_beam.vTPar(0),run, Rmax, Ymax, tgt_zmin, tgt_zmax, RmaxMC) - cut 4 starts here
    if (!(PaAlgo::CrossCells(beam_track.vTPar(0),Run, 1.9, 1.2, -318.5, -78.5, 2))) continue;
    fluxFlags.passTarget_flag = PaAlgo::CrossCells(beam_track.vTPar(0), Run, 1.9, 1.2, -318.5, -78.5, 2);
    if (fluxFlags.passTarget_flag && !cutCounts["Cut 12"].second) {
      cutCounts["Cut 12"].first++;
      cutCounts["Cut 12"].second = true;
    }

    double mean_time = beam_track.MeanTime();  
    if (std::fabs(mean_time) >= 2) continue; // track meantime for incoming muon is within flux requirements 
    fluxFlags.meantime_flag = (std::fabs(mean_time) < 2); 
    if (fluxFlags.meantime_flag && !cutCounts["Cut 13"].second) {
      cutCounts["Cut 13"].first++;
      cutCounts["Cut 13"].second = true;
     }

    if (!fluxFlags.TiS_flag) continue; 
    if (fluxFlags.TiS_flag && !cutCounts["Cut 14"].second) {
      cutCounts["Cut 14"].first++;
      cutCounts["Cut 14"].second = true;
    }

    // Check if all requirements for proper beam flux are fulfilled
    if (fluxFlags.allVtx_flag && fluxFlags.pVtx_flag && fluxFlags.inMu_flag &&
      fluxFlags.inMuTrack_flag && fluxFlags.inMuPar_flag && fluxFlags.passTarget_flag &&
      fluxFlags.zFirst_flag && fluxFlags.FI_flag && fluxFlags.SI_flag 
      && fluxFlags.BMS_flag && fluxFlags.momRange_flag && fluxFlags.momErr_flag &&
      fluxFlags.meantime_flag && fluxFlags.TiS_flag) {eventFlags.flux_flag = true;}

    const PaTPar & Par_beam = beam.ParInVtx(iv); // beam parameters at the vertex
    bool isInTarget = PaAlgo::InTarget(Par_beam, 'O', Run, 1.9, 1.2, -318.5, -78.5, 2);
    if (!isInTarget) continue; // vertex is in the target  
    eventFlags.vertex_flag = isInTarget;
    if (eventFlags.vertex_flag && !cutCounts["Cut 15"].second) {
      cutCounts["Cut 15"].first++;
      cutCounts["Cut 15"].second = true;
    }

    if (!eventFlags.trig_flag) continue; 
    if (eventFlags.trig_flag && !cutCounts["Cut 16"].second) {
      cutCounts["Cut 16"].first++;
      cutCounts["Cut 16"].second = true;
    }

    //*******************************************
    // Store info about scattered muon (outMu) 
    // HodoHelper->iMuPrim(v, checkYokeSM2, reject2muEvents, checkCanBeMuon, true, minXX0muPr, true, true) 
    int i_omu = -1; 
		i_omu = HodoHelper->iMuPrim(v,false,false,true,false,15); // index of the scattered muon WITHOUT CHECKING IF IT PASSES the hodoscope check 
		if (i_omu == -1) continue;
    int i_omu_check_hodo = HodoHelper->iMuPrim(v,false,true,true,true,15,true,true); // index for the scattered muon IF IT PASSES the hodoscope check 
    // if scattered muon passed the hodoscope use the corresponding index, if not proceed with other index 
    if (i_omu_check_hodo != -1) {
      eventFlags.hodo_flag = true;
			i_omu = i_omu_check_hodo;
		}

    if (!eventFlags.hodo_flag) continue; 
    if (eventFlags.hodo_flag && !cutCounts["Cut 17"].second) {
      cutCounts["Cut 17"].first++;
      cutCounts["Cut 17"].second = true;
    }


    const PaParticle & outMu = e.vParticle(i_omu);
    int outMu_itrack = outMu.iTrack();
    if (outMu_itrack == -1) continue; // outgoing muon has a track associated with it 

    const PaTrack & outMu_track = e.vTrack(outMu_itrack); 
    const PaTPar& Par_outMu = outMu.ParInVtx(iv); // scattered muon parameters at the vertex
    double outMu_mom = outMu_track.vTPar(0).Mom();    

    if (outMu.Q() != beam.Q()) continue; // scattered muon has the same charge as the beam
    eventFlags.charge_flag = (outMu.Q() == beam.Q());
    if (eventFlags.charge_flag && !cutCounts["Cut 18"].second) {
      cutCounts["Cut 18"].first++;
      cutCounts["Cut 18"].second = true;
    }
    
    if (!(outMu_track.ZFirst() < 350 && outMu_track.ZLast() > 350)) continue; // first and last z coordinates are measured before and after SM1
    eventFlags.zFirstLast_flag = (outMu_track.ZFirst() < 350 && outMu_track.ZLast() > 350);
    if (eventFlags.zFirstLast_flag && !cutCounts["Cut 19"].second) {
      cutCounts["Cut 19"].first++;
      cutCounts["Cut 19"].second = true;
    }

    //*******************************************
    // Kinematic variables ... (1/2)
    TLorentzVector inMu_TL  = Par_beam.LzVec(M_mu); 
	  TLorentzVector outMu_TL = Par_outMu.LzVec(M_mu); 
    TLorentzVector targ_TL(0,0,0,M_p);
    TLorentzVector q = (inMu_TL - outMu_TL); // four momentum of the virtual photon

    Q2  = PaAlgo::Q2 (inMu_TL, outMu_TL); //Q2  = -(inMu_TL - outMu_TL).M2();
    y   = (inMu_TL.E() - outMu_TL.E()) / inMu_TL.E();
    nu  = (inMu_TL.E() - outMu_TL.E());
    W2  = PaAlgo::W2 (inMu_TL, outMu_TL);
    xbj = PaAlgo::xbj (inMu_TL, outMu_TL); //xbj = Q2/(2*q*targ_TL);

    // Current kinematic cuts are wide - they will be tightened after the kinematically constrained fit is applied 
    if (Q2 < 0.8 || Q2 > 10) continue;
    eventFlags.Q2_flag = (Q2 > 0.8 && Q2 < 10);
    if (eventFlags.Q2_flag && !cutCounts["Cut 20"].second) {
      cutCounts["Cut 20"].first++;
      cutCounts["Cut 20"].second = true;
    }

    if (y < 0.05 || y > 0.9) continue; 
    eventFlags.y_flag = (y > 0.05 && y < 0.9);
    if (eventFlags.y_flag && !cutCounts["Cut 21"].second) {
      cutCounts["Cut 21"].first++;
      cutCounts["Cut 21"].second = true;
    }

    // Check that DIS flags are satisfied so far 
    if (eventFlags.flux_flag && eventFlags.hodo_flag && eventFlags.vertex_flag && eventFlags.charge_flag 
      && eventFlags.zFirstLast_flag && eventFlags.trig_flag && eventFlags.Q2_flag && eventFlags.y_flag) 
      {eventFlags.DIS_flag = true;}

    //*******************************************
    // Exclusive selection starts here  ... 
    // Only one outgoing particle (scattered proton and photon are detected using ECals and CAMERA so will not be found here)
    if (v.NOutParticles()!= 1) continue; 
    eventFlags.singleTrack_flag = (v.NOutParticles() == 1); 
    if (eventFlags.singleTrack_flag && !cutCounts["Cut 22"].second) {
      cutCounts["Cut 22"].first++;
      cutCounts["Cut 22"].second = true;
    }

    if (!eventFlags.flux_flag || !eventFlags.DIS_flag) continue;

    //*******************************************
    // Store info about real photon - check ECals for single neutral cluster
    // ! CURRENTLY LOOKING FOR ANY NEUTRAL CLUSTERS, NOT SPECIFICALLY FOR ONE !
    int ecal0id = PaSetup::Ref().iCalorim("EC00P1__");
    int ecal1id = PaSetup::Ref().iCalorim("EC01P1__");
    int ecal2id = PaSetup::Ref().iCalorim("EC02P1__");

    std::vector<PaCaloClus> Clus_gamma;  // Store valid clusters dynamically
    std::vector<int> cl_id;  // Store calorimeter IDs for valid clusters

    int clusterCount = 0;

    for (int iclus = 0; iclus < e.NCaloClus(); iclus++) { // Begin loop over neutral clusters
      const PaCaloClus & cl = e.vCaloClus(iclus);
      int nam = cl.iCalorim();

      if ((nam != ecal0id) && (nam != ecal1id) && (nam != ecal2id)) continue; // cluster needs to be detected in one of the ECals
      if (cl.iTrack() != -1) continue; // cluster needs to have an associated track  

      if (Year == 2016 && !EcalTimeCut(beam_track, cl, Run)) continue; // check timing of the clusters

      if ((nam == ecal0id && cl.E() >= 4) || // Energy threshold check (above DVCS energy threshold)
          (nam == ecal1id && cl.E() >= 5) ||
          (nam == ecal2id && cl.E() >= 10)) {
          Clus_gamma.push_back(cl); // Store the valid cluster
          cl_id.push_back(nam);     // Store the ECal ID 
          clusterCount++; 
      }

    } // End loop over neutral clusters 
    eventFlags.checkCls_flag = true;  // Set flag to true after processing all clusters 

    // ! NOW CHECK THAT THERE IS ONLY ONE NEUTRAL CLUSTER !
    if(clusterCount != 1) continue;  
    if (clusterCount == 1) {eventFlags.singleCl_flag = true;} 
    if (eventFlags.singleCl_flag && !cutCounts["Cut 23"].second) {
      cutCounts["Cut 23"].first++;
      cutCounts["Cut 23"].second = true;
    }

    // Calculate TL vector for all photons in event (possible exclusive ones and low energy background) 
    std::vector<TLorentzVector> gamma_TLs; // Store TL vectors dynamically
    for (size_t i = 0; i < Clus_gamma.size(); ++i) { // Begin loop over cluster
      const auto& cluster = Clus_gamma[i];
      int ecalid = cl_id[i]; // Retrieve cluster ID from cl_id

      double phdz  = cluster.Z() - v.Z();
      double phdy  = cluster.Y() - v.Y();
      double phdx  = cluster.X() - v.X();
      double phdxy = TMath::Sqrt(phdx * phdx + phdy * phdy);
      double phL   = TMath::Sqrt(phdxy * phdxy + phdz * phdz);
      double cl_E  = cluster.E();

      TLorentzVector gamma_TL(cl_E*phdx/phL, cl_E*phdy/phL, cl_E*phdz/phL, cl_E);
      gamma_TLs.emplace_back(gamma_TL); // Store the TL vectors

      // Fill the reconstructed energy distributions for each ECal
      if (ecalid == 0) h97_gamma_E_EC0->Fill(cl_E);
      if (ecalid == 1) h97_gamma_E_EC1->Fill(cl_E);
      if (ecalid == 2) h97_gamma_E_EC2->Fill(cl_E);
    } // End loop over cluster

    //*******************************************
    // Kinematic variables ... (2/2)
    t = (q - gamma_TLs[0]).M2();
    E_miss = nu - gamma_TLs[0].E() + t/(2 * M_p);

    // Check that exclusive cut conditions are satisfied before checking camera
    // TODO: add exclusiveCl_flag check here as well 
    if (eventFlags.singleTrack_flag && eventFlags.singleCl_flag) {eventFlags.exclEvent_flag = true;}

    //*******************************************
    printDebug("      ");
    printDebug("*** Run: " + std::to_string(Run) + ", spill: " + std::to_string(Spill) + ", event: " + std::to_string(EvtInSpill) + " ***");
    printDebug("    Event Count: 1: " + std::to_string(cutCounts["Cut 01"].first) + ", 2: " + std::to_string(cutCounts["Cut 02"].first) + 
        ", 3: " + std::to_string(cutCounts["Cut 03"].first) + ", 4: " + std::to_string(cutCounts["Cut 04"].first) + 
        ", 5: " + std::to_string(cutCounts["Cut 05"].first) + ", 6: " + std::to_string(cutCounts["Cut 06"].first) +
        ", 7: " + std::to_string(cutCounts["Cut 07"].first) + ", 8: " + std::to_string(cutCounts["Cut 08"].first) +
        ", 9: " + std::to_string(cutCounts["Cut 09"].first) + ", 10: " + std::to_string(cutCounts["Cut 10"].first) +
        ", 11: " + std::to_string(cutCounts["Cut 11"].first) + ", 12: " + std::to_string(cutCounts["Cut 12"].first) +
        ", 13: " + std::to_string(cutCounts["Cut 13"].first) + ", 14: " + std::to_string(cutCounts["Cut 14"].first) +
        ", 15: " + std::to_string(cutCounts["Cut 15"].first) + ", 16: " + std::to_string(cutCounts["Cut 16"].first) +
        ", 17: " + std::to_string(cutCounts["Cut 17"].first) + ", 18: " + std::to_string(cutCounts["Cut 18"].first) +
        ", 19: " + std::to_string(cutCounts["Cut 19"].first) + ", 20: " + std::to_string(cutCounts["Cut 20"].first) +
        ", 21: " + std::to_string(cutCounts["Cut 21"].first) + ", 22: " + std::to_string(cutCounts["Cut 22"].first) +
        ", 23: " + std::to_string(cutCounts["Cut 23"].first) + ", 24: " + std::to_string(cutCounts["Cut 24"].first));

    //*******************************************
    // Fill the tree/histograms with the extracted event information 
    int bin = 1;
    for (const auto &cut : cutCounts) {
        h97_cutStats->SetBinContent(bin, cut.second.first);  // Use cut.second.first to get the count
        h97_cutStats->GetXaxis()->SetBinLabel(bin, cut.first.c_str());  // Set cut labels
        bin++;
    }

    inMu_pz = beam_track.vTPar(0).Pz();
    inMu_py = beam_track.vTPar(0).Py();
    inMu_px = beam_track.vTPar(0).Px();
    inMu_p  = beam_track.vTPar(0).Mom();
    inMu_E  = inMu_TL.E();
    //inMu_theta = acos(inMu_pz/inMu_p);
    //inMu_phi   = atan2(inMu_py, inMu_px);

    outMu_pz = outMu_track.vTPar(0).Pz();
    outMu_py = outMu_track.vTPar(0).Py();
    outMu_px = outMu_track.vTPar(0).Px();
    outMu_p  = outMu_track.vTPar(0).Mom();
    outMu_E  = outMu_TL.E();
    //outMu_theta = acos(outMu_pz/outMu_p);
    //outMu_phi   = atan2(outMu_py, outMu_px);

    tree->Fill();
    h97_Zprim->Fill(Zprim); 
    h97_Yprim->Fill(Yprim);
    h97_Xprim->Fill(Xprim);
    h97_XYprim->Fill(Xprim, Yprim);

    h97_inMu_p->Fill(inMu_p);
    h97_inMu_px->Fill(inMu_px);
    h97_inMu_py->Fill(inMu_py);
    //h97_inMu_theta->Fill(inMu_theta);
    //h97_inMu_phi->Fill(inMu_phi);

    h97_outMu_p->Fill(outMu_p);
    h97_outMu_px->Fill(outMu_px);
    h97_outMu_py->Fill(outMu_py);
    //h97_outMu_theta->Fill(outMu_theta);
    //h97_outMu_phi->Fill(outMu_phi);

    h97_y->Fill(y);
    h97_Q2->Fill(Q2);
    h97_W2->Fill(W2);
    h97_xbj->Fill(xbj); 
    h97_Q2xbj->Fill(xbj,Q2);

    h97_E_miss->Fill(E_miss); 

  } // end of loop over vertices

  //if (saveEvent_flag) {e.TagToSave();}
  e.TagToSave();

}


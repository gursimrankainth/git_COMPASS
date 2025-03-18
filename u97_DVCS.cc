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

// ************************************************************************** //
// UserEvent for preselecting exclusive photon events	(DVCS)			            //
// In the selection all possible combinations of:					                    //
// Vertices (incoming and outgoing muons), exclusive photon and recoil proton	//
// ************************************************************************** //

extern "C" float prob_(float&, int&);

// Global flag for verbose mode 
bool verbose_mode = false; // Set to true for verbose output, false to suppress
void printDebug(const std::string &message, bool forcePrint = false) {
    if (verbose_mode || forcePrint) {
        std::cout << message << std::endl;
    }
}

//*****************************************************************************
void UserEvent97(PaEvent & e) {

  // Define constants
  static PaHodoHelper* HodoHelper = NULL;
  static TiSRange*     tis_range  = NULL; 
  //static PaCamera*     cam_inst;

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
  
  static int    Hodo_i_omu;
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
  struct EventFlags {
    //DIS cuts
    bool flux_flag       = false;
    bool hodo_flag       = false;
    bool outMuTrack_flag = false;
    bool vertex_flag     = false;
    bool realoutMu_flag  = false;
    bool charge_flag     = false;
    bool zFirstLast_flag = false;
    bool trig_flag       = false;
    bool Q2_flag         = false;
    bool y_flag          = false;
    bool DIS_flag        = false; // true if all DIS cuts have been applied 
    //Exclusive cuts 
    bool singleTrack_flag = false;
    bool checkCls_flag    = false; // true if all found clusters have a track, and pass the timing and ECal energy cuts 
    bool singleCl_flag    = false; // true after a single neutral cluster has been found in the ECals
    bool exclEvent_flag   = false; // true if all exclusive cuts have ben applied 
    bool saveEvent_flag   = false;
  };

  struct FluxFlags {
    bool allVtx_flag     = false;
    bool pVtx_flag       = false;
    bool inMu_flag       = false;
    bool inMuTrack_flag  = false;
    bool inMuPar_flag    = false;
    bool passTarget_flag = false;
    bool zFirst_flag     = false;
    bool BMS_flag        = false;
    bool FI_flag         = false;
    bool SI_flag         = false;
    bool momRange_flag   = false;
    bool momErr_flag     = false;
    bool meantime_flag   = false;
    bool TiS_flag        = false;
  };

  EventFlags eventFlags;   // Create an instance of event flags
  FluxFlags  fluxFlags;    // Create an instance of event flags

  // Register XCheck flags to track cut statistics 
  statis bool xcheck_mode = false; 
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
  XCHECK_REGISTER_FLAG(Event_Hodo, "No. of events where scattered muon passes Hodoscope check");
  XCHECK_REGISTER_FLAG(Event_InTarget, "No. of events where the vertex is in the target");
  XCHECK_REGISTER_FLAG(Event_Trigger, "No. of events with MT, LT, OT or LAST physics triggers");
  XCHECK_REGISTER_FLAG(Event_RealoutMu, "No. of events where the scattered muon actually exists");
  XCHECK_REGISTER_FLAG(Event_Charge, "No. of events where scattered muon has the same charge as the beam");
  XCHECK_REGISTER_FLAG(Event_ZFirstLast, "No. of events where first and last scattered muon z coord. are measured before and after SM1");
  XCHECK_REGISTER_FLAG(Event_Q2, "No. of events where 1 < Q2 < 10");
  XCHECK_REGISTER_FLAG(Event_Y, "No. of events where 0.05 < y < 0.9");
  XCHECK_REGISTER_FLAG(Event_TrackMult, "No. of events where primary vertex only has one outgoing track");
  XCHECK_REGISTER_FLAG(Event_ClustMult, "No. of events where there is only a single neutral clsuter in the ECals");

  //*****************************************************************************

  static bool first(true);
  if (first) { // histograms and Ntuples booking block
    Phast::Ref().HistFileDir("UserEvent97");
   
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
  if (trig_mask & LAST) {
    trigCheck += "LAST ";
    eventFlags.trig_flag = true;
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

  fluxFlags.TiS_flag = tis_range->CheckWindow(Run, Spill, TimeInSpill); // check the time in spill 

  //*******************************************  
  XCHECK_COUNT_FLAG(Event_AllEvents, "Total no. of events processed by PHAST user script"); 

  // Loop over reconstructed vertices in the event 
  //std::cout << std::endl << "DEBUG :: " << Evt << ", " << e.NVertex() << std::endl;
  for (int iv = 0; iv < e.NVertex(); iv++) { // begin loop over vertices
  fluxFlags.allVtx_flag = true;
    //******************************************* 
    // Store info about primary vertex (if found) 
    const PaVertex & v = e.vVertex(iv);
    if (!v.IsPrimary()) continue;
    //std::cout << std::endl << "DEBUG :: " << Evt << std::endl;
    fluxFlags.pVtx_flag = v.IsPrimary();
    XCHECK_COUNT_FLAG(Event_PVtx, "No. of events with a primary vertex");
    Zprim   = v.Pos(2);
    Yprim   = v.Pos(1);
    Xprim   = v.Pos(0);
    Nprim   = v.NOutParticles(); // number of tracks in vertex 

    //*******************************************
    // Store info about incoming muon beam (inMu)
    int i_beam = v.InParticle(); 
    if (i_beam == -1) continue; // there is incoming particle associated with the primary vertex
    fluxFlags.inMu_flag = (i_beam != -1);  

    const PaParticle & beam = e.vParticle(i_beam);
		int it_beam = beam.iTrack();
		if (it_beam == -1) continue; // the incoming particle has a track associated with it
    fluxFlags.inMuTrack_flag = (it_beam != -1);

    const PaTrack & beam_track = e.vTrack(i_beam);
		if (beam_track.NTPar() == 0) continue; // the track has parameters
    fluxFlags.inMuPar_flag = (beam_track.NTPar() != 0); 
    XCHECK_COUNT_FLAG(Event_BeamTrack, "No. of events where beam has a track with parameters");

    if (beam_track.ZFirst() >= -78.5) continue; // incoming muon was first measured before the target
    fluxFlags.zFirst_flag = (beam_track.ZFirst() < -78.5);
    XCHECK_COUNT_FLAG(Event_ZFirst, "No. of events where beam was first measured before the target");

    double inMu_mom = beam_track.vTPar(0).Mom();
    if (inMu_mom < 140.0 || inMu_mom > 180.0) continue; // momentum falls within acceptable range
    fluxFlags.momRange_flag = (inMu_mom >= 140.0 && inMu_mom <= 180.0);
    XCHECK_COUNT_FLAG(Event_BeamMom, "No. of events where beam momentum falls within acceptable range");
    
    double inMu_momErr = sqrt(beam_track.vTPar(0)(5,5))/(beam_track.vTPar(0)(5)*beam_track.vTPar(0)(5));
    if (inMu_momErr > 0.025*inMu_mom) continue; // momentum error falls within acceptable range  
    fluxFlags.momErr_flag = (inMu_momErr <= 0.025*inMu_mom);
    XCHECK_COUNT_FLAG(Event_BeamMomErr, "No. of events where beam momentum error falls within acceptable range");

    // incoming muon is detected by detectors along the beamline  
    int nhits_BMS = beam_track.NHitsFoundInDetect("BM");
    int nhits_FI  = beam_track.NHitsFoundInDetect("FI"); 
		int nhits_SI  = beam_track.NHitsFoundInDetect("SI"); 

    if ((nhits_BMS < 3)) continue; 
    fluxFlags.BMS_flag = (nhits_BMS >= 3); 
    XCHECK_COUNT_FLAG(Event_BMS, "No. of events where beam is detected by BMS");

    if ((nhits_FI < 2)) continue; 
    fluxFlags.FI_flag = (nhits_FI >= 2);
    XCHECK_COUNT_FLAG(Event_FI, "No. of events where beam is detected by SCIFI");

		if ((nhits_SI < 3)) continue;
    fluxFlags.SI_flag = (nhits_SI >= 3);
    XCHECK_COUNT_FLAG(Event_SI, "No. of events where beam is detected by SI");

    // PaAlgo::CrossCells(t_beam.vTPar(0),run, Rmax, Ymax, tgt_zmin, tgt_zmax, RmaxMC) 
    const PaTPar & Par_beam = beam.ParInVtx(iv); // beam parameters at the vertex
    if (!(PaAlgo::CrossCells(beam_track.vTPar(0), Run, 1.9, 1.2, -318.5, -78.5, 2.0))) continue;
    fluxFlags.passTarget_flag = PaAlgo::CrossCells(beam_track.vTPar(0), Run, 1.9, 1.2, -318.5, -78.5, 2.0);
    XCHECK_COUNT_FLAG(Event_CrossCells, "No. of events where the beam crosses full target length");

    double mean_time = beam_track.MeanTime();  
    if (std::fabs(mean_time) >= 2) continue; // track meantime for incoming muon is within flux requirements 
    fluxFlags.meantime_flag = (std::fabs(mean_time) < 2); 
    XCHECK_COUNT_FLAG(Event_Meantime, "No. of events where beam track meantime is within acceptable flux requirements");

    if (!fluxFlags.TiS_flag) continue; 
    XCHECK_COUNT_FLAG(Event_TiS, "No. of events where time in spill is within acceptable flux requirements");

    // Check if all requirements for proper beam flux are fulfilled
    if (fluxFlags.allVtx_flag && fluxFlags.pVtx_flag && fluxFlags.inMu_flag &&
      fluxFlags.inMuTrack_flag && fluxFlags.inMuPar_flag && fluxFlags.passTarget_flag &&
      fluxFlags.zFirst_flag && fluxFlags.FI_flag && fluxFlags.SI_flag 
      && fluxFlags.BMS_flag && fluxFlags.momRange_flag && fluxFlags.momErr_flag &&
      fluxFlags.meantime_flag && fluxFlags.TiS_flag) {eventFlags.flux_flag = true;}

    //*******************************************
    // Store info about scattered muon (outMu) 
    // HodoHelper->iMuPrim(v, checkYokeSM2, reject2muEvents, checkCanBeMuon, true, minXX0muPr, true, true) 
    int i_omu = -1; 
		i_omu = HodoHelper->iMuPrim(v,false,false,true,false,15); // index of the scattered muon WITHOUT CHECKING IF IT PASSES the hodoscope check 
    if (i_omu == -1) continue;  

    if(!(PaAlgo::InTarget(Par_beam,'O',Run, 1.9, 1.2, -318.5, -78.5, 2.0))) continue;    
    eventFlags.vertex_flag = PaAlgo::InTarget(Par_beam,'O',Run, 1.9, 1.2, -318.5, -78.5, 2.0);
    XCHECK_COUNT_FLAG(Event_InTarget, "No. of events where the vertex is in the target");

    if (!eventFlags.trig_flag) continue; 
    XCHECK_COUNT_FLAG(Event_Trigger, "No. of events with MT, LT, OT or LAST physics triggers");

    int i_omu_check_hodo = HodoHelper->iMuPrim(v,false,true,true,true,15,true,true); // index for the scattered muon IF IT PASSES the hodoscope check 
    // if scattered muon passed the hodoscope use the corresponding index, if not proceed with other index 
    if (i_omu_check_hodo != -1) {
      eventFlags.hodo_flag = true;
			i_omu = i_omu_check_hodo;
		}
    XCHECK_COUNT_FLAG(Event_Hodo, "No. of events where scattered muon passes Hodoscope check");

    const PaParticle & outMu = e.vParticle(i_omu);
    // if outMu.Q or beam.Q return -777 it means the assocaited track was reconstructed in a field free region (charge is unkown)
    if (outMu.Q() != beam.Q() || outMu.Q() == -777 || beam.Q() == -777) continue; // scattered muon has the same charge as the beam
    eventFlags.charge_flag = (outMu.Q() == beam.Q() && outMu.Q() != -777 && beam.Q()!= -777);
    XCHECK_COUNT_FLAG(Event_Charge, "No. of events where scattered muon has the same charge as the beam");

    int outMu_itrack = outMu.iTrack();
    outMu_track = e.vTrack(outMu_itrack);
    Par_outMu = outMu.ParInVtx(vertexIndex);
    if (!(outMu_track.ZFirst() < 350 && outMu_track.ZLast() > 350)) continue; // first and last z coordinates are measured before and after SM1
    eventFlags.zFirstLast_flag = (outMu_track.ZFirst() < 350 && outMu_track.ZLast() > 350);
    XCHECK_COUNT_FLAG(Event_ZFirstLast, "No. of events where first and last scattered muon z coord. are measured before and after SM1");

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
    if (Q2 < 1 || Q2 > 10) continue;
    eventFlags.Q2_flag = (Q2 > 1 && Q2 < 10);
    XCHECK_COUNT_FLAG(Event_Q2, "No. of events where 1 < Q2 < 10");

    if (y < 0.05 || y > 0.9) continue; 
    eventFlags.y_flag = (y > 0.05 && y < 0.9);
    XCHECK_COUNT_FLAG(Event_Y, "No. of events where 0.05 < y < 0.9");

    // Check that DIS flags are satisfied so far 
    if (eventFlags.flux_flag && eventFlags.hodo_flag && eventFlags.vertex_flag && eventFlags.charge_flag 
      && eventFlags.zFirstLast_flag && eventFlags.trig_flag && eventFlags.Q2_flag && eventFlags.y_flag) 
      {eventFlags.DIS_flag = true;}

    //*******************************************
    // Exclusive selection starts here  ... 
    // Only one outgoing particle (scattered proton and photon are detected using ECals and CAMERA so will not be found here)
    if (v.NOutParticles()!= 1) continue; 
    eventFlags.singleTrack_flag = (v.NOutParticles() == 1); 
    XCHECK_COUNT_FLAG(Event_TrackMult, "No. of events where primary vertex only has one outgoing track");

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
    XCHECK_COUNT_FLAG(Event_ClustMult, "No. of events where there is only a single neutral clsuter in the ECals");

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
    t       = (q - gamma_TLs[0]).M2();
    E_miss  = nu - gamma_TLs[0].E() + t/(2 * M_p); // Missing energy (assuming proton)
    M2_miss = (2 * M_p * (nu - gamma_TLs[0].E())) + (M_p * M_p) + t; // Missing mass squared (assuming proton)

    //*******************************************
    // Store info about proton using Camera and calculate exclusivity variables (used for kinematic fit)
    //vector <CameraProton> protons = cam_inst->GetGoodCandidates(v);

    //*******************************************
    // Check that exclusive cut conditions are satisfied after checking camera
    // add additional flags here 
    if (eventFlags.singleTrack_flag && eventFlags.singleCl_flag) {eventFlags.exclEvent_flag = true;}

    //*******************************************
    printDebug("      ");
    printDebug("*** Run: " + std::to_string(Run) + ", spill: " + std::to_string(Spill) + ", event: " + std::to_string(Evt) + " ***");
    printDebug("    Vertex: (" + std::to_string(Xprim) + ", " + std::to_string(Yprim) + ", " + std::to_string(Zprim) + ")");
    printDebug("    mu: P: " + std::to_string(inMu_mom) + " GeV/c, Charge: " + std::to_string(beam.Q()));
    printDebug("    mu': P: " + std::to_string(outMu_mom) + " GeV/c, Charge: " + std::to_string(outMu.Q()));
    printDebug("    Kinematics: Q2: " + std::to_string(Q2) +  " GeV2, y: " + std::to_string(y) + ", W2: " + std::to_string(W2) + " GeV2, x: " + std::to_string(xbj));

    //*******************************************
    inMu_pz = beam_track.vTPar(0).Pz();
    inMu_py = beam_track.vTPar(0).Py();
    inMu_px = beam_track.vTPar(0).Px();
    inMu_p  = beam_track.vTPar(0).Mom();
    inMu_E  = inMu_TL.E();
    //inMu_theta = acos(inMu_pz/inMu_p); // how to get theta 
    //inMu_phi   = atan2(inMu_py, inMu_px); // how to get phi 

    outMu_pz = outMu_track.vTPar(0).Pz();
    outMu_py = outMu_track.vTPar(0).Py();
    outMu_px = outMu_track.vTPar(0).Px();
    outMu_p  = outMu_track.vTPar(0).Mom();
    outMu_E  = outMu_TL.E();

    tree->Fill();
    h97_Zprim->Fill(Zprim); 
    h97_Yprim->Fill(Yprim);
    h97_Xprim->Fill(Xprim);
    h97_XYprim->Fill(Xprim, Yprim);

    h97_inMu_p->Fill(inMu_p);
    h97_inMu_px->Fill(inMu_px);
    h97_inMu_py->Fill(inMu_py);

    h97_outMu_p->Fill(outMu_p);
    h97_outMu_px->Fill(outMu_px);
    h97_outMu_py->Fill(outMu_py);

    h97_y->Fill(y);
    h97_Q2->Fill(Q2);
    h97_W2->Fill(W2);
    h97_xbj->Fill(xbj); 
    h97_Q2xbj->Fill(xbj,Q2);

    h97_E_miss->Fill(E_miss); 
    h97_M2_miss->Fill(M2_miss);

  } // end of loop over vertices

  //if (saveEvent_flag) {e.TagToSave();}
  e.TagToSave();

}

void UserJobEnd97() {
  XCHECK_JOB_END(); // Print to output stream
}


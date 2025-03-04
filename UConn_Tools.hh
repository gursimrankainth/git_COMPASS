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

// ************************************************************** //
// Header script containing functions used for event selection.   //
// Notes for each functon:  					                  //
// 1. Name and use case                                           //
// 2. Inputs (NEED TO BE DEFINED BY USER)                         //
// 3. Output                                                      //
// 4. Example usage                                               //
// ************************************************************** //

// *************************  BEAMFLUX  *************************** 
// 1. Function: beamFlux -> check that all of the flux requirements are satisfied by the incoming muon beam
// 2. Input: PaEvent object, PaVertex object, Run, TiS_flag, BeamFluxParams object 
// 3. Output: flag with a boolean value that is true for events that pass the check and false for events that do not 
// 4. Example usage: ... 
"""
# include ... 
extern "C" float prob_(float&, int&);

void UserEventxxx(PaEvent & e) {

    static int    Run;         
    static int    LastRun = -1;     
    static int    Evt;         
    static int    Spill;       
    static double TimeInSpill; 
    static bool   TiS_flag = false;  

    static BeamFluxParams params; // Create an instance of BeamFluxParams
    Run = e.RunNum();
    if (Run != LastRun) {  
        // Set_TiSrange("PATH_TO_FLUX_FILES", runMin, runMax);
        Set_TiSrange("PATH_TO_FLUX_FILES", Run, Run);
        LastRun = Run;  // Update LastRun to the current run number
    }

    Evt          = e.UniqueEvNum();
    Spill        = e.SpillNum(); 
    TimeInSpill  = e.TimeInSpill(); 
    TiS_flag     = Check_TiS_window(Run,Spill,TimeInSpill);

    for (int iv = 0; iv < e.NVertex(); iv++) { // loop over all vertices 
        const PaVertex & v = e.vVertex(iv);
        if (!v.IsPrimary()) continue; // skip vertices that are not primary
        flux_flag = beamFlux(e, v, Run, TiS_flag, params); 
    }

    if (flux_flag) {e.TagToSave();} // only save events that pass the beamFlux check
}
"""

// Define a struct to hold all the parameters
struct BeamFluxParams {
    int Run;
    double Rmax;
    double Ymax;
    double tgt_zmin;
    double tgt_zmax;
    int RmaxMC;
    double zfirst;  
    int minFI;
    int minSI;
    int minBMS;
    double minMom; 
    double maxMom; 
    double percent; // use this to determine if error in momentum is within acceptable range 

    // Constructor with default values
    BeamFluxParams(double rmax = 1.9, double ymax = 1.2, double zmin = -318.5, 
                   double zmax = -78.5, int rmaxMC = 2, double zfirst = -78.5, 
                   int minFI = 2, int minSI = 3, int minBMS = 3, double minMom = 14.0,
                   double maxMom = 180.0, double percent = 0.025) 
        : Rmax(rmax), Ymax(ymax), tgt_zmin(zmin), tgt_zmax(zmax), RmaxMC(rmaxMC),
          zfirst(zfirst), minFI(minFI), minSI(minSI), minBMS(minBMS), minBMS(minMom),
          maxMom(maxMom), percent(percent) {}
};

// Define the function 
bool beamFlux(const PaEvent &e, const PaVertex &v, int Run , bool TiS_flag, const BeamFluxParams &params) { // beamFlux loop begins 
    // Initialize fluxFlag to true (assuming all conditions will be met unless one fails)
    bool fluxFlag = true;

    // Check that there is an incoming particle associated with the vertex
    int muBeam = v.InParticle();
    if (muBeam == -1) {
        fluxFlag = false;
    }

    // Check that the incoming particle has a track associated with it
    int it_beam = beam.iTrack();
    if (it_beam == -1) {
        fluxFlag = false;
    }

    // Check that the track has parameters
    const PaTrack &beam_track = e.vTrack(muBeam);
    if (beam_track.NTPar() == 0) {
        fluxFlag = false;
    }

    // Check that beam crosses full length of target
    if (!(PaAlgo::CrossCells(beam_track.vTPar(0), Run, params.Rmax, params.Ymax, params.tgt_zmin, params.tgt_zmax, params.RmaxMC))) {
        fluxFlag = false;
    }

    // Check that the beam was first measured before the target 
    if (beam_track.ZFirst() >= params.zfirst) {
        fluxFlag = false;
    }

    // Check that beam is detected by beamline detectors   
    int nhits_FI  = beam_track.NHitsFoundInDetect("FI"); 
	int nhits_SI  = beam_track.NHitsFoundInDetect("SI"); 
    int nhits_BMS = beam_track.NHitsFoundInDetect("BM");
    if (nhits_FI < params.minFI && nhits_SI < params.minSI && nhits_BMS < params.BMS) {
        fluxFlag = false;
    }

    // Check that the momentum falls within an acceptable range 
    double inMu_mom = beam_track.vTPar(0).Mom();
    if (inMu_mom < params.minMom && inMu_mom > params.maxMom) {
        fluxFlag = false;
    }

    // Check that the momentum error falls within an acceptable range 
    double inMu_momErr = sqrt(beam_track.vTPar(0)(5,5))/(beam_track.vTPar(0)(5)*beam_track.vTPar(0)(5));
    if (inMu_momErr > params.percent*inMu_mom) {
        fluxFlag = false;
    }

    // Passes time in spill window requirements 
    if (!TiS_flag) {
        fluxFlag = false;
    }

    // Vertex is in the target 
    const PaTPar & Par_beam = beam.ParInVtx(iv);
    bool isInTarget = PaAlgo::InTarget(Par_beam, 'O', Run, params.Rmax, params.Ymax, params.tgt_zmin, params.tgt_zmax, params.RmaxMC);
    if (!isInTarget) {
        fluxFlag = false;
    }

    // Return the value of fluxFlag
    return fluxFlag;
} // beamFlux loop begins


// *************************  ...................  *************************** 
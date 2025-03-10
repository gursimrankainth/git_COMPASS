#pragma once
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

// *************************  PRINTDEBUG  *************************** 
// 1. Function: printDebug -> print debug statements if verbose_mode is set to true 
// 2. Input: N/A, see below how to use  
// 3. Output: prints statements to console while PHAST us running for involved debugging 
// 4. Example usage: ... 

//verbose_mode = true; // set this manually in your own script for ease of use 
//printDebug("*** Run: " + std::to_string(Run) + ", spill: " + std::to_string(Spill) + ", event: " + std::to_string(EvtInSpill) + " ***");
// OR
//printDebug("test", true); // if second arguement is set to true it will print even if verbose_mode is off 

// Global flag for verbose mode 
bool verbose_mode = false; // Set to true for verbose output, false to suppress

// Define the function (to print statement even if verbose mode is off set second arguement to true)
void printDebug(const std::string &message, bool forcePrint = false) {
    if (verbose_mode || forcePrint) {
        std::cout << message << std::endl;
    }
}

// *************************  BEAMFLUXCHECK  *************************** 
// 1. Function: beamFlux -> check that all of the flux requirements are satisfied by the incoming muon beam
// 2. Input: PaEvent object, PaVertex object, int vertex index, Run, TiS_flag, BeamFluxParams object, xcheck_mode
// 3. Output: flag with a boolean value that is true for events that pass the check and false for events that do not 
// 4. Example usage: ... 

//TODO: example user event here  

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
          zfirst(zfirst), minFI(minFI), minSI(minSI), minBMS(minBMS), minMom(minMom),
          maxMom(maxMom), percent(percent) {}
};

// Define the function 
bool beamFluxCheck(const PaEvent & e, const PaVertex &v, int vertexIndex, int Run , bool TiS_flag, const BeamFluxParams &params) { // beamFlux loop begins 
	
	// Check that there is an incoming muon associated with the vertex
	int i_beam = v.InParticle(); 
	if (i_beam == -1) {
			return false;
	}

	// Check that the beam has a track associated with it
	const PaParticle & beam = e.vParticle(i_beam);
	int it_beam = beam.iTrack();
	if (it_beam == -1) {
			return false;
	}

	// Check that the track has parameters
	const PaTrack & beam_track = e.vTrack(i_beam);
	if (beam_track.NTPar() == 0) {
			return false;
	}

	// Check that the beam was first measured before the target
	if (beam_track.ZFirst() >= -78.5) {
			return false;
	}

	// Check that the beam momentum falls within acceptable range
	double inMu_mom = beam_track.vTPar(0).Mom();
	if (inMu_mom < 140.0 || inMu_mom > 180.0) {
		return false;
	}

	// Check that the beam momentum error falls within acceptable range
	double inMu_momErr = sqrt(beam_track.vTPar(0)(5,5))/(beam_track.vTPar(0)(5)*beam_track.vTPar(0)(5));
  if (inMu_momErr > 0.025*inMu_mom) {
		return false;
	}

	// Check that the beam is detected by detectors along the beamline
	int nhits_BMS = beam_track.NHitsFoundInDetect("BM");
	int nhits_FI  = beam_track.NHitsFoundInDetect("FI"); 
	int nhits_SI  = beam_track.NHitsFoundInDetect("SI"); 
	if (nhits_BMS < 3 || nhits_FI < 2 || nhits_SI < 3) {
    return false;
	}

	// Check that the beam crosses the full target length 
	// PaAlgo::CrossCells(t_beam.vTPar(0),run, Rmax, Ymax, tgt_zmin, tgt_zmax, RmaxMC) 
	const PaTPar & Par_beam = beam.ParInVtx(vertexIndex); // beam parameters at the vertex
  	if (!(PaAlgo::CrossCells(beam_track.vTPar(0), Run, 1.9, 1.2, -318.5, -78.5, 2.0))) {
		return false;
	}

	// Check that the track meantime is within flux requirements
	double mean_time = beam_track.MeanTime();  
    if (std::fabs(mean_time) >= 2) {
		return false;
	}

	// Time in spill is within flux requirements
	if (!TiS_flag) {
		return false; 
	}

	// Return true if all conditions are met
	return true;
} 


// *************************  ...................  *************************** 
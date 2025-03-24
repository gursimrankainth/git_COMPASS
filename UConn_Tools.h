#pragma once
#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <cmath>
#include <tuple>
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
// 5. Notes                                                       //
// ************************************************************** //

// *************************  EVENTFLAGS (EVT. STATISTICS STRUCTURE WITH MEMBER FUNCS)  ***************************
// 1. EventFlags: Structure (with member funcs) that holds all flags and an associated counter for event statistics
// 2. Inputs: N/A, you just need to create an isntance of this outside of your UserEvent loop 
// 3. Output: N/A, you print all cut statistics if wanted 
// 4. Example usage: https://github.com/gursimrankainth/git_COMPASS/blob/main/u970_DVCS.cc  
//    ... 
//    EventFlags eventFlags; // Create an instance of flags and counters 
//    eventFlags.createFlag("flag name", "description"); // create a new flag (needs to be in event loop above resetFlags)
//    eventFlags.resetFlags(); // reset all flags to false at the start of each event loop
//    eventFlags.setFlagByName("allEvts_flag") // set flag to true or false by its name 
//    eventFlags.incrementCounters() // increment all counters with true flags in main event loop (outside of vertex loop)
//    eventFlags.printFlags() // print cut statistics
// 5. Can also create a new flag usng setFlagByName though the reccomendation is to use createFlag and setFlagByName


// Define a structure and its member functions 
struct EventFlags {
    // Vector of tuples: (flag, counter, flag name, description)
    std::vector<std::tuple<bool, int, std::string, std::string>> flags;

    // Default constructor to initialize flags
    EventFlags() {
		// Event flags 
        flags.push_back({false, 0, "allEvts_flag", "Total no. of events processed by PHAST user script"});
        flags.push_back({false, 0, "pVtx_flag", "No. of events with a primary vertex"});
		// inMu flags
        flags.push_back({false, 0, "inMuTrack_flag", "No. of events where beam has a track with parameters"});
        flags.push_back({false, 0, "zFirst_flag", "No. of events where beam was first measured before the target"});
        flags.push_back({false, 0, "momRange_flag", "No. of events where beam momentum falls within flux requirements"});
        flags.push_back({false, 0, "momErr_flag", "No. of events where beam momentum error falls within flux requirements"});
        flags.push_back({false, 0, "BMS_flag", "No. of events where beam is detected by BMS"});
        flags.push_back({false, 0, "FI_flag", "No. of events where beam is detected by SCIFI"});
        flags.push_back({false, 0, "SI_flag", "No. of events where beam is detected by SI"});
        flags.push_back({false, 0, "crossCells_flag", "No. of events where beam crosses full target length"});
        flags.push_back({false, 0, "meantime_flag", "No. of events where beam track meantime is within flux requirements"});
        flags.push_back({false, 0, "timeInSpill_flag", "No. of events where time in spill is within flux requirements"});
		// outMu flags 
        flags.push_back({false, 0, "vtxInTarget_flag", "No. of events where the vertex is in the target"});
        flags.push_back({false, 0, "trigger_flag", "No. of events with MT, LT, OT or LAST physics triggers"});
        flags.push_back({false, 0, "passHodo_flag", "No. of events where scattered muon passes Hodoscope check"});
        flags.push_back({false, 0, "charge_flag", "No. of events where scattered muon has the same charge as the beam"});
        flags.push_back({false, 0, "zFirstLast_flag", "No. of events where first and last scattered muon z coord. are measured before and after SM1"});
    }

    // Function to create a new flag dynamically
    void createFlag(const std::string& flagName, const std::string& description) {
        // Check if flag already exists
        for (const auto& flag : flags) {
            if (std::get<2>(flag) == flagName) return; // Don't add duplicate flags
        }
        flags.push_back({false, 0, flagName, description});
    }

    // Function to reset all flags
    void resetFlags() {
        for (auto& flag : flags) {
            std::get<0>(flag) = false;
        }
    }

    // Function to set a flag by name (creates if not found)
    void setFlagByName(const std::string& flagName, bool value, const std::string& description = "") {
        for (auto& flag : flags) {
            if (std::get<2>(flag) == flagName) {
                std::get<0>(flag) = value;
                return;
            }
        }
        // If not found, create a new flag dynamically
        if (!description.empty()) {
            createFlag(flagName, description);
            std::get<0>(flags.back()) = value; // Set newly added flag
        }
    }

    // Function to increment counters of active flags
    void incrementCounters() {
        for (auto& flag : flags) {
            if (std::get<0>(flag)) {
                std::get<1>(flag)++;
            }
        }
    }

    // Function to print flags
    void printFlags() const {
        for (const auto& flag : flags) {
            std::cout << std::get<1>(flag) << " : " << std::get<3>(flag) << std::endl;
        }
    }
};


// *************************  PRINTDEBUG  *************************** 
// 1. Function: printDebug -> print debug statements if verbose_mode is set to true 
// 2. Input: N/A, see below how to use  
// 3. Output: prints statements to console while PHAST us running for involved debugging 
// 4. Example usage: https://github.com/gursimrankainth/git_COMPASS/blob/main/u970_DVCS.cc  
//   ... 
//   verbose_mode = true; // set this manually in your own script for ease of use 
//   printDebug("*** Run: " + std::to_string(Run) + ", spill: " + std::to_string(Spill) + ", event: " + std::to_string(EvtInSpill) + " ***");
//   OR
//   printDebug("test", true); // if second arguement is set to true it will print even if verbose_mode is off 

// Global flag for verbose mode 
extern bool verbose_mode; // Set to true for verbose output, false to suppress

// Define the function (to print statement even if verbose mode is off set second arguement to true)
void printDebug(const std::string &message, bool forcePrint = false) {
    if (verbose_mode || forcePrint) {
        std::cout << message << std::endl;
    }
}


// *************************  BEAMFLUXCHECK  *************************** 
// 1. Function: beamFluxCheck -> check that all of the flux requirements are satisfied by the incoming muon beam
// 2. Input (8): PaEvent object, PaVertex object, int vertex index, Run, bool TiS_flag, BeamFluxParams object,
//               PaParticle (beam), PaTrack (beam_track), PaTPar (Par_beam), EventFlags object 
// 3. Output (1): boolean value that is true for events that pass the check and false for events that do not 
// 4. Example usage: https://github.com/gursimrankainth/git_COMPASS/blob/main/u970_DVCS.cc 

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
bool beamFluxCheck(const PaEvent &e, const PaVertex &v, int vertexIndex, int Run, bool TiS_flag, 
				const BeamFluxParams &params, PaParticle &beam, PaTrack &beam_track, PaTPar &Par_beam, 
				EventFlags &flags) { // beamFlux loop begins 
	
	// Check that there is an incoming muon associated with the vertex
	int i_beam = v.InParticle(); 
	if (i_beam == -1) { 
		return false;
	} 

	// Check that the beam has a track associated with it
	beam = e.vParticle(i_beam);
	int it_beam = beam.iTrack();
	if (it_beam == -1) { 
		return false;
	} 

	// Check that the track has parameters
	beam_track = e.vTrack(i_beam);
	if (beam_track.NTPar() == 0) { 
		return false;
	} else {flags.setFlagByName("inMuTrack_flag", true);}

	// Check that the beam was first measured before the target
	if (beam_track.ZFirst() >= -78.5) {
		return false;
	} else {flags.setFlagByName("zFirst_flag", true);}

	// Check that the beam momentum falls within acceptable range
	double inMu_mom = beam_track.vTPar(0).Mom();
	if (inMu_mom < 140.0 || inMu_mom > 180.0) {
		return false;
	} else {flags.setFlagByName("momRange_flag", true);}

	// Check that the beam momentum error falls within acceptable range
	double inMu_momErr = sqrt(beam_track.vTPar(0)(5,5))/(beam_track.vTPar(0)(5)*beam_track.vTPar(0)(5));
  	if (inMu_momErr > 0.025*inMu_mom) {
		return false;
	} else {flags.setFlagByName("momErr_flag", true);}

	// Check that the beam is detected by detectors along the beamline
	int nhits_BMS = beam_track.NHitsFoundInDetect("BM");
	int nhits_FI  = beam_track.NHitsFoundInDetect("FI"); 
	int nhits_SI  = beam_track.NHitsFoundInDetect("SI"); 
	if (nhits_BMS < 3) {
    	return false;
	} else {flags.setFlagByName("BMS_flag", true);}
	if (nhits_FI < 2) {
    	return false;
	} else {flags.setFlagByName("FI_flag", true);}
	if (nhits_SI < 3) {
    	return false;
	} else {flags.setFlagByName("SI_flag", true);}

	// Check that the beam crosses the full target length 
	// PaAlgo::CrossCells(t_beam.vTPar(0),run, Rmax, Ymax, tgt_zmin, tgt_zmax, RmaxMC) 
	Par_beam = beam.ParInVtx(vertexIndex); // beam parameters at the vertex
  	if (!(PaAlgo::CrossCells(beam_track.vTPar(0), Run, params.Rmax, params.Ymax, params.tgt_zmin, params.tgt_zmax, params.RmaxMC))) {
		return false;
	} else {flags.setFlagByName("crossCells_flag", true);}

	// Check that the track meantime is within flux requirements
	double mean_time = beam_track.MeanTime();  
    if (std::fabs(mean_time) >= 2) {
		return false;
	} else {flags.setFlagByName("meantime_flag", true);}

	// Time in spill is within flux requirements
	if (!TiS_flag) {
		return false; 
	} else {flags.setFlagByName("timeInSpill_flag", true);}

	// Return true if all conditions are met
	return true;
} 


// *************************  OUTMUCHECK  *************************** 
// 1. Function: outMuCheck -> check that all requirements are satisfied by the scattered muon
// 2. Input (11): PaEvent object, PaVertex object, int vertex index, int Run, PaParticle object (beam), 
//                PaHodoHelper object, bool trig_flag, OutMuParams object, PaParticle object (outMu), 
//                PaTrack object (outMu_track), PaTPar object (Par_outMu), EventFlags object
// 3. Output (1): boolean value that is true for events that pass the check and false for events that do not 
// 4. Example usage: https://github.com/gursimrankainth/git_COMPASS/blob/main/u970_DVCS.cc 

// Define a struct to hold all the parameters
struct OutMuParams {
    int Run;
    double Rmax;
    double Ymax;
    double tgt_zmin;
    double tgt_zmax;
    int RmaxMC;
    double zfirstlast;   

    // Constructor with default values
    OutMuParams(double rmax = 1.9, double ymax = 1.2, double zmin = -318.5, 
                   double zmax = -78.5, int rmaxMC = 2, double zfirstlast = 350) 
        : Rmax(rmax), Ymax(ymax), tgt_zmin(zmin), tgt_zmax(zmax), RmaxMC(rmaxMC),
          zfirstlast(zfirstlast) {}
};

// Define the function 
bool outMuCheck(const PaEvent &e, const PaVertex &v, int vertexIndex, int Run, const PaParticle &beam, 
							PaHodoHelper* HodoHelper, bool trig_flag, const OutMuParams &params,
							PaParticle &outMu, PaTrack &outMu_track, PaTPar &Par_outMu, EventFlags &flags) {

	// Get the index of the scattered muon WITHOUT CHECKING IF IT PASSES the hodoscope check
	// HodoHelper->iMuPrim(v, checkYokeSM2, reject2muEvents, checkCanBeMuon, true, minXX0muPr, true, true) 
	int i_omu = -1; 
	i_omu = HodoHelper->iMuPrim(v,false,false,true,false,15);  
	if (i_omu == -1) {
		return false; 
	}

	// Check that the vertex is in the target (same whether you use beam or scattered muon since they have the same vertex)
	const PaTPar& Par_beam = beam.ParInVtx(vertexIndex); // beam parameters at the vertex
	const PaParticle & outMu_noHodo = e.vParticle(i_omu); 
	const PaTPar& Par_outMu_noHodo = outMu_noHodo.ParInVtx(vertexIndex); // scattered muon parameters at the vertex 
	if(!(PaAlgo::InTarget(Par_beam,'O',Run, params.Rmax, params.Ymax, params.tgt_zmin, params.tgt_zmax, params.RmaxMC))) {
		return false; 
	} else {flags.setFlagByName("vtxInTarget_flag", true);}

	// Check that there is a physics trigger for this event (MT, LT, OT or LAST)
	if (trig_flag == 0) {
		return false; 
	} else {flags.setFlagByName("trigger_flag", true);}

	// Get the index for the scattered muon IF IT PASSES the hodoscope check 
	int i_omu_check_hodo = HodoHelper->iMuPrim(v,false,false,true,true,15,true,true);  
	//int i_omu_check_hodo = HodoHelper->iMuPrim(v,false,true,true,true,15,true,true);
	// if scattered muon passed the hodoscope use the corresponding index, if not proceed with other index 
	if (i_omu_check_hodo == -1) {
		return false; 
	} else {flags.setFlagByName("passHodo_flag", true);}
	i_omu = i_omu_check_hodo;

	// Check that the scattered muon has the same charge as the beam  
	outMu = e.vParticle(i_omu);
	//outMu = e.vParticle(i_omu_check_hodo);
	// if outMu.Q or beam.Q return -777 it means the assocaited track was reconstructed in a field free region (charge is unkown)
	if (outMu.Q() != beam.Q() || outMu.Q() == -777 || beam.Q() == -777) {
		return false; 
	} else {flags.setFlagByName("charge_flag", true);}

	int outMu_itrack = outMu.iTrack();
	outMu_track = e.vTrack(outMu_itrack);
	Par_outMu = outMu.ParInVtx(vertexIndex);
	// Check that the first and last z coordinates are measured before and after SM1
	if (!(outMu_track.ZFirst() < params.zfirstlast && outMu_track.ZLast() > params.zfirstlast)) {
		return false; 
	} else {flags.setFlagByName("zFirstLast_flag", true);}

	// If all checks pass, return true and the relevant objects
	return true;

}
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
#include "UConn_Tools.h"

// ************************************************************** //
// .CC script containing functions used for event selection.      //
// Notes for each functon:  					                  //
// 1. Name and use case                                           //
// 2. Inputs (NEED TO BE DEFINED BY USER)                         //
// 3. Output                                                      //
// 4. Example usage                                               //
// 5. Notes                                                       //
// ************************************************************** //

// *************************  EVENTFLAGS (EVT. STATISTICS STRUCTURE WITH MEMBER FUNCS)  ***************************
// 1. EventFlags: Structure (with member funcs) that holds all flags and associated counters for event statistics
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

// Constructor with mode selector 
EventFlags::EventFlags(Mode mode) {
	switch (mode) {
		case DVCS:  
			// Event flags 
			flags.push_back({false, 0, 0, "allEvts_flag", "Total no. of events processed by PHAST user script"});
			flags.push_back({false, 0, 0, "pVtx_flag", "No. of events with a primary vertex"});
			// inMu flags
			flags.push_back({false, 0, 0, "inMuTrack_flag", "No. of events where beam has a track with parameters"});
			flags.push_back({false, 0, 0, "zFirst_flag", "No. of events where beam was first measured before the target"});
			flags.push_back({false, 0, 0, "momRange_flag", "No. of events where beam momentum falls within flux requirements"});
			flags.push_back({false, 0, 0, "momErr_flag", "No. of events where beam momentum error falls within flux requirements"});
			flags.push_back({false, 0, 0, "BMS_flag", "No. of events where beam is detected by BMS"});
			flags.push_back({false, 0, 0, "FI_flag", "No. of events where beam is detected by SCIFI"});
			flags.push_back({false, 0, 0, "SI_flag", "No. of events where beam is detected by SI"});
			flags.push_back({false, 0, 0, "crossCells_flag", "No. of events where beam crosses full target length"});
			flags.push_back({false, 0, 0, "meantime_flag", "No. of events where beam track meantime is within flux requirements"});
			flags.push_back({false, 0, 0, "timeInSpill_flag", "No. of events where time in spill is within flux requirements"});
			// outMu flags 
			flags.push_back({false, 0, 0, "outMuFound_flag", "No. of events with a scattered muon (before full Hodo check)"});
			flags.push_back({false, 0, 0, "passHodo_flag", "No. of events where scattered muon passes full Hodoscope check"});
			flags.push_back({false, 0, 0, "vtxInTarget_flag", "No. of events where the vertex is in the target"});
			flags.push_back({false, 0, 0, "trigger_flag", "No. of events with MT, LT, OT or LAST physics triggers"});
			flags.push_back({false, 0, 0, "passHodoStat_flag", "No. of events where scattered muon passes Hodoscope check (DVCS stats)"});
			flags.push_back({false, 0, 0, "charge_flag", "No. of events where scattered muon has the same charge as the beam"});
			flags.push_back({false, 0, 0, "zFirstLast_flag", "No. of events where first and last scattered muon z coord. are measured before and after SM1"});
			break; 

		case RHO0:
			flags.push_back({false, 0, 0, "allEvts_flag", "Total no. of events processed by PHAST user script"});
			flags.push_back({false, 0, 0, "pVtx_flag", "No. of events with a primary vertex"});
			flags.push_back({false, 0, 0, "inMuTrack_flag", "No. of events where beam has a track with parameters"});
			flags.push_back({false, 0, 0, "passHodo_flag", "No. of events where scattered muon passes Hodoscope check"});
			flags.push_back({false, 0, 0, "charge_flag", "No. of events where scattered muon has the same charge as the beam"});
			flags.push_back({false, 0, 0, "zFirst_flag", "No. of events where beam was first measured before the target"});
			flags.push_back({false, 0, 0, "zFirstLast_flag", "No. of events where first and last scattered muon z coord. are measured before and after SM1"});
			flags.push_back({false, 0, 0, "BMS_flag", "No. of events where beam is detected by BMS"});
			flags.push_back({false, 0, 0, "mu0LH_flag", "No. of events where beam is probability of back propagation is bigger than 0.01"});
			flags.push_back({false, 0, 0, "muChi2_flag", "Beam and mu chi2 check"});
			flags.push_back({false, 0, 0, "momRange_flag", "No. of events where beam momentum falls within flux requirements"});
			flags.push_back({false, 0, 0, "momErr_flag", "No. of events where beam momentum error falls within flux requirements"});
			flags.push_back({false, 0, 0, "radLen_flag", "No. of events where radiation len of mu > 15"});
			flags.push_back({false, 0, 0, "crossCells_flag", "No. of events where beam crosses full target length"});
			flags.push_back({false, 0, 0, "vtxInTarget_flag", "No. of events where the vertex is in the target"});
			flags.push_back({false, 0, 0, "trigger_flag", "No. of events with MT, LT, OT or LAST physics triggers"});
			flags.push_back({false, 0, 0, "timeInSpill_flag", "No. of events where time in spill is within flux requirements"});
			break; 
	}
}

// Create a new flag 
void EventFlags::createFlag(const std::string& flagName, const std::string& description) {
    for (const auto& flag : flags) {
        if (std::get<3>(flag) == flagName) return; // avoid duplicates
    }
    flags.push_back({false, 0, 0, flagName, description});
}

void EventFlags::resetFlags() {
    for (auto& flag : flags) {
        std::get<0>(flag) = false;
    }
}

// Set the flag to true by name (default is false for all flags)
void EventFlags::setFlagByName(const std::string& flagName, bool value, const std::string& description) {
    for (auto& flag : flags) {
        if (std::get<3>(flag) == flagName) {
            std::get<0>(flag) = value;
            if (value) {
                std::get<2>(flag)++;
            }
            return;
        }
    }
    if (!description.empty()) {
        createFlag(flagName, description);
        std::get<0>(flags.back()) = value;
        if (value) {
            std::get<2>(flags.back())++;
        }
    }
}

// Increments counters for all flags (use once at the end of event loop)
void EventFlags::incrementCounters() {
    for (auto& flag : flags) {
        if (std::get<0>(flag)) {
            std::get<1>(flag)++;
        }
    }
} 

// Prints the values of all flags (place in the UserJobEnd block)
void EventFlags::printFlags() const {
    for (const auto& flag : flags) {
        std::cout << std::get<1>(flag) << " : "      // event-level counter
                  << std::get<2>(flag) << " : "      // multi counter
                  << std::get<4>(flag) << std::endl; // description
    }
}

// Returns the current boolean value of the flag with the given name 
bool EventFlags::getFlag(const std::string& flagName) const {
    for (const auto& flag : flags) {
        if (std::get<3>(flag) == flagName)
            return std::get<0>(flag);  // Return the flag's boolean value
    }
    // If flag not found, optionally print warning or just return false
    // std::cerr << "Warning: Flag '" << flagName << "' not found." << std::endl;
    return false;
}

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
void printDebug(const std::string &message, bool forcePrint) {
    if (verbose_mode || forcePrint) {
        std::cout << message << std::endl;
    }
}


// *************************  BEAMFLUXCHECK DVCS *************************** 
// 1. Function: beamFluxCheck -> check that all of the flux requirements are satisfied by the incoming muon beam
// 2. Input (8): PaEvent object, PaVertex object, int vertex index, Run, bool TiS_flag, BeamFluxParams object,
//               PaParticle (beam), PaTrack (beam_track), PaTPar (Par_beam), EventFlags object 
// 3. Output (1): boolean value that is true for events that pass the check and false for events that do not 
// 4. Example usage: https://github.com/gursimrankainth/git_COMPASS/blob/main/u970_DVCS.cc 

// Constructor with default values
BeamFluxParams::BeamFluxParams(double rmax, double ymax, double zmin, 
                               double zmax, int rmaxMC, double zfirst, 
                               int minFI, int minSI, int minBMS, double minMom,
                               double maxMom, double percent)
    : Rmax(rmax), Ymax(ymax), tgt_zmin(zmin), tgt_zmax(zmax), RmaxMC(rmaxMC),
      zfirst(zfirst), minFI(minFI), minSI(minSI), minBMS(minBMS), minMom(minMom),
      maxMom(maxMom), percent(percent)
{}

// Define the function 
bool beamFluxCheck(const PaEvent &e, const PaVertex &v, int vertexIndex, int Run,
				const BeamFluxParams &params, PaParticle &beam, PaTrack &beam_track, 
				PaTPar &Par_beam, EventFlags &flags) { // beamFlux loop begins 

	// Keep just the signatures in a header file and move the body to a cc file 
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
	beam_track = e.vTrack(it_beam);
	if (beam_track.NTPar() == 0) { 
		return false;
	} else {flags.setFlagByName("inMuTrack_flag", true);}

	// Check that the beam was first measured before the target
	if (beam_track.ZFirst() >= -78.5) {
		return false;
	} else {flags.setFlagByName("zFirst_flag", true);}

	// Check that the beam momentum falls within acceptable range
	//double inMu_mom = beam.ParInVtx(vertexIndex).Mom();
	double inMu_mom = beam_track.vTPar(0).Mom();
	if (inMu_mom < 140.0 || inMu_mom > 180.0) {
		return false;
	} else {flags.setFlagByName("momRange_flag", true);}

	// Check that the beam momentum error falls within acceptable range
	//double inMu_momErr = sqrt(beam.ParInVtx(vertexIndex)(5,5))/(beam.ParInVtx(vertexIndex)(5)*beam.ParInVtx(vertexIndex)(5));
	double inMu_momErr = sqrt(beam_track.vTPar(0)(5,5))/(beam_track.vTPar(0)(5)*beam_track.vTPar(0)(5));
  	if (inMu_momErr > 0.025*inMu_mom) {
		return false;
	} else {flags.setFlagByName("momErr_flag", true);}

	// Check that the beam is detected by detectors along the beamline
	int nhits_BMS = beam_track.NHitsFoundInDetect("BM");
	int nhits_FI  = beam_track.NHitsFoundInDetect("FI"); 
	int nhits_SI  = beam_track.NHitsFoundInDetect("SI"); 
	
	if (!e.IsMC()) { // skip BMS check for MC events 
		if (nhits_BMS < 3) {
			return false; 
		} else {flags.setFlagByName("BMS_flag", true);}
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
  	if (!PaAlgo::CrossCells(beam_track.vTPar(0), Run, params.Rmax, params.Ymax, params.tgt_zmin, params.tgt_zmax, params.RmaxMC)) {
		return false; 
	} else {flags.setFlagByName("crossCells_flag", true);}

	// Check that the track meantime is within flux requirements
	double mean_time = beam_track.MeanTime();  
	if (!e.IsMC()) { // ignore meantime check for MC 
		if (std::fabs(mean_time) >= 2) {
			return false;
		} else {flags.setFlagByName("meantime_flag", true);}
	} else {flags.setFlagByName("meantime_flag", true);}

	// Return true if all conditions are met
	return true;
} 


// *************************  OUTMUCHECK DVCS *************************** 
// 1. Function: outMuCheck -> check that all requirements are satisfied by the scattered muon
// 2. Input (11): PaEvent object, PaVertex object, int vertex index, int Run, PaParticle object (beam), 
//                PaHodoHelper object, bool trig_flag, OutMuParams object, PaParticle object (outMu), 
//                PaTrack object (outMu_track), PaTPar object (Par_outMu), EventFlags object
// 3. Output (1): boolean value that is true for events that pass the check and false for events that do not 
// 4. Example usage: https://github.com/gursimrankainth/git_COMPASS/blob/main/u970_DVCS.cc 

// Constructor definition
OutMuParams::OutMuParams(double rmax, double ymax, double zmin, double zmax,
                         double rmaxMC, double zfirstlast)
    : Rmax(rmax), Ymax(ymax), tgt_zmin(zmin), tgt_zmax(zmax), RmaxMC(rmaxMC),
      zfirstlast(zfirstlast) {}

// Define the function 
bool outMuCheck(const PaEvent &e, const PaVertex &v, int vertexIndex, int Run, const PaParticle &beam, 
				PaHodoHelper* HodoHelper, bool trig_flag, bool TiS_flag, const OutMuParams &params,
				PaParticle &outMu, PaTrack &outMu_track, PaTPar &Par_outMu, EventFlags &flags) {

	// HodoHelper->iMuPrim(v, checkYokeSM2, reject2muEvents, checkCanBeMuon, true, minXX0muPr, true, true) 
	// Loose check: no hodo requirement
	int i_omu = HodoHelper->iMuPrim(v, false, false, true, false, 15);
	if (i_omu == -1) {
		return false; 
	} else {flags.setFlagByName("outMuFound_flag", true);}

	// Check that the vertex is in the target
	const PaTPar& Par_beam = beam.ParInVtx(vertexIndex); // beam parameters at the vertex
	if(!PaAlgo::InTarget(Par_beam,'O', Run, params.Rmax, params.Ymax, params.tgt_zmin, params.tgt_zmax, params.RmaxMC)) {
		return false; 
	} else {
		if (TiS_flag) {flags.setFlagByName("vtxInTarget_flag", true);} 
	}

	// Check that there is a physics trigger for this event (MT, LT, OT or LAST)
	if (trig_flag == 0) {
		return false; 
	} else {
		if (TiS_flag) {flags.setFlagByName("trigger_flag", true);}
	}
 
	// Strict check: with hodo and full cuts
	bool hodo_passed_DVCS = false; 
	int i_omu_check_hodo = HodoHelper->iMuPrim(v, false, false, true, true, 15, true, true);
	if (i_omu_check_hodo != -1) {
		i_omu = i_omu_check_hodo;  // Use better muon if found
		if (TiS_flag) {flags.setFlagByName("passHodoStat_flag", true);}
		hodo_passed_DVCS = true;
		flags.setFlagByName("passHodo_flag", true);
	}

	// Check that the scattered muon has the same charge as the beam  
	bool charge = false; 
	outMu = e.vParticle(i_omu);
	TVector3 SM2center(0, 0, 1825);
    const TVector3 SM2field = PaSetup::Ref().MagField(SM2center);
    const int inMuQ = SM2field(1) < 0 ? 1 : -1;
	if (outMu.Q() != inMuQ) {
		return false; 
	} else {
		if (hodo_passed_DVCS && TiS_flag) {
			flags.setFlagByName("charge_flag", true);
			charge = true; 
		}
	} 

	int outMu_itrack = outMu.iTrack();
	outMu_track = e.vTrack(outMu_itrack);
	Par_outMu = outMu.ParInVtx(vertexIndex);
	// Check that the first and last z coordinates are measured before and after SM1
	if (!(outMu_track.ZFirst() < params.zfirstlast && outMu_track.ZLast() > params.zfirstlast)) {
		return false; 
	} else {
		if (charge && TiS_flag && hodo_passed_DVCS) {flags.setFlagByName("zFirstLast_flag", true);}
	}

	// If all checks pass, return true and the relevant objects
	return true;

}


// *************************  REMOVE EXCLUSIVE LEPTO EVTS (DVCS)  ***************************
// This funtion is an adaptation of LeptoPi0Filter.cc (https://gitlab.cern.ch/compass/dvcs_sidis/2016-dvcs-analysis/-/blob/master/UserEvent/LeptoPi0Filter.cc?ref_type=heads)
// 1. Function: exclLepto() -> check if the event is exclusive (these need to be removed from the sample)
// 2. Input (2): PaEvent object, bool leptoMC 
// 3. Output (1): boolean value that is true for events that pass the check and false for events that do not 
// 4. Example usage: https://github.com/gursimrankainth/git_COMPASS/blob/main/u970_DVCS.cc 

// Global flag for LEPTO MC  
extern bool leptoMC; // Set to true when processing LEPTO data, false otherwise 

bool exclLepto(const PaEvent &e) {
	// Define constants 
	const int kLujetFlavourCodeMuMinus = +13;
	const int kLujetFlavourCodeMuPlus  = -13;
	const int kLujetFlavourCodeProton  = +2212;
	const int kLujetFlavourCodeGamma   = +22;
	const int kLujetStatusCodeUnfragmentedParticle = +1;

	// Particle counters
	unsigned int topology_count = 0U;
	unsigned int muon_count     = 0U;
	unsigned int proton_count   = 0U;
	unsigned int gamma_count    = 0U;
	unsigned int other_count    = 0U;

	int Nparticles;
	std::vector<LUJET> vlj;
	e.MCgen(Nparticles, vlj);

	for (const auto &v : vlj) { // Begin loop over particles in lujet
		if (v.k[0] != kLujetStatusCodeUnfragmentedParticle) continue;
		++topology_count;

		const int flavour_code = v.k[1];
		if (flavour_code == kLujetFlavourCodeMuMinus || flavour_code == kLujetFlavourCodeMuPlus) {
			++muon_count;
		} else if (flavour_code == kLujetFlavourCodeProton) {
			++proton_count;
		} else if (flavour_code == kLujetFlavourCodeGamma) {
			++gamma_count;
		} else {
			++other_count;
		}

	} // end loop over particles in lujet 

	bool is_pi0 = (muon_count > 0U) && (proton_count > 0U) && (gamma_count > 1U) && (topology_count == 4U);

	// Return true if this is NOT a pi0-exclusive event
	return !is_pi0;
}


// *************************  ANGLE BETWEEN THE REAL AND VIRTUAL PHOTON SCATTERING PLANE (DVCS)  ***************************
// This funtion is an adaptation of calc_phi_2gamma.cc (https://gitlab.cern.ch/compass/dvcs_sidis/2016-dvcs-analysis/-/blob/master/UserEvent/calc_phi_2gamma.cc?ref_type=heads)
// 1. Function: phiRV() -> find the azimuthal angle between the lepton scattering plane and the real photon production plane
// 2. Input (3): TLorentz vectors for the inMU, outMu, and real photon  
// 3. Output (1): double value for the angle in radians
// 4. Example usage: https://github.com/gursimrankainth/git_COMPASS/blob/main/u970_DVCS.cc 

double phiRV(TLorentzVector inMu_TL, TLorentzVector outMu_TL, TLorentzVector proton_TL, TLorentzVector gamma_TL, bool eIsMC) {
	TLorentzVector q = inMu_TL - outMu_TL; // four-vector of the virtual photon 
	TVector3 q3 = q.Vect(); // three-vector of the virtual photon 
	TVector3 g3 = gamma_TL.Vect(); // three-vector of the real photon 
	TVector3 p3 = proton_TL.Vect(); // three-vector of the recoil proton

	TVector3 l3; // three-vector of the incoming lepton
	if (eIsMC) {
		l3 = -inMu_TL.Vect(); // flip the sign for MC to align with real data 
	} else {
		l3 = inMu_TL.Vect(); 
	}

	TVector3 norm_lep = q3.Cross(l3).Unit(); // unit normal vector to the lepton plane 
	TVector3 norm_had = g3.Cross(p3).Unit(); // unit normal vector to the real photon plane 

	double sign_lab = norm_lep.Dot(g3); // Get direction of phi
	double phi = norm_lep.Angle(norm_had);
	if (sign_lab < 0) phi = -phi; // correct direction of phi if needed

	return phi; 
}


// *************************  BUILD PHOTON VECTORS (DVCS)  ***************************
// This funtion is used to build TLorentz vectors containing information about the cluster and reconstructed real photon 
// 1. Function: buildClusterVecs() -> build the TLorentz vectors for the clusters 
// 2. Input (2): PaEvent object, PaVertex object, int cluster index, TLorentz vectors for the photon and cluster 
// 3. Output (1): nothing -> it will jsut modify the already existing vectors it takes as an input  
// 4. Example usage: https://github.com/gursimrankainth/git_COMPASS/blob/main/u970_DVCS.cc 
void buildClusterVecs(const PaEvent &e, const PaVertex &v, int cl_id, TLorentzVector &gamma_TL, TLorentzVector &cluster_TL) {
	const auto& cluster = e.vCaloClus(cl_id);
  double phdz = cluster.Z() - v.Z();
  double phdy = cluster.Y() - v.Y();
  double phdx = cluster.X() - v.X();
  double phL  = TMath::Sqrt(phdx * phdx + phdy * phdy + phdz * phdz);
  double cl_E = cluster.E();
  gamma_TL.SetPxPyPzE(cl_E * phdx / phL, cl_E * phdy / phL, cl_E * phdz / phL, cl_E);
  cluster_TL.SetXYZT(cluster.X(), cluster.Y(), cluster.Z(), cluster.E());
}

 // *************************  X-CHECK 1 RHO0 *************************** 
bool crossCheck(const PaEvent &e, const PaVertex & v,int iv, int Run,const BeamFluxParams &params, PaParticle &beam, PaTrack &beam_track, PaTPar &Par_beam, PaHodoHelper* HodoHelper, const OutMuParams &outparams, PaParticle &outMu, PaTrack &outMu_track,  PaTPar &Par_outMu,EventFlags &flags){
	
	flags.setFlagByName("TEST","TEST");
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


	int i_omu = -1;
        i_omu = HodoHelper->iMuPrim(v,false,false,true,false,15);
        if (i_omu == -1) {
                return false;
        }

        const PaParticle & outMu_noHodo = e.vParticle(i_omu);
        const PaTPar& Par_outMu_noHodo = outMu_noHodo.ParInVtx(iv); // scattered muon parameters at the vertex
        /*if(!(PaAlgo::InTarget(Par_beam,'O',Run, params.Rmax, params.Ymax, params.tgt_zmin, params.tgt_zmax, params.RmaxMC))) {
                return false;
        } else {flags.setFlagByName("vtxInTarget_flag", true);}*/

        // Check that there is a physics trigger for this event (MT, LT, OT or LAST)
        /*if (trig_flag == 0) {
                return false;
        } else {flags.setFlagByName("trigger_flag", true);}*/

        // Get the index for the scattered muon IF IT PASSES the hodoscope check
        int i_omu_check_hodo = HodoHelper->iMuPrim(v,false,false,true,true,15,true,true);
        //int i_omu_check_hodo = HodoHelper->iMuPrim(v,false,true,true,true,15,true,true);
        // if scattered muon passed the hodoscope use the corresponding index, if not proceed with other index
        if (i_omu_check_hodo == -1) {
                return false;
        } else {flags.setFlagByName("passHodo_flag", true);}
        i_omu = i_omu_check_hodo;
	
	if (beam_track.ZFirst() > -78.5) {
                return false;
        } else {flags.setFlagByName("zFirst_flag", true);}
	// Check that the scattered muon has the same charge as the beam
        outMu = e.vParticle(i_omu);
        // if outMu.Q or beam.Q return -777 it means the assocaited track was reconstructed in a field free region (charge is unkown)
        if (outMu.Q() != beam.Q() || outMu.Q() == -777 || beam.Q() == -777) {
                return false;
        } else {flags.setFlagByName("charge_flag", true);}
	
	int outMu_itrack = outMu.iTrack();
        outMu_track = e.vTrack(outMu_itrack);
        Par_outMu = outMu.ParInVtx(iv);
        // Check that the first and last z coordinates are measured before and after SM1
        if (!(outMu_track.ZFirst() < outparams.zfirstlast && outMu_track.ZLast() > outparams.zfirstlast)) {
                return false;
        } else {flags.setFlagByName("zFirstLast_flag", true);}

	  double mean_time = beam_track.MeanTime();
    	if (std::fabs(mean_time) >= 2) {
                return false;
        } else {flags.setFlagByName("meantime_flag", true);}

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
        Par_beam = beam.ParInVtx(iv); // beam parameters at the vertex
        if (!(PaAlgo::CrossCells(beam_track.vTPar(0), Run, params.Rmax, params.Ymax, params.tgt_zmin, params.tgt_zmax, params.RmaxMC))) {
                return false;
        } else {flags.setFlagByName("crossCells_flag", true);}

        if(!(PaAlgo::InTarget(Par_beam,'O',Run, params.Rmax, params.Ymax, params.tgt_zmin, params.tgt_zmax, params.RmaxMC))) {
                return false;
        } else {flags.setFlagByName("vtxInTarget_flag", true);}	


	// Check that the beam momentum falls within acceptable range
        double inMu_mom = beam.ParInVtx(iv).Mom();
        if (inMu_mom < 140.0 || inMu_mom > 180.0) {
                return false;
        } else {flags.setFlagByName("momRange_flag", true);}

        // Check that the beam momentum error falls within acceptable range
        double inMu_momErr = sqrt(beam.ParInVtx(iv)(5,5))/(beam.ParInVtx(iv)(5)*beam.ParInVtx(iv)(5));
        if (inMu_momErr > 0.025*inMu_mom) {
                return false;
        } else {flags.setFlagByName("momErr_flag", true);}

	return true;

}


// *************************  X-CHECK 2 RHO0 *************************** 
bool crossCheck2(const PaEvent &e, const PaVertex & v,int iv, int Run,const BeamFluxParams &params, PaParticle &beam, PaTrack &beam_track, PaTPar &Par_beam, PaHodoHelper* HodoHelper, const OutMuParams &outparams, PaParticle &outMu, PaTrack &outMu_track,  PaTPar &Par_outMu,EventFlags &flags,bool isMC){
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


        int i_omu = -1;
        i_omu = HodoHelper->iMuPrim(v,false,false,true,false,15);
        if (i_omu == -1) {
                return false;
        }

        const PaParticle & outMu_noHodo = e.vParticle(i_omu);
        const PaTPar& Par_outMu_noHodo = outMu_noHodo.ParInVtx(iv); // scattered muon parameters at the vertex
        /*if(!(PaAlgo::InTarget(Par_beam,'O',Run, params.Rmax, params.Ymax, params.tgt_zmin, params.tgt_zmax, params.RmaxMC))) {
                return false;
        } else {flags.setFlagByName("vtxInTarget_flag", true);}*/

        // Check that there is a physics trigger for this event (MT, LT, OT or LAST)
        /*if (trig_flag == 0) {
                return false;
        } else {flags.setFlagByName("trigger_flag", true);}*/

        // Get the index for the scattered muon IF IT PASSES the hodoscope check
        int i_omu_check_hodo = HodoHelper->iMuPrim(v,false,false,true,true,15,true,true);
	if (i_omu_check_hodo == -1) {
		 return false;
	} else {flags.setFlagByName("passHodo_flag", true);}
	i_omu = i_omu_check_hodo;
	
	outMu = e.vParticle(i_omu);
	// if outMu.Q or beam.Q return -777 it means the assocaited track was reconstructed in a field free region (charge is unkown)
	if (outMu.Q() != beam.Q() || outMu.Q() == -777 || beam.Q() == -777) {
		return false;
	} else {flags.setFlagByName("charge_flag", true);}
	int outMu_itrack = outMu.iTrack();
	outMu_track = e.vTrack(outMu_itrack);
	Par_outMu = outMu.ParInVtx(iv);
	
	double Zprim = v.Pos(2);
	if (Zprim<-318.5 || Zprim>-78.5){
		return false;

	} else { flags.setFlagByName("zFirst_flag", true);}
	

        // Check that the first and last z coordinates are measured before and after SM1
        if (!(outMu_track.ZFirst() < outparams.zfirstlast && outMu_track.ZLast() > outparams.zfirstlast)) {
                return false;
        } else {flags.setFlagByName("zFirstLast_flag", true);}	
	
	int nhits_BMS = beam_track.NHitsFoundInDetect("BM");
	if (nhits_BMS < 3 and !(isMC)) {
		   return false;
	} else {flags.setFlagByName("BMS_flag", true);}
	
	 if (beam.BackPropLH() < 0.01) {
                return false;
        }else {flags.setFlagByName("mu0LH_flag", true);}



	if (beam_track.Chi2tot()/beam_track.Ndf() > 10){
		 return false;
	}

	if (outMu_track.Chi2tot()/outMu_track.Ndf() > 10){
                 return false;
        } else 
	flags.setFlagByName("muChi2_flag", true);
	//beam momentum falls within acceptable range
        /*double inMu_mom = beam.ParInVtx(iv).Mom();
        if (inMu_mom < 140.0 || inMu_mom > 180.0) {
                return false;
        } else {flags.setFlagByName("momRange_flag", true);}

        // Check that the beam momentum error falls within acceptable range
        double inMu_momErr = sqrt(beam.ParInVtx(iv)(5,5))/(beam.ParInVtx(iv)(5)*beam.ParInVtx(iv)(5));
        if (inMu_momErr > 0.025*inMu_mom) {
                return false; 
        } else {flags.setFlagByName("momErr_flag", true);}
	*/

	if(outMu_track.XX0() < 15){
		return false;
	} else {flags.setFlagByName("radLen_flag", true);}

	// Check that the beam crosses the full target length 
	// PaAlgo::CrossCells(t_beam.vTPar(0),run, Rmax, Ymax, tgt_zmin, tgt_zmax, RmaxMC) 
	Par_beam = beam.ParInVtx(iv); // beam parameters at the vertex
  	/*if (!(PaAlgo::CrossCells(beam_track.vTPar(0), Run, params.Rmax, params.Ymax, params.tgt_zmin, params.tgt_zmax, params.RmaxMC))) {
		return false;
	} else {flags.setFlagByName("crossCells_flag", true);}
	*/
	/*if(!(PaAlgo::InTarget(Par_beam,'O',Run, params.Rmax, params.Ymax, params.tgt_zmin, params.tgt_zmax, params.RmaxMC))) {
		return false; 
	} else {flags.setFlagByName("vtxInTarget_flag", true);}*/

/*	PaVertex tempVertex = v;
	PaTrack tempTrack = e.vTrack(i_beam);
	if(!(target.InTarget(&tempVertex, 1.9))) {
		 return false;
	}
	if(!target.InTarget(&tempTrack,-318.5,-78.5,1.9)){
		return false;
        }
	else{flags.setFlagByName("vtxInTarget_flag", true);}	

*/


	return true;
} 

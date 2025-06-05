#ifndef UCONN_TOOLS_H
#define UCONN_TOOLS_H

#pragma once
#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <tuple>
#include "Phast.h"
#include "PaSetup.h"
#include "PaAlgo.h"
#include "PaEvent.h"
#include "TLorentzVector.h"
#include "PaHodoHelper.h"

// ************************************************************** //
// .H script containing functions used for event selection.       //
// See .CC script for function def. and additional details.       //
// ************************************************************** //

// Flags and an associated counters for event statistics
struct EventFlags { 
    std::vector<std::tuple<bool, int, int, std::string, std::string>> flags;

    enum Mode { DVCS, RHO0 };

    EventFlags(Mode mode = DVCS);

    void createFlag(const std::string& flagName, const std::string& description);
    void resetFlags();
    void setFlagByName(const std::string& flagName, bool value, const std::string& description = "");
    void incrementCounters();
    void printFlags() const;
};

// Print debug statements if verbose_mode is true 
extern bool verbose_mode;
void printDebug(const std::string &message, bool forcePrint = false); 


// *************************  DVCS FUNCTIONS FOR PHAST EVENT SELECTION  ***************************
//Structure for storing default values for beamFluxCheck()
struct BeamFluxParams {
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
    double percent;

    // Constructor declaration
    BeamFluxParams(double rmax = 1.9, double ymax = 1.2, double zmin = -318.5, 
                   double zmax = -78.5, int rmaxMC = 2, double zfirst = -78.5, 
                   int minFI = 2, int minSI = 3, int minBMS = 3, double minMom = 14.0,
                   double maxMom = 180.0, double percent = 0.025);
};


bool beamFluxCheck(const PaEvent &e, const PaVertex &v, int vertexIndex, int Run, bool TiS_flag, 
                const BeamFluxParams &params, PaParticle &beam, PaTrack &beam_track, 
                PaTPar &Par_beam, EventFlags &flags);

//Structure for storing default values for outMuCheck()
struct OutMuParams {
    int Run;
    double Rmax;
    double Ymax;
    double tgt_zmin;
    double tgt_zmax;
    int RmaxMC;
    double zfirstlast;

    // Constructor declaration only
    OutMuParams(double rmax = 1.9, double ymax = 1.2, double zmin = -318.5, 
                double zmax = -78.5, int rmaxMC = 2, double zfirstlast = 350);
};

bool outMuCheck(const PaEvent &e, const PaVertex &v, int vertexIndex, int Run, const PaParticle &beam, 
				PaHodoHelper* HodoHelper, bool trig_flag, const OutMuParams &params,
				PaParticle &outMu, PaTrack &outMu_track, PaTPar &Par_outMu, EventFlags &flags); 

bool exclLepto (const PaEvent &e, bool leptoMC); 

double phiRV(TLorentzVector inMu_TL, TLorentzVector outMu_TL, TLorentzVector proton_TL, TLorentzVector gamma_TL, bool eIsMC = false);

// *************************  RHO0 FUNCTIONS FOR PHAST EVENT SELECTION  ***************************
bool crossCheck(const PaEvent &e, const PaVertex & v,int iv, int Run,const BeamFluxParams &params, 
                PaParticle &beam, PaTrack &beam_track, PaTPar &Par_beam, PaHodoHelper* HodoHelper, const OutMuParams &outparams, 
                PaParticle &outMu, PaTrack &outMu_track,  PaTPar &Par_outMu,EventFlags &flags); 

bool crossCheck2(const PaEvent &e, const PaVertex & v,int iv, int Run,const BeamFluxParams &params, 
                PaParticle &beam, PaTrack &beam_track, PaTPar &Par_beam, PaHodoHelper* HodoHelper, const OutMuParams &outparams, 
                PaParticle &outMu, PaTrack &outMu_track,  PaTPar &Par_outMu,EventFlags &flags,bool isMC = false); 

#endif // UCONN_TOOLS_H
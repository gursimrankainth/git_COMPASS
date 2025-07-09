#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <cmath>

#include "RooRealVar.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "RooPlot.h"
#include "RooAddPdf.h"
#include "RooFitResult.h"
#include "RooCategory.h"
#include "RooSimultaneous.h"

#include "TLorentzVector.h"

using namespace RooFit;

bool isZeroVector(const TLorentzVector& vec) {
    return vec.E() == 0 && vec.Px() == 0 && vec.Py() == 0 && vec.Pz() == 0;
}


int make_plots() {

    // Open the ROOT file
    TFile *file = TFile::Open("/Users/gursimran/Desktop/merged_P09_626.root");
    if (!file || file->IsZombie()) {
        std::cerr << "Failed to open file!" << std::endl;
        return 1;
    }

    // Get the TTree from the file
    TTree *tree = (TTree*)file->Get("USR970"); // Replace "tree_name" with the actual tree name
    if (!tree) {
        std::cerr << "Failed to get TTree!" << std::endl;
        return 1;
    }

    // Set Global style options 
    //gStyle->SetOptStat(1110);              // Show entries, mean, RMS (adjust flags if needed)
    gStyle->SetStatColor(0);               // Transparent background
    gStyle->SetStatStyle(0);
    gStyle->SetStatTextColor(kBlack);      // Text color
    gStyle->SetStatX(0.9);                 // X position (NDC)
    gStyle->SetStatY(0.9);                 // Y position (NDC)
    gStyle->SetStatFontSize(0.03);         // Font size
    gStyle->SetLabelSize(0.04, "XYZ");     // Axis label size

    gStyle->SetPadLeftMargin(0.15);
    gStyle->SetPadBottomMargin(0.15);
    gStyle->SetTitleOffset(1.3, "X");

    // Variable to hold the value you want to read
    TLorentzVector *inMu_TL = nullptr, *outMu_TL = nullptr, *gamma_TL = nullptr, *p_camera_TL = nullptr;
    TLorentzVector *gammaLow_TL = nullptr, targ_TL(0, 0, 0, 0.93827231);
    tree->SetBranchAddress("inMu_TL", &inMu_TL);
    tree->SetBranchAddress("outMu_TL", &outMu_TL);
    tree->SetBranchAddress("gamma_TL", &gamma_TL);
    tree->SetBranchAddress("p_camera_TL", &p_camera_TL);
    tree->SetBranchAddress("gammaLow_TL", &gammaLow_TL);

    double y, t, nu, Q2, xbj, W2;
    tree->SetBranchAddress("y", &y);
    tree->SetBranchAddress("t", &t);
    tree->SetBranchAddress("nu", &nu);
    tree->SetBranchAddress("Q2", &Q2);
    tree->SetBranchAddress("W2", &W2);
    tree->SetBranchAddress("xbj", &xbj);

    int low_calo, excl_calo; 
    double M_pi0, E_gammaLow; 
    tree->SetBranchAddress("M_pi0", &M_pi0);
    tree->SetBranchAddress("low_calo", &low_calo);
    tree->SetBranchAddress("excl_calo", &excl_calo);
    tree->SetBranchAddress("E_gammaLow", &E_gammaLow); 

    double delta_phi, delta_pt, delta_Z, M2x;
    tree->SetBranchAddress("delta_phi", &delta_phi);
    tree->SetBranchAddress("delta_pt", &delta_pt);
    tree->SetBranchAddress("delta_Z", &delta_Z);
    tree->SetBranchAddress("M2x", &M2x);

    // Create histogramS
    const int nBins = 100;      // Number of bins
    double xMin = 1e-3;         // Minimum x value (avoid 0 because log(0) is undefined)
    double xMax = 1.0;          // Maximum x value
    double binEdges[nBins + 1]; // Bin edges array
    for (int i = 0; i <= nBins; ++i) {
        binEdges[i] = xMin * pow(xMax / xMin, double(i) / nBins);
    }
    TH1F *h97_xbj   = new TH1F("h97_xbj", "Elasticity of the Scattering Process (x_{bj}); x_{bj}; Events", nBins, binEdges);
    TH1F *h97_y     = new TH1F("h97_y", "Fractional Energy Loss of Incoming Muon (y); y; Events", 100, 0, 1);
    TH1F *h97_nu    = new TH1F("h97_nu", "Energy of the virtual photon (#nu); #nu [GeV]; Events", 100, 0, 180);
    TH1F *h97_Q2    = new TH1F("h97_Q2", "Four-momentum Transfer Squared (Lepton) (Q^{2}); Q^{2} [GeV^{2}]; Events", 100, 0, 10);
    TH1F *h97_W2    = new TH1F("h97_W2", "Effective Mass of final state hadrons Squared (W2); W^{2} [GeV^{2}]; Events", 100, 0, 325);
    //TH2F *h97_Q2xbj = new TH2F("h97_Q2xbj", "Kinematic Coverage of Dataset; x_{bj}; Q^{2} [GeV^{2}]", 150, 0, 1, 150, 0, 100);
    TH1F *h97_t     = new TH1F("h97_t", "Four-momentum Transfer Squared (Nucleon) (t); t [GeV^{2}]; Events", 100, -4, 1);

    TH1F *h97_delta_phi = new TH1F("h97_delta_phi", "#Delta#varphi = #varphi^{cam} - #varphi^{spec}; #Delta#phi [rad]; Counts", 100, -0.5, 0.5);
    TH1F *h97_delta_pt  = new TH1F("h97_delta_pt", "#DeltaP_{t} = P_{t}^{cam} - P_{t}^{spec}; #DeltaP_{t} [GeV/c]; Counts", 100, -0.4, 0.4);
    TH1F *h97_delta_Z   = new TH1F("h97_delta_Z", "#DeltaZ_{A} = Z_{A}^{cam} - Z_{A}^{inter}; #DeltaZ_{A} [cm]; Counts", 100, -20, 20);
    TH1F *h97_M2x       = new TH1F("h97_M2x", "M^{2}_{undet} = (k + p - k'- q'- p')^{2}; M^{2}_{undet} [rad]; Counts", 100, -0.5, 0.5);

    TH1F *h97_pi0_M_Comb = new TH1F("h97_pi0_M_comb", "Invariant Mass #pi^{0}; M [GeV/c^{2}]; Counts", 50, 0, 1.6);
    TH1F *h97_pi0_M_EC0 = new TH1F("h97_pi0_M_EC0", "Invariant Mass #pi^{0}; M [GeV/c^{2}]; Counts", 50, 0, 0.3);
    TH1F *h97_pi0_M_EC1 = new TH1F("h97_pi0_M_EC1", "Invariant Mass #pi^{0}; M [GeV/c^{2}]; Counts", 50, 0, 0.3);
    TH1F *h97_pi0_M_EC0_EC1 = new TH1F("h97_pi0_M_EC0_EC1", "Invariant Mass #pi^{0}; M [GeV/c^{2}]; Counts", 50, 0, 0.3);
    //TH1F *h97_E_gammaLow = new TH1F("h97_E_gammaLow", "Low-energy Photon Energy - EC0 & EC1; E [GeV]; Counts", 100, 0, 5);
    //TH1F *h97_M_miss_gamma = new TH1F("h97_M_miss_gamma", "Missing Mass: #mup#rightarrow#mupX; M_{X} [GeV/c^{2}]; Counts", 100, 0, 17);
    TH1F *h97_M2Miss_pR1 = new TH1F("h97_M2Miss_pR1", "Missing Mass Squared (R1): #mup#rightarrow#mupX; M^{2}_{X} [GeV^{2}/c^{4}]; Counts", 50, -30, 30);
    TH1F *h97_M2Miss_pR2 = new TH1F("h97_M2Miss_pR2", "Missing Mass Squared (R2): #mup#rightarrow#mupX; M^{2}_{X} [GeV^{2}/c^{4}]; Counts", 50, -30, 40);
    TH1F *h97_M2Miss_pi0 = new TH1F("h97_M2Miss_pi0", "Missing Mass Squared (R1): #mup#rightarrow#mu#pi^{0}X; M^{2}_{X} [GeV^{2}/c^{4}]; Counts", 50, -20, 80);

    // Loop through entries and fill histogram
    Long64_t nEntries = tree->GetEntries();
    for (Long64_t i = 0; i < nEntries; ++i) {
        tree->GetEntry(i);

        if (isZeroVector(*inMu_TL) || isZeroVector(*outMu_TL) || isZeroVector(*p_camera_TL)) {
            continue; // Skip event if any 4-vector is zero
        }

        bool excl_kin = false; 
        if ((Q2 < 10 && Q2 > 1) && (y > 0.05 && y < 0.9) && (nu > 10 && nu < 144) && (t > -0.64 && t < -0.08)) {
            excl_kin = true;
        }

        h97_xbj->Fill(xbj);
        h97_y->Fill(y);
        h97_Q2->Fill(Q2);
        h97_nu->Fill(nu);
        h97_t->Fill(t);
        h97_W2->Fill(W2);

        h97_delta_phi->Fill(delta_phi); 
        h97_delta_pt->Fill(delta_pt);
        h97_delta_Z->Fill(delta_Z);
        h97_M2x->Fill(M2x);

        if ((low_calo != -999) && (low_calo != 2) && (excl_calo != 2) && excl_kin) {
            // Fill the pi0 invariant mass distributions 
            h97_pi0_M_Comb->Fill(M_pi0); 
            if ((low_calo == 0) && (excl_calo == 0)) {
                h97_pi0_M_EC0->Fill(M_pi0); 
            }
            else if ((low_calo == 1) && (excl_calo == 1)) {
                h97_pi0_M_EC1->Fill(M_pi0); 
            }
            else {
                h97_pi0_M_EC0_EC1->Fill(M_pi0); 
            }

            // Missing Mass distributions 
            TLorentzVector pi0_TL = *gammaLow_TL + *gamma_TL; 
            if (M_pi0 < 0.2) { // region 1 (R1)
                // mup -> mupX (assume the X is missing -> measure the proton) 
                TLorentzVector pMiss_pR1 = *inMu_TL + targ_TL - *outMu_TL - *p_camera_TL; 
                double M2Miss_pR1 = pMiss_pR1.M2(); 
                h97_M2Miss_pR1->Fill(M2Miss_pR1); 
                // mup -> mupi0X (assume the proton is missing -> measure the pi0)
                TLorentzVector pMiss_pi0 = *inMu_TL + targ_TL - *outMu_TL - pi0_TL;
                double M2Miss_pi0 = pMiss_pi0.M2();
                h97_M2Miss_pi0->Fill(M2Miss_pi0);
            }
            else { // region 2 (R2)
                // mup -> mupX (assume the X is missing -> measure the proton) 
                TLorentzVector pMiss_pR2 = *inMu_TL + targ_TL - *outMu_TL - *p_camera_TL; 
                double M2Miss_pR2 = pMiss_pR2.M2(); 
                h97_M2Miss_pR2->Fill(M2Miss_pR2);
            }

        }


/*         if (!inMu_TL || !outMu_TL || !gamma_TL || !p_camera_TL) continue;
        TLorentzVector p_miss_p = *inMu_TL + targ_TL - *outMu_TL - *gamma_TL;
        double M_miss_p = p_miss_p.M();
        TLorentzVector p_miss_gamma = *inMu_TL + targ_TL - *outMu_TL - *p_camera_TL;
        double M_miss_gamma = p_miss_gamma.M(); 

        if ((M_pi0 < 0.0 || M_pi0 > 0.2) && excl_kin) {
            h97_M_miss_p->Fill(M_miss_p);
            h97_M_miss_gamma->Fill(M_miss_gamma);
        } */

    }

/*     TFile myfile("pi0_DVCS.root","RECREATE");
    h97_pi0_M_EC0->Write();
    myfile.Close();  */

    // First canvas: first 5 histograms
    TCanvas *c1 = new TCanvas("c1", "Kinematic Distributions Page 1", 1200, 800);
    c1->Divide(3, 2);
    c1->cd(1);
    gPad->SetLogx();
    h97_xbj->Draw();
    c1->cd(2); h97_y->Draw();
    c1->cd(3); h97_Q2->Draw();
    c1->cd(4); h97_nu->Draw();
    c1->cd(5); h97_t->Draw();
    c1->cd(6); h97_W2->Draw();
    c1->SaveAs("plots.pdf(");  

    TCanvas *c2 = new TCanvas("c2", "Exclusivity Variables Page 2", 1200, 800);
    c2->Divide(3, 2);
    c2->cd(1); h97_delta_phi->Draw();
    c2->cd(2); h97_delta_pt->Draw();
    c2->cd(3); h97_delta_Z->Draw();
    c2->cd(4); h97_M2x->Draw();
    c2->SaveAs("plots.pdf");  

    TCanvas *c3 = new TCanvas("c3", "Pi0 Invariant Mass Page 3", 1200, 800);
    c3->Divide(3, 2);

    gStyle->SetOptFit(1111);
    c3->cd(1); h97_pi0_M_Comb->Draw();
    /*     TF1 *gausFitComb = new TF1("gausFit", "gaus", 0.11, 0.16);
    gausFitComb->SetParameter(1, 0.135); // mean 
    gausFitComb->SetParameter(2, 0.01); // sigma 
    h97_pi0_M_Comb->Fit(gausFitComb, "R"); 
    gausFitComb->SetRange(0.0, 0.3); */
/*     TF1 *fitFunc = new TF1("fitFunc", "[0]*exp(-0.5*((x-[1])/[2])^2) + [3] + [4]*x", 0.1, 0.17);
    fitFunc->SetParameters(140, 0.135, 0.01, 10, 0.15);
    c3->cd(1); h97_pi0_M_Comb->Draw();
    h97_pi0_M_Comb->Fit(fitFunc, "R"); */
    //fitFunc->Draw("same");

/*     c3->cd(1); 
    h97_pi0_M_EC0->SetFillColorAlpha(kRed, 0.3);   // Red, 30% opacity
    h97_pi0_M_EC0->SetLineColor(kRed+1);
    h97_pi0_M_EC1->SetFillColorAlpha(kAzure-5, 0.3);  // Blue, 30% opacity
    h97_pi0_M_EC1->SetLineColor(kAzure-5);
    h97_pi0_M_EC0->Draw("HIST");
    h97_pi0_M_EC1->Draw("HIST SAME"); */
    c3->cd(2); h97_pi0_M_EC0->Draw();
    c3->cd(3); h97_pi0_M_EC1->Draw();
    c3->cd(4); 
    h97_pi0_M_EC0_EC1->Draw();
    h97_pi0_M_EC0_EC1->SetFillColorAlpha(kViolet+5, 0.3);
    h97_pi0_M_EC0_EC1->SetLineColor(kViolet+5);
    c3->SaveAs("plots.pdf");

    TCanvas *c4 = new TCanvas("c4", "Missing Mass Page 4", 1200, 800);
    c4->Divide(3, 2);
    c4->cd(1); h97_M2Miss_pR1->Draw();
    c4->cd(2); h97_M2Miss_pR2->Draw();
    c4->cd(3); h97_M2Miss_pi0->Draw();
    c4->SaveAs("plots.pdf)");  

    // Cleanup
    file->Close();
    delete file;
    return 0;
}

#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TStyle.h>

#include "TLorentzVector.h"

using namespace std;

bool isZeroVector(const TLorentzVector& vec) {
    return vec.E() == 0 && vec.Px() == 0 && vec.Py() == 0 && vec.Pz() == 0;
}

int kinFitPerformance() {
    // Open the ROOT file
    TFile *file = TFile::Open("/Users/gursimran/Desktop/merged_P09_626.root");
    if (!file || file->IsZombie()) {
        std::cerr << "Failed to open file!" << std::endl;
        return 1;
    }

    // Get the TTree
    TTree *tree = (TTree*)file->Get("USR970");
    if (!tree) {
        std::cerr << "Failed to get TTree!" << std::endl;
        return 1;
    }

    // Style options
    gStyle->SetStatColor(0);
    gStyle->SetStatStyle(0);
    gStyle->SetStatTextColor(kBlack);
    gStyle->SetStatX(0.9);
    gStyle->SetStatY(0.9);
    gStyle->SetStatFontSize(0.03);
    gStyle->SetLabelSize(0.04, "XYZ");
    gStyle->SetPadLeftMargin(0.15);
    gStyle->SetPadBottomMargin(0.15);
    gStyle->SetTitleOffset(1.3, "X");

    // Define variables
    TLorentzVector *inMu_TL = nullptr, *outMu_TL = nullptr, *gamma_TL = nullptr, *p_camera_TL = nullptr;
    TLorentzVector *inMuFit_TL = nullptr, *outMuFit_TL = nullptr, *gammaFit_TL = nullptr, *protonFit_TL = nullptr;
    TLorentzVector *cluster_TL = nullptr; 
    TVector3 *pVtx_vec = nullptr, *posRingA_vec = nullptr, *posRingB_vec = nullptr;
    TVector3 *pVtxFit_vec = nullptr, *posRingAFit_vec = nullptr, *posRingBFit_vec = nullptr, *clusterFit_vec = nullptr;

    double y, t, nu, Q2;
    int ndf_fit;
    double chi2_fit;
    bool fit_conv;

    double inMu_sigmaX, inMu_sigmaY, inMu_sigmaPx, inMu_sigmaPy, inMu_sigmaPz; 
    double outMu_sigmaX, outMu_sigmaY, outMu_sigmaPx, outMu_sigmaPy, outMu_sigmaPz;
    double gamma_sigmaX, gamma_sigmaY, gamma_sigmaE, proton_sigmaP; 
    double ringA_sigmaR, ringA_sigmaPhi, ringA_sigmaZ;
    double ringB_sigmaR, ringB_sigmaPhi, ringB_sigmaZ;

    // Set branch addresses
    tree->SetBranchAddress("inMu_TL", &inMu_TL);
    tree->SetBranchAddress("outMu_TL", &outMu_TL);
    tree->SetBranchAddress("gamma_TL", &gamma_TL);
    tree->SetBranchAddress("cluster_TL", &cluster_TL);
    tree->SetBranchAddress("p_camera_TL", &p_camera_TL);
    tree->SetBranchAddress("pVtx_vec", &pVtx_vec);
    tree->SetBranchAddress("posRingA_vec", &posRingA_vec);
    tree->SetBranchAddress("posRingB_vec", &posRingB_vec);

    tree->SetBranchAddress("inMuFit_TL", &inMuFit_TL);
    tree->SetBranchAddress("outMuFit_TL", &outMuFit_TL);
    tree->SetBranchAddress("gammaFit_TL", &gammaFit_TL);
    tree->SetBranchAddress("protonFit_TL", &protonFit_TL);
    tree->SetBranchAddress("pVtxFit_vec", &pVtxFit_vec);
    tree->SetBranchAddress("clusterFit_vec", &clusterFit_vec);
    tree->SetBranchAddress("posRingAFit_vec", &posRingAFit_vec);
    tree->SetBranchAddress("posRingBFit_vec", &posRingBFit_vec);

    tree->SetBranchAddress("y", &y);
    tree->SetBranchAddress("t", &t);
    tree->SetBranchAddress("nu", &nu);
    tree->SetBranchAddress("Q2", &Q2);
    tree->SetBranchAddress("ndf_fit", &ndf_fit);
    tree->SetBranchAddress("chi2_fit", &chi2_fit);
    tree->SetBranchAddress("fit_conv", &fit_conv);

    tree->SetBranchAddress("inMu_sigmaX", &inMu_sigmaX);
    tree->SetBranchAddress("inMu_sigmaY", &inMu_sigmaY);
    tree->SetBranchAddress("inMu_sigmaPx", &inMu_sigmaPx);
    tree->SetBranchAddress("inMu_sigmaPy", &inMu_sigmaPy);
    tree->SetBranchAddress("inMu_sigmaPz", &inMu_sigmaPz);

    tree->SetBranchAddress("outMu_sigmaX", &outMu_sigmaX);
    tree->SetBranchAddress("outMu_sigmaY", &outMu_sigmaY);
    tree->SetBranchAddress("outMu_sigmaPx", &outMu_sigmaPx);
    tree->SetBranchAddress("outMu_sigmaPy", &outMu_sigmaPy);
    tree->SetBranchAddress("outMu_sigmaPz", &outMu_sigmaPz);

    tree->SetBranchAddress("gamma_sigmaX", &gamma_sigmaX);
    tree->SetBranchAddress("gamma_sigmaY", &gamma_sigmaY);
    tree->SetBranchAddress("gamma_sigmaE", &gamma_sigmaE);
    tree->SetBranchAddress("proton_sigmaP", &proton_sigmaP);

    tree->SetBranchAddress("ringA_sigmaR", &ringA_sigmaR);
    tree->SetBranchAddress("ringA_sigmaPhi", &ringA_sigmaPhi);
    tree->SetBranchAddress("ringA_sigmaZ", &ringA_sigmaZ);

    tree->SetBranchAddress("ringB_sigmaR", &ringB_sigmaR);
    tree->SetBranchAddress("ringB_sigmaPhi", &ringB_sigmaPhi);
    tree->SetBranchAddress("ringB_sigmaZ", &ringB_sigmaZ);

    // Histograms

    TH1F *h97_CL = new TH1F("h97_CL", "Confidence Level; CL; Counts", 100, 0, 10);
    TH1F *h97_pull_inMuX = new TH1F("h97_pull_inMuX", "Pull: inMu X; #Delta X_{#mu} / #sigma_{X_{#mu}} ; Counts", 100, -10, 10);
    TH1F *h97_pull_inMuY = new TH1F("h97_pull_inMuY", "Pull: inMu Y; #Delta Y_{#mu} / #sigma_{Y_{#mu}} ; Counts", 100, -10, 10);
    TH1F *h97_pull_inMuPx = new TH1F("h97_pull_inMuPx", "Pull: inMu Px; #Delta Px_{#mu} / #sigma_{Px_{#mu}} ; Counts", 100, -10, 10);
    TH1F *h97_pull_inMuPy = new TH1F("h97_pull_inMuPy", "Pull: inMu Py; #Delta Py_{#mu} / #sigma_{Py_{#mu}} ; Counts", 100, -10, 10);
    TH1F *h97_pull_inMuPz = new TH1F("h97_pull_inMuPz", "Pull: inMu Pz; #Delta Pz_{#mu} / #sigma_{Pz_{#mu}} ; Counts", 100, -10, 10);

    TH1F *h97_pull_outMuX = new TH1F("h97_pull_outMuX", "Pull: outMu X; #Delta X_{#mu'} / #sigma_{X_{#mu'}} ; Counts", 100, -10, 10);
    TH1F *h97_pull_outMuY = new TH1F("h97_pull_outMuY", "Pull: outMu Y; #Delta Y_{#mu} / #sigma_{Y_{#mu}} ; Counts", 100, -10, 10);
    TH1F *h97_pull_outMuPx = new TH1F("h97_pull_outMuPx", "Pull: outMu Px; #Delta Px_{#mu} / #sigma_{Px_{#mu}} ; Counts", 100, -10, 10);
    TH1F *h97_pull_outMuPy = new TH1F("h97_pull_outMuPy", "Pull: outMu Py; #Delta Py_{#mu'} / #sigma_{Py_{#mu'}} ; Counts", 100, -10, 10);
    TH1F *h97_pull_outMuPz = new TH1F("h97_pull_outMuPz", "Pull: outMu Pz; #Delta Pz_{#mu'} / #sigma_{Pz_{#mu'}} ; Counts", 100, -10, 10);

    TH1F *h97_pull_gammaX = new TH1F("h97_pull_gammaX", "Pull: gamma X; #Delta X_{#gamma} / #sigma_{X_{#gamma}} ; Counts", 100, -10, 10);
    TH1F *h97_pull_gammaY = new TH1F("h97_pull_gammaY", "Pull: gamma Y; #Delta Y_{#gamma} / #sigma_{Y_{#gamma}} ; Counts", 100, -10, 10);
    TH1F *h97_pull_gammaE = new TH1F("h97_pull_gammaE", "Pull: gamma E; #Delta E_{#gamma} / #sigma_{E_{#gamma}} ; Counts", 100, -10, 10);
    TH1F *h97_pull_protonP = new TH1F("h97_pull_protonP", "Pull: proton P; #Delta P_{p} / #sigma_{P_{p}} ; Counts", 100, -10, 10);

    TH1F *h97_pull_ringAR = new TH1F("h97_pull_ringAR", "Pull: ring_A R; #Delta R_{A} / #sigma_{R_{A}} ; Counts", 100, -10, 10);
    TH1F *h97_pull_ringAPhi = new TH1F("h97_pull_ringAPhi", "Pull: ring_A #phi; #Delta #phi_{A} / #sigma_{#phi_{A}} ; Counts", 100, -10, 10);
    TH1F *h97_pull_ringAZ = new TH1F("h97_pull_ringAZ", "Pull: ring_A Z; #Delta Z_{A} / #sigma_{Z_{A}} ; Counts", 100, -10, 10);

    TH1F *h97_pull_ringBR = new TH1F("h97_pull_ringBR", "Pull: ring_B R; #Delta R_{B} / #sigma_{R_{B}} ; Counts", 100, -10, 10);
    TH1F *h97_pull_ringBPhi = new TH1F("h97_pull_ringBPhi", "Pull: ring_B #phi; #Delta #phi_{B} / #sigma_{#phi_{B}} ; Counts", 100, -10, 10);
    TH1F *h97_pull_ringBZ = new TH1F("h97_pull_ringBZ", "Pull: ring_B Z; #Delta Z_{B} / #sigma_{Z_{B}} ; Counts", 100, -10, 10);

    // Loop over entries
    Long64_t nEntries = tree->GetEntries();
    for (Long64_t i = 0; i < nEntries; ++i) {
        tree->GetEntry(i);

        if (isZeroVector(*inMu_TL) || isZeroVector(*outMu_TL) || isZeroVector(*p_camera_TL))
            continue;

        bool excl_kin = (Q2 > 1 && Q2 < 10) && (y > 0.05 && y < 0.9) && (nu > 10 && nu < 144) && (t > -0.64 && t < -0.08);
        if (!(excl_kin && fit_conv && chi2_fit > 10))
            continue;

        double conf_level = TMath::Prob(chi2_fit, ndf_fit);
        h97_CL->Fill(conf_level); 


        // Fill pull histograms (example for inMu X)
        //double pull_inMuX = (pVtx_vec->X() - pVtx_vec->X()) / inMu_sigmaX;
        //double pull_inMuY = (pVtx_vec->Y() - pVtx_vec->Y()) / inMu_sigmaY;
        double pull_inMuPx = (inMu_TL->Px() - inMuFit_TL->Px()) / inMu_sigmaPx;
        double pull_inMuPy = (inMu_TL->Py() - inMuFit_TL->Py()) / inMu_sigmaPy;
        double pull_inMuPz = (inMu_TL->Pz() - inMuFit_TL->Pz()) / inMu_sigmaPz;

        double pull_outMuPx = (outMu_TL->Px() - outMuFit_TL->Px()) / outMu_sigmaPx;
        double pull_outMuPy = (outMu_TL->Py() - outMuFit_TL->Py()) / outMu_sigmaPy;
        double pull_outMuPz = (outMu_TL->Pz() - outMuFit_TL->Pz()) / outMu_sigmaPz;

        //h97_pull_inMuX->Fill(pull_inMuX);
        //h97_pull_inMuY->Fill(pull_inMuY);
        h97_pull_inMuPx->Fill(pull_inMuPx);
        h97_pull_inMuPy->Fill(pull_inMuPy);
        h97_pull_inMuPz->Fill(pull_inMuPz);

        //h97_pull_outMuX->Fill(pull_outMuX);
        //h97_pull_outMuY->Fill(pull_outMuY);
        h97_pull_outMuPx->Fill(pull_outMuPx);
        h97_pull_outMuPy->Fill(pull_outMuPy);
        h97_pull_outMuPz->Fill(pull_outMuPz);

        double pull_gammaX = (cluster_TL->X() - clusterFit_vec->X()) / gamma_sigmaX; 
        double pull_gammaY = (cluster_TL->Y() - clusterFit_vec->Y()) / gamma_sigmaY;
        double pull_gammaE = (gamma_TL->E() - gammaFit_TL->E()) / gamma_sigmaE;
        double pull_protonP = (p_camera_TL->P() - protonFit_TL->P()) / proton_sigmaP; 
        h97_pull_gammaX->Fill(pull_gammaX); 
        h97_pull_gammaY->Fill(pull_gammaY); 
        h97_pull_gammaE->Fill(pull_gammaE); 
        h97_pull_protonP->Fill(pull_protonP);


        //double pull_ringAZ = (posRingA_vec->Z() - posRingAFit_vec->Z()) / ringA_sigmaZ; 
        //h97_pull_ringAZ->Fill(pull_ringAZ); 

    }

    // First canvas: first 5 histograms
    TCanvas *c1 = new TCanvas("c1", "Pull Distributions Page 1", 1200, 800);
    c1->Divide(3, 2);
    c1->cd(1); h97_CL->Draw(); 
    c1->cd(2); h97_pull_inMuX->Draw();
    c1->cd(3); h97_pull_inMuY->Draw();
    c1->cd(4); h97_pull_inMuPx->Draw();
    c1->cd(5); h97_pull_inMuPy->Draw();
    c1->cd(6); h97_pull_inMuPz->Draw();
    c1->SaveAs("pulls.pdf(");  // Open multi-page PDF

    // Second canvas: remaining 5 histograms
    TCanvas *c2 = new TCanvas("c2", "Pull Distributions Page 2", 1200, 800);
    c2->Divide(3, 2);
    c2->cd(1); h97_pull_outMuX->Draw();
    c2->cd(2); h97_pull_outMuY->Draw();
    c2->cd(3); h97_pull_outMuPx->Draw();
    c2->cd(4); h97_pull_outMuPy->Draw();
    c2->cd(5); h97_pull_outMuPz->Draw();
    c2->SaveAs("pulls.pdf");

    TCanvas *c3 = new TCanvas("c3", "Pull Distributions Page 3", 1200, 800);
    c3->Divide(3, 2);
    c3->cd(1); h97_pull_gammaX->Draw();
    c3->cd(2); h97_pull_gammaY->Draw();
    c3->cd(3); h97_pull_gammaE->Draw();
    c3->cd(4); h97_pull_protonP->Draw();
    c3->SaveAs("pulls.pdf");

    TCanvas *c4 = new TCanvas("c4", "Pull Distributions Page 4", 1200, 800);
    c4->Divide(3, 2);
    c4->cd(1); h97_pull_ringAR->Draw();
    c4->cd(2); h97_pull_ringAPhi->Draw();
    c4->cd(3); h97_pull_ringAZ->Draw();
    c4->cd(4); h97_pull_ringBR->Draw();
    c4->cd(5); h97_pull_ringBPhi->Draw();
    c4->cd(6); h97_pull_ringBZ->Draw();
    c4->SaveAs("pulls.pdf)");  // Close multi-page PDF

    file->Close();
    return 0;
}

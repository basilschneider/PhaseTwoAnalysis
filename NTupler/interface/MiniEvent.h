#ifndef _minievent_h_
#define _minievent_h_
// -*- C++ -*-
//
// Package:     PhaseTwoAnalysis/NTupler
// Class:       MiniEvent
// Description: Define the structure of ntuples

#include "TTree.h"
#include "TH2D.h"

struct MiniEvent_t
{
    MiniEvent_t()
    {
        ;
    }

    Int_t run, event, lumi;

    // Cutflow control
    //static constexpr bool fill_rle = false;
    static constexpr bool fill_vld = false;
    static constexpr bool event_by_event_comparison_primary = false;
    static constexpr bool event_by_event_comparison_secondary = false;

    // Cut variables
    static constexpr double el_pt_lo = 5.;
    static constexpr double el_pt_hi = 30.;
    static constexpr double mu_pt_lo = 5.;
    static constexpr double mu_pt_hi = 30.;
    static constexpr double jet_pt_lo = 25.;
    static constexpr double mass_el = .000511;
    static constexpr double mass_mu = .105658;
    static constexpr double iso_cut_rel = .5;
    static constexpr double iso_cut_abs = 5.;
    static constexpr double truth_match_diff_pt_rel = .5;
    static constexpr double truth_match_diff_dr = .1;

    double genWeight;
    Int_t nLep, nMu/*, nEl*/;
    Int_t nSoftLep, nSoftMu/*, nSoftEl*/;
    Int_t nJet25, nJet40, nJet60, nJet100, nJet150, nBJet;

    //std::vector<double> el1_pt, el1_eta, el1_phi, el2_pt, el2_eta, el2_phi;
    //std::vector<int> el1_q, el2_q;
    std::vector<double> mu_pt, mu_eta, mu_phi;
    std::vector<int> mu_q, mu_mother, mu_directmother;
    std::vector<bool> mu_matched, mu_st20to30;
    std::vector<double> mu_woIso_pt, mu_woIso_eta, mu_woIso_phi;
    std::vector<double> mu_pt_truth, mu_eta_truth, mu_phi_truth;
    //std::vector<double> lep1_pt, lep1_eta, lep1_phi, lep2_pt, lep2_eta, lep2_phi;
    //std::vector<int> lep1_mass, lep2_mass;

    //std::vector<double> el1_pt_truth, el1_eta_truth, el1_phi_truth, el2_pt_truth, el2_eta_truth, el2_phi_truth;
    //std::vector<int> el1_q_truth, el2_q_truth;
    //std::vector<double> mu1_pt_truth, mu1_eta_truth, mu1_phi_truth, mu2_pt_truth, mu2_eta_truth, mu2_phi_truth;
    //std::vector<int> mu1_q_truth, mu2_q_truth;
    //std::vector<double> lep1_pt_truth, lep1_eta_truth, lep1_phi_truth, lep2_pt_truth, lep2_eta_truth, lep2_phi_truth, lep1_mass_truth, lep2_mass_truth;

    std::vector<double> jet_pt, jet_eta, jet_phi, jet_mass;
    //std::vector<double> jet_pt_truth, jet_eta_truth, jet_phi_truth, jet_mass_truth;

    double met, met_eta, met_phi, genmet, genmet_eta, genmet_phi;
    double ht25, ht40, ht60, ht100, ht150;
    double genht25, genht40;
    double mlt, mht25, mht40, mhlt25, mhlt40;
    std::vector<double> mllMin, mllMax, mt1, mt2, pt2l, mtautau;

    unsigned int ZtoLL;
    bool crazyMuon50, crazyMuon200, crazyMuon500;

    //std::vector<double> vld_el_pt, vld_el_iso_pt, vld_el_good_pt, vld_el_is_tight;
    //std::vector<double> vld_el_hs_pt, vld_el_iso_hs_pt, vld_el_good_hs_pt;
    //std::vector<double> vld_el_py8_pt, vld_el_iso_py8_pt, vld_el_good_py8_pt;
    //std::vector<double> vld_el_others_pt, vld_el_iso_others_pt, vld_el_good_others_pt;
    //std::vector<double> vld_el_eta, vld_el_iso_eta, vld_el_good_eta;
    //std::vector<double> vld_el_hs_eta, vld_el_iso_hs_eta, vld_el_good_hs_eta;
    //std::vector<double> vld_el_py8_eta, vld_el_iso_py8_eta, vld_el_good_py8_eta;
    //std::vector<double> vld_el_others_eta, vld_el_iso_others_eta, vld_el_good_others_eta;

    //std::vector<double> vld_mu_pt, vld_mu_iso_pt, vld_mu_good_pt, vld_mu_is_tight;
    //std::vector<double> vld_mu_hs_pt, vld_mu_iso_hs_pt, vld_mu_good_hs_pt;
    //std::vector<double> vld_mu_py8_pt, vld_mu_iso_py8_pt, vld_mu_good_py8_pt;
    //std::vector<double> vld_mu_others_pt, vld_mu_iso_others_pt, vld_mu_good_others_pt;
    //std::vector<double> vld_mu_eta, vld_mu_iso_eta, vld_mu_good_eta;
    //std::vector<double> vld_mu_hs_eta, vld_mu_iso_hs_eta, vld_mu_good_hs_eta;
    //std::vector<double> vld_mu_py8_eta, vld_mu_iso_py8_eta, vld_mu_good_py8_eta;
    //std::vector<double> vld_mu_others_eta, vld_mu_iso_others_eta, vld_mu_good_others_eta;
    std::vector<double> vld_mu_hs_pt, vld_mu_py8_pt, vld_mu_others_pt;
    std::vector<double> vld_mu_iso_hs_pt, vld_mu_iso_py8_pt, vld_mu_iso_others_pt;
    std::vector<double> vld_mu_good_hs_pt, vld_mu_good_py8_pt, vld_mu_good_others_pt;
    std::vector<double> vld_mu_hs_eta, vld_mu_py8_eta, vld_mu_others_eta;
    std::vector<double> vld_mu_iso_hs_eta, vld_mu_iso_py8_eta, vld_mu_iso_others_eta;
    std::vector<double> vld_mu_good_hs_eta, vld_mu_good_py8_eta, vld_mu_good_others_eta;

    //std::vector<double> vld_genel_pt, vld_genel_hs_pt, vld_genel_py8_pt;
    //std::vector<double> vld_genel_eta, vld_genel_hs_eta, vld_genel_py8_eta;
    //std::vector<double> vld_genmu_pt, vld_genmu_hs_pt, vld_genmu_py8_pt;
    //std::vector<double> vld_genmu_eta, vld_genmu_hs_eta, vld_genmu_py8_eta;
    std::vector<double> vld_genmu_hs_pt, vld_genmu_py8_pt;
    std::vector<double> vld_genmu_hs_eta, vld_genmu_py8_eta;

    //std::vector<double> vld_el_absiso, vld_el_iso_absiso, vld_el_good_absiso;
    //std::vector<double> vld_el_hs_absiso, vld_el_iso_hs_absiso, vld_el_good_hs_absiso;
    //std::vector<double> vld_el_py8_absiso, vld_el_iso_py8_absiso, vld_el_good_py8_absiso;
    //std::vector<double> vld_el_others_absiso, vld_el_iso_others_absiso, vld_el_good_others_absiso;
    //std::vector<double> vld_el_reliso, vld_el_iso_reliso, vld_el_good_reliso;
    //std::vector<double> vld_el_hs_reliso, vld_el_iso_hs_reliso, vld_el_good_hs_reliso;
    //std::vector<double> vld_el_py8_reliso, vld_el_iso_py8_reliso, vld_el_good_py8_reliso;
    //std::vector<double> vld_el_others_reliso, vld_el_iso_others_reliso, vld_el_good_others_reliso;
    //std::vector<double> vld_el_dxy, vld_el_iso_dxy, vld_el_good_dxy;
    //std::vector<double> vld_el_hs_dxy, vld_el_iso_hs_dxy, vld_el_good_hs_dxy;
    //std::vector<double> vld_el_py8_dxy, vld_el_iso_py8_dxy, vld_el_good_py8_dxy;
    //std::vector<double> vld_el_others_dxy, vld_el_iso_others_dxy, vld_el_good_others_dxy;
    //std::vector<double> vld_el_dz, vld_el_iso_dz, vld_el_good_dz;
    //std::vector<double> vld_el_hs_dz, vld_el_iso_hs_dz, vld_el_good_hs_dz;
    //std::vector<double> vld_el_py8_dz, vld_el_iso_py8_dz, vld_el_good_py8_dz;
    //std::vector<double> vld_el_others_dz, vld_el_iso_others_dz, vld_el_good_others_dz;

    //std::vector<double> vld_mu_absiso, vld_mu_iso_absiso, vld_mu_good_absiso;
    //std::vector<double> vld_mu_hs_absiso, vld_mu_iso_hs_absiso, vld_mu_good_hs_absiso;
    //std::vector<double> vld_mu_py8_absiso, vld_mu_iso_py8_absiso, vld_mu_good_py8_absiso;
    //std::vector<double> vld_mu_others_absiso, vld_mu_iso_others_absiso, vld_mu_good_others_absiso;
    //std::vector<double> vld_mu_reliso, vld_mu_iso_reliso, vld_mu_good_reliso;
    //std::vector<double> vld_mu_hs_reliso, vld_mu_iso_hs_reliso, vld_mu_good_hs_reliso;
    //std::vector<double> vld_mu_py8_reliso, vld_mu_iso_py8_reliso, vld_mu_good_py8_reliso;
    //std::vector<double> vld_mu_others_reliso, vld_mu_iso_others_reliso, vld_mu_good_others_reliso;
    //std::vector<double> vld_mu_dxy, vld_mu_iso_dxy, vld_mu_good_dxy;
    //std::vector<double> vld_mu_hs_dxy, vld_mu_iso_hs_dxy, vld_mu_good_hs_dxy;
    //std::vector<double> vld_mu_py8_dxy, vld_mu_iso_py8_dxy, vld_mu_good_py8_dxy;
    //std::vector<double> vld_mu_others_dxy, vld_mu_iso_others_dxy, vld_mu_good_others_dxy;
    //std::vector<double> vld_mu_dz, vld_mu_iso_dz, vld_mu_good_dz;
    //std::vector<double> vld_mu_hs_dz, vld_mu_iso_hs_dz, vld_mu_good_hs_dz;
    //std::vector<double> vld_mu_py8_dz, vld_mu_iso_py8_dz, vld_mu_good_py8_dz;
    //std::vector<double> vld_mu_others_dz, vld_mu_iso_others_dz, vld_mu_good_others_dz;
    std::vector<double> vld_mu_hs_reliso, vld_mu_py8_reliso, vld_mu_others_reliso;
    std::vector<double> vld_mu_iso_hs_reliso, vld_mu_iso_py8_reliso, vld_mu_iso_others_reliso;
    std::vector<double> vld_mu_good_hs_reliso, vld_mu_good_py8_reliso, vld_mu_good_others_reliso;

    //TH2D* vld_el_hs_pt_eta = new TH2D("vld_el_hs_pt_eta", "vld_el_hs_pt_eta", 36, 0., 90., 40, 0., 4.);
    //TH2D* vld_el_iso_hs_pt_eta = new TH2D("vld_el_iso_hs_pt_eta", "vld_el_iso_hs_pt_eta", 36, 0., 90., 40, 0., 4.);
    //TH2D* vld_el_good_hs_pt_eta = new TH2D("vld_el_good_hs_pt_eta", "vld_el_good_hs_pt_eta", 36, 0., 90., 40, 0., 4.);
    //TH2D* vld_el_py8_pt_eta = new TH2D("vld_el_py8_pt_eta", "vld_el_py8_pt_eta", 36, 0., 90., 40, 0., 4.);
    //TH2D* vld_el_iso_py8_pt_eta = new TH2D("vld_el_iso_py8_pt_eta", "vld_el_iso_py8_pt_eta", 36, 0., 90., 40, 0., 4.);
    //TH2D* vld_el_good_py8_pt_eta = new TH2D("vld_el_good_py8_pt_eta", "vld_el_good_py8_pt_eta", 36, 0., 90., 40, 0., 4.);
    //TH2D* vld_el_others_pt_eta = new TH2D("vld_el_others_pt_eta", "vld_el_others_pt_eta", 36, 0., 90., 40, 0., 4.);
    //TH2D* vld_el_iso_others_pt_eta = new TH2D("vld_el_iso_others_pt_eta", "vld_el_iso_others_pt_eta", 36, 0., 90., 40, 0., 4.);
    //TH2D* vld_el_good_others_pt_eta = new TH2D("vld_el_good_others_pt_eta", "vld_el_good_others_pt_eta", 36, 0., 90., 40, 0., 4.);
    //TH2D* vld_mu_hs_pt_eta = new TH2D("vld_mu_hs_pt_eta", "vld_mu_hs_pt_eta", 36, 0., 90., 40, 0., 4.);
    //TH2D* vld_mu_iso_hs_pt_eta = new TH2D("vld_mu_iso_hs_pt_eta", "vld_mu_iso_hs_pt_eta", 36, 0., 90., 40, 0., 4.);
    //TH2D* vld_mu_good_hs_pt_eta = new TH2D("vld_mu_good_hs_pt_eta", "vld_mu_good_hs_pt_eta", 36, 0., 90., 40, 0., 4.);
    //TH2D* vld_mu_py8_pt_eta = new TH2D("vld_mu_py8_pt_eta", "vld_mu_py8_pt_eta", 36, 0., 90., 40, 0., 4.);
    //TH2D* vld_mu_iso_py8_pt_eta = new TH2D("vld_mu_iso_py8_pt_eta", "vld_mu_iso_py8_pt_eta", 36, 0., 90., 40, 0., 4.);
    //TH2D* vld_mu_good_py8_pt_eta = new TH2D("vld_mu_good_py8_pt_eta", "vld_mu_good_py8_pt_eta", 36, 0., 90., 40, 0., 4.);
    //TH2D* vld_mu_others_pt_eta = new TH2D("vld_mu_others_pt_eta", "vld_mu_others_pt_eta", 36, 0., 90., 40, 0., 4.);
    //TH2D* vld_mu_iso_others_pt_eta = new TH2D("vld_mu_iso_others_pt_eta", "vld_mu_iso_others_pt_eta", 36, 0., 90., 40, 0., 4.);
    //TH2D* vld_mu_good_others_pt_eta = new TH2D("vld_mu_good_others_pt_eta", "vld_mu_good_others_pt_eta", 36, 0., 90., 40, 0., 4.);

    ////TH2D* vld_el_pt_iso_abs = new TH2D("vld_el_pt_iso_abs", "vld_el_pt_iso_abs", 12, 0., 30., 40, 0., 20.);
    ////TH2D* vld_el_iso_pt_iso_abs = new TH2D("vld_el_iso_pt_iso_abs", "vld_el_iso_pt_iso_abs", 12, 0., 30., 40, 0., 20.);
    ////TH2D* vld_el_pt_iso_rel = new TH2D("vld_el_pt_iso_rel", "vld_el_pt_iso_rel", 12, 0., 30., 40, 0., 2.);
    ////TH2D* vld_el_iso_pt_iso_rel = new TH2D("vld_el_iso_pt_iso_rel", "vld_el_iso_pt_iso_rel", 12, 0., 30., 40, 0., 2.);
    ////TH2D* vld_el_pt_dxy = new TH2D("vld_el_pt_dxy", "vld_el_pt_dxy", 12, 0., 30., 40 ,0., .1);
    ////TH2D* vld_el_iso_pt_dxy = new TH2D("vld_el_iso_pt_dxy", "vld_el_iso_pt_dxy", 12, 0., 30., 40 ,0., .1);
    ////TH2D* vld_el_pt_dz = new TH2D("vld_el_pt_dz", "vld_el_pt_dz", 12, 0., 30., 40 ,0., .1);
    ////TH2D* vld_el_iso_pt_dz = new TH2D("vld_el_iso_pt_dz", "vld_el_iso_pt_dz", 12, 0., 30., 40 ,0., .1);
    ////TH2D* vld_mu_pt_iso_abs = new TH2D("vld_mu_pt_iso_abs", "vld_mu_pt_iso_abs", 12, 0., 30., 40, 0., 20.);
    ////TH2D* vld_mu_iso_pt_iso_abs = new TH2D("vld_mu_iso_pt_iso_abs", "vld_mu_iso_pt_iso_abs", 12, 0., 30., 40, 0., 20.);
    ////TH2D* vld_mu_pt_iso_rel = new TH2D("vld_mu_pt_iso_rel", "vld_mu_pt_iso_rel", 12, 0., 30., 40, 0., 2.);
    ////TH2D* vld_mu_iso_pt_iso_rel = new TH2D("vld_mu_iso_pt_iso_rel", "vld_mu_iso_pt_iso_rel", 12, 0., 30., 40, 0., 2.);
    ////TH2D* vld_mu_pt_dxy = new TH2D("vld_mu_pt_dxy", "vld_mu_pt_dxy", 12, 0., 30., 40 ,0., .1);
    ////TH2D* vld_mu_iso_pt_dxy = new TH2D("vld_mu_iso_pt_dxy", "vld_mu_iso_pt_dxy", 12, 0., 30., 40 ,0., .1);
    ////TH2D* vld_mu_pt_dz = new TH2D("vld_mu_pt_dz", "vld_mu_pt_dz", 12, 0., 30., 40 ,0., .1);
    ////TH2D* vld_mu_iso_pt_dz = new TH2D("vld_mu_iso_pt_dz", "vld_mu_iso_pt_dz", 12, 0., 30., 40 ,0., .1);

    // Real lepton efficiency histograms
    //TH2D* rle_el_num = new TH2D("rle_el_num", "rle_el_num", 6, 0., 30., 8, 0., 4.);
    //TH2D* rle_el_den = new TH2D("rle_el_den", "rle_el_den", 6, 0., 30., 8, 0., 4.);
    //TH2D* rle_mu_num = new TH2D("rle_mu_num", "rle_mu_num", 6, 0., 30., 8, 0., 4.);
    //TH2D* rle_mu_den = new TH2D("rle_mu_den", "rle_mu_den", 6, 0., 30., 8, 0., 4.);
};

//void createMiniEventTree(TTree *t_event_, TTree *t_genParts_, TTree *t_vertices_, TTree *t_genJets_, TTree *t_looseElecs_, TTree *t_tightElecs_, TTree *t_looseMuons_, TTree *t_tightMuons_, TTree *t_puppiJets_, TTree *t_puppiMET_,MiniEvent_t &ev);
void createMiniEventTree(TTree *t_tree_, MiniEvent_t &ev);

#endif

#include "PhaseTwoAnalysis/NTupler/interface/MiniEvent.h"

void createMiniEventTree(TTree* t_tree_, MiniEvent_t &ev)
{
    t_tree_->Branch("Run", &ev.run);
    t_tree_->Branch("Event", &ev.event);
    t_tree_->Branch("Lumi", &ev.lumi);

    t_tree_->Branch("genWeight", &ev.genWeight);
    // ntot
    // xs
    t_tree_->Branch("el1_pt", &ev.el1_pt);
    t_tree_->Branch("el1_eta", &ev.el1_eta);
    t_tree_->Branch("el1_phi", &ev.el1_phi);
    t_tree_->Branch("el1_q", &ev.el1_q);
    t_tree_->Branch("el2_pt", &ev.el2_pt);
    t_tree_->Branch("el2_eta", &ev.el2_eta);
    t_tree_->Branch("el2_phi", &ev.el2_phi);
    t_tree_->Branch("el2_q", &ev.el2_q);
    t_tree_->Branch("mu1_pt", &ev.mu1_pt);
    t_tree_->Branch("mu1_eta", &ev.mu1_eta);
    t_tree_->Branch("mu1_phi", &ev.mu1_phi);
    t_tree_->Branch("mu1_q", &ev.mu1_q);
    t_tree_->Branch("mu1_mother", &ev.mu1_mother);
    t_tree_->Branch("mu2_pt", &ev.mu2_pt);
    t_tree_->Branch("mu2_eta", &ev.mu2_eta);
    t_tree_->Branch("mu2_phi", &ev.mu2_phi);
    t_tree_->Branch("mu2_q", &ev.mu2_q);
    t_tree_->Branch("mu2_mother", &ev.mu2_mother);
    t_tree_->Branch("lep1_pt", &ev.lep1_pt);
    t_tree_->Branch("lep1_eta", &ev.lep1_eta);
    t_tree_->Branch("lep1_phi", &ev.lep1_phi);
    t_tree_->Branch("lep1_mass", &ev.lep1_mass);
    t_tree_->Branch("lep2_pt", &ev.lep2_pt);
    t_tree_->Branch("lep2_eta", &ev.lep2_eta);
    t_tree_->Branch("lep2_phi", &ev.lep2_phi);
    t_tree_->Branch("lep2_mass", &ev.lep2_mass);

    //t_tree_->Branch("el1_pt_truth", &ev.el1_pt_truth);
    //t_tree_->Branch("el1_eta_truth", &ev.el1_eta_truth);
    //t_tree_->Branch("el1_phi_truth", &ev.el1_phi_truth);
    //t_tree_->Branch("el1_q_truth", &ev.el1_q_truth);
    //t_tree_->Branch("el2_pt_truth", &ev.el2_pt_truth);
    //t_tree_->Branch("el2_eta_truth", &ev.el2_eta_truth);
    //t_tree_->Branch("el2_phi_truth", &ev.el2_phi_truth);
    //t_tree_->Branch("el2_q_truth", &ev.el2_q_truth);
    //t_tree_->Branch("mu1_pt_truth", &ev.mu1_pt_truth);
    //t_tree_->Branch("mu1_eta_truth", &ev.mu1_eta_truth);
    //t_tree_->Branch("mu1_phi_truth", &ev.mu1_phi_truth);
    //t_tree_->Branch("mu1_q_truth", &ev.mu1_q_truth);
    //t_tree_->Branch("mu2_pt_truth", &ev.mu2_pt_truth);
    //t_tree_->Branch("mu2_eta_truth", &ev.mu2_eta_truth);
    //t_tree_->Branch("mu2_phi_truth", &ev.mu2_phi_truth);
    //t_tree_->Branch("mu2_q_truth", &ev.mu2_q_truth);
    //t_tree_->Branch("lep1_pt_truth", &ev.lep1_pt_truth);
    //t_tree_->Branch("lep1_eta_truth", &ev.lep1_eta_truth);
    //t_tree_->Branch("lep1_phi_truth", &ev.lep1_phi_truth);
    //t_tree_->Branch("lep1_mass_truth", &ev.lep1_mass_truth);
    //t_tree_->Branch("lep2_pt_truth", &ev.lep2_pt_truth);
    //t_tree_->Branch("lep2_eta_truth", &ev.lep2_eta_truth);
    //t_tree_->Branch("lep2_phi_truth", &ev.lep2_phi_truth);
    //t_tree_->Branch("lep2_mass_truth", &ev.lep2_mass_truth);

    t_tree_->Branch("jet1_pt", &ev.jet1_pt);
    t_tree_->Branch("jet1_eta", &ev.jet1_eta);
    t_tree_->Branch("jet1_phi", &ev.jet1_phi);
    t_tree_->Branch("jet1_mass", &ev.jet1_mass);

    t_tree_->Branch("jet1_pt_truth", &ev.jet1_pt_truth);
    t_tree_->Branch("jet1_eta_truth", &ev.jet1_eta_truth);
    t_tree_->Branch("jet1_phi_truth", &ev.jet1_phi_truth);
    t_tree_->Branch("jet1_mass_truth", &ev.jet1_mass_truth);

    t_tree_->Branch("nLep", &ev.nLep);
    t_tree_->Branch("nEl", &ev.nEl);
    t_tree_->Branch("nMu", &ev.nMu);
    t_tree_->Branch("nSoftLep", &ev.nSoftLep);
    t_tree_->Branch("nSoftEl", &ev.nSoftEl);
    t_tree_->Branch("nSoftMu", &ev.nSoftMu);
    t_tree_->Branch("nJet", &ev.nJet);
    t_tree_->Branch("nBJet", &ev.nBJet);

    t_tree_->Branch("met", &ev.met);
    t_tree_->Branch("genmet", &ev.genmet);
    t_tree_->Branch("ht", &ev.ht);
    t_tree_->Branch("mllMin", &ev.mllMin);
    t_tree_->Branch("mllMax", &ev.mllMax);
    t_tree_->Branch("mt1", &ev.mt1);
    t_tree_->Branch("mt2", &ev.mt2);
    t_tree_->Branch("pt2l", &ev.pt2l);
    t_tree_->Branch("mtautau", &ev.mtautau);

    if (ev.fill_vld){
        //t_tree_->Branch("vld_el_pt", &ev.vld_el_pt);
        //t_tree_->Branch("vld_el_iso_pt", &ev.vld_el_iso_pt);
        //t_tree_->Branch("vld_el_good_pt", &ev.vld_el_good_pt);
        //t_tree_->Branch("vld_el_is_tight", &ev.vld_el_is_tight);
        //t_tree_->Branch("vld_el_hs_pt", &ev.vld_el_hs_pt);
        //t_tree_->Branch("vld_el_iso_hs_pt", &ev.vld_el_iso_hs_pt);
        //t_tree_->Branch("vld_el_good_hs_pt", &ev.vld_el_good_hs_pt);
        //t_tree_->Branch("vld_el_py8_pt", &ev.vld_el_py8_pt);
        //t_tree_->Branch("vld_el_iso_py8_pt", &ev.vld_el_iso_py8_pt);
        //t_tree_->Branch("vld_el_good_py8_pt", &ev.vld_el_good_py8_pt);
        //t_tree_->Branch("vld_el_others_pt", &ev.vld_el_others_pt);
        //t_tree_->Branch("vld_el_iso_others_pt", &ev.vld_el_iso_others_pt);
        //t_tree_->Branch("vld_el_good_others_pt", &ev.vld_el_good_others_pt);
        //t_tree_->Branch("vld_el_eta", &ev.vld_el_eta);
        //t_tree_->Branch("vld_el_iso_eta", &ev.vld_el_iso_eta);
        //t_tree_->Branch("vld_el_good_eta", &ev.vld_el_good_eta);
        //t_tree_->Branch("vld_el_hs_eta", &ev.vld_el_hs_eta);
        //t_tree_->Branch("vld_el_iso_hs_eta", &ev.vld_el_iso_hs_eta);
        //t_tree_->Branch("vld_el_good_hs_eta", &ev.vld_el_good_hs_eta);
        //t_tree_->Branch("vld_el_py8_eta", &ev.vld_el_py8_eta);
        //t_tree_->Branch("vld_el_iso_py8_eta", &ev.vld_el_iso_py8_eta);
        //t_tree_->Branch("vld_el_good_py8_eta", &ev.vld_el_good_py8_eta);
        //t_tree_->Branch("vld_el_others_eta", &ev.vld_el_others_eta);
        //t_tree_->Branch("vld_el_iso_others_eta", &ev.vld_el_iso_others_eta);
        //t_tree_->Branch("vld_el_good_others_eta", &ev.vld_el_good_others_eta);

        //t_tree_->Branch("vld_mu_pt", &ev.vld_mu_pt);
        //t_tree_->Branch("vld_mu_iso_pt", &ev.vld_mu_iso_pt);
        //t_tree_->Branch("vld_mu_good_pt", &ev.vld_mu_good_pt);
        //t_tree_->Branch("vld_mu_is_tight", &ev.vld_mu_is_tight);
        t_tree_->Branch("vld_mu_hs_pt", &ev.vld_mu_hs_pt);
        t_tree_->Branch("vld_mu_iso_hs_pt", &ev.vld_mu_iso_hs_pt);
        t_tree_->Branch("vld_mu_good_hs_pt", &ev.vld_mu_good_hs_pt);
        t_tree_->Branch("vld_mu_py8_pt", &ev.vld_mu_py8_pt);
        t_tree_->Branch("vld_mu_iso_py8_pt", &ev.vld_mu_iso_py8_pt);
        t_tree_->Branch("vld_mu_good_py8_pt", &ev.vld_mu_good_py8_pt);
        t_tree_->Branch("vld_mu_others_pt", &ev.vld_mu_others_pt);
        t_tree_->Branch("vld_mu_iso_others_pt", &ev.vld_mu_iso_others_pt);
        t_tree_->Branch("vld_mu_good_others_pt", &ev.vld_mu_good_others_pt);
        //t_tree_->Branch("vld_mu_eta", &ev.vld_mu_eta);
        //t_tree_->Branch("vld_mu_iso_eta", &ev.vld_mu_iso_eta);
        //t_tree_->Branch("vld_mu_good_eta", &ev.vld_mu_good_eta);
        t_tree_->Branch("vld_mu_hs_eta", &ev.vld_mu_hs_eta);
        t_tree_->Branch("vld_mu_iso_hs_eta", &ev.vld_mu_iso_hs_eta);
        t_tree_->Branch("vld_mu_good_hs_eta", &ev.vld_mu_good_hs_eta);
        t_tree_->Branch("vld_mu_py8_eta", &ev.vld_mu_py8_eta);
        t_tree_->Branch("vld_mu_iso_py8_eta", &ev.vld_mu_iso_py8_eta);
        t_tree_->Branch("vld_mu_good_py8_eta", &ev.vld_mu_good_py8_eta);
        t_tree_->Branch("vld_mu_others_eta", &ev.vld_mu_others_eta);
        t_tree_->Branch("vld_mu_iso_others_eta", &ev.vld_mu_iso_others_eta);
        t_tree_->Branch("vld_mu_good_others_eta", &ev.vld_mu_good_others_eta);

        //t_tree_->Branch("vld_genel_pt", &ev.vld_genel_pt);
        //t_tree_->Branch("vld_genel_hs_pt", &ev.vld_genel_hs_pt);
        //t_tree_->Branch("vld_genel_py8_pt", &ev.vld_genel_py8_pt);
        //t_tree_->Branch("vld_genel_eta", &ev.vld_genel_eta);
        //t_tree_->Branch("vld_genel_hs_eta", &ev.vld_genel_hs_eta);
        //t_tree_->Branch("vld_genel_py8_eta", &ev.vld_genel_py8_eta);
        //t_tree_->Branch("vld_genmu_pt", &ev.vld_genmu_pt);
        t_tree_->Branch("vld_genmu_hs_pt", &ev.vld_genmu_hs_pt);
        t_tree_->Branch("vld_genmu_py8_pt", &ev.vld_genmu_py8_pt);
        //t_tree_->Branch("vld_genmu_eta", &ev.vld_genmu_eta);
        t_tree_->Branch("vld_genmu_hs_eta", &ev.vld_genmu_hs_eta);
        t_tree_->Branch("vld_genmu_py8_eta", &ev.vld_genmu_py8_eta);

        //t_tree_->Branch("vld_el_absiso", &ev.vld_el_absiso);
        //t_tree_->Branch("vld_el_iso_absiso", &ev.vld_el_iso_absiso);
        //t_tree_->Branch("vld_el_good_absiso", &ev.vld_el_good_absiso);
        //t_tree_->Branch("vld_el_hs_absiso", &ev.vld_el_hs_absiso);
        //t_tree_->Branch("vld_el_iso_hs_absiso", &ev.vld_el_iso_hs_absiso);
        //t_tree_->Branch("vld_el_good_hs_absiso", &ev.vld_el_good_hs_absiso);
        //t_tree_->Branch("vld_el_py8_absiso", &ev.vld_el_py8_absiso);
        //t_tree_->Branch("vld_el_iso_py8_absiso", &ev.vld_el_iso_py8_absiso);
        //t_tree_->Branch("vld_el_good_py8_absiso", &ev.vld_el_good_py8_absiso);
        //t_tree_->Branch("vld_el_others_absiso", &ev.vld_el_others_absiso);
        //t_tree_->Branch("vld_el_iso_others_absiso", &ev.vld_el_iso_others_absiso);
        //t_tree_->Branch("vld_el_good_others_absiso", &ev.vld_el_good_others_absiso);
        //t_tree_->Branch("vld_el_reliso", &ev.vld_el_reliso);
        //t_tree_->Branch("vld_el_iso_reliso", &ev.vld_el_iso_reliso);
        //t_tree_->Branch("vld_el_good_reliso", &ev.vld_el_good_reliso);
        //t_tree_->Branch("vld_el_hs_reliso", &ev.vld_el_hs_reliso);
        //t_tree_->Branch("vld_el_iso_hs_reliso", &ev.vld_el_iso_hs_reliso);
        //t_tree_->Branch("vld_el_good_hs_reliso", &ev.vld_el_good_hs_reliso);
        //t_tree_->Branch("vld_el_py8_reliso", &ev.vld_el_py8_reliso);
        //t_tree_->Branch("vld_el_iso_py8_reliso", &ev.vld_el_iso_py8_reliso);
        //t_tree_->Branch("vld_el_good_py8_reliso", &ev.vld_el_good_py8_reliso);
        //t_tree_->Branch("vld_el_others_reliso", &ev.vld_el_others_reliso);
        //t_tree_->Branch("vld_el_iso_others_reliso", &ev.vld_el_iso_others_reliso);
        //t_tree_->Branch("vld_el_good_others_reliso", &ev.vld_el_good_others_reliso);
        //t_tree_->Branch("vld_el_dxy", &ev.vld_el_dxy);
        //t_tree_->Branch("vld_el_iso_dxy", &ev.vld_el_iso_dxy);
        //t_tree_->Branch("vld_el_good_dxy", &ev.vld_el_good_dxy);
        //t_tree_->Branch("vld_el_hs_dxy", &ev.vld_el_hs_dxy);
        //t_tree_->Branch("vld_el_iso_hs_dxy", &ev.vld_el_iso_hs_dxy);
        //t_tree_->Branch("vld_el_good_hs_dxy", &ev.vld_el_good_hs_dxy);
        //t_tree_->Branch("vld_el_py8_dxy", &ev.vld_el_py8_dxy);
        //t_tree_->Branch("vld_el_iso_py8_dxy", &ev.vld_el_iso_py8_dxy);
        //t_tree_->Branch("vld_el_good_py8_dxy", &ev.vld_el_good_py8_dxy);
        //t_tree_->Branch("vld_el_others_dxy", &ev.vld_el_others_dxy);
        //t_tree_->Branch("vld_el_iso_others_dxy", &ev.vld_el_iso_others_dxy);
        //t_tree_->Branch("vld_el_good_others_dxy", &ev.vld_el_good_others_dxy);
        //t_tree_->Branch("vld_el_dz", &ev.vld_el_dz);
        //t_tree_->Branch("vld_el_iso_dz", &ev.vld_el_iso_dz);
        //t_tree_->Branch("vld_el_good_dz", &ev.vld_el_good_dz);
        //t_tree_->Branch("vld_el_hs_dz", &ev.vld_el_hs_dz);
        //t_tree_->Branch("vld_el_iso_hs_dz", &ev.vld_el_iso_hs_dz);
        //t_tree_->Branch("vld_el_good_hs_dz", &ev.vld_el_good_hs_dz);
        //t_tree_->Branch("vld_el_py8_dz", &ev.vld_el_py8_dz);
        //t_tree_->Branch("vld_el_iso_py8_dz", &ev.vld_el_iso_py8_dz);
        //t_tree_->Branch("vld_el_good_py8_dz", &ev.vld_el_good_py8_dz);
        //t_tree_->Branch("vld_el_others_dz", &ev.vld_el_others_dz);
        //t_tree_->Branch("vld_el_iso_others_dz", &ev.vld_el_iso_others_dz);
        //t_tree_->Branch("vld_el_good_others_dz", &ev.vld_el_good_others_dz);

        //t_tree_->Branch("vld_mu_absiso", &ev.vld_mu_absiso);
        //t_tree_->Branch("vld_mu_iso_absiso", &ev.vld_mu_iso_absiso);
        //t_tree_->Branch("vld_mu_good_absiso", &ev.vld_mu_good_absiso);
        //t_tree_->Branch("vld_mu_hs_absiso", &ev.vld_mu_hs_absiso);
        //t_tree_->Branch("vld_mu_iso_hs_absiso", &ev.vld_mu_iso_hs_absiso);
        //t_tree_->Branch("vld_mu_good_hs_absiso", &ev.vld_mu_good_hs_absiso);
        //t_tree_->Branch("vld_mu_py8_absiso", &ev.vld_mu_py8_absiso);
        //t_tree_->Branch("vld_mu_iso_py8_absiso", &ev.vld_mu_iso_py8_absiso);
        //t_tree_->Branch("vld_mu_good_py8_absiso", &ev.vld_mu_good_py8_absiso);
        //t_tree_->Branch("vld_mu_others_absiso", &ev.vld_mu_others_absiso);
        //t_tree_->Branch("vld_mu_iso_others_absiso", &ev.vld_mu_iso_others_absiso);
        //t_tree_->Branch("vld_mu_good_others_absiso", &ev.vld_mu_good_others_absiso);
        //t_tree_->Branch("vld_mu_reliso", &ev.vld_mu_reliso);
        //t_tree_->Branch("vld_mu_iso_reliso", &ev.vld_mu_iso_reliso);
        //t_tree_->Branch("vld_mu_good_reliso", &ev.vld_mu_good_reliso);
        t_tree_->Branch("vld_mu_hs_reliso", &ev.vld_mu_hs_reliso);
        t_tree_->Branch("vld_mu_iso_hs_reliso", &ev.vld_mu_iso_hs_reliso);
        t_tree_->Branch("vld_mu_good_hs_reliso", &ev.vld_mu_good_hs_reliso);
        t_tree_->Branch("vld_mu_py8_reliso", &ev.vld_mu_py8_reliso);
        t_tree_->Branch("vld_mu_iso_py8_reliso", &ev.vld_mu_iso_py8_reliso);
        t_tree_->Branch("vld_mu_good_py8_reliso", &ev.vld_mu_good_py8_reliso);
        t_tree_->Branch("vld_mu_others_reliso", &ev.vld_mu_others_reliso);
        t_tree_->Branch("vld_mu_iso_others_reliso", &ev.vld_mu_iso_others_reliso);
        t_tree_->Branch("vld_mu_good_others_reliso", &ev.vld_mu_good_others_reliso);
        //t_tree_->Branch("vld_mu_dxy", &ev.vld_mu_dxy);
        //t_tree_->Branch("vld_mu_iso_dxy", &ev.vld_mu_iso_dxy);
        //t_tree_->Branch("vld_mu_good_dxy", &ev.vld_mu_good_dxy);
        //t_tree_->Branch("vld_mu_hs_dxy", &ev.vld_mu_hs_dxy);
        //t_tree_->Branch("vld_mu_iso_hs_dxy", &ev.vld_mu_iso_hs_dxy);
        //t_tree_->Branch("vld_mu_good_hs_dxy", &ev.vld_mu_good_hs_dxy);
        //t_tree_->Branch("vld_mu_py8_dxy", &ev.vld_mu_py8_dxy);
        //t_tree_->Branch("vld_mu_iso_py8_dxy", &ev.vld_mu_iso_py8_dxy);
        //t_tree_->Branch("vld_mu_good_py8_dxy", &ev.vld_mu_good_py8_dxy);
        //t_tree_->Branch("vld_mu_others_dxy", &ev.vld_mu_others_dxy);
        //t_tree_->Branch("vld_mu_iso_others_dxy", &ev.vld_mu_iso_others_dxy);
        //t_tree_->Branch("vld_mu_good_others_dxy", &ev.vld_mu_good_others_dxy);
        //t_tree_->Branch("vld_mu_dz", &ev.vld_mu_dz);
        //t_tree_->Branch("vld_mu_iso_dz", &ev.vld_mu_iso_dz);
        //t_tree_->Branch("vld_mu_good_dz", &ev.vld_mu_good_dz);
        //t_tree_->Branch("vld_mu_hs_dz", &ev.vld_mu_hs_dz);
        //t_tree_->Branch("vld_mu_iso_hs_dz", &ev.vld_mu_iso_hs_dz);
        //t_tree_->Branch("vld_mu_good_hs_dz", &ev.vld_mu_good_hs_dz);
        //t_tree_->Branch("vld_mu_py8_dz", &ev.vld_mu_py8_dz);
        //t_tree_->Branch("vld_mu_iso_py8_dz", &ev.vld_mu_iso_py8_dz);
        //t_tree_->Branch("vld_mu_good_py8_dz", &ev.vld_mu_good_py8_dz);
        //t_tree_->Branch("vld_mu_others_dz", &ev.vld_mu_others_dz);
        //t_tree_->Branch("vld_mu_iso_others_dz", &ev.vld_mu_iso_others_dz);
        //t_tree_->Branch("vld_mu_good_others_dz", &ev.vld_mu_good_others_dz);
    }

    ////event header
    //t_tree_->Branch("Run",               &ev.run,        "Run/I");
    //t_tree_->Branch("Event",             &ev.event,      "Event/I");
    //t_tree_->Branch("Lumi",              &ev.lumi,       "Lumi/I");

    ////gen level event
    //t_tree_->Branch("Particle_size",  &ev.ngl,        "Particle_size/I");
    //t_tree_->Branch("PID",            ev.gl_pid,      "PID[Particle_size]/I");
    //t_tree_->Branch("Charge",         ev.gl_ch,       "Charge[Particle_size]/I");
    //t_tree_->Branch("Status",         ev.gl_st,       "Status[Particle_size]/I");
    //t_tree_->Branch("P",              ev.gl_p,        "P[Particle_size]/F");
    //t_tree_->Branch("Px",             ev.gl_px,       "Px[Particle_size]/F");
    //t_tree_->Branch("Py",             ev.gl_py,       "Py[Particle_size]/F");
    //t_tree_->Branch("Pz",             ev.gl_pz,       "Pz[Particle_size]/F");
    //t_tree_->Branch("E",              ev.gl_nrj,      "E[Particle_size]/F");
    //t_tree_->Branch("PT",             ev.gl_pt,       "PT[Particle_size]/F");
    //t_tree_->Branch("Eta",            ev.gl_eta,      "Eta[Particle_size]/F");
    //t_tree_->Branch("Phi",            ev.gl_phi,      "Phi[Particle_size]/F");
    //t_tree_->Branch("Mass",           ev.gl_mass,     "Mass[Particle_size]/F");
    //t_tree_->Branch("IsolationVar",   ev.gl_relIso,   "IsolationVar/F");

    //t_tree_->Branch("GenJet_size",     &ev.ngj,        "GenJet_size/I");
    //t_tree_->Branch("PT",              ev.gj_pt,       "PT[GenJet_size]/F");
    //t_tree_->Branch("Eta",             ev.gj_eta,      "Eta[GenJet_size]/F");
    //t_tree_->Branch("Phi",             ev.gj_phi,      "Phi[GenJet_size]/F");
    //t_tree_->Branch("Mass",            ev.gj_mass,     "Mass[GenJet_size]/F");

    ////reco level event
    //t_tree_->Branch("Vertex_size",    &ev.nvtx,       "Vertex_size/I");
    //t_tree_->Branch("SumPT2",         &ev.v_pt2,      "SumPT2[Vertex_size]/F");

    //t_tree_->Branch("ElectronLoose_size", &ev.nle,  "ElectronLoose_size/I");
    //t_tree_->Branch("Charge",       ev.le_ch,       "Charge[ElectronLoose_size]/I");
    //t_tree_->Branch("Particle",     ev.le_g,        "Particle[ElectronLoose_size]/I");
    //t_tree_->Branch("PT",           ev.le_pt,       "PT[ElectronLoose_size]/F");
    //t_tree_->Branch("Eta",          ev.le_eta,      "Eta[ElectronLoose_size]/F");
    //t_tree_->Branch("Phi",          ev.le_phi,      "Phi[ElectronLoose_size]/F");
    //t_tree_->Branch("Mass",         ev.le_mass,     "Mass[ElectronLoose_size]/F");
    //t_tree_->Branch("IsolationVar", ev.le_relIso,   "IsolationVar[ElectronLoose_size]/F");

    //t_tree_->Branch("ElectronTight_size", &ev.nte,  "ElectronTight_size/I");
    //t_tree_->Branch("Charge",       ev.te_ch,       "Charge[ElectronTight_size]/I");
    //t_tree_->Branch("Particle",     ev.te_g,        "Particle[ElectronTight_size]/I");
    //t_tree_->Branch("PT",           ev.te_pt,       "PT[ElectronTight_size]/F");
    //t_tree_->Branch("Eta",          ev.te_eta,      "Eta[ElectronTight_size]/F");
    //t_tree_->Branch("Phi",          ev.te_phi,      "Phi[ElectronTight_size]/F");
    //t_tree_->Branch("Mass",         ev.te_mass,     "Mass[ElectronTight_size]/F");
    //t_tree_->Branch("IsolationVar", ev.te_relIso,   "IsolationVar[ElectronTight_size]/F");

    //t_tree_->Branch("MuonLoose_size", &ev.nlm,      "MuonLoose_size/I");
    //t_tree_->Branch("Charge",       ev.lm_ch,       "Charge[MuonLoose_size]/I");
    //t_tree_->Branch("Particle",     ev.lm_g,        "Particle[MuonLoose_size]/I");
    //t_tree_->Branch("PT",           ev.lm_pt,       "PT[MuonLoose_size]/F");
    //t_tree_->Branch("Eta",          ev.lm_eta,      "Eta[MuonLoose_size]/F");
    //t_tree_->Branch("Phi",          ev.lm_phi,      "Phi[MuonLoose_size]/F");
    //t_tree_->Branch("Mass",         ev.lm_mass,     "Mass[MuonLoose_size]/F");
    //t_tree_->Branch("IsolationVar", ev.lm_relIso,   "IsolationVar[MuonLoose_size]/F");

    //t_tree_->Branch("MuonTight_size", &ev.ntm,      "MuonTight_size/I");
    //t_tree_->Branch("Charge",       ev.tm_ch,       "Charge[MuonTight_size]/I");
    //t_tree_->Branch("Particle",     ev.tm_g,        "Particle[MuonTight_size]/I");
    //t_tree_->Branch("PT",           ev.tm_pt,       "PT[MuonTight_size]/F");
    //t_tree_->Branch("Eta",          ev.tm_eta,      "Eta[MuonTight_size]/F");
    //t_tree_->Branch("Phi",          ev.tm_phi,      "Phi[MuonTight_size]/F");
    //t_tree_->Branch("Mass",         ev.tm_mass,     "Mass[MuonTight_size]/F");
    //t_tree_->Branch("IsolationVar", ev.tm_relIso,   "IsolationVar[MuonTight_size]/F");

    //t_tree_->Branch("JetPUPPI_size", &ev.nj,         "JetPUPPI_size/I");
    //t_tree_->Branch("ID",            ev.j_id,        "ID[JetPUPPI_size]/I");
    //t_tree_->Branch("GenJet",        ev.j_g,         "GenJet[JetPUPPI_size]/I");
    //t_tree_->Branch("PT",            ev.j_pt,        "PT[JetPUPPI_size]/F");
    //t_tree_->Branch("Eta",           ev.j_eta,       "Eta[JetPUPPI_size]/F");
    //t_tree_->Branch("Phi",           ev.j_phi,       "Phi[JetPUPPI_size]/F");
    //t_tree_->Branch("Mass",          ev.j_mass,      "Mass[JetPUPPI_size]/F");
    //t_tree_->Branch("MVAv2",         ev.j_mvav2,     "MVAv2[JetPUPPI_size]/I");
    //t_tree_->Branch("DeepCSV",       ev.j_deepcsv,   "DeepCSV[JetPUPPI_size]/I");
    //t_tree_->Branch("PartonFlavor",  ev.j_flav,      "PartonFlavor[JetPUPPI_size]/I");
    //t_tree_->Branch("HadronFlavor",  ev.j_hadflav,   "HadronFlavor[JetPUPPI_size]/I");
    //t_tree_->Branch("GenPartonPID",  ev.j_pid,       "GenPartonPID[JetPUPPI_size]/I");

    //t_tree_->Branch("PuppiMissingET_size", &ev.nmet,  "PuppiMissingET_size/I");
    //t_tree_->Branch("MET",            ev.met_pt,      "MET[PuppiMissingET_size]/F");
    //t_tree_->Branch("Phi",            ev.met_phi,     "Phi[PuppiMissingET_size]/F");
    //t_tree_->Branch("Eta",            ev.met_eta,     "Eta[PuppiMissingET_size]/F");
}


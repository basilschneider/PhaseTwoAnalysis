#include "PhaseTwoAnalysis/NTupler/interface/MiniEvent.h"

void createMiniEventTree(TTree* t_tree_, MiniEvent_t &ev)
{
  //event header
  t_tree_->Branch("Run",               &ev.run,        "Run/I");
  t_tree_->Branch("Event",             &ev.event,      "Event/I");
  t_tree_->Branch("Lumi",              &ev.lumi,       "Lumi/I");

  //gen level event
  t_tree_->Branch("Particle_size",  &ev.ngl,        "Particle_size/I");
  t_tree_->Branch("PID",            ev.gl_pid,      "PID[Particle_size]/I");
  t_tree_->Branch("Charge",         ev.gl_ch,       "Charge[Particle_size]/I");
  t_tree_->Branch("Status",         ev.gl_st,       "Status[Particle_size]/I");
  t_tree_->Branch("P",              ev.gl_p,        "P[Particle_size]/F");
  t_tree_->Branch("Px",             ev.gl_px,       "Px[Particle_size]/F");
  t_tree_->Branch("Py",             ev.gl_py,       "Py[Particle_size]/F");
  t_tree_->Branch("Pz",             ev.gl_pz,       "Pz[Particle_size]/F");
  t_tree_->Branch("E",              ev.gl_nrj,      "E[Particle_size]/F");
  t_tree_->Branch("PT",             ev.gl_pt,       "PT[Particle_size]/F");
  t_tree_->Branch("Eta",            ev.gl_eta,      "Eta[Particle_size]/F");
  t_tree_->Branch("Phi",            ev.gl_phi,      "Phi[Particle_size]/F");
  t_tree_->Branch("Mass",           ev.gl_mass,     "Mass[Particle_size]/F");
  t_tree_->Branch("IsolationVar",   ev.gl_relIso,   "IsolationVar/F");

  t_tree_->Branch("GenJet_size",     &ev.ngj,        "GenJet_size/I");
  t_tree_->Branch("PT",              ev.gj_pt,       "PT[GenJet_size]/F");
  t_tree_->Branch("Eta",             ev.gj_eta,      "Eta[GenJet_size]/F");
  t_tree_->Branch("Phi",             ev.gj_phi,      "Phi[GenJet_size]/F");
  t_tree_->Branch("Mass",            ev.gj_mass,     "Mass[GenJet_size]/F");

  //reco level event
  t_tree_->Branch("Vertex_size",    &ev.nvtx,       "Vertex_size/I");
  t_tree_->Branch("SumPT2",         &ev.v_pt2,      "SumPT2[Vertex_size]/F");

  t_tree_->Branch("ElectronLoose_size", &ev.nle,  "ElectronLoose_size/I");
  t_tree_->Branch("Charge",       ev.le_ch,       "Charge[ElectronLoose_size]/I");
  t_tree_->Branch("Particle",     ev.le_g,        "Particle[ElectronLoose_size]/I");
  t_tree_->Branch("PT",           ev.le_pt,       "PT[ElectronLoose_size]/F");
  t_tree_->Branch("Eta",          ev.le_eta,      "Eta[ElectronLoose_size]/F");
  t_tree_->Branch("Phi",          ev.le_phi,      "Phi[ElectronLoose_size]/F");
  t_tree_->Branch("Mass",         ev.le_mass,     "Mass[ElectronLoose_size]/F");
  t_tree_->Branch("IsolationVar", ev.le_relIso,   "IsolationVar[ElectronLoose_size]/F");

  t_tree_->Branch("ElectronTight_size", &ev.nte,  "ElectronTight_size/I");
  t_tree_->Branch("Charge",       ev.te_ch,       "Charge[ElectronTight_size]/I");
  t_tree_->Branch("Particle",     ev.te_g,        "Particle[ElectronTight_size]/I");
  t_tree_->Branch("PT",           ev.te_pt,       "PT[ElectronTight_size]/F");
  t_tree_->Branch("Eta",          ev.te_eta,      "Eta[ElectronTight_size]/F");
  t_tree_->Branch("Phi",          ev.te_phi,      "Phi[ElectronTight_size]/F");
  t_tree_->Branch("Mass",         ev.te_mass,     "Mass[ElectronTight_size]/F");
  t_tree_->Branch("IsolationVar", ev.te_relIso,   "IsolationVar[ElectronTight_size]/F");

  t_tree_->Branch("MuonLoose_size", &ev.nlm,      "MuonLoose_size/I");
  t_tree_->Branch("Charge",       ev.lm_ch,       "Charge[MuonLoose_size]/I");
  t_tree_->Branch("Particle",     ev.lm_g,        "Particle[MuonLoose_size]/I");
  t_tree_->Branch("PT",           ev.lm_pt,       "PT[MuonLoose_size]/F");
  t_tree_->Branch("Eta",          ev.lm_eta,      "Eta[MuonLoose_size]/F");
  t_tree_->Branch("Phi",          ev.lm_phi,      "Phi[MuonLoose_size]/F");
  t_tree_->Branch("Mass",         ev.lm_mass,     "Mass[MuonLoose_size]/F");
  t_tree_->Branch("IsolationVar", ev.lm_relIso,   "IsolationVar[MuonLoose_size]/F");

  t_tree_->Branch("MuonTight_size", &ev.ntm,      "MuonTight_size/I");
  t_tree_->Branch("Charge",       ev.tm_ch,       "Charge[MuonTight_size]/I");
  t_tree_->Branch("Particle",     ev.tm_g,        "Particle[MuonTight_size]/I");
  t_tree_->Branch("PT",           ev.tm_pt,       "PT[MuonTight_size]/F");
  t_tree_->Branch("Eta",          ev.tm_eta,      "Eta[MuonTight_size]/F");
  t_tree_->Branch("Phi",          ev.tm_phi,      "Phi[MuonTight_size]/F");
  t_tree_->Branch("Mass",         ev.tm_mass,     "Mass[MuonTight_size]/F");
  t_tree_->Branch("IsolationVar", ev.tm_relIso,   "IsolationVar[MuonTight_size]/F");

  t_tree_->Branch("JetPUPPI_size", &ev.nj,         "JetPUPPI_size/I");
  t_tree_->Branch("ID",            ev.j_id,        "ID[JetPUPPI_size]/I");
  t_tree_->Branch("GenJet",        ev.j_g,         "GenJet[JetPUPPI_size]/I");
  t_tree_->Branch("PT",            ev.j_pt,        "PT[JetPUPPI_size]/F");
  t_tree_->Branch("Eta",           ev.j_eta,       "Eta[JetPUPPI_size]/F");
  t_tree_->Branch("Phi",           ev.j_phi,       "Phi[JetPUPPI_size]/F");
  t_tree_->Branch("Mass",          ev.j_mass,      "Mass[JetPUPPI_size]/F");
  t_tree_->Branch("MVAv2",         ev.j_mvav2,     "MVAv2[JetPUPPI_size]/I");
  t_tree_->Branch("DeepCSV",       ev.j_deepcsv,   "DeepCSV[JetPUPPI_size]/I");
  t_tree_->Branch("PartonFlavor",  ev.j_flav,      "PartonFlavor[JetPUPPI_size]/I");
  t_tree_->Branch("HadronFlavor",  ev.j_hadflav,   "HadronFlavor[JetPUPPI_size]/I");
  t_tree_->Branch("GenPartonPID",  ev.j_pid,       "GenPartonPID[JetPUPPI_size]/I");

  t_tree_->Branch("PuppiMissingET_size", &ev.nmet,  "PuppiMissingET_size/I");
  t_tree_->Branch("MET",            ev.met_pt,      "MET[PuppiMissingET_size]/F");
  t_tree_->Branch("Phi",            ev.met_phi,     "Phi[PuppiMissingET_size]/F");
  t_tree_->Branch("Eta",            ev.met_eta,     "Eta[PuppiMissingET_size]/F");
}


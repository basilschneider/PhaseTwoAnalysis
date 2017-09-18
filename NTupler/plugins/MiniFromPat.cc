// -*- C++ -*-
//
// Package:    PhaseTwoAnalysis/NTupler
// Class:      MiniFromPat
//
/**\class MiniFromPat MiniFromPat.cc PhaseTwoAnalysis/NTupler/plugins/MiniFromPat.cc

Description: produces flat ntuples from PAT collections
- storing gen, reco, and pf leptons with pT > 10 GeV and |eta| < 3
- storing gen and reco jets with pT > 20 GeV and |eta| < 5

Implementation:
- muon isolation comes from https://twiki.cern.ch/twiki/bin/viewauth/CMS/Phase2MuonBarrelRecipes#Muon_isolation
- muon ID comes from https://twiki.cern.ch/twiki/bin/viewauth/CMS/Phase2MuonBarrelRecipes#Muon_identification
- electron isolation might need to be refined
- electron ID comes from https://indico.cern.ch/event/623893/contributions/2531742/attachments/1436144/2208665/UPSG_EGM_Workshop_Mar29.pdf
/!\ no ID is implemented for forward electrons as:
- PFClusterProducer does not run on miniAOD
- jurassic isolation needs tracks
- PF jet ID comes from Run-2 https://github.com/cms-sw/cmssw/blob/CMSSW_9_1_1_patch1/PhysicsTools/SelectorUtils/interface/PFJetIDSelectionFunctor.h
- no JEC applied
- b-tagging WPs come from https://twiki.cern.ch/twiki/bin/viewauth/CMS/Phase2MuonBarrelRecipes#B_tagging
*/

//
// Original Author:  Elvire Bouvier
//         Created:  Tue, 20 Jun 2017 11:27:12 GMT
//
//


// system include files
#include <memory>
#include <cmath>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"//
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "PhysicsTools/SelectorUtils/interface/PFJetIDSelectionFunctor.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/JetReco/interface/GenJet.h"

#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexFitter.h"
#include "SimTracker/Records/interface/TrackAssociatorRecord.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "Geometry/GEMGeometry/interface/ME0EtaPartitionSpecs.h"
#include "Geometry/GEMGeometry/interface/ME0Geometry.h"

#include "FWCore/Framework/interface/ESHandle.h"

#include "RecoVertex/KinematicFitPrimitives/interface/ParticleMass.h"
#include <RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h>
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicVertex.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "PhysicsTools/CandUtils/interface/AddFourMomenta.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "PhaseTwoAnalysis/NTupler/interface/MiniEvent.h"

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "Math/GenVector/VectorUtil.h"

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class MiniFromPat : public edm::one::EDAnalyzer<edm::one::SharedResources, edm::one::WatchRuns>  {
    public:
        explicit MiniFromPat(const edm::ParameterSet&);
        ~MiniFromPat();

        static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

        enum ElectronMatchType {UNMATCHED = 0,
            TRUE_PROMPT_ELECTRON,
            TRUE_ELECTRON_FROM_TAU,
            TRUE_NON_PROMPT_ELECTRON};


    private:
        virtual void beginJob() override;
        virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
        void genAnalysis(const edm::Event& iEvent, const edm::EventSetup& iSetup);
        void recoAnalysis(const edm::Event& iEvent, const edm::EventSetup& iSetup);
        virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
        virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
        virtual void endJob() override;

        bool isLooseElec(const pat::Electron & patEl, edm::Handle<reco::ConversionCollection> conversions, const reco::BeamSpot beamspot);
        bool isMediumElec(const pat::Electron & patEl, edm::Handle<reco::ConversionCollection> conversions, const reco::BeamSpot beamspot);
        bool isTightElec(const pat::Electron & patEl, edm::Handle<reco::ConversionCollection> conversions, const reco::BeamSpot beamspot);
        bool isGoodElecSOS(const pat::Electron & patEl, edm::Handle<reco::ConversionCollection> conversions, const reco::BeamSpot beamspot, edm::Handle<std::vector<reco::Vertex>> vertices);
        bool isGoodMuonSOS(const pat::Muon & patMu, edm::Handle<std::vector<reco::Vertex>> vertices, int prVtx);
        bool isGoodJetSOS(const pat::Jet & patJet);
        bool isGoodElecTruthSOS(const pat::PackedGenParticle truthEl, const std::vector<size_t> jGenJets, const edm::Handle<std::vector<reco::GenJet>> genJets);
        bool isGoodMuonTruthSOS(const pat::PackedGenParticle truthMu, const std::vector<size_t> jGenJets, const edm::Handle<std::vector<reco::GenJet>> genJets);
        bool isME0MuonSel(reco::Muon, double pullXCut, double dXCut, double pullYCut, double dYCut, double dPhi);
        bool isME0MuonSelNew(reco::Muon, double, double, double);

        // ----------member data ---------------------------
        edm::Service<TFileService> fs_;

        unsigned int pileup_;
        edm::EDGetTokenT<std::vector<reco::Vertex>> verticesToken_;
        edm::EDGetTokenT<std::vector<pat::Electron>> elecsToken_;
        edm::EDGetTokenT<reco::BeamSpot> bsToken_;
        edm::EDGetTokenT<std::vector<reco::Conversion>> convToken_;
        edm::EDGetTokenT<std::vector<pat::Muon>> muonsToken_;
        edm::EDGetTokenT<std::vector<pat::Jet>> jetsToken_;
        PFJetIDSelectionFunctor jetIDLoose_;
        PFJetIDSelectionFunctor jetIDTight_;
        edm::EDGetTokenT<std::vector<pat::MET>> metsToken_;
        edm::EDGetTokenT<std::vector<reco::GenJet>> genJetsToken_;
        edm::EDGetTokenT<std::vector<pat::PackedGenParticle>> genPartsToken_;
        const ME0Geometry* ME0Geometry_;
        double mvaThres_[3];
        double deepThres_[3];

        //TTree *t_event_, *t_genParts_, *t_vertices_, *t_genJets_, *t_looseElecs_, *t_tightElecs_, *t_looseMuons_, *t_tightMuons_, *t_puppiJets_, *t_puppiMET_;
        TTree *t_tree_;

        MiniEvent_t ev_;
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
MiniFromPat::MiniFromPat(const edm::ParameterSet& iConfig):
    pileup_(iConfig.getParameter<unsigned int>("pileup")),
    verticesToken_(consumes<std::vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("vertices"))),
    elecsToken_(consumes<std::vector<pat::Electron>>(iConfig.getParameter<edm::InputTag>("electrons"))),
    bsToken_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("beamspot"))),
    convToken_(consumes<std::vector<reco::Conversion>>(iConfig.getParameter<edm::InputTag>("conversions"))),
    muonsToken_(consumes<std::vector<pat::Muon>>(iConfig.getParameter<edm::InputTag>("muons"))),
    jetsToken_(consumes<std::vector<pat::Jet>>(iConfig.getParameter<edm::InputTag>("jets"))),
    jetIDLoose_(PFJetIDSelectionFunctor::FIRSTDATA, PFJetIDSelectionFunctor::LOOSE),
    jetIDTight_(PFJetIDSelectionFunctor::FIRSTDATA, PFJetIDSelectionFunctor::TIGHT),
    metsToken_(consumes<std::vector<pat::MET>>(iConfig.getParameter<edm::InputTag>("mets"))),
    genJetsToken_(consumes<std::vector<reco::GenJet>>(iConfig.getParameter<edm::InputTag>("genJets"))),
    genPartsToken_(consumes<std::vector<pat::PackedGenParticle>>(iConfig.getParameter<edm::InputTag>("genParts")))
{
    //now do what ever initialization is needed
    if (pileup_ == 0) {
        mvaThres_[0] = -0.694;
        mvaThres_[1] = 0.128;
        mvaThres_[2] = 0.822;
        deepThres_[0] = 0.131;
        deepThres_[1] = 0.432;
        deepThres_[2] = 0.741;
    } else if (pileup_ == 140) {
        mvaThres_[0] = -0.654;
        mvaThres_[1] = 0.214;
        mvaThres_[2] = 0.864;
        deepThres_[0] = 0.159;
        deepThres_[1] = 0.507;
        deepThres_[2] = 0.799;
    } else if (pileup_ == 200) {
        mvaThres_[0] = -0.642;
        mvaThres_[1] = 0.236;
        mvaThres_[2] = 0.878;
        deepThres_[0] = 0.170;
        deepThres_[1] = 0.527;
        deepThres_[2] = 0.821;
    } else {
        mvaThres_[0] = -1.;
        mvaThres_[1] = -1.;
        mvaThres_[2] = -1.;
        deepThres_[0] = 0.;
        deepThres_[1] = 0.;
        deepThres_[2] = 0.;
    }

    //t_event_      = fs_->make<TTree>("Event","Event");
    //t_genParts_   = fs_->make<TTree>("Particle","Particle");
    //t_vertices_   = fs_->make<TTree>("Vertex","Vertex");
    //t_genJets_    = fs_->make<TTree>("GenJet","GenJet");
    //t_looseElecs_ = fs_->make<TTree>("ElectronLoose","ElectronLoose");
    //t_tightElecs_ = fs_->make<TTree>("ElectronTight","ElectronTight");
    //t_looseMuons_ = fs_->make<TTree>("MuonLoose","MuonLoose");
    //t_tightMuons_ = fs_->make<TTree>("MuonTight","MuonTight");
    //t_puppiJets_  = fs_->make<TTree>("JetPUPPI","JetPUPPI");
    //t_puppiMET_   = fs_->make<TTree>("PuppiMissingET","PuppiMissingET");
    t_tree_ = fs_->make<TTree>("Delphes", "Delphes");
    //createMiniEventTree(t_event_, t_genParts_, t_vertices_, t_genJets_, t_looseElecs_, t_tightElecs_, t_looseMuons_, t_tightMuons_, t_puppiJets_, t_puppiMET_, ev_);
    createMiniEventTree(t_tree_, ev_);

}


MiniFromPat::~MiniFromPat()
{

    // do anything here that needs to be done at desctruction time
    // (e.g. close files, deallocate resources etc.)

    if (ev_.fill_rle){
        ev_.rle_el_num->Write();
        ev_.rle_el_den->Write();
        ev_.rle_mu_num->Write();
        ev_.rle_mu_den->Write();
        TH2D* rle_el = (TH2D*)ev_.rle_el_num->Clone();
        TH2D* rle_mu = (TH2D*)ev_.rle_mu_num->Clone();
        rle_el->SetNameTitle("rle_el", "rle_el");
        rle_mu->SetNameTitle("rle_mu", "rle_mu");
        rle_el->Divide(ev_.rle_el_den);
        rle_mu->Divide(ev_.rle_mu_den);
        rle_el->Write();
        rle_mu->Write();
    }
    if (ev_.fill_vld){
        ev_.vld_el_pt_iso_abs->Write();
        ev_.vld_el_tight_pt_iso_abs->Write();
        ev_.vld_el_pt_iso_rel->Write();
        ev_.vld_el_tight_pt_iso_rel->Write();
        ev_.vld_el_pt_dxy->Write();
        ev_.vld_el_tight_pt_dxy->Write();
        ev_.vld_el_pt_dz->Write();
        ev_.vld_el_tight_pt_dz->Write();
        ev_.vld_mu_pt_iso_abs->Write();
        ev_.vld_mu_tight_pt_iso_abs->Write();
        ev_.vld_mu_pt_iso_rel->Write();
        ev_.vld_mu_tight_pt_iso_rel->Write();
        ev_.vld_mu_pt_dxy->Write();
        ev_.vld_mu_tight_pt_dxy->Write();
        ev_.vld_mu_pt_dz->Write();
        ev_.vld_mu_tight_pt_dz->Write();
    }
}


//
// member functions
//

// ------------ method to fill gen level pat -------------
    void
MiniFromPat::genAnalysis(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    using namespace edm;

    Handle<std::vector<pat::PackedGenParticle>> genParts;
    iEvent.getByToken(genPartsToken_, genParts);
    Handle<std::vector<reco::GenJet>> genJets;
    iEvent.getByToken(genJetsToken_, genJets);
    Handle<std::vector<pat::Electron>> elecs;
    iEvent.getByToken(elecsToken_, elecs);
    Handle<reco::ConversionCollection> conversions;
    iEvent.getByToken(convToken_, conversions);
    Handle<reco::BeamSpot> bsHandle;
    iEvent.getByToken(bsToken_, bsHandle);
    const reco::BeamSpot &beamspot = *bsHandle.product();
    Handle<std::vector<pat::Muon>> muons;
    iEvent.getByToken(muonsToken_, muons);
    Handle<std::vector<reco::Vertex>> vertices;
    iEvent.getByToken(verticesToken_, vertices);

    // Vertices
    int prVtx = -1;
    ev_.nvtx = 0;
    for (size_t i=0; i<vertices->size(); ++i){
        if (vertices->at(i).isFake()){ continue; }
        if (vertices->at(i).ndof() <= 4){ continue; }
        if (prVtx < 0){ prVtx = i; }
        ev_.v_pt2[ev_.nvtx] = vertices->at(i).p4().pt();
        ev_.nvtx++;
    }
    if (prVtx < 0) return;

    // Truth Jets
    std::vector<size_t> jGenJets;
    for (size_t i = 0; i < genJets->size(); i++) {
        if (genJets->at(i).pt() < 20.){ continue; }
        if (fabs(genJets->at(i).eta()) > 5){ continue; }

        bool overlaps = false;
        for (size_t j = 0; j < genParts->size(); j++) {
            if (abs(genParts->at(j).pdgId()) != 11 && abs(genParts->at(j).pdgId()) != 13) continue;
            if (fabs(genJets->at(i).pt()-genParts->at(j).pt()) < .01*genParts->at(j).pt() && ROOT::Math::VectorUtil::DeltaR(genParts->at(j).p4(),genJets->at(i).p4()) < .01) {
                overlaps = true;
                break;
            }
        }
        if (overlaps){ continue; }
        jGenJets.push_back(i);

        if (ev_.jet1_pt_truth.size() == 0){
            ev_.jet1_pt_truth.push_back(genJets->at(i).pt());
            ev_.jet1_eta_truth.push_back(genJets->at(i).eta());
            ev_.jet1_phi_truth.push_back(genJets->at(i).phi());
            ev_.jet1_mass_truth.push_back(genJets->at(i).mass());
        }
    }

    // Count particles with status 23 for real lepton efficiency histograms
    // This is not working, as status 23 particles are not stored. Only select
    // the hardest lepton (for semileptonic ttbar); if you want to store more
    // than one lepton, the logic with maxLepPt needs to be rewritten
    unsigned int nGenElStatus23 = 0;
    unsigned int nGenMuStatus23 = 0;
    double maxLepPt = -1.;
    if (ev_.fill_rle){
        for (size_t i=0; i<genParts->size(); ++i){
            if (genParts->at(i).pt() > 2 && genParts->at(i).pt() > maxLepPt){
                if (fabs(genParts->at(i).pdgId()) == 11){
                    if (!isGoodElecTruthSOS(genParts->at(i), jGenJets, genJets)){ continue; }
                    nGenElStatus23 = 1;
                    nGenMuStatus23 = 0;
                    maxLepPt = genParts->at(i).pt();
                }else if (fabs(genParts->at(i).pdgId()) == 13){
                    if (!isGoodMuonTruthSOS(genParts->at(i), jGenJets, genJets)){ continue; }
                    nGenElStatus23 = 0;
                    nGenMuStatus23 = 1;
                    maxLepPt = genParts->at(i).pt();
                }
            }
            //if (genParts->at(i).status() == 23){
            //    if (fabs(genParts->at(i).pdgId()) == 11){
            //        nGenElStatus23++;
            //    }else if (fabs(genParts->at(i).pdgId()) == 13){
            //        nGenMuStatus23++;
            //    }
            //}
        }
    }

    // Truth leptons
    for (size_t i = 0; i < genParts->size(); i++) {
        if (genParts->at(i).status() != 1){ continue; }

        if (ev_.fill_rle){
            // Fill real lepton efficiency histograms as many times as there are
            // Status 23 particles, but use Status 1 particles for this
            if (fabs(genParts->at(i).pdgId()) == 11){
                if (!nGenElStatus23){ continue; }
                if (!isGoodElecTruthSOS(genParts->at(i), jGenJets, genJets)){ continue; }
                if (genParts->at(i).pt() != maxLepPt){ continue; }

                // Fill denominator histogram
                ev_.rle_el_den->Fill(genParts->at(i).pt(), fabs(genParts->at(i).eta()));
                nGenElStatus23--;

                // Check if we can match that particle
                for (size_t j=0; j<elecs->size(); ++j){
                    if (!isGoodElecSOS(elecs->at(j), conversions, beamspot, vertices)){ continue; }
                    // If we make it here, this is a proper reco electron
                    // Fill numerator histogram if it can be matched
                    if (fabs(genParts->at(i).pt() - elecs->at(j).pt()) < ev_.truth_match_diff_pt \
                            && fabs(genParts->at(i).eta() - elecs->at(j).eta()) < ev_.truth_match_diff_eta){
                        ev_.rle_el_num->Fill(genParts->at(i).pt(), fabs(genParts->at(i).eta()));
                        break;
                    }
                }

            }else if (fabs(genParts->at(i).pdgId()) == 13){
                if (!nGenMuStatus23){ continue; }
                if (!isGoodMuonTruthSOS(genParts->at(i), jGenJets, genJets)){ continue; }
                if (genParts->at(i).pt() != maxLepPt){ continue; }

                // Fill denominator histogram
                ev_.rle_mu_den->Fill(genParts->at(i).pt(), fabs(genParts->at(i).eta()));
                nGenMuStatus23--;

                // Check if we can match that particle
                for (size_t j=0; j<muons->size(); ++j){
                    if (!isGoodMuonSOS(muons->at(j), vertices, prVtx)){ continue; }
                    // If we make it here, this is a proper reco muon
                    // Fill numerator histogram if it can be matched
                    if (fabs(genParts->at(i).pt() - muons->at(j).pt()) < ev_.truth_match_diff_pt \
                            && fabs(genParts->at(i).eta() - muons->at(j).eta()) < ev_.truth_match_diff_eta){
                        ev_.rle_mu_num->Fill(genParts->at(i).pt(), fabs(genParts->at(i).eta()));
                        break;
                    }
                }
            }
        }

        // FIXME: truth leptons are not ordered
        //if (abs(genParts->at(i).pdgId()) == 11){
        //    // Only select good truth electrons
        //    if (!isGoodElecTruthSOS(genParts->at(i), jGenJets, genJets)){ continue; }
        //    if (ev_.el1_pt_truth.size() == 0){
        //        ev_.el1_pt_truth.push_back(genParts->at(i).pt());
        //        ev_.el1_eta_truth.push_back(genParts->at(i).eta());
        //        ev_.el1_phi_truth.push_back(genParts->at(i).phi());
        //        ev_.el1_q_truth.push_back(genParts->at(i).charge());
        //    }else if (ev_.el2_pt_truth.size() == 0){
        //        ev_.el2_pt_truth.push_back(genParts->at(i).pt());
        //        ev_.el2_eta_truth.push_back(genParts->at(i).eta());
        //        ev_.el2_phi_truth.push_back(genParts->at(i).phi());
        //        ev_.el2_q_truth.push_back(genParts->at(i).charge());
        //    }
        //}else if (abs(genParts->at(i).pdgId()) == 13){
        //    if (!isGoodMuonTruthSOS(genParts->at(i), jGenJets, genJets)){ continue; }
        //    if (ev_.mu1_pt_truth.size() == 0){
        //        ev_.mu1_pt_truth.push_back(genParts->at(i).pt());
        //        ev_.mu1_eta_truth.push_back(genParts->at(i).eta());
        //        ev_.mu1_phi_truth.push_back(genParts->at(i).phi());
        //        ev_.mu1_q_truth.push_back(genParts->at(i).charge());
        //    }else if (ev_.mu2_pt_truth.size() == 0){
        //        ev_.mu2_pt_truth.push_back(genParts->at(i).pt());
        //        ev_.mu2_eta_truth.push_back(genParts->at(i).eta());
        //        ev_.mu2_phi_truth.push_back(genParts->at(i).phi());
        //        ev_.mu2_q_truth.push_back(genParts->at(i).charge());
        //    }
        //}
    }

    // FIXME: Once the truth lepton ordering above is fixed, you can put them
    // back in
    //// Fill truth leptons
    //// Put pT and eta into vector of vector for sorting
    //std::vector<std::vector<double>> lepvec_truth;
    //if (ev_.el1_pt_truth.size() != 0){
    //    lepvec_truth.push_back({ev_.el1_pt_truth.at(0), ev_.el1_eta_truth.at(0), ev_.el1_phi_truth.at(0), ev_.mass_el});
    //}
    //if (ev_.el2_pt_truth.size() != 0){
    //    lepvec_truth.push_back({ev_.el2_pt_truth.at(0), ev_.el2_eta_truth.at(0), ev_.el2_phi_truth.at(0), ev_.mass_el});
    //}
    //if (ev_.mu1_pt_truth.size() != 0){
    //    lepvec_truth.push_back({ev_.mu1_pt_truth.at(0), ev_.mu1_eta_truth.at(0), ev_.mu1_phi_truth.at(0), ev_.mass_mu});
    //}
    //if (ev_.mu2_pt_truth.size() != 0){
    //    lepvec_truth.push_back({ev_.mu2_pt_truth.at(0), ev_.mu2_eta_truth.at(0), ev_.mu2_phi_truth.at(0), ev_.mass_mu});
    //}
    //// By definition, this sorts by the first element of the vector (in this case pT)
    //if (lepvec_truth.size() > 1){
    //    std::sort(begin(lepvec_truth), end(lepvec_truth));
    //    std::reverse(begin(lepvec_truth), end(lepvec_truth));
    //}
    //if (lepvec_truth.size() >= 1){
    //    ev_.lep1_pt_truth.push_back(lepvec_truth[0][0]);
    //    ev_.lep1_eta_truth.push_back(lepvec_truth[0][1]);
    //    ev_.lep1_phi_truth.push_back(lepvec_truth[0][2]);
    //    ev_.lep1_mass_truth.push_back(lepvec_truth[0][3]);
    //}
    //if (lepvec_truth.size() >= 2){
    //    ev_.lep2_pt_truth.push_back(lepvec_truth[1][0]);
    //    ev_.lep2_eta_truth.push_back(lepvec_truth[1][1]);
    //    ev_.lep2_phi_truth.push_back(lepvec_truth[1][2]);
    //    ev_.lep2_mass_truth.push_back(lepvec_truth[1][3]);
    //}
    //lepvec_truth.clear();


    //// Jets
    //std::vector<size_t> jGenJets;
    //ev_.ngj = 0;
    //for (size_t i = 0; i < genJets->size(); i++) {
    //  if (genJets->at(i).pt() < 20.) continue;
    //  if (fabs(genJets->at(i).eta()) > 5) continue;

    //  bool overlaps = false;
    //  for (size_t j = 0; j < genParts->size(); j++) {
    //    if (abs(genParts->at(j).pdgId()) != 11 && abs(genParts->at(j).pdgId()) != 13) continue;
    //    if (fabs(genJets->at(i).pt()-genParts->at(j).pt()) < 0.01*genParts->at(j).pt() && ROOT::Math::VectorUtil::DeltaR(genParts->at(j).p4(),genJets->at(i).p4()) < 0.01) {
    //      overlaps = true;
    //      break;
    //    }
    //  }
    //  if (overlaps) continue;
    //  jGenJets.push_back(i);

    //  ev_.gj_pt[ev_.ngj]   = genJets->at(i).pt();
    //  ev_.gj_phi[ev_.ngj]  = genJets->at(i).phi();
    //  ev_.gj_eta[ev_.ngj]  = genJets->at(i).eta();
    //  ev_.gj_mass[ev_.ngj] = genJets->at(i).mass();
    //  ev_.ngj++;
    //}

    //// Leptons
    //ev_.ngl = 0;
    //for (size_t i = 0; i < genParts->size(); i++) {
    //  if (abs(genParts->at(i).pdgId()) != 11 && abs(genParts->at(i).pdgId()) != 13) continue;
    //  if (genParts->at(i).pt() < 10.) continue;
    //  if (fabs(genParts->at(i).eta()) > 3.) continue;
    //  double genIso = 0.;
    //  for (size_t j = 0; j < jGenJets.size(); j++) {
    //    if (ROOT::Math::VectorUtil::DeltaR(genParts->at(i).p4(),genJets->at(jGenJets[j]).p4()) > 0.7) continue;
    //    std::vector<const reco::Candidate *> jconst = genJets->at(jGenJets[j]).getJetConstituentsQuick();
    //    for (size_t k = 0; k < jconst.size(); k++) {
    //      double deltaR = ROOT::Math::VectorUtil::DeltaR(genParts->at(i).p4(),jconst[k]->p4());
    //      if (deltaR < 0.01 || deltaR > 0.4) continue;
    //      genIso = genIso + jconst[k]->pt();
    //    }
    //  }
    //  genIso = genIso / genParts->at(i).pt();
    //  ev_.gl_pid[ev_.ngl]    = genParts->at(i).pdgId();
    //  ev_.gl_ch[ev_.ngl]     = genParts->at(i).charge();
    //  ev_.gl_st[ev_.ngl]     = genParts->at(i).status();
    //  ev_.gl_p[ev_.ngl]      = genParts->at(i).p();
    //  ev_.gl_px[ev_.ngl]     = genParts->at(i).px();
    //  ev_.gl_py[ev_.ngl]     = genParts->at(i).py();
    //  ev_.gl_pz[ev_.ngl]     = genParts->at(i).pz();
    //  ev_.gl_nrj[ev_.ngl]    = genParts->at(i).energy();
    //  ev_.gl_pt[ev_.ngl]     = genParts->at(i).pt();
    //  ev_.gl_phi[ev_.ngl]    = genParts->at(i).phi();
    //  ev_.gl_eta[ev_.ngl]    = genParts->at(i).eta();
    //  ev_.gl_mass[ev_.ngl]   = genParts->at(i).mass();
    //  ev_.gl_relIso[ev_.ngl] = genIso;
    //  ev_.ngl++;
    //}
}

// ------------ method to fill reco level pat -------------
    void
MiniFromPat::recoAnalysis(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    using namespace edm;

    Handle<std::vector<reco::Vertex>> vertices;
    iEvent.getByToken(verticesToken_, vertices);
    Handle<std::vector<pat::Electron>> elecs;
    iEvent.getByToken(elecsToken_, elecs);
    Handle<reco::ConversionCollection> conversions;
    iEvent.getByToken(convToken_, conversions);
    Handle<reco::BeamSpot> bsHandle;
    iEvent.getByToken(bsToken_, bsHandle);
    const reco::BeamSpot &beamspot = *bsHandle.product();
    Handle<std::vector<pat::Muon>> muons;
    iEvent.getByToken(muonsToken_, muons);
    Handle<std::vector<pat::MET>> mets;
    iEvent.getByToken(metsToken_, mets);
    Handle<std::vector<pat::Jet>> jets;
    iEvent.getByToken(jetsToken_, jets);

    // Vertices
    int prVtx = -1;
    ev_.nvtx = 0;
    for (size_t i=0; i<vertices->size(); ++i){
        if (vertices->at(i).isFake()){ continue; }
        if (vertices->at(i).ndof() <= 4){ continue; }
        if (prVtx < 0){ prVtx = i; }
        ev_.v_pt2[ev_.nvtx] = vertices->at(i).p4().pt();
        ev_.nvtx++;
    }
    if (prVtx < 0) return;

    // Validation
    if (ev_.fill_vld){

        // Electrons
        for (size_t i=0; i<elecs->size(); ++i) {
            double pt = elecs->at(i).pt();
            double iso_abs = elecs->at(i).pfIsolationVariables().sumChargedHadronPt;
            double iso_rel = iso_abs/pt;
            const reco::Vertex &pv = vertices->front();
            double dxy = std::abs(elecs->at(i).gsfTrack()->dxy(pv.position()));
            double dz = std::abs(elecs->at(i).gsfTrack()->dz(pv.position()));
            ev_.vld_el_pt.push_back(pt);
            ev_.vld_el_is_tight.push_back(isTightElec(elecs->at(i), conversions, beamspot));
            ev_.vld_el_iso_abs.push_back(iso_abs);
            ev_.vld_el_iso_rel.push_back(iso_rel);
            ev_.vld_el_dxy.push_back(dxy);
            ev_.vld_el_dz.push_back(dz);
            ev_.vld_el_pt_iso_abs->Fill(pt, iso_abs);
            ev_.vld_el_pt_iso_rel->Fill(pt, iso_rel);
            ev_.vld_el_pt_dxy->Fill(pt, dxy);
            ev_.vld_el_pt_dz->Fill(pt, dz);
            if (isGoodElecSOS(elecs->at(i), conversions, beamspot, vertices)){
                ev_.vld_el_tight_pt.push_back(pt);
                ev_.vld_el_tight_iso_abs.push_back(iso_abs);
                ev_.vld_el_tight_iso_rel.push_back(iso_rel);
                ev_.vld_el_tight_dxy.push_back(dxy);
                ev_.vld_el_tight_dz.push_back(dz);
                ev_.vld_el_tight_pt_iso_abs->Fill(pt, iso_abs);
                ev_.vld_el_tight_pt_iso_rel->Fill(pt, iso_rel);
                ev_.vld_el_tight_pt_dxy->Fill(pt, dxy);
                ev_.vld_el_tight_pt_dz->Fill(pt, dz);
            }
        }

        // Muons
        for (size_t i=0; i<muons->size(); ++i){
            double pt = muons->at(i).pt();
            double iso_abs = muons->at(i).puppiChargedHadronIso() + muons->at(i).puppiNeutralHadronIso() + muons->at(i).puppiPhotonIso();
            double iso_rel = iso_abs/pt;
            double dxy = std::abs(muons->at(i).muonBestTrack()->dxy(vertices->at(prVtx).position()));
            double dz = std::abs(muons->at(i).muonBestTrack()->dz(vertices->at(prVtx).position()));
            ev_.vld_mu_pt.push_back(pt);
            ev_.vld_mu_is_tight.push_back(muon::isTightMuon(muons->at(i), vertices->at(prVtx)));
            ev_.vld_mu_iso_abs.push_back(iso_abs);
            ev_.vld_mu_iso_rel.push_back(iso_rel);
            ev_.vld_mu_dxy.push_back(dxy);
            ev_.vld_mu_dz.push_back(dz);
            ev_.vld_mu_pt_iso_abs->Fill(pt, iso_abs);
            ev_.vld_mu_pt_iso_rel->Fill(pt, iso_rel);
            ev_.vld_mu_pt_dxy->Fill(pt, dxy);
            ev_.vld_mu_pt_dz->Fill(pt, dz);
            if (isGoodMuonSOS(muons->at(i), vertices, prVtx)){
                ev_.vld_mu_tight_pt.push_back(pt);
                ev_.vld_mu_tight_iso_abs.push_back(iso_abs);
                ev_.vld_mu_tight_iso_rel.push_back(iso_rel);
                ev_.vld_mu_tight_dxy.push_back(dxy);
                ev_.vld_mu_tight_dz.push_back(dz);
                ev_.vld_mu_tight_pt_iso_abs->Fill(pt, iso_abs);
                ev_.vld_mu_tight_pt_iso_rel->Fill(pt, iso_rel);
                ev_.vld_mu_tight_pt_dxy->Fill(pt, dxy);
                ev_.vld_mu_tight_pt_dz->Fill(pt, dz);
            }
        }
    }

    // Electrons
    for (size_t i=0; i<elecs->size(); ++i) {

        // Only select good electrons
        if (!isGoodElecSOS(elecs->at(i), conversions, beamspot, vertices)){ continue; }

        // Only select electrons above certain pT
        //if (elecs->at(i).pt() < ev_.el_pt_lo){ continue; }

        ev_.nLep++;
        ev_.nEl++;
        if (elecs->at(i).pt() < ev_.el_pt_hi){
            ev_.nSoftLep++;
            ev_.nSoftEl++;
        }

        // Fill electron variables
        if (ev_.el1_pt.size() == 0){
            ev_.el1_pt.push_back(elecs->at(i).pt());
            ev_.el1_eta.push_back(elecs->at(i).eta());
            ev_.el1_phi.push_back(elecs->at(i).phi());
            ev_.el1_q.push_back(elecs->at(i).charge());
        }else if (ev_.el2_pt.size() == 0){
            ev_.el2_pt.push_back(elecs->at(i).pt());
            ev_.el2_eta.push_back(elecs->at(i).eta());
            ev_.el2_phi.push_back(elecs->at(i).phi());
            ev_.el2_q.push_back(elecs->at(i).charge());
        }
    }

    // Muons
    for (size_t i=0; i<muons->size(); ++i) {

        // Only select good muons
        if (!isGoodMuonSOS(muons->at(i), vertices, prVtx)){ continue; }

        // Only select muons above certain pT
        //if (muons->at(i).pt() < ev_.mu_pt_lo){ continue; }

        ev_.nLep++;
        ev_.nMu++;
        if (muons->at(i).pt() < ev_.mu_pt_hi){
            ev_.nSoftLep++;
            ev_.nSoftMu++;
        }

        // Fill muon variables
        if (ev_.mu1_pt.size() == 0){
            ev_.mu1_pt.push_back(muons->at(i).pt());
            ev_.mu1_eta.push_back(muons->at(i).eta());
            ev_.mu1_phi.push_back(muons->at(i).phi());
            ev_.mu1_q.push_back(muons->at(i).charge());
        }else if (ev_.mu2_pt.size() == 0){
            ev_.mu2_pt.push_back(muons->at(i).pt());
            ev_.mu2_eta.push_back(muons->at(i).eta());
            ev_.mu2_phi.push_back(muons->at(i).phi());
            ev_.mu2_q.push_back(muons->at(i).charge());
        }
    }

    // Fill leptons
    // Put pT and eta into vector of vector for sorting
    std::vector<std::vector<double>> lepvec;
    if (ev_.el1_pt.size() != 0){
        lepvec.push_back({ev_.el1_pt.at(0), ev_.el1_eta.at(0), ev_.el1_phi.at(0), ev_.mass_el});
    }
    if (ev_.el2_pt.size() != 0){
        lepvec.push_back({ev_.el2_pt.at(0), ev_.el2_eta.at(0), ev_.el2_phi.at(0), ev_.mass_el});
    }
    if (ev_.mu1_pt.size() != 0){
        lepvec.push_back({ev_.mu1_pt.at(0), ev_.mu1_eta.at(0), ev_.mu1_phi.at(0), ev_.mass_mu});
    }
    if (ev_.mu2_pt.size() != 0){
        lepvec.push_back({ev_.mu2_pt.at(0), ev_.mu2_eta.at(0), ev_.mu2_phi.at(0), ev_.mass_mu});
    }
    // By definition, this sorts by the first element of the vector (in this case pT)
    if (lepvec.size() > 1){
        std::sort(begin(lepvec), end(lepvec));
        std::reverse(begin(lepvec), end(lepvec));
    }
    if (lepvec.size() >= 1){
        ev_.lep1_pt.push_back(lepvec[0][0]);
        ev_.lep1_eta.push_back(lepvec[0][1]);
        ev_.lep1_phi.push_back(lepvec[0][2]);
        ev_.lep1_mass.push_back(lepvec[0][3]);
    }
    if (lepvec.size() >= 2){
        ev_.lep2_pt.push_back(lepvec[1][0]);
        ev_.lep2_eta.push_back(lepvec[1][1]);
        ev_.lep2_phi.push_back(lepvec[1][2]);
        ev_.lep2_mass.push_back(lepvec[1][3]);
    }
    lepvec.clear();

    // Fill pt of two leptons
    if (ev_.lep2_pt.size() != 0){
        TLorentzVector l1, l2;
        l1.SetPtEtaPhiM(ev_.lep1_pt[0], ev_.lep1_eta[0], ev_.lep1_phi[0], ev_.lep1_mass[0]);
        l2.SetPtEtaPhiM(ev_.lep2_pt[0], ev_.lep2_eta[0], ev_.lep2_phi[0], ev_.lep2_mass[0]);
        ev_.pt2l.push_back((l1 + l2).Pt());
    }

    // Is a same flavour opposite sign lepton pair present?
    for (size_t i=0; i<elecs->size(); ++i){
        // Only select good electrons
        if (!isGoodElecSOS(elecs->at(i), conversions, beamspot, vertices)){ continue; }

        for (size_t j=i+1; j<elecs->size(); ++j){
            // Only select good electrons
            if (!isGoodElecSOS(elecs->at(j), conversions, beamspot, vertices)){ continue; }

            // is SFOS?
            if (elecs->at(i).charge()*elecs->at(j).charge() > 0){ continue; }

            // Are leptons soft?
            if (elecs->at(i).pt() > ev_.el_pt_hi || elecs->at(j).pt() > ev_.el_pt_hi){ continue; }

            // Mll for soft SFOS
            TLorentzVector l1, l2;
            l1.SetPtEtaPhiM(elecs->at(i).pt(), elecs->at(i).eta(), elecs->at(i).phi(), ev_.mass_el);
            l2.SetPtEtaPhiM(elecs->at(j).pt(), elecs->at(j).eta(), elecs->at(j).phi(), ev_.mass_el);
            double mll = (l1+l2).M();
            if (ev_.mllMin.size() == 0){
                ev_.mllMin.push_back(mll);
            }else if (ev_.mllMin.at(0) > mll){
                ev_.mllMin.at(0) = mll;
            }
            if (ev_.mllMax.size() == 0){
                ev_.mllMax.push_back(mll);
            }else if (ev_.mllMax.at(0) < mll){
                ev_.mllMax.at(0) = mll;
            }
        }
    }

    for (size_t i=0; i<muons->size(); ++i){
        // Only select good muons
        if (!isGoodMuonSOS(muons->at(i), vertices, prVtx)){ continue; }

        for (size_t j=i+1; j<muons->size(); ++j){
            // Only select good muons
            if (!isGoodMuonSOS(muons->at(j), vertices, prVtx)){ continue; }

            // is SFOS?
            if (muons->at(i).charge()*muons->at(j).charge() > 0){ continue; }

            // Are leptons soft?
            if (muons->at(i).pt() > ev_.mu_pt_hi || muons->at(j).pt() > ev_.mu_pt_hi){ continue; }

            // Mll for soft SFOS
            TLorentzVector l1, l2;
            l1.SetPtEtaPhiM(muons->at(i).pt(), muons->at(i).eta(), muons->at(i).phi(), ev_.mass_mu);
            l2.SetPtEtaPhiM(muons->at(j).pt(), muons->at(j).eta(), muons->at(j).phi(), ev_.mass_mu);
            double mll = (l1+l2).M();
            if (ev_.mllMin.size() == 0){
                ev_.mllMin.push_back(mll);
            }else if (ev_.mllMin.at(0) > mll){
                ev_.mllMin.at(0) = mll;
            }
            if (ev_.mllMax.size() == 0){
                ev_.mllMax.push_back(mll);
            }else if (ev_.mllMax.at(0) < mll){
                ev_.mllMax.at(0) = mll;
            }
        }
    }

    // Jets
    for (size_t i=0; i<jets->size(); ++i) {
        if (jets->at(i).pt() < ev_.jet_pt_lo) continue;
        //if (fabs(jets->at(i).eta()) > 5) continue;

        bool overlaps = false;
        for (size_t j = 0; j < elecs->size(); j++) {
            if (fabs(jets->at(i).pt()-elecs->at(j).pt()) < 0.01*elecs->at(j).pt() && ROOT::Math::VectorUtil::DeltaR(elecs->at(j).p4(),jets->at(i).p4()) < 0.01) {
                overlaps = true;
                break;
            }
        }
        if (overlaps) continue;
        for (size_t j = 0; j < muons->size(); j++) {
            if (fabs(jets->at(i).pt()-muons->at(j).pt()) < 0.01*muons->at(j).pt() && ROOT::Math::VectorUtil::DeltaR(muons->at(j).p4(),jets->at(i).p4()) < 0.01) {
                overlaps = true;
                break;
            }
        }
        if (overlaps) continue;

        // Only select good jets
        if (!isGoodJetSOS(jets->at(i))){ continue; }

        ev_.nJet++;
        if (ev_.jet1_pt.size() == 0){
            ev_.jet1_pt.push_back(jets->at(i).pt());
            ev_.jet1_eta.push_back(jets->at(i).eta());
            ev_.jet1_phi.push_back(jets->at(i).phi());
            ev_.jet1_mass.push_back(jets->at(i).mass());
        }
        ev_.ht += jets->at(i).pt();

        double deepcsv = jets->at(i).bDiscriminator("pfDeepCSVJetTags:probb") +
            jets->at(i).bDiscriminator("pfDeepCSVJetTags:probbb");
        bool isMediumDeepCSV = deepcsv > deepThres_[1];
        if (!isMediumDeepCSV){ continue; }

        ev_.nBJet++;
    }

    // MET
    if (mets->size() > 0) {
        ev_.met = mets->at(0).pt();
    }

    //// Muons
    //ev_.nlm = 0;
    //ev_.ntm = 0;

    //for (size_t i = 0; i < muons->size(); i++) {
    //  if (muons->at(i).pt() < 2.) continue;
    //  if (fabs(muons->at(i).eta()) > 2.8) continue;

    //  // Loose ID
    //  double dPhiCut = std::min(std::max(1.2/muons->at(i).p(),1.2/100),0.056);
    //  double dPhiBendCut = std::min(std::max(0.2/muons->at(i).p(),0.2/100),0.0096);
    //  bool isLoose = (fabs(muons->at(i).eta()) < 2.4 && muon::isLooseMuon(muons->at(i))) || (fabs(muons->at(i).eta()) > 2.4 && isME0MuonSelNew(muons->at(i), 0.077, dPhiCut, dPhiBendCut));

    //  // Medium ID -- needs to be updated
    //  bool ipxy = false, ipz = false, validPxlHit = false, highPurity = false;
    //  if (muons->at(i).innerTrack().isNonnull()){
    //      ipxy = std::abs(muons->at(i).muonBestTrack()->dxy(vertices->at(prVtx).position())) < 0.2;
    //      ipz = std::abs(muons->at(i).muonBestTrack()->dz(vertices->at(prVtx).position())) < 0.5;
    //      validPxlHit = muons->at(i).innerTrack()->hitPattern().numberOfValidPixelHits() > 0;
    //      highPurity = muons->at(i).innerTrack()->quality(reco::Track::highPurity);
    //  }
    //  // bool isMedium = (fabs(muons->at(i).eta()) < 2.4 && muon::isMediumMuon(muons->at(i))) || (fabs(muons->at(i).eta()) > 2.4 && isME0MuonSelNew(muons->at(i), 0.077, dPhiCut, dPhiBendCut) && ipxy && ipz && validPxlHit && highPurity);

    //  // Tight ID
    //  dPhiCut = std::min(std::max(1.2/muons->at(i).p(),1.2/100),0.032);
    //  dPhiBendCut = std::min(std::max(0.2/muons->at(i).p(),0.2/100),0.0041);
    //  bool isTight = (fabs(muons->at(i).eta()) < 2.4 && vertices->size() > 0 && muon::isTightMuon(muons->at(i),vertices->at(prVtx))) || (fabs(muons->at(i).eta()) > 2.4 && isME0MuonSelNew(muons->at(i), 0.048, dPhiCut, dPhiBendCut) && ipxy && ipz && validPxlHit && highPurity);

    //  if (!isLoose) continue;

    //  ev_.lm_ch[ev_.nlm]     = muons->at(i).charge();
    //  ev_.lm_pt[ev_.nlm]     = muons->at(i).pt();
    //  ev_.lm_phi[ev_.nlm]    = muons->at(i).phi();
    //  ev_.lm_eta[ev_.nlm]    = muons->at(i).eta();
    //  ev_.lm_mass[ev_.nlm]   = muons->at(i).mass();
    //  ev_.lm_relIso[ev_.nlm] = (muons->at(i).puppiNoLeptonsChargedHadronIso() + muons->at(i).puppiNoLeptonsNeutralHadronIso() + muons->at(i).puppiNoLeptonsPhotonIso()) / muons->at(i).pt();
    //  ev_.lm_g[ev_.nlm] = -1;
    //  for (int ig = 0; ig < ev_.ngl; ig++) {
    //    if (abs(ev_.gl_pid[ig]) != 13) continue;
    //    if (reco::deltaR(ev_.gl_eta[ig],ev_.gl_phi[ig],ev_.lm_eta[ev_.nlm],ev_.lm_phi[ev_.nlm]) > 0.4) continue;
    //    ev_.lm_g[ev_.nlm]    = ig;
    //  }
    //  ev_.nlm++;

    //  if (!isTight) continue;

    //  ev_.tm_ch[ev_.ntm]     = muons->at(i).charge();
    //  ev_.tm_pt[ev_.ntm]     = muons->at(i).pt();
    //  ev_.tm_phi[ev_.ntm]    = muons->at(i).phi();
    //  ev_.tm_eta[ev_.ntm]    = muons->at(i).eta();
    //  ev_.tm_mass[ev_.ntm]   = muons->at(i).mass();
    //  ev_.tm_relIso[ev_.ntm] = (muons->at(i).puppiNoLeptonsChargedHadronIso() + muons->at(i).puppiNoLeptonsNeutralHadronIso() + muons->at(i).puppiNoLeptonsPhotonIso()) / muons->at(i).pt();
    //  ev_.tm_g[ev_.ntm] = -1;
    //  for (int ig = 0; ig < ev_.ngl; ig++) {
    //    if (abs(ev_.gl_pid[ig]) != 13) continue;
    //    if (reco::deltaR(ev_.gl_eta[ig],ev_.gl_phi[ig],ev_.tm_eta[ev_.ntm],ev_.tm_phi[ev_.ntm]) > 0.4) continue;
    //    ev_.tm_g[ev_.ntm]    = ig;
    //  }
    //  ev_.ntm++;
    //}

    // Electrons

    //ev_.nle = 0;
    //ev_.nte = 0;

    //for (size_t i = 0; i < elecs->size(); i++) {
    //  if (elecs->at(i).pt() < 10.) continue;
    //  if (fabs(elecs->at(i).eta()) > 3.) continue;

    //  bool isLoose = isLooseElec(elecs->at(i),conversions,beamspot);
    //  // bool isMedium = isMediumElec(elecs->at(i),conversions,beamspot);
    //  bool isTight = isTightElec(elecs->at(i),conversions,beamspot);

    //  if (!isLoose) continue;

    //  ev_.le_ch[ev_.nle]     = elecs->at(i).charge();
    //  ev_.le_pt[ev_.nle]     = elecs->at(i).pt();
    //  ev_.le_phi[ev_.nle]    = elecs->at(i).phi();
    //  ev_.le_eta[ev_.nle]    = elecs->at(i).eta();
    //  ev_.le_mass[ev_.nle]   = elecs->at(i).mass();
    //  ev_.le_relIso[ev_.nle] = (elecs->at(i).puppiNoLeptonsChargedHadronIso() + elecs->at(i).puppiNoLeptonsNeutralHadronIso() + elecs->at(i).puppiNoLeptonsPhotonIso()) / elecs->at(i).pt();
    //  ev_.le_g[ev_.nle] = -1;
    //  for (int ig = 0; ig < ev_.ngl; ig++) {
    //    if (abs(ev_.gl_pid[ig]) != 11) continue;
    //    if (reco::deltaR(ev_.gl_eta[ig],ev_.gl_phi[ig],ev_.le_eta[ev_.nle],ev_.le_phi[ev_.nle]) > 0.4) continue;
    //    ev_.le_g[ev_.nle]    = ig;
    //  }
    //  ev_.nle++;

    //  if (!isTight) continue;

    //  if (ev_.el1_pt.size() == 0){
    //      ev_.el1_pt.push_back(elecs->at(i).pt());
    //  }
    //  ev_.te_ch[ev_.nte]     = elecs->at(i).charge();
    //  ev_.te_pt[ev_.nte]     = elecs->at(i).pt();
    //  ev_.te_phi[ev_.nte]    = elecs->at(i).phi();
    //  ev_.te_eta[ev_.nte]    = elecs->at(i).eta();
    //  ev_.te_mass[ev_.nte]   = elecs->at(i).mass();
    //  ev_.te_relIso[ev_.nte] = (elecs->at(i).puppiNoLeptonsChargedHadronIso() + elecs->at(i).puppiNoLeptonsNeutralHadronIso() + elecs->at(i).puppiNoLeptonsPhotonIso()) / elecs->at(i).pt();
    //  ev_.te_g[ev_.nte] = -1;
    //  for (int ig = 0; ig < ev_.ngl; ig++) {
    //    if (abs(ev_.gl_pid[ig]) != 11) continue;
    //    if (reco::deltaR(ev_.gl_eta[ig],ev_.gl_phi[ig],ev_.te_eta[ev_.nte],ev_.te_phi[ev_.nte]) > 0.4) continue;
    //    ev_.te_g[ev_.nte]    = ig;
    //  }
    //  ev_.nte++;
    //}

    //// Jets
    //ev_.nj = 0;
    //for (size_t i =0; i < jets->size(); i++) {
    //  if (jets->at(i).pt() < 20.) continue;
    //  if (fabs(jets->at(i).eta()) > 5) continue;

    //  bool overlaps = false;
    //  for (size_t j = 0; j < elecs->size(); j++) {
    //    if (fabs(jets->at(i).pt()-elecs->at(j).pt()) < 0.01*elecs->at(j).pt() && ROOT::Math::VectorUtil::DeltaR(elecs->at(j).p4(),jets->at(i).p4()) < 0.01) {
    //      overlaps = true;
    //      break;
    //    }
    //  }
    //  if (overlaps) continue;
    //  for (size_t j = 0; j < muons->size(); j++) {
    //    if (fabs(jets->at(i).pt()-muons->at(j).pt()) < 0.01*muons->at(j).pt() && ROOT::Math::VectorUtil::DeltaR(muons->at(j).p4(),jets->at(i).p4()) < 0.01) {
    //      overlaps = true;
    //      break;
    //    }
    //  }
    //  if (overlaps) continue;

    //  pat::strbitset retLoose = jetIDLoose_.getBitTemplate();
    //  retLoose.set(false);
    //  bool isLoose = jetIDLoose_(jets->at(i), retLoose);
    //  pat::strbitset retTight = jetIDTight_.getBitTemplate();
    //  retTight.set(false);
    //  bool isTight = jetIDTight_(jets->at(i), retTight);

    //  double mvav2   = jets->at(i).bDiscriminator("pfCombinedMVAV2BJetTags");
    //  bool isLooseMVAv2  = mvav2 > mvaThres_[0];
    //  bool isMediumMVAv2 = mvav2 > mvaThres_[1];
    //  bool isTightMVAv2  = mvav2 > mvaThres_[2];
    //  double deepcsv = jets->at(i).bDiscriminator("pfDeepCSVJetTags:probb") +
    //                          jets->at(i).bDiscriminator("pfDeepCSVJetTags:probbb");
    //  bool isLooseDeepCSV  = deepcsv > deepThres_[0];
    //  bool isMediumDeepCSV = deepcsv > deepThres_[1];
    //  bool isTightDeepCSV  = deepcsv > deepThres_[2];

    //  ev_.j_id[ev_.nj]      = (isTight | (isLoose<<1));
    //  ev_.j_pt[ev_.nj]      = jets->at(i).pt();
    //  ev_.j_phi[ev_.nj]     = jets->at(i).phi();
    //  ev_.j_eta[ev_.nj]     = jets->at(i).eta();
    //  ev_.j_mass[ev_.nj]    = jets->at(i).mass();
    //  ev_.j_mvav2[ev_.nj]   = (isTightMVAv2 | (isMediumMVAv2<<1) | (isLooseMVAv2<<2));
    //  ev_.j_deepcsv[ev_.nj] = (isTightDeepCSV | (isMediumDeepCSV<<1) | (isLooseDeepCSV<<2));
    //  ev_.j_flav[ev_.nj]    = jets->at(i).partonFlavour();
    //  ev_.j_hadflav[ev_.nj] = jets->at(i).hadronFlavour();
    //  ev_.j_pid[ev_.nj]     = (jets->at(i).genParton() ? jets->at(i).genParton()->pdgId() : 0);
    //  ev_.j_g[ev_.nj] = -1;
    //  for (int ig = 0; ig < ev_.ngj; ig++) {
    //    if (reco::deltaR(ev_.gj_eta[ig],ev_.gj_phi[ig],ev_.j_eta[ev_.nj],ev_.j_phi[ev_.nj]) > 0.4) continue;
    //    ev_.j_g[ev_.nj]     = ig;
    //    break;
    //  }
    //  ev_.nj++;

    //}
    //
    //// MET
    //ev_.nmet = 0;
    //if (mets->size() > 0) {
    //  ev_.met_pt[ev_.nmet]  = mets->at(0).pt();
    //  ev_.met_eta[ev_.nmet] = mets->at(0).eta();
    //  ev_.met_phi[ev_.nmet] = mets->at(0).phi();
    //  ev_.nmet++;
    //}

}

// ------------ method called for each event  ------------
    void
MiniFromPat::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

    // Clear vectors
    ev_.el1_pt.clear();
    ev_.el1_eta.clear();
    ev_.el1_phi.clear();
    ev_.el1_q.clear();
    ev_.el2_pt.clear();
    ev_.el2_eta.clear();
    ev_.el2_phi.clear();
    ev_.el2_q.clear();
    ev_.mu1_pt.clear();
    ev_.mu1_eta.clear();
    ev_.mu1_phi.clear();
    ev_.mu1_q.clear();
    ev_.mu2_pt.clear();
    ev_.mu2_eta.clear();
    ev_.mu2_phi.clear();
    ev_.mu2_q.clear();
    ev_.lep1_pt.clear();
    ev_.lep1_eta.clear();
    ev_.lep1_phi.clear();
    ev_.lep1_mass.clear();
    ev_.lep2_pt.clear();
    ev_.lep2_eta.clear();
    ev_.lep2_phi.clear();
    ev_.lep2_mass.clear();
    //ev_.el1_pt_truth.clear();
    //ev_.el1_eta_truth.clear();
    //ev_.el1_phi_truth.clear();
    //ev_.el1_q_truth.clear();
    //ev_.el2_pt_truth.clear();
    //ev_.el2_eta_truth.clear();
    //ev_.el2_phi_truth.clear();
    //ev_.el2_q_truth.clear();
    //ev_.mu1_pt_truth.clear();
    //ev_.mu1_eta_truth.clear();
    //ev_.mu1_phi_truth.clear();
    //ev_.mu1_q_truth.clear();
    //ev_.mu2_pt_truth.clear();
    //ev_.mu2_eta_truth.clear();
    //ev_.mu2_phi_truth.clear();
    //ev_.mu2_q_truth.clear();
    //ev_.lep1_pt_truth.clear();
    //ev_.lep1_eta_truth.clear();
    //ev_.lep1_phi_truth.clear();
    //ev_.lep1_mass_truth.clear();
    //ev_.lep2_pt_truth.clear();
    //ev_.lep2_eta_truth.clear();
    //ev_.lep2_phi_truth.clear();
    //ev_.lep2_mass_truth.clear();
    ev_.jet1_pt.clear();
    ev_.jet1_eta.clear();
    ev_.jet1_phi.clear();
    ev_.jet1_mass.clear();
    ev_.jet1_pt_truth.clear();
    ev_.jet1_eta_truth.clear();
    ev_.jet1_phi_truth.clear();
    ev_.jet1_mass_truth.clear();
    ev_.mllMin.clear();
    ev_.mllMax.clear();
    ev_.mt1.clear();
    ev_.mt2.clear();
    ev_.pt2l.clear();
    ev_.vld_el_pt.clear();
    ev_.vld_el_tight_pt.clear();
    ev_.vld_el_is_tight.clear();
    ev_.vld_el_iso_abs.clear();
    ev_.vld_el_tight_iso_abs.clear();
    ev_.vld_el_iso_rel.clear();
    ev_.vld_el_tight_iso_rel.clear();
    ev_.vld_el_dxy.clear();
    ev_.vld_el_tight_dxy.clear();
    ev_.vld_el_dz.clear();
    ev_.vld_el_tight_dz.clear();
    ev_.vld_mu_pt.clear();
    ev_.vld_mu_tight_pt.clear();
    ev_.vld_mu_is_tight.clear();
    ev_.vld_mu_iso_abs.clear();
    ev_.vld_mu_tight_iso_abs.clear();
    ev_.vld_mu_iso_rel.clear();
    ev_.vld_mu_tight_iso_rel.clear();
    ev_.vld_mu_dxy.clear();
    ev_.vld_mu_tight_dxy.clear();
    ev_.vld_mu_dz.clear();
    ev_.vld_mu_tight_dz.clear();

    ev_.nLep = ev_.nEl = ev_.nMu = 0;
    ev_.nSoftLep = ev_.nSoftEl = ev_.nSoftMu = 0;
    ev_.nJet = ev_.nBJet = 0;
    ev_.met = ev_.ht = 0.;

    //analyze the event
    if(!iEvent.isRealData()) genAnalysis(iEvent, iSetup);
    recoAnalysis(iEvent, iSetup);

    //save event if at least one lepton at gen or reco level
    ev_.run     = iEvent.id().run();
    ev_.lumi    = iEvent.luminosityBlock();
    ev_.event   = iEvent.id().event();
    ev_.genWeight = 1.;
    //t_event_->Fill();
    //t_genParts_->Fill();
    //t_vertices_->Fill();
    //t_genJets_->Fill();
    //t_looseElecs_->Fill();
    //t_tightElecs_->Fill();
    //t_looseMuons_->Fill();
    //t_tightMuons_->Fill();
    //t_puppiJets_->Fill();
    //t_puppiMET_->Fill();
    t_tree_->Fill();

}


// ------------ method check that an e passes loose ID ----------------------------------
    bool
MiniFromPat::isLooseElec(const pat::Electron & patEl, edm::Handle<reco::ConversionCollection> conversions, const reco::BeamSpot beamspot)
{
    if (fabs(patEl.superCluster()->eta()) > 1.479 && fabs(patEl.superCluster()->eta()) < 1.556) return false;
    if (patEl.full5x5_sigmaIetaIeta() > 0.02992) return false;
    if (fabs(patEl.deltaEtaSuperClusterTrackAtVtx()) > 0.004119) return false;
    if (fabs(patEl.deltaPhiSuperClusterTrackAtVtx()) > 0.05176) return false;
    if (patEl.hcalOverEcal() > 6.741) return false;
    if (patEl.pfIsolationVariables().sumChargedHadronPt / patEl.pt() > 2.5) return false;
    double Ooemoop = 999.;
    if (patEl.ecalEnergy() == 0) Ooemoop = 0.;
    else if (!std::isfinite(patEl.ecalEnergy())) Ooemoop = 998.;
    else Ooemoop = fabs(1./patEl.ecalEnergy() - patEl.eSuperClusterOverP()/patEl.ecalEnergy());
    if (Ooemoop > 73.76) return false;
    if (ConversionTools::hasMatchedConversion(patEl, conversions, beamspot.position())) return false;
    return true;
}

// ------------ method check that an e passes medium ID ----------------------------------
    bool
MiniFromPat::isMediumElec(const pat::Electron & patEl, edm::Handle<reco::ConversionCollection> conversions, const reco::BeamSpot beamspot)
{
    if (fabs(patEl.superCluster()->eta()) > 1.479 && fabs(patEl.superCluster()->eta()) < 1.556) return false;
    if (patEl.full5x5_sigmaIetaIeta() > 0.01609) return false;
    if (fabs(patEl.deltaEtaSuperClusterTrackAtVtx()) > 0.001766) return false;
    if (fabs(patEl.deltaPhiSuperClusterTrackAtVtx()) > 0.03130) return false;
    if (patEl.hcalOverEcal() > 7.371) return false;
    if (patEl.pfIsolationVariables().sumChargedHadronPt / patEl.pt() > 1.325) return false;
    double Ooemoop = 999.;
    if (patEl.ecalEnergy() == 0) Ooemoop = 0.;
    else if (!std::isfinite(patEl.ecalEnergy())) Ooemoop = 998.;
    else Ooemoop = fabs(1./patEl.ecalEnergy() - patEl.eSuperClusterOverP()/patEl.ecalEnergy());
    if (Ooemoop > 22.6) return false;
    if (ConversionTools::hasMatchedConversion(patEl, conversions, beamspot.position())) return false;
    return true;
}

// ------------ method check that an e passes tight ID ----------------------------------
    bool
MiniFromPat::isTightElec(const pat::Electron & patEl, edm::Handle<reco::ConversionCollection> conversions, const reco::BeamSpot beamspot)
{
    if (fabs(patEl.superCluster()->eta()) > 1.479 && fabs(patEl.superCluster()->eta()) < 1.556) return false;
    if (patEl.full5x5_sigmaIetaIeta() > 0.01614) return false;
    if (fabs(patEl.deltaEtaSuperClusterTrackAtVtx()) > 0.001322) return false;
    if (fabs(patEl.deltaPhiSuperClusterTrackAtVtx()) > 0.06129) return false;
    if (patEl.hcalOverEcal() > 4.492) return false;
    if (patEl.pfIsolationVariables().sumChargedHadronPt / patEl.pt() > 1.255) return false;
    double Ooemoop = 999.;
    if (patEl.ecalEnergy() == 0) Ooemoop = 0.;
    else if (!std::isfinite(patEl.ecalEnergy())) Ooemoop = 998.;
    else Ooemoop = fabs(1./patEl.ecalEnergy() - patEl.eSuperClusterOverP()/patEl.ecalEnergy());
    if (Ooemoop > 18.26) return false;
    if (ConversionTools::hasMatchedConversion(patEl, conversions, beamspot.position())) return false;
    return true;
}

bool MiniFromPat::isGoodElecSOS(const pat::Electron & patEl, edm::Handle<reco::ConversionCollection> conversions, const reco::BeamSpot beamspot, edm::Handle<std::vector<reco::Vertex>> vertices){
    // Default cuts
    if (!isTightElec(patEl, conversions, beamspot)){ return false; }

    // Isolation
    if (patEl.pfIsolationVariables().sumChargedHadronPt / patEl.pt() > 0.5){ return false; }
    if (patEl.pfIsolationVariables().sumChargedHadronPt > 5.){ return false; }

    // IP3D cuts
    const reco::Vertex &pv = vertices->front();
    double dxy = std::abs(patEl.gsfTrack()->dxy(pv.position()));
    double dz = std::abs(patEl.gsfTrack()->dz(pv.position()));
    if (sqrt(dxy*dz) > .01){ return false; }

    return true;
}

bool MiniFromPat::isGoodMuonSOS(const pat::Muon & patMu, edm::Handle<std::vector<reco::Vertex>> vertices, int prVtx){

    // Default cuts
    double dPhiCut = std::min(std::max(1.2/patMu.p(),1.2/100),0.032);
    double dPhiBendCut = std::min(std::max(0.2/patMu.p(),0.2/100),0.0041);
    bool ipxy = false, ipz = false, validPxlHit = false, highPurity = false;
    if (patMu.innerTrack().isNonnull()){
        ipxy = std::abs(patMu.muonBestTrack()->dxy(vertices->at(prVtx).position())) < 0.2;
        ipz = std::abs(patMu.muonBestTrack()->dz(vertices->at(prVtx).position())) < 0.5;
        validPxlHit = patMu.innerTrack()->hitPattern().numberOfValidPixelHits() > 0;
        highPurity = patMu.innerTrack()->quality(reco::Track::highPurity);
    }
    if (!(fabs(patMu.eta()) < 2.4 &&
                vertices->size() > 0 &&
                muon::isTightMuon(patMu,vertices->at(prVtx))) ||
            (fabs(patMu.eta()) > 2.4
             && isME0MuonSelNew(patMu, 0.048, dPhiCut, dPhiBendCut)
             && ipxy && ipz && validPxlHit && highPurity)){ return false; }

    // Isolation (might be revised)
    double muIsolation = patMu.puppiChargedHadronIso() + patMu.puppiNeutralHadronIso() + patMu.puppiPhotonIso();
    if (muIsolation / patMu.pt() > 0.5){ return false; }
    if (muIsolation > 5.){ return false; }

    // IP3D
    double dxy = std::abs(patMu.muonBestTrack()->dxy(vertices->at(prVtx).position()));
    double dz = std::abs(patMu.muonBestTrack()->dz(vertices->at(prVtx).position()));
    if (sqrt(dxy*dz) > .01){ return false; }

    return true;
}

bool MiniFromPat::isGoodJetSOS(const pat::Jet& patJet){
    // Only select tight jets
    pat::strbitset retTight = jetIDTight_.getBitTemplate();
    retTight.set(false);
    if (!jetIDTight_(patJet, retTight)){ return false; }

    // Additional requirements for ISR jets
    // (from PhysicsTools/Heppy/python/physicsobjects/Jet.py)
    double eta = abs(patJet.eta());
    double chf = patJet.chargedHadronEnergyFraction();
    double nhf = patJet.neutralHadronEnergyFraction();
    double phf = patJet.neutralEmEnergyFraction();
    double elf = patJet.chargedEmEnergyFraction();
    double chm = patJet.chargedHadronMultiplicity();
    double npr = patJet.chargedMultiplicity() + patJet.neutralMultiplicity();
    double npn = patJet.neutralMultiplicity();
    if (!(eta<3. and chf>.2 and nhf<.7 and phf<.7)){ return false; }
    if (!((eta<2.7 and ((npr>1 and phf<0.90 and nhf<0.90) and (eta>2.4 or (elf<0.99 and chf>0 and chm>0)))) or ((eta>2.7 and eta<3.0) and (nhf<0.98 and phf>0.01 and npn>2)) or (eta>3.0 and (phf<0.90 and npn>10)))){ return false; }
    return true;
}

bool MiniFromPat::isGoodElecTruthSOS(const pat::PackedGenParticle truthEl, const std::vector<size_t> jGenJets, const edm::Handle<std::vector<reco::GenJet>> genJets){
    // Isolation
    double absIso = 0.;
    for (size_t i = 0; i < jGenJets.size(); ++i) {
        if (ROOT::Math::VectorUtil::DeltaR(truthEl.p4(), genJets->at(jGenJets[i]).p4()) > .7){ continue; }
        std::vector<const reco::Candidate*> jconst = genJets->at(jGenJets[i]).getJetConstituentsQuick();
        for (size_t j = 0; j < jconst.size(); j++) {
            double deltaR = ROOT::Math::VectorUtil::DeltaR(truthEl.p4(),jconst[j]->p4());
            if (deltaR < .01 || deltaR > .4){ continue; }
            absIso = absIso + jconst[j]->pt();
        }
    }
    double relIso = absIso / truthEl.pt();

    if (relIso > .5){ return false; }
    if (absIso > 5.){ return false; }

    return true;
}

bool MiniFromPat::isGoodMuonTruthSOS(const pat::PackedGenParticle truthMu, const std::vector<size_t> jGenJets, const edm::Handle<std::vector<reco::GenJet>> genJets){
    // Isolation
    double absIso = 0.;
    for (size_t i = 0; i < jGenJets.size(); ++i) {
        if (ROOT::Math::VectorUtil::DeltaR(truthMu.p4(), genJets->at(jGenJets[i]).p4()) > .7){ continue; }
        std::vector<const reco::Candidate*> jconst = genJets->at(jGenJets[i]).getJetConstituentsQuick();
        for (size_t j = 0; j < jconst.size(); j++) {
            double deltaR = ROOT::Math::VectorUtil::DeltaR(truthMu.p4(),jconst[j]->p4());
            if (deltaR < .01 || deltaR > .4){ continue; }
            absIso = absIso + jconst[j]->pt();
        }
    }
    double relIso = absIso / truthMu.pt();

    if (relIso > .5){ return false; }
    if (absIso > 5.){ return false; }

    return true;
}

// ------------ method to improve ME0 muon ID ----------------
bool MiniFromPat::isME0MuonSel(reco::Muon muon, double pullXCut, double dXCut, double pullYCut, double dYCut, double dPhi)
{

    bool result = false;
    bool isME0 = muon.isME0Muon();

    if(isME0){

        double deltaX = 999;
        double deltaY = 999;
        double pullX = 999;
        double pullY = 999;
        double deltaPhi = 999;

        bool X_MatchFound = false, Y_MatchFound = false, Dir_MatchFound = false;

        const std::vector<reco::MuonChamberMatch>& chambers = muon.matches();
        for(std::vector<reco::MuonChamberMatch>::const_iterator chamber = chambers.begin(); chamber != chambers.end(); ++chamber){

            for (std::vector<reco::MuonSegmentMatch>::const_iterator segment = chamber->me0Matches.begin(); segment != chamber->me0Matches.end(); ++segment){

                if (chamber->detector() == 5){

                    deltaX   = std::abs(chamber->x - segment->x);
                    deltaY   = std::abs(chamber->y - segment->y);
                    pullX    = std::abs(chamber->x - segment->x) / std::sqrt(chamber->xErr + segment->xErr);
                    pullY    = std::abs(chamber->y - segment->y) / std::sqrt(chamber->yErr + segment->yErr);
                    deltaPhi = std::abs(atan(chamber->dXdZ) - atan(segment->dXdZ));

                }
            }
        }

        if ((pullX < pullXCut) || (deltaX < dXCut)) X_MatchFound = true;
        if ((pullY < pullYCut) || (deltaY < dYCut)) Y_MatchFound = true;
        if (deltaPhi < dPhi) Dir_MatchFound = true;

        result = X_MatchFound && Y_MatchFound && Dir_MatchFound;

    }

    return result;

}

    bool
MiniFromPat::isME0MuonSelNew(reco::Muon muon, double dEtaCut, double dPhiCut, double dPhiBendCut)
{

    bool result = false;
    bool isME0 = muon.isME0Muon();

    if(isME0){

        double deltaEta = 999;
        double deltaPhi = 999;
        double deltaPhiBend = 999;

        const std::vector<reco::MuonChamberMatch>& chambers = muon.matches();
        for( std::vector<reco::MuonChamberMatch>::const_iterator chamber = chambers.begin(); chamber != chambers.end(); ++chamber ){

            if (chamber->detector() == 5){

                for ( std::vector<reco::MuonSegmentMatch>::const_iterator segment = chamber->me0Matches.begin(); segment != chamber->me0Matches.end(); ++segment ){

                    LocalPoint trk_loc_coord(chamber->x, chamber->y, 0);
                    LocalPoint seg_loc_coord(segment->x, segment->y, 0);
                    LocalVector trk_loc_vec(chamber->dXdZ, chamber->dYdZ, 1);
                    LocalVector seg_loc_vec(segment->dXdZ, segment->dYdZ, 1);

                    const ME0Chamber * me0chamber = ME0Geometry_->chamber(chamber->id);

                    GlobalPoint trk_glb_coord = me0chamber->toGlobal(trk_loc_coord);
                    GlobalPoint seg_glb_coord = me0chamber->toGlobal(seg_loc_coord);

                    //double segDPhi = segment->me0SegmentRef->deltaPhi();
                    // need to check if this works
                    double segDPhi = me0chamber->computeDeltaPhi(seg_loc_coord, seg_loc_vec);
                    double trackDPhi = me0chamber->computeDeltaPhi(trk_loc_coord, trk_loc_vec);

                    deltaEta = std::abs(trk_glb_coord.eta() - seg_glb_coord.eta() );
                    deltaPhi = std::abs(trk_glb_coord.phi() - seg_glb_coord.phi() );
                    deltaPhiBend = std::abs(segDPhi - trackDPhi);

                    if (deltaEta < dEtaCut && deltaPhi < dPhiCut && deltaPhiBend < dPhiBendCut) result = true;

                }
            }
        }

    }

    return result;

}

// ------------ method called once each job just before starting event loop  ------------
    void
MiniFromPat::beginJob()
{
}

// ------------ method called once each run ----------------
    void
MiniFromPat::beginRun(edm::Run const& iRun, edm::EventSetup const& iSetup)
{
    edm::ESHandle<ME0Geometry> hGeom;
    iSetup.get<MuonGeometryRecord>().get(hGeom);
    ME0Geometry_ =( &*hGeom);
}

// ------------ method called when ending the processing of a run  ------------
    void
MiniFromPat::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called once each job just after ending the event loop  ------------
    void
MiniFromPat::endJob()
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MiniFromPat::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    //The following says we do not know what parameters are allowed so do no validation
    // Please change this to state exactly what you do use, even if it is no parameters
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MiniFromPat);

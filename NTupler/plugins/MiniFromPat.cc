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
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "DataFormats/PatCandidates/interface/Photon.h"

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

        bool isLooseElec(const pat::Electron & patEl);
        bool isMediumElec(const pat::Electron & patEl);
        bool isTightElec(const pat::Electron & patEl);
        bool isIsolatedElec(const pat::Electron & patEl);
        bool isGoodElecSOS(const pat::Electron & patEl, reco::Vertex primaryVertex);
        bool isLooseMuon(const pat::Muon & patMu, const edm::EventSetup& iSetup);
        bool isMediumMuon(const pat::Muon & patMu, reco::Vertex primaryVertex);
        bool isTightMuon(const pat::Muon & patMu, reco::Vertex primaryVertex, const edm::EventSetup& iSetup);
        bool isIsolatedMuon(const pat::Muon & patMu);
        bool isGoodMuonSOS(const pat::Muon & patMu, reco::Vertex primaryVertex, edm::EventSetup const& iSetup);
        bool isGoodJetSOS(const pat::Jet & patJet);
        bool isGoodElecTruthSOS(const pat::PackedGenParticle truthEl, const std::vector<size_t> jGenJets, const edm::Handle<std::vector<reco::GenJet>> genJets);
        bool isGoodMuonTruthSOS(const pat::PackedGenParticle truthMu, const std::vector<size_t> jGenJets, const edm::Handle<std::vector<reco::GenJet>> genJets);
        double getMTauTau(double met_pt, double met_phi, double l1_pt, double l1_eta, double l1_phi, double l2_pt, double l2_eta, double l2_phi);
        double DeltaR(double eta1, double eta2, double phi1, double phi2);
        double DeltaPhi(double phi1, double phi2);
        template <typename T> bool isMatched(const pat::PackedGenParticle truthEl, T particle);
        template <typename T> bool matchAny(const edm::Handle<std::vector<pat::PackedGenParticle>> genParts, T particle, bool hs);
        bool isHs(const pat::PackedGenParticle truthParticle, int pdgId);
        template <typename T> void pppWpidWstatus(const char* text, const size_t idx, const size_t noParticles, const T particle, const char* addText="") const ;
        template <typename T> void pppWisoWpassid(const char* text, const size_t idx, const size_t noParticles, const T particle, const int pid, const int status, const short int passID, const char* addText="") const ;
        template <typename T> void ppp(const char* text, const size_t idx, const size_t noParticles, const T particle, const int pid=-1, const int status=-1, const double iso=-1., const short int passID=-1., const char* addText="") const ;

        bool isME0MuonSelNew(const reco::Muon&, double, double, double, edm::EventSetup const& );

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
        edm::EDGetTokenT<GenEventInfoProduct> generatorToken_;
        edm::EDGetTokenT<LHEEventProduct> generatorlheToken_;
        //edm::EDGetTokenT<std::vector<pat::Photon>> photonsToken_;
        const ME0Geometry* ME0Geometry_;
        double mvaThres_[3];
        double deepThres_[3];

        //TTree *t_event_, *t_genParts_, *t_vertices_, *t_genJets_, *t_genPhotons_, *t_looseElecs_, *t_mediumElecs_, *t_tightElecs_, *t_looseMuons_, *t_tightMuons_, *t_puppiJets_, *t_puppiMET_, *t_loosePhotons_, *t_tightPhotons_;
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
    genPartsToken_(consumes<std::vector<pat::PackedGenParticle>>(iConfig.getParameter<edm::InputTag>("genParts"))),
    generatorToken_(consumes<GenEventInfoProduct>(edm::InputTag("generator"))),
    generatorlheToken_(consumes<LHEEventProduct>(edm::InputTag("externalLHEProducer","")))//,
    //photonsToken_(consumes<std::vector<pat::Photon>>(iConfig.getParameter<edm::InputTag>("photons")))
{
    ME0Geometry_=0;
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

    usesResource("TFileService");

    //t_event_      = fs_->make<TTree>("Event","Event");
    //t_genParts_   = fs_->make<TTree>("Particle","Particle");
    //t_genPhotons_ = fs_->make<TTree>("GenPhoton","GenPhoton");
    //t_vertices_   = fs_->make<TTree>("Vertex","Vertex");
    //t_genJets_    = fs_->make<TTree>("GenJet","GenJet");
    //t_looseElecs_ = fs_->make<TTree>("ElectronLoose","ElectronLoose");
    //t_mediumElecs_ = fs_->make<TTree>("ElectronMedium","ElectronMedium");
    //t_tightElecs_ = fs_->make<TTree>("ElectronTight","ElectronTight");
    //t_looseMuons_ = fs_->make<TTree>("MuonLoose","MuonLoose");
    //t_tightMuons_ = fs_->make<TTree>("MuonTight","MuonTight");
    //t_puppiJets_  = fs_->make<TTree>("JetPUPPI","JetPUPPI");
    //t_puppiMET_   = fs_->make<TTree>("PuppiMissingET","PuppiMissingET");
    //t_loosePhotons_ = fs_->make<TTree>("PhotonLoose","PhotonLoose");
    //t_tightPhotons_ = fs_->make<TTree>("PhotonTight","PhotonTight");
    //createMiniEventTree(t_event_, t_genParts_, t_vertices_, t_genJets_, t_genPhotons_, t_looseElecs_,
    //        t_mediumElecs_,t_tightElecs_, t_looseMuons_, t_tightMuons_, t_puppiJets_, t_puppiMET_,
    //        t_loosePhotons_, t_tightPhotons_, ev_);
    t_tree_ = fs_->make<TTree>("Delphes", "Delphes");
    createMiniEventTree(t_tree_, ev_);

}


MiniFromPat::~MiniFromPat()
{

    // do anything here that needs to be done at desctruction time
    // (e.g. close files, deallocate resources etc.)

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

    //// Jets
    //std::vector<size_t> jGenJets;
    //ev_.ngj = 0;
    //for (size_t i = 0; i < genJets->size(); i++) {
    //    if (ev_.ngj>=MiniEvent_t::maxjets) break;
    //    if (genJets->at(i).pt() < 20.) continue;
    //    if (fabs(genJets->at(i).eta()) > 5) continue;

    //    bool overlaps = false;
    //    for (size_t j = 0; j < genParts->size(); j++) {
    //        if (abs(genParts->at(j).pdgId()) != 11 && abs(genParts->at(j).pdgId()) != 13) continue;
    //        if (fabs(genJets->at(i).pt()-genParts->at(j).pt()) < 0.01*genParts->at(j).pt() && ROOT::Math::VectorUtil::DeltaR(genParts->at(j).p4(),genJets->at(i).p4()) < 0.01) {
    //            overlaps = true;
    //            break;
    //        }
    //    }
    //    if (overlaps) continue;
    //    jGenJets.push_back(i);

    //    ev_.gj_pt[ev_.ngj]   = genJets->at(i).pt();
    //    ev_.gj_phi[ev_.ngj]  = genJets->at(i).phi();
    //    ev_.gj_eta[ev_.ngj]  = genJets->at(i).eta();
    //    ev_.gj_mass[ev_.ngj] = genJets->at(i).mass();
    //    ev_.ngj++;
    //}

    //// Leptons
    //ev_.ngl = 0;
    //for (size_t i = 0; i < genParts->size(); i++) {
    //    if (ev_.ngl>=MiniEvent_t::maxpart) break;
    //    if (abs(genParts->at(i).pdgId()) != 11 && abs(genParts->at(i).pdgId()) != 13) continue;
    //    if (genParts->at(i).pt() < 10.) continue;
    //    if (fabs(genParts->at(i).eta()) > 3.) continue;
    //    double genIso = 0.;
    //    for (size_t j = 0; j < jGenJets.size(); j++) {
    //        if (ROOT::Math::VectorUtil::DeltaR(genParts->at(i).p4(),genJets->at(jGenJets[j]).p4()) > 0.7) continue;
    //        std::vector<const reco::Candidate *> jconst = genJets->at(jGenJets[j]).getJetConstituentsQuick();
    //        for (size_t k = 0; k < jconst.size(); k++) {
    //            double deltaR = ROOT::Math::VectorUtil::DeltaR(genParts->at(i).p4(),jconst[k]->p4());
    //            if (deltaR < 0.01 || deltaR > 0.4) continue;
    //            genIso = genIso + jconst[k]->pt();
    //        }
    //    }
    //    genIso = genIso / genParts->at(i).pt();
    //    ev_.gl_pid[ev_.ngl]    = genParts->at(i).pdgId();
    //    ev_.gl_ch[ev_.ngl]     = genParts->at(i).charge();
    //    ev_.gl_st[ev_.ngl]     = genParts->at(i).status();
    //    ev_.gl_p[ev_.ngl]      = genParts->at(i).p();
    //    ev_.gl_px[ev_.ngl]     = genParts->at(i).px();
    //    ev_.gl_py[ev_.ngl]     = genParts->at(i).py();
    //    ev_.gl_pz[ev_.ngl]     = genParts->at(i).pz();
    //    ev_.gl_nrj[ev_.ngl]    = genParts->at(i).energy();
    //    ev_.gl_pt[ev_.ngl]     = genParts->at(i).pt();
    //    ev_.gl_phi[ev_.ngl]    = genParts->at(i).phi();
    //    ev_.gl_eta[ev_.ngl]    = genParts->at(i).eta();
    //    ev_.gl_mass[ev_.ngl]   = genParts->at(i).mass();
    //    ev_.gl_relIso[ev_.ngl] = genIso;
    //    ev_.ngl++;
    //}

    //// Photons
    //ev_.ngp = 0;
    //for (size_t i = 0; i < genParts->size(); i++) {
    //    if (ev_.ngp>=MiniEvent_t::maxpart) break;
    //    if (abs(genParts->at(i).pdgId()) != 22) continue;
    //    if (genParts->at(i).pt() < 10.) continue;
    //    if (fabs(genParts->at(i).eta()) > 3.) continue;

    //    ev_.gp_st[ev_.ngp]     = genParts->at(i).status();
    //    ev_.gp_p[ev_.ngp]      = genParts->at(i).p();
    //    ev_.gp_px[ev_.ngp]     = genParts->at(i).px();
    //    ev_.gp_py[ev_.ngp]     = genParts->at(i).py();
    //    ev_.gp_pz[ev_.ngp]     = genParts->at(i).pz();
    //    ev_.gp_nrj[ev_.ngp]    = genParts->at(i).energy();
    //    ev_.gp_pt[ev_.ngp]     = genParts->at(i).pt();
    //    ev_.gp_phi[ev_.ngp]    = genParts->at(i).phi();
    //    ev_.gp_eta[ev_.ngp]    = genParts->at(i).eta();
    //    ev_.ngp++;
    //}

    //// Generator weights
    //ev_.g_nw = 0; ev_.g_w[0] = 1.0;
    //edm::Handle<GenEventInfoProduct> evt;
    //iEvent.getByToken( generatorToken_,evt);
    //if(evt.isValid()) {
    //    ev_.g_w[0] = evt->weight();
    //    ev_.g_nw++;
    //    for (unsigned int i = 1; i<evt->weights().size(); i++) {
    //        if (ev_.g_nw>=MiniEvent_t::maxweights) break;
    //        ev_.g_w[ev_.g_nw]=evt->weights()[i];
    //        ev_.g_nw++;
    //    }
    //}
    //// LHE weights
    //edm::Handle<LHEEventProduct> evet;
    //iEvent.getByToken(generatorlheToken_, evet);
    //if(evet.isValid()) {
    //    double asdd=evet->originalXWGTUP();
    //    for(unsigned int i=0  ; i<evet->weights().size();i++) {
    //        if (ev_.g_nw>=MiniEvent_t::maxweights) break;
    //        double asdde=evet->weights()[i].wgt;
    //        ev_.g_w[ev_.g_nw]=ev_.g_w[0]*asdde/asdd;
    //        ev_.g_nw++;
    //    }
    //}
}

// ------------ method to fill reco level pat -------------
    void
MiniFromPat::recoAnalysis(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    using namespace edm;

    Handle<std::vector<pat::PackedGenParticle>> genParts;
    iEvent.getByToken(genPartsToken_, genParts);
    Handle<std::vector<reco::GenJet>> genJets;
    iEvent.getByToken(genJetsToken_, genJets);

    Handle<std::vector<reco::Vertex>> vertices;
    iEvent.getByToken(verticesToken_, vertices);

    Handle<std::vector<pat::Electron>> elecs;
    iEvent.getByToken(elecsToken_, elecs);
    Handle<reco::ConversionCollection> conversions;
    iEvent.getByToken(convToken_, conversions);
    Handle<reco::BeamSpot> bsHandle;
    //iEvent.getByToken(bsToken_, bsHandle);
    //const reco::BeamSpot &beamspot = *bsHandle.product();

    Handle<std::vector<pat::Muon>> muons;
    iEvent.getByToken(muonsToken_, muons);

    Handle<std::vector<pat::MET>> mets;
    iEvent.getByToken(metsToken_, mets);

    Handle<std::vector<pat::Jet>> jets;
    iEvent.getByToken(jetsToken_, jets);

    //Handle<std::vector<pat::Photon>> photons;
    //iEvent.getByToken(photonsToken_, photons);

    // Vertices
    int prVtx = -1;
    //ev_.nvtx = 0;
    for (size_t i = 0; i < vertices->size(); i++) {
        if (vertices->at(i).isFake()) continue;
        if (vertices->at(i).ndof() <= 4) continue;
        if (prVtx < 0) prVtx = i;
        //ev_.v_pt2[ev_.nvtx] = vertices->at(i).p4().pt();
        //ev_.nvtx++;
    }
    if (prVtx < 0) return;
    reco::Vertex primaryVertex = vertices->at(prVtx);

    // Validation
    if (ev_.fill_vld){

        //// Count leptons coming from SUSY (pdgId >= 1000000)
        //unsigned int nGenElFromSusy = 0;
        //unsigned int nGenMuFromSusy = 0;
        //std::vector<double> nGenLepPts;
        //if (ev_.fill_vld){
        //    for (size_t i=0; i<genParts->size(); ++i){
        //        if (fabs(genParts->at(i).pdgId()) == 11){
        //            // Loop over all mothers to find SUSY mother particle (or not)
        //            const reco::Candidate* mom = genParts->at(i).mother(0);
        //            while (mom->numberOfMothers() != 0){
        //                mom = mom->mother(0);
        //                if (mom->pdgId() >= 1000000){
        //                    nGenElFromSusy++;
        //                    nGenLepPts.push_back(genParts->at(i).pt());
        //                    break;
        //                }
        //            }
        //        }else if (fabs(genParts->at(i).pdgId()) == 13){
        //            // Loop over all mothers to find SUSY mother particle (or not)
        //            const reco::Candidate* mom = genParts->at(i).mother(0);
        //            while (mom->numberOfMothers() != 0){
        //                mom = mom->mother(0);
        //                if (mom->pdgId() >= 1000000){
        //                    nGenMuFromSusy++;
        //                    nGenLepPts.push_back(genParts->at(i).pt());
        //                    break;
        //                }
        //            }
        //        }
        //    }
        //}

        // Electrons (TURNED OFF FOR NOW)
        //for (size_t i=0; i<elecs->size(); ++i){
        //    double pt = elecs->at(i).pt();
        //    double eta = fabs(elecs->at(i).eta());
        //    double absiso = elecs->at(i).pfIsolationVariables().sumChargedHadronPt;
        //    double reliso = absiso/pt;
        //    const reco::Vertex &pv = vertices->front();
        //    double dxy = std::abs(elecs->at(i).gsfTrack()->dxy(pv.position()));
        //    double dz = std::abs(elecs->at(i).gsfTrack()->dz(pv.position()));

        //    //ev_.vld_el_pt.push_back(pt);
        //    //ev_.vld_el_eta.push_back(eta);
        //    //ev_.vld_el_is_tight.push_back(isTightElec(elecs->at(i), conversions, beamspot, vertices));

        //    //ev_.vld_el_absiso.push_back(absiso);
        //    //ev_.vld_el_reliso.push_back(reliso);
        //    //ev_.vld_el_dxy.push_back(dxy);
        //    //ev_.vld_el_dz.push_back(dz);

        //    if (matchAny(genParts, elecs->at(i), true)){
        //        // Hard scattering electrons
        //        //ev_.vld_el_hs_pt.push_back(pt);
        //        //ev_.vld_el_hs_eta.push_back(eta);
        //        //ev_.vld_el_hs_pt_eta->Fill(pt, eta);
        //        //ev_.vld_el_hs_absiso.push_back(absiso);
        //        //ev_.vld_el_hs_reliso.push_back(reliso);
        //        //ev_.vld_el_hs_dxy.push_back(dxy);
        //        //ev_.vld_el_hs_dz.push_back(dz);
        //    }else if (matchAny(genParts, elecs->at(i), false)){
        //        // py8 electrons
        //        //ev_.vld_el_py8_pt.push_back(pt);
        //        //ev_.vld_el_py8_eta.push_back(eta);
        //        //ev_.vld_el_py8_pt_eta->Fill(pt, eta);
        //        //ev_.vld_el_py8_absiso.push_back(absiso);
        //        //ev_.vld_el_py8_reliso.push_back(reliso);
        //        //ev_.vld_el_py8_dxy.push_back(dxy);
        //        //ev_.vld_el_py8_dz.push_back(dz);
        //    }else{
        //        // Other electrons
        //        //ev_.vld_el_others_pt.push_back(pt);
        //        //ev_.vld_el_others_eta.push_back(eta);
        //        //ev_.vld_el_others_pt_eta->Fill(pt, eta);
        //        //ev_.vld_el_others_absiso.push_back(absiso);
        //        //ev_.vld_el_others_reliso.push_back(reliso);
        //        //ev_.vld_el_others_dxy.push_back(dxy);
        //        //ev_.vld_el_others_dz.push_back(dz);
        //    }

        //    if (isIsolatedElec(elecs->at(i))){

        //        //ev_.vld_el_iso_pt.push_back(pt);
        //        //ev_.vld_el_iso_eta.push_back(eta);

        //        //ev_.vld_el_iso_absiso.push_back(absiso);
        //        //ev_.vld_el_iso_reliso.push_back(reliso);
        //        //ev_.vld_el_iso_dxy.push_back(dxy);
        //        //ev_.vld_el_iso_dz.push_back(dz);

        //        if (matchAny(genParts, elecs->at(i), true)){
        //            // Hard scattering electrons
        //            //ev_.vld_el_iso_hs_pt.push_back(pt);
        //            //ev_.vld_el_iso_hs_eta.push_back(eta);
        //            //ev_.vld_el_iso_hs_pt_eta->Fill(pt, eta);
        //            //ev_.vld_el_iso_hs_absiso.push_back(absiso);
        //            //ev_.vld_el_iso_hs_reliso.push_back(reliso);
        //            //ev_.vld_el_iso_hs_dxy.push_back(dxy);
        //            //ev_.vld_el_iso_hs_dz.push_back(dz);
        //        }else if (matchAny(genParts, elecs->at(i), false)){
        //            // py8 electrons
        //            //ev_.vld_el_iso_py8_pt.push_back(pt);
        //            //ev_.vld_el_iso_py8_eta.push_back(eta);
        //            //ev_.vld_el_iso_py8_pt_eta->Fill(pt, eta);
        //            //ev_.vld_el_iso_py8_absiso.push_back(absiso);
        //            //ev_.vld_el_iso_py8_reliso.push_back(reliso);
        //            //ev_.vld_el_iso_py8_dxy.push_back(dxy);
        //            //ev_.vld_el_iso_py8_dz.push_back(dz);
        //        }else{
        //            // Other electrons
        //            //ev_.vld_el_iso_others_pt.push_back(pt);
        //            //ev_.vld_el_iso_others_eta.push_back(eta);
        //            //ev_.vld_el_iso_others_pt_eta->Fill(pt, eta);
        //            //ev_.vld_el_iso_others_absiso.push_back(absiso);
        //            //ev_.vld_el_iso_others_reliso.push_back(reliso);
        //            //ev_.vld_el_iso_others_dxy.push_back(dxy);
        //            //ev_.vld_el_iso_others_dz.push_back(dz);
        //        }
        //    }

        //    if (isGoodElecSOS(elecs->at(i), conversions, beamspot, vertices)){

        //        //ev_.vld_el_good_pt.push_back(pt);
        //        //ev_.vld_el_good_eta.push_back(eta);

        //        //ev_.vld_el_good_absiso.push_back(absiso);
        //        //ev_.vld_el_good_reliso.push_back(reliso);
        //        //ev_.vld_el_good_dxy.push_back(dxy);
        //        //ev_.vld_el_good_dz.push_back(dz);

        //        if (matchAny(genParts, elecs->at(i), true)){
        //            // Hard scattering electrons
        //            //ev_.vld_el_good_hs_pt.push_back(pt);
        //            //ev_.vld_el_good_hs_eta.push_back(eta);
        //            //ev_.vld_el_good_hs_pt_eta->Fill(pt, eta);
        //            //ev_.vld_el_good_hs_absiso.push_back(absiso);
        //            //ev_.vld_el_good_hs_reliso.push_back(reliso);
        //            //ev_.vld_el_good_hs_dxy.push_back(dxy);
        //            //ev_.vld_el_good_hs_dz.push_back(dz);
        //        }else if (matchAny(genParts, elecs->at(i), false)){
        //            // py8 electrons
        //            //ev_.vld_el_good_py8_pt.push_back(pt);
        //            //ev_.vld_el_good_py8_eta.push_back(eta);
        //            //ev_.vld_el_good_py8_pt_eta->Fill(pt, eta);
        //            //ev_.vld_el_good_py8_absiso.push_back(absiso);
        //            //ev_.vld_el_good_py8_reliso.push_back(reliso);
        //            //ev_.vld_el_good_py8_dxy.push_back(dxy);
        //            //ev_.vld_el_good_py8_dz.push_back(dz);
        //        }else{
        //            // Other electrons
        //            //ev_.vld_el_good_others_pt.push_back(pt);
        //            //ev_.vld_el_good_others_eta.push_back(eta);
        //            //ev_.vld_el_good_others_pt_eta->Fill(pt, eta);
        //            //ev_.vld_el_good_others_absiso.push_back(absiso);
        //            //ev_.vld_el_good_others_reliso.push_back(reliso);
        //            //ev_.vld_el_good_others_dxy.push_back(dxy);
        //            //ev_.vld_el_good_others_dz.push_back(dz);
        //        }
        //    }
        //}

        //// Truth electron validation
        //for (size_t i=0; i<genParts->size(); ++i){
        //    if (genParts->at(i).status() != 1){ continue; }
        //    if (fabs(genParts->at(i).pdgId()) != 11){ continue; }
        //    double pt = genParts->at(i).pt();
        //    double eta = fabs(genParts->at(i).eta());
        //    //ev_.vld_genel_pt.push_back(pt);
        //    //ev_.vld_genel_eta.push_back(eta);
        //    if (isHs(genParts->at(i), genParts->at(i).pdgId())){
        //        //ev_.vld_genel_hs_pt.push_back(pt);
        //        //ev_.vld_genel_hs_eta.push_back(eta);
        //    }else{
        //        //ev_.vld_genel_py8_pt.push_back(pt);
        //        //ev_.vld_genel_py8_eta.push_back(eta);
        //    }
        //}

        // Muons
        for (size_t i=0; i<muons->size(); ++i){
            double pt = muons->at(i).pt();
            double eta = fabs(muons->at(i).eta());
            double absiso = muons->at(i).isolationR03().sumPt;
            double reliso = absiso/pt;
            //double dxy = std::abs(muons->at(i).muonBestTrack()->dxy(vertices->at(prVtx).position()));
            //double dz = std::abs(muons->at(i).muonBestTrack()->dz(vertices->at(prVtx).position()));

            //ev_.vld_mu_pt.push_back(pt);
            //ev_.vld_mu_eta.push_back(eta);
            //ev_.vld_mu_is_tight.push_back(isTightMuon(muons->at(i), vertices, prVtx));

            //ev_.vld_mu_absiso.push_back(absiso);
            //ev_.vld_mu_reliso.push_back(reliso);
            //ev_.vld_mu_dxy.push_back(dxy);
            //ev_.vld_mu_dz.push_back(dz);

            if (matchAny(genParts, muons->at(i), true)){
                // Hard scattering muons
                ev_.vld_mu_hs_pt.push_back(pt);
                ev_.vld_mu_hs_eta.push_back(eta);
                //ev_.vld_mu_hs_pt_eta->Fill(pt, eta);
                //ev_.vld_mu_hs_absiso.push_back(absiso);
                ev_.vld_mu_hs_reliso.push_back(reliso);
                //ev_.vld_mu_hs_dxy.push_back(dxy);
                //ev_.vld_mu_hs_dz.push_back(dz);
            }else if (matchAny(genParts, muons->at(i), false)){
                // py8 muons
                ev_.vld_mu_py8_pt.push_back(pt);
                ev_.vld_mu_py8_eta.push_back(eta);
                //ev_.vld_mu_py8_pt_eta->Fill(pt, eta);
                //ev_.vld_mu_py8_absiso.push_back(absiso);
                ev_.vld_mu_py8_reliso.push_back(reliso);
                //ev_.vld_mu_py8_dxy.push_back(dxy);
                //ev_.vld_mu_py8_dz.push_back(dz);
            }else{
                // Other muons
                ev_.vld_mu_others_pt.push_back(pt);
                ev_.vld_mu_others_eta.push_back(eta);
                //ev_.vld_mu_others_pt_eta->Fill(pt, eta);
                //ev_.vld_mu_others_absiso.push_back(absiso);
                ev_.vld_mu_others_reliso.push_back(reliso);
                //ev_.vld_mu_others_dxy.push_back(dxy);
                //ev_.vld_mu_others_dz.push_back(dz);
            }

            if (isIsolatedMuon(muons->at(i))){

                //ev_.vld_mu_iso_pt.push_back(pt);
                //ev_.vld_mu_iso_eta.push_back(eta);

                //ev_.vld_mu_iso_absiso.push_back(absiso);
                //ev_.vld_mu_iso_reliso.push_back(reliso);
                //ev_.vld_mu_iso_dxy.push_back(dxy);
                //ev_.vld_mu_iso_dz.push_back(dz);

                if (matchAny(genParts, muons->at(i), true)){
                    // Hard scattering muons
                    ev_.vld_mu_iso_hs_pt.push_back(pt);
                    ev_.vld_mu_iso_hs_eta.push_back(eta);
                    //ev_.vld_mu_iso_hs_pt_eta->Fill(pt, eta);
                    //ev_.vld_mu_iso_hs_absiso.push_back(absiso);
                    ev_.vld_mu_iso_hs_reliso.push_back(reliso);
                    //ev_.vld_mu_iso_hs_dxy.push_back(dxy);
                    //ev_.vld_mu_iso_hs_dz.push_back(dz);
                }else if (matchAny(genParts, muons->at(i), false)){
                    // py8 muons
                    ev_.vld_mu_iso_py8_pt.push_back(pt);
                    ev_.vld_mu_iso_py8_eta.push_back(eta);
                    //ev_.vld_mu_iso_py8_pt_eta->Fill(pt, eta);
                    //ev_.vld_mu_iso_py8_absiso.push_back(absiso);
                    ev_.vld_mu_iso_py8_reliso.push_back(reliso);
                    //ev_.vld_mu_iso_py8_dxy.push_back(dxy);
                    //ev_.vld_mu_iso_py8_dz.push_back(dz);
                }else{
                    // Other muons
                    ev_.vld_mu_iso_others_pt.push_back(pt);
                    ev_.vld_mu_iso_others_eta.push_back(eta);
                    //ev_.vld_mu_iso_others_pt_eta->Fill(pt, eta);
                    //ev_.vld_mu_iso_others_absiso.push_back(absiso);
                    ev_.vld_mu_iso_others_reliso.push_back(reliso);
                    //ev_.vld_mu_iso_others_dxy.push_back(dxy);
                    //ev_.vld_mu_iso_others_dz.push_back(dz);
                }
            }

            if (isGoodMuonSOS(muons->at(i), primaryVertex, iSetup)){

                //ev_.vld_mu_good_pt.push_back(pt);
                //ev_.vld_mu_good_eta.push_back(eta);

                //ev_.vld_mu_good_absiso.push_back(absiso);
                //ev_.vld_mu_good_reliso.push_back(reliso);
                //ev_.vld_mu_good_dxy.push_back(dxy);
                //ev_.vld_mu_good_dz.push_back(dz);

                if (matchAny(genParts, muons->at(i), true)){
                    // Hard scattering muons
                    ev_.vld_mu_good_hs_pt.push_back(pt);
                    ev_.vld_mu_good_hs_eta.push_back(eta);
                    //ev_.vld_mu_good_hs_pt_eta->Fill(pt, eta);
                    //ev_.vld_mu_good_hs_absiso.push_back(absiso);
                    ev_.vld_mu_good_hs_reliso.push_back(reliso);
                    //ev_.vld_mu_good_hs_dxy.push_back(dxy);
                    //ev_.vld_mu_good_hs_dz.push_back(dz);
                }else if (matchAny(genParts, muons->at(i), false)){
                    // py8 muons
                    ev_.vld_mu_good_py8_pt.push_back(pt);
                    ev_.vld_mu_good_py8_eta.push_back(eta);
                    //ev_.vld_mu_good_py8_pt_eta->Fill(pt, eta);
                    //ev_.vld_mu_good_py8_absiso.push_back(absiso);
                    ev_.vld_mu_good_py8_reliso.push_back(reliso);
                    //ev_.vld_mu_good_py8_dxy.push_back(dxy);
                    //ev_.vld_mu_good_py8_dz.push_back(dz);
                }else{
                    // Other muons
                    ev_.vld_mu_good_others_pt.push_back(pt);
                    ev_.vld_mu_good_others_eta.push_back(eta);
                    //ev_.vld_mu_good_others_pt_eta->Fill(pt, eta);
                    //ev_.vld_mu_good_others_absiso.push_back(absiso);
                    ev_.vld_mu_good_others_reliso.push_back(reliso);
                    //ev_.vld_mu_good_others_dxy.push_back(dxy);
                    //ev_.vld_mu_good_others_dz.push_back(dz);
                }
            }
        }

        // Truth muon validation
        for (size_t i=0; i<genParts->size(); ++i){
            if (genParts->at(i).status() != 1){ continue; }
            if (fabs(genParts->at(i).pdgId()) != 13){ continue; }
            double pt = genParts->at(i).pt();
            double eta = fabs(genParts->at(i).eta());
            //ev_.vld_genmu_pt.push_back(pt);
            //ev_.vld_genmu_eta.push_back(eta);
            if (isHs(genParts->at(i), genParts->at(i).pdgId())){
                ev_.vld_genmu_hs_pt.push_back(pt);
                ev_.vld_genmu_hs_eta.push_back(eta);
            }else{
                ev_.vld_genmu_py8_pt.push_back(pt);
                ev_.vld_genmu_py8_eta.push_back(eta);
            }
        }
    }

    // In DYtoLL events, figure out what LL is
    unsigned int n11From1to4 = 0;
    unsigned int n13From1to4 = 0;
    unsigned int n16From1to4 = 0;
    for (size_t i=0; i<genParts->size(); ++i){
        // Count e's (11), mu's (13) and tau neutrinos (16); these are all
        // status one particles
        if (fabs(genParts->at(i).pdgId()) != 11 && fabs(genParts->at(i).pdgId()) != 13 && fabs(genParts->at(i).pdgId()) != 16){ continue; }
        // Only consider particles with a pT of at least 2 GeV
        if (genParts->at(i).pt() < 2){ continue; }
        // Check what the mother ID of that particle is, to figure out if it is
        // from the hard scattering event
        const reco::Candidate* mom = genParts->at(i).mother(0);
        int momId = mom->pdgId();
        while (mom->numberOfMothers() != 0){
            momId = mom->pdgId();
            mom = mom->mother(0);
        }
        // If second to last mother is between 1 to 4, then it is from the hard
        // scattering (at least that's my empiric assumption)
        if (fabs(momId) >= 1 && fabs(momId) <= 4){
            if (fabs(genParts->at(i).pdgId()) == 11){
                n11From1to4++;
            }else if (fabs(genParts->at(i).pdgId()) == 13){
                n13From1to4++;
            }else if (fabs(genParts->at(i).pdgId()) == 16){
                n16From1to4++;
            }
        }
    }
    // If there's at least one tau neutrino from hard scattering, it's a
    // Z --> tautau event
    // If there are at least two electrons from hard scattering or if there is
    // exactly one electron and zero muons, it's a
    // Z --> ee event
    // If there are at least two muons from hard scattering or if there is
    // exactly one muon and zero electrons, it's a
    // Z --> mumu event
    // All other events are unknown;
    if (n16From1to4 >= 1){
        ev_.ZtoLL = 15;
    }else if (n11From1to4 >= 2 || (n11From1to4 == 1 && n13From1to4 == 0)){
        ev_.ZtoLL = 11;
    }else if (n13From1to4 >= 2 || (n13From1to4 == 1 && n11From1to4 == 0)){
        ev_.ZtoLL = 13;
    }else{
        ev_.ZtoLL = 999;
    }

    // Check if there is a "crazy muon"
    // Definition of a "crazy muon":
    // pT above x, but not passing tight requirement
    ev_.crazyMuon50 = ev_.crazyMuon200 = ev_.crazyMuon500 = false;
    for (size_t i=0; i<muons->size(); ++i){
        if (isTightMuon(muons->at(i), primaryVertex, iSetup)){ continue; }
        if (muons->at(i).pt() > 50){
            ev_.crazyMuon50 = true;
        }
        if (muons->at(i).pt() > 200){
            ev_.crazyMuon200 = true;
        }
        if (muons->at(i).pt() > 500){
            ev_.crazyMuon500 = true;
        }
    }

    // Electrons
    for (size_t i=0; i<elecs->size(); ++i) {

        // Only select good electrons
        if (!isGoodElecSOS(elecs->at(i), primaryVertex)){ continue; }

        // Only select electrons above certain pT
        if (elecs->at(i).pt() < ev_.el_pt_lo){ continue; }

        ev_.nLep++;
        ev_.nEl++;
        if (elecs->at(i).pt() < ev_.el_pt_hi){
            ev_.nSoftLep++;
            ev_.nSoftEl++;
        }

        // Fill electron variables
        ev_.el_pt.push_back(elecs->at(i).pt());
        ev_.el_eta.push_back(elecs->at(i).eta());
        ev_.el_phi.push_back(elecs->at(i).phi());
        ev_.el_q.push_back(elecs->at(i).charge());
    }

    // Muons
    for (size_t i=0; i<muons->size(); ++i){

        // Only select muons above certain pT
        if (muons->at(i).pt() < ev_.mu_pt_lo){ continue; }

        // Fill muon without isolation, if not yet done so
        if (isTightMuon(muons->at(i), primaryVertex, iSetup)){
            ev_.mu_woIso_pt.push_back(muons->at(i).pt());
            ev_.mu_woIso_eta.push_back(muons->at(i).eta());
            ev_.mu_woIso_phi.push_back(muons->at(i).phi());
        }

        // Only select good muons
        if (!isGoodMuonSOS(muons->at(i), primaryVertex, iSetup)){ continue; }

        ev_.nLep++;
        ev_.nMu++;
        if (muons->at(i).pt() < ev_.mu_pt_hi){
            ev_.nSoftLep++;
            ev_.nSoftMu++;
        }

        // Skip the rest if we already have two good muons
        //if (ev_.mu2_pt.size() != 0){ continue; }

        // Check where the muon comes from
        int mother = -1;
        for (size_t j=0; j<genParts->size(); ++j){

            // Only consider muons
            if (fabs(genParts->at(j).pdgId()) != 13){ continue; }

            if (!isMatched(genParts->at(j), muons->at(i))){ continue; }

            // Loop over all mothers to find first non muon
            const reco::Candidate* mom = genParts->at(j).mother(0);
            while (mom->numberOfMothers() != 0){
                //printf("Mom: %d\n", mom->pdgId());
                if (fabs(mom->pdgId()) != 13){
                    //printf("MomFinal: %d\n", mom->pdgId());
                    mother = mom->pdgId();
                    break;
                }
                mom = mom->mother(0);
            }

            //printf("GEN  Muon: pt: %f, eta: %f, mother: %d\n", genParts->at(j).pt(), genParts->at(j).eta(), mother);
            //printf("RECO Muon: pt: %f, eta: %f, mother: %d\n", muons->at(i).pt(), muons->at(i).eta(), mother);
        }
        // Fill muon variables
        ev_.mu_pt.push_back(muons->at(i).pt());
        ev_.mu_eta.push_back(muons->at(i).eta());
        ev_.mu_phi.push_back(muons->at(i).phi());
        ev_.mu_q.push_back(muons->at(i).charge());
        ev_.mu_mother.push_back(23);
        ev_.mu_directmother.push_back(mother);
        ev_.mu_matched.push_back(true);
        ev_.mu_st20to30.push_back(true);
    }

    // Truth muons
    for (size_t i=0; i<genParts->size(); ++i){
        if (genParts->at(i).status() != 1){ continue; }
        if (fabs(genParts->at(i).pdgId()) != 13){ continue; }
        if (genParts->at(i).pt() < 10.){ continue; }
        ev_.mu_pt_truth.push_back(genParts->at(i).pt());
        ev_.mu_eta_truth.push_back(genParts->at(i).eta());
        ev_.mu_phi_truth.push_back(genParts->at(i).phi());
    }

    //// Fill leptons
    //// Put pT and eta into vector of vector for sorting
    //std::vector<std::vector<double>> lepvec;
    //if (ev_.el1_pt.size() != 0){
    //    lepvec.push_back({ev_.el1_pt.at(0), ev_.el1_eta.at(0), ev_.el1_phi.at(0), ev_.mass_el});
    //}
    //if (ev_.el2_pt.size() != 0){
    //    lepvec.push_back({ev_.el2_pt.at(0), ev_.el2_eta.at(0), ev_.el2_phi.at(0), ev_.mass_el});
    //}
    //if (ev_.mu1_pt.size() != 0){
    //    lepvec.push_back({ev_.mu1_pt.at(0), ev_.mu1_eta.at(0), ev_.mu1_phi.at(0), ev_.mass_mu});
    //}
    //if (ev_.mu2_pt.size() != 0){
    //    lepvec.push_back({ev_.mu2_pt.at(0), ev_.mu2_eta.at(0), ev_.mu2_phi.at(0), ev_.mass_mu});
    //}
    //// By definition, this sorts by the first element of the vector (in this case pT)
    //if (lepvec.size() > 1){
    //    std::sort(begin(lepvec), end(lepvec));
    //    std::reverse(begin(lepvec), end(lepvec));
    //}
    //if (lepvec.size() >= 1){
    //    ev_.lep1_pt.push_back(lepvec[0][0]);
    //    ev_.lep1_eta.push_back(lepvec[0][1]);
    //    ev_.lep1_phi.push_back(lepvec[0][2]);
    //    ev_.lep1_mass.push_back(lepvec[0][3]);
    //}
    //if (lepvec.size() >= 2){
    //    ev_.lep2_pt.push_back(lepvec[1][0]);
    //    ev_.lep2_eta.push_back(lepvec[1][1]);
    //    ev_.lep2_phi.push_back(lepvec[1][2]);
    //    ev_.lep2_mass.push_back(lepvec[1][3]);
    //}
    //lepvec.clear();

    //// Fill pt of two leptons
    //if (ev_.lep2_pt.size() != 0){
    //    TLorentzVector l1, l2;
    //    l1.SetPtEtaPhiM(ev_.lep1_pt[0], ev_.lep1_eta[0], ev_.lep1_phi[0], ev_.lep1_mass[0]);
    //    l2.SetPtEtaPhiM(ev_.lep2_pt[0], ev_.lep2_eta[0], ev_.lep2_phi[0], ev_.lep2_mass[0]);
    //    ev_.pt2l.push_back((l1 + l2).Pt());
    //}

    // MET
    if (mets->size() > 0) {
        ev_.met = mets->at(0).pt();
        ev_.met_eta = mets->at(0).eta();
        ev_.met_phi = mets->at(0).phi();
        ev_.genmet = mets->at(0).genMET()->pt();
        ev_.genmet_eta = mets->at(0).genMET()->eta();
        ev_.genmet_phi = mets->at(0).genMET()->phi();
    }

    //// Transverse mass
    //if (ev_.lep1_pt.size() != 0){
    //    ev_.mt1.push_back(std::sqrt(2*ev_.lep1_pt.at(0)*ev_.met*(1-std::cos(ev_.lep1_phi.at(0)-mets->at(0).phi()))));
    //}
    //if (ev_.lep2_pt.size() != 0){
    //    ev_.mt2.push_back(std::sqrt(2*ev_.lep2_pt.at(0)*ev_.met*(1-std::cos(ev_.lep2_phi.at(0)-mets->at(0).phi()))));
    //}

    //// mtautau
    //if (ev_.lep2_pt.size() != 0){
    //    ev_.mtautau.push_back(getMTauTau(mets->at(0).pt(), mets->at(0).phi(),
    //                ev_.lep1_pt.at(0), ev_.lep1_eta.at(0), ev_.lep1_phi.at(0),
    //                ev_.lep2_pt.at(0), ev_.lep2_eta.at(0), ev_.lep2_phi.at(0)));
    //}

    // Is a same flavour opposite sign lepton pair present?
    for (size_t i=0; i<elecs->size(); ++i){
        // Only select good electrons
        if (!isGoodElecSOS(elecs->at(i), primaryVertex)){ continue; }

        for (size_t j=i+1; j<elecs->size(); ++j){
            // Only select good electrons
            if (!isGoodElecSOS(elecs->at(j), primaryVertex)){ continue; }

            // is SFOS?
            if (elecs->at(i).charge()*elecs->at(j).charge() > 0){ continue; }

            // Are electrons soft?
            if (elecs->at(i).pt() > ev_.el_pt_hi || elecs->at(i).pt() < ev_.el_pt_lo){ continue; }
            if (elecs->at(j).pt() > ev_.el_pt_hi || elecs->at(j).pt() < ev_.el_pt_lo){ continue; }

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
        if (!isGoodMuonSOS(muons->at(i), primaryVertex, iSetup)){ continue; }

        for (size_t j=i+1; j<muons->size(); ++j){
            // Only select good muons
            if (!isGoodMuonSOS(muons->at(j), primaryVertex, iSetup)){ continue; }

            // is SFOS?
            if (muons->at(i).charge()*muons->at(j).charge() > 0){ continue; }

            // Are muons soft?
            if (muons->at(i).pt() > ev_.mu_pt_hi || muons->at(i).pt() < ev_.mu_pt_lo){ continue; }
            if (muons->at(j).pt() > ev_.mu_pt_hi || muons->at(j).pt() < ev_.mu_pt_lo){ continue; }

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
    for (size_t i=0; i<jets->size(); ++i){
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

        ev_.jet_pt.push_back(jets->at(i).pt());
        ev_.jet_eta.push_back(jets->at(i).eta());
        ev_.jet_phi.push_back(jets->at(i).phi());
        ev_.jet_mass.push_back(jets->at(i).mass());

        if (jets->at(i).pt() > 25.){
            ev_.ht25 += jets->at(i).pt();
            ev_.nJet25++;
        }
        if (jets->at(i).pt() > 40.){
            ev_.ht40 += jets->at(i).pt();
            ev_.nJet40++;
        }
        if (jets->at(i).pt() > 60.){
            ev_.ht60 += jets->at(i).pt();
            ev_.nJet60++;
        }
        if (jets->at(i).pt() > 100.){
            ev_.ht100 += jets->at(i).pt();
            ev_.nJet100++;
        }
        if (jets->at(i).pt() > 150.){
            ev_.ht150 += jets->at(i).pt();
            ev_.nJet150++;
        }

        double deepcsv = jets->at(i).bDiscriminator("pfDeepCSVJetTags:probb") +
            jets->at(i).bDiscriminator("pfDeepCSVJetTags:probbb");
        bool isMediumDeepCSV = deepcsv > deepThres_[1];
        if (!isMediumDeepCSV){ continue; }

        ev_.nBJet++;
    }

    // GenJets
    for (size_t i=0; i<genJets->size(); ++i){
        if (genJets->at(i).pt() > 25.){
            ev_.genht25 += genJets->at(i).pt();
        }
        if (genJets->at(i).pt() > 40.){
            ev_.genht40 += genJets->at(i).pt();
        }
    }

    // Fill poor man's MET
    TLorentzVector mlt4, mht4v25, mht4v40, mhlt4v25, mhlt4v40;
    for (size_t i=0; i<muons->size(); ++i) {
        if (!isGoodMuonSOS(muons->at(i), primaryVertex, iSetup)){ continue; }
        if (muons->at(i).pt() < ev_.mu_pt_lo){ continue; }
        TLorentzVector m4;
        m4.SetPtEtaPhiM(muons->at(i).pt(), muons->at(i).eta(), muons->at(i).phi(), ev_.mass_mu);
        mlt4 += m4;
        mhlt4v25 += m4;
        mhlt4v40 += m4;
    }
    for (size_t i=0; i<elecs->size(); ++i) {
        if (!isGoodElecSOS(elecs->at(i), primaryVertex)){ continue; }
        if (elecs->at(i).pt() < ev_.el_pt_lo){ continue; }
        TLorentzVector e4;
        e4.SetPtEtaPhiM(elecs->at(i).pt(), elecs->at(i).eta(), elecs->at(i).phi(), ev_.mass_el);
        mlt4 += e4;
        mhlt4v25 += e4;
        mhlt4v40 += e4;
    }
    for (size_t i=0; i<jets->size(); ++i){
        if (!isGoodJetSOS(jets->at(i))){ continue; }
        if (jets->at(i).pt() < ev_.jet_pt_lo) { continue; }
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
        TLorentzVector j4;
        j4.SetPtEtaPhiM(jets->at(i).pt(), jets->at(i).eta(), jets->at(i).phi(), jets->at(i).mass());
        if (jets->at(i).pt() > 25.){ mht4v25 +=j4; mhlt4v25 += j4; }
        if (jets->at(i).pt() > 40.){ mht4v40 +=j4; mhlt4v40 += j4; }
    }
    ev_.mlt = mlt4.Pt();
    ev_.mht25 = mht4v25.Pt();
    ev_.mht40 = mht4v40.Pt();
    ev_.mhlt25 = mhlt4v25.Pt();
    ev_.mhlt40 = mhlt4v40.Pt();

    bool evtFound = false;
    // Secondary event by event comparison: Choose events in Delphes and find them here
    if (ev_.event_by_event_comparison_secondary){
        for (size_t i=0; i<genParts->size(); ++i){
            if (genParts->at(i).status() != 1){ continue; }
            if (fabs(genParts->at(i).pdgId()) != 13){ continue; }
            // It's a final state muon!
            double mupt = genParts->at(i).pt();
            double mueta = genParts->at(i).eta();
            double muphi = genParts->at(i).phi();
            if ((fabs(mupt-6.882438) < 1.e-3) && (fabs(mueta-1.476477)) < 1.e-3 && (fabs(muphi-0.460743) < 1.e-3)){ evtFound = true; }
            if ((fabs(mupt-6.789154) < 1.e-3) && (fabs(mueta-1.970843)) < 1.e-3 && (fabs(muphi+2.242026) < 1.e-3)){ evtFound = true; }
            if ((fabs(mupt-8.218305) < 1.e-3) && (fabs(mueta-2.705792)) < 1.e-3 && (fabs(muphi-0.584903) < 1.e-3)){ evtFound = true; }
            std::cout << "Foo02: pT: " << mupt << "; eta: " << mueta << "; phi: " << muphi << std::endl;
        }
    }
    if (evtFound){
        std::cout << "Foo01" << std::endl;
    }

    // Primary event by event comparison: Choose events here and find them in Delphes
    if (ev_.event_by_event_comparison_primary || evtFound){

        //// Find specific events to compare with Delphes ...
        //if (iEvent.id().run() != 1 || iEvent.luminosityBlock() != 680){ return; }
        //unsigned long long int e = iEvent.id().event();
        //if (e!=964815 && e!=965098 && e!=965115 && e!=965602){ return; }

        // ... or find event with kinematic cuts
        //if (mets->at(0).pt() < 300.){ return; }

        // Found event
        printf("Event by event comparison. Compare event:\n");
        printf("Run Number: %d; Lumi Section: %d; Event Number: %lld\n",
                iEvent.id().run(), iEvent.luminosityBlock(), iEvent.id().event());
        printf("Has crazy muon: %d\n", ev_.crazyMuon200);

        // Truth objects
        printf("\n%20s\n", "Truth objects");
        // Truth electrons
        for (size_t i=0; i<genParts->size(); ++i){
            if (fabs(genParts->at(i).pdgId()) != 11){ continue; }
            pppWpidWstatus("Truth electrons", i, genParts->size(), genParts->at(i));
        }
        // Truth muons
        for (size_t i=0; i<genParts->size(); ++i){
            if (fabs(genParts->at(i).pdgId()) != 13){ continue; }
            pppWpidWstatus("Truth muons", i, genParts->size(), genParts->at(i));
        }
        // Truth particles (all)
        for (size_t i=0; i<genParts->size(); ++i){
            pppWpidWstatus("Truth particles", i, genParts->size(), genParts->at(i));
        }
        // Truth jets
        for (size_t i=0; i<genJets->size(); ++i){
            ppp("Truth jets", i, genJets->size(), genJets->at(i));
        }
        // Truth MET
        printf("%20s: Idx: %3d/%3d; ID: %8s; Status: %3s; pt: %8.3f; eta: %6.3f; phi: %6.3f\n",
                "Truth MET", 1, 1, "-1", "-1", mets->at(0).genMET()->pt(), mets->at(0).genMET()->eta(), mets->at(0).genMET()->phi());

        // Reco objects
        printf("\n%20s\n", "Reco objects");
        // Reco electrons
        for (size_t i=0; i<elecs->size(); ++i){
            pppWisoWpassid("Reco electrons", i, elecs->size(), elecs->at(i), elecs->at(i).charge()>0 ? -11: 11, 1, isGoodElecSOS(elecs->at(i), primaryVertex));
        }
        // Reco muons
        for (size_t i=0; i<muons->size(); ++i){
            pppWisoWpassid("Reco muons", i, muons->size(), muons->at(i), muons->at(i).charge()>0 ? -13: 13, 1, isGoodMuonSOS(muons->at(i), primaryVertex, iSetup));
        }
        // Reco jets
        for (size_t i=0; i<jets->size(); ++i){
            ppp("Reco jets", i, jets->size(), jets->at(i));
        }
        // Reco MET
        printf("%20s: Idx: %3d/%3d; ID: %8s; Status: %3s; pt: %8.3f; eta: %6.3f; phi: %6.3f\n",
                "Reco MET", 1, 1, "-", "-", mets->at(0).pt(), mets->at(0).eta(), mets->at(0).phi());


        ///// MAYBE YOU WANT TO ADD IF THE PARTICLES PASS ID SELECTION

        //for (size_t i=0; i<muons->size(); ++i) {
        //    const pat::Muon& m = muons->at(i);
        //    printf("Foo02: pT: %f; eta: %f; phi: %f; goodMuon: %d; loose: %d; medium: %d; tight: %d\n",
        //            m.pt(), m.eta(), m.phi(), isGoodMuonSOS(m, primaryVertex, iSetup), isLooseMuon(m, iSetup),
        //            isMediumMuon(m, primaryVertex), isTightMuon(m, primaryVertex, iSetup));
        //}
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

    //// Photons

    //ev_.nlp = 0;
    //ev_.ntp = 0;

    //for (size_t i = 0; i < photons->size(); i++) {
    //    if (photons->at(i).pt() < 10.) continue;
    //    if (fabs(photons->at(i).eta()) > 3.) continue;

    //    float mvaValue = photons->at(i).userFloat("mvaValue");
    //    bool isEB = photons->at(i).isEB();

    //    bool isLoose = 0;
    //    bool isTight = 0;

    //    if( isEB )
    //    {
    //        isLoose = (mvaValue > 0.00);
    //        isTight = (mvaValue > 0.56);
    //    }
    //    else
    //    {
    //        isLoose = (mvaValue > 0.20);
    //        isTight = (mvaValue > 0.68);
    //    }

    //    if (!isLoose) continue;

    //    ev_.lp_pt[ev_.nlp]     = photons->at(i).pt();
    //    ev_.lp_phi[ev_.nlp]    = photons->at(i).phi();
    //    ev_.lp_eta[ev_.nlp]    = photons->at(i).eta();
    //    ev_.lp_nrj[ev_.nlp]    = photons->at(i).energy();
    //    ev_.lp_g[ev_.nlp] = -1;
    //    for (int ig = 0; ig < ev_.ngp; ig++) {
    //        if (reco::deltaR(ev_.gp_eta[ig],ev_.gp_phi[ig],ev_.lp_eta[ev_.nlp],ev_.lp_phi[ev_.nlp]) > 0.4) continue;
    //        ev_.lp_g[ev_.nlp]    = ig;
    //    }
    //    ev_.nlp++;

    //    if (!isTight) continue;

    //    ev_.tp_pt[ev_.ntp]     = photons->at(i).pt();
    //    ev_.tp_phi[ev_.ntp]    = photons->at(i).phi();
    //    ev_.tp_eta[ev_.ntp]    = photons->at(i).eta();
    //    ev_.tp_nrj[ev_.ntp]    = photons->at(i).energy();
    //    ev_.tp_g[ev_.ntp] = -1;
    //    for (int ig = 0; ig < ev_.ngp; ig++) {
    //        if (reco::deltaR(ev_.gp_eta[ig],ev_.gp_phi[ig],ev_.tp_eta[ev_.ntp],ev_.tp_phi[ev_.ntp]) > 0.4) continue;
    //        ev_.tp_g[ev_.ntp]    = ig;
    //    }
    //    ev_.ntp++;
    //}

}

// ------------ method called for each event  ------------
    void
MiniFromPat::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

    // Clear vectors
    ev_.el_pt.clear();
    ev_.el_eta.clear();
    ev_.el_phi.clear();
    ev_.el_q.clear();
    ev_.mu_pt.clear();
    ev_.mu_eta.clear();
    ev_.mu_phi.clear();
    ev_.mu_q.clear();
    ev_.mu_mother.clear();
    ev_.mu_directmother.clear();
    ev_.mu_matched.clear();
    ev_.mu_st20to30.clear();
    ev_.mu_woIso_pt.clear();
    ev_.mu_woIso_eta.clear();
    ev_.mu_woIso_phi.clear();
    ev_.mu_pt_truth.clear();
    ev_.mu_eta_truth.clear();
    ev_.mu_phi_truth.clear();
    //ev_.lep1_pt.clear();
    //ev_.lep1_eta.clear();
    //ev_.lep1_phi.clear();
    //ev_.lep1_mass.clear();
    //ev_.lep2_pt.clear();
    //ev_.lep2_eta.clear();
    //ev_.lep2_phi.clear();
    //ev_.lep2_mass.clear();
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
    ev_.jet_pt.clear();
    ev_.jet_eta.clear();
    ev_.jet_phi.clear();
    ev_.jet_mass.clear();
    //ev_.jet_pt_truth.clear();
    //ev_.jet_eta_truth.clear();
    //ev_.jet_phi_truth.clear();
    //ev_.jet_mass_truth.clear();
    ev_.mllMin.clear();
    ev_.mllMax.clear();
    ev_.mt1.clear();
    ev_.mt2.clear();
    ev_.pt2l.clear();
    ev_.mtautau.clear();

    //ev_.vld_el_iso_abs.clear();
    //ev_.vld_el_good_abs.clear();
    //ev_.vld_el_iso_rel.clear();
    //ev_.vld_el_good_rel.clear();
    //ev_.vld_el_dxy.clear();
    //ev_.vld_el_iso_dxy.clear();
    //ev_.vld_el_dz.clear();
    //ev_.vld_el_iso_dz.clear();
    //ev_.vld_mu_iso_abs.clear();
    //ev_.vld_mu_good_abs.clear();
    //ev_.vld_mu_iso_rel.clear();
    //ev_.vld_mu_good_rel.clear();
    //ev_.vld_mu_dxy.clear();
    //ev_.vld_mu_iso_dxy.clear();
    //ev_.vld_mu_dz.clear();
    //ev_.vld_mu_iso_dz.clear();

    //ev_.vld_el_pt.clear();
    //ev_.vld_el_iso_pt.clear();
    //ev_.vld_el_good_pt.clear();
    //ev_.vld_el_is_tight.clear();
    //ev_.vld_el_hs_pt.clear();
    //ev_.vld_el_iso_hs_pt.clear();
    //ev_.vld_el_good_hs_pt.clear();
    //ev_.vld_el_py8_pt.clear();
    //ev_.vld_el_iso_py8_pt.clear();
    //ev_.vld_el_good_py8_pt.clear();
    //ev_.vld_el_others_pt.clear();
    //ev_.vld_el_iso_others_pt.clear();
    //ev_.vld_el_good_others_pt.clear();
    //ev_.vld_el_eta.clear();
    //ev_.vld_el_iso_eta.clear();
    //ev_.vld_el_good_eta.clear();
    //ev_.vld_el_hs_eta.clear();
    //ev_.vld_el_iso_hs_eta.clear();
    //ev_.vld_el_good_hs_eta.clear();
    //ev_.vld_el_py8_eta.clear();
    //ev_.vld_el_iso_py8_eta.clear();
    //ev_.vld_el_good_py8_eta.clear();
    //ev_.vld_el_others_eta.clear();
    //ev_.vld_el_iso_others_eta.clear();
    //ev_.vld_el_good_others_eta.clear();

    //ev_.vld_mu_pt.clear();
    //ev_.vld_mu_iso_pt.clear();
    //ev_.vld_mu_good_pt.clear();
    //ev_.vld_mu_is_tight.clear();
    ev_.vld_mu_hs_pt.clear();
    ev_.vld_mu_iso_hs_pt.clear();
    ev_.vld_mu_good_hs_pt.clear();
    ev_.vld_mu_py8_pt.clear();
    ev_.vld_mu_iso_py8_pt.clear();
    ev_.vld_mu_good_py8_pt.clear();
    ev_.vld_mu_others_pt.clear();
    ev_.vld_mu_iso_others_pt.clear();
    ev_.vld_mu_good_others_pt.clear();
    //ev_.vld_mu_eta.clear();
    //ev_.vld_mu_iso_eta.clear();
    //ev_.vld_mu_good_eta.clear();
    ev_.vld_mu_hs_eta.clear();
    ev_.vld_mu_iso_hs_eta.clear();
    ev_.vld_mu_good_hs_eta.clear();
    ev_.vld_mu_py8_eta.clear();
    ev_.vld_mu_iso_py8_eta.clear();
    ev_.vld_mu_good_py8_eta.clear();
    ev_.vld_mu_others_eta.clear();
    ev_.vld_mu_iso_others_eta.clear();
    ev_.vld_mu_good_others_eta.clear();

    //ev_.vld_genel_pt.clear();
    //ev_.vld_genel_hs_pt.clear();
    //ev_.vld_genel_py8_pt.clear();
    //ev_.vld_genel_eta.clear();
    //ev_.vld_genel_hs_eta.clear();
    //ev_.vld_genel_py8_eta.clear();
    //ev_.vld_genmu_pt.clear();
    ev_.vld_genmu_hs_pt.clear();
    ev_.vld_genmu_py8_pt.clear();
    //ev_.vld_genmu_eta.clear();
    ev_.vld_genmu_hs_eta.clear();
    ev_.vld_genmu_py8_eta.clear();

    //ev_.vld_el_absiso.clear();
    //ev_.vld_el_iso_absiso.clear();
    //ev_.vld_el_good_absiso.clear();
    //ev_.vld_el_hs_absiso.clear();
    //ev_.vld_el_iso_hs_absiso.clear();
    //ev_.vld_el_good_hs_absiso.clear();
    //ev_.vld_el_py8_absiso.clear();
    //ev_.vld_el_iso_py8_absiso.clear();
    //ev_.vld_el_good_py8_absiso.clear();
    //ev_.vld_el_others_absiso.clear();
    //ev_.vld_el_iso_others_absiso.clear();
    //ev_.vld_el_good_others_absiso.clear();
    //ev_.vld_el_reliso.clear();
    //ev_.vld_el_iso_reliso.clear();
    //ev_.vld_el_good_reliso.clear();
    //ev_.vld_el_hs_reliso.clear();
    //ev_.vld_el_iso_hs_reliso.clear();
    //ev_.vld_el_good_hs_reliso.clear();
    //ev_.vld_el_py8_reliso.clear();
    //ev_.vld_el_iso_py8_reliso.clear();
    //ev_.vld_el_good_py8_reliso.clear();
    //ev_.vld_el_others_reliso.clear();
    //ev_.vld_el_iso_others_reliso.clear();
    //ev_.vld_el_good_others_reliso.clear();
    //ev_.vld_el_dxy.clear();
    //ev_.vld_el_iso_dxy.clear();
    //ev_.vld_el_good_dxy.clear();
    //ev_.vld_el_hs_dxy.clear();
    //ev_.vld_el_iso_hs_dxy.clear();
    //ev_.vld_el_good_hs_dxy.clear();
    //ev_.vld_el_py8_dxy.clear();
    //ev_.vld_el_iso_py8_dxy.clear();
    //ev_.vld_el_good_py8_dxy.clear();
    //ev_.vld_el_others_dxy.clear();
    //ev_.vld_el_iso_others_dxy.clear();
    //ev_.vld_el_good_others_dxy.clear();
    //ev_.vld_el_dz.clear();
    //ev_.vld_el_iso_dz.clear();
    //ev_.vld_el_good_dz.clear();
    //ev_.vld_el_hs_dz.clear();
    //ev_.vld_el_iso_hs_dz.clear();
    //ev_.vld_el_good_hs_dz.clear();
    //ev_.vld_el_py8_dz.clear();
    //ev_.vld_el_iso_py8_dz.clear();
    //ev_.vld_el_good_py8_dz.clear();
    //ev_.vld_el_others_dz.clear();
    //ev_.vld_el_iso_others_dz.clear();
    //ev_.vld_el_good_others_dz.clear();

    //ev_.vld_mu_absiso.clear();
    //ev_.vld_mu_iso_absiso.clear();
    //ev_.vld_mu_good_absiso.clear();
    //ev_.vld_mu_hs_absiso.clear();
    //ev_.vld_mu_iso_hs_absiso.clear();
    //ev_.vld_mu_good_hs_absiso.clear();
    //ev_.vld_mu_py8_absiso.clear();
    //ev_.vld_mu_iso_py8_absiso.clear();
    //ev_.vld_mu_good_py8_absiso.clear();
    //ev_.vld_mu_others_absiso.clear();
    //ev_.vld_mu_iso_others_absiso.clear();
    //ev_.vld_mu_good_others_absiso.clear();
    //ev_.vld_mu_reliso.clear();
    //ev_.vld_mu_iso_reliso.clear();
    //ev_.vld_mu_good_reliso.clear();
    ev_.vld_mu_hs_reliso.clear();
    ev_.vld_mu_iso_hs_reliso.clear();
    ev_.vld_mu_good_hs_reliso.clear();
    ev_.vld_mu_py8_reliso.clear();
    ev_.vld_mu_iso_py8_reliso.clear();
    ev_.vld_mu_good_py8_reliso.clear();
    ev_.vld_mu_others_reliso.clear();
    ev_.vld_mu_iso_others_reliso.clear();
    ev_.vld_mu_good_others_reliso.clear();
    //ev_.vld_mu_dxy.clear();
    //ev_.vld_mu_iso_dxy.clear();
    //ev_.vld_mu_good_dxy.clear();
    //ev_.vld_mu_hs_dxy.clear();
    //ev_.vld_mu_iso_hs_dxy.clear();
    //ev_.vld_mu_good_hs_dxy.clear();
    //ev_.vld_mu_py8_dxy.clear();
    //ev_.vld_mu_iso_py8_dxy.clear();
    //ev_.vld_mu_good_py8_dxy.clear();
    //ev_.vld_mu_others_dxy.clear();
    //ev_.vld_mu_iso_others_dxy.clear();
    //ev_.vld_mu_good_others_dxy.clear();
    //ev_.vld_mu_dz.clear();
    //ev_.vld_mu_iso_dz.clear();
    //ev_.vld_mu_good_dz.clear();
    //ev_.vld_mu_hs_dz.clear();
    //ev_.vld_mu_iso_hs_dz.clear();
    //ev_.vld_mu_good_hs_dz.clear();
    //ev_.vld_mu_py8_dz.clear();
    //ev_.vld_mu_iso_py8_dz.clear();
    //ev_.vld_mu_good_py8_dz.clear();
    //ev_.vld_mu_others_dz.clear();
    //ev_.vld_mu_iso_others_dz.clear();
    //ev_.vld_mu_good_others_dz.clear();

    ev_.nLep = ev_.nEl = ev_.nMu = 0;
    ev_.nSoftLep = ev_.nSoftEl = ev_.nSoftMu = 0;
    ev_.nJet25 = ev_.nJet40 = ev_.nJet60 = ev_.nJet100 = ev_.nJet150 = ev_.nBJet = 0;
    ev_.met = ev_.met_eta = ev_.met_phi = ev_.genmet = ev_.genmet_eta = ev_.genmet_phi = 0.;
    ev_.ht25 = ev_.ht40 = ev_.ht60 = ev_.ht100 = ev_.ht150 = 0.;
    ev_.genht25 = ev_.genht40 = 0.;

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

bool MiniFromPat::isLooseElec(const pat::Electron & patEl){

    if (patEl.pt() < 5.){ return false; }
    if (fabs(patEl.eta()) > 3.){ return false; }

    float mvaValue = patEl.userFloat("mvaValue");
    bool isEB = patEl.isEB();

    bool isLoose = 0;

    if( isEB ) {
        if (patEl.pt() < 20.) {
            isLoose = (mvaValue > -0.661);
        }
        else {
            isLoose = (mvaValue > -0.797);
        }
    }
    else {
        if (not (patEl.userFloat("hgcElectronID:ecEnergy") > 0)){ return false; }
        if (not (patEl.userFloat("hgcElectronID:sigmaUU") > 0)){ return false; }
        if (not (patEl.fbrem() > -1)){ return false; }
        if (not (patEl.userFloat("hgcElectronID:measuredDepth") < 40)){ return false; }
        if (not (patEl.userFloat("hgcElectronID:nLayers") > 20)){ return false; }
        if (patEl.pt() < 20.) {
            isLoose = (mvaValue > -0.320);
        }
        else {
            isLoose = (mvaValue > -0.919);
        }
    }
    return isLoose;
}

bool MiniFromPat::isMediumElec(const pat::Electron & patEl){

    if (patEl.pt() < 5.){ return false; }
    if (fabs(patEl.eta()) > 3.){ return false; }

    float mvaValue = patEl.userFloat("mvaValue");
    bool isEB = patEl.isEB();

    bool isMedium = 0;

    if( isEB ) {
        if (patEl.pt() < 20.) {
            isMedium = mvaValue > 0.855;
        }
        else {
            isMedium = mvaValue > 0.723;
        }
    }
    else {
        if (not (patEl.userFloat("hgcElectronID:ecEnergy") > 0)){ return false; }
        if (not (patEl.userFloat("hgcElectronID:sigmaUU") > 0)){ return false; }
        if (not (patEl.fbrem() > -1)){ return false; }
        if (not (patEl.userFloat("hgcElectronID:measuredDepth") < 40)){ return false; }
        if (not (patEl.userFloat("hgcElectronID:nLayers") > 20)){ return false; }
        if (patEl.pt() < 20.) {
            isMedium = mvaValue > 0.777;
        }
        else {
            isMedium = mvaValue > 0.591;
        }
    }
    return isMedium;
}

bool MiniFromPat::isTightElec(const pat::Electron & patEl){

    if (patEl.pt() < 5.){ return false; }
    if (fabs(patEl.eta()) > 3.){ return false; }

    float mvaValue = patEl.userFloat("mvaValue");
    bool isEB = patEl.isEB();

    bool isTight = 0;

    if( isEB ) {
        if (patEl.pt() < 20.) {
            isTight = (mvaValue > 0.986);
        }
        else {
            isTight = (mvaValue > 0.988);
        }
    }
    else {
        if (not (patEl.userFloat("hgcElectronID:ecEnergy") > 0)){ return false; }
        if (not (patEl.userFloat("hgcElectronID:sigmaUU") > 0)){ return false; }
        if (not (patEl.fbrem() > -1)){ return false; }
        if (not (patEl.userFloat("hgcElectronID:measuredDepth") < 40)){ return false; }
        if (not (patEl.userFloat("hgcElectronID:nLayers") > 20)){ return false; }
        if (patEl.pt() < 20.) {
            isTight = (mvaValue > 0.969);
        }
        else {
            isTight = (mvaValue > 0.983);
        }
    }
    return isTight;
}

bool MiniFromPat::isIsolatedElec(const pat::Electron & patEl){
    double elIso;
    if( patEl.isEB() )
        elIso = (patEl.puppiNoLeptonsChargedHadronIso() + patEl.puppiNoLeptonsNeutralHadronIso() + patEl.puppiNoLeptonsPhotonIso()) / patEl.pt();
    else
        elIso = (patEl.userFloat("hgcElectronID:caloIsoRing1") + patEl.userFloat("hgcElectronID:caloIsoRing2") + patEl.userFloat("hgcElectronID:caloIsoRing3") + patEl.userFloat("hgcElectronID:caloIsoRing4")) / patEl.energy();
    if (elIso/patEl.pt() > .5){ return false; }
    if (elIso > 5.){ return false; }
    return true;
}

bool MiniFromPat::isGoodElecSOS(const pat::Electron & patEl, reco::Vertex primaryVertex){
    // Isolation
    if (!isIsolatedElec(patEl)){ return false; }

    // ID cuts
    if (!isTightElec(patEl)){ return false; }

    float dxy=0;
    float dz=0;
    if(patEl.gsfTrack().isNonnull()){
        dxy=std::abs(patEl.gsfTrack()->dxy(primaryVertex.position()));
        dz=std::abs(patEl.gsfTrack()->dz(primaryVertex.position()));
    }
    if (sqrt(dxy*dz) > .01){ return false; }

    return true;
}

bool MiniFromPat::isLooseMuon(const pat::Muon & patMu, const edm::EventSetup& iSetup){

    if (patMu.pt() < 2.) return false;
    if (fabs(patMu.eta()) > 2.8) return false;

    // Loose ID
    double dPhiCut = std::min(std::max(1.2/patMu.p(),1.2/100),0.056);
    double dPhiBendCut = std::min(std::max(0.2/patMu.p(),0.2/100),0.0096);
    bool isLoose = (fabs(patMu.eta()) < 2.4 && muon::isLooseMuon(patMu))
        || (fabs(patMu.eta()) > 2.4 && isME0MuonSelNew(patMu, .077, dPhiCut, dPhiBendCut, iSetup));
    return isLoose;
}

bool MiniFromPat::isMediumMuon(const pat::Muon & patMu, reco::Vertex primaryVertex){

    if (patMu.pt() < 2.) return false;
    if (fabs(patMu.eta()) > 2.8) return false;

    // Medium ID -- needs to be updated
    // bool isMedium = (fabs(patMu.eta()) < 2.4 && muon::isMediumMuon(patMu)) || (fabs(patMu.eta()) > 2.4 && isME0MuonSelNew(patMu, 0.077, dPhiCut, dPhiBendCut) && ipxy && ipz && validPxlHit && highPurity);
    return false;
}

bool MiniFromPat::isTightMuon(const pat::Muon & patMu, reco::Vertex primaryVertex, const edm::EventSetup& iSetup){

    if (patMu.pt() < 2.) return false;
    if (fabs(patMu.eta()) > 2.8) return false;

    bool ipxy = false, ipz = false, validPxlHit = false, highPurity = false;
    float dz=0, dxy=0;
    if (patMu.innerTrack().isNonnull()){
        dxy=std::abs(patMu.muonBestTrack()->dxy(primaryVertex.position()));
        dz= std::abs(patMu.muonBestTrack()->dz(primaryVertex.position()));
        ipxy = dxy < 0.2;
        ipz = dz < 0.5;
        validPxlHit = patMu.innerTrack()->hitPattern().numberOfValidPixelHits() > 0;
        highPurity = patMu.innerTrack()->quality(reco::Track::highPurity);
    }

    // Tight ID
    double dPhiCut = std::min(std::max(1.2/patMu.p(),1.2/100),0.032);
    double dPhiBendCut = std::min(std::max(0.2/patMu.p(),0.2/100),0.0041);
    bool isTight = (fabs(patMu.eta()) < 2.4 && muon::isTightMuon(patMu, primaryVertex))
        || (fabs(patMu.eta()) > 2.4 && isME0MuonSelNew(patMu, 0.048, dPhiCut, dPhiBendCut, iSetup) && ipxy && ipz && validPxlHit && highPurity);
    return isTight;
}

bool MiniFromPat::isIsolatedMuon(const pat::Muon & patMu){

    // PUPPI isolation (apparently buggy)
    //double muIsolation = patMu.puppiChargedHadronIso() + patMu.puppiNeutralHadronIso() + patMu.puppiPhotonIso();
    // PF isolation (not recommended for some reason)
    //double muIsolation = patMu.pfIsolationR04().sumChargedHadronPt +
    //                         std::max(0., patMu.pfIsolationR04().sumNeutralHadronEt + patMu.pfIsolationR04().sumPhotonEt - 0.5*patMu.pfIsolationR04().sumPUPt);

    // Tracker based isolation (best what we have for now)
    //double muIsolation = patMu.isolationR03().sumPt;
    double muIsolation = patMu.trackIso();
    if (muIsolation / patMu.pt() > 0.5){ return false; }
    if (muIsolation > 5.){ return false; }
    return true;
}

bool MiniFromPat::isGoodMuonSOS(const pat::Muon & patMu, reco::Vertex primaryVertex, edm::EventSetup const& iSetup){
    // Isolation
    if (!isIsolatedMuon(patMu)){ return false; }

    // ID cuts
    if (!isTightMuon(patMu, primaryVertex, iSetup)){ return false; }

    // IP3D cuts
    double dxy = std::abs(patMu.muonBestTrack()->dxy(primaryVertex.position()));
    double dz = std::abs(patMu.muonBestTrack()->dz(primaryVertex.position()));
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

double MiniFromPat::getMTauTau(double met_pt, double met_phi, double l1_pt, double l1_eta, double l1_phi, double l2_pt, double l2_eta, double l2_phi){
    typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<double> > PtEtaPhiMVector;
    typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzM4D<double>   > PxPyPzMVector;
    PtEtaPhiMVector Met( met_pt, 0.     , met_phi , 0.   );
    PtEtaPhiMVector L1(  l1_pt , l1_eta , l1_phi  , 0.106 );
    PtEtaPhiMVector L2(  l2_pt , l2_eta , l2_phi  , 0.106 );   // 0.106 mu mass
    float A00,A01,A10,A11,  C0,C1,  X0,X1,  inv_det;     // Define A:2x2 matrix, C,X 2x1 vectors & det[A]^-1
    inv_det = 1./( L1.Px()*L2.Py() - L2.Px()*L1.Py() );
    A00 = inv_det*L2.Py();     A01 =-inv_det*L2.Px();
    A10 =-inv_det*L1.Py();     A11 = inv_det*L1.Px();
    C0  = (Met+L1+L2).Px();    C1  = (Met+L1+L2).Py();
    X0  = A00*C0 + A01*C1;     X1  = A10*C0 + A11*C1;
    PxPyPzMVector T1( L1.Px()*X0 , L1.Py()*X0 , L1.Pz()*X0 , 1.777 );    // 1.777 tau mass
    PxPyPzMVector T2( L2.Px()*X1 , L2.Py()*X1 , L2.Pz()*X1 , 1.777 );
    if(X0>0.&&X1>0.)return  (T1+T2).M();
    else            return -(T1+T2).M();
}

double MiniFromPat::DeltaR(double eta1, double eta2, double phi1, double phi2){
    double dEta = eta1 - eta2;
    double dPhi = DeltaPhi(phi1, phi2);
    return TMath::Sqrt(dEta*dEta + dPhi*dPhi);
}

double MiniFromPat::DeltaPhi(double phi1, double phi2){
    double dPhi = phi1 - phi2;
    while (dPhi  >  TMath::Pi()) dPhi -= 2*TMath::Pi();
    while (dPhi <= -TMath::Pi()) dPhi += 2*TMath::Pi();
    return fabs(dPhi);
}

template <typename T> bool MiniFromPat::isMatched(const pat::PackedGenParticle truthParticle, const T particle){
    if (fabs((truthParticle.pt() - particle.pt())/truthParticle.pt()) > ev_.truth_match_diff_pt_rel){ return false; }
    if (DeltaR(truthParticle.eta(), particle.eta(), truthParticle.phi(), particle.phi()) > ev_.truth_match_diff_dr){ return false; }
    return true;
}

template <typename T> bool MiniFromPat::matchAny(const edm::Handle<std::vector<pat::PackedGenParticle>> genParts, T particle, bool hs){
    for (size_t i=0; i<genParts->size(); ++i){
        if (genParts->at(i).pdgId() != particle.pdgId()){ continue; }
        // Only consider particles that are from hard scattering when hs is
        // true, or that are *not* from hard scattering when hs is false
        if (hs ^ isHs(genParts->at(i), particle.pdgId())){ continue; }
        if (isMatched(genParts->at(i), particle)){
            return true;
        }
    }
    return false;
}

// Is truth particle from hard scattering?
bool MiniFromPat::isHs(const pat::PackedGenParticle truthParticle, int pdgId){
    if (truthParticle.pdgId() == pdgId){
        // Loop over all mothers to find SUSY mother particle (or not)
        const reco::Candidate* mom = truthParticle.mother(0);
        if (mom->pdgId() >= 1000000){ return true; }
        while (mom->numberOfMothers() != 0){
            mom = mom->mother(0);
            if (mom->pdgId() >= 1000000){
                return true;
            }
        }
    }
    return false;
}


// Print properties of particle (ppp = print particle properties)
// ppp with PID and with status (for truth leptons)
template <typename T> void MiniFromPat::pppWpidWstatus(const char* text, const size_t idx, const size_t noParticles, const T particle, const char* addText) const {
    ppp(text, idx, noParticles, particle, particle.pdgId(), particle.status(), -1., -1., addText);
}

// ppp with isolation and boolean if ID is passed (for reco leptons)
template <typename T> void MiniFromPat::pppWisoWpassid(const char* text, const size_t idx, const size_t noParticles, const T particle, const int pid, const int status, const short int passID, const char* addText) const {
    // Isolation depends on particle type
    double iso = -1.;
    if (fabs(particle.pdgId()) == 13){ // muon
        iso = particle.trackIso();
    }
    ppp(text, idx, noParticles, particle, particle.pdgId(), particle.status(), iso, passID, addText);
}

// Generic base ppp method
template <typename T> void MiniFromPat::ppp(const char* text, const size_t idx, const size_t noParticles, const T particle, const int pid, const int status, const double iso, const short int passID, const char* addText) const {
    printf("%20s: Idx: %3lu/%3lu; ID: %8d; Status: %3d; pt: %8.3f; eta: %6.3f; phi: %6.3f; iso: %8.3f; passID: %1d; %s\n",
            text, idx, noParticles, pid, status, particle.pt(), particle.eta(), particle.phi(), iso, passID, addText);
    fflush(stdout);
}


bool MiniFromPat::isME0MuonSelNew(const reco::Muon& muon, double dEtaCut, double dPhiCut, double dPhiBendCut, edm::EventSetup const& iSetup){

    bool result = false;
    bool isME0 = muon.isME0Muon();

    if(isME0){

        double deltaEta = 999;
        double deltaPhi = 999;
        double deltaPhiBend = 999;

        if(!ME0Geometry_){
            edm::ESHandle<ME0Geometry> hGeom;
            iSetup.get<MuonGeometryRecord>().get(hGeom);
            ME0Geometry_ =( &*hGeom);
            if(!ME0Geometry_)
                return false;
        }

        const std::vector<reco::MuonChamberMatch>& chambers = muon.matches();
        for( std::vector<reco::MuonChamberMatch>::const_iterator chamber = chambers.begin(); chamber != chambers.end(); ++chamber ){

            if (chamber->detector() == 5){

                for ( std::vector<reco::MuonSegmentMatch>::const_iterator segment = chamber->me0Matches.begin(); segment != chamber->me0Matches.end(); ++segment ){

                    LocalPoint trk_loc_coord(chamber->x, chamber->y, 0);
                    LocalPoint seg_loc_coord(segment->x, segment->y, 0);
                    LocalVector trk_loc_vec(chamber->dXdZ, chamber->dYdZ, 1);
                    LocalVector seg_loc_vec(segment->dXdZ, segment->dYdZ, 1);

                    const ME0Chamber * me0chamber = ME0Geometry_->chamber(chamber->id);
                    if(!me0chamber)continue;

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

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Utilities/interface/EDMException.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Common/interface/Handle.h"
#include "CommonTools/Utils/interface/StringObjectFunction.h"
#include "CommonTools/CandAlgos/interface/ModifyObjectValueBase.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalTools.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "DataFormats/EgammaCandidates/interface/PhotonCore.h"
#include "DataFormats/EgammaCandidates/interface/Conversion.h"
#include "DataFormats/EgammaCandidates/interface/ConversionFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronCore.h"
#include "DataFormats/ParticleFlowReco/interface/GsfPFRecTrack.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Photon.h"
#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "Geometry/CaloTopology/interface/CaloSubdetectorTopology.h"
#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/Records/interface/CaloTopologyRecord.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "FWCore/Utilities/interface/isFinite.h"
#include "TLorentzVector.h"
#include "TMath.h"
#include "TTree.h"
#include <vector>
#include <map>
#include <algorithm>

using namespace edm;
using namespace std;
using namespace reco;

class ElePhoDumper : public edm::EDAnalyzer
{
      public:
         explicit ElePhoDumper(const edm::ParameterSet&);
	 ~ElePhoDumper();
  
  
      private:
	 void beginJob() override;
	 void analyze(const edm::Event&, const edm::EventSetup&) override;
         void endJob() override;
        
      // ----------additional functions-------------------
      void setTree(TTree* tree);
      void setVectors(int genParticle_size, int patElectron_size, int patPhoton_size);
      template <typename T> void setDefaultValues(std::vector<T>& vec, const T& value);
      const reco::GenParticle* getFirstMother( const reco::GenParticle* particle, int pdgId, int status ); 
      int getGenIndex( edm::Handle<std::vector<reco::GenParticle> > genParticles, const reco::GenParticle* particle ); 
      void printAllDaughters( const reco::GenParticle* particle, int depth ); 
      math::XYZTLorentzVector setP4(const reco::SuperClusterRef scRef, const reco::Vertex& vtx, double energy);
      int passPreselections(const pat::Photon* photon);
      
      // ----------collection tokens-------------------
      edm::EDGetTokenT<double> rhoToken_;
      edm::EDGetTokenT<std::vector<PileupSummaryInfo> > pileupSummaryToken_;
      edm::EDGetTokenT<reco::VertexCollection> vtxToken_; 
      edm::EDGetTokenT<std::vector<reco::GenParticle> > genToken_; 
      edm::EDGetTokenT<LHEEventProduct> lheEventsToken_;
      edm::EDGetTokenT<EcalRecHitCollection> ebRechitToken_; 
      edm::EDGetTokenT<EcalRecHitCollection> eeRechitToken_; 
      edm::EDGetTokenT<std::vector<pat::Electron> > patElectronToken_;
      edm::EDGetTokenT<std::vector<pat::Photon> > patPhotonToken_;
      
      edm::Service<TFileService> iFile;
      
      // ----------config inputs-------------------
      bool debug_;
      bool isMC_;
      bool isParticleGun_;
      std::vector<int> motherPdgId_;
      std::vector<int> motherStatus_;
      std::vector<int> pdgIdBeforeFSR_;
      std::vector<int> pdgId_;
      std::vector<int> status_;
      bool doLHEMatching_;
      double dRMaxLHEMatch_;
      double dEFracMaxLHEMatch_;
      bool ignoreTauDecays_;
      double dRMin_;
      bool savePhotons_; 
      bool saveElectrons_; 
      
      std::string egmCutBasedElectronIDVeto_;
      std::string egmCutBasedElectronIDloose_;
      std::string egmCutBasedElectronIDmedium_;
      std::string egmCutBasedElectronIDtight_; 
      std::string egmMVAElectronIDloose_;
      std::string egmMVAElectronIDmedium_;
      std::string egmMVAElectronIDtight_; 
      std::string egmMVAElectronIDlooseNoIso_;
      std::string egmMVAElectronIDmediumNoIso_;
      std::string egmMVAElectronIDtightNoIso_; 
      std::string heepElectronID_; 
      std::string egmCutBasedPhotonIDloose_;
      std::string egmCutBasedPhotonIDmedium_;
      std::string egmCutBasedPhotonIDtight_;  
      std::string egmMVAPhotonIDmedium_;
      std::string egmMVAPhotonIDtight_; 
      
      std::vector<int> genIndices; 
      std::map<int,int> matchedPhotons;
      const reco::Vertex* vtx0;
      math::XYZTLorentzVector p4_vtx0;
      TLorentzVector lvec;
      std::vector<TLorentzVector> lheParticles;
      std::vector<int> lheParticlesID;
      std::vector<int> lheParticlesStatus;
      std::map<int,double> fsrEnergyMap;
      
      // ----------histograms & trees & branches-------------------
      TTree* tree;
      
      long int eventId;
      int lumiId;
      int runId; 
      double truePU;
      double obsPU;
      int nVtx;
      double rho; 
      int genParticle_size; 
      std::vector<int> genParticle_index;
      std::vector<int> genParticle_pdgId;
      std::vector<int> genParticle_status; 
      std::vector<double> genParticle_energy;
      std::vector<double> genParticle_pt;
      std::vector<double> genParticle_eta;
      std::vector<double> genParticle_phi;
      std::vector<double> genParticle_mass;
      std::vector<int> genParticle_charge;
      std::vector<int> genParticle_isLHEMatched;
      int patElectron_size;
      std::vector<int> patElectron_index; 
      std::vector<int> patElectron_genIndex; 
      std::vector<double> patElectron_genDR; 
      std::vector<double> patElectron_genFSREnergy;
      std::vector<double> patElectron_genFSRNergy;
      std::vector<int> patElectron_classification;
      std::vector<int> patElectron_nrSatCrys;
      std::vector<int> patElectron_scNPFClusters;
      std::vector<int> patElectron_charge;
      std::vector<int> patElectron_isEB;
      std::vector<int> patElectron_isEE;
      std::vector<int> patElectron_isEBEEGap;
      std::vector<int> patElectron_isEBEtaGap;
      std::vector<int> patElectron_isEBPhiGap;
      std::vector<int> patElectron_isEEDeeGap;
      std::vector<int> patElectron_isEERingGap;
      std::vector<int> patElectron_isEcalDriven;
      std::vector<int> patElectron_isTrackerDriven;
      std::vector<int> patElectron_passConversionVeto;
      std::vector<double> patElectron_eta;
      std::vector<double> patElectron_phi;
      std::vector<double> patElectron_p;
      std::vector<double> patElectron_pt;
      std::vector<double> patElectron_pIn;
      std::vector<double> patElectron_pOut;
      std::vector<double> patElectron_trackFbrem;
      std::vector<double> patElectron_superClusterFbrem;
      std::vector<double> patElectron_energy;
      std::vector<double> patElectron_energyErr;
      std::vector<double> patElectron_ecalEnergy;
      std::vector<double> patElectron_ecalEnergyErr;
      std::vector<double> patElectron_ecalEnergyPhoReg;
      std::vector<double> patElectron_ecalEnergyPhoRegErr;
      std::vector<double> patElectron_et;
      std::vector<double> patElectron_HoE;
      std::vector<int> patElectron_scIsEB;
      std::vector<double> patElectron_scEnergy;
      std::vector<double> patElectron_scRawEnergy;
      std::vector<double> patElectron_scRawESEnergy;
      std::vector<double> patElectron_scSeedEnergy;
      std::vector<double> patElectron_scDPhiSeed;
      std::vector<double> patElectron_scDEtaSeed;
      std::vector<int> patElectron_scIEtaOrX;
      std::vector<int> patElectron_scIPhiOrY;
      std::vector<int> patElectron_scIEtaMod5;
      std::vector<int> patElectron_scIPhiMod2;
      std::vector<int> patElectron_scIEtaMod20;
      std::vector<int> patElectron_scIPhiMod20;
      std::vector<double> patElectron_scEt;
      std::vector<double> patElectron_scPhiWidth;
      std::vector<double> patElectron_scEtaWidth;
      std::vector<double> patElectron_scEoP;
      std::vector<double> patElectron_scEta;
      std::vector<double> patElectron_scPhi;
      std::vector<double> patElectron_scEta_vtx0;
      std::vector<double> patElectron_scPhi_vtx0;
      std::vector<double> patElectron_SwissCross;
      std::vector<double> patElectron_r9;
      std::vector<double> patElectron_full5x5_r9;
      std::vector<double> patElectron_full5x5_e3x3;
      std::vector<double> patElectron_full5x5_e5x5;
      std::vector<double> patElectron_full5x5_eMax;
      std::vector<double> patElectron_full5x5_e2nd;
      std::vector<double> patElectron_full5x5_eTop;
      std::vector<double> patElectron_full5x5_eBottom;
      std::vector<double> patElectron_full5x5_eLeft;
      std::vector<double> patElectron_full5x5_eRight;
      std::vector<double> patElectron_full5x5_e2x5Max;
      std::vector<double> patElectron_full5x5_e2x5Left;
      std::vector<double> patElectron_full5x5_e2x5Right;
      std::vector<double> patElectron_full5x5_e2x5Top;
      std::vector<double> patElectron_full5x5_e2x5Bottom;
      std::vector<double> patElectron_full5x5_sigmaIEtaIEta;
      std::vector<double> patElectron_full5x5_sigmaIEtaIPhi;
      std::vector<double> patElectron_full5x5_sigmaIPhiIPhi;
      std::vector<double> patElectron_egmMVAElectronIDScore;
      std::vector<double> patElectron_egmMVAElectronIDNoIsoScore;
      std::vector<int> patElectron_egmCutBasedElectronIDVeto;
      std::vector<int> patElectron_egmCutBasedElectronIDloose;
      std::vector<int> patElectron_egmCutBasedElectronIDmedium;
      std::vector<int> patElectron_egmCutBasedElectronIDtight;
      std::vector<int> patElectron_egmMVAElectronIDloose;
      std::vector<int> patElectron_egmMVAElectronIDmedium;
      std::vector<int> patElectron_egmMVAElectronIDtight;
      std::vector<int> patElectron_egmMVAElectronIDlooseNoIso;
      std::vector<int> patElectron_egmMVAElectronIDmediumNoIso;
      std::vector<int> patElectron_egmMVAElectronIDtightNoIso;
      std::vector<int> patElectron_heepElectronID;
      int patPhoton_size;
      std::vector<int> patPhoton_index; 
      std::vector<int> patPhoton_genIndex; 
      std::vector<double> patPhoton_genDR; 
      std::vector<int> patPhoton_nrSatCrys;
      std::vector<int> patPhoton_scNPFClusters;
      std::vector<int> patPhoton_isEB;
      std::vector<int> patPhoton_isEE;
      std::vector<int> patPhoton_isEBEEGap;
      std::vector<int> patPhoton_isEBEtaGap;
      std::vector<int> patPhoton_isEBPhiGap;
      std::vector<int> patPhoton_isEEDeeGap;
      std::vector<int> patPhoton_isEERingGap;
      std::vector<int> patPhoton_passElectronVeto;
      std::vector<int> patPhoton_hasPixelSeed;
      std::vector<int> patPhoton_hasConversionTracks;
      std::vector<int> patPhoton_nConversions;
      std::vector<int> patPhoton_nConversionsOneLeg;  
      std::vector<double> patPhoton_eta;
      std::vector<double> patPhoton_phi;
      std::vector<double> patPhoton_energy; 
      std::vector<double> patPhoton_energyErr;
      std::vector<double> patPhoton_ecalEnergy;
      std::vector<double> patPhoton_ecalEnergyErr;
      std::vector<double> patPhoton_et;
      std::vector<double> patPhoton_pt;
      std::vector<double> patPhoton_mt;
      std::vector<double> patPhoton_dphiMET;   
      std::vector<double> patPhoton_HoE;
      std::vector<int> patPhoton_scIsEB;  
      std::vector<double> patPhoton_scEnergy;  
      std::vector<double> patPhoton_scRawEnergy;  
      std::vector<double> patPhoton_scRawESEnergy;
      std::vector<double> patPhoton_scSeedEnergy;
      std::vector<double> patPhoton_scDPhiSeed;
      std::vector<double> patPhoton_scDEtaSeed;
      std::vector<int> patPhoton_scIEtaOrX;
      std::vector<int> patPhoton_scIPhiOrY;
      std::vector<int> patPhoton_scIEtaMod5;
      std::vector<int> patPhoton_scIPhiMod2;
      std::vector<int> patPhoton_scIEtaMod20;
      std::vector<int> patPhoton_scIPhiMod20;
      std::vector<double> patPhoton_scEt;
      std::vector<double> patPhoton_scEtaWidth;
      std::vector<double> patPhoton_scPhiWidth;    
      std::vector<double> patPhoton_scEta;
      std::vector<double> patPhoton_scPhi;
      std::vector<double> patPhoton_scEta_vtx0;
      std::vector<double> patPhoton_scPhi_vtx0;
      std::vector<double> patPhoton_SwissCross;
      std::vector<double> patPhoton_r9;
      std::vector<double> patPhoton_full5x5_r9;
      std::vector<double> patPhoton_full5x5_e3x3;
      std::vector<double> patPhoton_full5x5_e5x5;
      std::vector<double> patPhoton_full5x5_eMax;
      std::vector<double> patPhoton_full5x5_e2nd;
      std::vector<double> patPhoton_full5x5_eTop;
      std::vector<double> patPhoton_full5x5_eBottom;
      std::vector<double> patPhoton_full5x5_eLeft;
      std::vector<double> patPhoton_full5x5_eRight;
      std::vector<double> patPhoton_full5x5_e2x5Max;
      std::vector<double> patPhoton_full5x5_e2x5Left;
      std::vector<double> patPhoton_full5x5_e2x5Right;
      std::vector<double> patPhoton_full5x5_e2x5Top;
      std::vector<double> patPhoton_full5x5_e2x5Bottom;
      std::vector<double> patPhoton_full5x5_sigmaIEtaIEta;
      std::vector<double> patPhoton_full5x5_sigmaIEtaIPhi;
      std::vector<double> patPhoton_full5x5_sigmaIPhiIPhi;
      std::vector<double> patPhoton_full5x5_effSigmaRR;
      std::vector<double> patPhoton_ecalPFClusterIso;
      std::vector<double> patPhoton_hcalPFClusterIso;
      std::vector<double> patPhoton_trkSumPtHollowConeDR03;
      std::vector<double> patPhoton_trkSumPtSolidConeDR04;
      std::vector<double> patPhoton_chargedHadronIso;
      std::vector<double> patPhoton_chargedHadronWorstVtxIso;
      std::vector<double> patPhoton_egmMVAPhotonIDScore;
      std::vector<int> patPhoton_egmCutBasedPhotonIDloose;
      std::vector<int> patPhoton_egmCutBasedPhotonIDmedium;
      std::vector<int> patPhoton_egmCutBasedPhotonIDtight;
      std::vector<int> patPhoton_egmMVAPhotonIDmedium;
      std::vector<int> patPhoton_egmMVAPhotonIDtight;  
      std::vector<int> patPhoton_passPreselections;
  
};

//
// constructors and destructor
//
ElePhoDumper::ElePhoDumper(const edm::ParameterSet& iConfig)
{
   pileupSummaryToken_            = consumes<std::vector<PileupSummaryInfo> >(iConfig.getParameter<edm::InputTag>("pileupSummary"));
   vtxToken_                      = consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertexCollection"));
   rhoToken_                      = consumes<double>(iConfig.getParameter<edm::InputTag>("rhoCollection"));
   genToken_                      = consumes<std::vector<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("genParticleCollection"));
   lheEventsToken_                = consumes<LHEEventProduct>(iConfig.getParameter<edm::InputTag>("lheEvents"));
   ebRechitToken_                 = consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("ebRechitCollection"));
   eeRechitToken_                 = consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("eeRechitCollection"));
   patElectronToken_              = consumes<std::vector<pat::Electron> >(iConfig.getParameter<edm::InputTag>("patElectronCollection"));
   patPhotonToken_                = consumes<std::vector<pat::Photon> >(iConfig.getParameter<edm::InputTag>("patPhotonCollection")); 
   
   debug_                         = iConfig.getParameter<bool>("debug"); 
   isMC_                          = iConfig.getParameter<bool>("isMC"); 
   isParticleGun_                 = iConfig.getParameter<bool>("isParticleGun"); 
   motherPdgId_                   = iConfig.getParameter<std::vector<int> >( "motherPdgId" );
   motherStatus_                  = iConfig.getParameter<std::vector<int> >( "motherStatus" );
   pdgIdBeforeFSR_                = iConfig.getParameter<std::vector<int> >( "pdgIdBeforeFSR" );
   pdgId_                         = iConfig.getParameter<std::vector<int> >( "pdgId" );
   status_                        = iConfig.getParameter<std::vector<int> >( "status" );
   ignoreTauDecays_               = iConfig.getParameter<bool>( "ignoreTauDecays" );
   dRMin_                         = iConfig.getParameter<double>( "dRMin" );
   doLHEMatching_                 = iConfig.getParameter<bool>( "doLHEMatching" );
   dRMaxLHEMatch_                 = iConfig.getParameter<double>( "dRMaxLHEMatch" );
   dEFracMaxLHEMatch_             = iConfig.getParameter<double>( "dEFracMaxLHEMatch" );
    
   savePhotons_                   = iConfig.getParameter<bool>("savePhotons"); 
   saveElectrons_                 = iConfig.getParameter<bool>("saveElectrons");

   egmCutBasedElectronIDVeto_     = iConfig.getParameter<std::string>("egmCutBasedElectronIDVeto"); 
   egmCutBasedElectronIDloose_    = iConfig.getParameter<std::string>("egmCutBasedElectronIDloose");  
   egmCutBasedElectronIDmedium_   = iConfig.getParameter<std::string>("egmCutBasedElectronIDmedium"); 
   egmCutBasedElectronIDtight_    = iConfig.getParameter<std::string>("egmCutBasedElectronIDtight");   
   egmMVAElectronIDloose_         = iConfig.getParameter<std::string>("egmMVAElectronIDloose");  
   egmMVAElectronIDmedium_        = iConfig.getParameter<std::string>("egmMVAElectronIDmedium"); 
   egmMVAElectronIDtight_         = iConfig.getParameter<std::string>("egmMVAElectronIDtight");   
   egmMVAElectronIDlooseNoIso_    = iConfig.getParameter<std::string>("egmMVAElectronIDlooseNoIso");  
   egmMVAElectronIDmediumNoIso_   = iConfig.getParameter<std::string>("egmMVAElectronIDmediumNoIso"); 
   egmMVAElectronIDtightNoIso_    = iConfig.getParameter<std::string>("egmMVAElectronIDtightNoIso");
   heepElectronID_                = iConfig.getParameter<std::string>("heepElectronID");   
   egmCutBasedPhotonIDloose_      = iConfig.getParameter<std::string>("egmCutBasedPhotonIDloose");  
   egmCutBasedPhotonIDmedium_     = iConfig.getParameter<std::string>("egmCutBasedPhotonIDmedium"); 
   egmCutBasedPhotonIDtight_      = iConfig.getParameter<std::string>("egmCutBasedPhotonIDtight"); 
   egmMVAPhotonIDmedium_          = iConfig.getParameter<std::string>("egmMVAPhotonIDmedium"); 
   egmMVAPhotonIDtight_           = iConfig.getParameter<std::string>("egmMVAPhotonIDtight");
  
   //output file, historgrams and trees
   tree = iFile->make<TTree>("caloTree","caloTree"); 
   setTree(tree);
}

ElePhoDumper::~ElePhoDumper()
{
        // do anything here that needs to be done at desctruction time
        // (e.g. close files, deallocate resources etc.)
}


//
// member functions
//

// ------------ method called to for each event  ------------
void ElePhoDumper::analyze(const edm::Event& ev, const edm::EventSetup& iSetup)
{
   //MC-only info and collections
   truePU=-1.;
   obsPU=-1.;
   edm::Handle<std::vector<PileupSummaryInfo> > PupInfo;
   edm::Handle<std::vector<reco::GenParticle> > genParticles; 

   if(isMC_){    
      ev.getByToken(pileupSummaryToken_, PupInfo);
      if(PupInfo.isValid()) 
      {
         for(auto &pu : *PupInfo){
           if(pu.getBunchCrossing() == 0 ){
              truePU = pu.getTrueNumInteractions();
              obsPU = pu.getPU_NumInteractions();
              break;
           } 
         } 
      }else{
         std::cerr << "Analyze --> PupInfo not found" << std::endl;
      }

      ev.getByToken(genToken_,genParticles);
      if(!genParticles.isValid()) {
         std::cerr << "Analyze --> genParticles not found" << std::endl; 
         return;
      }
    
   }

   //Other collections
   edm::Handle<double> rhos;
   ev.getByToken(rhoToken_,rhos);
   if (!rhos.isValid()) {
       std::cerr << "Analyze --> rhos not found" << std::endl; 
       return;
   }

   edm::Handle<reco::VertexCollection> vertices;
   ev.getByToken(vtxToken_,vertices);
   if (!vertices.isValid()) {
       std::cerr << "Analyze --> vertices not found" << std::endl; 
       return;
   }
   
   edm::Handle<EcalRecHitCollection> recHitsEB;
   ev.getByToken(ebRechitToken_, recHitsEB);
   if (!recHitsEB.isValid()) {
       std::cerr << "Analyze --> recHitsEB not found" << std::endl; 
       return;
   }

   edm::Handle<EcalRecHitCollection> recHitsEE;
   ev.getByToken(eeRechitToken_, recHitsEE);
   if (!recHitsEE.isValid()) {
       std::cerr << "Analyze --> recHitsEE not found" << std::endl; 
       return;
   } 
    
   edm::Handle<std::vector<pat::Photon> > patPhoton;
   ev.getByToken(patPhotonToken_,patPhoton);
   if(!patPhoton.isValid()) {
      std::cerr << "Analyze --> patPhotons not found" << std::endl; 
      return;
   }
   
   edm::Handle<std::vector<pat::Electron> > patElectron;
   ev.getByToken(patElectronToken_,patElectron);
   if (!patElectron.isValid()) {
       std::cerr << "Analyze --> patElectrons not found" << std::endl; 
       return;
   }
   
   edm::Handle<LHEEventProduct> lheEvents;
   if(doLHEMatching_){
      ev.getByToken(lheEventsToken_, lheEvents);
      if(!lheEvents.isValid()) {
         std::cerr << "Analyze --> LHEEventProduct not found" << std::endl; 
         return;
      }
   }
   
   lheParticles.clear();
   lheParticlesID.clear();
   lheParticlesStatus.clear();
   if(doLHEMatching_){
      const auto& hepeup = lheEvents->hepeup();
      const auto& particles = hepeup.PUP;
      if(debug_) std::cout << "LHE-particles: " << particles.size() << std::endl;
         for(size_t i = 0; i < particles.size(); ++i){
             const auto& p = particles[i];
             lvec.SetPxPyPzE(p[0],p[1],p[2],p[3]);
             if(debug_){
                std::cout<< "Particle " << i
                         << ": pdgId = " << hepeup.IDUP[i]
        	         << ", status = " << hepeup.ISTUP[i]
        	         << ", px = " << p[0]
        	         << ", py = " << p[1]
        	         << ", pz = " << p[2]
        	         << ", eta = " << lvec.Eta()
        	         << ", phi = " << lvec.Phi()
        	         << ", pt = " << lvec.Pt()
        	         << ", E = " << p[3]
        	         << ", m = " << p[4] << std::endl;
             }            
             lheParticles.push_back(lvec);
             lheParticlesID.push_back(hepeup.IDUP[i]);
             lheParticlesStatus.push_back(hepeup.ISTUP[i]);
         }    	 
   }
   
   edm::ESHandle<CaloTopology> caloTopology;
   iSetup.get<CaloTopologyRecord>().get(caloTopology);
   const CaloTopology* topology = caloTopology.product();
   
   runId = ev.id().run();
   lumiId = ev.luminosityBlock();
   eventId = ev.id().event();
   nVtx = vertices->size();
   rho = *(rhos.product());
   
   vtx0 = &(*(vertices.product())).at(0);  

   if(isMC_) genParticle_size = (*(genParticles.product())).size();
   patElectron_size = (*(patElectron.product())).size();
   patPhoton_size = (*(patPhoton.product())).size();
   if(saveElectrons_) patPhoton_size = patElectron_size;
   setVectors(genParticle_size, patElectron_size,patPhoton_size);
   
   genIndices.clear();
   fsrEnergyMap.clear();
   if(isMC_){ 
      int iGen=-1;
      for(const auto& iGenPart : *(genParticles.product()))
      {
          iGen++; 
          genParticle_pdgId[iGen] = iGenPart.pdgId(); 
          genParticle_status[iGen] = iGenPart.status(); 
          genParticle_energy[iGen] = iGenPart.energy(); 
          genParticle_pt[iGen] = iGenPart.pt(); 
          genParticle_eta[iGen] = iGenPart.eta(); 
          genParticle_phi[iGen] = iGenPart.phi();   
          genParticle_mass[iGen] = iGenPart.mass();  
          genParticle_charge[iGen] = iGenPart.charge(); 
          genParticle_isLHEMatched[iGen] = false;  
          if(isParticleGun_){  
             genIndices.push_back(iGen);
             continue;
          } 
          
          const reco::GenParticle* motherPart = nullptr;
          if(saveElectrons_)
          { 
             if(abs(iGenPart.pdgId())!=11 || iGenPart.status()!=1) continue;
             motherPart = getFirstMother( &iGenPart, 11, -999 ); 
             if(motherPart==nullptr){ 
                //if(debug_) std::cout << "WARNING: no good mother!" << std::endl; 
                continue;
             }
             if(abs(motherPart->pdgId())!=11) continue;
             if(debug_ && motherPart!=nullptr) std::cout << "good mother found! - pdgId = " << motherPart->pdgId() << " - status = " << motherPart->status() << std::endl;
          } 
          if(savePhotons_)
          { 
             if(abs(iGenPart.pdgId())!=22 || iGenPart.status()!=1) continue;
             motherPart = nullptr;
            
             if(iGenPart.numberOfMothers()==0) motherPart = &iGenPart;
             else{
                motherPart = getFirstMother( &iGenPart, 11, -999 ); 
                if(abs(motherPart->pdgId())!=11) motherPart = iGenPart.motherRef().get();
             }   
             if(motherPart==nullptr){ 
                if(debug_) std::cout << "WARNING: no good mother!" << std::endl; 
                continue;
             }
             if(abs(motherPart->pdgId())>99) continue;
             if(abs(motherPart->pdgId())==23) motherPart = &iGenPart;
          }
         
          if(doLHEMatching_)
          {
             int lheIndex = -1;
             double dR_tmp=dRMaxLHEMatch_; 
             for(unsigned int j=0; j<lheParticles.size(); j++)
             {
        	 if(abs(lheParticlesID.at(j))!=11) continue;
        	 if(lheParticlesStatus.at(j)!=1) continue;
                 double dR = deltaR(motherPart->eta(),motherPart->phi(),lheParticles.at(j).Eta(),lheParticles.at(j).Phi()); 
                 if(dR<dR_tmp){
                    dR_tmp = dR;
                    lheIndex = j; 
                 }
             }
             if(lheIndex<0) continue;
             double dE_fraction = abs(lheParticles.at(lheIndex).Energy()-iGenPart.energy())/lheParticles.at(lheIndex).Energy(); 
             if(dE_fraction>dEFracMaxLHEMatch_) continue;
             genParticle_isLHEMatched[iGen] = true;  
             if(debug_){ 
                std::cout << "Particle before FSR: pdgId = " << motherPart->pdgId() << " - status = " <<  motherPart->status() << " -  eta = " << motherPart->eta() << " - phi = " << motherPart->phi() << " - pt = " << motherPart->pt() << " -  energy = " << motherPart->energy() << std::endl; 
                std::cout << " --> Matched lhe-particle: dR = " << dR_tmp << " - lheIndex = " << lheIndex << std::endl;  
             } 
          }      
          genIndices.push_back(iGen); 
      }
      
      if(genIndices.size()==0){ 
         std::cout << "WARINING: genIndices-size = " << genIndices.size() << " --> Skipping event! " << std::endl;
         return;
      }   
      
      //Remove duplicates
      std::sort( genIndices.begin(), genIndices.end() );
      genIndices.erase( unique( genIndices.begin(), genIndices.end() ), genIndices.end() );
      
      if(debug_){
         if(!isParticleGun_){ 
            std::cout << "GenParticles trees... " << std::endl;
            for(unsigned int i=0; i<genIndices.size(); i++)
            {
                const reco::GenParticle* motherPart = dynamic_cast<const reco::GenParticle*>((&genParticles->at(genIndices.at(i)))->motherRef().get()); 
                printAllDaughters( motherPart, 0 );      
            }   
         }
      }
    
      //check FSR
      if(saveElectrons_){
         for(const auto& iGenPart : *(genParticles.product()))
         {     
             if(abs(iGenPart.pdgId())!=22 || iGenPart.status()!=1) continue;
         
             const reco::GenParticle* mother = getFirstMother( &iGenPart, 11, -999 ); 
             //if(abs(mother->pdgId())==11) std::cout << "mother: pdgId = " << mother->pdgId() << " - status = " << mother->status() << std::endl; 
             if(abs(mother->pdgId())!=11 || mother->status()==1) continue;
           
             if(mother==nullptr){ 
                if(debug_) std::cout << "WARNING: no good particle before FSR!" << std::endl; 
                continue;
             }
             
             int genIndex = getGenIndex(genParticles,mother); 
             if(debug_) std::cout << "mother: pdgId = " << mother->pdgId() << " -  status = " << mother->status() << " - index = " << genIndex << " - fsrEnrgy = " << iGenPart.energy() << std::endl; 
             if(fsrEnergyMap.find(genIndex) == fsrEnergyMap.end()) fsrEnergyMap[genIndex] = iGenPart.energy();
             else fsrEnergyMap[genIndex] += iGenPart.energy();
         }
         if(debug_){
            for(auto &myPair : fsrEnergyMap){
                std::cout << "FsrMap: " << myPair.first << " - " << myPair.second << std::endl;
            }
         }
      }
   }
           
   //save electron info
   if(saveElectrons_)
   { 
      int iEle=0;
      matchedPhotons.clear();
      for(const auto& iElectron : *(patElectron.product())){ 
 
          reco::SuperClusterRef scRef = iElectron.superCluster();
          double R  = TMath::Sqrt(scRef->x()*scRef->x() + scRef->y()*scRef->y() +scRef->z()*scRef->z());
          double Rt = TMath::Sqrt(scRef->x()*scRef->x() + scRef->y()*scRef->y());
          bool scIsEB = (scRef->seed()->seed().subdetId()==EcalBarrel);
          
          double swissCross = -999.;
          double e3x3 = -999.; 
          double full5x5_e3x3 = -999.; 
          double full5x5_e5x5 = -999.; 
          double full5x5_eMax = -999.; 
          double full5x5_e2nd  = -999.;
          double full5x5_eTop  = -999.;
          double full5x5_eBottom  = -999.;
          double full5x5_eLeft  = -999.;
          double full5x5_eRight  = -999.;
          double full5x5_e2x5Max  = -999.;
          double full5x5_e2x5Left  = -999.;
          double full5x5_e2x5Right  = -999.;
          double full5x5_e2x5Top  = -999.;
          double full5x5_e2x5Bottom  = -999.;
          double full5x5_sigmaIEtaIEta = -999.;
          double full5x5_sigmaIEtaIPhi = -999.;
          double full5x5_sigmaIPhiIPhi = -999.;
         
          reco::GsfElectron::ShowerShape eleSS = iElectron.showerShape();
          reco::GsfElectron::ShowerShape full5x5_eleSS = iElectron.full5x5_showerShape();
          const std::vector<std::pair<DetId,float> > &hits= iElectron.superCluster()->hitsAndFractions();
          if(scIsEB)
          {
             std::pair<DetId, float> id = EcalClusterTools::getMaximum(hits,&(*(recHitsEB.product())));   
             swissCross = EcalTools::swissCross(id.first,*(recHitsEB.product()),0.);     
             e3x3 = EcalClusterTools::e3x3( *(scRef->seed()), &(*(recHitsEB.product())), topology);
             full5x5_e3x3 = noZS::EcalClusterTools::e3x3( *(scRef->seed()), &(*(recHitsEB.product())), topology);
             full5x5_e5x5 = noZS::EcalClusterTools::e5x5( *(scRef->seed()), &(*(recHitsEB.product())), topology);
             full5x5_eMax = noZS::EcalClusterTools::eMax( *(scRef->seed()), &(*(recHitsEB.product())));
             full5x5_e2nd = noZS::EcalClusterTools::e2nd( *(scRef->seed()), &(*(recHitsEB.product())));
             full5x5_eTop = noZS::EcalClusterTools::eTop( *(scRef->seed()), &(*(recHitsEB.product())), topology);
             full5x5_eBottom = noZS::EcalClusterTools::eBottom( *(scRef->seed()), &(*(recHitsEB.product())), topology);
             full5x5_eLeft = noZS::EcalClusterTools::eLeft( *(scRef->seed()), &(*(recHitsEB.product())), topology);
             full5x5_eRight = noZS::EcalClusterTools::eRight( *(scRef->seed()), &(*(recHitsEB.product())), topology);
             full5x5_e2x5Max = noZS::EcalClusterTools::e2x5Max( *(scRef->seed()), &(*(recHitsEB.product())), topology);
             full5x5_e2x5Left = noZS::EcalClusterTools::e2x5Left( *(scRef->seed()), &(*(recHitsEB.product())), topology);
             full5x5_e2x5Right = noZS::EcalClusterTools::e2x5Right( *(scRef->seed()), &(*(recHitsEB.product())), topology);
             full5x5_e2x5Top = noZS::EcalClusterTools::e2x5Top( *(scRef->seed()), &(*(recHitsEB.product())), topology);
             full5x5_e2x5Bottom = noZS::EcalClusterTools::e2x5Bottom( *(scRef->seed()), &(*(recHitsEB.product())), topology);
             
             const auto localCovs =  noZS::EcalClusterTools::localCovariances(*(scRef->seed()), &(*(recHitsEB.product())), topology);
             full5x5_sigmaIEtaIEta =std::sqrt(localCovs[0]);
             full5x5_sigmaIEtaIPhi = std::numeric_limits<float>::max();
             full5x5_sigmaIPhiIPhi = std::numeric_limits<float>::max();
             if(!edm::isNotFinite(localCovs[2])) full5x5_sigmaIPhiIPhi = std::sqrt(localCovs[2]) ;
             const bool applySPPBug = false;
             const float seeBySpp = applySPPBug ? full5x5_sigmaIEtaIEta*std::numeric_limits<float>::max() : full5x5_sigmaIEtaIEta*full5x5_sigmaIPhiIPhi;  
             if( seeBySpp > 0 ) { 
                 full5x5_sigmaIEtaIPhi = localCovs[1] / seeBySpp;
             } else if ( localCovs[1] > 0 ) {
                 full5x5_sigmaIEtaIPhi = 1.f;
             } else {
                 full5x5_sigmaIEtaIPhi = -1.f;
             }
          }else{
             std::pair<DetId, float> id = EcalClusterTools::getMaximum(hits,&(*(recHitsEE.product())));   
             swissCross = EcalTools::swissCross(id.first,*(recHitsEE.product()),0.);
             e3x3 = EcalClusterTools::e3x3( *(scRef->seed()), &(*(recHitsEE.product())), topology);
             full5x5_e3x3 = noZS::EcalClusterTools::e3x3( *(scRef->seed()), &(*(recHitsEE.product())), topology);
             full5x5_e5x5 = noZS::EcalClusterTools::e5x5( *(scRef->seed()), &(*(recHitsEE.product())), topology);
             full5x5_eMax = noZS::EcalClusterTools::eMax( *(scRef->seed()), &(*(recHitsEE.product())));
             full5x5_e2nd = noZS::EcalClusterTools::e2nd( *(scRef->seed()), &(*(recHitsEE.product())));
             full5x5_eTop = noZS::EcalClusterTools::eTop( *(scRef->seed()), &(*(recHitsEE.product())), topology);
             full5x5_eBottom = noZS::EcalClusterTools::eBottom( *(scRef->seed()), &(*(recHitsEE.product())), topology);
             full5x5_eLeft = noZS::EcalClusterTools::eLeft( *(scRef->seed()), &(*(recHitsEE.product())), topology);
             full5x5_eRight = noZS::EcalClusterTools::eRight( *(scRef->seed()), &(*(recHitsEE.product())), topology);
             full5x5_e2x5Max = noZS::EcalClusterTools::e2x5Max( *(scRef->seed()), &(*(recHitsEE.product())), topology);
             full5x5_e2x5Left = noZS::EcalClusterTools::e2x5Left( *(scRef->seed()), &(*(recHitsEE.product())), topology);
             full5x5_e2x5Right = noZS::EcalClusterTools::e2x5Right( *(scRef->seed()), &(*(recHitsEE.product())), topology);
             full5x5_e2x5Top = noZS::EcalClusterTools::e2x5Top( *(scRef->seed()), &(*(recHitsEE.product())), topology);
             full5x5_e2x5Bottom = noZS::EcalClusterTools::e2x5Bottom( *(scRef->seed()), &(*(recHitsEE.product())), topology);
             
             const auto localCovs =  noZS::EcalClusterTools::localCovariances(*(scRef->seed()), &(*(recHitsEE.product())), topology);
             full5x5_sigmaIEtaIEta =std::sqrt(localCovs[0]);
             full5x5_sigmaIEtaIPhi = std::numeric_limits<float>::max();
             full5x5_sigmaIPhiIPhi = std::numeric_limits<float>::max();
             if(!edm::isNotFinite(localCovs[2])) full5x5_sigmaIPhiIPhi = std::sqrt(localCovs[2]) ;
             const bool applySPPBug = false;
             const float seeBySpp = applySPPBug ? full5x5_sigmaIEtaIEta*std::numeric_limits<float>::max() : full5x5_sigmaIEtaIEta*full5x5_sigmaIPhiIPhi;  
             if( seeBySpp > 0 ) { 
                 full5x5_sigmaIEtaIPhi = localCovs[1] / seeBySpp;
             } else if ( localCovs[1] > 0 ) {
                 full5x5_sigmaIEtaIPhi = 1.f;
             } else {
                 full5x5_sigmaIEtaIPhi = -1.f;
             }
          }
          
          p4_vtx0 = setP4(scRef, *vtx0, iElectron.correctedEcalEnergy());
          
          double dR_tmp = dRMin_;
          double dR = -999.;
          int genIndex = -1;
          if(isMC_){
             for(unsigned int i=0; i<genIndices.size(); i++){
                 int genIndex_tmp = genIndices.at(i);
                 dR = deltaR((&genParticles->at(genIndex_tmp))->eta(),(&genParticles->at(genIndex_tmp))->phi(),iElectron.eta(),iElectron.phi());
                 if(dR<=dR_tmp){
                    dR_tmp=dR;
                    genIndex = genIndex_tmp;
                 }
             }
             if(genIndex>=0) dR = deltaR((&genParticles->at(genIndex))->eta(),(&genParticles->at(genIndex))->phi(),iElectron.eta(),iElectron.phi());
             else dR = 999.;
          }
          
          if(genIndex>=0 && !isParticleGun_ && debug_) std::cout << "genIndex = " << genIndex << " - dR = " << dR << " - pdgId = " << (&genParticles->at(genIndex))->pdgId() << " - status = " << (&genParticles->at(genIndex))->status() << " - eta = " << (&genParticles->at(genIndex))->eta() << " - phi = " << (&genParticles->at(genIndex))->phi() << " - pt = " << (&genParticles->at(genIndex))->pt() << " - genEnergy = " << (&genParticles->at(genIndex))->energy() << " - electronEnergy = " << iElectron.energy() << std::endl; 
          
          const reco::GenParticle* mother = getFirstMother( &genParticles->at(genIndex), 11, -999 ); 
          if(mother==nullptr){ 
              if(debug_) std::cout << "WARNING: no good mother!" << std::endl; 
              continue;
          }
          int genMotherIndex = getGenIndex(genParticles,mother); 
          
          //compute photon regression
          if(patPhoton.isValid()) {
             int iPho = -1;
             auto eleSeedId = iElectron.superCluster()->seed()->seed().rawId();
             for(const auto& iPhoton : *(patPhoton.product()))
             {
                 iPho++;
                 auto phoSeedId = iPhoton.superCluster()->seed()->seed().rawId();
                 if(phoSeedId==eleSeedId){
                    patElectron_ecalEnergyPhoReg[iEle] = iPhoton.energy();
                    patElectron_ecalEnergyPhoRegErr[iEle] = iPhoton.getCorrectedEnergyError(iPhoton.getCandidateP4type());
                    matchedPhotons[iPho] = iEle;
                    break;
                 }
	     }
	  }
         
          const int iEtaCorr = scIsEB ? (EBDetId(scRef->seed()->seed().rawId()).ieta() - (EBDetId(scRef->seed()->seed().rawId()).ieta() > 0 ? +1 : -1)) : -999;
          const int iEtaCorr26 = scIsEB ? (EBDetId(scRef->seed()->seed().rawId()).ieta() - (EBDetId(scRef->seed()->seed().rawId()).ieta() > 0 ? +26 : -26)) : -999;
          
          patElectron_index[iEle] = iEle;
          patElectron_genIndex[iEle] = genIndex;
          patElectron_genDR[iEle] = dR;
          patElectron_genFSREnergy[iEle] = fsrEnergyMap[genMotherIndex];
          patElectron_classification[iEle] = iElectron.classification();
          patElectron_nrSatCrys[iEle] = iElectron.nSaturatedXtals();
          patElectron_scNPFClusters[iEle] = scRef->clusters().size();
          patElectron_charge[iEle] = iElectron.charge(); 
          patElectron_isEB[iEle] = iElectron.isEB(); 
          patElectron_isEE[iEle] = iElectron.isEE();  
          patElectron_isEBEEGap[iEle] = iElectron.isEBEEGap();  
          patElectron_isEBEtaGap[iEle] = iElectron.isEBEtaGap();  
          patElectron_isEBPhiGap[iEle] = iElectron.isEBPhiGap();  
          patElectron_isEEDeeGap[iEle] = iElectron.isEEDeeGap();  
          patElectron_isEERingGap[iEle] = iElectron.isEERingGap();   
          patElectron_isEcalDriven[iEle] = iElectron.ecalDrivenSeed(); 
          patElectron_isTrackerDriven[iEle] = iElectron.trackerDrivenSeed(); 
          patElectron_passConversionVeto[iEle] = iElectron.passConversionVeto(); 
          patElectron_eta[iEle] = iElectron.p4().eta();
          patElectron_phi[iEle] = iElectron.p4().phi();
          patElectron_p[iEle] = iElectron.trackMomentumAtVtx().R();
          patElectron_pt[iEle] = TMath::Sqrt(iElectron.trackMomentumAtVtx().Perp2());
          patElectron_pIn[iEle] = iElectron.trackMomentumAtVtx().R();
          patElectron_pOut[iEle] = iElectron.trackMomentumOut().R();
          patElectron_trackFbrem[iEle] = iElectron.trackFbrem();
          patElectron_superClusterFbrem[iEle] = iElectron.superClusterFbrem(); 
          patElectron_energy[iEle] = iElectron.energy();
          patElectron_energyErr[iEle] = iElectron.p4Error(reco::GsfElectron::P4_COMBINATION);
          patElectron_ecalEnergy[iEle] = iElectron.ecalEnergy();
          patElectron_ecalEnergyErr[iEle] = iElectron.ecalEnergyError();
          patElectron_et[iEle] = iElectron.energy()/TMath::CosH(iElectron.eta());
          patElectron_HoE[iEle] = iElectron.hcalOverEcalBc();
          patElectron_scIsEB[iEle] = scIsEB;
          patElectron_scEta[iEle] = scRef->eta();
          patElectron_scPhi[iEle] = scRef->phi();
          patElectron_scEta_vtx0[iEle] = p4_vtx0.Eta();
          patElectron_scPhi_vtx0[iEle] = p4_vtx0.Phi();
          patElectron_scEnergy[iEle] = scRef->energy(); 
          patElectron_scRawEnergy[iEle] = scRef->rawEnergy(); 
          patElectron_scRawESEnergy[iEle] = scRef->preshowerEnergy(); 
          patElectron_scSeedEnergy[iEle] = scRef->seed()->energy(); 
          patElectron_scDPhiSeed[iEle] = reco::deltaPhi(scRef->seed()->phi(),scRef->position().Phi());
          patElectron_scDEtaSeed[iEle] = scRef->seed()->eta() - scRef->position().Eta();
          patElectron_scIEtaOrX[iEle] = scIsEB ? EBDetId(scRef->seed()->seed().rawId()).ieta() : EEDetId(scRef->seed()->seed().rawId()).ix(); 
          patElectron_scIPhiOrY[iEle] = scIsEB ? EBDetId(scRef->seed()->seed().rawId()).iphi() : EEDetId(scRef->seed()->seed().rawId()).iy(); 
          patElectron_scIEtaMod5[iEle] = scIsEB ? iEtaCorr%5 : -999;
          patElectron_scIPhiMod2[iEle] = scIsEB ? (EBDetId(scRef->seed()->seed().rawId()).iphi()-1)%2 : -999;
          patElectron_scIEtaMod20[iEle] = scIsEB ? (std::abs(EBDetId(scRef->seed()->seed().rawId()).ieta())<=25 ? iEtaCorr%20 : iEtaCorr26%20) : -999;
          patElectron_scIPhiMod20[iEle] = scIsEB ? (EBDetId(scRef->seed()->seed().rawId()).iphi()-1)%20 : -999;
          patElectron_scEt[iEle] = scRef->energy()*(Rt/R);
          patElectron_scPhiWidth[iEle] = scRef->phiWidth();
          patElectron_scEtaWidth[iEle] = scRef->etaWidth(); 
          patElectron_scEoP[iEle] = scRef->energy()/iElectron.trackMomentumAtVtx().R(); 
          patElectron_SwissCross[iEle] = swissCross; 
          patElectron_r9[iEle] = e3x3/scRef->rawEnergy(); 
          patElectron_full5x5_r9[iEle] = full5x5_e3x3/scRef->rawEnergy();
          patElectron_full5x5_e3x3[iEle] = full5x5_e3x3;
          patElectron_full5x5_e5x5[iEle] = full5x5_e5x5;
          patElectron_full5x5_eMax[iEle] = full5x5_eMax;
          patElectron_full5x5_e2nd[iEle] = full5x5_e2nd;
          patElectron_full5x5_eTop[iEle] = full5x5_eTop;
          patElectron_full5x5_eBottom[iEle] = full5x5_eBottom;
          patElectron_full5x5_eLeft[iEle] = full5x5_eLeft;
          patElectron_full5x5_eRight[iEle] = full5x5_eRight;
          patElectron_full5x5_e2x5Max[iEle] = full5x5_e2x5Max;
          patElectron_full5x5_e2x5Left[iEle] = full5x5_e2x5Left;
          patElectron_full5x5_e2x5Right[iEle] = full5x5_e2x5Right;
          patElectron_full5x5_e2x5Top[iEle] = full5x5_e2x5Top;
          patElectron_full5x5_e2x5Bottom[iEle] = full5x5_e2x5Bottom;
          patElectron_full5x5_sigmaIEtaIEta[iEle] = full5x5_sigmaIEtaIEta;
          patElectron_full5x5_sigmaIEtaIPhi[iEle] = full5x5_sigmaIEtaIPhi;
          patElectron_full5x5_sigmaIPhiIPhi[iEle] = full5x5_sigmaIPhiIPhi; 
          patElectron_egmMVAElectronIDScore[iEle]  = iElectron.userFloat("ElectronMVAEstimatorRun2Fall17IsoV2Values"); 
          patElectron_egmMVAElectronIDNoIsoScore[iEle]  = iElectron.userFloat("ElectronMVAEstimatorRun2Fall17NoIsoV2Values"); 
          patElectron_egmCutBasedElectronIDVeto[iEle] = iElectron.electronID(egmCutBasedElectronIDVeto_.c_str());
          patElectron_egmCutBasedElectronIDloose[iEle] = iElectron.electronID(egmCutBasedElectronIDloose_.c_str());
          patElectron_egmCutBasedElectronIDmedium[iEle] = iElectron.electronID(egmCutBasedElectronIDmedium_.c_str());
          patElectron_egmCutBasedElectronIDtight[iEle] = iElectron.electronID(egmCutBasedElectronIDtight_.c_str());
          patElectron_egmMVAElectronIDloose[iEle] = iElectron.electronID(egmMVAElectronIDloose_.c_str());
          patElectron_egmMVAElectronIDmedium[iEle] = iElectron.electronID(egmMVAElectronIDmedium_.c_str());
          patElectron_egmMVAElectronIDtight[iEle] = iElectron.electronID(egmMVAElectronIDtight_.c_str());
          patElectron_egmMVAElectronIDlooseNoIso[iEle] = iElectron.electronID(egmMVAElectronIDlooseNoIso_.c_str());
          patElectron_egmMVAElectronIDmediumNoIso[iEle] = iElectron.electronID(egmMVAElectronIDmediumNoIso_.c_str());
          patElectron_egmMVAElectronIDtightNoIso[iEle] = iElectron.electronID(egmMVAElectronIDtightNoIso_.c_str());
          patElectron_heepElectronID[iEle] = iElectron.electronID(heepElectronID_.c_str());
          
          iEle++; 
      }
   }    
   
   //save photon info
   if(savePhotons_ || saveElectrons_)
   {    
      int iPho=0;
      int iPho_tmp=-1;
      for(const auto& iPhoton : *(patPhoton.product())){ 
      
          if(saveElectrons_){
             iPho_tmp++;
             if(matchedPhotons.find(iPho_tmp)==matchedPhotons.end()) continue; 
             int iEle = matchedPhotons[iPho_tmp];
             iPho = iEle;
          } 
         
          reco::SuperClusterRef scRef = iPhoton.superCluster();
          reco::PhotonCoreRef phoCoreRef = iPhoton.photonCore();
          double R  = TMath::Sqrt(scRef->x()*scRef->x() + scRef->y()*scRef->y() +scRef->z()*scRef->z());
          double Rt = TMath::Sqrt(scRef->x()*scRef->x() + scRef->y()*scRef->y());
          bool scIsEB = (scRef->seed()->seed().subdetId()==EcalBarrel);

          double swissCross = -999.;
          double e3x3 = -999.; 
          double full5x5_e3x3 = -999.; 
          double full5x5_e5x5 = -999.; 
          double full5x5_eMax = -999.; 
          double full5x5_e2nd  = -999.;
          double full5x5_eTop  = -999.;
          double full5x5_eBottom  = -999.;
          double full5x5_eLeft  = -999.;
          double full5x5_eRight  = -999.;
          double full5x5_e2x5Max  = -999.;
          double full5x5_e2x5Left  = -999.;
          double full5x5_e2x5Right  = -999.;
          double full5x5_e2x5Top  = -999.;
          double full5x5_e2x5Bottom  = -999.;
          double full5x5_sigmaIEtaIEta = -999.;
          double full5x5_sigmaIEtaIPhi = -999.;
          double full5x5_sigmaIPhiIPhi = -999.;
          
          reco::Photon::ShowerShape phoSS = iPhoton.showerShapeVariables(); 
          reco::Photon::ShowerShape full5x5_phoSS = iPhoton.full5x5_showerShapeVariables();    
          const std::vector<std::pair<DetId,float> > &hits= iPhoton.superCluster()->hitsAndFractions();
          if(scIsEB)
          {
             std::pair<DetId, float> id = EcalClusterTools::getMaximum(hits,&(*(recHitsEB.product())));   
             swissCross = EcalTools::swissCross(id.first,*(recHitsEB.product()),0.);     
             e3x3 = EcalClusterTools::e3x3( *(scRef->seed()), &(*(recHitsEB.product())), topology);
             full5x5_e3x3 = noZS::EcalClusterTools::e3x3( *(scRef->seed()), &(*(recHitsEB.product())), topology);
             full5x5_e5x5 = noZS::EcalClusterTools::e5x5( *(scRef->seed()), &(*(recHitsEB.product())), topology);
             full5x5_eMax = noZS::EcalClusterTools::eMax( *(scRef->seed()), &(*(recHitsEB.product())));
             full5x5_e2nd = noZS::EcalClusterTools::e2nd( *(scRef->seed()), &(*(recHitsEB.product())));
             full5x5_eTop = noZS::EcalClusterTools::eTop( *(scRef->seed()), &(*(recHitsEB.product())), topology);
             full5x5_eBottom = noZS::EcalClusterTools::eBottom( *(scRef->seed()), &(*(recHitsEB.product())), topology);
             full5x5_eLeft = noZS::EcalClusterTools::eLeft( *(scRef->seed()), &(*(recHitsEB.product())), topology);
             full5x5_eRight = noZS::EcalClusterTools::eRight( *(scRef->seed()), &(*(recHitsEB.product())), topology);
             full5x5_e2x5Max = noZS::EcalClusterTools::e2x5Max( *(scRef->seed()), &(*(recHitsEB.product())), topology);
             full5x5_e2x5Left = noZS::EcalClusterTools::e2x5Left( *(scRef->seed()), &(*(recHitsEB.product())), topology);
             full5x5_e2x5Right = noZS::EcalClusterTools::e2x5Right( *(scRef->seed()), &(*(recHitsEB.product())), topology);
             full5x5_e2x5Top = noZS::EcalClusterTools::e2x5Top( *(scRef->seed()), &(*(recHitsEB.product())), topology);
             full5x5_e2x5Bottom = noZS::EcalClusterTools::e2x5Bottom( *(scRef->seed()), &(*(recHitsEB.product())), topology);
             
             const auto localCovs =  noZS::EcalClusterTools::localCovariances(*(scRef->seed()), &(*(recHitsEB.product())), topology);
             full5x5_sigmaIEtaIEta =std::sqrt(localCovs[0]);
             full5x5_sigmaIEtaIPhi = std::numeric_limits<float>::max();
             full5x5_sigmaIPhiIPhi = std::numeric_limits<float>::max();
             if(!edm::isNotFinite(localCovs[2])) full5x5_sigmaIPhiIPhi = std::sqrt(localCovs[2]) ;
             const bool applySPPBug = false;
             const float seeBySpp = applySPPBug ? full5x5_sigmaIEtaIEta*std::numeric_limits<float>::max() : full5x5_sigmaIEtaIEta*full5x5_sigmaIPhiIPhi;  
             if( seeBySpp > 0 ) { 
                 full5x5_sigmaIEtaIPhi = localCovs[1] / seeBySpp;
             } else if ( localCovs[1] > 0 ) {
                 full5x5_sigmaIEtaIPhi = 1.f;
             } else {
                 full5x5_sigmaIEtaIPhi = -1.f;
             }
          }else{
             std::pair<DetId, float> id = EcalClusterTools::getMaximum(hits,&(*(recHitsEE.product())));   
             swissCross = EcalTools::swissCross(id.first,*(recHitsEE.product()),0.);
             e3x3 = EcalClusterTools::e3x3( *(scRef->seed()), &(*(recHitsEE.product())), topology);
             full5x5_e3x3 = noZS::EcalClusterTools::e3x3( *(scRef->seed()), &(*(recHitsEE.product())), topology);
             full5x5_e5x5 = noZS::EcalClusterTools::e5x5( *(scRef->seed()), &(*(recHitsEE.product())), topology);
             full5x5_eMax = noZS::EcalClusterTools::eMax( *(scRef->seed()), &(*(recHitsEE.product())));
             full5x5_e2nd = noZS::EcalClusterTools::e2nd( *(scRef->seed()), &(*(recHitsEE.product())));
             full5x5_eTop = noZS::EcalClusterTools::eTop( *(scRef->seed()), &(*(recHitsEE.product())), topology);
             full5x5_eBottom = noZS::EcalClusterTools::eBottom( *(scRef->seed()), &(*(recHitsEE.product())), topology);
             full5x5_eLeft = noZS::EcalClusterTools::eLeft( *(scRef->seed()), &(*(recHitsEE.product())), topology);
             full5x5_eRight = noZS::EcalClusterTools::eRight( *(scRef->seed()), &(*(recHitsEE.product())), topology);
             full5x5_e2x5Max = noZS::EcalClusterTools::e2x5Max( *(scRef->seed()), &(*(recHitsEE.product())), topology);
             full5x5_e2x5Left = noZS::EcalClusterTools::e2x5Left( *(scRef->seed()), &(*(recHitsEE.product())), topology);
             full5x5_e2x5Right = noZS::EcalClusterTools::e2x5Right( *(scRef->seed()), &(*(recHitsEE.product())), topology);
             full5x5_e2x5Top = noZS::EcalClusterTools::e2x5Top( *(scRef->seed()), &(*(recHitsEE.product())), topology);
             full5x5_e2x5Bottom = noZS::EcalClusterTools::e2x5Bottom( *(scRef->seed()), &(*(recHitsEE.product())), topology);
             
             const auto localCovs =  noZS::EcalClusterTools::localCovariances(*(scRef->seed()), &(*(recHitsEE.product())), topology);
             full5x5_sigmaIEtaIEta =std::sqrt(localCovs[0]);
             full5x5_sigmaIEtaIPhi = std::numeric_limits<float>::max();
             full5x5_sigmaIPhiIPhi = std::numeric_limits<float>::max();
             if(!edm::isNotFinite(localCovs[2])) full5x5_sigmaIPhiIPhi = std::sqrt(localCovs[2]) ;
             const bool applySPPBug = false;
             const float seeBySpp = applySPPBug ? full5x5_sigmaIEtaIEta*std::numeric_limits<float>::max() : full5x5_sigmaIEtaIEta*full5x5_sigmaIPhiIPhi;  
             if( seeBySpp > 0 ) { 
                 full5x5_sigmaIEtaIPhi = localCovs[1] / seeBySpp;
             } else if ( localCovs[1] > 0 ) {
                 full5x5_sigmaIEtaIPhi = 1.f;
             } else {
                 full5x5_sigmaIEtaIPhi = -1.f;
             }
          }
          
          p4_vtx0 = setP4(scRef, *vtx0, iPhoton.energy());

          double dR_tmp = dRMin_;
          double dR = 999.;
          int genIndex = -1;
          if(isMC_){
             for(unsigned int i=0; i<genIndices.size(); i++){
                 int genIndex_tmp = genIndices.at(i);
                 dR = deltaR((&genParticles->at(genIndex_tmp))->eta(),(&genParticles->at(genIndex_tmp))->phi(),iPhoton.eta(),iPhoton.phi());
                 if(dR<=dR_tmp){
                    dR_tmp=dR;
                    genIndex = genIndex_tmp;
                 }
             }
             if(genIndex>=0) dR = deltaR((&genParticles->at(genIndex))->eta(),(&genParticles->at(genIndex))->phi(),iPhoton.eta(),iPhoton.phi());
             else dR = 999.;
          }
          
          if(genIndex>=0 && !isParticleGun_ && debug_) std::cout << "genIndex = " << genIndex << " - dR = " << dR << " - pdgId = " << (&genParticles->at(genIndex))->pdgId() << " - status = " << (&genParticles->at(genIndex))->status() << " - eta = " << (&genParticles->at(genIndex))->eta() << " - phi = " << (&genParticles->at(genIndex))->phi() << " - pt = " << (&genParticles->at(genIndex))->pt() << " - energy = " << (&genParticles->at(genIndex))->energy() << " - photonEnergy = " << iPhoton.energy() << std::endl; 
          
          const int iEtaCorr = scIsEB ? (EBDetId(scRef->seed()->seed().rawId()).ieta() - (EBDetId(scRef->seed()->seed().rawId()).ieta() > 0 ? +1 : -1)) : -999;
          const int iEtaCorr26 = scIsEB ? (EBDetId(scRef->seed()->seed().rawId()).ieta() - (EBDetId(scRef->seed()->seed().rawId()).ieta() > 0 ? +26 : -26)) : -999;
          
          patPhoton_index[iPho] = iPho;
          patPhoton_genIndex[iPho] = genIndex;
          patPhoton_genDR[iPho] = dR;
          patPhoton_nrSatCrys[iPho] = iPhoton.nSaturatedXtals();
          patPhoton_scNPFClusters[iPho] = scRef->clusters().size();
          patPhoton_isEB[iPho] = iPhoton.isEB();
          patPhoton_isEE[iPho] = iPhoton.isEE();  
          patPhoton_isEBEEGap[iPho] = iPhoton.isEBEEGap();  
          patPhoton_isEBEtaGap[iPho] = iPhoton.isEBEtaGap();  
          patPhoton_isEBPhiGap[iPho] = iPhoton.isEBPhiGap();  
          patPhoton_isEEDeeGap[iPho] = iPhoton.isEEDeeGap();  
          patPhoton_isEERingGap[iPho] = iPhoton.isEERingGap();   
          patPhoton_passElectronVeto[iPho] = iPhoton.passElectronVeto(); 
          patPhoton_hasPixelSeed[iPho] = iPhoton.hasPixelSeed();  
          patPhoton_hasConversionTracks[iPho] = iPhoton.hasConversionTracks();
          patPhoton_nConversions[iPho] = phoCoreRef->conversions().size();
          patPhoton_nConversionsOneLeg[iPho] = phoCoreRef->conversionsOneLeg().size();
          patPhoton_eta[iPho] = iPhoton.p4().eta();
          patPhoton_phi[iPho] = iPhoton.p4().phi();
          patPhoton_energy[iPho] = iPhoton.energy();
          patPhoton_energyErr[iPho] = iPhoton.getCorrectedEnergyError(reco::Photon::regression2);
          patPhoton_ecalEnergy[iPho] = iPhoton.energyCorrections().phoEcalEnergy;
          patPhoton_ecalEnergyErr[iPho] = iPhoton.energyCorrections().phoEcalEnergyError; 
          patPhoton_et[iPho] = iPhoton.p4().Et();
          patPhoton_pt[iPho] = iPhoton.pt();
          patPhoton_HoE[iPho] = iPhoton.hadronicOverEm(); 
          patPhoton_scIsEB[iPho] = (scRef->seed()->seed().subdetId()==EcalBarrel);
          patPhoton_scEta[iPho] = scRef->eta();
          patPhoton_scPhi[iPho] = scRef->phi();
          patPhoton_scEta_vtx0[iPho] = p4_vtx0.Eta();
          patPhoton_scPhi_vtx0[iPho] = p4_vtx0.Phi();
          patPhoton_scEnergy[iPho] = scRef->energy(); 
          patPhoton_scRawEnergy[iPho] = scRef->rawEnergy(); 
          patPhoton_scRawESEnergy[iPho] = scRef->preshowerEnergy(); 
          patPhoton_scSeedEnergy[iPho] = scRef->seed()->energy(); 
          patPhoton_scDPhiSeed[iPho] = reco::deltaPhi(scRef->seed()->phi(),scRef->position().Phi());
          patPhoton_scDEtaSeed[iPho] = scRef->seed()->eta() - scRef->position().Eta();
          patPhoton_scIEtaOrX[iPho] = scIsEB ? EBDetId(scRef->seed()->seed().rawId()).ieta() : EEDetId(scRef->seed()->seed().rawId()).ix(); 
          patPhoton_scIPhiOrY[iPho] = scIsEB ? EBDetId(scRef->seed()->seed().rawId()).iphi() : EEDetId(scRef->seed()->seed().rawId()).iy(); 
          patPhoton_scIEtaMod5[iPho] = scIsEB ? iEtaCorr%5 : -999;
          patPhoton_scIPhiMod2[iPho] = scIsEB ? (EBDetId(scRef->seed()->seed().rawId()).iphi()-1)%2 : -999;
          patPhoton_scIEtaMod20[iPho] = scIsEB ? (std::abs(EBDetId(scRef->seed()->seed().rawId()).ieta())<=25 ? iEtaCorr%20 : iEtaCorr26%20) : -999;
          patPhoton_scIPhiMod20[iPho] = scIsEB ? (EBDetId(scRef->seed()->seed().rawId()).iphi()-1)%20 : -999;
          patPhoton_scEt[iPho] = scRef->energy()*(Rt/R);
          patPhoton_scPhiWidth[iPho] = scRef->phiWidth();
          patPhoton_scEtaWidth[iPho] = scRef->etaWidth(); 
          patPhoton_SwissCross[iPho] = swissCross; 
          patPhoton_r9[iPho] = e3x3/scRef->rawEnergy(); 
          patPhoton_full5x5_r9[iPho] = full5x5_e3x3/scRef->rawEnergy();
          patPhoton_full5x5_e3x3[iPho] = full5x5_e3x3;
          patPhoton_full5x5_e5x5[iPho] = full5x5_e5x5;
          patPhoton_full5x5_eMax[iPho] = full5x5_eMax;
          patPhoton_full5x5_e2nd[iPho] = full5x5_e2nd;
          patPhoton_full5x5_eTop[iPho] = full5x5_eTop;
          patPhoton_full5x5_eBottom[iPho] = full5x5_eBottom;
          patPhoton_full5x5_eLeft[iPho] = full5x5_eLeft;
          patPhoton_full5x5_eRight[iPho] = full5x5_eRight;
          patPhoton_full5x5_e2x5Max[iPho] = full5x5_e2x5Max;
          patPhoton_full5x5_e2x5Left[iPho] = full5x5_e2x5Left;
          patPhoton_full5x5_e2x5Right[iPho] = full5x5_e2x5Right;
          patPhoton_full5x5_e2x5Top[iPho] = full5x5_e2x5Top;
          patPhoton_full5x5_e2x5Bottom[iPho] = full5x5_e2x5Bottom;
          patPhoton_full5x5_sigmaIEtaIEta[iPho] = full5x5_sigmaIEtaIEta;
          patPhoton_full5x5_sigmaIEtaIPhi[iPho] = full5x5_sigmaIEtaIPhi;
          patPhoton_full5x5_sigmaIEtaIPhi[iPho] = full5x5_phoSS.sigmaIetaIphi;
          patPhoton_full5x5_sigmaIPhiIPhi[iPho] = full5x5_sigmaIPhiIPhi;   
          patPhoton_full5x5_effSigmaRR[iPho] = full5x5_phoSS.effSigmaRR;
          patPhoton_ecalPFClusterIso[iPho] = iPhoton.ecalPFClusterIso();
          patPhoton_hcalPFClusterIso[iPho] = iPhoton.hcalPFClusterIso();
          patPhoton_trkSumPtHollowConeDR03[iPho] = iPhoton.trkSumPtHollowConeDR03(); 
          patPhoton_trkSumPtSolidConeDR04[iPho] = iPhoton.trkSumPtSolidConeDR04();
          patPhoton_chargedHadronIso[iPho] = iPhoton.chargedHadronIso();
          patPhoton_chargedHadronWorstVtxIso[iPho] = iPhoton.chargedHadronWorstVtxIso(); 
          patPhoton_egmMVAPhotonIDScore[iPho] = iPhoton.userFloat("PhotonMVAEstimatorRunIIFall17v2Values"); 
          patPhoton_egmCutBasedPhotonIDloose[iPho] = iPhoton.photonID(egmCutBasedPhotonIDloose_.c_str());
          patPhoton_egmCutBasedPhotonIDmedium[iPho] = iPhoton.photonID(egmCutBasedPhotonIDmedium_.c_str());
          patPhoton_egmCutBasedPhotonIDtight[iPho] = iPhoton.photonID(egmCutBasedPhotonIDtight_.c_str());
          patPhoton_egmMVAPhotonIDmedium[iPho] = iPhoton.photonID(egmMVAPhotonIDmedium_.c_str());
          patPhoton_egmMVAPhotonIDtight[iPho] = iPhoton.photonID(egmMVAPhotonIDtight_.c_str());    
          patPhoton_passPreselections[iPho] = passPreselections(&iPhoton);
          
          iPho++; 
      }
   }
   
   //fill tree for each event
   tree->Fill();
}

void ElePhoDumper::beginJob()
{

}

void ElePhoDumper::endJob() 
{
    

}

///------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
void ElePhoDumper::setTree(TTree* tree)
{
   tree->Branch("eventId", &eventId, "eventId/L");
   tree->Branch("lumiId", &lumiId, "lumiId/I");
   tree->Branch("runId", &runId, "runId/I");
   tree->Branch("rho", &rho, "rho/D"); 
   tree->Branch("nVtx", &nVtx, "nVtx/I");
   if(isMC_){ 
      tree->Branch("truePU", &truePU, "truePU/D");
      tree->Branch("obsPU", &obsPU, "obsPU/D");
      tree->Branch("genParticle_size", &genParticle_size, "genParticle_size/I");  
      tree->Branch("genParticle_index","std::vector<int>",&genParticle_index);
      tree->Branch("genParticle_pdgId","std::vector<int>",&genParticle_pdgId);
      tree->Branch("genParticle_status","std::vector<int>",&genParticle_status); 
      tree->Branch("genParticle_energy","std::vector<double>",&genParticle_energy);
      tree->Branch("genParticle_pt","std::vector<double>",&genParticle_pt);
      tree->Branch("genParticle_eta","std::vector<double>",&genParticle_eta);
      tree->Branch("genParticle_phi","std::vector<double>",&genParticle_phi);
      tree->Branch("genParticle_mass","std::vector<double>",&genParticle_mass);
      tree->Branch("genParticle_charge","std::vector<int>",&genParticle_charge);
      tree->Branch("genParticle_isLHEMatched","std::vector<int>",&genParticle_isLHEMatched);
   }
   if(saveElectrons_){
      tree->Branch("patElectron_size", &patElectron_size, "patElectron_size/I");  
      tree->Branch("patElectron_index","std::vector<int>",&patElectron_index); 
      if(isMC_){ 
         tree->Branch("patElectron_genIndex","std::vector<int>",&patElectron_genIndex); 
         tree->Branch("patElectron_genDR","std::vector<double>",&patElectron_genDR); 
         tree->Branch("patElectron_genFSREnergy","std::vector<double>",&patElectron_genFSREnergy); 
      }   
      tree->Branch("patElectron_classification","std::vector<int>",&patElectron_classification); 
      tree->Branch("patElectron_nrSatCrys","std::vector<int>",&patElectron_nrSatCrys); 
      tree->Branch("patElectron_scNPFClusters","std::vector<int>",&patElectron_scNPFClusters); 
      tree->Branch("patElectron_isEB","std::vector<int> ",&patElectron_isEB); 
      tree->Branch("patElectron_isEE","std::vector<int> ",&patElectron_isEE); 
      tree->Branch("patElectron_isEBEEGap","std::vector<int> ",&patElectron_isEBEEGap);
      tree->Branch("patElectron_isEBEtaGap","std::vector<int> ",&patElectron_isEBEtaGap);
      tree->Branch("patElectron_isEBPhiGap","std::vector<int> ",&patElectron_isEBPhiGap);
      tree->Branch("patElectron_isEEDeeGap","std::vector<int> ",&patElectron_isEEDeeGap);
      tree->Branch("patElectron_isEERingGap","std::vector<int> ",&patElectron_isEERingGap);
      tree->Branch("patElectron_isEcalDriven","std::vector<int> ",&patElectron_isEcalDriven);
      tree->Branch("patElectron_isTrackerDriven","std::vector<int> ",&patElectron_isTrackerDriven);
      tree->Branch("patElectron_passConversionVeto","std::vector<int> ",&patElectron_passConversionVeto);
      tree->Branch("patElectron_eta","std::vector<double>",&patElectron_eta); 
      tree->Branch("patElectron_phi","std::vector<double>",&patElectron_phi); 
      tree->Branch("patElectron_p","std::vector<double>",&patElectron_p); 
      tree->Branch("patElectron_pt","std::vector<double>",&patElectron_pt); 
      tree->Branch("patElectron_pIn","std::vector<double>",&patElectron_pIn); 
      tree->Branch("patElectron_pOut","std::vector<double>",&patElectron_pOut); 
      tree->Branch("patElectron_trackFbrem","std::vector<double>",&patElectron_trackFbrem);
      tree->Branch("patElectron_superClusterFbrem","std::vector<double>",&patElectron_superClusterFbrem);
      tree->Branch("patElectron_energy","std::vector<double>",&patElectron_energy);   
      tree->Branch("patElectron_energyErr","std::vector<double>",&patElectron_energyErr);   
      tree->Branch("patElectron_ecalEnergy","std::vector<double>",&patElectron_ecalEnergy); 
      tree->Branch("patElectron_ecalEnergyErr","std::vector<double>",&patElectron_ecalEnergyErr);  
      tree->Branch("patElectron_ecalEnergyPhoReg","std::vector<double>",&patElectron_ecalEnergyPhoReg); 
      tree->Branch("patElectron_ecalEnergyPhoRegErr","std::vector<double>",&patElectron_ecalEnergyPhoRegErr);  
      tree->Branch("patElectron_et","std::vector<double>",&patElectron_et); 
      tree->Branch("patElectron_HoE","std::vector<double>",&patElectron_HoE);    
      tree->Branch("patElectron_scIsEB","std::vector<int>",&patElectron_scIsEB); 
      tree->Branch("patElectron_scEta","std::vector<double>",&patElectron_scEta); 
      tree->Branch("patElectron_scPhi","std::vector<double>",&patElectron_scPhi);
      tree->Branch("patElectron_scEta_vtx0","std::vector<double>",&patElectron_scEta_vtx0); 
      tree->Branch("patElectron_scPhi_vtx0","std::vector<double>",&patElectron_scPhi_vtx0);
      tree->Branch("patElectron_scEnergy","std::vector<double>",&patElectron_scEnergy); 
      tree->Branch("patElectron_scRawEnergy","std::vector<double>",&patElectron_scRawEnergy); 
      tree->Branch("patElectron_scRawESEnergy","std::vector<double>",&patElectron_scRawESEnergy); 
      tree->Branch("patElectron_scSeedEnergy","std::vector<double>",&patElectron_scSeedEnergy); 
      tree->Branch("patElectron_scDPhiSeed","std::vector<double>",&patElectron_scDPhiSeed); 
      tree->Branch("patElectron_scDEtaSeed","std::vector<double>",&patElectron_scDEtaSeed); 
      tree->Branch("patElectron_scIEtaOrX","std::vector<int>",&patElectron_scIEtaOrX); 
      tree->Branch("patElectron_scIPhiOrY","std::vector<int>",&patElectron_scIPhiOrY); 
      tree->Branch("patElectron_scIEtaMod5","std::vector<int>",&patElectron_scIEtaMod5); 
      tree->Branch("patElectron_scIPhiMod2","std::vector<int>",&patElectron_scIPhiMod2); 
      tree->Branch("patElectron_scIEtaMod20","std::vector<int>",&patElectron_scIEtaMod20); 
      tree->Branch("patElectron_scIPhiMod20","std::vector<int>",&patElectron_scIPhiMod20); 
      tree->Branch("patElectron_scEt","std::vector<double>",&patElectron_scEt); 
      tree->Branch("patElectron_scPhiWidth","std::vector<double>",&patElectron_scPhiWidth);  
      tree->Branch("patElectron_scEtaWidth","std::vector<double>",&patElectron_scEtaWidth);  
      tree->Branch("patElectron_scEoP","std::vector<double>",&patElectron_scEoP);  
      tree->Branch("patElectron_SwissCross","std::vector<double>",&patElectron_SwissCross);  
      tree->Branch("patElectron_r9","std::vector<double>",&patElectron_r9);  
      tree->Branch("patElectron_full5x5_r9","std::vector<double>",&patElectron_full5x5_r9);  
      tree->Branch("patElectron_full5x5_e3x3","std::vector<double>",&patElectron_full5x5_e3x3);  
      tree->Branch("patElectron_full5x5_e5x5","std::vector<double>",&patElectron_full5x5_e5x5);  
      tree->Branch("patElectron_full5x5_eMax","std::vector<double>",&patElectron_full5x5_eMax);  
      tree->Branch("patElectron_full5x5_e2nd","std::vector<double>",&patElectron_full5x5_e2nd);  
      tree->Branch("patElectron_full5x5_eTop","std::vector<double>",&patElectron_full5x5_eTop);  
      tree->Branch("patElectron_full5x5_eBottom","std::vector<double>",&patElectron_full5x5_eBottom);  
      tree->Branch("patElectron_full5x5_eLeft","std::vector<double>",&patElectron_full5x5_eLeft);  
      tree->Branch("patElectron_full5x5_eRight","std::vector<double>",&patElectron_full5x5_eRight);  
      tree->Branch("patElectron_full5x5_e2x5Max","std::vector<double>",&patElectron_full5x5_e2x5Max); 
      tree->Branch("patElectron_full5x5_e2x5Top","std::vector<double>",&patElectron_full5x5_e2x5Top);  
      tree->Branch("patElectron_full5x5_e2x5Bottom[","std::vector<double>",&patElectron_full5x5_e2x5Bottom);  
      tree->Branch("patElectron_full5x5_e2x5Left","std::vector<double>",&patElectron_full5x5_e2x5Left);  
      tree->Branch("patElectron_full5x5_e2x5Right","std::vector<double>",&patElectron_full5x5_e2x5Right);  
      tree->Branch("patElectron_full5x5_sigmaIEtaIEta","std::vector<double>",&patElectron_full5x5_sigmaIEtaIEta);  
      tree->Branch("patElectron_full5x5_sigmaIEtaIPhi","std::vector<double>",&patElectron_full5x5_sigmaIEtaIPhi);  
      tree->Branch("patElectron_full5x5_sigmaIPhiIPhi","std::vector<double>",&patElectron_full5x5_sigmaIPhiIPhi);  
      tree->Branch("patElectron_egmCutBasedElectronIDVeto","std::vector<int>",&patElectron_egmCutBasedElectronIDVeto);  
      tree->Branch("patElectron_egmCutBasedElectronIDloose","std::vector<int>",&patElectron_egmCutBasedElectronIDloose);  
      tree->Branch("patElectron_egmCutBasedElectronIDmedium","std::vector<int>",&patElectron_egmCutBasedElectronIDmedium);
      tree->Branch("patElectron_egmCutBasedElectronIDtight","std::vector<int>",&patElectron_egmCutBasedElectronIDtight); 
      tree->Branch("patElectron_egmMVAElectronIDScore","std::vector<double>",&patElectron_egmMVAElectronIDScore);   
      tree->Branch("patElectron_egmMVAElectronIDNoIsoScore","std::vector<double>",&patElectron_egmMVAElectronIDNoIsoScore);   
      tree->Branch("patElectron_egmMVAElectronIDloose","std::vector<int>",&patElectron_egmMVAElectronIDloose);  
      tree->Branch("patElectron_egmMVAElectronIDmedium","std::vector<int>",&patElectron_egmMVAElectronIDmedium);
      tree->Branch("patElectron_egmMVAElectronIDtight","std::vector<int>",&patElectron_egmMVAElectronIDtight);  
      tree->Branch("patElectron_egmMVAElectronIDlooseNoIso","std::vector<int>",&patElectron_egmMVAElectronIDlooseNoIso);  
      tree->Branch("patElectron_egmMVAElectronIDmediumNoIso","std::vector<int>",&patElectron_egmMVAElectronIDmediumNoIso);
      tree->Branch("patElectron_egmMVAElectronIDtightNoIso","std::vector<int>",&patElectron_egmMVAElectronIDtightNoIso);  
      tree->Branch("patElectron_heepElectronID","std::vector<int>",&patElectron_heepElectronID);  
   }
   if(savePhotons_ || saveElectrons_){
      tree->Branch("patPhoton_size", &patPhoton_size, "patPhoton_size/I");  
      tree->Branch("patPhoton_index","std::vector<int>",&patPhoton_index); 
      if(isMC_){ 
         tree->Branch("patPhoton_genIndex","std::vector<int>",&patPhoton_genIndex); 
         tree->Branch("patPhoton_genDR","std::vector<double>",&patPhoton_genDR); 
      } 
      tree->Branch("patPhoton_nrSatCrys","std::vector<int>",&patPhoton_nrSatCrys); 
      tree->Branch("patPhoton_scNPFClusters","std::vector<int>",&patPhoton_scNPFClusters);   
      tree->Branch("patPhoton_passElectronVeto","std::vector<int> ",&patPhoton_passElectronVeto); 
      tree->Branch("patPhoton_hasPixelSeed","std::vector<int> ",&patPhoton_hasPixelSeed); 
      tree->Branch("patPhoton_hasConversionTracks","std::vector<int> ",&patPhoton_hasConversionTracks); 
      tree->Branch("patPhoton_nConversions","std::vector<int>",&patPhoton_nConversions);     
      tree->Branch("patPhoton_nConversionsOneLeg","std::vector<int>",&patPhoton_nConversionsOneLeg);     
      tree->Branch("patPhoton_isEB","std::vector<int> ",&patPhoton_isEB); 
      tree->Branch("patPhoton_isEE","std::vector<int> ",&patPhoton_isEE); 
      tree->Branch("patPhoton_isEBEEGap","std::vector<int> ",&patPhoton_isEBEEGap);
      tree->Branch("patPhoton_isEBEtaGap","std::vector<int> ",&patPhoton_isEBEtaGap);
      tree->Branch("patPhoton_isEBPhiGap","std::vector<int> ",&patPhoton_isEBPhiGap);
      tree->Branch("patPhoton_isEEDeeGap","std::vector<int> ",&patPhoton_isEEDeeGap);
      tree->Branch("patPhoton_isEERingGap","std::vector<int> ",&patPhoton_isEERingGap);  
      tree->Branch("patPhoton_eta","std::vector<double>",&patPhoton_eta); 
      tree->Branch("patPhoton_phi","std::vector<double>",&patPhoton_phi); 
      tree->Branch("patPhoton_energy","std::vector<double>",&patPhoton_energy);   
      tree->Branch("patPhoton_energyErr","std::vector<double>",&patPhoton_energyErr);   
      tree->Branch("patPhoton_ecalEnergy","std::vector<double>",&patPhoton_ecalEnergy); 
      tree->Branch("patPhoton_ecalEnergyErr","std::vector<double>",&patPhoton_ecalEnergyErr);  
      tree->Branch("patPhoton_et","std::vector<double>",&patPhoton_et); 
      tree->Branch("patPhoton_pt","std::vector<double>",&patPhoton_pt); 
      tree->Branch("patPhoton_HoE","std::vector<double>",&patPhoton_HoE);  
      tree->Branch("patPhoton_scIsEB","std::vector<int>",&patPhoton_scIsEB); 
      tree->Branch("patPhoton_scEta","std::vector<double>",&patPhoton_scEta); 
      tree->Branch("patPhoton_scPhi","std::vector<double>",&patPhoton_scPhi);  
      tree->Branch("patPhoton_scEta_vtx0","std::vector<double>",&patPhoton_scEta_vtx0); 
      tree->Branch("patPhoton_scPhi_vtx0","std::vector<double>",&patPhoton_scPhi_vtx0);  
      tree->Branch("patPhoton_scEnergy","std::vector<double>",&patPhoton_scEnergy); 
      tree->Branch("patPhoton_scRawEnergy","std::vector<double>",&patPhoton_scRawEnergy); 
      tree->Branch("patPhoton_scRawESEnergy","std::vector<double>",&patPhoton_scRawESEnergy); 
      tree->Branch("patPhoton_scSeedEnergy","std::vector<double>",&patPhoton_scSeedEnergy); 
      tree->Branch("patPhoton_scDPhiSeed","std::vector<double>",&patPhoton_scDPhiSeed); 
      tree->Branch("patPhoton_scDEtaSeed","std::vector<double>",&patPhoton_scDEtaSeed); 
      tree->Branch("patPhoton_scIEtaOrX","std::vector<int>",&patPhoton_scIEtaOrX); 
      tree->Branch("patPhoton_scIPhiOrY","std::vector<int>",&patPhoton_scIPhiOrY); 
      tree->Branch("patPhoton_scIEtaMod5","std::vector<int>",&patPhoton_scIEtaMod5); 
      tree->Branch("patPhoton_scIPhiMod2","std::vector<int>",&patPhoton_scIPhiMod2); 
      tree->Branch("patPhoton_scIEtaMod20","std::vector<int>",&patPhoton_scIEtaMod20); 
      tree->Branch("patPhoton_scIPhiMod20","std::vector<int>",&patPhoton_scIPhiMod20); 
      tree->Branch("patPhoton_scEt","std::vector<double>",&patPhoton_scEt); 
      tree->Branch("patPhoton_scPhiWidth","std::vector<double>",&patPhoton_scPhiWidth);  
      tree->Branch("patPhoton_scEtaWidth","std::vector<double>",&patPhoton_scEtaWidth);  
      tree->Branch("patPhoton_SwissCross","std::vector<double>",&patPhoton_SwissCross);  
      tree->Branch("patPhoton_r9","std::vector<double>",&patPhoton_r9);  
      tree->Branch("patPhoton_full5x5_r9","std::vector<double>",&patPhoton_full5x5_r9);  
      tree->Branch("patPhoton_full5x5_e3x3","std::vector<double>",&patPhoton_full5x5_e3x3);  
      tree->Branch("patPhoton_full5x5_e5x5","std::vector<double>",&patPhoton_full5x5_e5x5);  
      tree->Branch("patPhoton_full5x5_eMax","std::vector<double>",&patPhoton_full5x5_eMax);  
      tree->Branch("patPhoton_full5x5_e2nd","std::vector<double>",&patPhoton_full5x5_e2nd);  
      tree->Branch("patPhoton_full5x5_eTop","std::vector<double>",&patPhoton_full5x5_eTop);  
      tree->Branch("patPhoton_full5x5_eBottom","std::vector<double>",&patPhoton_full5x5_eBottom);  
      tree->Branch("patPhoton_full5x5_eLeft","std::vector<double>",&patPhoton_full5x5_eLeft);  
      tree->Branch("patPhoton_full5x5_eRight","std::vector<double>",&patPhoton_full5x5_eRight);  
      tree->Branch("patPhoton_full5x5_e2x5Max","std::vector<double>",&patPhoton_full5x5_e2x5Max); 
      tree->Branch("patPhoton_full5x5_e2x5Top","std::vector<double>",&patPhoton_full5x5_e2x5Top);  
      tree->Branch("patPhoton_full5x5_e2x5Bottom[","std::vector<double>",&patPhoton_full5x5_e2x5Bottom);  
      tree->Branch("patPhoton_full5x5_e2x5Left","std::vector<double>",&patPhoton_full5x5_e2x5Left);  
      tree->Branch("patPhoton_full5x5_e2x5Right","std::vector<double>",&patPhoton_full5x5_e2x5Right);  
      tree->Branch("patPhoton_full5x5_sigmaIEtaIEta","std::vector<double>",&patPhoton_full5x5_sigmaIEtaIEta);  
      tree->Branch("patPhoton_full5x5_sigmaIEtaIPhi","std::vector<double>",&patPhoton_full5x5_sigmaIEtaIPhi);  
      tree->Branch("patPhoton_full5x5_sigmaIPhiIPhi","std::vector<double>",&patPhoton_full5x5_sigmaIPhiIPhi);   
      tree->Branch("patPhoton_full5x5_effSigmaRR","std::vector<double>",&patPhoton_full5x5_effSigmaRR);   
      tree->Branch("patPhoton_ecalPFClusterIso","std::vector<double>",&patPhoton_ecalPFClusterIso);   
      tree->Branch("patPhoton_hcalPFClusterIso","std::vector<double>",&patPhoton_hcalPFClusterIso);   
      tree->Branch("patPhoton_trkSumPtHollowConeDR03","std::vector<double>",&patPhoton_trkSumPtHollowConeDR03);   
      tree->Branch("patPhoton_trkSumPtSolidConeDR04","std::vector<double>",&patPhoton_trkSumPtSolidConeDR04);   
      tree->Branch("patPhoton_chargedHadronIso","std::vector<double>",&patPhoton_chargedHadronIso);   
      tree->Branch("patPhoton_chargedHadronWorstVtxIso","std::vector<double>",&patPhoton_chargedHadronWorstVtxIso);   
      tree->Branch("patPhoton_egmCutBasedPhotonIDloose","std::vector<int>",&patPhoton_egmCutBasedPhotonIDloose);  
      tree->Branch("patPhoton_egmCutBasedPhotonIDmedium","std::vector<int>",&patPhoton_egmCutBasedPhotonIDmedium);
      tree->Branch("patPhoton_egmCutBasedPhotonIDtight","std::vector<int>",&patPhoton_egmCutBasedPhotonIDtight);  
      tree->Branch("patPhoton_egmMVAPhotonIDScore","std::vector<double>",&patPhoton_egmMVAPhotonIDScore); 
      tree->Branch("patPhoton_egmMVAPhotonIDmedium","std::vector<int>",&patPhoton_egmMVAPhotonIDmedium);
      tree->Branch("patPhoton_egmMVAPhotonIDtight","std::vector<int>",&patPhoton_egmMVAPhotonIDtight);  
      tree->Branch("patPhoton_passPreselections","std::vector<int>",&patPhoton_passPreselections);
   }
}

void ElePhoDumper::setVectors(int nGenParts, int nElectrons, int nPhotons)
{
   genParticle_index.clear();
   genParticle_pdgId.clear();
   genParticle_status.clear();
   genParticle_energy.clear();
   genParticle_pt.clear();
   genParticle_eta.clear();
   genParticle_phi.clear(); 
   genParticle_mass.clear();
   genParticle_charge.clear();
   genParticle_isLHEMatched.clear();
   genParticle_index.resize(nGenParts);
   setDefaultValues(genParticle_index,-1);
   genParticle_pdgId.resize(nGenParts);
   setDefaultValues(genParticle_pdgId,-999);
   genParticle_status.resize(nGenParts);
   setDefaultValues(genParticle_status,-999);
   genParticle_energy.resize(nGenParts);
   setDefaultValues(genParticle_energy,-999.);
   genParticle_pt.resize(nGenParts);
   setDefaultValues(genParticle_pt,-999.);
   genParticle_eta.resize(nGenParts);
   setDefaultValues(genParticle_eta,-999.);
   genParticle_phi.resize(nGenParts);
   setDefaultValues(genParticle_phi,-999.);
   genParticle_mass.resize(nGenParts);
   setDefaultValues(genParticle_mass,-999.);
   genParticle_charge.resize(nGenParts);
   setDefaultValues(genParticle_charge,-999);
   genParticle_isLHEMatched.resize(nGenParts);
   setDefaultValues(genParticle_isLHEMatched,-999);
    
   patElectron_index.clear();
   patElectron_genIndex.clear();
   patElectron_genDR.clear();
   patElectron_genFSREnergy.clear();
   patElectron_classification.clear();
   patElectron_nrSatCrys.clear();
   patElectron_scNPFClusters.clear();
   patElectron_charge.clear();
   patElectron_isEB.clear();
   patElectron_isEE.clear();  
   patElectron_isEBEEGap.clear();
   patElectron_isEBEtaGap.clear();
   patElectron_isEBPhiGap.clear();
   patElectron_isEEDeeGap.clear();
   patElectron_isEERingGap.clear(); 
   patElectron_isEcalDriven.clear();
   patElectron_isTrackerDriven.clear();
   patElectron_passConversionVeto.clear();
   patElectron_eta.clear();
   patElectron_phi.clear();
   patElectron_p.clear();
   patElectron_pt.clear();
   patElectron_pIn.clear();
   patElectron_pOut.clear();
   patElectron_trackFbrem.clear();
   patElectron_superClusterFbrem.clear();
   patElectron_energy.clear();
   patElectron_energyErr.clear();
   patElectron_ecalEnergy.clear();
   patElectron_ecalEnergyErr.clear();
   patElectron_ecalEnergyPhoReg.clear();
   patElectron_ecalEnergyPhoRegErr.clear();
   patElectron_et.clear();
   patElectron_HoE.clear();
   patElectron_scIsEB.clear();
   patElectron_scEta.clear();
   patElectron_scPhi.clear();
   patElectron_scEta_vtx0.clear();
   patElectron_scPhi_vtx0.clear();
   patElectron_scEnergy.clear();
   patElectron_scRawEnergy.clear();
   patElectron_scRawESEnergy.clear();
   patElectron_scSeedEnergy.clear();
   patElectron_scDPhiSeed.clear();
   patElectron_scDEtaSeed.clear(); 
   patElectron_scIEtaOrX.clear(); 
   patElectron_scIPhiOrY.clear(); 
   patElectron_scIEtaMod5.clear(); 
   patElectron_scIPhiMod2.clear(); 
   patElectron_scIEtaMod20.clear(); 
   patElectron_scIPhiMod20.clear(); 
   patElectron_scEt.clear();
   patElectron_scPhiWidth.clear();
   patElectron_scEtaWidth.clear();
   patElectron_scEoP.clear();
   patElectron_SwissCross.clear();
   patElectron_r9.clear();
   patElectron_full5x5_r9.clear();
   patElectron_full5x5_e3x3.clear();
   patElectron_full5x5_e5x5.clear();
   patElectron_full5x5_eMax.clear();
   patElectron_full5x5_e2nd.clear();
   patElectron_full5x5_eTop.clear();
   patElectron_full5x5_eBottom.clear();
   patElectron_full5x5_eLeft.clear();
   patElectron_full5x5_eRight.clear();
   patElectron_full5x5_e2x5Max.clear();
   patElectron_full5x5_e2x5Left.clear();
   patElectron_full5x5_e2x5Right.clear();
   patElectron_full5x5_e2x5Top.clear();
   patElectron_full5x5_e2x5Bottom.clear();
   patElectron_full5x5_sigmaIEtaIEta.clear();
   patElectron_full5x5_sigmaIEtaIPhi.clear();
   patElectron_full5x5_sigmaIPhiIPhi.clear();
   patElectron_egmCutBasedElectronIDVeto.clear();
   patElectron_egmCutBasedElectronIDloose.clear();
   patElectron_egmCutBasedElectronIDmedium.clear();
   patElectron_egmCutBasedElectronIDtight.clear();
   patElectron_egmMVAElectronIDScore.clear();
   patElectron_egmMVAElectronIDNoIsoScore.clear();
   patElectron_egmMVAElectronIDloose.clear();
   patElectron_egmMVAElectronIDmedium.clear();
   patElectron_egmMVAElectronIDtight.clear();
   patElectron_egmMVAElectronIDlooseNoIso.clear();
   patElectron_egmMVAElectronIDmediumNoIso.clear();
   patElectron_egmMVAElectronIDtightNoIso.clear();
   patElectron_heepElectronID.clear();
   if(saveElectrons_){
      patElectron_index.resize(nElectrons);
      setDefaultValues(patElectron_index,-1);
      patElectron_genIndex.resize(nElectrons);
      setDefaultValues(patElectron_genIndex,-1);
      patElectron_genDR.resize(nElectrons);
      setDefaultValues(patElectron_genDR,-999.);
      patElectron_genFSREnergy.resize(nElectrons);
      setDefaultValues(patElectron_genFSREnergy,-999.);
      patElectron_classification.resize(nElectrons);
      setDefaultValues(patElectron_classification,-999);
      patElectron_nrSatCrys.resize(nElectrons);
      setDefaultValues(patElectron_nrSatCrys,-1);
      patElectron_scNPFClusters.resize(nElectrons);
      setDefaultValues(patElectron_scNPFClusters,-1);
      patElectron_charge.resize(nElectrons);
      setDefaultValues(patElectron_charge,-999);
      patElectron_isEB.resize(nElectrons);
      setDefaultValues(patElectron_isEB,-1);
      patElectron_isEE.resize(nElectrons);  
      setDefaultValues(patElectron_isEE,-1);
      patElectron_isEBEEGap.resize(nElectrons);
      setDefaultValues(patElectron_isEBEEGap,-1);
      patElectron_isEBEtaGap.resize(nElectrons);
      setDefaultValues(patElectron_isEBEtaGap,-1);
      patElectron_isEBPhiGap.resize(nElectrons);
      setDefaultValues(patElectron_isEBPhiGap,-1);
      patElectron_isEEDeeGap.resize(nElectrons);
      setDefaultValues(patElectron_isEEDeeGap,-1);
      patElectron_isEERingGap.resize(nElectrons); 
      setDefaultValues(patElectron_isEERingGap,-1);
      patElectron_isEcalDriven.resize(nElectrons);
      setDefaultValues(patElectron_isEcalDriven,-1);
      patElectron_isTrackerDriven.resize(nElectrons);
      setDefaultValues(patElectron_isTrackerDriven,-1);
      patElectron_passConversionVeto.resize(nElectrons);
      setDefaultValues(patElectron_passConversionVeto,-1);
      patElectron_eta.resize(nElectrons);
      setDefaultValues(patElectron_eta,-999.);
      patElectron_phi.resize(nElectrons);
      setDefaultValues(patElectron_phi,-999.);
      patElectron_p.resize(nElectrons);
      setDefaultValues(patElectron_p,-999.);
      patElectron_pt.resize(nElectrons);
      setDefaultValues(patElectron_pt,-999.);
      patElectron_pIn.resize(nElectrons);
      setDefaultValues(patElectron_pIn,-999.);
      patElectron_pOut.resize(nElectrons);
      setDefaultValues(patElectron_pOut,-999.);
      patElectron_trackFbrem.resize(nElectrons);
      setDefaultValues(patElectron_trackFbrem,-999.);
      patElectron_superClusterFbrem.resize(nElectrons);
      setDefaultValues(patElectron_superClusterFbrem,-999.);
      patElectron_energy.resize(nElectrons);
      setDefaultValues(patElectron_energy,-999.);
      patElectron_energyErr.resize(nElectrons);
      setDefaultValues(patElectron_energyErr,-999.);
      patElectron_ecalEnergy.resize(nElectrons);
      setDefaultValues(patElectron_ecalEnergy,-999.);
      patElectron_ecalEnergyErr.resize(nElectrons);
      setDefaultValues(patElectron_ecalEnergyErr,-999.);
      patElectron_ecalEnergyPhoReg.resize(nElectrons);
      setDefaultValues(patElectron_ecalEnergyPhoReg,-999.);
      patElectron_ecalEnergyPhoRegErr.resize(nElectrons);
      setDefaultValues(patElectron_ecalEnergyPhoRegErr,-999.);
      patElectron_et.resize(nElectrons);
      setDefaultValues(patElectron_et,-999.);
      patElectron_HoE.resize(nElectrons);
      setDefaultValues(patElectron_HoE,-999.);
      patElectron_scIsEB.resize(nElectrons);
      setDefaultValues(patElectron_scIsEB,-1);
      patElectron_scEta.resize(nElectrons);
      setDefaultValues(patElectron_scEta,-999.);
      patElectron_scPhi.resize(nElectrons);
      setDefaultValues(patElectron_scPhi,-999.);
      patElectron_scEta_vtx0.resize(nElectrons);
      setDefaultValues(patElectron_scEta_vtx0,-999.); 
      patElectron_scPhi_vtx0.resize(nElectrons);
      setDefaultValues(patElectron_scPhi_vtx0,-999.); 
      patElectron_scEnergy.resize(nElectrons);
      setDefaultValues(patElectron_scEnergy,-999.); 
      patElectron_scRawEnergy.resize(nElectrons);
      setDefaultValues(patElectron_scRawEnergy,-999.); 
      patElectron_scRawESEnergy.resize(nElectrons);
      setDefaultValues(patElectron_scRawESEnergy,-999.);
      patElectron_scSeedEnergy.resize(nElectrons);
      setDefaultValues(patElectron_scSeedEnergy,-999.);
      patElectron_scDPhiSeed.resize(nElectrons);
      setDefaultValues(patElectron_scDPhiSeed,-999.);
      patElectron_scDEtaSeed.resize(nElectrons);
      setDefaultValues(patElectron_scDEtaSeed,-999.);
      patElectron_scIEtaOrX.resize(nElectrons);
      setDefaultValues(patElectron_scIEtaOrX,-999);
      patElectron_scIPhiOrY.resize(nElectrons);
      setDefaultValues(patElectron_scIPhiOrY,-999);
      patElectron_scIEtaMod5.resize(nElectrons);
      setDefaultValues(patElectron_scIEtaMod5,-999);
      patElectron_scIPhiMod2.resize(nElectrons);
      setDefaultValues(patElectron_scIPhiMod2,-999);
      patElectron_scIEtaMod20.resize(nElectrons);
      setDefaultValues(patElectron_scIEtaMod20,-999);
      patElectron_scIPhiMod20.resize(nElectrons);
      setDefaultValues(patElectron_scIPhiMod20,-999);
      patElectron_scEt.resize(nElectrons);
      setDefaultValues(patElectron_scEt,-999.);
      patElectron_scPhiWidth.resize(nElectrons);
      setDefaultValues(patElectron_scPhiWidth,-999.);
      patElectron_scEtaWidth.resize(nElectrons);
      setDefaultValues(patElectron_scEtaWidth,-999.);
      patElectron_scEoP.resize(nElectrons);
      setDefaultValues(patElectron_scEoP,-999.);
      patElectron_SwissCross.resize(nElectrons);
      setDefaultValues(patElectron_SwissCross,-999.);
      patElectron_r9.resize(nElectrons);
      setDefaultValues(patElectron_r9,-999.);
      patElectron_full5x5_r9.resize(nElectrons);
      setDefaultValues(patElectron_full5x5_r9,-999.);
      patElectron_full5x5_e3x3.resize(nElectrons);
      setDefaultValues(patElectron_full5x5_e3x3,-999.);
      patElectron_full5x5_e5x5.resize(nElectrons);
      setDefaultValues(patElectron_full5x5_e5x5,-999.);
      patElectron_full5x5_eMax.resize(nElectrons);
      setDefaultValues(patElectron_full5x5_eMax,-999.);
      patElectron_full5x5_e2nd.resize(nElectrons);
      setDefaultValues(patElectron_full5x5_e2nd,-999.);
      patElectron_full5x5_eTop.resize(nElectrons);
      setDefaultValues(patElectron_full5x5_eTop,-999.);
      patElectron_full5x5_eBottom.resize(nElectrons);
      setDefaultValues(patElectron_full5x5_eBottom,-999.);
      patElectron_full5x5_eLeft.resize(nElectrons);
      setDefaultValues(patElectron_full5x5_eLeft,-999.);
      patElectron_full5x5_eRight.resize(nElectrons);
      setDefaultValues(patElectron_full5x5_eRight,-999.);
      patElectron_full5x5_e2x5Max.resize(nElectrons);
      setDefaultValues(patElectron_full5x5_e2x5Max,-999.);
      patElectron_full5x5_e2x5Left.resize(nElectrons);
      setDefaultValues(patElectron_full5x5_e2x5Left,-999.);
      patElectron_full5x5_e2x5Right.resize(nElectrons);
      setDefaultValues(patElectron_full5x5_e2x5Right,-999.);
      patElectron_full5x5_e2x5Top.resize(nElectrons);
      setDefaultValues(patElectron_full5x5_e2x5Top,-999.);
      patElectron_full5x5_e2x5Bottom.resize(nElectrons);
      setDefaultValues(patElectron_full5x5_e2x5Bottom,-999.);
      patElectron_full5x5_sigmaIEtaIEta.resize(nElectrons);
      setDefaultValues(patElectron_full5x5_sigmaIEtaIEta,-999.);
      patElectron_full5x5_sigmaIEtaIPhi.resize(nElectrons);
      setDefaultValues(patElectron_full5x5_sigmaIEtaIPhi,-999.);
      patElectron_full5x5_sigmaIPhiIPhi.resize(nElectrons);
      setDefaultValues(patElectron_full5x5_sigmaIPhiIPhi,-999.);
      patElectron_egmCutBasedElectronIDVeto.resize(nElectrons);
      setDefaultValues(patElectron_egmCutBasedElectronIDVeto,-1);
      patElectron_egmCutBasedElectronIDloose.resize(nElectrons);
      setDefaultValues(patElectron_egmCutBasedElectronIDloose,-1);
      patElectron_egmCutBasedElectronIDmedium.resize(nElectrons);
      setDefaultValues(patElectron_egmCutBasedElectronIDmedium,-1);
      patElectron_egmCutBasedElectronIDtight.resize(nElectrons);
      setDefaultValues(patElectron_egmCutBasedElectronIDtight,-1);
      patElectron_egmMVAElectronIDScore.resize(nElectrons);
      setDefaultValues(patElectron_egmMVAElectronIDScore,-999.);
      patElectron_egmMVAElectronIDNoIsoScore.resize(nElectrons);
      setDefaultValues(patElectron_egmMVAElectronIDNoIsoScore,-999.);
      patElectron_egmMVAElectronIDloose.resize(nElectrons);
      setDefaultValues(patElectron_egmMVAElectronIDloose,-1);
      patElectron_egmMVAElectronIDmedium.resize(nElectrons);
      setDefaultValues(patElectron_egmMVAElectronIDmedium,-1);
      patElectron_egmMVAElectronIDtight.resize(nElectrons);
      setDefaultValues(patElectron_egmMVAElectronIDtight,-1);
      patElectron_egmMVAElectronIDlooseNoIso.resize(nElectrons);
      setDefaultValues(patElectron_egmMVAElectronIDlooseNoIso,-1);
      patElectron_egmMVAElectronIDmediumNoIso.resize(nElectrons);
      setDefaultValues(patElectron_egmMVAElectronIDmediumNoIso,-1);
      patElectron_egmMVAElectronIDtightNoIso.resize(nElectrons);
      setDefaultValues(patElectron_egmMVAElectronIDtightNoIso,-1);
      patElectron_heepElectronID.resize(nElectrons);
      setDefaultValues(patElectron_heepElectronID,-1);
   }

   patPhoton_index.clear();
   patPhoton_genIndex.clear();
   patPhoton_genDR.clear();
   patPhoton_nrSatCrys.clear();
   patPhoton_scNPFClusters.clear();
   patPhoton_isEB.clear();
   patPhoton_isEE.clear();
   patPhoton_isEBEEGap.clear(); 
   patPhoton_isEBEtaGap.clear();
   patPhoton_isEBPhiGap.clear();
   patPhoton_isEEDeeGap.clear();
   patPhoton_isEERingGap.clear(); 
   patPhoton_passElectronVeto.clear();
   patPhoton_hasPixelSeed.clear(); 
   patPhoton_hasConversionTracks.clear();
   patPhoton_nConversions.clear();
   patPhoton_nConversionsOneLeg.clear();
   patPhoton_eta.clear();
   patPhoton_phi.clear();
   patPhoton_energy.clear();
   patPhoton_energyErr.clear();
   patPhoton_ecalEnergy.clear();
   patPhoton_ecalEnergyErr.clear();
   patPhoton_et.clear();
   patPhoton_pt.clear();
   patPhoton_HoE.clear();
   patPhoton_scIsEB.clear();
   patPhoton_scEta.clear();
   patPhoton_scPhi.clear();
   patPhoton_scEta_vtx0.clear();
   patPhoton_scPhi_vtx0.clear();
   patPhoton_scEnergy.clear();
   patPhoton_scRawEnergy.clear();
   patPhoton_scRawESEnergy.clear();
   patPhoton_scSeedEnergy.clear();
   patPhoton_scDPhiSeed.clear();
   patPhoton_scDEtaSeed.clear();
   patPhoton_scIEtaOrX.clear();
   patPhoton_scIPhiOrY.clear();
   patPhoton_scIEtaMod5.clear();
   patPhoton_scIPhiMod2.clear();
   patPhoton_scIEtaMod20.clear();
   patPhoton_scIPhiMod20.clear();
   patPhoton_scEt.clear();
   patPhoton_scPhiWidth.clear();
   patPhoton_scEtaWidth.clear();
   patPhoton_SwissCross.clear();
   patPhoton_r9.clear();
   patPhoton_full5x5_r9.clear();
   patPhoton_full5x5_e3x3.clear();
   patPhoton_full5x5_e5x5.clear();
   patPhoton_full5x5_eMax.clear();
   patPhoton_full5x5_e2nd.clear();
   patPhoton_full5x5_eTop.clear();
   patPhoton_full5x5_eBottom.clear();
   patPhoton_full5x5_eLeft.clear();
   patPhoton_full5x5_eRight.clear();
   patPhoton_full5x5_e2x5Max.clear();
   patPhoton_full5x5_e2x5Left.clear();
   patPhoton_full5x5_e2x5Right.clear();
   patPhoton_full5x5_e2x5Top.clear();
   patPhoton_full5x5_e2x5Bottom.clear();
   patPhoton_full5x5_sigmaIEtaIEta.clear();
   patPhoton_full5x5_sigmaIEtaIPhi.clear();
   patPhoton_full5x5_sigmaIPhiIPhi.clear();
   patPhoton_full5x5_effSigmaRR.clear();
   patPhoton_ecalPFClusterIso.clear();
   patPhoton_hcalPFClusterIso.clear();
   patPhoton_trkSumPtHollowConeDR03.clear();
   patPhoton_trkSumPtSolidConeDR04.clear();
   patPhoton_chargedHadronIso.clear();
   patPhoton_chargedHadronWorstVtxIso.clear();
   patPhoton_egmCutBasedPhotonIDloose.clear();
   patPhoton_egmCutBasedPhotonIDmedium.clear();
   patPhoton_egmCutBasedPhotonIDtight.clear();
   patPhoton_egmMVAPhotonIDScore.clear();
   patPhoton_egmMVAPhotonIDmedium.clear();
   patPhoton_egmMVAPhotonIDtight.clear();
   patPhoton_passPreselections.clear();
   if(savePhotons_ || saveElectrons_){
      patPhoton_index.resize(nPhotons);
      setDefaultValues(patPhoton_index,-1);
      patPhoton_genIndex.resize(nPhotons);
      setDefaultValues(patPhoton_genIndex,-1);
      patPhoton_genDR.resize(nPhotons);
      setDefaultValues(patPhoton_genDR,-999.);
      patPhoton_nrSatCrys.resize(nPhotons);
      setDefaultValues(patPhoton_nrSatCrys,-1);
      patPhoton_scNPFClusters.resize(nPhotons);
      setDefaultValues(patPhoton_scNPFClusters,-1);
      patPhoton_isEB.resize(nPhotons);
      setDefaultValues(patPhoton_isEB,-1);
      patPhoton_isEE.resize(nPhotons);
      setDefaultValues(patPhoton_isEE,-1);
      patPhoton_isEBEEGap.resize(nPhotons); 
      setDefaultValues(patPhoton_isEBEEGap,-1);
      patPhoton_isEBEtaGap.resize(nPhotons);
      setDefaultValues(patPhoton_isEBEtaGap,-1);
      patPhoton_isEBPhiGap.resize(nPhotons);
      setDefaultValues(patPhoton_isEBPhiGap,-1);
      patPhoton_isEEDeeGap.resize(nPhotons);
      setDefaultValues(patPhoton_isEEDeeGap,-1);
      patPhoton_isEERingGap.resize(nPhotons); 
      setDefaultValues(patPhoton_isEERingGap,-1);
      patPhoton_passElectronVeto.resize(nPhotons);
      setDefaultValues(patPhoton_passElectronVeto,-1);
      patPhoton_hasPixelSeed.resize(nPhotons); 
      setDefaultValues(patPhoton_hasPixelSeed,-1);
      patPhoton_hasConversionTracks.resize(nPhotons);
      setDefaultValues(patPhoton_hasConversionTracks,-1);
      patPhoton_nConversions.resize(nPhotons);
      setDefaultValues(patPhoton_nConversions,-1);
      patPhoton_nConversionsOneLeg.resize(nPhotons);
      setDefaultValues(patPhoton_nConversionsOneLeg,-1);
      patPhoton_eta.resize(nPhotons);
      setDefaultValues(patPhoton_eta,-999.);
      patPhoton_phi.resize(nPhotons);
      setDefaultValues(patPhoton_phi,-999.);
      patPhoton_energy.resize(nPhotons);
      setDefaultValues(patPhoton_energy,-999.);
      patPhoton_energyErr.resize(nPhotons);
      setDefaultValues(patPhoton_energyErr,-999.);
      patPhoton_ecalEnergy.resize(nPhotons);
      setDefaultValues(patPhoton_ecalEnergy,-999.);
      patPhoton_ecalEnergyErr.resize(nPhotons);
      setDefaultValues(patPhoton_ecalEnergyErr,-999.);
      patPhoton_et.resize(nPhotons);
      setDefaultValues(patPhoton_et,-999.);
      patPhoton_pt.resize(nPhotons);
      setDefaultValues(patPhoton_pt,-999.);
      patPhoton_HoE.resize(nPhotons);
      setDefaultValues(patPhoton_HoE,-999.);
      patPhoton_scIsEB.resize(nPhotons);
      setDefaultValues(patPhoton_scIsEB,-1);
      patPhoton_scEta.resize(nPhotons);
      setDefaultValues(patPhoton_scEta,-999.);
      patPhoton_scPhi.resize(nPhotons);
      setDefaultValues(patPhoton_scPhi,-999.);
      patPhoton_scEta_vtx0.resize(nPhotons);
      setDefaultValues(patPhoton_scEta_vtx0,-999.);
      patPhoton_scPhi_vtx0.resize(nPhotons);
      setDefaultValues(patPhoton_scPhi_vtx0,-999.);
      patPhoton_scEnergy.resize(nPhotons);
      setDefaultValues(patPhoton_scEnergy,-999.);
      patPhoton_scRawEnergy.resize(nPhotons);
      setDefaultValues(patPhoton_scRawEnergy,-999.);
      patPhoton_scRawESEnergy.resize(nPhotons);
      setDefaultValues(patPhoton_scRawESEnergy,-999.);
      patPhoton_scSeedEnergy.resize(nPhotons);
      setDefaultValues(patPhoton_scRawESEnergy,-999.);
      patPhoton_scDPhiSeed.resize(nPhotons);
      setDefaultValues(patPhoton_scDPhiSeed,-999.);
      patPhoton_scDEtaSeed.resize(nPhotons);
      setDefaultValues(patPhoton_scDEtaSeed,-999.);
      patPhoton_scIEtaOrX.resize(nPhotons);
      setDefaultValues(patPhoton_scIEtaOrX,-999);
      patPhoton_scIPhiOrY.resize(nPhotons);
      setDefaultValues(patPhoton_scIPhiOrY,-999);
      patPhoton_scIEtaMod5.resize(nPhotons);
      setDefaultValues(patPhoton_scIEtaMod5,-999);
      patPhoton_scIPhiMod2.resize(nPhotons);
      setDefaultValues(patPhoton_scIPhiMod2,-999);
      patPhoton_scIEtaMod20.resize(nPhotons);
      setDefaultValues(patPhoton_scIEtaMod20,-999);
      patPhoton_scIPhiMod20.resize(nPhotons);
      setDefaultValues(patPhoton_scIPhiMod20,-999);
      patPhoton_scEt.resize(nPhotons);
      setDefaultValues(patPhoton_scEt,-999.);
      patPhoton_scPhiWidth.resize(nPhotons);
      setDefaultValues(patPhoton_scPhiWidth,-999.);
      patPhoton_scEtaWidth.resize(nPhotons);
      setDefaultValues(patPhoton_scEtaWidth,-999.);
      patPhoton_SwissCross.resize(nPhotons);
      setDefaultValues(patPhoton_SwissCross,-999.);
      patPhoton_r9.resize(nPhotons);
      setDefaultValues(patPhoton_r9,-999.);
      patPhoton_full5x5_r9.resize(nPhotons);
      setDefaultValues(patPhoton_full5x5_r9,-999.);
      patPhoton_full5x5_e3x3.resize(nPhotons);
      setDefaultValues(patPhoton_full5x5_e3x3,-999.);
      patPhoton_full5x5_e5x5.resize(nPhotons);
      setDefaultValues(patPhoton_full5x5_e5x5,-999.);
      patPhoton_full5x5_eMax.resize(nPhotons);
      setDefaultValues(patPhoton_full5x5_eMax,-999.);
      patPhoton_full5x5_e2nd.resize(nPhotons);
      setDefaultValues(patPhoton_full5x5_e2nd,-999.);
      patPhoton_full5x5_eTop.resize(nPhotons);
      setDefaultValues(patPhoton_full5x5_eTop,-999.);
      patPhoton_full5x5_eBottom.resize(nPhotons);
      setDefaultValues(patPhoton_full5x5_eBottom,-999.);
      patPhoton_full5x5_eLeft.resize(nPhotons);
      setDefaultValues(patPhoton_full5x5_eLeft,-999.);
      patPhoton_full5x5_eRight.resize(nPhotons);
      setDefaultValues(patPhoton_full5x5_eRight,-999.);
      patPhoton_full5x5_e2x5Max.resize(nPhotons);
      setDefaultValues(patPhoton_full5x5_e2x5Max,-999.);
      patPhoton_full5x5_e2x5Left.resize(nPhotons);
      setDefaultValues(patPhoton_full5x5_e2x5Left,-999.);
      patPhoton_full5x5_e2x5Right.resize(nPhotons);
      setDefaultValues(patPhoton_full5x5_e2x5Right,-999.);
      patPhoton_full5x5_e2x5Top.resize(nPhotons);
      setDefaultValues(patPhoton_full5x5_e2x5Top,-999.);
      patPhoton_full5x5_e2x5Bottom.resize(nPhotons);
      setDefaultValues(patPhoton_full5x5_e2x5Bottom,-999.);
      patPhoton_full5x5_sigmaIEtaIEta.resize(nPhotons);
      setDefaultValues(patPhoton_full5x5_sigmaIEtaIEta,-999.);
      patPhoton_full5x5_sigmaIEtaIPhi.resize(nPhotons);
      setDefaultValues(patPhoton_full5x5_sigmaIEtaIPhi,-999.);
      patPhoton_full5x5_sigmaIPhiIPhi.resize(nPhotons);
      setDefaultValues(patPhoton_full5x5_sigmaIPhiIPhi,-999.);
      patPhoton_full5x5_effSigmaRR.resize(nPhotons);
      setDefaultValues(patPhoton_full5x5_effSigmaRR,-999.);
      patPhoton_ecalPFClusterIso.resize(nPhotons);
      setDefaultValues(patPhoton_ecalPFClusterIso,-999.);
      patPhoton_hcalPFClusterIso.resize(nPhotons);
      setDefaultValues(patPhoton_hcalPFClusterIso,-999.);
      patPhoton_trkSumPtHollowConeDR03.resize(nPhotons);
      setDefaultValues(patPhoton_trkSumPtHollowConeDR03,-999.);
      patPhoton_trkSumPtSolidConeDR04.resize(nPhotons);
      setDefaultValues(patPhoton_trkSumPtSolidConeDR04,-999.);
      patPhoton_chargedHadronIso.resize(nPhotons);
      setDefaultValues(patPhoton_chargedHadronIso,-999.);
      patPhoton_chargedHadronWorstVtxIso.resize(nPhotons);
      setDefaultValues(patPhoton_chargedHadronWorstVtxIso,-999.);
      patPhoton_egmCutBasedPhotonIDloose.resize(nPhotons);
      setDefaultValues(patPhoton_egmCutBasedPhotonIDloose,-1);
      patPhoton_egmCutBasedPhotonIDmedium.resize(nPhotons);
      setDefaultValues(patPhoton_egmCutBasedPhotonIDmedium,-1);
      patPhoton_egmCutBasedPhotonIDtight.resize(nPhotons);
      setDefaultValues(patPhoton_egmCutBasedPhotonIDtight,-1);
      patPhoton_egmMVAPhotonIDScore.resize(nPhotons);
      setDefaultValues(patPhoton_egmMVAPhotonIDScore,-999.);
      patPhoton_egmMVAPhotonIDmedium.resize(nPhotons);
      setDefaultValues(patPhoton_egmMVAPhotonIDmedium,-1);
      patPhoton_egmMVAPhotonIDtight.resize(nPhotons);
      setDefaultValues(patPhoton_egmMVAPhotonIDtight,-1);
      patPhoton_passPreselections.resize(nPhotons);
      setDefaultValues(patPhoton_passPreselections,-1);
   }
   
}

template <typename T>
void ElePhoDumper::setDefaultValues(std::vector<T>& vec, const T& value)
{
    std::fill(vec.begin(), vec.end(), value);
}

const reco::GenParticle* ElePhoDumper::getFirstMother( const reco::GenParticle* particle, int pdgId, int status ) 
{
    if (!particle) return nullptr;
    
    const reco::GenParticle* current = particle;
    const reco::GenParticle* first = particle;
    while (current->numberOfMothers() > 0) 
    {
           const reco::GenParticle* mother = dynamic_cast<const reco::GenParticle*>(current->motherRef().get());
           if(!mother){ 
              return nullptr;
           }   
           if(status==-999){
              //if(debug_) std::cout << "getFirstMother - Mother: pdgId = " << mother->pdgId() << " - status = " << mother->status() << std::endl;
              if(abs(mother->pdgId())==pdgId) first = mother;
           }else{
              //if(debug_) std::cout << "getFirstMother - Mother: pdgId = " << mother->pdgId() << " - status = " << mother->status() << std::endl;
              if(mother->status()==status && abs(mother->pdgId())==pdgId) first = mother;
           }
           current = mother;  // Move up the ancestry
    } 

    // return first matched mother
    return first;
}

int ElePhoDumper::getGenIndex( edm::Handle<std::vector<reco::GenParticle> > genParticles, const reco::GenParticle* particle ) 
{
    int index = -1;
    if(!particle || particle==nullptr) return -1;
    for(size_t i = 0; i < genParticles->size(); ++i) 
    {
        if(&genParticles->at(i) == particle) 
        {
           index = i;
           break;
        }
    }   
    return index;  
}

void ElePhoDumper::printAllDaughters( const reco::GenParticle* particle, int depth ) 
{
   // Base case: safety check
   if(!particle) return;
   if(depth==0){
      std::cout << "Mother PDG ID: " << particle->pdgId()
                << ", status: " << particle->status()
                << ", energy: " << particle->energy() 
                << ", pt: " << particle->pt() 
                << ", eta: " << particle->eta() 
                << ", phi: " << particle->phi() << std::endl;
   }
   
   // Loop over direct daughters
   for(unsigned int i = 0; i < particle->numberOfDaughters(); ++i) 
   {
       const reco::Candidate* dau = particle->daughterRef(i).get();

       // Indentation for visual hierarchy
       std::cout << std::string(depth * 2, ' ')
                 << "Daughter PDG ID: " << dau->pdgId()
                 << ", status: " << dau->status()
                 << ", energy: " << dau->energy() 
                 << ", pt: " << dau->pt() 
                 << ", eta: " << dau->eta() 
                 << ", phi: " << dau->phi() << std::endl;

       // Try to cast to reco::GenParticle for recursive call
       const reco::GenParticle* genDau = dynamic_cast<const reco::GenParticle*>(dau);
       if(genDau) {
          printAllDaughters(genDau, depth + 1);  // recursive call
       }
   }
}

math::XYZTLorentzVector ElePhoDumper::setP4(const reco::SuperClusterRef scRef, const reco::Vertex& vtx, double energy)
{
    double sc_X = scRef->x();
    double sc_Y = scRef->y();
    double sc_Z = scRef->z();
    
    double vtx_X = vtx.x();
    double vtx_Y = vtx.y();
    double vtx_Z = vtx.z();
    
    //std::cout << "vtx = ("<< vtx_X << "," << vtx_Y << "," << vtx_Z << ")" << std::endl;
    
    math::XYZVector vtx_Pos( vtx_X, vtx_Y, vtx_Z );
    math::XYZVector sc_Pos( sc_X, sc_Y, sc_Z );

    math::XYZVector direction = sc_Pos - vtx_Pos;
    math::XYZVector p = ( direction.Unit() ) * ( energy );
    math::XYZTLorentzVector corrected_p4( p.x(), p.y(), p.z(), energy );
    
    return corrected_p4;
}

int ElePhoDumper::passPreselections(const pat::Photon* photon)
{
   int pass = 0;

   double pt = photon->pt();
   double scEta = photon->superCluster()->eta(); 
   double r9 = photon->full5x5_r9();
   double hoe = photon->hadronicOverEm(); 
   double sieie = photon->full5x5_sigmaIetaIeta();
   double pfPhoIso03 = photon->photonIso(); 
   double trackIso03 = photon->trkSumPtHollowConeDR03(); 
   double egCHIso = photon->chargedHadronIso();
   bool passEleVeto = photon->passElectronVeto();

   if( (( abs(scEta)<1.479 && r9>0.5 && r9>0.85 ) || 
        ( abs(scEta)<1.479 && r9>0.5 && r9<=0.85 && sieie<0.015 && pfPhoIso03<4.0 && trackIso03<6.0 ) ||
        ( abs(scEta)>1.479 && r9>0.8 && r9>0.90 ) ||
        ( abs(scEta)>1.479 && r9>0.8 && r9<=0.90 && sieie<0.035 && pfPhoIso03<4.0 && trackIso03<6.0 )) 
       && (( r9>0.8 && egCHIso<20. ) || ( egCHIso/pt<0.3 && pt>14. && hoe<0.15 ))
       && ( passEleVeto>0 && hoe<0.08 ) 
   ) pass = 1;
   
   return pass;
}   

///------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"
DEFINE_FWK_MODULE(ElePhoDumper);

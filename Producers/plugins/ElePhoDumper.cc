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
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
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
#include "TLorentzVector.h"
#include "TMath.h"
#include "TTree.h"

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
      const reco::GenParticle* getFirstMother( const reco::GenParticle* particle, int pdgId, int status );
      void printAllDaughters( const reco::GenParticle* particle, int depth ); 
      math::XYZTLorentzVector setP4(const reco::SuperClusterRef scRef, const reco::Vertex& vtx, double energy);
      int passPreselections(const pat::Photon* photon);
      
      // ----------collection tokens-------------------
      edm::EDGetTokenT<double> rhoToken_;
      edm::EDGetTokenT<std::vector<PileupSummaryInfo> > pileupSummaryToken_;
      edm::EDGetTokenT<reco::VertexCollection> vtxToken_; 
      edm::EDGetTokenT<std::vector<reco::GenParticle> > genToken_; 
      edm::EDGetTokenT<EcalRecHitCollection> ebRechitToken_; 
      edm::EDGetTokenT<EcalRecHitCollection> eeRechitToken_; 
      edm::EDGetTokenT<std::vector<pat::Electron> > patElectronToken_;
      edm::EDGetTokenT<std::vector<pat::Photon> > patPhotonToken_;
      
      edm::Service<TFileService> iFile;
      
      // ----------config inputs-------------------
      bool isMC_;
      std::vector<int> motherPdgId_;
      std::vector<int> motherStatus_;
      std::vector<int> pdgIdBeforeFSR_;
      std::vector<int> pdgId_;
      std::vector<int> status_;
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
      const reco::Vertex* vtx0;
      math::XYZTLorentzVector p4_vtx0;
      
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
      int patElectron_size;
      std::vector<int> patElectron_index; 
      std::vector<int> patElectron_genIndex; 
      std::vector<double> patElectron_genDR; 
      std::vector<double> patElectron_genFSREnergy;
      std::vector<double> patElectron_genFSRNergy;
      std::vector<int> patElectron_classification;
      std::vector<int> patElectron_scNPFClusters;
      std::vector<int> patElectron_charge;
      std::vector<bool> patElectron_isEB;
      std::vector<bool> patElectron_isEE;
      std::vector<bool> patElectron_isEBEEGap;
      std::vector<bool> patElectron_isEBEtaGap;
      std::vector<bool> patElectron_isEBPhiGap;
      std::vector<bool> patElectron_isEEDeeGap;
      std::vector<bool> patElectron_isEERingGap;
      std::vector<bool> patElectron_isEcalDriven;
      std::vector<bool> patElectron_isTrackerDriven;
      std::vector<bool> patElectron_passConversionVeto;
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
      std::vector<double> patElectron_et;
      std::vector<double> patElectron_scEnergy;
      std::vector<double> patElectron_scRawEnergy;
      std::vector<double> patElectron_scRawESEnergy;
      std::vector<double> patElectron_scEt;
      std::vector<double> patElectron_scPhiWidth;
      std::vector<double> patElectron_scEtaWidth;
      std::vector<double> patElectron_scEoP;
      std::vector<double> patElectron_scEta;
      std::vector<double> patElectron_scPhi;
      std::vector<double> patElectron_scEta_vtx0;
      std::vector<double> patElectron_scPhi_vtx0;
      std::vector<double> patElectron_scSwissCross;
      std::vector<double> patElectron_scEMax;
      std::vector<double> patElectron_scR9; 
      std::vector<double> patElectron_scSigmaIEtaIEta;
      std::vector<double> patElectron_scSigmaIEtaIPhi;
      std::vector<double> patElectron_scSigmaIPhiIPhi;
      std::vector<double> patElectron_full5x5_scEMax;
      std::vector<double> patElectron_full5x5_scR9;  
      std::vector<double> patElectron_full5x5_scSigmaIEtaIEta;
      std::vector<double> patElectron_full5x5_scSigmaIEtaIPhi;
      std::vector<double> patElectron_full5x5_scSigmaIPhiIPhi;
      std::vector<double> patElectron_HoE;
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
      std::vector<int> patPhoton_scNPFClusters;
      std::vector<bool> patPhoton_isEB;
      std::vector<bool> patPhoton_isEE;
      std::vector<bool> patPhoton_isEBEEGap;
      std::vector<bool> patPhoton_isEBEtaGap;
      std::vector<bool> patPhoton_isEBPhiGap;
      std::vector<bool> patPhoton_isEEDeeGap;
      std::vector<bool> patPhoton_isEERingGap;
      std::vector<bool> patPhoton_passElectronVeto;
      std::vector<bool> patPhoton_hasPixelSeed;
      std::vector<bool> patPhoton_hasConversionTracks;
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
      std::vector<double> patPhoton_scEnergy;  
      std::vector<double> patPhoton_scRawEnergy;  
      std::vector<double> patPhoton_scRawESEnergy;
      std::vector<double> patPhoton_scEt;
      std::vector<double> patPhoton_scEtaWidth;
      std::vector<double> patPhoton_scPhiWidth;    
      std::vector<double> patPhoton_scEta;
      std::vector<double> patPhoton_scPhi;
      std::vector<double> patPhoton_scEta_vtx0;
      std::vector<double> patPhoton_scPhi_vtx0;
      std::vector<double> patPhoton_scSwissCross;
      std::vector<double> patPhoton_scE2x2;
      std::vector<double> patPhoton_scE3x3; 
      std::vector<double> patPhoton_scE5x5; 
      std::vector<double> patPhoton_scEMax;
      std::vector<double> patPhoton_scR9;
      std::vector<double> patPhoton_scSigmaIEtaIEta;
      std::vector<double> patPhoton_scSigmaIEtaIPhi;
      std::vector<double> patPhoton_scSigmaIPhiIPhi;
      std::vector<double> patPhoton_full5x5_scEMax;
      std::vector<double> patPhoton_full5x5_scR9;
      std::vector<double> patPhoton_full5x5_scSigmaIEtaIEta;
      std::vector<double> patPhoton_full5x5_scSigmaIEtaIPhi;
      std::vector<double> patPhoton_full5x5_scSigmaIPhiIPhi;
      std::vector<double> patPhoton_HoE;
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
   ebRechitToken_                 = consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("ebRechitCollection"));
   eeRechitToken_                 = consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("eeRechitCollection"));
   patElectronToken_              = consumes<std::vector<pat::Electron> >(iConfig.getParameter<edm::InputTag>("patElectronCollection"));
   patPhotonToken_                = consumes<std::vector<pat::Photon> >(iConfig.getParameter<edm::InputTag>("patPhotonCollection")); 
   
   isMC_                          = iConfig.getParameter<bool>("isMC"); 
   motherPdgId_                   = iConfig.getParameter<std::vector<int> >( "motherPdgId" );
   motherStatus_                  = iConfig.getParameter<std::vector<int> >( "motherStatus" );
   pdgIdBeforeFSR_                = iConfig.getParameter<std::vector<int> >( "pdgIdBeforeFSR" );
   pdgId_                         = iConfig.getParameter<std::vector<int> >( "pdgId" );
   status_                        = iConfig.getParameter<std::vector<int> >( "status" );
   ignoreTauDecays_               = iConfig.getParameter<bool>( "ignoreTauDecays" );
   dRMin_                         = iConfig.getParameter<double>( "dRMin" );
    
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
   setVectors(genParticle_size, patElectron_size,patPhoton_size);
   
   genIndices.clear();
   double FSREnergy_tot = -999;
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
          
          for(unsigned int i = 0; i <pdgId_.size(); ++i)
          {
              if(iGenPart.pdgId()==pdgId_.at(i) && iGenPart.status()==status_.at(i))
              {   
                 if(ignoreTauDecays_ && iGenPart.isDirectPromptTauDecayProductFinalState()==1) continue; 
               
                  for(unsigned int j = 0; j <pdgIdBeforeFSR_.size(); ++j)
                  {        
                      const reco::GenParticle* motherPart = getFirstMother(&iGenPart,pdgIdBeforeFSR_.at(j),-999); //get the particle before the FSR
                      if(motherPart==nullptr) continue;
                      
                      bool goodMotherID=false;
                      bool goodMotherStatus=false;
                      for(unsigned int k=0; k<motherPdgId_.size(); k++){
                          if(motherPdgId_.at(k)!=-999 && motherPart->motherRef().get()->pdgId()==motherPdgId_.at(k)) goodMotherID=true; 
                          else goodMotherID=true;     
                          if(motherStatus_.at(k)!=-999 && motherPart->motherRef().get()->status()==motherStatus_.at(k)) goodMotherStatus=true;   
                          else goodMotherStatus=true;  
                      }    
                      if(!(goodMotherID && goodMotherStatus)) continue;
                      std::cout << "genIndex: " << iGen << " - pdgId: " << iGenPart.pdgId() << " - status = " << iGenPart.status() << " - mother-pdgId = " << motherPart->pdgId() << " - mother-status = " << motherPart->status() << std::endl;
                      genIndices.push_back(iGen);
                  }    
              }   
          } 
      }
      
      //Remove duplicates
      std::sort( genIndices.begin(), genIndices.end() );
      genIndices.erase( unique( genIndices.begin(), genIndices.end() ), genIndices.end() );
      
      std::cout << "GenParticles trees... " << std::endl;
      for(unsigned int i=0; i<genIndices.size(); i++){
          const reco::GenParticle* motherPart = dynamic_cast<const reco::GenParticle*>((&genParticles->at(genIndices.at(i)))->motherRef().get()); 
          printAllDaughters( motherPart, 0 );    
      }
    
      //check FSR
      iGen=-1;
      FSREnergy_tot = 0.;
      for(const auto& iGenPart : *(genParticles.product()))
      {
        iGen++;  
        for(unsigned int i = 0; i <pdgIdBeforeFSR_.size(); ++i)
        {
            if(iGenPart.pdgId()==22 && iGenPart.status()==1)
            { 
               const reco::GenParticle* motherPart = getFirstMother(&iGenPart,pdgIdBeforeFSR_.at(i),-999);
               if(motherPart==nullptr) continue;
               
               bool goodMotherID=false;
               bool goodMotherStatus=false;
               for(unsigned int k=0; k<motherPdgId_.size(); k++){
                   if(motherPdgId_.at(k)!=-999 && motherPart->motherRef().get()->pdgId()==motherPdgId_.at(k)) goodMotherID=true; 
                   else goodMotherID=true;     
                   if(motherStatus_.at(k)!=-999 && motherPart->motherRef().get()->status()==motherStatus_.at(k)) goodMotherStatus=true;   
                   else goodMotherStatus=true;  
               }    
               if(!(goodMotherID && goodMotherStatus)) continue;
               
               FSREnergy_tot += iGenPart.energy();
               if(saveElectrons_) std::cout << "FSR photon: p4 = (" << iGenPart.pt() << "," << iGenPart.eta() << "," << iGenPart.phi() << "," << iGenPart.energy() << ") - last-mother: pdgId = " << motherPart->pdgId() << " , status = " << motherPart->status() << " , p4 = (" << motherPart->pt() << "," << motherPart->eta() << "," << motherPart->phi() << "," << motherPart->energy() << ")" << std::endl;           
              
            }
        }
      } 
   }
   
   
          
   //save electron info
   if(saveElectrons_)
   { 
      int iEle=0;
      for(const auto& iElectron : *(patElectron.product())){ 
 
          reco::SuperClusterRef scRef = iElectron.superCluster();
          double R  = TMath::Sqrt(scRef->x()*scRef->x() + scRef->y()*scRef->y() +scRef->z()*scRef->z());
          double Rt = TMath::Sqrt(scRef->x()*scRef->x() + scRef->y()*scRef->y());
          
          double swissCross = -999.;
          double e3x3 = -999.;
          double eMax = -999.;
          double full5x5_e3x3 = -999.; 
          double full5x5_eMax = -999.; 
         
          reco::GsfElectron::ShowerShape eleSS = iElectron.showerShape();
          reco::GsfElectron::ShowerShape full5x5_eleSS = iElectron.full5x5_showerShape();
          const std::vector<std::pair<DetId,float> > &hits= iElectron.superCluster()->hitsAndFractions();
          if(iElectron.isEB())
          {
             std::pair<DetId, float> id = EcalClusterTools::getMaximum(hits,&(*(recHitsEB.product())));   
             swissCross = EcalTools::swissCross(id.first,*(recHitsEB.product()),0.);
             e3x3 = EcalClusterTools::e3x3( *(scRef->seed()), &(*(recHitsEB.product())), topology);
             eMax = EcalClusterTools::eMax( *(scRef->seed()), &(*(recHitsEB.product())));             
             full5x5_e3x3 = noZS::EcalClusterTools::e3x3( *(scRef->seed()), &(*(recHitsEB.product())), topology);
             full5x5_eMax = noZS::EcalClusterTools::eMax( *(scRef->seed()), &(*(recHitsEB.product())));
          }
          if(iElectron.isEE())
          {
             std::pair<DetId, float> id = EcalClusterTools::getMaximum(hits,&(*(recHitsEE.product())));   
             swissCross = EcalTools::swissCross(id.first,*(recHitsEE.product()),0.);
             e3x3 = EcalClusterTools::e3x3( *(scRef->seed()), &(*(recHitsEE.product())), topology);
             eMax = EcalClusterTools::eMax( *(scRef->seed()), &(*(recHitsEE.product())));
             full5x5_e3x3 = noZS::EcalClusterTools::e3x3( *(scRef->seed()), &(*(recHitsEE.product())), topology);
             full5x5_eMax = noZS::EcalClusterTools::eMax( *(scRef->seed()), &(*(recHitsEE.product())));
          }
          
          p4_vtx0 = setP4(scRef, *vtx0, iElectron.correctedEcalEnergy());
          
          double dR_tmp = dRMin_;
          double dR = 999.;
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
          
         if(genIndex>=0) std::cout << "genIndex = " << genIndex << " - dR = " << dR << " - pdgId = " << (&genParticles->at(genIndex))->pdgId() << " - status = " << (&genParticles->at(genIndex))->status() << " - eta = " << (&genParticles->at(genIndex))->eta() << " - phi = " << (&genParticles->at(genIndex))->phi() << " - pt = " << (&genParticles->at(genIndex))->pt() << " - energy = " << (&genParticles->at(genIndex))->energy() << " - electronEnergy = " << iElectron.energy() << std::endl; 
          
          patElectron_index[iEle] = iEle;
          patElectron_genIndex[iEle] = genIndex;
          patElectron_genDR[iEle] = dR;
          patElectron_genFSREnergy[iEle] = FSREnergy_tot;
          patElectron_classification[iEle] = iElectron.classification();
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
          patElectron_et[iEle] = iElectron.p4().Et();
          patElectron_HoE[iEle] = iElectron.hadronicOverEm(); 
          patElectron_scEta[iEle] = scRef->eta();
          patElectron_scPhi[iEle] = scRef->phi();
          patElectron_scEta_vtx0[iEle] = p4_vtx0.Eta();
          patElectron_scPhi_vtx0[iEle] = p4_vtx0.Phi();
          patElectron_scEnergy[iEle] = scRef->energy(); 
          patElectron_scRawEnergy[iEle] = scRef->rawEnergy(); 
          patElectron_scRawESEnergy[iEle] = scRef->preshowerEnergy(); 
          patElectron_scEt[iEle] = scRef->energy()*(Rt/R);
          patElectron_scPhiWidth[iEle] = scRef->phiWidth();
          patElectron_scEtaWidth[iEle] = scRef->etaWidth(); 
          patElectron_scEoP[iEle] = scRef->energy()/iElectron.trackMomentumAtVtx().R(); 
          patElectron_scSwissCross[iEle] = swissCross; 
          patElectron_scR9[iEle] = e3x3/scRef->rawEnergy(); 
          patElectron_scEMax[iEle] = eMax; 
          patElectron_scSigmaIEtaIEta[iEle] = eleSS.sigmaIetaIeta;
          patElectron_scSigmaIEtaIPhi[iEle] = eleSS.sigmaIetaIphi;
          patElectron_scSigmaIPhiIPhi[iEle] = eleSS.sigmaIphiIphi;  
          patElectron_full5x5_scR9[iEle] = full5x5_e3x3/scRef->rawEnergy();
          patElectron_full5x5_scEMax[iEle] = full5x5_eMax;
          patElectron_full5x5_scSigmaIEtaIEta[iEle] = full5x5_eleSS.sigmaIetaIeta;
          patElectron_full5x5_scSigmaIEtaIPhi[iEle] = full5x5_eleSS.sigmaIetaIphi;
          patElectron_full5x5_scSigmaIPhiIPhi[iEle] = full5x5_eleSS.sigmaIphiIphi;  
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
   if(savePhotons_)
   {    
      int iPho=0;
      for(const auto& iPhoton : *(patPhoton.product())){ 
 
          reco::SuperClusterRef scRef = iPhoton.superCluster();
          reco::PhotonCoreRef phoCoreRef = iPhoton.photonCore();
          double R  = TMath::Sqrt(scRef->x()*scRef->x() + scRef->y()*scRef->y() +scRef->z()*scRef->z());
          double Rt = TMath::Sqrt(scRef->x()*scRef->x() + scRef->y()*scRef->y());

          double swissCross = -999.;
          double e3x3 = -999.;
          double eMax = -999.;
          double full5x5_e3x3 = -999.; 
          double full5x5_eMax = -999.; 
         
          reco::Photon::ShowerShape phoSS = iPhoton.showerShapeVariables(); 
          reco::Photon::ShowerShape full5x5_phoSS = iPhoton.full5x5_showerShapeVariables();    
          const std::vector<std::pair<DetId,float> > &hits= iPhoton.superCluster()->hitsAndFractions();
          if(iPhoton.isEB())
          {
             std::pair<DetId, float> id = EcalClusterTools::getMaximum(hits,&(*(recHitsEB.product())));   

             swissCross = EcalTools::swissCross(id.first,*(recHitsEB.product()),0.);
             e3x3 = EcalClusterTools::e3x3( *(scRef->seed()), &(*(recHitsEB.product())), topology);
             eMax = EcalClusterTools::eMax( *(scRef->seed()), &(*(recHitsEB.product())));
             full5x5_e3x3 = noZS::EcalClusterTools::e3x3( *(scRef->seed()), &(*(recHitsEB.product())), topology);
             full5x5_eMax = noZS::EcalClusterTools::eMax( *(scRef->seed()), &(*(recHitsEB.product())));
          }
          if(iPhoton.isEE())
          {
             std::pair<DetId, float> id = EcalClusterTools::getMaximum(hits,&(*(recHitsEE.product())));   

             swissCross = EcalTools::swissCross(id.first,*(recHitsEE.product()),0.);
             e3x3 = EcalClusterTools::e3x3( *(scRef->seed()), &(*(recHitsEE.product())), topology);
             eMax = EcalClusterTools::eMax( *(scRef->seed()), &(*(recHitsEE.product())));
             full5x5_e3x3 = noZS::EcalClusterTools::e3x3( *(scRef->seed()), &(*(recHitsEE.product())), topology);
             full5x5_eMax = noZS::EcalClusterTools::eMax( *(scRef->seed()), &(*(recHitsEE.product())));
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
          
          if(genIndex>=0) std::cout << "genIndex = " << genIndex << " - dR = " << dR << " - pdgId = " << (&genParticles->at(genIndex))->pdgId() << " - status = " << (&genParticles->at(genIndex))->status() << " - eta = " << (&genParticles->at(genIndex))->eta() << " - phi = " << (&genParticles->at(genIndex))->phi() << " - pt = " << (&genParticles->at(genIndex))->pt() << " - energy = " << (&genParticles->at(genIndex))->energy() << " - photonEnergy = " << iPhoton.energy() << std::endl; 
          
          patPhoton_index[iPho] = iPho;
          patPhoton_genIndex[iPho] = genIndex;
          patPhoton_genDR[iPho] = dR;
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
          patPhoton_scEta[iPho] = scRef->eta();
          patPhoton_scPhi[iPho] = scRef->phi();
          patPhoton_scEta_vtx0[iPho] = p4_vtx0.Eta();
          patPhoton_scPhi_vtx0[iPho] = p4_vtx0.Phi();
          patPhoton_scEnergy[iPho] = scRef->energy(); 
          patPhoton_scRawEnergy[iPho] = scRef->rawEnergy(); 
          patPhoton_scRawESEnergy[iPho] = scRef->preshowerEnergy(); 
          patPhoton_scEt[iPho] = scRef->energy()*(Rt/R);
          patPhoton_scPhiWidth[iPho] = scRef->phiWidth();
          patPhoton_scEtaWidth[iPho] = scRef->etaWidth(); 
          patPhoton_scSwissCross[iPho] = swissCross; 
          patPhoton_scR9[iPho] = e3x3/scRef->rawEnergy(); 
          patPhoton_scEMax[iPho] = eMax; 
          patPhoton_scSigmaIEtaIEta[iPho] = phoSS.sigmaIetaIeta;
          patPhoton_scSigmaIEtaIPhi[iPho] = phoSS.sigmaIetaIphi;
          patPhoton_scSigmaIPhiIPhi[iPho] = phoSS.sigmaIphiIphi;   
          patPhoton_full5x5_scR9[iPho] = full5x5_e3x3/scRef->rawEnergy();  
          patPhoton_full5x5_scEMax[iPho] = full5x5_eMax; 
          patPhoton_full5x5_scSigmaIEtaIEta[iPho] = full5x5_phoSS.sigmaIetaIeta;
          patPhoton_full5x5_scSigmaIEtaIPhi[iPho] = full5x5_phoSS.sigmaIetaIphi;
          patPhoton_full5x5_scSigmaIPhiIPhi[iPho] = full5x5_phoSS.sigmaIphiIphi;  
          patPhoton_egmMVAPhotonIDScore[iPho]  = iPhoton.userFloat("PhotonMVAEstimatorRunIIFall17v2Values"); 
          patPhoton_egmCutBasedPhotonIDloose[iPho] = iPhoton.photonID(egmCutBasedPhotonIDloose_.c_str());
          patPhoton_egmCutBasedPhotonIDmedium[iPho] = iPhoton.photonID(egmCutBasedPhotonIDmedium_.c_str());
          patPhoton_egmCutBasedPhotonIDtight[iPho] = iPhoton.photonID(egmCutBasedPhotonIDtight_.c_str());
          patPhoton_egmMVAPhotonIDmedium[iPho] = iPhoton.photonID(egmMVAPhotonIDmedium_.c_str());
          patPhoton_egmMVAPhotonIDtight[iPho] = iPhoton.photonID(egmMVAPhotonIDtight_.c_str());    
          patPhoton_passPreselections[iPho] = passPreselections(&iPhoton);
          
          iPho++; 
      }
   }
   
   std::cout << "" << std::endl;
   
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
      tree->Branch("patElectron_scNPFClusters","std::vector<int>",&patElectron_scNPFClusters); 
      tree->Branch("patElectron_isEB","std::vector<bool> ",&patElectron_isEB); 
      tree->Branch("patElectron_isEE","std::vector<bool> ",&patElectron_isEE); 
      tree->Branch("patElectron_isEBEEGap","std::vector<bool> ",&patElectron_isEBEEGap);
      tree->Branch("patElectron_isEBEtaGap","std::vector<bool> ",&patElectron_isEBEtaGap);
      tree->Branch("patElectron_isEBPhiGap","std::vector<bool> ",&patElectron_isEBPhiGap);
      tree->Branch("patElectron_isEEDeeGap","std::vector<bool> ",&patElectron_isEEDeeGap);
      tree->Branch("patElectron_isEERingGap","std::vector<bool> ",&patElectron_isEERingGap);
      tree->Branch("patElectron_isEcalDriven","std::vector<bool> ",&patElectron_isEcalDriven);
      tree->Branch("patElectron_isTrackerDriven","std::vector<bool> ",&patElectron_isTrackerDriven);
      tree->Branch("patElectron_passConversionVeto","std::vector<bool> ",&patElectron_passConversionVeto);
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
      tree->Branch("patElectron_et","std::vector<double>",&patElectron_et); 
      tree->Branch("patElectron_HoE","std::vector<double>",&patElectron_HoE);    
      tree->Branch("patElectron_scEta","std::vector<double>",&patElectron_scEta); 
      tree->Branch("patElectron_scPhi","std::vector<double>",&patElectron_scPhi);
      tree->Branch("patElectron_scEta_vtx0","std::vector<double>",&patElectron_scEta_vtx0); 
      tree->Branch("patElectron_scPhi_vtx0","std::vector<double>",&patElectron_scPhi_vtx0);
      tree->Branch("patElectron_scEnergy","std::vector<double>",&patElectron_scEnergy); 
      tree->Branch("patElectron_scRawEnergy","std::vector<double>",&patElectron_scRawEnergy); 
      tree->Branch("patElectron_scRawESEnergy","std::vector<double>",&patElectron_scRawESEnergy); 
      tree->Branch("patElectron_scEt","std::vector<double>",&patElectron_scEt); 
      tree->Branch("patElectron_scPhiWidth","std::vector<double>",&patElectron_scPhiWidth);  
      tree->Branch("patElectron_scEtaWidth","std::vector<double>",&patElectron_scEtaWidth);  
      tree->Branch("patElectron_scEoP","std::vector<double>",&patElectron_scEoP);  
      tree->Branch("patElectron_scSwissCross","std::vector<double>",&patElectron_scSwissCross);  
      tree->Branch("patElectron_scR9","std::vector<double>",&patElectron_scR9); 
      tree->Branch("patElectron_scEMax","std::vector<double>",&patElectron_scEMax); 
      tree->Branch("patElectron_scSigmaIEtaIEta","std::vector<double>",&patElectron_scSigmaIEtaIEta);  
      tree->Branch("patElectron_scSigmaIEtaIPhi","std::vector<double>",&patElectron_scSigmaIEtaIPhi);  
      tree->Branch("patElectron_scSigmaIPhiIPhi","std::vector<double>",&patElectron_scSigmaIPhiIPhi);  
      tree->Branch("patElectron_full5x5_scR9","std::vector<double>",&patElectron_full5x5_scR9);  
      tree->Branch("patElectron_full5x5_scEMax","std::vector<double>",&patElectron_full5x5_scEMax);  
      tree->Branch("patElectron_full5x5_scSigmaIEtaIEta","std::vector<double>",&patElectron_full5x5_scSigmaIEtaIEta);  
      tree->Branch("patElectron_full5x5_scSigmaIEtaIPhi","std::vector<double>",&patElectron_full5x5_scSigmaIEtaIPhi);  
      tree->Branch("patElectron_full5x5_scSigmaIPhiIPhi","std::vector<double>",&patElectron_full5x5_scSigmaIPhiIPhi);  
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
   if(savePhotons_){
      tree->Branch("patPhoton_size", &patPhoton_size, "patPhoton_size/I");  
      tree->Branch("patPhoton_index","std::vector<int>",&patPhoton_index); 
      if(isMC_){ 
         tree->Branch("patPhoton_genIndex","std::vector<int>",&patPhoton_genIndex); 
         tree->Branch("patPhoton_genDR","std::vector<double>",&patPhoton_genDR); 
      } 
      tree->Branch("patPhoton_scNPFClusters","std::vector<int>",&patPhoton_scNPFClusters);   
      tree->Branch("patPhoton_passElectronVeto","std::vector<bool> ",&patPhoton_passElectronVeto); 
      tree->Branch("patPhoton_hasPixelSeed","std::vector<bool> ",&patPhoton_hasPixelSeed); 
      tree->Branch("patPhoton_hasConversionTracks","std::vector<bool> ",&patPhoton_hasConversionTracks); 
      tree->Branch("patPhoton_nConversions","std::vector<int>",&patPhoton_nConversions);     
      tree->Branch("patPhoton_nConversionsOneLeg","std::vector<int>",&patPhoton_nConversionsOneLeg);     
      tree->Branch("patPhoton_isEB","std::vector<bool> ",&patPhoton_isEB); 
      tree->Branch("patPhoton_isEE","std::vector<bool> ",&patPhoton_isEE); 
      tree->Branch("patPhoton_isEBEEGap","std::vector<bool> ",&patPhoton_isEBEEGap);
      tree->Branch("patPhoton_isEBEtaGap","std::vector<bool> ",&patPhoton_isEBEtaGap);
      tree->Branch("patPhoton_isEBPhiGap","std::vector<bool> ",&patPhoton_isEBPhiGap);
      tree->Branch("patPhoton_isEEDeeGap","std::vector<bool> ",&patPhoton_isEEDeeGap);
      tree->Branch("patPhoton_isEERingGap","std::vector<bool> ",&patPhoton_isEERingGap);  
      tree->Branch("patPhoton_eta","std::vector<double>",&patPhoton_eta); 
      tree->Branch("patPhoton_phi","std::vector<double>",&patPhoton_phi); 
      tree->Branch("patPhoton_energy","std::vector<double>",&patPhoton_energy);   
      tree->Branch("patPhoton_energyErr","std::vector<double>",&patPhoton_energyErr);   
      tree->Branch("patPhoton_ecalEnergy","std::vector<double>",&patPhoton_ecalEnergy); 
      tree->Branch("patPhoton_ecalEnergyErr","std::vector<double>",&patPhoton_ecalEnergyErr);  
      tree->Branch("patPhoton_et","std::vector<double>",&patPhoton_et); 
      tree->Branch("patPhoton_pt","std::vector<double>",&patPhoton_pt); 
      tree->Branch("patPhoton_HoE","std::vector<double>",&patPhoton_HoE);  
      tree->Branch("patPhoton_scEta","std::vector<double>",&patPhoton_scEta); 
      tree->Branch("patPhoton_scPhi","std::vector<double>",&patPhoton_scPhi);  
      tree->Branch("patPhoton_scEta_vtx0","std::vector<double>",&patPhoton_scEta_vtx0); 
      tree->Branch("patPhoton_scPhi_vtx0","std::vector<double>",&patPhoton_scPhi_vtx0);  
      tree->Branch("patPhoton_scEnergy","std::vector<double>",&patPhoton_scEnergy); 
      tree->Branch("patPhoton_scRawEnergy","std::vector<double>",&patPhoton_scRawEnergy); 
      tree->Branch("patPhoton_scRawESEnergy","std::vector<double>",&patPhoton_scRawESEnergy); 
      tree->Branch("patPhoton_scEt","std::vector<double>",&patPhoton_scEt); 
      tree->Branch("patPhoton_scPhiWidth","std::vector<double>",&patPhoton_scPhiWidth);  
      tree->Branch("patPhoton_scEtaWidth","std::vector<double>",&patPhoton_scEtaWidth);  
      tree->Branch("patPhoton_scSwissCross","std::vector<double>",&patPhoton_scSwissCross);  
      tree->Branch("patPhoton_scR9","std::vector<double>",&patPhoton_scR9); 
      tree->Branch("patPhoton_scEMax","std::vector<double>",&patPhoton_scEMax); 
      tree->Branch("patPhoton_scSigmaIEtaIEta","std::vector<double>",&patPhoton_scSigmaIEtaIEta);  
      tree->Branch("patPhoton_scSigmaIEtaIPhi","std::vector<double>",&patPhoton_scSigmaIEtaIPhi);  
      tree->Branch("patPhoton_scSigmaIPhiIPhi","std::vector<double>",&patPhoton_scSigmaIPhiIPhi);  
      tree->Branch("patPhoton_full5x5_scR9","std::vector<double>",&patPhoton_full5x5_scR9);  
      tree->Branch("patPhoton_full5x5_scEMax","std::vector<double>",&patPhoton_full5x5_scEMax);  
      tree->Branch("patPhoton_full5x5_scSigmaIEtaIEta","std::vector<double>",&patPhoton_full5x5_scSigmaIEtaIEta);  
      tree->Branch("patPhoton_full5x5_scSigmaIEtaIPhi","std::vector<double>",&patPhoton_full5x5_scSigmaIEtaIPhi);  
      tree->Branch("patPhoton_full5x5_scSigmaIPhiIPhi","std::vector<double>",&patPhoton_full5x5_scSigmaIPhiIPhi);  
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
   genParticle_index.resize(nGenParts);
   genParticle_pdgId.resize(nGenParts);
   genParticle_status.resize(nGenParts);
   genParticle_energy.resize(nGenParts);
   genParticle_pt.resize(nGenParts);
   genParticle_eta.resize(nGenParts);
   genParticle_phi.resize(nGenParts);
   genParticle_mass.resize(nGenParts);
   genParticle_charge.resize(nGenParts);
    
   patElectron_index.clear();
   patElectron_genIndex.clear();
   patElectron_genDR.clear();
   patElectron_genFSREnergy.clear();
   patElectron_classification.clear();
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
   patElectron_et.clear();
   patElectron_HoE.clear();
   patElectron_scEta.clear();
   patElectron_scPhi.clear();
   patElectron_scEta_vtx0.clear();
   patElectron_scPhi_vtx0.clear();
   patElectron_scEnergy.clear();
   patElectron_scRawEnergy.clear();
   patElectron_scRawESEnergy.clear();
   patElectron_scEt.clear();
   patElectron_scPhiWidth.clear();
   patElectron_scEtaWidth.clear();
   patElectron_scEoP.clear();
   patElectron_scSwissCross.clear();
   patElectron_scR9.clear();
   patElectron_scEMax.clear();
   patElectron_scSigmaIEtaIEta.clear();
   patElectron_scSigmaIEtaIPhi.clear();
   patElectron_scSigmaIPhiIPhi.clear();
   patElectron_full5x5_scR9.clear();
   patElectron_full5x5_scEMax.clear();
   patElectron_full5x5_scSigmaIEtaIEta.clear();
   patElectron_full5x5_scSigmaIEtaIPhi.clear();
   patElectron_full5x5_scSigmaIPhiIPhi.clear();
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
      patElectron_genIndex.resize(nElectrons);
      patElectron_genDR.resize(nElectrons);
      patElectron_genFSREnergy.resize(nElectrons);
      patElectron_classification.resize(nElectrons);
      patElectron_scNPFClusters.resize(nElectrons);
      patElectron_charge.resize(nElectrons);
      patElectron_isEB.resize(nElectrons);
      patElectron_isEE.resize(nElectrons);  
      patElectron_isEBEEGap.resize(nElectrons);
      patElectron_isEBEtaGap.resize(nElectrons);
      patElectron_isEBPhiGap.resize(nElectrons);
      patElectron_isEEDeeGap.resize(nElectrons);
      patElectron_isEERingGap.resize(nElectrons); 
      patElectron_isEcalDriven.resize(nElectrons);
      patElectron_isTrackerDriven.resize(nElectrons);
      patElectron_passConversionVeto.resize(nElectrons);
      patElectron_eta.resize(nElectrons);
      patElectron_phi.resize(nElectrons);
      patElectron_p.resize(nElectrons);
      patElectron_pt.resize(nElectrons);
      patElectron_pIn.resize(nElectrons);
      patElectron_pOut.resize(nElectrons);
      patElectron_trackFbrem.resize(nElectrons);
      patElectron_superClusterFbrem.resize(nElectrons);
      patElectron_energy.resize(nElectrons);
      patElectron_energyErr.resize(nElectrons);
      patElectron_ecalEnergy.resize(nElectrons);
      patElectron_ecalEnergyErr.resize(nElectrons);
      patElectron_et.resize(nElectrons);
      patElectron_HoE.resize(nElectrons);
      patElectron_scEta.resize(nElectrons);
      patElectron_scPhi.resize(nElectrons);
      patElectron_scEta_vtx0.resize(nElectrons);
      patElectron_scPhi_vtx0.resize(nElectrons);
      patElectron_scEnergy.resize(nElectrons);
      patElectron_scRawEnergy.resize(nElectrons);
      patElectron_scRawESEnergy.resize(nElectrons);
      patElectron_scEt.resize(nElectrons);
      patElectron_scPhiWidth.resize(nElectrons);
      patElectron_scEtaWidth.resize(nElectrons);
      patElectron_scEoP.resize(nElectrons);
      patElectron_scSwissCross.resize(nElectrons);
      patElectron_scR9.resize(nElectrons);
      patElectron_scEMax.resize(nElectrons);
      patElectron_scSigmaIEtaIEta.resize(nElectrons);
      patElectron_scSigmaIEtaIPhi.resize(nElectrons);
      patElectron_scSigmaIPhiIPhi.resize(nElectrons);
      patElectron_full5x5_scR9.resize(nElectrons);
      patElectron_full5x5_scEMax.resize(nElectrons);
      patElectron_full5x5_scSigmaIEtaIEta.resize(nElectrons);
      patElectron_full5x5_scSigmaIEtaIPhi.resize(nElectrons);
      patElectron_full5x5_scSigmaIPhiIPhi.resize(nElectrons);
      patElectron_egmCutBasedElectronIDVeto.resize(nElectrons);
      patElectron_egmCutBasedElectronIDloose.resize(nElectrons);
      patElectron_egmCutBasedElectronIDmedium.resize(nElectrons);
      patElectron_egmCutBasedElectronIDtight.resize(nElectrons);
      patElectron_egmMVAElectronIDScore.resize(nElectrons);
      patElectron_egmMVAElectronIDNoIsoScore.resize(nElectrons);
      patElectron_egmMVAElectronIDloose.resize(nElectrons);
      patElectron_egmMVAElectronIDmedium.resize(nElectrons);
      patElectron_egmMVAElectronIDtight.resize(nElectrons);
      patElectron_egmMVAElectronIDlooseNoIso.resize(nElectrons);
      patElectron_egmMVAElectronIDmediumNoIso.resize(nElectrons);
      patElectron_egmMVAElectronIDtightNoIso.resize(nElectrons);
      patElectron_heepElectronID.resize(nElectrons);
   }

   patPhoton_index.clear();
   patPhoton_genIndex.clear();
   patPhoton_genDR.clear();
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
   patPhoton_scEta.clear();
   patPhoton_scPhi.clear();
   patPhoton_scEta_vtx0.clear();
   patPhoton_scPhi_vtx0.clear();
   patPhoton_scEnergy.clear();
   patPhoton_scRawEnergy.clear();
   patPhoton_scRawESEnergy.clear();
   patPhoton_scEt.clear();
   patPhoton_scPhiWidth.clear();
   patPhoton_scEtaWidth.clear();
   patPhoton_scSwissCross.clear();
   patPhoton_scR9.clear();
   patPhoton_scEMax.clear();
   patPhoton_scSigmaIEtaIEta.clear();
   patPhoton_scSigmaIEtaIPhi.clear();
   patPhoton_scSigmaIPhiIPhi.clear(); 
   patPhoton_full5x5_scR9.clear();
   patPhoton_full5x5_scEMax.clear();
   patPhoton_full5x5_scSigmaIEtaIEta.clear();
   patPhoton_full5x5_scSigmaIEtaIPhi.clear();
   patPhoton_full5x5_scSigmaIPhiIPhi.clear();
   patPhoton_egmCutBasedPhotonIDloose.clear();
   patPhoton_egmCutBasedPhotonIDmedium.clear();
   patPhoton_egmCutBasedPhotonIDtight.clear();
   patPhoton_egmMVAPhotonIDScore.clear();
   patPhoton_egmMVAPhotonIDmedium.clear();
   patPhoton_egmMVAPhotonIDtight.clear();
   patPhoton_passPreselections.clear();
   if(savePhotons_){
      patPhoton_index.resize(nPhotons);
      patPhoton_genIndex.resize(nPhotons);
      patPhoton_genDR.resize(nPhotons);
      patPhoton_scNPFClusters.resize(nPhotons);
      patPhoton_isEB.resize(nPhotons);
      patPhoton_isEE.resize(nPhotons);
      patPhoton_isEBEEGap.resize(nPhotons); 
      patPhoton_isEBEtaGap.resize(nPhotons);
      patPhoton_isEBPhiGap.resize(nPhotons);
      patPhoton_isEEDeeGap.resize(nPhotons);
      patPhoton_isEERingGap.resize(nPhotons); 
      patPhoton_passElectronVeto.resize(nPhotons);
      patPhoton_hasPixelSeed.resize(nPhotons); 
      patPhoton_hasConversionTracks.resize(nPhotons);
      patPhoton_nConversions.resize(nPhotons);
      patPhoton_nConversionsOneLeg.resize(nPhotons);
      patPhoton_eta.resize(nPhotons);
      patPhoton_phi.resize(nPhotons);
      patPhoton_energy.resize(nPhotons);
      patPhoton_energyErr.resize(nPhotons);
      patPhoton_ecalEnergy.resize(nPhotons);
      patPhoton_ecalEnergyErr.resize(nPhotons);
      patPhoton_et.resize(nPhotons);
      patPhoton_pt.resize(nPhotons);
      patPhoton_HoE.resize(nPhotons);
      patPhoton_scEta.resize(nPhotons);
      patPhoton_scPhi.resize(nPhotons);
      patPhoton_scEta_vtx0.resize(nPhotons);
      patPhoton_scPhi_vtx0.resize(nPhotons);
      patPhoton_scEnergy.resize(nPhotons);
      patPhoton_scRawEnergy.resize(nPhotons);
      patPhoton_scRawESEnergy.resize(nPhotons);
      patPhoton_scEt.resize(nPhotons);
      patPhoton_scPhiWidth.resize(nPhotons);
      patPhoton_scEtaWidth.resize(nPhotons);
      patPhoton_scSwissCross.resize(nPhotons);
      patPhoton_scR9.resize(nPhotons);
      patPhoton_scEMax.resize(nPhotons);
      patPhoton_scSigmaIEtaIEta.resize(nPhotons);
      patPhoton_scSigmaIEtaIPhi.resize(nPhotons);
      patPhoton_scSigmaIPhiIPhi.resize(nPhotons); 
      patPhoton_full5x5_scR9.resize(nPhotons);
      patPhoton_full5x5_scEMax.resize(nPhotons);
      patPhoton_full5x5_scSigmaIEtaIEta.resize(nPhotons);
      patPhoton_full5x5_scSigmaIEtaIPhi.resize(nPhotons);
      patPhoton_full5x5_scSigmaIPhiIPhi.resize(nPhotons);
      patPhoton_egmCutBasedPhotonIDloose.resize(nPhotons);
      patPhoton_egmCutBasedPhotonIDmedium.resize(nPhotons);
      patPhoton_egmCutBasedPhotonIDtight.resize(nPhotons);
      patPhoton_egmMVAPhotonIDScore.resize(nPhotons);
      patPhoton_egmMVAPhotonIDmedium.resize(nPhotons);
      patPhoton_egmMVAPhotonIDtight.resize(nPhotons);
      patPhoton_passPreselections.resize(nPhotons);
   }
   
}

const reco::GenParticle* ElePhoDumper::getFirstMother( const reco::GenParticle* particle, int pdgId, int status ) 
{
    if (!particle) return nullptr;
    
    const reco::GenParticle* current = particle;
    const reco::GenParticle* first = nullptr;
    while (current->numberOfMothers() > 0) 
    {
           const reco::GenParticle* mother = dynamic_cast<const reco::GenParticle*>(current->motherRef().get());
           if(!mother){ 
              return nullptr;
           }   
           if(status==-999){
              //if(debug_) std::cout << "getFirstMother - Mother: pdgId = " << mother->pdgId() << " - status = " << mother->status() << std::endl;
              if(mother->pdgId()==pdgId) first = mother;
           }else{
              //if(debug_) std::cout << "getFirstMother - Mother: pdgId = " << mother->pdgId() << " - status = " << mother->status() << std::endl;
              if(mother->status()==status && mother->pdgId()==pdgId) first = mother;
           }
           current = mother;  // Move up the ancestry
    } 

    // return first matched mother
    return first;
}

void ElePhoDumper::printAllDaughters( const reco::GenParticle* particle, int depth ) 
{
   // Base case: safety check
   if(!particle) return;
   if(depth==0){
      std::cout << "Mother PDG ID: " << particle->pdgId()
                << ", status: " << particle->status()
                << ", energy: " << particle->energy() << std::endl;
   }
   
   // Loop over direct daughters
   for(unsigned int i = 0; i < particle->numberOfDaughters(); ++i) 
   {
       const reco::Candidate* dau = particle->daughterRef(i).get();

       // Indentation for visual hierarchy
       std::cout << std::string(depth * 2, ' ')
                 << "Daughter PDG ID: " << dau->pdgId()
                 << ", status: " << dau->status()
                 << ", energy: " << dau->energy() << std::endl;

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

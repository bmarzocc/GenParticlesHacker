#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "TLorentzVector.h"
#include <iostream>

class GenParticleAnalyzer : public edm::EDAnalyzer {
public:
    explicit GenParticleAnalyzer(const edm::ParameterSet&);
    ~GenParticleAnalyzer() override = default;

    void analyze(const edm::Event&, const edm::EventSetup&) override;
    int getGenIndex( edm::Handle<std::vector<reco::GenParticle> > genParticles, const reco::GenParticle* particle ); 
    const reco::GenParticle* getLastDaughterSamePdgID( const reco::GenParticle* particle ); 
    const reco::GenParticle* getFirstMother( const reco::GenParticle* particle, int pdgId, int status ); 
    const reco::GenParticle* getLastMother( const reco::GenParticle* particle, int pdgId, int status); 
    void printAllDaughters( const reco::GenParticle* particle, int depth, bool printAll ); 

private:
    edm::EDGetTokenT<reco::GenParticleCollection> genParticlesToken_;
    edm::EDGetTokenT<LHEEventProduct> lheEventsToken_;
    bool debug_;
    bool isParticleGun_;
    std::vector<int> motherPdgId_;
    std::vector<int> motherStatus_;
    std::vector<int> excludeMomPdgId_;
    std::vector<int> pdgIdBeforeFSR_;
    std::vector<int> pdgIdIn_;
    std::vector<int> statusIn_;
    bool ignoreTauDecays_;
    bool doLHEMatching_;
    double dRMaxLHEMatch_;
    const reco::GenParticle* genPart;
    const reco::GenParticle* motherPart;
    const reco::GenParticle* mother;
    std::vector<const reco::GenParticle*> mothers;
    TLorentzVector vec;
    std::vector<TLorentzVector> lheParticles;
    std::vector<int> lheParticlesID;
    std::vector<int> lheParticlesStatus;
};

GenParticleAnalyzer::GenParticleAnalyzer(const edm::ParameterSet& iConfig): 
    genParticlesToken_(consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genParticles"))),
    lheEventsToken_(consumes<LHEEventProduct>(iConfig.getParameter<edm::InputTag>("lheEvents")))
{
    debug_           = iConfig.getParameter<bool>( "debug" ); 
    isParticleGun_   = iConfig.getParameter<bool>( "isParticleGun" ); 
    motherPdgId_     = iConfig.getParameter<std::vector<int> >( "motherPdgId" );
    motherStatus_    = iConfig.getParameter<std::vector<int> >( "motherStatus" );
    excludeMomPdgId_ = iConfig.getParameter<std::vector<int> >( "excludeMomPdgId" );
    pdgIdBeforeFSR_  = iConfig.getParameter<std::vector<int> >( "pdgIdBeforeFSR" );
    pdgIdIn_         = iConfig.getParameter<std::vector<int> >( "pdgIdIn" );
    statusIn_        = iConfig.getParameter<std::vector<int> >( "statusIn" );
    ignoreTauDecays_ = iConfig.getParameter<bool>( "ignoreTauDecays" );
    doLHEMatching_   = iConfig.getParameter<bool>( "doLHEMatching" );
    dRMaxLHEMatch_   = iConfig.getParameter<double>( "dRMaxLHEMatch" );
    
    assert(motherPdgId_.size()==motherStatus_.size());
    assert(pdgIdIn_.size()==statusIn_.size());
}

void GenParticleAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
    edm::Handle<std::vector<reco::GenParticle>> genParticles;
    iEvent.getByToken(genParticlesToken_, genParticles);

    if (!genParticles.isValid()) {
        edm::LogWarning("GenParticleAnalyzer") << "GenParticles not found!";
        return;
    }
    
    edm::Handle<LHEEventProduct> lheEvents;
    if(doLHEMatching_){
       iEvent.getByToken(lheEventsToken_, lheEvents);
    }
    
    if(doLHEMatching_){
       if(!lheEvents.isValid()) {
          edm::LogWarning("GenParticleAnalyzer") << "No LHEEventProduct found!";
          return;
       }
    }
    
    if(isParticleGun_){
       for(unsigned int i = 0; i <genParticles->size(); ++i)
       {
           genPart = &genParticles->at(i);
           for(unsigned int j = 0; j <pdgIdIn_.size(); ++j)
           {
               if(genPart->pdgId()==pdgIdIn_.at(j) && genPart->status()==statusIn_.at(j))
                  printAllDaughters( genPart, 0, true ); 
           }
       } 
       return;   
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
           vec.SetPxPyPzE(p[0],p[1],p[2],p[3]);
           if(debug_){
              std::cout<< "Particle " << i
                       << ": pdgId = " << hepeup.IDUP[i]
        	       << ", status = " << hepeup.ISTUP[i]
        	       << ", px = " << p[0]
        	       << ", py = " << p[1]
        	       << ", pz = " << p[2]
        	       << ", eta = " << vec.Eta()
        	       << ", phi = " << vec.Phi()
        	       << ", pt = " << vec.Pt()
        	       << ", E = " << p[3]
        	       << ", m = " << p[4] << std::endl;
           }          	       
           lheParticles.push_back(vec);
           lheParticlesID.push_back(hepeup.IDUP[i]);
           lheParticlesStatus.push_back(hepeup.ISTUP[i]);
       }    	 
    }

    //Look for the mothers
    mothers.clear();
    motherPart = nullptr;                 
    for(unsigned int i = 0; i <genParticles->size(); ++i)
    {
        genPart = &genParticles->at(i);
        for(unsigned int j = 0; j <pdgIdIn_.size(); ++j)
        {
            if(genPart->pdgId()==pdgIdIn_.at(j) && genPart->status()==statusIn_.at(j))
            { 
               if(ignoreTauDecays_ && genPart->isDirectPromptTauDecayProductFinalState()==1) continue; 
               
               for(unsigned int k = 0; k <pdgIdBeforeFSR_.size(); ++k)
               {
                   motherPart = getFirstMother(&(*genPart),pdgIdBeforeFSR_.at(k),-999); //get the particle before the FSR
                   if(motherPart==nullptr){ 
                      if(debug_) std::cout << "WARNING: no good particle before FSR!" << std::endl; 
                      continue;
                   }
                   if(motherPart->motherRef().get()->pdgId()==motherPart->pdgId()) continue;
                   
                   if(doLHEMatching_){
                      double dR_tmp=0.8;
                      int lheIndex = -1;
                      for(unsigned int k=0; k<lheParticles.size(); k++){
                          double dR = deltaR(motherPart->eta(),motherPart->phi(),lheParticles.at(k).Eta(),lheParticles.at(k).Phi()); 
                          if(dR<dR_tmp){
                             dR_tmp = dR;
                             lheIndex = k; 
                          }
                      }
                      if(lheIndex<0) continue;    
                      if(debug_){ 
                         std::cout << "Particle before FSR: pdgId = " << motherPart->pdgId() << " - status = " <<  motherPart->status() << " -  eta = " << motherPart->eta() << " - phi = " << motherPart->phi() << " - pt = " << motherPart->pt() << " -  energy = " << motherPart->energy() << std::endl; 
                         std::cout << " --> Matched lhe-particle: dR = " << dR_tmp << " - lheIndex = " << lheIndex << std::endl;  
                      }   
                   }
                   
                   bool goodMotherID=false;
                   bool goodMotherStatus=false;
                   for(unsigned int k=0; k<motherPdgId_.size(); k++){
                       if(motherPdgId_.at(k)!=-999 && motherPart->motherRef().get()->pdgId()==motherPdgId_.at(k)) goodMotherID=true; 
                       else goodMotherID=true;     
                       if(motherStatus_.at(k)!=-999 && motherPart->motherRef().get()->status()==motherStatus_.at(k)) goodMotherStatus=true;   
                       else goodMotherStatus=true;  
                   }    
                   if(!(goodMotherID && goodMotherStatus)) continue;
                   
                   bool badMotherID=false;
                   for(unsigned int k=0; k<excludeMomPdgId_.size(); k++){
                       if(excludeMomPdgId_.at(k)!=-999 && motherPart->motherRef().get()->pdgId()==excludeMomPdgId_.at(k)) badMotherID=true;    
                   }
                   if(badMotherID) continue; //exclude first mothers
                   
                   mothers.push_back(motherPart->motherRef().get()); 
               }       
            }           
        }
    }
 
    if(mothers.size()==0){ 
       std::cout << "WARNING: Bad event --> no good genParticles..." << std::endl;
       return;
    }
    
    //Remove duplicates
    std::sort( mothers.begin(), mothers.end() );
    mothers.erase( unique( mothers.begin(), mothers.end() ), mothers.end() );
    
    std::cout << "GenParticles trees... " << std::endl;
    for(unsigned int i=0; i<mothers.size(); i++){
        printAllDaughters( mothers.at(i), 0, false );    
    }
    
    //Check FSR photons
    if(mother!=nullptr) std::cout << "--> CHECK FSR-PHOTONS " << std::endl;
    for(unsigned int i = 0; i <genParticles->size(); ++i)
    {
        genPart = &genParticles->at(i);
        for(unsigned int j = 0; j <pdgIdBeforeFSR_.size(); ++j)
        {
            if(genPart->pdgId()==22 && genPart->status()==1)
            { 
               if(genPart->mother()->pdgId()!=pdgIdBeforeFSR_.at(j)) continue;
               motherPart = getFirstMother(&(*genPart),pdgIdBeforeFSR_.at(j),-999);
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
               
               std::cout << "FSR photon: p4 = (" << genPart->pt() << "," << genPart->eta() << "," << genPart->phi() << "," << genPart->energy() << ") - last-mother: pdgId = " << motherPart->pdgId() << " , status = " << motherPart->status() << " , p4 = (" << motherPart->pt() << "," << motherPart->eta() << "," << motherPart->phi() << "," << motherPart->energy() << ")" << std::endl;             
            }
        }
    }
    std::cout << "" << std::endl;
}

int GenParticleAnalyzer::getGenIndex( edm::Handle<std::vector<reco::GenParticle> > genParticles, const reco::GenParticle* particle ) 
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

const reco::GenParticle* GenParticleAnalyzer::getLastDaughterSamePdgID( const reco::GenParticle* particle ) 
{
    if (!particle) return nullptr;

    int pdgId = particle->pdgId();
    for (size_t i = 0; i < particle->numberOfDaughters(); ++i) {
        const reco::GenParticle* dau = dynamic_cast<const reco::GenParticle*>(particle->daughter(i));
        if(!dau) continue; // Safety check

        // Check if daughter has same pdgId
        if(dau->pdgId() == pdgId){
           // Recursive call to go deeper
           const reco::GenParticle* lastDaughter = getLastDaughterSamePdgID(dau);
           if(lastDaughter) return lastDaughter;
           else return dau;  // If no further daughters, return this one
        }
    }

    // No daughter with same pdgId found
    return particle;
}

const reco::GenParticle* GenParticleAnalyzer::getLastMother( const reco::GenParticle* particle, int pdgId, int status ) 
{
    if (!particle) return nullptr;
    
    const reco::GenParticle* current = particle;
    while (current->numberOfMothers() > 0) 
    {
           const reco::GenParticle* mother = dynamic_cast<const reco::GenParticle*>(current->motherRef().get());
           if(!mother){ 
              return nullptr;
           }   
           if(status==-999){
              //if(debug_) std::cout << "getLastMother - Mother: pdgId = " << mother->pdgId() << " - status = " << mother->status() << std::endl;
              if(mother->pdgId()==pdgId) return mother;
           }else{
              //if(debug_) std::cout << "getLastMother - Mother: pdgId = " << mother->pdgId() << " - status = " << mother->status() << std::endl;
              if(mother->status()==status && mother->pdgId()==pdgId) return mother;
           }
           current = mother;  // Move up the ancestry
    } 
    
    // If no match found
    return nullptr;
}

const reco::GenParticle* GenParticleAnalyzer::getFirstMother( const reco::GenParticle* particle, int pdgId, int status ) 
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

void GenParticleAnalyzer::printAllDaughters( const reco::GenParticle* particle, int depth, bool printAll=false ) 
{
   // Base case: safety check
   if(!particle) return;
   if(depth==0){
      if(!printAll){ 
         std::cout << "Mother PDG ID: " << particle->pdgId()
                << ", status: " << particle->status()
                << ", energy: " << particle->energy() << std::endl;
      }else{
         std::cout << "Mother PDG ID: " << particle->pdgId()
                << ", status: " << particle->status()
                << ", energy: " << particle->energy()  
                << ", pt: " << particle->pt()      
                << ", eta: " << particle->eta()
                << ", phi: " << particle->phi() << std::endl;  
      }          
   }
   
   // Loop over direct daughters
   for(unsigned int i = 0; i < particle->numberOfDaughters(); ++i) 
   {
       const reco::Candidate* dau = particle->daughterRef(i).get();

       // Indentation for visual hierarchy
       if(!printAll){ 
          std::cout << std::string(depth * 2, ' ')
                 << "Daughter PDG ID: " << dau->pdgId()
                 << ", status: " << dau->status()
                 << ", energy: " << dau->energy() << std::endl;
       }else{
          std::cout << std::string(depth * 2, ' ')
                 << "Daughter PDG ID: " << dau->pdgId()
                 << ", status: " << dau->status()
                 << ", energy: " << dau->energy() 
                 << ", pt: " << dau->pt()
                 << ", eta: " << dau->eta()  
                 << ", phi: " << dau->phi() << std::endl;   
       } 
       
       // Try to cast to reco::GenParticle for recursive call
       const reco::GenParticle* genDau = dynamic_cast<const reco::GenParticle*>(dau);
       if(genDau) {
          printAllDaughters(genDau, depth + 1);  // recursive call
       }
   }
}

DEFINE_FWK_MODULE(GenParticleAnalyzer);


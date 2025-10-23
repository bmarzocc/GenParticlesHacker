#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Utilities/interface/EDMException.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

using namespace edm;
using namespace std;

class GenParticlesLabelChanger : public edm::EDProducer 
{
    public:
        explicit GenParticlesLabelChanger( const ParameterSet &  );
        void produce( Event &evt, const EventSetup & ) override;
    private:
        edm::EDGetTokenT<std::vector<int>> genIntToken_;
        edm::EDGetTokenT<std::vector<reco::GenParticle>> genToken_;
        std::unordered_map<const reco::GenParticle*, size_t> origToNewIndex;
};

GenParticlesLabelChanger::GenParticlesLabelChanger(const edm::ParameterSet& iConfig): 
    genIntToken_( consumes<std::vector<int> >( iConfig.getParameter<InputTag> ( "genInts" ) ) ),
    genToken_( consumes<std::vector<reco::GenParticle> >( iConfig.getParameter<InputTag> ( "genCollection" ) ) )
{
    produces<vector<int>>("");
    produces<vector<reco::GenParticle>>("");
}

void GenParticlesLabelChanger::produce( Event &evt, const EventSetup & ) 
{
    edm::Handle<std::vector<int> > genParticleInts; 
    evt.getByToken(genIntToken_,genParticleInts);
    if(!genParticleInts.isValid()) {
       std::cerr << "Analyze --> genParticle Ints not found" << std::endl; 
       return;
    }
        
    edm::Handle<std::vector<reco::GenParticle> > genParticles; 
    evt.getByToken(genToken_,genParticles);
    if(!genParticles.isValid()) {
       std::cerr << "Analyze --> genParticles not found" << std::endl; 
       return;
    }
   
     unique_ptr<vector<int> > newGenInt( new vector<int> ); 
     unique_ptr<vector<reco::GenParticle> > newGenColl( new vector<reco::GenParticle> ); 
     const edm::RefProd<vector<reco::GenParticle> > refProd = evt.getRefBeforePut<vector<reco::GenParticle> >("");
     
     //fill the new collections
     origToNewIndex.clear();
     const reco::GenParticle* genPart;
     for(unsigned int i = 0; i <genParticles->size(); ++i){
         int genInt = genParticleInts->at(i);
         genPart = &genParticles->at(i);
         
         newGenInt->push_back( genInt );
         newGenColl->push_back( *genPart );
         origToNewIndex[genPart] = newGenInt->size()-1;
    }
    
    //fix mother/daughter relationships
    for(size_t i = 0; i < newGenColl->size(); ++i) {
        const reco::GenParticle* original = &genParticles->at(i);
        reco::GenParticle& clone = newGenColl->at(i);

        // Store new mother/daughter refs temporarily
        std::vector<edm::Ref<reco::GenParticleCollection>> newMotherRefs;
        std::vector<edm::Ref<reco::GenParticleCollection>> newDaughterRefs;

        for(size_t m = 0; m < original->numberOfMothers(); ++m) {
            const reco::GenParticle* oldMomPtr = dynamic_cast<const reco::GenParticle*>(original->motherRef(m).get());
            if(!oldMomPtr) continue;

            auto it = origToNewIndex.find(oldMomPtr);
            if(it == origToNewIndex.end()) continue;

            edm::Ref<std::vector<reco::GenParticle>> momRef(refProd, it->second);
            newMotherRefs.push_back(momRef);
        }
            
        for(size_t d = 0; d <original->numberOfDaughters(); ++d) {
            const reco::GenParticle* oldDauPtr = dynamic_cast<const reco::GenParticle*>(original->daughterRef(d).get());
            if(!oldDauPtr) continue;

            auto it = origToNewIndex.find(oldDauPtr);
            if(it == origToNewIndex.end()) continue;

            edm::Ref<std::vector<reco::GenParticle>> dauRef(refProd, it->second);
            newDaughterRefs.push_back(dauRef);
        }

        // Clear old mothers and add new edm::Refs
        clone.clearMothers();
        for(auto& momRef : newMotherRefs) {
            clone.addMother(momRef);
        }
            
        // Clear old daughters and add new edm::Refs
        clone.clearDaughters();
        for(auto& dauRef : newDaughterRefs) {
            clone.addDaughter(dauRef);
        }
    }
            
    assert(newGenInt->size()==newGenColl->size());
    evt.put(std::move(newGenInt),"");
    evt.put(std::move(newGenColl),"");
}

DEFINE_FWK_MODULE(GenParticlesLabelChanger);


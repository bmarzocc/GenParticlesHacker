#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Utilities/interface/EDMException.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "TLorentzVector.h"
#include "TMath.h"

using namespace edm;
using namespace std;

class NewGenParticlesWithFSRProducer : public EDProducer
{

    public:
        NewGenParticlesWithFSRProducer( const ParameterSet & );
    private:
        void produce( Event &, const EventSetup & ) override;
        int getGenIndex( edm::Handle<std::vector<reco::GenParticle> > genParticles, const reco::GenParticle* particle ); 
        const reco::GenParticle* getLastMother( const reco::GenParticle* particle, int pdgId, int status );
        const reco::GenParticle* getFirstMother( const reco::GenParticle* particle, int pdgId, int status );
        void printAllDaughters( const reco::GenParticle* particle, int depth ); 
        edm::EDGetTokenT<std::vector<int> > genIntToken_; 
        edm::EDGetTokenT<std::vector<reco::GenParticle> > genToken_; 
        edm::EDGetTokenT<LHEEventProduct> lheEventsToken_;

        std::vector<int> motherPdgId_;
        std::vector<int> motherStatus_;
        std::vector<int> excludeMomPdgId_;
        std::vector<int> pdgIdIn_;
        std::vector<int> pdgIdOut_;
        std::vector<int> statusIn_;
        std::vector<int> statusOut_;
        std::vector<int> chargeOut_;
        std::vector<double> massOut_;
        bool doLHEMatching_;
        double dRMaxLHEMatch_;
        bool ignoreTauDecays_;
        bool suppressFSR_;
        bool preserveAngles_;
        bool preserveMotherMass_;
        bool debug_;
        
        int genIndex;
        const reco::GenParticle* genPart;
        reco::GenParticle* genPart_clone;
        const reco::GenParticle* motherPart;
        std::vector<const reco::GenParticle*> mothers;
        std::map<int,reco::GenParticle*> modified_genParticles;
        std::vector<int> fsrIndices;
        std::vector<int> fsrMothers;
        std::vector<int> newDaughterRefs;
        std::map<int,std::vector<double>> fsrEnergyMap;
        std::unordered_map<const reco::GenParticle*, size_t> origToNewIndex;
        std::unordered_map<size_t, size_t> newToOrigIndex;
        std::map<int,std::vector<int> > mothersMap;
        std::map<int,std::vector<int> > daughtersMap;
        TLorentzVector lvec;
        math::XYZTLorentzVector vec;
        std::vector<TLorentzVector> lheParticles;
        std::vector<int> lheParticlesID;
        std::vector<int> lheParticlesStatus;

};

NewGenParticlesWithFSRProducer::NewGenParticlesWithFSRProducer( const ParameterSet &iConfig ) :
        genIntToken_( consumes<std::vector<int> >( iConfig.getParameter<InputTag> ( "genInts" ) ) ),
        genToken_( consumes<std::vector<reco::GenParticle> >( iConfig.getParameter<InputTag> ( "genCollection" ) ) ),
        lheEventsToken_(consumes<LHEEventProduct>(iConfig.getParameter<edm::InputTag>("lheEvents") ) )
{
        motherPdgId_     = iConfig.getParameter<std::vector<int> >( "motherPdgId" );
        motherStatus_    = iConfig.getParameter<std::vector<int> >( "motherStatus" );
        excludeMomPdgId_ = iConfig.getParameter<std::vector<int> >( "excludeMomPdgId" );
        pdgIdIn_         = iConfig.getParameter<std::vector<int> >( "pdgIdIn" );
        pdgIdOut_        = iConfig.getParameter<std::vector<int> >( "pdgIdOut" );
        statusIn_        = iConfig.getParameter<std::vector<int> >( "statusIn" );
        statusOut_       = iConfig.getParameter<std::vector<int> >( "statusOut" );
        chargeOut_       = iConfig.getParameter<std::vector<int> >( "chargeOut" );
        massOut_         = iConfig.getParameter<std::vector<double> >( "massOut" );
        doLHEMatching_   = iConfig.getParameter<bool>( "doLHEMatching" );
        dRMaxLHEMatch_   = iConfig.getParameter<double>( "dRMaxLHEMatch" );
        ignoreTauDecays_ = iConfig.getParameter<bool>( "ignoreTauDecays" );
        suppressFSR_     = iConfig.getParameter<bool>( "suppressFSR" );
        preserveAngles_  = iConfig.getParameter<bool>( "preserveAngles" );
        debug_           = iConfig.getParameter<bool>( "debug" );
        
        produces<vector<int> >("modified");
        produces<vector<reco::GenParticle> >("modified");
        
        size_t ref_size = pdgIdIn_.size();
        assert(pdgIdOut_.size()==ref_size && statusIn_.size()==ref_size && statusOut_.size()==ref_size && chargeOut_.size()==ref_size && massOut_.size()==ref_size);
        assert(motherPdgId_.size()==motherStatus_.size());
}

void NewGenParticlesWithFSRProducer::produce( Event &evt, const EventSetup & )
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
        
        edm::Handle<LHEEventProduct> lheEvents;
        if(doLHEMatching_){
           evt.getByToken(lheEventsToken_, lheEvents);
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
        
        unique_ptr<vector<int> > newGenInt( new vector<int> ); 
        unique_ptr<vector<reco::GenParticle> > newGenColl( new vector<reco::GenParticle> ); 
        const edm::RefProd<vector<reco::GenParticle> > refProd = evt.getRefBeforePut<vector<reco::GenParticle> >("modified");
       
        fsrIndices.clear();
        fsrMothers.clear();
        fsrEnergyMap.clear();
        modified_genParticles.clear();
        if(suppressFSR_)
        {
            //find FSR photons indices
            for(unsigned int i = 0; i <genParticles->size(); ++i){
                genPart = &genParticles->at(i);
                genPart_clone = new reco::GenParticle(genParticles->at(i));
                for(unsigned int j = 0; j <pdgIdIn_.size(); ++j){
                    if(genPart->pdgId()==22 && genPart->status()==1){ 
                    
                       if(genPart->mother()->pdgId()!=pdgIdIn_.at(j)) continue;
                       motherPart = getFirstMother(&(*genPart),pdgIdIn_.at(j),-999);
                       if(motherPart==nullptr){ 
                          if(debug_) std::cout << "WARNING: no good particle before FSR!" << std::endl; 
                          continue;
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
                       
                       fsrIndices.push_back(i);
                       genIndex = getGenIndex(genParticles,motherPart);    
                       if(std::find(fsrMothers.begin(), fsrMothers.end(), genIndex) == fsrMothers.end()) fsrMothers.push_back(genIndex);  
                       fsrEnergyMap[genIndex].push_back(genPart->energy()); 
                       lvec.SetPtEtaPhiE(0.,0.,0.,0.);
                       vec.SetPxPyPzE(lvec.Px(),lvec.Py(),lvec.Pz(),lvec.Energy());
                       genPart_clone->setP4(vec);
                       genPart_clone->setStatus(-999); 
                       modified_genParticles[i] = genPart_clone;             
                    }
                }
            }            
        }
        
        //modify the desired genParticles
        mothers.clear();
        motherPart=nullptr;
        for(unsigned int i = 0; i <genParticles->size(); ++i){
            genPart = &genParticles->at(i);
            genPart_clone = new reco::GenParticle(genParticles->at(i));
            for(unsigned int j = 0; j <pdgIdIn_.size(); ++j){
                if(genPart->pdgId()==pdgIdIn_.at(j) && genPart->status()==statusIn_.at(j))
                { 
                   motherPart = getFirstMother(&(*genPart),genPart->pdgId(),-999); //get the particle before the FSR
                   if(motherPart==nullptr){ 
                       std::cout << "motherPart: nullptr" << std::endl;
                      continue;
                   }   
                   if(motherPart->numberOfMothers()>0 && motherPart->motherRef().isNonnull() && motherPart->motherRef().isAvailable()){
                      if(motherPart->motherRef().get()->pdgId()==motherPart->pdgId()) continue;
                   }   
                   
                   if(doLHEMatching_){
                      double dR_tmp=dRMaxLHEMatch_;
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
                       if(motherPart->numberOfMothers()>0 && motherPart->motherRef().isNonnull() && motherPart->motherRef().isAvailable()){
                          if(motherPdgId_.at(k)!=-999 && motherPart->motherRef().get()->pdgId()==motherPdgId_.at(k)) goodMotherID=true; 
                          else goodMotherID=true;    
                          if(motherStatus_.at(k)!=-999 && motherPart->motherRef().get()->status()==motherStatus_.at(k)) goodMotherStatus=true;  
                          else goodMotherStatus=true;     
                       }else{
                          goodMotherID=true; 
                          goodMotherStatus=true;  
                       }
                   } 
                   if(!(goodMotherID && goodMotherStatus)) continue;
                   
                   bool badMotherID=false;
                   for(unsigned int k=0; k<excludeMomPdgId_.size(); k++){
                       if(motherPart->numberOfMothers()>0 && motherPart->motherRef().isNonnull() && motherPart->motherRef().isAvailable()){
                          if(excludeMomPdgId_.at(k)!=-999 && motherPart->motherRef().get()->pdgId()==excludeMomPdgId_.at(k)) badMotherID=true;  
                       }  
                   }
                   
                   if(badMotherID) continue; //exclude first mothers
                   if(motherPart->numberOfMothers()>0 && motherPart->motherRef().isNonnull() && motherPart->motherRef().isAvailable()){
                      mothers.push_back(motherPart->motherRef().get());
                   }    
                   genPart_clone->setPdgId(pdgIdOut_.at(j));
                   genPart_clone->setStatus(statusOut_.at(j)); 
                   genPart_clone->setCharge(chargeOut_.at(j));
                   genPart_clone->setMass(massOut_.at(j));
                   modified_genParticles[i] = genPart_clone;
                   
                   //modify the 4-momenta accordingly if you suppress FSR photons
                   if(suppressFSR_)
                   {
                      if(preserveAngles_){
                         //preserve the angles: don't change eta-phi and update only the energy and pt, summing up all FSR photons energies 
                         genIndex = getGenIndex(genParticles,motherPart);
                         if(fsrEnergyMap.find(genIndex) == fsrEnergyMap.end()) continue;
                      
                         double fsrEnergyTot = 0.;
                         for(unsigned int iFSR=0; iFSR<fsrEnergyMap[genIndex].size(); iFSR++){
                             fsrEnergyTot += fsrEnergyMap[genIndex].at(iFSR);    
                         }
                         double energy = genPart->energy()+fsrEnergyTot; 
                         double pt = sqrt(energy*energy - massOut_.at(j)*massOut_.at(j))/TMath::CosH(genPart->eta());  
                         lvec.SetPtEtaPhiE(pt,genPart->eta(),genPart->phi(),energy);
                      }else{
                         //preserve the mother mass (such as the Z), and therefore change both energy and angles
                         lvec.SetPtEtaPhiE(motherPart->pt(),motherPart->eta(),motherPart->phi(),motherPart->energy());
                      }   
                      vec.SetPxPyPzE(lvec.Px(),lvec.Py(),lvec.Pz(),lvec.Energy());
                      genPart_clone->setP4(vec);
                      modified_genParticles[i] = genPart_clone;
                   }
                }        
            }
        }
        
        if(debug_){ 
           //Remove duplicates
           std::sort( mothers.begin(), mothers.end() );
           mothers.erase( unique( mothers.begin(), mothers.end() ), mothers.end() );
    
           std::cout << "GenParticles trees... " << std::endl;
           for(unsigned int i=0; i<mothers.size(); i++){
               printAllDaughters( mothers.at(i), 0 );    
           } 
           std::cout << "" << std::endl; 
        }
        
        //fill the modified collections and remove FSR photons, if suppressFSR == true
        origToNewIndex.clear();
        newToOrigIndex.clear();
        mothersMap.clear();
        daughtersMap.clear();
        for(unsigned int i = 0; i <genParticles->size(); ++i){
            int genInt = genParticleInts->at(i);
            genPart = &genParticles->at(i);
            if(suppressFSR_){
               if(modified_genParticles.find(i)!=modified_genParticles.end()){
                  if(modified_genParticles[i]->status()==-999){ 
                     if(debug_) std::cout << "FSR photon: p4 = (" << genPart->pt() << "," << genPart->eta() << "," << genPart->phi() << "," << genPart->energy() << std::endl;    
                     continue;
                  }   
               }   
            }
            newGenInt->push_back( genInt );
            if(modified_genParticles.find(i)==modified_genParticles.end()) newGenColl->push_back( *genPart );
            else newGenColl->push_back( *(modified_genParticles[i]) ); 
            origToNewIndex[genPart] = newGenInt->size()-1;
            newToOrigIndex[newGenInt->size()-1] = i;
        }
       
        //fix mother/daughter relationships
        for(size_t i = 0; i < newGenColl->size(); ++i){
            const reco::GenParticle* original = &genParticles->at(newToOrigIndex[i]);
            reco::GenParticle& clone = newGenColl->at(i);
            
            if(suppressFSR_){
               genIndex = getGenIndex(genParticles,original);  
               if(std::find(fsrIndices.begin(), fsrIndices.end(), genIndex) != fsrIndices.end())  continue; 
            }

            // Store new mother/daughter refs temporarily
            std::vector<edm::Ref<reco::GenParticleCollection>> newMotherRefs;
            std::vector<edm::Ref<reco::GenParticleCollection>> newDaughterRefs;

            for(size_t m = 0; m < original->numberOfMothers(); ++m){
                const reco::GenParticle* oldMomPtr = dynamic_cast<const reco::GenParticle*>(original->motherRef(m).get());
                if(!oldMomPtr) continue;
                
                if(suppressFSR_){
                   genIndex = getGenIndex(genParticles,oldMomPtr);  
                   if(std::find(fsrIndices.begin(), fsrIndices.end(), genIndex) != fsrIndices.end())  continue; 
                }

                auto it = origToNewIndex.find(oldMomPtr);
                if(it == origToNewIndex.end()) continue;

                edm::Ref<std::vector<reco::GenParticle>> momRef(refProd, it->second);
                newMotherRefs.push_back(momRef);
            }
            
            for(size_t d = 0; d <original->numberOfDaughters(); ++d){
                const reco::GenParticle* oldDauPtr = dynamic_cast<const reco::GenParticle*>(original->daughterRef(d).get());
                if(!oldDauPtr) continue;
                
                if(suppressFSR_){
                   genIndex = getGenIndex(genParticles,oldDauPtr);  
                   if(std::find(fsrIndices.begin(), fsrIndices.end(), genIndex) != fsrIndices.end())  continue; 
                }
                
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
        evt.put( std::move( newGenInt ) , "modified" );
        evt.put( std::move( newGenColl ) , "modified" );

}

int NewGenParticlesWithFSRProducer::getGenIndex( edm::Handle<std::vector<reco::GenParticle> > genParticles, const reco::GenParticle* particle ) 
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

const reco::GenParticle* NewGenParticlesWithFSRProducer::getLastMother( const reco::GenParticle* particle, int pdgId, int status ) 
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

const reco::GenParticle* NewGenParticlesWithFSRProducer::getFirstMother( const reco::GenParticle* particle, int pdgId, int status ) 
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

void NewGenParticlesWithFSRProducer::printAllDaughters( const reco::GenParticle* particle, int depth ) 
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

//typedef NewGenParticlesWithFSRProducer NewGenParticlesWithFSRProducer;
DEFINE_FWK_MODULE( NewGenParticlesWithFSRProducer );

// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4

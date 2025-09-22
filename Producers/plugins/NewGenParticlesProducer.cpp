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
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/Math/interface/Point3D.h"
#include "TLorentzVector.h"
#include "TMath.h"

using namespace edm;
using namespace std;

class NewGenParticlesProducer : public EDProducer
{

    public:
        NewGenParticlesProducer( const ParameterSet & );
    private:
        void produce( Event &, const EventSetup & ) override;
        int getGenIndex( edm::Handle<std::vector<reco::GenParticle> > genParticles, const reco::GenParticle* particle ); 
        const reco::GenParticle* getFirstMother( const reco::GenParticle* particle, int pdgId, int status); 
        const reco::GenParticle* getLastMother( const reco::GenParticle* particle, int pdgId, int status ); 
        void printAllDaughters( const reco::GenParticle* particle, int depth ); 
        edm::EDGetTokenT<std::vector<int> > genIntToken_; 
        edm::EDGetTokenT<std::vector<reco::GenParticle> > genToken_; 

        int motherPdgId_;
        int motherStatus_;
        std::vector<int> pdgIdIn_;
        std::vector<int> pdgIdOut_;
        std::vector<int> statusIn_;
        std::vector<int> statusOut_;
        std::vector<int> chargeOut_;
        std::vector<double> massOut_;
        bool ignoreTauDecays_;
        bool suppressFSR_;
        bool preserveAngles_;
        bool preserveMotherMass_;
        
        int genIndex;
        const reco::GenParticle* genPart;
        reco::GenParticle* genPart_clone;
        const reco::GenParticle* motherPart;
        std::map<int,reco::GenParticle*> modified_genParticles;
        std::vector<int> fsrMothers;
        std::vector<int> fsrIndex;
        std::vector<int> newDaughterRefs;
        std::map<int,std::vector<double>> fsrEnergyMap;
        TLorentzVector lvec;
        math::XYZTLorentzVector vec;

};

NewGenParticlesProducer::NewGenParticlesProducer( const ParameterSet &iConfig ) :
        genIntToken_( consumes<std::vector<int> >( iConfig.getParameter<InputTag> ( "genInts" ) ) ),
        genToken_( consumes<std::vector<reco::GenParticle> >( iConfig.getParameter<InputTag> ( "genCollection" ) ) )
{
        motherPdgId_     = iConfig.getParameter<int>( "motherPdgId" );
        motherStatus_    = iConfig.getParameter<int>( "motherStatus" );
        pdgIdIn_         = iConfig.getParameter<std::vector<int> >( "pdgIdIn" );
        pdgIdOut_        = iConfig.getParameter<std::vector<int> >( "pdgIdOut" );
        statusIn_        = iConfig.getParameter<std::vector<int> >( "statusIn" );
        statusOut_       = iConfig.getParameter<std::vector<int> >( "statusOut" );
        chargeOut_       = iConfig.getParameter<std::vector<int> >( "chargeOut" );
        massOut_         = iConfig.getParameter<std::vector<double> >( "massOut" );
        ignoreTauDecays_ = iConfig.getParameter<bool>( "ignoreTauDecays" );
        suppressFSR_     = iConfig.getParameter<bool>( "suppressFSR" );
        preserveAngles_  = iConfig.getParameter<bool>( "preserveAngles" );
        
        produces<vector<int> >("modified");
        produces<vector<reco::GenParticle> >("modified");
        
        size_t ref_size = pdgIdIn_.size();
        assert(pdgIdOut_.size()==ref_size && statusIn_.size()==ref_size && statusOut_.size()==ref_size && chargeOut_.size()==ref_size && massOut_.size()==ref_size);
}

void NewGenParticlesProducer::produce( Event &evt, const EventSetup & )
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
       
        fsrIndex.clear();
        fsrMothers.clear();
        fsrEnergyMap.clear();
        modified_genParticles.clear();
        if(suppressFSR_)
        {
            //find FSR photons indices
            for(unsigned int i = 0; i <genParticles->size(); ++i){
                genPart = &genParticles->at(i);
                for(unsigned int j = 0; j <pdgIdIn_.size(); ++j){
                    if(genPart->pdgId()==22 && genPart->status()==1){ 
                       motherPart = getLastMother(&(*genPart),pdgIdIn_.at(j),-999);
                       if(motherPart==nullptr) continue;
                       if(motherPart->mother()->pdgId()!=motherPdgId_) continue;
                       if(motherStatus_!=-999 && motherPart->mother()->status()!=motherStatus_) continue;
                       
                       genIndex = getGenIndex(genParticles,motherPart);    
                       if(std::find(fsrMothers.begin(), fsrMothers.end(), genIndex) == fsrMothers.end()) fsrMothers.push_back(genIndex);  
                       fsrEnergyMap[genIndex].push_back(genPart->energy());
                       
                       genIndex = getGenIndex(genParticles,genPart);   
                       if(genIndex>=0) fsrIndex.push_back(genIndex);
                    }
                }
            }
            
            //clear the refs of FSR photons from the mother
            for(unsigned int j=0; j<fsrIndex.size(); j++){
                newDaughterRefs.clear();
                genPart = &genParticles->at(fsrIndex.at(j));
                const reco::GenParticle* fsrMother = dynamic_cast<const reco::GenParticle*>(genPart->mother());
                genIndex = getGenIndex(genParticles,fsrMother); 
                genPart_clone = new reco::GenParticle(*fsrMother);
                for(unsigned int i = 0; i < genPart_clone->numberOfDaughters(); ++i) 
                {
                    const reco::GenParticle* genDaughter = dynamic_cast<const reco::GenParticle*>(genPart_clone->daughterRef(i).get());
                    int daughterIndex = getGenIndex( genParticles, genDaughter );
                    if(daughterIndex!=fsrIndex.at(j)) newDaughterRefs.push_back(daughterIndex);
                }
                     
                genPart_clone->clearDaughters();
                for(unsigned int i=0; i<newDaughterRefs.size(); i++){
                    edm::Ref<reco::GenParticleCollection> genRef(genParticles, newDaughterRefs.at(i));
                    if(genRef.isNonnull()) {
                       genPart_clone->addDaughter(genRef);
                    }else{
                       std::cerr << "Null edm::Ref at daughter-index " << newDaughterRefs.at(i) << std::endl;
                    }
                }
                modified_genParticles[genIndex] = genPart_clone;
            }           
              
        }
        
        //modify the wanted genParticles
        for(unsigned int i = 0; i <genParticles->size(); ++i){
            genPart = &genParticles->at(i);
            genPart_clone = new reco::GenParticle(genParticles->at(i));
            for(unsigned int j = 0; j <pdgIdIn_.size(); ++j){
                if(genPart->pdgId()==pdgIdIn_.at(j) && genPart->status()==statusIn_.at(j))
                { 
                   if(ignoreTauDecays_ && genPart->isDirectPromptTauDecayProductFinalState()==1) continue; 
                   
                   motherPart = getFirstMother(&(*genPart),motherPdgId_,motherStatus_); 
                   if(motherPart==nullptr) continue;
                   if(motherPart->pdgId()!=motherPdgId_) continue;
                   if(motherStatus_!=-999 && motherPart->status()!=motherStatus_) continue;
                   
                   genPart_clone->setPdgId(pdgIdOut_.at(j));
                   genPart_clone->setStatus(statusOut_.at(j)); 
                   genPart_clone->setCharge(chargeOut_.at(j));
                   genPart_clone->setMass(massOut_.at(j));
                   modified_genParticles[i] = genPart_clone;
                   
                   //modify the 4-momenta accordingly if you suppress FSR photons
                   if(suppressFSR_)
                   {
                      motherPart = getLastMother(&(*genPart),genPart->pdgId(),-999); //get the particle before the FSR
                      if(motherPart==nullptr) continue;
                      if(motherPart->mother()->pdgId()!=motherPdgId_) continue;
                      
                      if(preserveAngles_){
                         //preserve the angles: don't change eta-phi and update only the energy summing up all FSR photons energies 
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
             
        //fill the modified collections and skip FSR photons, if suppressFSR     
        for(unsigned int i = 0; i <genParticles->size(); ++i){
            int genInt = genParticleInts->at(i);
            genPart = &genParticles->at(i);
            genPart_clone = new reco::GenParticle(genParticles->at(i));
            if(suppressFSR_){
               for(unsigned int j=0; j<fsrIndex.size(); j++){
                   genIndex = fsrIndex.at(j);
                   if(int(i)==genIndex) continue;  
               }    
            }
            
            newGenInt->push_back( genInt );
            if(modified_genParticles.find(i)==modified_genParticles.end()) newGenColl->push_back( *genPart_clone );
            else newGenColl->push_back( *(modified_genParticles[i]) ); 
        }
        
        assert(newGenInt->size()==newGenColl->size());
        evt.put( std::move( newGenInt ) , "modified" );
        evt.put( std::move( newGenColl ) , "modified" );

}

int NewGenParticlesProducer::getGenIndex( edm::Handle<std::vector<reco::GenParticle> > genParticles, const reco::GenParticle* particle ) 
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

const reco::GenParticle* NewGenParticlesProducer::getFirstMother( const reco::GenParticle* particle, int pdgId, int status ) 
{
    if (!particle) return nullptr;
    
    const reco::GenParticle* current = particle;
    while (current->numberOfMothers() > 0) 
    {
           const reco::GenParticle* mother = dynamic_cast<const reco::GenParticle*>(current->mother());
           if(!mother){ 
              return nullptr;
           }   
           if(status==-999){
              //std::cout << "Mother: pdgId = " << mother->pdgId() << " - status = " << mother->status() << std::endl;
              if(mother->pdgId()==pdgId) return mother;
           }else{
              //std::cout << "Mother: pdgId = " << mother->pdgId() << " - status = " << mother->status() << std::endl; 
              if(mother->status()==status && mother->pdgId()==pdgId) return mother;
           }
           
           current = mother;  // Move up the ancestry
    } 
    
    // If no match found
    return nullptr;
}

const reco::GenParticle* NewGenParticlesProducer::getLastMother( const reco::GenParticle* particle, int pdgId, int status ) 
{
    if (!particle) return nullptr;
    
    const reco::GenParticle* current = particle;
    const reco::GenParticle* last = nullptr;
    while (current->numberOfMothers() > 0) 
    {
           const reco::GenParticle* mother = dynamic_cast<const reco::GenParticle*>(current->mother());
           if(!mother){ 
              return nullptr;
              break;
           }   
           if(status==-999){
              //std::cout << "Mother: pdgId = " << mother->pdgId() << " - status = " << mother->status() << std::endl;
              if(mother->pdgId()==pdgId) last = mother;
           }else{
              //std::cout << "Mother: pdgId = " << mother->pdgId() << " - status = " << mother->status() << std::endl; 
              if(mother->status()==status && mother->pdgId()==pdgId) last = mother;
           }
           
           current = mother;  // Move up the ancestry
    } 

    // return last matched mother
    return last;
}

void NewGenParticlesProducer::printAllDaughters( const reco::GenParticle* particle, int depth ) 
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

//typedef NewGenParticlesProducer NewGenParticlesProducer;
DEFINE_FWK_MODULE( NewGenParticlesProducer );

// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4

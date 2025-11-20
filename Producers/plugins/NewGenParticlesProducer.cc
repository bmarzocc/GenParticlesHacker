#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Utilities/interface/EDMException.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
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

class NewGenParticlesProducer : public EDProducer
{

    public:
        NewGenParticlesProducer( const ParameterSet & );
    private:
        void produce( Event &, const EventSetup & ) override;
        edm::EDGetTokenT<std::vector<int> > genIntToken_; 
        edm::EDGetTokenT<std::vector<reco::GenParticle> > genToken_; 

        std::vector<int> pdgIdIn_;
        std::vector<int> pdgIdOut_;
        std::vector<int> statusIn_;
        std::vector<int> statusOut_;
        std::vector<int> chargeOut_;
        std::vector<double> massOut_;
        bool conserveEnergy_;
        bool debug_;
        
        HepMC::FourVector p4;
        TLorentzVector lvec;
        math::XYZTLorentzVector vec;
        const reco::GenParticle* genPart;
        reco::GenParticle* genPart_clone;

};

NewGenParticlesProducer::NewGenParticlesProducer( const ParameterSet &iConfig ) :
        genIntToken_( consumes<std::vector<int> >( iConfig.getParameter<InputTag> ( "genInts" ) ) ),
        genToken_( consumes<std::vector<reco::GenParticle> >( iConfig.getParameter<InputTag> ( "genCollection" ) ) )
{
        pdgIdIn_         = iConfig.getParameter<std::vector<int> >( "pdgIdIn" );
        pdgIdOut_        = iConfig.getParameter<std::vector<int> >( "pdgIdOut" );
        statusIn_        = iConfig.getParameter<std::vector<int> >( "statusIn" );
        statusOut_       = iConfig.getParameter<std::vector<int> >( "statusOut" );
        chargeOut_       = iConfig.getParameter<std::vector<int> >( "chargeOut" );
        massOut_         = iConfig.getParameter<std::vector<double> >( "massOut" );
        conserveEnergy_  = iConfig.getParameter<bool>( "conserveEnergy" );
        debug_           = iConfig.getParameter<bool>( "debug" );
        
        produces<vector<int>>("modified");
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
        
        //modify genParticles
        if(debug_) std::cout << "Modifying genParticles..." << std::endl;
        for(unsigned int i = 0; i <genParticles->size(); ++i)
        {
            int genInt = genParticleInts->at(i);
            genPart = &genParticles->at(i);
            genPart_clone = new reco::GenParticle(genParticles->at(i));
            for(unsigned j=0; j<pdgIdIn_.size(); j++)
            {
                if(genPart_clone->pdgId()==pdgIdIn_.at(j) && genPart_clone->status()==statusIn_.at(j))
                {
                   double energy = genPart_clone->energy();
                   double pt = genPart_clone->pt();
                   double eta = genPart_clone->eta();
                   double phi = genPart_clone->phi();
                   if(debug_){
                      std::cout << "gen-before --> PDG ID: " << genPart_clone->pdgId() 
                                << "  pT: " << pt
                                << "  eta: " << eta
                                << "  phi: " << phi
                                << "  E: " << energy
                                << "  status: " << genPart_clone->status()
                                << "  charge: " << genPart_clone->charge()
                                << std::endl;
                   }
                   
                   if(conserveEnergy_){ 
                      pt = sqrt(energy*energy - massOut_.at(j)*massOut_.at(j))/std::cosh(eta);
                   }else{ 
                      energy = sqrt(massOut_.at(j)*massOut_.at(j) + pt*pt*std::cosh(eta)*std::cosh(eta));
                   }   
                   lvec.SetPtEtaPhiE(pt,eta,phi,energy);
                   vec.SetPxPyPzE(lvec.Px(),lvec.Py(),lvec.Pz(),lvec.Energy());
                   genPart_clone->setP4(vec);
                   genPart_clone->setPdgId(pdgIdOut_.at(j));
                   genPart_clone->setStatus(statusOut_.at(j)); 
                   genPart_clone->setCharge(chargeOut_.at(j));
                   genPart_clone->setMass(massOut_.at(j));
                   
                   if(debug_){
                      std::cout << "gen-after  --> PDG ID: " << genPart_clone->pdgId() 
                                << "  pT: " << pt
                                << "  eta: " << eta
                                << "  phi: " << phi
                                << "  E: " << energy
                                << "  status: " << genPart_clone->status()
                                << "  charge: " << genPart_clone->charge()
                                << std::endl;
                   }
                }   
            }
            newGenInt->push_back( genInt );
            newGenColl->push_back( *genPart_clone );
        }

        assert(newGenInt->size()==newGenColl->size());
        evt.put( std::move( newGenInt ) , "modified" );
        evt.put( std::move( newGenColl ) , "modified" );

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

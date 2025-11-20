#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/Utilities/interface/EDMException.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Common/interface/Handle.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"

using namespace edm;
using namespace std;

class HepMCLabelChanger : public edm::EDProducer 
{
    public:
        explicit HepMCLabelChanger( const ParameterSet &  );
        void produce( Event &evt, const EventSetup & ) override;
    private:
        edm::EDGetTokenT<edm::HepMCProduct> hepMCToken_; 
};

HepMCLabelChanger::HepMCLabelChanger(const edm::ParameterSet& iConfig): 
    hepMCToken_(consumes<edm::HepMCProduct>(iConfig.getParameter<edm::InputTag>("hepMC") ) )
{
    produces<edm::HepMCProduct>("");
}

void HepMCLabelChanger::produce( Event &evt, const EventSetup & ) 
{
    edm::Handle<edm::HepMCProduct> hepMC;
    evt.getByToken(hepMCToken_, hepMC);
    if(!hepMC.isValid()) {
       std::cerr << "Analyze --> HepMCProduct not found" << std::endl; 
       return;
    }
   
    const HepMC::GenEvent* hepEvent = hepMC->GetEvent();
    std::unique_ptr<HepMC::GenEvent> newHepMC(new HepMC::GenEvent(*hepEvent));
     
    auto outputHepMC = std::make_unique<edm::HepMCProduct>();
    outputHepMC->addHepMCData(newHepMC.release()); 
    evt.put( std::move( outputHepMC ) , "" );
}

DEFINE_FWK_MODULE(HepMCLabelChanger);


#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ParameterSetReader/interface/ParameterSetReader.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/FWLite/interface/ChainEvent.h"
#include "FWCore/FWLite/interface/FWLiteEnabler.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "TLorentzVector.h"
#include "TFile.h"
#include "TTree.h"
#include "TROOT.h"
#include "TChain.h"
#include "TGraphErrors.h"
#include "TH1F.h"
#include "TH1D.h"
#include "TF1.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TVector2.h"
#include "TMath.h"
#include "TLegend.h"
#include "TEfficiency.h"
#include "TProfile.h"
#include "TStyle.h"
#include "TTreeReader.h"
#include "TLatex.h"
#include "TSystem.h"
#include <algorithm> 
#include <iostream>
#include <utility>
#include <string>
#include <vector>

using namespace std;

const reco::GenParticle* getFirstMother( const reco::GenParticle* particle, int pdgId, int status ) 
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

int getGenIndex( edm::Handle<std::vector<reco::GenParticle> > genParticles, const reco::GenParticle* particle ) 
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

void printAllDaughters( const reco::GenParticle* particle, int depth ) 
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

int main(int argc, char** argv)
{

   const edm::ParameterSet &process  = edm::readPSetsFrom( argv[1] )->getParameter<edm::ParameterSet>( "process" );
   const edm::ParameterSet &filesOpt = process.getParameter<edm::ParameterSet>( "ioFilesOpt" );
   
   int skipEvents                    = filesOpt.getParameter<int>( "skipEvents" );
   int maxEvents                     = filesOpt.getParameter<int>( "maxEvents" );
   std::string hepMC                 = filesOpt.getParameter<std::string>( "hepMC" );
   std::string genCollection         = filesOpt.getParameter<std::string>( "genCollection" );
   bool debug                        = filesOpt.getParameter<bool>( "debug" ); 
   bool isParticleGun                = filesOpt.getParameter<bool>( "isParticleGun" ); 
   bool checkHepMC                   = filesOpt.getParameter<bool>( "checkHepMC" ); 
   bool doLHEMatching                = filesOpt.getParameter<bool>( "doLHEMatching" ); 
   double dRMaxLHEMatch              = filesOpt.getParameter<double>( "dRMaxLHEMatch" ); 
   double dEFracMaxLHEMatch          = filesOpt.getParameter<double>( "dEFracMaxLHEMatch" ); 
   std::vector<string> phoInputFiles = filesOpt.getParameter<std::vector<string>>( "phoInputFiles" ); 
   std::vector<string> eleInputFiles = filesOpt.getParameter<std::vector<string>>( "eleInputFiles" ); 
   
   
   gROOT->SetBatch(kTRUE);
   //gStyle->SetOptStat(0);
   
   // Enable auto-loading of CMSSW libraries for FWLite
   gSystem->Load("libFWCoreFWLite");
   FWLiteEnabler::enable();

   // Create a ChainEvent (works like a TChain but for FWLite/EDM)
   fwlite::ChainEvent phoEvents(phoInputFiles);
   fwlite::ChainEvent eleEvents(eleInputFiles);
   
   TLorentzVector vec;
   std::vector<TLorentzVector> lheParticles;
   std::vector<int> lheParticlesID;
   std::vector<int> lheParticlesStatus;
   
   std::map<std::tuple<int, int, long int>,int> eventMap;
   std::map<std::tuple<int, int, long int>,std::vector<double>> hepEtaMap;
   std::map<std::tuple<int, int, long int>,std::vector<double>> hepPhiMap;
   std::map<std::tuple<int, int, long int>,std::vector<double>> hepPtMap;
   std::map<std::tuple<int, int, long int>,std::vector<double>> hepEnergyMap;
   std::map<std::tuple<int, int, long int>,std::vector<double>> hepPdgIDMap;
   std::map<std::tuple<int, int, long int>,std::vector<double>> hepStatusMap;
   std::map<std::tuple<int, int, long int>,std::vector<double>> genEtaMap;
   std::map<std::tuple<int, int, long int>,std::vector<double>> genPhiMap;
   std::map<std::tuple<int, int, long int>,std::vector<double>> genPtMap;
   std::map<std::tuple<int, int, long int>,std::vector<double>> genEnergyMap;
   std::map<std::tuple<int, int, long int>,std::vector<double>> genPdgIDMap;
   std::map<std::tuple<int, int, long int>,std::vector<double>> genStatusMap;
   std::map<std::tuple<int, int, long int>,std::vector<double>> genChargeMap;
   
   std::vector<double> hepEtaMatch;
   std::vector<double> hepPhiMatch;
   std::vector<double> hepPtMatch;
   std::vector<double> hepEnergyMatch;
   std::vector<double> hepPdgIDMatch;
   std::vector<double> hepStatusMatch;
   std::vector<double> genEtaMatch;
   std::vector<double> genPhiMatch;
   std::vector<double> genPtMatch;
   std::vector<double> genEnergyMatch;
   std::vector<double> genPdgIDMatch;
   std::vector<double> genStatusMatch;
   std::vector<double> genChargeMatch;
  
   // Run on the photon samples
   std::cout << "" << std::endl;
   long int num_events = 0;
   for(phoEvents.toBegin(); !phoEvents.atEnd(); ++phoEvents, ++num_events) 
   {
       std::cout << "" << std::endl;
       if((num_events % 1) == 0) std::cout << "Processing event " << num_events << std::endl;
       //if(skipEvents>=0 && num_events<skipEvents) continue;
       if(maxEvents>=0 && num_events>maxEvents) break; 
       
       edm::Handle<edm::HepMCProduct> hepMCHandle;
       if(checkHepMC){
          phoEvents.getByLabel(hepMC, hepMCHandle);
          if(!hepMCHandle.isValid()) {
             std::cerr << "Analyze --> HepMCProduct not found" << std::endl; 
             continue;
          }
       }
       
       edm::Handle<LHEEventProduct> lheHandle;
       if(doLHEMatching){
          phoEvents.getByLabel(std::string("externalLHEProducer"), lheHandle);
          if(!lheHandle.isValid()){
             std::cout << "Invalid lheProduct-handle!" << std::endl;
             continue;
          }
       }
        
       edm::Handle<reco::GenParticleCollection> genHandle;
       phoEvents.getByLabel(genCollection, genHandle);
       if(!genHandle.isValid()){
          std::cout << "Invalid genParticles-handle!" << std::endl;
          continue;
       }
       
       long int eventId = phoEvents.id().run();
       int lumiId = phoEvents.id().luminosityBlock();
       int runId = phoEvents.id().event();
       
       lheParticles.clear();
       lheParticlesID.clear();
       lheParticlesStatus.clear();
       const LHEEventProduct& lheEvents = *lheHandle;
       if(doLHEMatching)
       {
          const auto& hepeup = lheEvents.hepeup();
          const auto& lheparticles = hepeup.PUP;
          for(size_t j = 0; j < lheparticles.size(); ++j)
          {
              const auto& l = lheparticles[j];
              vec.SetPxPyPzE(l[0],l[1],l[2],l[3]);
              if(debug){
                 std::cout << "Photon-sample: particle " << j
                           << ": pdgId = " << hepeup.IDUP[j]
        	           << ", status = " << hepeup.ISTUP[j]
        	           << ", px = " << l[0]
        	           << ", py = " << l[1]
        	           << ", pz = " << l[2]
        	           << ", eta = " << vec.Eta()
        	           << ", phi = " << vec.Phi()
        	           << ", pt = " << vec.Pt()
        	           << ", E = " << l[3]
        	           << ", m = " << l[4] << std::endl;
              }
              lheParticles.push_back(vec);
              lheParticlesID.push_back(hepeup.IDUP[j]);
              lheParticlesStatus.push_back(hepeup.ISTUP[j]);
          }
       }
       
       hepEtaMatch.clear();
       hepPhiMatch.clear();
       hepPtMatch.clear();
       hepEnergyMatch.clear();  
       hepPdgIDMatch.clear();
       hepStatusMatch.clear();
       
       if(checkHepMC){
          const HepMC::GenEvent* hepEvent = hepMCHandle->GetEvent();
          for(auto p = hepEvent->particles_begin(); p != hepEvent->particles_end(); ++p){
              const HepMC::FourVector& mom = (*p)->momentum();
              double px = mom.px();
              double py = mom.py();
              double pz = mom.pz();
              double energy  = mom.e();

              double pt  = std::sqrt(px*px + py*py);
              double eta = 0.;
              if(pt > 0.0){
                 double theta = std::atan2(pt, pz);
                 eta = -std::log(std::tan(theta / 2.0));
              }
              double phi = std::atan2(py, px);
              if(isParticleGun){
                 hepEtaMatch.push_back(eta);
                 hepPhiMatch.push_back(phi);
                 hepPtMatch.push_back(pt);
                 hepEnergyMatch.push_back(energy);
                 hepPdgIDMatch.push_back((*p)->pdg_id());
                 hepStatusMatch.push_back((*p)->status());
                 continue;
              }
          }      
       }
          	  
       genEtaMatch.clear();
       genPhiMatch.clear();
       genPtMatch.clear();
       genEnergyMatch.clear();  
       genPdgIDMatch.clear();
       genStatusMatch.clear();
       genChargeMatch.clear();        	  
       const reco::GenParticleCollection& genParticles = *genHandle;
       for(size_t i = 0; i < genParticles.size(); ++i) 
       {
           const reco::GenParticle& genPart = genParticles[i];
           if(abs(genPart.pdgId())!=22) continue;
           if(genPart.status()!=1) continue;
           
           if(isParticleGun){
              genEtaMatch.push_back(genPart.eta());
              genPhiMatch.push_back(genPart.phi());
              genPtMatch.push_back(genPart.pt());
              genEnergyMatch.push_back(genPart.energy());
              genPdgIDMatch.push_back(genPart.pdgId());
              genStatusMatch.push_back(genPart.status());
              genChargeMatch.push_back(genPart.charge());
              continue;
           }
           
           const reco::GenParticle* mother = nullptr;
           if(genPart.numberOfMothers()==0) mother = &(genHandle->at(i));
           else{
              mother = getFirstMother( &(genHandle->at(i)), 11, -999 ); 
              if(abs(mother->pdgId())!=11) mother = genPart.motherRef().get();
           }   
           if(mother==nullptr){ 
              if(debug) std::cout << "WARNING: no good mother!" << std::endl; 
              continue;
           }
           if(abs(mother->pdgId())>99) continue;
           if(abs(mother->pdgId())==23) mother = &(genHandle->at(i));
           
           //if(abs(mother->pdgId())==11) std::cout << "number of daughters = " << mother->numberOfMothers() << std::endl;
           
           if(doLHEMatching)
           {
              int lheIndex = -1;
              double dR_tmp=dRMaxLHEMatch; 
              for(size_t j = 0; j < lheParticles.size(); ++j)
              {
        	  if(abs(lheParticlesID.at(j))!=11) continue;
        	  if(lheParticlesStatus.at(j)!=1) continue;
                  double dR = deltaR(mother->eta(),mother->phi(),lheParticles.at(j).Eta(),lheParticles.at(j).Phi()); 
                  if(dR<dR_tmp){
                     dR_tmp = dR;
                     lheIndex = j; 
                  }
              }   
             
              if(lheIndex<0) continue;
              double dE_fraction = abs(lheParticles.at(lheIndex).Energy()-genPart.energy())/lheParticles.at(lheIndex).Energy(); 
              if(dE_fraction>dEFracMaxLHEMatch) continue;
              if(debug) std::cout << "Mathced-LHE: " << lheIndex << " - pdgId = " << lheParticlesID.at(lheIndex) << " - dR = " << dR_tmp << " - lheEnergy = " << lheParticles.at(lheIndex).Energy() << " - dE_fraction = " << dE_fraction << std::endl; 
           }
           
           genEtaMatch.push_back(genPart.eta());
           genPhiMatch.push_back(genPart.phi());
           genPtMatch.push_back(genPart.pt());
           genEnergyMatch.push_back(genPart.energy());
           genPdgIDMatch.push_back(genPart.pdgId());
           genStatusMatch.push_back(genPart.status());
           genChargeMatch.push_back(genPart.charge());
       }   
       
       eventMap.insert(std::pair<typename std::tuple<int, int, long int>, int>(std::make_tuple(runId, lumiId, eventId), 1)); 
       hepEtaMap.insert(std::pair<typename std::tuple<int, int, long int>, std::vector<double>>(std::make_tuple(runId, lumiId, eventId), hepEtaMatch)); 
       hepPhiMap.insert(std::pair<typename std::tuple<int, int, long int>, std::vector<double>>(std::make_tuple(runId, lumiId, eventId), hepPhiMatch)); 
       hepPtMap.insert(std::pair<typename std::tuple<int, int, long int>, std::vector<double>>(std::make_tuple(runId, lumiId, eventId), hepPtMatch)); 
       hepEnergyMap.insert(std::pair<typename std::tuple<int, int, long int>, std::vector<double>>(std::make_tuple(runId, lumiId, eventId), hepEnergyMatch)); 
       hepPdgIDMap.insert(std::pair<typename std::tuple<int, int, long int>, std::vector<double>>(std::make_tuple(runId, lumiId, eventId), hepPdgIDMatch)); 
       hepStatusMap.insert(std::pair<typename std::tuple<int, int, long int>, std::vector<double>>(std::make_tuple(runId, lumiId, eventId), hepStatusMatch));
       genEtaMap.insert(std::pair<typename std::tuple<int, int, long int>, std::vector<double>>(std::make_tuple(runId, lumiId, eventId), genEtaMatch)); 
       genPhiMap.insert(std::pair<typename std::tuple<int, int, long int>, std::vector<double>>(std::make_tuple(runId, lumiId, eventId), genPhiMatch)); 
       genPtMap.insert(std::pair<typename std::tuple<int, int, long int>, std::vector<double>>(std::make_tuple(runId, lumiId, eventId), genPtMatch)); 
       genEnergyMap.insert(std::pair<typename std::tuple<int, int, long int>, std::vector<double>>(std::make_tuple(runId, lumiId, eventId), genEnergyMatch)); 
       genPdgIDMap.insert(std::pair<typename std::tuple<int, int, long int>, std::vector<double>>(std::make_tuple(runId, lumiId, eventId), genPdgIDMatch)); 
       genStatusMap.insert(std::pair<typename std::tuple<int, int, long int>, std::vector<double>>(std::make_tuple(runId, lumiId, eventId), genStatusMatch)); 
       genChargeMap.insert(std::pair<typename std::tuple<int, int, long int>, std::vector<double>>(std::make_tuple(runId, lumiId, eventId), genChargeMatch)); 
   }
   
   // Run on the electron samples
   std::cout << "" << std::endl;
   num_events = 0;
   int nBugged = 0;
   std::map<int,double> fsrEnergyMap;
   std::map<int,double> fsrNumberMap;
   for(eleEvents.toBegin(); !eleEvents.atEnd(); ++eleEvents, ++num_events) 
   {
       std::cout << "" << std::endl;
       if((num_events % 1) == 0) std::cout << "Processing event " << num_events << std::endl;
       //if(skipEvents>=0 && num_events<skipEvents) continue;
       if(maxEvents>=0 && num_events>maxEvents) break; 
       
       edm::Handle<edm::HepMCProduct> hepMCHandle;
       if(checkHepMC){
          eleEvents.getByLabel(hepMC, hepMCHandle);
          if(!hepMCHandle.isValid()) {
             std::cerr << "Analyze --> HepMCProduct not found" << std::endl; 
             continue;
          }
       }
       
       edm::Handle<LHEEventProduct> lheHandle;
       if(doLHEMatching){
          eleEvents.getByLabel(std::string("externalLHEProducer"), lheHandle);
          if(!lheHandle.isValid()){
             std::cout << "Invalid lheProduct-handle!" << std::endl;
             continue;
          }
       }
       
       edm::Handle<reco::GenParticleCollection> genHandle;
       eleEvents.getByLabel(genCollection, genHandle);
       if(!genHandle.isValid()){
          std::cout << "Invalid genParticles-handle!" << std::endl;
          continue;
       }
       
       long int eventId = eleEvents.id().run();
       int lumiId = eleEvents.id().luminosityBlock();
       int runId = eleEvents.id().event();
       if((eventMap.find(std::make_tuple(runId, lumiId, eventId)))->second!=1) continue;
       
       lheParticles.clear();
       lheParticlesID.clear();
       lheParticlesStatus.clear();
       const LHEEventProduct& lheEvents = *lheHandle;
       if(doLHEMatching)
       {
          const auto& hepeup = lheEvents.hepeup();
          const auto& lheparticles = hepeup.PUP;
          for(size_t j = 0; j < lheparticles.size(); ++j)
          {
              const auto& l = lheparticles[j];
              vec.SetPxPyPzE(l[0],l[1],l[2],l[3]);
              if(debug){
                 std::cout << "Electron-sample: particle " << j
                           << ": pdgId = " << hepeup.IDUP[j]
        	           << ", status = " << hepeup.ISTUP[j]
        	           << ", px = " << l[0]
        	           << ", py = " << l[1]
        	           << ", pz = " << l[2]
        	           << ", eta = " << vec.Eta()
        	           << ", phi = " << vec.Phi()
        	           << ", pt = " << vec.Pt()
        	           << ", E = " << l[3]
        	           << ", m = " << l[4] << std::endl;
              }
              lheParticles.push_back(vec);
              lheParticlesID.push_back(hepeup.IDUP[j]);
              lheParticlesStatus.push_back(hepeup.ISTUP[j]);
          }
       }
          
       hepEtaMatch = (hepEtaMap.find(std::make_tuple(runId, lumiId, eventId)))->second;
       hepPhiMatch = (hepPhiMap.find(std::make_tuple(runId, lumiId, eventId)))->second;
       hepPtMatch = (hepPtMap.find(std::make_tuple(runId, lumiId, eventId)))->second;
       hepEnergyMatch = (hepEnergyMap.find(std::make_tuple(runId, lumiId, eventId)))->second; 
       hepPdgIDMatch = (hepPdgIDMap.find(std::make_tuple(runId, lumiId, eventId)))->second; 
       hepStatusMatch = (hepStatusMap.find(std::make_tuple(runId, lumiId, eventId)))->second;
       genEtaMatch = (genEtaMap.find(std::make_tuple(runId, lumiId, eventId)))->second;
       genPhiMatch = (genPhiMap.find(std::make_tuple(runId, lumiId, eventId)))->second;
       genPtMatch = (genPtMap.find(std::make_tuple(runId, lumiId, eventId)))->second;
       genEnergyMatch = (genEnergyMap.find(std::make_tuple(runId, lumiId, eventId)))->second; 
       genPdgIDMatch = (genPdgIDMap.find(std::make_tuple(runId, lumiId, eventId)))->second; 
       genStatusMatch = (genStatusMap.find(std::make_tuple(runId, lumiId, eventId)))->second; 
       genChargeMatch = (genChargeMap.find(std::make_tuple(runId, lumiId, eventId)))->second; 
       
       const reco::GenParticleCollection& genParticles = *genHandle;
       
       // get FSR photons
       fsrNumberMap.clear();
       fsrEnergyMap.clear();
       for(size_t i = 0; i < genParticles.size(); ++i) 
       {
           if(isParticleGun) break;
           
           const reco::GenParticle& genPart = genParticles[i];
           if(abs(genPart.pdgId())!=22) continue;
           if(genPart.status()!=1) continue;
           
           const reco::GenParticle* mother = getFirstMother( &(genHandle->at(i)), 11, -999 ); 
           //if(abs(mother->pdgId())==11) std::cout << "mother: pdgId = " << mother->pdgId() << " - status = " << mother->status() << std::endl; 
           if(abs(mother->pdgId())!=11 || mother->status()==1) continue;
           
           if(mother==nullptr){ 
              if(debug) std::cout << "WARNING: no good particle before FSR!" << std::endl; 
              continue;
           }
           
           int genIndex = getGenIndex(genHandle,mother); 
           if(debug) std::cout << "mother: pdgId = " << mother->pdgId() << " -  status = " << mother->status() << " - index = " << genIndex << " - fsrEnrgy = " << genPart.energy() << std::endl; 
           if(fsrNumberMap.find(genIndex) == fsrNumberMap.end()) fsrNumberMap[genIndex] = 1;
           else fsrNumberMap[genIndex] += 1;
           if(fsrEnergyMap.find(genIndex) == fsrEnergyMap.end()) fsrEnergyMap[genIndex] = genPart.energy();
           else fsrEnergyMap[genIndex] += genPart.energy(); 
       }
       if(debug){
          for(auto &myPair : fsrNumberMap){
              std::cout << "FsrNum: " << myPair.first << " - " << myPair.second << std::endl;
          }
          for(auto &myPair : fsrEnergyMap){
              std::cout << "FsrMap: " << myPair.first << " - " << myPair.second << std::endl;
          }
       }
       
       std::cout << "Running on hepMC" << std::endl; 
       if(checkHepMC){
          const HepMC::GenEvent* hepEvent = hepMCHandle->GetEvent();
          for(auto p = hepEvent->particles_begin(); p != hepEvent->particles_end(); ++p){
       
              if(abs((*p)->pdg_id())!=11) continue;
              if((*p)->status()!=1) continue;
            
              const HepMC::FourVector& mom = (*p)->momentum();
              double px = mom.px();
              double py = mom.py();
              double pz = mom.pz();
              double energy  = mom.e();

              double pt  = std::sqrt(px*px + py*py);
              double eta = 0.;
              if(pt > 0.0){
                 double theta = std::atan2(pt, pz);
                 eta = -std::log(std::tan(theta / 2.0));
              }
              double phi = std::atan2(py, px);
           
              if(isParticleGun){
                 int genIndex = -1;
                 double dR_tmp = 0.1;
                 for(unsigned int iGen=0; iGen<hepEtaMatch.size(); iGen++){
                     double dR = deltaR(hepEtaMatch.at(iGen),hepPhiMatch.at(iGen),eta,phi); 
                     if(dR<dR_tmp){
                        dR_tmp = dR;
                        genIndex = iGen;               
                     }
                 }    
                 //std::cout << "selected: dR = " << dR_tmp << " - genIndex = " << genIndex << std::endl;
                 if(genIndex<0) continue;
              
                 std::cout << "pho-hepEnergy = " << hepEnergyMatch.at(genIndex) << " - ele-hepEnergy = " << energy << std::endl;  
                 std::cout << "pho-hepPt     = " << hepPtMatch.at(genIndex) << " - ele-hepPt     = " << pt << std::endl;     
                 std::cout << "pho-hepEta    = " << hepEtaMatch.at(genIndex) << " - ele-hepEta    = " << eta << std::endl;    
                 std::cout << "pho-hepPhi    = " << hepPhiMatch.at(genIndex) << " - ele-hepPhi    = " << phi << std::endl;  
                 std::cout << "pho-pdgId     = " << hepPdgIDMatch.at(genIndex) << " - ele-pdgId    = " << (*p)->pdg_id() << std::endl;    
                 std::cout << "pho-status    = " << hepStatusMatch.at(genIndex) << " - ele-status    = " << (*p)->status() << std::endl;  
                 continue;
              }
          } 
       }
       
       std::cout << "" << std::endl;  
       std::cout << "genParticles.size() = " << genParticles.size() << std::endl;       	       	 
       for(size_t i = 0; i < genParticles.size(); ++i) 
       {
           const reco::GenParticle& genPart = genParticles[i];
           if(abs(genPart.pdgId())!=11) continue;
           if(genPart.status()!=1) continue;
           
           if(isParticleGun){
              int genIndex = -1;
              double dR_tmp = 0.1;
              for(unsigned int iGen=0; iGen<genEtaMatch.size(); iGen++){
                  double dR = deltaR(genEtaMatch.at(iGen),genPhiMatch.at(iGen),genPart.eta(),genPart.phi()); 
                  if(dR<dR_tmp){
                     dR_tmp = dR;
                     genIndex = iGen;               
                  }
              }    
              //std::cout << "selected: dR = " << dR_tmp << " - genIndex = " << genIndex << std::endl;
              if(genIndex<0) continue;
              
              std::cout << "pho-GenEnergy = " << genEnergyMatch.at(genIndex) << " - ele-GenEnergy = " << genPart.energy() << std::endl;  
              std::cout << "pho-GenPt     = " << genPtMatch.at(genIndex) << " - ele-GenPt     = " << genPart.pt() << std::endl;     
              std::cout << "pho-GenEta    = " << genEtaMatch.at(genIndex) << " - ele-GenEta    = " << genPart.eta() << std::endl;    
              std::cout << "pho-GenPhi    = " << genPhiMatch.at(genIndex) << " - ele-GenPhi    = " << genPart.phi() << std::endl;  
              std::cout << "pho-pdgId     = " << genPdgIDMatch.at(genIndex) << " - ele-pdgId    = " << genPart.pdgId() << std::endl;    
              std::cout << "pho-charge    = " << genChargeMatch.at(genIndex) << " - ele-charge    = " << genPart.charge() << std::endl;  
              std::cout << "pho-status    = " << genStatusMatch.at(genIndex) << " - ele-status    = " << genPart.status() << std::endl;  
              continue;
           }
          
           
           const reco::GenParticle* mother = getFirstMother( &(genHandle->at(i)), 11, -999 ); 
           if(mother==nullptr){ 
              if(debug) std::cout << "WARNING: no good mother!" << std::endl; 
              continue;
           }
           if(abs(mother->pdgId())!=11) continue;
           int genMotherIndex = getGenIndex(genHandle,mother); 
           
           if(doLHEMatching)
           {
              int lheIndex = -1;
              double dR_tmp=dRMaxLHEMatch; 
              for(size_t j = 0; j < lheParticles.size(); ++j)
              {
        	  if(abs(lheParticlesID.at(j))!=11) continue;
        	  if(lheParticlesStatus.at(j)!=1) continue;
                  double dR = deltaR(mother->eta(),mother->phi(),lheParticles.at(j).Eta(),lheParticles.at(j).Phi()); 
                  if(dR<dR_tmp){
                     dR_tmp = dR;
                     lheIndex = j; 
                  }
              }   
             
              if(lheIndex<0) continue;
              double dE_fraction = abs(lheParticles.at(lheIndex).Energy()-genPart.energy())/lheParticles.at(lheIndex).Energy(); 
              if(dE_fraction>dEFracMaxLHEMatch) continue;
              if(debug) std::cout << "Mathced-LHE: " << lheIndex << " - pdgId = " << lheParticlesID.at(lheIndex) << " - dR = " << dR_tmp << " - lheEnergy = " << lheParticles.at(lheIndex).Energy() << " - dE_fraction = " << dE_fraction << std::endl; 
           }
              
           int genIndex = -1;
           double dR_tmp = 0.1;
           for(unsigned int iGen=0; iGen<genEtaMatch.size(); iGen++){
               double dR = deltaR(genEtaMatch.at(iGen),genPhiMatch.at(iGen),genPart.eta(),genPart.phi()); 
               if(dR<dR_tmp){
                  dR_tmp = dR;
                  genIndex = iGen;               
               }
           }
           if(genIndex<0) continue;
           //printAllDaughters( mother, 0 ); 
           
           if(abs(genEnergyMatch.at(genIndex)-genPart.energy())<1.e-6 && fsrEnergyMap[genMotherIndex]>0.){ 
              nBugged++;
              continue;
           }    
           if(abs(genEnergyMatch.at(genIndex)-(genPart.energy()+fsrEnergyMap[genMotherIndex]))>0.001) std::cout << "WARNING no match between electron and photon: " << i << " - matchedGenPhoton: " << genIndex << " - genEnergyPho = " << genEnergyMatch.at(genIndex) << " -  genEnergyEle = " << genPart.energy() << " -  fsrEnergy = " << fsrEnergyMap[genMotherIndex] << " - dE = " << abs(genEnergyMatch.at(genIndex)-(genPart.energy()+fsrEnergyMap[genMotherIndex])) << std::endl; 
       }   
   }
   std::cout << "nBugged events = " << nBugged << std::endl;
} 



# GenParticlesHacker

1) Install:

    * cmssw-el7
    * scram project CMSSW_10_6_19_patch3
    * cd CMSSW_10_6_19_patch3/src/
    * cmsenv
    * git cms-init
    * git clone git@github.com:bmarzocc/GenParticlesHacker.git
    * scram b -j 5

2) Produce gen-step input:

    * cd CMSSW_10_6_19_patch3/src/GenParticlesHacker/Producers/test/
    * cmsRun DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIISummer20UL18MiniAODv2_GENSTEP.py
    
3) Produce modified genParticles collection:

    * #Set all the options in CMSSW_10_6_19_patch3/src/GenParticlesHacker/python/NewGenParticlesProducer_cfi.py
    * cd CMSSW_10_6_19_patch3/src/GenParticlesHacker/
    * cmsRun python/NewGenParticlesProducer_cfg.py
    
4) Produce sim-step output using the modified genParticles:

    * cd CMSSW_10_6_19_patch3/src/GenParticlesHacker/Producers/test/
    * cmsRun DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIISummer20UL18MiniAODv2_SIMSTEP.py

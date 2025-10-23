#!/bin/bash    

SEED=$(( $1 + 1 ))
                                                                                                                                                                                   
set -e
BASE=$PWD
RELEASE_BASE=$CMSSW_BASE

export SCRAM_ARCH=slc7_amd64_gcc700
source /cvmfs/cms.cern.ch/cmsset_default.sh

eval `scram project CMSSW_10_6_48`
cd CMSSW_10_6_48/src/
eval `scram runtime -sh`
cp ../../DYJetsToEE_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIISummer20UL18MiniAODv2_stdGEN_cfg.py ./

echo "cmsRun -e -j FrameworkJobReport.xml DYJetsToEE_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIISummer20UL18MiniAODv2_stdGEN_cfg.py runNumber="$1" maxEvents=10000"
cmsRun -e -j FrameworkJobReport.xml DYJetsToEE_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8_RunIISummer20UL18MiniAODv2_stdGEN_cfg.py runNumber=$SEED

cp step1_stdGEN*.root $BASE/
cp FrameworkJobReport.xml $BASE/
cd $RELEASE_BASE
eval `scram runtime -sh`
cd $BASE
echo $BASE

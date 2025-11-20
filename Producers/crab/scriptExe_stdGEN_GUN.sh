#!/bin/bash    

SEED=$(( $1 + 1 ))
                                                                                                                                                                                   
set -e
BASE=$PWD
RELEASE_BASE=$CMSSW_BASE

export SCRAM_ARCH=slc7_amd64_gcc700
source /cvmfs/cms.cern.ch/cmsset_default.sh

echo "cmsRun -e -j FrameworkJobReport.xml ManyElectron_Pt-15To250-gun-RunIISummer20UL18MiniAODv2_stdGEN_cfg.py runNumber="$SEED
cmsRun -e -j FrameworkJobReport.xml ManyElectron_Pt-15To250-gun-RunIISummer20UL18MiniAODv2_stdGEN_cfg.py runNumber=$SEED


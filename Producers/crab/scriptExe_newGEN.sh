#!/bin/bash    

SEED=$(( $1 + 1 ))
                                                                                                                                                                                   
set -e
BASE=$PWD
RELEASE_BASE=$CMSSW_BASE

export SCRAM_ARCH=slc7_amd64_gcc700
source /cvmfs/cms.cern.ch/cmsset_default.sh

echo "Assign input files from the job number: $1"
inFiles="input_files.dat"
inputFile="root://cms-xrd-global.cern.ch/$(grep -E "_$1\.root$" "$inFiles")"
outputName="step1_newGEN.root"
outputNameTmp="step1_newGEN_tmp.root"

echo "Input  file: $inputFile"
echo "Output file1: $outputName"
echo "Output file2: $outputNameTmp"
echo "cmsRun -e -j FrameworkJobReport.xml NewGenParticlesProducer_cfg.py inputFile=$inputFile outputFile=$outputName"
cmsRun -e -j FrameworkJobReport.xml NewGenParticlesProducer_cfg.py inputFile=$inputFile outputFile=$outputName

mv $outputName $outputNameTmp 

echo "cmsRun -e -j FrameworkJobReport.xml GenParticlesLabelChanger_cfg.py inputFile=$outputNameTmp outputFile=$outputName"
cmsRun -e -j FrameworkJobReport.xml GenParticlesLabelChanger_cfg.py inputFile=$outputNameTmp outputFile=$outputName

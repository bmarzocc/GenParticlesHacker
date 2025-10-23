#!/bin/bash    

SEED=$(( $1 + 1 ))
                                                                                                                                                                                   
set -e
BASE=$PWD
RELEASE_BASE=$CMSSW_BASE

export SCRAM_ARCH=slc7_amd64_gcc700
source /cvmfs/cms.cern.ch/cmsset_default.sh

inFiles="input_files_newGEN_preserveAngles.dat"
inputFile="root://cms-xrd-global.cern.ch/$(grep -E "_$1\.root$" "$inFiles")"
pileupFiles="input_pileup_unique_list.txt"
inputPremix="root://cms-xrd-global.cern.ch/$(sed -n "$1p" "$pileupFiles")"

outputNameMINIAOD="step6_stdMINIAOD.root"
outputNameMINIAODTmp="step6_stdMINIAOD_tmp.root"

echo "Input  file: $inputFile"
echo "Input Premix file: $inputPremix"
echo "Output file: $outputNameMINIAOD"

eval `scram project CMSSW_10_6_47`
cd CMSSW_10_6_47/src/
eval `scram runtime -sh`
cp ../../step2_SIM_cfg.py ./
cp ../../step3_GEN-SIM-DIGI_inputNEWGEN_cfg.py ./

echo "cmsRun -e -j FrameworkJobReport.xml step2_SIM_cfg.py jobId=$1 inputFile=$inputFile outputFile=$outputNameMINIAOD"
cmsRun -e -j FrameworkJobReport.xml step2_SIM_cfg.py jobId=$1 inputFile=$inputFile outputFile=$outputNameMINIAOD

mv $outputNameMINIAOD $outputNameMINIAODTmp 

outputNameRAW="step3_stdRAW.root"
echo "cmsRun -e -j FrameworkJobReport.xml step3_GEN-SIM-DIGI_cfg.py jobId=$1 inputFile=$outputNameMINIAODTmp outputFile=$outputNameRAW inputPremix=$inputPremix"
cmsRun -e -j FrameworkJobReport.xml step3_GEN-SIM-DIGI_inputNEWGEN_cfg.py jobId=$1 inputFile=$outputNameMINIAODTmp outputFile=$outputNameRAW inputPremix=$inputPremix


cp $outputNameRAW $BASE/
cp FrameworkJobReport.xml $BASE/
cd $BASE

eval `scram project  CMSSW_10_2_20_UL`
cd  CMSSW_10_2_20_UL/src/
eval `scram runtime -sh`
cp $BASE/step4_HLT_cfg.py ./
cp $BASE/FrameworkJobReport.xml ./
cp $BASE/$outputNameRAW ./

outputNameHLT="step4_stdHLT.root"
echo "cmsRun -e -j FrameworkJobReport.xml step4_HLT_cfg.py jobId=$1 inputFile=$outputNameRAW output=$outputNameHLT"
cmsRun -e -j FrameworkJobReport.xml step4_HLT_cfg.py jobId=$1 inputFile=$outputNameRAW outputFile=$outputNameHLT

cp $outputNameHLT $BASE/
cp FrameworkJobReport.xml $BASE/
cd $BASE

cd CMSSW_10_6_47/src/
eval `scram runtime -sh`
cp $BASE/step5_RECO_cfg.py ./
cp $BASE/FrameworkJobReport.xml ./
cp $BASE/$outputNameHLT ./

outputNameRECO="step5_stdRECO.root"
echo "cmsRun -e -j FrameworkJobReport.xml step5_RECO_cfg.py jobId=$1 inputFile=$outputNameHLT output=$outputNameRECO"
cmsRun -e -j FrameworkJobReport.xml step5_RECO_cfg.py jobId=$1 inputFile=$outputNameHLT outputFile=$outputNameRECO

cp $outputNameRECO $BASE/
cp FrameworkJobReport.xml $BASE/
cd $BASE

eval `scram project  CMSSW_10_6_47_patch1`
cd CMSSW_10_6_47_patch1/src/
eval `scram runtime -sh`
cp $BASE/step6_MINIAOD_cfg.py ./
cp $BASE/FrameworkJobReport.xml ./
cp $BASE/$outputNameRECO ./

echo "cmsRun -e -j FrameworkJobReport.xml step6_MINIAOD_cfg.py jobId=$1 inputFile=$outputNameRECO output=$outputNameMINIAOD"
cmsRun -e -j FrameworkJobReport.xml step6_MINIAOD_cfg.py jobId=$1 inputFile=$outputNameRECO outputFile=$outputNameMINIAOD

cp $outputNameMINIAOD $BASE/
cp FrameworkJobReport.xml $BASE/
cd $BASE











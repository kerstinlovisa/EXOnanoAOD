EGamma_data='/store/data/Run2023D/EGamma0/AOD/PromptReco-v1/000/370/580/00001/ff70112c-b5a0-4011-b054-a423a3cb1464.root'

#cmsRun $CMSSW_BASE/src/PhysicsTools/EXOnanoAOD/test/Run3_2023_PAT_EXONANO_template_data.py maxEvents=1000 inputFile=$EGamma_data outputTag="EGamma"
cmsRun $CMSSW_BASE/src/PhysicsTools/EXOnanoAOD/test/Run3_2023_PAT_EXONANO_template_data.py maxEvents=1000 inputFile=$EGamma_data outputTag="EGamma_v1"
# cmsRun $CMSSW_BASE/src/PhysicsTools/EXOnanoAOD/test/Run3_2023_PAT_EXONANO_template_data.py maxEvents=1000 inputFile=$EGamma_data outputTag="EGamma_ref" EXO="FALSE"

$CMSSW_BASE/src/PhysicsTools/NanoAOD/test/inspectNanoFile.py EXONANO_EGamma_ref.root --size size_EGamma_ref.html -j size_EGamma_ref.json
$CMSSW_BASE/src/PhysicsTools/NanoAOD/test/inspectNanoFile.py EXONANO_EGamma_v1.root --size size_EGamma_v1.html -j size_EGamma_v1.json

# TT_MC='/store/mc/Run3Summer23DRPremix/TTto4Q_TuneCP5_13p6TeV_powheg-pythia8/AODSIM/130X_mcRun3_2023_realistic_v14-v2/2550000/0377b97b-92c5-4fad-9f8c-a5e65e528810.root'

# cmsRun $CMSSW_BASE/src/PhysicsTools/EXOnanoAOD/test/Run3_2023_PAT_EXONANO_MC.py maxEvents=1000 inputFile=$TT_MC outputTag="TT_4Q_v1" 
# cmsRun $CMSSW_BASE/src/PhysicsTools/EXOnanoAOD/test/Run3_2023_PAT_EXONANO_MC.py maxEvents=1000 inputFile=$TT_MC outputTag="TT_4Q_ref"  EXO="FALSE" 

# #$CMSSW_BASE/src/PhysicsTools/NanoAOD/test/inspectNanoFile.py EXONANO_MC_TT_4Q.root --size size_TT.html -j size_TT.json
# $CMSSW_BASE/src/PhysicsTools/NanoAOD/test/inspectNanoFile.py EXONANO_MC_TT_4Q_v1.root --size size_TT_v1.html -j size_TT_v1.json
# $CMSSW_BASE/src/PhysicsTools/NanoAOD/test/inspectNanoFile.py EXONANO_MC_TT_4Q_ref.root --size size_TT_ref.html -j size_TT_ref.json

# MDSskim="file:/eos/cms/store/user/kakwok/HLT/Commissioning2024/mds_nano/Muon0_EXOskim_0.root"

# #cmsRun $CMSSW_BASE/src/PhysicsTools/EXOnanoAOD/test/Run3_2023_PAT_EXONANO_template_data.py maxEvents=1000 inputFile=$MDSskim outputTag="MDSskim" 
# cmsRun $CMSSW_BASE/src/PhysicsTools/EXOnanoAOD/test/Run3_2023_PAT_EXONANO_template_data.py maxEvents=1000 inputFile=$MDSskim outputTag="MDSskim_v1" 
# cmsRun $CMSSW_BASE/src/PhysicsTools/EXOnanoAOD/test/Run3_2023_PAT_EXONANO_template_data.py maxEvents=1000 inputFile=$MDSskim outputTag="MDSskim_ref" EXO="FALSE" 
# #$CMSSW_BASE/src/PhysicsTools/NanoAOD/test/inspectNanoFile.py EXONANO_MDSskim.root --size size_MDSskim.html -j size_MDSskim.json
# $CMSSW_BASE/src/PhysicsTools/NanoAOD/test/inspectNanoFile.py EXONANO_MDSskim_v1.root --size size_MDSskim_v1.html -j size_MDSskim_v1.json
# $CMSSW_BASE/src/PhysicsTools/NanoAOD/test/inspectNanoFile.py EXONANO_MDSskim_ref.root --size size_MDSskim_ref.html -j size_MDSskim_ref.json

# Muon_data='/store/data/Run2024I/Muon1/AOD/PromptReco-v2/000/386/694/00000/00d3c6c0-84fa-44e3-89c6-55aee9ad8bf1.root'

# #cmsRun $CMSSW_BASE/src/PhysicsTools/EXOnanoAOD/test/Run3_2023_PAT_EXONANO_template_data.py maxEvents=1000 inputFile=$Muon_data outputTag="Muon" 
# cmsRun $CMSSW_BASE/src/PhysicsTools/EXOnanoAOD/test/Run3_2023_PAT_EXONANO_template_data.py maxEvents=1000 inputFile=$Muon_data outputTag="Muon_v1" 
# cmsRun $CMSSW_BASE/src/PhysicsTools/EXOnanoAOD/test/Run3_2023_PAT_EXONANO_template_data.py maxEvents=1000 inputFile=$Muon_data outputTag="Muon_ref" EXO="FALSE" 

# #$CMSSW_BASE/src/PhysicsTools/NanoAOD/test/inspectNanoFile.py EXONANO_Muon.root --size size_Muon.html -j size_Muon.json
# $CMSSW_BASE/src/PhysicsTools/NanoAOD/test/inspectNanoFile.py EXONANO_Muon_v1.root --size size_Muon_v1.html -j size_Muon_v1.json
# $CMSSW_BASE/src/PhysicsTools/NanoAOD/test/inspectNanoFile.py EXONANO_Muon_ref.root --size size_Muon_ref.html -j size_Muon_ref.json

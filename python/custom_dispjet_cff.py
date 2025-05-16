import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import *

process.load('PhysicsTools.EXOnanoAOD.displacedInclusiveVertexing_cff')

outputTable = cms.EDProducer("DispJetTableProducer",
                             rho=cms.InputTag("fixedGridRhoFastjetAll"),
                             electrons=cms.InputTag("linkedObjects","electrons"),
                             muons=cms.InputTag("linkedObjects","muons"),
                             jets=cms.InputTag("linkedObjects","jets"),
                             jetsFat=cms.InputTag("slimmedJetsAK8"),
                             jetsSub=cms.InputTag("slimmedJetsAK8PFPuppiSoftDropPacked", "SubJets"),
                             primaryVertex = cms.InputTag("offlineSlimmedPrimaryVertices"),
                             secondaryVertex = cms.InputTag("displacedInclusiveSecondaryVertices")
)

def add_dispJetTables(process):

    process.outputTable = outputTable
    process.exonanoaodTask = cms.Task(process.outputTable)
    process.nanoTableTaskCommon.add(process.exonanoaodTask)
    
    return process

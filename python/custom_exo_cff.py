import FWCore.ParameterSet.Config as cms
from PhysicsTools.NanoAOD.common_cff import *

from RecoMuon.MuonRechitClusterProducer.cscRechitClusterProducer_cfi import cscRechitClusterProducer
from RecoMuon.MuonRechitClusterProducer.dtRechitClusterProducer_cfi import dtRechitClusterProducer 

# MDSnano tables 
cscMDSClusterTable = cms.EDProducer("SimpleMuonRecHitClusterFlatTableProducer",
    src  = cms.InputTag("ca4CSCrechitClusters"),
    name = cms.string("cscMDSHLTCluster"),
    doc  = cms.string("MDS cluster at HLT"),
    variables = cms.PSet(
        eta = Var("eta", float, doc="cluster eta"),
        phi = Var("phi", float, doc="cluster phi"),
        x   = Var("x"  , float, doc="cluster x"),
        y   = Var("y"  , float, doc="cluster y"),
        z   = Var("z"  , float, doc="cluster z"),
        r   = Var("r"  , float, doc="cluster r"),
        size   = Var("size"  , int, doc="cluster size"),
        nStation   = Var("nStation"  , int, doc="cluster nStation"),
        avgStation   = Var("avgStation"  , float, doc="cluster avgStation"),
        nMB1   = Var("nMB1"  , int, doc="cluster nMB1"),
        nMB2   = Var("nMB2"  , int, doc="cluster nMB2"),
        nME11   = Var("nME11"  , int, doc="cluster nME11"),
        nME12   = Var("nME12"  , int, doc="cluster nME12"),
        nME41   = Var("nME41"  , int, doc="cluster nME41"),
        nME42   = Var("nME42"  , int, doc="cluster nME42"),
        time   = Var("time"  , float, doc="cluster time = avg cathode and anode time"),
        timeSpread   = Var("timeSpread"  , float, doc="cluster timeSpread")
    )
)

dtMDSClusterTable = cscMDSClusterTable.clone(
    src = cms.InputTag("ca4DTrechitClusters"),
    name= cms.string("dtMDSHLTCluster")
)

#DSA muon tables
DSAmuonsTable = cms.EDProducer("DSAMuonTableProducer",
    dsaMuons=cms.InputTag("displacedStandAloneMuons"),
    muons=cms.InputTag("linkedObjects","muons"),
    primaryVertex = cms.InputTag("offlineSlimmedPrimaryVertices"),
    beamspot = cms.InputTag("offlineBeamSpot")
)

DSAmuonVertexTable = cms.EDProducer("MuonVertexTableProducer",
    dsaMuons=cms.InputTag("displacedStandAloneMuons"),
    patMuons=cms.InputTag("linkedObjects","muons"),
    beamspot=cms.InputTag("offlineBeamSpot"),
    generalTracks=cms.InputTag("generalTracks"),
    primaryVertex=cms.InputTag("offlineSlimmedPrimaryVertices")
)

PATmuonExtendedTable = cms.EDProducer("MuonExtendedTableProducer",
    name=cms.string("Muon"),
    muons=cms.InputTag("linkedObjects","muons"),
    dsaMuons=cms.InputTag("displacedStandAloneMuons"),
    primaryVertex=cms.InputTag("offlineSlimmedPrimaryVertices"),
    beamspot=cms.InputTag("offlineBeamSpot"),
    generalTracks=cms.InputTag("generalTracks")
)

electronVertexTable = cms.EDProducer("ElectronVertexTableProducer",
    electrons=cms.InputTag("linkedObjects","electrons"),
    beamspot=cms.InputTag("offlineBeamSpot"),
    generalTracks=cms.InputTag("generalTracks"),
    primaryVertex=cms.InputTag("offlineSlimmedPrimaryVertices")
)

electronExtendedTable = cms.EDProducer("ElectronExtendedTableProducer",
    name=cms.string("Electron"),
    electrons=cms.InputTag("linkedObjects","electrons"),
)

def add_mdsTables(process):
    process.ca4CSCrechitClusters = cscRechitClusterProducer    
    process.ca4DTrechitClusters = dtRechitClusterProducer
    process.cscMDSClusterTable = cscMDSClusterTable
    process.dtMDSClusterTable = dtMDSClusterTable 

    process.MDSTask = cms.Task(process.ca4CSCrechitClusters)
    process.MDSTask.add(process.ca4DTrechitClusters)
    process.MDSTask.add(process.cscMDSClusterTable)
    process.MDSTask.add(process.dtMDSClusterTable)

    process.nanoTableTaskCommon.add(process.MDSTask)

    return process

def add_dsamuonTables(process):
    process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
    process.load('TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorAny_cfi')

    process.DSAmuonsTable = DSAmuonsTable
    process.DSAmuonVertexTable = DSAmuonVertexTable
    process.PATmuonExtendedTable = PATmuonExtendedTable

    process.dsamuonTask = cms.Task(process.DSAmuonsTable)
    process.dsamuonVertexTask = cms.Task(process.DSAmuonVertexTable)
    process.patmuonTask = cms.Task(process.PATmuonExtendedTable)

    process.nanoTableTaskCommon.add(process.dsamuonTask)
    process.nanoTableTaskCommon.add(process.dsamuonVertexTask)
    process.nanoTableTaskCommon.add(process.patmuonTask)

    return process

def add_electronVertexTables(process):

    process.load("TrackingTools.TransientTrack.TransientTrackBuilder_cfi")
    process.load('TrackPropagation.SteppingHelixPropagator.SteppingHelixPropagatorAny_cfi')

    process.electronVertexTable = electronVertexTable
    process.electronExtendedTable = electronExtendedTable

    process.electronVertexTask = cms.Task(process.electronVertexTable)
    process.electronExtendedTask = cms.Task(process.electronExtendedTable)

    process.nanoTableTaskCommon.add(process.electronVertexTask)
    process.nanoTableTaskCommon.add(process.electronExtendedTask)

    return process

def update_genParticleTable(process):

    process.genParticleTable.variables.vx = Var("vx",float, doc = "gen particle production vertex x coordinate (cm)", precision=8)
    process.genParticleTable.variables.vy = Var("vy",float, doc = "gen particle production vertex y coordinate (cm)", precision=8)
    process.genParticleTable.variables.vz = Var("vz",float, doc = "gen particle production vertex z coordinate (cm)", precision=8)

    process.genParticleTable.variables.px = Var("px",float, doc = "gen particle momentum x coordinate", precision=8)
    process.genParticleTable.variables.py = Var("py",float, doc = "gen particle momentum y coordinate", precision=8)
    process.genParticleTable.variables.pz = Var("pz",float, doc = "gen particle momentum z coordinate", precision=8)

    return process

def add_exonanoTables(process):

    process = add_mdsTables(process)
    process = add_dsamuonTables(process)
    process = add_electronVertexTables(process)
    process = update_genParticleTable(process)

    return process

def add_exonanoMCTables(process):

    process = update_genParticleTable(process)

    return process

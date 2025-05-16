import FWCore.ParameterSet.Config as cms

unpackedTracksAndVertices = cms.EDProducer('PATTrackAndVertexUnpacker',
                                           slimmedVertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
                                           slimmedSecondaryVertices = cms.InputTag("slimmedSecondaryVertices"),
                                           additionalTracks= cms.InputTag("lostTracks"),
                                           packedCandidates = cms.InputTag("packedPFCandidates"))

# ref: https://github.com/cms-sw/cmssw/blob/2bed69b1658e4deeaef914e462741919e9183be3/RecoVertex/AdaptiveVertexFinder/plugins/InclusiveVertexFinder.h#L48

displacedInclusiveVertexFinder  = cms.EDProducer("InclusiveVertexFinder",
                                                 beamSpot = cms.InputTag("offlineBeamSpot"),
                                                 primaryVertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
                                                 tracks = cms.InputTag("unpackedTracksAndVertices"),
                                                 minHits = cms.uint32(6), # 8
                                                 maximumLongitudinalImpactParameter = cms.double(20), # 0.3
                                                 minPt = cms.double(0.8), # 0.8
                                                 maxNTracks = cms.uint32(100), # 30
                                                 
                                                 clusterizer = cms.PSet(
                                                     seedMax3DIPSignificance = cms.double(9999.), # 9999.
                                                     seedMax3DIPValue = cms.double(9999.), # 9999.
                                                     seedMin3DIPSignificance = cms.double(1.2), # 1.2
                                                     seedMin3DIPValue = cms.double(0.005), # 0.005
                                                     clusterMaxDistance = cms.double(0.4), # 0.05
                                                     clusterMaxSignificance = cms.double(4.5), # 4.5
                                                     distanceRatio = cms.double(20), # 20
                                                     clusterMinAngleCosine = cms.double(0.5), # 0.5
                                                 ),
                                                     
                                                 vertexMinAngleCosine = cms.double(0.95), # 0.95
                                                 vertexMinDLen2DSig = cms.double(2.5), # 2.5
                                                 vertexMinDLenSig = cms.double(0.5), # 0.5
                                                 fitterSigmacut =  cms.double(3), # 3
                                                 fitterTini = cms.double(256), # 256
                                                 fitterRatio = cms.double(0.25), # 0.25
                                                 useDirectVertexFitter = cms.bool(True), # True
                                                 useVertexReco  = cms.bool(True), # True
                                                 vertexReco = cms.PSet(
                                                     finder = cms.string('avr'),
                                                     primcut = cms.double(1.0), # 1.0
                                                     seccut = cms.double(3), # 3
                                                     smoothing = cms.bool(True)) # True
                                            )
                                            
displacedVertexMerger = cms.EDProducer("VertexMerger",
                                       secondaryVertices = cms.InputTag("displacedInclusiveVertexFinder"),
                                       maxFraction = cms.double(0.7),
                                       minSignificance = cms.double(2))

# ref: https://github.com/cms-sw/cmssw/blob/2bed69b1658e4deeaef914e462741919e9183be3/RecoVertex/AdaptiveVertexFinder/plugins/VertexArbitrators.cc#L54

displacedTrackVertexArbitrator = cms.EDProducer("TrackVertexArbitrator",
                                                beamSpot = cms.InputTag("offlineBeamSpot"),
                                                primaryVertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
                                                tracks = cms.InputTag("unpackedTracksAndVertices"),
                                                secondaryVertices = cms.InputTag("displacedVertexMerger"),
                                                dLenFraction = cms.double(0.333), # 0.333
                                                dRCut = cms.double(1), # 0.4
                                                distCut = cms.double(0.1), # 0.04
                                                sigCut = cms.double(5), # 5
                                                fitterSigmacut =  cms.double(3), # 3
                                                fitterTini = cms.double(256), # 256
                                                fitterRatio = cms.double(0.25), # 0.25
                                                trackMinLayers = cms.int32(4), # 4
                                                trackMinPt = cms.double(.4), # 0.4
                                                trackMinPixels = cms.int32(0) # 1
)
    
displacedInclusiveSecondaryVertices = displacedVertexMerger.clone()
displacedInclusiveSecondaryVertices.secondaryVertices = cms.InputTag("displacedTrackVertexArbitrator")
displacedInclusiveSecondaryVertices.maxFraction = 0.2
displacedInclusiveSecondaryVertices.minSignificance = 10
    
displacedInclusiveVertexing = cms.Sequence(unpackedTracksAndVertices * displacedInclusiveVertexFinder  * displacedVertexMerger * displacedTrackVertexArbitrator * displacedInclusiveSecondaryVertices)

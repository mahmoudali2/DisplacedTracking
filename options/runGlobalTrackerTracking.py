# Vertex (barrel + endcap) + Si wrapper (barrel + disks)
#
# CombinatorialTrackFinder on the GLOBAL TRACKER readout of IDEA_o1_v04.
# Pattern recognition AND fitting run in ONE pass over all input subdetectors:
# MergeTrackerHitPlanes combines their hit and reco->sim link collections first,
# so a track crossing several subdetectors is found as a single track.
#
# Readout differences vs the muon system:
#   muon   : system:5,type:-2,layer:2,chamber:13,slice:1,y:-10,z:-10
#   global : system:5,side:-2,layer:5,module:12,sensor:8
# "type" and "side" are both signed 2-bit with the same 0/+1/-1 barrel/+endcap/
# -endcap meaning, so only the NAME differs -> RegionFieldName.
import os
from Gaudi.Configuration import INFO
from k4FWCore import ApplicationMgr, IOSvc
from Configurables import EventDataSvc, GeoSvc
from Configurables import CombinatorialTrackFinder, MergeTrackerHitPlanes

geoservice = GeoSvc("GeoSvc")
geoservice.detectors = [os.path.join(os.environ["K4GEO"],
                        "FCCee/IDEA/compact/IDEA_o1_v04/IDEA_o1_v04.xml")]
geoservice.OutputLevel = INFO
geoservice.EnableGeant4Geo = False

merge = MergeTrackerHitPlanes("MergeGlobalTrackerHits")
merge.InputHitCollections         = ["VTXBDigis", "VTXDDigis", "SiWrBDigis", "SiWrDDigis"]
merge.InputRecoSimLinkCollections = ["VTXBSimDigiLinks", "VTXDSimDigiLinks", "SiWrBSimDigiLinks", "SiWrDSimDigiLinks"]
merge.OutputHitCollection         = ["GlobalTrackerHits"]
merge.OutputRecoSimLinkCollection = ["GlobalTrackerHitSimLinks"]
merge.OutputLevel = INFO

t = CombinatorialTrackFinder()
# -- readout adoption ---------------------------------------------------------
t.EncodingStringParameterName = "GlobalTrackerReadoutID"
t.RegionFieldName = "side"          # "type" for MuonSystemReadoutID
t.LayerFieldName  = "layer"
t.DetectorNames   = ["Vertex", "SiWrB", "SiWrD"]        # surface maps of all contributing DetElements
t.UseSystemInCompositeID = 1        # else VTXB layer 0 and SiWrB layer 0 collide
t.ParticleType    = "muon"
t.OutputLevel     = INFO

t.UseGenFit          = False
t.DoInnerPropagation = False        # hits are already inside the solenoid

# -- detector-scale retuning (vertex is um-scale, the muon system is m-scale) --
t.SigmaHitDefault           = 0.0003   # cm; per-hit du/dv from the digitiser wins
t.MaxPairHitDistance        = 400.0    # cm; wrapper disks reach |z| ~ 2.3 m
t.MaxConsecutiveHitDistance = -1
t.NeighbourTrackMaxDist     = 1.0
t.PairCoincidenceMaxDistMM  = 20.0     # 250 mm is muon-scale, meaningless here

t.MinTrackHits          = 3
t.MinLayerSpan          = 3
t.MaxCombinatorialHits  = 8
t.NHitMaxChi2NDF        = 10.0
t.MaxComboPT            = 200.0
t.OutlierSigma          = 5.0
t.MaxOutlierIterations  = 10
t.MaxHitsPerCompositeID = 0
t.MaxHitEdepKeV         = 0
t.MaxConsecDeltaPhi     = -1           # road cuts are tuned for muon radii
t.MaxComboPhiSpread     = -1

t.InputHitCollection          = ["GlobalTrackerHits"]
t.InputRecoSimLinkCollection  = ["GlobalTrackerHitSimLinks"]
t.OutputTrackCollection       = ["GlobalTrackerTracks"]
t.OutputHitCollection         = ["GlobalTrackerTrackHits"]
t.OutputRecoSimLinkCollection = ["GlobalTrackerTrackHitSimLinks"]

iosvc = IOSvc()
iosvc.Input  = "/eos/user/m/maali/output_files/global_tracker/o1v04_gun_digi_10k.root"
iosvc.Output = "/eos/user/m/maali/output_files/global_tracker/o1v04_gun_tracks_all_10k.root"

ApplicationMgr(TopAlg=[merge, t], EvtSel="NONE", EvtMax=-1,
               ExtSvc=[EventDataSvc("EventDataSvc")], OutputLevel=INFO)

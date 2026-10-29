# CombinatorialTrackFinder on an inner tracker: IDEA_o1_v04 vertex detector (barrel + disks)
# and silicon wrapper (barrel + disks), which share the GlobalTrackerReadoutID readout.
# MergeTrackerHitPlanes first combines their hit and reco->sim link collections, so a
# track crossing several subdetectors is found as one track.
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
# -- readout -------------------------------------------------------------------
t.EncodingStringParameterName = "GlobalTrackerReadoutID"
t.RegionFieldName = "side"          # barrel / endcap field of this readout
t.LayerFieldName  = "layer"
t.DetectorNames   = ["Vertex", "SiWrB", "SiWrD"]        # surface maps of all contributing DetElements
t.UseSystemInCompositeID = 1        # several subdetectors share the readout
t.ParticleType    = "muon"
t.OutputLevel     = INFO

t.UseGenFit          = False
t.DoInnerPropagation = False        # hits are inside the solenoid

# -- stored track states: the fitted helix AtFirstHit and AtLastHit, plus AtIP ----
t.IPStateMode        = "Helix"      # uniform field: the same helix at its closest approach to the origin
t.StoreLastHitState  = True

# -- scale of a silicon tracker ---------------------------------------------------
t.SigmaHitDefault           = 0.0003   # cm; per-hit du/dv from the digitiser wins
t.MaxPairHitDistance        = 400.0    # cm; wrapper disks reach |z| ~ 2.3 m
t.MaxConsecutiveHitDistance = -1
t.NeighbourTrackMaxDist     = 1.0
t.PairCoincidenceMaxDistMM  = 20.0

t.MinTrackHits          = 3
t.MinLayerSpan          = 3
t.MaxCombinatorialHits  = 8
t.NHitMaxChi2NDF        = 10.0
t.MaxComboPT            = 200.0
t.OutlierSigma          = 5.0
t.MaxOutlierIterations  = 10
t.MaxHitsPerCompositeID = 0
t.MaxHitEdepKeV         = 0
t.MaxConsecDeltaPhi     = -1           # no phi road cuts
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

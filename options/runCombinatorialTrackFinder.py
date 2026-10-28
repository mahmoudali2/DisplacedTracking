# Author : Mahmoud Althakeel (mahmoud.althakeel@cern.ch)
#
# Steering file for running the standalone muon-system tracking algorithm
# (pattern recognition + optional GenFit fit over IDEA Muon-System hits).
# This script configures and runs the CombinatorialTrackFinder Gaudi component.
#

import os
from Gaudi.Configuration import INFO
from k4FWCore import ApplicationMgr, IOSvc
from Configurables import EventDataSvc
from Configurables import CombinatorialTrackFinder
from Configurables import GeoSvc

# Geometry service - provides access to detector geometry.
# IMPORTANT: the compact MUST match the geometry the input file was digitised
# with. The muon-system configurations use dedicated compacts, NOT the nominal
# IDEA_o1_v03.xml; feeding a 3L/4L digi file to the nominal compact still
# "works" but wrecks everything downstream of pattern recognition (analytical
# inner propagation 94% -> 35%, charge-sign correctness 77% -> 45%).
#   3L_30cm -> IDEA_3L_30cm.xml      4L_50cm -> IDEA_MS_4L_50cm.xml
# (see IDEA_muon_optimization/gen_tracking_steerings.py for the full mapping)
geoservice = GeoSvc("GeoSvc")
# This compact matches iosvc.Input below (the 3L_30cm muon-gun digi file).
# Change BOTH together if you switch configuration.
geoservice.detectors = [os.environ["K4GEO"] + "/FCCee/IDEA/compact/IDEA_o1_v03/IDEA_3L_30cm.xml"]
geoservice.OutputLevel = INFO
geoservice.EnableGeant4Geo = False

# Combinatorial track finder algorithm
combinatorialTrackFinder = CombinatorialTrackFinder()
combinatorialTrackFinder.DetectorName = "Muon-System"           # Name of the detector to process
combinatorialTrackFinder.ParticleType = "muon"                  # Particle type for material effects
combinatorialTrackFinder.EncodingStringParameterName = "MuonSystemReadoutID"
combinatorialTrackFinder.OutputLevel = INFO

# ── Readout schema ───────────────────────────────────────────────────────────
# The component is readout-agnostic: only the cellID FIELD NAMES are assumed.
#   MuonSystemReadoutID   : system:5,type:-2,layer:2,chamber:13,slice:1,y:-10,z:-10
#   GlobalTrackerReadoutID: system:5,side:-2,layer:5,module:12,sensor:8
# "type" and "side" are both signed 2-bit with the same 0 / +1 / -1 meaning
# barrel / +endcap / -endcap, so only the name differs. initialize() validates
# these against the actual encoding and fails loudly on a mismatch.
combinatorialTrackFinder.RegionFieldName = "type"               # "side" for GlobalTrackerReadoutID
combinatorialTrackFinder.LayerFieldName  = "layer"
# Fold the cellID "system" field into the per-layer key. -1 = auto (on when
# DetectorNames lists more than one detector), 0 = off, 1 = on. Required when a
# readout spans several subdetectors, else e.g. VertexBarrel layer 0 and SiWrB
# layer 0 collapse onto the same key. Irrelevant for the single-DetElement
# muon system, where it is a constant offset.
combinatorialTrackFinder.UseSystemInCompositeID = -1
# DetectorNames merges the surface maps of several DetElements sharing one
# readout (GlobalTrackerReadoutID covers Vertex, SiWrB, SiWrD). Leave empty to
# use DetectorName alone, as here.
# combinatorialTrackFinder.DetectorNames = ["Vertex", "SiWrB", "SiWrD"]

# See runGlobalTrackerTracking.py in this directory for the full
# GlobalTrackerReadoutID configuration, including MergeTrackerHitPlanes.

# GenFit configuration
combinatorialTrackFinder.UseGenFit = False                      # Enable GenFit track fitting
combinatorialTrackFinder.MaxFitIterations = 100                 # Maximum iterations for the Kalman fit
# Inner propagation (disabled: solenoid boundary extrapolation not yet calibrated)
combinatorialTrackFinder.DoInnerPropagation = False
combinatorialTrackFinder.InnerPropTargetRadius = 100.0          # target radius in cm if enabled

# ── N-hit combinatorial strategy ─────────────────────────────────────────────
combinatorialTrackFinder.MinTrackHits = 3                       # minimum hits to form a track
combinatorialTrackFinder.MaxCombinatorialHits = 8               # cap on k in C(N,k) enumeration
combinatorialTrackFinder.NHitMaxChi2NDF = 10.0                  # max chi2/NDF of N-hit circle fit
combinatorialTrackFinder.HitIsolationCut = 0                    # min cross-region distance (cm); 0 to disable
combinatorialTrackFinder.MinLayerSpan = 3                       # min distinct detector layers per combo
combinatorialTrackFinder.MaxComboPT = 200.0                     # GeV: maximum combo pT (excludes combos above this threshold during search)
combinatorialTrackFinder.OutlierSigma = 5.0                     # outlier removal threshold (σ units)
combinatorialTrackFinder.SigmaHitDefault = 0.04                 # default hit resolution (cm = 0.4 mm)
combinatorialTrackFinder.MaxOutlierIterations = 10              # max rounds of outlier removal
# Road cuts: guard against cross-track hit mixing
combinatorialTrackFinder.MaxConsecDeltaPhi = 0.15               # max |Δφ| between consecutive layers (rad); -1 to disable
combinatorialTrackFinder.MaxComboPhiSpread = 0.5                # max total φ spread across combo (rad); -1 to disable
# Distance cuts: guard against ghost combos spanning different detector regions
combinatorialTrackFinder.MaxConsecutiveHitDistance = -1         # max 3D dist between adjacent (inner->outer) hits (cm); -1 to disable
combinatorialTrackFinder.MaxPairHitDistance = 350.0             # max 3D dist between any two hits in a combo (cm); 300 accommodates barrel->endcap
# Hit quality pre-filter and proximity-based combo selection
combinatorialTrackFinder.MaxHitEdepKeV = 0                      # reject hits with edep > 10 keV (delta-rays / hadronic secondaries)
combinatorialTrackFinder.ProximityScoreWeight = 1.0             # weight for edep-homogeneity in combo ranking (0 to disable)
combinatorialTrackFinder.EdepNormKeV = 3.0                      # normalisation: 2 keV avg deviation ~ +1 chi2 unit (for k>=4 combined score)
combinatorialTrackFinder.MaxHitsPerCompositeID = 3              # 0 = use all hits; >0 caps per layer to top-N by proximity score
# conditions for second track in the events
combinatorialTrackFinder.NeighbourTrackMaxDist = 10.0           # cm (= 50 mm)
combinatorialTrackFinder.NeighbourTrackMinHits = 5              # must have >=4 hits
combinatorialTrackFinder.NeighbourTrackMaxChi2NDF = 3.0         # must have chi2/NDF <= 3
# Max 3D distance (mm) between two hits treated as a paired-layer doublet.
# 250 mm suits muon-system paired layers (12-24 cm apart); scale it down for
# inner-tracker readouts where layers sit millimetres apart.
combinatorialTrackFinder.PairCoincidenceMaxDistMM = 250.0
# Guided crowded-layer hit selection
combinatorialTrackFinder.doCrowdedLayerHitSelection = False     # enable guided crowded-layer hit selection

# ── Collections ──────────────────────────────────────────────────────────────
combinatorialTrackFinder.InputHitCollection = ["MSTrackerHits"]
combinatorialTrackFinder.OutputTrackCollection = ["StandaloneMuonTracks"]
combinatorialTrackFinder.OutputHitCollection = ["StandaloneMuonTrackHits"]
combinatorialTrackFinder.InputRecoSimLinkCollection = ["MSTrackerHitRelations"]
combinatorialTrackFinder.OutputRecoSimLinkCollection = ["StandaloneMuonTrackHitSimLinks"]

# Input/Output service
iosvc = IOSvc()
iosvc.Input = "/eos/user/m/maali/output_files/performance/3L_uGun_IDEA_MS_digi_100k.root"
iosvc.Output = "test_stdMuonTracks.root"

# Application manager configuration
ApplicationMgr(
    TopAlg=[combinatorialTrackFinder],
    EvtSel="NONE",
    EvtMax=-1,
    ExtSvc=[EventDataSvc("EventDataSvc")],
    OutputLevel=INFO
)

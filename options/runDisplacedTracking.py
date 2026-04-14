# Author : Mahmoud Al-Thakeel
#
# Steering file for running the displaced tracking algorithm
# This script configures and runs the DisplacedTracking Gaudi component
#

import os
from Gaudi.Configuration import INFO
from k4FWCore import ApplicationMgr, IOSvc
from Configurables import EventDataSvc
from Configurables import DisplacedTracking
from Configurables import GeoSvc

# Geometry service - provides access to detector geometry
geoservice = GeoSvc("GeoSvc")
geoservice.detectors = ["./k4geo/FCCee/IDEA/compact/IDEA_o1_v03/IDEA_o1_v03+.xml"]
geoservice.OutputLevel = INFO
geoservice.EnableGeant4Geo = False

# displaced tracker algorithm
displaced = DisplacedTracking()
displaced.DetectorName = "Muon-System"                   # Name of the detector to process
displaced.ParticleType = "muon"                          # Particle type for material effects
displaced.EncodingStringParameterName = "MuonSystemReadoutID"
displaced.OutputLevel = 0                                # DEBUG level -> 1 , Info level -> 0
# GenFit configuration
displaced.UseGenFit = False                               # Enable GenFit track fitting
displaced.MaxFitIterations = 100                         # Maximum iterations for the Kalman fit
#displaced.DebugLevel = 0                                # GenFit internal debug level
# Inner propagation (disabled: solenoid boundary extrapolation not yet calibrated)
displaced.DoInnerPropagation = False
displaced.InnerPropTargetRadius = 100.0                  # target radius in cm if enabled

# ── N-hit combinatorial strategy ─────────────────────────────────────────────
displaced.MinTrackHits = 4                               # minimum hits to form a track
displaced.MaxCombinatorialHits = 8                       # cap on k in C(N,k) enumeration
displaced.NHitMaxChi2NDF = 10.0                           # max chi2/NDF of N-hit circle fit
displaced.HitIsolationCut = 0                          # min cross-region distance (cm); 0 to disable
displaced.MinLayerSpan = 4                               # min distinct detector layers per combo
displaced.MaxComboPT = 200.0                             # GeV: maximum combo pT (excludes combos above this threshold during search)
displaced.OutlierSigma = 5.0                             # outlier removal threshold (σ units)
displaced.SigmaHitDefault = 0.04                         # default hit resolution (cm = 0.4 mm)
displaced.MaxOutlierIterations = 10                       # max rounds of outlier removal
# Road cuts: guard against cross-track hit mixing
displaced.MaxConsecDeltaPhi = 0.15                        # max |Δφ| between consecutive layers (rad); -1 to disable
displaced.MaxComboPhiSpread = 0.5                        # max total φ spread across combo (rad); -1 to disable
# Distance cuts: guard against ghost combos spanning different detector regions
displaced.MaxConsecutiveHitDistance = -1            # max 3D dist between adjacent (inner→outer) hits (cm); -1 to disable
displaced.MaxPairHitDistance = 350.0                     # max 3D dist between any two hits in a combo (cm); 300 accommodates barrel→endcap
# Hit quality pre-filter and proximity-based combo selection
displaced.MaxHitEdepKeV = 0                           # reject hits with edep > 10 keV (delta-rays / hadronic secondaries)
displaced.ProximityScoreWeight = 1.0                     # weight for edep-homogeneity in combo ranking (0 to disable)
displaced.EdepNormKeV = 3.0                              # normalisation: 2 keV avg deviation ≈ +1 chi2 unit (for k>=4 combined score)
displaced.MaxHitsPerCompositeID = 3                      # 0 = use all hits; >0 caps per layer to top-N by proximity score
#conditions for second track in the events
displaced.NeighbourTrackMaxDist    = 10.0                 # cm (= 50 mm, as requested)
displaced.NeighbourTrackMinHits    = 5                   # must have ≥4 hits
displaced.NeighbourTrackMaxChi2NDF = 3.0                 # must have chi2/NDF ≤ 3
#Guided crowded-layer hit selection
displaced.doCrowdedLayerHitSelection = False               # enable Guided crowded-layer hit selection
# ── Collections ──────────────────────────────────────────────────────────────
displaced.InputHitCollection = ["MSTrackerHits"]
displaced.OutputTrackCollection = ["DisplacedTracks"]
displaced.OutputHitCollection = ["DisplacedTrackHits"]
displaced.InputRecoSimLinkCollection = ["MSTrackerHitRelations"]
displaced.OutputRecoSimLinkCollection = ["DisplacedTrackHitSimLinks"]

# Input/Output service
iosvc = IOSvc()
iosvc.Input = "output_files/fake_rate/6L_o2_50cm_10K_digi.root"
iosvc.Output = "output_files/fake_rate/6L_o2_50cm_10K_tracks_test.root"

# Application manager configuration
ApplicationMgr(
    TopAlg=[displaced],
    EvtSel="NONE",
    EvtMax=1,
    ExtSvc=[EventDataSvc("EventDataSvc")],
    OutputLevel=INFO
)

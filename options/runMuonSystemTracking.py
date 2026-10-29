# Author : Mahmoud Althakeel (mahmoud.althakeel@cern.ch)
#
# CombinatorialTrackFinder on the IDEA muon system: tracks are found from the muon-chamber
# hits alone (no inner tracker), so muons from decays far from the IP are found too.
# Settings are those of the HNL muon-system studies (IDEA_o1_v05, 4-layer muon system).
#
# Stored track states (all EDM4hep convention):
#   AtFirstHit  the fitted helix at the innermost muon hit
#   AtLastHit   the same helix at the outermost muon hit
#   AtIP        the track at its closest approach to the origin: the outer helix is followed
#               back to the solenoid and matched there to a helix in the inner field
#               (IPStateMode = SolenoidCrossing)
# Note: the muon chambers sit in the return field (B < 0), so the sign of omega is MINUS
# the charge of the muon (omega > 0 means clockwise seen from +z).
#
# usage:
#   k4run runMuonSystemTracking.py --inputFile digi.root --outputFile tracks.root \
#         --compact <compact the input was digitised with> [--numEvents N] [--minHits 3] [--genfit]
#   compact examples: <k4geo>/FCCee/IDEA/compact/IDEA_o1_v05/IDEA_o1_v05.xml
#                     <k4geo>/FCCee/IDEA/compact/IDEA_o1_v03/IDEA_3L_30cm.xml

from Gaudi.Configuration import INFO
from k4FWCore import ApplicationMgr, IOSvc
from k4FWCore.parseArgs import parser
from Configurables import EventDataSvc, GeoSvc
from Configurables import CombinatorialTrackFinder

parser.add_argument("--inputFile", required=True, help="digitised file with MSTrackerHits")
parser.add_argument("--outputFile", required=True)
parser.add_argument("--compact", required=True,
                    help="geometry compact; MUST be the one the input file was digitised with")
parser.add_argument("--numEvents", type=int, default=-1, help="-1 = all events")
parser.add_argument("--minHits", type=int, default=3,
                    help="minimum hits (and layers) per track; 3 for a 3-4 layer system, 4 for 5-6 layers")
parser.add_argument("--genfit", action="store_true", help="also fit the tracks with GenFit (state AtOther)")
args = parser.parse_args()

# ── Geometry ─────────────────────────────────────────────────────────────────
geoservice = GeoSvc("GeoSvc")
geoservice.detectors = [args.compact]
geoservice.OutputLevel = INFO
geoservice.EnableGeant4Geo = False

t = CombinatorialTrackFinder("MuonSystemTracking")
t.OutputLevel = INFO

# ── Detector and readout ─────────────────────────────────────────────────────
t.DetectorName = "Muon-System"
t.EncodingStringParameterName = "MuonSystemReadoutID"
t.RegionFieldName = "type"           # 0 barrel, +1/-1 endcaps
t.LayerFieldName = "layer"
t.UseSystemInCompositeID = -1        # one detector: no "system" in the layer key
t.ParticleType = "muon"              # mass hypothesis for material effects
t.SolenoidRadius = 2300.0            # mm, start of the search for the field sign change
t.SolenoidHalfZ = 2180.0             # mm

# ── Stored track states ──────────────────────────────────────────────────────
t.StoreLastHitState = True
t.IPStateMode = "SolenoidCrossing"   # AtIP from the helix matched at the solenoid boundary

# ── Optional fits and propagation ────────────────────────────────────────────
t.UseGenFit = args.genfit
t.MaxFitIterations = 100
t.DoInnerPropagation = False         # RK4 through the solenoid (printed only)
t.InnerPropTargetRadius = 100.0      # cm

# ── Track finding: all k-hit combinations over the layers ────────────────────
t.Strategy = "Exhaustive"
t.MinTrackHits = args.minHits
t.MinLayerSpan = args.minHits        # hits must be on this many different layers
t.MaxCombinatorialHits = 8
t.NHitMaxChi2NDF = 10.0              # circle-fit chi2/ndf cut
t.MaxComboPT = 200.0                 # GeV; higher-pT combinations only as fallback
t.SigmaHitDefault = 0.04             # cm, used when a hit has no resolution
t.OutlierSigma = 5.0                 # outlier removal at 5 x SigmaHitDefault
t.MaxOutlierIterations = 10
t.HitIsolationCut = 0                # off

# road cuts: hits of one track stay in a narrow azimuth window
t.MaxConsecDeltaPhi = 0.15           # rad, between consecutive layers
t.MaxComboPhiSpread = 0.5            # rad, over the whole track
t.MaxConsecutiveHitDistance = -1     # off
t.MaxPairHitDistance = 350.0         # cm, allows barrel -> endcap tracks

# energy-deposit cuts against delta rays and hadronic secondaries
t.MaxHitEdepKeV = 10                 # hits above 10 keV are not used
t.ProximityScoreWeight = 1.0         # weight of the energy-deposit score in the ranking
t.EdepNormKeV = 3.0
t.MaxHitsPerCompositeID = 0          # all hits of a layer are tried

# stricter cuts for a track next to an accepted one (delta-ray shadows)
t.NeighbourTrackMaxDist = 10.0       # cm
t.NeighbourTrackMinHits = 4
t.NeighbourTrackMaxChi2NDF = 3.0
t.doCrowdedLayerHitSelection = False

# ── Collections ──────────────────────────────────────────────────────────────
t.InputHitCollection = ["MSTrackerHits"]
t.InputRecoSimLinkCollection = ["MSTrackerHitRelations"]
t.OutputTrackCollection = ["StandaloneMuonTracks"]
t.OutputHitCollection = ["StandaloneMuonTrackHits"]
t.OutputRecoSimLinkCollection = ["StandaloneMuonTrackHitSimLinks"]

iosvc = IOSvc()
iosvc.Input = args.inputFile
iosvc.Output = args.outputFile

ApplicationMgr(
    TopAlg=[t],
    EvtSel="NONE",
    EvtMax=args.numEvents,
    ExtSvc=[EventDataSvc("EventDataSvc")],
    OutputLevel=INFO,
)

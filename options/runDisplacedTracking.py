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
geoservice.detectors = ["./k4geo/FCCee/IDEA/compact/IDEA_o1_v03/IDEA_o1_v03.xml"]
geoservice.OutputLevel = INFO
geoservice.EnableGeant4Geo = False

# displaced tracker algorithm
displaced = DisplacedTracking()
displaced.DetectorName = "Muon-System"                   # Name of the detector to process      
displaced.MaxChi2 = 100.0                                # Maximum chi-square for hit acceptance
displaced.MaxDist = 160.0                                # Maximum distance between two hits (cm)
displaced.MinCosAngle2d = 0.5                            # Minimum cos(angle) between consecutive hit-segment vectors in the transverse (xy) plane for triplet seeding.
displaced.MaxSeedPT = 200.0                              # GeV Maximum Pt for the seed 
displaced.ParticleType = "muon"                          # Particle type for material effects
displaced.EncodingStringParameterName = "MuonSystemReadoutID"
displaced.OutputLevel = 0                               # DEBUG level -> 1 , Info level -> 0
# GenFit configuration options
displaced.UseGenFit = True                                  # Enable GenFit track fitting
displaced.MaxFitIterations = 100                            # Maximum iterations for the displaced fit
#displaced.DebugLevel = 0                                    # DEBUG level of GenFit
#Inner propagation settings (Inner track segment)
displaced.DoInnerPropagation = False                         # disable inner propagation
displaced.InnerPropTargetRadius = 100.0                       # change target radius in cm

displaced.InputHitCollection = ["MSTrackerHits"]                     # Input digitized hits from DDPlanarDigi
displaced.OutputTrackCollection = ["DisplacedTracks"]                   # Output reconstructed tracks
displaced.OutputHitCollection = ["DisplacedTrackHits"]                  # Output hits used in tracks (for PODIO relations)
displaced.InputRecoSimLinkCollection = ["MSTrackerHitRelations"] # Input reco->sim links from digitizer
displaced.OutputRecoSimLinkCollection = ["DisplacedTrackHitSimLinks"]    # Output reco->sim links for track hits

# Input/Output service
iosvc = IOSvc()
iosvc.Input = "output_files/performance/3L_uGun_IDEA_MS_digi_100k.root"            # Input file with digitized hits
iosvc.Output = "output_files/pt_res/3L_30cm_uGun_100k_tracks.root"              # Output file for reconstructed tracks

# Application manager configuration
ApplicationMgr(
    TopAlg=[displaced],                         # Algorithms to run
    EvtSel="NONE",                              # Event selection policy
    EvtMax=-1,                                   # Process all events (-1)
    ExtSvc=[EventDataSvc("EventDataSvc")],      # Event data service
    OutputLevel=INFO                            # Logging level
)
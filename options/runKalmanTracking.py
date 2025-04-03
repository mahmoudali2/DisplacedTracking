# Author : Mahmoud Al-Thakeel
# 
# Steering file for running the Kalman tracking algorithm
# This script configures and runs the KalmanTracking Gaudi component
#

import os
from Gaudi.Configuration import INFO
from k4FWCore import ApplicationMgr, IOSvc
from Configurables import EventDataSvc
from Configurables import KalmanTracking
from Configurables import GeoSvc

# Geometry service - provides access to detector geometry
geoservice = GeoSvc("GeoSvc")
geoservice.detectors = ["./k4geo/FCCee/IDEA/compact/IDEA_o1_v03/IDEA_o1_v03.xml"]
geoservice.OutputLevel = INFO
geoservice.EnableGeant4Geo = False

# Kalman tracker algorithm
kalman = KalmanTracking()
kalman.DetectorName = "Muon-System"                   # Name of the detector to process
kalman.InputHitCollection = "MSTrackerHits"           # Input digitized hits from DDPlanarDigi
kalman.OutputTrackCollection = "KalmanTracks"         # Output reconstructed tracks
kalman.MaxChi2 = 15.0                                 # Maximum chi-square for hit acceptance
kalman.MaxRadius = 100000000000.0                     # Maximum radius of track helix circle in x-y
kalman.ParticleType = "muon"                          # Particle type for material effects
kalman.InitialMomentum = 1.0                          # Initial momentum estimate for seeding (GeV)
kalman.MaxDistanceToSurface = 10.0                    # Maximum distance to consider surface intersections (mm)
kalman.EncodingStringParameterName = "MuonSystemReadoutID"
kalman.OutputLevel = 1  # DEBUG level

# Input/Output service
iosvc = IOSvc()
iosvc.Input = "output_digi.root"                # Input file with digitized hits
iosvc.Output = "output_tracks.root"             # Output file for reconstructed tracks

# Application manager configuration
ApplicationMgr(
    TopAlg=[kalman],                            # Algorithms to run
    EvtSel="NONE",                              # Event selection policy
    EvtMax=-1,                                  # Process all events (-1)
    ExtSvc=[EventDataSvc("EventDataSvc")],      # Event data service
    OutputLevel=INFO                            # Logging level
)
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
displaced.MaxDist = 150.0                                # Maximum distance between two hits (cm)
displaced.ParticleType = "muon"                          # Particle type for material effects
displaced.EncodingStringParameterName = "MuonSystemReadoutID"
displaced.OutputLevel = 1                                # DEBUG level
# GenFit configuration options
displaced.UseGenFit = True                                # Enable GenFit track fitting
displaced.MaxFitIterations = 8                            # Maximum iterations for the displaced fit
displaced.DebugLevel = 0                                  # DEBUG level of GenFit

displaced.InputHitCollection = ["MSTrackerHits"]         # Input digitized hits from DDPlanarDigi
displaced.OutputTrackCollection = ["DisplacedTracks"]       # Output reconstructed tracks

# Input/Output service
iosvc = IOSvc()
iosvc.Input = "output_digi_10k_10GeV.root"            # Input file with digitized hits
iosvc.Output = "output_tracks_test.root"              # Output file for reconstructed tracks

# Application manager configuration
ApplicationMgr(
    TopAlg=[displaced],                         # Algorithms to run
    EvtSel="NONE",                              # Event selection policy
    EvtMax=-1,                                  # Process all events (-1)
    ExtSvc=[EventDataSvc("EventDataSvc")],      # Event data service
    OutputLevel=INFO                            # Logging level
)

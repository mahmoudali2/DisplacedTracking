#include <iostream>
#include <iomanip>
#include "DisplacedTracking.h"
#include "DDRec/SurfaceManager.h"
#include "DDRec/SurfaceHelper.h"
#include "DD4hep/DD4hepUnits.h"
#include "DDRec/MaterialManager.h"
#include "podio/Frame.h"
#include "k4FWCore/PodioDataSvc.h"
#include "k4FWCore/DataHandle.h"

// Some physical constants
constexpr double electron_mass = 0.000511;  // GeV
constexpr double muon_mass = 0.10566;       // GeV
constexpr double pion_mass = 0.13957;       // GeV
constexpr double kaon_mass = 0.49368;       // GeV
constexpr double proton_mass = 0.93827;     // GeV
constexpr double c_light = 299792458.0;     // m/s
constexpr double qe = 1.602176634e-19;      // coulombs

//------------------------------------------------------------------------------
// DisplacedTracking Implementation
//------------------------------------------------------------------------------

// Constructor
DisplacedTracking::DisplacedTracking(const std::string& name, ISvcLocator* pSvcLocator)
    : MultiTransformer(name, pSvcLocator,
        {
            KeyValues{"InputHitCollection",       {"TrackerHits"}},
            KeyValues{"HeaderCollectionName",     {"EventHeader"}},
            KeyValues{"InputRecoSimLinkCollection", {"TrackerHitSimTrackerHitLink"}}
        },
        {
            KeyValues{"OutputTrackCollection",    {"KalmanTracks"}},
            KeyValues{"OutputHitCollection",      {"TrackHits"}},
            KeyValues{"OutputRecoSimLinkCollection", {"DisplacedTrackHitSimLinks"}}
        }) {
    
    // Initialize particle map with common particles
    m_particleMap["electron"] = ParticleProperties(electron_mass, -1.0, "electron");
    m_particleMap["positron"] = ParticleProperties(electron_mass, 1.0, "positron");
    m_particleMap["muon"] = ParticleProperties(muon_mass, -1.0, "muon");
    m_particleMap["antimuon"] = ParticleProperties(muon_mass, 1.0, "antimuon");
    m_particleMap["pion"] = ParticleProperties(pion_mass, 1.0, "pion");
    m_particleMap["kaon"] = ParticleProperties(kaon_mass, 1.0, "kaon");
    m_particleMap["proton"] = ParticleProperties(proton_mass, 1.0, "proton");
    m_particleMap["antiproton"] = ParticleProperties(proton_mass, -1.0, "antiproton");
}

// Initialize method
StatusCode DisplacedTracking::initialize() {
    StatusCode sc = Algorithm::initialize();
    if (!sc.isSuccess()) return sc;

    // Get detector service
    if (!m_geoSvc.retrieve()) {
        error() << "Failed to retrieve GeoSvc" << endmsg;
        return StatusCode::FAILURE;
    }

    // Get detector
    m_detector = m_geoSvc->getDetector();
    if (!m_detector) {
        error() << "Failed to access detector" << endmsg;
        return StatusCode::FAILURE;
    }

    // Get BitFieldCoder for cell ID decoding
    std::string cellIDEncodingString = m_geoSvc->constantAsString(m_encodingStringParameter);
    m_bitFieldCoder = std::make_unique<dd4hep::DDSegmentation::BitFieldCoder>(cellIDEncodingString);
    
    // Get surface map
    auto surfaceMan = m_detector->extension<dd4hep::rec::SurfaceManager>();
    m_surfaceMap = surfaceMan->map(m_detectorName);
    if (!m_surfaceMap) {
        error() << "Could not find surface map for detector: " << m_detectorName << endmsg;
        return StatusCode::FAILURE;
    }

    // Get magnetic field
    m_field = m_detector->field();
   
    dd4hep::Position center(0, 0, 0);
    // check the field strength at different point at the detector
    dd4hep::Direction bfield = m_field.magneticField(center);
    info() << "Magnetic field at detector center: ("
        << bfield.x() / dd4hep::tesla << ", "
        << bfield.y() / dd4hep::tesla << ", "
        << bfield.z() / dd4hep::tesla << ") Tesla" << endmsg;

    // Check field at several points
    const double solenoidRadius = 230.0; // cm
    const std::vector<dd4hep::Position> testPoints = { // in cm
        dd4hep::Position(0, 0, 0),  
        dd4hep::Position(100, 100, 100), 
        dd4hep::Position(solenoidRadius, 0, 0),
        dd4hep::Position(0, solenoidRadius, 0),
        dd4hep::Position(solenoidRadius/sqrt(2), solenoidRadius/sqrt(2), 0),           
        dd4hep::Position(229, 0, 0),
        dd4hep::Position(0, 231, 0),
        dd4hep::Position(0, 0, 219),
        dd4hep::Position(0, 0, 218),
        dd4hep::Position(0, 0, -217),
        dd4hep::Position(600, 500, 200),
        dd4hep::Position(500, 500, 250),
        dd4hep::Position(500, 500, -250),
        dd4hep::Position(100, 300, 350),
        dd4hep::Position(500, 500, 500),
        dd4hep::Position(-500, -500, -500)     
    };

    info() << "Magnetic field at different positions:" << endmsg;
    for (const auto& point : testPoints) {
        dd4hep::Direction bfieldVec = m_field.magneticField(point);
        info() << "  Position (" << point.x() << ", " << point.y() << ", " << point.z() 
            << ") cm: B = (" << bfieldVec.x() / dd4hep::tesla << ", "
            << bfieldVec.y() / dd4hep::tesla << ", "
            << bfieldVec.z() / dd4hep::tesla << ") Tesla" << endmsg;
    }
    
    // Create material manager
    dd4hep::Volume worldVol = m_detector->worldVolume();
    m_materialManager = new dd4hep::rec::MaterialManager(worldVol);
    //m_materialManager = new dd4hep::rec::MaterialManager(*m_detector);

    // Get all surfaces and group them by layer 
    // Another possible and simpler solution is : to call surface by layer directly from surfaceMap,
    // but for the moment we keep the following effeicient way of surfaceList.
    dd4hep::DetElement world = m_detector->world();
    dd4hep::rec::SurfaceHelper surfHelper(world);
    //dd4hep::rec::SurfaceHelper surfHelper(*m_detector);

    const dd4hep::rec::SurfaceList& surfList = surfHelper.surfaceList();

    // Get all surfaces and group them by layer 
    // Use only Muon-System surfaces from the surface map, not all world surfaces
    m_surfaces.clear();
    for (const auto& [cellID, surface] : *m_surfaceMap) {
        // Cast from ISurface* to Surface*
        const dd4hep::rec::Surface* concreteSurface = 
            dynamic_cast<const dd4hep::rec::Surface*>(surface);
        if (concreteSurface) {
            m_surfaces.push_back(concreteSurface);
        }
    }
    
    // Group surfaces by layer
    m_surfacesByLayer = getSurfacesByLayer();
    
    // Set particle properties based on configuration
    auto it = m_particleMap.find(m_particleType);
    if (it != m_particleMap.end()) {
        m_particleProperties = it->second;
        info() << "Using particle type: " << m_particleType 
               << " (mass=" << m_particleProperties.mass 
               << " GeV, charge=" << m_particleProperties.charge << ")" << endmsg;
    } else {
        warning() << "Unknown particle type: " << m_particleType 
                 << ", defaulting to pion" << endmsg;
    }
    
    info() << "Found " << m_surfaces.size() << " detector surfaces in " 
           << m_surfacesByLayer.size() << " layers" << endmsg;

    // Initialize GenFit if enabled
    if (m_useGenFit) {
        try {   
            // ---- Field ----
            // Use DD4hepFieldAdapter which wraps the DD4hep OverlayedField directly.
            // This is the authoritative field and properly handles the field flip at R=230 cm.
            // We store the adapter as a member so it lives for the algorithm lifetime.
            m_genFitField = std::make_unique<DD4hepFieldAdapter>(m_field);
            genfit::FieldManager::getInstance()->init(m_genFitField.get());

            // ---- Material ----
            // Use DD4hepMaterialAdapter wrapping the DD4hep MaterialManager.
            // This is consistent with the DD4hep geometry used for surface lookup.
            m_genFitMaterial = std::make_unique<DD4hepMaterialAdapter>(m_materialManager);
            genfit::MaterialEffects* matEff = genfit::MaterialEffects::getInstance();
            matEff->init(m_genFitMaterial.get());

            // Set material effects options
            matEff->setMscModel("GEANE");          // Multiple-scattering model
            matEff->setEnergyLossBetheBloch(true);
            matEff->setNoiseBetheBloch(true);
            matEff->setDebugLvl(0);

            info() << "GenFit initialized with DD4hep field+material adapters" << endmsg;
        } catch (const genfit::Exception& e) {
            error() << "Error initializing GenFit: " << e.what() << endmsg;
            return StatusCode::FAILURE;
        }
    }

    return StatusCode::SUCCESS;
}

// Operator method
std::tuple<edm4hep::TrackCollection,
               edm4hep::TrackerHitPlaneCollection,
               edm4hep::TrackerHitSimTrackerHitLinkCollection> DisplacedTracking::operator()(
    const edm4hep::TrackerHitPlaneCollection& hits,
    const edm4hep::EventHeaderCollection& headers,
    const edm4hep::TrackerHitSimTrackerHitLinkCollection& recoSimLinks) const {
    
    // Find tracks using the direct EDM4hep approach
    edm4hep::TrackCollection trackCollection;
    edm4hep::TrackerHitPlaneCollection outputHits;
    edm4hep::TrackerHitSimTrackerHitLinkCollection outputLinks;

    m_statTotalEvents++;

    // Per-event disambiguation flags (set inside findTracks, counted once per event)
    bool evtUsedPtInner    = false;
    bool evtUsedPtInner2   = false;
    bool evtUsedPtFallback = false;

    // Print event separator and basic info
    info() << "\n" << std::string(80, '=') << endmsg;
    
    // Get event number from the header
    unsigned int eventNumber = 0;
    if (!headers.empty()) {
        eventNumber = headers[0].getEventNumber();
        info() << "Processing Event #" << eventNumber << endmsg;
    } else {
        info() << "Processing new event (unknown number)" << endmsg;
    }
    
    info() << std::string(80, '-') << endmsg;
    info() << "Processing " << hits.size() << " tracker hits" << endmsg;
    
    // Print some details about hits
    if (msgLevel(MSG::DEBUG)) {
        debug() << "Hit details:" << endmsg;
        for (size_t i = 0; i < hits.size(); ++i) {
            const auto& hit = hits[i];
            const auto& pos = hit.getPosition();
            debug() << "  Hit " << i << " at ("
                    << pos[0] << ", " << pos[1] << ", " << pos[2] << ") mm" << endmsg;
        }
    }
    
    // Check if we have enough hits for tracking
    if (hits.size() < 3) {
        warning() << "Not enough hits to create tracks. Need at least 3." << endmsg;
        info() << std::string(80, '=') << "\n" << endmsg; // Bottom separator
        return std::make_tuple(std::move(trackCollection), std::move(outputHits), std::move(outputLinks)); // Return empty collections
    }
    
    try {
        findTracks(&hits, trackCollection, outputHits, recoSimLinks, outputLinks,
                   evtUsedPtInner, evtUsedPtInner2, evtUsedPtFallback);
    } catch (const std::exception& ex) {
        error() << "Exception during track finding: " << ex.what() << endmsg;
        info() << std::string(80, '=') << "\n" << endmsg; // Bottom separator
       return std::make_tuple(std::move(trackCollection), std::move(outputHits), std::move(outputLinks)); // Return empty collections
    }
    
    // Count per-event disambiguation (once per event, not per track)
    if (evtUsedPtInner)    m_statHighestPtInner++;
    if (evtUsedPtInner2)   m_statHighestPtInner2++;
    if (evtUsedPtFallback) m_statHighestPtFallback++;

    info() << "Found " << trackCollection.size() << " tracks" << endmsg;
    
    // Print track details
    for (size_t i = 0; i < trackCollection.size(); ++i) {
        const auto& track = trackCollection[i];
        
        // Get track state at IP if available
        edm4hep::TrackState state;
        bool foundState = false;
        
        for (int j = 0; j < track.trackStates_size(); ++j) {
            if (track.getTrackStates(j).location == edm4hep::TrackState::AtFirstHit) {
                state = track.getTrackStates(j);
                foundState = true;
                break;
            }
        }
        
        if (!foundState && track.trackStates_size() > 0) {
            state = track.getTrackStates(0);
            foundState = true;
        }
        
        if (foundState) {
            // Get magnetic field at track position
            dd4hep::Position fieldPos(0, 0, 0);  // Default center position
            if (track.trackerHits_size() > 0) {
                // Use average position of hits for better field estimate
                double sumX = 0, sumY = 0, sumZ = 0;
                for (int j = 0; j < track.trackerHits_size(); j++) {
                    auto hit = track.getTrackerHits(j);
                    sumX += hit.getPosition()[0] / 10.0;  // mm to cm
                    sumY += hit.getPosition()[1] / 10.0;
                    sumZ += hit.getPosition()[2] / 10.0;
                }
                fieldPos = dd4hep::Position(
                    sumX / track.trackerHits_size(),
                    sumY / track.trackerHits_size(),
                    sumZ / track.trackerHits_size()
                );
            }
            
            // Extract field value in Tesla
            double bField = 0.0;
            try {
                bField = m_field.magneticField(fieldPos).z() / dd4hep::tesla;
            } catch (...) {
                bField = 0.0;  // Handle any field access errors
            }
            
            // Calculate pT from omega and magnetic field
            double omega = state.omega;
            double pT = 0.3 * std::abs(bField) / std::abs(omega) * 0.001; // GeV/c
            
            // Calculate eta from tanLambda
            double tanLambda = state.tanLambda;
            double theta = std::atan2(1.0, tanLambda);
            double eta = -std::log(std::tan(theta/2.0));
            
            // Calculate d0 and z0 in cm for display
            double d0 = state.D0 / 10.0;  // mm to cm
            double z0 = state.Z0 / 10.0;  // mm to cm
            
            // Display track information
            info() << "Track " << i << ":" << endmsg;
            info() << "  pT = " << pT << " GeV/c" << endmsg;
            info() << "  phi = " << state.phi << " rad" << endmsg;
            info() << "  eta = " << eta << endmsg;
            info() << "  d0 = " << d0 << " cm" << endmsg;
            info() << "  z0 = " << z0 << " cm" << endmsg;
            info() << "  chi2/ndof = " << track.getChi2() / track.getNdf() << endmsg;
            info() << "  # hits = " << track.trackerHits_size() << endmsg;
        }
    }
    
    // Bottom separator
    info() << std::string(80, '=') << "\n" << endmsg;
    
    return std::make_tuple(std::move(trackCollection), std::move(outputHits), std::move(outputLinks));
}

// Finalize method
StatusCode DisplacedTracking::finalize() {

    // ================================================================
    //                   RUN STATISTICS SUMMARY
    // ================================================================
    int nEvt    = m_statTotalEvents.load();
    int nTrk    = m_statTracksReconstructed.load();
    int nPos    = m_statPositiveCharge.load();
    int nNeg    = m_statNegativeCharge.load();
    int nInner  = m_statInnerSeedUsed.load();
    int nFall   = m_statFallbackSeedUsed.load();
    int nPtI    = m_statHighestPtInner.load();
    int nPtI2   = m_statHighestPtInner2.load();
    int nPtF    = m_statHighestPtFallback.load();
    int nCombos = m_statTotalTripletCombos.load();
    int nValid  = m_statValidTriplets.load();
    int n4hit   = m_statFourHitTracks.load();
    int n3hit   = m_statThreeHitTracks.load();
    int nSF     = m_statStateAtFirstHit.load();
    int nSL     = m_statStateAtLastHit.load();
    int nSC     = m_statStateAtCalorimeter.load();
    int nSO     = m_statStateAtOther.load();
    int nSV     = m_statStateAtVertex.load();
    int nProp   = m_statInnerPropSuccess.load();
    int nGF     = m_statGenFitSuccess.load();
    int nGFFail = m_statGenFitFailed.load();
    int nGFIll  = m_statGenFitIllCond.load();
    int nGFExc  = m_statGenFitException.load();
    int nGFChi2N = m_statGenFitChi2Count.load();
    double avgChi2ndf = (nGFChi2N > 0) 
                        ? (m_statGenFitChi2Sum.load() / 1000.0) / nGFChi2N 
                        : 0.0;
    int nBadGeom = m_statTripletBadGeom.load();
    int nPTCut   = m_statTripletsCutByPT.load();
    double trkPerEvt = (nEvt>0)   ? double(nTrk)/double(nEvt) : 0.0;
    double seedEff   = (nCombos>0)? 100.0*double(nValid)/double(nCombos) : 0.0;
    double innerFrac = (nTrk>0)   ? 100.0*double(nInner)/double(nTrk) : 0.0;
    double propFrac  = (nTrk>0)   ? 100.0*double(nProp)/double(nTrk) : 0.0;
    double gfFrac    = (nTrk>0)   ? 100.0*double(nGF)/double(nTrk) : 0.0;
    info() << "\n"
        << "╔══════════════════════════════════════════════════════════╗\n"
        << "║          DISPLACED TRACKING  -  RUN SUMMARY             ║\n"
        << "╠══════════════════════════════════════════════════════════╣\n"
        << "║  Events processed             : " << std::setw(7) << nEvt    << "                   ║\n"
        << "║  Tracks reconstructed         : " << std::setw(7) << nTrk    << "                   ║\n"
        << "║  Avg tracks / event           : " << std::setw(7) << std::fixed << std::setprecision(2) << trkPerEvt << "                   ║\n"
        << "╠══════════════════════════════════════════════════════════╣\n"
        << "║  CURVATURE / CHARGE SIGN                                 ║\n"
        << "║    Positive particles (w < 0) : " << std::setw(7) << nPos    << "                   ║\n"
        << "║    Negative particles (w > 0) : " << std::setw(7) << nNeg    << "                   ║\n"
        << "╠══════════════════════════════════════════════════════════╣\n"
        << "║  SEEDING STRATEGY                                        ║\n"
        << "║    Inner-layer seed (0,1,2)   : " << std::setw(7) << nInner  << "  (" << std::setw(5) << std::fixed << std::setprecision(1) << innerFrac << "% of tracks)  ║\n"
        << "║    Fallback Phase 1/2 seed    : " << std::setw(7) << nFall   << "                   ║\n"
        << "╠══════════════════════════════════════════════════════════╣\n"
        << "║  HIGHEST-pT DISAMBIGUATION (# events with >1 candidate) ║\n"
        << "║    Best-pT - inner 0,1,2      : " << std::setw(7) << nPtI   << "                   ║\n"
        << "║    Best-pT - inner+layer3     : " << std::setw(7) << nPtI2  << "                   ║\n"
        << "║    Best-pT - fallback Ph.1/2  : " << std::setw(7) << nPtF   << "                   ║\n"
        << "╠══════════════════════════════════════════════════════════╣\n"
        << "║  TRIPLET SEEDING                                         ║\n"
        << "║    Combinations tried         : " << std::setw(7) << nCombos << "                   ║\n"
        << "║    Valid triplets formed      : " << std::setw(7) << nValid  << "  (" << std::setw(5) << std::fixed << std::setprecision(1) << seedEff  << "% efficiency)  ║\n"
        << "║    Marginal geometry (warned) : " << std::setw(7) << nBadGeom << "                   ║\n"
        << "║    Rejected by MaxSeedPT cut  : " << std::setw(7) << nPTCut  << "                   ║\n"
        << "╠══════════════════════════════════════════════════════════╣\n"
        << "║  HIT MULTIPLICITY                                        ║\n"
        << "║    Tracks with 4 hits         : " << std::setw(7) << n4hit   << "                   ║\n"
        << "║    Tracks with 3 hits         : " << std::setw(7) << n3hit   << "                   ║\n"
        << "╠══════════════════════════════════════════════════════════╣\n"
        << "║  TRACK STATES STORED                                     ║\n"
        << "║    AtFirstHit  (seed state)   : " << std::setw(7) << nSF     << "                   ║\n"
        << "║    AtLastHit   (4-hit circle) : " << std::setw(7) << nSL     << "                   ║\n"
        << "║    AtCalorimeter              : " << std::setw(7) << nSC     << "                   ║\n"
        << "║    AtOther     (GenFit fit)   : " << std::setw(7) << nSO     << "                   ║\n"
        << "║    AtVertex    (inner prop.)  : " << std::setw(7) << nSV     << "                   ║\n"
        << "╠══════════════════════════════════════════════════════════╣\n"
        << "║  DOWNSTREAM PROCESSING                                   ║\n"
        << "║    Inner solenoid prop. OK    : " << std::setw(7) << nProp   << "  (" << std::setw(5) << std::fixed << std::setprecision(1) << propFrac << "% of tracks)  ║\n"
        << "║    GenFit fit succeeded       : " << std::setw(7) << nGF     << "  (" << std::setw(5) << std::fixed << std::setprecision(1) << gfFrac   << "% of tracks)  ║\n"
        << "╠══════════════════════════════════════════════════════════╣\n"
        << "║  GENFIT DETAILS                                          ║\n"
        << "║    Fit succeeded              : " << std::setw(7) << nGF     << "                   ║\n"
        << "║    Failed (isFitted=false)    : " << std::setw(7) << nGFFail << "                   ║\n"
        << "║    Failed (ill-conditioned)   : " << std::setw(7) << nGFIll  << "                   ║\n"
        << "║    Failed (other exception)   : " << std::setw(7) << nGFExc  << "                   ║\n"
        << "║    Fits w/ non-physical chi2  : " << std::setw(7) << m_statGenFitChi2BadCount.load() << "  (chi2/ndf>1000, excl.)  ║\n"
        << "║    Avg chi2/ndf (sane fitted) : " << std::setw(7) << std::fixed << std::setprecision(3) << avgChi2ndf << "                   ║\n"
        << "╚══════════════════════════════════════════════════════════╝\n"
        << endmsg;
    // ================================================================

    // Clean up material manager
    if (m_materialManager) {
        delete m_materialManager;
        m_materialManager = nullptr;
    }
    
    // Clean up GenFit resources if needed
    if (m_useGenFit) {
        // Field manager and material effects are singletons
        // They will be cleaned up automatically at program exit
        
        // Release our adapters
        m_genFitField.reset();
        m_genFitMaterial.reset();
    }
    
    return Algorithm::finalize();
}

// Find surface for a hit
const dd4hep::rec::Surface* DisplacedTracking::findSurface(const edm4hep::TrackerHitPlane& hit) const {
    return findSurfaceByID(hit.getCellID());
}

// Find surface by cell ID
const dd4hep::rec::Surface* DisplacedTracking::findSurfaceByID(uint64_t cellID) const {
    // Find surface for this cell ID using the surface map
    dd4hep::rec::SurfaceMap::const_iterator sI = m_surfaceMap->find(cellID);
    
    if (sI != m_surfaceMap->end()) {
        return dynamic_cast<const dd4hep::rec::Surface*>(sI->second);
    }
    
    return nullptr;
}

// Extract Type (Barrel, Endcap) ID from cell ID
int DisplacedTracking::getTypeID(uint64_t cellID) const {
    // Use BitFieldCoder to extract type ID
    return m_bitFieldCoder->get(cellID, "type");
}

// Extract layer ID from cell ID
int DisplacedTracking::getLayerID(uint64_t cellID) const {
    // Use BitFieldCoder to extract layer ID
    return m_bitFieldCoder->get(cellID, "layer");
}

// Compute composite layer ID that encodes both region and layer number:
//   barrel:    compositeID = layerID
//   +endcap:   compositeID = 1000 + layerID
//   -endcap:   compositeID = -1000 - layerID
int DisplacedTracking::getCompositeID(uint64_t cellID) const {
    int layerID = getLayerID(cellID);
    int typeID  = getTypeID(cellID);
    if      (typeID ==  0) return layerID;
    else if (typeID ==  1) return  1000 + layerID;
    else                   return -1000 - layerID;
}

// Group surfaces by detector layer
std::map<int, std::vector<const dd4hep::rec::Surface*>> DisplacedTracking::getSurfacesByLayer() const {
    std::map<int, std::vector<const dd4hep::rec::Surface*>> result;
    
    // Process all surfaces and group them by layer ID
    for (const auto& surface : m_surfaces) {
        // Get detector element
        dd4hep::DetElement det = surface->detElement();
        if (!det.isValid()) continue;
        
        // Get volume ID
        uint64_t volID = det.volumeID();
        int layerID = 0;
        int type = 0;  // 0=barrel, 1=positive endcap, -1=negative endcap
        

        // Get the layer ID
        layerID = m_bitFieldCoder->get(volID, "layer");
            
        // Try to get the type field
        type = m_bitFieldCoder->get(volID, "type");
            
        // Create a composite ID that distinguishes barrel and endcaps
        // For barrel: layerID
        // For endcaps: 1000 * type + layerID (1000+layerID for positive, -1000+layerID for negative)
        int compositeID;
        if (type == 0) {
                // Barrel
                compositeID = layerID;
        } else if (type == 1) {
                // Endcaps
                compositeID = 1000 * type + layerID;
        } else {
                compositeID = 1000 * type - layerID;
        }
            
        // Add surface to the appropriate layer group
        result[compositeID].push_back(surface);
            
        // Debug output
        if (msgLevel(MSG::DEBUG)) {
                dd4hep::rec::Vector3D pos = surface->origin();
                std::string typeStr;
                if (type == 0) typeStr = "Barrel";
                else if (type == 1) typeStr = "Positive Endcap";
                else if (type == -1) typeStr = "Negative Endcap";
                else typeStr = "Unknown";
                
               /* debug() << "Surface at (" << pos.x() << ", " << pos.y() << ", " 
                       << pos.z() << ") - " << typeStr << " Layer " << layerID 
                       << " (ID=" << compositeID << ")" << endmsg; */
        }
    } 
    
    // Output summary of layer distribution
    if (msgLevel(MSG::INFO)) {
        info() << "Layer distribution:" << endmsg;
        for (const auto& [layer, surfaces] : result) {
            std::string typeStr;
            int actualLayer;
            
            if (layer >= 1000) {
                // Positive endcap
                typeStr = "Positive Endcap";
                actualLayer = layer - 1000;
            } else if (layer <= -1000) {
                // Negative endcap
                typeStr = "Negative Endcap";
                actualLayer = -(layer + 1000);  // Make layer number positive for display
            } else {
                // Barrel
                typeStr = "Barrel";
                actualLayer = layer;
            }
            
            info() << "  " << typeStr << " Layer " << actualLayer 
                  << " (ID=" << layer << "): " << surfaces.size() << " surfaces" << endmsg;
        }
    }
    
    return result;
}


double DisplacedTracking::getRadiationLength(const dd4hep::rec::Vector3D& start, 
                                        const dd4hep::rec::Vector3D& end) const {
    if (!m_materialManager) {
        warning() << "Material manager not initialized for radiation length calculation" << endmsg;
        return 0.0;
    }
    
    // Get materials along the path
    const dd4hep::rec::MaterialVec& materials = 
        m_materialManager->materialsBetween(start, end);
    
    double totalRadiationLength = 0.0;
    
    // Process all materials along the path
    for (const auto& material_pair : materials) {
        // Each entry is a pair of (Material, thickness)
        double thickness = material_pair.second;
        const dd4hep::Material& material = material_pair.first;
        
        // Get radiation length from the Material object
        double radLength = material.radLength(); // Or use material.radiationLength()
        
        if (radLength > 0) {
            totalRadiationLength += thickness / radLength;
        }
    }
    
    return totalRadiationLength;
}

// Get Transverse Momentum
double DisplacedTracking::getPT(const edm4hep::TrackState& state) const {
    // Get parameters
    double omega = state.omega;  // 1/mm

    // Use reference point (first hit position)
    double x = state.referencePoint[0];
    double y = state.referencePoint[1];
    double z = state.referencePoint[2];
    
    // Get field at this position to calculate momentum
    dd4hep::Position fieldPos(x/10.0, y/10.0, z/10.0); // Convert mm to cm for DD4hep
    double bField = m_field.magneticField(fieldPos).z() / dd4hep::tesla;

    // Calculate pT from omega
    double pT = 0.3 * std::abs(bField) / std::abs(omega) * 0.001; // GeV/c

    return pT;
}
// Get Position
Eigen::Vector3d DisplacedTracking::getPosition(const edm4hep::TrackState& state) const {
    // Get parameters
    double d0 = state.D0;      // mm
    double phi = state.phi;    // rad
    double z0 = state.Z0;      // mm
    
    // Position at closest approach to origin
    double x0 = -d0 * std::sin(phi);
    double y0 = d0 * std::cos(phi);
    
    return Eigen::Vector3d(x0, y0, z0);
}
// EDM4HEP Track state
edm4hep::TrackState DisplacedTracking::createTrackState(
    double d0, double phi, double omega, double z0, double tanLambda,
    int location) const {
    
    edm4hep::TrackState state;
    
    // Set parameters
    state.D0 = d0;
    state.phi = phi;
    state.omega = omega;
    state.Z0 = z0;
    state.tanLambda = tanLambda;
    state.location = location;
    
    // Zero the full 6×6 covariance (21 lower-triangle elements).
    // Callers that know their fit uncertainties should set the diagonal
    // explicitly after calling this function.  The defaults below are
    // conservative fall-backs for the muon system (poor resolution).
    for (int i = 0; i < 6; ++i)
        for (int j = 0; j <= i; ++j)
            state.setCovMatrix(0.0,
                static_cast<edm4hep::TrackParams>(i),
                static_cast<edm4hep::TrackParams>(j));

    // Default diagonal: generous seed uncertainties for the outer muon system.
    // Units: D0/Z0 in mm², phi/tanLambda in rad², omega in (1/mm)².
    state.setCovMatrix(1.0,   edm4hep::TrackParams::d0,        edm4hep::TrackParams::d0);
    state.setCovMatrix(0.01,  edm4hep::TrackParams::phi,       edm4hep::TrackParams::phi);
    state.setCovMatrix(1e-10,  edm4hep::TrackParams::omega,     edm4hep::TrackParams::omega);
    state.setCovMatrix(1.0,   edm4hep::TrackParams::z0,        edm4hep::TrackParams::z0);
    state.setCovMatrix(0.01,  edm4hep::TrackParams::tanLambda, edm4hep::TrackParams::tanLambda);

    return state;
}

//------------------------------------------------------------------------------
// Core tracking methods
//------------------------------------------------------------------------------

void DisplacedTracking::findTracks(
    const edm4hep::TrackerHitPlaneCollection* hits,
    edm4hep::TrackCollection& finalTracks,
    edm4hep::TrackerHitPlaneCollection& outputHits,
    const edm4hep::TrackerHitSimTrackerHitLinkCollection& recoSimLinks,
    edm4hep::TrackerHitSimTrackerHitLinkCollection& outputLinks,
    bool& evtUsedPtInner,
    bool& evtUsedPtInner2,
    bool& evtUsedPtFallback) const {

    if (hits->size() < 3) {
        info() << "Not enough hits for tracking. Need at least 3, found " << hits->size() << endmsg;
        return;
    }

    // Build a map from input TrackerHitPlane objectID -> SimTrackerHit
    // using the existing RecoSim link collection from the digitizer
    std::unordered_map<uint32_t, edm4hep::SimTrackerHit> recoToSimMap;
    for (const auto& link : recoSimLinks) {
        // link.getFrom() is a TrackerHit interface — get the underlying object ID
        auto recoHit = link.getFrom();
        recoToSimMap[recoHit.id().index] = link.getTo();
    }

    // Helper lambda: given a newly created output hit whose data was copied from
    // an input hit, propagate the RecoSim link to the output link collection
    auto propagateLink = [&](const edm4hep::MutableTrackerHitPlane& outHit,
                              const edm4hep::TrackerHitPlane& inHit) {
        auto it = recoToSimMap.find(inHit.id().index);
        if (it != recoToSimMap.end()) {
            auto outLink = outputLinks.create();
            outLink.setFrom(outHit);
            outLink.setTo(it->second);
            outLink.setWeight(1.0f);
        }
    };

    // Global hit usage tracking across all track-finding iterations
    std::vector<bool> globalUsedHits(hits->size(), false);
    int trackNumber = 0;
    
    // MAIN LOOP: Continue finding tracks until no more can be formed
    while (true) {
        trackNumber++;
        info() << "\n=== TRACK FINDING ITERATION " << trackNumber << " ===" << endmsg;
        
        // Count remaining unused hits
        int remainingHits = 0;
        std::map<int, int> unusedHitsPerLayer;
        for (size_t i = 0; i < globalUsedHits.size(); ++i) {
            if (!globalUsedHits[i]) {
                remainingHits++;
                int layerID = getLayerID((*hits)[i].getCellID());
                unusedHitsPerLayer[layerID]++;
            }
        }
        
        info() << "Remaining unused hits: " << remainingHits << endmsg;
        for (const auto& [layer, count] : unusedHitsPerLayer) {
            info() << "  Layer " << layer << ": " << count << " unused hits" << endmsg;
        }
        
        if (remainingHits < 3) {
            info() << "Not enough unused hits remaining (" << remainingHits << ") to form another track" << endmsg;
            break;
        }

        // Structure to hold hit information for this iteration (only unused hits)
        struct HitInfo {
            size_t index;
            edm4hep::TrackerHitPlane hit;
            int layerID;
            int typeID;
            int compositeID;

            HitInfo(size_t idx, const edm4hep::TrackerHitPlane& h, int layer, int type, int composite)
                : index(idx), hit(h), layerID(layer), typeID(type), compositeID(composite) {}
        };

        std::vector<HitInfo> allHitInfo;
        std::map<int, std::vector<size_t>> hitIndicesByCompositeLayer; 

        // Reserve space to prevent reallocation during population
        allHitInfo.reserve(remainingHits);

        // Build hit info for ONLY unused hits
        for (size_t i = 0; i < hits->size(); ++i) {
            if (globalUsedHits[i]) continue;  // Skip already used hits
            
            const auto& hit = (*hits)[i];
            const dd4hep::rec::Surface* surface = findSurface(hit);

            if (surface) {
                int layerID     = getLayerID(hit.getCellID());
                int typeID      = getTypeID(hit.getCellID());
                int compositeID = getCompositeID(hit.getCellID());

                allHitInfo.emplace_back(i, hit, layerID, typeID, compositeID);
                hitIndicesByCompositeLayer[compositeID].push_back(allHitInfo.size() - 1);

                debug() << "Available hit " << i << ": CellID=" << hit.getCellID() 
                        << ", Layer=" << layerID << ", Type=" << typeID 
                        << ", Composite=" << compositeID << endmsg;
            } else {
                debug() << "Could not find surface for unused hit with cellID " << hit.getCellID() << endmsg;
            }
        }

        if (allHitInfo.size() < 3) {
            info() << "Not enough unused hits with valid surfaces (" << allHitInfo.size() << ") for iteration " << trackNumber << endmsg;
            break;
        }

        info() << "Available hits for iteration " << trackNumber << ": " << allHitInfo.size() 
               << " in " << hitIndicesByCompositeLayer.size() << " detector different layers/regions" << endmsg;

        // Get composite layer IDs for available hits
        std::vector<int> compositeLayerIDs;
        for (const auto& [compositeID, indices] : hitIndicesByCompositeLayer) {
            compositeLayerIDs.push_back(compositeID);
        }

        // Sort compositeLayerIDs by layer number (nearest to IP first), then by compositeID as tie-breaker.
        auto layerNumber = [](int compositeID) {
            return std::abs(compositeID) % 1000;
        };
        std::sort(compositeLayerIDs.begin(), compositeLayerIDs.end(),
            [&](int a, int b) {
                int la = layerNumber(a);
                int lb = layerNumber(b);
                if (la != lb) return la < lb;            // smaller layer number => nearer to IP
                return a < b;                            // tie-break by compositeID numeric order
            });

        if (compositeLayerIDs.size() < 3) {
            info() << "Not enough detector regions with unused hits for triplet seeding" << endmsg;
            break;
        }

        // Build mapping layerNumber -> vector of HitInfo indices (aggregating across composite types)
        std::map<int, std::vector<size_t>> hitsByLayerNumber;
        for (size_t idx = 0; idx < allHitInfo.size(); ++idx) {
            int lyr = allHitInfo[idx].layerID;
            hitsByLayerNumber[lyr].push_back(idx);
        }

        // Create local used hits array for this iteration (starts all false for available hits)
        std::vector<bool> usedHits(hits->size(), false);
        
        // Copy global usage to local array
        for (size_t i = 0; i < globalUsedHits.size(); ++i) {
            usedHits[i] = globalUsedHits[i];
        }

        struct TrackCandidate {
            edm4hep::Track track;
            double pT;
            std::vector<size_t> hitIndices;
            std::vector<int> usedCompositeIDs;

            TrackCandidate(edm4hep::Track t, double pt, const std::vector<size_t>& indices, 
                          const std::vector<int>& composites)
                : track(t), pT(pt), hitIndices(indices), usedCompositeIDs(composites) {}
        };

        std::vector<TrackCandidate> trackCandidates;
        int tripletCandidates = 0;
        int validTriplets = 0;

        // Helper to extract pT from a track (using getPT method)
        auto extractPTFromTrack = [&](const edm4hep::Track& t) {
            double pT = 0.0;
            for (int l = 0; l < t.trackStates_size(); ++l) {
                auto state = t.getTrackStates(l);
                if (state.location == edm4hep::TrackState::AtFirstHit) {
                    pT = getPT(state);
                    break;
                }
            }
            return pT;
        };

        edm4hep::TrackCollection candidateTracks;  // Temporary collection for this iteration

        bool innerSeedChosen = false;

        // -------------------------------------------------------------------------
        // Try to seed from the INNER layers (0,1,2) first (across all composites)
        // -------------------------------------------------------------------------
        info() << "Phase A: Trying inner-layer seeding (layers 0,1,2) for iteration " << trackNumber << "..." << endmsg;

        // Determine availability in layers 0,1,2 and layer3
        int layer0Count = hitsByLayerNumber.count(0) ? hitsByLayerNumber[0].size() : 0;
        int layer1Count = hitsByLayerNumber.count(1) ? hitsByLayerNumber[1].size() : 0;
        int layer2Count = hitsByLayerNumber.count(2) ? hitsByLayerNumber[2].size() : 0;
        int layer3Count = hitsByLayerNumber.count(3) ? hitsByLayerNumber[3].size() : 0;

        // Case 1: we have hits in all inner three layers -> form all combinations [layer0 x layer1 x layer2],
        // pick the single best triplet (highest pT), mark its hits used and add it as a candidate.
        if (layer0Count > 0 && layer1Count > 0 && layer2Count > 0) {
            info() << "Inner layers 0,1,2 have unused hits. Forming triplets from these layers and selecting best pT." << endmsg;

            std::vector<TrackCandidate> innerCandidates;
            // iterate all combinations using aggregated layer lists (which reference allHitInfo indices)
            for (size_t idx0 : hitsByLayerNumber[0]) {
                for (size_t idx1 : hitsByLayerNumber[1]) {
                    for (size_t idx2 : hitsByLayerNumber[2]) {
                        const HitInfo& h0 = allHitInfo[idx0];
                        const HitInfo& h1 = allHitInfo[idx1];
                        const HitInfo& h2 = allHitInfo[idx2];

                        // check bounds & used flags
                        if (h0.index >= usedHits.size() || h1.index >= usedHits.size() || h2.index >= usedHits.size()) continue;
                        if (usedHits[h0.index] || usedHits[h1.index] || usedHits[h2.index]) continue;

                        tripletCandidates++;

                        size_t prevSize = candidateTracks.size();
                        bool seedValid = false;
                        try {
                            seedValid = createTripletSeed(h0.hit, h1.hit, h2.hit, &candidateTracks, outputHits, propagateLink, usedHits, h0.index, h1.index, h2.index);
                        } catch (const std::exception& ex) {
                            warning() << "Exception in createTripletSeed (inner): " << ex.what() << endmsg;
                            continue;
                        } catch (...) {
                            continue;
                        }

                        if (seedValid && candidateTracks.size() > prevSize) {
                            // newly added track
                            validTriplets++;

                            auto newTrack = candidateTracks[candidateTracks.size() - 1];
                            double pT = extractPTFromTrack(newTrack);
                            std::vector<size_t> hitIndices = {h0.index, h1.index, h2.index};
                            std::vector<int> usedComposites = {h0.compositeID, h1.compositeID, h2.compositeID};
                            innerCandidates.emplace_back(newTrack, pT, hitIndices, usedComposites);

                            // Undo marking used for selection stage
                            usedHits[h0.index] = false;
                            usedHits[h1.index] = false;
                            usedHits[h2.index] = false;
                        }
                    }
                }
            }

            if (!innerCandidates.empty()) {
                // choose best by pT (highest)
                auto bestIt = std::max_element(innerCandidates.begin(), innerCandidates.end(),
                                               [](const TrackCandidate& a, const TrackCandidate& b){ return a.pT < b.pT; });
                if (bestIt != innerCandidates.end()) {
                    // Add the selected best candidate (and mark its hits used)
                    trackCandidates.push_back(*bestIt);
                    for (size_t hi : bestIt->hitIndices) {
                        if (hi < usedHits.size()) usedHits[hi] = true;
                    }
                    innerSeedChosen = true;
                    info() << "Chosen best inner-layer triplet for iteration " << trackNumber << " with pT=" << bestIt->pT << " GeV/c" << endmsg;
                    m_statInnerSeedUsed++;
                    if (innerCandidates.size() > 1) evtUsedPtInner = true;
                }
            } else {
                info() << "No valid inner triplet seeds found in layers 0,1,2 for iteration " << trackNumber << endmsg;
            }
        }
        // Case 2: fewer than 3 hits across layers 0..2 -> allow third hit from layer 3 (if present)
        else {
            int totalInnerHits = layer0Count + layer1Count + layer2Count;
            if (totalInnerHits >= 2 && layer3Count > 0) {
                info() << "Less than full inner set; forming triplets using any two hits from layers 0..2 plus one hit from layer 3 for iteration " << trackNumber << endmsg;

                // build vector of all inner indices
                std::vector<size_t> innerIndices;
                for (int l = 0; l <= 2; ++l) {
                    if (hitsByLayerNumber.count(l)) {
                        for (size_t idx : hitsByLayerNumber[l]) innerIndices.push_back(idx);
                    }
                }

                // layer3 indices
                std::vector<size_t> layer3Indices;
                if (hitsByLayerNumber.count(3)) layer3Indices = hitsByLayerNumber[3];

                // combine: two from innerIndices + one from layer3Indices
                for (size_t a = 0; a < innerIndices.size(); ++a) {
                    for (size_t b = a + 1; b < innerIndices.size(); ++b) {
                        for (size_t idx3 : layer3Indices) {
                            const HitInfo& hA = allHitInfo[innerIndices[a]];
                            const HitInfo& hB = allHitInfo[innerIndices[b]];
                            const HitInfo& hC = allHitInfo[idx3];

                            if (hA.index >= usedHits.size() || hB.index >= usedHits.size() || hC.index >= usedHits.size()) continue;
                            if (usedHits[hA.index] || usedHits[hB.index] || usedHits[hC.index]) continue;

                            tripletCandidates++;
                            size_t prevSize = candidateTracks.size();
                            bool seedValid = false;
                            try {
                                seedValid = createTripletSeed(hA.hit, hB.hit, hC.hit, &candidateTracks, outputHits, propagateLink, usedHits, hA.index, hB.index, hC.index);
                            } catch (const std::exception& ex) {
                                warning() << "Exception in createTripletSeed (inner+layer3): " << ex.what() << endmsg;
                                continue;
                            } catch (...) {
                                continue;
                            }

                            if (seedValid && candidateTracks.size() > prevSize) {
                                validTriplets++;
                                auto newTrack = candidateTracks[candidateTracks.size() - 1];
                                double pT = extractPTFromTrack(newTrack);
                                std::vector<size_t> hitIndices = {hA.index, hB.index, hC.index};
                                std::vector<int> usedComposites = {hA.compositeID, hB.compositeID, hC.compositeID};
                                trackCandidates.emplace_back(newTrack, pT, hitIndices, usedComposites);
                            }
                        }
                    }
                }

                if (!trackCandidates.empty()) {
                    // If multiple choices, keep the one with highest pT
                    std::sort(trackCandidates.begin(), trackCandidates.end(),
                              [](const TrackCandidate& a, const TrackCandidate& b){ return a.pT > b.pT; });
                    // Mark used hits for the chosen best and discard the rest
                    const auto& best = trackCandidates.front();
                    for (size_t hi : best.hitIndices) if (hi < usedHits.size()) usedHits[hi] = true;
                    // keep only best
                    std::vector<TrackCandidate> kept = {best};
                    trackCandidates.swap(kept);
                    innerSeedChosen = true;
                    info() << "Selected best combined inner+layer3 triplet for iteration " << trackNumber << " with pT=" << best.pT << " GeV/c" << endmsg;
                    m_statInnerSeedUsed++;
                    if (kept.size() > 1) evtUsedPtInner2 = true;
                } else {
                    info() << "No valid triplets formed from inner+layer3 strategy for iteration " << trackNumber << endmsg;
                }
            } else {
                info() << "Not enough inner hits (or no layer3) to form inner-preferred triplets for iteration " << trackNumber << "; will proceed to Phase 1/2 fallback." << endmsg;
            }
        }

        // If we successfully chose an inner seed, skip broad Phase 1/2 search for this iteration
        if (!innerSeedChosen) {
            // ---------------------------------------------------------------------
            // PHASE 1: Consecutive layer triplets (but now compositeLayerIDs are sorted nearest-to-IP)
            // ---------------------------------------------------------------------
            bool foundGoodConsecutiveTriplets = false;
            bool foundInnerLayerTriplets = false;  // Track if we found triplets in inner layers (0,1,2)
            int consecutiveTriplets = 0;
            int maxEarlyLayersToTry = std::min(size_t(5), compositeLayerIDs.size());

            info() << "Phase 1: Testing consecutive layer triplets for iteration " << trackNumber << "..." << endmsg;

            for (size_t i = 0; i + 2 < compositeLayerIDs.size() && i < (size_t)std::max(0, maxEarlyLayersToTry - 2); ++i) {
                int composite1 = compositeLayerIDs[i];
                int composite2 = compositeLayerIDs[i+1];
                int composite3 = compositeLayerIDs[i+2];

                debug() << "Testing consecutive triplet: regions " 
                        << composite1 << " -> " << composite2 << " -> " << composite3 << endmsg;

                // Check if this triplet uses the first inner layers (0,1,2)
                bool isInnerLayerTriplet = false;
                int layer1 = std::abs(composite1) % 1000;
                int layer2 = std::abs(composite2) % 1000;
                int layer3 = std::abs(composite3) % 1000;

                if ((layer1 == 0 && layer2 == 1 && layer3 == 2) ||
                    (layer1 <= 2 && layer2 <= 2 && layer3 <= 2)) {
                    isInnerLayerTriplet = true;
                    debug() << "Testing inner layer triplet: layers " << layer1 << "," << layer2 << "," << layer3 << endmsg;
                }

                int tripletsBefore = validTriplets;

                for (size_t idx1 : hitIndicesByCompositeLayer[composite1]) {
                    for (size_t idx2 : hitIndicesByCompositeLayer[composite2]) {
                        for (size_t idx3 : hitIndicesByCompositeLayer[composite3]) {

                            if (idx1 >= allHitInfo.size() || idx2 >= allHitInfo.size() || idx3 >= allHitInfo.size()) {
                                error() << "HitInfo index out of bounds: " << idx1 << ", " << idx2 << ", " << idx3 
                                        << " (allHitInfo.size: " << allHitInfo.size() << ")" << endmsg;
                                continue;
                            }

                            const HitInfo& hitInfo1 = allHitInfo[idx1];
                            const HitInfo& hitInfo2 = allHitInfo[idx2];
                            const HitInfo& hitInfo3 = allHitInfo[idx3];

                            if (hitInfo1.index >= usedHits.size() || 
                                hitInfo2.index >= usedHits.size() || 
                                hitInfo3.index >= usedHits.size()) {
                                error() << "Hit index out of bounds: " << hitInfo1.index << ", " 
                                        << hitInfo2.index << ", " << hitInfo3.index 
                                        << " (usedHits.size: " << usedHits.size() << ")" << endmsg;
                                continue;
                            }

                            if (usedHits[hitInfo1.index] || 
                                usedHits[hitInfo2.index] || 
                                usedHits[hitInfo3.index]) {
                                continue;
                            }

                            tripletCandidates++;

                            size_t prevSize = candidateTracks.size();

                            bool seedValid = false;
                            try {
                                seedValid = createTripletSeed(
                                    hitInfo1.hit, hitInfo2.hit, hitInfo3.hit, 
                                    &candidateTracks, outputHits, propagateLink, usedHits,
                                    hitInfo1.index, hitInfo2.index, hitInfo3.index);
                            } catch (const std::exception& ex) {
                                warning() << "Exception in createTripletSeed: " << ex.what() << endmsg;
                                continue;
                            }

                            if (seedValid && candidateTracks.size() > prevSize) {
                                validTriplets++;
                                consecutiveTriplets++;

                                if (!usedHits[hitInfo1.index] || !usedHits[hitInfo2.index] || !usedHits[hitInfo3.index]) {
                                    warning() << "createTripletSeed didn't mark hits as used - fixing" << endmsg;
                                    usedHits[hitInfo1.index] = true;
                                    usedHits[hitInfo2.index] = true;
                                    usedHits[hitInfo3.index] = true;
                                }

                                auto newTrack = candidateTracks[candidateTracks.size() - 1];
                                double pT = extractPTFromTrack(newTrack);

                                std::vector<size_t> hitIndices = {hitInfo1.index, hitInfo2.index, hitInfo3.index};
                                std::vector<int> usedComposites = {composite1, composite2, composite3};

                                trackCandidates.emplace_back(newTrack, pT, hitIndices, usedComposites);

                                if (msgLevel(MSG::INFO)) {
                                    edm4hep::TrackState trackState;
                                    bool foundTrackState = false;
                                    for (int l = 0; l < newTrack.trackStates_size(); ++l) {
                                        auto state = newTrack.getTrackStates(l);
                                        if (state.location == edm4hep::TrackState::AtFirstHit) {
                                            trackState = state;
                                            foundTrackState = true;
                                            break;
                                        }
                                    }

                                    if (foundTrackState) {
                                        double tanLambda = trackState.tanLambda;
                                        double theta = std::atan2(1.0, tanLambda);
                                        double eta = -std::log(std::tan(theta/2.0));

                                        info() << "✓ Valid triplet #" << validTriplets << " found for iteration " << trackNumber << ":" << endmsg;
                                        info() << "  Layers: " << layer1 << " -> " << layer2 << " -> " << layer3 
                                               << " (composites: " << composite1 << "," << composite2 << "," << composite3 << ")" << endmsg;
                                        info() << "  Track parameters: pT=" << pT << " GeV/c, eta=" << eta 
                                               << ", phi=" << trackState.phi << " rad" << endmsg;
                                        info() << "  Impact parameters: d0=" << trackState.D0/10.0 << " cm, z0=" << trackState.Z0/10.0 << " cm" << endmsg;
                                    }
                                }
                            }
                        }
                    }
                }

                if (validTriplets > tripletsBefore) {
                    foundGoodConsecutiveTriplets = true;
                    if (isInnerLayerTriplet) {
                        foundInnerLayerTriplets = true;
                        info() << "Found " << (validTriplets - tripletsBefore) 
                               << " good triplets in inner layers (0,1,2) for iteration " << trackNumber << " - stopping consecutive layer tests" << endmsg;
                        // We break the top loop because we prefer inner triplets if found here
                        break;
                    }
                }
            }

            // PHASE 2: All combinations if needed - SKIP if we found inner layer triplets above
            bool tryAllCombinations = !foundInnerLayerTriplets && (!foundGoodConsecutiveTriplets || compositeLayerIDs.size() <= 3);

            if (tryAllCombinations) {
                info() << "Phase 2: Trying all layer combinations for iteration " << trackNumber << " (no inner layer triplets found)..." << endmsg;

                for (size_t i = 0; i < compositeLayerIDs.size() - 2; ++i) {
                    for (size_t j = i + 1; j < compositeLayerIDs.size() - 1; ++j) {
                        for (size_t k = j + 1; k < compositeLayerIDs.size(); ++k) {

                            int composite1 = compositeLayerIDs[i];
                            int composite2 = compositeLayerIDs[j];  
                            int composite3 = compositeLayerIDs[k];

                            // Skip if already tested in Phase 1
                            if (foundGoodConsecutiveTriplets && 
                                i < maxEarlyLayersToTry - 2 && j == i + 1 && k == i + 2) {
                                continue;
                            }

                            for (size_t idx1 : hitIndicesByCompositeLayer[composite1]) {
                                for (size_t idx2 : hitIndicesByCompositeLayer[composite2]) {
                                    for (size_t idx3 : hitIndicesByCompositeLayer[composite3]) {

                                        if (idx1 >= allHitInfo.size() || idx2 >= allHitInfo.size() || idx3 >= allHitInfo.size()) {
                                            continue;
                                        }

                                        const HitInfo& hitInfo1 = allHitInfo[idx1];
                                        const HitInfo& hitInfo2 = allHitInfo[idx2];
                                        const HitInfo& hitInfo3 = allHitInfo[idx3];

                                        if (hitInfo1.index >= usedHits.size() || 
                                            hitInfo2.index >= usedHits.size() || 
                                            hitInfo3.index >= usedHits.size()) {
                                            continue;
                                        }

                                        if (usedHits[hitInfo1.index] || 
                                            usedHits[hitInfo2.index] || 
                                            usedHits[hitInfo3.index]) {
                                            continue;
                                        }

                                        tripletCandidates++;

                                        size_t prevSize = candidateTracks.size();

                                        try {
                                            bool seedValid = createTripletSeed(
                                                hitInfo1.hit, hitInfo2.hit, hitInfo3.hit, 
                                                &candidateTracks, outputHits, propagateLink, usedHits,
                                                hitInfo1.index, hitInfo2.index, hitInfo3.index);

                                            if (seedValid && candidateTracks.size() > prevSize) {
                                                validTriplets++;

                                                if (!usedHits[hitInfo1.index] || !usedHits[hitInfo2.index] || !usedHits[hitInfo3.index]) {
                                                    warning() << "createTripletSeed didn't mark hits as used - fixing" << endmsg;
                                                    usedHits[hitInfo1.index] = true;
                                                    usedHits[hitInfo2.index] = true;
                                                    usedHits[hitInfo3.index] = true;
                                                }

                                                auto newTrack = candidateTracks[candidateTracks.size() - 1];
                                                double pT = extractPTFromTrack(newTrack);

                                                std::vector<size_t> hitIndices = {hitInfo1.index, hitInfo2.index, hitInfo3.index};
                                                std::vector<int> usedComposites = {composite1, composite2, composite3};

                                                trackCandidates.emplace_back(newTrack, pT, hitIndices, usedComposites);

                                                if (msgLevel(MSG::INFO)) {
                                                    int layer1 = std::abs(composite1) % 1000;
                                                    int layer2 = std::abs(composite2) % 1000;
                                                    int layer3 = std::abs(composite3) % 1000;

                                                    info() << "✓ Valid triplet #" << validTriplets << " found (Phase 2) for iteration " << trackNumber << ":" << endmsg;
                                                    info() << "  Layers: " << layer1 << " -> " << layer2 << " -> " << layer3 
                                                           << " (composites: " << composite1 << "," << composite2 << "," << composite3 << ")" << endmsg;
                                                    info() << "  Track pT: " << pT << " GeV/c" << endmsg;
                                                }
                                            }
                                        } catch (...) {
                                            continue;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            } else if (foundInnerLayerTriplets) {
                info() << "Phase 2: SKIPPED - Found good triplets in inner layers (0,1,2) for iteration " << trackNumber << endmsg;
            }
        } // end fallback seeding

        info() << "Seed generation complete for iteration " << trackNumber << ": " << trackCandidates.size()
            << " triplet candidates from " << tripletCandidates
            << " combinations (" << validTriplets << " valid)" << endmsg;

        m_statTotalTripletCombos += tripletCandidates;
        m_statValidTriplets      += validTriplets;
        if (!innerSeedChosen && !trackCandidates.empty()) {
            m_statFallbackSeedUsed++;
            if (trackCandidates.size() > 1) evtUsedPtFallback = true;
        }

        // If no track candidates found in this iteration, break
        if (trackCandidates.empty()) {
            info() << "No track candidates found in iteration " << trackNumber << " - stopping track search" << endmsg;
            break;
        }

        // ── MaxSeedPT threshold filter ──────────────────────────────────────
        // Reject triplets whose pT exceeds m_maxSeedPT (e.g. numerically
        // pathological high-curvature seeds).  If ALL candidates are above
        // the threshold we keep the single lowest-pT one as a fallback so
        // we never silently drop a track-finding iteration.
        if (m_maxSeedPT < std::numeric_limits<double>::max()) {
            std::vector<TrackCandidate> below, above;
            for (auto& tc : trackCandidates) {
                if (tc.pT <= m_maxSeedPT) below.push_back(tc);
                else                       above.push_back(tc);
            }

            if (!below.empty()) {
                // Normal case: some candidates pass the cut
                int nCut = static_cast<int>(above.size());
                if (nCut > 0) {
                    info() << "MaxSeedPT filter: rejected " << nCut
                           << " triplet(s) with pT > " << m_maxSeedPT << " GeV/c for iteration " << trackNumber << endmsg;
                    m_statTripletsCutByPT += nCut;
                }
                trackCandidates.swap(below);
            } else {
                // All candidates exceeded threshold — fallback to lowest-pT one
                auto lowestIt = std::min_element(above.begin(), above.end(),
                    [](const TrackCandidate& a, const TrackCandidate& b){ return a.pT < b.pT; });
                warning() << "MaxSeedPT filter: all " << above.size()
                          << " triplet(s) exceed threshold (" << m_maxSeedPT
                          << " GeV/c) for iteration " << trackNumber
                          << ". Keeping lowest-pT fallback with pT=" << lowestIt->pT << " GeV/c" << endmsg;
                int nCut = static_cast<int>(above.size()) - 1;
                if (nCut > 0) m_statTripletsCutByPT += nCut;
                trackCandidates = { *lowestIt };
            }
        }
        // ────────────────────────────────────────────────────────────────────

        // Sort candidates by decreasing pT for better track selection
        std::sort(trackCandidates.begin(), trackCandidates.end(),
                  [](const TrackCandidate& a, const TrackCandidate& b) {
                      return a.pT > b.pT;
                  });

        // Process the best candidate from this iteration
        const auto& candidate = trackCandidates.front();
        
        // Mark the hits used by this track in the global array
        for (size_t idx : candidate.hitIndices) {
            globalUsedHits[idx] = true;
        }

        // All candidates should have their hits marked as used from seed generation
        // Just verify the indices are valid
        bool validIndices = true;
        for (size_t idx : candidate.hitIndices) {
            if (idx >= usedHits.size()) {
                error() << "Final processing: hit index " << idx 
                        << " out of bounds (usedHits size: " << usedHits.size() << ")" << endmsg;
                validIndices = false;
                break;
            }
        }

        if (!validIndices) {
            continue;
        }

        info() << "Processing track candidate " << trackNumber << " with pT=" << candidate.pT 
               << " GeV/c using hits: [" << candidate.hitIndices[0] 
               << "," << candidate.hitIndices[1] << "," << candidate.hitIndices[2] << "]" << endmsg;

        // Get the seed track and its state
        const auto& seedTrack = candidate.track;

        // Find a suitable seed state
        // Priority 1: 4-hit circle fit (AtLastHit) — most constrained analytical estimate
        // Priority 2: 3-hit triplet seed (AtFirstHit)
        edm4hep::TrackState seedState;
        bool foundState = false;

        for (int j = 0; j < seedTrack.trackStates_size(); ++j) {
            if (seedTrack.getTrackStates(j).location == edm4hep::TrackState::AtLastHit) {
                seedState = seedTrack.getTrackStates(j);
                foundState = true;
                debug() << "Using AtLastHit (4-hit circle fit) as GenFit seed for track " << trackNumber << endmsg;
                break;
            }
        }
        if (!foundState) {
            for (int j = 0; j < seedTrack.trackStates_size(); ++j) {
                if (seedTrack.getTrackStates(j).location == edm4hep::TrackState::AtFirstHit) {
                    seedState = seedTrack.getTrackStates(j);
                    foundState = true;
                    debug() << "Using AtFirstHit (3-hit seed) as GenFit seed for track " << trackNumber << endmsg;
                    break;
                }
            }
        }
        if (!foundState && seedTrack.trackStates_size() > 0) {
            seedState = seedTrack.getTrackStates(0);
            foundState = true;
        }

        if (!foundState) {
            warning() << "No seed state found for track " << trackNumber << ", skipping" << endmsg;
            continue;
        }

        debug() << "Found seed state for track candidate " << trackNumber << endmsg;

        // Build trackHits vector directly using hit indices (avoid cellID duplicate issues)
        std::vector<edm4hep::TrackerHitPlane> trackHits;
        trackHits.reserve(3); // Start with 3 hits from triplet

        // Use direct index access instead of cellID lookup to handle duplicate cellIDs
        for (size_t idx : candidate.hitIndices) {
            if (idx < hits->size()) {
                trackHits.push_back((*hits)[idx]); // Direct access by index
            } else {
                warning() << "Hit index " << idx << " out of bounds - skipping hit" << endmsg;
            }
        }

        if (trackHits.size() != candidate.hitIndices.size()) {
            warning() << "Could not retrieve all hits for track candidate " << trackNumber << " - expected " 
                      << candidate.hitIndices.size() << ", got " << trackHits.size() << endmsg;
            continue;
        }

        debug() << "Successfully built trackHits vector with " << trackHits.size() << " hits for track " << trackNumber << endmsg;

        // Extend the track by searching for extra hits in subsequent layers.
        // The seed direction is derived from the 3-hit seed (hit2 - hit0) in 3D.
        // The cone half-angle (35°) gates candidates; the loop continues until
        // no further compatible hit is found.
        bool foundExtraHit = false;
        if (trackHits.size() >= 3) {
            try {
                // Compute seed direction from first and last seed hit (mm -> cm, then normalize)
                const auto& p0 = trackHits.front().getPosition();
                const auto& p2 = trackHits.back().getPosition();
                Eigen::Vector3d seedDir(
                    (p2[0] - p0[0]) / 10.0,
                    (p2[1] - p0[1]) / 10.0,
                    (p2[2] - p0[2]) / 10.0);
                if (seedDir.norm() < 1e-6) {
                    warning() << "Degenerate seed direction for track " << trackNumber << " — skipping extension" << endmsg;
                } else {
                    seedDir.normalize();

                    // Count available unused hits per composite layer for diagnostics
                    int unusedHitCount = 0;
                    std::map<int, int> hitsPerComposite;
                    for (size_t i = 0; i < globalUsedHits.size(); ++i) {
                        if (!globalUsedHits[i]) {
                            unusedHitCount++;
                            hitsPerComposite[getCompositeID((*hits)[i].getCellID())]++;
                        }
                    }
                    info() << "Searching for extra hits for track " << trackNumber
                           << " - last compositeID=" << getCompositeID(trackHits.back().getCellID())
                           << " | unused hits: " << unusedHitCount << endmsg;
                    for (const auto& [composite, count] : hitsPerComposite) {
                        info() << "  CompositeID " << composite << ": " << count << " unused hits" << endmsg;
                    }

                    size_t hitsBefore = trackHits.size();
                    foundExtraHit = findCompatibleExtraHit(trackHits, hits, globalUsedHits, seedDir);

                    if (foundExtraHit) {
                        info() << "Extended track " << trackNumber
                               << " from " << hitsBefore << " to " << trackHits.size() << " hits" << endmsg;
                    }
                }
            } catch (const std::exception& ex) {
                warning() << "Exception in findCompatibleExtraHit for track " << trackNumber << ": " << ex.what() << endmsg;
            } catch (...) {
                warning() << "Unknown exception in findCompatibleExtraHit for track " << trackNumber << endmsg;
            }
        }

        if (!foundExtraHit) {
            debug() << "Keeping track " << trackNumber << " with " << trackHits.size() << " hits (no extra hits found)" << endmsg;
        }

        // Create a new track
        auto finalTrack = finalTracks.create();

        debug() << "Created new final track object for track " << trackNumber << endmsg;

        // ── Always preserve the analytical 3-hit seed state at AtFirstHit ──────────────
        // This is the baseline from circle fitting and must always be present for
        // comparison against the GenFit result and the 4-hit circle fit.
        // We add it now BEFORE the 4-hit extension and GenFit fitting so it is never lost.
        {
            edm4hep::TrackState seedAtFirst = seedState;
            seedAtFirst.location = edm4hep::TrackState::AtFirstHit;
            finalTrack.addToTrackStates(seedAtFirst);
            debug() << "Saved analytical seed state at AtFirstHit:"
                    << "  d0=" << seedAtFirst.D0 << " mm"
                    << "  phi=" << seedAtFirst.phi << " rad"
                    << "  omega=" << seedAtFirst.omega << " 1/mm"
                    << "  z0=" << seedAtFirst.Z0 << " mm"
                    << "  tanL=" << seedAtFirst.tanLambda << endmsg;
        }

        // Test 4-hit circle fitting if we have at least 4 hits.
        // Always uses the first 4 hits for the analytical circle fit — a quality
        // gate on the seed extension. If the 4th hit produces a poor chi2 it is
        // removed and any further hits that were appended by the extension loop
        // are also discarded (they were built on top of a bad 4th hit).
        if (trackHits.size() >= 4) {
            // Slice first 4 hits for the circle fitter
            std::vector<edm4hep::TrackerHitPlane> first4(trackHits.begin(), trackHits.begin() + 4);
            double x0, y0, radius, chi2;
            Eigen::Matrix3d fitCov3x3 = Eigen::Matrix3d::Zero();
            bool fit4Hits = false;
            try {
                fit4Hits = fitCircleToFourHits(first4, x0, y0, radius, chi2, fitCov3x3);
            } catch (const std::exception& ex) {
                warning() << "Exception in fitCircleToFourHits for track " << trackNumber << ": " << ex.what() << endmsg;
            } catch (...) {
                warning() << "Unknown exception in fitCircleToFourHits for track " << trackNumber << endmsg;
            }

            if (fit4Hits) {
                info() << "4-hit circle fit successful for track " << trackNumber
                       << " (total hits: " << trackHits.size() << ")" << endmsg;
                info() << "  Center: (" << x0 << ", " << y0 << ") cm" << endmsg;
                info() << "  Radius: " << radius << " cm" << endmsg;
                info() << "  Chi2/DOF: " << chi2 << endmsg;

                if (chi2 > 50.0) {
                    warning() << "Poor circle fit (chi2=" << chi2 << ") for track " << trackNumber
                              << " — removing 4th hit and all subsequent extension hits" << endmsg;
                    // Truncate back to 3 hits: the 4th hit (and anything built on top of it) is rejected
                    trackHits.resize(3);
                } else {
                    // Query actual B field at the average position of the 4 hits
                    double sumX = 0, sumY = 0, sumZ = 0;
                    for (const auto& h : trackHits) {
                        const auto& pos = h.getPosition();
                        sumX += pos[0];
                        sumY += pos[1];
                        sumZ += pos[2];
                    }
                    dd4hep::Position avgPos(sumX / trackHits.size(),
                                           sumY / trackHits.size(),
                                           sumZ / trackHits.size());
                    double bField4 = m_field.magneticField(avgPos).z() / dd4hep::tesla;
                    double pT = 0.3 * std::abs(bField4) * radius / 100.0; // GeV/c
                    double d0 = std::sqrt(x0*x0 + y0*y0) - radius;
                    double phi = std::atan2(y0, x0) - M_PI/2;
                    double omega = 1.0 / (radius * 10.0); // Assuming positive charge

                    edm4hep::TrackState circleFitState = createTrackState(
                        d0*10.0, phi, omega, seedState.Z0, seedState.tanLambda, 
                        edm4hep::TrackState::AtLastHit);

                    // Set reference point to position of last hit
                    const auto& firstHitPos = trackHits.front().getPosition();
                    circleFitState.referencePoint = edm4hep::Vector3f(
                        firstHitPos[0], firstHitPos[1], firstHitPos[2]);

                // ── Propagate fit covariance [x0,y0,R] → track parameters ──
                // All EDM4hep units: D0/Z0 in mm², phi/omega in rad²/(1/mm)².
                {
                    double cx = x0;
                    double cy = y0;

                    double r2 = cx*cx + cy*cy;
                    double dist0 = std::sqrt(r2);

                    // ---- Gradients ----------------------------------------------------

                    double dDist_dx0 = (dist0 > 1e-9) ? cx / dist0 : 0.0;
                    double dDist_dy0 = (dist0 > 1e-9) ? cy / dist0 : 0.0;

                    // d0 = sqrt(x0² + y0²) − R
                    Eigen::Vector3d grad_d0(dDist_dx0, dDist_dy0, -1.0);

                    // phi = atan2(y0,x0) − π/2
                    Eigen::Vector3d grad_phi(
                        (r2 > 1e-12) ? -cy / r2 : 0.0,
                        (r2 > 1e-12) ?  cx / r2 : 0.0,
                        0.0
                    );

                    // omega = charge/(R*10)   (R in cm → omega in 1/mm)
                    double chargeSign = (omega > 0) ? 1.0 : -1.0;
                    double dOmega_dR = -chargeSign / (radius * radius * 10.0);

                    Eigen::Vector3d grad_omega(0.0, 0.0, dOmega_dR);

                    // ---- Covariance propagation --------------------------------------

                    double var_d0_cm2    = grad_d0.dot(fitCov3x3 * grad_d0);
                    double var_phi_rad2  = grad_phi.dot(fitCov3x3 * grad_phi);
                    double var_omega_mm2 = grad_omega.dot(fitCov3x3 * grad_omega);

                    // Convert d0 variance from cm² → mm²
                    double var_d0_mm2 = var_d0_cm2 * 100.0;

                    // ---- Inherit longitudinal parameters from seed -------------------

                    double var_z0 = seedState.getCovMatrix(
                        edm4hep::TrackParams::z0,
                        edm4hep::TrackParams::z0);

                    double var_tanL = seedState.getCovMatrix(
                        edm4hep::TrackParams::tanLambda,
                        edm4hep::TrackParams::tanLambda);

                    if (var_z0   < 1e-30) var_z0   = 1.0;
                    if (var_tanL < 1e-30) var_tanL = 0.01;

                    // ---- Apply covariance floors -------------------------------------

                    var_d0_mm2   = std::max(var_d0_mm2,   0.01);   // 0.1 mm resolution floor
                    var_phi_rad2 = std::max(var_phi_rad2, 0.01);

                    double omega_rel_floor = std::pow(0.1 * std::abs(omega), 2);

                    var_omega_mm2 = std::max({var_omega_mm2, omega_rel_floor});

                    // ---- Reset covariance matrix -------------------------------------

                    for (int ii = 0; ii < 6; ++ii)
                        for (int jj = 0; jj <= ii; ++jj)
                            circleFitState.setCovMatrix(
                                0.0,
                                static_cast<edm4hep::TrackParams>(ii),
                                static_cast<edm4hep::TrackParams>(jj)
                            );

                    // ---- Fill diagonal terms -----------------------------------------

                    circleFitState.setCovMatrix(
                        static_cast<float>(var_d0_mm2),
                        edm4hep::TrackParams::d0,
                        edm4hep::TrackParams::d0);

                    circleFitState.setCovMatrix(
                        static_cast<float>(var_phi_rad2),
                        edm4hep::TrackParams::phi,
                        edm4hep::TrackParams::phi);

                    circleFitState.setCovMatrix(
                        static_cast<float>(var_omega_mm2),
                        edm4hep::TrackParams::omega,
                        edm4hep::TrackParams::omega);

                    circleFitState.setCovMatrix(
                        static_cast<float>(var_z0),
                        edm4hep::TrackParams::z0,
                        edm4hep::TrackParams::z0);

                    circleFitState.setCovMatrix(
                        static_cast<float>(var_tanL),
                        edm4hep::TrackParams::tanLambda,
                        edm4hep::TrackParams::tanLambda);

                    // ---- Debug print --------------------------------------------------

                    info() << "AtLastHit covariance from 4-hit circle fit:"
                        << "  σ_d0="   << std::sqrt(var_d0_mm2)    << " mm"
                        << "  σ_phi="  << std::sqrt(var_phi_rad2)  << " rad"
                        << "  σ_ω="    << std::sqrt(var_omega_mm2) << " 1/mm"
                        << "  σ_z0="   << std::sqrt(var_z0)        << " mm"
                        << "  σ_tanL=" << std::sqrt(var_tanL)
                        << endmsg;
                }
                    finalTrack.addToTrackStates(circleFitState);
                    finalTrack.setChi2(chi2);

                    // Promote to GenFit seed: 4-hit fit is a better starting point
                    seedState = circleFitState;
                    debug() << "Promoted AtLastHit (4-hit circle fit) to GenFit seed state" << endmsg;

                    info() << "Added 4-hit circle fit track state for track " << trackNumber << ": pT=" << pT 
                           << " GeV/c, d0=" << d0*10.0 << " mm" << endmsg;
                }
            } else {
                warning() << "4-hit circle fit failed for track " << trackNumber << endmsg;
            }
        }

        // Fit the track with GenFit (if enabled)
        bool fitSuccess = false;
        if (m_useGenFit) {
            try {
                fitSuccess = fitTrackWithGenFit(trackHits, seedState, finalTrack, outputHits, propagateLink, hits);
            } catch (const std::exception& ex) {
                error() << "Exception in fitTrackWithGenFit for track " << trackNumber << ": " << ex.what() << endmsg;
            } catch (...) {
                error() << "Unknown exception in fitTrackWithGenFit for track " << trackNumber << endmsg;
            }
        }
        if (fitSuccess) {
            m_statGenFitSuccess++;
            info() << "Successfully fitted track " << trackNumber << " with GenFit: " 
                << finalTrack.trackerHits_size() << " hits, chi2/ndf = " 
                << (finalTrack.getNdf() > 0 ? finalTrack.getChi2() / finalTrack.getNdf() : -1.0) << endmsg;
        } else {
            if (m_useGenFit) {
                warning() << "GenFit track fitting failed for track " << trackNumber << ", using seed parameters" << endmsg;
            } else {
                debug() << "Using analytical seed parameters for track " << trackNumber << " (GenFit disabled)" << endmsg;
            }

            finalTrack.setNdf(std::max(1, static_cast<int>(trackHits.size() * 2 - 5)));

            // Copy seed track states — but skip AtFirstHit since we already saved it above
            // to avoid duplicating the seed state.
            for (int j = 0; j < seedTrack.trackStates_size(); ++j) {
                auto st = seedTrack.getTrackStates(j);
                if (st.location != edm4hep::TrackState::AtFirstHit) {
                    finalTrack.addToTrackStates(st);
                }
            }

            // Copy hits into output collection so PODIO can resolve the relation
            for (const auto& hit : trackHits) {
                auto outHit = outputHits.create();
                outHit.setCellID(hit.getCellID());
                outHit.setTime(hit.getTime());
                outHit.setEDep(hit.getEDep());
                outHit.setEDepError(hit.getEDepError());
                outHit.setPosition(hit.getPosition());
                outHit.setCovMatrix(hit.getCovMatrix());
                outHit.setDu(hit.getDu());
                outHit.setDv(hit.getDv());
                propagateLink(outHit, hit);
                finalTrack.addToTrackerHits(outHit);
            }
        }

        // ======= Reconstruct Inner Track Segment =========
        // Controlled by the steering property DoInnerPropagation (default: true).
        // The result is saved as TrackState::AtVertex so it is distinct from all other states.
        if (m_doInnerPropagation) {
        edm4hep::TrackState bestState;
        bool foundStateForProp = false;
        // ===== EXTRACT BEST STATE FOR INNER PROPAGATION =====
        // Priority 1: GenFit-fitted state (AtOther) — most accurate
        // Priority 2: 4-hit circle fit (AtLastHit)
        // Priority 3: 3-hit seed (AtFirstHit)
        // The reference point in all these states is the first/last hit position
        // in the outer muon system, which is where RK4 propagation starts.
        for (int j = 0; j < finalTrack.trackStates_size(); ++j) {
            auto st = finalTrack.getTrackStates(j);
            if (st.location == edm4hep::TrackState::AtOther) {
                bestState = st;
                foundStateForProp = true;
                info() << "Using AtOther state (GenFit fitted) for inner propagation" << endmsg;
                break;
            }
        }
        if (!foundStateForProp) {
            for (int j = 0; j < finalTrack.trackStates_size(); ++j) {
                auto st = finalTrack.getTrackStates(j);
                if (st.location == edm4hep::TrackState::AtLastHit) {
                    bestState = st;
                    foundStateForProp = true;
                    info() << "Using AtLastHit state (4-hit circle) for inner propagation" << endmsg;
                    break;
                }
            }
        }
        if (!foundStateForProp) {
            for (int j = 0; j < finalTrack.trackStates_size(); ++j) {
                auto st = finalTrack.getTrackStates(j);
                if (st.location == edm4hep::TrackState::AtFirstHit) {
                    bestState = st;
                    foundStateForProp = true;
                    info() << "Using AtFirstHit state (3-hit seed) for inner propagation" << endmsg;
                    break;
                }
            }
        }
            
        // ===== INNER PROPAGATION =====
        if (foundStateForProp) {
            // Extract helix parameters from best state
            double d0 = bestState.D0 / 10.0;           // mm → cm
            double phi = bestState.phi;
            double omega = bestState.omega * 10.0;     // 1/mm → 1/cm
            double z0 = bestState.Z0 / 10.0;           // mm → cm
            double tanLambda = bestState.tanLambda;
            
            if (std::abs(omega) > 1e-10) {
                double radius = 1.0 / std::abs(omega);
                
                // Position at closest approach
                Eigen::Vector3d pos;
                if (finalTrack.trackerHits_size() > 0) {
                    auto firstHit = finalTrack.getTrackerHits(0);
                    pos = Eigen::Vector3d(
                        firstHit.getPosition()[0] / 10.0,  // mm to cm
                        firstHit.getPosition()[1] / 10.0,
                        firstHit.getPosition()[2] / 10.0
                    );
                }
                
                // Field at position
                dd4hep::Position fieldPos(pos.x(), pos.y(), pos.z());
                double bField = m_field.magneticField(fieldPos).z() / dd4hep::tesla;
                
                // Momentum
                double lambda = std::atan(tanLambda);
                double cosL = std::cos(lambda);
                double sinL = std::sin(lambda);
                double pT = 0.3 * std::abs(bField) * radius / 100.0;
                
                Eigen::Vector3d mom(
                    pT * std::cos(phi) * cosL,
                    pT * std::sin(phi) * cosL,
                    pT * sinL
                );
                
                info() << ">>> Attempting inner solenoid propagation for track " << trackNumber << "..." << endmsg;
                info() << "  Propagating track " << trackNumber << " to inner region..." << endmsg;
                info() << "    Outer R=" << std::sqrt(pos.x()*pos.x() + pos.y()*pos.y())
                    << " cm, |p|=" << mom.norm() << " GeV/c" << endmsg;
                
                // Propagate using RK4 to m_innerPropTargetRadius
                auto inner = propagateToInner(pos, mom, m_innerPropTargetRadius);
                
                if (inner.success) {
                    m_statInnerPropSuccess++;
                    double rxy = std::sqrt(
                        inner.finalPosition.x() * inner.finalPosition.x() +
                        inner.finalPosition.y() * inner.finalPosition.y()
                    );
                    double pMagFinal = inner.finalMomentum.norm();
                    double pRatio = pMagFinal / mom.norm();
                    
                    // ===== SAVE INNER STATE AS AtVertex =====
                    double d0_inner = std::sqrt(
                        inner.finalPosition.x() * inner.finalPosition.x() +
                        inner.finalPosition.y() * inner.finalPosition.y()
                    ) - (1.0 / std::abs(omega)) / 10.0;  // Approximation
                    
                    double phi_inner = std::atan2(
                        inner.finalPosition.y(),
                        inner.finalPosition.x()
                    );
                    
                    double omega_inner = omega * -1.176;  // Field flip: B_outer/B_inner ratio .. need to be fixed for automated ratio
                    
                    double z0_inner = inner.finalPosition.z();
                    
                    // Create inner track state at AtVertex (distinct from AtOther used by GenFit)
                    edm4hep::TrackState innerState = createTrackState(
                        d0_inner * 10.0,    // cm to mm
                        phi_inner,
                        omega_inner / 10.0, // 1/cm to 1/mm
                        z0_inner * 10.0,    // cm to mm
                        tanLambda,
                        edm4hep::TrackState::AtVertex  // ← Inner propagation result
                    );
                    
                    // Add to track
                    finalTrack.addToTrackStates(innerState);
                    m_statStateAtVertex++;
                    
                    info() << "  ✓ Inner propagation SUCCESS → saved as AtVertex:" << endmsg;
                    info() << "    Final R=" << rxy << " cm, |p|=" << pMagFinal 
                        << " GeV/c (ratio=" << pRatio << ")" << endmsg;
                    info() << "    Arc length: " << inner.arcLength << " cm, steps: "
                        << inner.numSteps << endmsg;
                    
                    // ===== TRUTH COMPARISON =====
                    try {
                        if (m_mcParticles.exist()) {
                            const auto* mcParts = m_mcParticles.get();
                            if (mcParts && !mcParts->empty()) {
                                edm4hep::Track trackView = finalTrack;
                                validateTrackWithTruth(trackView, *mcParts, trackNumber);
                            }
                        }
                    } catch (const std::exception& e) {
                        debug() << "Could not access MCParticles for truth comparison: " << e.what() << endmsg;
                    }
                } else {
                    debug() << "  ✗ Inner propagation failed: " << inner.message << endmsg;
                }
            }
        }
        } // end if (m_doInnerPropagation)
        // ===== END INNER PROPAGATION =====

        // ======================= Printing final track=============================
        info() << "Successfully created track " << trackNumber
            << " with " << finalTrack.trackerHits_size() << " hits"
            << " and " << finalTrack.trackStates_size() << " track states" << endmsg;

        // ---- per-track run statistics ----
        m_statTracksReconstructed++;
        bool hasFirstHit = false, hasLastHit = false, hasCalorimeter = false, hasOther = false, hasVertex = false;
        for (int sIdx = 0; sIdx < finalTrack.trackStates_size(); ++sIdx) {
            auto ts = finalTrack.getTrackStates(sIdx);
            switch (ts.location) {
                case edm4hep::TrackState::AtFirstHit:    hasFirstHit    = true; break;
                case edm4hep::TrackState::AtLastHit:     hasLastHit     = true; break;
                case edm4hep::TrackState::AtCalorimeter: hasCalorimeter = true; break;
                case edm4hep::TrackState::AtOther:       hasOther       = true; break;
                case edm4hep::TrackState::AtVertex:      hasVertex      = true; break;
                default: break;
            }
        }
        if (hasFirstHit)    m_statStateAtFirstHit++;
        if (hasLastHit)     m_statStateAtLastHit++;
        if (hasCalorimeter) m_statStateAtCalorimeter++;
        if (hasOther)       m_statStateAtOther++;
        if (hasVertex)      m_statStateAtVertex++;
        if (hasLastHit) m_statFourHitTracks++; else m_statThreeHitTracks++;
        // Count charge from the best available state
        // (prefer AtFirstHit, fall back to AtOther, then AtLastHit)
        {
            int chargePriority = 99;
            double chargeOmega = 0.0;
            for (int sIdx = 0; sIdx < finalTrack.trackStates_size(); ++sIdx) {
                auto ts = finalTrack.getTrackStates(sIdx);
                int prio = (ts.location == edm4hep::TrackState::AtFirstHit) ? 0 :
                           (ts.location == edm4hep::TrackState::AtOther)    ? 1 :
                           (ts.location == edm4hep::TrackState::AtLastHit)  ? 2 : 99;
                if (prio < chargePriority) {
                    chargePriority = prio;
                    chargeOmega = ts.omega;
                }
            }
            if (chargePriority < 99) {
                if (chargeOmega < 0.0) m_statPositiveCharge++;
                else                   m_statNegativeCharge++;
            }
        }
        // Continue to next iteration to find more tracks
    } // End of main while loop

    info() << "Track reconstruction complete: created " << finalTracks.size() 
           << " final tracks across " << (trackNumber-1) << " iterations" << endmsg;
}

// function to calculate circle center and radius using the Direct Formula Method
bool DisplacedTracking::calculateCircleCenterDirect(
    double x1, double y1, double x2, double y2, double x3, double y3,
    double& x0, double& y0, double& radius) const{
    
    // Calculate the determinant
    double D = 2 * (x1 * (y2 - y3) + x2 * (y3 - y1) + x3 * (y1 - y2));
    
    // Check if points are collinear (determinant would be zero)
    if (std::abs(D) < 1e-10) {
        return false;
    }
    
    // Calculate center coordinates
    x0 = ((x1*x1 + y1*y1) * (y2 - y3) + 
          (x2*x2 + y2*y2) * (y3 - y1) + 
          (x3*x3 + y3*y3) * (y1 - y2)) / D;
    
    y0 = ((x1*x1 + y1*y1) * (x3 - x2) + 
          (x2*x2 + y2*y2) * (x1 - x3) + 
          (x3*x3 + y3*y3) * (x2 - x1)) / D;
    
    // Calculate radius as distance from center to any of the points
    radius = std::sqrt((x1 - x0)*(x1 - x0) + (y1 - y0)*(y1 - y0));
    
    return true;
}

// Improved sagitta method with circle center calculation
double DisplacedTracking::calculateSagitta(const Eigen::Vector3d& p1, 
                                        const Eigen::Vector3d& p2, 
                                        const Eigen::Vector3d& p3) const{
    // Project points onto the xy plane (transverse plane)
    Eigen::Vector2d p1_2d(p1.x(), p1.y());
    Eigen::Vector2d p2_2d(p2.x(), p2.y());
    Eigen::Vector2d p3_2d(p3.x(), p3.y());
    
    // Vector from p1 to p3 (chord)
    Eigen::Vector2d chord = p3_2d - p1_2d;
    double chordLength = chord.norm();
    
    // Unit vector along the chord
    Eigen::Vector2d chordDir = chord / chordLength;
    
    // Vector from p1 to p2
    Eigen::Vector2d v1to2 = p2_2d - p1_2d;
    
    // Project v1to2 onto the chord direction
    double projection = v1to2.dot(chordDir);
    
    // Calculate the perpendicular distance (sagitta)
    Eigen::Vector2d projectionVec = projection * chordDir;
    Eigen::Vector2d perpVec = v1to2 - projectionVec;
    double sagitta = perpVec.norm();
    
    return sagitta;
}

// Function to calculate circle center using sagitta method
bool DisplacedTracking::calculateCircleCenterSagitta(
    const Eigen::Vector3d& p1, const Eigen::Vector3d& p2, const Eigen::Vector3d& p3,
    double& x0, double& y0, double& radius) const{
    
    // Project points onto the xy plane (transverse plane)
    Eigen::Vector2d p1_2d(p1.x(), p1.y());
    Eigen::Vector2d p2_2d(p2.x(), p2.y());
    Eigen::Vector2d p3_2d(p3.x(), p3.y());
    
    // Vector from p1 to p3 (chord)
    Eigen::Vector2d chord = p3_2d - p1_2d;
    double chordLength = chord.norm();
    
    // Check if points are too close
    if (chordLength < 1e-6) {
        return false;
    }
    
    // Vector from p1 to p2
    Eigen::Vector2d v1to2 = p2_2d - p1_2d;
    
    // Project v1to2 onto the chord direction
    Eigen::Vector2d chordDir = chord / chordLength;
    double projection = v1to2.dot(chordDir);
    
    // Calculate the perpendicular vector and sagitta
    Eigen::Vector2d projectionVec = projection * chordDir;
    Eigen::Vector2d perpVec = v1to2 - projectionVec;
    double sagitta = perpVec.norm();
    
    // If sagitta is too small, points are nearly collinear
    if (sagitta < 1e-6) {
        return false;
    }
    
    // Calculate radius using sagitta formula
    radius = (chordLength * chordLength) / (8 * sagitta) + (sagitta / 2);
    
    // Find the midpoint of chord p1-p3
    Eigen::Vector2d midpoint = (p1_2d + p3_2d) / 2.0;
    
    // Perpendicular direction to chord (normalized)
    Eigen::Vector2d perpDir(-chordDir.y(), chordDir.x());
    
    // Determine on which side of the chord the center lies
    // Cross product to check if p2 is "above" or "below" the chord
    double crossProduct = chordDir.x() * (p2_2d.y() - p1_2d.y()) - 
                         chordDir.y() * (p2_2d.x() - p1_2d.x());
    double directionFactor = (crossProduct > 0) ? 1.0 : -1.0;
    
    // Height from chord to center
    double height = std::sqrt(radius * radius - (chordLength / 2.0) * (chordLength / 2.0));
    
    // Calculate center coordinates
    x0 = midpoint.x() + directionFactor * height * perpDir.x();
    y0 = midpoint.y() + directionFactor * height * perpDir.y();
    
    return true;
}

bool DisplacedTracking::createTripletSeed(
    const edm4hep::TrackerHitPlane& hit1,
    const edm4hep::TrackerHitPlane& hit2,
    const edm4hep::TrackerHitPlane& hit3,
    edm4hep::TrackCollection* tracks,
    edm4hep::TrackerHitPlaneCollection& outputHits,
    const LinkPropagator& propagateLink,
    std::vector<bool>& usedHits,
    size_t idx1, size_t idx2, size_t idx3) const {
    
    // Get hit positions
    const auto& pos1 = hit1.getPosition();
    const auto& pos2 = hit2.getPosition();
    const auto& pos3 = hit3.getPosition();
    
    // Get surfaces for each hit
    const dd4hep::rec::Surface* surf1 = findSurface(hit1);
    const dd4hep::rec::Surface* surf2 = findSurface(hit2);
    const dd4hep::rec::Surface* surf3 = findSurface(hit3);
    
    if (!surf1 || !surf2 || !surf3) {
        debug() << "Could not find surfaces for triplet hits" << endmsg;
        return false; // Couldn't find surfaces
    }
    
    // Convert to Eigen vectors in cm
    Eigen::Vector3d p1(pos1[0] / 10.0, pos1[1] / 10.0, pos1[2] / 10.0);
    Eigen::Vector3d p2(pos2[0] / 10.0, pos2[1] / 10.0, pos2[2] / 10.0);
    Eigen::Vector3d p3(pos3[0] / 10.0, pos3[1] / 10.0, pos3[2] / 10.0);
    
    // Check if hits are spatially compatible
    double maxDist = m_maxDist; // cm
    if ((p2 - p1).norm() > maxDist || (p3 - p2).norm() > maxDist) {
        debug() << "Hits too far apart spatially, more than 1 m." << endmsg;
        return false;
    }
    
    // Check angle consistency
    Eigen::Vector3d v1 = p2 - p1;
    Eigen::Vector3d v2 = p3 - p2;
    v1.normalize();
    v2.normalize();
    /*
    double cosAngle = v1.dot(v2);
    if (cosAngle < 0.75) { // Allow up to about 25 degrees deviation
        debug() << "Hits not along a consistent path, angle too large" << endmsg;
        return false;
    }
    */
    debug() << "Fitting circle through points (cm): " 
            << "(" << p1.x() << "," << p1.y() << "), "
            << "(" << p2.x() << "," << p2.y() << "), "
            << "(" << p3.x() << "," << p3.y() << ")" << endmsg;
    
    // Calculate using Direct Formula Method
    double x0_direct, y0_direct, radius_direct;
    bool directValid = calculateCircleCenterDirect(
        p1.x(), p1.y(), p2.x(), p2.y(), p3.x(), p3.y(),
        x0_direct, y0_direct, radius_direct);
    
    if (!directValid) {
        debug() << "Direct formula method failed" << endmsg;
        return false;
    }
    
    debug() << "Direct formula method: center=(" << x0_direct << "," << y0_direct 
            << "), radius=" << radius_direct << " cm" << endmsg;
    
    // Calculate using sagitta method
    double sagitta = calculateSagitta(p1, p2, p3);
    
    // Calculate chord length
    Eigen::Vector2d p1_2d(p1.x(), p1.y());
    Eigen::Vector2d p3_2d(p3.x(), p3.y());
    double chordLength = (p3_2d - p1_2d).norm();
    
    // Simplified sagitta radius calculation
    double sagittaRadius = (chordLength * chordLength) / (8 * sagitta);
    
    debug() << "Sagitta method: sagitta = " << sagitta << " cm, chord = " 
            << chordLength << " cm, radius = " << sagittaRadius << " cm" << endmsg;
    
    // Calculate full sagitta center and radius
    double x0_sagitta, y0_sagitta, radius_sagitta;
    bool sagittaValid = calculateCircleCenterSagitta(
        p1, p2, p3, x0_sagitta, y0_sagitta, radius_sagitta);
    
    if (!sagittaValid) {
        debug() << "Sagitta center calculation failed" << endmsg;
        return false;
    }
    
    debug() << "Sagitta full method: center=(" << x0_sagitta << "," << y0_sagitta 
            << "), radius=" << radius_sagitta << " cm" << endmsg;
    
    // Calculate magnetic field
    dd4hep::Position fieldPos((p1.x() + p2.x() + p3.x())/3.0, 
                            (p1.y() + p2.y() + p3.y())/3.0, 
                            (p1.z() + p2.z() + p3.z())/3.0);
    double actualBz = m_field.magneticField(fieldPos).z() / dd4hep::tesla;
    
    debug() << "Magnetic field: actual=" << actualBz << " Tesla" << endmsg;
    
    // Calculate pT using all methods for comparison
    double pT_direct = 0.3 * std::abs(actualBz) * radius_direct / 100.0;
    double pT_sagitta = 0.3 * std::abs(actualBz) * sagittaRadius / 100.0;
    double pT_sagitta_full = 0.3 * std::abs(actualBz) * radius_sagitta / 100.0;
    
    debug() << "Comparison of pT estimates:" << endmsg;
    debug() << "  Direct formula pT: " << pT_direct << " GeV/c" << endmsg;
    debug() << "  Sagitta simple pT: " << pT_sagitta << " GeV/c" << endmsg;
    debug() << "  Sagitta full pT: " << pT_sagitta_full << " GeV/c" << endmsg;

    // --- Geometry sanity check: direct vs sagitta radius should agree within a factor of 2 ---
    // If they disagree wildly, the triplet is nearly collinear / numerically unstable.
    // We log a warning but do NOT reject — the sagitta method is more robust for high-pT tracks.
    {
        double rRatio = (radius_direct > 1e-3) ? radius_sagitta / radius_direct : 0.0;
        if (rRatio < 0.3 || rRatio > 3.0) {
            m_statTripletBadGeom++;
            debug() << "WARNING: radius disagreement — direct=" << radius_direct
                    << " cm  sagitta=" << radius_sagitta << " cm  ratio=" << rRatio
                    << "  (numerically marginal triplet — proceeding with sagitta)" << endmsg;
        }
    }
    /*
    // Use sagitta method for track parameters
    double radius = radius_direct;
    double x0 = x0_direct;
    double y0 = y0_direct;
    double pT = pT_direct;
    */
    // Use sagitta method for track parameters
    double radius = radius_sagitta;
    double x0 = x0_sagitta;
    double y0 = y0_sagitta;
    double pT = pT_sagitta_full;
    
    // Determine helix direction and charge
    double phi1 = std::atan2(p1.y() - y0, p1.x() - x0);
    double phi2 = std::atan2(p2.y() - y0, p2.x() - x0);
    double phi3 = std::atan2(p3.y() - y0, p3.x() - x0);
    
    // Unwrap angles
    if (phi2 - phi1 > M_PI) phi2 -= 2*M_PI;
    if (phi2 - phi1 < -M_PI) phi2 += 2*M_PI;
    if (phi3 - phi2 > M_PI) phi3 -= 2*M_PI;
    if (phi3 - phi2 < -M_PI) phi3 += 2*M_PI;
    
    bool clockwise = (phi3 < phi1);

    // Geometric rotation direction gives charge sign only when B > 0.
    // For B < 0 the Lorentz force reverses, so the rotation sense flips:
    //   q_geometric = clockwise ? +1 : -1   (assumes B > 0)
    //   true charge = q_geometric * sign(Bz)
    // We read actualBz directly from the detector so this is field-agnostic.
    double bSign  = (actualBz >= 0.0) ? 1.0 : -1.0;
    double charge = bSign * (clockwise ? 1.0 : -1.0);

    debug() << "Track direction: " << (clockwise ? "clockwise" : "counter-clockwise")
            << ", Bz=" << actualBz << " T, charge: " << charge << endmsg;
    
    // Fit z-component
    double s1 = 0;
    double s2 = radius * std::abs(phi2 - phi1);
    double s3 = radius * std::abs(phi3 - phi1);
    
    double a, b;
    fitLine(s1, p1.z(), s2, p2.z(), s3, p3.z(), a, b);
    
    debug() << "Z-fit: z = " << a << " * s + " << b << endmsg;
    
    double theta = std::atan2(1.0, a);
    double eta = -std::log(std::tan(theta/2.0));
    
    debug() << "Track angles: theta=" << theta << ", eta=" << eta << endmsg;
    
    // Calculate impact parameters using two-segment model
    // Query inner and outer field strengths from DD4hep at representative positions
    // Inner: a point well inside the solenoid boundary (R=230 cm)
    double innerFieldStrength = m_field.magneticField(dd4hep::Position(0.0, 0.0, 0.0)).z() / dd4hep::tesla;
    // Outer: a point well outside the solenoid boundary, near the muon system
    double outerFieldStrength = m_field.magneticField(dd4hep::Position(350.0, 0.0, 0.0)).z() / dd4hep::tesla;

    double d0 = std::sqrt(std::pow(x0, 2) + std::pow(y0, 2)) - radius;
    /*
    // Use the two-segment model for d0 calculation
    double d0 = calculateImpactParameter(x0, y0, radius, clockwise, 
                                   innerFieldStrength, outerFieldStrength,
                                   p1, p2, p3);
    */
    double z0 = b; // z0 calculation remains the same
    
    debug() << "Impact parameters: d0=" << d0 << " cm, z0=" << z0 << " cm" << endmsg;
    
    // Track parameters
    double qOverPt = charge / pT;
    double phi = std::atan2(y0, x0) + (clockwise ? -M_PI/2 : M_PI/2);
    
    // Normalize phi
    if (phi > M_PI) phi -= 2*M_PI;
    if (phi < -M_PI) phi += 2*M_PI;
    
    // Convert to EDM4hep parameters
    double d0_mm = d0 * 10.0;  // Convert from cm to mm
    double z0_mm = z0 * 10.0;  // Convert from cm to mm
    
    // omega = 1/R with correct sign (R in mm)
    double omega = charge / (radius * 10.0);  // 1/mm
    
    // tanLambda = tan of dip angle
    double tanLambda = std::sinh(eta);
    
    debug() << "Track parameters: (" 
            << qOverPt << ", " << phi << ", " << eta << ", " << d0 << ", " << z0 << ")" << endmsg;
    
    debug() << "EDM4hep parameters: d0=" << d0_mm << " mm, phi=" << phi 
            << ", omega=" << omega << " 1/mm, z0=" << z0_mm 
            << " mm, tanLambda=" << tanLambda << endmsg;

    // Create EDM4hep track
    auto edm_track = tracks->create();
    
    // Create track state at IP
    edm4hep::TrackState state = createTrackState(
        d0_mm, phi, omega, z0_mm, tanLambda, edm4hep::TrackState::AtFirstHit);

    // Set reference point to the first hit position
    const auto& firstHitPos = hit1.getPosition();
    state.referencePoint = edm4hep::Vector3f(firstHitPos[0], firstHitPos[1], firstHitPos[2]);

    // ── Propagate hit uncertainties to track-parameter covariance ──────────────
    // All EDM4hep cov units:
    //   D0, Z0  in mm²;  phi, tanLambda in rad²;  omega in (1/mm)².
    {
        // Hit position uncertainties in cm
        auto getSigmaXY = [](const edm4hep::TrackerHitPlane& h) -> double {
            double su = h.getDu(); 
            double sv = h.getDv();
            if (su <= 0) su = 0.4; // mm
            if (sv <= 0) sv = 0.4;
            return std::sqrt(su*su + sv*sv) / 10.0; // cm
        };

        auto getSigmaZ = [](const edm4hep::TrackerHitPlane& h) -> double {
            double sv = h.getDv();
            if (sv <= 0) sv = 0.4;
            return sv / 10.0; // cm
        };

        double sigma_xy1 = getSigmaXY(hit1);
        double sigma_xy2 = getSigmaXY(hit2);
        double sigma_xy3 = getSigmaXY(hit3);

        double sigma_z1  = getSigmaZ(hit1);
        double sigma_z3  = getSigmaZ(hit3);

        // lever arm
        double chordL = std::sqrt(std::pow(p3.x()-p1.x(),2) +
                                std::pow(p3.y()-p1.y(),2));

        // --- omega ---
        double sigma_R_cm = (chordL > 1e-3)
            ? 2.0 * sigma_xy2 * radius / chordL
            : 0.5 * radius;

        double dOmega_dR = -charge / (radius * radius * 10.0);
        double var_omega  = dOmega_dR * dOmega_dR *
                            sigma_R_cm * sigma_R_cm;

        double omega_abs_floor = std::pow(0.2 * std::abs(omega), 2);

        debug() << "Omega covariance values AtFirstHit:"
                << "  σ_ω="  << std::sqrt(var_omega)  << " 1/mm"
                << "  omega_abs_floor="   << std::sqrt(omega_abs_floor)<< " 1/mm"
                << endmsg;

        var_omega = std::max({var_omega,
                            omega_abs_floor});

        // --- phi ---
        double var_phi = (chordL > 1e-6)
            ? (sigma_xy1 * sigma_xy1) / (chordL * chordL)
            : 0.01;

         debug() << "Phi covariance values AtFirstHit:"
                << "  σ_phi=" << std::sqrt(var_phi)  << " rad"
                << "  phi_floor= 0.1 rad"
                << endmsg;

        var_phi = std::max(var_phi, 0.01);  // floor for sigma = 0.1 

        // --- d0 ---
        double avg_sigma_xy_mm =
            (sigma_xy1 + sigma_xy2 + sigma_xy3) / 3.0 * 10.0;

        double var_d0 = avg_sigma_xy_mm * avg_sigma_xy_mm;

        // --- z0 / tanLambda ---
        double avg_sigma_z_mm =
            (sigma_z1 + sigma_z3) / 2.0 * 10.0;

        double chordL_mm = chordL * 10.0;

        double var_tanL = (chordL_mm > 1e-6)
            ? (avg_sigma_z_mm * avg_sigma_z_mm) /
            (chordL_mm * chordL_mm)
            : 0.01;

        debug() << "Tan lamda covariance values AtFirstHit:"
                << "  σ_tanL="<< std::sqrt(var_tanL) 
                << "  tanL_floor= 0.1"
                << endmsg;

        var_tanL = std::max(var_tanL, 0.01); // floor for sigma = 0.1 

        double var_z0 = avg_sigma_z_mm * avg_sigma_z_mm;

        // reset covariance
        for (int i = 0; i < 6; ++i)
            for (int j = 0; j <= i; ++j)
                state.setCovMatrix(0.0,
                    static_cast<edm4hep::TrackParams>(i),
                    static_cast<edm4hep::TrackParams>(j));

        state.setCovMatrix(var_d0,
            edm4hep::TrackParams::d0,
            edm4hep::TrackParams::d0);

        state.setCovMatrix(var_phi,
            edm4hep::TrackParams::phi,
            edm4hep::TrackParams::phi);

        state.setCovMatrix(var_omega,
            edm4hep::TrackParams::omega,
            edm4hep::TrackParams::omega);

        state.setCovMatrix(var_z0,
            edm4hep::TrackParams::z0,
            edm4hep::TrackParams::z0);

        state.setCovMatrix(var_tanL,
            edm4hep::TrackParams::tanLambda,
            edm4hep::TrackParams::tanLambda);

        debug() << "AtFirstHit covariance from hit uncertainties:"
                << "  σ_d0="  << std::sqrt(var_d0)  << " mm"
                << "  σ_phi=" << std::sqrt(var_phi)  << " rad"
                << "  σ_ω="   << std::sqrt(var_omega)<< " 1/mm"
                << "  σ_z0="  << std::sqrt(var_z0)   << " mm"
                << "  σ_tanL="<< std::sqrt(var_tanL) << endmsg;
    }

    // Add track state to track
    edm_track.addToTrackStates(state);
    
    // Copy input hits into the output collection so PODIO can resolve the relation
    auto outHit1 = outputHits.create();
    outHit1.setCellID(hit1.getCellID()); outHit1.setTime(hit1.getTime());
    outHit1.setEDep(hit1.getEDep()); outHit1.setEDepError(hit1.getEDepError());
    outHit1.setPosition(hit1.getPosition()); outHit1.setCovMatrix(hit1.getCovMatrix());
    outHit1.setDu(hit1.getDu()); outHit1.setDv(hit1.getDv());
    auto outHit2 = outputHits.create();
    outHit2.setCellID(hit2.getCellID()); outHit2.setTime(hit2.getTime());
    outHit2.setEDep(hit2.getEDep()); outHit2.setEDepError(hit2.getEDepError());
    outHit2.setPosition(hit2.getPosition()); outHit2.setCovMatrix(hit2.getCovMatrix());
    outHit2.setDu(hit2.getDu()); outHit2.setDv(hit2.getDv());
    auto outHit3 = outputHits.create();
    outHit3.setCellID(hit3.getCellID()); outHit3.setTime(hit3.getTime());
    outHit3.setEDep(hit3.getEDep()); outHit3.setEDepError(hit3.getEDepError());
    outHit3.setPosition(hit3.getPosition()); outHit3.setCovMatrix(hit3.getCovMatrix());
    outHit3.setDu(hit3.getDu()); outHit3.setDv(hit3.getDv());
    propagateLink(outHit1, hit1);
    propagateLink(outHit2, hit2);
    propagateLink(outHit3, hit3);
    edm_track.addToTrackerHits(outHit1);
    edm_track.addToTrackerHits(outHit2);
    edm_track.addToTrackerHits(outHit3);
    
    // Set chi2 and ndf
    edm_track.setChi2(-1.0);  // Initial seed has no chi2 yet
   // edm_track.setNdf(-1);  // 
    
    // Mark hits as used
    usedHits[idx1] = true;
    usedHits[idx2] = true;
    usedHits[idx3] = true;
    
    debug() << "Created valid EDM4hep track with 3 hits" << endmsg;
    return true;
}
/*
// This calculates d0 using a two-segment track model for field transitions
double DisplacedTracking::calculateImpactParameter(
    double x0, double y0, double radius, bool clockwise,
    double innerFieldStrength, double outerFieldStrength,
    const Eigen::Vector3d& p1, const Eigen::Vector3d& p2, const Eigen::Vector3d& p3) const{
    
    // Constants
    const double solenoidRadius = 230.0; // cm - radius of the solenoid (matches detector geometry)
    double centerToOriginDistance = std::sqrt(x0*x0 + y0*y0);
    
    debug() << "---- Analytical Impact Parameter Calculation ----" << endmsg;
    debug() << "Outer circle center: (" << x0 << ", " << y0 << ") cm" << endmsg;
    debug() << "Outer circle radius: " << radius << " cm" << endmsg;
    debug() << "Track curvature direction: " << (clockwise ? "clockwise" : "counter-clockwise") << endmsg;
    debug() << "Center to origin distance: " << centerToOriginDistance << " cm" << endmsg;
    
    // Check if the particle would cross the solenoid boundary
    bool intersectsSolenoid = false;
    double intersectionX1 = 0, intersectionY1 = 0;
    double intersectionX2 = 0, intersectionY2 = 0;
    int numIntersections = 0;
    
    // Case 1: The circle intersects with the solenoid boundary
    if (std::abs(centerToOriginDistance - solenoidRadius) < radius && 
        radius < centerToOriginDistance + solenoidRadius) {
        // Circle intersects with solenoid boundary
        intersectsSolenoid = true;
        
        // Calculate the circle-circle intersection points
        double centerAngle = std::atan2(y0, x0);
        double distanceToChord = (centerToOriginDistance*centerToOriginDistance + 
                                 solenoidRadius*solenoidRadius - radius*radius) / 
                                 (2.0 * centerToOriginDistance);
        double halfChordLength = std::sqrt(solenoidRadius*solenoidRadius - distanceToChord*distanceToChord);
        
        // Calculate the position of the center of the chord
        double chordCenterX = distanceToChord * std::cos(centerAngle);
        double chordCenterY = distanceToChord * std::sin(centerAngle);
        
        // Calculate the angle perpendicular to the center angle
        double perpAngle = centerAngle + M_PI/2.0;
        
        // Calculate the two intersection points
        intersectionX1 = chordCenterX + halfChordLength * std::cos(perpAngle);
        intersectionY1 = chordCenterY + halfChordLength * std::sin(perpAngle);
        intersectionX2 = chordCenterX - halfChordLength * std::cos(perpAngle);
        intersectionY2 = chordCenterY - halfChordLength * std::sin(perpAngle);
        
        numIntersections = 2;
        
        debug() << "Track circle intersects with solenoid boundary" << endmsg;
        debug() << "Intersection point 1: (" << intersectionX1 << ", " << intersectionY1 << ") cm" << endmsg;
        debug() << "Intersection point 2: (" << intersectionX2 << ", " << intersectionY2 << ") cm" << endmsg;
    } 
    // Case 2: The circle is completely inside the solenoid
    else if (centerToOriginDistance + radius < solenoidRadius) {
        debug() << "Track circle is completely inside the solenoid - using inner circle directly" << endmsg;
        intersectsSolenoid = false;
        
        // For particles inside the solenoid, calculate the actual impact parameter
        // Impact parameter is the closest distance from the origin to the circle
        double d0 = std::abs(centerToOriginDistance - radius);
        
        // Apply sign convention based on curvature
        if (!clockwise) {
            d0 = -d0;
        }
        
        debug() << "Particle is fully inside the solenoid" << endmsg;
        debug() << "Distance from center to origin: " << centerToOriginDistance << " cm" << endmsg;
        debug() << "Track radius: " << radius << " cm" << endmsg;
        debug() << "Calculated d0 = " << d0 << " cm" << endmsg;
        debug() << "---- End Analytical Impact Parameter Calculation ----" << endmsg;
        
        return d0;
    }
    // Case 3: The circle is completely outside the solenoid
    else {
        debug() << "Track circle is completely outside the solenoid - checking if it crosses" << endmsg;
        
        // Calculate the closest approach of the circle to the origin
        double closestApproach = std::abs(centerToOriginDistance - radius);
        
        if (closestApproach < solenoidRadius) {
            // Track would pass through the solenoid when extrapolated
            debug() << "Extrapolated track would pass through the solenoid" << endmsg;
            
            // Calculate the intersection point with the solenoid
            double approachAngle = std::atan2(y0, x0);
            if (centerToOriginDistance < radius) {
                // If the origin is inside the circle, flip the approach angle
                approachAngle += M_PI;
            }
            
            // Vector from origin towards closest approach point
            intersectionX1 = solenoidRadius * std::cos(approachAngle);
            intersectionY1 = solenoidRadius * std::sin(approachAngle);
            numIntersections = 1;
            
            intersectsSolenoid = true;
            debug() << "Assuming intersection point: (" << intersectionX1 << ", " << intersectionY1 << ") cm" << endmsg;
        } else {
            // Track doesn't cross solenoid when extrapolated back
            debug() << "Track does not cross solenoid when extrapolated" << endmsg;
            debug() << "Closest approach to origin: " << closestApproach << " cm" << endmsg;
            debug() << "This indicates a displaced vertex outside the solenoid" << endmsg;
            
            // For particles not crossing the solenoid, use the original d0 calculation
            double d0 = centerToOriginDistance - radius;
            if (!clockwise) {
                d0 = -d0;  // Apply sign convention
            }
            
            debug() << "Using outer circle only for impact parameter" << endmsg;
            debug() << "Calculated d0 = " << d0 << " cm" << endmsg;
            debug() << "---- End Analytical Impact Parameter Calculation ----" << endmsg;
            
            return d0;
        }
    }
        
    // If we're here, the track intersects the solenoid boundary
    // We need to select the intersection point closest to the detector hits
    double intersectionX, intersectionY;
    if (numIntersections == 2) {
        // Calculate the angles from circle center to first hit and from center to both intersections
        double hitAngle = std::atan2(p1.y() - y0, p1.x() - x0);
        double intersection1Angle = std::atan2(intersectionY1 - y0, intersectionX1 - x0);
        double intersection2Angle = std::atan2(intersectionY2 - y0, intersectionX2 - x0);
        
        // Normalize angles to [0, 2π)
        hitAngle = (hitAngle < 0) ? hitAngle + 2 * M_PI : hitAngle;
        intersection1Angle = (intersection1Angle < 0) ? intersection1Angle + 2 * M_PI : intersection1Angle;
        intersection2Angle = (intersection2Angle < 0) ? intersection2Angle + 2 * M_PI : intersection2Angle;
        
        // Calculate the angle differences in the direction of motion
        double angleDiff1, angleDiff2;
        if (clockwise) {
            angleDiff1 = hitAngle - intersection1Angle;
            angleDiff2 = hitAngle - intersection2Angle;
        } else {
            angleDiff1 = intersection1Angle - hitAngle;
            angleDiff2 = intersection2Angle - hitAngle;
        }
        
        // Normalize angle differences to [-π, π]
        if (angleDiff1 > M_PI) angleDiff1 -= 2 * M_PI;
        if (angleDiff1 < -M_PI) angleDiff1 += 2 * M_PI;
        if (angleDiff2 > M_PI) angleDiff2 -= 2 * M_PI;
        if (angleDiff2 < -M_PI) angleDiff2 += 2 * M_PI;
        
        // Calculate cross products for additional info
        double cross1 = intersectionX1 * (y0 - intersectionY1) - intersectionY1 * (x0 - intersectionX1);
        double cross2 = intersectionX2 * (y0 - intersectionY2) - intersectionY2 * (x0 - intersectionX2);
        debug() << "Cross product 1: " << cross1 << ", Cross product 2: " << cross2 << endmsg;
        
        // The correct intersection should be in the opposite direction of motion from the hit
        // and have a positive angle difference
        if (angleDiff1 > 0 && (angleDiff2 < 0 || angleDiff1 < angleDiff2)) {
            intersectionX = intersectionX1;
            intersectionY = intersectionY1;
            debug() << "Chose intersection 1 (based on particle direction, angle diff=" << angleDiff1 << ")" << endmsg;
        } else {
            intersectionX = intersectionX1;
            intersectionY = intersectionY1;
            debug() << "Chose intersection 1 (based on particle direction, angle diff=" << angleDiff2 << ")" << endmsg;
        }
    } else {
        intersectionX = intersectionX1;
        intersectionY = intersectionY1;
    }
    
    // Calculate the tangent direction precisely at the intersection point
    // This is perpendicular to the radius vector from outer circle center
    double outerRadiusX = intersectionX - x0;
    double outerRadiusY = intersectionY - y0;
    double outerRadiusLength = std::sqrt(outerRadiusX*outerRadiusX + outerRadiusY*outerRadiusY);
    outerRadiusX /= outerRadiusLength; // Normalize
    outerRadiusY /= outerRadiusLength;
    
    // Calculate tangent vector (perpendicular to radius)
    double tangentX, tangentY;
    if (clockwise) {
        tangentX = -outerRadiusY;
        tangentY = outerRadiusX;
    } else {
        tangentX = outerRadiusY;
        tangentY = -outerRadiusX;
    }
    
    // Verify tangent is perpendicular to outer radius vector
    double outerDotProduct = tangentX * outerRadiusX + tangentY * outerRadiusY;
    debug() << "Outer radius-tangent dot product (should be ~0): " << outerDotProduct << endmsg;
    
    // Calculate inner circle radius based on field strength ratio
    double innerRadius = radius * std::abs(outerFieldStrength / innerFieldStrength);
    bool innerClockwise = !clockwise; // Curvature direction flips with field
    
    debug() << "Inner field: " << innerFieldStrength << " T, Outer field: " << outerFieldStrength << " T" << endmsg;
    debug() << "Inner circle radius: " << innerRadius << " cm" << endmsg;
    debug() << "Inner circle direction: " << (innerClockwise ? "clockwise" : "counter-clockwise") << endmsg;
    
    // Calculate the inner circle center
    // For fields of opposite signs, the center must be on the other side
    double innerX0, innerY0;
    bool oppositeSigns = (innerFieldStrength * outerFieldStrength < 0);
    
    // Direction factor for inner circle calculation
    double dirFactor = (innerClockwise) ? 1.0 : -1.0;
    if (oppositeSigns) dirFactor = -dirFactor;
    
    // Calculate inner circle center
    innerX0 = intersectionX + dirFactor * innerRadius * tangentY;
    innerY0 = intersectionY - dirFactor * innerRadius * tangentX;
    
    debug() << "Inner circle center: (" << innerX0 << ", " << innerY0 << ") cm" << endmsg;
    
    // Calculate and verify inner radius vector at intersection
    double innerRadiusX = intersectionX - innerX0;
    double innerRadiusY = intersectionY - innerY0;
    double innerRadiusLength = std::sqrt(innerRadiusX*innerRadiusX + innerRadiusY*innerRadiusY);
    innerRadiusX /= innerRadiusLength; // Normalize
    innerRadiusY /= innerRadiusLength;
    
    // Verify the inner radius is perpendicular to the tangent
    double innerDotProduct = tangentX * innerRadiusX + tangentY * innerRadiusY;
    debug() << "Inner radius-tangent dot product (should be ~0): " << innerDotProduct << endmsg;
    
    // Verify the calculated radius matches the expected inner radius
    debug() << "Calculated inner radius length: " << innerRadiusLength << " cm" << endmsg;
    debug() << "Expected inner radius: " << innerRadius << " cm" << endmsg;
    debug() << "Radius difference: " << std::abs(innerRadiusLength - innerRadius) << " cm" << endmsg;
    
    // Calculate distance from inner circle center to origin
    double innerCenterToOriginDistance = std::sqrt(innerX0*innerX0 + innerY0*innerY0);
    
    // Calculate how close the inner circle comes to the origin
    double distanceToOrigin = std::abs(innerCenterToOriginDistance - innerRadius);
    debug() << "Distance from inner circle to origin: " << distanceToOrigin << " cm" << endmsg;
    debug() << "Inner circle to origin distance: " << innerCenterToOriginDistance << " cm" << endmsg;
    debug() << "Inner circle radius: " << innerRadius << " cm" << endmsg;
    
    // Calculate the ratio between the distance to origin and the radius
    // For a perfect track from the origin, this ratio should be exactly 1.0
    double distanceRatio = innerCenterToOriginDistance / innerRadius;
    debug() << "Distance ratio (should be close to 1.0): " << distanceRatio << endmsg;
    
    // Calculate precision metrics
    debug() << "Precision error: " << std::abs(distanceRatio - 1.0) * 100.0 << "%" << endmsg;
    
    // Calculate the impact parameter (closest approach to origin)
    double d0 = distanceToOrigin;
    
    // Apply sign convention based on curvature
    if (!innerClockwise) {
        d0 = -d0;
    }
    
    debug() << "Final calculated d0 = " << d0 << " cm" << endmsg;
    debug() << "---- End Analytical Impact Parameter Calculation ----" << endmsg;
    
    return d0;
}
*/
void DisplacedTracking::fitLine(double x1, double y1, double x2, double y2, double x3, double y3,
                           double& slope, double& intercept) const{
    // Fit line using least squares method
    double sumx = x1 + x2 + x3;
    double sumy = y1 + y2 + y3;
    double sumxy = x1*y1 + x2*y2 + x3*y3;
    double sumx2 = x1*x1 + x2*x2 + x3*x3;
    
    // Number of points
    const int n = 3;
    
    // Calculate slope and intercept
    double denominator = n * sumx2 - sumx * sumx;
    if (std::abs(denominator) < 1e-6) {
        // Near-vertical line, use large slope
        slope = 1e6;
        intercept = sumy / n;
    } else {
        slope = (n * sumxy - sumx * sumy) / denominator;
        intercept = (sumy - slope * sumx) / n;
    }
}
// ------------------------------------ //
// -------------- GenFit -------------- //
// ------------------------------------ //
genfit::MeasuredStateOnPlane DisplacedTracking::convertToGenFitState(
    const edm4hep::TrackState& state,
    genfit::AbsTrackRep* rep) const {
    
    // ---------- Extract EDM4hep track parameters ----------
    // All EDM4hep track params: D0 [mm], phi [rad], omega [1/mm], Z0 [mm], tanLambda
    double d0        = state.D0;          // mm
    double phi       = state.phi;         // rad
    double omega     = state.omega;       // 1/mm, signed: omega = charge/R
    double tanLambda = state.tanLambda;

    // Reference point (hit position in mm → cm for GenFit)
    double refX = state.referencePoint[0] / 10.0;  // cm
    double refY = state.referencePoint[1] / 10.0;
    double refZ = state.referencePoint[2] / 10.0;

    // ---------- Magnetic field at reference point ----------
    dd4hep::Position fieldPos(refX, refY, refZ);  // cm
    double bField = m_field.magneticField(fieldPos).z() / dd4hep::tesla;  // Tesla

    if (std::abs(omega) < 1e-10 || std::abs(bField) < 1e-6) {
        throw genfit::Exception("convertToGenFitState: omega or B too small", __LINE__, __FILE__);
    }

    // ---------- pT and pz ----------
    // pT [GeV/c] = 0.3 * |Bz| [T] / (|omega| [1/mm] * 1000)
    double pT = 0.3 * std::abs(bField) / (std::abs(omega) * 1000.0);  // GeV/c
    double pz = pT * tanLambda;

    // ---------- Momentum direction at reference point ----------
    // Helix center (mm), derived from perigee parameters
    double R_mm = 1.0 / std::abs(omega);
    double centerX_mm, centerY_mm;
    if (omega < 0) {  // positive charge → clockwise
        centerX_mm = -d0 * std::sin(phi) - R_mm * std::cos(phi);
        centerY_mm =  d0 * std::cos(phi) - R_mm * std::sin(phi);
    } else {          // negative charge → counter-clockwise
        centerX_mm = -d0 * std::sin(phi) + R_mm * std::cos(phi);
        centerY_mm =  d0 * std::cos(phi) + R_mm * std::sin(phi);
    }

    // Unit radial vector from center to reference point (in cm)
    double vecX = refX - centerX_mm / 10.0;
    double vecY = refY - centerY_mm / 10.0;
    double vecMag = std::sqrt(vecX*vecX + vecY*vecY);
    if (vecMag < 1e-10) {
        throw genfit::Exception("convertToGenFitState: degenerate center-to-hit vector", __LINE__, __FILE__);
    }
    vecX /= vecMag;
    vecY /= vecMag;

    // Tangent (momentum direction in xy) = radial rotated 90°, sign by charge
    double momDirX, momDirY;
    if (omega < 0) {  // positive charge → clockwise
        momDirX = -vecY;  momDirY =  vecX;
    } else {          // negative charge → counter-clockwise
        momDirX =  vecY;  momDirY = -vecX;
    }

    double px = pT * momDirX;
    double py = pT * momDirY;

    TVector3 posVec(refX, refY, refZ);    // cm
    TVector3 momVec(px, py, pz);          // GeV/c

    // ---------- Build GenFit state with pos/mom only ----------
    // The covariance is set by the caller via genfit::Track(rep, stateVec, covInit),
    // which accepts a 6×6 Cartesian seed covariance.  Do not call setCov() here.
    genfit::MeasuredStateOnPlane gfState(rep);
    gfState.setPosMom(posVec, momVec);

    return gfState;
}


edm4hep::TrackState DisplacedTracking::convertToEDM4hepState(
    const genfit::MeasuredStateOnPlane& state,
    int location) const {
    
    // Get position [cm] and momentum [GeV/c] from GenFit state
    TVector3 pos = state.getPos();   // cm
    TVector3 mom = state.getMom();   // GeV/c
    double charge = (state.getQop() > 0) ? 1.0 : -1.0;
    
    double px = mom.X();
    double py = mom.Y();
    double pz = mom.Z();
    double pt = std::sqrt(px*px + py*py);
    double p  = mom.Mag();
    
    // ---------- EDM4hep 5-parameter track state ----------
    // phi: azimuthal angle of momentum at the track point
    double phi = std::atan2(py, px);

    // tanLambda: pz/pT
    double tanLambda = (pt > 1e-10) ? pz / pt : 0.0;

    // Magnetic field at fit position
    dd4hep::Position fieldPos(pos.X(), pos.Y(), pos.Z());  // cm
    double bField = m_field.magneticField(fieldPos).z() / dd4hep::tesla;
    
    // omega [1/mm] = charge * 0.3 * |B| / (pT [GeV/c] * 1000)
    double omega = 0.0;
    if (pt > 1e-10 && std::abs(bField) > 1e-6)
        omega = charge * 0.3 * std::abs(bField) / (pt * 1000.0);  // 1/mm

    // d0 [mm]: signed transverse impact parameter w.r.t. origin
    // d0 = -(y*cos(phi) - x*sin(phi)) evaluated at fit position (cm → mm)
    double d0 = (-pos.Y() * std::cos(phi) + pos.X() * std::sin(phi)) * 10.0;  // mm

    // z0 [mm]: z at closest approach
    double z0 = pos.Z() * 10.0;  // cm → mm

    edm4hep::TrackState outState;
    outState.D0        = d0;
    outState.phi       = phi;
    outState.omega     = omega;
    outState.Z0        = z0;
    outState.tanLambda = tanLambda;
    outState.location  = location;

    // Reference point: GenFit fit position converted mm
    outState.referencePoint = edm4hep::Vector3f(
        pos.X() * 10.0f,   // cm → mm
        pos.Y() * 10.0f,
        pos.Z() * 10.0f);

    // ---------- Covariance back-propagation 6D → 5D ----------
    // IMPORTANT: GenFit's MeasuredStateOnPlane::getCov() returns a 5×5 matrix in the
    // local track-representation basis [q/p, u', v', u, v].  It is NOT a 6×6 Cartesian
    // matrix.  Calling getCov() and treating the result as 6×6 causes the
    // "Request column(5)/row(5) outside matrix range of 0-5" errors because ROOT's
    // TMatrixDSym only has indices 0-4 for a 5×5 matrix.
    //
    // The correct approach is to use getPosMomCov() which explicitly constructs the
    // 6×6 Cartesian covariance [x,y,z,px,py,pz] from the internal representation.
    // We then propagate this back to EDM4hep's 5D parameter space via our Jacobian.
    {
        // --- Get full 6×6 Cartesian covariance from GenFit ---
        // getPosMomCov fills a 6×6 TMatrixDSym with the covariance in (x,y,z,px,py,pz)
        // ordering.  This is safe regardless of the internal representation dimension.
        TMatrixDSym C6sym(6);
        {
            TVector3 tmpPos, tmpMom;
            state.getPosMomCov(tmpPos, tmpMom, C6sym);  // fills 6×6 Cartesian cov (x,y,z,px,py,pz)
        }

        if (msgLevel(MSG::DEBUG)) {
            TVector3 fitPos = state.getPos();
            TVector3 fitMom = state.getMom();
            debug() << "  GenFit fitted pos (cm)    : (" << fitPos.X() << ", " << fitPos.Y() << ", " << fitPos.Z()
                    << ")  |p|=" << fitMom.Mag() << " GeV/c  pT=" << fitMom.Perp() << " GeV/c" << endmsg;
            debug() << "  GenFit fitted cov σ (C6)  :"
                    << "  x="  << std::sqrt(std::max(0., C6sym(0,0))) << " cm"
                    << "  y="  << std::sqrt(std::max(0., C6sym(1,1))) << " cm"
                    << "  z="  << std::sqrt(std::max(0., C6sym(2,2))) << " cm"
                    << "  px=" << std::sqrt(std::max(0., C6sym(3,3))) << " GeV/c"
                    << "  py=" << std::sqrt(std::max(0., C6sym(4,4))) << " GeV/c"
                    << "  pz=" << std::sqrt(std::max(0., C6sym(5,5))) << " GeV/c" << endmsg;
        }

        // Copy to plain TMatrixD for subsequent multiplications
        TMatrixD C6(6, 6);
        for (int i = 0; i < 6; ++i)
            for (int j = 0; j < 6; ++j)
                C6(i, j) = C6sym(i, j);

        // --- Jacobian ∂(x,y,z,px,py,pz)/∂(d0,phi,omega,z0,tanLambda) ---
        // GenFit stores position in cm and momentum in GeV/c.
        // EDM4hep stores d0,z0 in mm; omega in 1/mm; phi,tanLambda dimensionless.
        // We work in cm/GeV throughout and scale at the end.
        constexpr double mm2cm = 0.1;   // factor for mm → cm conversions
        double dpTdomega_inv = (std::abs(omega) > 1e-10) ? -pt / omega : 0.0;
        double momDirX = (pt > 1e-10) ? px / pt : 1.0;
        double momDirY = (pt > 1e-10) ? py / pt : 0.0;

        TMatrixD J(6, 5);
        J.Zero();
        // Row 0: ∂x/∂(d0,phi,omega,z0,tanL) — d0 in mm → multiply by mm2cm
        J(0, 0) = -std::sin(phi) * mm2cm;
        J(0, 1) = -d0 * std::cos(phi) * mm2cm;
        // Row 1: ∂y
        J(1, 0) =  std::cos(phi) * mm2cm;
        J(1, 1) = -d0 * std::sin(phi) * mm2cm;
        // Row 2: ∂z — z0 in mm
        J(2, 3) =  mm2cm;
        // Row 3: ∂px
        J(3, 1) = -py;
        J(3, 2) =  momDirX * dpTdomega_inv;
        // Row 4: ∂py
        J(4, 1) =  px;
        J(4, 2) =  momDirY * dpTdomega_inv;
        // Row 5: ∂pz
        J(5, 2) =  tanLambda * dpTdomega_inv;
        J(5, 4) =  pt;

        TMatrixD JT(TMatrixD::kTransposed, J);    // 5×6

    // Pseudo-inverse: C5 = (J^T J)^{-1} J^T C6 J (J^T J)^{-1}
    // All steps use explicit Mult() to avoid ROOT operator* dimension errors
    TMatrixD JtJ(5, 5);
    JtJ.Mult(JT, J);                           // (5×6)·(6×5) = 5×5
    TMatrixD JtJ_inv(TMatrixD::kInverted, JtJ);

    TMatrixD JtC6(5, 6);
    JtC6.Mult(JT, C6);                         // (5×6)·(6×6) = 5×6
    TMatrixD JtC6J(5, 5);
    JtC6J.Mult(JtC6, J);                       // (5×6)·(6×5) = 5×5

    TMatrixD tmpL(5, 5);
    tmpL.Mult(JtJ_inv, JtC6J);                 // 5×5
    TMatrixD C5mat(5, 5);
    C5mat.Mult(tmpL, JtJ_inv);                 // 5×5

    // Fill EDM4hep covariance; D0(col 0) and Z0(col 3) are in mm→ scale ×100 per axis
    for (int i = 0; i < 5; ++i) {
        for (int j = 0; j <= i; ++j) {
            double v = 0.5 * (C5mat(i,j) + C5mat(j,i));
            auto scaleRow = [](int k) { return (k == 0 || k == 3) ? 100.0 : 1.0; };
            v *= scaleRow(i) * scaleRow(j);
            outState.setCovMatrix(v,
                static_cast<edm4hep::TrackParams>(i),
                static_cast<edm4hep::TrackParams>(j));
        }
    }

    // Debug: print EDM4hep 5×5 covariance diagonal
    if (msgLevel(MSG::DEBUG)) {
        double v_d0    = outState.getCovMatrix(edm4hep::TrackParams::d0,        edm4hep::TrackParams::d0);
        double v_phi   = outState.getCovMatrix(edm4hep::TrackParams::phi,       edm4hep::TrackParams::phi);
        double v_omega = outState.getCovMatrix(edm4hep::TrackParams::omega,     edm4hep::TrackParams::omega);
        double v_z0    = outState.getCovMatrix(edm4hep::TrackParams::z0,        edm4hep::TrackParams::z0);
        double v_tanL  = outState.getCovMatrix(edm4hep::TrackParams::tanLambda, edm4hep::TrackParams::tanLambda);
        debug() << "  EDM4hep output cov σ (C5) :"
                << "  d0="    << std::sqrt(std::max(0., v_d0))    << " mm"
                << "  phi="   << std::sqrt(std::max(0., v_phi))   << " rad"
                << "  ω="     << std::sqrt(std::max(0., v_omega)) << " 1/mm"
                << "  z0="    << std::sqrt(std::max(0., v_z0))    << " mm"
                << "  tanL="  << std::sqrt(std::max(0., v_tanL))  << endmsg;
    }

    } // end covariance block

    debug() << "  EDM4hep params    :"
            << "  d0="     << d0      << " mm"
            << "  phi="    << phi     << " rad"
            << "  ω="      << omega   << " 1/mm"
            << "  z0="     << z0      << " mm"
            << "  tanL="   << tanLambda << endmsg;
    
    return outState;
}

genfit::AbsMeasurement* DisplacedTracking::createGenFitMeasurement(
    const edm4hep::TrackerHitPlane& hit,
    const dd4hep::rec::Surface* surface,
    int hitId,
    genfit::TrackPoint* trackPoint) const {
    
    // Get surface properties for measurement plane
    dd4hep::rec::Vector3D u = surface->u();
    dd4hep::rec::Vector3D v = surface->v();
    dd4hep::rec::Vector3D origin = surface->origin();
    
    // Convert to ROOT's TVector3
    TVector3 o(origin.x(), origin.y(), origin.z());  // Note: origin already in cm
    TVector3 uVec(u.x(), u.y(), u.z());
    TVector3 vVec(v.x(), v.y(), v.z());
    
    // Get hit position
    const auto& pos = hit.getPosition();  // mm
    
    // Convert global position to local position on the surface
    dd4hep::rec::Vector3D hitGlobal(pos[0], pos[1], pos[2]);  // mm
    dd4hep::rec::Vector2D hitLocal = surface->globalToLocal(dd4hep::mm * hitGlobal);  // cm
    
    // Create measurement coordinates
    TVectorD hitCoords(2);
    hitCoords[0] = hitLocal[0];  // u coordinate in cm
    hitCoords[1] = hitLocal[1];  // v coordinate in cm
    
    // Get measurement errors
    double sigma_u = hit.getDu();  // mm
    double sigma_v = hit.getDv();  // mm
    
    // If resolutions are zero or not set, use reasonable defaults
    if (sigma_u <= 0) sigma_u = 0.4;  // 400 microns default
    if (sigma_v <= 0) sigma_v = 0.4;  // 400 microns default
    
    // Create covariance matrix
    TMatrixDSym hitCov(2);
    hitCov(0,0) = std::pow(dd4hep::mm * sigma_u, 2);  // cm²
    hitCov(0,1) = 0;
    hitCov(1,0) = 0; 
    hitCov(1,1) = std::pow(dd4hep::mm * sigma_v, 2);  // cm²
    
    // Extract detector ID from cell ID
    int detId = hit.getCellID();
    
    // Create planar measurement
    genfit::PlanarMeasurement* meas = 
        new genfit::PlanarMeasurement(hitCoords, hitCov, detId, hitId, trackPoint);
    
    // Create and set measurement plane
    genfit::SharedPlanePtr plane(new genfit::DetPlane(o, uVec, vVec));
    meas->setPlane(plane, detId);  // Use cellID as planeID as in the example
    
    return meas;
}

bool DisplacedTracking::fitTrackWithGenFit(
    const std::vector<edm4hep::TrackerHitPlane>& hits,
    const edm4hep::TrackState& seedState,
    edm4hep::MutableTrack& finalTrack,
    edm4hep::TrackerHitPlaneCollection& outputHits,
    const LinkPropagator& propagateLink,
    const edm4hep::TrackerHitPlaneCollection* allHits) const {
    
    if (hits.empty()) {
        warning() << "Cannot fit track - no hits provided" << endmsg;
        return false;
    }
    
    // Determine particle type for the track representation
    int pdgCode = getPDGCode();
    
    debug() << "Using particle type: " << m_particleType << " (PDG code: " << pdgCode << ")" << endmsg;
    
    try {
        // Create a track representation
        genfit::AbsTrackRep* rep = new genfit::RKTrackRep(pdgCode);
        //genfit::AbsTrackRep* rep = new genfit::GeaneTrackRep();  // takes more time, but more efficient with material and particle independent.
        // Convert seed state (EDM4hep helix params) → Cartesian (pos, mom) for GenFit
        genfit::MeasuredStateOnPlane seedGFState = convertToGenFitState(seedState, rep);
        TVector3 pos = seedGFState.getPos();  // cm
        TVector3 mom = seedGFState.getMom();  // GeV/c

        // Build 6D state vector [x, y, z, px, py, pz]
        TVectorD stateVec(6);
        stateVec[0] = pos.X();  stateVec[1] = pos.Y();  stateVec[2] = pos.Z();
        stateVec[3] = mom.X();  stateVec[4] = mom.Y();  stateVec[5] = mom.Z();

        // ── Build 6×6 Cartesian seed covariance from seed state's EDM4hep covariance ──
        // We propagate the 5D EDM4hep cov C5 to Cartesian via the analytic Jacobian
        //   J = ∂(x,y,z,px,py,pz)/∂(d0,phi,ω,z0,tanλ)
        // so that covInit = J · C5 · J^T.
        // This replaces the old hardcoded (10% of |p|) momentum sigma.
        TMatrixDSym covInit(6);
        covInit.Zero();
        {
            // Read EDM4hep diagonal variances (units: mm², rad², (1/mm)², mm², dimensionless)
            double var_d0   = seedState.getCovMatrix(edm4hep::TrackParams::d0,        edm4hep::TrackParams::d0);
            double var_phi  = seedState.getCovMatrix(edm4hep::TrackParams::phi,       edm4hep::TrackParams::phi);
            double var_om   = seedState.getCovMatrix(edm4hep::TrackParams::omega,     edm4hep::TrackParams::omega);
            double var_z0   = seedState.getCovMatrix(edm4hep::TrackParams::z0,        edm4hep::TrackParams::z0);
            double var_tanL = seedState.getCovMatrix(edm4hep::TrackParams::tanLambda, edm4hep::TrackParams::tanLambda);

            // Apply floors (protect against unset or collapsed covariances)
            if (var_d0   < 1e-30) var_d0   = 1.0;  // (1 mm)²
            if (var_phi  < 1e-30) var_phi  = 0.1;
            //if (var_om   < 1e-30) var_om   = 1e-10;  // absolute minimum
            if (var_z0   < 1e-30) var_z0   = 1.0;
            if (var_tanL < 1e-30) var_tanL = 0.1;

            // Track parameters at ref point (mm, rad, 1/mm)
            double phi  = seedState.phi;
            double om   = seedState.omega;   // 1/mm, signed
            double tanL = seedState.tanLambda;
            double pT   = mom.Perp();        // GeV/c  (from Cartesian conversion)
            double pMag = mom.Mag();

            // ∂(px,py)/∂phi: tangent direction rotates with phi
            //   px = pT * sin(phi_track);  phi_track depends on helix geometry
            //   Proxy: use the actual Cartesian (px,py) directions
            double ux = (pMag > 1e-10) ? mom.X() / pMag : 0.0;
            double uy = (pMag > 1e-10) ? mom.Y() / pMag : 0.0;
            double uz = (pMag > 1e-10) ? mom.Z() / pMag : 0.0;

            // Position block (mm² → cm²):
            //   σ²(x) ≈ σ²(y) ≈ σ²(d0) * mm²→cm²
            //   σ²(z) ≈ σ²(z0) * mm²→cm²
            const double mm2cm2 = 0.01;
            double var_xy = std::max(var_d0 * mm2cm2, 0.1);   // floor: 0.1 cm²
            double var_z  = std::max(var_z0 * mm2cm2, 0.1);
            covInit(0,0) = var_xy;
            covInit(1,1) = var_xy;
            covInit(2,2) = var_z;

            // Momentum block via Jacobian:
            //   pT  = 0.3 * |Bz| / (|ω| * 1000)  → ∂pT/∂ω = -pT/ω
            //   pz  = pT * tanL
            //   σ²(pT) from σ²(ω):  var_pT = (pT/ω)² * var_om   [GeV² per (1/mm)²]
            //   σ²(pz) from σ²(tanL): var_pz = pT² * var_tanL
            //   σ²(phi) contributes rotation of (px,py) but not |pT|
            double dpT_dom  = (std::abs(om) > 1e-10) ? -pT / om : 0.0;  // GeV·mm
            double var_pT   = dpT_dom * dpT_dom * var_om;      // (GeV/c)²
            double var_pz   = pT * pT * var_tanL;              // (GeV/c)²
            // phi variance rotates (px,py) in the transverse plane:
            //   px = pT*cos(φ_track), py = pT*sin(φ_track) → ∂px/∂phi ≈ -py, ∂py/∂phi ≈ px
            double var_px_phi = mom.Y() * mom.Y() * var_phi;   // (GeV/c)²
            double var_py_phi = mom.X() * mom.X() * var_phi;

            // Combine pT magnitude uncertainty + angular uncertainty, then apply floor
            double pT_floor = std::max(pT * 0.1, 0.1);  // 1% of pT or 0.1 GeV/c
            double var_pT_floor = pT_floor * pT_floor;

            covInit(3,3) = std::max(var_pT * ux*ux + var_px_phi, var_pT_floor);
            covInit(4,4) = std::max(var_pT * uy*uy + var_py_phi, var_pT_floor);
            covInit(5,5) = std::max(var_pz,                       var_pT_floor);
        }

        debug() << "─── GenFit Input ──────────────────────────────────────────────────────" << endmsg;
        debug() << "  Seed state source : "
                << (seedState.location == edm4hep::TrackState::AtLastHit ? "AtLastHit (4-hit circle fit)" :
                    seedState.location == edm4hep::TrackState::AtFirstHit ? "AtFirstHit (3-hit seed)" : "other")
                << endmsg;
        debug() << "  Seed pos  (cm)    : (" << pos.X() << ", " << pos.Y() << ", " << pos.Z() << ")" << endmsg;
        debug() << "  Seed mom  (GeV/c) : (" << mom.X() << ", " << mom.Y() << ", " << mom.Z()
                << ")  |p|=" << mom.Mag() << "  pT=" << mom.Perp() << endmsg;
        debug() << "  Seed cov diag (σ) :"
                << "  x="   << std::sqrt(covInit(0,0)) << " cm"
                << "  y="   << std::sqrt(covInit(1,1)) << " cm"
                << "  z="   << std::sqrt(covInit(2,2)) << " cm"
                << "  px="  << std::sqrt(covInit(3,3)) << " GeV/c"
                << "  py="  << std::sqrt(covInit(4,4)) << " GeV/c"
                << "  pz="  << std::sqrt(covInit(5,5)) << " GeV/c" << endmsg;
        debug() << "  Hits to fit       : " << hits.size() << endmsg;
        debug() << "───────────────────────────────────────────────────────────────────────" << endmsg;

        genfit::Track finalGFTrack(rep, stateVec, covInit);
        
        // Add all seed hits to the track
        for (size_t i = 0; i < hits.size(); ++i) {
            const auto& hit = hits[i];
            
            // Find the surface for this hit
            const dd4hep::rec::Surface* surface = findSurface(hit);
            if (!surface) {
                warning() << "Could not find surface for hit " << i << endmsg;
                continue;
            }
            
            // Create a new track point
            genfit::TrackPoint* trackPoint = new genfit::TrackPoint();
            
            // Store the hit index if it's from the allHits collection
            int hitIdx = -1;
            if (allHits) {
                for (size_t j = 0; j < allHits->size(); ++j) {
                    if ((*allHits)[j].getCellID() == hit.getCellID()) {
                        hitIdx = j;
                        break;
                    }
                }
            }
            
            // Create a measurement for this hit
            genfit::AbsMeasurement* measurement = createGenFitMeasurement(hit, surface, hitIdx, trackPoint);
            if (!measurement) {
                warning() << "Could not create measurement for hit " << i << endmsg;
                delete trackPoint;
                continue;
            }
            
            // Add measurement to track point
            trackPoint->addRawMeasurement(measurement);
            
            // Add track point to track
            finalGFTrack.insertPoint(trackPoint);
        }
        
        //------------------------------------------------------------------
        // Fit the track with GenFit (KalmanFitterRefTrack)
        //------------------------------------------------------------------
        genfit::KalmanFitterRefTrack fitter;
        fitter.setMaxIterations(m_maxFitIterations);
        fitter.setMinIterations(3);

        debug() << "─── GenFit Kalman Fit ─────────────────────────────────────────────────" << endmsg;
        debug() << "  Fitter            : KalmanFitterRefTrack"
                << "  max_iter=" << m_maxFitIterations << "  min_iter=3" << endmsg;
        debug() << "  Points in track   : " << finalGFTrack.getNumPoints() << endmsg;
        std::string capturedStderr;
        {
            // Save original stderr fd
            int savedStderr = dup(STDERR_FILENO);
            // Create a pipe to capture stderr
            int pipefd[2];
            if (pipe(pipefd) == 0) {
                dup2(pipefd[1], STDERR_FILENO);
                close(pipefd[1]);

                fitter.processTrack(&finalGFTrack);

                // Restore stderr — flush both C and C++ stderr before restoring
                fflush(stderr);
                std::cerr.flush();
                dup2(savedStderr, STDERR_FILENO);
                close(savedStderr);

                // Read captured output
                // Set non-blocking and drain the pipe
                fcntl(pipefd[0], F_SETFL, O_NONBLOCK);
                char buf[4096];
                ssize_t n;
                while ((n = read(pipefd[0], buf, sizeof(buf)-1)) > 0) {
                    buf[n] = '\0';
                    capturedStderr += buf;
                }
                close(pipefd[0]);

                // If anything was captured, print it at warning level
                if (!capturedStderr.empty()) {
                    warning() << "GenFit stderr during fit: " << capturedStderr << endmsg;
                }
            } else {
                // pipe() failed — just run without capture
                close(savedStderr);
                fitter.processTrack(&finalGFTrack);
            }
        }

        // ---- Backward fit is commented out for now ----
        // The forward fit is sufficient for the muon system seed quality.
        // A proper bidirectional smooth requires combining forward+backward states
        // at each measurement site (smoothing step), which GenFit's KalmanFitterRefTrack
        // supports via processTrackWithRep with bidir=true.  Enable when convergence is stable.
        //
        // genfit::Track backwardTrack = finalGFTrack;
        // backwardTrack.reverseTrack();
        // fitter.processTrack(&backwardTrack);

        // Check fit quality
        // processTrack() swallows ill-conditioned exceptions internally and sets isFitted=false.
        // We detect Cholesky failures by capturing stderr during processTrack (GenFit prints to stderr).
        if (!finalGFTrack.getFitStatus()->isFitted()) {
            // Classify: ill-conditioned covariance vs other failure.
            // The stderr output from processTrack was already captured in capturedStderr above.
            if (capturedStderr.find("ill-conditioned") != std::string::npos ||
                capturedStderr.find("not positive definite") != std::string::npos ||
                capturedStderr.find("TDecompChol") != std::string::npos) {
                warning() << "GenFit fit failed: ill-conditioned covariance (Cholesky failure)" << endmsg;
                m_statGenFitIllCond++;
            } else {
                warning() << "GenFit track fitting failed (isFitted=false)" << endmsg;
                m_statGenFitFailed++;
            }
            return false;
        }
        
        // Get fit quality metrics
        double chi2 = finalGFTrack.getFitStatus()->getChi2();
        int ndf     = finalGFTrack.getFitStatus()->getNdf();
        double chi2ndf = (ndf > 0) ? chi2 / ndf : -1.0;

        debug() << "─── GenFit Result ─────────────────────────────────────────────────────" << endmsg;
        debug() << "  Fit status        : " << (finalGFTrack.getFitStatus()->isFitted() ? "FITTED" : "FAILED") << endmsg;
        debug() << "  chi2=" << std::fixed << std::setprecision(3) << chi2
                << "  ndf=" << ndf
                << "  chi2/ndf=" << chi2ndf << endmsg;

        // Accumulate chi2/ndf for run statistics (store as integer × 1000 to keep atomic)
        // Cap at chi2/ndf < 1000 to exclude numerically blown-up fits that are still
        // flagged isFitted=true but have non-physical covariances.
        // NOTE: GenFit chi2 is purely from position measurement residuals, NOT momentum.
        // A huge chi2/ndf indicates the track model doesn't describe the hits well.
        if (ndf > 0) {
            if (chi2ndf < 1000.0) {
                m_statGenFitChi2Count++;
                m_statGenFitChi2Sum += static_cast<long long>(chi2ndf * 1000.0);
                debug() << "  Running avg chi2/ndf=" << std::fixed << std::setprecision(3)
                        << (m_statGenFitChi2Sum.load() / 1000.0 / m_statGenFitChi2Count.load())
                        << "  (over " << m_statGenFitChi2Count.load() << " tracks)" << endmsg;
            } else {
                warning() << "  Non-physical chi2/ndf=" << std::fixed << std::setprecision(1) << chi2ndf
                          << "  (chi2=" << chi2 << "  ndf=" << ndf << ") — excluded from average" << endmsg;
                m_statGenFitChi2BadCount++;
            }
        }

        // ---- Update chi2/ndf on the EDM4hep track ----
        finalTrack.setChi2(chi2);
        finalTrack.setNdf(ndf);

        // Copy hits into output collection (only if not already filled by the seed fallback path)
        if (finalTrack.trackerHits_size() == 0) {
            for (const auto& hit : hits) {
                auto outHit = outputHits.create();
                outHit.setCellID(hit.getCellID()); outHit.setTime(hit.getTime());
                outHit.setEDep(hit.getEDep()); outHit.setEDepError(hit.getEDepError());
                outHit.setPosition(hit.getPosition()); outHit.setCovMatrix(hit.getCovMatrix());
                outHit.setDu(hit.getDu()); outHit.setDv(hit.getDv());
                propagateLink(outHit, hit);
                finalTrack.addToTrackerHits(outHit);
            }
        }

        // ----------------------------------------------------------------
        // Save GenFit fitted states.
        // Convention:
        //   • AtFirstHit  → fitted state at first measurement (innermost hit in outer muon system)
        //   • AtLastHit   → fitted state at last  measurement (used only if no 4-hit circle fit)
        //   • AtOther     → GenFit fitted state at first hit (labelled separately so we can compare
        //                    with the analytical seed which is also stored at AtFirstHit)
        //
        // We do NOT overwrite any state that was already added by the seeding / 4-hit circle fit
        // (those are preserved for comparison).  Instead we use AtOther for the GenFit output.
        // Inner-propagation result will be saved at AtVertex (see propagation section below).
        // ----------------------------------------------------------------

        // Determine which locations are already occupied in this track
        auto hasLocation = [&](int loc) {
            for (int j = 0; j < finalTrack.trackStates_size(); ++j)
                if (finalTrack.getTrackStates(j).location == loc) return true;
            return false;
        };

        // ---- State at first measurement ----
        if (finalGFTrack.getNumPoints() > 0) {
            try {
                genfit::MeasuredStateOnPlane stateFirst = finalGFTrack.getFittedState(0);
                // Save as AtOther so it does not overwrite the analytical AtFirstHit seed
                edm4hep::TrackState gfFirst = convertToEDM4hepState(stateFirst,
                                                    edm4hep::TrackState::AtOther);
                // Only add if AtOther is not yet used (inner propagation uses it later)
                if (!hasLocation(edm4hep::TrackState::AtOther)) {
                    finalTrack.addToTrackStates(gfFirst);
                    debug() << "  Saved fitted state : AtOther (first measurement)" << endmsg;
                } else {
                    finalTrack.addToTrackStates(gfFirst);
                    debug() << "  Saved fitted state : AtOther (additional)" << endmsg;
                }
            } catch (genfit::Exception& e) {
                warning() << "Could not extract fitted state at first hit: " << e.what() << endmsg;
            }
        }

        // ---- State at last measurement ----
        // AtLastHit is RESERVED for the 4-hit analytical circle fit (set earlier in findTracks).
        // GenFit's fitted state at the last measurement is NOT stored here to keep the semantics
        // clean: AtFirstHit = analytical seed, AtLastHit = 4-hit circle, AtOther = GenFit at first hit.
        // For 3-hit tracks the last-measurement GenFit state is redundant (very close to AtOther).
        // Only add it if the 4-hit circle fit already filled AtLastHit (i.e. it's a 4-hit track)
        // in which case we add the GenFit last-hit state at AtOther2 — but since EDM4hep has no
        // second AtOther slot, we simply skip it. Summary: for 3-hit tracks, no AtLastHit is stored.
        if (finalGFTrack.getNumPoints() > 1 && hasLocation(edm4hep::TrackState::AtLastHit)) {
            // 4-hit track: AtLastHit already has the circle fit; optionally compare with GenFit last.
            // Currently skipped — the circle fit is the preferred reference.
            debug() << "Skipping GenFit AtLastHit add (4-hit circle fit already occupies AtLastHit)" << endmsg;
        }

        debug() << "  Track states total : " << finalTrack.trackStates_size()
                << "  hits: " << finalTrack.trackerHits_size() << endmsg;
        debug() << "───────────────────────────────────────────────────────────────────────" << endmsg;
        
        return true;
        
    } catch (genfit::Exception& e) {
        std::string exMsg = e.what();
        // Distinguish ill-conditioned covariance (Cholesky failure) from other errors
        if (exMsg.find("ill-conditioned") != std::string::npos ||
            exMsg.find("not positive definite") != std::string::npos) {
            warning() << "GenFit: ill-conditioned covariance matrix (numerics): " << exMsg << endmsg;
            m_statGenFitIllCond++;
        } else {
            error() << "GenFit exception during track fitting: " << exMsg << endmsg;
            m_statGenFitException++;
        }
        return false;
    }
}

// Extend a track by searching for compatible extra hits in subsequent layers.
// Starts from the current last hit and loops until no more hits are found.
// 
// Layer acceptance: uses compositeID (barrel=layerID, +endcap=1000+layerID,
//   -endcap=-1000-layerID) so barrel and endcap layers are never conflated.
//   Same-composite-layer hits are NOT accepted (disabled — needs tighter spatial
//   cut to avoid fake combinatorics; TODO: to fix ..re-enable with ~1 cm threshold).
//   Cross-region transitions (e.g. barrel seed -> endcap continuation) are handled
//   naturally by the cone cut — no compositeID arithmetic needed across regions.
//
// Hit acceptance uses a forward cone from the seed direction:
//   - the delta vector from the last hit to the candidate must point within
//     coneHalfAngleDeg of seedDirection (rejects hits behind or off-axis)
//   - distance must be within m_maxDist
bool DisplacedTracking::findCompatibleExtraHit(
    std::vector<edm4hep::TrackerHitPlane>& trackHits,
    const edm4hep::TrackerHitPlaneCollection* allHits,
    std::vector<bool>& usedHits,
    const Eigen::Vector3d& seedDirection,
    double coneHalfAngleDeg) const {

    if (trackHits.size() < 3) {
        debug() << "Need at least 3 hits for extension, got " << trackHits.size() << endmsg;
        return false;
    }

    const double cosHalfAngle = std::cos(coneHalfAngleDeg * M_PI / 180.0);
    // Normalize seed direction once
    const Eigen::Vector3d dir = seedDirection.normalized();

    bool anyFound = false;

    // Loop: each iteration tries to add one more hit beyond the current last hit.
    // Stops when no compatible hit is found for the current last layer.
    while (true) {
        const auto& lastHit      = trackHits.back();
        const auto& lastPos      = lastHit.getPosition();
        Eigen::Vector3d lastHitPos(lastPos[0] / 10.0, lastPos[1] / 10.0, lastPos[2] / 10.0); // mm->cm

        int lastCompositeID = getCompositeID(lastHit.getCellID());

        debug() << "Extra-hit search loop: current last hit compositeID=" << lastCompositeID
                << " at (" << lastHitPos.x() << "," << lastHitPos.y() << "," << lastHitPos.z() << ") cm"
                << " | track now has " << trackHits.size() << " hits" << endmsg;

        double bestDistance  = 1e6;
        int    bestHitIndex  = -1;
        edm4hep::TrackerHitPlane bestHit;

        for (size_t i = 0; i < allHits->size(); ++i) {
            if (usedHits[i]) continue;

            const auto& candHit      = (*allHits)[i];
            int candCompositeID      = getCompositeID(candHit.getCellID());

            // --- Layer gate ---
            // Reject hits in the same composite layer (see comment above).
            if (candCompositeID == lastCompositeID) continue;

            // For same-region (barrel->barrel or endcap->endcap) we require the
            // composite layer number to be strictly greater (outward direction).
            // For cross-region transitions we rely purely on the cone cut below,
            // so we only hard-reject hits that are clearly inward in the same region.
            int lastType    = getTypeID(lastHit.getCellID());
            int candType    = getTypeID(candHit.getCellID());
            int lastLayerN  = getLayerID(lastHit.getCellID());
            int candLayerN  = getLayerID(candHit.getCellID());
            if (lastType == candType) {
                // Same region: require outward step in physical layer number.
                // For barrel and +endcap, layerID increases outward.
                // For -endcap, layerID also increases outward (layer 0 is closest to IP),
                // so comparing raw layerIDs is correct for all three regions.
                if (candLayerN <= lastLayerN) continue;
                // Allow at most a skip of 2 physical layers (no wild jumps)
                if (candLayerN > lastLayerN + 2) continue;
            }
            // Cross-region: no compositeID arithmetic constraint — cone decides

            // --- Cone gate ---
            const auto& candPos = candHit.getPosition();
            Eigen::Vector3d candHitPos(candPos[0] / 10.0, candPos[1] / 10.0, candPos[2] / 10.0);

            Eigen::Vector3d delta = candHitPos - lastHitPos;
            double distance = delta.norm();
            if (distance < 1e-6) continue; // degenerate

            double cosTheta = delta.dot(dir) / distance;

            // Must be within the forward cone and within max distance
            if (cosTheta < cosHalfAngle) {
                debug() << "Candidate hit " << i << " compositeID=" << candCompositeID
                        << " rejected by cone (cosTheta=" << cosTheta
                        << " < " << cosHalfAngle << ")" << endmsg;
                continue;
            }
            if (distance > m_maxDist.value()) continue;

            debug() << "Candidate hit " << i << " compositeID=" << candCompositeID
                    << " distance=" << distance << " cm cosTheta=" << cosTheta << endmsg;

            if (distance < bestDistance) {
                bestDistance = distance;
                bestHitIndex = static_cast<int>(i);
                bestHit      = candHit;
            }
        }

        if (bestHitIndex < 0) {
            debug() << "No compatible extra hit found beyond compositeID=" << lastCompositeID
                    << " — stopping extension loop" << endmsg;
            break;
        }

        int bestCompositeID = getCompositeID(bestHit.getCellID());
        info() << "Found extra hit: compositeID=" << bestCompositeID
               << " distance=" << bestDistance << " cm"
               << " (track now has " << trackHits.size() + 1 << " hits)" << endmsg;

        trackHits.push_back(bestHit);
        usedHits[bestHitIndex] = true;
        anyFound = true;
    }

    return anyFound;
}

bool DisplacedTracking::fitCircleToFourHits(
    const std::vector<edm4hep::TrackerHitPlane>& hits,
    double& x0, double& y0, double& radius, double& chi2,
    Eigen::Matrix3d& fitCov3x3) const {
    
    if (hits.size() != 4) {
        warning() << "fitCircleToFourHits expects exactly 4 hits, got " << hits.size() << endmsg;
        return false;
    }
    
    // Extract hit positions in cm
    std::vector<Eigen::Vector2d> points;
    std::vector<double> uncertainties;
    
    for (const auto& hit : hits) {
        const auto& pos = hit.getPosition();
        points.emplace_back(pos[0] / 10.0, pos[1] / 10.0); // Convert mm to cm
        
        // Get measurement uncertainties
        double sigma_u = hit.getDu(); // mm
        double sigma_v = hit.getDv(); // mm
        
        // Use defaults if not set
        if (sigma_u <= 0) sigma_u = 0.4; // 400 microns
        if (sigma_v <= 0) sigma_v = 0.4; // 400 microns
        
        // Use average uncertainty (could be improved with proper error propagation)
        double avgSigma = std::sqrt(sigma_u * sigma_u + sigma_v * sigma_v) / 10.0; // Convert to cm
        uncertainties.push_back(avgSigma);
    }
    
    debug() << "Fitting circle to 4 points:" << endmsg;
    for (size_t i = 0; i < points.size(); ++i) {
        debug() << "  Point " << i << ": (" << points[i].x() << ", " << points[i].y() 
                << ") cm, σ = " << uncertainties[i] << " cm" << endmsg;
    }
    
    // Method 1: Algebraic circle fit using least squares
    // Circle equation: x² + y² + Dx + Ey + F = 0
    // Where D = -2x0, E = -2y0, F = x0² + y0² - R²
    
    Eigen::MatrixXd A(4, 3);
    Eigen::VectorXd b(4);
    Eigen::VectorXd weights(4);
    
    // Setup the linear system
    for (int i = 0; i < 4; ++i) {
        double x = points[i].x();
        double y = points[i].y();
        double w = 1.0 / (uncertainties[i] * uncertainties[i]); // Weight = 1/σ²
        
        A(i, 0) = w * x;           // D coefficient
        A(i, 1) = w * y;           // E coefficient  
        A(i, 2) = w * 1.0;         // F coefficient
        b(i) = -w * (x*x + y*y);   // RHS
        weights(i) = w;
    }
    
    // Solve the weighted least squares system: A^T W A x = A^T W b
    Eigen::Vector3d solution = (A.transpose() * A).ldlt().solve(A.transpose() * b);
    
    double D = solution(0);
    double E = solution(1); 
    double F = solution(2);
    
    // Convert back to center and radius
    x0 = -D / 2.0;
    y0 = -E / 2.0;
    
    double discriminant = D*D + E*E - 4*F;
    if (discriminant <= 0) {
        debug() << "Invalid circle fit: discriminant = " << discriminant << endmsg;
        return false;
    }
    
    radius = std::sqrt(discriminant) / 2.0;
    
    debug() << "Algebraic fit result: center=(" << x0 << ", " << y0 
            << "), radius=" << radius << " cm" << endmsg;
    
    // Method 2: Geometric refinement using Gauss-Newton
    // This improves the fit by minimizing geometric distance rather than algebraic distance
    
    // Initial guess from algebraic fit
    Eigen::Vector3d params(x0, y0, radius);
    
    // Gauss-Newton iterations
    const int maxIterations = 10;
    const double tolerance = 1e-6;
    
    for (int iter = 0; iter < maxIterations; ++iter) {
        Eigen::MatrixXd J(4, 3); // Jacobian matrix
        Eigen::VectorXd residuals(4);
        
        double currentX0 = params(0);
        double currentY0 = params(1);
        double currentR = params(2);
        
        // Calculate residuals and Jacobian
        for (int i = 0; i < 4; ++i) {
            double x = points[i].x();
            double y = points[i].y();
            double w = std::sqrt(weights(i));
            
            // Distance from point to circle center
            double dx = x - currentX0;
            double dy = y - currentY0;
            double dist = std::sqrt(dx*dx + dy*dy);
            
            // Residual: weighted geometric distance from point to circle
            residuals(i) = w * (dist - currentR);
            
            // Jacobian elements (derivatives of residual w.r.t. parameters)
            if (dist > 1e-12) { // Avoid division by zero
                J(i, 0) = -w * dx / dist;  // ∂r/∂x0
                J(i, 1) = -w * dy / dist;  // ∂r/∂y0
                J(i, 2) = -w;              // ∂r/∂R
            } else {
                J(i, 0) = J(i, 1) = J(i, 2) = 0;
            }
        }
        
        // Gauss-Newton update: Δp = -(J^T J)^(-1) J^T r
        Eigen::Vector3d delta = -(J.transpose() * J).ldlt().solve(J.transpose() * residuals);
        params += delta;
        
        // Check convergence
        if (delta.norm() < tolerance) {
            debug() << "Geometric refinement converged after " << iter + 1 << " iterations" << endmsg;
            break;
        }
    }
    
    // Update final parameters
    x0 = params(0);
    y0 = params(1);
    radius = params(2);

    // ── Compute fit covariance on [x0, y0, R] from the final Gauss-Newton Jacobian ──
    // Cov(params) = (J^T W J)^{-1}  where W=diag(w_i) from the WLS weights.
    // Build the final J at the converged parameters.
    {
        Eigen::MatrixXd Jf(4, 3);
        for (int i = 0; i < 4; ++i) {
            double xi = points[i].x(), yi = points[i].y();
            double dx = xi - x0, dy = yi - y0;
            double dist = std::sqrt(dx*dx + dy*dy);
            double w = std::sqrt(weights(i));
            if (dist > 1e-12) {
                Jf(i,0) = -w * dx / dist;
                Jf(i,1) = -w * dy / dist;
                Jf(i,2) = -w;
            } else {
                Jf(i,0) = Jf(i,1) = Jf(i,2) = 0;
            }
        }
        Eigen::Matrix3d JtJ = Jf.transpose() * Jf;
        // Use pseudo-inverse via LDLT (should be full-rank for a good fit)
        fitCov3x3 = JtJ.ldlt().solve(Eigen::Matrix3d::Identity());
    }

    dd4hep::Position fieldPos(0, 0, 0);  // Default center position      
    // Use average position of hits for better field estimate
    double sumX = 0, sumY = 0, sumZ = 0;
    for (const auto& hit : hits) {
        sumX += hit.getPosition()[0] / 10.0;  // mm to cm
        sumY += hit.getPosition()[1] / 10.0;
        sumZ += hit.getPosition()[2] / 10.0;
    }
    fieldPos = dd4hep::Position(
                    sumX / hits.size(),
                    sumY / hits.size(),
                    sumZ / hits.size());

    double actualBz = m_field.magneticField(fieldPos).z() / dd4hep::tesla;
    double pT = 0.3 * std::abs(actualBz) * radius / 100.0;

    debug() << "Refined fit result: center=(" << x0 << ", " << y0 
            << "), radius=" << radius << " cm, B field=" << actualBz << " ,pT=" << pT << " GeV" << endmsg;
    
    // Calculate chi2
    chi2 = 0.0;
    for (int i = 0; i < 4; ++i) {
        double x = points[i].x();
        double y = points[i].y();
        double sigma = uncertainties[i];
        
        // Calculate distance from point to circle
        double dx = x - x0;
        double dy = y - y0;
        double distToCenter = std::sqrt(dx*dx + dy*dy);
        double residual = std::abs(distToCenter - radius);
        
        // Add to chi2
        chi2 += (residual * residual) / (sigma * sigma);
        
        debug() << "  Point " << i << ": residual = " << residual 
                << " cm, σ = " << sigma << " cm, contribution = " 
                << (residual * residual) / (sigma * sigma) << endmsg;
    }
    
    // Degrees of freedom = number of points - number of parameters
    int ndf = 3; // 3 degree of freedom for 4 hits (4 hits * 2 D measurements - 5 track parameters)
    chi2 = chi2/ndf;
    
    if (chi2 / ndf > 100) { // Very poor fit
        debug() << "Poor circle fit: chi2/ndf = " << chi2 / ndf << endmsg;
        return false;
    }
    
    debug() << "Circle fit successful!" << endmsg;
    return true;
}

// Helper to get PDG code
int DisplacedTracking::getPDGCode() const {
    // Map particle type to PDG code
    if (m_particleType == "electron") return 11;
    else if (m_particleType == "positron") return -11;
    else if (m_particleType == "muon") return 13;
    else if (m_particleType == "antimuon") return -13;
    else if (m_particleType == "pion") return 211;
    else if (m_particleType == "kaon") return 321;
    else if (m_particleType == "proton") return 2212;
    else if (m_particleType == "antiproton") return -2212;
    else return 13; // Default to muon
}

// ============================================================================
// INNER SOLENOID PROPAGATION IMPLEMENTATIONS
// ============================================================================
bool DisplacedTracking::isInsideSolenoid(const Eigen::Vector3d& pos) const {
    double r = std::sqrt(pos.x() * pos.x() + pos.y() * pos.y());
    double z = std::abs(pos.z());
    return (r < SOLENOID_RADIUS_CM) && (z < SOLENOID_Z_HALF_CM);
}

void DisplacedTracking::rk4PropagationStep(
    Eigen::Vector3d& pos,
    Eigen::Vector3d& mom,
    double stepSize) const {

    double pMagInitial = mom.norm();  // Save initial momentum magnitude
    const double kappa_factor = 0.3 * 1000.0;  // Conversion factor for B*0.3 in cm^-1
    
    // Lambda for field at position
    auto getFieldAt = [this](const Eigen::Vector3d& r) -> double {
        dd4hep::Position ddPos(r.x(), r.y(), r.z());
        dd4hep::Direction B = m_field.magneticField(ddPos);
        return B.z() / dd4hep::tesla;
    };
    
    // Lambda for derivatives dr/ds, dp/ds
    auto derivatives = [&](const Eigen::Vector3d& r, const Eigen::Vector3d& p)
        -> std::pair<Eigen::Vector3d, Eigen::Vector3d> {
        
        double pMag = p.norm();
        if (pMag < 1e-12) {
            return {Eigen::Vector3d::Zero(), Eigen::Vector3d::Zero()};
        }
        
        Eigen::Vector3d pHat = p / pMag;
        double B = getFieldAt(r);
        Eigen::Vector3d B_vec(0, 0, B);
        
        // dr/ds = p̂
        Eigen::Vector3d drds = pHat;
        
        // dp/ds = q(p̂ × B) where q absorbed in field
        Eigen::Vector3d dpds = kappa_factor * pHat.cross(B_vec);
        
        return {drds, dpds};
    };

    // RK4 stages
    auto [k1r, k1p] = derivatives(pos, mom);
    auto [k2r, k2p] = derivatives(pos + 0.5 * stepSize * k1r, mom + 0.5 * stepSize * k1p);
    auto [k3r, k3p] = derivatives(pos + 0.5 * stepSize * k2r, mom + 0.5 * stepSize * k2p);
    auto [k4r, k4p] = derivatives(pos + stepSize * k3r, mom + stepSize * k3p);
    
    // Update state
    pos += (stepSize / 6.0) * (k1r + 2.0*k2r + 2.0*k3r + k4r);
    mom += (stepSize / 6.0) * (k1p + 2.0*k2p + 2.0*k3p + k4p);

    if (mom.norm() > 1e-12) {
        mom = mom.normalized() * pMagInitial;
    }
}

InnerTrajectory DisplacedTracking::propagateToInner(
    const Eigen::Vector3d& outerPos,
    const Eigen::Vector3d& outerMom,
    double targetRadius) const {
    
    InnerTrajectory result;
    result.success = false;
    result.numSteps = 0;
    result.arcLength = 0.0;
    
    // Starting position must be outside
    if (isInsideSolenoid(outerPos)) {
        result.message = "Starting position already inside solenoid";
        return result;
    }
    
    Eigen::Vector3d r = outerPos;
    Eigen::Vector3d p = outerMom;
    double arcLength = 0.0;
    
    const double RK4_STEP = -10.0;  // Negative = inward
    const int MAX_STEPS = 10000;
    
    info() << "\n=== Starting Inner Propagation ===" << endmsg;
    info() << "From R=" << std::sqrt(r.x()*r.x() + r.y()*r.y()) 
           << " cm, z=" << r.z() << " cm" << endmsg;
    info() << "Target radius: " << targetRadius << " cm" << endmsg;
    
    // Propagate
    for (int step = 0; step < MAX_STEPS; ++step) {
        double rxy = std::sqrt(r.x() * r.x() + r.y() * r.y());
        
        // Check if reached target
        if (isInsideSolenoid(r) && rxy < targetRadius) {
            result.success = true;
            result.message = "Reached target radius";
            result.finalPosition = r;
            result.finalMomentum = p;
            result.arcLength = arcLength;
            result.numSteps = step;
            break;
        }
        
        // RK4 step
        rk4PropagationStep(r, p, RK4_STEP);
        arcLength += std::abs(RK4_STEP);
        result.numSteps = step + 1;
        
        // Print every 1000 steps
        if (step % 1000 == 0 && msgLevel(MSG::DEBUG)) {
            double B = m_field.magneticField(dd4hep::Position(r.x(), r.y(), r.z())).z() / dd4hep::tesla;
            debug() << "Step " << step << ": R=" << rxy << " cm, B=" << B << " T" << endmsg;
        }
        
        // Safety
        if (step >= MAX_STEPS - 1) {
            result.message = "Max steps reached";
            result.finalPosition = r;
            result.finalMomentum = p;
            result.arcLength = arcLength;
            return result;
        }
    }
    
    return result;
}

const edm4hep::MCParticle* DisplacedTracking::findTruthParticle(
    const edm4hep::Track& recoTrack,
    const edm4hep::MCParticleCollection& mcParticles) const {
    
    if (recoTrack.trackerHits_size() == 0) {
        debug() << "Track has no hits - cannot find truth" << endmsg;
        return nullptr;
    }
    
    // Count which MCParticle appears most in the track hits
    std::map<const edm4hep::MCParticle*, int> hitCounts;
    
    for (int i = 0; i < recoTrack.trackerHits_size(); ++i) {
        const auto& hit = recoTrack.getTrackerHits(i);
        
        // For each MCParticle, check if it could have created this hit
        // (This is a simple matching - in reality, you'd use SimTrackerHit links)
        for (const auto& mcPart : mcParticles) {
            // Get truth vertex and momentum
            const auto& vertex = mcPart.getVertex();
            const auto& momentum = mcPart.getMomentum();
            
            // Check if this particle's trajectory passes near this hit
            // Simple check: does hit lie roughly on particle's path?
            const auto& hitPos = hit.getPosition();
            double dist = std::sqrt(
                (hitPos[0]/10.0 - vertex.x) * (hitPos[0]/10.0 - vertex.x) +
                (hitPos[1]/10.0 - vertex.y) * (hitPos[1]/10.0 - vertex.y) +
                (hitPos[2]/10.0 - vertex.z) * (hitPos[2]/10.0 - vertex.z)
            );
            
            // If hit is within 50cm of vertex, count it (rough association)
            if (dist < 50.0) {
                hitCounts[&mcPart]++;
            }
        }
    }
    
    // Find MCParticle with most hits
    const edm4hep::MCParticle* bestMatch = nullptr;
    int maxHits = 0;
    for (const auto& [particle, count] : hitCounts) {
        if (count > maxHits) {
            maxHits = count;
            bestMatch = particle;
        }
    }
    
    if (bestMatch && maxHits > 0) {
        info() << "Truth match: " << maxHits << " / " << recoTrack.trackerHits_size() 
               << " hits matched to MCParticle with momentum " 
               << std::sqrt(bestMatch->getMomentum().x * bestMatch->getMomentum().x +
                           bestMatch->getMomentum().y * bestMatch->getMomentum().y +
                           bestMatch->getMomentum().z * bestMatch->getMomentum().z)
               << " GeV/c" << endmsg;
    }
    
    return bestMatch;
}

void DisplacedTracking::validateTrackWithTruth(
    const edm4hep::Track& recoTrack,
    const edm4hep::MCParticleCollection& mcParticles,
    int trackNumber) const {
    
    // Find truth match
    const auto* truthParticle = findTruthParticle(recoTrack, mcParticles);
    if (!truthParticle) {
        debug() << "No truth match found for track " << trackNumber << endmsg;
        return;
    }
    
    // Get truth information
    const auto& truthVertex = truthParticle->getVertex();
    const auto& truthMom = truthParticle->getMomentum();
    
    double truthMomMag = std::sqrt(
        truthMom.x * truthMom.x +
        truthMom.y * truthMom.y +
        truthMom.z * truthMom.z
    );
    
    // Get reconstructed state (AtFirstHit for outer, AtOther for inner)
    edm4hep::TrackState recoState;
    bool foundOuter = false, foundInner = false;
    
    for (int i = 0; i < recoTrack.trackStates_size(); ++i) {
        auto st = recoTrack.getTrackStates(i);
        if (st.location == edm4hep::TrackState::AtFirstHit) {
            recoState = st;
            foundOuter = true;
            break;
        }
    }
    
    edm4hep::TrackState innerState;
    for (int i = 0; i < recoTrack.trackStates_size(); ++i) {
        auto st = recoTrack.getTrackStates(i);
        if (st.location == edm4hep::TrackState::AtOther) {
            innerState = st;
            foundInner = true;
            break;
        }
    }
    
    // Extract reconstructed position and momentum
    double d0 = recoState.D0 / 10.0;  // mm to cm
    double phi = recoState.phi;
    
    Eigen::Vector3d recoPos(
        d0 * std::cos(phi),
        d0 * std::sin(phi),
        recoState.Z0 / 10.0
    );
    
    // Calculate reconstructed momentum magnitude
    double omega = recoState.omega * 10.0;  // 1/mm to 1/cm
    double radius = 1.0 / std::abs(omega);
    dd4hep::Position fieldPos(recoPos.x(), recoPos.y(), recoPos.z());
    double bField = m_field.magneticField(fieldPos).z() / dd4hep::tesla;
    double lambda = std::atan(recoState.tanLambda);
    double recoMomMag = 0.3 * std::abs(bField) * radius / (100.0 * std::cos(lambda));
    
    // Calculate residuals
    double deltaX = recoPos.x() - truthVertex.x;
    double deltaY = recoPos.y() - truthVertex.y;
    double deltaZ = recoPos.z() - truthVertex.z;
    double deltaRxy = std::sqrt(deltaX*deltaX + deltaY*deltaY);
    
    double deltaMom = std::abs(recoMomMag - truthMomMag);
    double momResolution = (truthMomMag > 0) ? deltaMom / truthMomMag : 0;
    
    // Log comparison
    info() << "\n=== TRUTH COMPARISON (Track " << trackNumber << ") ===" << endmsg;
    info() << "Truth vertex: (" << truthVertex.x << ", " 
           << truthVertex.y << ", " << truthVertex.z << ") cm" << endmsg;
    info() << "Reco position (outer): (" << recoPos.x() << ", " 
           << recoPos.y() << ", " << recoPos.z() << ") cm" << endmsg;
    info() << "Position error: ΔR_xy=" << deltaRxy << " cm, ΔZ=" << deltaZ << " cm" << endmsg;
    
    info() << "Truth momentum: " << truthMomMag << " GeV/c" << endmsg;
    info() << "Reco momentum (outer): " << recoMomMag << " GeV/c" << endmsg;
    info() << "Momentum resolution: " << momResolution * 100.0 << "%" << endmsg;
    
    // Check curvature flip if inner exists
    if (foundInner) {
        double omegaInner = innerState.omega * 10.0;
        if (std::abs(omega) > 1e-10 && std::abs(omegaInner) > 1e-10) {
            if (omega * omegaInner < 0) {
                info() << "✓ Curvature properly flipped at boundary:" << endmsg;
                info() << "  ω_outer = " << omega << " 1/cm (sign=" << (omega > 0 ? "+" : "-") << ")" << endmsg;
                info() << "  ω_inner = " << omegaInner << " 1/cm (sign=" << (omegaInner > 0 ? "+" : "-") << ")" << endmsg;
            } else {
                warning() << "✗ Curvature DID NOT flip at boundary!" << endmsg;
            }
        }
    }
    
    info() << "\n";
}
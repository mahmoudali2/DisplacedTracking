#include <iostream>
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
            KeyValues{"InputHitCollection", {"TrackerHits"}},
            KeyValues{"HeaderCollectionName", {"EventHeader"}}
        },
        {
            KeyValues{"OutputTrackCollection", {"KalmanTracks"}}
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
    const std::vector<dd4hep::Position> testPoints = {
        dd4hep::Position(0, 0, 0),              
        dd4hep::Position(450, 0, 0),           // it expects dimensions in cm
        dd4hep::Position(0, 450, 0),           
        dd4hep::Position(0, 0, 450),           
        dd4hep::Position(500, 500, 0),
        dd4hep::Position(500, 600, 100),
        dd4hep::Position(600, 500, 200),
        dd4hep::Position(500, 500, 250),
        dd4hep::Position(100, 300, 350),
        dd4hep::Position(500, 500, 500)    
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

    // Convert to your vector of Surface pointers
    m_surfaces.clear();
    for (const auto& surf : surfList) {
        // Cast to Surface* if necessary
        auto surface = dynamic_cast<const dd4hep::rec::Surface*>(surf);
        if (surface) {
            m_surfaces.push_back(surface);
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
            //materialEffects = genfit::MaterialEffects::getInstance();
            //materialEffects->init(new genfit::TGeoMaterialInterface());
                    
            m_detector = m_geoSvc->getDetector();
            m_field = m_detector->field();
            m_genfitField=new GenfitField(m_field);
            m_geoMaterial=GenfitMaterialInterface::getInstance(m_detector);

            fieldManager = genfit::FieldManager::getInstance();
            fieldManager->init(m_genfitField); // kGauss
            /*
            // Create and register the field adapter
            m_genFitField = std::unique_ptr<DD4hepFieldAdapter>(new DD4hepFieldAdapter(m_field));
            genfit::FieldManager::getInstance()->init(m_genFitField.get());
            
            // Create and register the material adapter
            m_genFitMaterial = std::unique_ptr<DD4hepMaterialAdapter>(new DD4hepMaterialAdapter(m_materialManager));
            genfit::MaterialEffects::getInstance()->init(m_genFitMaterial.get());
 */           
            // Optionally set the multiple scattering model
            genfit::MaterialEffects::getInstance()->setMscModel("GEANE");
            
            info() << "GenFit initialized successfully" << endmsg;
        } catch (const genfit::Exception& e) {
            error() << "Error initializing GenFit: " << e.what() << endmsg;
            return StatusCode::FAILURE;
        }
    }       

    return StatusCode::SUCCESS;
}

// Operator method
std::tuple<edm4hep::TrackCollection> DisplacedTracking::operator()(
    const edm4hep::TrackerHitPlaneCollection& hits,
    const edm4hep::EventHeaderCollection& headers) const {
    
    // Find tracks using the direct EDM4hep approach
    edm4hep::TrackCollection trackCollection;

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
        return std::make_tuple(std::move(trackCollection)); // Return empty collection
    }
    
    try {
        trackCollection = findTracks(&hits);
    } catch (const std::exception& ex) {
        error() << "Exception during track finding: " << ex.what() << endmsg;
        info() << std::string(80, '=') << "\n" << endmsg; // Bottom separator
       return std::make_tuple(std::move(trackCollection)); // Return empty collection
    }
    
    info() << "Found " << trackCollection.size() << " tracks" << endmsg;
    
    // Print track details
    for (size_t i = 0; i < trackCollection.size(); ++i) {
        const auto& track = trackCollection[i];
        
        // Get track state at IP if available
        edm4hep::TrackState state;
        bool foundState = false;
        
        for (int j = 0; j < track.trackStates_size(); ++j) {
            if (track.getTrackStates(j).location == edm4hep::TrackState::AtOther) {
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
            
            // Use default if field is too small
            if (std::abs(bField) < 0.1) {
                debug() << "Magnetic field too small (" << bField 
                        << " T), using default value -1.7 T" << endmsg;
                bField = -1.7;  // Default field value
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
    
    return std::make_tuple(std::move(trackCollection));
}

// Finalize method
StatusCode DisplacedTracking::finalize() {
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

    // Use default field if actual field is too small -- temporary workaround
    if (std::abs(bField) < 0.1) {
        bField = -1.7; // Tesla
    }

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
    
    // Create covariance matrix with 21 elements (6x6 symmetric)
    std::array<float, 21> covValues = {0};  // Initialize all to zero
    
    // Set diagonal elements for our 5 track parameters
    // The indices would be different in a 6x6 matrix, so we need to map them correctly
    covValues[0] = 1.0;   // d0 variance
    covValues[2] = 0.01;  // phi variance
    covValues[5] = 1e-6;  // omega variance 
    covValues[9] = 1.0;   // z0 variance
    covValues[14] = 0.01; // tanLambda variance
    
    // Create covariance matrix with individual elements
    // Note: EDM4hep expects us to set each element individually
    for (int i = 0; i < 6; ++i) {
        for (int j = 0; j <= i; ++j) {
            int index = i*(i+1)/2 + j;
            if (index < 21) {
                state.setCovMatrix(covValues[index], static_cast<edm4hep::TrackParams>(i), 
                                                    static_cast<edm4hep::TrackParams>(j));
            }
        }
    }
    
    return state;
}

//------------------------------------------------------------------------------
// Core tracking methods
//------------------------------------------------------------------------------

edm4hep::TrackCollection DisplacedTracking::findTracks(
    const edm4hep::TrackerHitPlaneCollection* hits) const {
    
    // Create output collection for track candidates
    edm4hep::TrackCollection candidateTracks;
    
    // Check if we have enough hits for tracking
    if (hits->size() < 3) {
        info() << "Not enough hits for tracking. Need at least 3, found " << hits->size() << endmsg;
        return std::move(candidateTracks);
    }
    
    // Group hits by layer
    std::map<int, std::vector<std::pair<size_t, edm4hep::TrackerHitPlane>>> hitsByLayer;
    
    // Find surfaces and organize hits by layer
    for (size_t i = 0; i < hits->size(); ++i) {
        const auto& hit = (*hits)[i];
        const dd4hep::rec::Surface* surface = findSurface(hit);
        
        if (surface) {
            int layerID = getLayerID(hit.getCellID());
            int typeID = getTypeID(hit.getCellID());
            hitsByLayer[layerID].push_back(std::make_pair(i, hit));
            debug() << "Found hit at CellId " << hit.getCellID() << " ,in Type: " << typeID << " , and layer: " << layerID << endmsg;
        } else {
            debug() << "Could not find surface for hit with cellID " << hit.getCellID() << endmsg;
        }
    }
    
    info() << "Found hits in " << hitsByLayer.size() << " layers" << endmsg;
    
    // Display hit distribution by layer
    if (msgLevel(MSG::DEBUG)) {
        debug() << "Hit distribution by layer:" << endmsg;
        for (const auto& [layer, layerHits] : hitsByLayer) {
            debug() << "  Layer " << layer << ": " << layerHits.size() << " hits" << endmsg;
        }
    }
    
    // Get sorted layer IDs
    std::vector<int> layerIDs;
    for (const auto& layer : hitsByLayer) {
        layerIDs.push_back(layer.first);
    }
    std::sort(layerIDs.begin(), layerIDs.end());
    
    // Report seeding configuration
    if (layerIDs.size() >= 3) {
        info() << "Seeding with " << std::min(layerIDs.size(), size_t(3)) << " layers" << endmsg;
        
        if (layerIDs.size() > 3) {
            info() << "Additional layers available for extension" << endmsg;
        }
    } else {
        info() << "Not enough layers with hits for triplet seeding (need 3, found " 
               << layerIDs.size() << ")" << endmsg;
        return std::move(candidateTracks);
    }
    
    // Keep track of which hits are used
    std::vector<bool> usedHits(hits->size(), false);
    
    // Structure to hold track candidates
    struct TrackCandidate {
        edm4hep::Track track;
        double pT;
        std::vector<size_t> hitIndices;
        int highestLayerID;
        
        TrackCandidate(edm4hep::Track t, double pt, const std::vector<size_t>& indices, int layerID)
            : track(t), pT(pt), hitIndices(indices), highestLayerID(layerID) {}
    };
    
    std::vector<TrackCandidate> trackCandidates;
    
    // Try consecutive layer triplets first
    for (size_t i = 0; i <= layerIDs.size() - 3; ++i) {
        int layer1 = layerIDs[i];
        int layer2 = layerIDs[i+1];
        int layer3 = layerIDs[i+2];

        int tripletCandidates = 0;
        int validTriplets = 0;

        if (layer3 == 3 && validTriplets >= 1){
            debug() << "Skipping due to found valid triplets in first 3 layers " << endmsg;
            break;
        } else {
            debug() << "Keep searching for valid triplets in +3 layers " << endmsg;
        }

        debug() << "Testing hit triplets from layers " << layer1 << ", " 
                << layer2 << ", " << layer3 << endmsg;
        
        // Try all hit combinations from these three layers
        for (const auto& [idx1, hit1] : hitsByLayer[layer1]) {
            for (const auto& [idx2, hit2] : hitsByLayer[layer2]) {
                for (const auto& [idx3, hit3] : hitsByLayer[layer3]) {
                    tripletCandidates++;
                    
                    // Create a temporary vector for tracking used hits
                    std::vector<bool> tempUsedHits = usedHits;
                    
                    // Try to create a triplet seed
                    size_t prevSize = candidateTracks.size();
                    bool seedValid = createTripletSeed(
                        hit1, hit2, hit3, &candidateTracks, tempUsedHits, idx1, idx2, idx3);
                    
                    if (seedValid && candidateTracks.size() > prevSize) {
                        validTriplets++;
                        
                        // Get the new track
                        auto newTrack = candidateTracks[candidateTracks.size() - 1];
                        
                        // Calculate pT
                        double pT = 0.;
                        bool foundPt = false;
                        for (int j = 0; j < newTrack.trackStates_size(); ++j) {
                            auto state = newTrack.getTrackStates(j);
                            if (state.location == edm4hep::TrackState::AtOther) {
                                pT = getPT(state);
                                foundPt = true;
                                break;
                            }
                        }

                        if (!foundPt) {
                            std::cerr << "Warning: No TrackState::AtOther found, pT = 0\n";
                        }
                        
                        // Store hit indices
                        std::vector<size_t> hitIndices = {idx1, idx2, idx3};
                        
                        // Add to track candidates
                        trackCandidates.emplace_back(newTrack, pT, hitIndices, layer3);
                    }
                }
            }
        }
        
        debug() << "Tested " << tripletCandidates << " triplet candidates, found " 
                << validTriplets << " valid triplets" << endmsg;
    }
    
    info() << "Found " << trackCandidates.size() << " triplet track candidates" << endmsg;
    
    // Create the final collection for high-quality tracks
    edm4hep::TrackCollection finalTracks;
    
    // Sort candidates by decreasing pT
    std::sort(trackCandidates.begin(), trackCandidates.end(),
              [](const TrackCandidate& a, const TrackCandidate& b) {
                  return a.pT > b.pT;
              });
    
    // Process candidates in order of pT, avoiding conflicts
    for (const auto& candidate : trackCandidates) {
        // Skip if any hits are already used
        bool hasConflict = false;
        for (size_t idx : candidate.hitIndices) {
            if (usedHits[idx]) {
                hasConflict = true;
                break;
            }
        }
        
        if (hasConflict) {
            continue;
        }
        
        // Get the seed track and its state
        const auto& seedTrack = candidate.track;
        
        // Find a suitable seed state
        edm4hep::TrackState seedState;
        bool foundState = false;
        
        for (int j = 0; j < seedTrack.trackStates_size(); ++j) {
            if (seedTrack.getTrackStates(j).location == edm4hep::TrackState::AtOther) {
                seedState = seedTrack.getTrackStates(j);
                foundState = true;
                break;
            }
        }
        
        if (!foundState && seedTrack.trackStates_size() > 0) {
            seedState = seedTrack.getTrackStates(0);
            foundState = true;
        }
        
        if (!foundState) {
            warning() << "No seed state found for track, skipping" << endmsg;
            continue;
        }
        
        // Get the original hits
        std::vector<edm4hep::TrackerHitPlane> trackHits;
        for (int j = 0; j < seedTrack.trackerHits_size(); ++j) {
            auto trackHit = seedTrack.getTrackerHits(j);
            
            // Find the matching hit in the collection
            for (size_t i = 0; i < hits->size(); ++i) {
                const auto& hit = (*hits)[i];
                if (hit.getCellID() == trackHit.getCellID()) {
                    trackHits.push_back(hit);
                    break;
                }
            }
        }
        
        // Create a new track
        auto finalTrack = finalTracks.create();
        
        // Fit the track with GenFit, passing the full hit collection
        bool fitSuccess = fitTrackWithGenFit(
            trackHits, seedState, finalTrack, hits);
        
        if (fitSuccess) {
            info() << "Successfully fitted track with " << finalTrack.trackerHits_size() 
                   << " hits, chi2/ndf = " << finalTrack.getChi2() / finalTrack.getNdf() << endmsg;
                   
            // Mark hits as used
            for (size_t idx : candidate.hitIndices) {
                usedHits[idx] = true;
            }
            
            // If the track was extended with an additional hit, mark that hit as used too
            if (finalTrack.trackerHits_size() > trackHits.size()) {
                info() << "Track was extended with " 
                       << (finalTrack.trackerHits_size() - trackHits.size()) 
                       << " additional hit(s)" << endmsg;
                       
                // We need to find which hit was added and mark it as used
                for (int j = 0; j < finalTrack.trackerHits_size(); ++j) {
                    auto finalHit = finalTrack.getTrackerHits(j);
                    
                    // Check if this hit is in the original trackHits
                    bool isOriginalHit = false;
                    for (const auto& origHit : trackHits) {
                        if (finalHit.getCellID() == origHit.getCellID()) {
                            isOriginalHit = true;
                            break;
                        }
                    }
                    
                    // If not in original hits, find and mark as used
                    if (!isOriginalHit) {
                        // Find the hit index in the original collection
                        for (size_t i = 0; i < hits->size(); ++i) {
                            if ((*hits)[i].getCellID() == finalHit.getCellID() && !usedHits[i]) {
                                usedHits[i] = true;
                                debug() << "Marked added hit at index " << i << " as used" << endmsg;
                                break;
                            }
                        }
                    }
                }
            }
        } else {
            warning() << "Track fitting failed, using seed parameters" << endmsg;
            
            // Fall back to seed parameters
            finalTrack.setChi2(seedTrack.getChi2());
            finalTrack.setNdf(seedTrack.getNdf());
            
            // Copy track states
            for (int j = 0; j < seedTrack.trackStates_size(); ++j) {
                finalTrack.addToTrackStates(seedTrack.getTrackStates(j));
            }
            
            // Add hits to track
            for (const auto& hit : trackHits) {
                finalTrack.addToTrackerHits(hit);
            }
            
            // Mark hits as used
            for (size_t idx : candidate.hitIndices) {
                usedHits[idx] = true;
            }
        }
    }
    
    info() << "Found " << finalTracks.size() << " tracks" << endmsg;
    
    // Print track details
    for (size_t i = 0; i < finalTracks.size(); ++i) {
        const auto& track = finalTracks[i];
        
        // Get track state at IP if available
        edm4hep::TrackState state;
        bool foundState = false;
        
        for (int j = 0; j < track.trackStates_size(); ++j) {
            if (track.getTrackStates(j).location == edm4hep::TrackState::AtOther) {
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
            // Calculate parameters for display
            double omega = state.omega;
            double bField = -1.7; // Tesla
            double pT = 0.3 * std::abs(bField) / std::abs(omega) * 0.001; // GeV/c
            
            double tanLambda = state.tanLambda;
            double theta = std::atan2(1.0, tanLambda);
            double eta = -std::log(std::tan(theta/2.0));
            
            double d0 = state.D0 / 10.0;  // mm to cm
            double z0 = state.Z0 / 10.0;  // mm to cm
        }
    }
    
    return std::move(finalTracks);
}

// New function to calculate circle center and radius using the Direct Formula Method
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
    const double estimatedBz = -1.7; // Tesla
    
    debug() << "Magnetic field: actual=" << actualBz << " Tesla, using estimated=" 
            << estimatedBz << " Tesla for calculation as a temporary solution, due to a problem in k4geo in retrieving the right value of Bz" << endmsg;
    
    // Calculate pT using all methods for comparison
    double pT_direct = 0.3 * std::abs(estimatedBz) * radius_direct / 100.0;
    double pT_sagitta = 0.3 * std::abs(estimatedBz) * sagittaRadius / 100.0;
    double pT_sagitta_full = 0.3 * std::abs(estimatedBz) * radius_sagitta / 100.0;
    
    debug() << "Comparison of pT estimates:" << endmsg;
    debug() << "  Direct formula pT: " << pT_direct << " GeV/c" << endmsg;
    debug() << "  Sagitta simple pT: " << pT_sagitta << " GeV/c" << endmsg;
    debug() << "  Sagitta full pT: " << pT_sagitta_full << " GeV/c" << endmsg;
    
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
    double charge = clockwise ? 1.0 : -1.0;
    
    debug() << "Track direction: " << (clockwise ? "clockwise" : "counter-clockwise") 
            << ", charge: " << charge << endmsg;
    
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
    // Determine inner and outer fields
    double innerFieldStrength = 2.0;   // 2T inside solenoid
    double outerFieldStrength = -1.7;  // -1.7T outside solenoid

    double d0 = std::sqrt(std::pow(x0, 2) + std::pow(y0, 2)) - radius;
    /*
    // Use the two-segment model for d0 calculation
    double d0 = calculateImpactParameter(x0, y0, radius, clockwise, 
                                   innerFieldStrength, outerFieldStrength,
                                   p1, p2, p3);
    */
    double z0 = b; // z0 calculation remains the same
    
    debug() << "Impact parameters: d0=" << d0 << " cm, z0=" << z0 << " cm" << endmsg;
    
    std::cout << "Debug d0: centerToOrigin=" << std::sqrt(std::pow(x0, 2) + std::pow(y0, 2)) 
          << "cm, radius=" << radius << "cm, raw d0=" 
          << (std::sqrt(std::pow(x0, 2) + std::pow(y0, 2)) - radius) << "cm" << std::endl;

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
        d0_mm, phi, omega, z0_mm, tanLambda, edm4hep::TrackState::AtOther);

    // Set reference point to the first hit position
    const auto& firstHitPos = hit1.getPosition();
    state.referencePoint = edm4hep::Vector3f(firstHitPos[0], firstHitPos[1], firstHitPos[2]);

    // Add track state to track
    edm_track.addToTrackStates(state);
    
    // Add hits to track
    edm_track.addToTrackerHits(hit1);
    edm_track.addToTrackerHits(hit2);
    edm_track.addToTrackerHits(hit3);
    
    // Set chi2 and ndf
    edm_track.setChi2(-1.0);  // Initial seed has no chi2 yet
    //edm_track.setNdf(2 * 3 - 5);  // 2 DOF per hit, 5 parameters
    
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
    const double solenoidRadius = 200.0; // cm - radius of the solenoid
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
    
    // Extract EDM4hep track parameters
    double d0 = state.D0;          // Impact parameter in mm
    double phi = state.phi;        // Azimuthal angle at PCA
    double omega = state.omega;    // Signed curvature (1/R) in 1/mm
    double z0 = state.Z0;          // Z position at PCA
    double tanLambda = state.tanLambda; // Tangent of dip angle
    
    // Use reference point (first hit position)
    double x = state.referencePoint[0];
    double y = state.referencePoint[1];
    double z = state.referencePoint[2];
    
    // Get field at this position to calculate momentum
    dd4hep::Position fieldPos(x/10.0, y/10.0, z/10.0); // Convert mm to cm for DD4hep
    double bField = m_field.magneticField(fieldPos).z() / dd4hep::tesla;
    
    // Use default field if actual field is too small
    if (std::abs(bField) < 0.1) {
        bField = -1.7; // Default field value in Tesla
    }
    
    // Calculate pT from curvature
    double pT = 0.3 * std::abs(bField) / (std::abs(omega) * 1000.0); // GeV/c
    
    // Calculate center of helix
    double radius = 1.0 / std::abs(omega); // mm
    double centerX, centerY;
    
    if (omega < 0) { // Positive charge
        centerX = -d0 * sin(phi) - radius * cos(phi);
        centerY = d0 * cos(phi) - radius * sin(phi);
    } else { // Negative charge
        centerX = -d0 * sin(phi) + radius * cos(phi);
        centerY = d0 * cos(phi) + radius * sin(phi);
    }
    debug() << "Circle Radius: " << radius;
    debug() << "Circle center, x: " << centerX;
    debug() << "Circle center, y: " << centerY;

    // Calculate momentum direction at the first hit
    double hitX = x/10.0; // convert mm to cm
    double hitY = y/10.0;
    
    // Vector from center to hit position
    double vecX = hitX - centerX/10.0; // convert center to cm
    double vecY = hitY - centerY/10.0;
    double vecMag = sqrt(vecX*vecX + vecY*vecY);
    
    // Normalize vector
    vecX /= vecMag;
    vecY /= vecMag;
    
    // Rotate 90 degrees to get tangent direction (momentum)
    double momDirX, momDirY;
    if (omega < 0) { // Positive charge -> clockwise
        momDirX = -vecY;
        momDirY = vecX;
    } else { // Negative charge -> counter-clockwise
        momDirX = vecY;
        momDirY = -vecX;
    }
    
    // Create momentum vector at the hit
    double px = pT * momDirX;
    double py = pT * momDirY;
    double pz = pT * tanLambda;
    double p = sqrt(px*px + py*py + pz*pz);
    
    // Charge from omega sign
    double charge = (omega < 0) ? 1.0 : -1.0;
    
    // Create GenFit state vectors
    TVector3 posVec(x/10.0, y/10.0, z/10.0); // Convert mm to cm
    TVector3 momVec(px, py, pz);
    
    // Create a new GenFit state on the reference plane
    genfit::MeasuredStateOnPlane state_gf(rep);
    
    // Set state parameters
    state_gf.setPosMom(posVec, momVec);
    state_gf.setQop(charge / p); // q/p = charge / momentum magnitude
    
    // Convert covariance matrix from EDM4hep to GenFit format
    // This is a simplified conversion - a full conversion would need careful parameter mapping
    TMatrixDSym covMat(6); // 6x6 symmetric matrix

    // This is a simplified mapping - a proper mapping would require coordinate transformation
    covMat(0, 0) = 0.004; //xx
    covMat(1, 1) = 0.004; //yy
    covMat(2, 2) = 0.004; //zz
    covMat(3, 3) = 0.01; //pxpx
    covMat(4, 4) = 0.01; //pypy
    covMat(5, 5) = 0.01; //pzpz
    
    state_gf.setCov(covMat);
    
    return state_gf;
}


edm4hep::TrackState DisplacedTracking::convertToEDM4hepState(
    const genfit::MeasuredStateOnPlane& state,
    int location) const {
    
    // Get position and momentum from GenFit state
    TVector3 pos = state.getPos(); // cm
    TVector3 mom = state.getMom(); // GeV/c
    double charge = (state.getQop() > 0) ? 1.0 : -1.0;
    
    // Calculate derived quantities
    double px = mom.X();
    double py = mom.Y();
    double pz = mom.Z();
    double pt = sqrt(px*px + py*py);
    double p = mom.Mag();
    
    // Calculate EDM4hep track parameters
    double phi = atan2(py, px);
    double tanLambda = pz / pt;
    
    // Get field at this position to calculate curvature
    dd4hep::Position fieldPos(pos.X(), pos.Y(), pos.Z());
    double bField = m_field.magneticField(fieldPos).z() / dd4hep::tesla;
    
    // Use default field if actual field is too small
    if (std::abs(bField) < 0.1) {
        bField = -1.7; // Default field value in Tesla
    }
    
    // Calculate curvature (1/R) from pT and B
    // pT [GeV/c] = 0.3 * |B| [T] * R [m]
    // omega [1/mm] = 1/R [1/m] * 0.001
    double omega = charge * 0.3 * std::abs(bField) / (pt * 1000.0); // 1/mm
    
    // Calculate impact parameter d0
    // For a helix, d0 is the distance from (0,0) to the center of the helix, minus the radius
    // In this simplified calculation, we project the position onto the x-y plane to get d0
    double d0 = -pos.Y() * cos(phi) + pos.X() * sin(phi);
    d0 *= 10.0; // Convert from cm to mm
    
    // Calculate z position at PCA
    double z0 = pos.Z() * 10.0; // Convert from cm to mm
    
    // Create EDM4hep track state with appropriate parameters
    edm4hep::TrackState outState;
    outState.D0 = d0;
    outState.phi = phi;
    outState.omega = omega;
    outState.Z0 = z0;
    outState.tanLambda = tanLambda;
    outState.location = location;
    
    // Get covariance matrix from GenFit state
    const TMatrixDSym& covMat = state.getCov();
    
    // Convert covariance matrix to EDM4hep format
    // Only access up to 5x5 dimensions to avoid out-of-bounds errors
    for (int i = 0; i < 5; ++i) {
        for (int j = 0; j <= i; ++j) {
            // Apply appropriate scaling based on parameter type
            double value = covMat(i, j);
            
            // Scale position elements from cm² to mm²
            if (i >= 3 && j >= 3) {
                value *= 100.0;  // cm² to mm²
            }
            // Scale mixed position-momentum elements
            else if (i >= 3 || j >= 3) {
                value *= 10.0;   // cm to mm for mixed terms
            }
            
            // Set covariance matrix element
            outState.setCovMatrix(value, 
                               static_cast<edm4hep::TrackParams>(i), 
                               static_cast<edm4hep::TrackParams>(j));
        }
    }
    
    debug() << "Converted track state: d0=" << d0 << "mm, phi=" << phi 
           << ", omega=" << omega << "1/mm, z0=" << z0 
           << "mm, tanLambda=" << tanLambda << endmsg;
    
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
    if (sigma_u <= 0) sigma_u = 0.1;  // 100 microns default
    if (sigma_v <= 0) sigma_v = 0.1;  // 100 microns default
    
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
        
        // Convert the seed state to GenFit format
        genfit::MeasuredStateOnPlane seedGFState = convertToGenFitState(seedState, rep);

        // Extract position and momentum from the MeasuredStateOnPlane
        TVector3 pos = seedGFState.getPos();
        TVector3 mom = seedGFState.getMom();

        debug() << "Pos initialization:" << endmsg;
        pos.Print();
        debug() << "Mom initialization:" << endmsg;
        mom.Print();
        
        //------------------------------------------------------------------
        // Create a GenFit track with the seed hits (no extrapolation)
        //------------------------------------------------------------------
        
        debug() << "Building final track with " << hits.size() << " seed hits (no extrapolation)" << endmsg;
        
        // Create a new GenFit track
        genfit::Track finalGFTrack(rep, pos, mom);
        
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
        // Fit the track with GenFit
        //------------------------------------------------------------------
        
        // Create and configure Kalman fitter
        genfit::KalmanFitterRefTrack fitter;
        fitter.setMaxIterations(m_maxFitIterations);
        fitter.setMinIterations(2);
        fitter.setDebugLvl(1); // Minimal output
        fitter.setBlowUpMaxVal(1000.0);
        
        // Perform the fit
        debug() << "Starting GenFit Kalman fitting with " << finalGFTrack.getNumPoints() << " hits" << endmsg;
        fitter.processTrack(&finalGFTrack);
        
        // Check fit quality
        if (!finalGFTrack.getFitStatus()->isFitted()) {
            warning() << "GenFit track fitting failed" << endmsg;
            return false;
        }
        
        // Get fit quality metrics
        double chi2 = finalGFTrack.getFitStatus()->getChi2();
        int ndf = finalGFTrack.getFitStatus()->getNdf();
        
        debug() << "Track fitted successfully: chi2=" << chi2 
                << ", ndf=" << ndf << ", chi2/ndf=" << (chi2/(ndf*1.0)) << endmsg;
        
        // Update the EDM4hep track with fit results
        finalTrack.setChi2(chi2);
        finalTrack.setNdf(ndf);

        // Add hits to the EDM4hep track
        for (const auto& hit : hits) {
            finalTrack.addToTrackerHits(hit);
        }
        // Add the original seed state to the track
        edm4hep::TrackState seedStateCopy = seedState;
        seedStateCopy.location = edm4hep::TrackState::AtOther;  // Mark it as "other" location
        finalTrack.addToTrackStates(seedStateCopy);

        // Get the track states at key positions and add them to EDM4hep track
        try {
            // State at first hit
            if (finalGFTrack.getNumPoints() > 0) {
                genfit::MeasuredStateOnPlane stateFirst = 
                    finalGFTrack.getFittedState(0);
                edm4hep::TrackState firstState = convertToEDM4hepState(stateFirst, 
                                                                      edm4hep::TrackState::AtFirstHit);
                finalTrack.addToTrackStates(firstState);
            }
/*           
            // State at last hit
            if (finalGFTrack.getNumPoints() > 1) {
                genfit::MeasuredStateOnPlane stateLast = 
                    finalGFTrack.getFittedState(finalGFTrack.getNumPoints() - 1);
                edm4hep::TrackState lastState = convertToEDM4hepState(stateLast, 
                                                                     edm4hep::TrackState::AtLastHit);
                finalTrack.addToTrackStates(lastState);
            }
 **/            
            // Add state at IP or other reference point
            //edm4hep::TrackState refState = convertToEDM4hepState(seedGFState, 
            //                                                    edm4hep::TrackState::AtOther);
            //finalTrack.addToTrackStates(refState);
            
        } catch (genfit::Exception& e) {
            warning() << "Error extracting track states: " << e.what() << endmsg;
            // Continue anyway - we'll use what we have
        }
        
        debug() << "Successfully created " << finalTrack.trackStates_size() 
                << " track states from GenFit fit" << endmsg;
        
        return true;
        
    } catch (genfit::Exception& e) {
        error() << "GenFit exception during track fitting: " << e.what() << endmsg;
        return false;
    }
}

// Extrapolate the track to +3 layers
genfit::MeasuredStateOnPlane DisplacedTracking::extrapolateToSurface(
    const genfit::MeasuredStateOnPlane& state,
    genfit::AbsTrackRep* rep,
    const dd4hep::rec::Surface* targetSurface) const {
    
    // Extract surface properties
    dd4hep::rec::Vector3D normal = targetSurface->normal();
    dd4hep::rec::Vector3D origin = targetSurface->origin();
    
    // Convert to ROOT's TVector3
    TVector3 planeOrigin(origin.x(), origin.y(), origin.z());  // Already in cm
    TVector3 planeNormal(normal.x(), normal.y(), normal.z());
    
    // Create destination plane
    genfit::SharedPlanePtr destPlane(new genfit::DetPlane(planeOrigin, planeNormal));
    
    // Clone state to avoid modifying the original
    genfit::MeasuredStateOnPlane newState(state);
    
    try {
        // Extrapolate to the destination plane
        rep->extrapolateToPlane(newState, destPlane);
        return newState;
    }
    catch(genfit::Exception& e) {
        warning() << "Extrapolation failed: " << e.what() << endmsg;
        throw; // Re-throw to allow caller to handle the error
    }
}

std::vector<std::pair<size_t, edm4hep::TrackerHitPlane>> DisplacedTracking::findCompatibleHitsInLayer(
    const genfit::MeasuredStateOnPlane& state,
    genfit::AbsTrackRep* rep,
    int layerID,
    const edm4hep::TrackerHitPlaneCollection& allHits,
    double maxChi2) const {
    
    std::vector<std::pair<size_t, edm4hep::TrackerHitPlane>> compatibleHits;
    
    // Get all surfaces in the target layer
    auto layerIter = m_surfacesByLayer.find(layerID);
    if (layerIter == m_surfacesByLayer.end()) {
        return compatibleHits; // No surfaces in this layer
    }
    
    const auto& layerSurfaces = layerIter->second;
    
    // For each hit in the collection
    for (size_t i = 0; i < allHits.size(); ++i) {
        const auto& hit = allHits[i];
        
        // Check if hit is in the target layer
        int hitLayerID = getLayerID(hit.getCellID());
        if (hitLayerID != layerID) continue;
        
        // Find surface for this hit
        const dd4hep::rec::Surface* surface = findSurface(hit);
        if (!surface) continue;
        
        try {
            // Extrapolate state to this surface
            genfit::MeasuredStateOnPlane extrapolatedState = 
                extrapolateToSurface(state, rep, surface);
            
            // Get hit position and error
            const auto& pos = hit.getPosition();
            dd4hep::rec::Vector3D hitGlobal(pos[0], pos[1], pos[2]);
            dd4hep::rec::Vector2D hitLocal = surface->globalToLocal(dd4hep::mm * hitGlobal);
            
            // Get extrapolated position on surface
            TVector3 extrapPos = extrapolatedState.getPos();
            dd4hep::rec::Vector3D extrapGlobal(extrapPos.X(), extrapPos.Y(), extrapPos.Z());
            dd4hep::rec::Vector2D extrapLocal = surface->globalToLocal(extrapGlobal);
            
            // Calculate chi2
            double du = hit.getDu() > 0 ? hit.getDu() : 0.1; // Default if not set
            double dv = hit.getDv() > 0 ? hit.getDv() : 0.1; // Default if not set
            
            // Convert from mm to cm for comparison
            du /= 10.0;
            dv /= 10.0;
            
            double residualU = (hitLocal[0] - extrapLocal[0]);
            double residualV = (hitLocal[1] - extrapLocal[1]);
            
            double chi2 = (residualU*residualU)/(du*du) + (residualV*residualV)/(dv*dv);
            
            if (chi2 < maxChi2) {
                compatibleHits.push_back(std::make_pair(i, hit));
            }
        }
        catch (genfit::Exception&) {
            // Extrapolation failed for this surface, continue to next hit
            continue;
        }
    }
    
    return compatibleHits;
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
// Helper to get seed state
edm4hep::TrackState DisplacedTracking::getSeedState(const edm4hep::Track& track) const {
    // Find a state location AtOther (our seed)
    for (int i = 0; i < track.trackStates_size(); ++i) {
        if (track.getTrackStates(i).location == edm4hep::TrackState::AtOther) {
            return track.getTrackStates(i);
        }
    }
    
    // Fall back to first state
    if (track.trackStates_size() > 0) {
        return track.getTrackStates(0);
    }
    
    // Create empty state if none found
    edm4hep::TrackState emptyState;
    return emptyState;
}

// Helper to get original hits
std::vector<edm4hep::TrackerHitPlane> DisplacedTracking::getOriginalHits(
    const edm4hep::Track& track,
    const edm4hep::TrackerHitPlaneCollection& allHits) const {
    
    std::vector<edm4hep::TrackerHitPlane> hits;
    
    for (int j = 0; j < track.trackerHits_size(); ++j) {
        auto trackHit = track.getTrackerHits(j);
        
        // Find matching hit in collection by ID
        for (size_t i = 0; i < allHits.size(); ++i) {
            if (allHits[i].getCellID() == trackHit.getCellID()) {
                hits.push_back(allHits[i]);
                break;
            }
        }
    }
    
    return hits;
}

// Helper to get highest layer
int DisplacedTracking::getHighestLayerID(const std::vector<edm4hep::TrackerHitPlane>& hits) const {
    int highest = -1;
    
    for (const auto& hit : hits) {
        int layer = getLayerID(hit.getCellID());
        highest = std::max(highest, layer);
    }
    return highest;
}

// Helper to copy parameters
void DisplacedTracking::copyTrackParameters(
    const edm4hep::Track& source, 
    edm4hep::MutableTrack& target) const {
    
    target.setChi2(source.getChi2());
    target.setNdf(source.getNdf());
    
    // Copy track states
    for (int i = 0; i < source.trackStates_size(); ++i) {
        target.addToTrackStates(source.getTrackStates(i));
    }
}
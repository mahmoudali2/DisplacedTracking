#ifndef KALMAN_TRACKING_H
#define KALMAN_TRACKING_H

#include <vector>
#include <map>
#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>

#include "GaudiKernel/Algorithm.h"
#include "Gaudi/Property.h"
#include "GaudiKernel/ServiceHandle.h"
#include "GAUDI_VERSION.h"

#include <Eigen/Dense>
#include "DD4hep/Detector.h"
#include "DD4hep/Fields.h"
#include "DD4hep/DD4hepUnits.h"
#include "DDRec/Surface.h"
#include "DDRec/SurfaceManager.h"
#include "DDRec/MaterialManager.h"
#include "DDSegmentation/BitFieldCoder.h"

#include "k4FWCore/DataHandle.h"
#include "k4FWCore/Transformer.h"
//#include "k4FWCore/MultiTransformer.h"
#include "k4FWCore/PodioDataSvc.h"
#include "k4Interface/IUniqueIDGenSvc.h"

// EDM4hep includes
#include "edm4hep/TrackerHitPlaneCollection.h"
#include "edm4hep/TrackCollection.h"
#include "edm4hep/EventHeaderCollection.h"
#include "edm4hep/TrackerHitSimTrackerHitLinkCollection.h"

// Interface includes
#include "k4Interface/IGeoSvc.h"

// Forward declarations
namespace dd4hep {
    class OverlayedField;
    namespace rec {
        class Surface;
        class MaterialManager;
    }
}

// Define the 5-parameter vector and matrix types
namespace Eigen {
    typedef Matrix<double, 5, 1> Vector5d;
    typedef Matrix<double, 5, 5> Matrix5d;
}

/**
 * Particle properties struct
 */
struct ParticleProperties {
    double mass;          // Mass in GeV
    double charge;        // Charge in e
    std::string name;     // Particle name

    // Add a default constructor
    ParticleProperties() : mass(0), charge(0), name("unknown") {}

    ParticleProperties(double m, double q, const std::string& n)
        : mass(m), charge(q), name(n) {}
};

/**
 * Track state at a specific surface
 * Uses 5 parameters: (q/pT, phi, eta, d0, z0)
 * This parametrization is suitable for both barrel and endcap regions
 */
class TrackState {
public:
    // Constructor for 5-parameter track state
    TrackState(const Eigen::Vector5d& params, const Eigen::Matrix5d& cov, 
               const dd4hep::rec::Surface* surface);
    
    // Predict state to a new surface, including material effects
    TrackState predictTo(const dd4hep::rec::Surface* newSurface, 
                         const dd4hep::OverlayedField& field,
                         const dd4hep::rec::MaterialManager& matMgr,
                         const ParticleProperties& particle) const;
    
    // Filter: update state with a measurement
    TrackState update(const edm4hep::TrackerHitPlane& hit, 
                     const dd4hep::rec::Surface* surface) const;
    
    // Get chi-square increment when adding this hit
    double getChi2Increment(const edm4hep::TrackerHitPlane& hit, 
                           const dd4hep::rec::Surface* surface) const;
    
    // Accessors
    const Eigen::Vector5d& parameters() const { return _params; }
    const Eigen::Matrix5d& covariance() const { return _cov; }
    const dd4hep::rec::Surface* surface() const { return _surface; }
    
    // Get 3D momentum vector from parameters
    Eigen::Vector3d momentum() const;
    
    // Get position at closest approach to origin
    Eigen::Vector3d positionAtOrigin() const;
    
private:
    // Helper functions for material effects
    void addMultipleScatteringNoise(Eigen::Matrix5d& cov, double theta0) const;
    double calculateEnergyLoss(double momentum, double beta, double radLength, 
                              const ParticleProperties& particle) const;
    
    // Helper for projecting track parameters to measurement space
    Eigen::Matrix<double, 2, 5> projectionMatrix(const dd4hep::rec::Surface* surface) const;
    
    Eigen::Vector5d _params;           // Track parameters: (q/pT, phi, eta, d0, z0)
    Eigen::Matrix5d _cov;              // Covariance matrix
    const dd4hep::rec::Surface* _surface; // Reference surface
};

/**
 * Track class representing a reconstructed particle trajectory
 */
class Track {
public:
    Track(const TrackState& initialState);
    
    // Add a hit to the track
    void addHit(const edm4hep::TrackerHitPlane& hit, 
               const dd4hep::rec::Surface* surface,
               const TrackState& state);
    
    // Get track states
    const std::vector<TrackState>& states() const { return _states; }
    
    // Get hits associated with this track
    const std::vector<edm4hep::TrackerHitPlane>& hits() const { return _hits; }
    
    // Get surfaces associated with the hits
    const std::vector<const dd4hep::rec::Surface*>& surfaces() const { return _surfaces; }
    
    // Get chi-square of the track fit
    double chi2() const { return _chi2; }
    
    // Get number of degrees of freedom
    int ndf() const { return 2 * _hits.size() - 5; }
    
    // Momentum at a given point (e.g., vertex)
    Eigen::Vector3d momentumAt(const dd4hep::Position& point) const;
    
    // Get track parameters at first state
    Eigen::Vector5d parameters() const { 
        return _states.front().parameters(); 
    }
    
private:
    std::vector<TrackState> _states;   // Track states at each hit
    std::vector<edm4hep::TrackerHitPlane> _hits;  // Hits used in the track
    std::vector<const dd4hep::rec::Surface*> _surfaces; // Surfaces for each hit
    double _chi2;                      // Total chi-square
};

/**
 * Kalman Filter based Track Finder/Fitter implemented as a Gaudi Algorithm
 */
class KalmanTracking : public Algorithm {
public:
    KalmanTracking(const std::string& name, ISvcLocator* pSvcLocator);
    virtual ~KalmanTracking() = default;
    
    // Gaudi Algorithm methods
    virtual StatusCode initialize() override;
    virtual StatusCode execute() override;
    virtual StatusCode finalize() override;
    
    // Kalman Tracker methods
    StatusCode configure() override;
    
    // Find tracks in a collection of hits
    std::vector<Track> findTracks(const edm4hep::TrackerHitPlaneCollection* hits);
    
    // Create a track seed from three hits (triplet seeding)
    bool createTripletSeed(const edm4hep::TrackerHitPlane& hit1,
                        const edm4hep::TrackerHitPlane& hit2,
                        const edm4hep::TrackerHitPlane& hit3,
                        std::vector<Track>& tracks,
                        std::vector<bool>& usedHits,
                        size_t idx1, size_t idx2, size_t idx3);

    // Helper functions for circle fitting
    bool fitCircle(double x1, double y1, double x2, double y2, double x3, double y3, 
                double& x0, double& y0, double& radius);

    // Helper function for line fitting
    void fitLine(double x1, double y1, double x2, double y2, double x3, double y3,
                double& slope, double& intercept);

    // Fit an existing track candidate
    Track fitTrack(const std::vector<edm4hep::TrackerHitPlane>& hits,
                  const std::vector<const dd4hep::rec::Surface*>& surfaces,
                  const TrackState& seedState);
    
    // Set maximum chi-square for hit acceptance
    void setMaxChi2(double chi2) { m_maxChi2 = chi2; }
    
    // Get maximum chi-square
    double getMaxChi2() const { return m_maxChi2; }
    
    // Calculate radiation length between two points
    double getRadiationLength(const dd4hep::rec::Vector3D& start, const dd4hep::rec::Vector3D& end) const;
    
    // Find surface for a hit
    const dd4hep::rec::Surface* findSurface(const edm4hep::TrackerHitPlane& hit) const;
    
private:
    // Helper method to find surface for a given cell ID
    const dd4hep::rec::Surface* findSurfaceByID(uint64_t cellID) const;
    
    // Extract layer ID from cell ID
    int getLayerID(uint64_t cellID) const;
    
    // Find potential intersecting surfaces for a track
    std::vector<const dd4hep::rec::Surface*> findIntersectingSurfaces(
        const TrackState& state, double maxDistance) const;
    
    // Kalman filter operations
    TrackState predictStep(const TrackState& state, const dd4hep::rec::Surface* nextSurface);
    TrackState filterStep(const TrackState& predicted, 
                         const edm4hep::TrackerHitPlane& hit,
                         const dd4hep::rec::Surface* surface);
    
    // Track building operations
    std::vector<std::tuple<size_t, edm4hep::TrackerHitPlane, const dd4hep::rec::Surface*>> findCompatibleHits(
        const TrackState& state, 
        const edm4hep::TrackerHitPlaneCollection* hits,
        const std::vector<bool>& usedHits);
    
    // Group surfaces by detector layer
    std::map<int, std::vector<const dd4hep::rec::Surface*>> getSurfacesByLayer() const;
    
    // Convert internal track representation to standard track format
    void createTrack(edm4hep::TrackCollection* trackCollection, const Track& track) const;
    
    // Properties
    Gaudi::Property<std::string> m_inputHitCollection{this, "InputHitCollection", "TrackerHits", "Input hit collection path"};
    Gaudi::Property<std::string> m_outputTrackCollection{this, "OutputTrackCollection", "KalmanTracks", "Output track collection path"};
    Gaudi::Property<std::string> m_detectorName{this, "DetectorName", "Tracker", "Name of detector to process"};
    Gaudi::Property<double> m_maxChi2{this, "MaxChi2", 10.0, "Maximum chi2 for hit-track compatibility"};
    Gaudi::Property<std::string> m_particleType{this, "ParticleType", "pion", "Particle type for material effects"};
    Gaudi::Property<double> m_initialMomentum{this, "InitialMomentum", 1.0, "Initial momentum estimate for seeding (GeV)"};
    Gaudi::Property<double> m_maxDistanceToSurface{this, "MaxDistanceToSurface", 10.0, "Maximum distance to consider surface intersections (mm)"};
    Gaudi::Property<std::string> m_encodingStringParameter{this, "EncodingStringParameterName", "GlobalTrackerReadoutID", "Name of DD4hep parameter with the encoding string"};
    Gaudi::Property<double> m_maxRadius{this, "MaxRadius", 100000.0, "Maximum radius for track curvature (mm)"};

    // Services
    ServiceHandle<IGeoSvc> m_geoSvc{this, "GeoSvc", "GeoSvc", "Detector geometry service"};
    
    // Member variables
    dd4hep::Detector* m_detector{nullptr};  // Detector instance
    dd4hep::OverlayedField m_field;         // Magnetic field
    dd4hep::rec::MaterialManager* m_materialManager{nullptr}; // Material manager
    const dd4hep::rec::SurfaceMap* m_surfaceMap{nullptr};    // Surface map from detector
    std::vector<const dd4hep::rec::Surface*> m_surfaces;      // All detector surfaces
    std::map<int, std::vector<const dd4hep::rec::Surface*>> m_surfacesByLayer; // Surfaces grouped by layer
    ParticleProperties m_particleProperties;  // Particle properties for material effects
    std::unique_ptr<dd4hep::DDSegmentation::BitFieldCoder> m_bitFieldCoder; // For cell ID decoding
    
    // Map of known particle types
    std::map<std::string, ParticleProperties> m_particleMap;
};
DECLARE_COMPONENT(KalmanTracking)
#endif // KALMAN_TRACKER_H
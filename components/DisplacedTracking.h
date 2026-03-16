#ifndef DISPLACED_TRACKING_H
#define DISPLACED_TRACKING_H

#include <vector>
#include <map>
#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <utility>

// POSIX headers for stderr capture (GenFit Cholesky failure detection)
#include <unistd.h>
#include <fcntl.h>

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
#include "edm4hep/MCParticleCollection.h"

// Interface includes
#include "k4Interface/IGeoSvc.h"

// GenFit includes
#include "Track.h"
#include "TrackCand.h"
#include "AbsTrackRep.h"
#include "RKTrackRep.h"
//#include "GeaneTrackRep.h"
#include "KalmanFitter.h"
#include "KalmanFitterRefTrack.h"
#include <DAF.h>
#include "KalmanFitterInfo.h"
#include "MeasuredStateOnPlane.h"
#include "PlanarMeasurement.h"
#include "SpacepointMeasurement.h"
#include "MaterialEffects.h"
#include "TGeoMaterialInterface.h"
#include "FieldManager.h"
#include "AbsBField.h"
#include "FitStatus.h"
#include <EventDisplay.h>
#include <unordered_map>
#include <unordered_set>
#include <atomic>
#include <iomanip>

// Include our adapter classes
#include "GenFitAdapters.h"
#include "GenfitField.hpp"
#include "GenfitMaterialInterface.h"

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
 * Result of inner solenoid propagation
 */
struct InnerTrajectory {
    bool success;
    Eigen::Vector3d finalPosition;      // cm
    Eigen::Vector3d finalMomentum;      // GeV/c
    double arcLength;                   // cm, total path length
    int numSteps;                       // total steps taken
    std::string message;                // status/error message
};

/**
 * Displaced Track Finder/Fitter implemented as a Key4hep Transformer
 */
// Convenience alias for the link-propagation callback
using LinkPropagator = std::function<void(const edm4hep::MutableTrackerHitPlane&, const edm4hep::TrackerHitPlane&)>;

class DisplacedTracking final
    : public k4FWCore::MultiTransformer<
        std::tuple<edm4hep::TrackCollection,
                   edm4hep::TrackerHitPlaneCollection,
                   edm4hep::TrackerHitSimTrackerHitLinkCollection>(
                       const edm4hep::TrackerHitPlaneCollection&,
                       const edm4hep::EventHeaderCollection&,
                       const edm4hep::TrackerHitSimTrackerHitLinkCollection&)
      > {
public:
    DisplacedTracking(const std::string& name, ISvcLocator* pSvcLocator);
    virtual ~DisplacedTracking() = default;

    StatusCode initialize() override;
    StatusCode finalize() override;

    std::tuple<edm4hep::TrackCollection,
               edm4hep::TrackerHitPlaneCollection,
               edm4hep::TrackerHitSimTrackerHitLinkCollection> operator()(
        const edm4hep::TrackerHitPlaneCollection& hits,
        const edm4hep::EventHeaderCollection& headers,
        const edm4hep::TrackerHitSimTrackerHitLinkCollection& recoSimLinks) const override;
    
    // Make sure findTracks is const
    void findTracks(
                            const edm4hep::TrackerHitPlaneCollection* hits,
                            edm4hep::TrackCollection& outputTracks,
                            edm4hep::TrackerHitPlaneCollection& outputHits,
                            const edm4hep::TrackerHitSimTrackerHitLinkCollection& recoSimLinks,
                            edm4hep::TrackerHitSimTrackerHitLinkCollection& outputLinks,
                            bool& evtUsedPtInner,
                            bool& evtUsedPtInner2,
                            bool& evtUsedPtFallback) const;

    // Calculate circle center and radius using the Direct Formula Method
    bool calculateCircleCenterDirect(
        double x1, double y1, double x2, double y2, double x3, double y3,
        double& x0, double& y0, double& radius) const;
    // calculate sagitta    
    double calculateSagitta(const Eigen::Vector3d& p1, 
                            const Eigen::Vector3d& p2, 
                            const Eigen::Vector3d& p3) const; 
    // Calculate circle center and radius using the Sagitta Method
    bool calculateCircleCenterSagitta(
        const Eigen::Vector3d& p1, const Eigen::Vector3d& p2, const Eigen::Vector3d& p3,
        double& x0, double& y0, double& radius) const;
    // Calculate impact parameter
    double calculateImpactParameter(
        double x0, double y0, double radius, bool clockwise,
        double innerFieldStrength, double outerFieldStrength,
        const Eigen::Vector3d& p1, const Eigen::Vector3d& p2, const Eigen::Vector3d& p3) const;    

    // Create a track seed from three hits (triplet seeding)
    bool createTripletSeed(const edm4hep::TrackerHitPlane& hit1,
                        const edm4hep::TrackerHitPlane& hit2,
                        const edm4hep::TrackerHitPlane& hit3,
                        edm4hep::TrackCollection* tracks,
                        edm4hep::TrackerHitPlaneCollection& outputHits,
                        const LinkPropagator& propagateLink,
                        std::vector<bool>& usedHits,
                        size_t idx1, size_t idx2, size_t idx3) const;

    // Helper functions for circle fitting
    bool fitCircle(double x1, double y1, double x2, double y2, double x3, double y3, 
                double& x0, double& y0, double& radius) const;

    // Helper function for line fitting
    void fitLine(double x1, double y1, double x2, double y2, double x3, double y3,
                double& slope, double& intercept) const;

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
    // Extract type ID from cell ID
    int getTypeID(uint64_t cellID) const;
    // Extract layer ID from cell ID
    int getLayerID(uint64_t cellID) const;
    // Compute composite layer ID: barrel->layerID, +endcap->1000+layerID, -endcap->-1000-layerID
    int getCompositeID(uint64_t cellID) const;
    // Get PT
    double getPT(const edm4hep::TrackState& state) const;
    // Get position vector from track state
    Eigen::Vector3d getPosition(const edm4hep::TrackState& state) const;
    
    // Create a track state with the given parameters
    edm4hep::TrackState createTrackState(
                    double d0, double phi, double omega, double z0, double tanLambda,
                    int location = edm4hep::TrackState::AtOther) const;
        
    // Predict track state to a new surface
    edm4hep::TrackState predictToSurface(
        const edm4hep::TrackState& state,
        const dd4hep::rec::Surface* surface) const;
        
    // Calculate chi-square between track state and hit
    double getChi2(
        const edm4hep::TrackState& state,
        const edm4hep::TrackerHitPlane& hit,
        const dd4hep::rec::Surface* surface) const;

    // Group surfaces by detector layer
    std::map<int, std::vector<const dd4hep::rec::Surface*>> getSurfacesByLayer() const;

    // Looking for compatible hits in the +3 layer (initially extra 4th layer)
    bool findCompatibleExtraHit(
        std::vector<edm4hep::TrackerHitPlane>& trackHits,
        const edm4hep::TrackerHitPlaneCollection* allHits,
        std::vector<bool>& usedHits,
        const Eigen::Vector3d& seedDirection,
        double coneHalfAngleDeg = 35.0) const;
        
    // Fit circle to 4 hits using least squares and calculate chi2
    // fitCov3x3 is the 3×3 covariance matrix on [x0, y0, R] in cm² from the WLS fit.
    // It is set on success and used by the caller to propagate uncertainties to EDM4hep track params.
    bool fitCircleToFourHits(const std::vector<edm4hep::TrackerHitPlane>& hits,
                            double& x0, double& y0, double& radius, double& chi2,
                            Eigen::Matrix3d& fitCov3x3) const;

    // fitTrackWithGenFit with extrapolation support
    bool fitTrackWithGenFit(
        const std::vector<edm4hep::TrackerHitPlane>& hits,
        const edm4hep::TrackState& seedState,
        edm4hep::MutableTrack& finalTrack,
        edm4hep::TrackerHitPlaneCollection& outputHits,
        const LinkPropagator& propagateLink,
        const edm4hep::TrackerHitPlaneCollection* allHits = nullptr) const;
    
    // Convert EDM4hep track state to GenFit track representation
    genfit::MeasuredStateOnPlane convertToGenFitState(
        const edm4hep::TrackState& state,
        genfit::AbsTrackRep* rep) const;
    
    // Convert GenFit state back to EDM4hep track state
    edm4hep::TrackState convertToEDM4hepState(
        const genfit::MeasuredStateOnPlane& state,
        int location = edm4hep::TrackState::AtOther) const;
    
    // Create a GenFit measurement from an EDM4hep hit
    genfit::AbsMeasurement* createGenFitMeasurement(
        const edm4hep::TrackerHitPlane& hit,
        const dd4hep::rec::Surface* surface,
        int hitId,
        genfit::TrackPoint* trackPoint) const;

    // Add field and material adapters as members
    std::unique_ptr<DD4hepFieldAdapter> m_genFitField;
    std::unique_ptr<DD4hepMaterialAdapter> m_genFitMaterial;
    
    // ==================== TRUTH MATCHING ====================
    // Data handle for truth information - use mutable to allow get() in const methods
    mutable k4FWCore::DataHandle<edm4hep::MCParticleCollection> m_mcParticles{
        "MCParticles", Gaudi::DataHandle::Reader, this};

    // Helper method for truth matching
    const edm4hep::MCParticle* findTruthParticle(
        const edm4hep::Track& recoTrack,
        const edm4hep::MCParticleCollection& mcParticles) const;

    // Validate track against truth
    void validateTrackWithTruth(
        const edm4hep::Track& recoTrack,
        const edm4hep::MCParticleCollection& mcParticles,
        int trackNumber) const;
    // ==================== INNER SOLENOID PROPAGATION ====================
    /**
     * Propagate outer track inward to solenoid using RK4
     */
    InnerTrajectory propagateToInner(
        const Eigen::Vector3d& outerPos,    // Starting position in cm
        const Eigen::Vector3d& outerMom,    // Starting momentum in GeV/c
        double targetRadius = 100.0         // Stop at this R in cm
    ) const;
    
    /**
     * Check if position is inside solenoid boundary
     */
    bool isInsideSolenoid(const Eigen::Vector3d& pos) const;
    
    /**
     * Single RK4 propagation step
     */
    void rk4PropagationStep(
        Eigen::Vector3d& pos,      // Position (modified in place)
        Eigen::Vector3d& mom,      // Momentum (modified in place)
        double stepSize            // Step size in cm (negative = inward)
    ) const;
    
    // Inner propagation constants
    static constexpr double SOLENOID_RADIUS_CM = 230.0;  // shouldn't be hardcoded fix
    static constexpr double SOLENOID_Z_HALF_CM = 218.0;  // shouldn't be hardcoded fix
    // =====================================================================
    
    // Properties
    Gaudi::Property<std::string> m_detectorName{this, "DetectorName", "Tracker", "Name of detector to process"};
    Gaudi::Property<double> m_maxChi2{this, "MaxChi2", 10.0, "Maximum chi2 for hit-track compatibility"};
    Gaudi::Property<std::string> m_particleType{this, "ParticleType", "pion", "Particle type for material effects"};
    Gaudi::Property<std::string> m_encodingStringParameter{this, "EncodingStringParameterName", "GlobalTrackerReadoutID", "Name of DD4hep parameter with the encoding string"};
    Gaudi::Property<double> m_maxDist{this, "MaxDist", 150.0, "Maximum distance between two hits (cm)"};
    // GenFit properties
    Gaudi::Property<int> m_maxFitIterations{this, "MaxFitIterations", 4, "Maximum iterations for track fitting"};
    Gaudi::Property<bool> m_useGenFit{this, "UseGenFit", true, "Use GenFit for track fitting"};
    Gaudi::Property<int> m_debugLevel{this, "DebugLevel", 0, "Debug level of GenFit"};
    // Inner propagation steering
    Gaudi::Property<bool> m_doInnerPropagation{this, "DoInnerPropagation", true,
        "Propagate outer track inward through the solenoid boundary using RK4 (saves state at AtVertex)"};
    Gaudi::Property<double> m_innerPropTargetRadius{this, "InnerPropTargetRadius", 100.0,
        "Target radius in cm for inner RK4 propagation"};
    Gaudi::Property<double> m_maxSeedPT{this, "MaxSeedPT", 200.0,
        "Maximum allowed seed pT in GeV/c. Triplets above this threshold are rejected unless they are the only option, in which case the lowest-pT candidate is kept."};

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

    int getPDGCode() const;
    
    // Map of known particle types
    std::map<std::string, ParticleProperties> m_particleMap;

    // ==================== RUN STATISTICS ====================
    mutable std::atomic<int> m_statTotalEvents{0};
    mutable std::atomic<int> m_statTracksReconstructed{0};
    mutable std::atomic<int> m_statPositiveCharge{0};
    mutable std::atomic<int> m_statNegativeCharge{0};
    mutable std::atomic<int> m_statInnerSeedUsed{0};
    mutable std::atomic<int> m_statFallbackSeedUsed{0};
    mutable std::atomic<int> m_statHighestPtInner{0};
    mutable std::atomic<int> m_statHighestPtInner2{0};
    mutable std::atomic<int> m_statHighestPtFallback{0};
    mutable std::atomic<int> m_statTotalTripletCombos{0};
    mutable std::atomic<int> m_statValidTriplets{0};
    mutable std::atomic<int> m_statFourHitTracks{0};
    mutable std::atomic<int> m_statThreeHitTracks{0};
    mutable std::atomic<int> m_statStateAtFirstHit{0};
    mutable std::atomic<int> m_statStateAtLastHit{0};
    mutable std::atomic<int> m_statStateAtCalorimeter{0};
    mutable std::atomic<int> m_statStateAtOther{0};
    mutable std::atomic<int> m_statStateAtVertex{0};
    mutable std::atomic<int> m_statInnerPropSuccess{0};
    mutable std::atomic<int> m_statGenFitSuccess{0};

    // GenFit detailed counters
    mutable std::atomic<int>   m_statGenFitFailed{0};        // isFitted==false after processTrack
    mutable std::atomic<int>   m_statGenFitIllCond{0};       // ill-conditioned covariance
    mutable std::atomic<int>   m_statGenFitException{0};     // genfit::Exception thrown
    mutable std::atomic<int>   m_statGenFitChi2Count{0};     // number of valid chi2/ndf values
    mutable std::atomic<int>   m_statGenFitChi2BadCount{0};  // fits with chi2/ndf > 1000 (non-physical)
    // chi2 sum stored as integer × 1000 to avoid non-atomic float (multiply by 1000 before storing)
    mutable std::atomic<long long> m_statGenFitChi2Sum{0};   // sum of (chi2/ndf * 1000)

    // Triplet quality: flag bad geometry (direct vs sagitta radius disagreement > 2×)
    mutable std::atomic<int> m_statTripletBadGeom{0};
    // Triplets rejected by MaxSeedPT threshold (fallback used when all above threshold)
    mutable std::atomic<int> m_statTripletsCutByPT{0};
    // =========================================================

};

DECLARE_COMPONENT(DisplacedTracking)
#endif // DISPLACED_TRACKING_H
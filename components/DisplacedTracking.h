#ifndef DISPLACED_TRACKING_H
#define DISPLACED_TRACKING_H

#include <vector>
#include <map>
#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <utility>

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

// GenFit includes
#include "Track.h"
#include "TrackCand.h"
#include "AbsTrackRep.h"
#include "RKTrackRep.h"
#include "KalmanFitter.h"
#include "KalmanFitterRefTrack.h"
#include "KalmanFitterInfo.h"
#include "MeasuredStateOnPlane.h"
#include "PlanarMeasurement.h"
#include "SpacepointMeasurement.h"
#include "MaterialEffects.h"
#include "TGeoMaterialInterface.h"
#include "FieldManager.h"
#include "AbsBField.h"
#include "FitStatus.h"

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
 * Kalman Filter based Track Finder/Fitter implemented as a Key4hep Transformer
 */
class DisplacedTracking final
    : public k4FWCore::MultiTransformer<
        std::tuple<edm4hep::TrackCollection>(const edm4hep::TrackerHitPlaneCollection&, const edm4hep::EventHeaderCollection&)
      > {
public:
    DisplacedTracking(const std::string& name, ISvcLocator* pSvcLocator);
    virtual ~DisplacedTracking() = default;

    StatusCode initialize() override;
    StatusCode finalize() override;

    std::tuple<edm4hep::TrackCollection> operator()(
        const edm4hep::TrackerHitPlaneCollection& hits,
        const edm4hep::EventHeaderCollection& headers) const override;
    
    // Make sure findTracks is const
    edm4hep::TrackCollection findTracks(
                            const edm4hep::TrackerHitPlaneCollection* hits) const;

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
    genfit::MaterialEffects* materialEffects;
    genfit::FieldManager* fieldManager;
    GenfitMaterialInterface* m_geoMaterial;
    GenfitField* m_genfitField;
    // Helper method to find surface for a given cell ID
    const dd4hep::rec::Surface* findSurfaceByID(uint64_t cellID) const;
    // Extract type ID from cell ID
    int getTypeID(uint64_t cellID) const;
    // Extract layer ID from cell ID
    int getLayerID(uint64_t cellID) const;
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

    // Extrapolate a state to a detector surface
    genfit::MeasuredStateOnPlane extrapolateToSurface(
        const genfit::MeasuredStateOnPlane& state,
        genfit::AbsTrackRep* rep,
        const dd4hep::rec::Surface* targetSurface) const;

    // Find hits compatible with a track state
    std::vector<std::pair<size_t, edm4hep::TrackerHitPlane>> findCompatibleHitsInLayer(
        const genfit::MeasuredStateOnPlane& state,
        genfit::AbsTrackRep* rep,
        int layerID,
        const edm4hep::TrackerHitPlaneCollection& allHits,
        double maxChi2) const;

    // fitTrackWithGenFit with extrapolation support
    bool fitTrackWithGenFit(
        const std::vector<edm4hep::TrackerHitPlane>& hits,
        const edm4hep::TrackState& seedState,
        edm4hep::MutableTrack& finalTrack,
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
    
    // Properties
    Gaudi::Property<std::string> m_detectorName{this, "DetectorName", "Tracker", "Name of detector to process"};
    Gaudi::Property<double> m_maxChi2{this, "MaxChi2", 10.0, "Maximum chi2 for hit-track compatibility"};
    Gaudi::Property<std::string> m_particleType{this, "ParticleType", "pion", "Particle type for material effects"};
    Gaudi::Property<std::string> m_encodingStringParameter{this, "EncodingStringParameterName", "GlobalTrackerReadoutID", "Name of DD4hep parameter with the encoding string"};
    Gaudi::Property<double> m_maxDist{this, "MaxDist", 150.0, "Maximum distance between two hits (cm)"};
    // GenFit properties
    Gaudi::Property<int> m_maxFitIterations{this, "MaxFitIterations", 4, "Maximum iterations for track fitting"};
    Gaudi::Property<bool> m_useGenFit{this, "UseGenFit", true, "Use GenFit for track fitting"};

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

    edm4hep::TrackState getSeedState(const edm4hep::Track& track) const;
    std::vector<edm4hep::TrackerHitPlane> getOriginalHits(
        const edm4hep::Track& track,
        const edm4hep::TrackerHitPlaneCollection& allHits) const;
    int getHighestLayerID(const std::vector<edm4hep::TrackerHitPlane>& hits) const;
    void copyTrackParameters(const edm4hep::Track& source, edm4hep::MutableTrack& target) const;
    int getPDGCode() const;
    
    // Map of known particle types
    std::map<std::string, ParticleProperties> m_particleMap;
};

DECLARE_COMPONENT(DisplacedTracking)
#endif // KALMAN_TRACKER_H
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

    // Calculate impact parameter
    double calculateImpactParameter(
        double x0, double y0, double radius, bool clockwise,
        double innerFieldStrength, double outerFieldStrength,
        const Eigen::Vector3d& p1, const Eigen::Vector3d& p2, const Eigen::Vector3d& p3) const;    

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

    // Fit circle to N≥3 hits using weighted least-squares + Gauss-Newton refinement.
    // fitCov3x3 is the 3×3 covariance on [x0,y0,R] (cm²).
    // residualsCm returns the per-hit geometric residual (cm) — used for outlier rejection.
    // Returns false if the fit is degenerate or NDF<1.
    bool fitCircleNHits(const std::vector<edm4hep::TrackerHitPlane>& hits,
                        double& x0, double& y0, double& radius, double& chi2ndf,
                        Eigen::Matrix3d& fitCov3x3,
                        std::vector<double>& residualsCm) const;

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
        double targetRadius = 100.0,        // Stop at this R in cm
        double chargeSign   = 1.0           // +1 or -1 particle charge sign
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
        double stepSize,           // Step size in cm
        double chargeSign = 1.0    // +1 or -1 particle charge sign
    ) const;

    /**
     * Analytical two-segment inner propagation.
     * Uses the outer circle (cx_o, cy_o, R_outer) and the field-reversal geometry
     * to construct the inner helix and save it as TrackState::AtIP.
     *
     * Parameters (all in cm unless noted):
     *   cx_o, cy_o, R_outer  – outer circle center and radius
     *   omega_outer_mm       – outer omega (1/mm, sign encodes curvature direction)
     *   pT_GeV               – transverse momentum (GeV/c, conserved across boundary)
     *   tanLambda            – dz/d(arc length), conserved
     *   x_ref, y_ref, z_ref  – innermost muon-system hit (cm); used for intersection
     *                          selection and z anchor
     *   finalTrack           – track to which the AtIP state is appended
     * Returns true on success.
     */
    bool analyticalInnerPropagation(
        double cx_o_cm, double cy_o_cm, double R_outer_cm,
        double omega_outer_mm,
        double pT_GeV,
        double tanLambda,
        double x_ref_cm, double y_ref_cm, double z_ref_cm,
        bool   isEndcap,             // true = endcap hit → use z-cap crossing
        edm4hep::MutableTrack& finalTrack) const;

    /**
     * Scan B_z radially to locate the exact solenoid boundary (sign-change radius).
     * Starts from the SolenoidRadius property, walks in 1-mm steps to bracket the
     * zero, then binary-searches to < 0.01 mm precision.
     * Returns the boundary radius in mm.
     */
    double findSolenoidBoundary() const;

    /**
     * Scan B_z along the z-axis (x=y=0) to locate the exact solenoid endcap
     * boundary (sign-change |z|). Starts from SolenoidHalfZ property, same
     * algorithm as findSolenoidBoundary(). Returns the half-length in mm.
     */
    double findSolenoidBoundaryZ() const;

    /**
     * Print a compact truth-comparison table for all track states
     * (AtLastHit, AtVertex/RK4, AtIP/analytical) vs. the best-matched MC particle.
     * For a prompt muon gun AtIP should show d0 ≈ 0, z0 ≈ 0.
     */
    void compareStatesWithTruth(
        const edm4hep::Track& track,
        const edm4hep::MCParticleCollection& mcParticles,
        int trackNumber) const;

    // Inner propagation constants
    static constexpr double SOLENOID_RADIUS_CM = 230.0;  // fallback default
    static constexpr double SOLENOID_Z_HALF_CM = 218.0;
    double m_solenoidBoundaryMm{2300.0};  // runtime value set by findSolenoidBoundary()
    double m_solenoidHalfZmm{2180.0};     // runtime value set by findSolenoidBoundaryZ()
    // =====================================================================
    
    // Properties
    Gaudi::Property<std::string> m_detectorName{this, "DetectorName", "Tracker", "Name of detector to process"};
    Gaudi::Property<std::string> m_particleType{this, "ParticleType", "pion", "Particle type for material effects"};
    Gaudi::Property<std::string> m_encodingStringParameter{this, "EncodingStringParameterName", "GlobalTrackerReadoutID", "Name of DD4hep parameter with the encoding string"};
    // GenFit properties
    Gaudi::Property<int> m_maxFitIterations{this, "MaxFitIterations", 4, "Maximum iterations for track fitting"};
    Gaudi::Property<bool> m_useGenFit{this, "UseGenFit", true, "Use GenFit for track fitting"};
    Gaudi::Property<int> m_debugLevel{this, "DebugLevel", 0, "Debug level of GenFit"};
    // Solenoid geometry
    Gaudi::Property<double> m_solenoidRadius{this, "SolenoidRadius", 2300.0,
        "Nominal solenoid inner radius in mm; start for Bz sign-change scan in initialize()"};
    Gaudi::Property<double> m_solenoidHalfZ{this, "SolenoidHalfZ", 2180.0,
        "Nominal solenoid half-length in mm; start for Bz sign-change scan along z-axis in initialize()"};
    // Inner propagation steering
    Gaudi::Property<bool> m_doInnerPropagation{this, "DoInnerPropagation", true,
        "Propagate outer track inward through the solenoid boundary using RK4 (saves state at AtVertex)"};
    Gaudi::Property<double> m_innerPropTargetRadius{this, "InnerPropTargetRadius", 100.0,
        "Target radius in cm for inner RK4 propagation"};
    Gaudi::Property<double> m_maxComboPT{this, "MaxComboPT", 200.0,
        "Maximum allowed combo pT in GeV/c. Combinations whose circle-fit radius implies pT above this threshold are excluded during the combinatorial search."};


    // ── N-hit combinatorial strategy properties ─────────────────────────────
    Gaudi::Property<int> m_minTrackHits{this, "MinTrackHits", 3,
        "Minimum number of hits required to form a track."};
    Gaudi::Property<int> m_maxCombinatorialHits{this, "MaxCombinatorialHits", 7,
        "Maximum number of hits considered per combination (caps C(N,k) combinatorics). "
        "Increase cautiously — C(10,7)=120 is fine, C(20,7)=77520 may be slow."};
    Gaudi::Property<double> m_nHitMaxChi2NDF{this, "NHitMaxChi2NDF", 5.0,
        "Maximum chi2/NDF of the N-hit circle fit to accept a combination."};
    Gaudi::Property<double> m_outlierSigma{this, "OutlierSigma", 5.0,
        "Outlier rejection threshold in units of SigmaHitDefault. "
        "Hits with residual > OutlierSigma * SigmaHitDefault are iteratively removed."};
    Gaudi::Property<double> m_sigmaHitDefault{this, "SigmaHitDefault", 0.04,
        "Default hit position resolution in cm (0.04 cm = 0.4 mm) used when hit du/dv is not set."};
    Gaudi::Property<double> m_hitIsolationCut{this, "HitIsolationCut", 2.0,
        "Minimum distance in cm between a candidate hit and any other unused hit not in the "
        "combination. Set to 0 to disable. Rejects hits inside secondary shower clusters."};
    Gaudi::Property<int> m_minLayerSpan{this, "MinLayerSpan", 3,
        "Minimum number of distinct composite detector layers that must be covered by a combination."};
    Gaudi::Property<int> m_maxOutlierIterations{this, "MaxOutlierIterations", 3,
        "Maximum rounds of outlier removal per track."};
    Gaudi::Property<double> m_maxConsecDeltaPhi{this, "MaxConsecDeltaPhi", 0.5,
        "Maximum |Δφ| in radians between consecutive (adjacent-layer) hits in a combination. "
        "Rejects cross-track combos where the azimuthal step between layers is too large. "
        "Set to -1 to disable."};
    Gaudi::Property<double> m_maxComboPhiSpread{this, "MaxComboPhiSpread", 1.0,
        "Maximum total φ spread in radians across all hits in a combination (max - min φ, "
        "relative to the innermost hit). Rejects combos spanning more than this azimuthal window. "
        "Set to -1 to disable."};
    Gaudi::Property<double> m_maxConsecHitDist{this, "MaxConsecutiveHitDistance", 100.0,
        "Maximum 3D distance in cm between consecutive (inner→outer ordered) hits in a combination. "
        "Rejects combos where adjacent-layer hits are implausibly far apart (e.g. hits from "
        "different physical tracks combined across a large gap). Set to -1 to disable."};
    Gaudi::Property<double> m_maxPairHitDist{this, "MaxPairHitDistance", 300.0,
        "Maximum 3D distance in cm between any pair of hits in a combination. "
        "Catches ghost combos that span widely separated detector regions. "
        "Set to -1 to disable. 300 cm accommodates barrel→endcap transition tracks."};
    Gaudi::Property<double> m_maxComboTanLambdaDev{this, "MaxComboTanLambdaDev", 0.5,
        "Maximum deviation of any consecutive-pair tanλ (= Δz / chord_xy) from the "
        "median tanλ of all pairs in a combination. Rejects combos where hits from "
        "different tracks or delta-ray loops produce wildly inconsistent dip angles. "
        "A real track has nearly constant tanλ; cross-particle fakes deviate by O(1). "
        "Default 0.5 (≈27° dip-angle spread). Set to -1 to disable."};
    Gaudi::Property<double> m_maxHitEdepKeV{this, "MaxHitEdepKeV", 10.0,
        "Maximum energy deposit in keV for a hit to be used in tracking. "
        "Hits above this threshold are excluded before combo enumeration. "
        "Rejects high-edep delta-rays and hadronic secondaries. Set to -1 to disable."};
    Gaudi::Property<double> m_proxScoreWeight{this, "ProximityScoreWeight", 1.0,
        "Weight for edep-homogeneity scoring when ranking same-k combos. "
        "avgEdepDev = mean |edep(hit) - meanEdepSingle| across combo hits, where "
        "meanEdepSingle is the event-level muon edep baseline from single-hit layers. "
        "For k=3: tiebreaker on equal radius. "
        "For k>=4: combined score = chi2/ndf + weight * avgEdepDev / EdepNormKeV "
        "(active criterion — not just tiebreaker). Lower combined score wins. "
        "Set to 0 to disable edep scoring."};
    Gaudi::Property<double> m_edepNormKeV{this, "EdepNormKeV", 2.0,
        "Normalisation scale (keV) for converting avgEdepDev into chi2-equivalent units "
        "when computing the combined combo score for k>=4. "
        "combined_score = chi2/ndf + ProximityScoreWeight * avgEdepDev / EdepNormKeV. "
        "A value of 2.0 means a 2 keV average edep deviation adds 1.0 to the effective chi2. "
        "Tune this to balance geometry (chi2) vs. edep-homogeneity discrimination."};
    Gaudi::Property<int> m_maxHitsPerComposite{this, "MaxHitsPerCompositeID", 0,
        "If > 0, at most this many hits per compositeID are kept before combo enumeration, "
        "ranked by |edep - mean_single_layer_edep| (lowest deviation first). "
        "Reduces combinatorial explosion in high-occupancy layers. "
        "Set to 0 to use all hits (default)."};
    Gaudi::Property<bool> m_doCrowdedLayerHitSelection{this, "doCrowdedLayerHitSelection", false,
        "Guided crowded-layer hit selection"};
    Gaudi::Property<bool> m_doHitClustering{this, "DoHitClustering", false,
        "Pre-cluster unused hits into spatially coherent road candidates before the combinatorial "
        "search. Uses Union-Find on adjacent compositeID layers with the same MaxConsecDeltaPhi and "
        "MaxConsecutiveHitDistance thresholds. Each cluster becomes an independent hit-pool, "
        "preventing hits from unrelated track regions from ever entering the same combo. "
        "Clusters spanning fewer than MinTrackHits distinct layers are skipped. "
        "Default false (original single-pool behaviour). Enable for high-occupancy events "
        "where delta-ray or multi-track hit pollution inflates the Cartesian product."};

    // ── Neighbour-track quality gate ──────────────────────────────────────
    Gaudi::Property<double> m_neighbourTrackMaxDist{this, "NeighbourTrackMaxDist", 5.0,
        "If the best-combo candidate has any hit within this distance (cm) of a hit from "
        "an already-accepted track, apply stricter quality requirements defined by "
        "NeighbourTrackMinHits and NeighbourTrackMaxChi2NDF. "
        "Rejects delta-ray shadow tracks that form near a parent muon track. "
        "Default 5.0 cm = 50 mm. Set to 0 to disable."};
    Gaudi::Property<int> m_neighbourTrackMinHits{this, "NeighbourTrackMinHits", 4,
        "Minimum hit count required for a combo flagged as a neighbour track "
        "(any hit lies within NeighbourTrackMaxDist of an accepted track hit). "
        "Stronger than MinTrackHits — forces ≥4 hits to reject 3-hit shadow fakes."};
    Gaudi::Property<double> m_neighbourTrackMaxChi2NDF{this, "NeighbourTrackMaxChi2NDF", 3.0,
        "Maximum chi2/NDF accepted for a neighbour-track combo. "
        "Applied together with NeighbourTrackMinHits when a combo is spatially near "
        "an already-accepted track."};
    // ────────────────────────────────────────────────────────────────────────

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
    // General
    mutable std::atomic<int> m_statTotalEvents{0};
    mutable std::atomic<int> m_statTracksReconstructed{0};

    // Charge sign
    mutable std::atomic<int> m_statPositiveCharge{0};
    mutable std::atomic<int> m_statNegativeCharge{0};
    mutable std::atomic<int> m_statChargeUndetermined{0};  // |sinBend| below reliability threshold

    // N-hit combinatorial search
    mutable std::atomic<int> m_statTotalTripletCombos{0};  // total C(N,k) combinations evaluated
    mutable std::atomic<int> m_statValidTriplets{0};        // combinations that produced a track
    mutable std::atomic<int> m_statDupCompositeRejected{0};  // rejected: ≥2 hits share a compositeID
    mutable std::atomic<int> m_statLayerSpanRejected{0};    // rejected by MinLayerSpan cut
    mutable std::atomic<int> m_statIsolationRejected{0};    // rejected by HitIsolationCut
    mutable std::atomic<int> m_statConsecPhiRejected{0};    // rejected by MaxConsecDeltaPhi
    mutable std::atomic<int> m_statPhiSpreadRejected{0};    // rejected by MaxComboPhiSpread
    mutable std::atomic<int> m_statConsecDistRejected{0};   // rejected by MaxConsecutiveHitDistance
    mutable std::atomic<int> m_statPairDistRejected{0};     // rejected by MaxPairHitDistance
    mutable std::atomic<int> m_statEdepRejected{0};         // hits removed by MaxHitEdepKeV pre-filter
    mutable std::atomic<int> m_statOutlierHitsRemoved{0};   // total hits removed by outlier rejection
    mutable std::atomic<int> m_statNeighbourRejected{0};    // combos rejected by neighbour-track gate
    mutable std::atomic<int> m_statComboPTRejected{0};      // combos rejected by MaxComboPT threshold
    mutable std::atomic<int> m_statTanLambdaRejected{0};   // combos rejected by MaxComboTanLambdaDev
    mutable std::atomic<int> m_statClustersBuilt{0};        // total clusters built across all iterations

    // Hit multiplicity distribution
    mutable std::atomic<int> m_statNhit3{0};
    mutable std::atomic<int> m_statNhit4{0};
    mutable std::atomic<int> m_statNhit5{0};
    mutable std::atomic<int> m_statNhit6{0};
    mutable std::atomic<int> m_statNhit7plus{0};

    // Track states
    mutable std::atomic<int> m_statStateAtLastHit{0};
    mutable std::atomic<int> m_statStateAtOther{0};
    mutable std::atomic<int> m_statStateAtVertex{0};
    mutable std::atomic<int> m_statStateAtIP{0};
    mutable std::atomic<int> m_statInnerPropSuccess{0};
    mutable std::atomic<int> m_statAnalyticalPropSuccess{0};
    mutable std::atomic<int> m_statGenFitSuccess{0};

    // GenFit detailed counters
    mutable std::atomic<int>       m_statGenFitFailed{0};
    mutable std::atomic<int>       m_statGenFitIllCond{0};
    mutable std::atomic<int>       m_statGenFitException{0};
    mutable std::atomic<int>       m_statGenFitChi2Count{0};
    mutable std::atomic<int>       m_statGenFitChi2BadCount{0};
    mutable std::atomic<long long> m_statGenFitChi2Sum{0};

    // Legacy counters (kept to avoid reference errors; no longer incremented)
    mutable std::atomic<int> m_statTripletBadGeom{0};
    mutable std::atomic<int> m_statAngleGuardRejected{0};
    mutable std::atomic<int> m_statHighestPtInner{0};
    mutable std::atomic<int> m_statHighestPtInner2{0};
    mutable std::atomic<int> m_statHighestPtFallback{0};
    mutable std::atomic<int> m_statFourHitTracks{0};
    mutable std::atomic<int> m_statThreeHitTracks{0};
    // =========================================================

};

DECLARE_COMPONENT(DisplacedTracking)
#endif // DISPLACED_TRACKING_H
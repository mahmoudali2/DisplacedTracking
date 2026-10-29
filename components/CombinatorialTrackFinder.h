#ifndef COMBINATORIAL_TRACK_FINDER_H
#define COMBINATORIAL_TRACK_FINDER_H

#include <vector>
#include <map>
#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <utility>

// stderr capture (detects GenFit Cholesky failures)
#include <unistd.h>
#include <fcntl.h>

#include "GaudiKernel/Algorithm.h"
#include "Gaudi/Property.h"
#include "GaudiKernel/ServiceHandle.h"
#include "GAUDI_VERSION.h"

#include <Eigen/Dense>
#include "LayerGeometry.h"
#include "SeedFinder.h"
#include "DD4hep/Detector.h"
#include "DD4hep/Fields.h"
#include "DD4hep/DD4hepUnits.h"
#include "DDRec/Surface.h"
#include "DDRec/SurfaceManager.h"
#include "DDSegmentation/BitFieldCoder.h"

#include "k4FWCore/DataHandle.h"
#include "k4FWCore/Transformer.h"
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

// DD4hep field and material for GenFit
#include "GenfitField.h"
#include "GenfitMaterialInterface.h"

// Forward declarations
namespace dd4hep {
    class OverlayedField;
    namespace rec {
        class Surface;
    }
}

// Define the 5-parameter vector and matrix types
namespace Eigen {
    typedef Matrix<double, 5, 1> Vector5d;
    typedef Matrix<double, 5, 5> Matrix5d;
}

// Mass and charge of the particle hypothesis used for material effects.
struct ParticleProperties {
    double mass;          // Mass in GeV
    double charge;        // Charge in e
    std::string name;     // Particle name

    ParticleProperties() : mass(0), charge(0), name("unknown") {}

    ParticleProperties(double m, double q, const std::string& n)
        : mass(m), charge(q), name(n) {}
};

// Result of the RK4 inner propagation.
struct InnerTrajectory {
    bool success;
    Eigen::Vector3d finalPosition;      // cm
    Eigen::Vector3d finalMomentum;      // GeV/c
    double arcLength;                   // cm, total path length
    int numSteps;                       // total steps taken
    std::string message;                // status/error message
};

/** @class CombinatorialTrackFinder
 *
 *  Track finding on the planar hits of one tracking system (inner tracker, silicon
 *  wrapper, muon system, ...), with no input from other detectors.
 *
 *  How it works, per event:
 *   1. Hits are grouped by layer (barrel / endcap layers from the cellID).
 *   2. Track candidates are found either
 *        - Exhaustive: all k-hit combinations over the layers, each tested with a
 *          circle fit in xy and a line fit of z versus arc length, or
 *        - SeedAndExtend: 3-5 hit seeds on neighbouring layers, extended layer by
 *          layer with a helix (BestFirst: best seeds first, kink chi2 with multiple
 *          scattering, hits claimed as tracks are accepted).
 *      Road, outlier and quality cuts reject wrong combinations.
 *   3. The best candidate is stored as a track; its hits are removed and the search
 *      repeats until no candidate is left.
 *   4. Each track stores the fitted helix as EDM4hep track states: AtFirstHit and
 *      AtLastHit, optionally AtIP (IPStateMode). Optional: a GenFit fit (UseGenFit,
 *      state AtOther).
 *
 *  @author Mahmoud Althakeel
 *  @email  mahmoud.althakeel@cern.ch
 */
// Callback that copies the reco->sim link of an input hit to its output copy
using LinkPropagator = std::function<void(const edm4hep::MutableTrackerHitPlane&, const edm4hep::TrackerHitPlane&)>;

class CombinatorialTrackFinder final
    : public k4FWCore::MultiTransformer<
        std::tuple<edm4hep::TrackCollection,
                   edm4hep::TrackerHitPlaneCollection,
                   edm4hep::TrackerHitSimTrackerHitLinkCollection>(
                       const edm4hep::TrackerHitPlaneCollection&,
                       const edm4hep::EventHeaderCollection&,
                       const edm4hep::TrackerHitSimTrackerHitLinkCollection&)
      > {
public:
    CombinatorialTrackFinder(const std::string& name, ISvcLocator* pSvcLocator);
    virtual ~CombinatorialTrackFinder() = default;

    StatusCode initialize() override;
    StatusCode finalize() override;

    std::tuple<edm4hep::TrackCollection,
               edm4hep::TrackerHitPlaneCollection,
               edm4hep::TrackerHitSimTrackerHitLinkCollection> operator()(
        const edm4hep::TrackerHitPlaneCollection& hits,
        const edm4hep::EventHeaderCollection& headers,
        const edm4hep::TrackerHitSimTrackerHitLinkCollection& recoSimLinks) const override;
    
    // Pattern recognition and track building for one event
    void findTracks(
                            const edm4hep::TrackerHitPlaneCollection* hits,
                            edm4hep::TrackCollection& outputTracks,
                            edm4hep::TrackerHitPlaneCollection& outputHits,
                            const edm4hep::TrackerHitSimTrackerHitLinkCollection& recoSimLinks,
                            edm4hep::TrackerHitSimTrackerHitLinkCollection& outputLinks,
                            bool& evtUsedPtInner,
                            bool& evtUsedPtInner2,
                            bool& evtUsedPtFallback) const;

    // Detector surface of a hit
    const dd4hep::rec::Surface* findSurface(const edm4hep::TrackerHitPlane& hit) const;
    
private:
    // Detector surface of a cellID
    const dd4hep::rec::Surface* findSurfaceByID(uint64_t cellID) const;
    // Extract type ID from cell ID
    int getTypeID(uint64_t cellID) const;
    // Extract layer ID from cell ID
    int getLayerID(uint64_t cellID) const;
    // Layer key: barrel -> layer, +endcap -> 1000 + layer, -endcap -> -1000 - layer
    int getCompositeID(uint64_t cellID) const;
    // pT from omega, with the field at the state's reference point
    double getPT(const edm4hep::TrackState& state) const;

    // Create a track state with the given parameters
    edm4hep::TrackState createTrackState(
                    double d0, double phi, double omega, double z0, double tanLambda,
                    int location = edm4hep::TrackState::AtOther) const;
        
    // Group surfaces by detector layer
    std::map<int, std::vector<const dd4hep::rec::Surface*>> getSurfacesByLayer() const;

    // Circle fit to N >= 3 hits (weighted least squares + Gauss-Newton).
    // fitCov3x3: covariance of (x0, y0, R) in cm^2; residualsCm: per-hit residual (cm).
    // Returns false if the fit is degenerate.
    bool fitCircleNHits(const std::vector<edm4hep::TrackerHitPlane>& hits,
                        double& x0, double& y0, double& radius, double& chi2ndf,
                        Eigen::Matrix3d& fitCov3x3,
                        std::vector<double>& residualsCm) const;

    // GenFit fit of the track hits, seeded by seedState; adds the fitted state (AtOther)
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

    // DD4hep field and material for GenFit
    std::unique_ptr<GenfitInterface::GenfitField> m_genFitField;
    GenfitInterface::GenfitMaterialInterface* m_genFitMaterial{nullptr}; // singleton, not owned
    
    // ==================== TRUTH (diagnostic printouts only) ====================
    mutable k4FWCore::DataHandle<edm4hep::MCParticleCollection> m_mcParticles{
        "MCParticles", Gaudi::DataHandle::Reader, this};

    const edm4hep::MCParticle* findTruthParticle(
        const edm4hep::Track& recoTrack,
        const edm4hep::MCParticleCollection& mcParticles) const;

    void validateTrackWithTruth(
        const edm4hep::Track& recoTrack,
        const edm4hep::MCParticleCollection& mcParticles,
        int trackNumber) const;
    // ==================== PROPAGATION THROUGH THE SOLENOID ====================
    // RK4 propagation of a track measured outside the solenoid, inward
    InnerTrajectory propagateToInner(
        const Eigen::Vector3d& outerPos,    // Starting position in cm
        const Eigen::Vector3d& outerMom,    // Starting momentum in GeV/c
        double targetRadius = 100.0,        // Stop at this R in cm
        double chargeSign   = 1.0           // +1 or -1 particle charge sign
    ) const;

    bool isInsideSolenoid(const Eigen::Vector3d& pos) const;
    bool allHitsOutsideSolenoid(const std::vector<edm4hep::TrackerHitPlane>& trackHits) const;

    // One RK4 step
    void rk4PropagationStep(
        Eigen::Vector3d& pos,      // Position (modified in place)
        Eigen::Vector3d& mom,      // Momentum (modified in place)
        double stepSize,           // Step size in cm
        double chargeSign = 1.0    // +1 or -1 particle charge sign
    ) const;

    // IPStateMode = SolenoidCrossing: the outer circle (centre, radius in cm) is matched
    // at the solenoid boundary to an inner helix with the inner field, which is stored as
    // AtIP. pT and tanLambda are kept across the boundary; the innermost hit (cm) selects
    // the crossing point and anchors z. Returns true on success.
    bool analyticalInnerPropagation(
        double cx_o_cm, double cy_o_cm, double R_outer_cm,
        double omega_outer_mm,
        double pT_GeV,
        double tanLambda,
        double x_ref_cm, double y_ref_cm, double z_ref_cm,
        bool   isEndcap,             // true = endcap hit → use z-cap crossing
        edm4hep::MutableTrack& finalTrack) const;

    // Radius (mm) where Bz changes sign, searched around SolenoidRadius
    double findSolenoidBoundary() const;

    // |z| (mm) on the axis where Bz changes sign, searched around SolenoidHalfZ
    double findSolenoidBoundaryZ() const;

    // Table of the track states versus the matched MC particle (debug printout)
    void compareStatesWithTruth(
        const edm4hep::Track& track,
        const edm4hep::MCParticleCollection& mcParticles,
        int trackNumber) const;

    // Solenoid boundary
    static constexpr double SOLENOID_RADIUS_CM = 230.0;  // fallback default
    static constexpr double SOLENOID_Z_HALF_CM = 218.0;
    double m_solenoidBoundaryMm{2300.0};  // set by findSolenoidBoundary()
    double m_solenoidHalfZmm{2180.0};     // set by findSolenoidBoundaryZ()

    // Properties
    Gaudi::Property<std::string> m_detectorName{this, "DetectorName", "Tracker", "Name of detector to process"};
    Gaudi::Property<std::vector<std::string>> m_detectorNames{this, "DetectorNames", {},
        "Optional list of detector names whose surface maps are merged. Use when one readout "
        "spans several DetElements (e.g. GlobalTrackerReadoutID covers VertexBarrel, VertexEndcap, "
        "SiWrB, SiWrD). If empty, DetectorName is used on its own."};
    Gaudi::Property<std::string> m_particleType{this, "ParticleType", "pion", "Particle type for material effects"};
    Gaudi::Property<std::string> m_encodingStringParameter{this, "EncodingStringParameterName", "GlobalTrackerReadoutID", "Name of DD4hep parameter with the encoding string"};
    // CellID field names: the region field is 0 / +1 / -1 for barrel / +endcap / -endcap
    Gaudi::Property<std::string> m_regionFieldName{this, "RegionFieldName", "type",
        "CellID field separating barrel (0) from +endcap (1) and -endcap (-1). "
        "\"type\" for MuonSystemReadoutID, \"side\" for GlobalTrackerReadoutID."};
    Gaudi::Property<std::string> m_layerFieldName{this, "LayerFieldName", "layer",
        "CellID field holding the layer index."};
    // ── Strategy ────────────────────────────────────────────────────────────────
    // Exhaustive   : all k-hit combinations over the layers (sparse events).
    // SeedAndExtend: seeds on NumSeedLayers layers, then helix extension (dense events).
    Gaudi::Property<std::string> m_strategy{this, "Strategy", "Exhaustive",
        "Exhaustive | SeedAndExtend"};
    Gaudi::Property<int> m_numSeedLayers{this, "NumSeedLayers", 3,
        "SeedAndExtend: hits (= layers) per seed, >= 3"};
    Gaudi::Property<std::string> m_seedLayerPolicy{this, "SeedLayerPolicy", "consecutive",
        "SeedAndExtend: consecutive (any run of NumSeedLayers consecutive layers; works for "
        "displaced tracks) | innermost | outermost | explicit (use SeedLayers)"};
    Gaudi::Property<std::vector<std::vector<int>>> m_seedLayers{this, "SeedLayers", {},
        "SeedAndExtend: explicit seed layer sets as lists of compositeIDs, inside-out"};
    Gaudi::Property<std::vector<int>> m_seedPasses{this, "SeedPasses", {},
        "SeedAndExtend: iterative seeding, hits per seed for each pass, e.g. [5,4,3]; later "
        "passes only use hits not claimed by earlier ones. Empty = one pass of NumSeedLayers"};
    Gaudi::Property<int> m_seedClaimMinHits{this, "SeedClaimMinHits", 4,
        "SeedAndExtend: a candidate with at least this many hits claims them from later passes"};
    Gaudi::Property<int> m_seedLayerSkip{this, "SeedLayerSkip", 0,
        "SeedAndExtend: layers that may be skipped inside one seed window (inefficiency)"};
    Gaudi::Property<double> m_seedMaxDeltaPhi{this, "SeedMaxDeltaPhi", 0.3,
        "SeedAndExtend: max azimuth difference (rad) between consecutive seed hits"};
    Gaudi::Property<double> m_seedMaxHitDist{this, "SeedMaxHitDistance", 100.0,
        "SeedAndExtend: max 3D distance (cm) between consecutive seed hits"};
    Gaudi::Property<double> m_seedMinRadius{this, "SeedMinRadius", 0.0,
        "SeedAndExtend: min seed circle radius (cm), i.e. a pT floor; 0 = off"};
    Gaudi::Property<double> m_seedMaxZResid{this, "SeedMaxZResidual", 0.5,
        "SeedAndExtend: max |z residual| (cm) of seed hits about their z-vs-s line"};
    Gaudi::Property<double> m_seedMaxChi2NDF{this, "SeedMaxChi2NDF", 10.0,
        "SeedAndExtend: seed circle chi2/ndf cut (only bites for NumSeedLayers > 3)"};
    Gaudi::Property<std::string> m_extDirection{this, "ExtrapolationDirection", "Outward",
        "SeedAndExtend: Outward | Inward | Bidirectional"};
    Gaudi::Property<double> m_extWindow{this, "ExtensionWindow", 1.0,
        "SeedAndExtend: hit search radius (cm) around the predicted crossing ..."};
    Gaudi::Property<double> m_extWindowSlope{this, "ExtensionWindowSlope", 0.02,
        "... growing by this many cm per cm extrapolated"};
    Gaudi::Property<double> m_extMaxChi2NDF{this, "ExtensionMaxChi2NDF", 10.0,
        "SeedAndExtend: refit chi2/ndf a newly attached hit must keep"};
    Gaudi::Property<int> m_maxMissedLayers{this, "MaxMissedLayers", 1,
        "SeedAndExtend: stop extending after this many consecutive layers without a hit"};
    Gaudi::Property<int> m_maxCandPerLayer{this, "MaxCandidatesPerLayer", 3,
        "SeedAndExtend: hits tried per crossed layer, nearest first"};
    Gaudi::Property<std::string> m_seMode{this, "SeedAndExtendMode", "Classic",
        "SeedAndExtend: Classic (all seeds extended, then ranked) | BestFirst (windows inside-out, "
        "best seed first, immediate hit claiming, multiple-scattering-aware kink chi2, confirmation "
        "pass, merging of track pieces)"};
    Gaudi::Property<double> m_kinkMaxChi2Triplet{this, "KinkMaxChi2Triplet", 16.0,
        "BestFirst: a hit is attached only if the kink chi2 (2 dof) of its triplet is below this"};
    Gaudi::Property<double> m_seedKinkMaxChi2NDF{this, "SeedKinkMaxChi2NDF", 9.0,
        "BestFirst: seed acceptance on kink chi2/ndf"};
    Gaudi::Property<double> m_trackKinkMaxChi2NDF{this, "TrackKinkMaxChi2NDF", 5.0,
        "BestFirst: a finished track is kept, and claims its hits, below this kink chi2/ndf"};
    Gaudi::Property<int> m_confirmIterations{this, "ConfirmationIterations", 2,
        "BestFirst: confirmation passes (fill crossed-but-empty layers, swap hits, drop outliers)"};
    Gaudi::Property<bool> m_mergePieces{this, "MergeTrackPieces", true,
        "BestFirst: merge sequential pieces of one trajectory; overlapping arcs are flagged as loopers (track type 1)"};
    Gaudi::Property<int> m_shortTrackMaxHits{this, "ShortTrackMaxHits", 4,
        "BestFirst: tracks with at most this many hits must pass the hole check (ShortTrackMaxHoles)"};
    Gaudi::Property<int> m_shortTrackMaxHoles{this, "ShortTrackMaxHoles", 1,
        "BestFirst: max layers crossed without a hit (inside the track's span and outward to the "
        "tracker exit; inward not counted) for a short track"};
    Gaudi::Property<double> m_shortTrackMinPt{this, "ShortTrackMinPt", 0.0,
        "BestFirst: tracks with at most ShortTrackMaxHits hits need fitted pT >= this (GeV); 0 = off"};
    Gaudi::Property<double> m_shadowHitWindow{this, "ShadowHitWindow", 0.1,
        "BestFirst: unclaimed hits within this distance (cm) of an accepted track's trajectory are "
        "claimed (overlapping-sensor double hits), so they cannot seed a duplicate; 0 = off"};
    Gaudi::Property<double> m_layerMatHalfDepth{this, "LayerMaterialHalfDepth", 150.0,
        "BestFirst: per-layer x/X0 is integrated from halfway to the inner neighbour layer to halfway "
        "to the outer one, each side capped at this depth (cm)"};
    Gaudi::Property<long> m_maxSeedsPerEvent{this, "MaxSeedsPerEvent", 200000,
        "SeedAndExtend: safety cap on seeds evaluated per event"};

    Gaudi::Property<std::string> m_layerGeoDump{this, "LayerGeometryDump", "",
        "If set, write the per-layer geometry table (radius/z extents, barrel/endcap) to "
        "this CSV path at initialize()."};
    Gaudi::Property<int> m_systemInCompositeID{this, "UseSystemInCompositeID", -1,
        "Fold the cellID 'system' field into the per-layer key. -1 = auto (on when "
        "DetectorNames lists more than one detector), 0 = off, 1 = on. Required when "
        "one readout spans several subdetectors, otherwise e.g. VertexBarrel layer 0 "
        "and SiWrB layer 0 collapse onto the same key."};
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
    Gaudi::Property<bool> m_doInnerPropagation{this, "DoInnerPropagation", false,
        "Propagate the track inward through the solenoid boundary with RK4 (tracks measured "
        "outside the solenoid). The result is printed in the log; no state is stored."};
    Gaudi::Property<double> m_innerPropTargetRadius{this, "InnerPropTargetRadius", 100.0,
        "Target radius in cm for inner RK4 propagation"};
    Gaudi::Property<std::string> m_ipStateMode{this, "IPStateMode", "None",
        "Extra track state AtIP: None (not stored), Helix (the fitted helix at its closest "
        "approach to the origin; for tracks in a uniform field, e.g. an inner tracker) or "
        "SolenoidCrossing (helix matched at the solenoid boundary; only for tracks measured "
        "outside the solenoid, e.g. the muon system)"};
    Gaudi::Property<bool> m_storeLastHitState{this, "StoreLastHitState", true,
        "Also store the fitted helix at the outermost hit (AtLastHit)"};
    Gaudi::Property<double> m_maxComboPT{this, "MaxComboPT", 200.0,
        "Maximum allowed combo pT in GeV/c. Combinations whose circle-fit radius implies pT above this threshold are excluded during the combinatorial search."};
    Gaudi::Property<double> m_minComboPT{this, "MinComboPT", -1.0,
        "Minimum allowed combo pT in GeV/c. Combinations whose circle-fit radius implies pT below this threshold are "
        "rejected during the combinatorial search (no fallback). A value <= 0 disables the cut (default: disabled)."};


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
        "Enable for high-occupancy events where delta-ray or multi-track hits inflate the "
        "number of combinations."};

    Gaudi::Property<bool> m_requirePairCoinc{this, "RequirePairCoincidence", false,
        "If true, exploit paired-layer geometry (doublets at layers (0,1), (2,3), (4,5)) "
        "as a hit-quality filter. Real muons traverse both layers of a doublet within a "
        "tight (Δφ, Δη) window; delta-rays/secondaries usually fire only one. "
        "Two-pass: pass 1 enumerates combos using only hits with a coincident partner in "
        "the paired layer; pass 2 (if pass 1 finds nothing) uses all hits. Either way, "
        "candidate combos are required to have ≥2 complete doublets — exactly 3 is best, "
        "2 is acceptable, <2 rejected."};

    Gaudi::Property<bool> m_pairCoincSoftScoring{this, "PairCoincidenceSoftScoring", false,
        "Only honoured when RequirePairCoincidence=true. Switches the doublet quality "
        "logic from HARD (reject combos with <2 complete doublets) to SOFT: keep all "
        "combos in the running but bias bestCombo selection to prefer higher doublet "
        "count (3>2>1>0) as the primary criterion ahead of hit count and chi². "
        "Soft mode disables the seed-eligible hit filter and the two-pass fallback — "
        "it runs a single pass over all hits, just with doublet-aware ranking."};

    Gaudi::Property<double> m_pairCoincMaxDistMM{this, "PairCoincidenceMaxDistMM", 250.0,
        "Max 3D distance (mm) between two hits treated as a paired-layer doublet. "
        "250 mm suits the muon system (paired layers 12-24 cm apart); scale it down "
        "for inner-tracker readouts where layers are far closer."};

    // ── Neighbour-track quality cuts ──────────────────────────────────────
    Gaudi::Property<double> m_neighbourTrackMaxDist{this, "NeighbourTrackMaxDist", 5.0,
        "If the best-combo candidate has any hit within this distance (cm) of a hit from "
        "an already-accepted track, apply stricter quality requirements defined by "
        "NeighbourTrackMinHits and NeighbourTrackMaxChi2NDF. "
        "Rejects delta-ray shadow tracks that form next to a real track. Set to 0 to disable."};
    Gaudi::Property<int> m_neighbourTrackMinHits{this, "NeighbourTrackMinHits", 4,
        "Minimum hit count required for a combo flagged as a neighbour track "
        "(any hit lies within NeighbourTrackMaxDist of an accepted track hit). "
        "Stronger than MinTrackHits, to reject 3-hit shadow fakes."};
    Gaudi::Property<double> m_neighbourTrackMaxChi2NDF{this, "NeighbourTrackMaxChi2NDF", 3.0,
        "Maximum chi2/NDF accepted for a neighbour-track combo. "
        "Applied together with NeighbourTrackMinHits when a combo is spatially near "
        "an already-accepted track."};
    // Services
    ServiceHandle<IGeoSvc> m_geoSvc{this, "GeoSvc", "GeoSvc", "Detector geometry service"};
    
    // Geometry and field
    dd4hep::Detector* m_detector{nullptr};
    dd4hep::OverlayedField m_field;
    const dd4hep::rec::SurfaceMap* m_surfaceMap{nullptr};
    std::vector<const dd4hep::rec::Surface*> m_surfaces;
    std::map<int, std::vector<const dd4hep::rec::Surface*>> m_surfacesByLayer;
    k4rt::LayerGeometry m_layerGeo;          // per-layer radius/z extents, built from m_surfacesByLayer
    std::unordered_map<int, int> m_cidRegion;  // compositeID -> readout region (0, +1, -1); see k4rt::regionOfCid
    bool m_useSeedExtend{false};             // Strategy == "SeedAndExtend"
    bool m_seBestFirst{false};               // SeedAndExtendMode == "BestFirst"
    std::unique_ptr<k4rt::SeedFinder> m_seedFinder;  // built in initialize() when m_useSeedExtend
    ParticleProperties m_particleProperties;  // hypothesis for material effects
    std::unique_ptr<dd4hep::DDSegmentation::BitFieldCoder> m_bitFieldCoder; // cellID decoding
    bool m_hasRegionField{true};   // false when the readout has no barrel/endcap field at all
    bool m_useSystemInCompositeID{false};  // fold "system" into the layer key (multi-subdetector readouts)
    // merged surface map when several detectors are combined (DetectorNames)
    dd4hep::rec::SurfaceMap m_mergedSurfaceMap;

    int getPDGCode() const;
    
    // known particle hypotheses
    std::map<std::string, ParticleProperties> m_particleMap;

    // ==================== RUN STATISTICS (printed in finalize) ====================
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
    mutable std::atomic<int> m_statComboMinPTRejected{0};   // combos rejected by MinComboPT threshold
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

    // SolenoidCrossing AtIP: failure reasons
    mutable std::atomic<int> m_statAnalBarrelNoIsect{0};    // outer circle doesn't cross solenoid barrel
    mutable std::atomic<int> m_statAnalBarrelBothZfail{0};  // both crossings have |z_sol|>|z_ref|
    mutable std::atomic<int> m_statAnalBarrelZswapped{0};   // swapped successfully to alt crossing
    mutable std::atomic<int> m_statAnalEndcapSmallTanL{0};  // |tanLambda| too small to extrapolate z
    mutable std::atomic<int> m_statAnalEndcapPastCap{0};    // ref hit already past endcap cap
    mutable std::atomic<int> m_statAnalEndcapArcLong{0};    // arc would exceed one loop

    // Pair-coincidence (doublet) counters
    mutable std::atomic<int> m_statHitsNonCoinc{0};       // hits with no partner in paired layer
    mutable std::atomic<int> m_statCombosRejByCoinc{0};   // combos rejected by <2 complete doublets
    mutable std::atomic<int> m_statCoincPass2Tracks{0};   // tracks recovered by pass-2 fallback

    // GenFit detailed counters
    mutable std::atomic<int>       m_statGenFitFailed{0};
    mutable std::atomic<int>       m_statGenFitIllCond{0};
    mutable std::atomic<int>       m_statGenFitException{0};
    mutable std::atomic<int>       m_statGenFitChi2Count{0};
    mutable std::atomic<int>       m_statGenFitChi2BadCount{0};
    mutable std::atomic<long long> m_statGenFitChi2Sum{0};

    // not filled by the current code; kept for the finalize printout
    mutable std::atomic<int> m_statTripletBadGeom{0};
    mutable std::atomic<int> m_statAngleGuardRejected{0};
    mutable std::atomic<int> m_statHighestPtInner{0};
    mutable std::atomic<int> m_statHighestPtInner2{0};
    mutable std::atomic<int> m_statHighestPtFallback{0};
    mutable std::atomic<int> m_statFourHitTracks{0};
    mutable std::atomic<int> m_statThreeHitTracks{0};
};

DECLARE_COMPONENT(CombinatorialTrackFinder)
#endif // COMBINATORIAL_TRACK_FINDER_H
/**
 * CombinatorialTrackFinder: track finding on the planar hits of one tracking system.
 * See CombinatorialTrackFinder.h for an overview of the algorithm.
 *
 * @author Mahmoud Althakeel
 * @email  mahmoud.althakeel@cern.ch
 */

#include <iostream>
#include <iomanip>
#include <chrono>
#include "CombinatorialTrackFinder.h"
#include "TrackingPrimitives.h"
#include "DDRec/SurfaceManager.h"
#include "DDRec/SurfaceHelper.h"
#include "DD4hep/DD4hepUnits.h"
#include "podio/Frame.h"
#include "k4FWCore/DataHandle.h"

// Particle masses for the material-effects hypothesis
constexpr double electron_mass = 0.000511;  // GeV
constexpr double muon_mass = 0.10566;       // GeV
constexpr double pion_mass = 0.13957;       // GeV
constexpr double kaon_mass = 0.49368;       // GeV
constexpr double proton_mass = 0.93827;     // GeV
constexpr double c_light = 299792458.0;     // m/s
constexpr double qe = 1.602176634e-19;      // coulombs

CombinatorialTrackFinder::CombinatorialTrackFinder(const std::string& name, ISvcLocator* pSvcLocator)
    : MultiTransformer(name, pSvcLocator,
        {
            KeyValues{"InputHitCollection",       {"TrackerHits"}},
            KeyValues{"HeaderCollectionName",     {"EventHeader"}},
            KeyValues{"InputRecoSimLinkCollection", {"TrackerHitSimTrackerHitLink"}}
        },
        {
            KeyValues{"OutputTrackCollection",    {"DisplacedTracks"}},
            KeyValues{"OutputHitCollection",      {"TrackHits"}},
            KeyValues{"OutputRecoSimLinkCollection", {"DisplacedTrackHitSimLinks"}}
        }) {

    m_particleMap["electron"] = ParticleProperties(electron_mass, -1.0, "electron");
    m_particleMap["positron"] = ParticleProperties(electron_mass, 1.0, "positron");
    m_particleMap["muon"] = ParticleProperties(muon_mass, -1.0, "muon");
    m_particleMap["antimuon"] = ParticleProperties(muon_mass, 1.0, "antimuon");
    m_particleMap["pion"] = ParticleProperties(pion_mass, 1.0, "pion");
    m_particleMap["kaon"] = ParticleProperties(kaon_mass, 1.0, "kaon");
    m_particleMap["proton"] = ParticleProperties(proton_mass, 1.0, "proton");
    m_particleMap["antiproton"] = ParticleProperties(proton_mass, -1.0, "antiproton");
}

StatusCode CombinatorialTrackFinder::initialize() {
    StatusCode sc = Algorithm::initialize();
    if (!sc.isSuccess()) return sc;

    if (!m_geoSvc.retrieve()) {
        error() << "Failed to retrieve GeoSvc" << endmsg;
        return StatusCode::FAILURE;
    }

    m_detector = m_geoSvc->getDetector();
    if (!m_detector) {
        error() << "Failed to access detector" << endmsg;
        return StatusCode::FAILURE;
    }

    // cellID decoder
    std::string cellIDEncodingString = m_geoSvc->constantAsString(m_encodingStringParameter);
    m_bitFieldCoder = std::make_unique<dd4hep::DDSegmentation::BitFieldCoder>(cellIDEncodingString);

    // Check that the configured cellID field names exist in this readout.
    {
        const auto& desc = m_bitFieldCoder->fieldDescription();
        auto hasField = [&](const std::string& n) {
            for (const auto& f : m_bitFieldCoder->fields())
                if (f.name() == n) return true;
            return false;
        };
        if (!hasField(m_layerFieldName.value())) {
            error() << "Readout '" << m_encodingStringParameter.value() << "' (" << desc
                    << ") has no field '" << m_layerFieldName.value()
                    << "'. Set LayerFieldName to match this readout." << endmsg;
            return StatusCode::FAILURE;
        }
        m_hasRegionField = hasField(m_regionFieldName.value());
        if (!m_hasRegionField) {
            warning() << "Readout '" << m_encodingStringParameter.value() << "' (" << desc
                      << ") has no field '" << m_regionFieldName.value()
                      << "'; treating every hit as barrel. Set RegionFieldName "
                      << "(\"type\" for MuonSystemReadoutID, \"side\" for GlobalTrackerReadoutID)."
                      << endmsg;
        }
        // When one readout spans several subdetectors, "system" goes into the layer key,
        // so that equally numbered layers of different subdetectors stay apart.
        const bool multiDet = m_detectorNames.value().size() > 1;
        if (m_systemInCompositeID.value() < 0)
            m_useSystemInCompositeID = multiDet && hasField("system");
        else
            m_useSystemInCompositeID = (m_systemInCompositeID.value() != 0);
        if (m_useSystemInCompositeID && !hasField("system")) {
            error() << "UseSystemInCompositeID requested but readout '"
                    << m_encodingStringParameter.value() << "' has no 'system' field." << endmsg;
            return StatusCode::FAILURE;
        }
        if (multiDet && !m_useSystemInCompositeID) {
            warning() << m_detectorNames.value().size() << " detectors share this readout but "
                      << "'system' is NOT folded into the layer key -- layers with the same "
                      << "index in different subdetectors will be merged." << endmsg;
        }

        info() << "CellID decoding: layer field '" << m_layerFieldName.value()
               << "', region field '" << m_regionFieldName.value()
               << (m_hasRegionField ? "'" : "' (absent, barrel-only)")
               << ", system in layer key: " << (m_useSystemInCompositeID ? "yes" : "no") << endmsg;
    }

    // Surface map: merged over DetectorNames when given, else the one of DetectorName.
    auto surfaceMan = m_detector->extension<dd4hep::rec::SurfaceManager>();
    std::vector<std::string> detNames = m_detectorNames.value();
    if (detNames.empty()) detNames.push_back(m_detectorName.value());

    if (detNames.size() == 1) {
        m_surfaceMap = surfaceMan->map(detNames.front());
        if (!m_surfaceMap) {
            error() << "Could not find surface map for detector: " << detNames.front() << endmsg;
            return StatusCode::FAILURE;
        }
    } else {
        m_mergedSurfaceMap.clear();
        for (const auto& name : detNames) {
            const auto* m = surfaceMan->map(name);
            if (!m) {
                error() << "Could not find surface map for detector: " << name << endmsg;
                return StatusCode::FAILURE;
            }
            m_mergedSurfaceMap.insert(m->begin(), m->end());
            info() << "  merged " << m->size() << " surfaces from '" << name << "'" << endmsg;
        }
        m_surfaceMap = &m_mergedSurfaceMap;
        info() << "Merged surface map over " << detNames.size() << " detectors: "
               << m_surfaceMap->size() << " surfaces" << endmsg;
    }

    // Magnetic field (printed at a few points for reference)
    m_field = m_detector->field();

    dd4hep::Position center(0, 0, 0);
    dd4hep::Direction bfield = m_field.magneticField(center);
    info() << "Magnetic field at detector center: ("
        << bfield.x() / dd4hep::tesla << ", "
        << bfield.y() / dd4hep::tesla << ", "
        << bfield.z() / dd4hep::tesla << ") Tesla" << endmsg;

    const double solenoidRadius = 230.0; // cm
    const std::vector<dd4hep::Position> testPoints = { // in cm
        dd4hep::Position(0, 0, 0),
        dd4hep::Position(100, 100, 100),
        dd4hep::Position(solenoidRadius, 0, 0),
        dd4hep::Position(0, solenoidRadius, 0),
        dd4hep::Position(solenoidRadius/sqrt(2), solenoidRadius/sqrt(2), 0),
        dd4hep::Position(249, 0, 0),
        dd4hep::Position(0, 251, 0),
        dd4hep::Position(0, 0, 249),
        dd4hep::Position(0, 0, 218),
        dd4hep::Position(0, 0, -217),
        dd4hep::Position(0, 500, 0),
        dd4hep::Position(460, 0, 0),
        dd4hep::Position(200, 200, 460),
        dd4hep::Position(460, 0, 460)
    };

    info() << "Magnetic field at different positions:" << endmsg;
    for (const auto& point : testPoints) {
        dd4hep::Direction bfieldVec = m_field.magneticField(point);
        info() << "  Position (" << point.x() << ", " << point.y() << ", " << point.z()
            << ") cm: B = (" << bfieldVec.x() / dd4hep::tesla << ", "
            << bfieldVec.y() / dd4hep::tesla << ", "
            << bfieldVec.z() / dd4hep::tesla << ") Tesla" << endmsg;
    }

    dd4hep::DetElement world = m_detector->world();
    dd4hep::rec::SurfaceHelper surfHelper(world);

    const dd4hep::rec::SurfaceList& surfList = surfHelper.surfaceList();

    // Surfaces of the selected detector(s), grouped by layer
    m_surfaces.clear();
    for (const auto& [cellID, surface] : *m_surfaceMap) {
        const dd4hep::rec::Surface* concreteSurface =
            dynamic_cast<const dd4hep::rec::Surface*>(surface);
        if (concreteSurface) {
            m_surfaces.push_back(concreteSurface);
        }
    }

    m_surfacesByLayer = getSurfacesByLayer();

    // Layer key -> region (barrel / endcap), taken from the readout's region field.
    // Two different layers must not share a key.
    m_cidRegion.clear();
    {
        int nMismatch = 0;
        for (const auto* s : m_surfaces) {
            dd4hep::DetElement det = s->detElement();
            if (!det.isValid()) continue;
            const uint64_t vid = det.volumeID();
            const int cid = getCompositeID(vid);
            const int reg = getTypeID(vid);
            auto [it, inserted] = m_cidRegion.emplace(cid, reg);
            if (!inserted && it->second != reg) {
                error() << "compositeID " << cid << " is produced by surfaces in different regions ("
                        << it->second << " and " << reg << "): the 100*system+layer encoding collides "
                        << "for this detector" << endmsg;
                return StatusCode::FAILURE;
            }
            if (inserted && reg != k4rt::regionOfCid(cid, nullptr)) ++nMismatch;
        }
        if (nMismatch > 0)
            info() << nMismatch << " compositeID(s) whose region differs from the cid-range rule "
                   << "(barrel system ID >= 10); using the readout region field" << endmsg;
    }

    // Per-layer geometry (radius / z extents, barrel or endcap from the surface normals).
    m_layerGeo.build(m_surfacesByLayer);
    info() << "Layer geometry (inside-out by closest approach to the IP):" << endmsg;
    for (const auto& d : m_layerGeo.layers()) {
        info() << "  cid " << std::setw(6) << d.compositeID << "  " << std::setw(7) << k4rt::kindName(d.kind)
               << "  r=[" << std::fixed << std::setprecision(2) << d.rMin_cm << "," << d.rMax_cm << "]"
               << "  z=[" << d.zMin_cm << "," << d.zMax_cm << "] cm"
               << "  dMin=" << d.dMin_cm << "  nSurf=" << d.nSurfaces
               << std::defaultfloat << std::setprecision(6) << endmsg;
    }
    for (const auto& p : m_layerGeo.validate()) warning() << "LayerGeometry: " << p << endmsg;

    if (m_ipStateMode.value() != "None" && m_ipStateMode.value() != "Helix" &&
        m_ipStateMode.value() != "SolenoidCrossing") {
        error() << "Unknown IPStateMode '" << m_ipStateMode.value()
                << "' (expected None, Helix or SolenoidCrossing)" << endmsg;
        return StatusCode::FAILURE;
    }
    // ── Strategy ────────────────────────────────────────────────────────────────
    if (m_strategy.value() == "SeedAndExtend") {
        m_useSeedExtend = true;
    } else if (m_strategy.value() != "Exhaustive") {
        error() << "Unknown Strategy '" << m_strategy.value()
                << "' (expected Exhaustive or SeedAndExtend)" << endmsg;
        return StatusCode::FAILURE;
    }
    if (m_useSeedExtend) {
        k4rt::SeedConfig sc;
        sc.numSeedLayers   = m_numSeedLayers.value();
        sc.passes          = m_seedPasses.value();
        sc.claimMinHits    = m_seedClaimMinHits.value();
        sc.policy          = m_seedLayerPolicy.value();
        sc.explicitSets    = m_seedLayers.value();
        sc.layerSkip       = m_seedLayerSkip.value();
        sc.maxDeltaPhi     = m_seedMaxDeltaPhi.value();
        sc.maxHitDist_cm   = m_seedMaxHitDist.value();
        sc.minRadius_cm    = m_seedMinRadius.value();
        sc.maxZResid_cm    = m_seedMaxZResid.value();
        sc.seedMaxChi2NDF  = m_seedMaxChi2NDF.value();
        sc.direction       = m_extDirection.value();
        sc.window_cm       = m_extWindow.value();
        sc.windowSlope     = m_extWindowSlope.value();
        sc.extMaxChi2NDF   = m_extMaxChi2NDF.value();
        sc.maxMissedLayers = m_maxMissedLayers.value();
        sc.maxCandPerLayer = m_maxCandPerLayer.value();
        sc.minTrackHits    = m_minTrackHits.value();
        sc.maxSeeds        = m_maxSeedsPerEvent.value();
        sc.mode                = m_seMode.value();
        sc.kinkMaxChi2Triplet  = m_kinkMaxChi2Triplet.value();
        sc.seedKinkMaxChi2NDF  = m_seedKinkMaxChi2NDF.value();
        sc.trackKinkMaxChi2NDF = m_trackKinkMaxChi2NDF.value();
        sc.confirmIterations   = m_confirmIterations.value();
        sc.mergePieces         = m_mergePieces.value();
        sc.shortTrackMaxHits   = m_shortTrackMaxHits.value();
        sc.shortTrackMaxHoles  = m_shortTrackMaxHoles.value();
        sc.shortTrackMinPt     = m_shortTrackMinPt.value();
        sc.shadowHitWindow_cm  = m_shadowHitWindow.value();
        sc.bField_T            = m_field.magneticField(dd4hep::Position(0, 0, 0)).z() / dd4hep::tesla;
        {   // mass of the particle hypothesis
            auto pit = m_particleMap.find(m_particleType);
            sc.mass_GeV = pit != m_particleMap.end() ? pit->second.mass : 0.13957;
        }
        m_seBestFirst = (sc.mode == "BestFirst");
        if (sc.mode != "Classic" && !m_seBestFirst) {
            error() << "Unknown SeedAndExtendMode '" << sc.mode << "' (expected Classic or BestFirst)" << endmsg;
            return StatusCode::FAILURE;
        }
        if (m_seBestFirst) {
            // material of each layer, for the multiple-scattering term of the kink chi2
            m_layerGeo.measureMaterial(*m_detector, m_surfacesByLayer, 24, m_layerMatHalfDepth.value());
            auto& m = info();
            m << "BestFirst: B = " << sc.bField_T << " T, mass " << sc.mass_GeV << " GeV; layer x/X0 [%]:";
            for (const auto& d : m_layerGeo.layers()) m << " " << d.compositeID << ":" << std::fixed << std::setprecision(2) << 100.0 * d.xOverX0;
            m << std::defaultfloat << std::setprecision(6) << endmsg;
        }
        if (sc.direction != "Outward" && sc.direction != "Inward" && sc.direction != "Bidirectional") {
            error() << "Unknown ExtrapolationDirection '" << sc.direction << "'" << endmsg;
            return StatusCode::FAILURE;
        }
        // the seed finder uses the same circle fit and detector surfaces
        auto fit = [this](const std::vector<edm4hep::TrackerHitPlane>& hv, double& x0, double& y0,
                          double& R, double& chi2, Eigen::Matrix3d& cov, std::vector<double>& res) {
            try { return fitCircleNHits(hv, x0, y0, R, chi2, cov, res); } catch (...) { return false; }
        };
        auto plane = [this](uint64_t cellID, Eigen::Vector3d& o, Eigen::Vector3d& n) {
            const dd4hep::rec::Surface* s = findSurfaceByID(cellID);
            if (!s) return false;
            const auto so = s->origin();
            const auto sn = s->normal(so);
            o = Eigen::Vector3d(so.x(), so.y(), so.z());
            n = Eigen::Vector3d(sn.x(), sn.y(), sn.z()).normalized();
            return true;
        };
        m_seedFinder = std::make_unique<k4rt::SeedFinder>(m_layerGeo, sc, fit, plane);
        for (const auto& p : m_seedFinder->prepare()) warning() << "SeedFinder: " << p << endmsg;
        {
            auto& m = info();
            m << "Strategy SeedAndExtend: passes";
            for (int n : m_seedFinder->passes()) m << " [" << n << " hits: " << m_seedFinder->numSets(n) << " layer sets]";
            m << " (policy " << sc.policy << ", skip " << sc.layerSkip << ", claim >= " << sc.claimMinHits
              << " hits), extension " << sc.direction << endmsg;
        }
        if (m_seedFinder->seedLayerSets().empty()) {
            error() << "SeedAndExtend: no seed layer sets — check SeedLayerPolicy/SeedLayers" << endmsg;
            return StatusCode::FAILURE;
        }
    }
    if (!m_layerGeoDump.value().empty()) {
        if (m_layerGeo.dumpCSV(m_layerGeoDump.value()))
            info() << "Layer geometry table written to " << m_layerGeoDump.value() << endmsg;
        else
            warning() << "Could not write layer geometry table to " << m_layerGeoDump.value() << endmsg;
    }

    // Particle hypothesis
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

    // GenFit: DD4hep field and material
    if (m_useGenFit) {
        try {
            m_genFitField = std::make_unique<GenfitInterface::GenfitField>(m_field);
            genfit::FieldManager::getInstance()->init(m_genFitField.get());

            // singleton; getInstance() also registers it with genfit::MaterialEffects
            m_genFitMaterial = GenfitInterface::GenfitMaterialInterface::getInstance(m_detector);
            genfit::MaterialEffects* matEff = genfit::MaterialEffects::getInstance();

            matEff->setMscModel("GEANE");          // multiple-scattering model
            matEff->setEnergyLossBetheBloch(true);
            matEff->setNoiseBetheBloch(true);
            matEff->setDebugLvl(0);

            info() << "GenFit initialized with DD4hep field+material adapters" << endmsg;
        } catch (const genfit::Exception& e) {
            error() << "Error initializing GenFit: " << e.what() << endmsg;
            return StatusCode::FAILURE;
        }
    }

    // Solenoid boundary, from where Bz changes sign
    m_solenoidBoundaryMm = findSolenoidBoundary();
    info() << "Solenoid Bz boundary (radial): R = " << m_solenoidBoundaryMm << " mm  ("
           << m_solenoidBoundaryMm / 10.0 << " cm)" << endmsg;

    m_solenoidHalfZmm = findSolenoidBoundaryZ();
    info() << "Solenoid Bz boundary (endcap): |Z| = " << m_solenoidHalfZmm << " mm  ("
           << m_solenoidHalfZmm / 10.0 << " cm)" << endmsg;

    return StatusCode::SUCCESS;
}

// Operator method
std::tuple<edm4hep::TrackCollection,
               edm4hep::TrackerHitPlaneCollection,
               edm4hep::TrackerHitSimTrackerHitLinkCollection> CombinatorialTrackFinder::operator()(
    const edm4hep::TrackerHitPlaneCollection& hits,
    const edm4hep::EventHeaderCollection& headers,
    const edm4hep::TrackerHitSimTrackerHitLinkCollection& recoSimLinks) const {

    // Find tracks using the direct EDM4hep approach
    edm4hep::TrackCollection trackCollection;
    edm4hep::TrackerHitPlaneCollection outputHits;
    edm4hep::TrackerHitSimTrackerHitLinkCollection outputLinks;

    m_statTotalEvents++;

    // per-event flags filled by findTracks (run statistics only)
    bool evtUsedPtInner    = false;
    bool evtUsedPtInner2   = false;
    bool evtUsedPtFallback = false;

    info() << "\n" << std::string(80, '=') << endmsg;

    unsigned int eventNumber = 0;
    if (!headers.empty()) {
        eventNumber = headers[0].getEventNumber();
        info() << "Processing Event #" << eventNumber << endmsg;
    } else {
        info() << "Processing new event (unknown number)" << endmsg;
    }

    info() << std::string(80, '-') << endmsg;
    info() << "Processing " << hits.size() << " tracker hits" << endmsg;

    if (msgLevel(MSG::DEBUG)) {
        debug() << "Hit details:" << endmsg;
        for (size_t i = 0; i < hits.size(); ++i) {
            const auto& hit = hits[i];
            const auto& pos = hit.getPosition();
            debug() << "  Hit " << i << " at ("
                    << pos[0] << ", " << pos[1] << ", " << pos[2] << ") mm" << endmsg;
        }
    }

    if (hits.size() < 3) {
        warning() << "Not enough hits to create tracks. Need at least 3." << endmsg;
        info() << std::string(80, '=') << "\n" << endmsg;
        return std::make_tuple(std::move(trackCollection), std::move(outputHits), std::move(outputLinks));
    }

    try {
        findTracks(&hits, trackCollection, outputHits, recoSimLinks, outputLinks,
                   evtUsedPtInner, evtUsedPtInner2, evtUsedPtFallback);
    } catch (const std::exception& ex) {
        error() << "Exception during track finding: " << ex.what() << endmsg;
        info() << std::string(80, '=') << "\n" << endmsg;
       return std::make_tuple(std::move(trackCollection), std::move(outputHits), std::move(outputLinks));
    }

    if (evtUsedPtInner)    m_statHighestPtInner++;
    if (evtUsedPtInner2)   m_statHighestPtInner2++;
    if (evtUsedPtFallback) m_statHighestPtFallback++;

    info() << "Found " << trackCollection.size() << " tracks" << endmsg;

    // Per-track printout, with the matched MC particle when the sim links are there
    std::unordered_map<uint32_t, edm4hep::SimTrackerHit> summarySimMap;
    for (const auto& link : recoSimLinks)
        summarySimMap[link.getFrom().id().index] = link.getTo();

    // ── Per-track summary ─────────────────────────────────────────────────────
    for (size_t i = 0; i < trackCollection.size(); ++i) {
        const auto& track = trackCollection[i];
        info() << "─── Track " << (i+1) << "  (" << track.trackerHits_size() << " hits) ───" << endmsg;

        // state at the first hit
        bool hasOuter = false;
        for (int j = 0; j < track.trackStates_size(); ++j) {
            auto ts = track.getTrackStates(j);
            if (ts.location != edm4hep::TrackState::AtFirstHit) continue;
            double B_outer = m_field.magneticField(dd4hep::Position(0,0,0)).z() / dd4hep::tesla;
            // field at the first hit
            if (track.trackerHits_size() > 0) {
                auto h0 = track.getTrackerHits(0);
                B_outer = m_field.magneticField(dd4hep::Position(
                    h0.getPosition()[0]/10.0,
                    h0.getPosition()[1]/10.0,
                    h0.getPosition()[2]/10.0)).z() / dd4hep::tesla;
            }
            double pT_out = (std::abs(ts.omega) > 1e-12)
                            ? 0.3 * std::abs(B_outer) * 1e-3 / std::abs(ts.omega) : 0.0;
            double eta_out = std::asinh(ts.tanLambda);
            info() << "  [AtFirstHit]       pT=" << std::fixed << std::setprecision(2) << pT_out
                   << " GeV  phi=" << std::setprecision(3) << ts.phi
                   << " rad  eta=" << std::setprecision(3) << eta_out
                   << "  d0=" << std::setprecision(2) << ts.D0/10.0 << " cm"
                   << "  z0=" << ts.Z0/10.0 << " cm  (w.r.t. the first hit)" << endmsg;
            hasOuter = true;
            break;
        }

        // state at the IP, if stored (IPStateMode)
        bool hasIP = false;
        for (int j = 0; j < track.trackStates_size(); ++j) {
            auto ts = track.getTrackStates(j);
            if (ts.location != edm4hep::TrackState::AtIP) continue;
            double B_inner = m_field.magneticField(dd4hep::Position(0,0,0)).z() / dd4hep::tesla;
            double pT_ip = (std::abs(ts.omega) > 1e-12)
                           ? 0.3 * std::abs(B_inner) * 1e-3 / std::abs(ts.omega) : 0.0;
            double eta_ip = std::asinh(ts.tanLambda);
            info() << "  [Inner/AtIP]       pT=" << std::fixed << std::setprecision(2) << pT_ip
                   << " GeV  phi=" << std::setprecision(3) << ts.phi
                   << " rad  eta=" << std::setprecision(3) << eta_ip
                   << "  d0=" << std::setprecision(2) << ts.D0/10.0 << " cm"
                   << "  z0=" << ts.Z0/10.0 << " cm" << endmsg;
            hasIP = true;
            break;
        }

        // MC particle with most hits on the track
        std::unordered_map<int, std::pair<edm4hep::MCParticle,int>> mcCount;
        for (int j = 0; j < track.trackerHits_size(); ++j) {
            auto hit = track.getTrackerHits(j);
            auto it  = summarySimMap.find(hit.id().index);
            if (it == summarySimMap.end()) continue;
            auto mc = it->second.getParticle();
            if (!mc.isAvailable()) continue;
            mcCount[mc.id().index].first  = mc;
            mcCount[mc.id().index].second++;
        }
        int bestN = 0;  edm4hep::MCParticle bestMC;  bool hasMC = false;
        for (const auto& [idx, pr] : mcCount)
            if (pr.second > bestN) { bestN = pr.second; bestMC = pr.first; hasMC = true; }

        if (hasMC) {
            const auto& mom = bestMC.getMomentum();
            double mcPT  = std::sqrt(mom.x*mom.x + mom.y*mom.y);
            double mcEta = (mcPT > 0) ? std::asinh(mom.z / mcPT) : 0.0;
            double mcPhi = std::atan2(mom.y, mom.x);
            const auto& vtx = bestMC.getVertex();
            double mcD0 = (-vtx.x * std::sin(mcPhi) + vtx.y * std::cos(mcPhi)) / 10.0; // mm→cm
            double mcZ0 = vtx.z / 10.0;  // mm → cm

            info() << "  [Truth PDG=" << bestMC.getPDG() << " q=" << bestMC.getCharge()
                   << "]        pT=" << std::fixed << std::setprecision(2) << mcPT
                   << " GeV  phi=" << std::setprecision(3) << mcPhi
                   << " rad  eta=" << std::setprecision(3) << mcEta
                   << "  d0=" << std::setprecision(2) << mcD0 << " cm"
                   << "  z0=" << mcZ0 << " cm" << endmsg;

            if (hasIP) {
                for (int j = 0; j < track.trackStates_size(); ++j) {
                    auto ts = track.getTrackStates(j);
                    if (ts.location != edm4hep::TrackState::AtIP) continue;
                    double B_inner = m_field.magneticField(dd4hep::Position(0,0,0)).z()
                                     / dd4hep::tesla;
                    double pT_ip = (std::abs(ts.omega)>1e-12)
                                   ? 0.3*std::abs(B_inner)*1e-3/std::abs(ts.omega) : 0.0;
                    double dphi = ts.phi - mcPhi;
                    while (dphi >  M_PI) dphi -= 2*M_PI;
                    while (dphi < -M_PI) dphi += 2*M_PI;
                    info() << "  [Resolution AtIP]  ΔpT/pT="
                           << std::fixed << std::setprecision(1)
                           << (mcPT>0 ? (pT_ip-mcPT)/mcPT*100.0 : 0.0) << "%"
                           << "  Δphi=" << std::setprecision(4) << dphi << " rad"
                           << "  Δeta=" << std::setprecision(4) << std::asinh(ts.tanLambda)-mcEta
                           << "  Δd0=" << std::setprecision(2) << ts.D0/10.0-mcD0 << " cm"
                           << "  Δz0=" << ts.Z0/10.0-mcZ0 << " cm" << endmsg;
                    break;
                }
            }
        } else {
            info() << "  [Truth] no MC link found" << endmsg;
        }
    }

    info() << std::string(80, '=') << "\n" << endmsg;

    return std::make_tuple(std::move(trackCollection), std::move(outputHits), std::move(outputLinks));
}

// Prints the run statistics.
StatusCode CombinatorialTrackFinder::finalize() {
    int nEvt      = m_statTotalEvents.load();
    int nTrk      = m_statTracksReconstructed.load();
    int nPos      = m_statPositiveCharge.load();
    int nNeg      = m_statNegativeCharge.load();
    int nUndet    = m_statChargeUndetermined.load();
    int nCombos   = m_statTotalTripletCombos.load();
    int nValid    = m_statValidTriplets.load();
    int nDupComp  = m_statDupCompositeRejected.load();
    int nLSRej    = m_statLayerSpanRejected.load();
    int nIsoRej        = m_statIsolationRejected.load();
    int nConsecPhiRej  = m_statConsecPhiRejected.load();
    int nPhiSpreadRej  = m_statPhiSpreadRejected.load();
    int nConsecDistRej = m_statConsecDistRejected.load();
    int nPairDistRej   = m_statPairDistRejected.load();
    int nEdepRej       = m_statEdepRejected.load();
    int nOutlier       = m_statOutlierHitsRemoved.load();
    int nNeighbRej     = m_statNeighbourRejected.load();
    int nPTRej         = m_statComboPTRejected.load();
    int nMinPTRej      = m_statComboMinPTRejected.load();
    int nTanLRej       = m_statTanLambdaRejected.load();
    int nClusters      = m_statClustersBuilt.load();
    int nH3       = m_statNhit3.load();
    int nH4       = m_statNhit4.load();
    int nH5       = m_statNhit5.load();
    int nH6       = m_statNhit6.load();
    int nH7p      = m_statNhit7plus.load();
    int nSL       = m_statStateAtLastHit.load();
    int nSO       = m_statStateAtOther.load();
    int nSV       = m_statStateAtVertex.load();
    int nSIP      = m_statStateAtIP.load();
    int nProp     = m_statInnerPropSuccess.load();
    int nAnalyt   = m_statAnalyticalPropSuccess.load();
    int nGF       = m_statGenFitSuccess.load();
    int nGFFail   = m_statGenFitFailed.load();
    int nGFIll    = m_statGenFitIllCond.load();
    int nGFExc    = m_statGenFitException.load();
    int nGFChi2N  = m_statGenFitChi2Count.load();
    double avgChi2ndf = (nGFChi2N > 0)
                        ? (m_statGenFitChi2Sum.load() / 1000.0) / nGFChi2N
                        : 0.0;
    double trkPerEvt  = (nEvt > 0)    ? double(nTrk)  / double(nEvt)   : 0.0;
    double comboEff   = (nCombos > 0) ? 100.0*double(nValid)/double(nCombos) : 0.0;
    double propFrac   = (nTrk > 0)    ? 100.0*double(nProp)   / double(nTrk) : 0.0;
    double analFrac   = (nTrk > 0)    ? 100.0*double(nAnalyt) / double(nTrk) : 0.0;
    double gfFrac     = (nTrk > 0)    ? 100.0*double(nGF)   / double(nTrk)  : 0.0;
    double avgOutlier = (nValid > 0)  ? double(nOutlier)/double(nValid) : 0.0;
    info() << "\n"
        << "╔══════════════════════════════════════════════════════════╗\n"
        << "║       COMBINATORIAL TRACK FINDER  -  RUN SUMMARY        ║\n"
        << "╠══════════════════════════════════════════════════════════╣\n"
        << "║  Events processed             : " << std::setw(7) << nEvt    << "                   ║\n"
        << "║  Tracks reconstructed         : " << std::setw(7) << nTrk    << "                   ║\n"
        << "║  Avg tracks / event           : " << std::setw(7) << std::fixed << std::setprecision(2) << trkPerEvt << "                   ║\n"
        << "╠══════════════════════════════════════════════════════════╣\n"
        << "║  CURVATURE / CHARGE SIGN                                 ║\n"
        << "║    Positive particles (w > 0) : " << std::setw(7) << nPos    << "                   ║\n"
        << "║    Negative particles (w < 0) : " << std::setw(7) << nNeg    << "                   ║\n"
        << "║    Unreliable sign (near-str) : " << std::setw(7) << nUndet  << "                   ║\n"
        << "╠══════════════════════════════════════════════════════════╣\n"
        << "║  N-HIT COMBINATORIAL SEARCH                              ║\n"
        << "║    C(N,k) combos evaluated    : " << std::setw(7) << nCombos  << "                   ║\n"
        << "║    Tracks found (best combo)  : " << std::setw(7) << nValid   << "  (" << std::setw(5) << std::fixed << std::setprecision(1) << comboEff << "% of combos)  ║\n"
        << "║    Rejected by dupComposite   : " << std::setw(7) << nDupComp << "                   ║\n"
        << "║    Rejected by LayerSpan cut  : " << std::setw(7) << nLSRej   << "                   ║\n"
        << "║    Rejected by consec Δφ cut  : " << std::setw(7) << nConsecPhiRej  << "                   ║\n"
        << "║    Rejected by φ spread cut   : " << std::setw(7) << nPhiSpreadRej  << "                   ║\n"
        << "║    Rejected by consec dist    : " << std::setw(7) << nConsecDistRej << "                   ║\n"
        << "║    Rejected by pair dist      : " << std::setw(7) << nPairDistRej   << "                   ║\n"
        << "║    Rejected by Isolation cut  : " << std::setw(7) << nIsoRej        << "                   ║\n"
        << "║    Hits removed by edep cut   : " << std::setw(7) << nEdepRej       << "                   ║\n"
        << "║    Outlier hits removed       : " << std::setw(7) << nOutlier<< "  (" << std::setw(5) << std::fixed << std::setprecision(2) << avgOutlier << " / track)      ║\n"
        << "║    Rejected by neighbour gate : " << std::setw(7) << nNeighbRej     << "                   ║\n"
        << "║    Rejected by MaxComboPT     : " << std::setw(7) << nPTRej         << "                   ║\n"
        << "║    Rejected by MinComboPT     : " << std::setw(7) << nMinPTRej      << "                   ║\n"
        << "║    Rejected by tanλ dev cut   : " << std::setw(7) << nTanLRej       << "                   ║\n"
        << "║    Hit clusters built (total) : " << std::setw(7) << nClusters      << "                   ║\n"
        << "╠══════════════════════════════════════════════════════════╣\n"
        << "║  HIT MULTIPLICITY PER TRACK                              ║\n"
        << "║    3 hits                     : " << std::setw(7) << nH3     << "                   ║\n"
        << "║    4 hits                     : " << std::setw(7) << nH4     << "                   ║\n"
        << "║    5 hits                     : " << std::setw(7) << nH5     << "                   ║\n"
        << "║    6 hits                     : " << std::setw(7) << nH6     << "                   ║\n"
        << "║    7+ hits                    : " << std::setw(7) << nH7p    << "                   ║\n"
        << "╠══════════════════════════════════════════════════════════╣\n"
        << "║  TRACK STATES STORED                                     ║\n"
        << "║    AtLastHit   (helix)        : " << std::setw(7) << nSL     << "                   ║\n"
        << "║    AtOther     (GenFit fit)   : " << std::setw(7) << nSO     << "                   ║\n"
        << "║    AtVertex    (RK4 inner)    : " << std::setw(7) << nSV     << "                   ║\n"
        << "║    AtIP        (IPStateMode)  : " << std::setw(7) << nSIP    << "                   ║\n"
        << "╠══════════════════════════════════════════════════════════╣\n"
        << "║  DOWNSTREAM PROCESSING                                   ║\n"
        << "║    RK4 inner prop. OK         : " << std::setw(7) << nProp   << "  (" << std::setw(5) << std::fixed << std::setprecision(1) << propFrac << "% of tracks)  ║\n"
        << "║    Analytical inner prop. OK  : " << std::setw(7) << nAnalyt << "  (" << std::setw(5) << std::fixed << std::setprecision(1) << analFrac << "% of tracks)  ║\n"
        << "╠══════════════════════════════════════════════════════════╣\n"
        << "║  SOLENOID-CROSSING AtIP: FAILURE REASONS                 ║\n"
        << "║    Barrel: no solenoid isect  : " << std::setw(7) << m_statAnalBarrelNoIsect.load()   << "                   ║\n"
        << "║    Barrel: both z fail        : " << std::setw(7) << m_statAnalBarrelBothZfail.load() << "                   ║\n"
        << "║    Barrel: z swapped (ok)     : " << std::setw(7) << m_statAnalBarrelZswapped.load()  << "                   ║\n"
        << "║    Endcap: |tanLam| too small : " << std::setw(7) << m_statAnalEndcapSmallTanL.load() << "                   ║\n"
        << "║    Endcap: already past cap   : " << std::setw(7) << m_statAnalEndcapPastCap.load()   << "                   ║\n"
        << "║    Endcap: arc > 2π R         : " << std::setw(7) << m_statAnalEndcapArcLong.load()   << "                   ║\n"
        << "╠══════════════════════════════════════════════════════════╣\n"
        << "║  PAIR-COINCIDENCE FILTER                                 ║\n"
        << "║    Hits demoted (no partner)  : " << std::setw(7) << m_statHitsNonCoinc.load()      << "                   ║\n"
        << "║    Combos rejected (<2 doubl) : " << std::setw(7) << m_statCombosRejByCoinc.load()   << "                   ║\n"
        << "║    Tracks recovered by pass-2 : " << std::setw(7) << m_statCoincPass2Tracks.load()   << "                   ║\n"
        << "╠══════════════════════════════════════════════════════════╣\n"
        << "║  GENFIT DOWNSTREAM                                       ║\n"
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

    // GenFit's field manager and material effects are singletons, freed at exit
    if (m_useGenFit) {
        m_genFitField.reset();
    }

    return Algorithm::finalize();
}

const dd4hep::rec::Surface* CombinatorialTrackFinder::findSurface(const edm4hep::TrackerHitPlane& hit) const {
    return findSurfaceByID(hit.getCellID());
}

const dd4hep::rec::Surface* CombinatorialTrackFinder::findSurfaceByID(uint64_t cellID) const {
    dd4hep::rec::SurfaceMap::const_iterator sI = m_surfaceMap->find(cellID);

    if (sI != m_surfaceMap->end()) {
        return dynamic_cast<const dd4hep::rec::Surface*>(sI->second);
    }

    return nullptr;
}

// Region of a cellID: 0 barrel, +1 / -1 endcaps (field RegionFieldName).
int CombinatorialTrackFinder::getTypeID(uint64_t cellID) const {
    if (!m_hasRegionField) return 0;   // single-region readout: treat everything as barrel
    return m_bitFieldCoder->get(cellID, m_regionFieldName.value());
}

int CombinatorialTrackFinder::getLayerID(uint64_t cellID) const {
    return m_bitFieldCoder->get(cellID, m_layerFieldName.value());
}

// Layer key of a cellID:
//   barrel:    100*system + layer
//   +endcap:   100*system + layer + 1000
//   -endcap: -(100*system + layer + 1000)
// The system term is only added when several subdetectors share one readout.
int CombinatorialTrackFinder::getCompositeID(uint64_t cellID) const {
    int layerID = getLayerID(cellID);
    int typeID  = getTypeID(cellID);
    int base    = layerID;
    if (m_useSystemInCompositeID) base += 100 * m_bitFieldCoder->get(cellID, "system");
    if      (typeID ==  0) return base;
    else if (typeID ==  1) return  1000 + base;
    else                   return -1000 - base;
}

// Detector surfaces grouped by layer key.
std::map<int, std::vector<const dd4hep::rec::Surface*>> CombinatorialTrackFinder::getSurfacesByLayer() const {
    std::map<int, std::vector<const dd4hep::rec::Surface*>> result;

    for (const auto& surface : m_surfaces) {
        dd4hep::DetElement det = surface->detElement();
        if (!det.isValid()) continue;

        uint64_t volID = det.volumeID();
        int layerID = 0;
        int type = 0;  // 0=barrel, 1=positive endcap, -1=negative endcap


        layerID = getLayerID(volID);
        type = getTypeID(volID);

        // same key as for the hits
        const int compositeID = getCompositeID(volID);
        result[compositeID].push_back(surface);

        if (msgLevel(MSG::DEBUG)) {
                dd4hep::rec::Vector3D pos = surface->origin();
                std::string typeStr;
                if (type == 0) typeStr = "Barrel";
                else if (type == 1) typeStr = "Positive Endcap";
                else if (type == -1) typeStr = "Negative Endcap";
                else typeStr = "Unknown";
        }
    }

    if (msgLevel(MSG::INFO)) {
        info() << "Layer distribution:" << endmsg;
        for (const auto& [layer, surfaces] : result) {
            std::string typeStr;
            int actualLayer;

            if (layer >= 1000) {
                typeStr = "Positive Endcap";
                actualLayer = layer - 1000;
            } else if (layer <= -1000) {
                typeStr = "Negative Endcap";
                actualLayer = -(layer + 1000);
            } else {
                typeStr = "Barrel";
                actualLayer = layer;
            }

            info() << "  " << typeStr << " Layer " << actualLayer
                  << " (ID=" << layer << "): " << surfaces.size() << " surfaces" << endmsg;
        }
    }

    return result;
}

// pT from omega, with Bz at the state's reference point.
double CombinatorialTrackFinder::getPT(const edm4hep::TrackState& state) const {
    double omega = state.omega;  // 1/mm

    double x = state.referencePoint[0];
    double y = state.referencePoint[1];
    double z = state.referencePoint[2];

    dd4hep::Position fieldPos(x/10.0, y/10.0, z/10.0); // mm -> cm
    double bField = m_field.magneticField(fieldPos).z() / dd4hep::tesla;

    double pT = 0.3 * std::abs(bField) / std::abs(omega) * 0.001; // GeV/c

    return pT;
}

// Where along the fitted helix a stored track state is placed.
enum class StatePoint { FirstHit, LastHit, Origin };

// Describe the fitted helix as an EDM4hep track state at one point of the track.
//
// The helix is the fitted circle (centre, radius) plus a straight-line fit of z versus
// the transverse arc length. Everything is geometric, so it works for either field sign:
//   - the particle moves outward, from the hit nearest the origin to the farthest one;
//   - omega > 0 means clockwise seen from +z (EDM4hep); in a field along +z its sign is
//     the charge;
//   - the reference point is the innermost hit, the outermost hit or the origin, and
//     d0, phi, z0 are given at the circle point closest to it (EDM4hep convention:
//     pca = ref + (-d0 sin(phi), d0 cos(phi), z0), centre = pca + (1/omega)(sin(phi), -cos(phi)),
//     z(s) = pca_z + tanLambda * s).
// Units: centre and radius in cm, circleCov is the 3x3 covariance of (x0, y0, R) in cm^2;
// the state is written in mm. Returns false (state unchanged) for degenerate input.
static bool makeEdm4hepStateFromCircle(double centreXcm, double centreYcm, double radiusCm,
                                       const Eigen::Matrix3d& circleCov,
                                       const std::vector<edm4hep::TrackerHitPlane>& trackHits,
                                       StatePoint where, edm4hep::TrackState& state) {
    const double centreX = centreXcm * 10.0;   // mm
    const double centreY = centreYcm * 10.0;
    const double radius  = radiusCm * 10.0;
    if (radius <= 0.0 || trackHits.size() < 2)
        return false;

    // 1. Hits sorted by distance from the origin: the order an outward-moving particle
    //    crosses them.
    struct HitPoint {
        double x, y, z;      // mm
        double distance;     // 3D distance from the origin, mm
        double varianceZ;    // mm^2
    };
    std::vector<HitPoint> points;
    for (const auto& hit : trackHits) {
        const auto& pos = hit.getPosition();
        HitPoint point;
        point.x = pos[0];
        point.y = pos[1];
        point.z = pos[2];
        point.distance = std::sqrt(point.x * point.x + point.y * point.y + point.z * point.z);
        point.varianceZ = hit.getCovMatrix(edm4hep::Cartesian::z, edm4hep::Cartesian::z);
        points.push_back(point);
    }
    std::sort(points.begin(), points.end(),
              [](const HitPoint& a, const HitPoint& b) { return a.distance < b.distance; });

    // 2. Sense of rotation seen from +z (+1 = counter-clockwise): sum of the cross
    //    products of consecutive hits, measured from the circle centre.
    double crossSum = 0.0;
    for (size_t i = 0; i + 1 < points.size(); ++i) {
        const double ax = points[i].x - centreX;
        const double ay = points[i].y - centreY;
        const double bx = points[i + 1].x - centreX;
        const double by = points[i + 1].y - centreY;
        crossSum += ax * by - ay * bx;
    }
    if (crossSum == 0.0)
        return false;
    const double sense = (crossSum > 0.0) ? 1.0 : -1.0;

    // 3. Reference point and the circle point closest to it (pca).
    double refX = 0.0, refY = 0.0, refZ = 0.0;
    size_t nearestHit = 0;   // hit next to the reference point, used to count the arc
    if (where == StatePoint::FirstHit) {
        refX = points.front().x;
        refY = points.front().y;
        refZ = points.front().z;
    } else if (where == StatePoint::LastHit) {
        nearestHit = points.size() - 1;
        refX = points.back().x;
        refY = points.back().y;
        refZ = points.back().z;
    }
    const double toRefX = refX - centreX;
    const double toRefY = refY - centreY;
    const double centreToRef = std::sqrt(toRefX * toRefX + toRefY * toRefY);
    if (centreToRef < 1e-6)
        return false;   // reference point at the circle centre: no closest point
    const double anglePca = std::atan2(toRefY, toRefX);   // angle of the pca around the centre
    const double pcaX = centreX + radius * toRefX / centreToRef;
    const double pcaY = centreY + radius * toRefY / centreToRef;

    // 4. phi = direction of travel at the pca (radius vector turned by +-90 degrees),
    //    omega > 0 for clockwise motion, d0 from pca - ref = (-d0 sin(phi), d0 cos(phi)).
    double phi = anglePca + sense * M_PI / 2.0;
    while (phi > M_PI) phi -= 2.0 * M_PI;
    while (phi <= -M_PI) phi += 2.0 * M_PI;
    const double omega = -sense / radius;   // 1/mm
    const double d0 = -(pcaX - refX) * std::sin(phi) + (pcaY - refY) * std::cos(phi);

    // 5. Angle turned along the direction of travel, hit by hit (no 2 pi jumps), then the
    //    arc length of every hit measured from the pca.
    std::vector<double> turned(points.size());
    double previousAngle = std::atan2(points[0].y - centreY, points[0].x - centreX);
    turned[0] = 0.0;
    for (size_t i = 1; i < points.size(); ++i) {
        const double angle = std::atan2(points[i].y - centreY, points[i].x - centreX);
        double step = sense * (angle - previousAngle);
        while (step > M_PI) step -= 2.0 * M_PI;
        while (step <= -M_PI) step += 2.0 * M_PI;
        turned[i] = turned[i - 1] + step;
        previousAngle = angle;
    }
    const double angleNearest = std::atan2(points[nearestHit].y - centreY, points[nearestHit].x - centreX);
    double pcaStep = sense * (anglePca - angleNearest);
    while (pcaStep > M_PI) pcaStep -= 2.0 * M_PI;
    while (pcaStep <= -M_PI) pcaStep += 2.0 * M_PI;
    const double turnedPca = turned[nearestHit] + pcaStep;

    // 6. Weighted straight-line fit z = zPca + tanLambda * s, weights 1/sigma_z^2 from
    //    the hit covariance (sigma_z = 1 mm for hits without a z error).
    double sumW = 0.0, sumWS = 0.0, sumWSS = 0.0, sumWZ = 0.0, sumWSZ = 0.0;
    for (size_t i = 0; i < points.size(); ++i) {
        double varianceZ = points[i].varianceZ;
        if (!(varianceZ > 0.0))
            varianceZ = 1.0;
        const double w = 1.0 / varianceZ;
        const double s = radius * (turned[i] - turnedPca);
        sumW   += w;
        sumWS  += w * s;
        sumWSS += w * s * s;
        sumWZ  += w * points[i].z;
        sumWSZ += w * s * points[i].z;
    }
    // det >= 0; it is (close to) zero when all hits sit at the same arc length
    const double det = sumW * sumWSS - sumWS * sumWS;
    if (!(det > 1e-12 * sumW * sumWSS))
        return false;
    const double zPca      = (sumWSS * sumWZ - sumWS * sumWSZ) / det;
    const double tanLambda = (sumW * sumWSZ - sumWS * sumWZ) / det;

    // 7. Covariance. d0, phi and omega from the circle-fit covariance (gradients with
    //    respect to x0, y0, R in cm), with the same floors as the GenFit seed;
    //    z0 and tanLambda from the line fit.
    const double D2 = centreToRef * centreToRef / 100.0;   // cm^2
    const Eigen::Vector3d gradD0(-toRefX / centreToRef, -toRefY / centreToRef, -1.0);   // cm/cm
    const Eigen::Vector3d gradPhi(toRefY / 10.0 / D2, -toRefX / 10.0 / D2, 0.0);        // rad/cm
    const Eigen::Vector3d gradOmega(0.0, 0.0, 1.0 / (radiusCm * radiusCm * 10.0));     // (1/mm)/cm
    const double varD0    = std::max(gradD0.dot(circleCov * gradD0) * 100.0, 0.01);    // mm^2
    const double varPhi   = std::max(gradPhi.dot(circleCov * gradPhi), 0.01);
    const double varOmega = std::max(gradOmega.dot(circleCov * gradOmega), std::pow(0.1 * omega, 2));

    using TP = edm4hep::TrackParams;
    for (int i = 0; i < 6; ++i)
        for (int j = 0; j <= i; ++j)
            state.setCovMatrix(0.0f, static_cast<TP>(i), static_cast<TP>(j));
    state.setCovMatrix(static_cast<float>(varD0), TP::d0, TP::d0);
    state.setCovMatrix(static_cast<float>(varPhi), TP::phi, TP::phi);
    state.setCovMatrix(static_cast<float>(varOmega), TP::omega, TP::omega);
    state.setCovMatrix(static_cast<float>(sumWSS / det), TP::z0, TP::z0);
    state.setCovMatrix(static_cast<float>(sumW / det), TP::tanLambda, TP::tanLambda);
    state.setCovMatrix(static_cast<float>(-sumWS / det), TP::z0, TP::tanLambda);

    state.D0        = static_cast<float>(d0);
    state.phi       = static_cast<float>(phi);
    state.omega     = static_cast<float>(omega);
    state.Z0        = static_cast<float>(zPca - refZ);
    state.tanLambda = static_cast<float>(tanLambda);
    state.referencePoint = edm4hep::Vector3f(refX, refY, refZ);
    return true;
}

// Track state with the given parameters and a default diagonal covariance.
edm4hep::TrackState CombinatorialTrackFinder::createTrackState(
    double d0, double phi, double omega, double z0, double tanLambda,
    int location) const {

    edm4hep::TrackState state;

    state.D0 = d0;
    state.phi = phi;
    state.omega = omega;
    state.Z0 = z0;
    state.tanLambda = tanLambda;
    state.location = location;

    // callers that know the fit uncertainties overwrite these
    for (int i = 0; i < 6; ++i)
        for (int j = 0; j <= i; ++j)
            state.setCovMatrix(0.0,
                static_cast<edm4hep::TrackParams>(i),
                static_cast<edm4hep::TrackParams>(j));

    // D0/Z0 in mm^2, phi/tanLambda in rad^2, omega in (1/mm)^2
    state.setCovMatrix(1.0,   edm4hep::TrackParams::d0,        edm4hep::TrackParams::d0);
    state.setCovMatrix(0.01,  edm4hep::TrackParams::phi,       edm4hep::TrackParams::phi);
    state.setCovMatrix(1e-10,  edm4hep::TrackParams::omega,     edm4hep::TrackParams::omega);
    state.setCovMatrix(1.0,   edm4hep::TrackParams::z0,        edm4hep::TrackParams::z0);
    state.setCovMatrix(0.01,  edm4hep::TrackParams::tanLambda, edm4hep::TrackParams::tanLambda);

    return state;
}

//------------------------------------------------------------------------------
// Track finding
//------------------------------------------------------------------------------

void CombinatorialTrackFinder::findTracks(
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

    if (msgLevel(MSG::DEBUG) && hits->size() > 0) {
        const auto& dbgHit = (*hits)[0];
    }

    // Input hit -> all SimTrackerHits that contributed to it (one digi hit can have
    // several). Keyed on (collectionID, index): merged input collections all start
    // their indices at 0.
    auto hitKey = [](const auto& h) -> uint64_t {
        const auto id = h.id();
        return (static_cast<uint64_t>(id.collectionID) << 32) | static_cast<uint32_t>(id.index);
    };
    std::unordered_map<uint64_t, std::vector<edm4hep::SimTrackerHit>> recoToSimMap;
    for (const auto& link : recoSimLinks) {
        auto recoHit = link.getFrom();
        recoToSimMap[hitKey(recoHit)].push_back(link.getTo());
    }

    // DEBUG: the MC particles behind this event's hits
    if (msgLevel(MSG::DEBUG)) {
        debug() << "MC truth particles for this event (from RecoSim links):" << endmsg;
        std::unordered_set<int> seenMCIdx;
        int nLinked = 0;
        for (size_t gi = 0; gi < hits->size(); ++gi) {
            auto it = recoToSimMap.find(hitKey((*hits)[gi]));
            if (it == recoToSimMap.end()) continue;
            nLinked++;
            for (const auto& sim : it->second) {
                auto mcPart = sim.getParticle();
                if (!mcPart.isAvailable()) continue;
                int mcIdx = mcPart.id().index;
                if (!seenMCIdx.insert(mcIdx).second) continue;  // already printed
                const auto& mom = mcPart.getMomentum();
                const auto& vtx = mcPart.getVertex();
                double mcPt  = std::sqrt(mom.x*mom.x + mom.y*mom.y);
                double mcP   = std::sqrt(mcPt*mcPt + mom.z*mom.z);
                double mcEta = (mcPt > 0.0) ? std::asinh(mom.z / mcPt) : 0.0;
                double mcPhi = std::atan2(mom.y, mom.x);
                double vtxR  = std::sqrt(vtx.x*vtx.x + vtx.y*vtx.y) / 10.0;  // mm→cm
                debug() << "  [MC truth] PDG=" << mcPart.getPDG()
                        << "  q=" << mcPart.getCharge()
                        << "  pT=" << std::fixed << std::setprecision(3) << mcPt << " GeV/c"
                        << "  pz=" << mom.z << " GeV/c"
                        << "  |p|=" << mcP << " GeV/c"
                        << "  phi=" << mcPhi << " rad"
                        << "  eta=" << mcEta
                        << "  vtx_R=" << vtxR << " cm"
                        << "  vtx_z=" << vtx.z/10.0 << " cm"
                        << endmsg;
            }
        }
        debug() << "  => " << nLinked << "/" << hits->size() << " hits linked to sim, "
                << seenMCIdx.size() << " unique MC particle(s)" << endmsg;
    }

    // Helper lambda: given a newly created output hit whose data was copied from
    // an input hit, propagate every RecoSim link of the input hit to the output
    // link collection. For cluster digi this preserves multi-SimHit links
    // (muon + delta-rays in the same readout cell) instead of collapsing them
    // to one; baseline digi yields exactly one link per hit as before.
    auto propagateLink = [&](const edm4hep::MutableTrackerHitPlane& outHit,
                              const edm4hep::TrackerHitPlane& inHit) {
        auto it = recoToSimMap.find(hitKey(inHit));
        if (it == recoToSimMap.end()) return;
        for (const auto& sim : it->second) {
            auto outLink = outputLinks.create();
            outLink.setFrom(outHit);
            outLink.setTo(sim);
            outLink.setWeight(1.0f);
        }
    };

    // hits already used by accepted tracks
    std::vector<bool> globalUsedHits(hits->size(), false);

    // SeedAndExtend: all candidates are built in the first iteration, then taken one per
    // iteration, best first, skipping those that share a hit with an accepted track.
    std::vector<k4rt::SECandidate> seCands;
    bool   seBuilt = false;
    size_t seNext  = 0;
    // hit positions (mm) of the accepted tracks, for the neighbour-track cuts
    std::vector<std::array<double, 3>> acceptedTrackHitPositions;
    int trackNumber = 0;

    // One track per iteration, until no candidate is left
    while (true) {
        trackNumber++;
        info() << "\n=== TRACK FINDING ITERATION " << trackNumber << " ===" << endmsg;

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
        if (msgLevel(MSG::DEBUG)) {
            for (const auto& [layer, count] : unusedHitsPerLayer) {
                debug() << "  Layer " << layer << ": " << count << " unused hits" << endmsg;
            }
        }

        if (remainingHits < 3) {
            info() << "Not enough unused hits remaining (" << remainingHits << ") to form another track" << endmsg;
            break;
        }

        // unused hits of this iteration
        using HitInfo = k4rt::HitInfo;

        std::vector<HitInfo> allHitInfo;
        std::map<int, std::vector<size_t>> hitIndicesByCompositeLayer;

        allHitInfo.reserve(remainingHits);

        // MaxHitEdepKeV: skip hits with a large energy deposit (EDM4hep edep is in GeV)
        const double edepThreshGeV = (m_maxHitEdepKeV.value() > 0.0)
                                     ? m_maxHitEdepKeV.value() * 1e-6   // keV → GeV
                                     : -1.0;                              // disabled
        int evtEdepRej = 0;
        for (size_t i = 0; i < hits->size(); ++i) {
            if (globalUsedHits[i]) continue;

            const auto& hit = (*hits)[i];

            if (edepThreshGeV > 0.0 && hit.getEDep() > edepThreshGeV) {
                debug() << "  Hit " << i << " excluded by edep cut: "
                        << hit.getEDep()*1e6 << " keV > " << m_maxHitEdepKeV.value()
                        << " keV" << endmsg;
                ++evtEdepRej;
                continue;
            }

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
        m_statEdepRejected += evtEdepRej;
        if (evtEdepRej > 0)
            debug() << "  edep pre-filter removed " << evtEdepRej << " hits this iteration" << endmsg;

        if (allHitInfo.size() < 3) {
            info() << "Not enough unused hits with valid surfaces (" << allHitInfo.size() << ") for iteration " << trackNumber << endmsg;
            break;
        }

        info() << "Available hits for iteration " << trackNumber << ": " << allHitInfo.size()
               << " in " << hitIndicesByCompositeLayer.size() << " detector different layers/regions" << endmsg;

        std::vector<int> compositeLayerIDs;
        for (const auto& [compositeID, indices] : hitIndicesByCompositeLayer) {
            compositeLayerIDs.push_back(compositeID);
        }

        // layer keys sorted by layer number (innermost first), then by key
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

        // Field at the centroid of the hits, to turn MaxComboPT into a radius cap for the
        // k = 3 score: R_cap [cm] = MaxComboPT [GeV] * 100 / (0.3 * |B| [T])
        double bRefForK3 = 1.7;  // used only if the field lookup fails
        {
            double sx = 0.0, sy = 0.0, sz = 0.0;
            for (const auto& hi : allHitInfo) {
                const auto& p = hi.hit.getPosition();  // mm
                sx += p.x; sy += p.y; sz += p.z;
            }
            const double n = static_cast<double>(allHitInfo.size());
            if (n > 0) {
                dd4hep::Position centroid(sx / n / 10.0, sy / n / 10.0, sz / n / 10.0);  // cm
                try {
                    bRefForK3 = m_field.magneticField(centroid).z() / dd4hep::tesla;
                } catch (...) {}
            }
        }
        const double absBRefForK3 = std::max(std::abs(bRefForK3), 1e-6);
        const double rCapK3 = m_maxComboPT.value() * 100.0 / (0.3 * absBRefForK3);

        // ── Energy-deposit score ──────────────────────────────────────────────────
        // In a layer with several hits, the right one usually has an energy deposit like
        // the hits in layers with a single hit (the event baseline); delta rays and
        // secondaries deposit much more or much less.
        //   hit score  = distance of edep from [baseline, 2 x baseline]  (keV, lower = better)
        //   combo score: k = 3 : w * avgScore / EdepNormKeV - min(R, R_cap) / R_cap
        //                k >= 4: chi2/ndf + w * avgScore / EdepNormKeV
        // MaxHitsPerCompositeID keeps only the best-scoring hits of crowded layers.
        std::unordered_map<size_t, double> hitEdepKeV;   // allHitInfo index → edep (keV)
        hitEdepKeV.reserve(allHitInfo.size());
        for (size_t idx = 0; idx < allHitInfo.size(); ++idx)
            hitEdepKeV[idx] = allHitInfo[idx].hit.getEDep() * 1e6;   // GeV → keV

        // baseline: mean edep of the layers with exactly one hit
        double meanEdepSingle = 0.0;
        int    nSingle        = 0;
        for (const auto& [cid, cidIndices] : hitIndicesByCompositeLayer) {
            if (cidIndices.size() == 1) {
                meanEdepSingle += hitEdepKeV.at(cidIndices[0]);
                ++nSingle;
            }
        }
        if (nSingle > 0) meanEdepSingle /= static_cast<double>(nSingle);

        // no single-hit layer: use the median edep of all hits
        if (nSingle == 0 && !hitEdepKeV.empty()) {
            std::vector<double> allEdeps;
            allEdeps.reserve(hitEdepKeV.size());
            for (const auto& [idx, e] : hitEdepKeV) allEdeps.push_back(e);
            std::sort(allEdeps.begin(), allEdeps.end());
            const size_t mid = allEdeps.size() / 2;
            meanEdepSingle = (allEdeps.size() % 2 == 0)
                           ? 0.5 * (allEdeps[mid-1] + allEdeps[mid])
                           : allEdeps[mid];
            debug() << "  No single-hit compositeIDs — using median edep as baseline: "
                    << meanEdepSingle << " keV" << endmsg;
        }

        std::unordered_map<size_t, double> hitEdepScore;
        hitEdepScore.reserve(allHitInfo.size());
        if (m_proxScoreWeight.value() > 0.0) {
            for (size_t idx = 0; idx < allHitInfo.size(); ++idx) {
                const double edep_    = hitEdepKeV.at(idx);
                const double lowDev_  = std::max(0.0, meanEdepSingle - edep_);
                const double highDev_ = std::max(0.0, edep_ - 2.0 * meanEdepSingle);
                hitEdepScore[idx]     = lowDev_ + highDev_;
                debug() << "  EdepScore hit[" << idx << "] compositeID="
                        << allHitInfo[idx].compositeID
                        << " edep=" << edep_ << " keV  baseline=" << meanEdepSingle << " keV"
                        << " score=" << hitEdepScore[idx] << " keV" << endmsg;
            }

            const int maxHPC = m_maxHitsPerComposite.value();
            if (maxHPC > 0) {
                for (auto& [cid, cidIndices] : hitIndicesByCompositeLayer) {
                    if (static_cast<int>(cidIndices.size()) <= maxHPC) continue;
                    std::sort(cidIndices.begin(), cidIndices.end(),
                              [&](size_t a, size_t b) {
                                  return hitEdepScore[a] < hitEdepScore[b];
                              });
                    cidIndices.resize(maxHPC);
                    debug() << "  compositeID=" << cid << " trimmed to top "
                            << maxHPC << " hits by edep similarity" << endmsg;
                }
            }
        }

        // hits per layer number (barrel and endcaps together)
        std::map<int, std::vector<size_t>> hitsByLayerNumber;
        for (size_t idx = 0; idx < allHitInfo.size(); ++idx) {
            int lyr = allHitInfo[idx].layerID;
            hitsByLayerNumber[lyr].push_back(idx);
        }

        std::vector<bool> usedHits(hits->size(), false);
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

        edm4hep::TrackCollection candidateTracks;  // this iteration only

        // ── N-hit combinatorial search (Strategy = Exhaustive) ────────────────────
        // All C(N,k) combinations for k = MaxCombinatorialHits down to MinTrackHits,
        // each scored by the chi2/ndf of its circle fit, after road, isolation and
        // layer-span cuts; outliers are then removed from the best combination.

        // hit indices ordered by layer, inner -> outer
        std::vector<size_t> orderedHII;
        orderedHII.reserve(allHitInfo.size());
        for (const auto& cid : compositeLayerIDs) {
            for (size_t idx : hitIndicesByCompositeLayer[cid])
                orderedHII.push_back(idx);
        }

        const int nAvail = static_cast<int>(orderedHII.size());
        const int maxK   = std::min(nAvail, m_maxCombinatorialHits.value());
        const int minK   = std::max(3, m_minTrackHits.value());

        using CombResult = k4rt::CombResult;
        CombResult bestCombo;
        CombResult fallbackCombo;  // best over-threshold combo; promoted if no valid combo exists
        bool foundCombination = false;

        // per-iteration counters (debug printout)
        int evtCombosEval    = 0;  // combinations that reached the fit stage
        int evtIsoRej        = 0;  // rejected by isolation cut
        int evtChi2Rej       = 0;  // rejected by chi2/ndf threshold
        int evtFitFail       = 0;  // fit itself returned false
        int evtConsecPhiRej  = 0;  // rejected by MaxConsecDeltaPhi
        int evtPhiSpreadRej  = 0;  // rejected by MaxComboPhiSpread
        int evtConsecDistRej = 0;  // rejected by MaxConsecutiveHitDistance
        int evtPairDistRej   = 0;  // rejected by MaxPairHitDistance
        int evtTanLambdaRej  = 0;  // rejected by MaxComboTanLambdaDev

        // ── Quality cuts ──────────────────────────────────────────────────────────

        // hit index -> layer key
        std::unordered_map<size_t, int> gIdxToComposite;
        gIdxToComposite.reserve(allHitInfo.size());
        for (const auto& hi : allHitInfo) {
            gIdxToComposite[hi.index] = hi.compositeID;
        }

        // Isolation: each hit of the combination must be HitIsolationCut (cm) away from
        // every other unused hit in a different layer. Hits in the same layer (secondaries
        // next to the track) are not counted.
        auto allIsolated = [&](const std::vector<size_t>& gIdxVec) -> bool {
            const double cut = m_hitIsolationCut.value();
            if (cut <= 0.0) return true;
            for (size_t gi : gIdxVec) {
                const auto& refPos = (*hits)[gi].getPosition();
                Eigen::Vector3d ref(refPos[0]/10.0, refPos[1]/10.0, refPos[2]/10.0);
                int giComposite = -999;
                auto cgit = gIdxToComposite.find(gi);
                if (cgit != gIdxToComposite.end()) giComposite = cgit->second;
                for (size_t j = 0; j < hits->size(); ++j) {
                    if (j == gi) continue;
                    if (globalUsedHits[j]) continue;
                    bool inCombo = false;
                    for (size_t si : gIdxVec) { if (si == j) { inCombo = true; break; } }
                    if (inCombo) continue;
                    int jComposite = -999;
                    auto cjit = gIdxToComposite.find(j);
                    if (cjit == gIdxToComposite.end()) continue; // not in this iteration's hit pool
                    jComposite = cjit->second;
                    if (jComposite == giComposite) continue;
                    const auto& p = (*hits)[j].getPosition();
                    Eigen::Vector3d cand(p[0]/10.0, p[1]/10.0, p[2]/10.0);
                    double dist = (cand - ref).norm();
                    if (dist < cut) {
                        debug() << "    isolation FAIL: selected hit gIdx=" << gi
                                << " compositeID=" << giComposite
                                << " too close to neighbour gIdx=" << j
                                << " compositeID=" << jComposite
                                << "  dist=" << dist << " cm  cut=" << cut << " cm" << endmsg;
                        return false;
                    }
                }
            }
            return true;
        };

        // Road cuts: consecutive-layer Δφ and total φ spread across the combo.
        // Both cuts reject combinations whose hits span an azimuthal window too wide
        // to belong to a single track.  Hits in gIdxVec are ordered inner→outer
        // (slot order matches compositeLayerIDs which is sorted by layer number).
        // Hit-set predicates live in TrackingPrimitives so the seed finder can reuse them.
        // (the debug stream is fetched per call: a cached one keeps the level of the last message)
        const k4rt::RoadCutConfig roadCfg{m_maxConsecDeltaPhi.value(), m_maxComboPhiSpread.value(),
                                          m_maxConsecHitDist.value(), m_maxPairHitDist.value(),
                                          m_maxComboTanLambdaDev.value()};
        auto roadCtx = [&]() {
            return k4rt::RoadCutContext{*hits, gIdxToComposite, roadCfg,
                                        msgLevel(MSG::DEBUG) ? &debug() : nullptr, &m_cidRegion};
        };
        auto passesRoadCuts = [&](const std::vector<size_t>& gIdxVec) -> bool {
            return k4rt::passesRoadCuts(gIdxVec, roadCtx());
        };
        auto passesDistanceCuts = [&](const std::vector<size_t>& gIdxVec) -> bool {
            return k4rt::passesDistanceCuts(gIdxVec, roadCtx());
        };
        auto passesTanLambdaCut = [&](const std::vector<size_t>& gIdxVec) -> bool {
            return k4rt::passesTanLambdaCut(gIdxVec, roadCtx());
        };


        // ── Pair coincidence (RequirePairCoincidence), for detectors built of layer
        // doublets (0,1), (2,3), (4,5): two hits of a doublet are coincident when
        // |dphi| < 0.15 rad, |deta| < 0.10 and their distance is below PairCoincidenceMaxDistMM.
        constexpr double kPairTolPhi    = 0.15;
        constexpr double kPairTolEta    = 0.10;
        const double kPairMaxDistMM = m_pairCoincMaxDistMM.value();

        auto coincidentHits = [&](const edm4hep::TrackerHitPlane& a,
                                  const edm4hep::TrackerHitPlane& b) -> bool {
            const auto& pa = a.getPosition();
            const auto& pb = b.getPosition();
            double dx = pb[0] - pa[0];
            double dy = pb[1] - pa[1];
            double dz = pb[2] - pa[2];
            if ((dx*dx + dy*dy + dz*dz) > (kPairMaxDistMM * kPairMaxDistMM))
                return false;
            double dPhi = std::atan2(pb[1], pb[0]) - std::atan2(pa[1], pa[0]);
            while (dPhi >  M_PI) dPhi -= 2.0*M_PI;
            while (dPhi < -M_PI) dPhi += 2.0*M_PI;
            if (std::abs(dPhi) > kPairTolPhi) return false;
            double rA = std::sqrt(pa[0]*pa[0] + pa[1]*pa[1]);
            double rB = std::sqrt(pb[0]*pb[0] + pb[1]*pb[1]);
            if (rA < 1e-6 || rB < 1e-6) return false;
            double etaA = std::asinh(pa[2] / rA);
            double etaB = std::asinh(pb[2] / rB);
            return std::abs(etaB - etaA) <= kPairTolEta;
        };

        // partner layer of a doublet: 0 <-> 1, 2 <-> 3, 4 <-> 5
        auto pairPartnerCid = [](int composite) -> int {
            int absC = std::abs(composite);
            int layer = absC % 1000;
            int partnerLayer = layer ^ 1;
            int absPartner = (absC >= 1000) ? (1000 + partnerLayer) : partnerLayer;
            return (composite < 0) ? -absPartner : absPartner;
        };

        // Does hit at allHitInfo[hii] have any coincident partner in its paired layer?
        auto hasPartnerInPair = [&](size_t hii) -> bool {
            int partnerCid = pairPartnerCid(allHitInfo[hii].compositeID);
            auto it = hitIndicesByCompositeLayer.find(partnerCid);
            if (it == hitIndicesByCompositeLayer.end()) return false;
            for (size_t partnerHii : it->second) {
                if (coincidentHits(allHitInfo[hii].hit, allHitInfo[partnerHii].hit))
                    return true;
            }
            return false;
        };

        // Count complete doublets in a combo: a doublet is complete when the combo
        // contains hits in both layers of one of (0,1)/(2,3)/(4,5) AND they are
        // mutually coincident, both in the same region.
        auto countCompleteDoublets = [&](const std::vector<size_t>& gIdxVec) -> int {
            int n_complete = 0;
            const int pairs[3][2] = {{0,1}, {2,3}, {4,5}};
            for (int p = 0; p < 3; ++p) {
                bool complete = false;
                for (size_t gi_a : gIdxVec) {
                    auto itA = gIdxToComposite.find(gi_a);
                    if (itA == gIdxToComposite.end()) continue;
                    int cidA = itA->second;
                    int absA = std::abs(cidA);
                    if ((absA % 1000) != pairs[p][0]) continue;
                    int regionA = (cidA < 0) ? -1 : (absA >= 1000 ? 1 : 0);
                    for (size_t gi_b : gIdxVec) {
                        if (gi_b == gi_a) continue;
                        auto itB = gIdxToComposite.find(gi_b);
                        if (itB == gIdxToComposite.end()) continue;
                        int cidB = itB->second;
                        int absB = std::abs(cidB);
                        if ((absB % 1000) != pairs[p][1]) continue;
                        int regionB = (cidB < 0) ? -1 : (absB >= 1000 ? 1 : 0);
                        if (regionA != regionB) continue;
                        if (coincidentHits((*hits)[gi_a], (*hits)[gi_b])) {
                            complete = true; break;
                        }
                    }
                    if (complete) break;
                }
                if (complete) ++n_complete;
            }
            return n_complete;
        };

        // unused hits of each layer: hitsForCompAll[i] belongs to compositeLayerIDs[i]
        const int nDistinctAll = static_cast<int>(compositeLayerIDs.size());
        std::vector<std::vector<size_t>> hitsForCompAll(nDistinctAll);
        for (int ci = 0; ci < nDistinctAll; ++ci) {
            int cid = compositeLayerIDs[ci];
            for (size_t hii : hitIndicesByCompositeLayer.at(cid)) {
                if (!usedHits[allHitInfo[hii].index]) {
                    hitsForCompAll[ci].push_back(hii);
                }
            }
        }

        // Pair coincidence, hard mode: first pass only with hits that have a coincident
        // partner (hitsForCompSeed). Soft mode keeps all hits and ranks by doublets.
        const bool pcHard = m_requirePairCoinc.value() && !m_pairCoincSoftScoring.value();
        const bool pcSoft = m_requirePairCoinc.value() &&  m_pairCoincSoftScoring.value();
        std::vector<std::vector<size_t>> hitsForCompSeed;
        if (pcHard) {
            hitsForCompSeed.assign(nDistinctAll, std::vector<size_t>{});
            int evtNonCoinc = 0;
            for (int ci = 0; ci < nDistinctAll; ++ci) {
                for (size_t hii : hitsForCompAll[ci]) {
                    if (hasPartnerInPair(hii)) {
                        hitsForCompSeed[ci].push_back(hii);
                    } else {
                        ++evtNonCoinc;
                    }
                }
            }
            m_statHitsNonCoinc += evtNonCoinc;
            debug() << "  PairCoinc: " << evtNonCoinc << " hits demoted (no partner in paired layer)" << endmsg;
        }

        // ── Hit clustering (DoHitClustering) ──────────────────────────────────────
        // Hits of neighbouring layers that pass the road cuts are joined (union-find);
        // each cluster is searched on its own. Without clustering: one cluster, all hits.
        struct ClusterDef {
            std::vector<int>                 cids;
            std::vector<std::vector<size_t>> hitsForComp;
        };
        std::vector<ClusterDef> activeClusters;

        if (m_doHitClustering.value() && nDistinctAll >= 2) {
            const size_t nhii = allHitInfo.size();
            std::vector<size_t> ufParent(nhii);
            for (size_t i = 0; i < nhii; ++i) ufParent[i] = i;

            auto ufFind = [&](size_t x) -> size_t {
                while (ufParent[x] != x) {
                    ufParent[x] = ufParent[ufParent[x]]; // path halving
                    x = ufParent[x];
                }
                return x;
            };
            auto ufUnite = [&](size_t a, size_t b) {
                a = ufFind(a); b = ufFind(b);
                if (a != b) ufParent[a] = b;
            };

            const double phiCut  = m_maxConsecDeltaPhi.value();
            const double distCut = m_maxConsecHitDist.value(); // cm

            // join hits of neighbouring layers
            for (int ci = 0; ci < nDistinctAll - 1; ++ci) {
                for (size_t hii_a : hitsForCompAll[ci]) {
                    const auto& pa   = (*hits)[allHitInfo[hii_a].index].getPosition();
                    const double phiA = std::atan2(pa[1], pa[0]);
                    for (size_t hii_b : hitsForCompAll[ci + 1]) {
                        const auto& pb = (*hits)[allHitInfo[hii_b].index].getPosition();
                        if (phiCut > 0.0) {
                            double dPhi = std::atan2(pb[1], pb[0]) - phiA;
                            while (dPhi >  M_PI) dPhi -= 2*M_PI;
                            while (dPhi < -M_PI) dPhi += 2*M_PI;
                            if (std::abs(dPhi) > phiCut) continue;
                        }
                        if (distCut > 0.0) {
                            double dx=pb[0]-pa[0], dy=pb[1]-pa[1], dz=pb[2]-pa[2];
                            if (std::sqrt(dx*dx+dy*dy+dz*dz)*0.1 > distCut) continue;
                        }
                        ufUnite(hii_a, hii_b);
                    }
                }
            }

            // also across one missing layer, with doubled cuts
            const double phiCut2  = (phiCut  > 0.0) ? phiCut  * 2.0 : -1.0;
            const double distCut2 = (distCut > 0.0) ? distCut * 2.0 : -1.0;
            for (int ci = 0; ci < nDistinctAll - 2; ++ci) {
                for (size_t hii_a : hitsForCompAll[ci]) {
                    const auto& pa    = (*hits)[allHitInfo[hii_a].index].getPosition();
                    const double phiA = std::atan2(pa[1], pa[0]);
                    for (size_t hii_b : hitsForCompAll[ci + 2]) {
                        const auto& pb = (*hits)[allHitInfo[hii_b].index].getPosition();
                        if (phiCut2 > 0.0) {
                            double dPhi = std::atan2(pb[1], pb[0]) - phiA;
                            while (dPhi >  M_PI) dPhi -= 2*M_PI;
                            while (dPhi < -M_PI) dPhi += 2*M_PI;
                            if (std::abs(dPhi) > phiCut2) continue;
                        }
                        if (distCut2 > 0.0) {
                            double dx=pb[0]-pa[0], dy=pb[1]-pa[1], dz=pb[2]-pa[2];
                            if (std::sqrt(dx*dx+dy*dy+dz*dz)*0.1 > distCut2) continue;
                        }
                        ufUnite(hii_a, hii_b);
                    }
                }
            }

            std::unordered_map<size_t, std::vector<size_t>> rootToHIIs;
            for (int ci = 0; ci < nDistinctAll; ++ci) {
                for (size_t hii : hitsForCompAll[ci]) {
                    rootToHIIs[ufFind(hii)].push_back(hii);
                }
            }

            // keep clusters with at least minK layers, layers in inner -> outer order
            for (auto& [root, hiis] : rootToHIIs) {
                std::unordered_map<int, std::vector<size_t>> cidMap;
                for (size_t hii : hiis) cidMap[allHitInfo[hii].compositeID].push_back(hii);
                if (static_cast<int>(cidMap.size()) < minK) {
                    debug() << "  cluster root=" << root << " spans only "
                            << cidMap.size() << " layer(s) — skipped (minK=" << minK << ")" << endmsg;
                    continue;
                }
                ClusterDef cd;
                for (int cid : compositeLayerIDs) {
                    auto it = cidMap.find(cid);
                    if (it == cidMap.end()) continue;
                    cd.cids.push_back(cid);
                    cd.hitsForComp.push_back(std::move(it->second));
                }
                debug() << "  cluster root=" << root
                        << " layers=" << cd.cids.size()
                        << " hits=" << hiis.size() << endmsg;
                activeClusters.push_back(std::move(cd));
            }

            m_statClustersBuilt += static_cast<int>(activeClusters.size());
            info() << "  Hit clustering [trk " << trackNumber << "]: "
                   << nDistinctAll << " layer(s) → "
                   << activeClusters.size() << " cluster(s)" << endmsg;
        } else {
            ClusterDef cd;
            cd.cids        = compositeLayerIDs;
            cd.hitsForComp = hitsForCompAll;
            activeClusters.push_back(std::move(cd));
        }

        // hard pair-coincidence mode: the first pass uses only the coincident hits
        std::vector<ClusterDef> activeClustersSeed;
        if (pcHard) {
            ClusterDef cd;
            cd.cids        = compositeLayerIDs;
            cd.hitsForComp = hitsForCompSeed;
            activeClustersSeed.push_back(std::move(cd));
        }

        // passes: hard pair coincidence = coincident hits first, then all hits if the best
        // combination has < 3 complete doublets; otherwise one pass over all hits
        const int nPairPasses = pcHard ? 2 : 1;

        // ── SeedAndExtend: the next candidate whose hits are all free ─────────────
        // Everything after the choice of hits is shared with the exhaustive search.
        if (m_useSeedExtend) {
            if (!seBuilt) {
                std::vector<std::pair<size_t, int>> pool;
                pool.reserve(allHitInfo.size());
                for (const auto& hi : allHitInfo) pool.emplace_back(hi.index, hi.compositeID);
                k4rt::SeedStats ss;
                const auto t0 = std::chrono::steady_clock::now();
                seCands = m_seedFinder->run(*hits, pool, ss);
                const double ms = std::chrono::duration<double, std::milli>(
                                      std::chrono::steady_clock::now() - t0).count();
                seBuilt = true;
                if (m_seBestFirst)
                    info() << "SeedAndExtend BestFirst: " << ss.extended << " seeds extended, " << ss.skippedClaimed
                           << " skipped (hits claimed), " << ss.rejectedTracks << " rejected (" << ss.rejectedHoles
                           << " for holes, " << ss.rejectedShortPt << " short with low pT); confirmation +" << ss.confirmAdded
                           << " hits, " << ss.confirmSwapped << " swapped, " << ss.confirmRemoved << " removed; " << ss.merged
                           << " pieces merged, " << ss.loopers << " loopers, " << ss.shadowClaimed << " shadow hits claimed" << endmsg;
                info() << "SeedAndExtend: " << pool.size() << " hits, " << ss.rawTriplets << " raw combinations, "
                       << ss.seeds << " seeds ("
                       << ss.seedsAfterZ << " after z check, " << ss.extended << " extended) -> "
                       << ss.candidates << " candidates in " << ms << " ms"
                       << (ss.hitSeedCap ? "  [MaxSeedsPerEvent reached]" : "") << endmsg;
            }
            std::unordered_map<size_t, size_t> gToHii;
            for (size_t i = 0; i < allHitInfo.size(); ++i) gToHii[allHitInfo[i].index] = i;
            while (seNext < seCands.size()) {
                const auto& c = seCands[seNext++];
                bool free = true;
                for (size_t g : c.gIdx)
                    if (usedHits[g] || !gToHii.count(g)) { free = false; break; }
                if (!free) continue;
                bestCombo = CombResult{};
                bestCombo.gIdxVec = c.gIdx;
                for (size_t g : c.gIdx) bestCombo.hiiVec.push_back(gToHii[g]);
                bestCombo.chi2ndf   = c.chi2ndf;
                bestCombo.avgProxCm = 0.0;
                bestCombo.x0 = c.x0; bestCombo.y0 = c.y0; bestCombo.radius = c.R;
                bestCombo.fitCov    = c.cov;
                bestCombo.residuals = c.residuals;
                bestCombo.nHits     = static_cast<int>(c.gIdx.size());
                bestCombo.kinkChi2  = c.kinkChi2;
                bestCombo.kinkNdf   = c.kinkNdf;
                bestCombo.looper    = c.looper;
                foundCombination = true;
                break;
            }
        }

        for (int pairPassIdx = 0; !m_useSeedExtend && pairPassIdx < nPairPasses; ++pairPassIdx) {
            // second pass only if the first one's best has fewer than 3 complete doublets
            if (pairPassIdx == 1 && foundCombination) {
                int nDoub = countCompleteDoublets(bestCombo.gIdxVec);
                if (nDoub >= 3) break;
                debug() << "  PairCoinc: pass-1 best has " << nDoub
                        << " complete doublets — running pass-2 over all hits"
                        << endmsg;
            }
            bool foundBeforeThisPass = foundCombination;
            int  doubletsBeforeThisPass = foundCombination
                ? countCompleteDoublets(bestCombo.gIdxVec) : -1;
            const std::vector<ClusterDef>& clustersThisPass =
                (pcHard && pairPassIdx == 0) ? activeClustersSeed : activeClusters;

        // the names below shadow the event-wide ones with this cluster's layers and hits
        for (const auto& cluster : clustersThisPass) {
        const std::vector<int>&                 compositeLayerIDs = cluster.cids;        // NOLINT(shadow)
        const std::vector<std::vector<size_t>>& hitsForComp       = cluster.hitsForComp; // NOLINT(shadow)
        const int nDistinct = static_cast<int>(compositeLayerIDs.size());
        const int maxKeff   = std::min(maxK, nDistinct);

        // For k = maxKeff ... minK: every choice of k layers, and for each every
        // combination of one hit per chosen layer.
        for (int k = maxKeff; k >= minK; --k) {
            // Enumerate C(nDistinct, k) — indices into compositeLayerIDs
            std::vector<int> cidxCombo(k);
            for (int i = 0; i < k; ++i) cidxCombo[i] = i;

            while (true) {
                // the hits of each chosen layer
                bool anySlotEmpty = false;
                std::vector<std::vector<size_t>> slotHits(k);
                for (int i = 0; i < k; ++i) {
                    slotHits[i] = hitsForComp[cidxCombo[i]];
                    if (slotHits[i].empty()) { anySlotEmpty = true; break; }
                }

                if (!anySlotEmpty) {
                    // every combination of one hit per layer
                    std::vector<int> pos(k, 0);
                    while (true) {
                        std::vector<size_t> hiiVec(k), gIdxVec(k);
                        for (int i = 0; i < k; ++i) {
                            hiiVec[i]  = slotHits[i][pos[i]];
                            gIdxVec[i] = allHitInfo[hiiVec[i]].index;
                        }

                        bool anyUsed = false;
                        for (size_t gi : gIdxVec) {
                            if (usedHits[gi]) { anyUsed = true; break; }
                        }

                        if (!anyUsed) {
                            auto comboLabel = [&]() -> std::string {
                                std::string s = "[";
                                for (int ii = 0; ii < k; ++ii) {
                                    if (ii) s += ",";
                                    s += std::to_string(gIdxVec[ii]);
                                }
                                return s + "]";
                            };

                            if (!passesRoadCuts(gIdxVec)) {
                                // which road cut failed (for the counters): repeat the
                                // consecutive-dphi cut, with the same exceptions as passesRoadCuts
                                bool consecFail = false;
                                {
                                    const int nH = static_cast<int>(gIdxVec.size());
                                    double sumR = 0.0, sumAbsZ = 0.0;
                                    for (int i = 0; i < nH; ++i) {
                                        const auto& p = (*hits)[gIdxVec[i]].getPosition();
                                        sumR    += std::sqrt(p[0]*p[0] + p[1]*p[1]);
                                        sumAbsZ += std::abs(p[2]);
                                    }
                                    const double thetaCombo = std::atan2(sumR / nH, sumAbsZ / nH);
                                    constexpr double kNearBeamLo = 15.0 * M_PI / 180.0;
                                    constexpr double kNearBeamHi = 160.0 * M_PI / 180.0;
                                    const bool nearBeamPipe = (thetaCombo < kNearBeamLo || thetaCombo > kNearBeamHi);
                                    const double maxConsec = m_maxConsecDeltaPhi.value();
                                    if (!nearBeamPipe && maxConsec > 0.0) {
                                        for (int ri = 0; ri < nH - 1 && !consecFail; ++ri) {
                                            int cidA = 0, cidB = 0;
                                            auto itA = gIdxToComposite.find(gIdxVec[ri]);
                                            auto itB = gIdxToComposite.find(gIdxVec[ri+1]);
                                            if (itA != gIdxToComposite.end()) cidA = itA->second;
                                            if (itB != gIdxToComposite.end()) cidB = itB->second;
                                            const bool aBarrel = (k4rt::regionOfCid(cidA, &m_cidRegion) == 0);
                                            const bool bBarrel = (k4rt::regionOfCid(cidB, &m_cidRegion) == 0);
                                            if (aBarrel != bBarrel) continue; // barrel-endcap step: no cut
                                            const auto& pa = (*hits)[gIdxVec[ri]].getPosition();
                                            const auto& pb = (*hits)[gIdxVec[ri+1]].getPosition();
                                            double dPhi = std::atan2(pb[1],pb[0]) - std::atan2(pa[1],pa[0]);
                                            while (dPhi >  M_PI) dPhi -= 2*M_PI;
                                            while (dPhi < -M_PI) dPhi += 2*M_PI;
                                            if (std::abs(dPhi) > maxConsec) consecFail = true;
                                        }
                                    }
                                }
                                if (consecFail) {
                                    m_statConsecPhiRejected++;
                                    evtConsecPhiRej++;
                                    debug() << "  k=" << k << " " << comboLabel()
                                            << " -> REJECT consec Δφ" << endmsg;
                                } else {
                                    m_statPhiSpreadRejected++;
                                    evtPhiSpreadRej++;
                                    debug() << "  k=" << k << " " << comboLabel()
                                            << " -> REJECT φ spread" << endmsg;
                                }
                            } else if (!passesDistanceCuts(gIdxVec)) {
                                // which distance cut failed (for the counters)
                                bool consecFail = false;
                                {
                                    const double mc = m_maxConsecHitDist.value();
                                    if (mc > 0.0) {
                                        const int nH = static_cast<int>(gIdxVec.size());
                                        for (int ri = 0; ri < nH - 1 && !consecFail; ++ri) {
                                            const auto& pa = (*hits)[gIdxVec[ri]].getPosition();
                                            const auto& pb = (*hits)[gIdxVec[ri+1]].getPosition();
                                            double dx=pb[0]-pa[0], dy=pb[1]-pa[1], dz=pb[2]-pa[2];
                                            if (std::sqrt(dx*dx+dy*dy+dz*dz)*0.1 > mc)
                                                consecFail = true;
                                        }
                                    }
                                }
                                if (consecFail) {
                                    m_statConsecDistRejected++;
                                    evtConsecDistRej++;
                                    debug() << "  k=" << k << " " << comboLabel()
                                            << " -> REJECT consec dist" << endmsg;
                                } else {
                                    m_statPairDistRejected++;
                                    evtPairDistRej++;
                                    debug() << "  k=" << k << " " << comboLabel()
                                            << " -> REJECT pair dist" << endmsg;
                                }
                            } else if (!allIsolated(gIdxVec)) {
                                m_statIsolationRejected++;
                                evtIsoRej++;
                                debug() << "  k=" << k << " " << comboLabel()
                                        << " -> REJECT isolation" << endmsg;
                            } else if (!passesTanLambdaCut(gIdxVec)) {
                                m_statTanLambdaRejected++;
                                evtTanLambdaRej++;
                                debug() << "  k=" << k << " " << comboLabel()
                                        << " -> REJECT tanLambda dev" << endmsg;
                            } else if (pcHard
                                       && countCompleteDoublets(gIdxVec) < 2) {
                                // hard pair coincidence: at least 2 complete doublets
                                m_statCombosRejByCoinc++;
                                debug() << "  k=" << k << " " << comboLabel()
                                        << " -> REJECT <2 complete doublets" << endmsg;
                            } else {
                                m_statTotalTripletCombos++;
                                evtCombosEval++;
                                std::vector<edm4hep::TrackerHitPlane> comboHits;
                                comboHits.reserve(k);
                                for (size_t gi : gIdxVec) comboHits.push_back((*hits)[gi]);

                                double cx0, cy0, cr, cchi2ndf;
                                Eigen::Matrix3d cCov = Eigen::Matrix3d::Zero();
                                std::vector<double> cRes;
                                bool fitOK = false;
                                try {
                                    fitOK = fitCircleNHits(comboHits, cx0, cy0, cr, cchi2ndf, cCov, cRes);
                                } catch (...) {}

                                if (!fitOK) {
                                    evtFitFail++;
                                    debug() << "  k=" << k << " " << comboLabel()
                                            << " -> REJECT fit failed" << endmsg;
                                } else if (cchi2ndf >= m_nHitMaxChi2NDF.value()) {
                                    evtChi2Rej++;
                                    debug() << "  k=" << k << " " << comboLabel()
                                            << " pT=" << std::fixed << std::setprecision(1)
                                            << (0.3 * absBRefForK3 * cr / 100.0) << " GeV"
                                            << " chi2/ndf=" << cchi2ndf
                                            << " -> REJECT chi2 threshold="
                                            << m_nHitMaxChi2NDF.value() << endmsg;
                                } else if ((0.3 * absBRefForK3 * cr / 100.0) > m_maxComboPT.value()) {
                                    // pT above MaxComboPT: not a candidate, but kept as fallback in
                                    // case every combination is above it (near-straight tracks)
                                    m_statComboPTRejected++;
                                    const double comboPTfb = 0.3 * absBRefForK3 * cr / 100.0;
                                    debug() << "  k=" << k << " " << comboLabel()
                                            << " -> REJECT pT=" << comboPTfb
                                            << " > MaxComboPT=" << m_maxComboPT.value()
                                            << " (stored as fallback)" << endmsg;
                                    // best fallback: more hits; then lower chi2/ndf (k >= 4) or
                                    // lower pT (k = 3)
                                    bool isBetterFallback = false;
                                    if (k != fallbackCombo.nHits) {
                                        isBetterFallback = (k > fallbackCombo.nHits);
                                    } else if (k == 3) {
                                        isBetterFallback = (cr < fallbackCombo.radius);
                                    } else {
                                        isBetterFallback = (cchi2ndf < fallbackCombo.chi2ndf);
                                    }
                                    if (isBetterFallback) {
                                        fallbackCombo.hiiVec    = hiiVec;
                                        fallbackCombo.gIdxVec   = gIdxVec;
                                        fallbackCombo.chi2ndf   = cchi2ndf;
                                        fallbackCombo.avgProxCm = 0.0;
                                        fallbackCombo.x0        = cx0;
                                        fallbackCombo.y0        = cy0;
                                        fallbackCombo.radius    = cr;
                                        fallbackCombo.fitCov    = cCov;
                                        fallbackCombo.residuals = cRes;
                                        fallbackCombo.nHits     = k;
                                    }
                                } else if (m_minComboPT.value() > 0.0
                                           && (0.3 * absBRefForK3 * cr / 100.0) < m_minComboPT.value()) {
                                    // pT below MinComboPT: rejected (no fallback)
                                    m_statComboMinPTRejected++;
                                    debug() << "  k=" << k << " " << comboLabel()
                                            << " -> REJECT pT=" << (0.3 * absBRefForK3 * cr / 100.0)
                                            << " < MinComboPT=" << m_minComboPT.value() << endmsg;
                                } else {
                                    // Ranking of the combinations (lower score is better):
                                    //  1. more hits always wins;
                                    //  2. k >= 4: chi2/ndf + w * edep score / EdepNormKeV;
                                    //  3. k = 3 (chi2 is always 0 for 3 points):
                                    //     w * edep score / EdepNormKeV - min(R, R_cap) / R_cap,
                                    //     i.e. straighter (higher pT) and MIP-like is better.
                                    // The edep score is the mean hit score against the event
                                    // baseline (see above).
                                    double cEdepDev = 0.0;
                                    const double w = m_proxScoreWeight.value();
                                    if (w > 0.0 && !hitEdepScore.empty()) {
                                        for (size_t hii_i : hiiVec)
                                            cEdepDev += hitEdepScore.at(hii_i);
                                        cEdepDev /= static_cast<double>(k);
                                    }
                                    const double edepNorm = std::max(m_edepNormKeV.value(), 1e-9);
                                    bool isBetter;
                                    // soft pair coincidence: more complete doublets wins first
                                    int cDoub = -1, bDoub = -1;
                                    if (pcSoft) {
                                        cDoub = countCompleteDoublets(gIdxVec);
                                        bDoub = (bestCombo.nHits > 0)
                                            ? countCompleteDoublets(bestCombo.gIdxVec) : -1;
                                    }
                                    if (pcSoft && cDoub != bDoub) {
                                        isBetter = (cDoub > bDoub);
                                    } else if (k != bestCombo.nHits) {
                                        isBetter = (k > bestCombo.nHits);
                                    } else if (k == 3) {
                                        // k = 3 (radius term saturates at R_cap; w = 0: radius only)
                                        if (w > 0.0) {
                                            const double cRadNorm = std::min(cr, rCapK3) / rCapK3;
                                            const double bRadNorm = std::min(bestCombo.radius, rCapK3) / rCapK3;
                                            const double cScore3  = w * cEdepDev / edepNorm - cRadNorm;
                                            const double bScore3  = w * bestCombo.avgProxCm / edepNorm - bRadNorm;
                                            isBetter = (cScore3 < bScore3);
                                        } else {
                                            isBetter = (cr > bestCombo.radius);
                                        }
                                    } else {
                                        const double cScore = cchi2ndf
                                            + (w > 0.0 ? w * cEdepDev / edepNorm : 0.0);
                                        const double bScore = bestCombo.chi2ndf
                                            + (w > 0.0 ? w * bestCombo.avgProxCm / edepNorm : 0.0);
                                        isBetter = (cScore < bScore);
                                    }
                                    const double comboPT = 0.3 * absBRefForK3 * cr / 100.0;
                                    if (isBetter) {
                                        bestCombo.hiiVec    = hiiVec;
                                        bestCombo.gIdxVec   = gIdxVec;
                                        bestCombo.chi2ndf   = cchi2ndf;
                                        bestCombo.avgProxCm = cEdepDev;
                                        bestCombo.x0        = cx0;
                                        bestCombo.y0        = cy0;
                                        bestCombo.radius    = cr;
                                        bestCombo.fitCov    = cCov;
                                        bestCombo.residuals = cRes;
                                        bestCombo.nHits     = k;
                                        foundCombination    = true;
                                        debug() << "  k=" << k << " " << comboLabel()
                                                << " pT=" << std::fixed << std::setprecision(1)
                                                << comboPT << " GeV"
                                                << " chi2/ndf=" << cchi2ndf
                                                << " -> NEW BEST" << endmsg;
                                    } else {
                                        debug() << "  k=" << k << " " << comboLabel()
                                                << " pT=" << std::fixed << std::setprecision(1)
                                                << comboPT << " GeV"
                                                << " chi2/ndf=" << cchi2ndf
                                                << " -> not better (best k="
                                                << bestCombo.nHits << " chi2/ndf="
                                                << bestCombo.chi2ndf << ")" << endmsg;
                                    }
                                }
                            }
                        }

                        // next combination of hits
                        int carry = k - 1;
                        while (carry >= 0) {
                            pos[carry]++;
                            if (pos[carry] < static_cast<int>(slotHits[carry].size())) break;
                            pos[carry] = 0;
                            carry--;
                        }
                        if (carry < 0) break;
                    }
                }

                // next choice of k layers
                int i = k - 1;
                while (i >= 0 && cidxCombo[i] == nDistinct - k + i) --i;
                if (i < 0) break;
                ++cidxCombo[i];
                for (int j = i+1; j < k; ++j) cidxCombo[j] = cidxCombo[j-1] + 1;
            }
        }
        } // end cluster loop
            // run statistics: tracks found or improved by the second pass
            if (pairPassIdx == 1 && foundCombination) {
                int nDoubAfter = countCompleteDoublets(bestCombo.gIdxVec);
                if (!foundBeforeThisPass) {
                    m_statCoincPass2Tracks++;
                } else if (nDoubAfter > doubletsBeforeThisPass) {
                    m_statCoincPass2Tracks++;
                }
            }
        } // end pair-coincidence two-pass

        const int maxKeffAll = std::min(maxK, nDistinctAll);
        info() << "  N-hit search [trk " << trackNumber << "]: "
               << nAvail << " avail hits in " << nDistinctAll << " regions, k=" << maxKeffAll << "→" << minK
               << " | evaluated=" << evtCombosEval
               << " consecPhiRej=" << evtConsecPhiRej
               << " phiSpreadRej=" << evtPhiSpreadRej
               << " consecDistRej=" << evtConsecDistRej
               << " pairDistRej=" << evtPairDistRej
               << " isoRej=" << evtIsoRej
               << " tanLRej=" << evtTanLambdaRej
               << " chi2Rej=" << evtChi2Rej
               << " fitFail=" << evtFitFail
               << " edepRej=" << evtEdepRej
               << (foundCombination
                   ? (" | BEST k=" + std::to_string(bestCombo.nHits)
                      + " chi2/ndf=" + std::to_string(bestCombo.chi2ndf)
                      + " edepDev=" + std::to_string(bestCombo.avgProxCm) + " keV")
                   : " | NO COMBINATION FOUND")
               << endmsg;

        if (!foundCombination) {
            if (fallbackCombo.nHits > 0) {
                // every combination was above MaxComboPT: use the fallback
                bestCombo        = fallbackCombo;
                foundCombination = true;
                info() << "  Promoted fallback combo: k=" << bestCombo.nHits
                       << " chi2/ndf=" << bestCombo.chi2ndf
                       << " pT=" << std::fixed << std::setprecision(1)
                       << (0.3 * absBRefForK3 * bestCombo.radius / 100.0)
                       << " GeV (all combos exceeded MaxComboPT="
                       << m_maxComboPT.value() << " GeV)" << endmsg;
            } else {
                info() << "No valid N-hit combination found for track " << trackNumber << " — stopping" << endmsg;
                break;
            }
        }


        // ── Guided hit choice in crowded layers (doCrowdedLayerHitSelection) ──────
        // Layers with one hit are "clean". In each crowded layer the hit is chosen by
        //   A (>= 3 clean layers): distance to the circle through the clean hits,
        //      score 0.75 * dr / SigmaHitDefault + 0.25 * edep score / EdepNormKeV;
        //   B (2 clean layers): kink angle to the line through the two clean hits,
        //      score 0.75 * kink / MaxConsecDeltaPhi + 0.25 * edep score / EdepNormKeV.
        // The result replaces the best combination if its refit has a lower chi2/ndf.
        if (m_doCrowdedLayerHitSelection && !m_useSeedExtend) {
            int nCleanCIDs   = 0;
            int nCrowdedCIDs = 0;
            for (const auto& [cid, cidHIIs] : hitIndicesByCompositeLayer) {
                if (static_cast<int>(cidHIIs.size()) == 1) ++nCleanCIDs;
                else if (static_cast<int>(cidHIIs.size()) > 1) ++nCrowdedCIDs;
            }
            const int totalCIDs = nCleanCIDs + nCrowdedCIDs;

            // clean hits, ordered by layer key
            auto collectCleanHits = [&](std::vector<edm4hep::TrackerHitPlane>& cleanHits,
                                        std::vector<size_t>& cleanHIIs,
                                        std::vector<size_t>& cleanGIdxs) {
                for (const auto& [cid, cidHIIs] : hitIndicesByCompositeLayer) {
                    if (static_cast<int>(cidHIIs.size()) != 1) continue;
                    size_t hii = cidHIIs[0];
                    if (globalUsedHits[allHitInfo[hii].index]) continue;
                    cleanHIIs.push_back(hii);
                    cleanGIdxs.push_back(allHitInfo[hii].index);
                    cleanHits.push_back(allHitInfo[hii].hit);
                }
            };

            // refit the guided choice; keep it if it is different and better
            auto applyGuidedResult = [&](std::vector<edm4hep::TrackerHitPlane>& guidedHits,
                                         std::vector<size_t>& guidedHIIs,
                                         std::vector<size_t>& guidedGIdxs,
                                         const char* branchName) {
                double gx0 = 0, gy0 = 0, gr = 0, gchi2 = 0;
                Eigen::Matrix3d gCov;
                std::vector<double> gRes;
                bool guidedOK = false;
                try { guidedOK = fitCircleNHits(guidedHits, gx0, gy0, gr, gchi2, gCov, gRes); }
                catch (...) {}

                if (!guidedOK) {
                    debug() << branchName << " fit failed — keeping general combo" << endmsg;
                    return;
                }
                if (gchi2 > m_nHitMaxChi2NDF.value()) {
                    debug() << branchName << " fit chi2/ndf=" << gchi2
                            << " > threshold — keeping general combo" << endmsg;
                    return;
                }
                bool differentChoice = false;
                for (size_t gHII : guidedHIIs) {
                    bool inGeneral = false;
                    for (size_t bHII : bestCombo.hiiVec)
                        if (bHII == gHII) { inGeneral = true; break; }
                    if (!inGeneral) { differentChoice = true; break; }
                }
                if (differentChoice && gchi2 < bestCombo.chi2ndf) {
                    info() << branchName << ": "
                           << (differentChoice ? "OVERRIDING" : "CONFIRMING")
                           << " (" << nCrowdedCIDs << " crowded CID(s))"
                           << " | guided chi2/ndf=" << gchi2
                           << " | general chi2/ndf=" << bestCombo.chi2ndf << endmsg;
                    bestCombo.hiiVec    = guidedHIIs;
                    bestCombo.gIdxVec   = guidedGIdxs;
                    bestCombo.chi2ndf   = gchi2;
                    bestCombo.x0        = gx0;
                    bestCombo.y0        = gy0;
                    bestCombo.radius    = gr;
                    bestCombo.fitCov    = gCov;
                    bestCombo.residuals = gRes;
                    bestCombo.nHits     = static_cast<int>(guidedHIIs.size());
                } else {
                    debug() << branchName << ": confirms general combo"
                            << " (all same hits, guided chi2/ndf=" << gchi2 << ")" << endmsg;
                }
            };

            // ── A: circle through the clean hits ────────────────────────────────
            if (nCleanCIDs >= 3 && nCrowdedCIDs >= 1) {
                debug() << "Guided crowded-layer Branch A: " << nCleanCIDs
                        << " clean CIDs + " << nCrowdedCIDs << " crowded CID(s)" << endmsg;

                std::vector<edm4hep::TrackerHitPlane> cleanHits;
                std::vector<size_t> cleanHIIs, cleanGIdxs;
                collectCleanHits(cleanHits, cleanHIIs, cleanGIdxs);

                if (static_cast<int>(cleanHits.size()) >= 3) {
                    double cx0 = 0, cy0 = 0, cr = 0, cchi2 = 0;
                    Eigen::Matrix3d cCov;
                    std::vector<double> cRes;
                    bool cleanOK = false;
                    try { cleanOK = fitCircleNHits(cleanHits, cx0, cy0, cr, cchi2, cCov, cRes); }
                    catch (...) {}

                    if (cleanOK && cr > 0) {
                        const double drNorm   = std::max(m_sigmaHitDefault.value(), 1e-9); // cm
                        const double edepNorm = std::max(m_edepNormKeV.value(),     1e-9); // keV

                        std::vector<edm4hep::TrackerHitPlane> guidedHits = cleanHits;
                        std::vector<size_t> guidedHIIs  = cleanHIIs;
                        std::vector<size_t> guidedGIdxs = cleanGIdxs;
                        bool allCrowdedResolved = true;

                        for (const auto& [cid, cidHIIs] : hitIndicesByCompositeLayer) {
                            if (static_cast<int>(cidHIIs.size()) <= 1) continue;
                            size_t bestHII   = SIZE_MAX;
                            double bestScore = std::numeric_limits<double>::max();
                            for (size_t hii : cidHIIs) {
                                if (globalUsedHits[allHitInfo[hii].index]) continue;
                                const auto& hp = allHitInfo[hii].hit.getPosition(); // mm
                                double px = hp.x * 0.1, py = hp.y * 0.1;           // → cm
                                double dr      = std::abs(
                                    std::sqrt((px-cx0)*(px-cx0) + (py-cy0)*(py-cy0)) - cr);
                                double edepDev = hitEdepScore.count(hii)
                                                 ? hitEdepScore.at(hii) : 0.0;
                                double score   = 0.75 * (dr / drNorm) + 0.25 * (edepDev / edepNorm);
                                debug() << "    BranchA CID=" << cid << " hii=" << hii
                                        << " dr=" << dr << " cm edepDev=" << edepDev
                                        << " keV score=" << score << endmsg;
                                if (score < bestScore) { bestScore = score; bestHII = hii; }
                            }
                            if (bestHII == SIZE_MAX) { allCrowdedResolved = false; break; }
                            guidedHIIs.push_back(bestHII);
                            guidedGIdxs.push_back(allHitInfo[bestHII].index);
                            guidedHits.push_back(allHitInfo[bestHII].hit);
                        }

                        if (allCrowdedResolved)
                            applyGuidedResult(guidedHits, guidedHIIs, guidedGIdxs,
                                              "Guided crowded-layer BranchA");
                    }
                }
            }
            // ── B: line through the two clean hits ──────────────────────────────
            else if (nCleanCIDs == 2 && nCrowdedCIDs >= 1 && totalCIDs >= 3) {
                debug() << "Guided crowded-layer Branch B: 2 clean CIDs + "
                        << nCrowdedCIDs << " crowded CID(s)" << endmsg;

                std::vector<edm4hep::TrackerHitPlane> cleanHits;
                std::vector<size_t> cleanHIIs, cleanGIdxs;
                collectCleanHits(cleanHits, cleanHIIs, cleanGIdxs);

                if (static_cast<int>(cleanHits.size()) == 2) {
                    // direction from the inner (A) to the outer (B) clean hit, cm
                    const auto& posA_mm = cleanHits[0].getPosition();
                    const auto& posB_mm = cleanHits[1].getPosition();
                    double ax = posA_mm.x * 0.1, ay = posA_mm.y * 0.1, az = posA_mm.z * 0.1;
                    double bx = posB_mm.x * 0.1, by = posB_mm.y * 0.1, bz = posB_mm.z * 0.1;
                    double dabx = bx - ax, daby = by - ay, dabz = bz - az;
                    double lenAB = std::sqrt(dabx*dabx + daby*daby + dabz*dabz);

                    const double sigmaKappa = std::max(m_maxConsecDeltaPhi.value(), 1e-4); // rad
                    const double edepNorm   = std::max(m_edepNormKeV.value(),       1e-9); // keV

                    std::vector<edm4hep::TrackerHitPlane> guidedHits = cleanHits;
                    std::vector<size_t> guidedHIIs  = cleanHIIs;
                    std::vector<size_t> guidedGIdxs = cleanGIdxs;
                    bool allCrowdedResolved = (lenAB > 1e-9); // degenerate if A==B

                    if (allCrowdedResolved) {
                        for (const auto& [cid, cidHIIs] : hitIndicesByCompositeLayer) {
                            if (static_cast<int>(cidHIIs.size()) <= 1) continue;
                            size_t bestHII   = SIZE_MAX;
                            double bestScore = std::numeric_limits<double>::max();
                            for (size_t hii : cidHIIs) {
                                if (globalUsedHits[allHitInfo[hii].index]) continue;
                                const auto& hp = allHitInfo[hii].hit.getPosition(); // mm
                                double cx = hp.x * 0.1, cy = hp.y * 0.1, cz = hp.z * 0.1;
                                double dacx = cx - ax, dacy = cy - ay, dacz = cz - az;
                                double lenAC = std::sqrt(dacx*dacx + dacy*dacy + dacz*dacz);
                                // angle at A between AB and AC
                                double kappa = 0.0;
                                if (lenAC > 1e-9) {
                                    double cosK = (dabx*dacx + daby*dacy + dabz*dacz)
                                                  / (lenAB * lenAC);
                                    cosK  = std::max(-1.0, std::min(1.0, cosK));
                                    kappa = std::acos(cosK); // rad
                                }
                                double edepDev = hitEdepScore.count(hii)
                                                 ? hitEdepScore.at(hii) : 0.0;
                                double score   = 0.75 * (kappa / sigmaKappa)
                                               + 0.25 * (edepDev / edepNorm);
                                debug() << "    BranchB CID=" << cid << " hii=" << hii
                                        << " kappa=" << kappa << " rad edepDev=" << edepDev
                                        << " keV score=" << score << endmsg;
                                if (score < bestScore) { bestScore = score; bestHII = hii; }
                            }
                            if (bestHII == SIZE_MAX) { allCrowdedResolved = false; break; }
                            guidedHIIs.push_back(bestHII);
                            guidedGIdxs.push_back(allHitInfo[bestHII].index);
                            guidedHits.push_back(allHitInfo[bestHII].hit);
                        }

                        if (allCrowdedResolved)
                            applyGuidedResult(guidedHits, guidedHIIs, guidedGIdxs,
                                              "Guided crowded-layer BranchB");
                    }
                }
            }
        }

        // ── Neighbour-track cuts ──────────────────────────────────────────────────
        // A combination with a hit within NeighbourTrackMaxDist (cm) of an accepted track
        // may be a delta-ray shadow. It is kept only with >= NeighbourTrackMinHits hits
        // and chi2/ndf <= NeighbourTrackMaxChi2NDF; otherwise its hits near the accepted
        // track are marked as used and the search goes on.
        // Not applied in BestFirst mode, whose tracks already own their hits.
        if (!m_seBestFirst && m_neighbourTrackMaxDist.value() > 0.0 && !acceptedTrackHitPositions.empty()) {
            const double distThreshMm  = m_neighbourTrackMaxDist.value() * 10.0; // cm → mm
            const double distThreshSq  = distThreshMm * distThreshMm;
            bool isNeighbour = false;
            for (size_t hii : bestCombo.hiiVec) {
                const auto& hp = allHitInfo[hii].hit.getPosition();  // mm
                for (const auto& ap : acceptedTrackHitPositions) {
                    double dx = hp.x - ap[0], dy = hp.y - ap[1], dz = hp.z - ap[2];
                    if (dx*dx + dy*dy + dz*dz < distThreshSq) {
                        isNeighbour = true;
                        break;
                    }
                }
                if (isNeighbour) break;
            }

            if (isNeighbour) {
                const int  nHits    = static_cast<int>(bestCombo.hiiVec.size());
                const bool passHits = (nHits    >= m_neighbourTrackMinHits.value());
                const bool passChi2 = (bestCombo.chi2ndf <= m_neighbourTrackMaxChi2NDF.value());
                if (!passHits || !passChi2) {
                    info() << "Neighbour-gate: iter " << trackNumber
                           << " REJECTED shadow combo"
                           << " hits=" << nHits
                           << " chi2/ndf=" << bestCombo.chi2ndf
                           << "  (required >=" << m_neighbourTrackMinHits.value()
                           << " hits AND <=" << m_neighbourTrackMaxChi2NDF.value()
                           << " chi2/ndf)" << endmsg;
                    ++m_statNeighbourRejected;
                    // only the hits near the accepted track are used up
                    int nConsumed = 0;
                    for (size_t hii : bestCombo.hiiVec) {
                        const auto& hp = allHitInfo[hii].hit.getPosition();  // mm
                        bool hitIsNeighbour = false;
                        for (const auto& ap : acceptedTrackHitPositions) {
                            double dx = hp.x - ap[0], dy = hp.y - ap[1], dz = hp.z - ap[2];
                            if (dx*dx + dy*dy + dz*dz < distThreshSq) {
                                hitIsNeighbour = true;
                                break;
                            }
                        }
                        if (hitIsNeighbour) {
                            globalUsedHits[allHitInfo[hii].index] = true;
                            ++nConsumed;
                        }
                    }
                    info() << "  Neighbour-gate: consumed " << nConsumed << " of "
                           << bestCombo.hiiVec.size() << " combo hits (within "
                           << m_neighbourTrackMaxDist.value() * 10.0 << " mm of accepted track)"
                           << endmsg;
                    continue;
                } else {
                    info() << "Neighbour-gate: iter " << trackNumber
                           << " passes stricter cuts"
                           << " (hits=" << nHits
                           << ", chi2/ndf=" << bestCombo.chi2ndf << ")" << endmsg;
                }
            }
        }

        {
            info() << "  Best combo: " << bestCombo.nHits << " hits, chi2/ndf="
                   << bestCombo.chi2ndf << endmsg;
            for (int bi = 0; bi < bestCombo.nHits; ++bi) {
                size_t gi  = bestCombo.gIdxVec[bi];
                size_t hii = bestCombo.hiiVec[bi];
                const auto& pos = (*hits)[gi].getPosition();
                double res = (bi < static_cast<int>(bestCombo.residuals.size()))
                             ? bestCombo.residuals[bi] : -1.0;
                info() << "    hit[" << bi << "] gIdx=" << gi
                       << " compositeID=" << allHitInfo[hii].compositeID
                       << " pos=(" << std::fixed << std::setprecision(1)
                       << pos[0]/10.0 << "," << pos[1]/10.0 << "," << pos[2]/10.0
                       << ") cm  residual=" << std::setprecision(4) << res << " cm"
                       << endmsg;
            }
        }

        // ── Outlier removal: drop the worst hit while its residual is above
        //    OutlierSigma * SigmaHitDefault, and refit (not in BestFirst mode) ─────────
        {
            const double outlierCutCm = m_outlierSigma.value() * m_sigmaHitDefault.value();
            for (int iter = 0; !m_seBestFirst && iter < m_maxOutlierIterations.value(); ++iter) {
                if (static_cast<int>(bestCombo.gIdxVec.size()) <= minK) break;

                size_t worstIdx = 0;
                double worstRes = 0.0;
                for (size_t ri = 0; ri < bestCombo.residuals.size(); ++ri) {
                    if (bestCombo.residuals[ri] > worstRes) {
                        worstRes = bestCombo.residuals[ri];
                        worstIdx = ri;
                    }
                }
                if (worstRes <= outlierCutCm) break;

                {
                    size_t wGi = bestCombo.gIdxVec[worstIdx];
                    size_t wHii = bestCombo.hiiVec[worstIdx];
                    const auto& wPos = (*hits)[wGi].getPosition();
                    info() << "  Outlier iter " << iter+1 << ": removing hit[" << worstIdx
                           << "] gIdx=" << wGi
                           << " compositeID=" << allHitInfo[wHii].compositeID
                           << " pos=(" << std::fixed << std::setprecision(1)
                           << wPos[0]/10.0 << "," << wPos[1]/10.0 << "," << wPos[2]/10.0
                           << ") cm  residual=" << std::setprecision(4) << worstRes
                           << " cm (cut=" << outlierCutCm << " cm)" << endmsg;
                }
                m_statOutlierHitsRemoved++;

                bestCombo.gIdxVec.erase(bestCombo.gIdxVec.begin() + worstIdx);
                bestCombo.hiiVec.erase(bestCombo.hiiVec.begin()   + worstIdx);
                bestCombo.residuals.erase(bestCombo.residuals.begin() + worstIdx);
                bestCombo.nHits = static_cast<int>(bestCombo.gIdxVec.size());

                if (bestCombo.nHits < minK) break;

                std::vector<edm4hep::TrackerHitPlane> refitHits;
                for (size_t gi : bestCombo.gIdxVec) refitHits.push_back((*hits)[gi]);
                double rx0, ry0, rr, rchi2;
                Eigen::Matrix3d rCov;
                std::vector<double> rRes;
                bool refitOK = false;
                try { refitOK = fitCircleNHits(refitHits, rx0, ry0, rr, rchi2, rCov, rRes); }
                catch (...) {}
                if (refitOK) {
                    bestCombo.x0        = rx0; bestCombo.y0      = ry0;
                    bestCombo.radius    = rr;  bestCombo.chi2ndf = rchi2;
                    bestCombo.fitCov    = rCov; bestCombo.residuals = rRes;
                }
            }
        }

        for (size_t gi : bestCombo.gIdxVec) usedHits[gi] = true;

        // ── Internal helix parameters of the chosen hits ──────────────────────────
        // Circle of bestCombo + z-vs-arc-length fit. These internal parameters seed
        // GenFit and the inner propagation; their sign convention ("charge" +1 turns
        // counter-clockwise) holds in the muon system's return field (B < 0). The
        // stored track states are built separately (makeEdm4hepStateFromCircle).

        const int  Nh   = bestCombo.nHits;
        const double cx = bestCombo.x0;     // cm
        const double cy = bestCombo.y0;     // cm
        const double cr = bestCombo.radius; // cm

        // charge from the bending of the hits and the field at their centroid
        {   double sumXf=0, sumYf=0, sumZf=0;
            for (size_t gi : bestCombo.gIdxVec) {
                const auto& p = (*hits)[gi].getPosition();
                sumXf += p[0]/10.0; sumYf += p[1]/10.0; sumZf += p[2]/10.0;
            }
            dd4hep::Position avgPosN(sumXf/Nh, sumYf/Nh, sumZf/Nh);
            double bFieldN = m_field.magneticField(avgPosN).z() / dd4hep::tesla;
            std::vector<edm4hep::TrackerHitPlane> orderedN;
            orderedN.reserve(Nh);
            for (int i = 0; i < Nh; ++i) orderedN.push_back((*hits)[bestCombo.gIdxVec[i]]);

            const k4rt::BendSense bsN = k4rt::bendSense(orderedN);
            const double sinBendSumN = bsN.sinBendMean;
            const double avgChordN   = bsN.avgChord;
            bool chargeReliableN = (avgChordN > 1e-6) &&
                                   (std::abs(sinBendSumN) >= 0.5/avgChordN);
            double chargeN = k4rt::chargeFromBend(sinBendSumN, bFieldN);

            const k4rt::ZFit zfN = k4rt::fitZvsArc(orderedN, cx, cy, cr, chargeN);
            double z0_cm = zfN.z0_cm, tanLambda = zfN.tanLambda;

            const k4rt::Perigee2D pgN = k4rt::perigeeFromCircle(cx, cy, cr, chargeN);
            double d0_cm  = pgN.d0_cm;
            double phi    = pgN.phi_rad;
            double omega  = pgN.omega_mm;   // 1/mm
            double pTN    = k4rt::ptFromRadius(bFieldN, cr); // GeV/c

            info() << "N-hit circle fit: " << Nh << " hits  pT=" << pTN << " GeV/c"
                   << "  charge=" << chargeN
                   << (chargeReliableN ? "" : " [UNRELIABLE]")
                   << "  chi2/ndf=" << bestCombo.chi2ndf << endmsg;

            bestCombo.charge    = chargeN;
            bestCombo.d0_cm     = d0_cm;
            bestCombo.phi_rad   = phi;
            bestCombo.omega_mm  = omega;
            bestCombo.z0_cm     = z0_cm;
            bestCombo.tanLambda = tanLambda;
            bestCombo.pT        = pTN;

            auto nhSeedTrack = candidateTracks.create();  // empty track (no states)
            tripletCandidates = 1; validTriplets = 1;
            m_statTotalTripletCombos++;
            m_statValidTriplets++;

            std::vector<int> usedComps;
            for (size_t hii : bestCombo.hiiVec) usedComps.push_back(allHitInfo[hii].compositeID);
            trackCandidates.emplace_back(nhSeedTrack, pTN, bestCombo.gIdxVec, usedComps);
        }

        // one candidate per iteration: the best combination
        if (trackCandidates.empty()) {
            info() << "No track candidates found in iteration " << trackNumber << " - stopping track search" << endmsg;
            break;
        }

        std::sort(trackCandidates.begin(), trackCandidates.end(),
                  [](const TrackCandidate& a, const TrackCandidate& b) {
                      return a.pT > b.pT;
                  });

        const auto& candidate = trackCandidates.front();

        for (size_t idx : candidate.hitIndices) {
            globalUsedHits[idx] = true;
        }

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
               << " GeV/c  nHits=" << candidate.hitIndices.size() << endmsg;

        // GenFit seed (internal convention); replaced by the circle-fit state below
        edm4hep::TrackState seedState = createTrackState(
            bestCombo.d0_cm * 10.0,   // cm → mm
            bestCombo.phi_rad,
            bestCombo.omega_mm,        // already in 1/mm
            bestCombo.z0_cm * 10.0,   // cm → mm
            bestCombo.tanLambda,
            edm4hep::TrackState::AtLastHit);
        {
            const auto& firstHitPos = (*hits)[bestCombo.gIdxVec[0]].getPosition();
            seedState.referencePoint = edm4hep::Vector3f(firstHitPos[0], firstHitPos[1], firstHitPos[2]);
        }
        debug() << "Built initial seed state from N-hit circle fit for track " << trackNumber << endmsg;

        // the track's hits, by index (cellIDs need not be unique)
        std::vector<edm4hep::TrackerHitPlane> trackHits;
        trackHits.reserve(3);

        for (size_t idx : candidate.hitIndices) {
            if (idx < hits->size()) {
                trackHits.push_back((*hits)[idx]);
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

        info() << "N-hit track " << trackNumber << " has " << trackHits.size()
               << " hits after combinatorial selection and outlier rejection" << endmsg;

        // for the neighbour-track cuts of the next iterations
        if (m_neighbourTrackMaxDist.value() > 0.0) {
            for (const auto& th : trackHits) {
                const auto& p = th.getPosition();  // mm
                acceptedTrackHitPositions.push_back({p.x, p.y, p.z});
            }
        }

        auto finalTrack = finalTracks.create();

        debug() << "Created new final track object for track " << trackNumber << endmsg;

        // circle-fit state in the internal convention: seeds GenFit and the inner
        // propagation, is not stored
        edm4hep::TrackState internalLastHitState;
        bool haveInternalLastHit = false;

        // circle of bestCombo (after outlier removal)
        if (trackHits.size() >= 3) {
            {
                double x0     = bestCombo.x0;
                double y0     = bestCombo.y0;
                double radius = bestCombo.radius;
                double chi2   = bestCombo.chi2ndf;
                Eigen::Matrix3d fitCov3x3 = bestCombo.fitCov;
                bool fit4Hits = (radius > 0);

                if (fit4Hits) {
                info() << "N-hit circle fit for track " << trackNumber
                       << " (" << trackHits.size() << " hits): chi2/ndf=" << chi2 << endmsg;
                info() << "  Center: (" << x0 << ", " << y0 << ") cm" << endmsg;
                info() << "  Radius: " << radius << " cm" << endmsg;

                if (false) {
                } else {
                    // field at the centroid of the hits (DD4hep positions in cm)
                    double sumX = 0, sumY = 0, sumZ = 0;
                    for (const auto& h : trackHits) {
                        const auto& pos = h.getPosition();
                        sumX += pos[0] / 10.0;  // mm → cm
                        sumY += pos[1] / 10.0;
                        sumZ += pos[2] / 10.0;
                    }
                    dd4hep::Position avgPos(sumX / trackHits.size(),
                                           sumY / trackHits.size(),
                                           sumZ / trackHits.size());
                    double bField4 = m_field.magneticField(avgPos).z() / dd4hep::tesla;
                    double pT = k4rt::ptFromRadius(bField4, radius);  // GeV/c
                    double bSign4 = (bField4 >= 0.0) ? 1.0 : -1.0;

                    // charge from the mean bending angle over consecutive hit triplets
                    const k4rt::BendSense bs4 = k4rt::bendSense(trackHits);
                    double sinBendSum4 = bs4.sinBendMean;
                    double avgChord4   = bs4.avgChord;

                    double sinThreshold4 = (avgChord4 > 1e-6) ? 0.5 / avgChord4 : 1.0;
                    bool chargeReliable4 = (std::abs(sinBendSum4) >= sinThreshold4);
                    double charge4 = bSign4 * (sinBendSum4 < 0.0 ? 1.0 : -1.0);

                    if (!chargeReliable4) {
                        debug() << "N-hit charge UNRELIABLE (near-straight track):"
                                << "  sinBendSum4=" << sinBendSum4
                                << "  threshold=" << sinThreshold4
                                << "  avgChord=" << avgChord4 << " mm" << endmsg;
                    }

                    double chargeSeed = (seedState.omega >= 0.0) ? 1.0 : -1.0;
                    if (charge4 != chargeSeed) {
                        debug() << "N-hit charge DIFFERS from 3-hit seed:"
                                << "  charge4=" << charge4
                                << (chargeReliable4 ? "" : " [UNRELIABLE]")
                                << "  chargeSeed=" << chargeSeed
                                << "  sinBendSum4=" << sinBendSum4 << endmsg;
                    }

                    // DEBUG: MC particle of the first hit (a muon if there is one)
                    {
                        auto mcIt = recoToSimMap.find(hitKey(trackHits[0]));
                        if (mcIt != recoToSimMap.end() && !mcIt->second.empty()) {
                            edm4hep::SimTrackerHit chosenSim = mcIt->second.front();
                            for (const auto& s : mcIt->second) {
                                auto p = s.getParticle();
                                if (p.isAvailable() && std::abs(p.getPDG()) == 13) {
                                    chosenSim = s;
                                    break;
                                }
                            }
                            auto mcPart = chosenSim.getParticle();
                            if (mcPart.isAvailable()) {
                                const auto& mom = mcPart.getMomentum();
                                double mcPt = std::sqrt(mom.x*mom.x + mom.y*mom.y);
                                debug() << "[MC-TRUTH 4-hit] charge=" << mcPart.getCharge()
                                        << "  PDG=" << mcPart.getPDG()
                                        << "  px=" << mom.x << " GeV/c"
                                        << "  py=" << mom.y << " GeV/c"
                                        << "  pz=" << mom.z << " GeV/c"
                                        << "  pT=" << mcPt << " GeV/c"
                                        << "  reco_chargeN=" << charge4
                                        << (chargeReliable4 ? "" : " [UNRELIABLE]")
                                        << endmsg;
                            } else {
                                debug() << "[MC-TRUTH 4-hit] MCParticle not available for hit[0]" << endmsg;
                            }
                        } else {
                            debug() << "[MC-TRUTH 4-hit] No sim link found for hit[0]" << endmsg;
                        }
                    }

                    const k4rt::Perigee2D pg4 = k4rt::perigeeFromCircle(x0, y0, radius, charge4);
                    double d0    = pg4.d0_cm;
                    double phi   = pg4.phi_rad;
                    double omega = pg4.omega_mm;  // 1/mm

                    debug() << "N-hit AtLastHit: charge=" << charge4
                            << (chargeReliable4 ? "" : " [UNRELIABLE]")
                            << "  pT=" << pT << " GeV/c"
                            << "  sinBendSum4=" << sinBendSum4
                            << "  omega=" << std::scientific << std::setprecision(4)
                            << omega << std::defaultfloat << " 1/mm" << endmsg;

                    edm4hep::TrackState circleFitState = createTrackState(
                        d0*10.0, phi, omega, seedState.Z0, seedState.tanLambda,
                        edm4hep::TrackState::AtLastHit);

                    const auto& firstHitPos = trackHits.front().getPosition();
                    circleFitState.referencePoint = edm4hep::Vector3f(
                        firstHitPos[0], firstHitPos[1], firstHitPos[2]);

                // covariance of (x0, y0, R) -> d0, phi, omega (mm, rad, 1/mm)
                {
                    double cx = x0;
                    double cy = y0;

                    double r2 = cx*cx + cy*cy;
                    double dist0 = std::sqrt(r2);

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

                    double var_d0_cm2    = grad_d0.dot(fitCov3x3 * grad_d0);
                    double var_phi_rad2  = grad_phi.dot(fitCov3x3 * grad_phi);
                    double var_omega_mm2 = grad_omega.dot(fitCov3x3 * grad_omega);

                    double var_d0_mm2 = var_d0_cm2 * 100.0;

                    // z0, tanLambda variances from the seed
                    double var_z0 = seedState.getCovMatrix(
                        edm4hep::TrackParams::z0,
                        edm4hep::TrackParams::z0);

                    double var_tanL = seedState.getCovMatrix(
                        edm4hep::TrackParams::tanLambda,
                        edm4hep::TrackParams::tanLambda);

                    if (var_z0   < 1e-30) var_z0   = 1.0;
                    if (var_tanL < 1e-30) var_tanL = 0.01;

                    // floors
                    var_d0_mm2   = std::max(var_d0_mm2,   0.01);   // 0.1 mm
                    var_phi_rad2 = std::max(var_phi_rad2, 0.01);

                    double omega_rel_floor = std::pow(0.1 * std::abs(omega), 2);

                    var_omega_mm2 = std::max({var_omega_mm2, omega_rel_floor});

                    for (int ii = 0; ii < 6; ++ii)
                        for (int jj = 0; jj <= ii; ++jj)
                            circleFitState.setCovMatrix(
                                0.0,
                                static_cast<edm4hep::TrackParams>(ii),
                                static_cast<edm4hep::TrackParams>(jj)
                            );

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

                    info() << "Circle-fit covariance:"
                        << "  σ_d0="   << std::sqrt(var_d0_mm2)    << " mm"
                        << "  σ_phi="  << std::sqrt(var_phi_rad2)  << " rad"
                        << "  σ_ω="    << std::sqrt(var_omega_mm2) << " 1/mm"
                        << "  σ_z0="   << std::sqrt(var_z0)        << " mm"
                        << "  σ_tanL=" << std::sqrt(var_tanL)
                        << endmsg;
                }
                    internalLastHitState = circleFitState;
                    haveInternalLastHit = true;
                    {
                        // Stored states: the fitted helix at the first and at the last hit
                        // and, with IPStateMode = Helix, at its closest approach to the origin.
                        struct StoredState { StatePoint where; int location; bool wanted; };
                        const StoredState stored[] = {
                            {StatePoint::FirstHit, edm4hep::TrackState::AtFirstHit, true},
                            {StatePoint::LastHit,  edm4hep::TrackState::AtLastHit,  m_storeLastHitState.value()},
                            {StatePoint::Origin,   edm4hep::TrackState::AtIP,       m_ipStateMode.value() == "Helix"},
                        };
                        for (const auto& entry : stored) {
                            if (!entry.wanted) continue;
                            edm4hep::TrackState state = circleFitState;
                            if (!makeEdm4hepStateFromCircle(x0, y0, radius, fitCov3x3, trackHits, entry.where, state)) {
                                warning() << "Could not build track state " << entry.location << " of track "
                                          << trackNumber << " (degenerate hits)" << endmsg;
                                continue;
                            }
                            state.location = entry.location;
                            finalTrack.addToTrackStates(state);
                        }
                    }
                    // the circle fit has one radial residual per hit and 3 parameters
                    const int fitNdf = std::max(1, static_cast<int>(trackHits.size()) - 3);
                    finalTrack.setChi2(static_cast<float>(chi2 * fitNdf));   // store raw chi2
                    finalTrack.setNdf(fitNdf);

                    seedState = circleFitState;
                    debug() << "Circle-fit state used as GenFit seed" << endmsg;

                    info() << "Added circle-fit track states for track " << trackNumber << ": pT=" << pT
                           << " GeV/c, d0=" << d0*10.0 << " mm" << endmsg;
                }
            } else {
                warning() << "Circle fit failed for track " << trackNumber << endmsg;
            }
            }
        }

        // GenFit fit (UseGenFit); it also copies the hits to the output
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

            // circle-fit ndf (n_hits - 3)
            finalTrack.setNdf(std::max(1, static_cast<int>(trackHits.size()) - 3));
            // BestFirst: the stored chi2 is the kink chi2 (with multiple scattering);
            // loopers get track type 1
            if (bestCombo.kinkNdf > 0) {
                finalTrack.setChi2(static_cast<float>(bestCombo.kinkChi2));
                finalTrack.setNdf(bestCombo.kinkNdf);
            }
            if (bestCombo.looper) finalTrack.setType(1);

            // copy the hits to the output collection
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
                outHit.setU(hit.getU());
                outHit.setV(hit.getV());
                propagateLink(outHit, hit);
                finalTrack.addToTrackerHits(outHit);
            }
        }

        // ======= RK4 propagation inward through the solenoid (DoInnerPropagation) =======
        // For tracks measured outside the solenoid. Starts at the first hit from the
        // GenFit state if there is one, else from the internal circle-fit state.
        // The result is only printed.
        if (m_doInnerPropagation) {
        edm4hep::TrackState bestState;
        bool foundStateForProp = false;
        for (int j = 0; j < finalTrack.trackStates_size(); ++j) {
            auto st = finalTrack.getTrackStates(j);
            if (st.location == edm4hep::TrackState::AtOther) {
                bestState = st;
                foundStateForProp = true;
                debug() << "Using AtOther state (GenFit fitted) for inner propagation" << endmsg;
                break;
            }
        }
        if (!foundStateForProp && haveInternalLastHit) {
            bestState = internalLastHitState;
            foundStateForProp = true;
            debug() << "Using the circle-fit state for inner propagation" << endmsg;
        }

        if (foundStateForProp) {
            double d0 = bestState.D0 / 10.0;           // mm → cm
            double phi = bestState.phi;
            double omega = bestState.omega * 10.0;     // 1/mm → 1/cm
            double z0 = bestState.Z0 / 10.0;           // mm → cm
            double tanLambda = bestState.tanLambda;

            if (std::abs(omega) > 1e-10) {
                double radius = 1.0 / std::abs(omega);

                // start at the first hit
                Eigen::Vector3d pos;
                if (finalTrack.trackerHits_size() > 0) {
                    auto firstHit = finalTrack.getTrackerHits(0);
                    pos = Eigen::Vector3d(
                        firstHit.getPosition()[0] / 10.0,  // mm to cm
                        firstHit.getPosition()[1] / 10.0,
                        firstHit.getPosition()[2] / 10.0
                    );
                }

                dd4hep::Position fieldPos(pos.x(), pos.y(), pos.z());
                double bField = m_field.magneticField(fieldPos).z() / dd4hep::tesla;

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

                double chargeSign = (omega > 0.0) ? 1.0 : -1.0;

                auto inner = propagateToInner(pos, mom, m_innerPropTargetRadius, chargeSign);

                if (inner.success) {
                    m_statInnerPropSuccess++;
                    double rxy = std::sqrt(
                        inner.finalPosition.x() * inner.finalPosition.x() +
                        inner.finalPosition.y() * inner.finalPosition.y()
                    );
                    double pMagFinal = inner.finalMomentum.norm();
                    double pRatio = pMagFinal / mom.norm();

                    info() << "  ✓ Inner propagation SUCCESS:" << endmsg;
                    info() << "    Final R=" << rxy << " cm, |p|=" << pMagFinal
                        << " GeV/c (ratio=" << pRatio << ")" << endmsg;
                    info() << "    Arc length: " << inner.arcLength << " cm, steps: "
                        << inner.numSteps << endmsg;
                } else {
                    debug() << "  ✗ Inner propagation failed: " << inner.message << endmsg;
                }
            }
        }
        } // end if (m_doInnerPropagation)

        // ===== IPStateMode = SolenoidCrossing: AtIP from the helix matched at the solenoid =====
        // Only for tracks with all hits outside the solenoid.
        if (m_ipStateMode.value() == "SolenoidCrossing" && allHitsOutsideSolenoid(trackHits)
            && foundCombination && std::abs(bestCombo.omega_mm) > 1e-10
            && finalTrack.trackerHits_size() > 0)
        {
            // outer circle (cm)
            double R_outer_anal = bestCombo.radius;  // cm
            double cx_o_anal    = bestCombo.x0;      // cm
            double cy_o_anal    = bestCombo.y0;      // cm

            // the innermost hit anchors the crossing
            auto refHit = finalTrack.getTrackerHits(0);
            double xr = refHit.getPosition()[0] / 10.0;  // mm → cm
            double yr = refHit.getPosition()[1] / 10.0;
            double zr = refHit.getPosition()[2] / 10.0;

            // barrel hit: crossing of the solenoid barrel; endcap hit: of its end face
            bool isEndcapHit = (getTypeID(refHit.getCellID()) != 0);
            debug() << "  [Analytical] Ref hit typeID=" << getTypeID(refHit.getCellID())
                    << "  -> " << (isEndcapHit ? "endcap" : "barrel") << " crossing" << endmsg;

            bool atIPok = analyticalInnerPropagation(
                cx_o_anal, cy_o_anal, R_outer_anal,
                bestCombo.omega_mm, bestCombo.pT, bestCombo.tanLambda,
                xr, yr, zr,
                isEndcapHit,
                finalTrack);

            if (atIPok) {
                m_statAnalyticalPropSuccess++;

                // printout: d0 of the AtIP state vs the true production radius
                try {
                    int bestMCCount = 0;
                    edm4hep::MCParticle bestMCpart;
                    bool hasMCpart = false;
                    for (const auto& th : trackHits) {
                        auto it = recoToSimMap.find(hitKey(th));
                        if (it == recoToSimMap.end()) continue;
                        // Prefer a muon contributor when the digi-hit has several.
                        for (const auto& sim : it->second) {
                            auto mcp = sim.getParticle();
                            if (!mcp.isAvailable()) continue;
                            if (!hasMCpart) {
                                bestMCpart = mcp; bestMCCount = 1; hasMCpart = true;
                            }
                            if (std::abs(mcp.getPDG()) == 13) {
                                bestMCpart = mcp;
                                break;
                            }
                        }
                    }
                    if (hasMCpart) {
                        const auto& vtx = bestMCpart.getVertex();  // mm (EDM4hep)
                        double rxy_truth = std::sqrt(vtx.x*vtx.x + vtx.y*vtx.y);  // mm
                        double d0_reco_mm = 0.0;
                        for (int si = 0; si < finalTrack.trackStates_size(); ++si) {
                            auto ts = finalTrack.getTrackStates(si);
                            if (ts.location == edm4hep::TrackState::AtIP) {
                                d0_reco_mm = ts.D0;
                                break;
                            }
                        }
                        bool suspicious = (std::abs(d0_reco_mm) > rxy_truth + 20.0);
                        info() << "  [TruthRxy] rxy_vtx=" << std::fixed << std::setprecision(2)
                               << rxy_truth << " mm  d0_reco=" << d0_reco_mm << " mm"
                               << "  PDG=" << bestMCpart.getPDG()
                               << (suspicious ? "  ** |d0|>rxy+20mm **" : "") << endmsg;
                    }
                } catch (const std::exception& e) {
                    debug() << "[TruthRxy] Exception: " << e.what() << endmsg;
                }

                // printout: all states vs the matched MC particle
                try {
                    if (m_mcParticles.exist()) {
                        const auto* mcParts = m_mcParticles.get();
                        if (mcParts && !mcParts->empty()) {
                            edm4hep::Track tv = finalTrack;
                            compareStatesWithTruth(tv, *mcParts, trackNumber);
                        }
                    }
                } catch (const std::exception& e) {
                    debug() << "[TruthCmp] Exception: " << e.what() << endmsg;
                }
            } else {
                debug() << "  [Analytical] AtIP propagation failed for track "
                        << trackNumber << endmsg;
            }
        }

        info() << "Successfully created track " << trackNumber
            << " with " << finalTrack.trackerHits_size() << " hits"
            << " and " << finalTrack.trackStates_size() << " track states" << endmsg;

        // ---- per-track run statistics ----
        m_statTracksReconstructed++;
        bool hasLastHit = false, hasOther = false, hasVertex = false, hasIP = false;
        for (int sIdx = 0; sIdx < finalTrack.trackStates_size(); ++sIdx) {
            auto ts = finalTrack.getTrackStates(sIdx);
            switch (ts.location) {
                case edm4hep::TrackState::AtLastHit:  hasLastHit = true; break;
                case edm4hep::TrackState::AtOther:    hasOther   = true; break;
                case edm4hep::TrackState::AtVertex:   hasVertex  = true; break;
                case edm4hep::TrackState::AtIP:       hasIP      = true; break;
                default: break;
            }
        }
        if (hasLastHit) m_statStateAtLastHit++;
        if (hasOther)   m_statStateAtOther++;
        if (hasVertex)  m_statStateAtVertex++;
        if (hasIP)      m_statStateAtIP++;
        if (hasLastHit) m_statFourHitTracks++; else m_statThreeHitTracks++;

        {
            int nFinalHits = static_cast<int>(trackHits.size());
            if      (nFinalHits <= 3) m_statNhit3++;
            else if (nFinalHits == 4) m_statNhit4++;
            else if (nFinalHits == 5) m_statNhit5++;
            else if (nFinalHits == 6) m_statNhit6++;
            else                      m_statNhit7plus++;
        }
        m_statValidTriplets++;

        // DEBUG: AtFirstHit state vs the MC particle with most hits on the track
        if (msgLevel(MSG::DEBUG)) {
            std::unordered_map<int, std::pair<edm4hep::MCParticle, int>> mcCount;
            for (const auto& th : trackHits) {
                auto it = recoToSimMap.find(hitKey(th));
                if (it == recoToSimMap.end()) continue;
                // each MC particle counts once per hit
                std::unordered_set<int> seenInThisHit;
                for (const auto& sim : it->second) {
                    auto mcPart = sim.getParticle();
                    if (!mcPart.isAvailable()) continue;
                    int mcIdx = mcPart.id().index;
                    if (!seenInThisHit.insert(mcIdx).second) continue;
                    mcCount[mcIdx].first  = mcPart;
                    mcCount[mcIdx].second++;
                }
            }
            int bestCount = 0;
            edm4hep::MCParticle bestMC;
            bool hasMC = false;
            for (const auto& [idx, pair] : mcCount) {
                if (pair.second > bestCount) {
                    bestCount = pair.second;
                    bestMC    = pair.first;
                    hasMC     = true;
                }
            }
            if (hasMC) {
                const auto& mom = bestMC.getMomentum();
                double mcPt  = std::sqrt(mom.x*mom.x + mom.y*mom.y);
                double mcP   = std::sqrt(mcPt*mcPt + mom.z*mom.z);
                double mcEta = (mcPt > 0.0) ? std::asinh(mom.z / mcPt) : 0.0;
                double mcPhi = std::atan2(mom.y, mom.x);

                double recoPt = 0.0, recoPhi = 0.0, recoD0mm = 0.0, recoZ0mm = 0.0, recoEta = 0.0;
                bool hasRecoState = false;
                for (int sIdx = 0; sIdx < finalTrack.trackStates_size(); ++sIdx) {
                    auto ts = finalTrack.getTrackStates(sIdx);
                    if (ts.location == edm4hep::TrackState::AtFirstHit) {
                        double omegaMm = ts.omega;           // 1/mm
                        if (std::abs(omegaMm) > 1e-12) {
                            recoPt = getPT(ts);   // field at the first hit
                        }
                        recoPhi      = ts.phi;
                        recoD0mm     = ts.D0;
                        recoZ0mm     = ts.Z0;
                        recoEta      = (std::abs(ts.tanLambda) > 0.0)
                                       ? std::asinh(ts.tanLambda) : 0.0;
                        hasRecoState = true;
                        break;
                    }
                }

                double dPt  = (mcPt  > 0) ? (recoPt  - mcPt)  / mcPt  * 100.0 : 0.0;
                double dPhi = recoPhi - mcPhi;
                while (dPhi >  M_PI) dPhi -= 2*M_PI;
                while (dPhi < -M_PI) dPhi += 2*M_PI;

                debug() << "[RECO vs TRUTH trk " << trackNumber << "]"
                        << "  PDG=" << bestMC.getPDG()
                        << "  q_truth=" << bestMC.getCharge()
                        << "  hits_matched=" << bestCount << "/" << trackHits.size()
                        << endmsg;
                debug() << "  Truth: pT=" << std::fixed << std::setprecision(3) << mcPt
                        << " GeV/c  phi=" << mcPhi << " rad  eta=" << mcEta
                        << "  |p|=" << mcP << " GeV/c" << endmsg;
                debug() << "  Reco : pT=" << recoPt
                        << " GeV/c  phi=" << recoPhi << " rad  eta=" << recoEta
                        << "  d0=" << recoD0mm << " mm  z0=" << recoZ0mm << " mm" << endmsg;
                debug() << "  Delta: dpT/pT=" << std::setprecision(1) << dPt << "%"
                        << "  dphi=" << std::setprecision(4) << dPhi << " rad"
                        << "  deta=" << recoEta - mcEta << endmsg;
            } else {
                debug() << "[RECO vs TRUTH trk " << trackNumber << "] no MC link found for track hits" << endmsg;
            }
        }

        // run statistics: sign of omega of the AtFirstHit state
        {
            double omegaLastHit = 0.0;  bool hasLastHitOmega = false;
            for (int sIdx = 0; sIdx < finalTrack.trackStates_size(); ++sIdx) {
                auto ts = finalTrack.getTrackStates(sIdx);
                if (ts.location == edm4hep::TrackState::AtFirstHit) {
                    omegaLastHit    = ts.omega;
                    hasLastHitOmega = true;
                    break;
                }
            }
            if (hasLastHitOmega) {
                // omega > 0: clockwise (positive charge in a field along +z)
                if (omegaLastHit > 0.0) m_statPositiveCharge++;
                else                     m_statNegativeCharge++;
                debug() << "Charge classification: omega=" << omegaLastHit
                        << " 1/mm (AtFirstHit)  →  " << (omegaLastHit > 0 ? "positive" : "negative") << endmsg;
            }
        }
    } // end of the track loop

    info() << "Track reconstruction complete: created " << finalTracks.size()
           << " final tracks across " << (trackNumber-1) << " iterations" << endmsg;
}

// ─────────────────────────────────────────────────────────────────────────────
// IPStateMode = SolenoidCrossing: AtIP for a track measured outside the solenoid.
//
// The solenoid boundary is a cylinder at R = R_sol closed by two disks at
// |z| = Z_sol. Going backward from the innermost hit, the outer circle is followed
// to the boundary crossing (barrel or end disk); there an inner helix with the
// inner field and the same pT and tanLambda is attached, and its perigee is stored.
// ─────────────────────────────────────────────────────────────────────────────
bool CombinatorialTrackFinder::analyticalInnerPropagation(
    double cx_o, double cy_o, double R_outer,
    double omega_outer_mm,
    double pT,
    double tanLambda,
    double x_ref, double y_ref, double z_ref,
    bool   isEndcap,
    edm4hep::MutableTrack& finalTrack) const {

    auto normAngle = [](double a) -> double {   // → (−π, π]
        while (a >  M_PI) a -= 2.0*M_PI;
        while (a < -M_PI) a += 2.0*M_PI;
        return a;
    };
    auto normPos = [](double a) -> double {     // → [0, 2π)
        a = std::fmod(a, 2.0*M_PI);
        if (a < 0.0) a += 2.0*M_PI;
        return a;
    };

    const double R_sol = m_solenoidBoundaryMm / 10.0;   // mm → cm
    const double Z_sol = m_solenoidHalfZmm    / 10.0;   // mm → cm

    double d_oc = std::sqrt(cx_o*cx_o + cy_o*cy_o);  // outer center dist from IP (cm)
    info() << "  [Anal.Entry] " << (isEndcap ? "endcap" : "barrel")
           << "  R_o=" << std::fixed << std::setprecision(2) << R_outer
           << " d_oc=" << d_oc
           << " outerPOCA=" << (d_oc - R_outer) << " cm"   // outer circle closest approach to IP
           << "  pT=" << std::setprecision(3) << pT << " GeV"
           << "  z_ref=" << std::setprecision(2) << z_ref << " cm" << endmsg;

    // Backward direction on the outer circle (+1: increasing angle, for omega < 0);
    // the same sign gives the reversed curvature of the inner helix.
    const double inner_sign = (omega_outer_mm < 0) ? 1.0 : -1.0;

    // ── Solenoid boundary crossing ────────────────────────────────────────────
    double xc_sol, yc_sol, z_sol;

    if (!isEndcap) {
        // ── barrel: outer circle x solenoid cylinder ────────────────────────
        if (d_oc < 1e-6) {
            debug() << "  [Analytical] Barrel: outer circle center at IP – degenerate" << endmsg;
            m_statAnalBarrelNoIsect++;
            return false;
        }
        if (R_outer <= std::abs(d_oc - R_sol) || R_outer >= d_oc + R_sol) {
            debug() << "  [Analytical] Barrel: outer circle does not cross solenoid"
                    << "  d_oc=" << d_oc << " R_outer=" << R_outer
                    << " R_sol=" << R_sol << endmsg;
            m_statAnalBarrelNoIsect++;
            return false;
        }

        double h   = (R_sol*R_sol + d_oc*d_oc - R_outer*R_outer) / (2.0 * d_oc);
        double l   = std::sqrt(std::max(0.0, R_sol*R_sol - h*h));
        double ex  = cx_o / d_oc,  ey = cy_o / d_oc;
        double px  = -ey,          py = ex;

        double P1x = h*ex + l*px,  P1y = h*ey + l*py;
        double P2x = h*ex - l*px,  P2y = h*ey - l*py;

        // the crossing reached first going backward from the hit
        double theta_hit = std::atan2(y_ref - cy_o, x_ref - cx_o);
        auto bwdArc = [&](double Px, double Py) -> double {
            double tc = std::atan2(Py - cy_o, Px - cx_o);
            return (omega_outer_mm < 0) ? normPos(tc - theta_hit)
                                        : normPos(theta_hit - tc);
        };
        double arc1 = bwdArc(P1x, P1y);
        double arc2 = bwdArc(P2x, P2y);
        double arc_back;
        if (arc1 <= arc2) { xc_sol = P1x; yc_sol = P1y; arc_back = arc1; }
        else              { xc_sol = P2x; yc_sol = P2y; arc_back = arc2; }

        double s_bwd = R_outer * arc_back;
        z_sol = z_ref - tanLambda * s_bwd;

        // The crossing must be on the IP side of the hit in z; if not, try the other one.
        if (std::abs(z_sol) > std::abs(z_ref) + 1e-3) {
            double arc_alt, s_bwd_alt, z_sol_alt;
            double xc_alt, yc_alt;
            if (arc1 <= arc2) { xc_alt = P2x; yc_alt = P2y; arc_alt = arc2; }
            else              { xc_alt = P1x; yc_alt = P1y; arc_alt = arc1; }
            s_bwd_alt = R_outer * arc_alt;
            z_sol_alt = z_ref - tanLambda * s_bwd_alt;
            if (std::abs(z_sol_alt) > std::abs(z_ref) + 1e-3) {
                debug() << "  [Analytical] Barrel: both crossings have |z_sol|>|z_ref|"
                        << " z_ref=" << z_ref << " z_sol=" << z_sol
                        << " z_sol_alt=" << z_sol_alt << " – skip" << endmsg;
                m_statAnalBarrelBothZfail++;
                return false;
            }
            debug() << "  [Analytical] Barrel: swapped to alt crossing"
                    << " (|z_sol|=" << std::abs(z_sol) << " > |z_ref|=" << std::abs(z_ref)
                    << ") -> z_sol_alt=" << z_sol_alt << " cm" << endmsg;
            m_statAnalBarrelZswapped++;
            xc_sol = xc_alt; yc_sol = yc_alt;
            s_bwd  = s_bwd_alt;
            z_sol  = z_sol_alt;
        }

        info() << "  [Anal.Brl] crossing=(" << std::fixed << std::setprecision(2) << xc_sol
               << ", " << yc_sol << ") cm  bwd_arc=" << s_bwd
               << " cm  z_sol=" << z_sol << " cm" << endmsg;

    } else {
        // ── end disk: arc s backward from the hit with z(s) = +-Z_sol ────────
        if (std::abs(tanLambda) < 1e-6) {
            debug() << "  [Analytical] Endcap: |tanLambda| too small to extrapolate z" << endmsg;
            m_statAnalEndcapSmallTanL++;
            return false;
        }

        double z_cap = (z_ref < 0.0) ? -Z_sol : +Z_sol;

        // z_ref - tanLambda * s = z_cap
        double s_endcap = (z_ref - z_cap) / tanLambda;
        if (s_endcap < 0.0) {
            // hit already at |z| < Z_sol: use the barrel crossing instead
            debug() << "  [Analytical] Endcap: z_ref=" << z_ref
                    << " cm already inside z_cap=" << z_cap
                    << " cm – falling back to barrel formula" << endmsg;
            m_statAnalEndcapPastCap++;
            if (d_oc < 1e-6 || R_outer <= std::abs(d_oc - R_sol) || R_outer >= d_oc + R_sol) {
                debug() << "  [Analytical] Endcap→Barrel fallback: no barrel crossing" << endmsg;
                return false;
            }
            double h   = (R_sol*R_sol + d_oc*d_oc - R_outer*R_outer) / (2.0 * d_oc);
            double l   = std::sqrt(std::max(0.0, R_sol*R_sol - h*h));
            double ex  = cx_o / d_oc,  ey = cy_o / d_oc;
            double px  = -ey,          py = ex;
            double P1x = h*ex + l*px,  P1y = h*ey + l*py;
            double P2x = h*ex - l*px,  P2y = h*ey - l*py;
            double theta_hit = std::atan2(y_ref - cy_o, x_ref - cx_o);
            auto bwdArc2 = [&](double Px, double Py) -> double {
                double tc = std::atan2(Py - cy_o, Px - cx_o);
                return (omega_outer_mm < 0) ? normPos(tc - theta_hit) : normPos(theta_hit - tc);
            };
            double arc1 = bwdArc2(P1x, P1y), arc2 = bwdArc2(P2x, P2y);
            double s_bwd_fb;
            if (arc1 <= arc2) { xc_sol = P1x; yc_sol = P1y; s_bwd_fb = arc1 * R_outer; }
            else               { xc_sol = P2x; yc_sol = P2y; s_bwd_fb = arc2 * R_outer; }
            z_sol = z_ref - tanLambda * s_bwd_fb;
            if (std::abs(z_sol) > std::abs(z_ref) + 1e-3) {
                debug() << "  [Analytical] Endcap→Barrel fallback: z_sol=" << z_sol << " cm out of range" << endmsg;
                return false;
            }
            info() << "  [Anal.Brl*] fallback crossing=(" << std::fixed << std::setprecision(2) << xc_sol
                   << ", " << yc_sol << ") cm  bwd_arc=" << s_bwd_fb
                   << " cm  z_sol=" << z_sol << " cm" << endmsg;
        } else {
        if (s_endcap > 2.0 * M_PI * R_outer) {
            debug() << "  [Analytical] Endcap: backward arc would exceed one loop – skip" << endmsg;
            m_statAnalEndcapArcLong++;
            return false;
        }

        double theta_ref = std::atan2(y_ref - cy_o, x_ref - cx_o);
        double d_theta   = s_endcap / R_outer;
        double theta_cap = theta_ref + inner_sign * d_theta;

        xc_sol = cx_o + R_outer * std::cos(theta_cap);
        yc_sol = cy_o + R_outer * std::sin(theta_cap);
        z_sol  = z_cap;

        info() << "  [Anal.End] crossing=(" << std::fixed << std::setprecision(2) << xc_sol
               << ", " << yc_sol << ", " << z_sol << ") cm  bwd_arc=" << s_endcap << " cm" << endmsg;
    }   // end if (s_endcap < 0.0)
    }   // end if (!isEndcap)

    // ── tangent and outward radial direction at the crossing ──────────────────
    double rx = (xc_sol - cx_o) / R_outer;
    double ry = (yc_sol - cy_o) / R_outer;
    double tx = (omega_outer_mm < 0) ? -ry :  ry;
    double ty = (omega_outer_mm < 0) ?  rx : -rx;
    debug() << "  [Analytical] Tangent: (" << tx << ", " << ty
            << ")  outward radial: (" << rx << ", " << ry << ")" << endmsg;

    // ── inner field, at the origin ────────────────────────────────────────────
    double B_inner = m_field.magneticField(dd4hep::Position(0.0, 0.0, 0.0)).z()
                     / dd4hep::tesla;
    if (std::abs(B_inner) < 1e-6) {
        debug() << "  [Analytical] Inner B field vanishes (" << B_inner << " T)" << endmsg;
        return false;
    }

    // ── inner circle: same pT, field reversed, so the centre is on the other side ─
    double R_inner = pT * 100.0 / (0.3 * std::abs(B_inner));   // cm
    double cx_i    = xc_sol + R_inner * rx;
    double cy_i    = yc_sol + R_inner * ry;

    double B_outer_approx = pT * 100.0 / (0.3 * R_outer);
    info() << "  [Anal.Inner] R_i=" << std::fixed << std::setprecision(2) << R_inner
           << " cx_i=" << cx_i << " cy_i=" << cy_i
           << " B_i=" << std::setprecision(3) << B_inner << " T"
           << "  d0=" << std::setprecision(2) << (std::sqrt(cx_i*cx_i+cy_i*cy_i)-R_inner) << " cm" << endmsg;

    // ── inner helix at its perigee ────────────────────────────────────────────
    double d0_inner_cm    = std::sqrt(cx_i*cx_i + cy_i*cy_i) - R_inner;
    double phi0_inner     = normAngle(std::atan2(cy_i, cx_i) - inner_sign * M_PI / 2.0);
    double omega_inner_mm = inner_sign / (R_inner * 10.0);   // 1/mm, signed

    // ── z at the perigee: backward arc from the crossing, the shorter way (< pi) ─
    double theta_c  = std::atan2(yc_sol - cy_i, xc_sol - cx_i);
    double theta_ip = std::atan2(-cy_i,  -cx_i);
    double bwd_arc_inner = (omega_inner_mm > 0) ? normPos(theta_c  - theta_ip)
                                                 : normPos(theta_ip - theta_c);
    if (bwd_arc_inner > M_PI)
        bwd_arc_inner = 2.0 * M_PI - bwd_arc_inner;

    double s_sol_to_ip  = R_inner * bwd_arc_inner;
    double z0_inner_cm  = z_sol - tanLambda * s_sol_to_ip;

    debug() << "  [Analytical] Inner perigee:"
            << "  d0=" << d0_inner_cm << " cm"
            << "  phi=" << phi0_inner << " rad"
            << "  omega=" << omega_inner_mm << " /mm"
            << "  z0=" << z0_inner_cm << " cm"
            << "  (bwd arc=" << s_sol_to_ip << " cm)" << endmsg;

    // ── stored as AtIP ────────────────────────────────────────────────────────
    edm4hep::TrackState atIPState = createTrackState(
        d0_inner_cm  * 10.0,    // cm → mm
        phi0_inner,
        omega_inner_mm,
        z0_inner_cm  * 10.0,    // cm → mm
        tanLambda,
        edm4hep::TrackState::AtIP
    );
    finalTrack.addToTrackStates(atIPState);

    info() << "  [Analytical] AtIP (" << (isEndcap ? "endcap" : "barrel") << ") saved:"
           << "  d0=" << d0_inner_cm << " cm"
           << "  phi=" << phi0_inner << " rad"
           << "  omega=" << omega_inner_mm << " /mm"
           << "  z0=" << z0_inner_cm << " cm" << endmsg;

    return true;
}

// ------------------------------------ //
// -------------- GenFit -------------- //
// ------------------------------------ //
genfit::MeasuredStateOnPlane CombinatorialTrackFinder::convertToGenFitState(
    const edm4hep::TrackState& state,
    genfit::AbsTrackRep* rep) const {

    // Internal circle-fit convention (see findTracks); mm, rad, 1/mm
    double d0        = state.D0;          // mm
    double phi       = state.phi;         // rad
    double omega     = state.omega;       // 1/mm
    double tanLambda = state.tanLambda;

    // reference point = first hit (cm for GenFit)
    double refX = state.referencePoint[0] / 10.0;  // cm
    double refY = state.referencePoint[1] / 10.0;
    double refZ = state.referencePoint[2] / 10.0;

    dd4hep::Position fieldPos(refX, refY, refZ);  // cm
    double bField = m_field.magneticField(fieldPos).z() / dd4hep::tesla;  // Tesla

    if (std::abs(omega) < 1e-10 || std::abs(bField) < 1e-6) {
        throw genfit::Exception("convertToGenFitState: omega or B too small", __LINE__, __FILE__);
    }

    // pT [GeV/c] = 0.3 * |Bz| [T] / (|omega| [1/mm] * 1000)
    double pT = 0.3 * std::abs(bField) / (std::abs(omega) * 1000.0);  // GeV/c
    double pz = pT * tanLambda;

    // momentum direction at the reference point: tangent of the circle
    double R_mm = 1.0 / std::abs(omega);
    double centerX_mm, centerY_mm;
    if (omega < 0) {  // positive charge → clockwise
        centerX_mm = -d0 * std::sin(phi) - R_mm * std::cos(phi);
        centerY_mm =  d0 * std::cos(phi) - R_mm * std::sin(phi);
    } else {          // negative charge → counter-clockwise
        centerX_mm = -d0 * std::sin(phi) + R_mm * std::cos(phi);
        centerY_mm =  d0 * std::cos(phi) + R_mm * std::sin(phi);
    }

    double vecX = refX - centerX_mm / 10.0;
    double vecY = refY - centerY_mm / 10.0;
    double vecMag = std::sqrt(vecX*vecX + vecY*vecY);
    if (vecMag < 1e-10) {
        throw genfit::Exception("convertToGenFitState: degenerate center-to-hit vector", __LINE__, __FILE__);
    }
    vecX /= vecMag;
    vecY /= vecMag;

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

    // position and momentum only; the caller sets the seed covariance
    genfit::MeasuredStateOnPlane gfState(rep);
    gfState.setPosMom(posVec, momVec);

    return gfState;
}


edm4hep::TrackState CombinatorialTrackFinder::convertToEDM4hepState(
    const genfit::MeasuredStateOnPlane& state,
    int location) const {

    TVector3 pos = state.getPos();   // cm
    TVector3 mom = state.getMom();   // GeV/c
    double charge = (state.getQop() > 0) ? 1.0 : -1.0;

    double px = mom.X();
    double py = mom.Y();
    double pz = mom.Z();
    double pt = std::sqrt(px*px + py*py);
    double p  = mom.Mag();

    double phi = std::atan2(py, px);

    double tanLambda = (pt > 1e-10) ? pz / pt : 0.0;

    dd4hep::Position fieldPos(pos.X(), pos.Y(), pos.Z());  // cm
    double bField = m_field.magneticField(fieldPos).z() / dd4hep::tesla;

    // omega [1/mm] = charge * 0.3 * |B| / (pT [GeV/c] * 1000)
    double omega = 0.0;
    if (pt > 1e-10 && std::abs(bField) > 1e-6)
        omega = charge * 0.3 * std::abs(bField) / (pt * 1000.0);  // 1/mm

    // d0, z0 of the straight line through the fit position (mm)
    double d0 = (-pos.Y() * std::cos(phi) + pos.X() * std::sin(phi)) * 10.0;  // mm

    double z0 = pos.Z() * 10.0;  // cm → mm

    edm4hep::TrackState outState;
    outState.D0        = d0;
    outState.phi       = phi;
    outState.omega     = omega;
    outState.Z0        = z0;
    outState.tanLambda = tanLambda;
    outState.location  = location;

    outState.referencePoint = edm4hep::Vector3f(
        pos.X() * 10.0f,   // cm → mm
        pos.Y() * 10.0f,
        pos.Z() * 10.0f);

    // ---------- covariance: (x, y, z, px, py, pz) -> 5 track parameters ----------
    // getPosMomCov gives the 6x6 Cartesian covariance (getCov would give the 5x5 one
    // of GenFit's local parameters).
    {
        TMatrixDSym C6sym(6);
        {
            TVector3 tmpPos, tmpMom;
            state.getPosMomCov(tmpPos, tmpMom, C6sym);
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

    // pseudo-inverse: C5 = (J^T J)^{-1} J^T C6 J (J^T J)^{-1}
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

    // EDM4hep covariance; d0 and z0 go from cm to mm
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

genfit::AbsMeasurement* CombinatorialTrackFinder::createGenFitMeasurement(
    const edm4hep::TrackerHitPlane& hit,
    const dd4hep::rec::Surface* surface,
    int hitId,
    genfit::TrackPoint* trackPoint) const {

    // measurement plane of the surface (DD4hep lengths in cm)
    dd4hep::rec::Vector3D u = surface->u();
    dd4hep::rec::Vector3D v = surface->v();
    dd4hep::rec::Vector3D origin = surface->origin();

    TVector3 o(origin.x(), origin.y(), origin.z());
    TVector3 uVec(u.x(), u.y(), u.z());
    TVector3 vVec(v.x(), v.y(), v.z());

    // hit position in the local (u, v) coordinates of the surface
    const auto& pos = hit.getPosition();  // mm
    dd4hep::rec::Vector3D hitGlobal(pos[0], pos[1], pos[2]);  // mm
    dd4hep::rec::Vector2D hitLocal = surface->globalToLocal(dd4hep::mm * hitGlobal);  // cm

    TVectorD hitCoords(2);
    hitCoords[0] = hitLocal[0];  // u coordinate in cm
    hitCoords[1] = hitLocal[1];  // v coordinate in cm

    // resolutions (0.4 mm if not set)
    double sigma_u = hit.getDu();  // mm
    double sigma_v = hit.getDv();  // mm
    if (sigma_u <= 0) sigma_u = 0.4;
    if (sigma_v <= 0) sigma_v = 0.4;

    TMatrixDSym hitCov(2);
    hitCov(0,0) = std::pow(dd4hep::mm * sigma_u, 2);  // cm²
    hitCov(0,1) = 0;
    hitCov(1,0) = 0;
    hitCov(1,1) = std::pow(dd4hep::mm * sigma_v, 2);  // cm²

    int detId = hit.getCellID();

    genfit::PlanarMeasurement* meas =
        new genfit::PlanarMeasurement(hitCoords, hitCov, detId, hitId, trackPoint);

    genfit::SharedPlanePtr plane(new genfit::DetPlane(o, uVec, vVec));
    meas->setPlane(plane, detId);  // cellID as plane ID

    return meas;
}

bool CombinatorialTrackFinder::fitTrackWithGenFit(
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

    int pdgCode = getPDGCode();

    debug() << "Using particle type: " << m_particleType << " (PDG code: " << pdgCode << ")" << endmsg;

    try {
        genfit::AbsTrackRep* rep = new genfit::RKTrackRep(pdgCode);
        // seed state -> position and momentum
        genfit::MeasuredStateOnPlane seedGFState = convertToGenFitState(seedState, rep);
        TVector3 pos = seedGFState.getPos();  // cm
        TVector3 mom = seedGFState.getMom();  // GeV/c

        TVectorD stateVec(6);   // x, y, z, px, py, pz
        stateVec[0] = pos.X();  stateVec[1] = pos.Y();  stateVec[2] = pos.Z();
        stateVec[3] = mom.X();  stateVec[4] = mom.Y();  stateVec[5] = mom.Z();

        // ── Cartesian seed covariance from the seed state's diagonal variances ──
        TMatrixDSym covInit(6);
        covInit.Zero();
        {
            // mm^2, rad^2, (1/mm)^2, mm^2, 1
            double var_d0   = seedState.getCovMatrix(edm4hep::TrackParams::d0,        edm4hep::TrackParams::d0);
            double var_phi  = seedState.getCovMatrix(edm4hep::TrackParams::phi,       edm4hep::TrackParams::phi);
            double var_om   = seedState.getCovMatrix(edm4hep::TrackParams::omega,     edm4hep::TrackParams::omega);
            double var_z0   = seedState.getCovMatrix(edm4hep::TrackParams::z0,        edm4hep::TrackParams::z0);
            double var_tanL = seedState.getCovMatrix(edm4hep::TrackParams::tanLambda, edm4hep::TrackParams::tanLambda);

            // defaults for unset variances
            if (var_d0   < 1e-30) var_d0   = 1.0;
            if (var_phi  < 1e-30) var_phi  = 0.1;
            if (var_z0   < 1e-30) var_z0   = 1.0;
            if (var_tanL < 1e-30) var_tanL = 0.1;

            double phi  = seedState.phi;
            double om   = seedState.omega;   // 1/mm, signed
            double tanL = seedState.tanLambda;
            double pT   = mom.Perp();        // GeV/c
            double pMag = mom.Mag();

            double ux = (pMag > 1e-10) ? mom.X() / pMag : 0.0;
            double uy = (pMag > 1e-10) ? mom.Y() / pMag : 0.0;
            double uz = (pMag > 1e-10) ? mom.Z() / pMag : 0.0;

            // position: var(x) = var(y) = var(d0), var(z) = var(z0), at least 0.1 cm^2
            const double mm2cm2 = 0.01;
            double var_xy = std::max(var_d0 * mm2cm2, 0.1);
            double var_z  = std::max(var_z0 * mm2cm2, 0.1);
            covInit(0,0) = var_xy;
            covInit(1,1) = var_xy;
            covInit(2,2) = var_z;

            // momentum: var(pT) = (pT/omega)^2 var(omega), var(pz) = pT^2 var(tanLambda),
            // and phi rotates (px, py): dpx/dphi = -py, dpy/dphi = px
            double dpT_dom  = (std::abs(om) > 1e-10) ? -pT / om : 0.0;  // GeV·mm
            double var_pT   = dpT_dom * dpT_dom * var_om;      // (GeV/c)²
            double var_pz   = pT * pT * var_tanL;              // (GeV/c)²
            double var_px_phi = mom.Y() * mom.Y() * var_phi;   // (GeV/c)²
            double var_py_phi = mom.X() * mom.X() * var_phi;

            // at least (10% of pT)^2 or (0.1 GeV)^2
            double pT_floor = std::max(pT * 0.1, 0.1);
            double var_pT_floor = pT_floor * pT_floor;

            covInit(3,3) = std::max(var_pT * ux*ux + var_px_phi, var_pT_floor);
            covInit(4,4) = std::max(var_pT * uy*uy + var_py_phi, var_pT_floor);
            covInit(5,5) = std::max(var_pz,                       var_pT_floor);
        }

        debug() << "─── GenFit Input ──────────────────────────────────────────────────────" << endmsg;
        debug() << "  Seed state source : "
                << (seedState.location == edm4hep::TrackState::AtLastHit ? "circle fit" :
                    seedState.location == edm4hep::TrackState::AtFirstHit ? "AtFirstHit" : "other")
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

        // one planar measurement per hit
        for (size_t i = 0; i < hits.size(); ++i) {
            const auto& hit = hits[i];

            const dd4hep::rec::Surface* surface = findSurface(hit);
            if (!surface) {
                warning() << "Could not find surface for hit " << i << endmsg;
                continue;
            }

            genfit::TrackPoint* trackPoint = new genfit::TrackPoint();

            // index of the hit in the input collection
            int hitIdx = -1;
            if (allHits) {
                for (size_t j = 0; j < allHits->size(); ++j) {
                    if ((*allHits)[j].getCellID() == hit.getCellID()) {
                        hitIdx = j;
                        break;
                    }
                }
            }

            genfit::AbsMeasurement* measurement = createGenFitMeasurement(hit, surface, hitIdx, trackPoint);
            if (!measurement) {
                warning() << "Could not create measurement for hit " << i << endmsg;
                delete trackPoint;
                continue;
            }

            trackPoint->addRawMeasurement(measurement);
            finalGFTrack.insertPoint(trackPoint);
        }

        // ── Kalman fit (KalmanFitterRefTrack); GenFit's stderr is captured to detect
        //    Cholesky failures ──────────────────────────────────────────────────
        genfit::KalmanFitterRefTrack fitter;
        fitter.setMaxIterations(m_maxFitIterations);
        fitter.setMinIterations(3);

        debug() << "─── GenFit Kalman Fit ─────────────────────────────────────────────────" << endmsg;
        debug() << "  Fitter            : KalmanFitterRefTrack"
                << "  max_iter=" << m_maxFitIterations << "  min_iter=3" << endmsg;
        debug() << "  Points in track   : " << finalGFTrack.getNumPoints() << endmsg;
        std::string capturedStderr;
        {
            int savedStderr = dup(STDERR_FILENO);
            int pipefd[2];
            if (pipe(pipefd) == 0) {
                dup2(pipefd[1], STDERR_FILENO);
                close(pipefd[1]);

                fitter.processTrack(&finalGFTrack);

                // restore stderr
                fflush(stderr);
                std::cerr.flush();
                dup2(savedStderr, STDERR_FILENO);
                close(savedStderr);

                // read what GenFit printed
                fcntl(pipefd[0], F_SETFL, O_NONBLOCK);
                char buf[4096];
                ssize_t n;
                while ((n = read(pipefd[0], buf, sizeof(buf)-1)) > 0) {
                    buf[n] = '\0';
                    capturedStderr += buf;
                }
                close(pipefd[0]);

                if (!capturedStderr.empty()) {
                    warning() << "GenFit stderr during fit: " << capturedStderr << endmsg;
                }
            } else {
                // no pipe: fit without capture
                close(savedStderr);
                fitter.processTrack(&finalGFTrack);
            }
        }

        // a failed fit: ill-conditioned covariance (seen in GenFit's stderr) or other
        if (!finalGFTrack.getFitStatus()->isFitted()) {
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

        double chi2 = finalGFTrack.getFitStatus()->getChi2();
        int ndf     = finalGFTrack.getFitStatus()->getNdf();
        double chi2ndf = (ndf > 0) ? chi2 / ndf : -1.0;

        debug() << "─── GenFit Result ─────────────────────────────────────────────────────" << endmsg;
        debug() << "  Fit status        : " << (finalGFTrack.getFitStatus()->isFitted() ? "FITTED" : "FAILED") << endmsg;
        debug() << "  chi2=" << std::fixed << std::setprecision(3) << chi2
                << "  ndf=" << ndf
                << "  chi2/ndf=" << chi2ndf << endmsg;

        // run statistics: mean chi2/ndf (stored x 1000 as an integer), fits with
        // chi2/ndf >= 1000 counted separately
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

        finalTrack.setChi2(chi2);
        finalTrack.setNdf(ndf);

        // copy the hits to the output collection
        if (finalTrack.trackerHits_size() == 0) {
            for (const auto& hit : hits) {
                auto outHit = outputHits.create();
                outHit.setCellID(hit.getCellID()); outHit.setTime(hit.getTime());
                outHit.setEDep(hit.getEDep()); outHit.setEDepError(hit.getEDepError());
                outHit.setPosition(hit.getPosition()); outHit.setCovMatrix(hit.getCovMatrix());
                outHit.setDu(hit.getDu()); outHit.setDv(hit.getDv());
                outHit.setU(hit.getU()); outHit.setV(hit.getV());
                propagateLink(outHit, hit);
                finalTrack.addToTrackerHits(outHit);
            }
        }

        // The GenFit state at the first measurement is stored as AtOther; the helix
        // states (AtFirstHit, AtLastHit, AtIP) are not touched.
        auto hasLocation = [&](int loc) {
            for (int j = 0; j < finalTrack.trackStates_size(); ++j)
                if (finalTrack.getTrackStates(j).location == loc) return true;
            return false;
        };

        if (finalGFTrack.getNumPoints() > 0) {
            try {
                genfit::MeasuredStateOnPlane stateFirst = finalGFTrack.getFittedState(0);
                edm4hep::TrackState gfFirst = convertToEDM4hepState(stateFirst,
                                                    edm4hep::TrackState::AtOther);
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

        // the GenFit state at the last measurement is not stored
        if (finalGFTrack.getNumPoints() > 1 && hasLocation(edm4hep::TrackState::AtLastHit)) {
            debug() << "GenFit state at the last hit not stored (AtLastHit holds the helix)" << endmsg;
        }

        debug() << "  Track states total : " << finalTrack.trackStates_size()
                << "  hits: " << finalTrack.trackerHits_size() << endmsg;
        debug() << "───────────────────────────────────────────────────────────────────────" << endmsg;

        return true;

    } catch (genfit::Exception& e) {
        std::string exMsg = e.what();
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

bool CombinatorialTrackFinder::fitCircleNHits(
    const std::vector<edm4hep::TrackerHitPlane>& hits,
    double& x0, double& y0, double& radius, double& chi2ndf,
    Eigen::Matrix3d& fitCov3x3,
    std::vector<double>& residualsCm) const {

    const int N = static_cast<int>(hits.size());
    if (N < 3) {
        warning() << "fitCircleNHits: need at least 3 hits, got " << N << endmsg;
        return false;
    }

    // positions (cm) and weights: one radial residual per hit, so one resolution (du)
    std::vector<Eigen::Vector2d> pts(N);
    std::vector<double> sigma(N), w(N);
    for (int i = 0; i < N; ++i) {
        const auto& pos = hits[i].getPosition();
        pts[i] = { pos[0] / 10.0, pos[1] / 10.0 };
        double su = hits[i].getDu();
        if (su <= 0) su = m_sigmaHitDefault.value() * 10.0;  // property is in cm, du is in mm
        sigma[i] = su / 10.0;  // mm -> cm
        w[i]     = 1.0 / (sigma[i] * sigma[i]);
    }

    // ── Step 1: Algebraic WLS circle fit ──────────────────────────────────
    // Circle equation: x² + y² + D·x + E·y + F = 0
    Eigen::MatrixXd A(N, 3);
    Eigen::VectorXd b(N);
    for (int i = 0; i < N; ++i) {
        double x = pts[i].x(), y = pts[i].y();
        A(i,0) = w[i] * x;
        A(i,1) = w[i] * y;
        A(i,2) = w[i] * 1.0;
        b(i)   = -w[i] * (x*x + y*y);
    }
    Eigen::Vector3d sol = (A.transpose() * A).ldlt().solve(A.transpose() * b);
    double D = sol(0), E = sol(1), F = sol(2);
    x0 = -D / 2.0;
    y0 = -E / 2.0;
    double disc = D*D + E*E - 4*F;
    if (disc <= 0) {
        debug() << "fitCircleNHits: degenerate algebraic fit (disc=" << disc << ")" << endmsg;
        return false;
    }
    radius = std::sqrt(disc) / 2.0;

    // ── Step 2: Gauss-Newton geometric refinement ─────────────────────────
    Eigen::Vector3d params(x0, y0, radius);
    const int maxIter = 10;
    const double tol  = 1e-6;
    for (int iter = 0; iter < maxIter; ++iter) {
        double cx = params(0), cy = params(1), r = params(2);
        Eigen::MatrixXd J(N, 3);
        Eigen::VectorXd res(N);
        for (int i = 0; i < N; ++i) {
            double dx = pts[i].x() - cx, dy = pts[i].y() - cy;
            double dist = std::sqrt(dx*dx + dy*dy);
            double sqw  = std::sqrt(w[i]);
            res(i) = sqw * (dist - r);
            if (dist > 1e-12) {
                J(i,0) = -sqw * dx / dist;
                J(i,1) = -sqw * dy / dist;
                J(i,2) = -sqw;
            } else {
                J(i,0) = J(i,1) = J(i,2) = 0.0;
            }
        }
        Eigen::Vector3d delta = -(J.transpose() * J).ldlt().solve(J.transpose() * res);
        params += delta;
        if (delta.norm() < tol) break;
    }
    x0     = params(0);
    y0     = params(1);
    radius = params(2);
    if (radius <= 0) {
        debug() << "fitCircleNHits: non-positive radius after refinement" << endmsg;
        return false;
    }

    // ── Step 3: Covariance on [x0, y0, R] from final Jacobian ────────────
    {
        Eigen::MatrixXd Jf(N, 3);
        for (int i = 0; i < N; ++i) {
            double dx = pts[i].x() - x0, dy = pts[i].y() - y0;
            double dist = std::sqrt(dx*dx + dy*dy);
            double sqw  = std::sqrt(w[i]);
            if (dist > 1e-12) {
                Jf(i,0) = -sqw * dx / dist;
                Jf(i,1) = -sqw * dy / dist;
                Jf(i,2) = -sqw;
            } else {
                Jf(i,0) = Jf(i,1) = Jf(i,2) = 0.0;
            }
        }
        fitCov3x3 = (Jf.transpose() * Jf).ldlt().solve(Eigen::Matrix3d::Identity());
    }

    // ── Step 4: Per-hit residuals and chi2/NDF ────────────────────────────
    residualsCm.resize(N);
    double chi2 = 0.0;
    for (int i = 0; i < N; ++i) {
        double dx   = pts[i].x() - x0, dy = pts[i].y() - y0;
        double res  = std::abs(std::sqrt(dx*dx + dy*dy) - radius);
        residualsCm[i] = res;
        chi2 += (res / sigma[i]) * (res / sigma[i]);
    }
    int ndf = N - 3;
    if (ndf < 1) ndf = 1;
    chi2ndf = chi2 / ndf;

    debug() << "fitCircleNHits (" << N << " hits): center=(" << x0 << "," << y0
            << ") cm  R=" << radius << " cm  chi2/ndf=" << chi2ndf << endmsg;
    return true;
}
// PDG code of the particle hypothesis
int CombinatorialTrackFinder::getPDGCode() const {
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
// SOLENOID BOUNDARY SCAN
// ============================================================================
double CombinatorialTrackFinder::findSolenoidBoundary() const {
    // Bz along +x in 1 mm steps around SolenoidRadius to bracket the sign change,
    // then a binary search to 0.01 mm.

    const double nominalMm  = m_solenoidRadius.value();
    const double stepMm     = 1.0;
    const double scanRangeMm = 600.0;   // look up to ±600 mm from nominal

    // B_z query helper (radius in mm, z=0)
    auto getBz = [&](double R_mm) -> double {
        double R_cm = R_mm / 10.0;
        return m_field.magneticField(dd4hep::Position(R_cm, 0.0, 0.0)).z() / dd4hep::tesla;
    };

    double Bz0 = getBz(nominalMm);
    info() << "  [SolScan] B_z at nominal R=" << nominalMm << " mm: " << Bz0 << " T" << endmsg;

    // Find bracket in both directions and pick the one that finds it first
    double R_lo = -1.0, R_hi = -1.0;

    // Try outward scan
    {
        double Bz_prev = Bz0, R_prev = nominalMm;
        for (double R = nominalMm + stepMm; R <= nominalMm + scanRangeMm; R += stepMm) {
            double Bz = getBz(R);
            if (Bz * Bz_prev < 0.0) { R_lo = R_prev; R_hi = R; break; }
            Bz_prev = Bz; R_prev = R;
        }
    }
    // Try inward scan if outward failed or to compare
    if (R_lo < 0.0) {
        double Bz_prev = Bz0, R_prev = nominalMm;
        for (double R = nominalMm - stepMm;
             R >= std::max(10.0, nominalMm - scanRangeMm); R -= stepMm) {
            double Bz = getBz(R);
            if (Bz * Bz_prev < 0.0) { R_lo = R; R_hi = R_prev; break; }
            Bz_prev = Bz; R_prev = R;
        }
    }

    if (R_lo < 0.0) {
        warning() << "  [SolScan] No B_z sign change found within ±" << scanRangeMm
                  << " mm of nominal R=" << nominalMm << " mm — using nominal." << endmsg;
        return nominalMm;
    }

    // Log the 1-mm bracket
    info() << "  [SolScan] Coarse bracket: R_lo=" << R_lo << " mm (Bz=" << getBz(R_lo)
           << " T)  R_hi=" << R_hi << " mm (Bz=" << getBz(R_hi) << " T)" << endmsg;

    // Print field at every 1 mm within the bracket ± 5 mm for context
    info() << "  [SolScan] Fine field profile near boundary:" << endmsg;
    for (double R = std::max(10.0, R_lo - 5.0); R <= R_hi + 5.0; R += 1.0) {
        info() << "    R=" << std::setw(7) << std::fixed << std::setprecision(1) << R
               << " mm  Bz=" << std::setw(9) << std::setprecision(4) << getBz(R) << " T" << endmsg;
    }

    // Binary search for sub-mm precision
    for (int iter = 0; iter < 40; ++iter) {
        double R_mid = 0.5 * (R_lo + R_hi);
        double Bz_mid = getBz(R_mid);
        double Bz_lo  = getBz(R_lo);
        if (Bz_mid * Bz_lo < 0.0) R_hi = R_mid; else R_lo = R_mid;
        if (R_hi - R_lo < 0.01) break;   // converged to 0.01 mm
    }

    double boundary = 0.5 * (R_lo + R_hi);
    info() << "  [SolScan] Boundary (Bz=0): R = " << std::fixed << std::setprecision(2)
           << boundary << " mm  (Bz_lo=" << std::setprecision(4) << getBz(R_lo)
           << " T, Bz_hi=" << getBz(R_hi) << " T)" << endmsg;
    return boundary;
}

// Same as findSolenoidBoundary(), along the z axis.
double CombinatorialTrackFinder::findSolenoidBoundaryZ() const {
    const double nominalMm   = m_solenoidHalfZ.value();
    const double stepMm      = 1.0;
    const double scanRangeMm = 600.0;

    // B_z query along z-axis (x=y=0)
    auto getBz = [&](double Z_mm) -> double {
        double Z_cm = Z_mm / 10.0;
        return m_field.magneticField(dd4hep::Position(0.0, 0.0, Z_cm)).z() / dd4hep::tesla;
    };

    double Bz0 = getBz(nominalMm);
    info() << "  [SolScanZ] B_z at nominal |Z|=" << nominalMm << " mm: " << Bz0 << " T" << endmsg;

    // If B_z is already ≈ 0 at the nominal, we are exactly at the boundary — no scan needed.
    if (std::abs(Bz0) < 0.01) {
        info() << "  [SolScanZ] Nominal |Z|=" << nominalMm
               << " mm is already at Bz=0 — using nominal as endcap boundary." << endmsg;
        return nominalMm;
    }

    double Z_lo = -1.0, Z_hi = -1.0;

    // Outward scan (increasing |z|)
    {
        double Bz_prev = Bz0, Z_prev = nominalMm;
        for (double Z = nominalMm + stepMm; Z <= nominalMm + scanRangeMm; Z += stepMm) {
            double Bz = getBz(Z);
            if (Bz * Bz_prev < 0.0) { Z_lo = Z_prev; Z_hi = Z; break; }
            Bz_prev = Bz; Z_prev = Z;
        }
    }
    // Inward scan (decreasing |z|) if outward failed
    if (Z_lo < 0.0) {
        double Bz_prev = Bz0, Z_prev = nominalMm;
        for (double Z = nominalMm - stepMm;
             Z >= std::max(10.0, nominalMm - scanRangeMm); Z -= stepMm) {
            double Bz = getBz(Z);
            if (Bz * Bz_prev < 0.0) { Z_lo = Z; Z_hi = Z_prev; break; }
            Bz_prev = Bz; Z_prev = Z;
        }
    }

    if (Z_lo < 0.0) {
        warning() << "  [SolScanZ] No B_z sign change found within ±" << scanRangeMm
                  << " mm of nominal |Z|=" << nominalMm << " mm — using nominal." << endmsg;
        return nominalMm;
    }

    info() << "  [SolScanZ] Coarse bracket: Z_lo=" << Z_lo << " mm (Bz=" << getBz(Z_lo)
           << " T)  Z_hi=" << Z_hi << " mm (Bz=" << getBz(Z_hi) << " T)" << endmsg;

    // Fine profile ±5 mm around bracket
    info() << "  [SolScanZ] Fine field profile near endcap boundary:" << endmsg;
    for (double Z = std::max(10.0, Z_lo - 5.0); Z <= Z_hi + 5.0; Z += 1.0) {
        info() << "    |Z|=" << std::setw(7) << std::fixed << std::setprecision(1) << Z
               << " mm  Bz=" << std::setw(9) << std::setprecision(4) << getBz(Z) << " T" << endmsg;
    }

    // Binary search to 0.01 mm
    for (int iter = 0; iter < 40; ++iter) {
        double Z_mid  = 0.5 * (Z_lo + Z_hi);
        double Bz_mid = getBz(Z_mid);
        double Bz_lo  = getBz(Z_lo);
        if (Bz_mid * Bz_lo < 0.0) Z_hi = Z_mid; else Z_lo = Z_mid;
        if (Z_hi - Z_lo < 0.01) break;
    }

    double boundary = 0.5 * (Z_lo + Z_hi);
    info() << "  [SolScanZ] Boundary (Bz=0): |Z| = " << std::fixed << std::setprecision(2)
           << boundary << " mm  (Bz_lo=" << std::setprecision(4) << getBz(Z_lo)
           << " T, Bz_hi=" << getBz(Z_hi) << " T)" << endmsg;
    return boundary;
}

// ============================================================================
// TRUTH COMPARISON
// ============================================================================
void CombinatorialTrackFinder::compareStatesWithTruth(
    const edm4hep::Track& track,
    const edm4hep::MCParticleCollection& mcParticles,
    int trackNumber) const {

    // MC particle with the highest pT (single-particle samples)
    const edm4hep::MCParticle* mu = nullptr;
    double bestPT2 = -1.0;
    for (const auto& p : mcParticles) {
        const auto& m = p.getMomentum();
        double pT2 = m.x*m.x + m.y*m.y;
        if (pT2 > bestPT2) { bestPT2 = pT2; mu = &p; }
    }
    if (!mu) {
        debug() << "[TruthCmp] No MCParticle for track " << trackNumber << endmsg;
        return;
    }

    const auto& vtx = mu->getVertex();    // mm (EDM4hep convention)
    const auto& mom = mu->getMomentum(); // GeV/c
    double pT_true    = std::sqrt(mom.x*mom.x + mom.y*mom.y);
    double phi0_true  = std::atan2(mom.y, mom.x);
    double tanL_true  = (pT_true > 1e-9) ? mom.z / pT_true : 0.0;
    // straight-line d0 and z0 of the production vertex
    double d0_true_mm = -vtx.x * std::sin(phi0_true) + vtx.y * std::cos(phi0_true);
    double z0_true_mm = vtx.z;

    // fields inside (origin) and outside (50 cm beyond the boundary) the solenoid
    double B_ip = m_field.magneticField(dd4hep::Position(0.0, 0.0, 0.0)).z() / dd4hep::tesla;
    double R_outer_cm = (m_solenoidBoundaryMm / 10.0) + 50.0;
    double B_outer = m_field.magneticField(dd4hep::Position(R_outer_cm, 0.0, 0.0)).z()
                     / dd4hep::tesla;

    // ── Print header ─────────────────────────────────────────────────────────
    info() << "\n"
           << "╔══════════════════════════════════════════════════════════════════════════╗\n"
           << "║  TRUTH COMPARISON  Track " << std::setw(3) << trackNumber
           << "   PDG=" << std::setw(5) << mu->getPDG()
           << "   pT=" << std::fixed << std::setprecision(2) << std::setw(7) << pT_true << " GeV/c"
           << "   tanλ=" << std::setw(6) << std::setprecision(3) << tanL_true << "            ║\n"
           << "║  MC vertex: (" << std::setw(7) << std::setprecision(2) << vtx.x
           << ", " << std::setw(7) << vtx.y
           << ", " << std::setw(7) << vtx.z << ") mm"
           << "   B_ip=" << std::setw(6) << std::setprecision(3) << B_ip
           << " T   B_outer=" << std::setw(6) << B_outer << " T           ║\n"
           << "╠═════════════════╦═════════════╦═════════════╦════════════╦═════════════╣\n"
           << "║  State          ║  D0 (mm)    ║  Z0 (mm)    ║  phi (rad) ║  pT (GeV/c) ║\n"
           << "╠═════════════════╬═════════════╬═════════════╬════════════╬═════════════╣\n";

    // Truth row
    info() << "║  Truth (MC IP)  ║"
           << std::setw(11) << std::setprecision(3) << d0_true_mm << "  ║"
           << std::setw(11) << z0_true_mm << "  ║"
           << std::setw(10) << std::setprecision(4) << phi0_true << "  ║"
           << std::setw(11) << std::setprecision(2) << pT_true << "  ║" << endmsg;
    info() << "╠═════════════════╬═════════════╬═════════════╬════════════╬═════════════╣\n";

    // Helper: print one state row + residual row
    struct StateRow { const char* label; int loc; double B_for_pT; };
    StateRow rows[] = {
        { "AtFirstHit",   edm4hep::TrackState::AtFirstHit, B_outer },
        { "AtVertex(RK4)",edm4hep::TrackState::AtVertex,   B_ip    },
        { "AtIP(analyt)", edm4hep::TrackState::AtIP,       B_ip    },
    };
    for (auto& row : rows) {
        bool found = false;
        for (int i = 0; i < track.trackStates_size(); ++i) {
            auto ts = track.getTrackStates(i);
            if (ts.location != row.loc) continue;
            found = true;
            double pT_reco = (std::abs(ts.omega) > 1e-12 && std::abs(row.B_for_pT) > 1e-6)
                ? 0.3 * std::abs(row.B_for_pT) / (std::abs(ts.omega) * 1000.0) : 0.0;
            double resD0   = ts.D0 - d0_true_mm;
            double resZ0   = ts.Z0 - z0_true_mm;
            double resPhi  = ts.phi - phi0_true;
            double resPT   = pT_reco - pT_true;
            // phi residual: wrap to (−π, π]
            while (resPhi >  M_PI) resPhi -= 2*M_PI;
            while (resPhi < -M_PI) resPhi += 2*M_PI;
            info() << "║  " << std::left << std::setw(15) << row.label << "║"
                   << std::right << std::setw(11) << std::setprecision(3) << ts.D0 << "  ║"
                   << std::setw(11) << ts.Z0 << "  ║"
                   << std::setw(10) << std::setprecision(4) << ts.phi << "  ║"
                   << std::setw(11) << std::setprecision(2) << pT_reco << "  ║" << endmsg;
            info() << "║  " << std::left << std::setw(15) << "  Δ (reco−truth)" << "║"
                   << std::right << std::setw(9) << std::setprecision(3)
                   << (resD0 >= 0 ? "+" : "") << resD0 << "   ║"
                   << std::setw(9) << (resZ0 >= 0 ? "+" : "") << resZ0 << "   ║"
                   << std::setw(8) << std::setprecision(4)
                   << (resPhi >= 0 ? "+" : "") << resPhi << "    ║"
                   << std::setw(9) << std::setprecision(2)
                   << (resPT >= 0 ? "+" : "") << resPT << "   ║" << endmsg;
            break;
        }
        if (!found)
            info() << "║  " << std::left << std::setw(15) << row.label
                   << "║  (not saved)                                              ║" << endmsg;
    }
    info() << "╚═════════════════╩═════════════╩═════════════╩════════════╩═════════════╝\n"
           << endmsg;
}

// ============================================================================
// INNER SOLENOID PROPAGATION IMPLEMENTATIONS
// ============================================================================
bool CombinatorialTrackFinder::isInsideSolenoid(const Eigen::Vector3d& pos) const {
    double r = std::sqrt(pos.x() * pos.x() + pos.y() * pos.y());
    double z = std::abs(pos.z());
    return (r < m_solenoidBoundaryMm / 10.0) && (z < SOLENOID_Z_HALF_CM);
}

// True when no hit of the track lies inside the solenoid.
bool CombinatorialTrackFinder::allHitsOutsideSolenoid(const std::vector<edm4hep::TrackerHitPlane>& trackHits) const {
    for (const auto& hit : trackHits) {
        const auto& pos = hit.getPosition();
        const Eigen::Vector3d posCm(pos[0] / 10.0, pos[1] / 10.0, pos[2] / 10.0);
        if (isInsideSolenoid(posCm))
            return false;
    }
    return true;
}

void CombinatorialTrackFinder::rk4PropagationStep(
    Eigen::Vector3d& pos,
    Eigen::Vector3d& mom,
    double stepSize,
    double chargeSign) const {

    double pMagInitial = mom.norm();
    // dp/ds [GeV/c per cm] = charge * 0.003 * (p_hat x B [T])
    const double kappa_factor = 0.3 / 100.0;

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

        // dp/ds = q * kappa * (p̂ × B)  where q = ±1 particle charge
        Eigen::Vector3d dpds = chargeSign * kappa_factor * pHat.cross(B_vec);

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

InnerTrajectory CombinatorialTrackFinder::propagateToInner(
    const Eigen::Vector3d& outerPos,
    const Eigen::Vector3d& outerMom,
    double targetRadius,
    double chargeSign) const {

    InnerTrajectory result;
    result.success   = false;
    result.numSteps  = 0;
    result.arcLength = 0.0;

    // Starting position must be outside the solenoid
    if (isInsideSolenoid(outerPos)) {
        result.message = "Starting position already inside solenoid";
        return result;
    }

    // Backward propagation: negate the momentum so the track retraces its path
    // toward the IP.  After convergence the final momentum is negated again to
    // restore the physical (outward) direction at the endpoint.
    Eigen::Vector3d r = outerPos;
    Eigen::Vector3d p = -outerMom;          // backward = reversed momentum
    const double FWD_STEP = 5.0;            // cm; positive forward steps along -p
    const int    MAX_STEPS = 5000;

    debug() << "[RK4] Starting backward inner propagation"
            << "  R=" << std::sqrt(r.x()*r.x() + r.y()*r.y()) << " cm"
            << "  |p|=" << outerMom.norm() << " GeV/c"
            << "  charge=" << chargeSign << endmsg;

    for (int step = 0; step < MAX_STEPS; ++step) {
        double rxy = std::sqrt(r.x()*r.x() + r.y()*r.y());

        // Success: inside solenoid and within target radius
        if (isInsideSolenoid(r) && rxy < targetRadius) {
            result.success        = true;
            result.message        = "Reached target radius";
            result.finalPosition  = r;
            result.finalMomentum  = -p;     // restore physical direction
            result.arcLength      = result.arcLength;
            result.numSteps       = step;
            debug() << "[RK4] Success at step " << step
                    << "  R=" << rxy << " cm" << endmsg;
            break;
        }

        rk4PropagationStep(r, p, FWD_STEP, chargeSign);
        result.arcLength += FWD_STEP;
        result.numSteps   = step + 1;

        if (step % 200 == 0 && msgLevel(MSG::DEBUG)) {
            double B = m_field.magneticField(dd4hep::Position(r.x(), r.y(), r.z())).z()
                       / dd4hep::tesla;
            debug() << "[RK4] step " << step << "  R=" << rxy << " cm  B=" << B << " T" << endmsg;
        }
    }

    if (!result.success) {
        result.message        = "Max steps reached";
        result.finalPosition  = r;
        result.finalMomentum  = -p;
    }
    return result;
}

const edm4hep::MCParticle* CombinatorialTrackFinder::findTruthParticle(
    const edm4hep::Track& recoTrack,
    const edm4hep::MCParticleCollection& mcParticles) const {

    if (recoTrack.trackerHits_size() == 0) {
        debug() << "Track has no hits - cannot find truth" << endmsg;
        return nullptr;
    }

    // Rough match without sim links: MC particle whose production vertex is within
    // 50 cm of most hits.
    std::map<const edm4hep::MCParticle*, int> hitCounts;

    for (int i = 0; i < recoTrack.trackerHits_size(); ++i) {
        const auto& hit = recoTrack.getTrackerHits(i);

        for (const auto& mcPart : mcParticles) {
            const auto& vertex = mcPart.getVertex();
            const auto& momentum = mcPart.getMomentum();

            const auto& hitPos = hit.getPosition();
            double dist = std::sqrt(
                (hitPos[0]/10.0 - vertex.x) * (hitPos[0]/10.0 - vertex.x) +
                (hitPos[1]/10.0 - vertex.y) * (hitPos[1]/10.0 - vertex.y) +
                (hitPos[2]/10.0 - vertex.z) * (hitPos[2]/10.0 - vertex.z)
            );

            if (dist < 50.0) {
                hitCounts[&mcPart]++;
            }
        }
    }

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

void CombinatorialTrackFinder::validateTrackWithTruth(
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

    // AtFirstHit state and, if present, the GenFit state (AtOther)
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
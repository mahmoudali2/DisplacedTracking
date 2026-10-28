# CombinatorialTrackFinder (standalone)

A Gaudi/Key4hep algorithm that reconstructs tracks from the hits of a single tracking
system, without needing hits from any other sub-detector. It was written for the IDEA
muon system at FCC-ee, where muons from long-lived particle decays (e.g. HNL → μμ, decay
vertex at O(1–5 m)) leave hits only in the muon chambers. The same code runs on silicon
readouts (IDEA vertex + Si wrapper, CLD-like inner tracker) for prompt and displaced
tracks.

This repository is the standalone build of the component. The same sources live in
k4RecTracker under `Tracking/`, where the component was first called
`StandaloneMuonTracking`. The older `DisplacedTracking` component that used to be here was
replaced by it. That version is still in the git history (commit `8a9cb10` and earlier).

## Algorithm overview

The `Strategy` property chooses how candidates are found:

* **`Exhaustive`** (default): N-hit combinatorics. Hits are grouped per layer, k-hit
  combinations are formed for k from `MaxCombinatorialHits` down to `MinTrackHits`, and
  each one is tested with a weighted circle fit plus an sz line fit. Road cuts in φ,
  distance and tanλ keep the search small. Outliers are then removed and a neighbour-track
  check is applied. This mode works well on sparse muon-system events but scales as C(N,k),
  so it is not usable on dense silicon events.
* **`SeedAndExtend`**: forms triplet seeds on the inner layers (`SeedFinder`) and extends
  them layer by layer with an analytic helix (`HelixExtrapolator`) through the layer
  geometry read from DD4hep (`LayerGeometry`). `SeedAndExtendMode` selects one of two
  variants:
  * `Classic`: all seeds are built first, then extended.
  * `BestFirst`: seeds are taken in order of best kink-χ² first and claim their hits at
    once. The kink χ² accounts for multiple scattering, using the per-layer x/X0 measured
    from the DD4hep material (`TrackQuality`). This mode adds a confirmation pass, merges
    track pieces, flags loopers and checks short tracks for holes and minimum pT.
    `SeedPasses` (e.g. `[4, 3]`) runs iterative passes with fewer seed hits each time.

The found tracks can then be refitted with GenFit (Kalman/DAF, with DD4hep field and
material through `GenfitField` / `GenfitMaterialInterface`). With `DoInnerPropagation`,
tracks are propagated inward through the solenoid boundary, either analytically or with
RK4 steps.

`MergeTrackerHitPlanes` is a small helper. It concatenates several hit collections and
their reco→sim link collections, for example vertex barrel/disk + SiWrB/SiWrD, into one
collection that a single `CombinatorialTrackFinder` instance can read.

## Layout

```
components/CombinatorialTrackFinder.{h,cpp}   Gaudi algorithm
components/MergeTrackerHitPlanes.cpp          hit-collection merger
include/, src/
  TrackingPrimitives                          circle/helix fit primitives, hit helpers
  LayerGeometry                               layer surfaces + per-layer material from DD4hep
  HelixExtrapolator                           analytic helix–layer crossings
  SeedFinder                                  triplet seeding, seed-and-extend
  TrackQuality                                multiple-scattering-aware kink chi2
  genfit_interfaces/GenfitField, GenfitMaterialInterface   DD4hep → GenFit adapters
test/testHelixExtrapolator, test/testTrackQuality          standalone unit tests (ctest)
options/runCombinatorialTrackFinder.py        muon-system steering example
options/runGlobalTrackerTracking.py           IDEA vertex + Si wrapper example (with merger)
```

## Dependencies

* Key4hep stack: DD4hep, EDM4hep, Gaudi, k4FWCore, podio, ROOT, Eigen3
* GenFit2. `cmake/FindGenFit.cmake` looks for it on `CMAKE_PREFIX_PATH`, then through the
  environment hints `GENFIT_INSTALL_DIR` (install prefix) or `GENFIT_INC` / `GENFIT_LIB`
  (include directory / full library path).

Some recent nightlies (e.g. 2026-10-01) do not ship GenFit. In that case take it from an
older release and pass it through the hints. `env1001.sh` in the parent directory does this.

## Build

```bash
source /cvmfs/sw-nightlies.hsf.org/key4hep/setup.sh   # or env1001.sh
cmake -S . -B build -DCMAKE_INSTALL_PREFIX=install
cmake --build build -j16 --target install
(cd build && ctest --output-on-failure)               # unit tests
```

## Running

```bash
source /cvmfs/sw-nightlies.hsf.org/key4hep/setup.sh   # same environment as the build
source setup.sh                                        # puts install/ on LD_LIBRARY_PATH / PYTHONPATH
k4run options/runCombinatorialTrackFinder.py
```

The compact in `GeoSvc` **must** match the geometry the input file was digitised with.
With the wrong compact the job still runs, but every step after pattern recognition goes
wrong.

## Main configuration parameters

| Parameter | Default | Description |
|---|---|---|
| `DetectorName` / `DetectorNames` | `"Tracker"` / `[]` | DetElement(s) whose surfaces are used; list several when they share one readout |
| `EncodingStringParameterName` | `"GlobalTrackerReadoutID"` | readout whose cellID encoding is decoded |
| `RegionFieldName`, `LayerFieldName` | `"type"`, `"layer"` | cellID field names (`"side"` for GlobalTrackerReadoutID) |
| `UseSystemInCompositeID` | `-1` | fold `system` into the layer key (-1 = auto when several detectors are used) |
| `Strategy` | `"Exhaustive"` | `Exhaustive` or `SeedAndExtend` |
| `SeedAndExtendMode` | `"Classic"` | `Classic` or `BestFirst` |
| `NumSeedLayers`, `SeedPasses` | `3`, `[]` | seed size; iterative passes (e.g. `[4, 3]`) |
| `SeedMaxDeltaPhi`, `SeedMaxHitDistance` | `0.3`, `100.0` | seed road cuts (displaced tracks need loose φ) |
| `MaxMissedLayers`, `ExtensionWindow` | `1`, `1.0` | extension tolerance |
| `ShortTrackMaxHits`, `ShortTrackMinPt` | `4`, `0.0` | rules for short tracks; the minimum pT removes most random triplets |
| `MinTrackHits`, `MaxCombinatorialHits` | `3`, `7` | k range of the exhaustive search |
| `NHitMaxChi2NDF` | `5.0` | χ²/ndf cut on the N-hit fit |
| `UseGenFit`, `MaxFitIterations` | `true`, `4` | GenFit refit |
| `ParticleType` | `"pion"` | mass hypothesis for material effects |
| `DoInnerPropagation`, `DoAnalyticalPropagation` | `false`, `true` | inward propagation through the solenoid |

The header `components/CombinatorialTrackFinder.h` documents every property.

## Input / output collections

| Property | Type | Example |
|---|---|---|
| `InputHitCollection` | `TrackerHitPlaneCollection` | `MSTrackerHits` |
| `InputRecoSimLinkCollection` | `TrackerHitSimTrackerHitLinkCollection` | `MSTrackerHitRelations` |
| `OutputTrackCollection` | `TrackCollection` | `StandaloneMuonTracks` |
| `OutputHitCollection` | `TrackerHitPlaneCollection` | `StandaloneMuonTrackHits` |
| `OutputRecoSimLinkCollection` | `TrackerHitSimTrackerHitLinkCollection` | `StandaloneMuonTrackHitSimLinks` |

## Author

Mahmoud Al-Thakeel

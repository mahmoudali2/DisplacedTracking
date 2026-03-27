# DisplacedTracking

A Gaudi/Key4hep algorithm for reconstructing displaced tracks in the IDEA muon system at FCC-ee. Designed for long-lived particle (LLP) searches where signal tracks originate from a displaced vertex and leave hits only in the muon system, bypassing the inner tracker entirely.

## Physics motivation

Targets signatures such as HNL → μμ where the heavy neutral lepton decays at a displaced vertex (O(200–500 cm) from the IP). The decay products appear as muon-system-only tracks with large impact parameters (d0, z0) relative to the interaction point.

## Algorithm overview

1. **Triplet seeding** — forms seeds from hits in the three innermost muon layers using a circle fit and cross-product charge determination
2. **Hit extension** — extends seeds to a 4th layer hit using a distance+angle compatibility check
3. **4-hit circle fit** — weighted least-squares circle fit to all 4 hits for improved pT and charge
4. **GenFit Kalman fit** — full track fit with material effects using the GenFit library (DD4hep material interface)
5. **Inner propagation** — optional extrapolation of the fitted track back toward the solenoid

Track states stored per track:
- `AtFirstHit` — 3-hit analytical seed
- `AtLastHit` — 4-hit WLS circle fit
- `AtOther` — GenFit Kalman fit result
- `AtVertex` — inner propagation result (when enabled)

## Dependencies

- [Key4hep](https://github.com/key4hep/key4hep-spack) nightly stack (DD4hep, EDM4hep, Gaudi, k4FWCore, podio, ROOT, Eigen3)
- [GenFit2](https://github.com/GenFit/GenFit) — built and installed separately

## Build

```bash
source /cvmfs/sw-nightlies.hsf.org/key4hep/setup.sh

mkdir build && cd build
cmake .. -DCMAKE_INSTALL_PREFIX=../install
make install -j10
```

> **Note:** The GenFit install path is currently hardcoded in `CMakeLists.txt`. Update `GENFIT_INSTALL_DIR` to match your local GenFit installation before building.

## Running

```bash
source /cvmfs/sw-nightlies.hsf.org/key4hep/setup.sh
source build/displacedtrackingenv.sh

k4run options/runDisplacedTracking.py
```

Edit `options/runDisplacedTracking.py` to set the input file, geometry XML, and algorithm parameters before running.

## Key configuration parameters

| Parameter | Default | Description |
|---|---|---|
| `DetectorName` | `"Muon-System"` | Detector readout name |
| `MaxChi2` | `10.0` | Maximum chi2 for hit acceptance |
| `MaxDist` | `160.0 cm` | Maximum hit-to-hit distance in extension loop |
| `MinCosAngle2d` | `0.7` | Minimum cos(angle) between chord vectors for seeding |
| `MaxSeedPT` | `200.0 GeV` | Upper pT cut on triplet seeds |
| `ParticleType` | `"muon"` | Particle hypothesis for GenFit material effects |
| `UseGenFit` | `true` | Enable/disable GenFit Kalman fit |
| `MaxFitIterations` | `100` | GenFit Kalman filter max iterations |
| `DoInnerPropagation` | `false` | Enable extrapolation toward solenoid |
| `InnerPropTargetRadius` | `100.0 cm` | Target radius for inner propagation |

## Output collections

| Collection | Type | Contents |
|---|---|---|
| `DisplacedTracks` | `TrackCollection` | Reconstructed tracks with all states |
| `DisplacedTrackHits` | `TrackerHitPlaneCollection` | Hits associated to tracks |
| `DisplacedTrackHitSimLinks` | `TrackerHitSimTrackerHitLinkCollection` | Reco–sim hit links |

## Author

Mahmoud Al-Thakeel

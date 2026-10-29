# CombinatorialTrackFinder

A Gaudi algorithm for Key4hep that finds tracks in the planar hits of one tracking system
(an inner tracker, a silicon wrapper, a muon system, ...) without using any other detector.
It is built for displaced tracks, for example from long-lived particle decays, so it never
assumes that a track comes from the interaction point. The detector, field and cuts are set
in the steering file.

## How it works

1. The hits are grouped by layer (barrel and endcap layers, from the cellID).
2. Track candidates are found with one of two strategies (`Strategy`):
   - `Exhaustive`: every combination of k hits on k different layers is fitted (a circle
     in xy, a line of z versus arc length). Good for sparse events.
   - `SeedAndExtend`: seeds of 3 to 5 hits on neighbouring layers are extended layer by
     layer along their helix. Good for dense events. With `SeedAndExtendMode = BestFirst`,
     the best seeds go first, take their hits and are judged with a kink χ² that includes
     multiple scattering.
3. Road, outlier and quality cuts remove wrong combinations. The best candidate becomes a
   track, its hits are removed, and the search repeats.
4. Each track stores its fitted helix as EDM4hep track states (below). A GenFit fit can be
   added (`UseGenFit`).

## Track states

| Location | What it is | When |
|---|---|---|
| `AtFirstHit` | the fitted helix at the innermost hit (`referencePoint` = that hit) | always |
| `AtLastHit` | the same helix at the outermost hit | `StoreLastHitState = true` (default) |
| `AtIP` | the track at its closest approach to the origin | `IPStateMode` (below) |
| `AtOther` | GenFit fit at the first hit | `UseGenFit = true` |

All states follow the EDM4hep convention: the five parameters are relative to the
`referencePoint`, the track moves outward, and `omega > 0` means clockwise seen from +z.

`IPStateMode` says whether, and how, an `AtIP` state is made. The track is never forced to
come from the IP. The state only describes it at the point closest to the origin.

| `IPStateMode` | Use |
|---|---|
| `None` (default) | no `AtIP` state |
| `Helix` | the same helix at its closest approach to the origin. For tracks in a uniform field (inner trackers). |
| `SolenoidCrossing` | the helix is followed back to the solenoid boundary and matched there to a helix in the inner field. For tracks measured outside the solenoid (muon system). |

## Build

```bash
source /cvmfs/sw-nightlies.hsf.org/key4hep/setup.sh
cmake -S . -B build -DCMAKE_INSTALL_PREFIX=$PWD/install
cmake --build build -j8 --target install
(cd build && ctest)          # unit tests of the helix extrapolation and the kink chi2
```

GenFit is found on `CMAKE_PREFIX_PATH` or through `GENFIT_INSTALL_DIR`, or through
`GENFIT_INC` / `GENFIT_LIB` (see `cmake/FindGenFit.cmake`).

## Set up (in every new shell)

```bash
source /cvmfs/sw-nightlies.hsf.org/key4hep/setup.sh
export LD_LIBRARY_PATH=<path to>/install/lib64:$LD_LIBRARY_PATH
export PYTHONPATH=<path to>/install/python:$PYTHONPATH
```

## Run

```bash
k4run options/runGlobalTrackerTracking.py
```

| Steering | Detector |
|---|---|
| `runGlobalTrackerTracking.py` | IDEA vertex detector + silicon wrapper, merged with `MergeTrackerHitPlanes`; `IPStateMode = Helix` |
| `runMuonSystemTracking.py` | IDEA muon system, settings of the HNL muon-system studies; input, output and compact on the command line; `IPStateMode = SolenoidCrossing` |
| `runCombinatorialTrackFinder.py` | IDEA muon system, fixed example files; `IPStateMode = SolenoidCrossing` |

The geometry compact in the steering must be the one the input file was digitised with.
Every property is documented in `components/CombinatorialTrackFinder.h`.

## Files

| File | Content |
|---|---|
| `components/CombinatorialTrackFinder.{h,cpp}` | the algorithm |
| `components/MergeTrackerHitPlanes.cpp` | merges hit collections of several subdetectors |
| `src/`, `include/` | circle and helix tools, layer geometry, seed finder, kink χ², DD4hep field and material for GenFit |
| `test/` | unit tests |
| `options/` | example steering files |

## Author

Mahmoud Al-Thakeel

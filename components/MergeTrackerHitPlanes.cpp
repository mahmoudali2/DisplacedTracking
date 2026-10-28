/*
 * Merge several TrackerHitPlane collections (and their reco->sim link collections)
 * into one, so a single tracking pass can do pattern recognition AND fitting over
 * every subdetector sharing a readout at once.
 *
 * Motivation: GlobalTrackerReadoutID is shared by four DetElements — VertexBarrel,
 * VertexEndcap, SiWrB and SiWrD — each of which is digitised into its OWN
 * collection (VTXBDigis, VTXDDigis, SiWrBDigis, SiWrDDigis). CombinatorialTrackFinder
 * consumes a single TrackerHitPlaneCollection, so without merging, a track that
 * crosses the vertex barrel and then the silicon wrapper could never be found as
 * one track: each subdetector would have to be tracked separately and the segments
 * stitched afterwards.
 *
 * The outputs are podio SUBSET collections: they reference the original hits rather
 * than copying them, so every hit keeps its original ObjectID (collectionID, index).
 * That matters because the reco->sim links point at the original hits — copies would
 * silently break the truth association. It also means CombinatorialTrackFinder must
 * key its reco->sim map on (collectionID, index) and not on index alone, since the
 * per-collection indices all start at 0; that is done in CombinatorialTrackFinder.cpp.
 */
#include "edm4hep/TrackerHitPlaneCollection.h"
#include "edm4hep/TrackerHitSimTrackerHitLinkCollection.h"

#include "k4FWCore/Transformer.h"

#include <string>
#include <tuple>
#include <vector>

struct MergeTrackerHitPlanes final
    : k4FWCore::MultiTransformer<
          std::tuple<edm4hep::TrackerHitPlaneCollection,
                     edm4hep::TrackerHitSimTrackerHitLinkCollection>(
              const std::vector<const edm4hep::TrackerHitPlaneCollection*>&,
              const std::vector<const edm4hep::TrackerHitSimTrackerHitLinkCollection*>&)> {

  MergeTrackerHitPlanes(const std::string& name, ISvcLocator* svcLoc)
      : MultiTransformer(
            name, svcLoc,
            {KeyValues{"InputHitCollections", {"VTXBDigis"}},
             KeyValues{"InputRecoSimLinkCollections", {"VTXBSimDigiLinks"}}},
            {KeyValues{"OutputHitCollection", {"MergedTrackerHits"}},
             KeyValues{"OutputRecoSimLinkCollection", {"MergedTrackerHitSimLinks"}}}) {}

  std::tuple<edm4hep::TrackerHitPlaneCollection, edm4hep::TrackerHitSimTrackerHitLinkCollection>
  operator()(const std::vector<const edm4hep::TrackerHitPlaneCollection*>& hitColls,
             const std::vector<const edm4hep::TrackerHitSimTrackerHitLinkCollection*>& linkColls)
      const override {

    edm4hep::TrackerHitPlaneCollection hits;
    hits.setSubsetCollection(true);
    for (const auto* coll : hitColls) {
      if (!coll) continue;
      for (const auto hit : *coll) hits.push_back(hit);
    }

    edm4hep::TrackerHitSimTrackerHitLinkCollection links;
    links.setSubsetCollection(true);
    for (const auto* coll : linkColls) {
      if (!coll) continue;
      for (const auto link : *coll) links.push_back(link);
    }

    debug() << "Merged " << hits.size() << " hits from " << hitColls.size()
            << " collections and " << links.size() << " links from " << linkColls.size()
            << " collections" << endmsg;
    return std::make_tuple(std::move(hits), std::move(links));
  }
};

DECLARE_COMPONENT(MergeTrackerHitPlanes)

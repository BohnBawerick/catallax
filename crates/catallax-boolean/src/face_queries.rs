//! Repeated exact face membership without rebuilding a surface section per point.

use crate::{face_contains_with, region};
use catallax_numbers::{Decider, Refusal};
use catallax_shapes::Point3;
use catallax_topology::{FaceId, Model};
use std::collections::BTreeMap;

/// A batch of exact trimmed-face queries against one append-only model.
/// Boundary traces, projected enclosures and spherical reference windings are
/// reused. Queries first use certified chart homotopies, retaining exact surface
/// probes when no chart settles the winding. Keep the batch within one Decider,
/// numeric policy and operation-budget scope; drop it before leaving that scope.
/// Existing geometry and topology must remain unchanged until this value drops.
pub struct PreparedFaceQueries<'a> {
    model: &'a mut Model,
    regions: region::RegionCache,
    targets: BTreeMap<FaceId, Vec<i64>>,
}
impl<'a> PreparedFaceQueries<'a> {
    /// Borrows the model for a batch; preparing an empty batch does no geometry work.
    pub fn new(model: &'a mut Model) -> Self {
        Self {
            model,
            regions: region::RegionCache::default(),
            targets: BTreeMap::new(),
        }
    }

    /// Reads the model whose unchanged boundary definitions this batch retains.
    pub fn model(&self) -> &Model {
        self.model
    }

    /// Runs a query that may append auxiliary geometry, but must not alter any
    /// existing geometry or topology. Cached boundary proofs remain valid only
    /// under that restriction. Drop the batch before editing a face or edge.
    pub fn with_geometry<T>(&mut self, query: impl FnOnce(&mut Model) -> T) -> T {
        query(self.model)
    }

    /// Tests exact support and closed trimmed-face membership. A missing proof
    /// returns its refusal; it never becomes membership or a larger allowance.
    pub fn contains<P: Clone + Into<Point3>>(
        &mut self,
        face: FaceId,
        point: &P,
        dec: &mut Decider,
    ) -> Result<bool, Refusal> {
        let regions = &mut self.regions;
        let targets = &mut self.targets;
        face_contains_with(
            self.model,
            face,
            point,
            dec,
            |model, face, point, dec| match regions
                .locate_prepared(model, face, point, dec, targets)
            {
                Err(error) if error.code == "boolean.classify-ambiguous" => {
                    regions.locate(model, face, point, dec)
                }
                result => result,
            },
        )
    }
}

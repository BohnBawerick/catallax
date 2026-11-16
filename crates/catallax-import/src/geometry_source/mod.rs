//! Retained STEP geometry for certified and tolerant import views.
//!
//! Source numbers keep their written units. Reading a declaration does not
//! establish incidence, a material region or an approximation tolerance.

mod decode;
mod definitions;
pub use definitions::*;

use crate::ImportRefusal;
use crate::part21::{Entity, Graph};
use catallax_numbers::Decider;

/// One immutable Part 21 graph shared by geometry readers and the certified view.
#[derive(Clone, Debug)]
pub struct ParsedStep {
    graph: Graph,
}

impl ParsedStep {
    /// Read the Part 21 graph, charging source bytes before parsing.
    pub fn read(bytes: &[u8], dec: &mut Decider) -> Result<Self, ImportRefusal> {
        dec.charge_work(bytes.len() as u64)
            .map_err(|_| ImportRefusal::ResourceLimit)?;
        dec.check_deadline()
            .map_err(|_| ImportRefusal::ResourceLimit)?;
        let graph = crate::part21::parse(bytes)?;
        dec.check_deadline()
            .map_err(|_| ImportRefusal::ResourceLimit)?;
        Ok(Self { graph })
    }

    /// Original entity numbers in ascending order.
    pub fn entities(&self) -> impl Iterator<Item = EntityId> + '_ {
        self.graph.keys().copied()
    }

    /// Source inheritance components, retaining complex instances.
    pub fn entity_types(&self, id: EntityId) -> Result<Vec<&str>, ImportRefusal> {
        Ok(self
            .entity(id)?
            .parts
            .iter()
            .map(|(name, _)| name.as_str())
            .collect())
    }

    /// Read one geometry declaration without replacing its wrappers or units.
    pub fn geometry(&self, id: EntityId, dec: &mut Decider) -> Result<Geometry, ImportRefusal> {
        dec.charge_work(1)
            .map_err(|_| ImportRefusal::ResourceLimit)?;
        decode::geometry(self, id, dec)
    }

    /// Existing exact-source recognition over the same parsed graph.
    pub fn certified_view(&self, dec: &Decider) -> Result<crate::StepFile, ImportRefusal> {
        crate::recognize::read_graph(&self.graph, dec)
    }

    /// Source representation items and their context reference.
    pub fn representation(&self, id: EntityId) -> Result<Representation, ImportRefusal> {
        let e = self.entity(id)?;
        let args = e
            .part("REPRESENTATION")
            .or_else(|| {
                let [(kind, args)] = e.parts.as_slice() else {
                    return None;
                };
                (kind.ends_with("REPRESENTATION")).then_some(args.as_slice())
            })
            .ok_or_else(|| malformed(id, "expected a representation"))?;
        Ok(Representation {
            items: decode::references(args, 1, id)?,
            context: decode::reference(args, 2, id)?,
        })
    }

    /// Context dimensions, units and uncertainties without a global-unit guess.
    pub fn representation_context(
        &self,
        id: EntityId,
    ) -> Result<RepresentationContext, ImportRefusal> {
        let e = self.entity(id)?;
        let dimension = e
            .part("GEOMETRIC_REPRESENTATION_CONTEXT")
            .map(|a| {
                let index = usize::from(a.len() != 1) * 2;
                decode::integer(a, index, id)
            })
            .transpose()?;
        let units = e
            .part("GLOBAL_UNIT_ASSIGNED_CONTEXT")
            .map(|a| decode::references(a, 0, id))
            .transpose()?
            .unwrap_or_default();
        let uncertainties = e
            .part("GLOBAL_UNCERTAINTY_ASSIGNED_CONTEXT")
            .map(|a| decode::references(a, 0, id))
            .transpose()?
            .unwrap_or_default();
        Ok(RepresentationContext {
            dimension,
            units,
            uncertainties,
        })
    }

    fn entity(&self, id: EntityId) -> Result<&Entity, ImportRefusal> {
        self.graph
            .get(&id)
            .ok_or_else(|| malformed(id, "missing referenced entity"))
    }
}

fn malformed(id: EntityId, detail: impl std::fmt::Display) -> ImportRefusal {
    ImportRefusal::Malformed {
        detail: format!("#{id}: {detail}"),
    }
}

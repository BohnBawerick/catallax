//! Retained STEP geometry for certified and tolerant import views.
//!
//! Source numbers keep their written units. Reading a declaration does not
//! establish incidence, a material region or an approximation tolerance.

pub(crate) mod adapter;
mod decode;
mod definitions;
mod expression;
mod fitting;
mod mapping;
mod parametric;
mod topology;
mod units;
pub use definitions::*;
pub use fitting::{FitOptions, FittedCurve, FittedSurface};
pub use mapping::{Frame, MappingError};
pub use parametric::{AngleUnit, CoordinateSpace, ParameterChart};
pub use units::ContextUnits;

use crate::ImportRefusal;
use crate::part21::{Entity, Graph};
use catallax_numbers::Decider;
use std::sync::Arc;

/// One immutable Part 21 graph shared by geometry readers and the certified view.
#[derive(Clone, PartialEq, Eq, Debug)]
pub struct ParsedStep {
    graph: Arc<Graph>,
    source_sha256: String,
}

impl ParsedStep {
    /// Read the Part 21 graph, charging source bytes before parsing.
    pub fn read(bytes: &[u8], dec: &mut Decider) -> Result<Self, ImportRefusal> {
        dec.charge_work(bytes.len() as u64)
            .map_err(|_| ImportRefusal::ResourceLimit)?;
        dec.check_deadline()
            .map_err(|_| ImportRefusal::ResourceLimit)?;
        let graph = Arc::new(crate::part21::parse(bytes)?);
        dec.check_deadline()
            .map_err(|_| ImportRefusal::ResourceLimit)?;
        Ok(Self {
            graph,
            source_sha256: crate::sha256_hex(bytes),
        })
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

    /// Read topology declarations without inferring a material region.
    pub fn topology(&self, id: EntityId, dec: &mut Decider) -> Result<Topology, ImportRefusal> {
        let e = self.entity(id)?;
        decode::charge(e, dec)?;
        topology::read(e, id)
    }

    /// Read representation maps and transformation relationships.
    pub fn placement_link(
        &self,
        id: EntityId,
        dec: &mut Decider,
    ) -> Result<PlacementLink, ImportRefusal> {
        let e = self.entity(id)?;
        decode::charge(e, dec)?;
        topology::placement_link(e, id)
    }

    /// Model declarations in source order. Assembly representations can refer
    /// to these as prototypes; resolve mapping links before placing bodies.
    pub fn model_declarations(
        &self,
        dec: &mut Decider,
    ) -> Result<Vec<ModelDeclaration>, ImportRefusal> {
        let mut out = Vec::new();
        for (&id, e) in self.graph.iter() {
            dec.charge_work(1)
                .map_err(|_| ImportRefusal::ResourceLimit)?;
            let kind = if e.part("MANIFOLD_SOLID_BREP").is_some()
                || e.part("FACETED_BREP").is_some()
                || e.part("BREP_WITH_VOIDS").is_some()
            {
                Some(ModelDeclarationKind::Solid)
            } else if e.part("SHELL_BASED_SURFACE_MODEL").is_some() {
                Some(ModelDeclarationKind::Shells)
            } else if e.part("FACE_BASED_SURFACE_MODEL").is_some() {
                Some(ModelDeclarationKind::Faces)
            } else if e.part("GEOMETRIC_SET").is_some() || e.part("GEOMETRIC_CURVE_SET").is_some() {
                Some(ModelDeclarationKind::Geometry)
            } else {
                None
            };
            if let Some(kind) = kind {
                out.push(ModelDeclaration { entity: id, kind });
            }
        }
        Ok(out)
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
        if !e.parts.iter().any(|(kind, _)| {
            matches!(
                kind.as_str(),
                "REPRESENTATION_CONTEXT"
                    | "PARAMETRIC_REPRESENTATION_CONTEXT"
                    | "GEOMETRIC_REPRESENTATION_CONTEXT"
                    | "GLOBAL_UNIT_ASSIGNED_CONTEXT"
                    | "GLOBAL_UNCERTAINTY_ASSIGNED_CONTEXT"
            )
        }) {
            return Err(malformed(id, "expected a representation context"));
        }
        let dimension = e
            .part("GEOMETRIC_REPRESENTATION_CONTEXT")
            .map(|a| decode::integer(a, context_field(a.len(), id)?, id))
            .transpose()?;
        let units = e
            .part("GLOBAL_UNIT_ASSIGNED_CONTEXT")
            .map(|a| decode::references(a, context_field(a.len(), id)?, id))
            .transpose()?
            .unwrap_or_default();
        let uncertainties = e
            .part("GLOBAL_UNCERTAINTY_ASSIGNED_CONTEXT")
            .map(|a| decode::references(a, context_field(a.len(), id)?, id))
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

fn context_field(count: usize, id: EntityId) -> Result<usize, ImportRefusal> {
    match count {
        1 => Ok(0),
        3 => Ok(2),
        _ => Err(malformed(
            id,
            "unrecognized representation context field layout",
        )),
    }
}

/// An unmapped source carrier, retaining its complete graph and length units.
/// It supplies no native incidence or material certificate.
#[derive(Clone, PartialEq, Eq, Debug)]
pub struct RetainedGeometry {
    source: ParsedStep,
    entity: EntityId,
    length_mm: catallax_numbers::Rational,
}
impl RetainedGeometry {
    pub(crate) fn identity(&self) -> String {
        format!(
            "{} #{} {}",
            self.source.source_sha256, self.entity, self.length_mm
        )
    }
    /// Immutable original declarations, including wrappers and associated geometry.
    pub fn source(&self) -> &ParsedStep {
        &self.source
    }
    /// Source instance number within this graph.
    pub fn entity(&self) -> EntityId {
        self.entity
    }
    /// Millimetres per written Cartesian length unit.
    pub fn length_mm(&self) -> &catallax_numbers::Rational {
        &self.length_mm
    }
    pub(crate) fn admission_refusal(&self) -> ImportRefusal {
        ImportRefusal::NotSolid {
            detail: format!(
                "#{}: import.geometry-mapping-required: retained source carrier has no admitted native mapping",
                self.entity
            ),
        }
    }
}

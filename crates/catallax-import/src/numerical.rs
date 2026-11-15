//! A separate numerical result for source-incidence repair. Exact bodies never
//! receive a numerical cast, and numerical bodies expose no exact model.

use std::collections::BTreeMap;

use catallax_numbers::{
    Decider, Rational, Refusal, Sign, operation_budget_active, with_operation_budget,
};
use catallax_recipe::StepId;
use catallax_shapes::{Curve, Line, Plane, RPoint3, RVector3, Surface};
use catallax_topology::{BodyId, Model, Sense, numerical as n};

use crate::{ImportRefusal, RepairRecord, Solid, Vec3};

/// A source refusal or a numerical admission refusal. Neither contains partial geometry.
#[derive(Clone, Debug)]
pub enum ImportError {
    /// Recognition, source limits, repair or certified construction refused.
    Source(ImportRefusal),
    /// Numerical definition, tolerance or checked admission refused.
    Numerical(Refusal),
}

impl From<ImportRefusal> for ImportError {
    fn from(value: ImportRefusal) -> Self {
        Self::Source(value)
    }
}
impl From<Refusal> for ImportError {
    fn from(value: Refusal) -> Self {
        Self::Numerical(value)
    }
}

/// An import admitted by the existing certified path.
#[derive(Clone, Debug)]
pub struct CertifiedImport {
    model: Model,
    bodies: Vec<BodyId>,
    source_sha256: String,
}
impl CertifiedImport {
    /// The existing certified model; body marks retain their existing meaning.
    pub fn model(&self) -> &Model {
        &self.model
    }
    /// Source bodies, in source order.
    pub fn bodies(&self) -> &[BodyId] {
        &self.bodies
    }
    /// Actual immutable source digest.
    pub fn source_sha256(&self) -> &str {
        &self.source_sha256
    }
}

/// Numerical source bodies and the actual source-repair record.
/// Source repair history is separate from the numerical builder's own repairs.
#[derive(Clone, Debug)]
pub struct NumericalImport {
    bodies: Vec<(String, n::NumericalBody)>,
    record: RepairRecord,
}
impl NumericalImport {
    /// Numerical bodies keyed by source solid. Zero residuals remain numerical.
    pub fn bodies(&self) -> &[(String, n::NumericalBody)] {
        &self.bodies
    }
    /// Complete source interpretation, source limits and source repair history.
    pub fn record(&self) -> &RepairRecord {
        &self.record
    }
}

/// The result category is explicit and must be handled by the consumer.
#[derive(Clone, Debug)]
pub enum StepImport {
    /// An input already handled by the certified import path.
    Certified(Box<CertifiedImport>),
    /// A separately checked numerical B-rep, with no exact-model extraction.
    Numerical(Box<NumericalImport>),
}

fn unsupported(detail: &str) -> ImportError {
    Refusal::new("import.numerical-unsupported-family", detail).into()
}

/// Try the existing certified path, then numerically admit the bounded axis-
/// aligned planar source-incidence family. The raw source stays immutable.
///
/// The numerical path retains the recorded producer's actual source bounds and
/// repairs. It constructs a separate N01/N02 definition and rechecks it at the
/// caller's fixed tolerance. Every returned fallback stays numerical, including
/// zero sampled residuals and derivatives with cheap exact construction proofs.
/// No result is inserted into an exact caller model. The active request allowance
/// applies; without one, the whole query shares the unchanged 20M default.
pub fn import_step(
    bytes: &[u8],
    expected_sha256: &str,
    step: &StepId,
    options: &n::NumericalOptions,
    dec: &mut Decider,
) -> Result<StepImport, ImportError> {
    if operation_budget_active() {
        return import_inner(bytes, expected_sha256, step, options, dec);
    }
    let (result, work) = with_operation_budget(n::DEFAULT_WORK, || {
        import_inner(bytes, expected_sha256, step, options, dec)
    });
    if work.exhausted {
        return Err(ImportRefusal::ResourceLimit.into());
    }
    result
}

fn import_inner(
    bytes: &[u8],
    expected_sha256: &str,
    step: &StepId,
    options: &n::NumericalOptions,
    dec: &mut Decider,
) -> Result<StepImport, ImportError> {
    let mut model = Model::new();
    match crate::import_step_bodies(&mut model, bytes, expected_sha256, step) {
        Ok(bodies) => {
            return Ok(StepImport::Certified(Box::new(CertifiedImport {
                model,
                bodies,
                source_sha256: crate::sha256_hex(bytes),
            })));
        }
        Err(ImportRefusal::NotSolid { detail }) if detail.contains("exactly") => {}
        Err(error) => return Err(error.into()),
    }
    if options.tolerance_mm.sign() != Sign::Positive
        || options.tolerance_mm > crate::gap::auto_cap()
    {
        return Err(Refusal::new(
            "import.numerical-tolerance",
            "source tolerance must be positive and at most 0.01 mm",
        )
        .into());
    }
    let imported = crate::repair::apply_file(
        &mut Model::new(),
        crate::read_step(bytes)?,
        &crate::sha256_hex(bytes),
        "source.step",
        &Default::default(),
        &[],
        Some(step),
        None,
    )?;
    if !imported.left_open.is_empty() {
        return Err(unsupported("the source repair left an open solid"));
    }
    for bound in &imported.record.face_tolerances {
        let upper = Rational::parse(&bound.deviation_upper_mm).expect("producer bound");
        if upper > options.tolerance_mm {
            return Err(Refusal::new(
                "import.numerical-source-limit",
                format!(
                    "source bound {} mm exceeds the requested {} mm",
                    bound.deviation_upper_mm, options.tolerance_mm
                ),
            )
            .into());
        }
    }
    let definition = imported.comparison.complete_definition().ok_or_else(|| {
        unsupported("the builder did not retain the complete repaired definition")
    })?;
    let mut bodies = Vec::new();
    for solid in &definition.solids {
        let definition = planar_definition(solid, dec)?;
        bodies.push((
            solid.name.clone(),
            n::build(definition, options.clone(), dec)?,
        ));
    }
    Ok(StepImport::Numerical(Box::new(NumericalImport {
        bodies,
        record: imported.record,
    })))
}

fn point(value: &Vec3) -> catallax_shapes::Point3 {
    RPoint3::new(value[0].clone(), value[1].clone(), value[2].clone()).to_real()
}
fn vector(value: &Vec3) -> catallax_shapes::Vector3 {
    RVector3::new(value[0].clone(), value[1].clone(), value[2].clone()).to_real()
}

fn vertex(
    vertices: &mut Vec<n::NumericalVertex>,
    indices: &mut BTreeMap<Vec3, usize>,
    source: &crate::Vertex,
    dec: &mut Decider,
) -> Result<usize, ImportError> {
    dec.charge_work(1)?;
    if let Some(index) = indices.get(&source.point) {
        return Ok(*index);
    }
    let index = vertices.len();
    indices.insert(source.point.clone(), index);
    vertices.push(n::NumericalVertex {
        key: source.name.clone(),
        position: source.point.clone(),
    });
    Ok(index)
}

fn planar_definition(
    solid: &Solid,
    dec: &mut Decider,
) -> Result<n::NumericalDefinition, ImportError> {
    let mut definition = n::NumericalDefinition {
        name: solid.name.clone(),
        vertices: Vec::new(),
        edges: Vec::new(),
        faces: Vec::new(),
    };
    let mut vertices = BTreeMap::new();
    let mut edges = BTreeMap::new();
    for face in &solid.faces {
        dec.charge_work(1)?;
        let crate::Surface::Plane(plane) = &face.surface else {
            return Err(unsupported("this numerical family requires planar faces"));
        };
        let axes: Vec<_> = (0..3).filter(|i| !plane.axis[*i].is_zero()).collect();
        if axes.len() != 1 {
            return Err(unsupported(
                "this numerical family requires axis-aligned repaired planes",
            ));
        }
        let axis = axes[0];
        let sign = if plane.axis[axis].sign() == Sign::Positive {
            Rational::one()
        } else {
            -Rational::one()
        };
        let (u, v) = ((axis + 1) % 3, (axis + 2) % 3);
        let mut normal = std::array::from_fn(|_| Rational::zero());
        normal[axis] = sign.clone();
        let mut reference = std::array::from_fn(|_| Rational::zero());
        reference[u] = Rational::one();
        let uv = |p: &Vec3, dec: &mut Decider| -> Result<_, ImportError> {
            if p[axis] != plane.origin[axis] {
                return Err(unsupported(
                    "the retained repaired boundary is not on its plane",
                ));
            }
            Ok([
                crate::planar_arithmetic::arithmetic(&p[u], &plane.origin[u], '-', dec)?,
                crate::planar_arithmetic::arithmetic(
                    &crate::planar_arithmetic::arithmetic(&p[v], &plane.origin[v], '-', dec)?,
                    &sign,
                    '*',
                    dec,
                )?,
            ])
        };
        let mut loops = Vec::new();
        let mut boundaries: Vec<_> = face.loops.iter().collect();
        boundaries.sort_by_key(|bound| !bound.outer);
        for bound in boundaries {
            if bound.vertex.is_some() {
                return Err(unsupported(
                    "vertex-loop faces are outside this numerical family",
                ));
            }
            for chain in crate::build::chains(&face.name, bound.orientation, &bound.edges)? {
                let mut trims = Vec::new();
                for (edge, forward) in chain {
                    if !matches!(edge.curve, crate::Curve::Line(..)) {
                        return Err(unsupported(
                            "this numerical family requires finite straight boundaries",
                        ));
                    }
                    let (a, b) = if forward {
                        (&edge.start, &edge.end)
                    } else {
                        (&edge.end, &edge.start)
                    };
                    let ai = vertex(&mut definition.vertices, &mut vertices, a, dec)?;
                    let bi = vertex(&mut definition.vertices, &mut vertices, b, dec)?;
                    let pair = [ai.min(bi), ai.max(bi)];
                    let index = if let Some(index) = edges.get(&pair) {
                        *index
                    } else {
                        let index = definition.edges.len();
                        let start = &definition.vertices[pair[0]].position;
                        let end = &definition.vertices[pair[1]].position;
                        let direction = crate::planar_arithmetic::sub(end, start, dec)?;
                        definition.edges.push(n::NumericalEdge {
                            key: edge.edge.clone(),
                            curve: Curve::Line(Line {
                                origin: point(start),
                                dir: vector(&direction),
                            }),
                            vertices: pair,
                            parameters: [Rational::zero(), Rational::one()],
                        });
                        edges.insert(pair, index);
                        index
                    };
                    trims.push(n::NumericalTrim {
                        edge: index,
                        sense: if ai == pair[0] {
                            Sense::Same
                        } else {
                            Sense::Reversed
                        },
                        uv: [uv(&a.point, dec)?, uv(&b.point, dec)?],
                    });
                }
                loops.push(trims);
            }
        }
        definition.faces.push(n::NumericalFace {
            key: face.name.clone(),
            surface: Surface::Plane(Plane {
                origin: point(&plane.origin),
                normal: vector(&normal),
                ref_dir: vector(&reference),
            }),
            sense: if face.same_sense {
                Sense::Same
            } else {
                Sense::Reversed
            },
            loops,
        });
    }
    Ok(definition)
}

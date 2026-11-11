//! Recognition of supported B-rep geometry as exact source data.
//!
//! Decimal text becomes an exact rational. Lengths are converted to millimetres
//! exactly, directions and angles are not scaled. Anything outside the subset is
//! refused with its entity name and STEP type.

use std::collections::BTreeSet;

use catallax_numbers::{Decider, Rational};
use catallax_shapes::spline_tools::{RationalSplineCurve, RationalSplineSurface};

use crate::ImportRefusal;
use crate::part21::{Entity, Graph, Param, parse};

#[path = "spline.rs"]
mod spline;

/// Source spline curve. STEP logical flags are declarations, not proof.
#[derive(Clone, PartialEq, Eq, Debug)]
pub struct StepSplineCurve {
    /// Exact Euclidean control data and separate positive weights.
    pub geometry: RationalSplineCurve,
    /// Declared closure; `None` is STEP `.U.`.
    pub closed: Option<bool>,
    /// Declared self-intersection; never accepted as a validity certificate.
    pub self_intersect: Option<bool>,
}

/// Source spline surface. Reading it does not certify an imported solid.
#[derive(Clone, PartialEq, Eq, Debug)]
pub struct StepSplineSurface {
    /// Exact u-major, v-minor data.
    pub geometry: RationalSplineSurface,
    /// Declared U and V closure, not certified periodic transitions.
    pub closed: [Option<bool>; 2],
    /// Declared self-intersection.
    pub self_intersect: Option<bool>,
}

/// An exact point or vector, in millimetres for points.
pub type Vec3 = [Rational; 3];

/// An exact right-handed placement with its reference direction projected.
#[derive(Clone, PartialEq, Eq, Debug)]
pub struct Placement {
    /// Origin.
    pub origin: Vec3,
    /// Axis direction, not normalised.
    pub axis: Vec3,
    /// Projected x direction, perpendicular to the axis and not normalised.
    pub ref_dir: Vec3,
}

/// A supported surface.
#[derive(Clone, PartialEq, Eq, Debug)]
pub enum Surface {
    /// `PLANE`; the normal is the placement axis.
    Plane(Placement),
    /// `CYLINDRICAL_SURFACE`.
    Cylinder(Placement, Rational),
    /// `CONICAL_SURFACE` with semi-angle exactly 45 degrees (slope 1). Other
    /// angles have an irrational slope and are refused.
    Cone45(Placement, Rational),
    /// `CONICAL_SURFACE` at any other semi-angle, with its radius at the
    /// placement. The importer builds it through the face's rims; see
    /// [`SemiAngle`].
    Cone(Placement, Rational, SemiAngle),
    /// `SPHERICAL_SURFACE`.
    Sphere(Placement, Rational),
    /// `TOROIDAL_SURFACE` with major and minor radius.
    Torus(Placement, Rational, Rational),
    /// Exact finite B-spline source data, without a solid certificate.
    BSpline(StepSplineSurface),
}

/// A cone semi-angle exactly as the file writes it.
#[derive(Clone, PartialEq, Eq, Debug)]
pub struct SemiAngle {
    /// The written value.
    pub value: Rational,
    /// Degrees when the file declares a degree unit, otherwise radians.
    pub degrees: bool,
}

/// A supported curve.
#[derive(Clone, PartialEq, Eq, Debug)]
pub enum Curve {
    /// `LINE` through a point along a direction.
    Line(Vec3, Vec3),
    /// `CIRCLE`.
    Circle(Placement, Rational),
    /// `ELLIPSE` with semi-axes along the reference direction and perpendicular to it.
    Ellipse(Placement, Rational, Rational),
    /// Exact finite B-spline source data.
    BSpline(StepSplineCurve),
}

/// A `VERTEX_POINT`.
#[derive(Clone, PartialEq, Eq, Debug)]
pub struct Vertex {
    /// Entity name, such as `"#7"`.
    pub name: String,
    /// Position in millimetres.
    pub point: Vec3,
}

/// An `ORIENTED_EDGE` over an `EDGE_CURVE`.
#[derive(Clone, PartialEq, Eq, Debug)]
pub struct OrientedEdge {
    /// Name of the `EDGE_CURVE`.
    pub edge: String,
    /// Start of the `EDGE_CURVE` along its own direction.
    pub start: Vertex,
    /// End of the `EDGE_CURVE` along its own direction.
    pub end: Vertex,
    /// Edge direction against the curve's natural direction.
    pub curve_same_sense: bool,
    /// The curve.
    pub curve: Curve,
    /// Direction of travel in the loop against the `EDGE_CURVE` direction.
    pub forward: bool,
}

/// A face boundary.
#[derive(Clone, PartialEq, Eq, Debug)]
pub struct FaceLoop {
    /// True for `FACE_OUTER_BOUND`.
    pub outer: bool,
    /// The bound's orientation against its loop.
    pub orientation: bool,
    /// The loop's oriented edges in order. Empty for a `VERTEX_LOOP`.
    pub edges: Vec<OrientedEdge>,
    /// The single vertex of a `VERTEX_LOOP`, which has no edges.
    pub vertex: Option<Vertex>,
}

/// An `ADVANCED_FACE`.
#[derive(Clone, PartialEq, Eq, Debug)]
pub struct Face {
    /// Entity name, the recipe face key.
    pub name: String,
    /// The surface.
    pub surface: Surface,
    /// Face normal against the surface normal.
    pub same_sense: bool,
    /// Boundaries in file order.
    pub loops: Vec<FaceLoop>,
}

/// One `MANIFOLD_SOLID_BREP` with its closed shell.
#[derive(Clone, PartialEq, Eq, Debug)]
pub struct Solid {
    /// Entity name of the solid.
    pub name: String,
    /// Faces of the shell in file order.
    pub faces: Vec<Face>,
}

/// What a file contains, in millimetres.
#[derive(Clone, PartialEq, Eq, Debug)]
pub struct StepFile {
    /// Exact millimetres per file length unit.
    pub mm_per_unit: Rational,
    /// The solids in file order.
    pub solids: Vec<Solid>,
}

/// Reads and recognizes a Part 21 file.
pub fn read_step(bytes: &[u8]) -> Result<StepFile, ImportRefusal> {
    read_step_with_decider(bytes, &Decider::default())
}

/// Reads exact source geometry using caller-owned cancellation and work policy.
/// Spline data do not acquire incidence or solid proofs through recognition.
pub fn read_step_with_decider(bytes: &[u8], dec: &Decider) -> Result<StepFile, ImportRefusal> {
    let graph = parse(bytes)?;
    Reader::new(&graph, dec)?.file()
}

fn unsupported(id: u64, e: &Entity) -> ImportRefusal {
    ImportRefusal::Unsupported {
        entity: format!("#{id}"),
        type_name: e.type_name().to_string(),
    }
}

fn malformed(detail: impl Into<String>) -> ImportRefusal {
    ImportRefusal::Malformed {
        detail: detail.into(),
    }
}

struct Reader<'a> {
    g: &'a Graph,
    scale: Rational,
    dec: &'a Decider,
}

impl<'a> Reader<'a> {
    fn new(g: &'a Graph, dec: &'a Decider) -> Result<Self, ImportRefusal> {
        Ok(Self {
            g,
            scale: length_scale(g)?,
            dec,
        })
    }

    fn entity(&self, id: u64) -> Result<&'a Entity, ImportRefusal> {
        self.dec.check_deadline().map_err(spline::source_error)?;
        self.g
            .get(&id)
            .ok_or_else(|| malformed(format!("reference to missing entity #{id}")))
    }

    /// Arguments of a simple entity of the given type, or an unsupported refusal.
    fn args(&self, id: u64, ty: &str) -> Result<&'a [Param], ImportRefusal> {
        let e = self.entity(id)?;
        match e.parts.as_slice() {
            [(n, args)] if n == ty => Ok(args),
            _ => Err(unsupported(id, e)),
        }
    }

    fn file(&self) -> Result<StepFile, ImportRefusal> {
        let mut solids = Vec::new();
        for (&id, e) in self.g {
            match e.type_name() {
                "MANIFOLD_SOLID_BREP" => solids.push(self.solid(id)?),
                "BREP_WITH_VOIDS"
                | "BREP_WITH_VOIDS_AND_BOUNDS"
                | "SHELL_BASED_SURFACE_MODEL"
                | "NEXT_ASSEMBLY_USAGE_OCCURRENCE"
                | "FACETED_BREP" => return Err(unsupported(id, e)),
                _ => {}
            }
        }
        // Multiple direct solids are supported, but occurrence placement is not.
        // Never omit unsupported solid models or treat mapped instances as
        // additional unplaced source solids. This guard is deliberately limited
        // to the new multiple-direct-solid route; single-solid admission is unchanged.
        if solids.len() > 1 {
            for (&id, entity) in self.g {
                if entity.parts.iter().any(|(name, _)| {
                    matches!(
                        name.as_str(),
                        "MAPPED_ITEM"
                            | "REPRESENTATION_MAP"
                            | "REPRESENTATION_RELATIONSHIP_WITH_TRANSFORMATION"
                            | "ITEM_DEFINED_TRANSFORMATION"
                    ) || unsupported_solid_model(name)
                }) {
                    return Err(unsupported(id, entity));
                }
            }
        }
        if solids.is_empty() {
            return Err(malformed("no MANIFOLD_SOLID_BREP in the file"));
        }
        Ok(StepFile {
            mm_per_unit: self.scale.clone(),
            solids,
        })
    }

    fn solid(&self, id: u64) -> Result<Solid, ImportRefusal> {
        let args = self.args(id, "MANIFOLD_SOLID_BREP")?;
        let shell = reference(args, 1)?;
        let faces = self.args(shell, "CLOSED_SHELL")?;
        let mut out = Vec::new();
        for face in ref_list(faces, 1)? {
            out.push(self.face(face)?);
        }
        Ok(Solid {
            name: format!("#{id}"),
            faces: out,
        })
    }

    fn face(&self, id: u64) -> Result<Face, ImportRefusal> {
        let a = self.args(id, "ADVANCED_FACE")?;
        let mut loops = Vec::new();
        for b in ref_list(a, 1)? {
            let e = self.entity(b)?;
            let outer = match e.type_name() {
                "FACE_OUTER_BOUND" => true,
                "FACE_BOUND" => false,
                _ => return Err(unsupported(b, e)),
            };
            let ba = e.parts[0].1.as_slice();
            let lp = reference(ba, 1)?;
            let orientation = boolean(ba, 2)?;
            let loop_entity = self.entity(lp)?;
            if loop_entity.type_name() == "VERTEX_LOOP" {
                loops.push(FaceLoop {
                    outer,
                    orientation,
                    edges: Vec::new(),
                    vertex: Some(self.vertex(reference(self.args(lp, "VERTEX_LOOP")?, 1)?)?),
                });
                continue;
            }
            let edges = ref_list(self.args(lp, "EDGE_LOOP")?, 1)?;
            loops.push(FaceLoop {
                outer,
                orientation,
                edges: edges
                    .into_iter()
                    .map(|o| self.oriented_edge(o))
                    .collect::<Result<_, _>>()?,
                vertex: None,
            });
        }
        Ok(Face {
            name: format!("#{id}"),
            surface: self.surface(reference(a, 2)?)?,
            same_sense: boolean(a, 3)?,
            loops,
        })
    }

    fn oriented_edge(&self, id: u64) -> Result<OrientedEdge, ImportRefusal> {
        let a = self.args(id, "ORIENTED_EDGE")?;
        let ec = reference(a, 3)?;
        let c = self.args(ec, "EDGE_CURVE")?;
        Ok(OrientedEdge {
            edge: format!("#{ec}"),
            start: self.vertex(reference(c, 1)?)?,
            end: self.vertex(reference(c, 2)?)?,
            curve_same_sense: boolean(c, 4)?,
            curve: self.curve(reference(c, 3)?)?,
            forward: boolean(a, 4)?,
        })
    }

    fn vertex(&self, id: u64) -> Result<Vertex, ImportRefusal> {
        let a = self.args(id, "VERTEX_POINT")?;
        Ok(Vertex {
            name: format!("#{id}"),
            point: self.point(reference(a, 1)?)?,
        })
    }

    fn point(&self, id: u64) -> Result<Vec3, ImportRefusal> {
        let a = self.args(id, "CARTESIAN_POINT")?;
        Ok(self.triple(a)?.map(|c| c * self.scale.clone()))
    }

    fn direction(&self, id: u64) -> Result<Vec3, ImportRefusal> {
        let a = self.args(id, "DIRECTION")?;
        let v = self.triple(a)?;
        if v.iter().all(Rational::is_zero) {
            return Err(malformed(format!("zero DIRECTION #{id}")));
        }
        Ok(v)
    }

    fn triple(&self, a: &[Param]) -> Result<Vec3, ImportRefusal> {
        match a.get(1) {
            Some(Param::List(v)) if v.len() == 3 => {
                let n = |p: &Param| number(p);
                Ok([n(&v[0])?, n(&v[1])?, n(&v[2])?])
            }
            _ => Err(malformed("coordinates must be a list of three numbers")),
        }
    }

    fn length(&self, a: &[Param], i: usize) -> Result<Rational, ImportRefusal> {
        let v = number(a.get(i).ok_or_else(|| malformed("missing length"))?)?;
        Ok(&v * &self.scale)
    }

    fn positive(&self, a: &[Param], i: usize) -> Result<Rational, ImportRefusal> {
        let v = self.length(a, i)?;
        if v.sign() != catallax_numbers::Sign::Positive {
            return Err(malformed("radius must be positive"));
        }
        Ok(v)
    }

    fn placement(&self, id: u64) -> Result<Placement, ImportRefusal> {
        let a = self.args(id, "AXIS2_PLACEMENT_3D")?;
        let origin = self.point(reference(a, 1)?)?;
        let axis = match optional_ref(a, 2)? {
            Some(ax) => self.direction(ax)?,
            None => [Rational::zero(), Rational::zero(), Rational::one()],
        };
        let ref_dir = match optional_ref(a, 3)? {
            Some(rf) => self.direction(rf)?,
            None if axis[1].is_zero() && axis[2].is_zero() => {
                [Rational::zero(), Rational::one(), Rational::zero()]
            }
            None => [Rational::one(), Rational::zero(), Rational::zero()],
        };
        // ISO 10303-42 first_proj_axis; Placement carries ratios, not unit vectors.
        let dot: Rational = (0..3).fold(Rational::zero(), |s, i| s + &axis[i] * &ref_dir[i]);
        let norm_sq: Rational = (0..3).fold(Rational::zero(), |s, i| s + &axis[i] * &axis[i]);
        let along = &dot / &norm_sq;
        let ref_dir: Vec3 = std::array::from_fn(|i| &ref_dir[i] - &(&along * &axis[i]));
        if ref_dir.iter().all(Rational::is_zero) {
            return Err(malformed(format!(
                "placement #{id} axis and reference direction are parallel"
            )));
        }
        Ok(Placement {
            origin,
            axis,
            ref_dir,
        })
    }

    fn surface(&self, id: u64) -> Result<Surface, ImportRefusal> {
        let e = self.entity(id)?;
        if let Some(surface) = self.spline_surface(id, e)? {
            return Ok(Surface::BSpline(surface));
        }
        let [(ty, a)] = e.parts.as_slice() else {
            return Err(unsupported(id, e));
        };
        let pl = || self.placement(reference(a, 1)?);
        Ok(match ty.as_str() {
            "PLANE" => Surface::Plane(pl()?),
            "CYLINDRICAL_SURFACE" => Surface::Cylinder(pl()?, self.positive(a, 2)?),
            "SPHERICAL_SURFACE" => Surface::Sphere(pl()?, self.positive(a, 2)?),
            "TOROIDAL_SURFACE" => {
                let (major, minor) = (self.positive(a, 2)?, self.positive(a, 3)?);
                if minor >= major {
                    return Err(unsupported(id, e));
                }
                Surface::Torus(pl()?, major, minor)
            }
            "CONICAL_SURFACE" => {
                if is_45_degrees(a.get(3), self.g) {
                    Surface::Cone45(pl()?, self.length(a, 2)?)
                } else {
                    let radius = self.length(a, 2)?;
                    let value = number(a.get(3).ok_or_else(|| malformed("missing semi-angle"))?)?;
                    if radius.sign() == catallax_numbers::Sign::Negative {
                        return Err(malformed("cone radius must not be negative"));
                    }
                    Surface::Cone(
                        pl()?,
                        radius,
                        SemiAngle {
                            value,
                            degrees: degree_context(self.g),
                        },
                    )
                }
            }
            _ => return Err(unsupported(id, e)),
        })
    }

    fn curve(&self, id: u64) -> Result<Curve, ImportRefusal> {
        let e = self.entity(id)?;
        if let Some(curve) = self.spline_curve(id, e)? {
            return Ok(Curve::BSpline(curve));
        }
        let [(ty, a)] = e.parts.as_slice() else {
            return Err(unsupported(id, e));
        };
        Ok(match ty.as_str() {
            // The 3D curve decides; the parameter-space curves are not read.
            "SURFACE_CURVE" | "SEAM_CURVE" => self.curve(reference(a, 1)?)?,
            "LINE" => {
                let origin = self.point(reference(a, 1)?)?;
                let v = self.args(reference(a, 2)?, "VECTOR")?;
                Curve::Line(origin, self.direction(reference(v, 1)?)?)
            }
            "CIRCLE" => Curve::Circle(self.placement(reference(a, 1)?)?, self.positive(a, 2)?),
            "ELLIPSE" => Curve::Ellipse(
                self.placement(reference(a, 1)?)?,
                self.positive(a, 2)?,
                self.positive(a, 3)?,
            ),
            _ => return Err(unsupported(id, e)),
        })
    }
}

fn reference(a: &[Param], i: usize) -> Result<u64, ImportRefusal> {
    match a.get(i) {
        Some(Param::Ref(r)) => Ok(*r),
        _ => Err(malformed(format!(
            "expected an entity reference at argument {i}"
        ))),
    }
}

fn optional_ref(a: &[Param], i: usize) -> Result<Option<u64>, ImportRefusal> {
    match a.get(i) {
        Some(Param::Unset) => Ok(None),
        _ => reference(a, i).map(Some),
    }
}

fn ref_list(a: &[Param], i: usize) -> Result<Vec<u64>, ImportRefusal> {
    match a.get(i) {
        Some(Param::List(l)) => l
            .iter()
            .map(|p| match p {
                Param::Ref(r) => Ok(*r),
                _ => Err(malformed("expected a list of entity references")),
            })
            .collect(),
        _ => Err(malformed(format!(
            "expected a reference list at argument {i}"
        ))),
    }
}

fn boolean(a: &[Param], i: usize) -> Result<bool, ImportRefusal> {
    match a.get(i) {
        Some(Param::Enum(e)) if e == "T" => Ok(true),
        Some(Param::Enum(e)) if e == "F" => Ok(false),
        _ => Err(malformed(format!("expected .T. or .F. at argument {i}"))),
    }
}

/// An exact decimal. `1.` and `.5`-free STEP reals are normalised first.
fn number(p: &Param) -> Result<Rational, ImportRefusal> {
    let Param::Number(text) = p else {
        return Err(malformed("expected a number"));
    };
    let fixed = text.replace(".E", ".0E").replace(".e", ".0e");
    let fixed = if fixed.ends_with('.') {
        format!("{fixed}0")
    } else {
        fixed
    };
    Rational::parse(&fixed).map_err(|_| malformed(format!("unreadable number `{text}`")))
}

/// Longest chain of conversion-based units followed before the file is refused. It
/// also stops a conversion that names itself, directly or through other units.
const MAX_CONVERSION_DEPTH: u32 = 8;

/// Millimetres per file length unit, from the one length unit the file declares.
fn length_scale(g: &Graph) -> Result<Rational, ImportRefusal> {
    let mut scales = BTreeSet::new();
    // Only units a representation context assigns count; a conversion's base unit does not.
    let assigned: BTreeSet<u64> = g
        .values()
        .filter_map(|e| e.part("GLOBAL_UNIT_ASSIGNED_CONTEXT"))
        .filter_map(|a| match a.first() {
            Some(Param::List(l)) => Some(l),
            _ => None,
        })
        .flatten()
        .filter_map(|p| match p {
            Param::Ref(r) => Some(*r),
            _ => None,
        })
        .collect();
    for (id, e) in assigned.iter().filter_map(|id| Some((id, g.get(id)?))) {
        if e.part("LENGTH_UNIT").is_some() {
            scales.insert(unit_millimetres(g, *id, 0)?);
        }
    }
    let mut it = scales.into_iter();
    match (it.next(), it.next()) {
        (Some(s), None) => Ok(s),
        (None, _) => Err(unknown_unit("")),
        _ => Err(unknown_unit("conflicting length units")),
    }
}

fn unknown_unit(found: &str) -> ImportRefusal {
    ImportRefusal::UnknownUnit {
        found: found.to_string(),
    }
}

/// Millimetres in one of the length unit `id`: an SI metre with any prefix, or a
/// conversion-based unit worth its length measure of its own base unit. The
/// factors stay exact; a chain deeper than `MAX_CONVERSION_DEPTH`, which includes a
/// cycle, is refused.
fn unit_millimetres(g: &Graph, id: u64, depth: u32) -> Result<Rational, ImportRefusal> {
    let Some(e) = g.get(&id) else {
        return Err(unknown_unit("length unit"));
    };
    if let Some(si) = e.part("SI_UNIT") {
        let prefix = match si.first() {
            Some(Param::Enum(p)) => p.as_str(),
            _ => "",
        };
        let unit = match si.get(1) {
            Some(Param::Enum(u)) => u.as_str(),
            _ => "",
        };
        if unit != "METRE" {
            return Err(unknown_unit(unit));
        }
        // Powers of ten in metres, so a bare metre is 10^3 millimetres.
        let exponent: i32 = match prefix {
            "EXA" => 18,
            "PETA" => 15,
            "TERA" => 12,
            "GIGA" => 9,
            "MEGA" => 6,
            "KILO" => 3,
            "HECTO" => 2,
            "DECA" => 1,
            "" => 0,
            "DECI" => -1,
            "CENTI" => -2,
            "MILLI" => -3,
            "MICRO" => -6,
            "NANO" => -9,
            "PICO" => -12,
            "FEMTO" => -15,
            "ATTO" => -18,
            other => return Err(unknown_unit(other)),
        };
        return Rational::from_i64(10)
            .pow(exponent + 3)
            .ok_or_else(|| unknown_unit(prefix));
    }
    let Some(c) = e.part("CONVERSION_BASED_UNIT") else {
        return Err(unknown_unit("length unit"));
    };
    let name = match c.first() {
        Some(Param::Str(s)) => s.to_ascii_uppercase(),
        _ => String::new(),
    };
    if depth >= MAX_CONVERSION_DEPTH {
        return Err(unknown_unit(&name));
    }
    let measure = match c.get(1) {
        Some(Param::Ref(r)) => g.get(r).and_then(|e| e.part("LENGTH_MEASURE_WITH_UNIT")),
        _ => None,
    }
    .ok_or_else(|| unknown_unit(&name))?;
    let value = match measure.first() {
        Some(Param::Typed(_, v)) => v.first().map(number).transpose()?,
        _ => None,
    }
    .filter(|v| *v > Rational::zero())
    .ok_or_else(|| unknown_unit(&name))?;
    let base = match measure.get(1) {
        Some(Param::Ref(r)) => *r,
        _ => return Err(unknown_unit(&name)),
    };
    Ok(&value * &unit_millimetres(g, base, depth + 1)?)
}

/// True when a semi-angle parameter is exactly 45 in a degree context, or pi/4 is
/// never written exactly, so only degree files qualify.
fn is_45_degrees(p: Option<&Param>, g: &Graph) -> bool {
    degree_context(g) && p.and_then(|p| number(p).ok()) == Some(Rational::from_i64(45))
}

/// Whether the file declares a degree plane-angle unit.
fn degree_context(g: &Graph) -> bool {
    g.values().any(|e| {
        e.part("CONVERSION_BASED_UNIT").is_some_and(
            |c| matches!(c.first(), Some(Param::Str(s)) if s.eq_ignore_ascii_case("DEGREE")),
        )
    })
}

// Unsupported solid_model branches from the merged AP schema, including concrete
// subtypes (simple Part 21 entities do not spell out their supertypes). Check all
// parts of complex entities too. MANIFOLD_SOLID_BREP is the sole admitted branch.
// Schema audit: https://www.steptools.com/stds/stp_aim/html/t_solid_model.html
fn unsupported_solid_model(name: &str) -> bool {
    matches!(
        name,
        "SOLID_MODEL"
            | "BREP_WITH_VOIDS"
            | "BREP_WITH_VOIDS_AND_BOUNDS"
            | "FACETED_BREP"
            | "COMPLEX_SHELLED_SOLID"
            | "CSG_PRIMITIVE_SOLID_2D"
            | "CSG_SOLID"
            | "CSG_SOLID_2D"
            | "DOUBLE_OFFSET_SHELLED_SOLID"
            | "EDGE_BLENDED_SOLID"
            | "EXTRUDED_AREA_SOLID"
            | "EXTRUDED_FACE_SOLID"
            | "EXTRUDED_FACE_SOLID_WITH_DRAFT_ANGLE"
            | "EXTRUDED_FACE_SOLID_WITH_MULTIPLE_DRAFT_ANGLES"
            | "EXTRUDED_FACE_SOLID_WITH_TRIM_CONDITIONS"
            | "MODIFIED_SOLID"
            | "MODIFIED_SOLID_WITH_PLACED_CONFIGURATION"
            | "REVOLVED_AREA_SOLID"
            | "REVOLVED_FACE_SOLID"
            | "REVOLVED_FACE_SOLID_WITH_TRIM_CONDITIONS"
            | "RULED_SURFACE_SWEPT_AREA_SOLID"
            | "SCULPTURED_SOLID"
            | "SHELLED_SOLID"
            | "SOLID_REPLICA"
            | "SOLID_WITH_ANGLE_BASED_CHAMFER"
            | "SOLID_WITH_CHAMFERED_EDGES"
            | "SOLID_WITH_CIRCULAR_PATTERN"
            | "SOLID_WITH_CIRCULAR_POCKET"
            | "SOLID_WITH_CIRCULAR_PROTRUSION"
            | "SOLID_WITH_CONICAL_BOTTOM_ROUND_HOLE"
            | "SOLID_WITH_CONSTANT_RADIUS_EDGE_BLEND"
            | "SOLID_WITH_CURVED_SLOT"
            | "SOLID_WITH_DEPRESSION"
            | "SOLID_WITH_DOUBLE_OFFSET_CHAMFER"
            | "SOLID_WITH_FLAT_BOTTOM_ROUND_HOLE"
            | "SOLID_WITH_GENERAL_POCKET"
            | "SOLID_WITH_GENERAL_PROTRUSION"
            | "SOLID_WITH_GROOVE"
            | "SOLID_WITH_HOLE"
            | "SOLID_WITH_INCOMPLETE_CIRCULAR_PATTERN"
            | "SOLID_WITH_INCOMPLETE_RECTANGULAR_PATTERN"
            | "SOLID_WITH_POCKET"
            | "SOLID_WITH_PROTRUSION"
            | "SOLID_WITH_RECTANGULAR_PATTERN"
            | "SOLID_WITH_RECTANGULAR_POCKET"
            | "SOLID_WITH_RECTANGULAR_PROTRUSION"
            | "SOLID_WITH_SHAPE_ELEMENT_PATTERN"
            | "SOLID_WITH_SINGLE_OFFSET_CHAMFER"
            | "SOLID_WITH_SLOT"
            | "SOLID_WITH_SPHERICAL_BOTTOM_ROUND_HOLE"
            | "SOLID_WITH_STEPPED_ROUND_HOLE"
            | "SOLID_WITH_STEPPED_ROUND_HOLE_AND_CONICAL_TRANSITIONS"
            | "SOLID_WITH_STRAIGHT_SLOT"
            | "SOLID_WITH_TEE_SECTION_SLOT"
            | "SOLID_WITH_THROUGH_DEPRESSION"
            | "SOLID_WITH_TRAPEZOIDAL_SECTION_SLOT"
            | "SOLID_WITH_VARIABLE_RADIUS_EDGE_BLEND"
            | "SURFACE_CURVE_SWEPT_AREA_SOLID"
            | "SWEPT_AREA_SOLID"
            | "SWEPT_DISK_SOLID"
            | "SWEPT_FACE_SOLID"
            | "THICKENED_FACE_SOLID"
            | "TRACK_BLENDED_SOLID"
            | "TRACK_BLENDED_SOLID_WITH_END_CONDITIONS"
    )
}

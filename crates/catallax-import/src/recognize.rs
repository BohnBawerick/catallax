//! Recognition of the supported analytic B-rep subset as exact data.
//!
//! Decimal text becomes an exact rational. Lengths are converted to millimetres
//! exactly, directions and angles are not scaled. Anything outside the subset is
//! refused with its entity name and STEP type.

use std::collections::BTreeSet;

use catallax_numbers::Rational;

use crate::ImportRefusal;
use crate::part21::{Entity, Graph, Param, parse};

/// An exact point or vector, in millimetres for points.
pub type Vec3 = [Rational; 3];

/// A right-handed placement as written in the file.
#[derive(Clone, PartialEq, Eq, Debug)]
pub struct Placement {
    /// Origin.
    pub origin: Vec3,
    /// Axis direction, not normalised.
    pub axis: Vec3,
    /// Reference direction, not normalised.
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
    /// `SPHERICAL_SURFACE`.
    Sphere(Placement, Rational),
    /// `TOROIDAL_SURFACE` with major and minor radius.
    Torus(Placement, Rational, Rational),
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
    /// The loop's oriented edges in order.
    pub edges: Vec<OrientedEdge>,
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
    let graph = parse(bytes)?;
    Reader::new(&graph)?.file()
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
}

impl<'a> Reader<'a> {
    fn new(g: &'a Graph) -> Result<Self, ImportRefusal> {
        Ok(Self {
            g,
            scale: length_scale(g)?,
        })
    }

    fn entity(&self, id: u64) -> Result<&'a Entity, ImportRefusal> {
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
            let edges = ref_list(self.args(lp, "EDGE_LOOP")?, 1)?;
            loops.push(FaceLoop {
                outer,
                orientation: boolean(ba, 2)?,
                edges: edges
                    .into_iter()
                    .map(|o| self.oriented_edge(o))
                    .collect::<Result<_, _>>()?,
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
        let axis = optional_ref(a, 2)?;
        let ref_dir = optional_ref(a, 3)?;
        let z = || [Rational::zero(), Rational::zero(), Rational::one()];
        let x = || [Rational::one(), Rational::zero(), Rational::zero()];
        let (axis, ref_dir) = match (axis, ref_dir) {
            (None, None) => (z(), x()),
            (Some(ax), Some(rf)) => (self.direction(ax)?, self.direction(rf)?),
            _ => {
                return Err(malformed(format!(
                    "placement #{id} gives only one of axis and reference direction"
                )));
            }
        };
        let dot: Rational = (0..3).fold(Rational::zero(), |s, i| s + &axis[i] * &ref_dir[i]);
        if !dot.is_zero() {
            return Err(malformed(format!(
                "placement #{id} axis and reference direction are not perpendicular"
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
                if !is_45_degrees(a.get(3), self.g) {
                    return Err(unsupported(id, e));
                }
                Surface::Cone45(pl()?, self.length(a, 2)?)
            }
            _ => return Err(unsupported(id, e)),
        })
    }

    fn curve(&self, id: u64) -> Result<Curve, ImportRefusal> {
        let e = self.entity(id)?;
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

/// Millimetres per file length unit, from the one length unit the file declares.
fn length_scale(g: &Graph) -> Result<Rational, ImportRefusal> {
    let unknown = |found: &str| ImportRefusal::UnknownUnit {
        found: found.to_string(),
    };
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
    for e in assigned.iter().filter_map(|id| g.get(id)) {
        if e.part("LENGTH_UNIT").is_none() {
            continue;
        }
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
                return Err(unknown(unit));
            }
            let mm = match prefix {
                "MILLI" => Rational::one(),
                "CENTI" => Rational::from_i64(10),
                "" => Rational::from_i64(1000),
                other => return Err(unknown(other)),
            };
            scales.insert(mm);
        } else if let Some(c) = e.part("CONVERSION_BASED_UNIT") {
            let name = match c.first() {
                Some(Param::Str(s)) => s.to_ascii_uppercase(),
                _ => String::new(),
            };
            if name != "INCH" {
                return Err(unknown(&name));
            }
            let measure = match c.get(1) {
                Some(Param::Ref(r)) => *r,
                _ => return Err(unknown("INCH")),
            };
            let m = g
                .get(&measure)
                .and_then(|e| e.part("LENGTH_MEASURE_WITH_UNIT"))
                .ok_or_else(|| unknown("INCH"))?;
            let value = match m.first() {
                Some(Param::Typed(_, v)) => v.first().map(number).transpose()?,
                _ => None,
            }
            .ok_or_else(|| unknown("INCH"))?;
            // The conversion's own unit must be millimetres.
            let base = match m.get(1) {
                Some(Param::Ref(r)) => g.get(r).and_then(|e| e.part("SI_UNIT")),
                _ => None,
            };
            let ok = base.is_some_and(|b| {
                matches!(b.first(), Some(Param::Enum(p)) if p == "MILLI")
                    && matches!(b.get(1), Some(Param::Enum(u)) if u == "METRE")
            });
            if !ok {
                return Err(unknown("INCH"));
            }
            scales.insert(value);
        } else {
            return Err(unknown("length unit"));
        }
    }
    let mut it = scales.into_iter();
    match (it.next(), it.next()) {
        (Some(s), None) => Ok(s),
        (None, _) => Err(unknown("")),
        _ => Err(unknown("conflicting length units")),
    }
}

/// True when a semi-angle parameter is exactly 45 in a degree context, or pi/4 is
/// never written exactly, so only degree files qualify.
fn is_45_degrees(p: Option<&Param>, g: &Graph) -> bool {
    let degrees = g.values().any(|e| {
        e.part("CONVERSION_BASED_UNIT").is_some_and(
            |c| matches!(c.first(), Some(Param::Str(s)) if s.eq_ignore_ascii_case("DEGREE")),
        )
    });
    degrees && p.and_then(|p| number(p).ok()) == Some(Rational::from_i64(45))
}

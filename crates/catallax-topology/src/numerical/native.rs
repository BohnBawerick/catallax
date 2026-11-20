//! Retained finite parameter curves and sampled surface composition.

use super::*;
use catallax_numbers::Real;
use catallax_shapes::BSplineCurve;
use std::collections::BTreeMap;

mod domain;
mod joins;
pub(super) use joins::samples as join_samples;

/// A native parameter curve bound to one stable finite coedge use.
#[derive(Clone, Debug, catallax_numbers::Transfer)]
pub struct NumericalNativeTrim {
    /// Stable face key.
    pub face: String,
    /// Loop index, with zero naming the outer boundary.
    pub boundary: usize,
    /// Coedge index within that loop.
    pub use_index: usize,
    /// Stable edge key at that address.
    pub edge: String,
    /// X/y are surface parameters at the 3D edge's native parameter. Z is unused.
    /// The full carrier may extend beyond the edge's finite interval.
    pub curve: BSplineCurve,
    /// What this curve was made from.
    #[transfer(skip)]
    pub tag: CopyTag,
}

pub(super) struct Index<'a> {
    curves: BTreeMap<(usize, usize, usize), &'a BSplineCurve>,
    tags: BTreeMap<(usize, usize, usize), &'a CopyTag>,
    faces: BTreeSet<usize>,
    sampled_carriers: &'a BTreeSet<String>,
    numerical_edges: BTreeSet<usize>,
}

impl<'a> Index<'a> {
    pub fn new(
        definition: &NumericalDefinition,
        bindings: &'a [NumericalNativeTrim],
        sampled_carriers: &'a BTreeSet<String>,
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        let mut index = Self {
            curves: BTreeMap::new(),
            tags: BTreeMap::new(),
            faces: BTreeSet::new(),
            sampled_carriers,
            numerical_edges: BTreeSet::new(),
        };
        for face in &definition.faces {
            if sampled_carriers.contains(&face.key) {
                for trim in face.loops.iter().flatten() {
                    dec.charge_work(1)?;
                    index.numerical_edges.insert(trim.edge);
                }
            }
        }
        if bindings.is_empty() {
            return Ok(index);
        }
        let mut faces = BTreeMap::new();
        for (i, face) in definition.faces.iter().enumerate() {
            dec.charge_work(1)?;
            faces.insert(face.key.as_str(), i);
        }
        for binding in bindings {
            dec.charge_work(1)?;
            let face = *faces
                .get(binding.face.as_str())
                .ok_or_else(|| bad("trim", "native trim names an unknown face"))?;
            let trim = definition.faces[face]
                .loops
                .get(binding.boundary)
                .and_then(|row| row.get(binding.use_index))
                .ok_or_else(|| bad("trim", "native trim has a stale loop/use address"))?;
            let edge = definition
                .edges
                .get(trim.edge)
                .ok_or_else(|| bad("trim", "native trim names an unknown edge"))?;
            if edge.key != binding.edge {
                return Err(bad("trim", "native trim address names a different edge"));
            }
            let [lo, hi] = binding.curve.domain();
            if edge.parameters.iter().any(|t| t < lo || t > hi) {
                return Err(bad(
                    "trim",
                    "native trim does not cover the complete finite edge interval",
                ));
            }
            if index
                .curves
                .insert((face, binding.boundary, binding.use_index), &binding.curve)
                .is_some()
            {
                return Err(bad("trim", "duplicate native trim binding"));
            }
            index
                .tags
                .insert((face, binding.boundary, binding.use_index), &binding.tag);
            index.faces.insert(face);
        }
        Ok(index)
    }
    pub fn get(&self, face: usize, boundary: usize, use_index: usize) -> Option<&'a BSplineCurve> {
        self.curves.get(&(face, boundary, use_index)).copied()
    }
    /// The tag of the native curve bound to this use, if one is.
    pub fn tag(&self, face: usize, boundary: usize, use_index: usize) -> Option<&'a CopyTag> {
        self.tags.get(&(face, boundary, use_index)).copied()
    }
    pub fn has_face(&self, face: usize) -> bool {
        self.faces.contains(&face)
    }
    pub fn numerical_pair(&self, edge: usize) -> bool {
        self.numerical_edges.contains(&edge)
    }
}

fn parameters(
    edge: &NumericalEdge,
    curve: &BSplineCurve,
    sense: Sense,
    dec: &mut Decider,
) -> Result<Vec<Rational>, Refusal> {
    let mut interval = edge.parameters.clone();
    interval.sort();
    let mut breaks = BTreeSet::from(interval.clone());
    let xyz_knots = match &edge.curve {
        Curve::BSpline(c) => c.knots(),
        _ => &[],
    };
    for t in curve.knots().iter().chain(xyz_knots) {
        dec.charge_work(1)?;
        if t > &interval[0] && t < &interval[1] {
            breaks.insert(t.clone());
        }
    }
    let breaks: Vec<_> = breaks.into_iter().collect();
    let mut points = Vec::new();
    for span in breaks.windows(2) {
        for i in 0..SUBDIVISIONS {
            dec.charge_work(1)?;
            let share = Rational::from_ratio(i, SUBDIVISIONS).unwrap();
            points.push(&span[0] + &(&share * &(&span[1] - &span[0])));
        }
    }
    points.push(interval[1].clone());
    let first = if sense == Sense::Same {
        &edge.parameters[0]
    } else {
        &edge.parameters[1]
    };
    if first == &interval[1] {
        points.reverse();
    }
    Ok(points)
}

fn uv(curve: &BSplineCurve, t: &Rational, dec: &mut Decider) -> Result<[Real; 2], Refusal> {
    let point = curve.point_at(&Real::from(t.clone()), dec)?;
    Ok([point.x, point.y])
}

fn midpoint(value: &Real, dec: &mut Decider) -> Result<Rational, Refusal> {
    if let Some(q) = value.as_rational() {
        Ok(q)
    } else {
        Ok(dec.enclose_budgeted_at(value, BITS)?.mid())
    }
}

pub(super) fn point(
    surface: &Surface,
    uv: &[Real; 2],
    dec: &mut Decider,
) -> Result<[Ball; 3], Refusal> {
    if let Surface::BSpline(SplineSupport::RationalPatch(patch)) = surface {
        let p = patch.point_at(uv, dec)?;
        [&p.x, &p.y, &p.z]
            .map(|q| dec.enclose_budgeted_at(q, BITS))
            .into_iter()
            .collect::<Result<Vec<_>, _>>()?
            .try_into()
            .map_err(|_| unreachable!())
    } else {
        let u = dec.enclose_budgeted_at(&uv[0], BITS)?;
        let v = dec.enclose_budgeted_at(&uv[1], BITS)?;
        surface.evaluate(&u, &v, dec)
    }
}

/// Where a native u/v curve's two ends land on the surface against where the
/// trim's declared end values land, in traversal order.
pub(super) fn endpoints(
    face: &NumericalFace,
    trim: &NumericalTrim,
    edge: &NumericalEdge,
    curve: &BSplineCurve,
    dec: &mut Decider,
) -> Result<Vec<NumericalSample>, Refusal> {
    let ends = if trim.sense == Sense::Same {
        edge.parameters.each_ref()
    } else {
        [&edge.parameters[1], &edge.parameters[0]]
    };
    let mut samples = Vec::with_capacity(2);
    for (slot, t) in ends.into_iter().enumerate() {
        dec.charge_work(1)?;
        let at = uv(curve, t, dec)?;
        let found = point(&face.surface, &at, dec)?;
        let declared = trim.uv[slot].clone().map(Real::from);
        let expected = point(&face.surface, &declared, dec)?;
        samples.push(NumericalSample {
            parameters: vec![t.clone()],
            distance_upper_mm: distance(&found, &expected),
        });
    }
    Ok(samples)
}

pub(super) fn check_loops(
    face_index: usize,
    face: &NumericalFace,
    edges: &[NumericalEdge],
    index: &Index<'_>,
    tolerance: &Rational,
    dec: &mut Decider,
) -> Result<(), Refusal> {
    if index.has_face(face_index) && matches!(face.surface, Surface::Sphere(_) | Surface::Torus(_))
    {
        domain::check(face_index, face, edges, index, dec)?;
        analytic::check_periodic_domain(face, edges, dec)?;
    }
    if index.has_face(face_index)
        && (matches!(face.surface, Surface::Cone(_))
            || matches!(&face.surface, Surface::BSpline(source)
                if !matches!(source, SplineSupport::RationalPatch(_) | SplineSupport::RationalAtlas(_))))
    {
        domain::check(face_index, face, edges, index, dec)?;
    }
    if index.sampled_carriers.contains(&face.key) {
        return check_numerical_loops(face_index, face, edges, index, dec);
    }
    if !index.has_face(face_index) {
        return analytic::check_loops(face, edges, tolerance, dec);
    }
    let analytic_family = matches!(face.surface, Surface::Cylinder(_))
        || (matches!(face.surface, Surface::Plane(_))
            && face.loops.iter().flatten().any(|trim| {
                edges
                    .get(trim.edge)
                    .is_some_and(|edge| matches!(edge.curve, Curve::Circle(_)))
            }));
    if analytic_family {
        analytic::check_loops(face, edges, tolerance, dec)?;
    }
    let mut sampled = face.clone();
    sampled.loops.clear();
    for (row, boundary) in face.loops.iter().enumerate() {
        let mut polygon = Vec::new();
        for (use_index, trim) in boundary.iter().enumerate() {
            let edge = edges
                .get(trim.edge)
                .ok_or_else(|| bad("trim", "trim names an unknown edge"))?;
            if let Some(curve) = index.get(face_index, row, use_index) {
                let ts = parameters(edge, curve, trim.sense, dec)?;
                let mut previous = trim.uv[0].clone();
                for t in &ts[1..] {
                    let at = if t == ts.last().unwrap() {
                        trim.uv[1].clone()
                    } else {
                        let at = uv(curve, t, dec)?;
                        [midpoint(&at[0], dec)?, midpoint(&at[1], dec)?]
                    };
                    polygon.push(NumericalTrim {
                        tag: crate::numerical::CopyTag::Independent,
                        edge: trim.edge,
                        sense: trim.sense,
                        uv: [previous, at.clone()],
                    });
                    previous = at;
                }
            } else if matches!(edge.curve, Curve::Circle(_) | Curve::BSpline(_))
                && matches!(face.surface, Surface::Plane(_))
            {
                let [a, b] = if trim.sense == Sense::Same {
                    edge.parameters.each_ref()
                } else {
                    [&edge.parameters[1], &edge.parameters[0]]
                };
                let mut previous = trim.uv[0].clone();
                for i in 1..=SUBDIVISIONS {
                    dec.charge_work(1)?;
                    let at = if i == SUBDIVISIONS {
                        trim.uv[1].clone()
                    } else {
                        let share = Rational::from_ratio(i, SUBDIVISIONS).unwrap();
                        let t = a + &(&share * &(b - a));
                        let p = crate::numerical::curve_at(&edge.curve, &t, dec)?
                            .ok_or_else(|| bad("unsupported", "edge has no parameter evaluator"))?;
                        curved::parameters(&face.surface, &p, dec)?
                    };
                    polygon.push(NumericalTrim {
                        tag: crate::numerical::CopyTag::Independent,
                        edge: trim.edge,
                        sense: trim.sense,
                        uv: [previous, at.clone()],
                    });
                    previous = at;
                }
            } else {
                polygon.push(trim.clone());
            }
        }
        sampled.loops.push(polygon);
    }
    super::check_loops(&sampled, analytic_family.then_some((face, edges)), dec)
}

// Face-local UV topology uses nine locations per coedge for a declared
// numerical carrier. Its edge uses are measured like every other use, under
// labelled check names.
fn check_numerical_loops(
    face_index: usize,
    face: &NumericalFace,
    edges: &[NumericalEdge],
    index: &Index<'_>,
    dec: &mut Decider,
) -> Result<(), Refusal> {
    let mut loops = Vec::with_capacity(face.loops.len());
    for (row, boundary) in face.loops.iter().enumerate() {
        let mut polygon = Vec::with_capacity(boundary.len() * SUBDIVISIONS as usize);
        for (use_index, trim) in boundary.iter().enumerate() {
            let edge = edges
                .get(trim.edge)
                .ok_or_else(|| bad("trim", "trim names an unknown edge"))?;
            let curve = index.get(face_index, row, use_index).ok_or_else(|| {
                bad(
                    "trim",
                    "a numerical atlas face needs a native UV curve on every edge",
                )
            })?;
            let [a, b] = if trim.sense == Sense::Same {
                edge.parameters.each_ref()
            } else {
                [&edge.parameters[1], &edge.parameters[0]]
            };
            let step = (b - a) / Rational::from_i64(SUBDIVISIONS);
            let mut previous = trim.uv[0].clone();
            for i in 1..=SUBDIVISIONS {
                dec.charge_work(1)?;
                let at = if i == SUBDIVISIONS {
                    trim.uv[1].clone()
                } else {
                    let t = a + &(&step * &Rational::from_i64(i));
                    let at = uv(curve, &t, dec)?;
                    [midpoint(&at[0], dec)?, midpoint(&at[1], dec)?]
                };
                polygon.push(NumericalTrim {
                    tag: crate::numerical::CopyTag::Independent,
                    edge: trim.edge,
                    sense: trim.sense,
                    uv: [previous, at.clone()],
                });
                previous = at;
            }
        }
        loops.push(polygon);
    }
    super::check_boundaries(face, &loops, None, dec)
}

/// Builds with retained native UV curves at the 3D edges' finite parameters.
/// Every binding is checked; continuations retain and recheck these definitions.
/// Empty bindings preserve the ordinary numerical build path.
/// Cone parameter-curve hulls must retain a finite angular domain and stay
/// strictly above the apex.
pub fn build_with_native_trims(
    definition: NumericalDefinition,
    native_trims: Vec<NumericalNativeTrim>,
    options: NumericalOptions,
    dec: &mut Decider,
) -> Result<NumericalBody, Refusal> {
    guarded(|| {
        build_inner(
            definition,
            native_trims,
            BTreeSet::new(),
            options,
            Vec::new(),
            dec,
        )
    })
}

/// Admits explicitly numerical atlas carriers with sampled face-local checks.
/// Their UV walks use nine locations per coedge; chart coverage uses a 9 by 9
/// grid. Their edge uses are measured like every other use, under labelled
/// check names and the unchanged tolerance.
/// Undeclared face-local checks retain their original admission path.
pub fn build_with_sampled_carriers(
    definition: NumericalDefinition,
    native_trims: Vec<NumericalNativeTrim>,
    sampled_carriers: BTreeSet<String>,
    options: NumericalOptions,
    dec: &mut Decider,
) -> Result<NumericalBody, Refusal> {
    guarded(|| {
        build_inner(
            definition,
            native_trims,
            sampled_carriers,
            options,
            Vec::new(),
            dec,
        )
    })
}

impl NumericalBody {
    /// Original native parameter curves, including curves retained through edits.
    pub fn source_native_trims(&self) -> &[NumericalNativeTrim] {
        &self.source_native_trims
    }
    /// Original carriers declared numerical, retained through continuations.
    pub fn source_sampled_carriers(&self) -> &BTreeSet<String> {
        &self.source_sampled_carriers
    }
    /// Current declared numerical carriers using sampled face-local checks.
    pub fn sampled_carriers(&self) -> &BTreeSet<String> {
        &self.sampled_carriers
    }
    /// Current native parameter curves, bound to stable finite coedge uses.
    pub fn native_trims(&self) -> &[NumericalNativeTrim] {
        &self.native_trims
    }
    /// Rechecks an edited definition and explicitly replaced native trim bindings.
    ///
    /// A derived copy tag of the edit is kept only where this body shows both
    /// copies unchanged; every other one is reset and its use measured. The
    /// fillets, which make new blend faces and state their tags as their
    /// constructor, publish through
    /// [`NumericalBody::rebuild_with_sampled_carriers`] after resetting every
    /// use they take over.
    pub fn rebuild_with_native_trims(
        &self,
        mut definition: NumericalDefinition,
        mut native_trims: Vec<NumericalNativeTrim>,
        options: NumericalOptions,
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        self.keep_unchanged_copy_tags(&mut definition, &mut native_trims);
        let sampled_carriers = self
            .sampled_carriers
            .iter()
            .filter(|key| {
                definition.faces.iter().any(|face| {
                    &face.key == *key
                        && matches!(
                            face.surface,
                            Surface::BSpline(SplineSupport::RationalAtlas(_))
                        )
                })
            })
            .cloned()
            .collect();
        self.rebuild_with_sampled_carriers(definition, native_trims, sampled_carriers, options, dec)
    }

    /// Rechecks an edit with an explicit current numerical-carrier declaration.
    /// Source declarations and earlier check records remain historical.
    ///
    /// Copy tags are taken as stated, as [`build`](super::build) takes them:
    /// this is the producer path, for code that builds the copies it tags.
    /// Such code resets every use it takes over from this body
    /// ([`NumericalDefinition::reset_copy_tags`]) and tags only what it made,
    /// and the audit of section 4.4 checks each producer's rule. An edit that
    /// does not make both copies goes through
    /// [`NumericalBody::rebuild_with_native_trims`], which keeps only tags it
    /// can show still hold.
    pub fn rebuild_with_sampled_carriers(
        &self,
        definition: NumericalDefinition,
        native_trims: Vec<NumericalNativeTrim>,
        sampled_carriers: BTreeSet<String>,
        options: NumericalOptions,
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        guarded(|| {
            let mut next = build_inner_kind(
                self.kind(),
                definition,
                native_trims,
                sampled_carriers,
                options,
                self.evidence.repairs().to_vec(),
                dec,
            )?;
            next.inherit(self);
            Ok(next)
        })
    }
}

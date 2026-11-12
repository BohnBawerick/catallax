//! Finite boxes of trimmed faces and of finite supports. They are consulted
//! only after an exact pair or curve query has refused, to prove that the
//! refused meeting cannot reach a retained face. A box never claims a meeting
//! and never adds geometry to the model, so a query that answers keeps its
//! path and its work.

use catallax_numbers::{Ball, Decider, Rational, Real, Refusal, Sign};
use catallax_shapes::{Curve, Point3, SplineSupport, Surface, SurfaceId, Vector3};
use catallax_topology::{BodyId, BodyKind, EdgeId, FaceId, Loop, Model};
use std::collections::BTreeMap;

pub(crate) type Box3 = [Ball; 3];

fn enclose(x: &Real, dec: &mut Decider) -> Result<Ball, Refusal> {
    if let Some(q) = x.rational_ref() {
        Ok(Ball::from_rational(q, 128))
    } else {
        dec.enclose_budgeted_at(x, 128)
    }
}

fn point_box(p: &Point3, dec: &mut Decider) -> Result<Box3, Refusal> {
    Ok([
        enclose(&p.x, dec)?,
        enclose(&p.y, dec)?,
        enclose(&p.z, dec)?,
    ])
}

fn include(bounds: &mut Option<Box3>, next: &Box3) {
    if let Some(bounds) = bounds {
        for axis in 0..3 {
            bounds[axis] = bounds[axis].hull(&next[axis]);
        }
    } else {
        *bounds = Some(next.clone());
    }
}

fn around(center: &Box3, radius: &[Ball; 3]) -> Box3 {
    std::array::from_fn(|i| (&center[i] - &radius[i]).hull(&(&center[i] + &radius[i])))
}

/// `sqrt(1 - d_i^2 / |d|^2)` for each coordinate: the largest share of a unit
/// vector perpendicular to `d` that coordinate `i` can carry.
fn perpendicular_widths(d: &Vector3, dec: &mut Decider) -> Result<[Ball; 3], Refusal> {
    let norm = d.norm_squared();
    let parts = [&d.x, &d.y, &d.z];
    let mut widths = Vec::with_capacity(3);
    for i in 0..3 {
        let rest = &parts[(i + 1) % 3].square() + &parts[(i + 2) % 3].square();
        let share = &rest / &norm;
        widths.push(if share.rational_ref().is_some_and(Rational::is_zero) {
            Ball::from_rational(&Rational::zero(), 128)
        } else {
            enclose(&share.sqrt(), dec)?
        });
    }
    Ok(widths.try_into().expect("three coordinates"))
}

/// The unit direction of `d`, enclosed.
fn unit(d: &Vector3, dec: &mut Decider) -> Result<Option<[Ball; 3]>, Refusal> {
    let norm = d.norm_squared();
    let length = if norm.rational_ref() == Some(&Rational::one()) {
        Real::one()
    } else {
        norm.sqrt()
    };
    let length = enclose(&length, dec)?;
    let parts = [&d.x, &d.y, &d.z];
    let mut result = Vec::with_capacity(3);
    for part in parts {
        let Some(value) = enclose(part, dec)?.checked_div(&length) else {
            return Ok(None);
        };
        result.push(value);
    }
    Ok(Some(result.try_into().expect("three coordinates")))
}

fn control_box<'a>(
    controls: impl Iterator<Item = &'a Point3>,
    dec: &mut Decider,
) -> Result<Option<Box3>, Refusal> {
    let mut bounds = None;
    for control in controls {
        include(&mut bounds, &point_box(control, dec)?);
    }
    Ok(bounds)
}

fn enclosed_box(controls: Vec<([Ball; 3], Ball)>) -> Result<Option<Box3>, Refusal> {
    let mut bounds = None;
    for (control, _) in &controls {
        include(&mut bounds, control);
    }
    Ok(bounds)
}

/// A box around the complete support, when the support is bounded.
fn support_box(surface: &Surface, dec: &mut Decider) -> Result<Option<Box3>, Refusal> {
    match surface {
        // These families represent their geometry as the finite carrier over
        // its knot domain, and validated positive weights keep that carrier
        // inside the hull of its controls. An extrusion band can reach past
        // its carrier's domain, so other families give no box here.
        Surface::BSpline(SplineSupport::RollingEnvelope(band)) => {
            enclosed_box(band.enclosed_controls(dec)?)
        }
        Surface::BSpline(SplineSupport::RollingSaddle(band)) => {
            enclosed_box(band.enclosed_controls(dec)?)
        }
        Surface::BSpline(SplineSupport::RollingSphere(band)) => {
            enclosed_box(band.enclosed_controls(dec)?)
        }
        Surface::BSpline(SplineSupport::RationalPatch(patch)) => {
            enclosed_box(patch.enclosed_controls(dec)?)
        }
        Surface::BSpline(SplineSupport::RationalAtlas(atlas)) => {
            control_box(atlas.source().controls().iter(), dec)
        }
        Surface::Sphere(s) => {
            let radius = enclose(&s.radius, dec)?;
            Ok(Some(around(
                &point_box(&s.center, dec)?,
                &std::array::from_fn(|_| radius.clone()),
            )))
        }
        Surface::Torus(t) => {
            // Each point is centre + major e + minor n, with e a unit vector
            // perpendicular to the axis and n any unit vector.
            let major = enclose(&t.major, dec)?;
            let minor = enclose(&t.minor, dec)?;
            let widths = perpendicular_widths(&t.axis, dec)?;
            let radius = std::array::from_fn(|i| &(&major * &widths[i]) + &minor);
            Ok(Some(around(&point_box(&t.center, dec)?, &radius)))
        }
        _ => Ok(None),
    }
}

/// A box around an edge's retained arc, or around its whole finite carrier.
fn edge_box(
    model: &Model,
    edge: EdgeId,
    cache: &mut TrimBounds,
    dec: &mut Decider,
) -> Result<Option<Box3>, Refusal> {
    let edge_id = edge;
    let edge = model.edge(edge);
    Ok(match model.geom().curve(edge.curve) {
        Curve::Line(_) => {
            let (Some(start), Some(end)) = (edge.start, edge.end) else {
                return Ok(None);
            };
            let mut bounds = None;
            for vertex in [start, end] {
                let p = model.geom().position(model.vertex(vertex).point);
                include(&mut bounds, &point_box(&p, dec)?);
            }
            bounds
        }
        Curve::Circle(c) => {
            if let Some(bounds) = short_arc_box(model, edge_id, c, dec)? {
                return Ok(Some(bounds));
            }
            let radius = enclose(&c.radius_sq, dec)?;
            let Some(radius) = radius.sqrt() else {
                return Ok(None);
            };
            let widths = perpendicular_widths(&c.normal, dec)?;
            let radius = std::array::from_fn(|i| &radius * &widths[i]);
            Some(around(&point_box(&c.center, dec)?, &radius))
        }
        Curve::Ellipse(e) => {
            let reach = enclose(&e.major_sq, dec)?.hull(&enclose(&e.minor_sq, dec)?);
            let Some(reach) = reach.sqrt() else {
                return Ok(None);
            };
            Some(around(
                &point_box(&e.center, dec)?,
                &std::array::from_fn(|_| reach.clone()),
            ))
        }
        // Positive weights keep the complete finite curve inside its controls.
        Curve::BSpline(curve) => control_box(curve.controls(), dec)?,
        // A section component lies on both supports.
        Curve::Intersection(curve) => {
            let mut bounds: Option<Box3> = None;
            for surface in curve.surfaces {
                if let Some(next) = cache.support(model, surface, dec)? {
                    bounds = Some(match bounds {
                        None => next,
                        Some(old) => {
                            let meet: Option<Vec<Ball>> =
                                (0..3).map(|i| meet(&old[i], &next[i])).collect();
                            let Some(meet) = meet else {
                                // Disjoint supports carry no component.
                                return Ok(None);
                            };
                            meet.try_into().expect("three coordinates")
                        }
                    });
                }
            }
            bounds
        }
        _ => None,
    })
}

/// A circular arc turning less than half a circle lies in the triangle of its
/// two ends and the meeting point of the tangents there.
fn short_arc_box(
    model: &Model,
    edge: EdgeId,
    circle: &catallax_shapes::Circle,
    dec: &mut Decider,
) -> Result<Option<Box3>, Refusal> {
    let (Some(from), Some(to)) = crate::geometry::natural_ends(model, edge) else {
        return Ok(None);
    };
    if from == to {
        return Ok(None);
    }
    let a = model.geom().position(from);
    let b = model.geom().position(to);
    let ra = a.minus(&circle.center);
    let rb = b.minus(&circle.center);
    // A circle runs counter-clockwise about its normal from `from` to `to`.
    let turn = circle.normal.dot(&ra.cross(&rb));
    if crate::geometry::sign(&turn, dec, "boolean.bound-ambiguous")? != Sign::Positive {
        return Ok(None);
    }
    let half = Real::from(Rational::from_ratio(1, 2).expect("fixed denominator"));
    let middle = ra.plus(&rb).scaled(&half);
    let corner = circle
        .center
        .plus(&middle.scaled(&(&circle.radius_sq / &middle.norm_squared())));
    let mut bounds = None;
    for p in [&a, &b, &corner] {
        include(&mut bounds, &point_box(p, dec)?);
    }
    Ok(bounds)
}

fn meet(a: &Ball, b: &Ball) -> Option<Ball> {
    let lo = a.lo().max(b.lo());
    let hi = a.hi().min(b.hi());
    (lo <= hi).then(|| Ball::from_rational(&lo, 128).hull(&Ball::from_rational(&hi, 128)))
}

/// A box around every boundary edge and vertex of a face.
fn boundary_box(
    model: &Model,
    face: FaceId,
    cache: &mut TrimBounds,
    dec: &mut Decider,
) -> Result<Option<Box3>, Refusal> {
    let mut boundary = None;
    for lp in &model.face(face).loops {
        match model.loop_(*lp) {
            Loop::Vertex(vertex) => {
                let p = model.geom().position(model.vertex(*vertex).point);
                include(&mut boundary, &point_box(&p, dec)?);
            }
            Loop::Coedges(coedges) => {
                for coedge in coedges {
                    let Some(next) = edge_box(model, model.coedge(*coedge).edge, cache, dec)?
                    else {
                        return Ok(None);
                    };
                    include(&mut boundary, &next);
                }
            }
        }
    }
    Ok(boundary)
}

/// A box around a torus face that keeps to one side of the axis.
///
/// When the boundary's projection onto the plane across the axis misses the
/// axis, every boundary point lies in one pointed sector about the axis. The
/// angle about the axis has no critical point on a ring torus, so the face's
/// own range of angles ends at boundary points. It is either inside that
/// sector or covers the opposite direction, where no boundary point lies; one
/// point of the opposite meridian decides which.
fn torus_piece_box(
    model: &Model,
    face: FaceId,
    torus: &catallax_shapes::Torus,
    cache: &mut TrimBounds,
    dec: &mut Decider,
) -> Result<Option<Box3>, Refusal> {
    let Some(boundary) = boundary_box(model, face, cache, dec)? else {
        return Ok(None);
    };
    let Some(direction) = unit(&torus.axis, dec)? else {
        return Ok(None);
    };
    let center = point_box(&torus.center, dec)?;
    let relative: Box3 = std::array::from_fn(|i| &boundary[i] - &center[i]);
    let mut along = Ball::from_rational(&Rational::zero(), 128);
    for i in 0..3 {
        along = &along + &(&relative[i] * &direction[i]);
    }
    let across: Box3 = std::array::from_fn(|i| &relative[i] - &(&along * &direction[i]));
    // Squared distance from the axis, from below and above.
    let (mut near, mut far) = (Rational::zero(), Rational::zero());
    for value in &across {
        let (lo, hi) = (value.lo(), value.hi());
        let low = if lo > Rational::zero() {
            lo.clone()
        } else if hi < Rational::zero() {
            -hi.clone()
        } else {
            Rational::zero()
        };
        let high = lo.abs().max(hi.abs());
        near = near + low.clone() * low;
        far = far + high.clone() * high;
    }
    if near.is_zero() {
        return Ok(None);
    }
    let (Some(near), Some(far)) = (
        Ball::from_rational(&near, 128).sqrt(),
        Ball::from_rational(&far, 128).sqrt(),
    ) else {
        return Ok(None);
    };
    let (near, far) = (near.lo(), far.hi());
    if near <= Rational::zero() {
        return Ok(None);
    }
    let one = Rational::one();
    // Each unit direction q/|q| in the sector, coordinate by coordinate.
    let unit_range: Vec<Ball> = across
        .iter()
        .map(|value| {
            let (lo, hi) = (value.lo(), value.hi());
            let lower = if lo > Rational::zero() {
                lo / far.clone()
            } else {
                (lo / near.clone()).max(-one.clone())
            };
            let upper = if hi < Rational::zero() {
                hi / far.clone()
            } else {
                (hi / near.clone()).min(one.clone())
            };
            Ball::from_rational(&lower, 128).hull(&Ball::from_rational(&upper, 128))
        })
        .collect();
    // One boundary vertex fixes a direction inside the sector.
    let Some(vertex) = model
        .face(face)
        .loops
        .iter()
        .find_map(|lp| match model.loop_(*lp) {
            Loop::Vertex(vertex) => Some(*vertex),
            Loop::Coedges(coedges) => coedges.iter().find_map(|c| {
                let edge = model.edge(model.coedge(*c).edge);
                edge.start.or(edge.end)
            }),
        })
    else {
        return Ok(None);
    };
    let v = model.geom().position(model.vertex(vertex).point);
    let offset = v.minus(&torus.center);
    let inside = offset.minus(
        &torus
            .axis
            .scaled(&(&offset.dot(&torus.axis) / &torus.axis.norm_squared())),
    );
    let reach = &torus.major + &torus.minor;
    let opposite = torus
        .center
        .plus(&inside.scaled(&(&(-reach.clone()) / &inside.norm_squared().sqrt())));
    let probe = cache.probe.get_or_insert_with(|| model.clone());
    let point = probe
        .geom_mut()
        .add_point(catallax_shapes::PointDef::Explicit(opposite));
    if crate::region::locate(probe, face, point, dec)? != crate::region::Location::Outside {
        return Ok(None);
    }
    let low = enclose(&(&torus.major - &torus.minor), dec)?;
    let high = enclose(&reach, dec)?;
    let radius = low.hull(&high);
    let minor = enclose(&torus.minor, dec)?;
    Ok(Some(std::array::from_fn(|i| {
        let height = &minor * &direction[i];
        let tube = (-&height).hull(&height);
        &(&center[i] + &(&radius * &unit_range[i])) + &tube
    })))
}

/// A box around a bounded face of a solid.
fn face_box(
    model: &Model,
    face: FaceId,
    cache: &mut TrimBounds,
    dec: &mut Decider,
) -> Result<Option<Box3>, Refusal> {
    let support = cache.support(model, model.face(face).surface, dec)?;
    let surface = model.geom().surface(model.face(face).surface);
    if let Surface::Torus(torus) = surface {
        let support = support.expect("a torus support has a box");
        return Ok(Some(
            match torus_piece_box(model, face, torus, cache, dec)? {
                Some(piece) => {
                    let meet: Option<Vec<Ball>> =
                        (0..3).map(|i| meet(&support[i], &piece[i])).collect();
                    let Some(meet) = meet else {
                        return Ok(Some(piece));
                    };
                    meet.try_into().expect("three coordinates")
                }
                None => support,
            },
        ));
    }
    if let Some(bounds) = support {
        return Ok(Some(bounds));
    }
    let band = match surface {
        Surface::Plane(_) => None,
        Surface::Cylinder(c) => Some((&c.origin, &c.axis, None)),
        Surface::Cone(c) => Some((&c.apex, &c.axis, Some(&c.slope))),
        _ => return Ok(None),
    };
    let Some(boundary) = boundary_box(model, face, cache, dec)? else {
        return Ok(None);
    };
    let Some((origin, axis, slope)) = band else {
        // A bounded planar region lies in the convex hull of its boundary.
        return Ok(Some(boundary));
    };
    // The axial coordinate has no critical point on a cylinder, nor on a cone
    // nappe away from its apex, so over the face it peaks on the boundary. The
    // apex is the nappe's lowest point and is included for a cone.
    let Some(direction) = unit(axis, dec)? else {
        return Ok(None);
    };
    let origin = point_box(origin, dec)?;
    let mut axial = Ball::from_rational(&Rational::zero(), 128);
    for i in 0..3 {
        axial = &axial + &(&(&boundary[i] - &origin[i]) * &direction[i]);
    }
    let (lo, hi) = match slope {
        None => (axial.lo(), axial.hi()),
        Some(_) => (Rational::zero(), axial.hi().max(Rational::zero())),
    };
    let reach = match slope {
        None => {
            let Surface::Cylinder(c) = surface else {
                unreachable!()
            };
            let Some(radius) = enclose(&c.radius_sq, dec)?.sqrt() else {
                return Ok(None);
            };
            radius
        }
        Some(slope) => &enclose(slope, dec)? * &Ball::from_rational(&hi, 128),
    };
    let widths = perpendicular_widths(axis, dec)?;
    let span = Ball::from_rational(&lo, 128).hull(&Ball::from_rational(&hi, 128));
    let center: Box3 = std::array::from_fn(|i| &origin[i] + &(&direction[i] * &span));
    let radius = std::array::from_fn(|i| &reach * &widths[i]);
    Ok(Some(around(&center, &radius)))
}

fn apart(a: &Box3, b: &Box3) -> bool {
    (0..3).any(|i| a[i].hi() < b[i].lo() || b[i].hi() < a[i].lo())
}

/// Whether a face is one of the fitted rolling bands a fillet makes.
pub(crate) fn fitted_band(model: &Model, face: FaceId) -> bool {
    matches!(
        model.geom().surface(model.face(face).surface),
        Surface::BSpline(
            SplineSupport::RollingEnvelope(_)
                | SplineSupport::RollingSaddle(_)
                | SplineSupport::RollingSphere(_)
        )
    )
}

/// A complete, already charged separation proof for these unchanged bodies.
#[derive(Clone, Copy)]
pub(crate) struct BodyApartProof {
    a: BodyId,
    b: BodyId,
}

impl BodyApartProof {
    pub(crate) fn matches(&self, a: BodyId, b: BodyId) -> bool {
        self.a == a && self.b == b
    }
}

/// Memoized face and support boxes for one section pass.
///
/// The proofs run only for operands with a fitted rolling band. Everywhere
/// else a refusal stands at once, so callers that catch it and try another
/// route spend exactly the work they did before.
#[derive(Clone, Default)]
pub(crate) struct TrimBounds {
    enabled: bool,
    faces: BTreeMap<FaceId, Option<Box3>>,
    supports: BTreeMap<SurfaceId, Option<Box3>>,
    // Optional witnesses stay in this copy, outside the caller's model.
    probe: Option<Model>,
    // Work-span observation only: names this cache's entries.
    #[cfg(feature = "work-span")]
    span: catallax_numbers::work_span::Instance,
}

// The observer's resources for this cache's two maps. Each map has its own
// kind and every entry its own index, so no two entries share a resource.
#[cfg(feature = "work-span")]
impl TrimBounds {
    fn span_support(&self, surface: SurfaceId) -> (&'static str, u64, u64) {
        (
            "trim bounds cache: support",
            self.span.id(),
            u64::from(surface.index()),
        )
    }

    fn span_face(&self, face: FaceId) -> (&'static str, u64, u64) {
        (
            "trim bounds cache: face",
            self.span.id(),
            u64::from(face.index()),
        )
    }
}

impl TrimBounds {
    /// Box proofs for these faces, enabled when one of them is a fitted band.
    pub(crate) fn for_faces<'a>(
        model: &Model,
        faces: impl IntoIterator<Item = &'a FaceId>,
    ) -> Self {
        Self {
            enabled: faces.into_iter().any(|face| fitted_band(model, *face)),
            ..Self::default()
        }
    }

    /// Box proofs for two bodies, enabled when either has a fitted band.
    pub(crate) fn for_bodies(model: &Model, a: BodyId, b: BodyId) -> Self {
        let faces = [a, b]
            .into_iter()
            .flat_map(|body| model.faces_of_body(body))
            .collect::<Vec<_>>();
        Self::for_faces(model, &faces)
    }

    fn support(
        &mut self,
        model: &Model,
        surface: SurfaceId,
        dec: &mut Decider,
    ) -> Result<Option<Box3>, Refusal> {
        #[cfg(feature = "work-span")]
        {
            let (kind, owner, index) = self.span_support(surface);
            catallax_numbers::work_span::read_at(kind, owner, index);
        }
        if let Some(bounds) = self.supports.get(&surface) {
            return Ok(bounds.clone());
        }
        let bounds = support_box(model.geom().surface(surface), dec)?;
        self.supports.insert(surface, bounds.clone());
        #[cfg(feature = "work-span")]
        {
            let (kind, owner, index) = self.span_support(surface);
            catallax_numbers::work_span::write_at(kind, owner, index);
        }
        Ok(bounds)
    }

    fn face(&mut self, model: &Model, face: FaceId, dec: &mut Decider) -> Option<Box3> {
        #[cfg(feature = "work-span")]
        {
            let (kind, owner, index) = self.span_face(face);
            catallax_numbers::work_span::read_at(kind, owner, index);
        }
        if let Some(bounds) = self.faces.get(&face) {
            return bounds.clone();
        }
        // An unsettled box is no proof. The refusal it would excuse stands.
        let bounds = face_box(model, face, self, dec).ok().flatten();
        self.faces.insert(face, bounds.clone());
        #[cfg(feature = "work-span")]
        {
            let (kind, owner, index) = self.span_face(face);
            catallax_numbers::work_span::write_at(kind, owner, index);
        }
        bounds
    }

    /// Whether the two bounded faces are strictly apart.
    pub(crate) fn separated(
        &mut self,
        model: &Model,
        a: FaceId,
        b: FaceId,
        dec: &mut Decider,
    ) -> bool {
        if !self.enabled {
            return false;
        }
        let boxes = [self.face(model, a, dec), self.face(model, b, dec)];
        if let [Some(x), Some(y)] = &boxes
            && apart(x, y)
        {
            return true;
        }
        // A face's box on which the other support's equation keeps one
        // strict sign holds no point of that support.
        for (bounded, other) in [(&boxes[0], b), (&boxes[1], a)] {
            if let Some(bounds) = bounded
                && let Ok(Some(value)) = crate::bounds::support_box(
                    model.geom().surface(model.face(other).surface),
                    bounds,
                    dec,
                )
                && !value.contains_zero()
            {
                return true;
            }
        }
        false
    }

    /// Whether a carrier whose complete point set lies in `enclosure` misses
    /// the complete bounded support.
    pub(crate) fn carrier_misses_support(
        &mut self,
        model: &Model,
        enclosure: &[catallax_numbers::Interval; 3],
        surface: SurfaceId,
        dec: &mut Decider,
    ) -> bool {
        if !self.enabled {
            return false;
        }
        let Ok(Some(support)) = self.support(model, surface, dec) else {
            return false;
        };
        (0..3).any(|i| {
            let support = support[i].to_interval();
            enclosure[i].hi() < support.lo() || support.hi() < enclosure[i].lo()
        })
    }

    pub(crate) fn prove_bodies_apart(
        &mut self,
        model: &Model,
        a: BodyId,
        b: BodyId,
        dec: &mut Decider,
    ) -> Option<BodyApartProof> {
        self.bodies_apart(model, a, b, dec)
            .then_some(BodyApartProof { a, b })
    }

    /// Whether every lump of one solid is strictly apart from every lump of
    /// the other, by the boxes of their faces.
    pub(crate) fn bodies_apart(
        &mut self,
        model: &Model,
        a: BodyId,
        b: BodyId,
        dec: &mut Decider,
    ) -> bool {
        if !self.enabled {
            return false;
        }
        let mut lumps = |body: BodyId| -> Option<Vec<Box3>> {
            if model.body(body).kind != BodyKind::Solid {
                return None;
            }
            let mut boxes = Vec::new();
            for &lump in &model.body(body).lumps {
                let mut bounds = None;
                for &shell in &model.lump(lump).shells {
                    for &face in &model.shell(shell).faces {
                        include(&mut bounds, &self.face(model, face, dec)?);
                    }
                }
                boxes.push(bounds?);
            }
            (!boxes.is_empty()).then_some(boxes)
        };
        let (Some(a), Some(b)) = (lumps(a), lumps(b)) else {
            return false;
        };
        a.iter().all(|a| b.iter().all(|b| apart(a, b)))
    }
}

/// A short name for a surface kind, for refusals that name a missing meeting.
pub(crate) fn surface_name(surface: &Surface) -> &'static str {
    match surface {
        Surface::Plane(_) => "plane",
        Surface::Cylinder(_) => "cylinder",
        Surface::Cone(_) => "cone",
        Surface::Sphere(_) => "sphere",
        Surface::Torus(_) => "torus",
        Surface::TorusSheet(_) => "torus sheet",
        Surface::EllipticCylinder(_) => "elliptic cylinder",
        Surface::BSpline(spline) => match spline {
            SplineSupport::BezierExtrusion(_) => "Bezier extrusion spline",
            SplineSupport::EllipticOffset(_) => "fitted elliptic offset spline",
            SplineSupport::RollingEnvelope(_) => "fitted rolling-envelope spline",
            SplineSupport::RollingPipe(_) => "rolling-pipe spline",
            SplineSupport::RollingSaddle(_) => "fitted rolling-saddle spline",
            SplineSupport::RollingSphere(_) => "fitted rolling-sphere spline",
            SplineSupport::ObliqueDraft(_) => "oblique draft spline",
            SplineSupport::LinearFillet(_) => "linear fillet spline",
            SplineSupport::RollingMitre(_) => "rolling-mitre spline",
            SplineSupport::RationalPatch(_) => "rational patch",
            SplineSupport::RationalAtlas(_) => "rational atlas",
            _ => "spline",
        },
        _ => "surface",
    }
}

#[cfg(all(test, feature = "work-span"))]
mod span_tests {
    use super::TrimBounds;
    use catallax_numbers::{
        Decider, Rational,
        work_span::{parallel, read_at, with_work_span_on_this_thread, write_at},
    };
    use catallax_shapes::Pose;
    use catallax_topology::Model;

    // Adjacent faces once shared an observer key, so a task reading one face's
    // bounds waited for the task that wrote its neighbour's.
    #[test]
    fn adjacent_faces_and_supports_have_their_own_resources() {
        let mut model = Model::new();
        let body = catallax_topology::block(
            &mut model,
            [1i64, 1, 1].map(Rational::from),
            &Pose::identity(),
        )
        .unwrap();
        let faces = model.faces_of_body(body);
        let (first, second) = (faces[0], faces[1]);
        // An even index and its successor: the pair the old `index | 1` key merged.
        assert_eq!(first.index() % 2, 0);
        assert_eq!(first.index() + 1, second.index());
        let bounds = TrimBounds::default();
        let support = model.face(first).surface;
        let keys = [
            bounds.span_face(first),
            bounds.span_face(second),
            bounds.span_support(support),
        ];
        for (i, a) in keys.iter().enumerate() {
            for b in &keys[i + 1..] {
                assert_ne!(a, b);
            }
        }
        let ((), work) = catallax_numbers::with_operation_budget(1_000, || {
            let ((), report) = with_work_span_on_this_thread(|| {
                let mut dec = Decider::default();
                let tasks = parallel("faces");
                {
                    let _task = tasks.task();
                    dec.charge_work(10).unwrap();
                    let (kind, owner, index) = bounds.span_face(first);
                    write_at(kind, owner, index);
                }
                {
                    let _task = tasks.task();
                    let (kind, owner, index) = bounds.span_face(second);
                    read_at(kind, owner, index);
                    dec.charge_work(5).unwrap();
                }
            });
            // Each charge also pays one checkpoint unit. The second face does
            // not wait for the first: 11, not 17. Only the feature-level
            // model runs this loop's tasks in turn.
            let spans: Vec<_> = report.spans.iter().map(|row| row.span).collect();
            assert_eq!(spans, [17, 11, 11, 11, 11]);
        });
        assert_eq!(work.used, 17);
    }
}

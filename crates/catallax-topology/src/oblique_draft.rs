//! Complete finite boundary checks for obliquely drafted cylindrical cells.

use crate::{BodyId, BodyKind, EdgeId, FaceId, Loop, Model, Sense};
use catallax_numbers::{Decider, GEOMETRIC_GAP_UNPROVEN, Outcome, Real, Refusal, Sign, StrictSign};
use catallax_shapes::{
    BSplineCurve, Curve, ObliqueDraft, Plane, Point3, SplineBoundary, SplineSupport, Surface,
};

fn bad(detail: &str) -> Refusal {
    Refusal::new(GEOMETRIC_GAP_UNPROVEN, detail)
}
fn sign(x: &Real, dec: &mut Decider) -> Result<Sign, Refusal> {
    if let Some(q) = x.rational_ref() {
        return Ok(q.sign());
    }
    match dec.sign(x) {
        Outcome::Proven(StrictSign::Positive) => Ok(Sign::Positive),
        Outcome::Proven(StrictSign::Negative) => Ok(Sign::Negative),
        Outcome::ProvenCoincident(_) => Ok(Sign::Zero),
        Outcome::Refused(r) => Err(r),
        Outcome::NotExact(r) => Err(Refusal::from(r)),
    }
}
fn sense(s: Sense) -> i64 {
    if s == Sense::Same { 1 } else { -1 }
}
fn ring(model: &Model, face: FaceId, index: usize) -> Result<(EdgeId, Sense), Refusal> {
    let lp = model
        .face(face)
        .loops
        .get(index)
        .ok_or_else(|| bad("missing ruled-cell trim"))?;
    let Loop::Coedges(cs) = model.loop_(*lp) else {
        return Err(bad("ruled trim is not a physical ring"));
    };
    let [c] = cs.as_slice() else {
        return Err(bad("ruled trim contains additional edges"));
    };
    let c = model.coedge(*c);
    let e = model.edge(c.edge);
    if e.start.is_some() || e.end.is_some() {
        return Err(bad("closed ruled ring has seam vertices"));
    }
    Ok((c.edge, c.sense))
}

/// A complete ruled solid with two disks or two annular caps and a through hole.
/// The checked carrier proves star-shaped, convex horizontal sections; it does
/// not claim the entire solid is convex or replace a missing trim by its hull.
#[derive(Clone, Copy, Debug)]
pub struct DraftedCylinder<'a> {
    model: &'a Model,
    support: &'a ObliqueDraft,
    faces: [FaceId; 3],
    rings: [EdgeId; 2],
    hole: Option<(&'a ObliqueDraft, FaceId, [EdgeId; 2])>,
}
impl<'a> DraftedCylinder<'a> {
    /// Checks all shells, faces, exact shared rows and oriented cap domains.
    pub fn check(model: &'a Model, body: BodyId, dec: &mut Decider) -> Result<Self, Refusal> {
        dec.check_deadline()?;
        let entity = model.body(body);
        let [lump] = entity.lumps.as_slice() else {
            return Err(bad("ruled cell needs one material lump"));
        };
        let [shell] = model.lump(*lump).shells.as_slice() else {
            return Err(bad("ruled cell needs one complete shell"));
        };
        let fs = &model.shell(*shell).faces;
        if entity.kind != BodyKind::Solid
            || !(3..=4).contains(&fs.len())
            || !model.check(body).is_empty()
            || fs.iter().collect::<std::collections::BTreeSet<_>>().len() != fs.len()
        {
            return Err(bad("ruled cell is not one valid three- or four-face solid"));
        }
        let mut side = None;
        let mut hole = None;
        for &f in fs {
            match model.geom().surface(model.face(f).surface) {
                Surface::BSpline(SplineSupport::ObliqueDraft(s))
                    if model.face(f).sense == Sense::Same && side.is_none() =>
                {
                    side = Some((f, s))
                }
                Surface::BSpline(SplineSupport::ObliqueDraft(s))
                    if model.face(f).sense == Sense::Reversed && hole.is_none() =>
                {
                    hole = Some((f, s))
                }
                Surface::Plane(_) => (),
                _ => return Err(bad("ruled cell contains an extra or unsupported carrier")),
            }
        }
        let (side, support) = side.ok_or_else(|| bad("missing checked oblique draft support"))?;
        if model.face(side).sense != Sense::Same || model.face(side).loops.len() != 2 {
            return Err(bad("ruled side has incorrect sense or incomplete trims"));
        }
        if let Some((inner, s)) = hole {
            if model.face(inner).loops.len() != 2
                || sign(&s.origin().minus(support.origin()).norm_squared(), dec)? != Sign::Zero
                || sign(&(s.axis().dot(support.axis()) - Real::one()), dec)? != Sign::Zero
                || sign(&(s.across().dot(support.across()) - Real::one()), dec)? != Sign::Zero
            {
                return Err(bad(
                    "ruled hole does not share the complete outer frame and height",
                ));
            }
            for i in 0..2 {
                if sign(&(&s.heights()[i] - &support.heights()[i]), dec)? != Sign::Zero {
                    return Err(bad("ruled hole has a different finite height"));
                }
            }
            let lean = ((&support.lean()[0] - &s.lean()[0]).square()
                + (&support.lean()[1] - &s.lean()[1]).square())
            .sqrt();
            for height in support.heights() {
                if sign(
                    &(&(&support.amplitude(height) - &s.amplitude(height)) - &lean),
                    dec,
                )? != Sign::Positive
                {
                    return Err(bad("ruled hole touches or crosses the outer wall"));
                }
            }
        }
        let mut caps = [None, None];
        for &f in fs {
            if f == side || hole.is_some_and(|(inner, _)| f == inner) {
                continue;
            }
            let face = model.face(f);
            let Surface::Plane(p) = model.geom().surface(face.surface) else {
                return Err(bad("ruled cap is not planar"));
            };
            Plane::new(p.origin.clone(), p.normal.clone(), p.ref_dir.clone(), dec)?;
            if sign(&p.normal.cross(support.axis()).norm_squared(), dec)? != Sign::Zero {
                return Err(bad("ruled cap is not perpendicular to pull"));
            }
            let h = p.origin.minus(support.origin()).dot(support.axis());
            let slot = if sign(&(&h - &support.heights()[0]), dec)? == Sign::Zero {
                0
            } else if sign(&(&h - &support.heights()[1]), dec)? == Sign::Zero {
                1
            } else {
                return Err(bad("ruled cap does not meet its exact boundary row"));
            };
            let normal = match sign(&p.normal.dot(support.axis()), dec)? {
                Sign::Positive => 1,
                Sign::Negative => -1,
                Sign::Zero => return Err(bad("zero cap normal")),
            };
            if caps[slot].replace(f).is_some()
                || face.loops.len() != 1 + usize::from(hole.is_some())
                || normal * sense(face.sense) != [-1, 1][slot]
            {
                return Err(bad("ruled cap is duplicated, incomplete or reversed"));
            }
        }
        let [Some(lower), Some(upper)] = caps else {
            return Err(bad("ruled cell is missing a cap"));
        };
        let faces = [lower, upper, side];
        let rings = check_rings(model, faces, support, 0)?;
        let hole = hole
            .map(|(face, support)| {
                check_rings(model, [lower, upper, face], support, 1)
                    .map(|rings| (support, face, rings))
            })
            .transpose()?;
        if hole.is_some_and(|(_, _, inner)| rings.iter().any(|edge| inner.contains(edge))) {
            return Err(bad("ruled hole reuses an exterior edge"));
        }
        Ok(Self {
            model,
            support,
            faces,
            rings,
            hole,
        })
    }
    /// Returns a checked cell only when this body's own faces contain its side.
    pub fn of_body(
        model: &'a Model,
        body: BodyId,
        dec: &mut Decider,
    ) -> Result<Option<Self>, Refusal> {
        if !model.faces_of_body(body).iter().any(|f| {
            matches!(
                model.geom().surface(model.face(*f).surface),
                Surface::BSpline(SplineSupport::ObliqueDraft(_))
            )
        }) {
            return Ok(None);
        }
        Self::check(model, body, dec).map(Some)
    }
    /// Checks the complete solid that owns the queried face, without borrowing
    /// a neighbouring solid's cap or inferring a domain from a bare support.
    pub fn of_face(
        model: &'a Model,
        face: FaceId,
        dec: &mut Decider,
    ) -> Result<Option<Self>, Refusal> {
        for body in model.bodies_with_face(face) {
            if let Some(cell) = Self::of_body(model, body, dec)? {
                return Ok(Some(cell));
            }
        }
        Ok(None)
    }
    /// Negative inside, zero on boundary, positive outside the represented cell.
    pub fn locate(&self, p: &Point3, dec: &mut Decider) -> Result<Sign, Refusal> {
        let band = self.band(p, dec)?;
        if band == Sign::Positive {
            return Ok(band);
        }
        Ok(match self.radial_side(p, dec)? {
            Sign::Negative => band,
            s => s,
        })
    }
    fn radial_side(&self, p: &Point3, dec: &mut Decider) -> Result<Sign, Refusal> {
        let outer = self.support.side_at(p, dec)?;
        if outer != Sign::Negative {
            return Ok(outer);
        }
        if let Some((inner, _, _)) = self.hole {
            return Ok(match inner.side_at(p, dec)? {
                Sign::Negative => Sign::Positive,
                Sign::Zero => Sign::Zero,
                Sign::Positive => Sign::Negative,
            });
        }
        Ok(outer)
    }
    fn band(&self, p: &Point3, dec: &mut Decider) -> Result<Sign, Refusal> {
        let z = self.support.coordinates(p)[2].clone();
        Ok(
            match (
                sign(&(&z - &self.support.heights()[0]), dec)?,
                sign(&(&z - &self.support.heights()[1]), dec)?,
            ) {
                (Sign::Negative, _) | (_, Sign::Positive) => Sign::Positive,
                (Sign::Zero, _) | (_, Sign::Zero) => Sign::Zero,
                _ => Sign::Negative,
            },
        )
    }
    /// Trim membership for a point already proved on the selected support.
    pub fn locate_on_face(
        &self,
        f: FaceId,
        p: &Point3,
        dec: &mut Decider,
    ) -> Result<Sign, Refusal> {
        if f == self.faces[2] || self.hole.is_some_and(|(_, inner, _)| f == inner) {
            self.band(p, dec)
        } else {
            self.radial_side(p, dec)
        }
    }
    /// Exact retained exterior ruled support.
    pub fn support(&self) -> &'a ObliqueDraft {
        self.support
    }
    /// Inner support, reversed material face and complete lower/upper rings.
    pub fn hole(&self) -> Option<(&'a ObliqueDraft, FaceId, [EdgeId; 2])> {
        self.hole
    }
    /// Positive lower bound on the three-dimensional distance between the walls.
    /// Export and display must keep their combined displacement below this bound.
    pub fn wall_clearance(&self, dec: &mut Decider) -> Result<Option<Real>, Refusal> {
        let Some((inner, _, _)) = self.hole else {
            return Ok(None);
        };
        let minimum = |a: Real, b: Real, dec: &mut Decider| -> Result<Real, Refusal> {
            Ok(if sign(&(&a - &b), dec)? == Sign::Positive {
                b
            } else {
                a
            })
        };
        let lean = (inner.lean()[0].square() + inner.lean()[1].square()).sqrt();
        let radius = &minimum(
            inner.amplitude(&inner.heights()[0]),
            inner.amplitude(&inner.heights()[1]),
            dec,
        )? - &lean;
        let delta = ((&self.support.lean()[0] - &inner.lean()[0]).square()
            + (&self.support.lean()[1] - &inner.lean()[1]).square())
        .sqrt();
        let gap = minimum(
            self.support.amplitude(&inner.heights()[0]) - inner.amplitude(&inner.heights()[0]),
            self.support.amplitude(&inner.heights()[1]) - inner.amplitude(&inner.heights()[1]),
            dec,
        )? - delta;
        let axial = if sign(inner.tangent(), dec)? == Sign::Negative {
            -inner.tangent()
        } else {
            inner.tangent().clone()
        };
        // The inner residual is rho - A(z) - lean.dot(ray). Any segment
        // shorter than radius/2 stays at rho >= radius/2. Its transverse
        // derivative norm is at most 1 + 2*|lean|/radius and its axial
        // derivative is tangent. Their sum bounds the full gradient norm,
        // including wall points at different heights after XYZ rounding.
        minimum(
            &radius / &Real::from_i64(2),
            gap / (Real::one() + Real::from_i64(2) * lean / radius + axial),
            dec,
        )
        .map(Some)
    }
    /// Exact volume divided by pi, with the through hole removed.
    pub fn volume_over_pi(&self) -> Real {
        let outer = self.support.volume_over_pi();
        self.hole.map_or(outer.clone(), |(inner, _, _)| {
            outer - inner.volume_over_pi()
        })
    }
    /// Lower cap, upper cap and exterior ruled side, in material order.
    pub fn faces(&self) -> [FaceId; 3] {
        self.faces
    }
    /// Complete lower and upper exterior physical rings.
    pub fn rings(&self) -> [EdgeId; 2] {
        self.rings
    }
    /// Immutable checked model.
    pub fn model(&self) -> &'a Model {
        self.model
    }
}

fn check_rings(
    model: &Model,
    faces: [FaceId; 3],
    support: &ObliqueDraft,
    wall: usize,
) -> Result<[EdgeId; 2], Refusal> {
    let direction = if wall == 0 { 1 } else { -1 };
    let mut rings = Vec::with_capacity(2);
    for i in 0..2 {
        let (e, s) = ring(model, faces[i], wall)?;
        let edge = model.edge(e);
        let Curve::BSpline(curve) = model.geom().curve(edge.curve) else {
            return Err(bad("ruled cap ring is not the exact spline row"));
        };
        let expected = BSplineCurve::boundary(
            support.carrier(),
            [SplineBoundary::VMin, SplineBoundary::VMax][i],
        );
        if !curve.same_definition(&expected)
            || curve.is_closed() != Some(true)
            || sense(edge.sense) * sense(s) != direction * [-1, 1][i]
        {
            return Err(bad(
                "ruled cap ring differs from its complete oriented boundary row",
            ));
        }
        rings.push(e);
    }
    if rings[0] == rings[1] {
        return Err(bad("ruled cell uses one edge for different caps"));
    }
    let mut used = [false, false];
    for i in 0..2 {
        let (e, s) = ring(model, faces[2], i)?;
        let j = rings
            .iter()
            .position(|x| *x == e)
            .ok_or_else(|| bad("ruled side is detached from its cap"))?;
        if std::mem::replace(&mut used[j], true)
            || sense(model.edge(e).sense) * sense(s) != direction * [1, -1][j]
        {
            return Err(bad(
                "ruled side does not use both rings opposite their cap uses",
            ));
        }
    }
    Ok([rings[0], rings[1]])
}

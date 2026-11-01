//! Complete geometry checks for the first finite spline solid family.

use crate::{BodyId, BodyKind, EdgeId, FaceId, Loop, Model, Sense};
use catallax_numbers::{Decider, GEOMETRIC_GAP_UNPROVEN, Outcome, Real, Refusal, Sign, StrictSign};
use catallax_shapes::{
    BSplineCurve, Curve, Plane, Point3, SplineBoundary, SplineSupport, SplineSurface, Surface,
};

fn bad(detail: impl Into<String>) -> Refusal {
    Refusal::new(GEOMETRIC_GAP_UNPROVEN, detail)
}
fn sign(value: &Real, dec: &mut Decider) -> Result<Sign, Refusal> {
    match dec.sign(value) {
        Outcome::Proven(StrictSign::Positive) => Ok(Sign::Positive),
        Outcome::Proven(StrictSign::Negative) => Ok(Sign::Negative),
        Outcome::ProvenCoincident(_) => Ok(Sign::Zero),
        Outcome::NotExact(reason) => Err(Refusal::new(reason.code, reason.detail)),
        Outcome::Refused(reason) => Err(reason),
    }
}
fn orientation(sense: Sense) -> i64 {
    match sense {
        Sense::Same => 1,
        Sense::Reversed => -1,
    }
}
fn ring(model: &Model, face: FaceId, index: usize) -> Result<(EdgeId, Sense), Refusal> {
    let lp = model
        .face(face)
        .loops
        .get(index)
        .ok_or_else(|| bad("missing complete prism trim"))?;
    let Loop::Coedges(coedges) = model.loop_(*lp) else {
        return Err(bad("prism trim is not a closed ring"));
    };
    let [coedge] = coedges.as_slice() else {
        return Err(bad("prism trim contains extra physical seam edges"));
    };
    let coedge = model.coedge(*coedge);
    let edge = model.edge(coedge.edge);
    if edge.start.is_some() || edge.end.is_some() {
        return Err(bad("closed prism ring has artificial seam vertices"));
    }
    Ok((coedge.edge, coedge.sense))
}

// A prism side has a checked support; a cap is one loop of one fitted edge.
fn may_bound_prism(model: &Model, face: FaceId) -> bool {
    let entity = model.face(face);
    match model.geom().surface(entity.surface) {
        Surface::BSpline(SplineSupport::EllipticOffset(_)) => true,
        Surface::Plane(_) => {
            let [boundary] = entity.loops.as_slice() else {
                return false;
            };
            let Loop::Coedges(coedges) = model.loop_(*boundary) else {
                return false;
            };
            let [coedge] = coedges.as_slice() else {
                return false;
            };
            let edge = model.edge(model.coedge(*coedge).edge);
            matches!(model.geom().curve(edge.curve), Curve::BSpline(_))
        }
        _ => false,
    }
}

/// Borrowed valid represented solid: one checked convex periodic spline sweep,
/// two complete planar cap domains and their two shared physical rings.
///
/// This proves actual B-rep geometry, not correspondence to a recipe. Diagnostic
/// marks are neither consulted nor promoted to geometric authority. Features
/// must separately bind retained intent, local limits and whole-result evidence.
#[derive(Clone, Copy, Debug)]
pub struct SplinePrism<'a> {
    model: &'a Model,
    body: BodyId,
    spline: &'a SplineSurface,
    faces: [FaceId; 3],
    rings: [EdgeId; 2],
}
impl<'a> SplinePrism<'a> {
    /// Checks every reachable face, trim, incidence and material orientation.
    /// Partial covers, voids, detached caps and serialized seam assertions refuse.
    pub fn check(model: &'a Model, body: BodyId, dec: &mut Decider) -> Result<Self, Refusal> {
        dec.check_deadline()?;
        let entity = model.body(body);
        let [lump] = entity.lumps.as_slice() else {
            return Err(bad("spline prism requires one material lump"));
        };
        let [shell] = model.lump(*lump).shells.as_slice() else {
            return Err(bad("spline prism cannot hide extra shells or voids"));
        };
        let faces = &model.shell(*shell).faces;
        if entity.kind != BodyKind::Solid
            || faces.len() != 3
            || faces[0] == faces[1]
            || faces[0] == faces[2]
            || faces[1] == faces[2]
            || !model.check(body).is_empty()
        {
            return Err(bad("spline prism is not a valid complete three-face solid"));
        }
        let mut side = None;
        for &face in faces {
            match model.geom().surface(model.face(face).surface) {
                Surface::Plane(_) => (),
                Surface::BSpline(SplineSupport::EllipticOffset(spline)) => {
                    if side.replace((face, spline)).is_some() {
                        return Err(bad("spline prism has more than one side"));
                    }
                }
                _ => {
                    return Err(bad(
                        "spline prism has an unsupported or falsely labelled face",
                    ));
                }
            }
        }
        let (side, spline) =
            side.ok_or_else(|| bad("spline prism is missing its genuine checked side"))?;
        if model.face(side).sense != Sense::Same || model.face(side).loops.len() != 2 {
            return Err(bad("spline side has wrong material sense or trim count"));
        }
        let Surface::EllipticCylinder(original) = spline.definition().original() else {
            unreachable!("checked elliptic sweep source");
        };
        let mut caps = [None, None];
        for &face in faces {
            if face == side {
                continue;
            }
            let Surface::Plane(plane) = model.geom().surface(model.face(face).surface) else {
                return Err(bad("spline cap is not a plane"));
            };
            Plane::new(
                plane.origin.clone(),
                plane.normal.clone(),
                plane.ref_dir.clone(),
                dec,
            )?;
            if sign(&plane.normal.cross(spline.axis()).norm_squared(), dec)? != Sign::Zero {
                return Err(bad("cap plane is not perpendicular to the exact sweep"));
            }
            let normal = match sign(&plane.normal.dot(spline.axis()), dec)? {
                Sign::Positive => 1,
                Sign::Negative => -1,
                Sign::Zero => return Err(bad("cap normal vanished")),
            };
            let level = plane.origin.minus(&original.origin).dot(spline.axis());
            let slot = if sign(&(&level - spline.axial()[0]), dec)? == Sign::Zero {
                0
            } else if sign(&(&level - spline.axial()[1]), dec)? == Sign::Zero {
                1
            } else {
                return Err(bad("cap plane is detached from the finite spline sweep"));
            };
            if caps[slot].replace(face).is_some()
                || model.face(face).loops.len() != 1
                || normal * orientation(model.face(face).sense) != [-1, 1][slot]
            {
                return Err(bad("cap cover or material orientation is invalid"));
            }
        }
        let [Some(lower), Some(upper)] = caps else {
            return Err(bad("spline prism does not cover both complete caps"));
        };
        let ordered = [lower, upper, side];
        let mut rings = [None, None];
        for index in 0..2 {
            let (edge, travel) = ring(model, ordered[index], 0)?;
            let Curve::BSpline(curve) = model.geom().curve(model.edge(edge).curve) else {
                return Err(bad("spline cap was relabelled as an analytic curve"));
            };
            let expected = BSplineCurve::boundary(
                spline.carrier(),
                [SplineBoundary::VMin, SplineBoundary::VMax][index],
            );
            if !curve.same_definition(&expected)
                || curve.is_closed() != Some(true)
                || orientation(model.edge(edge).sense) * orientation(travel) != [-1, 1][index]
            {
                return Err(bad(
                    "cap trim is not the correctly oriented complete shared boundary",
                ));
            }
            rings[index] = Some(edge);
        }
        let [Some(lower), Some(upper)] = rings else {
            return Err(bad("incomplete physical ring inventory"));
        };
        if lower == upper {
            return Err(bad("distinct caps share one ring"));
        }
        let mut used = [false, false];
        for index in 0..2 {
            let (edge, travel) = ring(model, side, index)?;
            let slot = if edge == lower {
                0
            } else if edge == upper {
                1
            } else {
                return Err(bad("side trim is disconnected from its cap"));
            };
            if std::mem::replace(&mut used[slot], true)
                || orientation(model.edge(edge).sense) * orientation(travel) != [1, -1][slot]
            {
                return Err(bad(
                    "side and cap rings do not give one oriented material boundary",
                ));
            }
        }
        // The opaque support proves one regular convex closed profile. Complete
        // common cap loops bound its convex disks, not merely nearby planes.
        // Opposite coedge uses and one finite axial interval therefore give the
        // boundary of exactly one prism, with no branch jump or omitted region.
        Ok(Self {
            model,
            body,
            spline,
            faces: ordered,
            rings: [lower, upper],
        })
    }
    /// The checked prism of a body, when it has a checked side and is one.
    /// `None` claims nothing about the body; the caller keeps its other paths.
    /// Any refusal other than a failed geometry check is returned.
    pub fn of_body(
        model: &'a Model,
        body: BodyId,
        dec: &mut Decider,
    ) -> Result<Option<Self>, Refusal> {
        let has_side = model.faces_of_body(body).into_iter().any(|face| {
            matches!(
                model.geom().surface(model.face(face).surface),
                Surface::BSpline(SplineSupport::EllipticOffset(_))
            )
        });
        if !has_side {
            return Ok(None);
        }
        match Self::check(model, body, dec) {
            Ok(prism) => Ok(Some(prism)),
            Err(refusal) if refusal.code == GEOMETRIC_GAP_UNPROVEN => Ok(None),
            Err(refusal) => Err(refusal),
        }
    }
    /// The checked prism a face belongs to. A face with another trim, or one
    /// only sharing a ring with a prism, gets `None` and keeps its own domain.
    pub fn of_face(
        model: &'a Model,
        face: FaceId,
        dec: &mut Decider,
    ) -> Result<Option<Self>, Refusal> {
        if !may_bound_prism(model, face) {
            return Ok(None);
        }
        for body in model.bodies_with_face(face) {
            if let Some(prism) = Self::of_body(model, body, dec)? {
                return Ok(Some(prism));
            }
        }
        Ok(None)
    }
    /// Exact place of a point: negative inside the solid, zero on its boundary,
    /// positive outside. The side reads the represented section, not its target.
    pub fn locate(&self, point: &Point3, dec: &mut Decider) -> Result<Sign, Refusal> {
        let band = self.band(point, dec)?;
        if band == Sign::Positive {
            return Ok(Sign::Positive);
        }
        Ok(match self.spline.cross_section_side(point, dec)? {
            Sign::Negative => band,
            side => side,
        })
    }
    /// Exact place of a point on a boundary face's support within that closed
    /// face: negative in its interior, zero on its ring, positive off the face.
    pub fn locate_on_face(
        &self,
        face: FaceId,
        point: &Point3,
        dec: &mut Decider,
    ) -> Result<Sign, Refusal> {
        if face == self.faces[2] {
            self.band(point, dec)
        } else {
            self.spline.cross_section_side(point, dec)
        }
    }
    // Negative strictly between the caps, zero on a cap plane, positive beyond.
    fn band(&self, point: &Point3, dec: &mut Decider) -> Result<Sign, Refusal> {
        let axial = self.spline.axial_coordinate(point);
        let [lower, upper] = self.spline.axial();
        Ok(
            match (sign(&(&axial - lower), dec)?, sign(&(&axial - upper), dec)?) {
                (Sign::Negative, _) | (_, Sign::Positive) => Sign::Positive,
                (Sign::Zero, _) | (_, Sign::Zero) => Sign::Zero,
                _ => Sign::Negative,
            },
        )
    }
    /// Immutable model whose entities were checked.
    pub fn model(&self) -> &'a Model {
        self.model
    }
    /// Actual represented solid.
    pub fn body(&self) -> BodyId {
        self.body
    }
    /// Checked finite support and retained geometric definition.
    pub fn spline(&self) -> &'a SplineSurface {
        self.spline
    }
    /// Complete lower cap, upper cap and side face inventory.
    pub fn faces(&self) -> [FaceId; 3] {
        self.faces
    }
    /// Shared lower and upper physical ring edges.
    pub fn rings(&self) -> [EdgeId; 2] {
        self.rings
    }
}

//! Whether a stroke piece is hidden: whether the open ray from one of its
//! points toward the viewer meets the inside of the solid.
//!
//! OpenCASCADE's `HLRBRep_Data::Hide` hides an edge piece when a face covers
//! its middle point on the sheet and lies in front of it; SMLib's `SmHCR`
//! asks the same question with a ray from the piece. Here the ray is exact.
//! A crossing of a face's surface strictly inside the face, where the ray
//! is not tangent, passes between inside and outside, so the ray meets the
//! inside. Without one, the ray meets the boundary only at touches, at
//! edges and vertices, and along faces whose surface holds it. Between two
//! consecutive such parameters it is then wholly inside, wholly outside or
//! wholly on the boundary, and one exact point of each stretch decides, by
//! the Boolean engine's ray parity (`containment.rs`): a ray from the point
//! that crosses the boundary only transversally, strictly inside faces,
//! crosses it an odd number of times exactly when the point is inside
//! (the Jordan-Brouwer separation theorem). Here the parity rays skip faces
//! whose boxes they miss. A ray that only grazes the boundary, such as one
//! running along a face seen edge on, leaves its piece visible.

use catallax_boolean::face_contains;
use catallax_intersect::numerical_ray::{NumericalRayHit, patch_meetings};
use catallax_intersect::{Contact, CurveSurface, intersect_curve_surface_with_line_parameter};
use catallax_numbers::{Decider, MarkLog, Rational, Real, Refusal};
use catallax_shapes::{Curve, Line, Point3, PointId, RationalPatch, Vector3};
use catallax_topology::{FaceId, Loop, Model};

use super::geo::{Box3, V3, ray_meets};
use crate::queries::sign;
use crate::views::float::{double, rational};
use crate::views::{on_edge_span, same_point};

/// The solid seen by the rays: its faces with binary64 boxes widened by
/// `margin` (`None` where no cheap bound exists), and its vertices.
pub(super) struct Scene {
    pub faces: Vec<(FaceId, Option<Box3>)>,
    /// The rational patches of spline faces met by binary64 rays, with
    /// their control boxes.
    pub patches: std::collections::BTreeMap<FaceId, Vec<(RationalPatch, Box3)>>,
    pub vertices: Vec<(PointId, V3)>,
    pub margin: f64,
    /// The largest coordinate size of the body.
    pub size: f64,
}

/// Whether `point` lies on the boundary of `face`: on one of its edges'
/// spans or at an isolated vertex.
fn on_boundary(
    scratch: &mut Model,
    face: FaceId,
    point: PointId,
    dec: &mut Decider,
) -> Result<bool, Refusal> {
    for lp in scratch.face(face).loops.clone() {
        match scratch.loop_(lp).clone() {
            Loop::Vertex(vertex) => {
                if same_point(scratch, scratch.vertex(vertex).point, point, dec)? {
                    return Ok(true);
                }
            }
            Loop::Coedges(coedges) => {
                for coedge in coedges {
                    let edge = scratch.coedge(coedge).edge;
                    if on_edge_span(scratch, edge, point, dec)? {
                        return Ok(true);
                    }
                }
            }
        }
    }
    Ok(false)
}

/// What a binary64 ray meeting with a rational patch says, decided only
/// where it clears the tolerance: the ray's own start, a transversal
/// crossing strictly inside the patch, or a touch or boundary meeting at
/// this distance along the unit ray.
enum Meeting {
    Start,
    Crossing,
    Stop(Rational),
}

fn meeting(
    patch: &RationalPatch,
    hit: &NumericalRayHit,
    tolerance: &Rational,
    margin: f64,
    dec: &mut Decider,
) -> Result<Meeting, Refusal> {
    if hit.distance_mm <= *tolerance {
        return Ok(Meeting::Start);
    }
    let transversal = hit.normal_dot.lo().sign() == catallax_numbers::Sign::Positive
        || hit.normal_dot.hi().sign() == catallax_numbers::Sign::Negative;
    if !transversal {
        return Ok(Meeting::Stop(hit.distance_mm.clone()));
    }
    // Distance to the patch's boundary rows, to first order: the parameter
    // gap times the speed there.
    let at = hit.uv.clone().map(Real::from);
    let [du, dv] = patch.derivatives_at(&at, dec)?;
    let floor = Rational::from_f64_exact(margin).expect("a finite margin");
    let speed = |v: &Vector3, dec: &mut Decider| -> Result<f64, Refusal> {
        Ok(double(&v.norm_squared(), &floor, dec)?.sqrt())
    };
    let speeds = [speed(&du, dec)?, speed(&dv, dec)?];
    let domain = patch.domain();
    let mut gap = f64::INFINITY;
    for axis in 0..2 {
        let x = hit.uv[axis].to_f64_nearest();
        let [lo, hi] = [
            domain[axis][0].to_f64_nearest(),
            domain[axis][1].to_f64_nearest(),
        ];
        gap = gap.min((x - lo).min(hi - x) * speeds[axis]);
    }
    Ok(if gap > 2.0 * margin {
        Meeting::Crossing
    } else {
        Meeting::Stop(hit.distance_mm.clone())
    })
}

/// The binary64 meetings of a ray with a spline face's charts whose boxes
/// it meets, or `None` where the bounded solver refuses (a ray running along
/// the patch, or a meeting rounded off its domain); the exact intersector
/// then decides that face.
#[allow(clippy::too_many_arguments)]
fn meetings(
    charts: &[(RationalPatch, Box3)],
    origin: &Point3,
    o64: V3,
    dir: &Vector3,
    d64: V3,
    tolerance: &Rational,
    margin: f64,
    dec: &mut Decider,
) -> Result<Option<Vec<Meeting>>, Refusal> {
    let mut out = Vec::new();
    for (patch, bounds) in charts {
        if !ray_meets(bounds, o64, d64, margin) {
            continue;
        }
        let found = match patch_meetings(patch, origin, dir, tolerance, dec) {
            Ok(found) => found,
            Err(refusal) if refusal.code == "intersect.numerical-ray" => return Ok(None),
            Err(refusal) => return Err(refusal),
        };
        for hit in found.hits {
            out.push(meeting(patch, &hit, tolerance, margin, dec)?);
        }
    }
    Ok(Some(out))
}

/// Source: judgement. Directions of the parity rays: small primes, off the
/// axis-aligned and diagonal directions parts are built along. Any direction
/// off a set of measure zero decides; the next is tried when one grazes.
const PROBES: [[i64; 3]; 4] = [[2, 3, 5], [-3, 5, 7], [5, -7, 11], [7, 11, -13]];

/// Whether `point`, which lies on no face, is inside the solid: the parity
/// of a ray's transversal crossings of face interiors.
fn inside(
    scratch: &mut Model,
    scene: &Scene,
    point: &Point3,
    p64: V3,
    marks: &mut MarkLog,
    numeric: &mut bool,
    dec: &mut Decider,
) -> Result<bool, Refusal> {
    'probe: for probe in PROBES {
        dec.check_deadline()?;
        let dir = Vector3::new(
            Real::from_i64(probe[0]),
            Real::from_i64(probe[1]),
            Real::from_i64(probe[2]),
        );
        let d64 = probe.map(|x| x as f64);
        let ray = scratch.geom_mut().add_curve(Curve::Line(Line {
            origin: point.clone(),
            dir: dir.clone(),
        }));
        let mut crossings = 0_usize;
        for &(face, bounds) in &scene.faces {
            if let Some(b) = bounds
                && !ray_meets(&b, p64, d64, scene.margin)
            {
                continue;
            }
            if let Some(charts) = scene.patches.get(&face) {
                let tolerance = Rational::from_f64_exact(scene.margin).expect("a finite margin");
                if let Some(found) =
                    meetings(charts, point, p64, &dir, d64, &tolerance, scene.margin, dec)?
                {
                    *numeric = true;
                    for m in found {
                        match m {
                            Meeting::Crossing => crossings += 1,
                            // Within the tolerance of the boundary: on it.
                            Meeting::Start => return Ok(false),
                            Meeting::Stop(_) => continue 'probe,
                        }
                    }
                    continue;
                }
            }
            let surface = scratch.face(face).surface;
            let hits =
                intersect_curve_surface_with_line_parameter(scratch.geom_mut(), ray, surface, dec)?;
            marks.extend(&hits.marks);
            let hits = match hits.value {
                CurveSurface::Disjoint => continue,
                CurveSurface::Contained => continue 'probe,
                CurveSurface::Meets(hits) => hits,
            };
            for hit in hits {
                let at = scratch.geom().position(hit.point).minus(point).dot(&dir);
                if sign(&at, dec)? <= 0 {
                    continue;
                }
                let position = scratch.geom().position(hit.point);
                if !face_contains(scratch, face, &position, dec)? {
                    continue;
                }
                if hit.contact != Contact::Crossing || on_boundary(scratch, face, hit.point, dec)? {
                    continue 'probe;
                }
                crossings += 1;
            }
        }
        return Ok(crossings % 2 == 1);
    }
    Err(Refusal::new(
        "app.view-classify",
        "every parity ray grazed the boundary",
    ))
}

/// Sorted distinct values, by exact comparison.
fn insert(values: &mut Vec<Real>, value: Real, dec: &mut Decider) -> Result<(), Refusal> {
    let mut at = values.len();
    for (i, v) in values.iter().enumerate() {
        match sign(&(&value - v), dec)? {
            0 => return Ok(()),
            s if s < 0 => {
                at = i;
                break;
            }
            _ => {}
        }
    }
    values.insert(at, value);
    Ok(())
}

/// Whether the open ray from `p` along `d` meets the inside of the solid.
/// `p64` and `d64` are binary64 images for the box tests. Rational patch
/// faces are met by SMLib's bounded binary64 line/surface solver
/// (`patch_meetings`, as the Boolean engine's numerical loft
/// classification does), which sets `numeric`.
#[allow(clippy::too_many_arguments)]
pub(super) fn hidden(
    scratch: &mut Model,
    scene: &Scene,
    p: &Point3,
    p64: V3,
    d: &Vector3,
    d64: V3,
    marks: &mut MarkLog,
    numeric: &mut bool,
    dec: &mut Decider,
) -> Result<bool, Refusal> {
    let ray = scratch.geom_mut().add_curve(Curve::Line(Line {
        origin: p.clone(),
        dir: d.clone(),
    }));
    let mut stops: Vec<Real> = Vec::new();
    let mut along: Vec<FaceId> = Vec::new();
    let mut root_length: Option<Real> = None;
    for &(face, bounds) in &scene.faces {
        dec.check_deadline()?;
        if let Some(b) = bounds
            && !ray_meets(&b, p64, d64, scene.margin)
        {
            continue;
        }
        if let Some(charts) = scene.patches.get(&face) {
            let tolerance = Rational::from_f64_exact(scene.margin).expect("a finite margin");
            if let Some(found) = meetings(charts, p, p64, d, d64, &tolerance, scene.margin, dec)? {
                *numeric = true;
                let unit = root_length
                    .get_or_insert_with(|| d.norm_squared().sqrt())
                    .clone();
                for m in found {
                    match m {
                        Meeting::Start => {}
                        Meeting::Crossing => return Ok(true),
                        Meeting::Stop(distance) => {
                            insert(&mut stops, &Real::from(distance) * &unit, dec)?
                        }
                    }
                }
                continue;
            }
        }
        let surface = scratch.face(face).surface;
        let hits =
            intersect_curve_surface_with_line_parameter(scratch.geom_mut(), ray, surface, dec)?;
        marks.extend(&hits.marks);
        let hits = match hits.value {
            CurveSurface::Disjoint => continue,
            CurveSurface::Contained => {
                along.push(face);
                continue;
            }
            CurveSurface::Meets(hits) => hits,
        };
        for hit in hits {
            let at = scratch.geom().position(hit.point).minus(p).dot(d);
            if sign(&at, dec)? <= 0 {
                continue;
            }
            let position = scratch.geom().position(hit.point);
            if !face_contains(scratch, face, &position, dec)? {
                continue;
            }
            if hit.contact == Contact::Crossing && !on_boundary(scratch, face, hit.point, dec)? {
                return Ok(true);
            }
            insert(&mut stops, at, dec)?;
        }
    }
    for &(vertex, v64) in &scene.vertices {
        let r = [v64[0] - p64[0], v64[1] - p64[1], v64[2] - p64[2]];
        let ahead = r[0] * d64[0] + r[1] * d64[1] + r[2] * d64[2];
        let off = [
            r[0] - ahead * d64[0],
            r[1] - ahead * d64[1],
            r[2] - ahead * d64[2],
        ];
        if ahead < -scene.margin
            || (off[0] * off[0] + off[1] * off[1] + off[2] * off[2]).sqrt() > scene.margin
        {
            continue;
        }
        let r = scratch.geom().position(vertex).minus(p);
        if sign(&r.cross(d).norm_squared(), dec)? != 0 {
            continue;
        }
        let at = r.dot(d);
        if sign(&at, dec)? > 0 {
            insert(&mut stops, at, dec)?;
        }
    }
    let length = d.norm_squared();
    let two = Real::from_i64(2);
    let floor = Rational::from_f64_exact(scene.margin).expect("a finite margin");
    let mut previous = Real::zero();
    for stop in stops {
        dec.check_deadline()?;
        // A rational middle keeps the parity ray's numbers short. Binary64
        // images within half the floor of each end place it strictly
        // between them when they lie more than a floor apart.
        let ends = [double(&previous, &floor, dec)?, double(&stop, &floor, dec)?];
        let middle = if ends[1] - ends[0] > 2.0 * scene.margin {
            &Real::from(rational((ends[0] + ends[1]) / 2.0)?) / &length
        } else {
            &(&previous + &stop) / &(&two * &length)
        };
        previous = stop;
        let probe = p.plus(&d.scaled(&middle));
        let mut boundary = false;
        for &face in &along {
            if face_contains(scratch, face, &probe, dec)? {
                boundary = true;
                break;
            }
        }
        if boundary {
            continue;
        }
        let probe64 = [
            double(&probe.x, &floor, dec)?,
            double(&probe.y, &floor, dec)?,
            double(&probe.z, &floor, dec)?,
        ];
        if inside(scratch, scene, &probe, probe64, marks, numeric, dec)? {
            return Ok(true);
        }
    }
    Ok(false)
}

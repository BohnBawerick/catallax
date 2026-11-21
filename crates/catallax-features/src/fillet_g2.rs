//! Curvature-continuous (G2) fillets on straight edges between planes,
//! RECIPE.md section 4.5.3.
//!
//! This follows SMLib's blend-curve cross section,
//! `SmBlendCurveCrossSectionFSG::CreateSurfaceFromCurves` in usd-brep's
//! `SmFilletStandardSolver.cpp`. The rails are the constant-radius rolling
//! ball's contacts. Each section is the quintic Hermite curve through the two
//! contacts with the circular section's end tangent directions. Its tangent
//! lengths come from `SmBSplineCurve::FairBlendCurveDerivatives`, times the
//! thumbweight, and its second derivatives from the supports
//! (`SmSurface::ComputeHigherOrderDerivs`). A plane's second derivative is
//! zero, so the first three controls lie on the contact's tangent line and the
//! section's curvature at the rail is the plane's, zero. A straight centre
//! line makes SMLib's surface the ruled surface between the end sections
//! (`SmBSplineSurface::CreateRuledSurface`), here an exact degree-(5,1) patch.
//!
//! The applicable inputs and every proof but the section are the exact
//! variable-radius constructor's: one valid convex polyhedral shell, convex
//! straight edges whose ends are three-face vertices with plane end walls
//! normal to the edge, the tangent-triangle hull strictly clear of every face
//! but the edge's own and its end walls, and selected hulls apart.

use super::*;
use catallax_shapes::{BSplineCurve, BSplineSurface, LinearFillet, RationalPatch, SplineBoundary};
use catallax_topology::RationalPatchFace;

const UNSUPPORTED: &str = "fillet.unsupported-g2-contact";

fn unsupported(detail: &str) -> Refusal {
    Refusal::new(UNSUPPORTED, detail)
}

struct Band {
    edge: EdgeId,
    ends: [VertexId; 2],
    faces: [FaceId; 2],
    senses: [Sense; 2],
    caps: [FaceId; 2],
    /// The rolling ball's exact circular band: its rows are the contact
    /// foot, the edge station and the other foot, which bound every section.
    ball: LinearFillet,
    carrier: BSplineSurface,
    origin: Origin,
}

/// Binary64 holds 52 fraction bits, so its numbers just below 1 are
/// `2^-52` apart. SMLib evaluates its fairing rule at that precision.
const STEP_BITS: i64 = 52;

/// The ratio of the quintic's control step to the distance from a contact to
/// the edge station. `sine` is the sine of half the material angle at the
/// edge, the circular band's middle weight.
///
/// SMLib's fairing makes the cubic Hermite polygon's three legs equal. With
/// the contacts `s` from the station, the equal leg is `t = 2 s sine / (1 + 2
/// sine)`, the smallest positive root of its quadratic, and the end tangent
/// is `3 t`, times the thumbweight `w`. A quintic's first control step is a
/// fifth of its end tangent. So the step is `6 w sine / (5 (1 + 2 sine))` of
/// `s`. The recipe keeps `w` at most 5/4, so twice the step stays under the
/// whole distance: `12 sine / (5 (1 + 2 sine))` is below 4/5 for `sine <= 1`.
///
/// The step is rounded down to a multiple of `2^-52`, within SMLib's own
/// binary64 precision. No contact, tangent plane or rail curvature depends
/// on it, and it keeps the carrier rational wherever the contacts are, so
/// later exact decisions on the face stay in rational arithmetic.
fn control_step(sine: &Real, fullness: &Real, dec: &mut Decider) -> Result<Rational, Refusal> {
    let numerator = &Real::from_i64(6) * &(fullness * sine);
    let denominator = &Real::from_i64(5) * &(&Real::one() + &(sine + sine));
    let scaled = &(&numerator / &denominator) * &Real::from(Rational::one().mul_pow2(STEP_BITS));
    let ball = dec.enclose(&scaled, &Rational::from_ratio(1, 4).expect("nonzero"))?;
    if &ball.hi() - &ball.lo() >= Rational::one() {
        return Err(unsupported(
            "the curvature-continuous section step was not enclosed to one unit of 2^-52",
        ));
    }
    let mut floor = Rational::from_integer(ball.hi().floor());
    if Rational::from_integer(ball.lo().floor()) != floor {
        // An integer lies in the enclosure: decide its side exactly.
        if sign(&(&scaled - &Real::from(floor.clone())), dec)? == Sign::Negative {
            floor = &floor - &Rational::one();
        }
    }
    let step = floor.mul_pow2(-STEP_BITS);
    if step.sign() != Sign::Positive {
        return Err(unsupported(
            "the edge is too sharp for a positive curvature-continuous section step",
        ));
    }
    Ok(step)
}

/// The exact degree-(5,1) carrier, sections first and station second, as
/// the circular band orders its rows. Row `j` of a section is the contact
/// `a` moved `j` steps towards the station `p`, for `j` up to 2, and
/// mirrored from the other contact for the last three rows.
fn carrier(
    ball: &LinearFillet,
    fullness: &Real,
    dec: &mut Decider,
) -> Result<BSplineSurface, Refusal> {
    let rows = ball.carrier().controls();
    let sine = &ball.carrier().weights()[2];
    let step = Real::from(control_step(sine, fullness, dec)?);
    let mut controls = Vec::with_capacity(12);
    for (foot, steps) in [(0, 0), (0, 1), (0, 2), (2, 2), (2, 1), (2, 0)] {
        for station in 0..2 {
            let a = &rows[2 * foot + station];
            controls.push(if steps == 0 {
                a.clone()
            } else {
                let p = &rows[2 + station];
                a.plus(&p.minus(a).scaled(&(&step * &Real::from_i64(steps))))
            });
        }
    }
    let ends = |n: usize| {
        let mut knots = vec![Rational::zero(); n];
        knots.extend(vec![Rational::one(); n]);
        knots
    };
    BSplineSurface::new(
        [5, 1],
        [ends(6), ends(2)],
        [6, 2],
        controls,
        vec![Real::one(); 12],
        dec,
    )
}

/// Whether the chart normal points away from the material, towards the
/// removed edge. It is decided at the middle of the first section, against
/// the direction to its edge station.
fn outward_chart(band: &Band, dec: &mut Decider) -> Result<bool, Refusal> {
    let controls = band.carrier.controls();
    let row = |j: usize| &controls[2 * j];
    let binomial = [1, 5, 10, 10, 5, 1];
    let mut middle = Vector3::zero();
    for (j, weight) in binomial.iter().enumerate() {
        let scale = Real::from(Rational::from_ratio(*weight, 32).expect("nonzero"));
        middle = middle.plus(&row(j).minus(row(0)).scaled(&scale));
    }
    let middle = row(0).plus(&middle);
    let mut across = Vector3::zero();
    for (j, weight) in [1, 4, 6, 4, 1].iter().enumerate() {
        let scale = Real::from_i64(*weight);
        across = across.plus(&row(j + 1).minus(row(j)).scaled(&scale));
    }
    let station = &band.ball.carrier().controls()[2];
    let normal = across.cross(band.ball.along());
    match sign(&normal.dot(&station.minus(&middle)), dec)? {
        Sign::Positive => Ok(true),
        Sign::Negative => Ok(false),
        Sign::Zero => Err(bad(
            "a curvature-continuous section has no side towards its edge",
        )),
    }
}

pub(crate) fn build(
    model: &mut Model,
    body: BodyId,
    edges: &EdgeSelection,
    radius: &Real,
    fullness: &Real,
    step: &StepId,
    dec: &mut Decider,
) -> Result<BodyId, Refusal> {
    let uses = uses(model, body);
    let selected = resolve(model, edges, &uses, dec)?;
    if selected.is_empty() {
        return Err(unsupported("no input edge is selected"));
    }
    let planes = super::variable::convex_skin(model, body, unsupported, dec)?;
    let mut bands = Vec::<Band>::with_capacity(selected.len());
    for edge in selected {
        let [(a, sa), (b, sb)] = uses[&edge].as_slice() else {
            return Err(bad(
                "a curvature-continuous fillet edge needs two incident faces",
            ));
        };
        if !convex(model, edge, &uses[&edge], dec)? {
            return Err(unsupported(
                "a concave curvature-continuous band needs a nonconvex material-region proof",
            ));
        }
        let e = model.edge(edge);
        let (Some(start), Some(end)) = (e.start, e.end) else {
            return Err(unsupported(
                "a curvature-continuous band needs an open edge",
            ));
        };
        let ends = [start, end];
        let Curve::Line(_) = model.geom().curve(e.curve) else {
            return Err(unsupported(
                "a curvature-continuous band needs a straight edge",
            ));
        };
        let points = ends.map(|v| model.geom().position(model.vertex(v).point));
        let faces = [*a, *b];
        let supports = faces.map(|face| {
            planes
                .iter()
                .find(|(f, _, _)| *f == face)
                .expect("the convex skin lists every face")
                .1
                .clone()
        });
        let same_sense = faces.map(|face| model.face(face).sense == Sense::Same);
        let ball = LinearFillet::new(
            supports,
            same_sense,
            points,
            [radius.clone(), radius.clone()],
            true,
            dec,
        )?;
        let mut caps = [faces[0]; 2];
        for (end, old) in ends.iter().enumerate() {
            let incident: BTreeSet<_> = uses
                .iter()
                .filter(|(edge, _)| {
                    let e = model.edge(**edge);
                    e.start == Some(*old) || e.end == Some(*old)
                })
                .flat_map(|(_, faces)| faces.iter().map(|(face, _)| *face))
                .collect();
            if incident.len() != 3 {
                return Err(Refusal::new(
                    UNSUPPORTED_END,
                    "curvature-continuous band ends need a regular three-face input vertex",
                ));
            }
            caps[end] = *incident
                .iter()
                .find(|face| !faces.contains(face))
                .expect("a three-face vertex has a third face");
            let plane = &planes
                .iter()
                .find(|(face, _, _)| *face == caps[end])
                .expect("the convex skin lists every face")
                .1;
            if sign(&plane.normal.cross(ball.along()).norm_squared(), dec)? != Sign::Zero {
                return Err(Refusal::new(
                    UNSUPPORTED_END,
                    "an oblique curvature-continuous band end needs its end-section trimming construction",
                ));
            }
        }
        // Every section lies in its tangent triangle: its controls are the
        // two contacts and points on the segments from them to the station.
        // The triangle prism's six corners are the circular band's controls,
        // so the variable constructor's clearance proof covers the section.
        for (face, plane, normal) in &planes {
            for point in ball.carrier().controls() {
                let side = sign(&point.minus(&plane.origin).dot(normal), dec)?;
                if side == Sign::Positive
                    || (!faces.contains(face) && !caps.contains(face) && side != Sign::Negative)
                {
                    return Err(Refusal::new(
                        "fillet.unsupported-overflow",
                        "the curvature-continuous region reaches a neighbouring face and needs continuation",
                    ));
                }
            }
        }
        for previous in &bands {
            if !hulls_apart(
                ball.carrier().controls(),
                previous.ball.carrier().controls(),
                dec,
            )? {
                return Err(Refusal::new(
                    UNSUPPORTED_CORNER,
                    "curvature-continuous regions meet and need a junction construction",
                ));
            }
        }
        let name = model.edge_name(body, edge).ok_or_else(|| {
            bad("curvature-continuous fillet source edge has no unambiguous persistent name")
        })?;
        let carrier = carrier(&ball, fullness, dec)?;
        bands.push(Band {
            edge,
            ends,
            faces,
            senses: [*sa, *sb],
            caps,
            ball,
            carrier,
            origin: Origin {
                step: step.clone(),
                role: Role::Blend { edge: name },
            },
        });
    }
    let mut edit = Replacement::default();
    let mut end_changes = BTreeMap::new();
    let mut made = Vec::with_capacity(bands.len());
    for band in bands {
        let input = model.edge(band.edge).clone();
        let rows = band.ball.carrier().controls().to_vec();
        // feet[side][end]: rail 0 on the first face, rail 1 on the second.
        let mut feet = [[band.ends[0]; 2]; 2];
        for (side, rail) in feet.iter_mut().enumerate() {
            for (end, foot) in rail.iter_mut().enumerate() {
                *foot = vertex(model, rows[4 * side + end].clone());
            }
        }
        let mut path = Vec::with_capacity(4);
        for (side, rail) in feet.iter().enumerate() {
            let a = rows[4 * side].clone();
            let b = rows[4 * side + 1].clone();
            let curve = model.geom_mut().add_curve(Curve::Line(Line {
                dir: b.minus(&a),
                origin: a,
            }));
            let edge = model.add_edge(Edge {
                curve,
                sense: Sense::Same,
                start: Some(rail[0]),
                end: Some(rail[1]),
                tolerance: None,
            });
            edit.edges.insert((band.faces[side], band.edge), edge);
            path.push((edge, reverse(band.senses[side])));
        }
        for (end, old) in band.ends.iter().copied().enumerate() {
            let side = if end == 0 {
                SplineBoundary::VMin
            } else {
                SplineBoundary::VMax
            };
            let curve = model
                .geom_mut()
                .add_curve(Curve::BSpline(BSplineCurve::boundary(&band.carrier, side)));
            let cap = model.add_edge(Edge {
                curve,
                sense: Sense::Same,
                start: Some(feet[0][end]),
                end: Some(feet[1][end]),
                tolerance: None,
            });
            let first_finishes = if reverse(band.senses[0]) == Sense::Same {
                input.end == Some(old)
            } else {
                input.start == Some(old)
            };
            path.push((
                cap,
                if first_finishes {
                    Sense::Same
                } else {
                    Sense::Reversed
                },
            ));
            edit.gaps.insert((band.caps[end], old), Gap::Edge(cap));
            for (edge, adjacent) in &uses {
                if *edge == band.edge {
                    continue;
                }
                let e = model.edge(*edge);
                if e.start != Some(old) && e.end != Some(old) {
                    continue;
                }
                let side = usize::from(!adjacent.iter().any(|(face, _)| *face == band.faces[0]));
                end_changes.insert((*edge, old), feet[side][end]);
            }
        }
        let path = close(model, path)?;
        let sense = if outward_chart(&band, dec)? {
            Sense::Same
        } else {
            Sense::Reversed
        };
        let patch = RationalPatch::new(band.carrier.clone(), dec)?;
        let face = crate::extrude::face(
            model,
            Surface::BSpline(SplineSupport::RationalPatch(patch)),
            sense,
            vec![path],
            None,
            vec![band.origin],
        );
        edit.faces.push((band.faces[0], face));
        made.push(face);
    }
    for (old, incident) in &uses {
        let e = model.edge(*old).clone();
        let start = e
            .start
            .map(|v| end_changes.get(&(*old, v)).copied().unwrap_or(v));
        let end = e
            .end
            .map(|v| end_changes.get(&(*old, v)).copied().unwrap_or(v));
        if start == e.start && end == e.end {
            continue;
        }
        let edge = model.add_edge(Edge { start, end, ..e });
        for (face, _) in incident {
            edit.edges.insert((*face, *old), edge);
        }
    }
    let result = edge_blend::apply(model, body, &edit)?;
    for face in made {
        RationalPatchFace::of_face(model, face, dec)?
            .ok_or_else(|| bad("a curvature-continuous band lost its exact carrier"))?;
    }
    Ok(result)
}

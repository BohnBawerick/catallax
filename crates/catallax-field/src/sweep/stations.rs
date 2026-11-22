//! Independent coaxial loft sections in the saved uniform station coordinate.

use catallax_numbers::{Decider, Rational, Real, Refusal, Sign};
use catallax_recipe::{LoftSection, LoftShape, Recipe, Segment, StepId};
use catallax_shapes::spline_tools::RationalSplineSurface;
use catallax_shapes::{BSplineSurface, RPoint3, RVector3, RationalPatch};
use std::collections::{BTreeMap, BTreeSet};

use crate::{
    Plan,
    geometry::Leaf,
    patch_shell::{Face, PatchFace, Polygon, Shell},
    profile::{P, Profile, cross, dot, rational_vector, sub},
};

fn integer(n: i64) -> Rational {
    Rational::from_i64(n)
}

/// Natural cubic interpolation on the unit-spaced stations of RECIPE 4.12.
/// De Boor's second-derivative equations have diagonal 4 and off-diagonal 1;
/// the endpoint derivatives convert the solution to cubic Bernstein rows.
pub(super) fn station_controls(values: &[Rational], shape: LoftShape) -> Vec<Vec<Rational>> {
    if shape == LoftShape::Ruled || values.len() == 2 {
        return values.windows(2).map(<[_]>::to_vec).collect();
    }
    let count = values.len();
    let mut diagonal = vec![Rational::one(); count];
    let mut second = vec![Rational::zero(); count];
    for i in 1..count - 1 {
        let multiplier = if i == 1 {
            Rational::zero()
        } else {
            &Rational::one() / &diagonal[i - 1]
        };
        let rhs = ((&values[i + 1] - &values[i]) - (&values[i] - &values[i - 1])) * integer(6);
        diagonal[i] = &integer(4) - &multiplier;
        second[i] = rhs - &multiplier * &second[i - 1];
    }
    for i in (1..count - 1).rev() {
        second[i] = &(&second[i] - &second[i + 1]) / &diagonal[i];
    }
    (0..count - 1)
        .map(|i| {
            let a = &values[i];
            let b = &values[i + 1];
            let delta = (b - a) / integer(3);
            let m0 = &second[i];
            let m1 = &second[i + 1];
            vec![
                a.clone(),
                a + &delta - ((&(m0 * &integer(2)) + m1) / integer(18)),
                b - &delta - ((&(m1 * &integer(2)) + m0) / integer(18)),
                b.clone(),
            ]
        })
        .collect()
}

pub(super) fn knots(n: usize) -> Vec<Rational> {
    std::iter::repeat_n(Rational::zero(), n)
        .chain(std::iter::repeat_n(Rational::one(), n))
        .collect()
}

// Every section uses the same rational quarter-circle parameter. It is the
// recipe's equal-angle chart after applying one common Möbius map to u.
const QUARTERS: [[[i64; 2]; 3]; 4] = [
    [[1, 0], [1, 1], [0, 1]],
    [[0, 1], [-1, 1], [-1, 0]],
    [[-1, 0], [-1, -1], [0, -1]],
    [[0, -1], [1, -1], [1, 0]],
];

// Describe the authored chart cuts in coordinates relative to the saved start.
// Equal cuts and turns give all sections a common reparameterization by the
// quarters above, even when their starting directions differ.
fn circular_phase(
    recipe: &Recipe,
    section: &LoftSection,
    centre: &P,
    squared_radius: &Rational,
) -> Option<(P, Vec<P>)> {
    let [loop_] = super::authored_loops(recipe, &section.profile) else {
        return None;
    };
    let segments = &loop_.0;
    let start = segments.iter().position(|s| s.id() == &section.start[0])?;
    if let [Segment::Circle { radius, .. }] = segments.as_slice() {
        let cuts = QUARTERS.map(|[_, _, end]| end.map(integer)).into();
        return Some(([radius.value().clone(), Rational::zero()], cuts));
    }
    let endpoint = |i: usize| match &segments[i] {
        Segment::Arc { to, .. } => Some(sub(&to.0.each_ref().map(|n| n.value().clone()), centre)),
        _ => None,
    };
    let phase = endpoint((start + segments.len() - 1) % segments.len())?;
    let cuts = (start..start + segments.len())
        .map(|i| {
            let end = endpoint(i % segments.len())?;
            Some([
                &dot(&phase, &end) / squared_radius,
                &cross(&phase, &end) / squared_radius,
            ])
        })
        .collect::<Option<Vec<_>>>()?;
    Some((phase, cuts))
}

fn quarter_carrier(
    width: usize,
    positions: Vec<RPoint3>,
    dec: &Decider,
) -> Result<RationalSplineSurface, Refusal> {
    let half = Rational::one().mul_pow2(-1);
    let weights = [Rational::one(), half.clone(), half]
        .into_iter()
        .flat_map(|w| std::iter::repeat_n(w, width))
        .collect();
    RationalSplineSurface::new(
        [2, width - 1],
        [knots(3), knots(width)],
        [3, width],
        positions.into_iter().map(|p| [p.x, p.y, p.z]).collect(),
        weights,
        dec,
    )
}

fn quarter_patch(
    carrier: &RationalSplineSurface,
    outward: bool,
    dec: &mut Decider,
) -> Result<Face, Refusal> {
    let patch = RationalPatch::new(carrier.to_surface(dec)?, dec)?;
    Ok(Face::Patch(Box::new(PatchFace::new(&patch, outward, dec)?)))
}

impl Plan {
    pub(super) fn polygon_station_loft(
        &mut self,
        recipe: &Recipe,
        sections: &[LoftSection],
        shape: LoftShape,
        profiles: &BTreeMap<StepId, Profile>,
        dec: &mut Decider,
    ) -> Result<Option<usize>, Refusal> {
        let first = &profiles[&sections[0].profile];
        let inverse = first.pose.inverse();
        let mut heights = Vec::new();
        let mut sections_by_parameter: Vec<Vec<(Rational, P)>> = Vec::new();
        let mut breaks = BTreeSet::new();
        for section in sections {
            let profile = &profiles[&section.profile];
            let [edges] = super::authored_loops(recipe, &section.profile) else {
                return Ok(None);
            };
            let edges = &edges.0;
            if edges.len() < 3
                || edges.iter().any(|e| !matches!(e, Segment::Line { .. }))
                || !profile
                    .pose
                    .rotation
                    .same_rotation(&first.pose.rotation, dec)?
            {
                return Ok(None);
            }
            let Some(start) = edges.iter().position(|e| e.id() == &section.start[0]) else {
                return Ok(None);
            };
            let offset = profile
                .pose
                .then(&inverse)
                .translation
                .as_rational()
                .expect("rational sketch placements");
            heights.push(offset.z);
            let polygon = (0..edges.len())
                .map(|i| {
                    let Segment::Line { to, .. } =
                        &edges[(start + i + edges.len() - 1) % edges.len()]
                    else {
                        unreachable!("polygon section")
                    };
                    let parameter = Rational::from_ratio(i as i64, edges.len() as i64)
                        .expect("nonempty section");
                    breaks.insert(parameter.clone());
                    (
                        parameter,
                        [to.0[0].value() + &offset.x, to.0[1].value() + &offset.y],
                    )
                })
                .collect::<Vec<_>>();
            sections_by_parameter.push(polygon);
        }

        // Refine every authored loop on the union of its rational line
        // chart breaks. Merged collinear field edges do not change u.
        let mut polygons = Vec::new();
        let one = Rational::one();
        for section in &sections_by_parameter {
            let mut edge = 0;
            let mut polygon = Vec::with_capacity(breaks.len());
            for t in &breaks {
                dec.charge_work(1)?;
                if section.get(edge + 1).is_some_and(|(end, _)| t == end) {
                    edge += 1;
                }
                let (start, a) = &section[edge];
                let (end, b) = section
                    .get(edge + 1)
                    .map_or((&one, &section[0].1), |(u, p)| (u, p));
                let local = &(t - start) / &(end - start);
                polygon.push(std::array::from_fn(|k| {
                    &a[k] + &(&local * &(&b[k] - &a[k]))
                }));
            }
            polygons.push(polygon);
        }
        let count = breaks.len();
        let height_rows = station_controls(&heights, shape);
        let up = heights.last().expect("checked station count").sign();
        if up == Sign::Zero {
            return Ok(None);
        }
        let width = height_rows[0].len();
        let vertex_rows: Vec<_> = (0..count)
            .map(|vertex| {
                [0, 1].map(|axis| {
                    let values: Vec<_> = polygons.iter().map(|p| p[vertex][axis].clone()).collect();
                    station_controls(&values, shape)
                })
            })
            .collect();
        let mut faces = Vec::new();
        let mut winding = None;
        let world = |p: &[Rational; 2], h: &Rational| {
            crate::profile::rational_point(
                &first.pose,
                &RPoint3::new(p[0].clone(), p[1].clone(), h.clone()),
            )
        };
        for (span, z) in height_rows.iter().enumerate() {
            if z.windows(2).any(|p| (&p[1] - &p[0]).sign() != up) {
                return Ok(None);
            }
            let vertex =
                |i: usize, j: usize| vertex_rows[i].each_ref().map(|axis| axis[span][j].clone());
            // Cross products of every pair of incident edge controls
            // bound the turn by bilinearity and the positive Bernstein basis.
            for i in 0..count {
                for j in 0..width {
                    let incoming = sub(&vertex(i, j), &vertex((i + count - 1) % count, j));
                    for k in 0..width {
                        dec.charge_work(1)?;
                        let outgoing = sub(&vertex((i + 1) % count, k), &vertex(i, k));
                        let sign = cross(&incoming, &outgoing).sign();
                        if sign == Sign::Zero {
                            if dot(&incoming, &outgoing).sign() != Sign::Positive {
                                return Ok(None);
                            }
                        } else if sign != *winding.get_or_insert(sign) {
                            return Ok(None);
                        }
                    }
                }
            }
            let Some(winding) = winding else {
                return Ok(None);
            };
            let outward = winding == up;
            for i in 0..count {
                let points = [i, (i + 1) % count]
                    .into_iter()
                    .flat_map(|v| {
                        z.iter()
                            .enumerate()
                            .map(move |(j, h)| world(&vertex(v, j), h).to_real())
                    })
                    .collect();
                let carrier = BSplineSurface::new(
                    [1, width - 1],
                    [knots(2), knots(width)],
                    [2, width],
                    points,
                    vec![Real::one(); 2 * width],
                    dec,
                )?;
                let patch = RationalPatch::new(carrier, dec)?;
                faces.push(Face::Patch(Box::new(PatchFace::new(&patch, outward, dec)?)));
            }
        }
        let same = winding == Some(up);
        for (station, reverse) in [(0, same), (sections.len() - 1, !same)] {
            let mut cap: Vec<_> = polygons[station]
                .iter()
                .map(|p| world(p, &heights[station]))
                .collect();
            if reverse {
                cap.reverse();
            }
            faces.push(Face::Polygon(Box::new(Polygon::new(cap)?)));
        }
        let shell = Shell::new(faces)?;
        self.leaf(Leaf::Shell(Box::new(shell)), vec![], dec)
            .map(Some)
    }

    pub(super) fn circular_station_loft(
        &mut self,
        recipe: &Recipe,
        sections: &[LoftSection],
        shape: LoftShape,
        profiles: &BTreeMap<StepId, Profile>,
        dec: &mut Decider,
    ) -> Result<Option<usize>, Refusal> {
        let first = &profiles[&sections[0].profile];
        let basis = [[1, 0, 0], [0, 1, 0], [0, 0, 1]].map(|[x, y, z]| {
            rational_vector(
                &first.pose,
                &RVector3::new(integer(x), integer(y), integer(z)),
            )
        });
        let [x, y, normal] = basis;
        let mut centre: Option<RPoint3> = None;
        let mut phases = Vec::new();
        let mut chart_law = None;
        let mut heights = Vec::new();
        for section in sections {
            let profile = &profiles[&section.profile];
            let [loop_edges] = profile.loops.as_slice() else {
                return Ok(None);
            };
            let [edge] = loop_edges.as_slice() else {
                return Ok(None);
            };
            let Some((local, squared, turn)) = &edge.circle else {
                return Ok(None);
            };
            if !edge.full
                || !profile
                    .pose
                    .rotation
                    .same_rotation(&first.pose.rotation, dec)?
            {
                return Ok(None);
            }
            let Some((phase, cuts)) = circular_phase(recipe, section, local, squared) else {
                return Ok(None);
            };
            // Antipodal cuts do not distinguish clockwise and
            // counterclockwise traversal.
            let law = (*turn, cuts);
            if chart_law.get_or_insert_with(|| law.clone()) != &law {
                return Ok(None);
            }
            let current = profile.world(local);
            let origin = centre.get_or_insert(current.clone());
            let delta = current.minus(origin);
            if !delta.cross(&normal).is_zero() {
                return Ok(None);
            }
            heights.push(delta.dot(&normal));
            phases.push(phase);
        }

        let up = heights.last().expect("checked stations").sign();
        let [phase_x, phase_y] = [0, 1].map(|axis| {
            let values: Vec<_> = phases.iter().map(|phase| phase[axis].clone()).collect();
            station_controls(&values, shape)
        });
        let height_rows = station_controls(&heights, shape);
        // Strictly ordered height controls keep the sections disjoint over
        // the full cubic spans. A common open half-plane for the Bernstein
        // phase controls keeps every interpolated section away from zero.
        if up == Sign::Zero
            || height_rows
                .iter()
                .any(|row| row.windows(2).any(|p| (&p[1] - &p[0]).sign() != up))
        {
            return Ok(None);
        }
        let centre = centre.expect("checked stations");
        let world = |[u, v]: [i64; 2], phase: &P, height: &Rational| {
            let [u, v] = [u, v].map(integer);
            centre
                .plus(&normal.scaled(height))
                .plus(&x.scaled(&(&phase[0] * &u - &phase[1] * &v)))
                .plus(&y.scaled(&(&phase[1] * &u + &phase[0] * &v)))
        };
        let mut faces = Vec::new();
        for (span, heights) in height_rows.iter().enumerate() {
            let controls: Vec<P> = phase_x[span]
                .iter()
                .zip(&phase_y[span])
                .map(|(x, y)| [x.clone(), y.clone()])
                .collect();
            let across = std::array::from_fn(|i| &phases[span][i] + &phases[span + 1][i]);
            if controls
                .iter()
                .any(|p| dot(p, &across).sign() != Sign::Positive)
            {
                return Ok(None);
            }
            for quarter in QUARTERS {
                let positions = quarter
                    .into_iter()
                    .flat_map(|axes| {
                        controls
                            .iter()
                            .zip(heights)
                            .map(move |(phase, h)| world(axes, phase, h))
                    })
                    .collect();
                let carrier = quarter_carrier(controls.len(), positions, dec)?;
                let outward = up == Sign::Positive;
                match quarter_patch(&carrier, outward, dec) {
                    Ok(face) => faces.push(face),
                    Err(e) if e.code == "shapes.rational-patch-unproven" => {
                        // A rotating quarter can span a half turn in space.
                        // Exact knot insertion gives narrower projections
                        // while retaining the same correspondence and seams.
                        let mid = Rational::one().mul_pow2(-1);
                        for column in carrier.split(0, &mid, dec)? {
                            for piece in column.split(1, &mid, dec)? {
                                faces.push(quarter_patch(&piece, outward, dec)?);
                            }
                        }
                    }
                    Err(e) => return Err(e),
                }
            }
        }
        for (station, positive_normal) in [
            (0, up == Sign::Negative),
            (sections.len() - 1, up == Sign::Positive),
        ] {
            let phase = &phases[station];
            let h = &heights[station];
            // A central diamond and four annular pieces partition the disk.
            // Each piece interpolates between the chord and the quarter circle
            // on a common positive denominator.
            let inner = phase.each_ref().map(|v| v.mul_pow2(-1));
            let mut centre_cap = Vec::with_capacity(4);
            for [a, b, c] in QUARTERS {
                let first = world(a, &inner, h);
                let last = world(c, &inner, h);
                let midpoint =
                    first.plus(&last.minus(&first).scaled(&Rational::one().mul_pow2(-1)));
                centre_cap.push(first.clone());
                let carrier = quarter_carrier(
                    2,
                    vec![
                        first,
                        world(a, phase, h),
                        midpoint,
                        world(b, phase, h),
                        last,
                        world(c, phase, h),
                    ],
                    dec,
                )?;
                faces.push(quarter_patch(&carrier, !positive_normal, dec)?);
            }
            if !positive_normal {
                centre_cap.reverse();
            }
            faces.push(Face::Polygon(Box::new(Polygon::new(centre_cap)?)));
        }
        let shell = Shell::new(faces)?;
        self.leaf(Leaf::Shell(Box::new(shell)), vec![], dec)
            .map(Some)
    }
}

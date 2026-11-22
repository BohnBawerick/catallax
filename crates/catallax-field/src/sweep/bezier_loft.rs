//! Parallel lofts whose sections have a polynomial Bézier roof and a chord.

use catallax_numbers::{Decider, Rational, Refusal, Sign};
use catallax_recipe::{LoftParameterization, LoftSection, LoftShape, Operation, Recipe, Segment};
use catallax_shapes::spline_tools::{RationalPoint, RationalSplineCurve, RationalSplineSurface};
use catallax_shapes::{Pose, RPoint3, RationalPatch};

use crate::Plan;
use crate::geometry::Leaf;
use crate::patch_shell::{Face, PatchFace, Shell};
use crate::profile::{P, Profile, cross, dot, rational_point, sub};

use super::stations::{knots, station_controls};

fn unsupported() -> Refusal {
    Refusal::new(
        "field.loft-bezier-correspondence",
        "requires parallel sections with a corresponding Bézier roof over a shared straight chord",
    )
}

fn patch(
    surface: &RationalSplineSurface,
    outward: bool,
    dec: &mut Decider,
) -> Result<Face, Refusal> {
    let patch = RationalPatch::new(surface.to_surface(dec)?, dec)?;
    Ok(Face::Patch(Box::new(PatchFace::new(&patch, outward, dec)?)))
}

fn curve(controls: Vec<RationalPoint>, dec: &Decider) -> Result<RationalSplineCurve, Refusal> {
    let count = controls.len();
    RationalSplineCurve::new(
        count - 1,
        knots(count),
        controls,
        vec![Rational::one(); count],
        dec,
    )
}

fn cap(
    controls: Vec<RationalPoint>,
    outward: bool,
    dec: &mut Decider,
) -> Result<Vec<Face>, Refusal> {
    let start = controls[0].clone();
    let end = controls.last().expect("Bezier endpoint").clone();
    let half = Rational::one().mul_pow2(-1);
    let roof = curve(controls, dec)?;
    let mid = roof.point_at(&half, dec)?;
    let from = RPoint3::new(start[0].clone(), start[1].clone(), start[2].clone());
    let to = RPoint3::new(end[0].clone(), end[1].clone(), end[2].clone());
    let peak = RPoint3::new(mid[0].clone(), mid[1].clone(), mid[2].clone());
    let chord = to.minus(&from);
    let foot = from.plus(&chord.scaled(&(peak.minus(&from).dot(&chord) / chord.norm_squared())));
    let floor = [foot.x, foot.y, foot.z];
    let [left, right] = roof.split(&half, dec)?;
    let mut faces = Vec::new();
    for (arc, a, b) in [(left, start, mid.clone()), (right, mid, end)] {
        let centre =
            std::array::from_fn(|i| (&a[i] + &b[i] + floor[i].clone()) / Rational::from_i64(3));
        let boundary = [
            arc,
            curve(vec![b, floor.clone()], dec)?,
            curve(vec![floor.clone(), a], dec)?,
        ];
        // Splitting roof and chord at the same transverse line produces two
        // convex curved triangles. The N-sided Coons charts share straight
        // spokes, and each chart separately proves its planar embedding.
        for surface in RationalSplineSurface::n_sided_fill(&boundary, &centre, dec)? {
            faces.push(patch(&surface, outward, dec)?);
        }
    }
    Ok(faces)
}

impl Plan {
    pub(crate) fn bezier_loft(
        &mut self,
        recipe: &Recipe,
        sections: &[LoftSection],
        shape: LoftShape,
        dec: &mut Decider,
    ) -> Result<usize, Refusal> {
        if recipe.loft_parameterization().ok() != Some(LoftParameterization::RationalCharts) {
            return Err(Refusal::new(
                "loft.legacy-parameterization",
                "legacy arc-length lofts are not rational-chart lofts",
            ));
        }
        let mut base: Option<Pose> = None;
        let mut start_at_roof = None;
        let mut profiles: Vec<Vec<P>> = Vec::new();
        let mut heights = Vec::new();
        for section in sections {
            let Operation::Sketch { plane, loops } =
                &recipe.step(&section.profile).expect("checked profile").op
            else {
                unreachable!("checked loft profile")
            };
            let [boundary] = loops.as_slice() else {
                return Err(unsupported());
            };
            let (bezier, line) = match boundary.0.as_slice() {
                [bezier @ Segment::Bezier { .. }, line @ Segment::Line { .. }]
                | [line @ Segment::Line { .. }, bezier @ Segment::Bezier { .. }] => (bezier, line),
                _ => return Err(unsupported()),
            };
            let (Segment::Bezier { controls, to, .. }, Segment::Line { to: start, .. }) =
                (bezier, line)
            else {
                unreachable!("matched curve and chord")
            };
            let starts_at_roof = &section.start[0] == bezier.id();
            if *start_at_roof.get_or_insert(starts_at_roof) != starts_at_roof {
                return Err(unsupported());
            }
            let pose = Profile::new(plane, &[], dec)?.pose;
            let initial = base.get_or_insert_with(|| pose.clone());
            if !initial.rotation.same_rotation(&pose.rotation, dec)? {
                return Err(unsupported());
            }
            let shift = pose
                .then(&initial.inverse())
                .translation
                .as_rational()
                .expect("rational frames");
            let net: Vec<_> = std::iter::once(start)
                .chain(controls)
                .chain(std::iter::once(to))
                .map(|p| [p.0[0].value() + &shift.x, p.0[1].value() + &shift.y])
                .collect();
            if profiles.first().is_some_and(|first| {
                first.len() != net.len() || first[0] != net[0] || first.last() != net.last()
            }) {
                return Err(unsupported());
            }
            heights.push(shift.z);
            profiles.push(net);
        }
        let count = profiles[0].len();
        let a = &profiles[0][0];
        let b = &profiles[0][count - 1];
        let chord = sub(b, a);
        let sense = cross(&chord, &sub(&profiles[0][1], a)).sign();
        let h_sign = heights.last().expect("checked sections").sign();
        if heights.iter().all(Rational::is_zero) {
            return Ok(0);
        }
        if sense == Sign::Zero || h_sign == Sign::Zero {
            return Err(unsupported());
        }
        let rows: Vec<_> = (0..count)
            .map(|i| {
                [0, 1].map(|k| {
                    let values: Vec<_> = profiles.iter().map(|p| p[i][k].clone()).collect();
                    station_controls(&values, shape)
                })
            })
            .collect();
        let frame = base.expect("checked sections");
        let world = |p: &P, h: &Rational| {
            let p = rational_point(&frame, &RPoint3::new(p[0].clone(), p[1].clone(), h.clone()));
            [p.x, p.y, p.z]
        };
        let outward = sense != h_sign;
        let mut faces = Vec::new();
        let half = Rational::one().mul_pow2(-1);
        for (span, z) in station_controls(&heights, shape).iter().enumerate() {
            if z.windows(2)
                .any(|pair| (&pair[1] - &pair[0]).sign() != h_sign)
            {
                return Err(unsupported());
            }
            let width = z.len();
            let at = |u: usize, v: usize| rows[u].each_ref().map(|axis| axis[span][v].clone());
            for v in 0..width {
                // Bernstein positivity keeps every interpolated roof a graph
                // over the chord with its interior strictly on one side.
                for u in 0..count - 1 {
                    dec.charge_work(1)?;
                    let going = sub(&at(u + 1, v), &at(u, v));
                    if dot(&chord, &going).sign() == Sign::Negative
                        || (u + 1 < count - 1
                            && cross(&chord, &sub(&at(u + 1, v), a)).sign() != sense)
                    {
                        return Err(unsupported());
                    }
                }
            }
            let controls = (0..count)
                .flat_map(|u| z.iter().enumerate().map(move |(v, h)| world(&at(u, v), h)))
                .collect();
            let roof = RationalSplineSurface::new(
                [count - 1, width - 1],
                [knots(count), knots(width)],
                [count, width],
                controls,
                vec![Rational::one(); count * width],
                dec,
            )?;
            // Splitting at the interior midpoint lets both halves use their
            // own embedding projection when an endpoint tangent is vertical.
            for half_roof in roof.split(0, &half, dec)? {
                faces.push(patch(&half_roof, outward, dec)?);
            }
            let floor = RationalSplineSurface::new(
                [1, 1],
                [knots(2), knots(2)],
                [2, 2],
                [b, a]
                    .into_iter()
                    .flat_map(|p| [world(p, &z[0]), world(p, &z[width - 1])])
                    .collect(),
                vec![Rational::one(); 4],
                dec,
            )?;
            faces.push(patch(&floor, outward, dec)?);
        }
        for (station, outward) in [(0, !outward), (sections.len() - 1, outward)] {
            let p = &profiles[station];
            // The convex control polygon proves a convex roof and chord
            // region by the Bézier variation-diminishing property.
            for i in 0..count {
                let turn = cross(
                    &sub(&p[(i + 1) % count], &p[i]),
                    &sub(&p[(i + 2) % count], &p[(i + 1) % count]),
                );
                if turn.sign() == sense {
                    return Err(unsupported());
                }
            }
            faces.extend(cap(
                p.iter().map(|p| world(p, &heights[station])).collect(),
                outward,
                dec,
            )?);
        }
        self.leaf(Leaf::Shell(Box::new(Shell::new(faces)?)), vec![], dec)
    }
}

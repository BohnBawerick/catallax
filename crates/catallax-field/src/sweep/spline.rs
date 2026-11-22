//! Fixed-frame translations along a rational spline.

use catallax_numbers::{Decider, Rational, Real, Refusal, Sign};
use catallax_recipe::PathSegment;
use catallax_shapes::{BSplineCurve, BSplineSurface, Point3, RPoint3, RVector3, RationalPatch};

use crate::geometry::Leaf;
use crate::patch_shell::{Face, PatchFace, Polygon, Shell};
use crate::profile::{Profile, cross, loop_location, rational_vector, real};
use crate::{Location, Plan};

impl Plan {
    pub(super) fn spline_sweep(
        &mut self,
        profile: &Profile,
        segment: &PathSegment,
        dec: &mut Decider,
    ) -> Result<usize, Refusal> {
        let PathSegment::Spline {
            degree,
            control_points,
            weights,
            knots,
            ..
        } = segment
        else {
            unreachable!("fixed spline dispatch")
        };
        if profile.loops.is_empty() {
            return Ok(0);
        }
        if profile
            .loops
            .iter()
            .flatten()
            .any(|edge| edge.circle.is_some())
        {
            return Err(Refusal::new(
                "field.spline-sweep-profile",
                "spline sweep requires a polygonal profile",
            ));
        }

        let degree = *degree;
        let controls: Vec<_> = control_points.iter().map(super::point).collect();
        let origin = &controls[0];
        let end = controls.last().expect("checked spline").minus(origin);
        let curve = BSplineCurve::new(
            degree,
            knots.iter().map(|n| n.value().clone()).collect(),
            controls.iter().map(RPoint3::to_real).collect(),
            weights
                .iter()
                .map(|n| Real::from(n.value().clone()))
                .collect(),
            dec,
        )?;
        let normal = rational_vector(
            &profile.pose,
            &RVector3::new(Rational::zero(), Rational::zero(), Rational::one()),
        );
        let travel = end.dot(&normal).sign();
        let mut faces = Vec::new();
        let mut rings = Vec::new();
        for (i, edges) in profile.loops.iter().enumerate() {
            let signed_area = edges.iter().fold(Rational::zero(), |area, edge| {
                area + cross(&edge.start, &edge.end)
            });
            let mut hole = false;
            let witness = real(&edges[0].start);
            for (j, other) in profile.loops.iter().enumerate() {
                if i != j && loop_location(other, &witness, dec)? == Location::Inside {
                    hole = !hole;
                }
            }
            let positive = signed_area.sign() == travel;
            let vertices: Vec<_> = edges
                .iter()
                .map(|edge| profile.world(&edge.start))
                .collect();
            rings.push((vertices, positive, hole));
        }
        let line_knots = vec![
            Rational::zero(),
            Rational::zero(),
            Rational::one(),
            Rational::one(),
        ];
        let span_knots = vec![Rational::zero(); degree + 1]
            .into_iter()
            .chain(vec![Rational::one(); degree + 1])
            .collect::<Vec<_>>();

        // For positive weights, the derivative of the projected rational
        // Bézier curve has the sign of the sum of pairwise ordered control
        // differences times nonnegative Bernstein products. Strict ordering
        // keeps the transported copies in distinct parallel planes.
        for (span, _) in curve.spans() {
            let span_controls = curve.span_bezier_controls(span, dec)?;
            let mut shifts = Vec::with_capacity(span_controls.len());
            let mut weights = Vec::with_capacity(span_controls.len());
            for [x, y, z, w] in span_controls {
                let w = w.rational_ref().expect("rational spline weight");
                let coordinate = |q: &Real| q.rational_ref().expect("rational spline control") / w;
                let p = RPoint3::new(coordinate(&x), coordinate(&y), coordinate(&z));
                shifts.push(p.minus(origin));
                weights.push(Real::from(w.clone()));
            }
            if travel == Sign::Zero
                || shifts
                    .windows(2)
                    .any(|pair| pair[1].minus(&pair[0]).dot(&normal).sign() != travel)
            {
                return Err(Refusal::new(
                    "sweep.spline-admissibility-unproven",
                    format!("normal projections of spline span {span} are not strictly ordered"),
                ));
            }

            // SMLib N_CreateTransSweepSrf: the side is a tensor product
            // between the profile edge and the exact translation curve.
            for (ring, positive, hole) in &rings {
                for (a, b) in ring.iter().zip(ring.iter().cycle().skip(1)) {
                    let points: Vec<Point3> = [a, b]
                        .into_iter()
                        .flat_map(|p| shifts.iter().map(|shift| p.plus(shift).to_real()))
                        .collect();
                    let carrier = BSplineSurface::new(
                        [1, degree],
                        [line_knots.clone(), span_knots.clone()],
                        [2, degree + 1],
                        points,
                        [weights.as_slice(), weights.as_slice()].concat(),
                        dec,
                    )?;
                    let patch = RationalPatch::new(carrier, dec)?;
                    faces.push(Face::Patch(Box::new(PatchFace::new(
                        &patch,
                        positive != hole,
                        dec,
                    )?)));
                }
            }
        }
        // The cap uses the first ring for its normal. Use a material ring
        // there and orient every ring by nesting, regardless of authored order.
        let outer = rings
            .iter()
            .position(|(_, _, hole)| !hole)
            .expect("validated profile");
        rings.swap(0, outer);
        for at_end in [false, true] {
            let mut cap = rings.iter().map(|(ring, positive, hole)| {
                let mut vertices = if at_end {
                    ring.iter().map(|p| p.plus(&end)).collect()
                } else {
                    ring.clone()
                };
                if (*positive != *hole) != at_end {
                    vertices.reverse();
                }
                vertices
            });
            let first = cap.next().expect("nonempty profile");
            faces.push(Face::Polygon(Box::new(Polygon::with_rings(
                first,
                cap.collect(),
            )?)));
        }
        let shell = Shell::new(faces)?;
        self.leaf(Leaf::Shell(Box::new(shell)), vec![], dec)
    }
}

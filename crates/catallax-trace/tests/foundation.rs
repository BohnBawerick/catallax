#![allow(missing_docs)]
//! Differential gate for the port's foundation and step 1 (`SmGeomUtility`)
//! against saved outputs of the pinned C++ SMLib. See
//! `tests/fixtures/foundation/README.md` for provenance and regeneration.

mod common;

use catallax_trace::Fault;
use catallax_trace::extent::{Extent1, Extent2, Extent3};
use catallax_trace::geometry::{
    CurveJet, PlanePair, determinant_3_vectors, evaluate_curve_continuity, frenet_from_derivatives,
    intersect_three_planes, intersect_two_planes, line_closest_point,
    line_coplanar_circle_intersect, line_point_distance, plane_point_distance,
    point_project_to_plane, segment_point_angle_deg_to_perp, segment_point_distance,
    triangle_area_squared, triangle_point_distance, vec_lin_comb_two_vectors,
    vector_project_to_plane,
};
use catallax_trace::vector::Vector3;
use common::{Input, Output, run_gate};

fn ext2_box(out: &mut Output, e: Extent2) {
    out.v2(e.min);
    out.v2(e.max);
}

fn ext3_box(out: &mut Output, e: Extent3) {
    out.v3(e.min);
    out.v3(e.max);
}

fn evaluate(op: &str, input: &mut Input, out: &mut Output) {
    match op {
        "angle" => {
            let (a, b) = (input.v3(), input.v3());
            if let Some(angle) = out.status(a.angle_between(b)) {
                out.d(angle);
            }
        }
        "unitize" => {
            let mut a = input.v3();
            let status = a.unitize();
            out.status(status);
            out.v3(a);
        }
        "parallel" => {
            let (a, b, degrees) = (input.v3(), input.v3(), input.d());
            let (parallel, anti) = a.is_parallel_to(b, degrees);
            out.flag(parallel);
            out.i(anti.map_or(2, i64::from));
        }
        "perp" => {
            let (a, b, degrees) = (input.v3(), input.v3(), input.d());
            out.flag(a.is_perpendicular_to(b, degrees));
        }
        "ortho" => {
            let a = input.v3();
            let has = input.i() != 0;
            let reference = input.v3();
            if let Some([x, y, z]) = out.status(a.make_unit_ortho_vectors(has.then_some(reference)))
            {
                out.v3(x);
                out.v3(y);
                out.v3(z);
            }
        }
        "uderiv" => {
            let (a, d) = (input.v3(), input.v3());
            out.v3(a.unitized_derivative(d));
        }
        "divassign3" => {
            let mut a = input.v3();
            a /= input.d();
            out.v3(a);
        }
        "div2" => {
            let a = input.v2();
            out.v2(a / input.d());
        }
        "closer" => {
            let (a, b, maximum) = (input.v3(), input.v3(), input.d());
            out.flag(a.closer_than(maximum, b));
        }
        "differs2" => {
            let (a, b) = (input.v2(), input.v2());
            out.flag(a.differs(b));
        }
        "ext1" => {
            let [lo, hi, x, tol, olo, ohi, u] = std::array::from_fn(|_| input.d());
            let (e, o) = (Extent1::new(lo, hi), Extent1::new(olo, ohi));
            out.d(e.clamp_value(x));
            out.d(e.evaluate(u));
            out.i(e.classify_point(x, tol) as i64);
            let (on, other) = e.is_value_on_boundary(x, tol);
            out.flag(on);
            out.d(other);
            out.flag(e.are_disjoint(o, tol));
            out.flag(e.is_contained_by(o, tol));
            out.flag(e.contains_value(x, tol));
            out.flag(e.are_equal(o, tol));
            out.flag(e.equals_ends(olo, ohi, tol));
            out.flag(e.is_bounded());
            out.d(e.mid());
            out.d(e.length());
            let mut grown = e;
            grown.expand_absolute(tol);
            out.d(grown.min);
            out.d(grown.max);
            out.status(Extent1::default().set_min_max(olo, ohi));
        }
        "ext2" => {
            let (lo, hi, p, tol, olo, ohi) = (
                input.v2(),
                input.v2(),
                input.v2(),
                input.d(),
                input.v2(),
                input.v2(),
            );
            let (u, v) = (input.d(), input.d());
            let (e, o) = (Extent2::new(lo, hi), Extent2::new(olo, ohi));
            out.flag(e.contains_point(p, tol));
            let (inside, distance) = e.contains_point_distance(p, tol);
            out.flag(inside);
            out.d(distance);
            let (on, normal) = e.point_on_boundary(p, tol);
            out.flag(on);
            out.v2(normal);
            out.i(i64::from(e.point_boundaries(p, tol)));
            out.v2(e.clamp_point(p));
            out.v2(e.evaluate(u, v));
            ext2_box(out, e.union(o));
            out.flag(e.are_disjoint(o, tol));
            out.flag(e.is_contained_by(o, tol));
            out.v2(e.size());
            out.d(e.x_length());
            out.d(e.y_length());
            for interval in [e.u_interval(), e.v_interval()] {
                out.d(interval.min);
                out.d(interval.max);
            }
            if let Some(common) = out.status(e.intersect(o)) {
                ext2_box(out, common);
            }
        }
        "clip" => {
            let (lo, hi, point, vector) = (input.v2(), input.v2(), input.v2(), input.v2());
            out.d(Extent2::new(lo, hi).clip_line(point, vector, input.d()));
        }
        "ext3" => {
            let raw = input.i() != 0;
            let (lo, hi, p, tol) = (input.v3(), input.v3(), input.v3(), input.d());
            let (olo, ohi, n) = (input.v3(), input.v3(), input.v3());
            let e = if raw {
                Extent3::from_corners(lo, hi)
            } else {
                Extent3::new(lo, hi)
            };
            let o = Extent3::new(olo, ohi);
            out.flag(e.contains_point(p, tol));
            out.flag(e.are_disjoint(o, tol));
            ext3_box(out, e.union(o));
            out.flag(e.is_point_sized(tol));
            out.flag(e.has_negative_volume());
            let (center, radius, inner) = e.sphere_bound();
            out.v3(center);
            out.d(radius);
            out.d(inner);
            let mut grown = e;
            grown.expand_absolute(tol);
            ext3_box(out, grown);
            out.v3(e.evaluate(n));
            out.v3(e.size());
            out.d(e.max_dimension());
            if let Some((axis, angle)) = out.status(e.cone_bound_to(o)) {
                out.v3(axis);
                out.d(angle);
            }
            if let Some(common) = out.status(e.intersect(o)) {
                ext3_box(out, common);
            }
        }
        "triarea" => {
            let (a, b, c) = (input.v3(), input.v3(), input.v3());
            out.status(Ok(()));
            out.d(triangle_area_squared(a, b, c));
        }
        "tripoint" => {
            let (a, b, c, t) = (input.v3(), input.v3(), input.v3(), input.v3());
            let found = triangle_point_distance(a, b, c, t);
            out.status(Ok(()));
            out.d(found.distance);
            out.v3(found.closest);
            out.d(found.u);
            out.d(found.v);
        }
        "linepoint" => {
            let (p, d, t) = (input.v3(), input.v3(), input.v3());
            let (distance, parameter) = line_point_distance(p, d, t);
            out.status(Ok(()));
            out.d(distance);
            out.d(parameter);
        }
        "lineclosest" => {
            let (p, d, t) = (input.v3(), input.v3(), input.v3());
            out.status(Ok(()));
            out.d(line_closest_point(p, d, t));
        }
        "segpoint" => {
            let (s, e, t) = (input.v3(), input.v3(), input.v3());
            let (distance, parameter) = segment_point_distance(s, e, t);
            out.status(Ok(()));
            out.d(distance);
            out.d(parameter);
        }
        "segangle" => {
            let (s, e, t) = (input.v3(), input.v3(), input.v3());
            let (angle, parameter) = segment_point_angle_deg_to_perp(s, e, t);
            out.d(angle);
            out.d(parameter);
        }
        "linecircle" => {
            let (p, d, c, n) = (input.v3(), input.v3(), input.v3(), input.v3());
            let (radius, tol) = (input.d(), input.d());
            if let Some(hits) = out.status(line_coplanar_circle_intersect(p, d, c, n, radius, tol))
            {
                out.i(hits.count as i64);
                for parameter in &hits.parameters[..hits.count] {
                    out.d(*parameter);
                }
            }
        }
        "twoplanes" => {
            let (p1, n1, p2, n2) = (input.v3(), input.v3(), input.v3(), input.v3());
            match intersect_two_planes(p1, n1, p2, n2) {
                Ok(PlanePair::Line { point, direction }) => {
                    out.status(Ok(()));
                    out.flag(false);
                    out.v3(point);
                    out.v3(direction);
                }
                Ok(PlanePair::Parallel { coincident }) => {
                    out.status::<()>(Err(Fault::Degenerate));
                    out.flag(coincident);
                }
                Err(fault) => {
                    out.status::<()>(Err(fault));
                    out.flag(false);
                }
            }
        }
        "threeplanes" => {
            let [p1, n1, p2, n2, p3, n3] = std::array::from_fn(|_| input.v3());
            if let Some(point) = out.status(intersect_three_planes(p1, n1, p2, n2, p3, n3)) {
                out.v3(point);
            }
        }
        "planedist" => {
            let (p, n, t) = (input.v3(), input.v3(), input.v3());
            if let Some(distance) = out.status(plane_point_distance(p, n, t)) {
                out.d(distance);
            }
        }
        "ptproj" => {
            let (p, q, n) = (input.v3(), input.v3(), input.v3());
            out.status(Ok(()));
            out.v3(point_project_to_plane(p, q, n));
        }
        "vecproj" => {
            let (v, n) = (input.v3(), input.v3());
            out.status(Ok(()));
            out.v3(vector_project_to_plane(v, n));
        }
        "lincomb" => {
            let (a, b, v) = (input.v3(), input.v3(), input.v3());
            if let Some(uv) = out.status(vec_lin_comb_two_vectors(a, b, v)) {
                out.v2(uv);
            }
        }
        "det3" => {
            let (a, b, c) = (input.v3(), input.v3(), input.v3());
            out.d(determinant_3_vectors(a, b, c));
        }
        "frenet" => {
            let count = input.i() as usize;
            let values = std::array::from_fn(|_| input.v3());
            let mut frame = [Vector3::new(5.0, 5.0, 5.0); 4];
            out.status(frenet_from_derivatives(&values, count, &mut frame));
            for v in frame {
                out.v3(v);
            }
        }
        "continuity" => {
            let [p1, d1, s1, p2, d2, s2] = std::array::from_fn(|_| input.v3());
            let degrees = input.d();
            let has = input.i() != 0;
            let (t1, t2) = (input.v3(), input.v3());
            let a = CurveJet {
                position: p1,
                first: d1,
                second: s1,
                third: has.then_some(t1),
            };
            let b = CurveJet {
                position: p2,
                first: d2,
                second: s2,
                third: has.then_some(t2),
            };
            out.status(Ok(()));
            out.i(evaluate_curve_continuity(&a, &b, degrees) as i64);
        }
        other => panic!("unknown operation {other}"),
    }
}

#[test]
fn foundation_matches_smlib() {
    // Inverse sine and cosine come from libm in the port and from the C
    // library in the oracle, so their results may differ in the last places.
    let library = |op: &str| matches!(op, "angle" | "segangle" | "ext3" | "parallel" | "perp");
    let report = run_gate(
        include_str!("fixtures/foundation/cases.txt"),
        |op| if library(op) { 16 } else { 4 },
        evaluate,
    );
    report.assert_passed(&["tangent", "touching", "singular"]);
}

/// The one gate value that is not bit-identical (`angle`, fixture line 5)
/// comes from the arc cosine alone. Port and oracle form the same cosine
/// bit for bit. Its exact arc cosine, 1.168113111852943639..., lies 0.026
/// ulp from the midpoint between two binary64 values: the C library returns
/// the nearer one, `libm` the other. The port keeps `libm` because the same
/// candidate on every machine matters more than matching one C library.
#[test]
fn the_one_inexact_value_is_an_arc_cosine_near_tie() {
    let a = Vector3::new(
        f64::from_bits(0xc01ea3969132d455),
        f64::from_bits(0x4002025ff5719220),
        f64::from_bits(0x3ff1b3015fd2ffc0),
    );
    let b = Vector3::new(
        f64::from_bits(0xc020f9c1bcd4b3e9),
        f64::from_bits(0xc01810b916a9baa3),
        f64::from_bits(0xc020c35d46bdba00),
    );
    let dot = a.dot(b);
    let cosine = (dot * dot) / a.length_squared() / b.length_squared();
    let cosine = cosine.sqrt();
    assert_eq!(cosine.to_bits(), 0x3fd914b2d74833c6);
    let port = a.angle_between(b).unwrap();
    assert_eq!(port, catallax_trace::scalar::acos(cosine));
    let smlib = f64::from_bits(0x3ff2b0975fd6feb3);
    assert_eq!(port.to_bits() - smlib.to_bits(), 1);
}

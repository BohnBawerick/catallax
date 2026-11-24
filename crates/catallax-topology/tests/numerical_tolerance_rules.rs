//! Tolerance rules T0 to T8 (decision Q7, section 4.3 of the Checkpoint 2
//! decisions) on admitted bodies: the cap at its boundary, rule T1's size
//! warning, rule T3's refit before an opt-in, and rule T7 along a chain.
#[path = "common/analytic.rs"]
mod fixture;
#[path = "common/tent.rs"]
mod tent;
use catallax_numbers::{Decider, NumericalCheck, Rational};
use catallax_shapes::Surface;
use catallax_topology::numerical::{
    self as n, NumericalOptions, ToleranceWarning,
    tolerance::{self, EntityKind},
};
use fixture::{bowed_copy, q, quarter, scaled};

/// The quarter cylinder's curved face.
const WALL: usize = 4;

fn options(cap: &str) -> NumericalOptions {
    NumericalOptions {
        tolerance_mm: q(cap),
        vertex_tolerance_mm: None,
    }
}

fn with_radius(value: &str) -> n::NumericalDefinition {
    let mut source = quarter();
    set_radius(&mut source, value);
    source
}

fn set_radius(source: &mut n::NumericalDefinition, value: &str) {
    let Surface::Cylinder(cylinder) = &mut source.faces[WALL].surface else {
        unreachable!()
    };
    cylinder.radius_sq = (&q(value) * &q(value)).into();
}

fn wall(value: &str) -> Surface {
    with_radius(value).faces[WALL].surface.clone()
}

/// The largest measured edge-use record.
fn largest_record(body: &n::NumericalBody) -> Rational {
    body.evidence()
        .checks()
        .iter()
        .filter_map(|check| match check {
            NumericalCheck::Samples { name, samples, .. }
                if name == "trim.edge-surface-distance" =>
            {
                Some(samples[0].distance_upper_mm.clone())
            }
            _ => None,
        })
        .fold(Rational::zero(), Rational::max)
}

/// The near-cap boundary (reviewer note relayed in inbox 003): a raw gap
/// under the cap is admitted only when rule T4's full allowance, 1.1 times
/// the gap plus its measurement error, fits under it. Off the known-answer
/// fixtures: the quarter cylinder's wall is widened, so each wall edge sits
/// a known distance off it, measured by the checker's closed forms.
#[test]
fn a_raw_gap_under_the_cap_is_admitted_only_when_its_full_allowance_fits() {
    let cap = q("0.01");
    // 0.009 mm: 1.1 x 0.009 = 0.0099 mm, under the cap with its error.
    let fits = n::build(
        with_radius("1.009"),
        options("0.01"),
        &mut Decider::default(),
    )
    .expect("the full allowance fits");
    let stored = largest_record(&fits);
    assert!(stored <= cap, "{stored}");
    assert!(stored >= q("0.0098999"), "{stored}");
    // 0.0095 mm: the raw gap is under the cap, but 1.1 x 0.0095 = 0.01045
    // mm is not. The existing named refusal, with no threshold changed.
    let raw = q("0.0095");
    assert!(raw < cap);
    let refused = n::build(
        with_radius("1.0095"),
        options("0.01"),
        &mut Decider::default(),
    )
    .unwrap_err();
    assert_eq!(refused.code, "numbers.numerical-tolerance-exceeded");
}

/// Rule T1 on an admitted body: the floor of a body under 10 m is 1e-7 mm,
/// and a body over 1 km is admitted with a warning, its floor raised to
/// 1e-11 of its coordinate size.
#[test]
fn a_body_over_a_kilometre_is_admitted_with_a_warning_and_a_raised_floor() {
    let small = n::build(quarter(), options("0.01"), &mut Decider::default()).unwrap();
    assert_eq!(small.tolerances().floor_mm, q("0.0000001"));
    assert!(small.tolerances().warnings.is_empty());
    for entity in &small.tolerances().entities {
        assert!(entity.tolerance_mm >= small.tolerances().floor_mm);
    }
    let k = q("2000000");
    let big = n::build(
        scaled(quarter(), &k),
        options("0.01"),
        &mut Decider::default(),
    )
    .unwrap();
    let report = big.tolerances();
    // The frame's box holds each circle's centre widened by its radius, so
    // `R` is at least the 1 km of the quarter's own box.
    assert!(report.coordinate_size_mm >= q("1000000"));
    assert_eq!(
        report.floor_mm,
        tolerance::floor_mm(&report.coordinate_size_mm)
    );
    assert!(report.floor_mm >= q("0.00001"));
    assert!(report.warnings.iter().any(|warning| matches!(
        warning,
        ToleranceWarning::LargeBody { extent_mm } if *extent_mm > q("1000000")
    )));
}

/// Rule T7 along a chain of edits that starts at a near-zero input: the wall
/// of the quarter cylinder sits on its edges (each stores rule T1's floor),
/// then moves out in four steps. Each stage reports the absolute tolerance
/// before and after, the cap, and what is left under it. The warning fires
/// exactly where a tolerance rose more than ten times above its input, and
/// stays a warning: nothing is refused while the cap holds. The refit of a
/// grown edge cannot lower a gap that is in the geometry, so the warning of
/// the first step stands, small as its absolute value is.
#[test]
fn a_chain_reports_absolute_tolerances_and_headroom_at_every_stage() {
    let cap = q("0.01");
    let mut dec = Decider::default();
    let mut body = n::build(quarter(), options("0.01"), &mut dec).unwrap();
    let largest_edge = |body: &n::NumericalBody| {
        body.tolerances()
            .entities
            .iter()
            .filter(|entity| entity.kind == EntityKind::Edge)
            .map(|entity| entity.tolerance_mm.clone())
            .fold(Rational::zero(), Rational::max)
    };
    let start = largest_edge(&body);
    println!(
        "stage 0: gap 0 mm, largest edge {start} mm (near zero, at the floor), cap {cap} mm, headroom {} mm",
        body.tolerances().headroom_mm().to_f64_nearest()
    );
    assert_eq!(start, q("0.0000001"));
    let steps = [
        ("1.000001", true),
        ("1.000003", false),
        ("1.0005", true),
        ("1.0009", false),
    ];
    for (stage, (radius, warns)) in steps.into_iter().enumerate() {
        let input = largest_edge(&body);
        body = body.replace_surface("f4", wall(radius), &mut dec).unwrap();
        let gap = q(radius) - q("1");
        let output = largest_edge(&body);
        let report = body.tolerances();
        let growth: Vec<_> = report
            .warnings
            .iter()
            .filter_map(|warning| match warning {
                ToleranceWarning::Growth {
                    kind: EntityKind::Edge,
                    input_mm,
                    output_mm,
                    cap_mm,
                    headroom_mm,
                    same_entity,
                    ..
                } => Some((
                    input_mm.clone(),
                    output_mm.clone(),
                    cap_mm.clone(),
                    headroom_mm.clone(),
                    *same_entity,
                )),
                _ => None,
            })
            .collect();
        println!(
            "stage {}: radius {radius}, gap {} mm, largest edge {} mm -> {} mm (x{:.2}), cap {} mm, headroom {} mm, {} growth warnings",
            stage + 1,
            gap.to_f64_nearest(),
            input.to_f64_nearest(),
            output.to_f64_nearest(),
            (&output / &input).to_f64_nearest(),
            report.cap_mm,
            report.headroom_mm().to_f64_nearest(),
            growth.len()
        );
        assert_eq!(report.cap_mm, cap);
        assert!(output >= &gap * &q("1.1"));
        assert!(output <= cap);
        assert_eq!(report.headroom_mm(), &cap - &report.largest_mm());
        assert_eq!(!growth.is_empty(), warns, "stage {}", stage + 1);
        for (input_mm, output_mm, cap_mm, headroom_mm, same_entity) in growth {
            assert!(same_entity);
            assert!(tolerance::grew(&input_mm, &output_mm));
            assert_eq!(headroom_mm, &cap_mm - &output_mm);
            assert!(headroom_mm > Rational::zero());
        }
    }
}

/// Rule T7's refit: a continuation gives one circle use of the wall a face
/// copy bowed off its 3D curve (the curve itself still lies on the wall),
/// which raises that edge's tolerance more than ten times. The face copy is
/// refitted from the 3D curve, measured again, and the lower tolerance is
/// kept, so no growth is left to warn about.
#[test]
fn a_grown_face_copy_is_refitted_and_its_tolerance_lowered() {
    let mut dec = Decider::default();
    let body = n::build(quarter(), options("0.01"), &mut dec).unwrap();
    let (bowed, edge) = bowed_copy(body.definition(), "0.006");
    let before = body
        .tolerances()
        .of(EntityKind::Edge, &edge)
        .unwrap()
        .clone();
    // Without the refit the bowed copy stores about 1.1 x 0.0045 mm.
    let unrefitted = n::build_with_native_trims(
        body.definition().clone(),
        vec![bowed.clone()],
        options("0.01"),
        &mut dec,
    )
    .unwrap();
    let bowed_mm = unrefitted
        .tolerances()
        .of(EntityKind::Edge, &edge)
        .unwrap()
        .clone();
    assert!(bowed_mm > q("0.004"), "{bowed_mm}");
    let next = body
        .rebuild_with_native_trims(
            body.definition().clone(),
            vec![bowed],
            options("0.01"),
            &mut dec,
        )
        .unwrap();
    let after = next
        .tolerances()
        .of(EntityKind::Edge, &edge)
        .unwrap()
        .clone();
    println!("edge {edge}: {before} mm before, {bowed_mm} mm bowed, {after} mm after the refit");
    assert!(after < q("0.000001"), "{after}");
    assert!(!next.tolerances().warnings.iter().any(|warning| matches!(
        warning,
        ToleranceWarning::Growth { key, .. } if *key == edge
    )));
}

/// Rule T3: an import whose face copy puts a stored tolerance over the
/// automatic 0.01 mm cap is refitted first; where the refit brings it under,
/// the body is admitted at the automatic cap. A gap in the geometry itself
/// survives the refit and refuses, naming the opt-in; an opt-in up to 0.1
/// mm admits it, and one over 0.1 mm refuses.
#[test]
fn an_import_is_refitted_before_the_user_is_asked_to_opt_in() {
    let mut dec = Decider::default();
    let source = quarter();
    let (bowed, _) = bowed_copy(&source, "0.03");
    assert_eq!(
        n::build_with_native_trims(
            source.clone(),
            vec![bowed.clone()],
            options("0.01"),
            &mut dec
        )
        .unwrap_err()
        .code,
        "numbers.numerical-tolerance-exceeded"
    );
    let admitted = n::admit_import(source, vec![bowed], options("0.01"), &mut dec).unwrap();
    assert!(admitted.tolerances().largest_mm() <= q("0.01"));
    assert_eq!(admitted.options().tolerance_mm, q("0.01"));

    let refused =
        n::admit_import(with_radius("1.02"), Vec::new(), options("0.01"), &mut dec).unwrap_err();
    assert_eq!(refused.code, "numbers.numerical-tolerance-exceeded");
    assert!(refused.detail.contains("opt in"), "{}", refused.detail);
    let opted =
        n::admit_import(with_radius("1.02"), Vec::new(), options("0.05"), &mut dec).unwrap();
    assert!(opted.tolerances().largest_mm() <= q("0.05"));
    assert_eq!(
        n::admit_import(with_radius("1.02"), Vec::new(), options("0.2"), &mut dec)
            .unwrap_err()
            .code,
        "numerical.import-opt-in"
    );
}

/// The quarter cylinder with every height scaled by `h`: the top face moves
/// down to `z = h` and the bottom face stays where it is.
fn squashed(mut d: n::NumericalDefinition, h: &Rational) -> n::NumericalDefinition {
    let real = |r: &Rational| -> catallax_numbers::Real { r.clone().into() };
    for vertex in &mut d.vertices {
        vertex.position[2] = &vertex.position[2] * h;
    }
    for edge in &mut d.edges {
        match &mut edge.curve {
            catallax_shapes::Curve::Line(line) => {
                line.origin.z = line.origin.z.clone() * real(h);
                line.dir.z = line.dir.z.clone() * real(h);
            }
            catallax_shapes::Curve::Circle(circle) => {
                circle.center.z = circle.center.z.clone() * real(h);
            }
            _ => unreachable!(),
        }
    }
    for (index, face) in d.faces.iter_mut().enumerate() {
        if let Surface::Plane(plane) = &mut face.surface {
            plane.origin.z = plane.origin.z.clone() * real(h);
        }
        // The axis of each face's parameters that runs along z: v on the
        // two side planes and the wall, u on the plane x = 0.
        let along_z = match index {
            0 | 1 => None,
            3 => Some(0),
            _ => Some(1),
        };
        for trim in face.loops.iter_mut().flatten() {
            for uv in &mut trim.uv {
                if let Some(axis) = along_z {
                    uv[axis] = &uv[axis] * h;
                }
            }
        }
    }
    d
}

/// Rule T5 after an edit (relayed reviewer note, inbox 004): an edit moves
/// the top face down towards the bottom face, which it does not touch and
/// which the edit leaves unchanged. The two faces share no topology, and the
/// pair rule reads them through their edges: the top and bottom arcs of the
/// wall, each carrying the wall's 0.0044 mm. Their pair is checked at every
/// height: an error under half the distance, a warning under a fifth, and
/// nothing just beyond five times the tolerances' sum. No threshold changes.
#[test]
fn an_edit_that_brings_an_unchanged_face_near_is_checked_by_the_pair_rule() {
    let mut dec = Decider::default();
    let body = n::build(with_radius("1.004"), options("0.01"), &mut dec).unwrap();
    let edges_of = |face: usize| -> std::collections::BTreeSet<String> {
        let d = body.definition();
        d.faces[face].loops[0]
            .iter()
            .map(|trim| d.edges[trim.edge].key.clone())
            .collect()
    };
    let (bottom, top) = (edges_of(0), edges_of(1));
    assert!(bottom.is_disjoint(&top));
    let arc = |keys: &std::collections::BTreeSet<String>| {
        keys.iter()
            .find(|key| {
                let d = body.definition();
                d.edges.iter().any(|edge| {
                    &edge.key == *key && matches!(edge.curve, catallax_shapes::Curve::Circle(_))
                })
            })
            .cloned()
            .unwrap()
    };
    let (low_arc, high_arc) = (arc(&bottom), arc(&top));
    let report = body.tolerances();
    let sum = report.of(EntityKind::Edge, &low_arc).unwrap()
        + report.of(EntityKind::Edge, &high_arc).unwrap();
    println!(
        "arcs {low_arc} and {high_arc}: tolerances add up to {} mm",
        sum.to_f64_nearest()
    );

    let edit = |h: Rational, dec: &mut Decider| {
        body.rebuild(
            squashed(body.definition().clone(), &h),
            options("0.01"),
            dec,
        )
    };
    let names_the_arcs = |pair: &n::tolerance::FeaturePair| {
        let keys = [pair.a.1.as_str(), pair.b.1.as_str()];
        keys.contains(&low_arc.as_str()) && keys.contains(&high_arc.as_str())
    };
    // Between d/5 and d/2: the arcs' pair is found and warned about.
    let warned = edit(&sum * &q("3.5"), &mut dec).unwrap();
    assert!(warned.tolerances().warnings.iter().any(|warning| matches!(
        warning,
        ToleranceWarning::Pair(pair) if names_the_arcs(pair)
    )));
    // Under d/2: the same pair refuses the edit.
    let refused = edit(&sum * &q("1.5"), &mut dec).unwrap_err();
    assert_eq!(refused.code, "numerical.tolerance-pair");
    // The first pair over half its distance names an entity of each face:
    // the two arcs, or a corner of each.
    let named = |keys: &std::collections::BTreeSet<String>, vertices: [&str; 3]| {
        keys.iter()
            .any(|key| refused.detail.contains(&format!("{key:?}")))
            || vertices
                .iter()
                .any(|key| refused.detail.contains(&format!("{key:?}")))
    };
    assert!(
        named(&bottom, ["v0", "v1", "v2"]) && named(&top, ["v3", "v4", "v5"]),
        "{}",
        refused.detail
    );
    // Just beyond five times the sum: no pair warning at all.
    let clear = edit(&sum * &q("5.05"), &mut dec).unwrap();
    assert!(
        !clear
            .tolerances()
            .warnings
            .iter()
            .any(|warning| matches!(warning, ToleranceWarning::Pair(_))),
        "{:?}",
        clear.tolerances().warnings
    );
}

/// Rule T3 by entity class (reviewer hold on 83858ed4): an imported vertex
/// may store up to 0.02 mm while edges and faces keep 0.01 mm.
/// The quarter cylinder's corner `v0` is moved off the ends of its three
/// edges; the edges themselves are unchanged.
#[test]
fn an_imported_vertex_takes_its_own_cap_and_edges_keep_theirs() {
    let mut dec = Decider::default();
    let moved = |by: &str| {
        let mut d = quarter();
        d.vertices[0].position[2] = q(by);
        d
    };
    // 0.015 mm off: a vertex-end record over 0.01 mm and under 0.02 mm.
    let admitted = n::admit_import(moved("0.015"), Vec::new(), options("0.01"), &mut dec).unwrap();
    assert_eq!(admitted.options().vertex_tolerance_mm, Some(q("0.02")));
    assert_eq!(admitted.tolerances().vertex_cap_mm, Some(q("0.02")));
    let v0 = admitted.tolerances().of(EntityKind::Vertex, "v0").unwrap();
    assert!((v0 - &q("0.0165")).abs() < q("0.000000001"), "{v0}");
    for entity in &admitted.tolerances().entities {
        if entity.kind == EntityKind::Edge {
            assert!(entity.tolerance_mm <= q("0.01"), "{entity:?}");
        }
    }
    // Own geometry has one class cap: the same vertex refuses at 0.01 mm.
    let own = n::build(moved("0.015"), options("0.01"), &mut dec).unwrap_err();
    assert_eq!(own.code, "numbers.numerical-tolerance-exceeded");
    assert!(own.detail.contains("vertex-distance"), "{}", own.detail);
    // Past 0.02 mm the vertex refuses under its own cap.
    let far = n::admit_import(moved("0.025"), Vec::new(), options("0.01"), &mut dec).unwrap_err();
    assert_eq!(far.code, "numbers.numerical-tolerance-exceeded");
    assert!(
        far.detail.contains("vertex-distance") && far.detail.contains("1/50"),
        "{}",
        far.detail
    );
    // 0.019 mm fits under the cap before T4 and exceeds it as 0.0209 mm.
    let under_raw_cap =
        n::admit_import(moved("0.019"), Vec::new(), options("0.01"), &mut dec).unwrap_err();
    assert_eq!(under_raw_cap.code, "numbers.numerical-tolerance-exceeded");
    assert!(
        under_raw_cap.detail.contains("vertex-distance") && under_raw_cap.detail.contains("1/50"),
        "{}",
        under_raw_cap.detail
    );
    // An edge at 0.015 mm is under the vertex cap and still refuses: edges
    // and faces keep 0.01 mm.
    let edge =
        n::admit_import(with_radius("1.015"), Vec::new(), options("0.01"), &mut dec).unwrap_err();
    assert_eq!(edge.code, "numbers.numerical-tolerance-exceeded");
    assert!(
        edge.detail.contains("edge-surface-distance"),
        "{}",
        edge.detail
    );
}

/// Rule T5 on a face's interior (reviewer hold on 83858ed4): the apex of a
/// tent sits above the middle of its 100 mm base, which it does not touch,
/// while every edge and vertex of the base is more than 70 mm away. The apex
/// vertex is 0.002 mm off its edges' ends, so it carries 0.0022 mm; the base
/// carries none. The pair is found through the production admission: a
/// warning between `d/5` and `d/2`, a refusal under `d/2`, nothing just
/// beyond five times the tolerance.
#[test]
fn a_vertex_near_a_large_face_interior_is_a_pair() {
    let mut dec = Decider::default();
    let near = n::build(tent::tent("0.006", "0.002"), options("0.01"), &mut dec).unwrap();
    let apex = near
        .tolerances()
        .of(EntityKind::Vertex, tent::APEX)
        .unwrap()
        .clone();
    assert!(apex >= q("0.0022"), "{apex}");
    let face_pair = |body: &n::NumericalBody| {
        body.tolerances()
            .warnings
            .iter()
            .find_map(|warning| match warning {
                ToleranceWarning::Pair(pair)
                    if [&pair.a, &pair.b]
                        .contains(&&(EntityKind::Vertex, tent::APEX.to_string()))
                        && [&pair.a, &pair.b]
                            .contains(&&(EntityKind::Face, tent::BASE.to_string())) =>
                {
                    Some(pair.clone())
                }
                _ => None,
            })
    };
    // 0.008 mm above the base: 0.0022 > 0.008 / 5.
    let pair = face_pair(&near).expect("the apex and the base are a pair");
    assert_eq!(
        pair.tolerance_sum_mm.to_f64_nearest().to_bits(),
        apex.to_f64_nearest().to_bits()
    );
    assert!(
        (pair.distance_mm.to_f64_nearest() - 0.008).abs() < 1e-9,
        "{pair:?}"
    );
    // No edge or vertex of the base is near: the pair is the face's interior.
    assert!(
        near.tolerances()
            .warnings
            .iter()
            .all(|warning| match warning {
                ToleranceWarning::Pair(pair) =>
                    pair.a.0 == EntityKind::Face || pair.b.0 == EntityKind::Face,
                _ => true,
            })
    );
    // 0.004 mm above: 0.0022 > 0.004 / 2, an error that names the face.
    let refused = n::build(tent::tent("0.002", "0.002"), options("0.01"), &mut dec).unwrap_err();
    assert_eq!(refused.code, "numerical.tolerance-pair");
    assert!(
        refused.detail.contains("face \"base\""),
        "{}",
        refused.detail
    );
    // Just beyond five times the apex's tolerance: no pair.
    let beyond = n::build(tent::tent("0.0094", "0.002"), options("0.01"), &mut dec).unwrap();
    assert!(
        face_pair(&beyond).is_none(),
        "{:?}",
        beyond.tolerances().warnings
    );
}

#[test]
fn a_small_face_interior_feature_tightens_the_provisional_t2_fit() {
    let definition = tent::tent("0.00004", "0");
    let fit = tolerance::fitted_target_for_definition(
        &definition,
        &[],
        &q("0.0000001"),
        &mut Decider::default(),
    )
    .unwrap();
    assert!((fit.to_f64_nearest() - 0.0000004).abs() < 1e-12, "{fit}");
}

/// The pair search is charged as work, and an allowance it exhausts refuses
/// by the operation budget's name (reviewer hold on 83858ed4). The near tent
/// has one candidate pair to measure; the far one has none. The difference
/// in work is at least the 81 points of the face grid the search samples. An
/// allowance between the two admits the far tent and refuses the near one.
#[test]
fn the_pair_search_is_charged_and_its_exhaustion_refuses_by_name() {
    use catallax_numbers::with_operation_budget;
    let admit = |height: &str, limit: u64| {
        with_operation_budget(limit, || {
            n::build(
                tent::tent(height, "0.002"),
                options("0.01"),
                &mut Decider::default(),
            )
        })
    };
    let (near, near_work) = admit("0.006", 20_000_000);
    let (far, far_work) = admit("1", 20_000_000);
    near.unwrap();
    far.unwrap();
    println!(
        "pair search: near {} work, far {} work",
        near_work.used, far_work.used
    );
    assert!(
        near_work.used >= far_work.used + 81,
        "{near_work:?} {far_work:?}"
    );
    let between = far_work.used + (near_work.used - far_work.used) / 2;
    let (far, _) = admit("1", between);
    far.unwrap();
    let (refused, work) = admit("0.006", between);
    assert_eq!(refused.unwrap_err().code, "numbers.operation-budget");
    assert!(work.exhausted);
}

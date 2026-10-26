//! Public exporter contracts and independent-read-back fixture generation.
use catallax_iges::{IgesError, IgesOptions, Units, write_iges};
use catallax_numbers::{Decider, Rational};
use catallax_shapes::{Pose, RVector3, Rotation};
use catallax_topology::{BodyId, Model, block, cone, cylinder, empty, sphere, torus};

fn q(s: &str) -> Rational {
    Rational::parse(s).unwrap()
}
fn sector_recipe(kind: &str, angle: &str) -> String {
    let points: &[[&str; 2]] = match kind {
        "cylinder" => &[["3", "0"], ["3", "2"], ["1", "2"], ["1", "0"]],
        "cone" => &[["3", "0"], ["5", "2"], ["3", "2"], ["1", "0"]],
        "scaled-cone" => &[["1", "0"], ["5/3", "2/3"], ["1", "2/3"], ["1/3", "0"]],
        "reversed-cone" => &[["5", "0"], ["3", "2"], ["1", "2"], ["3", "0"]],
        "uncovered-cone" => &[["3", "0"], ["4", "2"], ["2", "2"], ["1", "0"]],
        "apex-cone" => &[["3", "0"], ["0", "3"], ["0", "0"]],
        _ => panic!("unknown sector"),
    };
    let segments = points
        .iter()
        .enumerate()
        .map(|(i, [x, y])| format!(r#"{{"kind":"line","id":"s{i}","to":["{x}","{y}"]}}"#))
        .collect::<Vec<_>>()
        .join(",");
    format!(
        r#"{{"format":"catallax-recipe/3","name":"public axial sector",
        "steps":[{{"id":"profile","op":"sketch",
          "plane":{{"rotation":["1","0","0","0"],"translation":["0","0","0"]}},
          "loops":[[{segments}]]}},
          {{"id":"body","op":"revolve","profile":"profile",
            "axis":{{"point":["0","0","0"],"direction":["0","1","0"]}},"angle":"{angle}"}}],
        "result":"body"}}"#
    )
}
fn sector(kind: &str, angle: &str) -> (Model, BodyId) {
    if let Some(kind) = kind.strip_prefix("seam-") {
        let (mut model, body) = sector(kind, angle);
        let body = rechart_and_reverse_edges(&mut model, body);
        assert!(model.check(body).is_empty(), "{:?}", model.check(body));
        return (model, body);
    }
    let recipe = catallax_recipe::Recipe::parse(&sector_recipe(kind, angle)).unwrap();
    let built = catallax_features::evaluate(&recipe, &mut Decider::default()).unwrap();
    assert!(built.model.check(built.result).is_empty());
    (built.model, built.result)
}

// Reparameterize the same exact supports by a quarter turn and reverse every
// topological edge. This changes no point set, trim, or material orientation.
fn rechart_and_reverse_edges(m: &mut Model, b: BodyId) -> BodyId {
    use catallax_numbers::{Fact, FactKind, FactOrigin};
    use catallax_shapes::{GeomKey, Surface};
    use catallax_topology::{Loop, Sense};
    let reverse = |s| {
        if s == Sense::Same {
            Sense::Reversed
        } else {
            Sense::Same
        }
    };
    let edges: std::collections::BTreeMap<_, _> = m
        .edges_of_body(b)
        .into_iter()
        .map(|id| {
            let mut e = m.edge(id).clone();
            if let catallax_shapes::Curve::Circle(mut c) = m.geom().curve(e.curve).clone() {
                c.ref_dir = c.normal.cross(&c.ref_dir);
                e.curve = m.geom_mut().add_curve(catallax_shapes::Curve::Circle(c));
                for v in [e.start, e.end].into_iter().flatten() {
                    let point = m.vertex(v).point;
                    m.geom_mut().record_fact(Fact {
                        kind: FactKind::Incident,
                        a: GeomKey::Point(point),
                        b: GeomKey::Curve(e.curve),
                        at: None,
                        origin: FactOrigin::Construction("test.same-circle"),
                    });
                }
            }
            std::mem::swap(&mut e.start, &mut e.end);
            e.sense = reverse(e.sense);
            (id, m.add_edge(e))
        })
        .collect();
    let mut body = m.body(b).clone();
    for lump in &mut body.lumps {
        let mut l = m.lump(*lump).clone();
        for shell in &mut l.shells {
            let mut s = m.shell(*shell).clone();
            for face in &mut s.faces {
                let mut f = m.face(*face).clone();
                let mut surface = m.geom().surface(f.surface).clone();
                let changed = match &mut surface {
                    Surface::Cylinder(c) => {
                        c.ref_dir = c.axis.cross(&c.ref_dir);
                        true
                    }
                    Surface::Cone(c) => {
                        c.ref_dir = c.axis.cross(&c.ref_dir);
                        true
                    }
                    _ => false,
                };
                if changed {
                    f.surface = m.geom_mut().add_surface(surface);
                }
                for lp in &mut f.loops {
                    let Loop::Coedges(cs) = m.loop_(*lp).clone() else {
                        panic!("finite sector")
                    };
                    let mut mapped = Vec::new();
                    for c in cs {
                        let mut c = m.coedge(c).clone();
                        c.edge = edges[&c.edge];
                        c.sense = reverse(c.sense);
                        {
                            let curve = m.edge(c.edge).curve;
                            m.geom_mut().record_fact(Fact {
                                kind: FactKind::Incident,
                                a: GeomKey::Curve(curve),
                                b: GeomKey::Surface(f.surface),
                                at: None,
                                origin: FactOrigin::Construction("test.same-axial-support"),
                            });
                        }
                        mapped.push(m.add_coedge(c));
                    }
                    *lp = m.add_loop(Loop::Coedges(mapped));
                }
                *face = m.add_face(f);
            }
            *shell = m.add_shell(s);
        }
        *lump = m.add_lump(l);
    }
    m.add_body(body)
}
fn fixture(name: &str) -> (Model, BodyId) {
    if let Some((kind, angle)) = name.split_once("-sector-") {
        return sector(kind, angle);
    }
    let mut m = Model::new();
    let b = match name {
        "box" => block(&mut m, [q("2"), q("3"), q("5")], &Pose::identity()).unwrap(),
        "cylinder" => cylinder(&mut m, q("2"), q("5"), &Pose::identity()).unwrap(),
        "cone" => cone(&mut m, q("3"), q("1"), q("2"), &Pose::identity()).unwrap(),
        "rounded-box" => block(&mut m, [q("1/3"), q("2/7"), q("5/11")], &Pose::identity()).unwrap(),
        "posed-cylinder" => {
            let pose = Pose::new(
                Rotation::from_quaternion(q("2"), q("1"), q("0"), q("0")).unwrap(),
                RVector3::new(q("1/3"), q("-7"), q("1000")).into(),
                &mut Decider::default(),
            )
            .unwrap();
            cylinder(&mut m, q("2/3"), q("5/7"), &pose).unwrap()
        }
        "P2-collar" => {
            let recipe =
                catallax_recipe::Recipe::parse(include_str!("../../../parts/P2-collar.json"))
                    .unwrap();
            let built = catallax_features::evaluate(&recipe, &mut Decider::default()).unwrap();
            return (built.model, built.result);
        }
        "bored-plate" => {
            let recipe=catallax_recipe::Recipe::parse(r#"{
              "format":"catallax-recipe/3","name":"bored plate",
              "steps":[{"id":"profile","op":"sketch",
              "plane":{"rotation":["1","0","0","0"],"translation":["0","0","0"]},
              "loops":[[{"kind":"line","id":"a","to":["10","0"]},
                        {"kind":"line","id":"b","to":["10","8"]},
                        {"kind":"line","id":"c","to":["0","8"]},
                        {"kind":"line","id":"d","to":["0","0"]}],
                       [{"kind":"circle","id":"bore","centre":["5","4"],"radius":"2"}]]},
                       {"id":"plate","op":"extrude","profile":"profile","distance":"2"}],"result":"plate"}"#).unwrap();
            let built = catallax_features::evaluate(&recipe, &mut Decider::default()).unwrap();
            return (built.model, built.result);
        }
        _ => panic!("unknown fixture"),
    };
    (m, b)
}
fn export(m: &Model, b: BodyId, o: &IgesOptions) -> (Vec<u8>, catallax_iges::IgesReport) {
    let mut out = vec![];
    let report = write_iges(m, &[b], o, &mut Decider::default(), &mut out).unwrap();
    (out, report)
}
#[test]
fn public_fixtures_are_deterministic_shared_breps() {
    for (name, faces) in [
        ("box", 6),
        ("bored-plate", 7),
        ("P2-collar", 6),
        ("cylinder", 3),
        ("cone", 3),
    ] {
        let (m, b) = fixture(name);
        let (a, r) = export(&m, b, &IgesOptions::default());
        assert_eq!(a, export(&m, b, &IgesOptions::default()).0);
        assert_eq!(r.faces, faces);
        assert_eq!(r.edges - r.added_edges, m.edges_of_body(b).len());
        assert!(r.max_deviation <= IgesOptions::default().tolerance);
        assert_eq!(r.scale, q("1"));
        let s = String::from_utf8(a).unwrap();
        assert!(s.lines().all(|l| l.len() == 80));
        assert!(s.contains("8HCatallax"));
        assert_eq!(
            s.lines()
                .filter(|l| &l[72..73] == "D" && l[..8].trim() == "186")
                .count(),
            2
        );
    }
}
#[test]
fn partial_axial_sectors_export_through_the_recipe_path() {
    for kind in [
        "cylinder",
        "cone",
        "reversed-cone",
        "seam-cylinder",
        "seam-cone",
    ] {
        for angle in ["90", "180", "270"] {
            let (m, b) = sector(kind, angle);
            let (bytes, report) = export(&m, b, &IgesOptions::default());
            assert_eq!(report.faces, 6);
            assert_eq!(report.edges, m.edges_of_body(b).len());
            assert_eq!(report.vertices, m.vertices_of_body(b).len());
            assert_eq!(report.added_edges, 0);
            assert_eq!(report.added_vertices, 0);
            assert_eq!((bytes, report), export(&m, b, &IgesOptions::default()));
        }
    }
}

fn entities(bytes: &[u8]) -> std::collections::BTreeMap<usize, Vec<String>> {
    let mut rows = std::collections::BTreeMap::<usize, String>::new();
    for line in std::str::from_utf8(bytes).unwrap().lines() {
        if &line[72..73] == "P" {
            rows.entry(line[64..72].trim().parse().unwrap())
                .or_default()
                .push_str(line[..64].trim_end());
        }
    }
    rows.into_iter()
        .map(|(id, row)| {
            (
                id,
                row.trim_end_matches(';')
                    .split(',')
                    .map(str::to_owned)
                    .collect(),
            )
        })
        .collect()
}

#[test]
fn sector_pcurves_close_in_one_finite_unwrapped_rectangle() {
    for kind in ["cylinder", "cone", "seam-cylinder", "seam-cone"] {
        for angle in ["90", "180", "270"] {
            let (m, b) = sector(kind, angle);
            let (bytes, _) = export(&m, b, &IgesOptions::default());
            let graph = entities(&bytes);
            let mut rectangles = 0;
            let mut beyond_seam = false;
            for wire in graph.values().filter(|p| p[0] == "508" && p[1] == "4") {
                if wire[6] == "0" {
                    continue;
                } // Planar wires have no pcurves.
                let pcs: Vec<_> = wire[2..]
                    .chunks(7)
                    .map(|entry| {
                        assert_eq!(entry[4], "1");
                        assert_eq!(entry[5], "1");
                        &graph[&entry[6].parse::<usize>().unwrap()]
                    })
                    .collect();
                assert_eq!(pcs.len(), 4);
                let mut arcs = 0;
                for i in 0..4 {
                    let pc = pcs[i];
                    assert_eq!(pc[0], "126");
                    assert_eq!(pc[1], "1");
                    assert_eq!(pc[2], "1");
                    assert_eq!(q(&pc[19]), q("0"));
                    assert!(q(&pc[20]) > q("0"));
                    assert_eq!(q(&pc[16]), q(&pcs[(i + 1) % 4][13]));
                    assert_eq!(q(&pc[17]), q(&pcs[(i + 1) % 4][14]));
                    let du = &q(&pc[16]) - &q(&pc[13]);
                    let dv = &q(&pc[17]) - &q(&pc[14]);
                    if du.is_zero() {
                        assert!(!dv.is_zero());
                        assert_eq!(q(&pc[20]), q("1"));
                    } else {
                        assert_eq!(du.abs(), q(angle));
                        assert!(dv.is_zero());
                        arcs += 1;
                    }
                    beyond_seam |= [q(&pc[13]), q(&pc[16])]
                        .iter()
                        .any(|u| u < &q("0") || u > &q("360"));
                }
                assert_eq!(arcs, 2);
                rectangles += 1;
            }
            assert_eq!(rectangles, 2);
            if kind.starts_with("seam-") && angle == "270" {
                assert!(beyond_seam, "{kind} {angle}");
            }
        }
    }
}

#[test]
fn sector_angular_slant_and_scale_errors_are_added_and_enforced() {
    use catallax_numbers::Ball;
    let (m, b) = sector("scaled-cone", "270");
    for digits in [4, 17] {
        let mut o = IgesOptions {
            scale_digits: digits,
            tolerance: q("1"),
            ..Default::default()
        };
        let (bytes, r) = export(&m, b, &o);
        assert_ne!(r.scale, q("1"));
        let graph = entities(&bytes);
        let pc = graph
            .values()
            .find(|p| p[0] == "126" && p[3] == "1" && p[13] != p[16])
            .unwrap();
        let t = q(&pc[20]);
        let radians = Ball::pi(160) * Ball::from_rational(&q("3/2"), 160);
        let angular =
            &(&t - &radians.lo()).abs().max((&t - &radians.hi()).abs()) * &r.coordinate_bound;
        assert!(r.trim_deviation > angular); // Includes the nonzero slant rounding too.
        assert_eq!(
            r.max_deviation,
            &(&r.scale - &q("1")).abs() * &r.coordinate_bound + &r.scale * &r.trim_deviation
        );
        o.tolerance = r.max_deviation.clone();
        assert_eq!(export(&m, b, &o).1, r);
        o.tolerance = &r.max_deviation / &q("2");
        let mut sink = b"retained artifact".to_vec();
        assert!(matches!(
            write_iges(&m, &[b], &o, &mut Decider::default(), &mut sink),
            Err(IgesError::ToleranceExceeded { .. })
        ));
        assert_eq!(sink, b"retained artifact");
    }
}

#[test]
fn unsupported_sector_definitions_refuse_atomically() {
    for (kind, angle) in [
        ("cylinder", "45"),
        ("uncovered-cone", "90"),
        ("apex-cone", "90"),
    ] {
        let (m, b) = sector(kind, angle);
        let mut sink = b"retained artifact".to_vec();
        assert!(matches!(
            write_iges(
                &m,
                &[b],
                &IgesOptions::default(),
                &mut Decider::default(),
                &mut sink
            ),
            Err(IgesError::UnrepresentableDefinition(_))
        ));
        assert_eq!(sink, b"retained artifact");
    }
}
#[test]
fn strict_bound_gate_is_inclusive_and_never_touches_sink_on_refusal() {
    let mut m = Model::new();
    let b = block(&mut m, [q("1/3"), q("2/7"), q("5/11")], &Pose::identity()).unwrap();
    let mut options = IgesOptions {
        scale_digits: 5,
        tolerance: q("1"),
        ..Default::default()
    };
    let (_, r) = export(&m, b, &options);
    assert!(r.max_deviation > q("0"));
    assert_eq!(
        r.max_deviation,
        (&r.scale - &q("1")).abs() * r.coordinate_bound.clone()
    );
    options.tolerance = r.max_deviation.clone();
    assert_eq!(export(&m, b, &options).1, r);
    options.tolerance = &r.max_deviation / &q("2");
    let mut out = b"existing artifact".to_vec();
    let before = out.clone();
    assert!(matches!(
        write_iges(&m, &[b], &options, &mut Decider::default(), &mut out),
        Err(IgesError::ToleranceExceeded { .. })
    ));
    assert_eq!(out, before);
}
#[test]
fn unsupported_kinds_and_options_are_typed_and_atomic() {
    let mut m = Model::new();
    let bodies = [
        sphere(&mut m, q("2"), &Pose::identity()).unwrap(),
        torus(&mut m, q("3"), q("1"), &Pose::identity()).unwrap(),
    ];
    for b in bodies {
        let mut out = vec![];
        assert!(matches!(
            write_iges(
                &m,
                &[b],
                &IgesOptions::default(),
                &mut Decider::default(),
                &mut out
            ),
            Err(IgesError::UnsupportedSurface(_))
        ));
        assert!(out.is_empty());
    }
    let (m, b) = fixture("box");
    let mut out = vec![];
    let o = IgesOptions {
        scale_digits: 0,
        ..Default::default()
    };
    assert!(matches!(
        write_iges(&m, &[b], &o, &mut Decider::default(), &mut out),
        Err(IgesError::InvalidOptions)
    ));
    assert!(out.is_empty());
}
#[test]
fn rational_poses_units_and_scale_are_certified() {
    let mut m = Model::new();
    let pose = Pose::new(
        Rotation::from_quaternion(q("2"), q("1"), q("0"), q("0")).unwrap(),
        RVector3::new(q("1/3"), q("-7"), q("1000")).into(),
        &mut Decider::default(),
    )
    .unwrap();
    for b in [
        block(&mut m, [q("2"), q("3"), q("5")], &pose).unwrap(),
        cylinder(&mut m, q("2/3"), q("5/7"), &pose).unwrap(),
    ] {
        {
            let units = Units::Millimetres;
            let o = IgesOptions {
                units,
                ..Default::default()
            };
            let (bytes, r) = export(&m, b, &o);
            assert!(r.max_deviation <= o.tolerance);
            assert_eq!(r.units, units);
            assert!(String::from_utf8(bytes).unwrap().contains(match units {
                Units::Millimetres => ",2,2HMM,",
            }));
        }
    }
}
#[test]
fn empty_and_multiple_solids() {
    let mut m = Model::new();
    let z = empty(&mut m);
    let (_, r) = export(&m, z, &IgesOptions::default());
    assert_eq!(r.solids, 0);
    let a = block(&mut m, [q("1"), q("2"), q("3")], &Pose::identity()).unwrap();
    let b = cylinder(&mut m, q("1"), q("2"), &Pose::identity()).unwrap();
    let mut out = vec![];
    let r = write_iges(
        &m,
        &[a, b],
        &IgesOptions::default(),
        &mut Decider::default(),
        &mut out,
    )
    .unwrap();
    assert_eq!(r.solids, 2);
}
#[test]
#[ignore = "writes independent OpenCASCADE read-back inputs; set CATALLAX_IGES_FIXTURES"]
fn write_readback_fixtures() {
    let out = std::path::PathBuf::from(std::env::var("CATALLAX_IGES_FIXTURES").unwrap());
    std::fs::create_dir_all(&out).unwrap();
    for name in [
        "box",
        "bored-plate",
        "P2-collar",
        "cylinder",
        "cone",
        "rounded-box",
        "posed-cylinder",
        "cylinder-sector-90",
        "cylinder-sector-180",
        "cylinder-sector-270",
        "cone-sector-90",
        "cone-sector-180",
        "cone-sector-270",
        "reversed-cone-sector-270",
        "seam-cylinder-sector-270",
        "seam-cone-sector-270",
        "scaled-cone-sector-270",
    ] {
        let (m, b) = fixture(name);
        {
            let (suffix, units) = ("mm", Units::Millimetres);
            let (bytes, r) = export(
                &m,
                b,
                &IgesOptions {
                    units,
                    ..Default::default()
                },
            );
            std::fs::write(out.join(format!("{name}-{suffix}.igs")), bytes).unwrap();
            let c = m.counts(b);
            std::fs::write(out.join(format!("{name}-{suffix}.json")),format!(
                "{{\"max_deviation\":\"{}\",\"trim_deviation\":\"{}\",\"scale\":\"{}\",\"coordinate_bound\":\"{}\",\"faces\":{},\"edges\":{},\"vertices\":{},\"added_edges\":{},\"added_vertices\":{},\"source_edges\":{},\"source_vertices\":{}}}\n",
                r.max_deviation,r.trim_deviation,r.scale,r.coordinate_bound,r.faces,r.edges,r.vertices,r.added_edges,r.added_vertices,c.edges,c.vertices)).unwrap();
            println!("{name}-{suffix}: {r:?}");
        }
    }
}

#[test]
fn emitted_vertices_match_the_proved_whole_body_similarity() {
    for digits in [3, 7, 17] {
        let mut m = Model::new();
        let b = block(&mut m, [q("1/3"), q("2/7"), q("5/11")], &Pose::identity()).unwrap();
        let (bytes, r) = export(
            &m,
            b,
            &IgesOptions {
                scale_digits: digits,
                tolerance: q("1"),
                ..Default::default()
            },
        );
        let file = String::from_utf8(bytes).unwrap();
        let params: String = file
            .lines()
            .filter(|l| &l[72..73] == "P")
            .map(|l| l[..64].trim_end())
            .collect();
        let vertices = params.split(';').find(|s| s.starts_with("502,")).unwrap();
        let fields: Vec<_> = vertices.split(',').collect();
        assert_eq!(fields[1], "8");
        let emitted: std::collections::BTreeSet<Vec<_>> = fields[2..]
            .chunks(3)
            .map(|p| p.iter().map(|s| q(s)).collect())
            .collect();
        for vertex in m.vertices_of_body(b) {
            let p = m.geom().position(m.vertex(vertex).point);
            let exact = [p.x, p.y, p.z].map(|x| x.as_rational().unwrap());
            let mapped: Vec<_> = exact.iter().map(|x| x * &r.scale).collect();
            assert!(emitted.contains(&mapped));
            let distance_sq = (0..3)
                .map(|i| {
                    let d = &mapped[i] - &exact[i];
                    &d * &d
                })
                .fold(q("0"), |a, b| a + b);
            assert!(distance_sq <= &r.max_deviation * &r.max_deviation);
        }
    }
}

#[test]
fn cone_trim_error_is_composed_and_strictly_enforced() {
    let (m, b) = fixture("P2-collar");
    let (_, r) = export(&m, b, &IgesOptions::default());
    assert!(r.trim_deviation > q("0"));
    assert_eq!(r.scale, q("1"));
    assert_eq!(r.max_deviation, r.trim_deviation);
    let mut out = vec![];
    let o = IgesOptions {
        tolerance: &r.max_deviation / &q("2"),
        ..Default::default()
    };
    assert!(matches!(
        write_iges(&m, &[b], &o, &mut Decider::default(), &mut out),
        Err(IgesError::ToleranceExceeded { .. })
    ));
    assert!(out.is_empty());
}

#[test]
fn unavailable_numeric_encodings_and_work_refuse_without_output() {
    let mut m = Model::new();
    let b = cone(&mut m, q("3"), q("1"), q("3"), &Pose::identity()).unwrap();
    let mut out = vec![];
    assert!(matches!(
        write_iges(
            &m,
            &[b],
            &IgesOptions::default(),
            &mut Decider::default(),
            &mut out
        ),
        Err(IgesError::UnrepresentableDefinition(_))
    ));
    assert!(out.is_empty());
    let (m, b) = fixture("box");
    let mut dec = Decider::new(catallax_numbers::BudgetPolicy {
        max_work: 0,
        ..Default::default()
    });
    assert!(matches!(
        write_iges(&m, &[b], &IgesOptions::default(), &mut dec, &mut out),
        Err(IgesError::Decision(_))
    ));
    assert!(out.is_empty());
}

#[test]
fn retained_step_bytes_and_source_model_are_unchanged() {
    for name in ["box", "bored-plate", "P2-collar"] {
        let (m, b) = fixture(name);
        let mut before = vec![];
        let mut after = vec![];
        let options = catallax_step::StepOptions::default();
        let report =
            catallax_step::write_step(&m, &[b], &options, &mut Decider::default(), &mut before)
                .unwrap();
        let _ = export(&m, b, &IgesOptions::default());
        let check =
            catallax_step::write_step(&m, &[b], &options, &mut Decider::default(), &mut after)
                .unwrap();
        assert_eq!(before, after);
        assert_eq!(report.value, check.value);
        assert_eq!(report.marks, check.marks);
    }
}

#[test]
fn marked_invalid_and_duplicate_bodies_are_rejected_atomically() {
    use catallax_numbers::{BoundKind, Mark};
    let (mut m, b) = fixture("box");
    let mut body = m.body(b).clone();
    body.marks.push(Mark::new(
        "test.bound",
        BoundKind::Distance,
        q("1/100"),
        "model bound",
    ));
    let marked = m.add_body(body);
    let mut invalid = m.body(b).clone();
    invalid.lumps.clear();
    let invalid = m.add_body(invalid);
    let mut out = b"retained".to_vec();
    assert!(matches!(
        write_iges(
            &m,
            &[marked],
            &IgesOptions::default(),
            &mut Decider::default(),
            &mut out
        ),
        Err(IgesError::UncertifiedModel)
    ));
    assert!(matches!(
        write_iges(
            &m,
            &[invalid],
            &IgesOptions::default(),
            &mut Decider::default(),
            &mut out
        ),
        Err(IgesError::InvalidTopology(_))
    ));
    assert!(matches!(
        write_iges(
            &m,
            &[b, b],
            &IgesOptions::default(),
            &mut Decider::default(),
            &mut out
        ),
        Err(IgesError::InvalidTopology(_))
    ));
    assert_eq!(out, b"retained");
}

#[test]
fn algebraic_positions_are_not_silently_rationalised() {
    use catallax_numbers::Real;
    use catallax_shapes::Vector3;
    let pose = Pose::new(
        Rotation::identity(),
        Vector3::new(Real::from_i64(2).sqrt(), Real::zero(), Real::zero()),
        &mut Decider::default(),
    )
    .unwrap();
    let mut m = Model::new();
    let b = block(&mut m, [q("1"), q("2"), q("3")], &pose).unwrap();
    let mut out = vec![];
    assert!(matches!(
        write_iges(
            &m,
            &[b],
            &IgesOptions::default(),
            &mut Decider::default(),
            &mut out
        ),
        Err(IgesError::UnrepresentableDefinition(_))
    ));
    assert!(out.is_empty());
}

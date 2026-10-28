use std::collections::BTreeMap;
use std::io;
use std::time::{Duration, Instant};

use catallax_intersect::{Component, SurfacePair, intersect_surfaces};
use catallax_numbers::{Answer, BudgetPolicy, Decider, MarkLog, Rational};
use catallax_shapes::{Cylinder, Pose, RPoint3, RVector3, Rotation, Sphere, Surface};
use catallax_topology::{
    self as topology, Body, BodyId, BodyKind, Coedge, Edge, Face, Loop, Lump, Model, Sense, Shell,
};

use crate::{StepOptions, StepReport, write_step};

mod certified_spline;

fn q(s: &str) -> Rational {
    Rational::parse(s).unwrap()
}
fn v(x: i64, y: i64, z: i64) -> RVector3 {
    RVector3::new(x.into(), y.into(), z.into())
}
fn point(x: i64, y: i64, z: i64) -> RPoint3 {
    RPoint3::new(x.into(), y.into(), z.into())
}

fn export(model: &Model, bodies: &[BodyId], options: &StepOptions) -> (String, Answer<StepReport>) {
    let mut bytes = Vec::new();
    let answer = write_step(model, bodies, options, &mut Decider::default(), &mut bytes).unwrap();
    (String::from_utf8(bytes).unwrap(), answer)
}

fn entities(file: &str) -> Vec<(&str, Vec<&str>)> {
    file.lines()
        .filter(|l| l.starts_with('#'))
        .enumerate()
        .map(|(i, line)| {
            let (number, entity) = line.split_once('=').unwrap();
            assert_eq!(number, format!("#{}", i + 1));
            let (name, rest) = entity.split_once('(').unwrap();
            (name, fields(rest.strip_suffix(");").unwrap()))
        })
        .collect()
}

fn fields(s: &str) -> Vec<&str> {
    let mut result = Vec::new();
    let (mut depth, mut quoted, mut start) = (0, false, 0);
    for (i, c) in s.char_indices() {
        match c {
            '\'' => quoted = !quoted,
            '(' if !quoted => depth += 1,
            ')' if !quoted => depth -= 1,
            ',' if !quoted && depth == 0 => {
                result.push(&s[start..i]);
                start = i + 1;
            }
            _ => {}
        }
    }
    result.push(&s[start..]);
    result
}

fn references(s: &str) -> Vec<usize> {
    let mut ids = Vec::new();
    let mut chars = s.chars().peekable();
    let mut quoted = false;
    while let Some(c) = chars.next() {
        if c == '\'' {
            quoted = !quoted;
        }
        if c == '#' && !quoted {
            let mut value = 0;
            while let Some(c) = chars.peek().copied().filter(char::is_ascii_digit) {
                value = value * 10 + c.to_digit(10).unwrap() as usize;
                chars.next();
            }
            ids.push(value);
        }
    }
    ids
}

fn step_real(s: &str) -> Rational {
    let s = s.replace(".E", "E");
    q(s.strip_suffix('.').unwrap_or(&s))
}

fn cartesian(fields: &[&str]) -> [Rational; 3] {
    let values = fields[1]
        .strip_prefix('(')
        .unwrap()
        .strip_suffix(')')
        .unwrap();
    let values = values.split(',').map(step_real).collect::<Vec<_>>();
    std::array::from_fn(|i| values[i].clone())
}

fn assert_walk(file: &str) {
    let es = entities(file);
    let mut pending = references(es[15].1[1]);
    pending.remove(0); // the fixed placement #17
    pending.reverse();
    let mut visited = vec![false; es.len() + 1];
    let mut next = 21;
    while let Some(id) = pending.pop() {
        assert!((1..=es.len()).contains(&id));
        if visited[id] {
            continue;
        }
        visited[id] = true;
        assert_eq!(id, next, "depth-first reference numbering");
        next += 1;
        let children = es[id - 1]
            .1
            .iter()
            .flat_map(|f| references(f))
            .collect::<Vec<_>>();
        pending.extend(children.into_iter().rev());
    }
    assert_eq!(next, es.len() + 1);
}

#[test]
fn primitive_entity_counts_and_added_topology() {
    let mut model = Model::new();
    let pose = Pose::identity();
    let block = topology::block(&mut model, [2.into(), 3.into(), 5.into()], &pose).unwrap();
    let cylinder = topology::cylinder(&mut model, 2.into(), 5.into(), &pose).unwrap();
    let sphere = topology::sphere(&mut model, 5.into(), &pose).unwrap();
    let torus = topology::torus(&mut model, 3.into(), 1.into(), &pose).unwrap();
    let cone = topology::cone(&mut model, 3.into(), 0.into(), 4.into(), &pose).unwrap();
    for (body, counts, added) in [
        (block, (6, 12, 8), (0, 0)),
        (cylinder, (3, 3, 2), (2, 1)),
        (sphere, (1, 0, 1), (1, 0)),
        (torus, (1, 2, 1), (1, 2)),
        (cone, (2, 2, 2), (1, 1)),
    ] {
        let (file, answer) = export(&model, &[body], &StepOptions::default());
        let report = answer.value;
        assert_eq!((report.faces, report.edges, report.vertices), counts);
        assert_eq!((report.added_vertices, report.added_edges), added);
        assert_eq!(report.solids, 1);
        assert_eq!(report.approximated_curves, 0);
        assert!(answer.marks.is_empty());
        if body == cone || body == torus || body == cylinder {
            assert!(report.max_deviation > Rational::zero());
        } else {
            assert_eq!(report.max_deviation, Rational::zero());
        }
        assert!(report.max_deviation <= StepOptions::default().curve_tolerance);
        assert_walk(&file);
        if body == block {
            let es = entities(&file);
            assert_eq!(es.len(), 170);
            let mut actual = BTreeMap::new();
            for (name, _) in es {
                *actual.entry(name).or_insert(0) += 1;
            }
            for (name, count) in [
                ("ADVANCED_FACE", 6),
                ("FACE_BOUND", 6),
                ("EDGE_LOOP", 6),
                ("ORIENTED_EDGE", 24),
                ("EDGE_CURVE", 12),
                ("VERTEX_POINT", 8),
                ("CARTESIAN_POINT", 27),
                ("LINE", 12),
                ("VECTOR", 12),
                ("DIRECTION", 26),
                ("AXIS2_PLACEMENT_3D", 7),
                ("PLANE", 6),
            ] {
                assert_eq!(actual[name], count, "{name}");
            }
            assert!(file.contains("LENGTH_MEASURE(1.E-7)"));
        }
        if body == cone {
            let es = entities(&file);
            let (_, fields) = es.iter().find(|(n, _)| *n == "CONICAL_SURFACE").unwrap();
            let angle = step_real(fields[3]);
            assert!((angle - q("0.6435011087932843868")).abs() < q("1e-17"));
        }
    }
}

#[test]
fn strings_empty_files_and_output_errors() {
    let mut model = Model::new();
    let empty = topology::empty(&mut model);
    let options = StepOptions {
        product_name: "a'b\\c é 😀\n".into(),
        ..StepOptions::default()
    };
    let (file, answer) = export(&model, &[empty], &options);
    assert_eq!(answer.value, StepReport::default());
    assert_eq!(entities(&file).len(), 20);
    assert!(file.contains("'a''b\\\\c \\X2\\00E9\\X0\\ \\X2\\D83DDE00\\X0\\\\X2\\000A\\X0\\'"));
    assert!(file.contains(",(#17),#14)"));
    assert_eq!(export(&model, &[], &options).0, file);
    struct Broken;
    impl io::Write for Broken {
        fn write(&mut self, _: &[u8]) -> io::Result<usize> {
            Err(io::Error::other("disk full"))
        }
        fn flush(&mut self) -> io::Result<()> {
            Ok(())
        }
    }
    let error =
        write_step(&model, &[], &options, &mut Decider::default(), &mut Broken).unwrap_err();
    assert_eq!(error.code, "step.io-error");
    assert_eq!(error.detail, "disk full");
}

#[test]
fn rounding_deviation_and_scale() {
    let mut model = Model::new();
    let third = topology::block(
        &mut model,
        [q("1/3"), q("1/3"), q("1/3")],
        &Pose::identity(),
    )
    .unwrap();
    let (file, answer) = export(
        &model,
        &[third],
        &StepOptions {
            significant_digits: 5,
            ..StepOptions::default()
        },
    );
    assert!(file.contains("(0.33333,0.33333,0.33333)"));
    assert!(&answer.value.max_deviation * &answer.value.max_deviation >= q("3/90000000000"));
    let mark = answer.marks.iter().next().unwrap();
    assert_eq!(mark.code, "step.tolerance-not-reached");
    assert_eq!(mark.bound, answer.value.max_deviation);
    assert_eq!(answer.marks.len(), 1);
    assert!(
        export(&model, &[third], &StepOptions::default())
            .1
            .marks
            .is_empty()
    );
    let huge = topology::block(
        &mut model,
        [q("2e30"), q("3e30"), q("5e30")],
        &Pose::identity(),
    )
    .unwrap();
    assert_eq!(
        export(&model, &[huge], &StepOptions::default())
            .1
            .value
            .max_deviation,
        Rational::zero()
    );
    let huge_third = topology::block(
        &mut model,
        [q("1e30") / q("3"), q("1"), q("1")],
        &Pose::identity(),
    )
    .unwrap();
    let answer = export(&model, &[huge_third], &StepOptions::default()).1;
    assert!(answer.value.max_deviation >= q("3e12"));
    assert_eq!(answer.marks.len(), 1);
}

#[test]
fn posed_blocks_bound_written_vertices_and_faces() {
    let mut seed = 0x76a9_72e5_0314_12f8u64;
    let mut next = || {
        seed ^= seed << 13;
        seed ^= seed >> 7;
        seed ^= seed << 17;
        (seed % 101) as i64 - 50
    };
    for case in 0..1000 {
        let mut model = Model::new();
        let rotation =
            Rotation::from_quaternion(1.into(), next().into(), next().into(), next().into())
                .unwrap();
        let translation = RVector3::new((next() * 20000).into(), next().into(), q("1/3"));
        let pose = Pose::new(
            rotation,
            (translation).into(),
            &mut catallax_numbers::Decider::default(),
        )
        .unwrap();
        let scale = q(&format!("1e{}", case % 61 - 30));
        let size = std::array::from_fn(|_| {
            Rational::from_ratio(next().abs() + 1, next().abs() + 1).unwrap() * scale.clone()
        });
        let body = topology::block(&mut model, size.clone(), &pose).unwrap();
        let (file, answer) = export(&model, &[body], &StepOptions::default());
        let es = entities(&file);
        let vertices = es
            .iter()
            .filter(|(n, _)| *n == "VERTEX_POINT")
            .map(|(_, f)| {
                let id = references(f[1])[0];
                cartesian(&es[id - 1].1)
            })
            .collect::<Vec<_>>();
        for corner in 0..8 {
            let p: [Rational; 3] = std::array::from_fn(|i| {
                if corner & (1 << i) == 0 {
                    Rational::zero()
                } else {
                    size[i].clone()
                }
            });
            let p = pose
                .apply_point(&RPoint3::new(p[0].clone(), p[1].clone(), p[2].clone()).into())
                .as_rational()
                .expect("rational STEP fixture");
            let xyz = [p.x, p.y, p.z];
            let distance_sq = vertices
                .iter()
                .map(|v| {
                    (0..3).fold(Rational::zero(), |s, i| {
                        s + (&xyz[i] - &v[i]).pow(2).unwrap()
                    })
                })
                .min()
                .unwrap();
            assert!(distance_sq <= answer.value.max_deviation.pow(2).unwrap());
        }
        let written_faces = es.iter().filter(|(name, _)| *name == "ADVANCED_FACE");
        for (face, (_, written)) in model.faces_of_body(body).into_iter().zip(written_faces) {
            let plane = references(written[2])[0];
            let placement = references(es[plane - 1].1[1])[0];
            let location = references(es[placement - 1].1[1])[0];
            let axis = references(es[placement - 1].1[2])[0];
            let location = cartesian(&es[location - 1].1);
            let axis = cartesian(&es[axis - 1].1);
            let norm_sq = axis.iter().fold(Rational::zero(), |sum, a| sum + a * a);
            for loop_ in &model.face(face).loops {
                let Loop::Coedges(coedges) = model.loop_(*loop_) else {
                    panic!("block loop");
                };
                for coedge in coedges {
                    let vertex = model.edge(model.coedge(*coedge).edge).start.unwrap();
                    let p = model.geom().position(model.vertex(vertex).point);
                    let p = [p.x, p.y, p.z].map(|r| r.as_rational().unwrap());
                    let signed = (0..3).fold(Rational::zero(), |s, i| {
                        s + (&p[i] - &location[i]) * axis[i].clone()
                    });
                    assert!(
                        signed.pow(2).unwrap()
                            <= answer.value.max_deviation.pow(2).unwrap() * norm_sq.clone()
                    );
                }
            }
        }
    }
}

#[test]
fn voids_multiple_lumps_and_repeated_bodies() {
    let mut model = Model::new();
    let outer = topology::block(
        &mut model,
        [10.into(), 10.into(), 10.into()],
        &Pose::identity(),
    )
    .unwrap();
    let inner = topology::block(
        &mut model,
        [2.into(), 2.into(), 2.into()],
        &Pose::new(
            Rotation::identity(),
            (v(4, 4, 4)).into(),
            &mut catallax_numbers::Decider::default(),
        )
        .unwrap(),
    )
    .unwrap();
    let outer_shell = model.lump(model.body(outer).lumps[0]).shells[0];
    let inner_shell = model.lump(model.body(inner).lumps[0]).shells[0];
    let faces = model
        .shell(inner_shell)
        .faces
        .clone()
        .into_iter()
        .map(|id| {
            let mut face = model.face(id).clone();
            face.sense = Sense::Reversed;
            face.loops = face
                .loops
                .iter()
                .map(|loop_| {
                    let Loop::Coedges(coedges) = model.loop_(*loop_).clone() else {
                        panic!("block loop");
                    };
                    let reversed = coedges
                        .into_iter()
                        .rev()
                        .map(|coedge| {
                            let mut coedge = model.coedge(coedge).clone();
                            coedge.sense = if coedge.sense == Sense::Same {
                                Sense::Reversed
                            } else {
                                Sense::Same
                            };
                            model.add_coedge(coedge)
                        })
                        .collect();
                    model.add_loop(Loop::Coedges(reversed))
                })
                .collect();
            model.add_face(face)
        })
        .collect();
    let reversed = model.add_shell(Shell { faces });
    let lump = model.add_lump(Lump {
        shells: vec![outer_shell, reversed],
    });
    let body = model.add_body(Body {
        kind: BodyKind::Solid,
        lumps: vec![lump],
        marks: MarkLog::new(),
    });
    assert!(model.check(body).is_empty());
    let (file, answer) = export(&model, &[body], &StepOptions::default());
    assert_eq!(
        (
            answer.value.solids,
            answer.value.faces,
            answer.value.edges,
            answer.value.vertices
        ),
        (1, 12, 24, 16)
    );
    let es = entities(&file);
    let void = es
        .iter()
        .find(|(n, _)| *n == "ORIENTED_CLOSED_SHELL")
        .unwrap();
    assert_eq!(void.1[3], ".F.");
    let shell = references(void.1[2])[0];
    for face in references(es[shell - 1].1[1]) {
        assert_eq!(es[face - 1].1[3], ".T.");
        for bound in references(es[face - 1].1[1]) {
            assert_eq!(es[bound - 1].1[2], ".F.");
        }
    }
    assert_walk(&file);
    let (file, answer) = export(&model, &[outer, outer], &StepOptions::default());
    assert_eq!(
        (
            answer.value.solids,
            answer.value.faces,
            answer.value.edges,
            answer.value.vertices
        ),
        (2, 12, 24, 16)
    );
    assert_eq!(
        entities(&file).iter().filter(|(n, _)| *n == "LINE").count(),
        12
    );
    assert_walk(&file);
    let second = topology::block(
        &mut model,
        [2.into(), 2.into(), 2.into()],
        &Pose::new(
            Rotation::identity(),
            (v(20, 0, 0)).into(),
            &mut catallax_numbers::Decider::default(),
        )
        .unwrap(),
    )
    .unwrap();
    let lumps = vec![model.body(outer).lumps[0], model.body(second).lumps[0]];
    let multiple = model.add_body(Body {
        kind: BodyKind::Solid,
        lumps,
        marks: MarkLog::new(),
    });
    assert!(model.check(multiple).is_empty());
    let (file, answer) = export(&model, &[multiple], &StepOptions::default());
    assert_eq!(
        (
            answer.value.solids,
            answer.value.faces,
            answer.value.edges,
            answer.value.vertices
        ),
        (2, 12, 24, 16)
    );
    assert_walk(&file);
}

fn branch_body() -> (Model, BodyId) {
    branch_body_with_sense(Sense::Same)
}

fn branch_body_with_sense(source_sense: Sense) -> (Model, BodyId) {
    let mut model = Model::new();
    let cylinder = model.geom_mut().add_surface(Surface::Cylinder(
        Cylinder::new(
            (point(0, 0, 0)).into(),
            (v(0, 0, 1)).into(),
            1.into(),
            (v(1, 0, 0)).into(),
            &mut catallax_numbers::Decider::default(),
        )
        .unwrap(),
    ));
    let sphere = model.geom_mut().add_surface(Surface::Sphere(
        Sphere::new(
            (point(2, 0, 0)).into(),
            2.into(),
            (v(0, 0, 1)).into(),
            (v(1, 0, 0)).into(),
            &mut catallax_numbers::Decider::default(),
        )
        .unwrap(),
    ));
    let pair = intersect_surfaces(model.geom_mut(), cylinder, sphere, &mut Decider::default())
        .unwrap()
        .value;
    let SurfacePair::Meets(components) = pair else {
        panic!("expected intersection");
    };
    let Component::Curve { curve, .. } = components[0] else {
        panic!("expected branch");
    };
    let edge = model.add_edge(Edge {
        curve,
        sense: source_sense,
        start: None,
        end: None,
    });
    let mut faces = Vec::new();
    for (surface, sense) in [(cylinder, Sense::Same), (sphere, Sense::Reversed)] {
        let sense = if source_sense == Sense::Same {
            sense
        } else if sense == Sense::Same {
            Sense::Reversed
        } else {
            Sense::Same
        };
        let coedge = model.add_coedge(Coedge { edge, sense });
        let loop_ = model.add_loop(Loop::Coedges(vec![coedge]));
        faces.push(model.add_face(Face {
            surface,
            sense: Sense::Same,
            loops: vec![loop_],
        }));
    }
    let shell = model.add_shell(Shell { faces });
    let lump = model.add_lump(Lump {
        shells: vec![shell],
    });
    let body = model.add_body(Body {
        kind: BodyKind::Solid,
        lumps: vec![lump],
        marks: MarkLog::new(),
    });
    (model, body)
}

fn check_spline(file: &str, report: &StepReport) {
    let es = entities(file);
    for (_, edge) in es.iter().filter(|(name, _)| *name == "EDGE_CURVE") {
        let mut curve = references(edge[3])[0];
        if matches!(es[curve - 1].0, "SURFACE_CURVE" | "SEAM_CURVE") {
            curve = references(es[curve - 1].1[1])[0];
        }
        let (name, spline) = &es[curve - 1];
        if *name != "B_SPLINE_CURVE_WITH_KNOTS" {
            continue;
        }
        let degree: usize = spline[1].parse().unwrap();
        let controls = references(spline[2])
            .iter()
            .map(|id| cartesian(&es[*id - 1].1))
            .collect::<Vec<_>>();
        let ends = if edge[4] == ".T." {
            [edge[1], edge[2]]
        } else {
            [edge[2], edge[1]]
        };
        for (index, end) in [(0, ends[0]), (controls.len() - 1, ends[1])] {
            let vertex = references(end)[0];
            let point = references(es[vertex - 1].1[1])[0];
            assert_eq!(controls[index], cartesian(&es[point - 1].1));
        }
        let multiplicities = fields(spline[6].trim_matches(['(', ')']));
        assert_eq!(multiplicities[0].parse::<usize>().unwrap(), degree + 1);
        assert_eq!(
            multiplicities.last().unwrap().parse::<usize>().unwrap(),
            degree + 1
        );
        for span in controls.windows(degree + 1).step_by(degree) {
            for j in 0..=8 {
                let t = Rational::from_ratio(j, 8).unwrap();
                let mut ps = span.to_vec();
                for level in (1..=degree).rev() {
                    for k in 0..level {
                        ps[k] = std::array::from_fn(|i| {
                            &ps[k][i] * &(Rational::one() - t.clone()) + &ps[k + 1][i] * &t
                        });
                    }
                }
                let p = &ps[0];
                let cylinder_residual = (&p[0] * &p[0] + &p[1] * &p[1] - Rational::one()).abs();
                let sphere_residual =
                    ((&p[0] - &q("2")).pow(2).unwrap() + &p[1] * &p[1] + &p[2] * &p[2] - q("4"))
                        .abs();
                assert!(cylinder_residual <= &report.max_deviation * &q("2"));
                assert!(sphere_residual <= &report.max_deviation * &q("4"));
            }
        }
    }
    assert_walk(file);
}

#[test]
fn branch_spline_has_certified_bounds_at_default_tight_and_small_budgets() {
    let (model, body) = branch_body();
    for tolerance in ["1e-7", "1e-12"] {
        let options = StepOptions {
            curve_tolerance: q(tolerance),
            ..StepOptions::default()
        };
        let (file, answer) = export(&model, &[body], &options);
        assert_eq!(answer.value.approximated_curves, 1);
        assert!(
            answer.value.max_deviation <= options.curve_tolerance,
            "{:?}",
            answer.value.max_deviation
        );
        assert!(answer.marks.is_empty());
        check_spline(&file, &answer.value);
    }
    let options = StepOptions::default();
    let mut bytes = Vec::new();
    let result = write_step(
        &model,
        &[body],
        &options,
        &mut Decider::new(BudgetPolicy {
            max_precision_bits: 64,
            max_work: 10_000,
        }),
        &mut bytes,
    );
    match result {
        Ok(answer) => {
            if answer.value.max_deviation > options.curve_tolerance {
                assert_eq!(answer.marks.len(), 1);
                assert_eq!(
                    answer.marks.iter().next().unwrap().bound,
                    answer.value.max_deviation
                );
            }
            check_spline(&String::from_utf8(bytes).unwrap(), &answer.value);
        }
        Err(_) => {
            assert!(bytes.is_empty());
        }
    }
}

#[test]
fn deterministic_threads_and_expired_deadline_write_nothing() {
    let mut model = Model::new();
    let body = topology::sphere(&mut model, 5.into(), &Pose::identity()).unwrap();
    let expected = export(&model, &[body], &StepOptions::default()).0;
    std::thread::scope(|scope| {
        for _ in 0..8 {
            let model = model.clone();
            let expected = &expected;
            scope.spawn(move || {
                assert_eq!(
                    &export(&model, &[body], &StepOptions::default()).0,
                    expected
                )
            });
        }
    });
    let mut dec = Decider::default();
    dec.set_deadline(Some(Instant::now() - Duration::from_secs(1)));
    let mut bytes = vec![0x7f];
    let error = write_step(
        &model,
        &[body],
        &StepOptions::default(),
        &mut dec,
        &mut bytes,
    )
    .unwrap_err();
    assert_eq!(error.code, "numbers.deadline");
    assert_eq!(bytes, vec![0x7f]);
    let empty = topology::empty(&mut model);
    let mut bytes = Vec::new();
    let answer = write_step(
        &model,
        &[empty],
        &StepOptions::default(),
        &mut dec,
        &mut bytes,
    )
    .unwrap();
    assert_eq!(answer.value, StepReport::default());
    assert_eq!(
        String::from_utf8(bytes).unwrap(),
        export(&model, &[empty], &StepOptions::default()).0
    );
}

#[test]
fn terminating_decimal_radii_and_rational_frames_have_zero_deviation() {
    let rotation = Rotation::from_quaternion(2.into(), 1.into(), 0.into(), 0.into()).unwrap();
    let pose = Pose::new(
        rotation,
        (RVector3::new(q("0.1"), q("-0.2"), q("2.5"))).into(),
        &mut catallax_numbers::Decider::default(),
    )
    .unwrap();
    let mut model = Model::new();
    let bodies = [topology::sphere(&mut model, q("0.6"), &pose).unwrap()];
    for body in bodies {
        let (file, answer) = export(&model, &[body], &StepOptions::default());
        assert_eq!(answer.value.max_deviation, Rational::zero());
        assert!(answer.marks.is_empty());
        assert_walk(&file);
    }
}

#[test]
fn irrational_intersection_circles_remain_analytic() {
    let mut model = Model::new();
    let cylinder = model.geom_mut().add_surface(Surface::Cylinder(
        Cylinder::new(
            (point(0, 0, 0)).into(),
            (v(0, 0, 1)).into(),
            5.into(),
            (v(1, 0, 0)).into(),
            &mut catallax_numbers::Decider::default(),
        )
        .unwrap(),
    ));
    let sphere = model.geom_mut().add_surface(Surface::Sphere(
        Sphere::new(
            (point(0, 0, 1)).into(),
            6.into(),
            (v(0, 0, 1)).into(),
            (v(1, 0, 0)).into(),
            &mut catallax_numbers::Decider::default(),
        )
        .unwrap(),
    ));
    let pair = intersect_surfaces(model.geom_mut(), cylinder, sphere, &mut Decider::default())
        .unwrap()
        .value;
    let SurfacePair::Meets(components) = pair else {
        panic!("expected circles");
    };
    let mut cylinder_loops = Vec::new();
    let mut sphere_faces = Vec::new();
    for component in components {
        let Component::Curve { curve, .. } = component else {
            panic!("expected circle");
        };
        let edge = model.add_edge(Edge {
            curve,
            sense: Sense::Same,
            start: None,
            end: None,
        });
        let coedge = model.add_coedge(Coedge {
            edge,
            sense: Sense::Same,
        });
        cylinder_loops.push(model.add_loop(Loop::Coedges(vec![coedge])));
        let coedge = model.add_coedge(Coedge {
            edge,
            sense: Sense::Reversed,
        });
        let loop_ = model.add_loop(Loop::Coedges(vec![coedge]));
        sphere_faces.push(model.add_face(Face {
            surface: sphere,
            sense: Sense::Same,
            loops: vec![loop_],
        }));
    }
    let cylinder_face = model.add_face(Face {
        surface: cylinder,
        sense: Sense::Same,
        loops: cylinder_loops,
    });
    let shell = model.add_shell(Shell {
        faces: std::iter::once(cylinder_face).chain(sphere_faces).collect(),
    });
    let lump = model.add_lump(Lump {
        shells: vec![shell],
    });
    let body = model.add_body(Body {
        kind: BodyKind::Solid,
        lumps: vec![lump],
        marks: MarkLog::new(),
    });
    assert!(model.check(body).is_empty());
    let (file, answer) = export(&model, &[body], &StepOptions::default());
    assert_eq!(answer.value.approximated_curves, 0);
    // Each cap's pole moves behind an export-only ring with one chart vertex.
    assert_eq!(answer.value.added_vertices, 4);
    assert!(answer.value.max_deviation > Rational::zero());
    assert!(answer.value.max_deviation <= StepOptions::default().curve_tolerance);
    let es = entities(&file);
    // Two trim circles, two pole rings and one meridian seam per cap chart.
    assert_eq!(es.iter().filter(|(name, _)| *name == "CIRCLE").count(), 6);
    // Every edge keeps an exact analytic 3D master: circles, and lines for
    // the cylinder's chart seam. Splines appear only as bounded pcurves.
    for (name, args) in &es {
        if *name != "EDGE_CURVE" {
            continue;
        }
        let mut id = references(args[3])[0];
        if matches!(es[id - 1].0, "SURFACE_CURVE" | "SEAM_CURVE") {
            id = references(es[id - 1].1[1])[0];
        }
        assert!(
            matches!(es[id - 1].0, "CIRCLE" | "LINE"),
            "{}",
            es[id - 1].0
        );
    }
    assert_walk(&file);
}

#[test]
fn low_digit_spline_knots_stay_distinct_and_the_bound_still_holds() {
    let (model, body) = branch_body();
    let (file, answer) = export(
        &model,
        &[body],
        &StepOptions {
            significant_digits: 1,
            ..StepOptions::default()
        },
    );
    let es = entities(&file);
    let (_, spline) = es
        .iter()
        .find(|(name, _)| *name == "B_SPLINE_CURVE_WITH_KNOTS")
        .unwrap();
    let knots = fields(spline[7].trim_matches(['(', ')']))
        .into_iter()
        .map(step_real)
        .collect::<Vec<_>>();
    assert!(knots.windows(2).all(|w| w[0] < w[1]));
    assert_eq!(answer.marks.len(), 1);
    check_spline(&file, &answer.value);
}

#[test]
fn reversed_branch_edge_fits_in_the_edges_direction() {
    let (forward, a) = branch_body();
    let (reverse, b) = branch_body_with_sense(Sense::Reversed);
    let (a, ar) = export(&forward, &[a], &StepOptions::default());
    let (b, br) = export(&reverse, &[b], &StepOptions::default());
    check_spline(&a, &ar.value);
    check_spline(&b, &br.value);
}

#[test]

fn irrational_intersection_lines_and_trimmed_seam_arcs_remain_analytic() {
    use catallax_intersect::{CurveSurface, intersect_curve_surface};
    use catallax_numbers::{Fact, FactKind, FactOrigin};
    use catallax_shapes::{Curve, GeomKey, Line, Plane};
    use catallax_topology::Vertex;
    let mut model = Model::new();
    let cylinder = model.geom_mut().add_surface(Surface::Cylinder(
        Cylinder::new(
            (point(0, 0, 0)).into(),
            (v(0, 0, 1)).into(),
            5.into(),
            (v(1, 0, 0)).into(),
            &mut catallax_numbers::Decider::default(),
        )
        .unwrap(),
    ));
    let cut = model.geom_mut().add_surface(Surface::Plane(
        Plane::new(
            (point(1, 0, 0)).into(),
            (v(1, 0, 0)).into(),
            (v(0, 1, 0)).into(),
            &mut catallax_numbers::Decider::default(),
        )
        .unwrap(),
    ));
    let bottom = model.geom_mut().add_surface(Surface::Plane(
        Plane::new(
            (point(0, 0, 0)).into(),
            (v(0, 0, -1)).into(),
            (v(1, 0, 0)).into(),
            &mut catallax_numbers::Decider::default(),
        )
        .unwrap(),
    ));
    let top = model.geom_mut().add_surface(Surface::Plane(
        Plane::new(
            (point(0, 0, 2)).into(),
            (v(0, 0, 1)).into(),
            (v(1, 0, 0)).into(),
            &mut catallax_numbers::Decider::default(),
        )
        .unwrap(),
    ));
    let mut intersect = |a, b| {
        let SurfacePair::Meets(components) =
            intersect_surfaces(model.geom_mut(), a, b, &mut Decider::default())
                .unwrap()
                .value
        else {
            panic!("expected curves");
        };
        components
            .into_iter()
            .map(|c| match c {
                Component::Curve { curve, .. } => curve,
                _ => panic!("expected curve"),
            })
            .collect::<Vec<_>>()
    };
    let lines = intersect(cut, cylinder);
    let bottom_circle = intersect(bottom, cylinder)[0];
    let top_circle = intersect(top, cylinder)[0];
    let mut vertices = Vec::new();
    for line in &lines {
        for cap in [bottom, top] {
            let CurveSurface::Meets(hits) =
                intersect_curve_surface(model.geom_mut(), *line, cap, &mut Decider::default())
                    .unwrap()
                    .value
            else {
                panic!("expected endpoint");
            };
            assert_eq!(hits.len(), 1);
            vertices.push(model.add_vertex(Vertex {
                point: hits[0].point,
            }));
        }
    }
    let [a, d, b, c] = vertices.try_into().unwrap();
    let bottom_chord = model.geom_mut().add_curve(Curve::Line(
        Line::new(
            (point(1, 0, 0)).into(),
            (v(0, 1, 0)).into(),
            &mut catallax_numbers::Decider::default(),
        )
        .unwrap(),
    ));
    let top_chord = model.geom_mut().add_curve(Curve::Line(
        Line::new(
            (point(1, 0, 2)).into(),
            (v(0, 1, 0)).into(),
            &mut catallax_numbers::Decider::default(),
        )
        .unwrap(),
    ));
    let mut edge = |curve, start, end, sense| {
        model.add_edge(Edge {
            curve,
            sense,
            start: Some(start),
            end: Some(end),
        })
    };
    let base_arc = edge(bottom_circle, a, b, Sense::Reversed);
    let top_arc = edge(top_circle, d, c, Sense::Same);
    let left = edge(lines[0], a, d, Sense::Reversed);
    let right = edge(lines[1], b, c, Sense::Same);
    let base_chord = edge(bottom_chord, a, b, Sense::Same);
    let top_chord = edge(top_chord, d, c, Sense::Same);
    let mut faces = Vec::new();
    for (surface, sense, boundary) in [
        (
            cylinder,
            Sense::Same,
            vec![
                (base_arc, Sense::Same),
                (right, Sense::Same),
                (top_arc, Sense::Reversed),
                (left, Sense::Reversed),
            ],
        ),
        (
            cut,
            Sense::Reversed,
            vec![
                (left, Sense::Same),
                (top_chord, Sense::Same),
                (right, Sense::Reversed),
                (base_chord, Sense::Reversed),
            ],
        ),
        (
            bottom,
            Sense::Same,
            vec![(base_arc, Sense::Reversed), (base_chord, Sense::Same)],
        ),
        (
            top,
            Sense::Same,
            vec![(top_arc, Sense::Same), (top_chord, Sense::Reversed)],
        ),
    ] {
        let mut coedges = Vec::new();
        for (edge, sense) in boundary {
            let source = model.edge(edge).clone();
            model.geom_mut().record_fact(Fact {
                kind: FactKind::Incident,
                a: GeomKey::Curve(source.curve),
                b: GeomKey::Surface(surface),
                at: None,
                origin: FactOrigin::Construction("step.test-clipped-cylinder"),
            });
            for vertex in [source.start.unwrap(), source.end.unwrap()] {
                let point = model.vertex(vertex).point;
                model.geom_mut().record_fact(Fact {
                    kind: FactKind::Incident,
                    a: GeomKey::Point(point),
                    b: GeomKey::Curve(source.curve),
                    at: None,
                    origin: FactOrigin::Construction("step.test-clipped-cylinder"),
                });
            }
            coedges.push(model.add_coedge(Coedge { edge, sense }));
        }
        let loop_ = model.add_loop(Loop::Coedges(coedges));
        faces.push(model.add_face(Face {
            surface,
            sense,
            loops: vec![loop_],
        }));
    }
    let shell = model.add_shell(Shell { faces });
    let lump = model.add_lump(Lump {
        shells: vec![shell],
    });
    let body = model.add_body(Body {
        kind: BodyKind::Solid,
        lumps: vec![lump],
        marks: MarkLog::new(),
    });
    assert!(model.check(body).is_empty());
    let (file, answer) = export(&model, &[body], &StepOptions::default());
    assert_eq!(
        (
            answer.value.faces,
            answer.value.edges,
            answer.value.vertices
        ),
        (4, 6, 4)
    );
    assert_eq!(answer.value.approximated_curves, 0);
    assert_eq!(answer.value.added_vertices, 0);
    assert!(answer.value.max_deviation > Rational::zero());
    assert!(answer.marks.is_empty());
    let es = entities(&file);
    for (_, edge) in es.iter().filter(|(name, _)| *name == "EDGE_CURVE") {
        let mut curve = references(edge[3])[0];
        if matches!(es[curve - 1].0, "SURFACE_CURVE" | "SEAM_CURVE") {
            curve = references(es[curve - 1].1[1])[0];
        }
        assert!(matches!(es[curve - 1].0, "LINE" | "CIRCLE"));
    }
    assert_walk(&file);
}

fn irrational_radius_cap() -> (Model, BodyId) {
    use catallax_numbers::{Fact, FactKind, FactOrigin};
    use catallax_shapes::{Circle, Curve, GeomKey, Plane};
    let mut model = Model::new();
    let sphere = model.geom_mut().add_surface(Surface::Sphere(
        Sphere::new(
            (point(0, 0, 0)).into(),
            5.into(),
            (v(0, 0, 1)).into(),
            (v(1, 0, 0)).into(),
            &mut catallax_numbers::Decider::default(),
        )
        .unwrap(),
    ));
    let plane = model.geom_mut().add_surface(Surface::Plane(
        Plane::new(
            (point(0, 0, 1)).into(),
            (v(0, 0, -1)).into(),
            (v(1, 0, 0)).into(),
            &mut catallax_numbers::Decider::default(),
        )
        .unwrap(),
    ));
    let curve = model.geom_mut().add_curve(Curve::Circle(
        Circle::new(
            (point(0, 0, 1)).into(),
            (v(0, 0, 1)).into(),
            24.into(),
            (v(1, 0, 0)).into(),
            &mut catallax_numbers::Decider::default(),
        )
        .unwrap(),
    ));
    let edge = model.add_edge(Edge {
        curve,
        sense: Sense::Same,
        start: None,
        end: None,
    });
    let mut faces = Vec::new();
    for (surface, sense) in [(sphere, Sense::Same), (plane, Sense::Reversed)] {
        model.geom_mut().record_fact(Fact {
            kind: FactKind::Incident,
            a: GeomKey::Curve(curve),
            b: GeomKey::Surface(surface),
            at: None,
            origin: FactOrigin::Construction("step.test-sphere-cap"),
        });
        let coedge = model.add_coedge(Coedge { edge, sense });
        let loop_ = model.add_loop(Loop::Coedges(vec![coedge]));
        faces.push(model.add_face(Face {
            surface,
            sense: Sense::Same,
            loops: vec![loop_],
        }));
    }
    let shell = model.add_shell(Shell { faces });
    let lump = model.add_lump(Lump {
        shells: vec![shell],
    });
    let body = model.add_body(Body {
        kind: BodyKind::Solid,
        lumps: vec![lump],
        marks: MarkLog::new(),
    });
    assert!(model.check(body).is_empty());
    (model, body)
}

#[test]
fn decimal_accuracy_is_certified_or_refused_without_writing() {
    let (model, body) = irrational_radius_cap();
    let options = StepOptions {
        significant_digits: 40,
        ..StepOptions::default()
    };
    let (file, answer) = export(&model, &[body], &options);
    assert!(answer.marks.is_empty());
    let es = entities(&file);
    // The trim circle has radius sqrt(24). The pole ring and the chart's
    // meridian circles have other radii.
    let (_, circle) = es
        .iter()
        .find(|(name, args)| {
            *name == "CIRCLE" && {
                let radius = step_real(args[2]);
                radius > q("4.8") && radius < q("4.95")
            }
        })
        .unwrap();
    let radius = step_real(circle[2]);
    let unit = q("1e-39");
    assert!((&radius - &unit).pow(2).unwrap() <= q("24"));
    assert!((&radius + &unit).pow(2).unwrap() >= q("24"));
    let mut bytes = vec![0x7f];
    let refusal = write_step(
        &model,
        &[body],
        &options,
        &mut Decider::new(BudgetPolicy {
            max_precision_bits: 64,
            max_work: 1_000_000_000,
        }),
        &mut bytes,
    )
    .unwrap_err();
    assert_eq!(refusal.code, "step.budget-exhausted");
    assert_eq!(bytes, vec![0x7f]);
}

#[test]
fn zero_work_answers_exact_rationals_or_refuses_before_writing() {
    let mut model = Model::new();
    let block = topology::block(
        &mut model,
        [2.into(), 3.into(), 5.into()],
        &Pose::identity(),
    )
    .unwrap();
    let cone = topology::cone(&mut model, 3.into(), 0.into(), 4.into(), &Pose::identity()).unwrap();
    let expected = export(&model, &[block], &StepOptions::default()).0;
    let mut bytes = Vec::new();
    let answer = write_step(
        &model,
        &[block],
        &StepOptions::default(),
        &mut Decider::new(BudgetPolicy {
            max_precision_bits: 0,
            max_work: 0,
        }),
        &mut bytes,
    )
    .unwrap();
    assert_eq!(answer.value.max_deviation, Rational::zero());
    assert_eq!(String::from_utf8(bytes).unwrap(), expected);
    let mut bytes = vec![0x7f];
    let refusal = write_step(
        &model,
        &[block, cone],
        &StepOptions::default(),
        &mut Decider::new(BudgetPolicy {
            max_precision_bits: 64,
            max_work: 0,
        }),
        &mut bytes,
    )
    .unwrap_err();
    assert_eq!(refusal.code, "step.budget-exhausted");
    assert_eq!(bytes, vec![0x7f]);
}

#[test]
fn a_trace_without_an_initial_fitting_budget_is_refused() {
    let (model, body) = branch_body();
    let edge = model.edges_of_body(body)[0];
    let curve = model.edge(edge).curve;
    let mut dec = Decider::default();
    let trace = catallax_intersect::trace(model.geom(), curve, None, None, &mut dec)
        .unwrap()
        .value;
    let origin = model.geom().curve(curve).parameter_origin(model.geom());
    let origin = crate::number::point(&origin, 17, &mut dec).unwrap();
    for max_precision_bits in [0, 64] {
        let refusal = crate::fit::fit(
            &trace,
            &origin,
            &origin,
            false,
            &StepOptions::default(),
            &Decider::new(BudgetPolicy {
                max_precision_bits,
                max_work: 0,
            }),
        )
        .err()
        .unwrap();
        assert_eq!(refusal.code, "step.budget-exhausted");
    }
}

#[test]
fn low_digit_cone_angles_stay_acute() {
    let mut model = Model::new();
    let body = topology::cone(
        &mut model,
        100.into(),
        0.into(),
        1.into(),
        &Pose::identity(),
    )
    .unwrap();
    let (file, answer) = export(
        &model,
        &[body],
        &StepOptions {
            significant_digits: 1,
            ..StepOptions::default()
        },
    );
    let es = entities(&file);
    let (_, cone) = es
        .iter()
        .find(|(name, _)| *name == "CONICAL_SURFACE")
        .unwrap();
    assert_eq!(step_real(cone[3]), Rational::one());
    assert_eq!(answer.marks.len(), 1);
    assert_eq!(
        answer.marks.iter().next().unwrap().bound,
        answer.value.max_deviation
    );
}

#[test]
fn steinmetz_ellipse_arcs_stay_analytic_with_a_valid_deviation_bound() {
    use catallax_numbers::Ball;
    use catallax_shapes::CurveId;
    use catallax_topology::Vertex;
    let mut model = Model::new();
    let a = model.geom_mut().add_surface(Surface::Cylinder(
        Cylinder::new(
            (point(0, 0, 0)).into(),
            (v(0, 0, 1)).into(),
            1.into(),
            (v(1, 0, 0)).into(),
            &mut catallax_numbers::Decider::default(),
        )
        .unwrap(),
    ));
    let b = model.geom_mut().add_surface(Surface::Cylinder(
        Cylinder::new(
            (point(0, 0, 0)).into(),
            (v(1, 0, 0)).into(),
            1.into(),
            (v(0, 0, 1)).into(),
            &mut catallax_numbers::Decider::default(),
        )
        .unwrap(),
    ));
    let SurfacePair::Meets(components) =
        intersect_surfaces(model.geom_mut(), a, b, &mut Decider::default())
            .unwrap()
            .value
    else {
        panic!("expected ellipses");
    };
    let curves = components
        .iter()
        .map(|component| match component {
            Component::Curve { curve, .. } => *curve,
            _ => panic!("expected curve"),
        })
        .collect::<Vec<CurveId>>();
    let Component::Curve { singular, .. } = &components[0] else {
        unreachable!()
    };
    let start = model.add_vertex(Vertex { point: singular[0] });
    let end = model.add_vertex(Vertex { point: singular[1] });
    let mut edges = Vec::new();
    for curve in curves {
        edges.push(model.add_edge(Edge {
            curve,
            sense: Sense::Same,
            start: Some(start),
            end: Some(end),
        }));
        edges.push(model.add_edge(Edge {
            curve,
            sense: Sense::Same,
            start: Some(end),
            end: Some(start),
        }));
    }
    let mut faces = Vec::new();
    for (surface, boundary) in [
        (a, [(edges[3], Sense::Same), (edges[1], Sense::Reversed)]),
        (a, [(edges[0], Sense::Same), (edges[2], Sense::Reversed)]),
        (b, [(edges[2], Sense::Same), (edges[1], Sense::Same)]),
        (
            b,
            [(edges[0], Sense::Reversed), (edges[3], Sense::Reversed)],
        ),
    ] {
        let coedges = boundary
            .into_iter()
            .map(|(edge, sense)| model.add_coedge(Coedge { edge, sense }))
            .collect();
        let loop_ = model.add_loop(Loop::Coedges(coedges));
        faces.push(model.add_face(Face {
            surface,
            sense: Sense::Same,
            loops: vec![loop_],
        }));
    }
    let shell = model.add_shell(Shell { faces });
    let lump = model.add_lump(Lump {
        shells: vec![shell],
    });
    let body = model.add_body(Body {
        kind: BodyKind::Solid,
        lumps: vec![lump],
        marks: MarkLog::new(),
    });
    assert!(model.check(body).is_empty());
    let (file, answer) = export(&model, &[body], &StepOptions::default());
    assert_eq!(
        (
            answer.value.faces,
            answer.value.edges,
            answer.value.vertices
        ),
        (4, 4, 2)
    );
    assert_eq!(answer.value.approximated_curves, 0);
    assert!(answer.marks.is_empty());
    let es = entities(&file);
    let ellipses = es
        .iter()
        .filter(|(name, _)| *name == "ELLIPSE")
        .collect::<Vec<_>>();
    assert_eq!(ellipses.len(), 2);
    let scalar = |q: &Rational| Ball::from_rational(q, 256);
    let length = |v: &[Ball; 3]| {
        (v[0].square() + v[1].square() + v[2].square())
            .sqrt()
            .unwrap()
    };
    for (_, ellipse) in ellipses {
        let placement = references(ellipse[1])[0];
        let origin = references(es[placement - 1].1[1])[0];
        let axis = references(es[placement - 1].1[2])[0];
        let reference = references(es[placement - 1].1[3])[0];
        let origin = cartesian(&es[origin - 1].1).each_ref().map(scalar);
        let axis = cartesian(&es[axis - 1].1).each_ref().map(scalar);
        let reference = cartesian(&es[reference - 1].1).each_ref().map(scalar);
        let norm = length(&axis);
        let z = axis.each_ref().map(|a| a.checked_div(&norm).unwrap());
        let dot = &reference[0] * &z[0] + &reference[1] * &z[1] + &reference[2] * &z[2];
        let x: [Ball; 3] = std::array::from_fn(|i| &reference[i] - &(&z[i] * &dot));
        let norm = length(&x);
        let x = x.each_ref().map(|a| a.checked_div(&norm).unwrap());
        let y: [Ball; 3] = std::array::from_fn(|i| {
            &z[(i + 1) % 3] * &x[(i + 2) % 3] - &z[(i + 2) % 3] * &x[(i + 1) % 3]
        });
        let major = scalar(&step_real(ellipse[2]));
        let minor = scalar(&step_real(ellipse[3]));
        for t in 0..=8 {
            let t = Rational::from_ratio(t, 4).unwrap();
            let d = Rational::one() + t.pow(2).unwrap();
            let cosine = scalar(&((Rational::one() - t.pow(2).unwrap()) / d.clone()));
            let sine = scalar(&((&t * &q("2")) / d));
            let radial = &major * &cosine;
            let transverse = &minor * &sine;
            let p: [Ball; 3] =
                std::array::from_fn(|i| &origin[i] + &(&x[i] * &radial) + &y[i] * &transverse);
            for pair in [[0, 1], [1, 2]] {
                let residual = p[pair[0]].square() + p[pair[1]].square() - scalar(&Rational::one());
                let magnitude = residual.lo().abs().max(residual.hi().abs());
                assert!(magnitude <= answer.value.max_deviation.clone() * q("2"));
            }
        }
    }
    assert_walk(&file);
}

#[test]
fn arithmetic_coordinates_can_cross_input_exponent_limits() {
    let huge = q("9e100000");
    for reversed in [false, true] {
        let mut model = Model::new();
        let rotation = if reversed {
            Rotation::from_quaternion(0.into(), 0.into(), 0.into(), 1.into()).unwrap()
        } else {
            Rotation::identity()
        };
        let translation = if reversed { -&huge } else { huge.clone() };
        let pose = Pose::new(
            rotation,
            (RVector3::new(translation, 0.into(), 0.into())).into(),
            &mut catallax_numbers::Decider::default(),
        )
        .unwrap();
        let body = topology::block(&mut model, [huge.clone(), 1.into(), 1.into()], &pose).unwrap();
        let mut bytes = Vec::new();
        let answer = write_step(
            &model,
            &[body],
            &StepOptions::default(),
            &mut Decider::new(BudgetPolicy {
                max_precision_bits: 64,
                max_work: 0,
            }),
            &mut bytes,
        )
        .unwrap();
        let file = String::from_utf8(bytes).unwrap();
        assert!(file.contains(if reversed {
            "-1.8E+100001"
        } else {
            "1.8E+100001"
        }));
        assert_eq!(answer.value.max_deviation, Rational::zero());
    }
    let mut model = Model::new();
    let body = topology::block(
        &mut model,
        [q("1e-100000") / q("3"), 1.into(), 1.into()],
        &Pose::identity(),
    )
    .unwrap();
    let mut bytes = Vec::new();
    let answer = write_step(
        &model,
        &[body],
        &StepOptions::default(),
        &mut Decider::new(BudgetPolicy {
            max_precision_bits: 64,
            max_work: 0,
        }),
        &mut bytes,
    )
    .unwrap();
    assert!(String::from_utf8(bytes).unwrap().contains("E-100001"));
    assert!(answer.value.max_deviation > Rational::zero());
    assert!(answer.value.max_deviation < q("1e-100000"));
}

fn algebraic_zero_with_biased_enclosures() -> catallax_numbers::Real {
    use catallax_numbers::{Outcome, Real, UPoly};
    let polynomial = UPoly::new(vec![q("-2"), q("0"), q("0"), q("1")]);
    let Outcome::Proven(mut roots) = Decider::default().real_roots(&polynomial) else {
        panic!("the cubic has one exact real root");
    };
    let alpha = Real::from_algebraic(roots.pop().unwrap().value);
    (alpha.pow_u32(3) - Real::from_i64(2)).square()
}

#[test]
fn stalled_algebraic_zero_decimal_uses_an_exact_proof() {
    let value = algebraic_zero_with_biased_enclosures();
    let mut dec = Decider::new(BudgetPolicy {
        max_precision_bits: 128,
        max_work: 1_000_000,
    });
    let number = crate::number::real(&value, 16, &mut dec).unwrap();
    assert_eq!(number.value, Rational::zero());
    assert_eq!(number.error, Rational::zero());
}

#[test]
fn stalled_tiny_positive_decimal_is_not_snapped_to_zero() {
    let value = algebraic_zero_with_biased_enclosures()
        + catallax_numbers::Real::from_rational(q("1e-200"));
    let mut dec = Decider::new(BudgetPolicy {
        max_precision_bits: 128,
        max_work: 1_000_000,
    });
    assert!(matches!(
        crate::number::real(&value, 16, &mut dec),
        Err(error) if error.code == "step.budget-exhausted"
    ));
}

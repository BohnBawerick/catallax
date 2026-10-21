use super::*;

#[derive(Clone)]
pub struct Chain {
    pub operands: Vec<Primitive>,
    pub steps: Vec<(usize, usize, BoolOp, Counts)>,
}
fn planes(f: usize, e: usize, v: usize) -> Counts {
    Counts {
        lumps: 1,
        shells: 1,
        faces: f,
        edges: e,
        rings: 0,
        vertices: v,
        surfaces: [f, 0, 0],
        curves: [e, 0, 0, 0],
    }
}
pub fn chain(number: usize) -> Chain {
    match number {
        41 => Chain {
            operands: vec![
                block_at([10; 3], [0; 3]),
                cyl_at(2, 12, [5, 5, -1], false),
                cyl_at(2, 12, [-1, 5, 5], true),
            ],
            steps: vec![
                (
                    0,
                    1,
                    BoolOp::Subtract,
                    Counts {
                        lumps: 1,
                        shells: 1,
                        faces: 7,
                        edges: 14,
                        rings: 2,
                        vertices: 8,
                        surfaces: [6, 1, 0],
                        curves: [12, 2, 0, 0],
                    },
                ),
                (
                    3,
                    2,
                    BoolOp::Subtract,
                    Counts {
                        lumps: 1,
                        shells: 1,
                        faces: 10,
                        edges: 20,
                        rings: 4,
                        vertices: 10,
                        surfaces: [6, 4, 0],
                        curves: [12, 4, 4, 0],
                    },
                ),
            ],
        },
        42 => Chain {
            operands: vec![
                block_at([3, 1, 1], [0; 3]),
                block_at([1, 3, 1], [2, 0, 0]),
                block_at([3, 1, 1], [0, 2, 0]),
                block_at([1, 3, 1], [0; 3]),
            ],
            steps: vec![
                (0, 1, BoolOp::Unite, planes(8, 18, 12)),
                (4, 2, BoolOp::Unite, planes(10, 24, 16)),
                (5, 3, BoolOp::Unite, planes(10, 24, 16)),
            ],
        },
        43 => Chain {
            operands: vec![
                ball_at(5, [0; 3]),
                ball_at(5, [6, 0, 0]),
                ball_at(5, [-6, 0, 0]),
            ],
            steps: vec![
                (
                    0,
                    1,
                    BoolOp::Subtract,
                    Counts {
                        lumps: 1,
                        shells: 1,
                        faces: 2,
                        edges: 1,
                        rings: 1,
                        vertices: 0,
                        surfaces: [0, 0, 2],
                        curves: [0, 1, 0, 0],
                    },
                ),
                (
                    3,
                    2,
                    BoolOp::Subtract,
                    Counts {
                        lumps: 1,
                        shells: 1,
                        faces: 3,
                        edges: 2,
                        rings: 2,
                        vertices: 0,
                        surfaces: [0, 0, 3],
                        curves: [0, 2, 0, 0],
                    },
                ),
            ],
        },
        44 => Chain {
            operands: (0..6).map(|k| block_at([2, 1, 1], [k, 0, 0])).collect(),
            steps: (1..6)
                .map(|k| {
                    (
                        if k == 1 { 0 } else { k + 4 },
                        k,
                        BoolOp::Unite,
                        planes(6, 12, 8),
                    )
                })
                .collect(),
        },
        _ => panic!("unknown chain"),
    }
}
impl Chain {
    pub fn transformed(&self, scale: &Rational, pose: &Pose) -> Self {
        Self {
            operands: self
                .operands
                .iter()
                .map(|s| s.transformed(scale, pose))
                .collect(),
            steps: self.steps.clone(),
        }
    }
    pub fn locate(&self, point: &RPoint3, step: usize) -> Option<Location> {
        let mut values: Vec<_> = self.operands.iter().map(|s| s.locate(point)).collect();
        if values.contains(&Location::On) {
            return None;
        }
        for &(a, b, op, _) in self.steps.iter().take(step + 1) {
            values.push(set_rule(op, values[a], values[b]).unwrap());
        }
        values.last().copied()
    }
}

fn chain_location(spec: &Chain, point: &RPoint3, step: usize) -> Location {
    let mut values: Vec<_> = spec.operands.iter().map(|p| p.locate(point)).collect();
    // Three-valued set logic certifies strict membership even on a discarded
    // operand boundary when another operand supplies a neighbourhood.
    for &(a, b, op, _) in spec.steps.iter().take(step + 1) {
        let (a, b) = (values[a], values[b]);
        let value = match op {
            BoolOp::Unite => {
                if a == Location::Inside || b == Location::Inside {
                    Location::Inside
                } else if a == Location::Outside && b == Location::Outside {
                    Location::Outside
                } else {
                    Location::On
                }
            }
            BoolOp::Intersect => {
                if a == Location::Outside || b == Location::Outside {
                    Location::Outside
                } else if a == Location::Inside && b == Location::Inside {
                    Location::Inside
                } else {
                    Location::On
                }
            }
            BoolOp::Subtract => {
                if a == Location::Outside || b == Location::Inside {
                    Location::Outside
                } else if a == Location::Inside && b == Location::Outside {
                    Location::Inside
                } else {
                    Location::On
                }
            }
            _ => unreachable!(),
        };
        values.push(value);
    }
    *values.last().unwrap()
}

fn chain_faces(m: &mut Model, body: BodyId, operands: &[BodyId], spec: &Chain, step: usize) {
    let surfaces: Vec<_> = operands
        .iter()
        .flat_map(|&b| m.faces_of_body(b))
        .map(|f| m.face(f).surface)
        .collect();
    let mut dec = Decider::default();
    for f in m.faces_of_body(body) {
        let surface = m.face(f).surface;
        assert!(surfaces.iter().any(|&s| matches!(
            m.geom().same_surface(surface, s, &mut dec),
            Outcome::ProvenCoincident(_)
        )));
        let point = exact(face_point(m, f, &mut dec));
        let (value, normal) = implicit_normal(m.geom().surface(surface), &point);
        assert_eq!(value, r(0));
        let normal = if m.face(f).sense == Sense::Same {
            normal
        } else {
            normal.scaled(&r(-1))
        };
        let t = spec
            .operands
            .windows(2)
            .map(|p| safe_step(&p[0], &p[1], &point, &normal))
            .min()
            .unwrap();
        assert_eq!(
            chain_location(spec, &point.plus(&normal.scaled(&t)), step),
            Location::Outside
        );
        assert_eq!(
            chain_location(spec, &point.plus(&normal.scaled(&(-t))), step),
            Location::Inside
        );
    }
}
pub fn run_chain(
    number: usize,
    scale: &Rational,
    pose: &Pose,
) -> (Counts, Vec<MarkLog>, Vec<DecisionStats>, Vec<String>) {
    let spec = chain(number).transformed(scale, pose);
    let mut m = Model::new();
    let mut bodies: Vec<_> = spec.operands.iter().map(|s| s.build(&mut m)).collect();
    let mut marks = Vec::new();
    let mut stats = Vec::new();
    let mut last = Counts::empty();
    for (step, (a, b, op, want)) in spec.steps.iter().enumerate() {
        let before = topology_snapshot(&m, &bodies);
        let mut dec = Decider::default();
        let answer = boolean(&mut m, bodies[*a], bodies[*b], *op, &mut dec).unwrap();
        assert_eq!(before, topology_snapshot(&m, &bodies));
        check_report(&answer, &DecisionStats::default(), dec.stats());
        assert!(answer.marks.is_empty());
        check_counts(&m, answer.body, want);
        topology_rules(&m, answer.body);
        stats.push(dec.stats().clone());
        marks.push(answer.marks.clone());
        last = counts(&m, answer.body);
        chain_faces(
            &mut m,
            answer.body,
            &bodies[..spec.operands.len()],
            &spec,
            step,
        );
        let mut rng = Rng(0xc4a105);
        let mut checked = 0;
        while checked < 32 {
            let point = moved_point(&rng.point(), scale, pose);
            if let Some(location) = spec.locate(&point, step) {
                query(&mut m, answer.body, &point, location, &mut dec);
                checked += 1;
            }
        }
        if step + 1 == spec.steps.len() {
            let explicit = match number {
                41 => vec![
                    (p("5", "5", "5"), Location::Outside),
                    (p("1", "1", "1"), Location::Inside),
                ],
                42 => vec![(p("3/2", "3/2", "1/2"), Location::Outside)],
                43 => vec![
                    (p("0", "0", "0"), Location::Inside),
                    (p("4", "0", "0"), Location::Outside),
                    (p("-4", "0", "0"), Location::Outside),
                ],
                44 => vec![
                    (p("13/2", "1/2", "1/2"), Location::Inside),
                    (p("7", "1/2", "1/2"), Location::On),
                    (p("15/2", "1/2", "1/2"), Location::Outside),
                ],
                _ => unreachable!(),
            };
            for (point, location) in explicit {
                query(
                    &mut m,
                    answer.body,
                    &moved_point(&point, scale, pose),
                    location,
                    &mut dec,
                );
            }
        }
        bodies.push(answer.body);
    }
    (last, marks, stats, topology_snapshot(&m, &bodies))
}

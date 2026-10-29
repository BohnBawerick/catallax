//! The final lattice stage: request building, refusals and rebuild from an edited body.
use catallax_field::{
    Location, Workspace,
    lattice_stage::{LatticeSpec, Settings, Stage, check_no_cad_after_lattice},
};
use catallax_numbers::{Decider, Rational};
use catallax_recipe::{Grade, LatticeCell, Number, Pose, Recipe};
use catallax_shapes::RPoint3;

fn q(value: &str) -> Rational {
    Rational::parse(value).unwrap()
}
fn body(id: &str, size: &str) -> Recipe {
    Recipe::parse(&format!(
        r#"{{"format":"catallax-recipe/3","name":"body","steps":[{{"id":"{id}","op":"block","size":["{size}","4","4"],"pose":{{"rotation":["1","0","0","0"],"translation":["0","0","0"]}}}}],"result":"{id}"}}"#
    ))
    .unwrap()
}
fn spec() -> LatticeSpec {
    LatticeSpec {
        cell: LatticeCell::Cubic,
        size: Number::parse("2").unwrap(),
        frame: serde_json::from_str::<Pose>(
            r#"{"rotation":["1","0","0","0"],"translation":["0","0","0"]}"#,
        )
        .unwrap(),
        radius: serde_json::from_str::<Grade>(r#"{"constant":"0.25"}"#).unwrap(),
    }
}
fn settings() -> Settings {
    Settings {
        accuracy: q("0.25"),
        minimum_detail: q("0.5"),
        max_cells: 2_000_000,
        max_vertices: 1_000_000,
    }
}
fn location(stage: &Stage, x: &str, y: &str, z: &str) -> Location {
    let point = RPoint3::new(q(x), q(y), q(z));
    stage
        .plan()
        .proven_checked(&point, &mut Workspace::default(), &mut Decider::default())
        .unwrap()
        .location
}

#[test]
fn the_stage_makes_the_lattice_the_result_of_the_body() {
    let stage = Stage::new(
        &body("block", "4"),
        spec(),
        settings(),
        &mut Decider::default(),
    )
    .unwrap();
    assert_eq!(stage.lattice_step().as_str(), "lattice");
    assert_eq!(stage.recipe().result, *stage.lattice_step());
    assert_eq!(location(&stage, "2", "1", "2"), Location::Inside);
    assert_eq!(location(&stage, "1", "1", "1"), Location::Outside);
}

fn first_again() -> Stage {
    Stage::new(
        &body("block", "4"),
        spec(),
        settings(),
        &mut Decider::default(),
    )
    .unwrap()
}

#[test]
fn an_edited_body_rebuilds_the_lattice_from_the_edited_body() {
    let first = Stage::new(
        &body("block", "4"),
        spec(),
        settings(),
        &mut Decider::default(),
    )
    .unwrap();
    let (earlier, kept_spec, kept_settings) = first.return_before_lattice();
    assert_eq!(kept_settings, settings());
    assert_eq!(kept_spec, spec());
    assert_eq!(earlier, body("block", "4"));
    assert_eq!(
        location(&first_again(), "3", "0.125", "2"),
        Location::Inside
    );
    // The block is shortened along x; the same lattice parameters now stop at x = 2.
    let edited = Stage::new(
        &body("block", "2"),
        kept_spec,
        kept_settings,
        &mut Decider::default(),
    )
    .unwrap();
    assert_eq!(location(&edited, "3", "0.125", "2"), Location::Outside);
    assert_eq!(location(&edited, "1", "0.125", "2"), Location::Inside);
}

#[test]
fn an_existing_lattice_in_the_body_is_cad_after_lattice() {
    let first = Stage::new(
        &body("block", "4"),
        spec(),
        settings(),
        &mut Decider::default(),
    )
    .unwrap();
    let error =
        Stage::new(first.recipe(), spec(), settings(), &mut Decider::default()).unwrap_err();
    assert_eq!(error.code, "lattice.cad-after-lattice");
}

#[test]
fn a_step_that_reads_a_lattice_is_cad_after_lattice() {
    let text = r#"{"format":"catallax-recipe/3","name":"x","steps":[{"id":"region","op":"block","size":["4","4","4"],"pose":{"rotation":["1","0","0","0"],"translation":["0","0","0"]}},{"id":"beams","op":"lattice","region":"region","cell":"cubic","size":"2","frame":{"rotation":["1","0","0","0"],"translation":["0","0","0"]},"radius":{"constant":"0.25"}},{"id":"after","op":"offset","body":"beams","distance":"0.1"}],"result":"after"}"#;
    let recipe = Recipe::parse(text).unwrap();
    assert_eq!(
        check_no_cad_after_lattice(&recipe).unwrap_err().code,
        "lattice.cad-after-lattice"
    );
    assert!(check_no_cad_after_lattice(&body("block", "4")).is_ok());
}

#[test]
fn impossible_settings_and_lattices_refuse_with_their_own_codes() {
    let mut dec = Decider::default();
    for bad in [
        Settings {
            accuracy: q("0"),
            ..settings()
        },
        Settings {
            minimum_detail: q("-1"),
            ..settings()
        },
        Settings {
            max_cells: 0,
            ..settings()
        },
        Settings {
            max_vertices: 3,
            ..settings()
        },
    ] {
        let error = Stage::new(&body("block", "4"), spec(), bad, &mut dec).unwrap_err();
        assert_eq!(error.code, "lattice.settings");
    }
    let thin = LatticeSpec {
        radius: serde_json::from_str::<Grade>(r#"{"constant":"0"}"#).unwrap(),
        ..spec()
    };
    let error = Stage::new(&body("block", "4"), thin, settings(), &mut dec).unwrap_err();
    assert_eq!(error.code, "lattice.impossible");
}

#[test]
fn the_lattice_step_id_avoids_a_taken_name() {
    let stage = Stage::new(
        &body("lattice", "4"),
        spec(),
        settings(),
        &mut Decider::default(),
    )
    .unwrap();
    assert_eq!(stage.lattice_step().as_str(), "lattice-1");
}

fn thin(radius: &str) -> LatticeSpec {
    LatticeSpec {
        radius: serde_json::from_str::<Grade>(&format!(r#"{{"constant":"{radius}"}}"#)).unwrap(),
        ..spec()
    }
}

fn thick() -> LatticeSpec {
    thin("0.3")
}
fn loose() -> Settings {
    Settings {
        accuracy: q("1"),
        ..settings()
    }
}

#[test]
fn a_thin_beam_below_the_minimum_detail_refuses_instead_of_refining() {
    let coarse = Settings {
        accuracy: q("1"),
        minimum_detail: q("1"),
        ..settings()
    };
    let stage = Stage::new(
        &body("block", "2"),
        thin("0.05"),
        coarse,
        &mut Decider::default(),
    )
    .unwrap();
    let error = stage.export(&mut Decider::default()).unwrap_err();
    assert_eq!(error.code, "field.print-unresolved");
    assert!(error.detail.contains("minimum detail"), "{}", error.detail);
}

#[test]
fn the_same_thin_beam_exports_when_the_minimum_detail_allows_the_refinement() {
    let fine = Settings {
        accuracy: q("1"),
        minimum_detail: q("0.01"),
        ..settings()
    };
    let stage = Stage::new(
        &body("block", "2"),
        thin("0.05"),
        fine,
        &mut Decider::default(),
    )
    .unwrap();
    let output = stage.export(&mut Decider::default()).unwrap();
    assert_eq!(output.verdict.components, 1);
    assert!(output.verdict.max_deviation <= q("1"));
}

#[test]
fn a_connected_lattice_exports_with_an_independently_checked_verdict() {
    let stage = Stage::new(
        &body("block", "2"),
        thick(),
        loose(),
        &mut Decider::default(),
    )
    .unwrap();
    let output = stage.export(&mut Decider::default()).unwrap();
    let v = &output.verdict;
    assert!(v.closed && v.non_self_intersecting && v.connected);
    assert!(v.max_deviation <= loose().accuracy);
    assert_eq!(output.mesh.component_count(), 1);
    assert_eq!(self_intersections(&output.mesh), 0);
    let mut stl = Vec::new();
    output.mesh.write_stl(&mut stl).unwrap();
    assert_eq!(stl.len(), 84 + 50 * output.mesh.triangles.len());
}

#[test]
fn separate_struts_refuse_as_disconnected() {
    let cell = serde_json::from_str::<LatticeCell>(
        r#"{"custom":{"nodes":[["0.25","0.25","0.25"],["0.75","0.75","0.75"]],"beams":[[0,1]]}}"#,
    )
    .unwrap();
    let stage = Stage::new(
        &body("block", "2"),
        LatticeSpec { cell, ..spec() },
        loose(),
        &mut Decider::default(),
    )
    .unwrap();
    let error = stage.export(&mut Decider::default()).unwrap_err();
    assert_eq!(error.code, "field.print-disconnected");
}

/// Counts triangle pairs that cross properly: an edge of one passes through the
/// interior of the other, tested in binary64 with a bounding-box prefilter. Pairs that
/// share an edge are skipped; pairs sharing one vertex test only their opposite edges.
fn self_intersections(mesh: &catallax_field::printing::Mesh) -> usize {
    type V = [f64; 3];
    let p = |i: usize| mesh.positions[i].map(f64::from);
    let sub = |a: V, b: V| [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
    let cross = |a: V, b: V| {
        [
            a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0],
        ]
    };
    let dot = |a: V, b: V| a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
    // Does segment (s, e) pass through the interior of triangle (a, b, c)?
    let crosses = |s: V, e: V, t: [V; 3]| {
        let n = cross(sub(t[1], t[0]), sub(t[2], t[0]));
        let (ds, de) = (dot(n, sub(s, t[0])), dot(n, sub(e, t[0])));
        if ds * de >= 0.0 {
            return false;
        }
        let u = ds / (ds - de);
        let x = [
            s[0] + u * (e[0] - s[0]),
            s[1] + u * (e[1] - s[1]),
            s[2] + u * (e[2] - s[2]),
        ];
        (0..3).all(|k| {
            let (a, b) = (t[k], t[(k + 1) % 3]);
            dot(n, cross(sub(b, a), sub(x, a))) > 0.0
        })
    };
    let boxes: Vec<(V, V)> = mesh
        .triangles
        .iter()
        .map(|t| {
            let ps = t.map(p);
            (
                std::array::from_fn(|k| ps.iter().map(|v| v[k]).fold(f64::MAX, f64::min)),
                std::array::from_fn(|k| ps.iter().map(|v| v[k]).fold(f64::MIN, f64::max)),
            )
        })
        .collect();
    let mut count = 0;
    for i in 0..mesh.triangles.len() {
        for j in i + 1..mesh.triangles.len() {
            if (0..3).any(|k| boxes[i].1[k] < boxes[j].0[k] || boxes[j].1[k] < boxes[i].0[k]) {
                continue;
            }
            let (a, b) = (mesh.triangles[i], mesh.triangles[j]);
            let shared = a.iter().filter(|v| b.contains(v)).count();
            if shared >= 2 {
                continue;
            }
            let (ta, tb) = (a.map(p), b.map(p));
            let opposite = |t: [usize; 3], u: [usize; 3]| -> Vec<(usize, usize)> {
                if shared == 1 {
                    let k = (0..3).find(|k| u.contains(&t[*k])).unwrap();
                    vec![((k + 1) % 3, (k + 2) % 3)]
                } else {
                    (0..3).map(|k| (k, (k + 1) % 3)).collect()
                }
            };
            let hit = opposite(a, b)
                .into_iter()
                .any(|(x, y)| crosses(ta[x], ta[y], tb))
                || opposite(b, a)
                    .into_iter()
                    .any(|(x, y)| crosses(tb[x], tb[y], ta));
            count += usize::from(hit);
        }
    }
    count
}

#[test]
fn the_preview_is_coarse_unproved_and_independent_of_the_export_settings() {
    use catallax_field::{display::PreviewOptions, lattice_stage::Proof};
    let stage = Stage::new(
        &body("block", "2"),
        thick(),
        loose(),
        &mut Decider::default(),
    )
    .unwrap();
    let options = PreviewOptions {
        tolerance: q("2"),
        max_vertices: 100_000,
        max_cells: 100_000,
    };
    let preview = stage.preview(&options, &mut Decider::default()).unwrap();
    assert_eq!(preview.proof, Proof::Unproved);
    assert!(!preview.mesh.triangles.is_empty());
    let output = stage.export(&mut Decider::default()).unwrap();
    assert!(preview.mesh.max_deviation >= output.verdict.max_deviation);
    assert!(preview.mesh.triangles.len() <= output.mesh.triangles.len());
}

#[test]
fn a_preview_over_its_cell_limit_refuses_without_a_partial_mesh() {
    use catallax_field::display::PreviewOptions;
    let stage = Stage::new(
        &body("block", "2"),
        thick(),
        loose(),
        &mut Decider::default(),
    )
    .unwrap();
    let options = PreviewOptions {
        tolerance: q("2"),
        max_vertices: 100_000,
        max_cells: 1,
    };
    let error = stage
        .preview(&options, &mut Decider::default())
        .unwrap_err();
    assert_eq!(error.code, "field.print-budget");
}

fn acceptance_stage(size: [&str; 3], translation: [&str; 3]) -> Stage {
    let recipe = Recipe::parse(&format!(
        r#"{{"format":"catallax-recipe/3","name":"V7","steps":[{{"id":"box","op":"block","size":{size:?},"pose":{{"rotation":["1","0","0","0"],"translation":{translation:?}}}}}],"result":"box"}}"#
    )).unwrap();
    Stage::new(
        &recipe,
        LatticeSpec {
            cell: LatticeCell::Cubic,
            size: Number::parse("5").unwrap(),
            radius: serde_json::from_str(r#"{"constant":"1"}"#).unwrap(),
            frame: spec().frame,
        },
        Settings {
            accuracy: q("0.1"),
            minimum_detail: q("0.5"),
            ..settings()
        },
        &mut Decider::default(),
    )
    .unwrap()
}

fn check_acceptance_mesh(
    positions: &[[f32; 3]],
    triangles: &[[usize; 3]],
    bounds: [[&str; 2]; 3],
    deviation: &Rational,
    euler: i64,
) {
    use std::collections::{BTreeMap, BTreeSet};
    let mut welded = BTreeMap::new();
    let ids: Vec<usize> = positions
        .iter()
        .map(|p| {
            let key = p.map(|x| {
                if x.to_bits() & 0x7fff_ffff == 0 {
                    0
                } else {
                    x.to_bits()
                }
            });
            let next = welded.len();
            *welded.entry(key).or_insert(next)
        })
        .collect();
    let mut edges = BTreeSet::new();
    let mut adjacency = vec![Vec::new(); welded.len()];
    for t in triangles {
        let [a, b, c] = t.map(|i| ids[i]);
        assert!(a != b && b != c && a != c);
        for (a, b) in [(a, b), (b, c), (c, a)] {
            edges.insert((a.min(b), a.max(b)));
            adjacency[a].push(b);
            adjacency[b].push(a);
        }
    }
    assert_eq!(
        welded.len() as i64 - edges.len() as i64 + triangles.len() as i64,
        euler
    );
    let mut seen = vec![false; welded.len()];
    let mut pending = vec![0];
    while let Some(i) = pending.pop() {
        if seen[i] {
            continue;
        }
        seen[i] = true;
        pending.extend(&adjacency[i]);
    }
    assert!(
        seen.into_iter().all(|v| v),
        "welded mesh must have one component"
    );
    for p in positions {
        for i in 0..3 {
            let coordinate = Rational::from_f64_exact(p[i] as f64).unwrap();
            assert!(coordinate >= &q(bounds[i][0]) - deviation);
            assert!(coordinate <= &q(bounds[i][1]) + deviation);
        }
    }
}

fn check_acceptance_verdict(output: &catallax_field::lattice_stage::Output) {
    let v = &output.verdict;
    assert!(v.closed && v.non_self_intersecting && v.connected);
    assert_eq!(v.components, 1);
    assert!(v.max_deviation > q("0") && v.max_deviation <= q("0.1"));
    assert_eq!(output.mesh.max_deviation, v.max_deviation);
    assert!(output.mesh.cells <= 2_000_000);
    assert!(output.mesh.positions.len() <= 1_000_000);
}

/// K1, fixed by the independent V7 designer: four nodes and one handle.
#[test]
fn k1_labelled_stage_file_contains_its_proven_four_cell_mesh() {
    let stage = acceptance_stage(["10", "10", "5"], ["-2.5", "-2.5", "-2.5"]);
    let output = stage.export(&mut Decider::default()).unwrap();
    check_acceptance_verdict(&output);
    let before = output.verdict.max_deviation.clone();
    let (bytes, output) = output.labelled_stl().unwrap();
    check_acceptance_verdict(&output);
    assert_eq!(&bytes[..8], b"CXPRINT1");
    assert_eq!(&bytes[8..12], &[2, 1, 1, 1]);
    assert_eq!(i64::from_le_bytes(bytes[12..20].try_into().unwrap()), 1);
    assert_eq!(u64::from_le_bytes(bytes[20..28].try_into().unwrap()), 10);
    assert_eq!(
        u64::from_le_bytes(bytes[28..36].try_into().unwrap()),
        2_000_000
    );
    assert_eq!(
        u64::from_le_bytes(bytes[36..44].try_into().unwrap()),
        1_000_000
    );
    let contribution = |start| {
        Rational::from_f64_exact(f64::from_le_bytes(
            bytes[start..start + 8].try_into().unwrap(),
        ))
        .unwrap()
    };
    let bound = contribution(68);
    assert!(bound >= before && bound <= q("0.1"));
    assert_eq!(bound, output.verdict.max_deviation);
    assert_eq!(contribution(44), q("0"));
    assert_eq!(contribution(60), q("0"));
    assert!(&(&contribution(44) + &contribution(52)) + &contribution(60) <= bound);
    assert_eq!(&bytes[76..80], &[0; 4]);
    let count = u32::from_le_bytes(bytes[80..84].try_into().unwrap()) as usize;
    assert_eq!(count, output.mesh.triangles.len());
    assert_eq!(bytes.len(), 84 + 50 * count);
    let mut positions = Vec::new();
    let mut triangles = Vec::new();
    for facet in bytes[84..].as_chunks::<50>().0 {
        let base = positions.len();
        for v in 0..3 {
            positions.push(std::array::from_fn(|i| {
                let start = 12 + v * 12 + i * 4;
                f32::from_le_bytes(facet[start..start + 4].try_into().unwrap())
            }));
        }
        triangles.push([base, base + 1, base + 2]);
    }
    check_acceptance_mesh(
        &positions,
        &triangles,
        [["-2.5", "7.5"], ["-2.5", "7.5"], ["-2.5", "2.5"]],
        &bound,
        0,
    );
}

/// K2, fixed by the independent V7 designer: S1 has 30 nodes and 30 handles.
#[test]
#[cfg_attr(
    debug_assertions,
    ignore = "V7 K2 is a release-mode proof-cost acceptance"
)]
fn k2_thirty_cell_lattice_proves_at_point_one_with_default_caps() {
    let stage = acceptance_stage(["25", "15", "10"], ["-12.5", "-7.5", "-2.5"]);
    let output = stage.export(&mut Decider::default()).unwrap();
    check_acceptance_verdict(&output);
    check_acceptance_mesh(
        &output.mesh.positions,
        &output.mesh.triangles,
        [["-12.5", "12.5"], ["-7.5", "7.5"], ["-2.5", "7.5"]],
        &output.verdict.max_deviation,
        -58,
    );
}

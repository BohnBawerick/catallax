//! The final lattice stage: request building, refusals and rebuild from an edited body.
use catallax_field::{
    Location, Workspace,
    lattice_stage::{GyroidSpec, LatticeSpec, Settings, Stage, check_no_cad_after_lattice},
};
use catallax_numbers::{Decider, Rational};
use catallax_recipe::{Grade, GyroidKind, LatticeCell, Number, Pose, Recipe};
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
    assert_eq!(output.verdict().components, 1);
    assert!(output.verdict().max_deviation <= q("1"));
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
    let v = output.verdict();
    assert!(v.closed && v.non_self_intersecting && v.connected);
    assert!(v.max_deviation <= loose().accuracy);
    assert_eq!(output.mesh().component_count(), 1);
    assert_eq!(self_intersections(output.mesh()), 0);
    let mut stl = Vec::new();
    output.mesh().write_stl(&mut stl).unwrap();
    assert_eq!(stl.len(), 84 + 50 * output.mesh().triangles.len());
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
    assert!(preview.mesh.max_deviation >= output.verdict().max_deviation);
    assert!(preview.mesh.triangles.len() <= output.mesh().triangles.len());
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
    let v = output.verdict();
    assert!(v.closed && v.non_self_intersecting && v.connected);
    assert_eq!(v.components, 1);
    assert!(v.max_deviation > q("0") && v.max_deviation <= q("0.1"));
    assert_eq!(output.mesh().max_deviation, v.max_deviation);
    assert!(output.mesh().cells <= 2_000_000);
    assert!(output.mesh().positions.len() <= 1_000_000);
}

/// K1, fixed by the independent V7 designer: four nodes and one handle.
#[test]
fn k1_labelled_stage_file_contains_its_proven_four_cell_mesh() {
    let stage = acceptance_stage(["10", "10", "5"], ["-2.5", "-2.5", "-2.5"]);
    let output = stage.export(&mut Decider::default()).unwrap();
    check_acceptance_verdict(&output);
    let before = output.verdict().max_deviation.clone();
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
    assert_eq!(bound, output.verdict().max_deviation);
    assert_eq!(contribution(44), q("0"));
    assert_eq!(contribution(60), q("0"));
    assert!(&(&contribution(44) + &contribution(52)) + &contribution(60) <= bound);
    assert_eq!(&bytes[76..80], &[0; 4]);
    let count = u32::from_le_bytes(bytes[80..84].try_into().unwrap()) as usize;
    assert_eq!(count, output.mesh().triangles.len());
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
        &output.mesh().positions,
        &output.mesh().triangles,
        [["-12.5", "12.5"], ["-7.5", "7.5"], ["-2.5", "7.5"]],
        &output.verdict().max_deviation,
        -58,
    );
}

fn profile_body(depth: &str) -> Recipe {
    Recipe::parse(&format!(
        r#"{{"format":"catallax-recipe/3","name":"S1 profile","steps":[{{"id":"profile","op":"sketch","plane":{{"rotation":["1","0","0","0"],"translation":["0","0","-2.5"]}},"loops":[[{{"kind":"line","id":"a","to":["12.5","-7.5"]}},{{"kind":"line","id":"b","to":["12.5","7.5"]}},{{"kind":"line","id":"c","to":["-12.5","7.5"]}},{{"kind":"line","id":"d","to":["-12.5","-7.5"]}}]]}},{{"id":"boss","op":"extrude","profile":"profile","distance":"{depth}"}}],"result":"boss"}}"#
    ))
    .unwrap()
}

/// Reads the labelled header fields the V7 tests fix: requested accuracy and caps.
fn labelled_header(bytes: &[u8]) -> (i64, u64, u64, u64) {
    let word = |start: usize| u64::from_le_bytes(bytes[start..start + 8].try_into().unwrap());
    (word(12) as i64, word(20), word(28), word(36))
}

/// V7 tenth-millimetre increment, fixed by the independent designer: S1's three profile
/// states certify at 1/10 mm within the default caps and the 100M work budget.
#[test]
#[cfg_attr(debug_assertions, ignore = "V7 S1 profile exports run in release mode")]
fn s1_profile_extrusion_rebuild_and_density_edit_certify_at_point_one() {
    use catallax_field::lattice_stage::{Proof, default_preview_options};
    for (depth, radius, top, euler) in [
        ("10", "1", "7.5", -58),
        ("15", "1", "12.5", -102),
        ("15", "1.5", "12.5", -102),
    ] {
        let (_, mut spec, settings) =
            acceptance_stage(["25", "15", "10"], ["-12.5", "-7.5", "-2.5"]).return_before_lattice();
        spec.radius = serde_json::from_str(&format!(r#"{{"constant":"{radius}"}}"#)).unwrap();
        assert_eq!(settings.accuracy, q("1/10"));
        let (output, work) = catallax_numbers::with_operation_budget(100_000_000, || {
            let mut dec = Decider::default();
            let stage = Stage::new(&profile_body(depth), spec, settings, &mut dec).unwrap();
            let controls = default_preview_options(stage.plan(), &mut dec).unwrap();
            assert_eq!(
                stage.preview(&controls, &mut dec).unwrap().proof,
                Proof::Unproved
            );
            stage.export(&mut dec).unwrap()
        });
        assert!(!work.exhausted, "{depth} {radius}: {work:?}");
        check_acceptance_verdict(&output);
        let (bytes, output) = output.labelled_stl().unwrap();
        check_acceptance_verdict(&output);
        assert_eq!(labelled_header(&bytes), (1, 10, 2_000_000, 1_000_000));
        check_acceptance_mesh(
            &output.mesh().positions,
            &output.mesh().triangles,
            [["-12.5", "12.5"], ["-7.5", "7.5"], ["-2.5", top]],
            &output.verdict().max_deviation,
            euler,
        );
    }
}

/// S1 export contract, fixed by the independent designer: each profile state keeps its
/// baseline mesh and labelled bytes within the app's 20M work allowance and default caps.
#[test]
#[cfg_attr(debug_assertions, ignore = "V7 S1 profile exports run in release mode")]
fn s1_profile_exports_keep_their_bytes_within_the_app_work_allowance() {
    for (depth, radius, top, euler, vertices, triangles, length, sha256) in [
        (
            "10",
            "1",
            "7.5",
            -58,
            149_522,
            299_160,
            14_958_084,
            "79790ed211669b811cd3b0086f24dab3a333ea779bdecf9000441ca138de3fd0",
        ),
        (
            "15",
            "1",
            "12.5",
            -102,
            217_518,
            435_240,
            21_762_084,
            "a6cc0e31ed490b193a0cb45ea3e5f6e149003c8ba71fd5e49fdf1d09d742a386",
        ),
        (
            "15",
            "1.5",
            "12.5",
            -102,
            263_346,
            526_896,
            26_344_884,
            "dedc73ca8eb8d1288b2658f1406296badbe06d202526b9b61985c968170cfe30",
        ),
    ] {
        let (_, mut spec, settings) =
            acceptance_stage(["25", "15", "10"], ["-12.5", "-7.5", "-2.5"]).return_before_lattice();
        spec.radius = serde_json::from_str(&format!(r#"{{"constant":"{radius}"}}"#)).unwrap();
        let mut dec = Decider::default();
        let stage = Stage::new(&profile_body(depth), spec, settings, &mut dec).unwrap();
        let (exported, work) = catallax_numbers::with_operation_budget(20_000_000, || {
            stage.export(&mut dec)?.labelled_stl()
        });
        assert!(!work.exhausted, "{depth} {radius}: {work:?}");
        let (bytes, output) = exported.unwrap();
        check_acceptance_verdict(&output);
        assert_eq!(labelled_header(&bytes), (1, 10, 2_000_000, 1_000_000));
        check_acceptance_mesh(
            &output.mesh().positions,
            &output.mesh().triangles,
            [["-12.5", "12.5"], ["-7.5", "7.5"], ["-2.5", top]],
            &output.verdict().max_deviation,
            euler,
        );
        assert_eq!(output.mesh().positions.len(), vertices);
        assert_eq!(output.mesh().triangles.len(), triangles);
        assert_eq!(bytes.len(), length);
        assert_eq!(catallax_import::sha256_hex(&bytes), sha256);
    }
}

/// Every acceptable S1 depth-15 mesh is one closed surface with Euler characteristic
/// -102, so E = 3(V + 102) <= V(V - 1)/2 and it needs at least 29 vertices. A 28-vertex
/// cap leaves no acceptable output: the stage must refuse, whatever its proof method.
#[test]
#[cfg_attr(debug_assertions, ignore = "V7 S1 profile exports run in release mode")]
fn s1_depth_fifteen_under_a_provably_too_small_vertex_cap_refuses() {
    let (_, mut spec, mut settings) =
        acceptance_stage(["25", "15", "10"], ["-12.5", "-7.5", "-2.5"]).return_before_lattice();
    spec.radius = serde_json::from_str(r#"{"constant":"1.5"}"#).unwrap();
    settings.max_vertices = 28;
    let stage = Stage::new(&profile_body("15"), spec, settings, &mut Decider::default()).unwrap();
    let error = stage.export(&mut Decider::default()).unwrap_err();
    assert_eq!(error.code, "field.print-budget", "{}", error.detail);
}

/// Under a small cell cap the stage either certifies within that cap at the requested
/// accuracy and records it, or refuses with the budget code. It never answers with a
/// coarser mesh or a larger cap.
#[test]
#[cfg_attr(debug_assertions, ignore = "V7 S1 profile exports run in release mode")]
fn s1_depth_fifteen_under_a_small_cell_cap_certifies_within_it_or_refuses() {
    let (_, mut spec, mut settings) =
        acceptance_stage(["25", "15", "10"], ["-12.5", "-7.5", "-2.5"]).return_before_lattice();
    spec.radius = serde_json::from_str(r#"{"constant":"1.5"}"#).unwrap();
    settings.max_cells = 1_000;
    let stage = Stage::new(&profile_body("15"), spec, settings, &mut Decider::default()).unwrap();
    match stage.export(&mut Decider::default()) {
        Err(error) => assert_eq!(error.code, "field.print-budget", "{}", error.detail),
        Ok(output) => {
            let v = output.verdict();
            assert!(v.closed && v.non_self_intersecting && v.connected && v.components == 1);
            assert!(v.max_deviation > q("0") && v.max_deviation <= q("1/10"));
            assert!(output.mesh().cells <= 1_000);
            let (bytes, _) = output.labelled_stl().unwrap();
            assert_eq!(labelled_header(&bytes), (1, 10, 1_000, 1_000_000));
        }
    }
}

fn gyroid_spec(level: &str) -> GyroidSpec {
    GyroidSpec {
        period: Number::parse("10").unwrap(),
        level: serde_json::from_str(&format!(r#"{{"constant":"{level}"}}"#)).unwrap(),
        kind: GyroidKind::Network,
        frame: spec().frame,
    }
}
fn gyroid_settings() -> Settings {
    Settings {
        accuracy: q("1/10"),
        minimum_detail: q("1/2"),
        max_cells: 2_000_000,
        max_vertices: 1_000_000,
    }
}
/// The sketch-plus-extrude body of the V7c scenario: a 20 by 20 rectangle at z = -5.
fn gyroid_body(depth: &str) -> Recipe {
    Recipe::parse(&format!(
        r#"{{"format":"catallax-recipe/3","name":"V7c body","steps":[{{"id":"profile","op":"sketch","plane":{{"rotation":["1","0","0","0"],"translation":["0","0","-5"]}},"loops":[[{{"kind":"line","id":"a","to":["10","-10"]}},{{"kind":"line","id":"b","to":["10","10"]}},{{"kind":"line","id":"c","to":["-10","10"]}},{{"kind":"line","id":"d","to":["-10","-10"]}}]]}},{{"id":"boss","op":"extrude","profile":"profile","distance":"{depth}"}}],"result":"boss"}}"#
    ))
    .unwrap()
}

/// An independent reading of one `catallax-print-stl/1` gyroid file, from the V7c
/// design's D1 and D3-D6. It uses no kernel value but the file's own bytes.
mod gyroid_oracle {
    use super::*;
    use std::collections::{BTreeMap, BTreeSet};

    type V = [f64; 3];
    const PERIOD: f64 = 10.0;
    const LIPSCHITZ: f64 = 2.176_559_237_081_062;

    fn g(p: V) -> f64 {
        let [x, y, z] = p.map(|v| v * std::f64::consts::TAU / PERIOD);
        x.sin() * y.cos() + y.sin() * z.cos() + z.sin() * x.cos()
    }
    fn sub(a: V, b: V) -> V {
        [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
    }
    fn dot(a: V, b: V) -> f64 {
        a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
    }
    fn cross(a: V, b: V) -> V {
        [
            a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0],
        ]
    }
    fn distance(p: V, t: [V; 3]) -> f64 {
        let along = |a: V, b: V| {
            let e = sub(b, a);
            let u = (dot(sub(p, a), e) / dot(e, e)).clamp(0.0, 1.0);
            let d = sub(p, [a[0] + u * e[0], a[1] + u * e[1], a[2] + u * e[2]]);
            dot(d, d)
        };
        let n = cross(sub(t[1], t[0]), sub(t[2], t[0]));
        let inside = (0..3).all(|k| dot(n, cross(sub(t[(k + 1) % 3], t[k]), sub(p, t[k]))) >= 0.0);
        let squared = if inside {
            dot(sub(p, t[0]), n).powi(2) / dot(n, n)
        } else {
            along(t[0], t[1])
                .min(along(t[1], t[2]))
                .min(along(t[2], t[0]))
        };
        squared.sqrt()
    }

    /// D4's fixed boundary samples: cyclic permutations of the level's seeds, moved by
    /// whole periods within the box. The negative seeds of the corrected G2 also need
    /// their positive copies, so its translations run one period further in x and y.
    fn samples(level: &Rational, top: f64) -> Vec<V> {
        let ([a, b], shifts): ([f64; 2], &[f64]) = if level.is_zero() {
            ([0.0, 5.0], &[-10.0, 0.0])
        } else {
            assert_eq!(*level, q("-1/2"));
            ([-5.0 / 6.0, -25.0 / 6.0], &[-10.0, 0.0, 10.0])
        };
        let mut points = Vec::new();
        for seed in [a, b] {
            for shift in 0..3 {
                let mut p = [0.0; 3];
                p[shift] = seed;
                for &i in shifts {
                    for &j in shifts {
                        for k in [0.0, 10.0] {
                            let t: V = [p[0] + i, p[1] + j, p[2] + k];
                            let inside = t[0].abs() <= 10.0
                                && t[1].abs() <= 10.0
                                && t[2] >= -5.0
                                && t[2] <= top;
                            if inside && !points.contains(&t) {
                                points.push(t);
                            }
                        }
                    }
                }
            }
        }
        points
    }

    /// D1's open material patches on each face of the box, with their outward normals.
    fn faces(top: f64) -> Vec<(V, V)> {
        let mut points = vec![
            ([10.0, -2.5, 0.0], [1.0, 0.0, 0.0]),
            ([-10.0, -2.5, 0.0], [-1.0, 0.0, 0.0]),
            ([0.0, 10.0, -2.5], [0.0, 1.0, 0.0]),
            ([0.0, -10.0, -2.5], [0.0, -1.0, 0.0]),
            ([-2.5, 0.0, -5.0], [0.0, 0.0, -1.0]),
            ([-2.5, 0.0, top], [0.0, 0.0, 1.0]),
        ];
        for (p, _) in &mut points {
            assert!(g(*p) < -0.999, "a face patch of D1 is material");
        }
        points
    }

    /// Checks D1 and D3-D6 on written bytes, for a state's level and upper face.
    pub(super) fn check(bytes: &[u8], level: &Rational, top: f64, caps: (u64, u64)) {
        let word = |at: usize| u64::from_le_bytes(bytes[at..at + 8].try_into().unwrap());
        let real = |at: usize| f64::from_le_bytes(bytes[at..at + 8].try_into().unwrap());
        assert_eq!(&bytes[..8], b"CXPRINT1");
        assert_eq!(&bytes[8..12], &[2, 1, 1, 1]);
        assert_eq!((word(12), word(20)), (1, 10));
        assert_eq!((word(28), word(36)), caps);
        let parts = [real(44), real(52), real(60), real(68)];
        assert!(parts.iter().all(|v| v.is_finite() && *v >= 0.0));
        let exact = parts.map(|v| Rational::from_f64_exact(v).unwrap());
        let w = &exact[3];
        assert!(*w > q("0") && *w <= q("1/10"));
        assert!(&(&exact[0] + &exact[1]) + &exact[2] <= *w);
        assert_eq!(&bytes[76..80], &[0; 4]);
        let count = u32::from_le_bytes(bytes[80..84].try_into().unwrap()) as usize;
        assert!(count > 0);
        assert_eq!(bytes.len(), 84 + 50 * count);

        let mut welded = BTreeMap::new();
        let mut points: Vec<V> = Vec::new();
        let mut facets = Vec::with_capacity(count);
        for facet in bytes[84..].as_chunks::<50>().0 {
            let mut ids = [0usize; 3];
            for (v, id) in ids.iter_mut().enumerate() {
                let c: [f32; 3] = std::array::from_fn(|i| {
                    let at = 12 + v * 12 + i * 4;
                    f32::from_le_bytes(facet[at..at + 4].try_into().unwrap())
                });
                assert!(c.iter().all(|x| x.is_finite()));
                let key = c.map(|x| if x == 0.0 { 0 } else { x.to_bits() });
                let next = points.len();
                *id = *welded.entry(key).or_insert(next);
                if *id == next {
                    points.push(c.map(f64::from));
                }
            }
            assert!(ids[0] != ids[1] && ids[1] != ids[2] && ids[2] != ids[0]);
            facets.push(ids);
        }
        assert!(points.len() as u64 <= caps.1);

        // Each undirected edge has exactly one facet in each direction.
        let mut directed = BTreeSet::new();
        let mut links: Vec<Vec<(usize, usize)>> = vec![Vec::new(); points.len()];
        for &[a, b, c] in &facets {
            for (p, x, y) in [(a, b, c), (b, c, a), (c, a, b)] {
                assert!(directed.insert((p, x)), "edge repeated in one direction");
                links[p].push((x, y));
            }
        }
        for &(a, b) in &directed {
            assert!(directed.contains(&(b, a)), "edge lacks its opposite");
        }
        // Each vertex link is a single cycle, and the surface is one component.
        for link in &links {
            let next: BTreeMap<usize, usize> = link.iter().copied().collect();
            assert_eq!(next.len(), link.len());
            let mut walk = link[0].0;
            for _ in 1..link.len() {
                walk = next[&walk];
                assert_ne!(walk, link[0].0, "vertex link has two cycles");
            }
            assert_eq!(next[&walk], link[0].0);
        }
        let mut seen = vec![false; points.len()];
        let mut pending = vec![0];
        while let Some(i) = pending.pop() {
            if !std::mem::replace(&mut seen[i], true) {
                pending.extend(links[i].iter().map(|l| l.0));
            }
        }
        assert!(seen.iter().all(|v| *v), "one edge-connected component");

        // D1: the box bound; D3: every vertex near the level or a clipping face.
        let w = w.to_f64_nearest();
        let c = level.to_f64_nearest();
        let reach = LIPSCHITZ * w;
        let (low, high) = ([-10.0, -10.0, -5.0], [10.0, 10.0, top]);
        for p in &points {
            for i in 0..3 {
                assert!(p[i] >= low[i] - w && p[i] <= high[i] + w, "{p:?}");
            }
            let value = g(*p);
            let near_face = (0..3).any(|i| p[i] <= low[i] + w || p[i] >= high[i] - w);
            assert!(
                (value - c).abs() <= reach + 1e-8 || (value <= c + reach + 1e-8 && near_face),
                "{p:?} has g = {value}"
            );
        }

        // D4: fixed boundary samples, by point-to-triangle distance.
        let triangles: Vec<[V; 3]> = facets.iter().map(|t| t.map(|v| points[v])).collect();
        let samples = samples(level, top);
        if !level.is_zero() {
            assert_eq!(samples.len(), 84);
        }
        for p in samples {
            let nearest = triangles
                .iter()
                .map(|t| distance(p, *t))
                .fold(f64::INFINITY, f64::min);
            assert!(nearest <= w + 1e-8, "{p:?} is {nearest} from the mesh");
        }
        // D1 and D5: the box faces are reached and face outward there.
        for (p, normal) in faces(top) {
            let (nearest, facet) = triangles
                .iter()
                .map(|t| (distance(p, *t), t))
                .min_by(|a, b| a.0.total_cmp(&b.0))
                .unwrap();
            assert!(nearest <= w + 1e-8, "{p:?} is {nearest} from the mesh");
            let n = cross(sub(facet[1], facet[0]), sub(facet[2], facet[0]));
            assert!(dot(n, normal) > 0.0, "{p:?} faces inward");
        }
    }
}

/// K-G2, fixed by the independent V7c designer: the rebuilt body and the level both reach
/// the geometry, at `1/10` mm within the default caps and one 100M allowance covering
/// compilation, preview, export and labelling. G2 uses level `-1/2`: at `1/2` this box
/// holds a separate exact corner piece, so that level has no one-component export.
#[test]
#[cfg_attr(debug_assertions, ignore = "V7c gyroid exports run in release mode")]
fn gyroid_profile_states_certify_at_point_one() {
    use catallax_field::lattice_stage::{Proof, default_preview_options};
    for (depth, level, top) in [("10", "0", 5.0), ("20", "0", 15.0), ("20", "-1/2", 15.0)] {
        let mut preview_time = std::time::Duration::ZERO;
        let started = std::time::Instant::now();
        let ((output, bytes), work) = catallax_numbers::with_operation_budget(100_000_000, || {
            let mut dec = Decider::default();
            let stage = Stage::new(
                &gyroid_body(depth),
                gyroid_spec(level),
                gyroid_settings(),
                &mut dec,
            )
            .unwrap();
            let controls = default_preview_options(stage.plan(), &mut dec).unwrap();
            let previewing = std::time::Instant::now();
            let preview = stage.preview(&controls, &mut dec).unwrap();
            preview_time = previewing.elapsed();
            assert_eq!(preview.proof, Proof::Unproved);
            assert!(!preview.mesh.triangles.is_empty());
            let output = stage.export(&mut dec).unwrap();
            check_acceptance_verdict(&output);
            let (bytes, output) = output.labelled_stl().unwrap();
            (output, bytes)
        });
        let total_time = started.elapsed();
        assert!(!work.exhausted, "{depth} {level}: {work:?}");
        check_acceptance_verdict(&output);
        assert_eq!(labelled_header(&bytes), (1, 10, 2_000_000, 1_000_000));
        gyroid_oracle::check(&bytes, &q(level), top, (2_000_000, 1_000_000));
        println!(
            "gyroid depth {depth} level {level}: cells {} vertices {} facets {} work {} preview {preview_time:?} total {total_time:?}",
            output.mesh().cells,
            output.mesh().positions.len(),
            output.mesh().triangles.len(),
            work.used,
        );
    }
}

/// K-G3, fixed by the independent V7c designer: a smaller cap either holds a fully
/// certified output or refuses with the budget code. It never coarsens the request.
#[test]
#[cfg_attr(debug_assertions, ignore = "V7c gyroid exports run in release mode")]
fn gyroid_small_caps_never_change_the_request() {
    for (max_cells, max_vertices) in [(1_000, 1_000_000), (2_000_000, 4)] {
        let settings = Settings {
            max_cells,
            max_vertices,
            ..gyroid_settings()
        };
        let stage = Stage::new(
            &gyroid_body("20"),
            gyroid_spec("-1/2"),
            settings,
            &mut Decider::default(),
        )
        .unwrap();
        match stage.export(&mut Decider::default()) {
            Err(error) => assert_eq!(error.code, "field.print-budget", "{}", error.detail),
            Ok(output) => {
                let (bytes, _) = output.labelled_stl().unwrap();
                gyroid_oracle::check(
                    &bytes,
                    &q("-1/2"),
                    15.0,
                    (max_cells as u64, max_vertices as u64),
                );
            }
        }
    }
    for (max_cells, max_vertices) in [(0, 1_000_000), (2_000_000, 3)] {
        let settings = Settings {
            max_cells,
            max_vertices,
            ..gyroid_settings()
        };
        let error = Stage::new(
            &gyroid_body("20"),
            gyroid_spec("-1/2"),
            settings,
            &mut Decider::default(),
        )
        .unwrap_err();
        assert_eq!(error.code, "lattice.settings");
    }
}

/// Terminal rule: no second field or CAD step may follow a gyroid, in either order.
#[test]
fn a_gyroid_is_terminal_like_a_lattice() {
    let mut dec = Decider::default();
    let gyroid = Stage::new(
        &gyroid_body("10"),
        gyroid_spec("0"),
        gyroid_settings(),
        &mut dec,
    )
    .unwrap();
    assert_eq!(gyroid.lattice_step().as_str(), "gyroid");
    assert_eq!(gyroid.recipe().result, *gyroid.lattice_step());
    let error = Stage::new(gyroid.recipe(), spec(), settings(), &mut dec).unwrap_err();
    assert_eq!(error.code, "lattice.cad-after-lattice");
    let error = Stage::new(
        gyroid.recipe(),
        gyroid_spec("0"),
        gyroid_settings(),
        &mut dec,
    )
    .unwrap_err();
    assert_eq!(error.code, "lattice.cad-after-lattice");
    let lattice = first_again();
    let error = Stage::new(
        lattice.recipe(),
        gyroid_spec("0"),
        gyroid_settings(),
        &mut dec,
    )
    .unwrap_err();
    assert_eq!(error.code, "lattice.cad-after-lattice");
    let (body, kept, _) = gyroid.return_before_lattice();
    assert_eq!(body, gyroid_body("10"));
    assert_eq!(kept, gyroid_spec("0"));

    let reading = Recipe::parse(
        r#"{"format":"catallax-recipe/3","name":"x","steps":[{"id":"region","op":"block","size":["4","4","4"],"pose":{"rotation":["1","0","0","0"],"translation":["0","0","0"]}},{"id":"net","op":"gyroid","region":"region","period":"4","level":{"constant":"0"},"kind":"network","frame":{"rotation":["1","0","0","0"],"translation":["0","0","0"]}},{"id":"after","op":"offset","body":"net","distance":"0.1"}],"result":"after"}"#,
    )
    .unwrap();
    assert_eq!(
        check_no_cad_after_lattice(&reading).unwrap_err().code,
        "lattice.cad-after-lattice"
    );
}

/// The design's public print input: a plain block region through the generic print route,
/// which needs no stage. The same certificate applies, and no B-rep builder is involved.
#[test]
#[cfg_attr(debug_assertions, ignore = "V7c gyroid exports run in release mode")]
fn a_block_region_certifies_through_the_generic_print_route() {
    let recipe = Recipe::parse(
        r#"{"format":"catallax-recipe/3","name":"V7c public gyroid diagnostic","steps":[{"id":"region","op":"block","size":["20","20","10"],"pose":{"rotation":["1","0","0","0"],"translation":["-10","-10","-5"]}},{"id":"gyroid","op":"gyroid","region":"region","period":"10","level":{"constant":"0"},"kind":"network","frame":{"rotation":["1","0","0","0"],"translation":["0","0","0"]}}],"result":"gyroid"}"#,
    )
    .unwrap();
    let mut dec = Decider::default();
    let plan = catallax_field::Plan::compile_checked(&recipe, &mut dec).unwrap();
    let options = catallax_field::printing::Options {
        tolerance: q("1/10"),
        max_cells: 2_000_000,
        max_vertices: 1_000_000,
    };
    let mesh = catallax_field::printing::mesh_plan(&plan, &options, &mut dec).unwrap();
    assert_eq!(mesh.component_count(), 1);
    assert!(mesh.max_deviation > q("0") && mesh.max_deviation <= q("1/10"));
    assert!(mesh.cells <= 2_000_000 && mesh.positions.len() <= 1_000_000);
}

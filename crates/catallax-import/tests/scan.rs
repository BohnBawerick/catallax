//! Open scans as reference data, RECIPE.md section 4.19.1: holes, flipped
//! triangles, crossings, degenerate triangles and isolated samples are kept and
//! counted, never repaired, and the placement is exact.

use catallax_import::mesh::ClosedMesh;
use catallax_import::scan::{ScanItem, ScanReference};
use catallax_numbers::{Decider, Rational, Refusal};
use catallax_recipe::{Number, Pose};

fn q(s: &str) -> Rational {
    Rational::parse(s).unwrap()
}

fn read(bytes: &[u8]) -> Result<ScanReference, Refusal> {
    ScanReference::read(
        bytes,
        &Rational::one(),
        &Pose::identity(),
        &mut Decider::default(),
    )
}

/// An OBJ file with these vertices and one-based triangles.
fn obj(points: &[[&str; 3]], triangles: &[[usize; 3]]) -> Vec<u8> {
    let mut s = String::new();
    for p in points {
        s.push_str(&format!("v {} {} {}\n", p[0], p[1], p[2]));
    }
    for t in triangles {
        s.push_str(&format!("f {} {} {}\n", t[0], t[1], t[2]));
    }
    s.into_bytes()
}

const CUBE: [[&str; 3]; 8] = [
    ["0", "0", "0"],
    ["2", "0", "0"],
    ["0", "2", "0"],
    ["2", "2", "0"],
    ["0", "0", "2"],
    ["2", "0", "2"],
    ["0", "2", "2"],
    ["2", "2", "2"],
];

/// The closed cube's twelve outward triangles, one-based; the last two are its top.
fn cube_triangles() -> Vec<[usize; 3]> {
    let faces: [[usize; 4]; 6] = [
        [0, 2, 3, 1],
        [0, 1, 5, 4],
        [2, 6, 7, 3],
        [0, 4, 6, 2],
        [1, 3, 7, 5],
        [4, 5, 7, 6],
    ];
    faces
        .iter()
        .flat_map(|f| [[f[0], f[1], f[2]], [f[0], f[2], f[3]]])
        .map(|t| t.map(|v| v + 1))
        .collect()
}

#[test]
fn a_cube_with_its_top_missing_is_a_reference_but_not_a_solid() {
    let mut triangles = cube_triangles();
    triangles.truncate(10);
    let bytes = obj(&CUBE, &triangles);
    let refusal = ClosedMesh::read(&bytes, &Rational::one(), &mut Decider::default()).unwrap_err();
    assert_eq!(refusal.code, "import.mesh-open");
    let scan = read(&bytes).unwrap();
    let report = scan.report();
    assert_eq!(report.triangles, 10);
    assert_eq!(report.vertices, 8);
    assert_eq!(
        report.boundary_edges, 4,
        "the four rim edges of the missing top"
    );
    assert_eq!(report.inconsistent_edges, 0);
    assert_eq!(report.nonmanifold_edges, 0);
    assert_eq!(scan.sha256(), catallax_import::sha256_hex(&bytes));
    // Over the hole the nearest surface is the rim or a side, not a lid.
    assert_eq!(scan.fast_distance([1.0, 1.0, 2.0]), 1.0);
    assert_eq!(scan.fast_distance([1.0, 1.0, 1.0]), 1.0);
    assert_eq!(scan.fast_distance([1.0, 1.0, 0.5]), 0.5);
    assert_eq!(scan.fast_distance([3.0, 1.0, 1.0]), 1.0);
}

#[test]
fn flipped_crossing_and_degenerate_triangles_are_kept_and_counted() {
    let mut triangles = cube_triangles();
    triangles[0] = [triangles[0][0], triangles[0][2], triangles[0][1]];
    let mut points = CUBE.to_vec();
    // A crossing triangle through the cube's middle, a collinear one and a point.
    points.extend([
        ["1", "1", "-1"],
        ["1", "-1", "3"],
        ["1", "3", "3"],
        ["5", "0", "0"],
        ["6", "0", "0"],
        ["7", "0", "0"],
    ]);
    triangles.extend([[9, 10, 11], [12, 14, 13], [12, 12, 12]]);
    let bytes = obj(&points, &triangles);
    let scan = read(&bytes).unwrap();
    let report = scan.report();
    assert_eq!(report.triangles, 15);
    assert_eq!(
        report.inconsistent_edges, 3,
        "the flipped triangle's three edges"
    );
    assert_eq!(
        report.boundary_edges, 3,
        "the crossing triangle's three edges"
    );
    assert_eq!(report.degenerate_segments, 1);
    assert_eq!(report.degenerate_points, 1);
    assert_eq!(
        scan.items()[13],
        ScanItem::Segment([11, 12]),
        "the longest side holds the middle corner"
    );
    assert_eq!(scan.items()[14], ScanItem::Point(11));
    assert_eq!(scan.fast_distance([6.0, 3.0, 0.0]), 3.0);
    assert_eq!(scan.fast_distance([4.0, 0.0, 0.0]), 1.0);
}

#[test]
fn a_point_cloud_is_a_set_of_isolated_samples() {
    let mut text = String::new();
    for i in 0..5 {
        text.push_str(&format!("v {i} 0 0.5\n"));
    }
    // An exact repeat is one sample; nothing is joined by distance.
    text.push_str("v 4 0 0.5\nv 4 0 0.50000000001\n");
    let scan = read(text.as_bytes()).unwrap();
    assert_eq!(scan.report().isolated_points, 6);
    assert_eq!(scan.report().triangles, 0);
    assert_eq!(scan.items().len(), 6);
    assert!(scan.items().iter().all(|i| matches!(i, ScanItem::Point(_))));
    assert_eq!(scan.fast_distance([2.0, 0.0, 0.5]), 0.0);
    assert_eq!(scan.fast_distance([2.5, 0.0, 0.5]), 0.5);
    // Unused vertices beside triangles are samples unless they equal a corner.
    let mut bytes = obj(&CUBE, &cube_triangles());
    bytes.extend_from_slice(b"v 2 2 2\nv 9 9 9\n");
    let scan = read(&bytes).unwrap();
    assert_eq!(scan.report().isolated_points, 1);
    assert_eq!(scan.points().last().unwrap(), &["9", "9", "9"].map(q));
}

#[test]
fn placement_is_exact_scale_then_rational_rotation_then_translation() {
    let bytes = obj(
        &[["2", "0", "0"], ["0", "4", "0"], ["0", "0", "6"]],
        &[[1, 2, 3]],
    );
    let n = |s: &str| Number::new(q(s));
    let pose = Pose {
        // A quarter turn about +z, written at length 2: x -> y, y -> -x.
        rotation: [n("1"), n("0"), n("0"), n("1")],
        translation: catallax_recipe::Vector3([n("10"), n("1/3"), n("0")]),
    };
    let scan = ScanReference::read(&bytes, &q("1/2"), &pose, &mut Decider::default()).unwrap();
    assert_eq!(
        scan.points(),
        &[
            ["10", "4/3", "0"].map(q),
            ["8", "1/3", "0"].map(q),
            ["10", "1/3", "3"].map(q),
        ]
    );
    let zero = Pose {
        rotation: [n("0"), n("0"), n("0"), n("0")],
        ..Pose::identity()
    };
    let code = |scale: &str, pose: &Pose| {
        ScanReference::read(&bytes, &q(scale), pose, &mut Decider::default())
            .unwrap_err()
            .code
    };
    assert_eq!(code("1", &zero), "shapes.zero-quaternion");
    assert_eq!(code("0", &Pose::identity()), "import.scan-scale");
    assert_eq!(code("-1", &Pose::identity()), "import.scan-scale");
}

#[test]
fn empty_data_is_an_empty_reference_and_limits_are_stated() {
    let mut empty = vec![b' '; 80];
    empty.extend_from_slice(&0u32.to_le_bytes());
    let scan = read(&empty).unwrap();
    assert!(scan.is_empty());
    assert_eq!(scan.fast_distance([0.0, 0.0, 0.0]), f64::INFINITY);
    let bytes = obj(&CUBE, &cube_triangles());
    let limited = |triangles, vertices| {
        ScanReference::read_with_limits(
            &bytes,
            &Rational::one(),
            &Pose::identity(),
            triangles,
            vertices,
            &mut Decider::default(),
        )
    };
    assert_eq!(limited(11, 8).unwrap_err().code, "import.scan-too-large");
    assert_eq!(limited(12, 7).unwrap_err().code, "import.scan-too-large");
    assert_eq!(limited(12, 8).unwrap().report().triangles, 12);
    assert_eq!(
        read(b"f 1 2 3\n").unwrap_err().code,
        "import.mesh-malformed",
        "unreadable bytes keep the reader's code"
    );
}

#[test]
fn reading_charges_the_callers_allowance() {
    let bytes = obj(&CUBE, &cube_triangles());
    let mut dec = Decider::default();
    let before = dec.stats().work;
    ScanReference::read(&bytes, &Rational::one(), &Pose::identity(), &mut dec).unwrap();
    assert!(dec.stats().work > before);
}

#[test]
fn a_recipe_step_binds_the_scan_to_its_hash() {
    let bytes = obj(&CUBE, &cube_triangles());
    let step = |sha256: String| catallax_recipe::Operation::ScanImport {
        file: catallax_recipe::FileRef {
            name: "cube.obj".into(),
            sha256,
        },
        pose: Pose::identity(),
        scale: Number::new(q("2")),
    };
    let mut dec = Decider::default();
    let good = step(catallax_import::sha256_hex(&bytes));
    let scan = ScanReference::from_step(&good, &bytes, &mut dec).unwrap();
    assert_eq!(scan.points()[7], ["4", "4", "4"].map(q));
    let other = step(catallax_import::sha256_hex(b"other bytes"));
    let refusal = ScanReference::from_step(&other, &bytes, &mut dec).unwrap_err();
    assert_eq!(refusal.code, "import.scan-hash");
    let block = catallax_recipe::Operation::Block {
        size: [1, 1, 1].map(Number::from),
        pose: Pose::identity(),
    };
    let refusal = ScanReference::from_step(&block, &bytes, &mut dec).unwrap_err();
    assert_eq!(refusal.code, "import.scan-step");
}

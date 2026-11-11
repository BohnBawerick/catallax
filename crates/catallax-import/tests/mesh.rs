//! Closed mesh reading and exact validity, RECIPE.md section 4.19.

use catallax_import::mesh::{ClosedMesh, MAX_TRIANGLES};
use catallax_numbers::{Decider, Rational, Refusal};

type Tri = [usize; 3];

fn read_mesh(bytes: &[u8], scale: &Rational) -> Result<ClosedMesh, Refusal> {
    ClosedMesh::read(bytes, scale, &mut Decider::default())
}

fn read_limited(bytes: &[u8], scale: &Rational, limit: usize) -> Result<ClosedMesh, Refusal> {
    ClosedMesh::read_with_limit(bytes, scale, limit, &mut Decider::default())
}

fn wind(mesh: &ClosedMesh, p: &[Rational; 3]) -> Option<i64> {
    mesh.winding(p, &mut Decider::default()).unwrap()
}

fn dist(mesh: &ClosedMesh, p: &[Rational; 3]) -> Option<Rational> {
    mesh.squared_distance(p, &mut Decider::default()).unwrap()
}

fn coplanar(mesh: &ClosedMesh, t: u32, slot: usize) -> bool {
    mesh.coplanar_across(t, slot, &mut Decider::default())
        .unwrap()
}

fn q(s: &str) -> Rational {
    Rational::parse(s).unwrap()
}

fn point(p: [&str; 3]) -> [Rational; 3] {
    p.map(q)
}

/// The eight corners of an axis box from `lo` to `hi`, and its twelve outward triangles.
fn cuboid(lo: [i32; 3], hi: [i32; 3], first: usize) -> (Vec<[i32; 3]>, Vec<Tri>) {
    let corners = (0..8)
        .map(|i| {
            [
                if i & 1 == 0 { lo[0] } else { hi[0] },
                if i & 2 == 0 { lo[1] } else { hi[1] },
                if i & 4 == 0 { lo[2] } else { hi[2] },
            ]
        })
        .collect();
    let faces: [[usize; 4]; 6] = [
        [0, 2, 3, 1], // z = lo
        [4, 5, 7, 6], // z = hi
        [0, 1, 5, 4], // y = lo
        [2, 6, 7, 3], // y = hi
        [0, 4, 6, 2], // x = lo
        [1, 3, 7, 5], // x = hi
    ];
    let triangles = faces
        .iter()
        .flat_map(|f| [[f[0], f[1], f[2]], [f[0], f[2], f[3]]])
        .map(|t| t.map(|v| v + first))
        .collect();
    (corners, triangles)
}

fn ascii(points: &[[i32; 3]], triangles: &[Tri]) -> Vec<u8> {
    let mut s = String::from("solid test\n");
    for t in triangles {
        s.push_str("  facet normal 0 0 0\n    outer loop\n");
        for &v in t {
            let [x, y, z] = points[v];
            s.push_str(&format!("      vertex {x} {y} {z}\n"));
        }
        s.push_str("    endloop\n  endfacet\n");
    }
    s.push_str("endsolid test\n");
    s.into_bytes()
}

fn binary(points: &[[f32; 3]], triangles: &[Tri]) -> Vec<u8> {
    let mut bytes = vec![b' '; 80];
    bytes.extend_from_slice(&(triangles.len() as u32).to_le_bytes());
    for t in triangles {
        bytes.extend_from_slice(&[0u8; 12]);
        for &v in t {
            for x in points[v] {
                bytes.extend_from_slice(&x.to_le_bytes());
            }
        }
        bytes.extend_from_slice(&[0u8; 2]);
    }
    bytes
}

fn read(bytes: &[u8]) -> Result<ClosedMesh, String> {
    read_mesh(bytes, &Rational::one()).map_err(|r| r.code.to_string())
}

fn code(bytes: &[u8]) -> String {
    read(bytes).expect_err("the mesh must be refused")
}

fn cube() -> (Vec<[i32; 3]>, Vec<Tri>) {
    cuboid([0, 0, 0], [2, 2, 2], 0)
}

#[test]
fn a_closed_cube_bounds_its_interior_by_winding_number() {
    let (points, triangles) = cube();
    let mesh = read(&ascii(&points, &triangles)).unwrap();
    assert_eq!(mesh.triangles().len(), 12);
    assert_eq!(mesh.vertex_count(), 8);
    assert_eq!(mesh.shells().len(), 1);
    assert!(mesh.shells()[0].outward);
    assert_eq!(mesh.lumps(), vec![(0, vec![])]);
    assert_eq!(wind(&mesh, &point(["1", "1", "1"])), Some(1));
    // Rays from these points pass exactly through edges and vertices.
    assert_eq!(wind(&mesh, &point(["-1", "0", "0"])), Some(0));
    assert_eq!(wind(&mesh, &point(["-1", "2", "2"])), Some(0));
    assert_eq!(wind(&mesh, &point(["1/2", "1", "1"])), Some(1));
    assert_eq!(wind(&mesh, &point(["1", "0", "0"])), None);
    assert_eq!(wind(&mesh, &point(["2", "1", "1/3"])), None);
    assert_eq!(wind(&mesh, &point(["3", "1", "1"])), Some(0));
    // No threshold: a point one part in 10^30 outside is outside.
    assert_eq!(
        wind(
            &mesh,
            &point(["2.000000000000000000000000000001", "1", "1"])
        ),
        Some(0)
    );
    assert_eq!(
        wind(
            &mesh,
            &point(["1.999999999999999999999999999999", "1", "1"])
        ),
        Some(1)
    );
    assert_eq!(dist(&mesh, &point(["1", "1", "1"])), Some(q("1")));
    assert_eq!(dist(&mesh, &point(["3", "3", "2"])), Some(q("2")));
    assert_eq!(
        mesh.bounds().unwrap(),
        [[q("0"), q("2")], [q("0"), q("2")], [q("0"), q("2")]]
    );
}

#[test]
fn binary_ascii_and_obj_files_give_the_same_exact_mesh() {
    let (points, triangles) = cube();
    let floats: Vec<[f32; 3]> = points.iter().map(|p| p.map(|x| x as f32)).collect();
    let mut obj = String::from("# cube\no cube\n");
    for p in &points {
        obj.push_str(&format!("v {}.0 {}.0 {}e0\n", p[0], p[1], p[2]));
    }
    for t in &triangles {
        obj.push_str(&format!(
            "f {}/1/1 {}//1 {}\n",
            t[0] + 1,
            t[1] + 1,
            t[2] + 1
        ));
    }
    let meshes = [
        read(&ascii(&points, &triangles)).unwrap(),
        read(&binary(&floats, &triangles)).unwrap(),
        read(obj.as_bytes()).unwrap(),
    ];
    for mesh in &meshes {
        assert_eq!(mesh.triangles(), meshes[0].triangles());
        for v in 0..8 {
            assert_eq!(mesh.file_vertex(v), meshes[0].file_vertex(v));
        }
    }
}

#[test]
fn scale_multiplies_exact_file_coordinates() {
    let (points, triangles) = cube();
    let mesh = read_mesh(&ascii(&points, &triangles), &q("25.4")).unwrap();
    assert!((0..8).any(|v| mesh.vertex(v) == point(["50.8", "50.8", "50.8"])));
    assert_eq!(mesh.bounds().unwrap()[0], [q("0"), q("50.8")]);
    assert_eq!(wind(&mesh, &point(["50.7", "1", "1"])), Some(1));
    assert_eq!(wind(&mesh, &point(["50.8", "1", "1"])), None);
    assert_eq!(dist(&mesh, &point(["60.8", "1", "1"])), Some(q("100")));
    for scale in ["0", "-1"] {
        let refused = read_mesh(&ascii(&points, &triangles), &q(scale)).unwrap_err();
        assert_eq!(refused.code, "import.mesh-scale");
    }
}

#[test]
fn binary_coordinates_are_the_exact_binary32_values() {
    let (points, triangles) = cube();
    let floats: Vec<[f32; 3]> = points
        .iter()
        .map(|p| p.map(|x| if x == 0 { 0.0 } else { 0.1 }))
        .collect();
    let mesh = read(&binary(&floats, &triangles)).unwrap();
    let tenth = Rational::from_f64_exact(f64::from(0.1f32)).unwrap();
    assert_ne!(tenth, q("0.1"));
    let corner = [tenth.clone(), tenth.clone(), tenth];
    assert!((0..8).any(|v| mesh.file_vertex(v) == &corner));
}

#[test]
fn an_open_mesh_is_impossible_and_names_the_open_edge() {
    let (points, mut triangles) = cube();
    triangles.pop();
    let refused = read_mesh(&ascii(&points, &triangles), &Rational::one()).unwrap_err();
    assert_eq!(refused.code, "import.mesh-open");
    assert!(refused.detail.contains("not closed"), "{}", refused.detail);
}

#[test]
fn a_reversed_triangle_is_an_orientation_error() {
    let (points, mut triangles) = cube();
    triangles[3].swap(1, 2);
    assert_eq!(code(&ascii(&points, &triangles)), "import.mesh-orientation");
}

#[test]
fn an_inside_out_surface_has_winding_minus_one_and_is_refused() {
    let (points, mut triangles) = cube();
    for t in &mut triangles {
        t.swap(1, 2);
    }
    assert_eq!(code(&ascii(&points, &triangles)), "import.mesh-winding");
}

#[test]
fn a_cavity_must_face_into_the_void() {
    let (mut points, mut triangles) = cuboid([0, 0, 0], [6, 6, 6], 0);
    let (inner, mut inner_triangles) = cuboid([2, 2, 2], [4, 4, 4], 8);
    points.extend(inner);
    // Facing outward, the inner box would make winding two inside it.
    let mut wrong = triangles.clone();
    wrong.extend(inner_triangles.clone());
    assert_eq!(code(&ascii(&points, &wrong)), "import.mesh-winding");
    for t in &mut inner_triangles {
        t.swap(1, 2);
    }
    triangles.extend(inner_triangles);
    let mesh = read(&ascii(&points, &triangles)).unwrap();
    assert_eq!(mesh.lumps(), vec![(0, vec![1])]);
    assert!(!mesh.shells()[1].outward);
    assert_eq!(mesh.shells()[1].parent, Some(0));
    assert_eq!(wind(&mesh, &point(["3", "3", "3"])), Some(0));
    assert_eq!(wind(&mesh, &point(["1", "3", "3"])), Some(1));
    assert_eq!(wind(&mesh, &point(["7", "3", "3"])), Some(0));
    assert_eq!(dist(&mesh, &point(["3", "3", "3"])), Some(q("1")));
}

#[test]
fn separate_bodies_and_an_island_inside_a_cavity_are_separate_lumps() {
    let (mut points, mut triangles) = cuboid([0, 0, 0], [10, 10, 10], 0);
    let (cavity, mut cavity_triangles) = cuboid([2, 2, 2], [8, 8, 8], 8);
    for t in &mut cavity_triangles {
        t.swap(1, 2);
    }
    let (island, island_triangles) = cuboid([4, 4, 4], [6, 6, 6], 16);
    let (apart, apart_triangles) = cuboid([20, 0, 0], [21, 1, 1], 24);
    points.extend(cavity.into_iter().chain(island).chain(apart));
    triangles.extend(
        cavity_triangles
            .into_iter()
            .chain(island_triangles)
            .chain(apart_triangles),
    );
    let mesh = read(&ascii(&points, &triangles)).unwrap();
    assert_eq!(mesh.lumps(), vec![(0, vec![1]), (2, vec![]), (3, vec![])]);
    assert_eq!(mesh.shells()[2].parent, Some(1));
    assert_eq!(wind(&mesh, &point(["5", "5", "5"])), Some(1));
    assert_eq!(wind(&mesh, &point(["3", "5", "5"])), Some(0));
}

#[test]
fn overlapping_bodies_cross_and_are_refused() {
    let (mut points, mut triangles) = cuboid([0, 0, 0], [2, 2, 2], 0);
    let (other, other_triangles) = cuboid([1, 1, 1], [3, 3, 3], 8);
    points.extend(other);
    triangles.extend(other_triangles);
    let refused = read_mesh(&ascii(&points, &triangles), &Rational::one()).unwrap_err();
    assert_eq!(refused.code, "import.mesh-crossing");
}

#[test]
fn bodies_touching_along_a_face_cross() {
    let (mut points, mut triangles) = cuboid([0, 0, 0], [2, 2, 2], 0);
    let (other, other_triangles) = cuboid([2, 0, 0], [4, 2, 2], 8);
    points.extend(other);
    triangles.extend(other_triangles);
    // Welding joins the four shared corners, so the shared face's edges have four uses.
    let refused = read(&ascii(&points, &triangles)).unwrap_err();
    assert!(
        refused == "import.mesh-nonmanifold-edge" || refused == "import.mesh-crossing",
        "{refused}"
    );
    let (mut points, mut triangles) = cuboid([0, 0, 0], [2, 2, 2], 0);
    let (other, other_triangles) = cuboid([2, 1, 1], [4, 3, 3], 8);
    points.extend(other);
    triangles.extend(other_triangles);
    assert_eq!(code(&ascii(&points, &triangles)), "import.mesh-crossing");
}

#[test]
fn a_folded_pair_of_triangles_crosses() {
    let points = [[0, 0, 0], [1, 0, 0], [0, 1, 0]];
    assert_eq!(
        code(&ascii(&points, &[[0, 1, 2], [0, 2, 1]])),
        "import.mesh-crossing"
    );
}

#[test]
fn a_surface_folded_through_itself_crosses() {
    // A cube with one corner pushed through the opposite face stays closed and
    // consistently oriented, but its triangles cross.
    let (mut points, triangles) = cube();
    points[7] = [1, 1, -1];
    assert_eq!(code(&ascii(&points, &triangles)), "import.mesh-crossing");
    // Pushed in only part of the way, the surface is valid and nonconvex.
    points[7] = [1, 1, 1];
    let mesh = read(&ascii(&points, &triangles)).unwrap();
    assert_eq!(wind(&mesh, &point(["1/2", "1/2", "1/2"])), Some(1));
    assert_eq!(wind(&mesh, &point(["3/2", "3/2", "3/2"])), Some(0));
}

#[test]
fn zero_area_triangles_have_no_orientation() {
    let (mut points, mut triangles) = cube();
    points.push([1, 0, 0]);
    // Split the bottom-front edge's two triangles through its midpoint,
    // leaving a zero-area sliver.
    triangles.push([0, 8, 1]);
    assert_eq!(code(&ascii(&points, &triangles)), "import.mesh-degenerate");
    let collapsed = [[0, 0, 0], [1, 0, 0], [1, 0, 0]];
    assert_eq!(
        code(&ascii(&collapsed, &[[0, 1, 2]])),
        "import.mesh-degenerate"
    );
}

#[test]
fn shared_edges_and_pinched_vertices_are_unsupported_not_impossible() {
    let (mut points, mut triangles) = cuboid([0, 0, 0], [2, 2, 2], 0);
    let (edge, edge_triangles) = cuboid([2, 2, 0], [4, 4, 2], 8);
    points.extend(edge);
    triangles.extend(edge_triangles);
    assert_eq!(
        code(&ascii(&points, &triangles)),
        "import.mesh-nonmanifold-edge"
    );
    let (mut points, mut triangles) = cuboid([0, 0, 0], [2, 2, 2], 0);
    let (corner, corner_triangles) = cuboid([2, 2, 2], [4, 4, 4], 8);
    points.extend(corner);
    triangles.extend(corner_triangles);
    assert_eq!(
        code(&ascii(&points, &triangles)),
        "import.mesh-pinched-vertex"
    );
}

#[test]
fn a_finely_split_flat_face_is_valid() {
    // Each cube face is split into four triangles around its centre, so many
    // neighbouring triangles are exactly coplanar.
    let (corners, _) = cube();
    let mut points = corners.clone();
    let faces: [[usize; 4]; 6] = [
        [0, 2, 3, 1],
        [4, 5, 7, 6],
        [0, 1, 5, 4],
        [2, 6, 7, 3],
        [0, 4, 6, 2],
        [1, 3, 7, 5],
    ];
    let mut triangles = Vec::new();
    for f in faces {
        let centre: [i32; 3] =
            std::array::from_fn(|i| f.iter().map(|&v| corners[v][i]).sum::<i32>() / 4);
        points.push(centre);
        let c = points.len() - 1;
        for i in 0..4 {
            triangles.push([f[i], f[(i + 1) % 4], c]);
        }
    }
    let mesh = read(&ascii(&points, &triangles)).unwrap();
    assert_eq!(mesh.triangles().len(), 24);
    assert_eq!(wind(&mesh, &point(["1", "1", "1"])), Some(1));
    // The ray from here runs along the split lines of two faces.
    assert_eq!(wind(&mesh, &point(["-5", "1", "1"])), Some(0));
    assert!(
        (0..24).all(|t| (0..3).all(|s| coplanar(&mesh, t, s) == (mesh.twin(t, s) / 4 == t / 4)))
    );
}

#[test]
fn the_size_limit_is_a_typed_refusal() {
    let (points, triangles) = cube();
    let floats: Vec<[f32; 3]> = points.iter().map(|p| p.map(|x| x as f32)).collect();
    let refused = read_limited(&binary(&floats, &triangles), &Rational::one(), 11).unwrap_err();
    assert_eq!(refused.code, "import.mesh-too-large");
    assert!(
        refused.detail.contains("12 triangles"),
        "{}",
        refused.detail
    );
    let refused = read_limited(&ascii(&points, &triangles), &Rational::one(), 11).unwrap_err();
    assert_eq!(refused.code, "import.mesh-too-large");
    assert!(read_limited(&ascii(&points, &triangles), &Rational::one(), 12).is_ok());
    // A binary header announcing more than the stated limit is refused before
    // any triangle is read.
    let mut header = vec![b' '; 80];
    let announced = MAX_TRIANGLES as u32 + 1;
    header.extend_from_slice(&announced.to_le_bytes());
    header.resize(84 + 50 * announced as usize, 0);
    let refused = read_mesh(&header, &Rational::one()).unwrap_err();
    assert_eq!(refused.code, "import.mesh-too-large");
    assert!(
        refused.detail.contains(&MAX_TRIANGLES.to_string()),
        "{}",
        refused.detail
    );
}

#[test]
fn malformed_files_say_why() {
    for bytes in [
        &b"solid x\n facet normal 0 0 0\n outer loop\n vertex 0 0\n"[..],
        b"solid x facet normal 0 0 0 outer loop vertex 0 0 nan",
        b"v 0 0 0\nv 1 0 0\nv 0 1 0\nv 0 0 1\nf 1 2 3 4\n",
        b"v 0 0 0\nf 1 2 3\n",
        b"v 0 0 1e9999\n",
        b"\xff\xfe\x00",
        b"",
    ] {
        assert_eq!(
            code(bytes),
            "import.mesh-malformed",
            "{:?}",
            String::from_utf8_lossy(bytes)
        );
    }
}

#[test]
fn an_empty_closed_mesh_bounds_nothing() {
    let mesh = read(&binary(&[], &[])).unwrap();
    assert!(mesh.triangles().is_empty());
    assert_eq!(wind(&mesh, &point(["0", "0", "0"])), Some(0));
    assert_eq!(mesh.bounds(), None);
    assert_eq!(dist(&mesh, &point(["0", "0", "0"])), None);
}

#[test]
fn validation_and_queries_charge_the_callers_allowance() {
    let (points, triangles) = cube();
    let bytes = ascii(&points, &triangles);
    let (_, used) = catallax_numbers::with_operation_budget(u64::MAX, || {
        read_mesh(&bytes, &Rational::one()).unwrap()
    });
    assert!(used.used > 12, "validation is metered: {used:?}");
    let (refused, work) = catallax_numbers::with_operation_budget(used.used / 2, || {
        read_mesh(&bytes, &Rational::one()).map(|_| ())
    });
    assert!(work.exhausted);
    assert_eq!(refused.unwrap_err().code, "numbers.operation-budget");
    let mesh = read_mesh(&bytes, &Rational::one()).unwrap();
    let centre = point(["1", "1", "1"]);
    let (answer, work) = catallax_numbers::with_operation_budget(1, || {
        mesh.squared_distance(&centre, &mut Decider::default())
    });
    assert!(work.exhausted);
    assert_eq!(answer.unwrap_err().code, "numbers.operation-budget");
    let mut late = Decider::default();
    late.set_deadline(Some(std::time::Instant::now()));
    assert_eq!(
        mesh.winding(&centre, &mut late).unwrap_err().code,
        "numbers.deadline"
    );
}

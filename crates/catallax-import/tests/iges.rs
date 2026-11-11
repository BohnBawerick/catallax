//! IGES 5.3 manifold solids read into the shared exact builder: files from
//! Catallax's own writer and small hand-made files, with typed refusals.

use catallax_boolean::{BoolOp, boolean};
use catallax_iges::{IgesOptions, write_iges};
use catallax_import::{
    Curve, ImportRefusal, Surface, build_file, import_iges_bodies, read_iges, sha256_hex,
};
use catallax_measure::{Accuracy, volume};
use catallax_numbers::{Decider, Rational};
use catallax_recipe::StepId;
use catallax_shapes::{Pose, RVector3, Rotation};
use catallax_topology::{BodyId, Model, block, cone, cylinder};

fn r(n: i64) -> Rational {
    Rational::from_i64(n)
}

fn posed(q: [i64; 4], t: [i64; 3]) -> Pose {
    let [w, x, y, z] = q.map(r);
    let [a, b, c] = t.map(r);
    Pose::from_rational(
        Rotation::from_quaternion(w, x, y, z).unwrap(),
        RVector3::new(a, b, c),
    )
}

fn written(model: &Model, body: BodyId) -> Vec<u8> {
    let mut out = Vec::new();
    write_iges(
        model,
        &[body],
        &IgesOptions::default(),
        &mut Decider::default(),
        &mut out,
    )
    .expect("IGES export");
    out
}

fn measured(model: &Model, body: BodyId) -> (Rational, Rational) {
    let v = volume(model, body, &Accuracy::default(), &mut Decider::default())
        .unwrap()
        .value;
    (v.lo, v.hi)
}

/// Reads `bytes` back into a fresh model and checks its volume overlaps the source's.
fn round_trip(model: &Model, body: BodyId, faces: usize) -> catallax_import::StepFile {
    let bytes = written(model, body);
    let file = read_iges(&bytes).unwrap_or_else(|e| panic!("{e:?}"));
    assert_eq!(file.solids.len(), 1);
    assert_eq!(file.solids[0].faces.len(), faces);
    let mut again = Model::new();
    let bodies = build_file(&mut again, &file, None).unwrap_or_else(|e| panic!("{e:?}"));
    assert_eq!(bodies.len(), 1);
    assert!(again.check(bodies[0]).is_empty());
    let (a, b) = (measured(model, body), measured(&again, bodies[0]));
    assert!(a.0 <= b.1 && b.0 <= a.1, "{a:?} {b:?}");
    file
}

#[test]
fn written_blocks_cylinders_cones_and_bores_read_back() {
    let mut m = Model::new();
    let b = block(&mut m, [2, 3, 5].map(r), &Pose::identity()).unwrap();
    let file = round_trip(&m, b, 6);
    assert!(
        file.solids[0]
            .faces
            .iter()
            .all(|f| matches!(f.surface, Surface::Plane(_)) && f.name.starts_with('D'))
    );

    let mut m = Model::new();
    let c = cylinder(&mut m, r(3), r(7), &Pose::identity()).unwrap();
    let file = round_trip(&m, c, 3);
    let faces = &file.solids[0].faces;
    assert!(
        faces
            .iter()
            .any(|f| matches!(&f.surface, Surface::Cylinder(_, radius) if *radius == r(3)))
    );
    assert!(
        faces
            .iter()
            .flat_map(|f| f.loops.iter().flat_map(|l| l.edges.iter()))
            .any(|e| matches!(e.curve, Curve::Circle(..)))
    );

    // A frustum with slope one is the exact 45 degree cone.
    let mut m = Model::new();
    let k = cone(&mut m, r(3), r(1), r(2), &Pose::identity()).unwrap();
    let file = round_trip(&m, k, 3);
    assert!(
        file.solids[0]
            .faces
            .iter()
            .any(|f| matches!(f.surface, Surface::Cone45(..)))
    );

    // A bored plate: a reversed cylindrical band with its seam used twice.
    let mut m = Model::new();
    let plate = block(&mut m, [4, 4, 2].map(r), &Pose::identity()).unwrap();
    let pose = posed([1, 0, 0, 0], [2, 2, -1]);
    let tool = cylinder(&mut m, r(1), r(4), &pose).unwrap();
    let bored = boolean(
        &mut m,
        plate,
        tool,
        BoolOp::Subtract,
        &mut Decider::default(),
    )
    .unwrap()
    .body;
    round_trip(&m, bored, 7);
}

#[test]
fn rotated_arcs_keep_their_exact_frames() {
    // A quarter turn about x is an exactly orthonormal transform.
    let mut m = Model::new();
    let pose = posed([1, 1, 0, 0], [1, 2, 3]);
    let c = cylinder(&mut m, r(2), r(5), &pose).unwrap();
    round_trip(&m, c, 3);
}

#[test]
fn the_recipe_entry_checks_the_hash_and_builds() {
    let mut m = Model::new();
    let b = block(&mut m, [1, 1, 1].map(r), &Pose::identity()).unwrap();
    let bytes = written(&m, b);
    let step = StepId::new("source").unwrap();
    let mut model = Model::new();
    let bodies = import_iges_bodies(&mut model, &bytes, &sha256_hex(&bytes), &step).unwrap();
    assert_eq!(model.faces_of_body(bodies[0]).len(), 6);
    let mut model = Model::new();
    assert!(matches!(
        import_iges_bodies(&mut model, &bytes, &"0".repeat(64), &step),
        Err(ImportRefusal::HashMismatch { .. })
    ));
}

/// A minimal IGES writer for hand-made fixtures: each entity is its type, form,
/// transformation matrix pointer and parameter text after the type number.
struct Fixture {
    global: String,
    entities: Vec<(u32, u32, usize, String)>,
}

impl Fixture {
    fn new(unit_flag: u32, unit_name: &str) -> Self {
        let name = format!("{}H{unit_name}", unit_name.len());
        Self {
            global: format!(
                "1H,,1H;,4Htest,7Hfixture,4Htest,4Htest,32,38,6,308,15,4Htest,1.0,{unit_flag},{name},1,1.0,15H19700101.000000,1.0E-9,1.0E3,0H,0H,11,0,15H19700101.000000;"
            ),
            entities: Vec::new(),
        }
    }

    /// Adds an entity and returns its directory sequence number.
    fn add(&mut self, kind: u32, form: u32, transform: usize, params: &str) -> usize {
        self.entities.push((kind, form, transform, params.into()));
        2 * self.entities.len() - 1
    }

    fn bytes(&self) -> Vec<u8> {
        let record = |out: &mut String, text: &str, section: char, n: usize| {
            out.push_str(&format!("{text:<72}{section}{n:>7}\n"));
        };
        let mut out = String::new();
        record(&mut out, "hand-made fixture", 'S', 1);
        let global: Vec<char> = self.global.chars().collect();
        let g = global.chunks(72).count();
        for (i, chunk) in global.chunks(72).enumerate() {
            record(&mut out, &chunk.iter().collect::<String>(), 'G', i + 1);
        }
        let mut directory = String::new();
        let mut parameters = String::new();
        let mut line = 1;
        for (i, (kind, form, transform, params)) in self.entities.iter().enumerate() {
            let sequence = 2 * i + 1;
            let text: Vec<char> = format!("{kind},{params};").chars().collect();
            let lines = text.chunks(64).count();
            directory.push_str(&format!(
                "{kind:>8}{line:>8}{:>8}{:>8}{:>8}{:>8}{transform:>8}{:>8}{:>8}D{sequence:>7}\n",
                0, 0, 0, 0, 0, "00000000"
            ));
            directory.push_str(&format!(
                "{kind:>8}{:>8}{:>8}{lines:>8}{form:>8}{:>8}{:>8}{:>8}{:>8}D{:>7}\n",
                0,
                0,
                "",
                "",
                "",
                0,
                sequence + 1
            ));
            for chunk in text.chunks(64) {
                let chunk: String = chunk.iter().collect();
                parameters.push_str(&format!("{chunk:<64} {sequence:>7}P{line:>7}\n"));
                line += 1;
            }
        }
        out.push_str(&directory);
        out.push_str(&parameters);
        let d = 2 * self.entities.len();
        record(
            &mut out,
            &format!("S{:>7}G{g:>7}D{d:>7}P{:>7}", 1, line - 1),
            'T',
            1,
        );
        out.into_bytes()
    }
}

/// A unit cube from planes and lines. `solid` receives the shell pointer and
/// returns the 186 parameters, so a test can change the solid entity.
fn cube(unit_flag: u32, unit_name: &str, solid: impl Fn(usize) -> String) -> Fixture {
    let mut f = Fixture::new(unit_flag, unit_name);
    let corners = [
        [0, 0, 0],
        [1, 0, 0],
        [1, 1, 0],
        [0, 1, 0],
        [0, 0, 1],
        [1, 0, 1],
        [1, 1, 1],
        [0, 1, 1],
    ];
    let coordinates: Vec<String> = corners
        .iter()
        .map(|c| format!("{}.,{}.,{}.", c[0], c[1], c[2]))
        .collect();
    let vertices = f.add(502, 1, 0, &format!("8,{}", coordinates.join(",")));
    // Faces as outward counter-clockwise corner cycles, with their outward normals.
    let faces: [([usize; 4], [i64; 3]); 6] = [
        ([0, 3, 2, 1], [0, 0, -1]),
        ([4, 5, 6, 7], [0, 0, 1]),
        ([0, 1, 5, 4], [0, -1, 0]),
        ([1, 2, 6, 5], [1, 0, 0]),
        ([2, 3, 7, 6], [0, 1, 0]),
        ([3, 0, 4, 7], [-1, 0, 0]),
    ];
    let mut edges: Vec<[usize; 2]> = Vec::new();
    for (cycle, _) in &faces {
        for i in 0..4 {
            let (a, b) = (cycle[i], cycle[(i + 1) % 4]);
            if !edges.contains(&[a, b]) && !edges.contains(&[b, a]) {
                edges.push([a, b]);
            }
        }
    }
    let mut rows = Vec::new();
    for [a, b] in &edges {
        let line = f.add(
            110,
            0,
            0,
            &format!("{},{}", coordinates[*a], coordinates[*b]),
        );
        rows.push(format!("{line},{vertices},{},{vertices},{}", a + 1, b + 1));
    }
    let edge_list = f.add(504, 1, 0, &format!("{},{}", edges.len(), rows.join(",")));
    let mut shell = Vec::new();
    for (cycle, normal) in &faces {
        let corner = &coordinates[cycle[0]];
        let point = f.add(116, 0, 0, corner);
        let direction = f.add(
            123,
            0,
            0,
            &format!("{}.,{}.,{}.", normal[0], normal[1], normal[2]),
        );
        let plane = f.add(190, 0, 0, &format!("{point},{direction}"));
        let mut uses = Vec::new();
        for i in 0..4 {
            let (a, b) = (cycle[i], cycle[(i + 1) % 4]);
            let (index, forward) = match edges.iter().position(|e| *e == [a, b]) {
                Some(index) => (index, 1),
                None => (edges.iter().position(|e| *e == [b, a]).unwrap(), 0),
            };
            uses.push(format!("0,{edge_list},{},{forward},0", index + 1));
        }
        let lp = f.add(508, 1, 0, &format!("4,{}", uses.join(",")));
        let face = f.add(510, 1, 0, &format!("{plane},1,1,{lp}"));
        shell.push(format!("{face},1"));
    }
    let shell = f.add(514, 1, 0, &format!("6,{}", shell.join(",")));
    let params = solid(shell);
    f.add(186, 0, 0, &params);
    f
}

fn ordinary(shell: usize) -> String {
    format!("{shell},1,0")
}

#[test]
fn a_hand_made_inch_cube_is_exactly_one_cubic_inch() {
    let bytes = cube(1, "IN", ordinary).bytes();
    let file = read_iges(&bytes).unwrap_or_else(|e| panic!("{e:?}"));
    assert_eq!(file.mm_per_unit, Rational::from_ratio(127, 5).unwrap());
    let mut m = Model::new();
    let bodies = build_file(&mut m, &file, None).unwrap_or_else(|e| panic!("{e:?}"));
    let inch = Rational::from_ratio(127, 5).unwrap();
    let (lo, hi) = measured(&m, bodies[0]);
    assert_eq!(
        (lo, hi),
        (&(&inch * &inch) * &inch, &(&inch * &inch) * &inch)
    );
    // Face keys are directory sequence numbers.
    assert!(file.solids[0].faces.iter().all(|f| f.name.starts_with('D')));
}

#[test]
fn unsupported_content_is_a_typed_refusal() {
    let refused = |bytes: Vec<u8>| read_iges(&bytes).unwrap_err();
    // Voids.
    let voids = cube(2, "MM", |shell| format!("{shell},1,1,{shell},0"));
    assert!(matches!(
        refused(voids.bytes()),
        ImportRefusal::Unsupported { type_name, .. } if type_name.contains("voids")
    ));
    // An unknown unit name.
    assert!(matches!(
        refused(cube(3, "FURLONG", ordinary).bytes()),
        ImportRefusal::UnknownUnit { .. }
    ));
    // A surface model has no solid.
    let mut surfaces = Fixture::new(2, "MM");
    let point = surfaces.add(116, 0, 0, "0.,0.,0.");
    let normal = surfaces.add(123, 0, 0, "0.,0.,1.");
    let plane = surfaces.add(190, 0, 0, &format!("{point},{normal}"));
    surfaces.add(144, 0, 0, &format!("{plane},1,0,0"));
    assert!(matches!(
        refused(surfaces.bytes()),
        ImportRefusal::Unsupported { type_name, .. } if type_name.contains("144")
    ));
    // A scaling transform on a curve is not a rigid motion.
    let mut scaled = cube(2, "MM", ordinary);
    let matrix = scaled.add(124, 0, 0, "2.,0.,0.,0.,0.,2.,0.,0.,0.,0.,2.,0.");
    let first_line = scaled.entities.iter().position(|e| e.0 == 110).unwrap();
    scaled.entities[first_line].2 = matrix;
    assert!(matches!(
        refused(scaled.bytes()),
        ImportRefusal::Unsupported { type_name, .. } if type_name.contains("orthonormal")
    ));
    // A record cut short of its section letter.
    assert!(matches!(
        refused(b"too short\n".to_vec()),
        ImportRefusal::Malformed { .. }
    ));
}

#[test]
fn an_arc_with_an_irrational_radius_is_refused() {
    let mut f = Fixture::new(2, "MM");
    let arc = f.add(100, 0, 0, "0.,0.,0.,1.,1.,-1.,1.");
    let start = f.add(502, 1, 0, "2,1.,1.,0.,-1.,1.,0.");
    let edges = f.add(504, 1, 0, &format!("1,{arc},{start},1,{start},2"));
    let lp = f.add(508, 1, 0, &format!("1,0,{edges},1,1,0"));
    let point = f.add(116, 0, 0, "0.,0.,0.");
    let normal = f.add(123, 0, 0, "0.,0.,1.");
    let plane = f.add(190, 0, 0, &format!("{point},{normal}"));
    let face = f.add(510, 1, 0, &format!("{plane},1,1,{lp}"));
    let shell = f.add(514, 1, 0, &format!("1,{face},1"));
    f.add(186, 0, 0, &format!("{shell},1,0"));
    assert!(matches!(
        read_iges(&f.bytes()).unwrap_err(),
        ImportRefusal::Unsupported { type_name, .. } if type_name.contains("radius")
    ));
}

/// Files from another writer: OpenCASCADE 8.0 (`cadquery-ocp==8.0.1.1.0`,
/// `IGESControl_Writer("MM", 1)`, B-rep mode), with the author field replaced.
fn occt(name: &str) -> Vec<u8> {
    let path = format!(
        "{}/tests/fixtures/iges/occt-{name}.igs",
        env!("CARGO_MANIFEST_DIR")
    );
    std::fs::read(path).unwrap()
}

#[test]
fn another_writers_box_and_bored_plate_build_exactly() {
    let file = read_iges(&occt("box")).unwrap();
    let mut m = Model::new();
    let bodies = build_file(&mut m, &file, None).unwrap();
    assert_eq!(measured(&m, bodies[0]), (r(30), r(30)));

    let file = read_iges(&occt("bored")).unwrap();
    assert_eq!(file.solids[0].faces.len(), 7);
    let mut m = Model::new();
    let bodies = build_file(&mut m, &file, None).unwrap();
    let mut source = Model::new();
    let plate = block(&mut source, [4, 4, 2].map(r), &Pose::identity()).unwrap();
    let tool = cylinder(&mut source, r(1), r(4), &posed([1, 0, 0, 0], [2, 2, -1])).unwrap();
    let bored = boolean(
        &mut source,
        plate,
        tool,
        BoolOp::Subtract,
        &mut Decider::default(),
    )
    .unwrap()
    .body;
    let (a, b) = (measured(&m, bodies[0]), measured(&source, bored));
    assert!(a.0 <= b.1 && b.0 <= a.1, "{a:?} {b:?}");
}

#[test]
fn another_writers_rounded_trig_values_are_refused_not_moved() {
    // The seam vertex sits at y = -7.3e-16, the binary sine of a full turn:
    // not exactly on the cap circle, so nothing is built until it is repaired.
    let mut m = Model::new();
    assert!(matches!(
        build_file(&mut m, &read_iges(&occt("cylinder")).unwrap(), None),
        Err(ImportRefusal::NotSolid { .. })
    ));
    // The sphere's seam arc is placed by a matrix with 2.4e-16 entries.
    assert!(matches!(
        read_iges(&occt("sphere")).unwrap_err(),
        ImportRefusal::Unsupported { type_name, .. } if type_name.contains("orthonormal")
    ));
}

#[test]
fn a_hand_made_sphere_closes_through_its_seam_and_poles() {
    let mut f = Fixture::new(2, "MM");
    // The seam meridian lies in y = 0, x > 0: a half circle turned a quarter
    // about x by an exactly orthonormal matrix.
    let turn = f.add(124, 0, 0, "1.,0.,0.,0.,0.,0.,-1.,0.,0.,1.,0.,0.");
    let seam = f.add(100, 0, turn, "0.,0.,0.,0.,-2.,0.,2.");
    let poles = f.add(502, 1, 0, "2,0.,0.,-2.,0.,0.,2.");
    let edges = f.add(504, 1, 0, &format!("1,{seam},{poles},1,{poles},2"));
    // Seam up, north pole, seam down, south pole: the poles are degenerate edges.
    let lp = f.add(
        508,
        1,
        0,
        &format!("4,0,{edges},1,1,0,1,{poles},2,1,0,0,{edges},1,0,0,1,{poles},1,1,0"),
    );
    let centre = f.add(116, 0, 0, "0.,0.,0.");
    let axis = f.add(123, 0, 0, "0.,0.,1.");
    let reference = f.add(123, 0, 0, "1.,0.,0.");
    let sphere = f.add(196, 1, 0, &format!("{centre},2.,{axis},{reference}"));
    let face = f.add(510, 1, 0, &format!("{sphere},1,1,{lp}"));
    let shell = f.add(514, 1, 0, &format!("1,{face},1"));
    f.add(186, 0, 0, &format!("{shell},1,0"));
    let file = read_iges(&f.bytes()).unwrap_or_else(|e| panic!("{e:?}"));
    let mut m = Model::new();
    let bodies = build_file(&mut m, &file, None).unwrap_or_else(|e| panic!("{e:?}"));
    let mut source = Model::new();
    let ball = catallax_topology::sphere(&mut source, r(2), &Pose::identity()).unwrap();
    let (a, b) = (measured(&m, bodies[0]), measured(&source, ball));
    assert!(a.0 <= b.1 && b.0 <= a.1, "{a:?} {b:?}");
}

#[test]
fn non_ascii_text_is_malformed_not_a_panic() {
    // The section letter is in column 73 by characters but not by bytes.
    let line = format!("{}é{}", " ".repeat(71), "G0000001\n");
    let mut bytes = cube(2, "MM", ordinary).bytes();
    bytes.extend_from_slice(line.as_bytes());
    assert!(matches!(
        read_iges(&bytes).unwrap_err(),
        ImportRefusal::Malformed { .. }
    ));
}

#[test]
fn huge_entity_counts_are_malformed_before_allocation() {
    let huge = "18446744073709551615";
    let mut shell = Fixture::new(2, "MM");
    let s = shell.add(514, 1, 0, huge);
    shell.add(186, 0, 0, &format!("{s},1,0"));
    let mut face = Fixture::new(2, "MM");
    let point = face.add(116, 0, 0, "0.,0.,0.");
    let normal = face.add(123, 0, 0, "0.,0.,1.");
    let plane = face.add(190, 0, 0, &format!("{point},{normal}"));
    let f = face.add(510, 1, 0, &format!("{plane},{huge},1"));
    let s = face.add(514, 1, 0, &format!("1,{f},1"));
    face.add(186, 0, 0, &format!("{s},1,0"));
    let mut lp = Fixture::new(2, "MM");
    let point = lp.add(116, 0, 0, "0.,0.,0.");
    let normal = lp.add(123, 0, 0, "0.,0.,1.");
    let plane = lp.add(190, 0, 0, &format!("{point},{normal}"));
    let l = lp.add(508, 1, 0, huge);
    let f = lp.add(510, 1, 0, &format!("{plane},1,1,{l}"));
    let s = lp.add(514, 1, 0, &format!("1,{f},1"));
    lp.add(186, 0, 0, &format!("{s},1,0"));
    let mut surface = Fixture::new(2, "MM");
    let spline = surface.add(128, 0, 0, &format!("{huge},{huge},1,1,0,0,1,0,0"));
    let f = surface.add(510, 1, 0, &format!("{spline},0,0"));
    let s = surface.add(514, 1, 0, &format!("1,{f},1"));
    surface.add(186, 0, 0, &format!("{s},1,0"));
    for (name, fixture) in [
        ("shell", shell),
        ("face", face),
        ("loop", lp),
        ("surface", surface),
    ] {
        assert!(
            matches!(
                read_iges(&fixture.bytes()),
                Err(ImportRefusal::Malformed { .. })
            ),
            "{name}"
        );
    }
}

//! Native faces whose coarse display needs its trim rebuilt from the source
//! geometry. Each fixture is a small public shape of the kind the display
//! used to leave out, with the reason in its test's comment.
use catallax_import::{sha256_hex, tolerant};
use catallax_mesh::{PreviewMesh, PreviewOptions, Triangle, preview_tolerant_body};
use catallax_numbers::{Decider, Rational, with_operation_budget};

type P = [f64; 3];

/// STEP text for a few entities, numbered in the order they are added.
#[derive(Default)]
struct Step {
    entities: Vec<String>,
}

fn real(value: f64) -> String {
    assert!(value.is_finite());
    let text = format!("{value}");
    if text.contains('.') {
        text
    } else {
        format!("{text}.")
    }
}

fn coordinates(p: P) -> String {
    format!("({},{},{})", real(p[0]), real(p[1]), real(p[2]))
}

fn references(ids: &[usize]) -> String {
    ids.iter()
        .map(|id| format!("#{id}"))
        .collect::<Vec<_>>()
        .join(",")
}

impl Step {
    fn add(&mut self, entity: impl Into<String>) -> usize {
        self.entities.push(entity.into());
        self.entities.len()
    }
    fn point(&mut self, p: P) -> usize {
        self.add(format!("CARTESIAN_POINT('',{})", coordinates(p)))
    }
    fn placement(&mut self, origin: P, axis: P, reference: P) -> usize {
        let origin = self.point(origin);
        let axis = self.add(format!("DIRECTION('',{})", coordinates(axis)));
        let reference = self.add(format!("DIRECTION('',{})", coordinates(reference)));
        self.add(format!(
            "AXIS2_PLACEMENT_3D('',#{origin},#{axis},#{reference})"
        ))
    }
    fn plane(&mut self, origin: P, normal: P) -> usize {
        let reference = if normal[0].abs() < 0.5 {
            [1.0, 0.0, 0.0]
        } else {
            [0.0, 0.0, 1.0]
        };
        let placement = self.placement(origin, normal, reference);
        self.add(format!("PLANE('',#{placement})"))
    }
    fn vertex(&mut self, p: P) -> usize {
        let point = self.point(p);
        self.add(format!("VERTEX_POINT('',#{point})"))
    }
    fn edge(&mut self, start: usize, end: usize, curve: usize) -> usize {
        self.add(format!("EDGE_CURVE('',#{start},#{end},#{curve},.T.)"))
    }
    fn line(&mut self, a: P, b: P) -> usize {
        let origin = self.point(a);
        let d: P = std::array::from_fn(|i| b[i] - a[i]);
        let length = d.iter().fold(0.0f64, |n, v| n.hypot(*v));
        let direction = self.add(format!(
            "DIRECTION('',{})",
            coordinates(d.map(|v| v / length))
        ));
        let vector = self.add(format!("VECTOR('',#{direction},{})", real(length)));
        self.add(format!("LINE('',#{origin},#{vector})"))
    }
    fn circle(&mut self, center: P, normal: P, radius: f64) -> usize {
        let placement = self.placement(center, normal, [1.0, 0.0, 0.0]);
        self.add(format!("CIRCLE('',#{placement},{})", real(radius)))
    }
    fn bound(&mut self, uses: &[(usize, bool)], outer: bool) -> usize {
        let oriented: Vec<_> = uses
            .iter()
            .map(|&(edge, forward)| {
                self.add(format!(
                    "ORIENTED_EDGE('',*,*,#{edge},{})",
                    if forward { ".T." } else { ".F." }
                ))
            })
            .collect();
        let edge_loop = self.add(format!("EDGE_LOOP('',({}))", references(&oriented)));
        self.add(format!(
            "{}('',#{edge_loop},.T.)",
            if outer {
                "FACE_OUTER_BOUND"
            } else {
                "FACE_BOUND"
            }
        ))
    }
    fn vertex_bound(&mut self, vertex: usize) -> usize {
        let vertex_loop = self.add(format!("VERTEX_LOOP('',#{vertex})"));
        self.add(format!("FACE_BOUND('',#{vertex_loop},.T.)"))
    }
    fn face(&mut self, bounds: &[usize], surface: usize, same_sense: bool) -> usize {
        self.add(format!(
            "ADVANCED_FACE('',({}),#{surface},{})",
            references(bounds),
            if same_sense { ".T." } else { ".F." }
        ))
    }
    fn solid(mut self, faces: &[usize], metre: bool) -> String {
        let shell = self.add(format!("CLOSED_SHELL('',({}))", references(faces)));
        let root = self.add(format!("MANIFOLD_SOLID_BREP('',#{shell})"));
        self.document(root, "ADVANCED_BREP_SHAPE_REPRESENTATION", metre)
    }
    fn surfaces(mut self, faces: &[usize], closed: bool) -> String {
        let shell = self.add(format!(
            "{}('',({}))",
            if closed { "CLOSED_SHELL" } else { "OPEN_SHELL" },
            references(faces)
        ));
        let root = self.add(format!("SHELL_BASED_SURFACE_MODEL('',(#{shell}))"));
        self.document(root, "MANIFOLD_SURFACE_SHAPE_REPRESENTATION", false)
    }
    fn document(mut self, root: usize, representation: &str, metre: bool) -> String {
        let unit = self.add(if metre {
            "( LENGTH_UNIT() NAMED_UNIT(*) SI_UNIT($,.METRE.) )"
        } else {
            "( LENGTH_UNIT() NAMED_UNIT(*) SI_UNIT(.MILLI.,.METRE.) )"
        });
        let context = self.add(format!(
            "( GEOMETRIC_REPRESENTATION_CONTEXT(3) GLOBAL_UNIT_ASSIGNED_CONTEXT((#{unit})) REPRESENTATION_CONTEXT('','') )"
        ));
        self.add(format!("{representation}('',(#{root}),#{context})"));
        let mut text = String::from(
            "ISO-10303-21;\nHEADER;\nFILE_DESCRIPTION(('native trim'),'2;1');\nFILE_NAME('trim.step','',(''),(''),'','','');\nFILE_SCHEMA(('AUTOMOTIVE_DESIGN'));\nENDSEC;\nDATA;\n",
        );
        for (index, entity) in self.entities.iter().enumerate() {
            text.push_str(&format!("#{}={entity};\n", index + 1));
        }
        text.push_str("ENDSEC;\nEND-ISO-10303-21;\n");
        text
    }
}

fn import(source: &str) -> tolerant::admission::TolerantImport {
    let bytes = source.as_bytes();
    let (import, work) = with_operation_budget(20_000_000, || {
        tolerant::admission::import_step(
            bytes,
            &sha256_hex(bytes),
            Default::default(),
            &mut Decider::default(),
        )
    });
    assert!(!work.exhausted);
    import.unwrap()
}

fn preview(
    import: &tolerant::admission::TolerantImport,
    options: &PreviewOptions,
) -> PreviewMesh<String, String> {
    let (preview, work) = with_operation_budget(100_000_000, || {
        preview_tolerant_body(import, 0, options, &mut Decider::default())
    });
    assert!(!work.exhausted);
    preview.unwrap()
}

/// Opens a file as the census does, and previews its one native body.
fn open_and_preview(source: &str, options: &PreviewOptions) -> PreviewMesh<String, String> {
    use catallax_import::tolerant::opened::{BodyView, open_step};
    let bytes = source.as_bytes();
    let (opened, work) = with_operation_budget(20_000_000, || {
        open_step(
            bytes,
            &sha256_hex(bytes),
            &catallax_recipe::StepId::new("import").unwrap(),
            Default::default(),
            &mut Decider::default(),
        )
    });
    assert!(!work.exhausted);
    let opened = opened.unwrap();
    let Some(BodyView::Tolerant { import, body }) = opened.body(0) else {
        panic!("the shape opens as a native body");
    };
    let (preview, work) = with_operation_budget(100_000_000, || {
        preview_tolerant_body(import, body, options, &mut Decider::default())
    });
    assert!(!work.exhausted);
    preview.unwrap()
}

fn corners(mesh: &PreviewMesh<String, String>, triangle: &Triangle<String>) -> [P; 3] {
    triangle.vertices.map(|index| mesh.positions[index])
}

fn area(mesh: &PreviewMesh<String, String>, face: &str) -> f64 {
    let face = mesh.faces.iter().find(|f| f.face == face).unwrap();
    mesh.triangles[face.triangles.clone()]
        .iter()
        .map(|triangle| {
            let [a, b, c] = corners(mesh, triangle);
            let u: P = std::array::from_fn(|i| b[i] - a[i]);
            let v: P = std::array::from_fn(|i| c[i] - a[i]);
            let n = [
                u[1] * v[2] - u[2] * v[1],
                u[2] * v[0] - u[0] * v[2],
                u[0] * v[1] - u[1] * v[0],
            ];
            n[0].hypot(n[1]).hypot(n[2]) * 0.5
        })
        .sum()
}

fn volume(mesh: &PreviewMesh<String, String>) -> f64 {
    mesh.triangles
        .iter()
        .map(|triangle| {
            let [a, b, c] = corners(mesh, triangle);
            (a[0] * (b[1] * c[2] - b[2] * c[1])
                + a[1] * (b[2] * c[0] - b[0] * c[2])
                + a[2] * (b[0] * c[1] - b[1] * c[0]))
                / 6.0
        })
        .sum()
}

/// A planar face bounded by two straight degree-one splines between the same
/// two vertices, 1.5 um apart. Admission welds them into one edge, so the
/// face's boundary runs out along it and back: it has no area at all.
fn welded_digon(step: &mut Step, at: P) -> usize {
    let plane = step.plane(at, [0.0, 0.0, 1.0]);
    let a = at;
    let b = [at[0] + 0.2, at[1] + 0.01, at[2]];
    let [start, end] = [step.vertex(a), step.vertex(b)];
    let there = [step.point(a), step.point(b)];
    let out = step.add(format!(
        "B_SPLINE_CURVE_WITH_KNOTS('',1,(#{},#{}),.UNSPECIFIED.,.F.,.F.,(2,2),(0.,1.),.UNSPECIFIED.)",
        there[0], there[1]
    ));
    let back = [
        step.point([b[0], b[1] + 0.0015, b[2]]),
        step.point([a[0], a[1] + 0.0015, a[2]]),
    ];
    let back = step.add(format!(
        "B_SPLINE_CURVE_WITH_KNOTS('',1,(#{},#{}),.UNSPECIFIED.,.F.,.F.,(2,2),(0.,1.),.UNSPECIFIED.)",
        back[0], back[1]
    ));
    let out = step.edge(start, end, out);
    let back = step.edge(end, start, back);
    let bound = step.bound(&[(out, true), (back, true)], false);
    step.add(format!("FACE_SURFACE('',(#{bound}),#{plane},.T.)"))
}

/// The face is left out of the body's face set before any sampling, so the
/// preview's first pass draws the body. It used to fail the face there and
/// at five coarse grids before leaving it out.
#[test]
fn a_welded_digon_with_no_area_draws_only_its_edge() {
    let mut step = Step::default();
    let digon = welded_digon(&mut step, [1.0, 2.0, 0.0]);
    let import = import(&step.surfaces(&[digon], true));
    let (view, work) = with_operation_budget(100_000_000, || {
        preview_tolerant_body(
            &import,
            0,
            &PreviewOptions::default(),
            &mut Decider::default(),
        )
    });
    let view = view.unwrap();
    eprintln!("digon preview work {}", work.used);
    assert!(work.used < 1_000, "preview work {}", work.used);
    assert!(view.faces.is_empty() && view.triangles.is_empty());
    let [edge] = view.edges.as_slice() else {
        panic!("one welded edge: {:?}", view.edges);
    };
    let length: f64 = edge
        .vertices
        .windows(2)
        .map(|pair| {
            let [a, b] = [view.positions[pair[0]], view.positions[pair[1]]];
            (a[0] - b[0]).hypot(a[1] - b[1]).hypot(a[2] - b[2])
        })
        .sum();
    assert!(
        (length - 0.2f64.hypot(0.01)).abs() < 1e-6,
        "length {length}"
    );
}

/// Control points of a closed curve: the corners of a regular pentagon.
fn pentagon(radius: f64, y: f64) -> [P; 5] {
    std::array::from_fn(|k| {
        let angle = std::f64::consts::TAU * k as f64 / 5.0;
        [radius * angle.cos(), y, radius * angle.sin()]
    })
}

/// Dense points of a uniform periodic cubic B-spline in the XZ plane.
fn periodic_cubic(controls: &[P; 5]) -> Vec<[f64; 2]> {
    let mut points = Vec::new();
    for span in 0..5 {
        let q = |k: usize| controls[(span + k) % 5];
        for step in 0..400 {
            let t = f64::from(step) / 400.0;
            let w = [
                (1.0 - t).powi(3),
                3.0 * t.powi(3) - 6.0 * t * t + 4.0,
                -3.0 * t.powi(3) + 3.0 * t * t + 3.0 * t + 1.0,
                t.powi(3),
            ];
            points.push([0, 2].map(|axis| (0..4).map(|k| w[k] * q(k)[axis]).sum::<f64>() / 6.0));
        }
    }
    points
}

/// The area inside a uniform periodic cubic B-spline in the XZ plane.
fn periodic_cubic_area(controls: &[P; 5]) -> f64 {
    let points = periodic_cubic(controls);
    (0..points.len())
        .map(|i| {
            let [a, b] = [points[i], points[(i + 1) % points.len()]];
            a[0] * b[1] - b[0] * a[1]
        })
        .sum::<f64>()
        .abs()
        * 0.5
}

/// The length of a uniform periodic cubic B-spline in the XZ plane.
fn periodic_cubic_length(controls: &[P; 5]) -> f64 {
    let points = periodic_cubic(controls);
    (0..points.len())
        .map(|i| {
            let [a, b] = [points[i], points[(i + 1) % points.len()]];
            (a[0] - b[0]).hypot(a[1] - b[1])
        })
        .sum()
}

/// A tube whose wall is one periodic cubic B-spline, closed along U with the
/// last three control rows repeating the first three, between two closed
/// spline rims. Its trims meet the parameter seam at both rims, and read as
/// open lines the wall has no inside to draw.
#[test]
fn a_closed_periodic_spline_band_draws_between_its_rims() {
    let height = 8.0;
    let mut step = Step::default();
    let knots = "(2,1,1,1,1,1,1,1,1,2)";
    let values = "(-0.4,-0.2,0.,0.2,0.4,0.6,0.8,1.,1.2,1.4)";
    let rows = [pentagon(10.0, 0.0), pentagon(10.0, height)];
    let controls: Vec<[usize; 2]> = (0..8)
        .map(|i| [step.point(rows[0][i % 5]), step.point(rows[1][i % 5])])
        .collect();
    let net = controls
        .iter()
        .map(|[a, b]| format!("(#{a},#{b})"))
        .collect::<Vec<_>>()
        .join(",");
    let wall = step.add(format!(
        "B_SPLINE_SURFACE_WITH_KNOTS('',3,1,({net}),.UNSPECIFIED.,.T.,.F.,.F.,{knots},(2,2),{values},(0.,1.),.UNSPECIFIED.)"
    ));
    let mut rims = Vec::new();
    for row in &rows {
        let points: Vec<usize> = (0..8).map(|i| step.point(row[i % 5])).collect();
        let curve = step.add(format!(
            "B_SPLINE_CURVE_WITH_KNOTS('',3,({}),.UNSPECIFIED.,.T.,.F.,{knots},{values},.UNSPECIFIED.)",
            references(&points)
        ));
        let start: P = std::array::from_fn(|i| (row[0][i] + 4.0 * row[1][i] + row[2][i]) / 6.0);
        let vertex = step.vertex(start);
        let edge = step.edge(vertex, vertex, curve);
        rims.push(edge);
    }
    let lower = step.bound(&[(rims[0], true)], true);
    let upper = step.bound(&[(rims[1], false)], false);
    let band = step.face(&[lower, upper], wall, false);
    let near_plane = step.plane([0.0; 3], [0.0, -1.0, 0.0]);
    let near_bound = step.bound(&[(rims[0], false)], true);
    let near = step.face(&[near_bound], near_plane, true);
    let far_plane = step.plane([0.0, height, 0.0], [0.0, 1.0, 0.0]);
    let far_bound = step.bound(&[(rims[1], true)], true);
    let far = step.face(&[far_bound], far_plane, true);
    let view = open_and_preview(
        &step.solid(&[band, near, far], false),
        &PreviewOptions::default(),
    );
    assert_eq!(view.faces.len(), 3, "{:?}", view.faces);
    let band_name = format!("#{band}");
    assert!(view.faces.iter().any(|face| face.face == band_name));
    // Every drawn point of the wall lies between the rims and inside the
    // control pentagon's circumcircle.
    let wall_face = view
        .faces
        .iter()
        .find(|face| face.face == band_name)
        .unwrap();
    for triangle in &view.triangles[wall_face.triangles.clone()] {
        for p in corners(&view, triangle) {
            assert!((-1e-9..=height + 1e-9).contains(&p[1]));
            assert!(p[0].hypot(p[2]) <= 10.0 + 1e-9);
        }
    }
    // Its chords lie inside the convex wall, and no fold adds area.
    let lateral = periodic_cubic_length(&rows[0]) * height;
    let drawn = area(&view, &band_name);
    assert!(
        drawn < lateral && drawn > lateral * 0.95,
        "wall area {drawn}, lateral area {lateral}"
    );
    let expected = periodic_cubic_area(&rows[0]) * height;
    let volume = volume(&view);
    assert!(
        (volume - expected).abs() < expected * 0.03,
        "volume {volume}, expected {expected}"
    );
}

/// A hemisphere written as one rational spline, closed along V, from its rim
/// at U = 0 to a pole at U = 1 that the face declares with a VERTEX_LOOP.
#[test]
fn a_closed_spline_cap_closes_at_its_declared_vertex() {
    let radius = 10.0;
    let s = std::f64::consts::FRAC_1_SQRT_2;
    let mut step = Step::default();
    // Meridian quarter circle as (distance from the axis, height, weight).
    let profile = [(radius, 0.0, 1.0), (radius, radius, s), (0.0, radius, 1.0)];
    let around: [(f64, f64, f64); 9] = [
        (1.0, 0.0, 1.0),
        (1.0, 1.0, s),
        (0.0, 1.0, 1.0),
        (-1.0, 1.0, s),
        (-1.0, 0.0, 1.0),
        (-1.0, -1.0, s),
        (0.0, -1.0, 1.0),
        (1.0, -1.0, s),
        (1.0, 0.0, 1.0),
    ];
    let mut rows = Vec::new();
    let mut weights = Vec::new();
    for &(rho, z, w) in &profile {
        let ids: Vec<usize> = around
            .iter()
            .map(|&(x, y, _)| step.point([rho * x, rho * y, z]))
            .collect();
        rows.push(format!("({})", references(&ids)));
        weights.push(format!(
            "({})",
            around
                .iter()
                .map(|&(_, _, v)| real(w * v))
                .collect::<Vec<_>>()
                .join(",")
        ));
    }
    let dome = step.add(format!(
        "( BOUNDED_SURFACE() B_SPLINE_SURFACE(2,2,({}),.UNSPECIFIED.,.F.,.T.,.F.) B_SPLINE_SURFACE_WITH_KNOTS((3,3),(3,2,2,2,3),(0.,1.),(0.,0.25,0.5,0.75,1.),.UNSPECIFIED.) GEOMETRIC_REPRESENTATION_ITEM() RATIONAL_B_SPLINE_SURFACE(({})) REPRESENTATION_ITEM('') SURFACE() )",
        rows.join(","),
        weights.join(",")
    ));
    let circle = step.circle([0.0; 3], [0.0, 0.0, 1.0], radius);
    let start = step.vertex([radius, 0.0, 0.0]);
    let rim = step.edge(start, start, circle);
    let pole = step.vertex([0.0, 0.0, radius]);
    let outer = step.bound(&[(rim, true)], true);
    let apex = step.vertex_bound(pole);
    let cap = step.face(&[outer, apex], dome, false);
    let plane = step.plane([0.0; 3], [0.0, 0.0, -1.0]);
    let disk_bound = step.bound(&[(rim, false)], true);
    let disk = step.face(&[disk_bound], plane, true);
    let view = open_and_preview(&step.solid(&[cap, disk], false), &PreviewOptions::default());
    assert_eq!(view.faces.len(), 2, "{:?}", view.faces);
    let cap_name = format!("#{cap}");
    let cap_face = view
        .faces
        .iter()
        .find(|face| face.face == cap_name)
        .unwrap();
    let mut reaches_pole = false;
    for triangle in &view.triangles[cap_face.triangles.clone()] {
        for p in corners(&view, triangle) {
            let r = p[0].hypot(p[1]).hypot(p[2]);
            assert!((r - radius).abs() < 1e-6, "radius {r}");
            assert!(p[2] >= -1e-9);
            reaches_pole |= (p[2] - radius).abs() < 1e-9;
        }
    }
    assert!(reaches_pole);
    // Every drawn point is on the sphere, so the drawing sits inside it.
    let hemisphere = std::f64::consts::TAU / 3.0 * radius.powi(3);
    let volume = volume(&view);
    assert!(
        volume < hemisphere && volume > hemisphere * 0.9,
        "volume {volume}, hemisphere {hemisphere}"
    );
    let surface = std::f64::consts::TAU * radius * radius;
    let drawn = area(&view, &cap_name);
    assert!(drawn < surface && drawn > surface * 0.9, "area {drawn}");
}

/// A cone in a metre file, 89.2 degrees from its axis with a 1 m radius at
/// its placement, cut 3 mm from the axis by a planar disk and closed at a
/// VERTEX_LOOP written to a nanometre near the origin. The placement is
/// 14 mm above the apex, and a meridian sampled one source unit away puts the
/// sampled pole farther off than the rounding of the small coordinates at
/// the apex allows.
#[test]
fn a_flat_metre_cone_closes_at_its_declared_vertex() {
    let (big, angle, rim) = (1.0f64, 1.5567785984f64, 0.003f64);
    let slope = angle.tan();
    let lift = (big / slope * 1e9).round() / 1e9;
    let rim_z = lift + (rim - big) / slope;
    let apex_z = ((lift - big / slope) * 1e9).round() / 1e9;
    let mut step = Step::default();
    let placement = step.placement([0.0, 0.0, lift], [0.0, 0.0, 1.0], [1.0, 0.0, 0.0]);
    let cone = step.add(format!(
        "CONICAL_SURFACE('',#{placement},{},{})",
        real(big),
        real(angle)
    ));
    let circle = step.circle([0.0, 0.0, rim_z], [0.0, 0.0, 1.0], rim);
    let start = step.vertex([rim, 0.0, rim_z]);
    let edge = step.edge(start, start, circle);
    let pole = step.vertex([0.0, 0.0, apex_z]);
    let outer = step.bound(&[(edge, true)], true);
    let apex = step.vertex_bound(pole);
    let lateral = step.face(&[outer, apex], cone, true);
    let plane = step.plane([0.0, 0.0, rim_z], [0.0, 0.0, 1.0]);
    let disk_bound = step.bound(&[(edge, false)], true);
    let disk = step.face(&[disk_bound], plane, true);
    let import = import(&step.solid(&[lateral, disk], true));
    let view = preview(&import, &PreviewOptions::default());
    assert_eq!(view.faces.len(), 2, "{:?}", view.faces);
    let name = format!("#{lateral}");
    let face = view.faces.iter().find(|face| face.face == name).unwrap();
    let apex_mm = apex_z * 1000.0;
    let mut reaches_apex = false;
    for triangle in &view.triangles[face.triangles.clone()] {
        for p in corners(&view, triangle) {
            assert!(p[0].hypot(p[1]) <= rim * 1000.0 + 1e-9);
            assert!((apex_mm - 1e-6..=rim_z * 1000.0 + 1e-9).contains(&p[2]));
            reaches_apex |= p[0].hypot(p[1]) < 1e-9 && (p[2] - apex_mm).abs() < 1e-9;
        }
    }
    assert!(reaches_apex);
    let cone_volume = std::f64::consts::PI * (rim * 1000.0).powi(2) * (rim / slope * 1000.0) / 3.0;
    let volume = volume(&view);
    assert!(
        (volume - cone_volume).abs() < cone_volume * 0.05,
        "volume {volume}, cone {cone_volume}"
    );
}

/// A plate's circular rim 12.5 mm out, round a plus-shaped hole whose corners
/// stand 0.04 mm inside it. At a coarse display the rim's chords cut the
/// hole's corners although the circle clears them all. Both sides of a
/// two-sided lamina share the rim and the hole.
#[test]
fn a_round_rim_close_to_a_straight_hole_draws_the_plate() {
    let (radius, reach, half) = (12.5, 12.45, 0.5);
    let mut step = Step::default();
    let plane = step.plane([0.0; 3], [0.0, 0.0, 1.0]);
    let circle = step.circle([0.0; 3], [0.0, 0.0, 1.0], radius);
    let start = step.vertex([-radius, 0.0, 0.0]);
    let rim = step.edge(start, start, circle);
    let corners_2d = [
        [-reach, half],
        [-half, half],
        [-half, reach],
        [half, reach],
        [half, half],
        [reach, half],
        [reach, -half],
        [half, -half],
        [half, -reach],
        [-half, -reach],
        [-half, -half],
        [-reach, -half],
    ];
    let hole_vertices: Vec<usize> = corners_2d
        .iter()
        .map(|[x, y]| step.vertex([*x, *y, 0.0]))
        .collect();
    let mut sides = Vec::new();
    for i in 0..corners_2d.len() {
        let j = (i + 1) % corners_2d.len();
        let [a, b] = [corners_2d[i], corners_2d[j]];
        let line = step.line([a[0], a[1], 0.0], [b[0], b[1], 0.0]);
        sides.push((step.edge(hole_vertices[i], hole_vertices[j], line), true));
    }
    let outer = step.bound(&[(rim, true)], true);
    let hole = step.bound(&sides, false);
    let plate = step.face(&[outer, hole], plane, true);
    // Its underside on the same plane closes the shell.
    let back: Vec<(usize, bool)> = sides.iter().rev().map(|&(e, f)| (e, !f)).collect();
    let under_outer = step.bound(&[(rim, false)], true);
    let under_hole = step.bound(&back, false);
    let under = step.face(&[under_outer, under_hole], plane, false);
    let import = import(&step.surfaces(&[plate, under], true));
    // A coarse display: chords of the rim at its default sampling sag 0.24 mm.
    let options = PreviewOptions {
        chord_tolerance: Rational::parse("0.5").unwrap(),
        ..PreviewOptions::default()
    };
    let view = preview(&import, &options);
    assert_eq!(view.faces.len(), 2, "{:?}", view.faces);
    let rim_points: std::collections::BTreeSet<usize> = view
        .edges
        .iter()
        .find(|edge| edge.edge == format!("#{rim}"))
        .unwrap()
        .vertices
        .iter()
        .copied()
        .collect();
    for triangle in &view.triangles {
        let points = corners(&view, triangle);
        for (&index, p) in triangle.vertices.iter().zip(points) {
            let r = p[0].hypot(p[1]);
            assert!(r <= radius + 1e-9);
            // The plate meets the rim only at the rim's own shared samples.
            if r > reach + half {
                assert!(rim_points.contains(&index), "{p:?}");
            }
        }
        let center: P = std::array::from_fn(|i| points.iter().map(|p| p[i]).sum::<f64>() / 3.0);
        let in_hole = (center[0].abs() < reach && center[1].abs() < half)
            || (center[0].abs() < half && center[1].abs() < reach);
        assert!(!in_hole, "triangle inside the hole at {center:?}");
    }
    let hole_area = 2.0 * (2.0 * reach) * (2.0 * half) - (2.0 * half).powi(2);
    let disk = std::f64::consts::PI * radius * radius;
    // Inside the circle, outside the hole, and no worse than the 16-gon.
    let sixteen = 8.0 * radius * radius * (std::f64::consts::TAU / 16.0).sin();
    for side in [plate, under] {
        let drawn = area(&view, &format!("#{side}"));
        assert!(
            drawn <= disk - hole_area + 1e-6 && drawn >= sixteen - hole_area - 1e-6,
            "plate area {drawn}"
        );
    }
}

/// The same closed cap as a polynomial spline: square-cornered rows from a
/// rim, which the face's first row also writes as its edge, to one pole
/// point. Read without its seam the rim's trim had no width in the chart.
#[test]
fn a_closed_polynomial_spline_cap_closes_at_its_declared_vertex() {
    let radius = 10.0;
    let around: [(f64, f64); 9] = [
        (1.0, 0.0),
        (1.0, 1.0),
        (0.0, 1.0),
        (-1.0, 1.0),
        (-1.0, 0.0),
        (-1.0, -1.0),
        (0.0, -1.0),
        (1.0, -1.0),
        (1.0, 0.0),
    ];
    let mut step = Step::default();
    let mut rows = Vec::new();
    let mut rim_controls = Vec::new();
    for (row, (rho, z)) in [(radius, 0.0), (radius, radius), (0.0, radius)]
        .into_iter()
        .enumerate()
    {
        let ids: Vec<usize> = around
            .iter()
            .map(|&(x, y)| step.point([rho * x, rho * y, z]))
            .collect();
        if row == 0 {
            rim_controls = ids.clone();
        }
        rows.push(format!("({})", references(&ids)));
    }
    let dome = step.add(format!(
        "B_SPLINE_SURFACE_WITH_KNOTS('',2,2,({}),.UNSPECIFIED.,.F.,.T.,.F.,(3,3),(3,2,2,2,3),(0.,1.),(0.,0.25,0.5,0.75,1.),.UNSPECIFIED.)",
        rows.join(",")
    ));
    let curve = step.add(format!(
        "B_SPLINE_CURVE_WITH_KNOTS('',2,({}),.UNSPECIFIED.,.T.,.F.,(3,2,2,2,3),(0.,0.25,0.5,0.75,1.),.UNSPECIFIED.)",
        references(&rim_controls)
    ));
    let start = step.vertex([radius, 0.0, 0.0]);
    let rim = step.edge(start, start, curve);
    let pole = step.vertex([0.0, 0.0, radius]);
    let outer = step.bound(&[(rim, true)], true);
    let apex = step.vertex_bound(pole);
    let cap = step.face(&[outer, apex], dome, false);
    let plane = step.plane([0.0; 3], [0.0, 0.0, -1.0]);
    let disk_bound = step.bound(&[(rim, false)], true);
    let disk = step.face(&[disk_bound], plane, true);
    let view = open_and_preview(&step.solid(&[cap, disk], false), &PreviewOptions::default());
    assert_eq!(view.faces.len(), 2, "{:?}", view.faces);
    let cap_name = format!("#{cap}");
    let cap_face = view
        .faces
        .iter()
        .find(|face| face.face == cap_name)
        .unwrap();
    let mut reaches_pole = false;
    for triangle in &view.triangles[cap_face.triangles.clone()] {
        for p in corners(&view, triangle) {
            assert!((-1e-9..=radius + 1e-9).contains(&p[2]));
            assert!(p[0].abs().max(p[1].abs()) <= radius + 1e-9);
            reaches_pole |= p[0].hypot(p[1]) < 1e-9 && (p[2] - radius).abs() < 1e-9;
        }
    }
    assert!(reaches_pole);
    // A closed dome over its own base: more than nothing, less than the box
    // round its control rows.
    let volume = volume(&view);
    let disk_area = area(&view, &format!("#{disk}"));
    assert!(
        volume > disk_area * radius / 4.0 && volume < disk_area * radius,
        "volume {volume} over a base of {disk_area}"
    );
}

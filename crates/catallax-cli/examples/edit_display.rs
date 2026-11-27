//! Verify a released edit's live analytic mesh against its exact definition.
use catallax_mesh::{MeshOptions, PreviewMesh, PreviewOptions, preview_mesh_body};
use catallax_numbers::{Decider, Interval, Rational, Real, with_operation_budget};
use catallax_shapes::Surface;
use serde::Deserialize;
use serde_json::json;
use std::{collections::BTreeMap, env, fs, path::Path};

#[derive(Deserialize)]
struct Capture {
    params: Tolerance,
    meshes: Vec<SceneMesh>,
}
#[derive(Deserialize)]
struct Tolerance {
    chord_height: f64,
    surface_angle_deg: f64,
    max_edge_length: f64,
    max_aspect_ratio: f64,
}
#[derive(Deserialize)]
struct SceneMesh {
    positions: Vec<[f32; 3]>,
    indices: Vec<usize>,
}

type I3 = [Interval; 3];
fn enclose(x: &Real, dec: &Decider) -> Interval {
    let ball = dec.enclose_at(x, 128).expect("analytic carrier enclosure");
    Interval::from_rational(&ball.lo()).hull(&Interval::from_rational(&ball.hi()))
}
fn dot(a: I3, b: I3) -> Interval {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
fn length(a: I3) -> Interval {
    (a[0].square() + a[1].square() + a[2].square())
        .sqrt()
        .unwrap()
}
fn delta(a: I3, b: I3) -> I3 {
    std::array::from_fn(|i| a[i] - b[i])
}
fn abs_upper(a: Interval) -> f64 {
    a.lo().abs().max(a.hi().abs())
}

struct Carrier {
    origin: I3,
    axis: I3,
    radius: Option<Interval>,
}
impl Carrier {
    fn new(surface: &Surface, dec: &Decider) -> Self {
        let (o, axis, radius) = match surface {
            Surface::Plane(p) => (&p.origin, &p.normal, None),
            Surface::Cylinder(c) => (
                &c.origin,
                &c.axis,
                Some(enclose(&c.radius_sq, dec).sqrt().unwrap()),
            ),
            _ => panic!("unproved live display carrier: {:?}", surface.kind()),
        };
        Self {
            origin: [&o.x, &o.y, &o.z].map(|v| enclose(v, dec)),
            axis: [&axis.x, &axis.y, &axis.z].map(|v| enclose(v, dec)),
            radius,
        }
    }

    fn bound(&self, points: [[f32; 3]; 3]) -> f64 {
        let n = self.axis;
        let d = points.map(|p| {
            delta(
                p.map(|x| Interval::point(f64::from(x)).unwrap()),
                self.origin,
            )
        });
        let Some(radius) = self.radius else {
            // Signed distance to the plane is affine over the whole facet.
            return d
                .map(|p| abs_upper(dot(p, n) / length(n)))
                .into_iter()
                .fold(0., f64::max);
        };
        let nn = dot(n, n);
        let radial = d.map(|p| {
            let h = dot(p, n) / nn;
            std::array::from_fn(|i| p[i] - h * n[i])
        });
        let (a, b) = [(0, 1), (1, 2), (0, 2)]
            .into_iter()
            .max_by_key(|&(a, b)| length(delta(radial[a], radial[b])).hi().to_bits())
            .unwrap();
        let w = std::array::from_fn(|i| radial[a][i] + radial[b][i]);
        let l = length(w);
        assert!(l.lo() > 0., "facet spans the cylinder axis");
        // The radial norm is convex. Its projection onto a fixed unit vector
        // is affine and bounds that norm below on the entire facet.
        let lo = radial
            .map(|p| (dot(p, w) / l).lo())
            .into_iter()
            .fold(f64::INFINITY, f64::min);
        let hi = radial
            .map(|p| length(p).hi())
            .into_iter()
            .fold(0., f64::max);
        abs_upper(Interval::new(lo, hi).unwrap() - radius)
    }
}

type Point = [u32; 3];
fn key(p: [f32; 3]) -> Point {
    p.map(|x| {
        assert!(x.is_finite());
        if x == 0. { 0 } else { x.to_bits() }
    })
}
fn check_area(points: [[f32; 3]; 3]) {
    let [a, b, c] = points.map(|p| p.map(|x| Interval::point(f64::from(x)).unwrap()));
    let u = delta(b, a);
    let v = delta(c, a);
    assert!(
        (0..3).any(|i| {
            let (j, k) = ((i + 1) % 3, (i + 2) % 3);
            !(u[j] * v[k] - u[k] * v[j]).contains_zero()
        }),
        "degenerate or unresolved facet"
    );
}
fn triangle(p: [[f32; 3]; 3]) -> [Point; 3] {
    let [a, b, c] = p.map(key);
    [a, b, c].min([b, c, a]).min([c, a, b])
}
fn closed(preview: &PreviewMesh, triangles: &[[Point; 3]]) {
    let mut uses = BTreeMap::new();
    for ([a, b, c], t) in triangles.iter().copied().zip(&preview.triangles) {
        for chord in [(a, b), (b, c), (c, a)] {
            assert!(
                uses.insert(chord, t.face).is_none(),
                "duplicate oriented chord"
            );
        }
    }
    assert!(
        uses.keys().all(|&(a, b)| uses.contains_key(&(b, a))),
        "open viewport mesh"
    );
    for e in &preview.edges {
        for p in e.vertices.windows(2) {
            let a = key(preview.positions[p[0]].map(|x| x as f32));
            let b = key(preview.positions[p[1]].map(|x| x as f32));
            assert_ne!(uses[&(a, b)], uses[&(b, a)], "unshared B-rep chord");
        }
    }
}

fn check(folder: &Path) -> Result<f64, Box<dyn std::error::Error>> {
    let recipe =
        catallax_app::soliddraft::convert(&fs::read_to_string(folder.join("recipe.json"))?)?.recipe;
    let capture: Capture = serde_json::from_slice(&fs::read(folder.join("viewport.mesh.json"))?)?;
    let [scene] = capture.meshes.as_slice() else {
        return Err("expected one live body".into());
    };
    let p = capture.params;
    assert!(p.chord_height.is_finite() && p.chord_height >= 0.0);
    assert!(p.surface_angle_deg.is_finite() && p.surface_angle_deg >= 0.0);
    assert_eq!((p.max_edge_length, p.max_aspect_ratio), (0.0, 0.0));
    let exact = |v: f64| Rational::parse(&v.to_string()).expect("finite positive tolerance");
    let mut options = MeshOptions::default();
    let bound = if p.chord_height == 0.0 {
        options.chord_tolerance.clone()
    } else {
        let requested = exact(p.chord_height);
        // The released adapter reserves half the bound for its f32 conversion.
        options.chord_tolerance = &requested / &Rational::from_i64(2);
        requested
    };
    if p.surface_angle_deg > 0.0 {
        let pi_below = Rational::from_ratio(314_159_265, 100_000_000).unwrap();
        options.angular_deviation =
            Some(&(&exact(p.surface_angle_deg) * &pi_below) / &Rational::from_i64(180));
    }
    let dec = &mut Decider::default();
    let app = catallax_app::build_incremental_without_display(&recipe, None, dec)?;
    let body = app.checkpoint.built();
    let preview = preview_mesh_body(
        &body.model,
        body.result,
        &PreviewOptions::from(&options),
        dec,
    )?;
    assert!(!preview.triangles.is_empty());
    assert_eq!(
        preview.faces.len(),
        body.model.faces_of_body(body.result).len()
    );
    let points = |t: usize| {
        preview.triangles[t]
            .vertices
            .map(|i| preview.positions[i].map(|x| x as f32))
    };
    let mut surface_bound: f64 = 0.;
    for face in &preview.faces {
        assert!(!face.triangles.is_empty());
        let surface = body
            .model
            .geom()
            .surface(body.model.face(face.face).surface);
        let carrier = Carrier::new(surface, dec);
        for t in face.triangles.clone() {
            let p = points(t);
            check_area(p);
            surface_bound = surface_bound.max(carrier.bound(p));
        }
    }
    assert!(
        surface_bound <= Interval::from_rational(&bound).lo(),
        "rendered facet exceeds chord bound"
    );
    let mut expected: Vec<_> = (0..preview.triangles.len())
        .map(|t| triangle(points(t)))
        .collect();
    closed(&preview, &expected);
    assert_eq!(
        scene.indices.len(),
        3 * expected.len(),
        "missing scene triangles"
    );
    let mut actual: Vec<_> = scene
        .indices
        .as_chunks::<3>()
        .0
        .iter()
        .map(|i| triangle(i.map(|v| scene.positions[v])))
        .collect();
    expected.sort_unstable();
    actual.sort_unstable();
    assert!(
        expected == actual,
        "live display differs from the edited recipe's checked preview"
    );
    Ok(surface_bound)
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    for folder in env::args_os().skip(1) {
        let folder = Path::new(&folder);
        let (result, work) = with_operation_budget(20_000_000, || check(folder));
        assert!(!work.exhausted);
        println!(
            "{}",
            json!({"step":folder.file_name().and_then(|x| x.to_str()),"carrier_bound_mm":result?})
        );
    }
    Ok(())
}

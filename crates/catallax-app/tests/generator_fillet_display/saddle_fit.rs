use catallax_mesh::Mesh;
use catallax_numbers::{Decider, Real};
use catallax_shapes::{SaddleRing, SplineSupport, Surface};
use catallax_topology::Model;

type Point = [f64; 3];

fn dot(a: Point, b: Point) -> f64 {
    a.into_iter().zip(b).map(|(a, b)| a * b).sum()
}

fn sub(a: Point, b: Point) -> Point {
    std::array::from_fn(|i| a[i] - b[i])
}

fn triangle_gap_sq(p: Point, [a, b, c]: [Point; 3]) -> f64 {
    let u = sub(b, a);
    let v = sub(c, a);
    let d = sub(p, a);
    let (uu, uv, vv) = (dot(u, u), dot(u, v), dot(v, v));
    let det = uu * vv - uv * uv;
    if det > 0.0 {
        let (du, dv) = (dot(d, u), dot(d, v));
        let s = (vv * du - uv * dv) / det;
        let t = (uu * dv - uv * du) / det;
        if s >= 0.0 && t >= 0.0 && s + t <= 1.0 {
            let delta = std::array::from_fn(|i| d[i] - s * u[i] - t * v[i]);
            return dot(delta, delta);
        }
    }
    [(a, b), (b, c), (c, a)]
        .into_iter()
        .map(|(a, b)| {
            let edge = sub(b, a);
            let d = sub(p, a);
            let length = dot(edge, edge);
            let t = if length > 0.0 {
                (dot(d, edge) / length).clamp(0.0, 1.0)
            } else {
                0.0
            };
            let delta = std::array::from_fn(|i| d[i] - t * edge[i]);
            dot(delta, delta)
        })
        .fold(f64::INFINITY, f64::min)
}

pub(super) fn check(case: &str, model: &Model, mesh: &Mesh) {
    let dec = Decider::default();
    let number = |x: Real| dec.enclose_at(&x, 80).unwrap().mid().to_f64_nearest();
    for face in &mesh.faces {
        let Surface::BSpline(SplineSupport::RollingSaddle(saddle)) =
            model.geom().surface(model.face(face.face).surface)
        else {
            continue;
        };
        let law = saddle.definition().law();
        let branches: &[i64] = match law.ring {
            SaddleRing::Through => std::slice::from_ref(&law.branch),
            SaddleRing::Overlap => &[-1, 1],
            SaddleRing::Crossing { .. } => continue,
        };
        let [narrow_axis, side, across, along] =
            [law.narrow_axis(), law.side_axis(), law.across, law.along]
                .map(|v| [v.x, v.y, v.z].map(number));
        let foot = [law.foot.x, law.foot.y, law.foot.z].map(number);
        let [wide, narrow, radius, gap, cosine, sine] = [
            law.wide, law.narrow, law.radius, law.gap, law.cosine, law.sine,
        ]
        .map(number);
        let [wide_side, narrow_side] = law.sides.map(|s| s as f64);
        let big = wide + wide_side * radius;
        let small = narrow + narrow_side * radius;
        let triangles: Vec<_> = mesh.triangles[face.triangles.clone()]
            .iter()
            .map(|t| t.vertices.map(|i| mesh.positions[i]))
            .collect();
        let limit = (&face.max_chord_error + &saddle.bound().upper).to_f64_nearest() + 1e-9;
        let mut samples = 0;
        for station in 0..48 {
            let theta = std::f64::consts::TAU * (f64::from(station) + 0.5) / 48.0;
            let (sin, cos) = theta.sin_cos();
            let radial: Point = std::array::from_fn(|i| cos * across[i] + sin * side[i]);
            let x = gap + small * cos;
            if x.abs() >= big {
                continue;
            }
            for &branch in branches {
                // Intersect the narrow offset cylinder's generator with the
                // wide offset cylinder to locate the rolling ball's centre.
                let y = branch as f64 * (big * big - x * x).sqrt();
                let height = (y - cosine * small * sin) / sine;
                let centre: Point =
                    std::array::from_fn(|i| foot[i] + small * radial[i] + height * narrow_axis[i]);
                let wide_ray: Point =
                    std::array::from_fn(|i| -wide_side * (x * across[i] + y * along[i]) / big);
                for t in [0.125, 0.5, 0.875] {
                    let ray = std::array::from_fn(|i| {
                        (1.0 - t) * wide_ray[i] - t * narrow_side * radial[i]
                    });
                    let length = dot(ray, ray).sqrt();
                    let p = std::array::from_fn(|i| centre[i] + radius * ray[i] / length);
                    let nearest = triangles
                        .iter()
                        .map(|&triangle| triangle_gap_sq(p, triangle))
                        .fold(f64::INFINITY, f64::min);
                    assert!(
                        nearest <= limit * limit,
                        "{case}, station {station}, branch {branch}, section {t}: {} mm > {limit} mm",
                        nearest.sqrt()
                    );
                    samples += 1;
                }
            }
        }
        assert!(samples >= 48, "{case}: too few offset-field samples");
    }
}

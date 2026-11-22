//! The patch-shell leaf: step 7.2 of 4.2a ([Checkpoint 2 decisions, section
//! 4.6]). A closed solid bounded by whole rational patches and planar
//! polygons, held as the region whose winding number is one, as the
//! closed-mesh leaf (`mesh.rs`) holds a mesh.
//!
//! The producer proves the shell closed, oriented outward and free of
//! crossings; the leaf keeps the faces, a rational rigid frame and nothing
//! of any B-rep. A point is classified by the signed count of a ray's
//! crossings with the faces: each patch's from [`PatchDistance::crossings`],
//! each polygon's in exact rationals. A ray that meets a face's border, or
//! touches a face, proves nothing, and the next direction of a fixed list is
//! tried. The value is the true signed distance: the bounds on the nearest
//! face's distance, negative inside. A point whose distance bounds reach zero
//! on a patch stays `Indeterminate`, since no exact zero test of a patch is
//! made; a point on a polygon is on the boundary exactly.
//!
//! Trimmed faces are not held (the decisions put them in 4.4).
//!
//! [Checkpoint 2 decisions, section 4.6]: ../../../docs/checkpoints/checkpoint-2-decisions.md

use catallax_numbers::{Ball, Decider, Interval, Rational, Refusal};
use catallax_recipe::{EdgeName, EdgeSelection, Origin, Role, StepId};
use catallax_shapes::{Isometry, RPoint3, RVector3, RationalPatch};

use crate::Location;
use crate::metric::PRECISION;
use crate::patch_distance::PatchDistance;
use crate::rolling_region::Range;
use crate::spline_spine::sweep::{Cap, loft::Round};

fn refusal(detail: &str) -> Refusal {
    Refusal::new("field.patch-shell", detail)
}

/// Ray directions tried in turn. Their coordinates are in no small ratio to
/// one another, so planes and lines of a recipe's own axes meet each ray in
/// isolated points; a ray that still grazes a face border is followed by the
/// next. Eight tries bound the work of a point that every ray grazes, which
/// then stays `Indeterminate`.
const DIRECTIONS: [[(i64, i64); 3]; 8] = [
    [(1, 1), (2, 7), (3, 11)],
    [(-2, 13), (1, 1), (5, 17)],
    [(3, 19), (-4, 23), (1, 1)],
    [(-1, 1), (-5, 29), (2, 31)],
    [(6, 37), (-1, 1), (-3, 41)],
    [(-5, 43), (7, 47), (-1, 1)],
    [(1, 1), (-8, 53), (-7, 59)],
    [(-9, 61), (1, 1), (10, 67)],
];

type V = [Rational; 3];

fn sub(a: &V, b: &V) -> V {
    std::array::from_fn(|i| &a[i] - &b[i])
}

fn dot(a: &V, b: &V) -> Rational {
    &(&(&a[0] * &b[0]) + &(&a[1] * &b[1])) + &(&a[2] * &b[2])
}

fn cross(a: &V, b: &V) -> V {
    [
        &(&a[1] * &b[2]) - &(&a[2] * &b[1]),
        &(&a[2] * &b[0]) - &(&a[0] * &b[2]),
        &(&a[0] * &b[1]) - &(&a[1] * &b[0]),
    ]
}

fn point(p: &RPoint3) -> V {
    [p.x.clone(), p.y.clone(), p.z.clone()]
}

/// A planar polygon face, its vertices in order counterclockwise seen from
/// outside the solid.
#[derive(Clone, Debug)]
pub struct Polygon {
    vertices: Vec<V>,
    other_rings: Vec<Vec<V>>,
    /// Outward, by Newell's sum of the vertices' cross products.
    normal: V,
    /// The coordinate axis dropped to test points in two dimensions: the one
    /// along which the normal is longest.
    dropped: usize,
}

/// Where a point of the polygon's plane lies.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Planar {
    Inside,
    Outside,
    Edge,
}

impl Polygon {
    /// A polygon from its vertices, counterclockwise seen from outside. The
    /// vertices must span a plane and lie in it exactly.
    pub fn new(vertices: Vec<RPoint3>) -> Result<Self, Refusal> {
        let vertices: Vec<V> = vertices.iter().map(point).collect();
        if vertices.len() < 3 {
            return Err(refusal("a polygon face needs three vertices"));
        }
        let mut normal: V = std::array::from_fn(|_| Rational::zero());
        for (i, a) in vertices.iter().enumerate() {
            let b = &vertices[(i + 1) % vertices.len()];
            let c = cross(a, b);
            normal = std::array::from_fn(|k| &normal[k] + &c[k]);
        }
        if normal.iter().all(Rational::is_zero) {
            return Err(refusal("a polygon face has no area"));
        }
        if vertices
            .iter()
            .any(|v| !dot(&normal, &sub(v, &vertices[0])).is_zero())
        {
            return Err(refusal("a polygon face's vertices are not in one plane"));
        }
        let dropped = (0..3)
            .max_by(|&i, &j| normal[i].abs().cmp(&normal[j].abs()))
            .expect("three axes");
        Ok(Self {
            vertices,
            other_rings: vec![],
            normal,
            dropped,
        })
    }

    /// Validated nonintersecting coplanar rings. The first ring determines
    /// the outward normal and even-odd nesting determines cap material.
    pub(crate) fn with_rings(
        first: Vec<RPoint3>,
        other_rings: Vec<Vec<RPoint3>>,
    ) -> Result<Self, Refusal> {
        let mut polygon = Self::new(first)?;
        polygon.other_rings = other_rings
            .iter()
            .map(|ring| ring.iter().map(point).collect())
            .collect();
        Ok(polygon)
    }

    fn rings(&self) -> impl Iterator<Item = &[V]> {
        std::iter::once(self.vertices.as_slice()).chain(self.other_rings.iter().map(Vec::as_slice))
    }

    fn segments(&self) -> impl Iterator<Item = (&V, &V)> {
        self.rings()
            .flat_map(|ring| ring.iter().zip(ring.iter().cycle().skip(1)))
    }

    fn flat(&self, v: &V) -> [Rational; 2] {
        let [a, b] = match self.dropped {
            0 => [1, 2],
            1 => [2, 0],
            _ => [0, 1],
        };
        [v[a].clone(), v[b].clone()]
    }

    /// Where a point of the plane lies, by the crossing number of a ray along
    /// the first flat axis, with exact rational tests.
    fn planar(&self, q: &V) -> Planar {
        let q = self.flat(q);
        let mut inside = false;
        for (a, b) in self.segments() {
            let a = self.flat(a);
            let b = self.flat(b);
            let side =
                &(&(&b[0] - &a[0]) * &(&q[1] - &a[1])) - &(&(&b[1] - &a[1]) * &(&q[0] - &a[0]));
            let within = |k: usize| {
                a[k].clone().min(b[k].clone()) <= q[k] && q[k] <= a[k].clone().max(b[k].clone())
            };
            if side.is_zero() && within(0) && within(1) {
                return Planar::Edge;
            }
            if (a[1] > q[1]) != (b[1] > q[1]) {
                // The edge's point at the height of `q`, against `q`.
                let x = &a[0] + &(&(&(&q[1] - &a[1]) * &(&b[0] - &a[0])) / &(&b[1] - &a[1]));
                if q[0] < x {
                    inside = !inside;
                }
            }
        }
        if inside {
            Planar::Inside
        } else {
            Planar::Outside
        }
    }

    /// The exact squared distance from `q` to the polygon.
    fn squared(&self, q: &V) -> Rational {
        let height = dot(&self.normal, &sub(q, &self.vertices[0]));
        let length = dot(&self.normal, &self.normal);
        let scale = &height / &length;
        let foot: V = std::array::from_fn(|k| &q[k] - &(&scale * &self.normal[k]));
        if self.planar(&foot) != Planar::Outside {
            return &(&height * &height) / &length;
        }
        let mut best: Option<Rational> = None;
        for (a, b) in self.segments() {
            let edge = sub(b, a);
            let t = (&dot(&sub(q, a), &edge) / &dot(&edge, &edge))
                .max(Rational::zero())
                .min(Rational::one());
            let near: V = std::array::from_fn(|k| &a[k] + &(&t * &edge[k]));
            let gap = sub(q, &near);
            let value = dot(&gap, &gap);
            best = Some(match best {
                Some(b) if b <= value => b,
                _ => value,
            });
        }
        best.expect("three vertices")
    }

    /// The signed crossing of the open ray `q + t d`, `t > 0`: one if it
    /// leaves through the polygon, minus one if it enters, zero if it misses.
    /// `None` where it runs in the plane or meets an edge.
    fn crossing(&self, q: &V, d: &V) -> Option<i64> {
        let along = dot(&self.normal, d);
        let height = dot(&self.normal, &sub(&self.vertices[0], q));
        if along.is_zero() {
            return if height.is_zero() { None } else { Some(0) };
        }
        let t = &height / &along;
        if t <= Rational::zero() {
            return Some(0);
        }
        let hit: V = std::array::from_fn(|k| &q[k] + &(&t * &d[k]));
        match self.planar(&hit) {
            Planar::Outside => Some(0),
            Planar::Edge => None,
            Planar::Inside => Some(if along > Rational::zero() { 1 } else { -1 }),
        }
    }

    /// Tentative squared distance and the side of the plane, from binary64
    /// arithmetic only.
    fn fast_sample(&self, q: [f64; 3]) -> (f64, f64, [f64; 3]) {
        let v0 = self.vertices[0].each_ref().map(Rational::to_f64_nearest);
        let n = self.normal.each_ref().map(Rational::to_f64_nearest);
        let length = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
        let side = (0..3).map(|k| (q[k] - v0[k]) * n[k]).sum::<f64>() / length;
        let foot: [f64; 3] = std::array::from_fn(|k| q[k] - side * n[k] / length);
        let [a, b] = match self.dropped {
            0 => [1, 2],
            1 => [2, 0],
            _ => [0, 1],
        };
        let mut inside = false;
        let mut edge = (f64::INFINITY, q);
        for (s, t) in self.segments() {
            let [s, t] = [s, t].map(|v| v.each_ref().map(Rational::to_f64_nearest));
            if (s[b] > foot[b]) != (t[b] > foot[b])
                && foot[a] < s[a] + (foot[b] - s[b]) * (t[a] - s[a]) / (t[b] - s[b])
            {
                inside = !inside;
            }
            let e: [f64; 3] = std::array::from_fn(|k| t[k] - s[k]);
            let ee = e[0] * e[0] + e[1] * e[1] + e[2] * e[2];
            let along = ((0..3).map(|k| (q[k] - s[k]) * e[k]).sum::<f64>() / ee).clamp(0.0, 1.0);
            let candidate = std::array::from_fn(|k| s[k] + along * e[k]);
            let distance = (0..3).map(|k| (q[k] - s[k] - along * e[k]).powi(2)).sum();
            if distance < edge.0 {
                edge = (distance, candidate);
            }
        }
        if inside {
            (side * side, side, foot)
        } else {
            (edge.0, side, edge.1)
        }
    }
}

/// A whole patch face of a shell.
#[derive(Clone, Debug)]
pub struct PatchFace {
    distance: PatchDistance,
    support: RationalPatch,
    /// Whether the patch's natural normal `S_u x S_v` points out of the solid.
    outward: bool,
}

impl PatchFace {
    /// A whole patch, with `outward` true where its natural normal points out
    /// of the solid.
    pub fn new(patch: &RationalPatch, outward: bool, dec: &mut Decider) -> Result<Self, Refusal> {
        Ok(Self {
            distance: PatchDistance::new(patch, dec)?,
            support: patch.clone(),
            outward,
        })
    }
}

/// One face of a shell.
#[derive(Clone, Debug)]
pub enum Face {
    /// A whole rational patch.
    Patch(Box<PatchFace>),
    /// A planar polygon.
    Polygon(Box<Polygon>),
}

/// A closed shell of whole patches and planar polygons, placed by a rational
/// rigid frame.
#[derive(Clone, Debug)]
pub struct Shell {
    faces: Vec<Face>,
    /// The faces' bounds in their own coordinates.
    local_bounds: [[Rational; 2]; 3],
    /// World point = rotation * local point + translation.
    rotation: [[Rational; 3]; 3],
    translation: V,
    inverse_rotation: [[Rational; 3]; 3],
    loft_sides: Vec<Origin>,
    round: Option<Round>,
}

fn apply(m: &[[Rational; 3]; 3], v: &V) -> V {
    std::array::from_fn(|i| dot(&m[i], v))
}

fn transpose(m: &[[Rational; 3]; 3]) -> [[Rational; 3]; 3] {
    std::array::from_fn(|i| std::array::from_fn(|j| m[j][i].clone()))
}

fn ball(interval: Interval) -> Ball {
    let end =
        |x: f64| Ball::from_rational(&Rational::from_f64_exact(x).expect("finite"), PRECISION);
    end(interval.lo()).hull(&end(interval.hi()))
}

impl Shell {
    /// Unsigned boundary distance for a proximity test with a positive limit.
    /// Once a witness is within the limit, zero remains a valid lower bound.
    pub fn boundary_distance_within(
        &self,
        p: &RPoint3,
        limit: f64,
        dec: &mut Decider,
    ) -> Result<Interval, Refusal> {
        let q = self.local(&point(p));
        let at = RPoint3::new(q[0].clone(), q[1].clone(), q[2].clone());
        let mut lo = f64::INFINITY;
        let mut hi = f64::INFINITY;
        for face in &self.faces {
            let distance = match face {
                Face::Polygon(polygon) => {
                    dec.charge_work(1)?;
                    Interval::from_rational(&polygon.squared(&q))
                        .sqrt()
                        .expect("non-negative")
                }
                // Leave half the requested limit for a nonzero lower bound.
                Face::Patch(patch) => {
                    patch
                        .distance
                        .distance_within(&at, limit * 0.5, dec)?
                        .distance
                }
            };
            lo = lo.min(distance.lo());
            hi = hi.min(distance.hi());
            if hi <= limit {
                return Ok(Interval::new(0.0, hi).expect("non-negative bounds"));
            }
        }
        Ok(Interval::new(lo, hi).expect("ordered bounds"))
    }

    /// A shell of `faces` in place. The caller proves them a closed,
    /// outward, crossing-free shell.
    pub fn new(faces: Vec<Face>) -> Result<Self, Refusal> {
        if faces.is_empty() {
            return Err(refusal("a shell needs faces"));
        }
        let mut bounds: Option<[[Rational; 2]; 3]> = None;
        let mut grow = |p: [Rational; 3]| {
            bounds = Some(match bounds.take() {
                None => p.map(|x| [x.clone(), x]),
                Some(b) => std::array::from_fn(|k| {
                    [
                        b[k][0].clone().min(p[k].clone()),
                        b[k][1].clone().max(p[k].clone()),
                    ]
                }),
            });
        };
        for face in &faces {
            match face {
                Face::Polygon(polygon) => polygon.rings().flatten().for_each(|v| grow(v.clone())),
                Face::Patch(patch) => {
                    let hull = patch.distance.hull();
                    let end = |x: f64| Rational::from_f64_exact(x).expect("finite control");
                    grow(hull.map(|k| end(k.lo())));
                    grow(hull.map(|k| end(k.hi())));
                }
            }
        }
        let identity: [[Rational; 3]; 3] = std::array::from_fn(|i| {
            std::array::from_fn(|j| {
                if i == j {
                    Rational::one()
                } else {
                    Rational::zero()
                }
            })
        });
        Ok(Self {
            faces,
            local_bounds: bounds.expect("faces"),
            inverse_rotation: identity.clone(),
            rotation: identity,
            translation: std::array::from_fn(|_| Rational::zero()),
            loft_sides: Vec::new(),
            round: None,
        })
    }

    pub(crate) fn named_loft(mut self, sides: Vec<Origin>) -> Self {
        self.loft_sides = sides;
        self
    }

    pub(crate) fn fillet(
        &self,
        edges: &EdgeSelection,
        step: &StepId,
        radius: Rational,
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        let unsupported = || {
            Refusal::new(
                "fillet.unsupported-spline",
                "a field spline fillet needs one named edge of an untrimmed convex ruled loft",
            )
        };
        let EdgeSelection::Refs(refs) = edges else {
            return Err(unsupported());
        };
        let [reference] = refs.as_slice() else {
            return Err(unsupported());
        };
        let sides = &self.loft_sides;
        if self.round.is_some()
            || sides.is_empty()
            || reference.near.is_some()
            || !reference.touching.is_empty()
        {
            return Err(unsupported());
        }
        let adjacent: Vec<_> = (0..sides.len())
            .map(|i| [i, (i + 1) % sides.len()])
            .filter(|&[i, j]| {
                let [a, b] = &reference.between;
                (a.matches(&sides[i]) && b.matches(&sides[j]))
                    || (a.matches(&sides[j]) && b.matches(&sides[i]))
            })
            .collect();
        let [selected] = adjacent.as_slice() else {
            return Err(Refusal::new(
                "fillet.lost-edge",
                "the side origins do not select one shared loft edge",
            ));
        };
        let edge = EdgeName::new(sides[selected[0]].clone(), sides[selected[1]].clone())
            .map_err(|e| refusal(&e.to_string()))?;
        let origin = Origin {
            step: step.clone(),
            role: Role::Blend { edge },
        };
        let patches = self.faces[..sides.len()]
            .iter()
            .map(|face| match face {
                Face::Patch(face) => Ok((&face.support, face.outward)),
                Face::Polygon(_) => Err(unsupported()),
            })
            .collect::<Result<Vec<_>, _>>()?;
        let cap = |face: &Face| -> Result<Cap, Refusal> {
            let Face::Polygon(face) = face else {
                return Err(unsupported());
            };
            let [x, y, z] = face.vertices[0].clone();
            let [nx, ny, nz] = face.normal.clone();
            Ok((RPoint3::new(x, y, z), RVector3::new(-nx, -ny, -nz)))
        };
        let caps = [
            cap(&self.faces[sides.len()])?,
            cap(&self.faces[sides.len() + 1])?,
        ];
        let mut shell = self.clone();
        shell.round = Some(Round::prove(
            origin, *selected, &patches, &caps, radius, dec,
        )?);
        Ok(shell)
    }

    pub(crate) fn blend_origin(&self) -> Option<&Origin> {
        self.round.as_ref().map(|r| &r.origin)
    }

    pub(crate) fn blend_distance(
        &self,
        p: &RPoint3,
        bound: &Rational,
        dec: &mut Decider,
    ) -> Result<Option<Range>, Refusal> {
        let Some(round) = &self.round else {
            return Ok(None);
        };
        let [x, y, z] = self.local(&point(p));
        round
            .distance_at_most(&RPoint3::new(x, y, z), bound, dec)
            .map(Some)
    }

    fn local_face_bound(
        &self,
        p: &RPoint3,
        round: &Round,
        width: &Rational,
        dec: &mut Decider,
    ) -> Result<Option<Rational>, Refusal> {
        let query = point(p);
        let caps = self.faces.iter().filter_map(|face| {
            let Face::Polygon(poly) = face else {
                return None;
            };
            let factor = &dot(&poly.normal, &sub(&query, &poly.vertices[0]))
                / &dot(&poly.normal, &poly.normal);
            let [x, y, z] = std::array::from_fn(|i| &query[i] - &(&factor * &poly.normal[i]));
            Some(RPoint3::new(x, y, z))
        });
        let square = width * width;
        let mut candidates: Vec<_> = caps
            .chain(round.side_witnesses(p))
            .filter_map(|q| {
                let d = p.minus(&q).norm_squared();
                (d <= square).then_some((d, q))
            })
            .collect();
        candidates.sort_by(|a, b| a.0.cmp(&b.0));
        for (d, candidate) in candidates {
            if round.location(&candidate, dec)? == Location::Boundary {
                let upper = Interval::from_rational(&d).sqrt().unwrap().hi();
                return Ok(Rational::from_f64_exact(upper));
            }
        }
        Ok(None)
    }

    fn local(&self, p: &V) -> V {
        apply(&self.inverse_rotation, &sub(p, &self.translation))
    }

    /// Composes a rational rigid motion, reflections included, in front of
    /// the frame. Queries map back into the faces' own coordinates, so a
    /// reflected copy keeps its winding numbers.
    pub fn transform(&mut self, map: &Isometry) -> Result<(), Refusal> {
        let frame = || {
            Refusal::new(
                "field.patch-shell",
                "a shell pose must be an exact rational rigid motion",
            )
        };
        let m = map.rational_matrix().ok_or_else(frame)?;
        let m: [[Rational; 3]; 3] =
            std::array::from_fn(|i| std::array::from_fn(|j| m[i][j].clone()));
        let t = map.translation.as_rational().ok_or_else(frame)?;
        let rotation = std::array::from_fn(|i| {
            std::array::from_fn(|j| {
                (0..3).fold(Rational::zero(), |s, k| s + &m[i][k] * &self.rotation[k][j])
            })
        });
        let shifted = apply(&m, &self.translation);
        self.translation = [&shifted[0] + &t.x, &shifted[1] + &t.y, &shifted[2] + &t.z];
        self.inverse_rotation = transpose(&rotation);
        self.rotation = rotation;
        Ok(())
    }

    /// Exact world bounds of the corners of the faces' own bounding box.
    pub fn bounds(&self) -> [[Rational; 2]; 3] {
        let mut out: Option<[[Rational; 2]; 3]> = None;
        for corner in 0..8 {
            let l: V = std::array::from_fn(|i| self.local_bounds[i][corner >> i & 1].clone());
            let r = apply(&self.rotation, &l);
            let w: V = std::array::from_fn(|i| &r[i] + &self.translation[i]);
            out = Some(match out {
                None => w.map(|x| [x.clone(), x]),
                Some(b) => std::array::from_fn(|i| {
                    [
                        b[i][0].clone().min(w[i].clone()),
                        b[i][1].clone().max(w[i].clone()),
                    ]
                }),
            });
        }
        out.expect("eight corners")
    }

    /// Outward world boxes whose union holds the shell's faces, so a box
    /// that meets none of them holds one winding number throughout: each
    /// polygon's vertex box, and each patch's [`PatchDistance::boxes`].
    pub fn surface_boxes(&self) -> Vec<[Interval; 3]> {
        let rotation = self
            .rotation
            .each_ref()
            .map(|row| row.each_ref().map(Interval::from_rational));
        let translation = self.translation.each_ref().map(Interval::from_rational);
        let world = |local: [Interval; 3]| -> [Interval; 3] {
            std::array::from_fn(|i| {
                rotation[i][0] * local[0]
                    + rotation[i][1] * local[1]
                    + rotation[i][2] * local[2]
                    + translation[i]
            })
        };
        let mut out = Vec::new();
        for face in &self.faces {
            match face {
                Face::Polygon(polygon) => {
                    let mut hull: Option<[Interval; 3]> = None;
                    for v in polygon.rings().flatten() {
                        let p = v.each_ref().map(Interval::from_rational);
                        hull = Some(match hull {
                            None => p,
                            Some(h) => std::array::from_fn(|k| h[k].hull(&p[k])),
                        });
                    }
                    out.push(world(hull.expect("three vertices")));
                }
                Face::Patch(patch) => out.extend(patch.distance.boxes().into_iter().map(world)),
            }
        }
        if let Some(round) = &self.round {
            out.extend(round.boxes().map(world));
        }
        out
    }

    /// Classifies a world point and bounds its distance to the shell.
    pub fn classify(
        &self,
        p: &RPoint3,
        dec: &mut Decider,
    ) -> Result<(Location, Interval), Refusal> {
        let q = self.local(&point(p));
        if let Some(round) = &self.round {
            let p = RPoint3::new(q[0].clone(), q[1].clone(), q[2].clone());
            let location = round.location(&p, dec)?;
            let distance = self.rounded_distance(&p, &q, round, location, dec)?;
            return Ok((location, distance));
        }
        let mut lo = f64::INFINITY;
        let mut hi = f64::INFINITY;
        for face in &self.faces {
            let distance = match face {
                Face::Polygon(polygon) => {
                    let squared = polygon.squared(&q);
                    if squared.is_zero() {
                        let zero = Interval::point(0.0).expect("zero");
                        return Ok((Location::Boundary, zero));
                    }
                    Interval::from_rational(&squared)
                        .sqrt()
                        .expect("non-negative")
                }
                Face::Patch(patch) => {
                    let at = RPoint3::new(q[0].clone(), q[1].clone(), q[2].clone());
                    patch.distance.distance(&at, dec)?.distance
                }
            };
            lo = lo.min(distance.lo());
            hi = hi.min(distance.hi());
        }
        let distance = Interval::new(lo, hi).expect("ordered bounds");
        if lo <= 0.0 {
            return Ok((Location::Indeterminate, distance));
        }
        if (0..3).any(|k| q[k] < self.local_bounds[k][0] || q[k] > self.local_bounds[k][1]) {
            return Ok((Location::Outside, distance));
        }
        for direction in DIRECTIONS {
            let d: V = direction.map(|(n, m)| Rational::from_ratio(n, m).expect("nonzero"));
            if let Some(winding) = self.winding(&q, &d, dec)? {
                let location = if winding == 0 {
                    Location::Outside
                } else {
                    Location::Inside
                };
                return Ok((location, distance));
            }
        }
        Ok((Location::Indeterminate, distance))
    }

    /// The winding number of the shell about `q` along `d`, or `None` where
    /// the ray grazes a face.
    fn winding(&self, q: &V, d: &V, dec: &mut Decider) -> Result<Option<i64>, Refusal> {
        let mut total = 0;
        for face in &self.faces {
            let count = match face {
                Face::Polygon(polygon) => polygon.crossing(q, d),
                Face::Patch(patch) => patch
                    .distance
                    .crossings(q, d, dec)?
                    .map(|n| if patch.outward { n } else { -n }),
            };
            match count {
                Some(n) => total += n,
                None => return Ok(None),
            }
        }
        Ok(Some(total))
    }

    /// Exact classification of a rational world point.
    pub fn location(&self, p: &[Rational; 3], dec: &mut Decider) -> Result<Location, Refusal> {
        if let Some(round) = &self.round {
            let [x, y, z] = self.local(p);
            return round.location(&RPoint3::new(x, y, z), dec);
        }
        let p = RPoint3::new(p[0].clone(), p[1].clone(), p[2].clone());
        Ok(self.classify(&p, dec)?.0)
    }

    /// Enclosure of the true signed distance, negative inside.
    pub fn proven_distance(
        &self,
        p: &RPoint3,
        location: Location,
        dec: &mut Decider,
    ) -> Result<Ball, Refusal> {
        if location == Location::Boundary {
            return Ok(Ball::from_rational(&Rational::zero(), PRECISION));
        }
        let distance = if let Some(round) = &self.round {
            let q = self.local(&point(p));
            let local = RPoint3::new(q[0].clone(), q[1].clone(), q[2].clone());
            self.rounded_distance(&local, &q, round, location, dec)?
        } else {
            self.classify(p, dec)?.1
        };
        let distance = ball(distance);
        Ok(match location {
            Location::Inside => -distance,
            Location::Outside => distance,
            _ => (-distance.clone()).hull(&distance),
        })
    }

    fn rounded_distance(
        &self,
        p: &RPoint3,
        q: &V,
        round: &Round,
        location: Location,
        dec: &mut Decider,
    ) -> Result<Interval, Refusal> {
        if location == Location::Boundary {
            return Ok(Interval::point(0.0).unwrap());
        }
        let width = Rational::from_ratio(1, 1_000_000).expect("one nanometre");
        if let Some(upper) = self.local_face_bound(p, round, &width, dec)? {
            return Ok(Interval::new(0.0, Interval::from_rational(&upper).hi()).unwrap());
        }
        let blend = round.distance(p, &width, dec)?;
        let mut lower = Interval::from_rational(&blend.lower).lo();
        let mut upper = blend
            .upper
            .as_ref()
            .map_or(f64::INFINITY, |v| Interval::from_rational(v).hi());
        for (i, face) in self.faces.iter().enumerate() {
            let distance = match face {
                Face::Polygon(polygon) => {
                    Interval::from_rational(&polygon.squared(q)).sqrt().unwrap()
                }
                Face::Patch(patch) => patch.distance.distance(p, dec)?.distance,
            };
            lower = lower.min(distance.lo());
            // Only an uncut face can furnish a witness from a whole-face
            // minimum. The blend itself already supplies an upper witness.
            if i < self.loft_sides.len() && !round.selected.contains(&i) {
                upper = upper.min(distance.hi());
            }
        }
        Interval::new(lower, upper)
            .ok_or_else(|| refusal("the rounded shell distance is not ordered"))
    }

    /// Tentative signed distance from binary64 searches: the nearest face's
    /// distance, signed by the side of that face's normal.
    pub fn fast_distance(&self, p: [f64; 3]) -> f64 {
        let t = self.translation.each_ref().map(Rational::to_f64_nearest);
        let m = self
            .inverse_rotation
            .each_ref()
            .map(|r| r.each_ref().map(Rational::to_f64_nearest));
        let s: [f64; 3] = std::array::from_fn(|k| p[k] - t[k]);
        let q: [f64; 3] = std::array::from_fn(|i| m[i][0] * s[0] + m[i][1] * s[1] + m[i][2] * s[2]);
        let mut best = (f64::INFINITY, 1.0);
        for (i, face) in self.faces.iter().enumerate() {
            let (squared, side, candidate) = match face {
                Face::Polygon(polygon) => polygon.fast_sample(q),
                Face::Patch(patch) => {
                    let (squared, side, candidate) = patch.distance.fast_sample(q);
                    (squared, if patch.outward { side } else { -side }, candidate)
                }
            };
            if squared < best.0
                && !self.round.as_ref().is_some_and(|round| {
                    (i >= self.loft_sides.len() || round.selected.contains(&i))
                        && round.preview(candidate).0
                })
            {
                best = (squared, side);
            }
        }
        let magnitude = crate::metric::Scalar::sqrt(&best.0);
        if let Some(round) = &self.round {
            let (_, blend) = round.preview(q);
            blend.abs().min(magnitude).copysign(blend)
        } else if best.1 < 0.0 {
            -magnitude
        } else {
            magnitude
        }
    }
}

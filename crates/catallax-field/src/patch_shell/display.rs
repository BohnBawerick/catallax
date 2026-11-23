use std::collections::BTreeMap;

use catallax_numbers::Sign;
use catallax_shapes::tessellation::{
    RationalSpan, checked_triangle_bound, chord_bound_interval, rational_spans,
};

use super::*;
use crate::{
    FaceId,
    display::{Mesh, Options, Triangle},
};

struct Sample {
    position: [f64; 3],
    error: Rational,
}

struct Builder<'a> {
    mesh: Mesh,
    options: &'a Options,
    rotation: [[Interval; 3]; 3],
    translation: [Interval; 3],
    reflected: bool,
    vertices: BTreeMap<[u64; 3], usize>,
}

fn exhausted() -> Refusal {
    Refusal::new(
        "field.display-budget",
        "patch shell exceeds its display limits",
    )
}

fn exact(value: f64) -> Rational {
    Rational::from_f64_exact(value).expect("a finite display bound")
}

fn singleton(value: f64) -> Interval {
    Interval::point(value).expect("a finite chart coordinate")
}

impl Shell {
    pub(crate) fn mesh(
        &self,
        faces: &[FaceId],
        options: &Options,
        dec: &mut Decider,
    ) -> Result<Mesh, Refusal> {
        assert_eq!(faces.len(), self.faces.len());
        let mut builder = Builder {
            mesh: Mesh {
                positions: Vec::new(),
                triangles: Vec::new(),
                max_deviation: Rational::zero(),
                cells: 0,
            },
            options,
            rotation: self
                .rotation
                .each_ref()
                .map(|row| row.each_ref().map(Interval::from_rational)),
            translation: self.translation.each_ref().map(Interval::from_rational),
            reflected: dot(
                &self.rotation[0],
                &cross(&self.rotation[1], &self.rotation[2]),
            )
            .sign()
                == Sign::Negative,
            vertices: BTreeMap::new(),
        };
        let mut sides = BTreeMap::new();
        for (face, &id) in self.faces.iter().zip(faces) {
            if let Face::Patch(patch) = face {
                let carrier = patch.support.carrier();
                let spans = rational_spans(carrier, dec)?;
                let [span]: [RationalSpan; 1] = spans.try_into().map_err(|_| unsupported_loft())?;
                if carrier.degrees() != [1, 1]
                    || carrier
                        .weights()
                        .iter()
                        .any(|w| w.as_rational() != Some(Rational::one()))
                {
                    return Err(unsupported_loft());
                }
                sides.insert(id.index(), span.affine_hessian());
            }
        }
        // Every side uses the same dyadic fractions along the two section
        // polygons. Its other two edges are straight and need only endpoints.
        let half_tolerance = Interval::from_rational(&options.tolerance.mul_pow2(-1));
        let mut columns = 1_usize;
        loop {
            dec.charge_work(sides.len() as u64)?;
            let width = [singleton(1.0 / columns as f64), singleton(1.0)];
            if sides
                .values()
                .all(|hessian| chord_bound_interval(width, hessian).hi() <= half_tolerance.lo())
            {
                break;
            }
            columns = columns
                .checked_mul(2)
                .filter(|&n| n <= options.max_cells)
                .ok_or_else(exhausted)?;
        }
        for (face, &id) in self.faces.iter().zip(faces) {
            match face {
                Face::Patch(patch) => {
                    let hessian = &sides[&id.index()];
                    let width = [singleton(1.0 / columns as f64), singleton(1.0)];
                    let chord = exact(chord_bound_interval(width, hessian).hi());
                    builder.patch(patch, id, columns, &chord, dec)?;
                }
                Face::Polygon(polygon) => builder.polygon(polygon, id, columns, dec)?,
            }
        }
        Ok(builder.mesh)
    }
}

fn unsupported_loft() -> Refusal {
    Refusal::new(
        "field.display-unsupported",
        "patch-shell loft side is not one bilinear span",
    )
}

impl Builder<'_> {
    fn visit(&mut self, dec: &mut Decider) -> Result<(), Refusal> {
        if self.mesh.cells == self.options.max_cells {
            return Err(exhausted());
        }
        dec.charge_work(1)?;
        self.mesh.cells += 1;
        Ok(())
    }

    fn sample(&self, local: [Interval; 3]) -> Result<Sample, Refusal> {
        let world: [Interval; 3] = std::array::from_fn(|i| {
            (0..3).fold(self.translation[i], |a, j| {
                a + self.rotation[i][j] * local[j]
            })
        });
        if world.iter().any(|x| !x.is_finite()) {
            return Err(exhausted());
        }
        let position = world.map(|bound| bound.mid());
        let error = world
            .into_iter()
            .zip(position)
            .fold(Rational::zero(), |sum, (x, value)| {
                let value = exact(value);
                sum + (&value - &exact(x.lo()))
                    .abs()
                    .max((&value - &exact(x.hi())).abs())
            });
        Ok(Sample { position, error })
    }

    fn triangle(
        &mut self,
        samples: [&Sample; 3],
        face: FaceId,
        outward: bool,
        dec: &mut Decider,
    ) -> Result<(), Refusal> {
        let [a, b, c] = samples.map(|sample| sample.position.map(singleton));
        dec.charge_work(6)?;
        let normal = cross_interval(
            std::array::from_fn(|i| b[i] - a[i]),
            std::array::from_fn(|i| c[i] - a[i]),
        );
        if !normal
            .iter()
            .any(|component| component.is_finite() && !component.contains_zero())
        {
            return Err(Refusal::new(
                "field.display-precision",
                "a rounded patch-shell triangle has no proved area",
            ));
        }
        let mut triangle = [0; 3];
        for (index, sample) in triangle.iter_mut().zip(samples) {
            let key = sample.position.map(f64::to_bits);
            let next = self.mesh.positions.len();
            if let Some(&known) = self.vertices.get(&key) {
                *index = known;
            } else {
                if next == self.options.max_vertices {
                    return Err(exhausted());
                }
                self.vertices.insert(key, next);
                self.mesh.positions.push(sample.position);
                *index = next;
            }
        }
        if triangle[0] == triangle[1] || triangle[1] == triangle[2] || triangle[0] == triangle[2] {
            return Err(exhausted());
        }
        if outward == self.reflected {
            triangle.swap(1, 2);
        }
        self.mesh.triangles.push(Triangle {
            vertices: triangle,
            face,
        });
        Ok(())
    }

    fn patch(
        &mut self,
        patch: &PatchFace,
        face: FaceId,
        columns: usize,
        chord: &Rational,
        dec: &mut Decider,
    ) -> Result<(), Refusal> {
        let [a, d, b, c] = patch.support.carrier().controls() else {
            return Err(unsupported_loft());
        };
        let local = |p: &catallax_shapes::Point3| {
            let [x, y, z] =
                [&p.x, &p.y, &p.z].map(|x| x.as_rational().ok_or_else(unsupported_loft));
            Ok::<_, Refusal>([x?, y?, z?])
        };
        let [a, b, c, d] = [local(a)?, local(b)?, local(c)?, local(d)?];
        let mut left = [
            self.section_point(&a, &b, 0, columns, dec)?,
            self.section_point(&d, &c, 0, columns, dec)?,
        ];
        for i in 1..=columns {
            self.visit(dec)?;
            let right = [
                self.section_point(&a, &b, i, columns, dec)?,
                self.section_point(&d, &c, i, columns, dec)?,
            ];
            let error = left
                .iter()
                .chain(&right)
                .map(|sample| &sample.error)
                .max()
                .unwrap();
            let bound = checked_triangle_bound(chord, error, &self.options.tolerance)
                .ok_or_else(exhausted)?;
            self.triangle([&left[0], &right[0], &right[1]], face, patch.outward, dec)?;
            self.triangle([&left[0], &right[1], &left[1]], face, patch.outward, dec)?;
            self.mesh.max_deviation = self.mesh.max_deviation.clone().max(bound);
            left = right;
        }
        Ok(())
    }

    fn section_point(
        &self,
        a: &V,
        b: &V,
        i: usize,
        columns: usize,
        dec: &mut Decider,
    ) -> Result<Sample, Refusal> {
        dec.charge_work(1)?;
        let t = Rational::from_ratio(i as i64, columns as i64).unwrap();
        let p: V = std::array::from_fn(|axis| &a[axis] + &(&t * &(&b[axis] - &a[axis])));
        self.sample(p.each_ref().map(Interval::from_rational))
    }

    fn polygon(
        &mut self,
        polygon: &Polygon,
        face: FaceId,
        columns: usize,
        dec: &mut Decider,
    ) -> Result<(), Refusal> {
        if !polygon.other_rings.is_empty() {
            return Err(Refusal::new(
                "field.display-unsupported",
                "patch-shell cap has holes",
            ));
        }
        let ring = &polygon.vertices;
        // A fan from an interior point keeps every subdivided rim segment.
        for (a, b) in polygon.segments() {
            for at in ring {
                dec.charge_work(1)?;
                if dot(&polygon.normal, &cross(&sub(b, a), &sub(at, a))).sign() == Sign::Negative {
                    return Err(Refusal::new(
                        "field.display-unsupported",
                        "patch-shell cap is not convex",
                    ));
                }
            }
        }
        dec.charge_work(ring.len() as u64)?;
        let count = Rational::from_i64(ring.len() as i64);
        let centre: V = std::array::from_fn(|axis| {
            let sum = ring.iter().fold(Rational::zero(), |sum, p| &sum + &p[axis]);
            &sum / &count
        });
        let centre = self.sample(centre.each_ref().map(Interval::from_rational))?;
        for (a, b) in polygon.segments() {
            let mut last = self.section_point(a, b, 0, columns, dec)?;
            for i in 1..=columns {
                self.visit(dec)?;
                let next = self.section_point(a, b, i, columns, dec)?;
                let error = (&centre.error).max(&last.error).max(&next.error);
                let bound =
                    checked_triangle_bound(&Rational::zero(), error, &self.options.tolerance)
                        .ok_or_else(exhausted)?;
                self.triangle([&centre, &last, &next], face, true, dec)?;
                self.mesh.max_deviation = self.mesh.max_deviation.clone().max(bound);
                last = next;
            }
        }
        Ok(())
    }
}

fn cross_interval(a: [Interval; 3], b: [Interval; 3]) -> [Interval; 3] {
    std::array::from_fn(|i| a[(i + 1) % 3] * b[(i + 2) % 3] - a[(i + 2) % 3] * b[(i + 1) % 3])
}

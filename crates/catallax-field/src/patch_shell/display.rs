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

type Cell = [[f64; 2]; 2];

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
        for (face, &id) in self.faces.iter().zip(faces) {
            match face {
                Face::Patch(patch) => {
                    for span in rational_spans(patch.support.carrier(), dec)? {
                        builder.patch(&span, id, patch.outward, dec)?;
                    }
                }
                Face::Polygon(polygon) => builder.polygon(polygon, id, dec)?,
            }
        }
        Ok(builder.mesh)
    }
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
        span: &RationalSpan,
        face: FaceId,
        outward: bool,
        dec: &mut Decider,
    ) -> Result<(), Refusal> {
        let mut pending: Vec<Cell> = vec![[[0.0, 1.0]; 2]];
        let hessian = span.affine_hessian();
        while let Some(cell) = pending.pop() {
            self.visit(dec)?;
            let [[u0, u1], [v0, v1]] = cell;
            let width = [singleton(u1) - singleton(u0), singleton(v1) - singleton(v0)];
            let chord = chord_bound_interval(width, &hessian);
            if !chord.is_finite() {
                return Err(exhausted());
            }
            if chord.hi() <= self.options.tolerance.to_f64_nearest() {
                let [a, b, c, d] = [[u0, v0], [u1, v0], [u1, v1], [u0, v1]]
                    .map(|uv| self.sample(span.evaluate(uv, dec)?));
                let [a, b, c, d] = [a?, b?, c?, d?];
                let replacement = [&a, &b, &c, &d]
                    .into_iter()
                    .map(|s| &s.error)
                    .max()
                    .unwrap();
                if let Some(bound) =
                    checked_triangle_bound(&exact(chord.hi()), replacement, &self.options.tolerance)
                {
                    self.triangle([&a, &b, &c], face, outward, dec)?;
                    self.triangle([&a, &c, &d], face, outward, dec)?;
                    self.mesh.max_deviation = self.mesh.max_deviation.clone().max(bound);
                    continue;
                }
            }
            let axis = (0..2)
                .max_by(|&a, &b| {
                    (hessian[2 * a] * width[a].square())
                        .hi()
                        .total_cmp(&(hessian[2 * b] * width[b].square()).hi())
                        .then_with(|| width[a].hi().total_cmp(&width[b].hi()))
                })
                .unwrap();
            let [lo, hi] = cell[axis];
            let mid = lo.midpoint(hi);
            if !(lo < mid && mid < hi) {
                return Err(exhausted());
            }
            let [mut first, mut second] = [cell; 2];
            first[axis][1] = mid;
            second[axis][0] = mid;
            pending.extend([second, first]);
        }
        Ok(())
    }

    fn polygon(
        &mut self,
        polygon: &Polygon,
        face: FaceId,
        dec: &mut Decider,
    ) -> Result<(), Refusal> {
        if !polygon.other_rings.is_empty() {
            return Err(Refusal::new(
                "field.display-unsupported",
                "patch-shell cap has holes",
            ));
        }
        let ring = &polygon.vertices;
        // Check both coverage and winding before using a fan from the first corner.
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
        let samples = ring
            .iter()
            .map(|p| self.sample(p.each_ref().map(Interval::from_rational)))
            .collect::<Result<Vec<_>, _>>()?;
        let error = samples.iter().map(|p| &p.error).max().unwrap();
        let bound = checked_triangle_bound(&Rational::zero(), error, &self.options.tolerance)
            .ok_or_else(exhausted)?;
        for i in 1..samples.len() - 1 {
            self.visit(dec)?;
            self.triangle([&samples[0], &samples[i], &samples[i + 1]], face, true, dec)?;
        }
        self.mesh.max_deviation = self.mesh.max_deviation.clone().max(bound);
        Ok(())
    }
}

fn cross_interval(a: [Interval; 3], b: [Interval; 3]) -> [Interval; 3] {
    std::array::from_fn(|i| a[(i + 1) % 3] * b[(i + 2) % 3] - a[(i + 2) % 3] * b[(i + 1) % 3])
}

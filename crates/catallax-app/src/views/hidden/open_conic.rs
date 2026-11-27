//! Rational charts of finite parabola and hyperbola arcs in a drawing view.

use catallax_intersect::CurveForm;
use catallax_numbers::{Decider, Rational, Real, RealPoly, Refusal, StrictSign};
use catallax_shapes::{Point3, Vector3};

use super::{Geo, geo::Span, point64, rational, root};
use crate::queries::sign;
use crate::views::float::double;

pub(super) struct OpenConic {
    pub(super) origin: Point3,
    pub(super) normal: Vector3,
    pub(super) edge_on: bool,
    u: Vector3,
    w: Vector3,
    ends: [Real; 2],
    hyperbola: bool,
}

impl OpenConic {
    pub(super) fn new(
        form: CurveForm,
        ends: [Point3; 2],
        d: &Vector3,
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        let unit = |v: &Vector3| v.scaled(&(&Real::one() / &root(&v.norm_squared())));
        let (origin, normal, u, w, hyperbola) = match form {
            CurveForm::Parabola {
                vertex,
                normal,
                axis_dir,
                focal,
            } => {
                let u = unit(&axis_dir).scaled(&(&Real::one() / &(&Real::from_i64(4) * &focal)));
                let w = unit(&normal.cross(&axis_dir));
                (vertex, normal, u, w, false)
            }
            CurveForm::Hyperbola {
                center,
                normal,
                axis_dir,
                transverse_sq,
                conjugate_sq,
                branch,
            } => {
                let axis = unit(&axis_dir);
                let u = axis.scaled(&match branch {
                    StrictSign::Positive => root(&transverse_sq),
                    StrictSign::Negative => -root(&transverse_sq),
                });
                let w = unit(&normal.cross(&axis_dir)).scaled(&root(&conjugate_sq));
                (center, normal, u, w, true)
            }
            _ => unreachable!("an open conic chart requires an open conic"),
        };
        // The parabola parameter is y. On the selected hyperbola branch it
        // is tanh(angle/2), or y/(1+x) in scaled branch coordinates.
        let ends = ends.map(|p| {
            let offset = p.minus(&origin);
            let y = offset.dot(&w) / w.norm_squared();
            if hyperbola {
                y / (Real::one() + offset.dot(&u) / u.norm_squared())
            } else {
                y
            }
        });
        let edge_on = sign(&normal.dot(d), dec)? == 0;
        Ok(Self {
            origin,
            normal,
            edge_on,
            u,
            w,
            ends,
            hyperbola,
        })
    }

    // A mixed degree-two homogeneous control, including the two endpoint
    // controls when a == b. All hyperbola weights are positive on (-1, 1).
    fn control(&self, a: &Real, b: &Real) -> (Point3, Real) {
        let product = a * b;
        let sum = a + b;
        let (x, y, weight) = if self.hyperbola {
            (&Real::one() + &product, sum, &Real::one() - &product)
        } else {
            (product, sum / Real::from_i64(2), Real::one())
        };
        let p = self.origin.plus(
            &self
                .u
                .scaled(&(&x / &weight))
                .plus(&self.w.scaled(&(&y / &weight))),
        );
        (p, weight)
    }

    pub(super) fn geo(&self, floor: &Rational, dec: &mut Decider) -> Result<Geo, Refusal> {
        let [a, b] = &self.ends;
        let mut ctrl = Vec::with_capacity(3);
        for (a, b) in [(a, a), (a, b), (b, b)] {
            let (p, weight) = self.control(a, b);
            let [x, y, z] = point64(&p, floor, dec)?;
            let w = double(&weight, floor, dec)?;
            if w <= 0.0 {
                return Err(Refusal::new(
                    "app.view-curve",
                    "open conic weight is too small for projection",
                ));
            }
            ctrl.push([x * w, y * w, z * w, w]);
        }
        Ok(Geo::Spline {
            spans: vec![Span {
                lo: 0.0,
                hi: 1.0,
                ctrl,
            }],
        })
    }

    pub(super) fn sample(&self, s: f64) -> Result<Point3, Refusal> {
        let [a, b] = &self.ends;
        let t = a + &((b - a) * Real::from(rational(s)?));
        Ok(self.control(&t, &t).0)
    }

    fn roots(
        &self,
        coefficients: [Real; 3],
        dec: &mut Decider,
    ) -> Result<Option<Vec<Point3>>, Refusal> {
        let polynomial = RealPoly::from_terms(
            1,
            coefficients
                .into_iter()
                .enumerate()
                .map(|(degree, c)| (vec![degree as u32], c))
                .collect(),
        );
        if polynomial.is_zero(dec)? {
            return Ok(None);
        }
        let mut points = Vec::new();
        for t in polynomial.real_roots(dec)? {
            if self.hyperbola && sign(&(Real::one() - t.square()), dec)? <= 0 {
                continue;
            }
            points.push(self.control(&t, &t).0);
        }
        Ok(Some(points))
    }

    pub(super) fn plane_points(
        &self,
        origin: &Point3,
        normal: &Vector3,
        dec: &mut Decider,
    ) -> Result<Option<Vec<Point3>>, Refusal> {
        let k = self.origin.minus(origin).dot(normal);
        let x = self.u.dot(normal);
        let y = self.w.dot(normal);
        self.roots(
            if self.hyperbola {
                [&k + &x, y * Real::from_i64(2), &x - &k]
            } else {
                [k, y, x]
            },
            dec,
        )
    }

    pub(super) fn folds(&self, view: &Vector3, dec: &mut Decider) -> Result<Vec<Point3>, Refusal> {
        let m = self.normal.cross(view);
        let x = self.u.dot(&m) * Real::from_i64(2);
        let y = self.w.dot(&m);
        // The numerator of the derivative of the projected coordinate.
        Ok(self
            .roots(
                if self.hyperbola {
                    [y.clone(), x, y]
                } else {
                    [y, x, Real::zero()]
                },
                dec,
            )?
            .unwrap_or_default())
    }
}

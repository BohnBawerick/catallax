//! Exact ruled cylindrical draft with a transverse, oblique neutral plane.
//!
//! Horizontal sections have polar radius `A(z) + b.cos(u) + c.sin(u)`.
//! The constructor proves a positive, strictly convex section throughout the
//! finite height. This proves an embedded ruled annulus, not a convex solid.

use std::sync::Arc;

use catallax_numbers::{Ball, Decider, Rational, Real, RealPoly, Refusal, Sign};

use crate::{BSplineSurface, Cylinder, Isometry, Plane, Point3, Pose, Vector3, math};

fn unsupported(detail: &str) -> Refusal {
    Refusal::new("draft.unsupported-oblique-arrangement", detail)
}
fn unit(v: &Vector3) -> Vector3 {
    v.scaled(&(&Real::one() / &math::exact_sqrt(&v.norm_squared())))
}

#[derive(Debug)]
struct Data {
    origin: Point3,
    axis: Vector3,
    across: Vector3,
    side: Vector3,
    heights: [Real; 2],
    intercept: Real,
    tangent: Real,
    lean: [Real; 2],
    carrier: BSplineSurface,
}

/// A finite exact ruled draft annulus with complete positive polar sections.
/// Construction retains the unsquared radial branch, including the axis's
/// interior side; the quartic equation alone is not its definition.
#[derive(Clone, Debug)]
pub struct ObliqueDraft {
    data: Arc<Data>,
}

impl ObliqueDraft {
    /// Constructs the neutral-curve generator law for a pull-parallel cylinder.
    /// `tangent` includes the material sense: negate it for a bore wall.
    /// Caps are perpendicular to the unit `pull`, at the supplied axial heights
    /// relative to the source cylinder's origin. Unproved folds remain refused.
    pub fn new(
        cylinder: &Cylinder,
        pull: Vector3,
        neutral: &Plane,
        tangent: Real,
        heights: [Real; 2],
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        math::validate_point(&cylinder.origin, dec)?;
        math::validate_vector(&cylinder.axis, dec)?;
        math::validate_vector(&cylinder.ref_dir, dec)?;
        math::directions(&cylinder.axis, &cylinder.ref_dir, dec)?;
        math::validate_vector(&pull, dec)?;
        math::validate_vector(&neutral.normal, dec)?;
        math::validate_point(&neutral.origin, dec)?;
        math::positive(&cylinder.radius_sq, "shapes.non-positive-radius", dec)?;
        math::validate(&tangent, dec)?;
        for h in &heights {
            math::validate(h, dec)?;
        }
        if math::sign(&(&pull.norm_squared() - &Real::one()), dec)? != Sign::Zero
            || math::sign(&pull.cross(&cylinder.axis).norm_squared(), dec)? != Sign::Zero
        {
            return Err(unsupported(
                "oblique draft requires a unit pull parallel to the cylinder axis",
            ));
        }
        if math::sign(&(&heights[1] - &heights[0]), dec)? != Sign::Positive {
            return Err(unsupported(
                "oblique draft requires a positive finite axial interval",
            ));
        }
        let rate = pull.dot(&neutral.normal);
        if math::sign(&rate, dec)? == Sign::Zero {
            return Err(unsupported(
                "an axis-parallel neutral plane has no closed transverse neutral curve",
            ));
        }
        let across = unit(&cylinder.ref_dir);
        let side = pull.cross(&across);
        if math::sign(&across.dot(&pull), dec)? != Sign::Zero {
            return Err(unsupported(
                "cylinder reference direction must be perpendicular to pull",
            ));
        }
        let radius = math::exact_sqrt(&cylinder.radius_sq);
        let neutral_height = &neutral.origin.minus(&cylinder.origin).dot(&neutral.normal) / &rate;
        let intercept = &radius + &(&tangent * &neutral_height);
        let factor = -(&(&radius * &tangent) / &rate);
        let lean = [neutral.normal.dot(&across), neutral.normal.dot(&side)].map(|x| &factor * &x);
        let lean_sq = &lean[0].square() + &lean[1].square();
        for height in &heights {
            let a = &intercept - &(&tangent * height);
            // A > 2|lean| gives rho > 0 and positive planar curvature:
            // rho^2 + 2 rho'^2 - rho rho'' >= (A-|lean|)(A-2|lean|).
            if math::sign(&a, dec)? != Sign::Positive
                || math::sign(&(&a.square() - &(&Real::from_i64(4) * &lean_sq)), dec)?
                    != Sign::Positive
            {
                return Err(unsupported(
                    "a section approaches the axis or loses the proved convex-section margin",
                ));
            }
        }
        let carrier = carrier(
            &cylinder.origin,
            &pull,
            &across,
            &side,
            &heights,
            &intercept,
            &tangent,
            &lean,
            dec,
        )?;
        Ok(Self {
            data: Arc::new(Data {
                origin: cylinder.origin.clone(),
                axis: pull,
                across,
                side,
                heights,
                intercept,
                tangent,
                lean,
                carrier,
            }),
        })
    }

    /// Exact rational degree-(4,1) carrier, on `[0,4] x [0,1]`.
    pub fn carrier(&self) -> &BSplineSurface {
        &self.data.carrier
    }
    /// Whole-domain two-way chord proof on a shared tensor grid. On each
    /// quarter, `phi=2 atan(tau)` has `|phi'|<=2`, `|phi''|<=4`.
    /// Thus `|S_uu|<=8 A_max+24|lean|`, `|S_uv|<=2|tangent| height`,
    /// and `S_vv=0`. Taylor interpolation bounds each complete triangle,
    /// including the ruled surface's mixed derivative, not only its rails.
    pub fn display_grid(
        &self,
        tolerance: &Rational,
        max_points: usize,
        dec: &Decider,
    ) -> Result<crate::EnvelopeGrid, Refusal> {
        if *tolerance <= Rational::zero() || max_points < 4 {
            return Err(Refusal::new(
                "mesh.invalid-options",
                "positive chord tolerance and at least four grid points required",
            ));
        }
        let bits = dec.policy().max_precision_bits.min(128);
        if bits == 0 {
            return Err(Refusal::new(
                "numbers.precision-limit",
                "no draft display precision available",
            ));
        }
        let lean = dec
            .enclose_at(
                &(&self.data.lean[0].square() + &self.data.lean[1].square()).sqrt(),
                bits,
            )?
            .hi();
        let amplitude = self
            .data
            .heights
            .each_ref()
            .map(|h| self.amplitude(h))
            .iter()
            .try_fold(Rational::zero(), |best, a| -> Result<_, Refusal> {
                Ok(best.max(dec.enclose_at(a, bits)?.hi()))
            })?;
        let tangent = dec.enclose_at(&self.data.tangent, bits)?;
        let height = dec
            .enclose_at(&(&self.data.heights[1] - &self.data.heights[0]), bits)?
            .hi();
        let uu = &(&Rational::from_i64(8) * &amplitude) + &(&Rational::from_i64(24) * &lean);
        let uv = &(&Rational::from_i64(2) * &tangent.lo().abs().max(tangent.hi().abs())) * &height;
        let (mut u, mut v) = (1usize, 1usize);
        let chord = loop {
            dec.check_deadline()?;
            if u > (i64::MAX / 4) as usize
                || v > (i64::MAX / 2) as usize
                || u.checked_mul(4)
                    .and_then(|n| n.checked_add(1))
                    .and_then(|n| v.checked_add(1).and_then(|rows| n.checked_mul(rows)))
                    .is_none_or(|n| n > max_points)
            {
                return Err(Refusal::new(
                    "mesh.resource-limit",
                    "oblique draft tensor grid exceeds the fixed vertex limit",
                ));
            }
            let du = Rational::from_ratio(1, u as i64).unwrap();
            let dv = Rational::from_ratio(1, v as i64).unwrap();
            let first = (&(&uu * &du) * &du).mul_pow2(-3);
            let mixed = (&(&uv * &du) * &dv).mul_pow2(-2);
            if &first + &mixed <= *tolerance {
                break first + mixed;
            }
            if first > tolerance.mul_pow2(-1) {
                u *= 2;
            } else {
                v *= 2;
            }
        };
        let (columns, rows) = (4 * u + 1, v + 1);
        let mut points = Vec::with_capacity(columns * rows);
        for i in 0..columns {
            dec.check_deadline()?;
            let s = Rational::from_ratio(i as i64, u as i64).unwrap();
            for j in 0..rows {
                let t = Rational::from_ratio(j as i64, v as i64).unwrap();
                let p = self.chart_point(&s, &t)?;
                points.push(math::point(&p, bits, dec)?);
            }
        }
        Ok(crate::EnvelopeGrid {
            columns,
            rows,
            points,
            chord,
        })
    }
    /// Original axis reference point.
    pub fn origin(&self) -> &Point3 {
        &self.data.origin
    }
    /// Unit pull direction.
    pub fn axis(&self) -> &Vector3 {
        &self.data.axis
    }
    /// Unit section reference direction.
    pub fn across(&self) -> &Vector3 {
        &self.data.across
    }
    /// Unit section direction completing the right-handed frame.
    pub fn side(&self) -> &Vector3 {
        &self.data.side
    }
    /// Exact finite cap heights relative to `origin`.
    pub fn heights(&self) -> &[Real; 2] {
        &self.data.heights
    }
    /// Constant angular coefficients in the polar section law.
    pub fn lean(&self) -> &[Real; 2] {
        &self.data.lean
    }
    /// Signed radial shrink per unit pull displacement.
    pub fn tangent(&self) -> &Real {
        &self.data.tangent
    }
    /// Angular-independent section radius at axial height `z`.
    pub fn amplitude(&self, z: &Real) -> Real {
        &self.data.intercept - &(&self.data.tangent * z)
    }
    /// Exact local coordinates, without floating-point frame reconstruction.
    pub fn coordinates(&self, p: &Point3) -> [Real; 3] {
        let d = p.minus(self.origin());
        [d.dot(self.across()), d.dot(self.side()), d.dot(self.axis())]
    }
    /// Candidate quartic in world coordinates. Its zero set also contains the
    /// entire axis: callers must retain the finite height and radial branch.
    /// Inside the finite height, its sign agrees with the radial side away
    /// from the axis, since `A > 2|lean|` makes the discarded factor positive.
    pub fn implicit_candidate(&self) -> RealPoly {
        let displacement = math::displacement(self.origin());
        let [x, y, z] =
            [self.across(), self.side(), self.axis()].map(|v| math::dot(&displacement, v));
        let squared = &(&x * &x) + &(&y * &y);
        let radial = &(&squared - &x.scale(&self.data.lean[0])) - &y.scale(&self.data.lean[1]);
        let a = RealPoly::constant(3, self.data.intercept.clone()) - z.scale(&self.data.tangent);
        &(&radial * &radial) - &(&(&a * &a) * &squared)
    }
    /// Exact point from a unit section ray and an axial height.
    pub fn point(&self, ray: [&Real; 2], z: &Real) -> Point3 {
        let radius =
            &self.amplitude(z) + &(&(&self.data.lean[0] * ray[0]) + &(&self.data.lean[1] * ray[1]));
        self.origin()
            .plus(&self.axis().scaled(z))
            .plus(&self.across().scaled(&(&radius * ray[0])))
            .plus(&self.side().scaled(&(&radius * ray[1])))
    }
    /// Exact chart point, evaluated from the generator law rather than control
    /// interpolation. Both finite parameter domains are checked.
    pub fn chart_point(&self, u: &Rational, v: &Rational) -> Result<Point3, Refusal> {
        if *u < Rational::zero()
            || *u > Rational::from_i64(4)
            || *v < Rational::zero()
            || *v > Rational::one()
        {
            return Err(Refusal::new(
                "shapes.bspline-domain",
                "oblique draft chart parameter leaves its finite domain",
            ));
        }
        let floor = u.floor();
        let tau = u - &Rational::from(floor.clone());
        let denominator = &Rational::one() + &(&tau * &tau);
        let cos = &(&Rational::one() - &(&tau * &tau)) / &denominator;
        let sin = &(&Rational::from_i64(2) * &tau) / &denominator;
        let [cos, sin] = match floor.to_i64().unwrap() % 4 {
            0 => [cos, sin],
            1 => [-sin, cos],
            2 => [-cos, -sin],
            _ => [sin, -cos],
        }
        .map(Real::from);
        let z = &self.data.heights[0]
            + &(&Real::from(v.clone()) * &(&self.data.heights[1] - &self.data.heights[0]));
        Ok(self.point([&cos, &sin], &z))
    }
    /// Radial inside/boundary/outside sign on the finite axial interval.
    pub fn side_at(&self, point: &Point3, dec: &mut Decider) -> Result<Sign, Refusal> {
        math::validate_point(point, dec)?;
        let [x, y, z] = self.coordinates(point);
        self.check_height(&z, dec)?;
        let squared = &x.square() + &y.square();
        if math::sign(&squared, dec)? == Sign::Zero {
            return Ok(Sign::Negative);
        }
        let radius = math::exact_sqrt(&squared);
        let residual = &(&(&squared - &(&self.data.lean[0] * &x)) - &(&self.data.lean[1] * &y))
            - &(&self.amplitude(&z) * &radius);
        math::sign(&residual, dec)
    }
    fn check_height(&self, z: &Real, dec: &mut Decider) -> Result<(), Refusal> {
        if math::sign(&(z - &self.data.heights[0]), dec)? == Sign::Negative
            || math::sign(&(&self.data.heights[1] - z), dec)? == Sign::Negative
        {
            return Err(Refusal::new(
                "shapes.bspline-domain",
                "point leaves the finite oblique draft height",
            ));
        }
        Ok(())
    }
    /// Outward normal at an actual point of the finite ruled support.
    pub fn normal_at(&self, point: &Point3, dec: &mut Decider) -> Result<Vector3, Refusal> {
        if self.side_at(point, dec)? != Sign::Zero {
            return Err(Refusal::new(
                "shapes.bspline-domain",
                "normal point is not on the oblique draft support",
            ));
        }
        let [x, y, z] = self.coordinates(point);
        let radius = math::exact_sqrt(&(&x.square() + &y.square()));
        let factor = &Real::from_i64(2) - &(&self.amplitude(&z) / &radius);
        Ok(self
            .across()
            .scaled(&(&(&factor * &x) - &self.data.lean[0]))
            .plus(&self.side().scaled(&(&(&factor * &y) - &self.data.lean[1])))
            .plus(&self.axis().scaled(&(&self.data.tangent * &radius))))
    }
    /// All exact line meetings, filtered against the unsquared finite branch.
    /// A line contained in the carrier refuses instead of publishing a miss.
    pub fn line_hits(
        &self,
        origin: &Point3,
        direction: &Vector3,
        dec: &mut Decider,
    ) -> Result<Vec<Point3>, Refusal> {
        math::validate_point(origin, dec)?;
        math::validate_vector(direction, dec)?;
        if math::sign(&direction.norm_squared(), dec)? == Sign::Zero {
            return Err(Refusal::new(
                "shapes.zero-direction",
                "line direction must be nonzero",
            ));
        }
        let start = self.coordinates(origin);
        let rate = [
            direction.dot(self.across()),
            direction.dot(self.side()),
            direction.dot(self.axis()),
        ];
        let parameter = RealPoly::var(1, 0);
        let c = |r: &Real| RealPoly::constant(1, r.clone());
        let [x, y, z] = std::array::from_fn(|i| &c(&start[i]) + &(&parameter * &c(&rate[i])));
        let squared = &(&x * &x) + &(&y * &y);
        let a = &c(&self.data.intercept) - &(&c(&self.data.tangent) * &z);
        let radial = &(&squared - &(&c(&self.data.lean[0]) * &x)) - &(&c(&self.data.lean[1]) * &y);
        let equation = &(&radial * &radial) - &(&(&a * &a) * &squared);
        if equation.is_zero(dec)? {
            return Err(unsupported(
                "line lies in the oblique draft algebraic carrier",
            ));
        }
        let mut hits = Vec::new();
        for t in equation.real_roots(dec)? {
            let point = origin.plus(&direction.scaled(&t));
            let z = self.coordinates(&point)[2].clone();
            if math::sign(&(&z - &self.data.heights[0]), dec)? == Sign::Negative
                || math::sign(&(&self.data.heights[1] - &z), dec)? == Sign::Negative
            {
                continue;
            }
            if self.side_at(&point, dec)? == Sign::Zero {
                hits.push(point);
            }
        }
        Ok(hits)
    }
    /// Volume divided by pi, from independently integrable polar sections.
    pub fn volume_over_pi(&self) -> Real {
        let [lo, hi] = &self.data.heights;
        let a = self.amplitude(lo);
        let b = self.amplitude(hi);
        let squared_lean = &self.data.lean[0].square() + &self.data.lean[1].square();
        (hi - lo)
            * (&(&(&a.square() + &(&a * &b)) + &b.square()) / &Real::from_i64(3)
                + &squared_lean / &Real::from_i64(2))
    }
    /// Exact rigid transport of the retained law and rational carrier.
    pub fn transformed(&self, pose: &Pose) -> Self {
        self.mapped(&Isometry::from(pose))
    }
    pub(crate) fn mapped(&self, pose: &Isometry) -> Self {
        let d = &self.data;
        Self {
            data: Arc::new(Data {
                origin: pose.apply_point(&d.origin),
                axis: pose.apply_vector(&d.axis),
                across: pose.apply_vector(&d.across),
                side: pose.apply_vector(&d.side),
                heights: d.heights.clone(),
                intercept: d.intercept.clone(),
                tangent: d.tangent.clone(),
                lean: d.lean.clone(),
                carrier: d.carrier.mapped(pose),
            }),
        }
    }
    /// Positive retained-definition identity for geometric-query caches.
    pub fn shares_definition_nodes(&self, other: &Self) -> bool {
        Arc::ptr_eq(&self.data, &other.data)
    }
    /// Outward enclosures of the actual carrier controls and positive weights.
    pub fn enclosed_controls(&self, dec: &Decider) -> Result<Vec<([Ball; 3], Ball)>, Refusal> {
        self.carrier()
            .controls()
            .iter()
            .zip(self.carrier().weights())
            .map(|(p, w)| {
                Ok((
                    [
                        dec.enclose_at(&p.x, 128)?,
                        dec.enclose_at(&p.y, 128)?,
                        dec.enclose_at(&p.z, 128)?,
                    ],
                    dec.enclose_at(w, 128)?,
                ))
            })
            .collect()
    }
}

fn product(a: &[Real; 3], b: &[Real; 3]) -> [Real; 5] {
    let choose = [1, 4, 6, 4, 1];
    let small = [1, 2, 1];
    std::array::from_fn(|k| {
        (0..3)
            .filter(|&i| k >= i && k - i < 3)
            .fold(Real::zero(), |sum, i| {
                let j = k - i;
                sum + &(&a[i] * &b[j])
                    * &Real::from(Rational::from_ratio(small[i] * small[j], choose[k]).unwrap())
            })
    })
}
#[allow(clippy::too_many_arguments)]
fn carrier(
    origin: &Point3,
    axis: &Vector3,
    across: &Vector3,
    side: &Vector3,
    heights: &[Real; 2],
    intercept: &Real,
    tangent: &Real,
    lean: &[Real; 2],
    dec: &mut Decider,
) -> Result<BSplineSurface, Refusal> {
    let wide = [1, 1, 2].map(Real::from_i64);
    let denominator = product(&wide, &wide);
    let mut controls = Vec::with_capacity(34);
    let mut weights = Vec::with_capacity(34);
    let mut scale = Real::one();
    for (span, (c, s)) in [(1, 0), (0, 1), (-1, 0), (0, -1)].into_iter().enumerate() {
        let flat = [1, 1, 0];
        let rise = [0, 1, 2];
        let cos = std::array::from_fn(|i| Real::from_i64(c * flat[i] - s * rise[i]));
        let sin = std::array::from_fn(|i| Real::from_i64(s * flat[i] + c * rise[i]));
        let rows: [[Point3; 5]; 2] = std::array::from_fn(|v| {
            let a = intercept - &(tangent * &heights[v]);
            let rho = std::array::from_fn(|i| {
                &(&a * &wide[i]) + &(&(&lean[0] * &cos[i]) + &(&lean[1] * &sin[i]))
            });
            let x = product(&rho, &cos);
            let y = product(&rho, &sin);
            std::array::from_fn(|i| {
                origin
                    .plus(&axis.scaled(&heights[v]))
                    .plus(&across.scaled(&(&x[i] / &denominator[i])))
                    .plus(&side.scaled(&(&y[i] / &denominator[i])))
            })
        });
        for i in usize::from(span != 0)..5 {
            for row in &rows {
                controls.push(row[i].clone());
                weights.push(&denominator[i] * &scale);
            }
        }
        scale = scale * Real::from_i64(4);
    }
    let knot = |x, n| vec![Rational::from_i64(x); n];
    let u = [knot(0, 5), knot(1, 4), knot(2, 4), knot(3, 4), knot(4, 5)].concat();
    BSplineSurface::new(
        [4, 1],
        [
            u,
            vec![
                Rational::zero(),
                Rational::zero(),
                Rational::one(),
                Rational::one(),
            ],
        ],
        [17, 2],
        controls,
        weights,
        dec,
    )
}

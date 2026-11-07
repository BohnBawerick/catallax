//! Exact contact-section fillets on a straight edge between two planes.
//!
//! Each section is a rational quadratic circular arc. Its endpoint feet and
//! tangent-intersection control are affine in edge station, so a linear radius
//! law has an exact degree-(2,1) carrier. A zero endpoint collapses one complete
//! boundary row to the original vertex, not to a positive-radius substitute.

use std::sync::Arc;

use catallax_numbers::{Ball, Decider, Rational, Real, RealPoly, Refusal, Sign};

use crate::{BSplineSurface, Circle, EnvelopeGrid, Isometry, Plane, Point3, Pose, Vector3, math};

fn bad(detail: &str) -> Refusal {
    Refusal::new("fillet.unsupported-variable-contact", detail)
}
fn domain(detail: &str) -> Refusal {
    Refusal::new("shapes.bspline-domain", detail)
}
fn unit(v: &Vector3) -> Vector3 {
    v.scaled(&(&Real::one() / &math::exact_sqrt(&v.norm_squared())))
}

#[derive(Debug)]
struct Data {
    supports: [Plane; 2],
    same_sense: [bool; 2],
    edge: [Point3; 2],
    radii: [Real; 2],
    zero: [bool; 2],
    convex: bool,
    along: Vector3,
    length_sq: Real,
    shift: Vector3,
    rays: [Vector3; 2],
    turn: Vector3,
    weight: Real,
    orientation: i64,
    carrier: BSplineSurface,
}

/// A checked exact variable-radius band, with its oriented input edge and law.
/// The interior chart is embedded and regular. At a zero endpoint its closed
/// parameter rectangle has the usual cone-apex quotient on that boundary row.
#[derive(Clone, Debug)]
pub struct LinearFillet {
    data: Arc<Data>,
}

impl LinearFillet {
    /// Constructs the shorter contact sections of two nonparallel planes.
    /// The input edge lies in both planes. Radius is affine from `edge[0]` to
    /// `edge[1]`, non-negative at both ends and positive in the open interval.
    /// Equal endpoints also define a regular exact cylindrical band.
    pub fn new(
        supports: [Plane; 2],
        same_sense: [bool; 2],
        edge: [Point3; 2],
        radii: [Real; 2],
        convex: bool,
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        for p in &edge {
            math::validate_point(p, dec)?;
        }
        let along = edge[1].minus(&edge[0]);
        let length_sq = along.norm_squared();
        math::positive(&length_sq, "fillet.invalid-edge", dec)?;
        for plane in &supports {
            math::validate_point(&plane.origin, dec)?;
            math::validate_vector(&plane.normal, dec)?;
            math::positive(&plane.normal.norm_squared(), "fillet.invalid-edge", dec)?;
            if math::sign(&plane.normal.dot(&along), dec)? != Sign::Zero
                || math::sign(&edge[0].minus(&plane.origin).dot(&plane.normal), dec)? != Sign::Zero
            {
                return Err(bad("the oriented straight edge must lie in both planes"));
            }
        }
        let signs = [math::sign(&radii[0], dec)?, math::sign(&radii[1], dec)?];
        if signs.contains(&Sign::Negative) || signs == [Sign::Zero, Sign::Zero] {
            return Err(Refusal::new(
                "fillet.non-positive-radius",
                "linear radius endpoints must be non-negative and not both zero",
            ));
        }
        let normals = std::array::from_fn::<_, 2, _>(|i| {
            unit(&supports[i].normal).scaled(&Real::from_i64(if same_sense[i] { 1 } else { -1 }))
        });
        let dot = normals[0].dot(&normals[1]);
        if math::sign(&(&Real::one() - &dot.square()), dec)? != Sign::Positive {
            return Err(bad("parallel planes have no regular contact section"));
        }
        let k = Real::from_i64(if convex { 1 } else { -1 });
        let shift = normals[0]
            .plus(&normals[1])
            .scaled(&(-(&k / &(&Real::one() + &dot))));
        let rays = normals.map(|normal| normal.scaled(&k));
        let turn = unit(&rays[0].cross(&rays[1]));
        let orientation = match math::sign(&turn.dot(&along), dec)? {
            Sign::Positive => 1,
            Sign::Negative => -1,
            Sign::Zero => return Err(bad("contact-section turn is not along the edge")),
        };
        let weight = math::exact_sqrt(&(&(&Real::one() + &dot) / &Real::from_i64(2)));
        let mut controls = Vec::with_capacity(6);
        for offset in [shift.plus(&rays[0]), Vector3::zero(), shift.plus(&rays[1])] {
            for i in 0..2 {
                controls.push(edge[i].plus(&offset.scaled(&radii[i])));
            }
        }
        let carrier = BSplineSurface::new(
            [2, 1],
            [
                vec![
                    Rational::zero(),
                    Rational::zero(),
                    Rational::zero(),
                    Rational::one(),
                    Rational::one(),
                    Rational::one(),
                ],
                vec![
                    Rational::zero(),
                    Rational::zero(),
                    Rational::one(),
                    Rational::one(),
                ],
            ],
            [3, 2],
            controls,
            vec![
                Real::one(),
                Real::one(),
                weight.clone(),
                weight.clone(),
                Real::one(),
                Real::one(),
            ],
            dec,
        )?;
        // In each row the middle control is the two contact tangents'
        // intersection, exactly the original edge station. The positive
        // quadratic weights therefore give the exact shorter circle arc for
        // every station, not only at sampled sections. Affine station and a
        // positive interior radius prove injectivity and nonzero chart rank.
        Ok(Self {
            data: Arc::new(Data {
                supports,
                same_sense,
                edge,
                radii,
                zero: signs.map(|s| s == Sign::Zero),
                convex,
                along,
                length_sq,
                shift,
                rays,
                turn,
                weight,
                orientation,
                carrier,
            }),
        })
    }

    /// Exact represented geometry on `[0,1] x [0,1]`, arc then station.
    pub fn carrier(&self) -> &BSplineSurface {
        &self.data.carrier
    }
    /// Original oriented input edge endpoints.
    pub fn edge(&self) -> &[Point3; 2] {
        &self.data.edge
    }
    /// Endpoint radii of the retained exact law.
    pub fn radii(&self) -> &[Real; 2] {
        &self.data.radii
    }
    /// Input supports, before any trim replacement.
    pub fn supports(&self) -> &[Plane; 2] {
        &self.data.supports
    }
    /// Input face senses relative to their original plane supports.
    pub fn same_sense(&self) -> [bool; 2] {
        self.data.same_sense
    }
    /// Whether the band removes convex material.
    pub fn convex(&self) -> bool {
        self.data.convex
    }
    /// Full directed input edge displacement.
    pub fn along(&self) -> &Vector3 {
        &self.data.along
    }
    /// Centre displacement per unit radius, perpendicular to the edge.
    pub fn shift(&self) -> &Vector3 {
        &self.data.shift
    }
    /// Unit contact rays, in arc-parameter order.
    pub fn rays(&self) -> &[Vector3; 2] {
        &self.data.rays
    }
    /// Unit normal of every section circle, from first foot to second.
    pub fn turn(&self) -> &Vector3 {
        &self.data.turn
    }
    /// Sign of the chart normal relative to the radial quadratic gradient.
    pub fn orientation(&self) -> i64 {
        self.data.orientation
    }
    /// Whether one endpoint is an exact point contact.
    pub fn zero_end(&self, end: usize) -> bool {
        self.data.zero[end]
    }
    /// One exact contact foot at an input endpoint.
    pub fn foot(&self, support: usize, end: usize) -> &Point3 {
        &self.carrier().controls()[support * 4 + end]
    }
    /// The original edge station coordinate of a world point.
    pub fn station(&self, point: &Point3) -> Real {
        &point.minus(&self.edge()[0]).dot(self.along()) / &self.data.length_sq
    }
    /// Radius at a station, with linear extension when the caller needs it.
    pub fn radius(&self, station: &Real) -> Real {
        &self.radii()[0] + &(&(&self.radii()[1] - &self.radii()[0]) * station)
    }
    fn centre(&self, station: &Real) -> Point3 {
        self.edge()[0]
            .plus(&self.along().scaled(station))
            .plus(&self.shift().scaled(&self.radius(station)))
    }
    /// Section circle at an endpoint; `None` is an exact collapsed row.
    pub fn section(&self, end: usize) -> Option<Circle> {
        (!self.zero_end(end)).then(|| Circle {
            center: self.edge()[end].plus(&self.shift().scaled(&self.radii()[end])),
            normal: self.turn().clone(),
            radius_sq: self.radii()[end].square(),
            ref_dir: self.rays()[0].clone(),
        })
    }
    /// Exact point from the retained circular-section formula.
    pub fn point(&self, arc: &Rational, station: &Rational) -> Result<Point3, Refusal> {
        if *arc < Rational::zero()
            || *arc > Rational::one()
            || *station < Rational::zero()
            || *station > Rational::one()
        {
            return Err(domain(
                "linear fillet point leaves its finite parameter domain",
            ));
        }
        let t = Real::from(arc.clone());
        let s = &Real::one() - &t;
        let middle = &(&Real::from_i64(2) * &self.data.weight) * &(&t * &s);
        let denominator = &(&s.square() + &middle) + &t.square();
        let first = self.shift().plus(&self.rays()[0]).scaled(&s.square());
        let last = self.shift().plus(&self.rays()[1]).scaled(&t.square());
        let offset = first.plus(&last).scaled(&(&Real::one() / &denominator));
        let v = Real::from(station.clone());
        Ok(self.edge()[0]
            .plus(&self.along().scaled(&v))
            .plus(&offset.scaled(&self.radius(&v))))
    }
    fn in_station(&self, station: &Real, dec: &mut Decider) -> Result<bool, Refusal> {
        Ok(math::sign(station, dec)? != Sign::Negative
            && math::sign(&(&Real::one() - station), dec)? != Sign::Negative)
    }
    fn in_sector(&self, ray: &Vector3, dec: &mut Decider) -> Result<bool, Refusal> {
        Ok(
            math::sign(&self.rays()[0].cross(ray).dot(self.turn()), dec)? != Sign::Negative
                && math::sign(&ray.cross(&self.rays()[1]).dot(self.turn()), dec)? != Sign::Negative,
        )
    }
    fn section_data(&self, point: &Point3) -> (Real, Real, Vector3) {
        let station = self.station(point);
        let radius = self.radius(&station);
        let centre = self.edge()[0]
            .plus(&self.along().scaled(&station))
            .plus(&self.shift().scaled(&radius));
        (station, radius, point.minus(&centre))
    }
    fn branch_side(
        &self,
        station: &Real,
        radius: &Real,
        ray: &Vector3,
        dec: &mut Decider,
    ) -> Result<Sign, Refusal> {
        if !self.in_station(station, dec)? || !self.in_sector(ray, dec)? {
            return Err(domain("point leaves the finite linear fillet branch"));
        }
        let sign = math::sign(&(&ray.norm_squared() - &radius.square()), dec)?;
        Ok(if self.orientation() > 0 {
            sign
        } else {
            match sign {
                Sign::Positive => Sign::Negative,
                Sign::Negative => Sign::Positive,
                Sign::Zero => Sign::Zero,
            }
        })
    }
    /// Complete incidence miss for the finite exact carrier, without assigning
    /// an artificial global material side to its excluded quadratic branches.
    pub fn misses_point(&self, point: &Point3, dec: &mut Decider) -> Result<bool, Refusal> {
        math::validate_point(point, dec)?;
        let (station, radius, ray) = self.section_data(point);
        Ok(!self.in_station(&station, dec)?
            || !self.in_sector(&ray, dec)?
            || math::sign(&(&ray.norm_squared() - &radius.square()), dec)? != Sign::Zero)
    }
    /// Radial side on the finite station interval and retained contact sector.
    pub fn side_at(&self, point: &Point3, dec: &mut Decider) -> Result<Sign, Refusal> {
        math::validate_point(point, dec)?;
        let (station, radius, ray) = self.section_data(point);
        self.branch_side(&station, &radius, &ray, dec)
    }
    /// Whether an incident point lies on a contact rail or an endpoint section.
    pub fn on_boundary(&self, point: &Point3, dec: &mut Decider) -> Result<bool, Refusal> {
        math::validate_point(point, dec)?;
        let (station, radius, ray) = self.section_data(point);
        if self.branch_side(&station, &radius, &ray, dec)? != Sign::Zero {
            return Err(domain("boundary point is not on the exact linear fillet"));
        }
        if math::sign(&station, dec)? == Sign::Zero
            || math::sign(&(&Real::one() - &station), dec)? == Sign::Zero
        {
            return Ok(true);
        }
        for contact in self.rays() {
            if math::sign(&contact.cross(&ray).norm_squared(), dec)? == Sign::Zero {
                return Ok(true);
            }
        }
        Ok(false)
    }
    /// Natural chart normal at an actual regular point, including contact rails.
    /// An apex has no single normal and returns a typed refusal.
    pub fn normal_at(&self, point: &Point3, dec: &mut Decider) -> Result<Vector3, Refusal> {
        math::validate_point(point, dec)?;
        let (station, radius, ray) = self.section_data(point);
        if self.branch_side(&station, &radius, &ray, dec)? != Sign::Zero {
            return Err(domain("normal point is not on the exact linear fillet"));
        }
        let rate = &self.radii()[1] - &self.radii()[0];
        let axial = &((&ray.dot(self.shift()) + &radius) * rate) / &self.data.length_sq;
        let normal = ray.minus(&self.along().scaled(&axial));
        if math::sign(&normal.norm_squared(), dec)? == Sign::Zero {
            return Err(Refusal::new(
                "shapes.singular-normal",
                "a zero-radius fillet endpoint has no unique normal",
            ));
        }
        Ok(normal.scaled(&Real::from_i64(self.orientation())))
    }
    /// Complete quadratic line intersections, retaining the finite contact sector.
    /// A contained generator refuses rather than being reported as a miss.
    pub fn line_hits(
        &self,
        origin: &Point3,
        direction: &Vector3,
        dec: &mut Decider,
    ) -> Result<Vec<Point3>, Refusal> {
        math::validate_point(origin, dec)?;
        math::validate_vector(direction, dec)?;
        math::positive(&direction.norm_squared(), "shapes.zero-direction", dec)?;
        let u = self.station(origin);
        let du = &direction.dot(self.along()) / &self.data.length_sq;
        let r = self.radius(&u);
        let dr = &(&self.radii()[1] - &self.radii()[0]) * &du;
        let q = origin.minus(&self.centre(&u));
        let dq = direction
            .minus(&self.along().scaled(&du))
            .minus(&self.shift().scaled(&dr));
        let t = RealPoly::var(1, 0);
        let c = |x: Real| RealPoly::constant(1, x);
        let base = [&q.x, &q.y, &q.z];
        let rate = [&dq.x, &dq.y, &dq.z];
        let components: [RealPoly; 3] =
            std::array::from_fn(|i| c(base[i].clone()) + t.scale(rate[i]));
        let radius = c(r) + t.scale(&dr);
        let xy = &(&components[0] * &components[0]) + &(&components[1] * &components[1]);
        let xyz = &xy + &(&components[2] * &components[2]);
        let equation = &xyz - &(&radius * &radius);
        if equation.is_zero(dec)? {
            return Err(bad("line is contained in the linear fillet carrier"));
        }
        let mut hits = Vec::new();
        for t in equation.real_roots(dec)? {
            let point = origin.plus(&direction.scaled(&t));
            let station = self.station(&point);
            if self.in_station(&station, dec)?
                && self.in_sector(&point.minus(&self.centre(&station)), dec)?
            {
                hits.push(point);
            }
        }
        Ok(hits)
    }
    /// Complete parameter enclosure for a spatial box. Apex boxes retain every
    /// arc parameter; that ambiguity is geometric, not a selected inverse branch.
    pub fn parameters_of_box(
        &self,
        point: &[Ball; 3],
        dec: &Decider,
    ) -> Result<Option<[Ball; 2]>, Refusal> {
        let bits = dec.policy().max_precision_bits.min(128);
        if bits == 0 {
            return Err(Refusal::new(
                "numbers.precision-limit",
                "no linear fillet query precision",
            ));
        }
        let origin = math::point(&self.edge()[0], bits, dec)?;
        let along = [
            dec.enclose_at(&self.along().x, bits)?,
            dec.enclose_at(&self.along().y, bits)?,
            dec.enclose_at(&self.along().z, bits)?,
        ];
        let delta: [Ball; 3] = std::array::from_fn(|i| &point[i] - &origin[i]);
        let numerator = (0..3).fold(Ball::from_rational(&Rational::zero(), bits), |sum, i| {
            &sum + &(&delta[i] * &along[i])
        });
        let u = numerator
            .checked_div(&dec.enclose_at(&self.data.length_sq, bits)?)
            .ok_or_else(|| domain("edge length lost its positive separation"))?;
        let lo = u.lo().max(Rational::zero());
        let hi = u.hi().min(Rational::one());
        if lo > hi {
            return Ok(None);
        }
        let unit = Ball::from_rational(&Rational::zero(), bits)
            .hull(&Ball::from_rational(&Rational::one(), bits));
        Ok(Some([
            unit,
            Ball::from_rational(&lo, bits).hull(&Ball::from_rational(&hi, bits)),
        ]))
    }
    /// Exact extrema candidates of the complete finite band. Each coordinate is
    /// affine in station, so only the two section arcs can contain its extrema.
    pub fn axis_extrema(&self, dec: &mut Decider) -> Result<Vec<Point3>, Refusal> {
        let axis = unit(self.along());
        let mut result = Vec::with_capacity(16);
        for end in 0..2 {
            if self.zero_end(end) {
                result.push(self.edge()[end].clone());
                continue;
            }
            result.extend([self.foot(0, end).clone(), self.foot(1, end).clone()]);
            for coordinate in 0..3 {
                let direction = Vector3::new(
                    Real::from_i64(i64::from(coordinate == 0)),
                    Real::from_i64(i64::from(coordinate == 1)),
                    Real::from_i64(i64::from(coordinate == 2)),
                );
                let radial = direction.minus(&axis.scaled(&direction.dot(&axis)));
                if math::sign(&radial.norm_squared(), dec)? == Sign::Zero {
                    continue;
                }
                let radial = unit(&radial);
                for sign in [1, -1] {
                    let ray = radial.scaled(&Real::from_i64(sign));
                    if self.in_sector(&ray, dec)? {
                        result.push(
                            self.edge()[end]
                                .plus(&self.shift().plus(&ray).scaled(&self.radii()[end])),
                        );
                    }
                }
            }
        }
        Ok(result)
    }
    /// Whole-cell chord proof for a shared display grid, station-major. The
    /// circular numerator has degree two, denominator at least one half, and
    /// derivative bounds `|h'| <= 12 H`, `|h''| <= 120 H`. Taylor interpolation
    /// bounds both complete triangles, including their mixed derivative.
    pub fn display_grid(
        &self,
        tolerance: &Rational,
        maximum: usize,
        dec: &Decider,
    ) -> Result<EnvelopeGrid, Refusal> {
        if *tolerance <= Rational::zero() || maximum < 4 {
            return Err(Refusal::new(
                "mesh.invalid-options",
                "positive tolerance and at least four grid points required",
            ));
        }
        let bits = dec.policy().max_precision_bits.min(128);
        if bits == 0 {
            return Err(Refusal::new(
                "numbers.precision-limit",
                "no linear fillet display precision",
            ));
        }
        let mut h = Rational::zero();
        for ray in self.rays() {
            h = h.max(
                dec.enclose_at(
                    &math::exact_sqrt(&self.shift().plus(ray).norm_squared()),
                    bits,
                )?
                .hi(),
            );
        }
        let radius = dec
            .enclose_at(&self.radii()[0], bits)?
            .hi()
            .max(dec.enclose_at(&self.radii()[1], bits)?.hi());
        let change = dec.enclose_at(&(&self.radii()[1] - &self.radii()[0]), bits)?;
        let uu = &(&Rational::from_i64(120) * &h) * &radius;
        let uv = &(&Rational::from_i64(12) * &h) * &change.lo().abs().max(change.hi().abs());
        let (mut arcs, mut stations) = (1usize, 1usize);
        let chord = loop {
            dec.check_deadline()?;
            if arcs > i64::MAX as usize
                || stations > i64::MAX as usize
                || arcs
                    .checked_add(1)
                    .and_then(|n| stations.checked_add(1).and_then(|m| n.checked_mul(m)))
                    .is_none_or(|n| n > maximum)
            {
                return Err(Refusal::new(
                    "mesh.resource-limit",
                    "linear fillet grid exceeds the requested vertex limit",
                ));
            }
            let du = Rational::from_ratio(1, arcs as i64).unwrap();
            let dv = Rational::from_ratio(1, stations as i64).unwrap();
            let first = (&(&uu * &du) * &du).mul_pow2(-3);
            let mixed = (&(&uv * &du) * &dv).mul_pow2(-2);
            if &first + &mixed <= *tolerance {
                break first + mixed;
            }
            if first > tolerance.mul_pow2(-1) {
                arcs *= 2;
            } else {
                stations *= 2;
            }
        };
        let mut points = Vec::with_capacity((arcs + 1) * (stations + 1));
        for i in 0..=stations {
            for j in 0..=arcs {
                dec.check_deadline()?;
                let p = self.point(
                    &Rational::from_ratio(j as i64, arcs as i64).unwrap(),
                    &Rational::from_ratio(i as i64, stations as i64).unwrap(),
                )?;
                points.push(math::point(&p, bits, dec)?);
            }
        }
        Ok(EnvelopeGrid {
            columns: stations + 1,
            rows: arcs + 1,
            points,
            chord,
        })
    }
    /// Exact transport of the retained law and carrier.
    pub fn transformed(&self, pose: &Pose) -> Self {
        self.mapped(&Isometry::from(pose))
    }
    pub(crate) fn mapped(&self, map: &Isometry) -> Self {
        let d = &self.data;
        let parity = if map.reversing() { -1 } else { 1 };
        Self {
            data: Arc::new(Data {
                supports: d.supports.each_ref().map(|p| Plane {
                    origin: map.apply_point(&p.origin),
                    normal: map.apply_vector(&p.normal),
                    ref_dir: map.apply_vector(&p.ref_dir),
                }),
                same_sense: d.same_sense,
                edge: d.edge.each_ref().map(|p| map.apply_point(p)),
                radii: d.radii.clone(),
                zero: d.zero,
                convex: d.convex,
                along: map.apply_vector(&d.along),
                length_sq: d.length_sq.clone(),
                shift: map.apply_vector(&d.shift),
                rays: d.rays.each_ref().map(|v| map.apply_vector(v)),
                turn: map.apply_vector(&d.turn).scaled(&Real::from_i64(parity)),
                weight: d.weight.clone(),
                orientation: d.orientation * parity,
                carrier: d.carrier.mapped(map),
            }),
        }
    }
    /// Positive retained-definition identity for geometric query caches.
    pub fn shares_definition_nodes(&self, other: &Self) -> bool {
        Arc::ptr_eq(&self.data, &other.data)
    }
    /// Outward enclosures of all represented controls and positive weights.
    pub fn enclosed_controls(&self, dec: &Decider) -> Result<Vec<([Ball; 3], Ball)>, Refusal> {
        let bits = dec.policy().max_precision_bits.min(128);
        if bits == 0 {
            return Err(Refusal::new(
                "numbers.precision-limit",
                "no linear fillet control precision",
            ));
        }
        self.carrier()
            .controls()
            .iter()
            .zip(self.carrier().weights())
            .map(|(p, w)| Ok((math::point(p, bits, dec)?, dec.enclose_at(w, bits)?)))
            .collect()
    }
}

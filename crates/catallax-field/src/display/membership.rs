//! Rational trim predicates for display witnesses. Unsupported profiles keep the
//! shared exact-real path. Orthogonal quarter-arc sections and complete ring
//! fillets expose their certified contact sectors to the display evaluator;
//! these caches never replace or widen the feature definition.
use crate::{Location, geometry, profile};
use catallax_numbers::{Interval, Rational, Sign};
use catallax_recipe::{ElementId, Role, Turn};
use catallax_shapes::{Isometry, RVector3, Vector3};
use std::collections::BTreeMap;

type P = [Rational; 2];
fn q(n: i64) -> Rational {
    Rational::from_i64(n)
}
fn cross(a: &P, b: &P) -> Rational {
    &a[0] * &b[1] - &a[1] * &b[0]
}
fn sub(a: &P, b: &P) -> P {
    [&a[0] - &b[0], &a[1] - &b[1]]
}
fn on_segment(p: &P, edge: &profile::Edge) -> bool {
    cross(&sub(p, &edge.start), &sub(&edge.end, &edge.start)).is_zero()
        && (0..2).all(|i| {
            p[i] >= *(&edge.start[i]).min(&edge.end[i])
                && p[i] <= *(&edge.start[i]).max(&edge.end[i])
        })
}
fn circle_value(p: &P, centre: &P, squared: &Rational) -> Rational {
    let d = sub(p, centre);
    &(&d[0] * &d[0]) + &(&d[1] * &d[1]) - squared.clone()
}
// Positive normalization preserves exact signs. Axis-aligned tests need no
// per-point multiplication, and all coefficients belong to the immutable face.
struct Linear2(P);
impl Linear2 {
    fn new(coefficients: P) -> Self {
        let scale = if coefficients[0].is_zero() {
            coefficients[1].abs()
        } else {
            coefficients[0].abs()
        };
        Self(coefficients.map(|q| q / scale.clone()))
    }
    fn sign(&self, p: &P) -> Sign {
        let oriented = |value: Sign, positive: bool| {
            if positive {
                value
            } else {
                match value {
                    Sign::Positive => Sign::Negative,
                    Sign::Negative => Sign::Positive,
                    Sign::Zero => Sign::Zero,
                }
            }
        };
        if self.0[0].is_zero() {
            return oriented(p[1].sign(), self.0[1].sign() == Sign::Positive);
        }
        if self.0[1].is_zero() {
            return oriented(p[0].sign(), self.0[0].sign() == Sign::Positive);
        }
        let term = |coefficient: &Rational, value: &Rational| {
            if *coefficient == q(1) {
                value.clone()
            } else if *coefficient == q(-1) {
                -value
            } else {
                coefficient * value
            }
        };
        (term(&self.0[0], &p[0]) + term(&self.0[1], &p[1])).sign()
    }
}
#[derive(Clone, Copy)]
enum SectorSide {
    Arc,
    Radial(usize),
}
struct Sector<'a> {
    centre: &'a P,
    squared: &'a Rational,
    major: bool,
    angle: [Linear2; 2],
    heading: [Linear2; 2],
    bounds: Option<[[Rational; 2]; 2]>,
    sides: BTreeMap<ElementId, SectorSide>,
}
impl<'a> Sector<'a> {
    fn new(edges: &'a [profile::Edge]) -> Option<Self> {
        if edges.len() != 3 || edges.iter().filter(|e| e.circle.is_some()).count() != 1 {
            return None;
        }
        let index = edges.iter().position(|e| e.circle.is_some())?;
        let arc = &edges[index];
        let (centre, squared, turn) = arc.circle.as_ref()?;
        let previous = &edges[(index + 2) % 3];
        let next = &edges[(index + 1) % 3];
        if arc.full
            || previous.start != *centre
            || next.end != *centre
            || previous.end != arc.start
            || next.start != arc.end
        {
            return None;
        }
        let forward = *turn == Turn::Ccw;
        let (start, end) = if forward {
            (&arc.start, &arc.end)
        } else {
            (&arc.end, &arc.start)
        };
        let start = sub(start, centre);
        let end = sub(end, centre);
        let mut sides = BTreeMap::new();
        for id in &arc.ids {
            sides.insert(id.clone(), SectorSide::Arc);
        }
        for id in &previous.ids {
            sides.insert(id.clone(), SectorSide::Radial(usize::from(!forward)));
        }
        for id in &next.ids {
            sides.insert(id.clone(), SectorSide::Radial(usize::from(forward)));
        }
        let mut result = Self {
            centre,
            squared,
            major: cross(&start, &end) < q(0),
            angle: [
                Linear2::new([-&start[1], start[0].clone()]),
                Linear2::new([end[1].clone(), -&end[0]]),
            ],
            heading: [Linear2::new(start), Linear2::new(end)],
            bounds: None,
            sides,
        };
        // Coordinate extrema occur at the centre, arc ends, or a cardinal
        // circle point in the angular domain. Irrational radii retain the
        // exact predicate without this optional rational bound.
        if let Some(radius) = profile::sqrt(squared) {
            let mut bounds = centre.each_ref().map(|c| [c.clone(), c.clone()]);
            let mut points = vec![arc.start.clone(), arc.end.clone()];
            for direction in [
                [radius.clone(), q(0)],
                [q(0), radius.clone()],
                [-&radius, q(0)],
                [q(0), -&radius],
            ] {
                if result.angle_contains(&direction) {
                    points.push([&centre[0] + &direction[0], &centre[1] + &direction[1]]);
                }
            }
            for point in points {
                for i in 0..2 {
                    bounds[i][0] = bounds[i][0].clone().min(point[i].clone());
                    bounds[i][1] = bounds[i][1].clone().max(point[i].clone());
                }
            }
            result.bounds = Some(bounds);
        }
        Some(result)
    }
    fn bounded_outside(&self, p: &P) -> bool {
        self.bounds
            .as_ref()
            .is_some_and(|b| (0..2).any(|i| p[i] < b[i][0] || p[i] > b[i][1]))
    }
    fn relative(&self, p: &P) -> P {
        std::array::from_fn(|i| {
            if self.centre[i].is_zero() {
                p[i].clone()
            } else {
                &p[i] - &self.centre[i]
            }
        })
    }
    fn accepted(&self, signs: [Sign; 2]) -> bool {
        if self.major {
            signs[0] != Sign::Negative || signs[1] != Sign::Negative
        } else {
            signs[0] != Sign::Negative && signs[1] != Sign::Negative
        }
    }
    fn angle_contains(&self, v: &P) -> bool {
        self.accepted(self.angle.each_ref().map(|test| test.sign(v)))
    }
    fn radial_value(&self, v: &P) -> Rational {
        &(&v[0] * &v[0] + &v[1] * &v[1]) - self.squared
    }
    fn location(&self, p: &P) -> Location {
        if self.bounded_outside(p) {
            return Location::Outside;
        }
        let v = self.relative(p);
        let signs = self.angle.each_ref().map(|test| test.sign(&v));
        if !self.accepted(signs) {
            return Location::Outside;
        }
        let value = self.radial_value(&v);
        if value > q(0) {
            return Location::Outside;
        }
        if value.is_zero()
            || (0..2).any(|i| signs[i] == Sign::Zero && self.heading[i].sign(&v) != Sign::Negative)
        {
            Location::Boundary
        } else {
            Location::Inside
        }
    }
    fn on_face(&self, id: &ElementId, p: &P) -> bool {
        if self.bounded_outside(p) {
            return false;
        }
        let Some(side) = self.sides.get(id) else {
            return false;
        };
        let v = self.relative(p);
        match side {
            SectorSide::Arc => self.radial_value(&v).is_zero() && self.angle_contains(&v),
            SectorSide::Radial(i) => {
                self.angle[*i].sign(&v) == Sign::Zero
                    && self.heading[*i].sign(&v) != Sign::Negative
                    && self.radial_value(&v) <= q(0)
            }
        }
    }
}

// These maps are rebuilt with the display classifier after every definition or
// placement edit; no point result is cached. A sparse identity projection is
// therefore exact and needs no per-cell zero arithmetic or division.
struct Affine3 {
    constant: Rational,
    terms: Vec<(usize, Rational)>,
}
impl Affine3 {
    fn new(coefficients: [Rational; 3], constant: Rational) -> Self {
        Self {
            constant,
            terms: coefficients
                .into_iter()
                .enumerate()
                .filter(|(_, c)| !c.is_zero())
                .collect(),
        }
    }
    fn evaluate(&self, p: &[Rational; 3]) -> Rational {
        let mut sum = if self.constant.is_zero() {
            None
        } else {
            Some(self.constant.clone())
        };
        for (i, c) in &self.terms {
            let value = if *c == q(1) {
                p[*i].clone()
            } else if *c == q(-1) {
                -&p[*i]
            } else {
                c * &p[*i]
            };
            sum = Some(if let Some(old) = sum {
                old + value
            } else {
                value
            });
        }
        sum.unwrap_or_else(|| q(0))
    }
}
struct SectorProjection {
    height: Affine3,
    point: [Affine3; 2],
}
impl SectorProjection {
    fn new(matrix: &[[Rational; 3]; 3], translation: &[Rational; 3], vector: &RVector3) -> Self {
        let factors = [&vector.x / &vector.z, &vector.y / &vector.z];
        Self {
            height: Affine3::new(matrix[2].clone(), translation[2].clone()),
            point: std::array::from_fn(|i| {
                Affine3::new(
                    std::array::from_fn(|j| &matrix[i][j] - &(&factors[i] * &matrix[2][j])),
                    &translation[i] - &(&factors[i] * &translation[2]),
                )
            }),
        }
    }
    fn project(&self, p: &[Rational; 3]) -> (Rational, P) {
        (
            self.height.evaluate(p),
            self.point.each_ref().map(|row| row.evaluate(p)),
        )
    }
}
fn profile_location(loops: &[Vec<profile::Edge>], p: &P) -> Location {
    let mut inside = false;
    for edges in loops {
        if let Some((centre, squared, _)) = &edges[0].circle {
            let value = circle_value(p, centre, squared);
            if value.is_zero() {
                return Location::Boundary;
            }
            inside ^= value < q(0);
            continue;
        }
        for edge in edges {
            if on_segment(p, edge) {
                return Location::Boundary;
            }
            if (edge.start[1] > p[1]) != (edge.end[1] > p[1]) {
                let sign = cross(&sub(&edge.end, &edge.start), &sub(p, &edge.start));
                inside ^= (sign > q(0)) == (edge.end[1] > edge.start[1]);
            }
        }
    }
    if inside {
        Location::Inside
    } else {
        Location::Outside
    }
}

// The section between two perpendicular contact lines and their quarter arc.
// Coordinates along the two lines give a unit square; the removed corner is
// exactly its portion outside the rolling circle centred at (1,1).
struct Corner {
    origin: P,
    axes: [P; 2],
    squared: Rational,
    centre: P,
}
impl Corner {
    fn new(loops: &[Vec<profile::Edge>]) -> Option<Self> {
        let [edges] = loops else {
            return None;
        };
        if edges.len() != 3 {
            return None;
        }
        let arc = edges.iter().find(|e| e.circle.is_some())?;
        if arc.full || edges.iter().filter(|e| e.circle.is_some()).count() != 1 {
            return None;
        }
        let (centre, squared, _) = arc.circle.as_ref()?;
        let origin: P = std::array::from_fn(|i| &arc.start[i] + &arc.end[i] - centre[i].clone());
        let axes = [sub(&arc.start, &origin), sub(&arc.end, &origin)];
        let dot = |a: &P, b: &P| &a[0] * &b[0] + &a[1] * &b[1];
        if squared <= &q(0)
            || !dot(&axes[0], &axes[1]).is_zero()
            || axes.iter().any(|a| dot(a, a) != *squared)
            || !edges.iter().filter(|e| e.circle.is_none()).all(|e| {
                (e.start == origin && (e.end == arc.start || e.end == arc.end))
                    || (e.end == origin && (e.start == arc.start || e.start == arc.end))
            })
        {
            return None;
        }
        // The profile validator fixes the shorter quarter arc. Verify its turn
        // here as well, so a complementary major arc cannot use this shortcut.
        let a = sub(&arc.start, centre);
        let b = sub(&arc.end, centre);
        if (cross(&a, &b) > q(0)) != (arc.circle.as_ref()?.2 == catallax_recipe::Turn::Ccw) {
            return None;
        }
        Some(Self {
            origin,
            axes,
            squared: squared.clone(),
            centre: centre.clone(),
        })
    }
    fn location(&self, p: &P) -> Location {
        let d = sub(p, &self.origin);
        let coordinates = self.axes.each_ref().map(|a| &d[0] * &a[0] + &d[1] * &a[1]);
        if coordinates.iter().any(|v| v < &q(0) || v > &self.squared) {
            return Location::Outside;
        }
        let circle = circle_value(p, &self.centre, &self.squared);
        if circle < q(0) {
            Location::Outside
        } else if circle.is_zero() || coordinates.iter().any(|v| v.is_zero()) {
            Location::Boundary
        } else {
            Location::Inside
        }
    }
}

// A complete elementary ring fillet is the quarter-section between its edge
// and two contact circles. For a ring torus (major > minor), the other factor
// of the squared implicit equation is positive everywhere: its polynomial has
// exactly the tube-distance sign, without taking a radial square root.
pub(crate) struct Ring {
    centre: [Rational; 3],
    axis: [Rational; 3],
    height: [Rational; 2],
    radial_squared: [Rational; 2],
    major_squared: Rational,
    minor_squared: Rational,
    polynomial: catallax_numbers::MPoly,
}
impl Ring {
    fn new(definition: &catallax_features::fillet::FilletDefinition) -> Option<Self> {
        use catallax_shapes::{Curve, Surface};
        let (Surface::Torus(torus), Curve::Circle(edge), [Curve::Circle(a), Curve::Circle(b)]) =
            (&definition.surface, &definition.edge, &definition.contacts)
        else {
            return None;
        };
        let major = torus.major.as_rational()?;
        let minor = torus.minor.as_rational()?;
        if major <= minor || minor <= q(0) {
            return None;
        }
        let raw_axis = torus.axis.rational_components()?.map(Clone::clone);
        let norm = crate::profile::sqrt(&raw_axis.iter().fold(q(0), |s, x| s + x * x))?;
        if norm.is_zero() {
            return None;
        }
        let axis = raw_axis.map(|x| x / norm.clone());
        let centre = [
            torus.center.x.as_rational()?,
            torus.center.y.as_rational()?,
            torus.center.z.as_rational()?,
        ];
        let offset = [
            edge.center.x.as_rational()?,
            edge.center.y.as_rational()?,
            edge.center.z.as_rational()?,
        ];
        let height = (0..3).fold(q(0), |s, i| s + (&offset[i] - &centre[i]) * axis[i].clone());
        if &height * &height != &minor * &minor {
            return None;
        }
        let radii = [a.radius_sq.as_rational()?, b.radius_sq.as_rational()?];
        let radial_squared = [
            radii[0].clone().min(radii[1].clone()),
            radii[0].clone().max(radii[1].clone()),
        ];
        Some(Self {
            centre,
            axis,
            height: [q(0).min(height.clone()), q(0).max(height)],
            radial_squared,
            major_squared: &major * &major,
            minor_squared: &minor * &minor,
            polynomial: definition.surface.implicit()?.rational()?.clone(),
        })
    }
    fn coordinates(&self, p: &[Rational; 3]) -> (Rational, Rational) {
        let d: [Rational; 3] = std::array::from_fn(|i| &p[i] - &self.centre[i]);
        let height = (0..3).fold(q(0), |s, i| s + &d[i] * &self.axis[i]);
        let radial = d.iter().fold(q(0), |s, x| s + x * x) - &height * &height;
        (height, radial)
    }
    fn query(&self, p: &[Rational; 3]) -> (Location, bool) {
        let (height, radial) = self.coordinates(p);
        if height < self.height[0]
            || height > self.height[1]
            || radial < self.radial_squared[0]
            || radial > self.radial_squared[1]
        {
            return (Location::Outside, false);
        }
        // Factored physical tube equation. Normalizing the source axis
        // preserves the expanded polynomial's sign and zero set.
        let s = &(&radial + &(&height * &height)) + &(&self.major_squared - &self.minor_squared);
        let tube = &s * &s - (&self.major_squared * &radial).mul_pow2(2);
        if tube < q(0) {
            (Location::Outside, false)
        } else if tube.is_zero() {
            (Location::Boundary, true)
        } else if self.height.contains(&height) || self.radial_squared.contains(&radial) {
            (Location::Boundary, false)
        } else {
            (Location::Inside, false)
        }
    }
    fn location(&self, p: &[Rational; 3]) -> Location {
        self.query(p).0
    }
    fn ray(&self, p: &[Rational; 3], direction: &Vector3) -> Option<Location> {
        let direction = direction.rational_components()?;
        let d: [Rational; 3] = std::array::from_fn(|i| &p[i] - &self.centre[i]);
        let (height, radial) = self.coordinates(p);
        let dh = (0..3).fold(q(0), |s, i| s + &self.axis[i] * direction[i]);
        let dot = (0..3).fold(q(0), |s, i| s + &d[i] * direction[i]);
        let norm = (0..3).fold(q(0), |s, i| s + direction[i] * direction[i]);
        let dr = (&dot - &(&height * &dh)).mul_pow2(1);
        let ddr = &norm - &(&dh * &dh);
        let inequalities = [
            vec![&self.height[0] - &height, -&dh],
            vec![&height - &self.height[1], dh],
            vec![&self.radial_squared[0] - &radial, -&dr, -&ddr],
            vec![&radial - &self.radial_squared[1], dr.clone(), ddr.clone()],
        ];
        let mut boundary = false;
        for coefficients in inequalities {
            match coefficients.iter().find(|c| !c.is_zero()).map(|c| c.sign()) {
                Some(catallax_numbers::Sign::Positive) => return Some(Location::Outside),
                None => boundary = true,
                _ => (),
            }
        }
        let s = &(&radial + &(&height * &height)) + &(&self.major_squared - &self.minor_squared);
        let tube = &s * &s - (&self.major_squared * &radial).mul_pow2(2);
        let first = if !tube.is_zero() {
            tube
        } else {
            let l = dot.mul_pow2(1);
            let first = (&s * &l).mul_pow2(1) - (&self.major_squared * &dr).mul_pow2(2);
            if !first.is_zero() {
                first
            } else {
                let second =
                    &l * &l + (&s * &norm).mul_pow2(1) - (&self.major_squared * &ddr).mul_pow2(2);
                if !second.is_zero() {
                    second
                } else {
                    let third = (&l * &norm).mul_pow2(1);
                    if !third.is_zero() {
                        third
                    } else {
                        &norm * &norm
                    }
                }
            }
        };
        // Material lies outside the rolling tube inside the contact sector.
        if first < q(0) {
            Some(Location::Outside)
        } else if first.is_zero() || boundary {
            Some(Location::Boundary)
        } else {
            Some(Location::Inside)
        }
    }

    fn excludes(&self, cover: &[Interval; 3]) -> bool {
        let d: [Interval; 3] =
            std::array::from_fn(|i| cover[i] - Interval::from_rational(&self.centre[i]));
        let height = (0..3).fold(Interval::point(0.0).unwrap(), |s, i| {
            s + d[i] * Interval::from_rational(&self.axis[i])
        });
        if height.hi() < Interval::from_rational(&self.height[0]).lo()
            || height.lo() > Interval::from_rational(&self.height[1]).hi()
        {
            return true;
        }
        let radial = d
            .iter()
            .fold(Interval::point(0.0).unwrap(), |s, x| s + x.square())
            - height.square();
        radial.hi() < Interval::from_rational(&self.radial_squared[0]).lo()
            || radial.lo() > Interval::from_rational(&self.radial_squared[1]).hi()
    }
}

pub(crate) enum Classifier<'a> {
    Extrusion(Box<Extrusion<'a>>),
    Ring(Box<Ring>),
}
impl<'a> Classifier<'a> {
    pub(crate) fn new(
        inverse: &Isometry,
        vector: &'a RVector3,
        through: bool,
        loops: &'a [Vec<profile::Edge>],
    ) -> Option<Self> {
        Extrusion::new(inverse, vector, through, loops).map(|e| Self::Extrusion(Box::new(e)))
    }
    pub(crate) fn fillet(definition: &catallax_features::fillet::FilletDefinition) -> Option<Self> {
        Ring::new(definition).map(|r| Self::Ring(Box::new(r)))
    }
    /// Exact closed intersections for the elementary fillet tool only.
    pub(crate) fn inequalities(&self) -> Option<Vec<catallax_numbers::MPoly>> {
        use catallax_numbers::MPoly;
        let constant = |q: Rational| MPoly::constant(3, q);
        let coordinates: [MPoly; 3] = std::array::from_fn(|i| MPoly::var(3, i));
        match self {
            Self::Extrusion(e) => {
                let corner = e.corner.as_ref()?;
                let local: [MPoly; 3] = std::array::from_fn(|i| {
                    (0..3).fold(constant(e.translation[i].clone()), |s, j| {
                        &s + &coordinates[j].scale(&e.matrix[i][j])
                    })
                });
                let t = local[2].scale(&(q(1) / e.vector.z.clone()));
                let p: [MPoly; 2] = std::array::from_fn(|i| {
                    &local[i] - &t.scale(if i == 0 { &e.vector.x } else { &e.vector.y })
                });
                let mut result = Vec::new();
                if !e.through {
                    result.extend([t.scale(&q(-1)), &t - &constant(q(1))]);
                }
                for axis in &corner.axes {
                    let u = (0..2).fold(MPoly::zero(3), |s, i| {
                        &s + &(&p[i] - &constant(corner.origin[i].clone())).scale(&axis[i])
                    });
                    result.extend([u.scale(&q(-1)), &u - &constant(corner.squared.clone())]);
                }
                let circle = (0..2).fold(constant(-&corner.squared), |s, i| {
                    let d = &p[i] - &constant(corner.centre[i].clone());
                    &s + &(&d * &d)
                });
                result.push(circle.scale(&q(-1)));
                Some(result)
            }
            Self::Ring(r) => {
                let d: [MPoly; 3] =
                    std::array::from_fn(|i| &coordinates[i] - &constant(r.centre[i].clone()));
                let h = (0..3).fold(MPoly::zero(3), |s, i| &s + &d[i].scale(&r.axis[i]));
                let radial = &d.iter().fold(MPoly::zero(3), |s, d| &s + &(d * d)) - &(&h * &h);
                Some(vec![
                    &constant(r.height[0].clone()) - &h,
                    &h - &constant(r.height[1].clone()),
                    &constant(r.radial_squared[0].clone()) - &radial,
                    &radial - &constant(r.radial_squared[1].clone()),
                    r.polynomial.scale(&q(-1)),
                ])
            }
        }
    }
    pub(crate) fn location(&self, p: &[Rational; 3]) -> Location {
        match self {
            Self::Extrusion(e) => e.location(p),
            Self::Ring(r) => r.location(p),
        }
    }
    pub(crate) fn on_face(&self, role: &Role, p: &[Rational; 3]) -> bool {
        match self {
            Self::Extrusion(e) => e.on_face(role, p),
            Self::Ring(r) => r.query(p).1,
        }
    }
    pub(crate) fn excludes(&self, cover: &[Interval; 3]) -> bool {
        match self {
            Self::Extrusion(e) => e.excludes(cover),
            Self::Ring(r) => r.excludes(cover),
        }
    }
    pub(crate) fn ray(
        &self,
        constraints: &[geometry::Constraint],
        p: &[Rational; 3],
        direction: &Vector3,
    ) -> Option<Location> {
        match self {
            Self::Extrusion(e) => e.ray(constraints, p, direction),
            Self::Ring(r) => r.ray(p, direction),
        }
    }
}

pub(crate) struct Extrusion<'a> {
    matrix: [[Rational; 3]; 3],
    translation: [Rational; 3],
    vector: &'a RVector3,
    through: bool,
    loops: &'a [Vec<profile::Edge>],
    corner: Option<Corner>,
    sector: Option<Sector<'a>>,
    sector_projection: Option<SectorProjection>,
    parameter_bounds: [Rational; 2],
}
impl<'a> Extrusion<'a> {
    pub(crate) fn new(
        inverse: &Isometry,
        vector: &'a RVector3,
        through: bool,
        loops: &'a [Vec<profile::Edge>],
    ) -> Option<Self> {
        let corner = Corner::new(loops);
        // Preserve the established fillet path when its exact corner proof applies.
        let sector = if corner.is_none() {
            loops
                .first()
                .filter(|_| loops.len() == 1)
                .and_then(|edges| Sector::new(edges))
        } else {
            None
        };
        if corner.is_none()
            && sector.is_none()
            && loops.iter().any(|edges| {
                edges.is_empty()
                    || !(edges.iter().all(|e| e.circle.is_none())
                        || (edges.len() == 1 && edges[0].full && edges[0].circle.is_some()))
            })
        {
            return None;
        }
        let mut matrix = std::array::from_fn(|_| std::array::from_fn(|_| q(0)));
        for (i, row) in inverse.matrix().iter().enumerate() {
            for (j, value) in row.iter().enumerate() {
                matrix[i][j] = value.as_rational()?;
            }
        }
        let translation = [
            inverse.translation.x.as_rational()?,
            inverse.translation.y.as_rational()?,
            inverse.translation.z.as_rational()?,
        ];
        let sector_projection = sector
            .as_ref()
            .map(|_| SectorProjection::new(&matrix, &translation, vector));
        let parameter_bounds = if sector_projection.is_some() {
            [q(0).min(vector.z.clone()), q(0).max(vector.z.clone())]
        } else {
            [q(0), q(1)]
        };
        Some(Self {
            matrix,
            translation,
            vector,
            through,
            loops,
            corner,
            sector,
            sector_projection,
            parameter_bounds,
        })
    }
    fn project(&self, p: &[Rational; 3]) -> (Rational, P) {
        if let Some(projector) = &self.sector_projection {
            return projector.project(p);
        }
        let local: [Rational; 3] = std::array::from_fn(|i| {
            (0..3).fold(self.translation[i].clone(), |s, j| {
                s + &self.matrix[i][j] * &p[j]
            })
        });
        let t = &local[2] / &self.vector.z;
        (
            t.clone(),
            [
                &local[0] - &(&t * &self.vector.x),
                &local[1] - &(&t * &self.vector.y),
            ],
        )
    }
    pub(crate) fn location(&self, p: &[Rational; 3]) -> Location {
        let (t, p) = self.project(p);
        if !self.through && (t < self.parameter_bounds[0] || t > self.parameter_bounds[1]) {
            return Location::Outside;
        }
        let location = self.profile_location(&p);
        if location == Location::Inside && !self.through && self.parameter_bounds.contains(&t) {
            Location::Boundary
        } else {
            location
        }
    }
    fn profile_location(&self, p: &P) -> Location {
        if let Some(corner) = &self.corner {
            return corner.location(p);
        }
        self.sector.as_ref().map_or_else(
            || profile_location(self.loops, p),
            |sector| sector.location(p),
        )
    }
    pub(crate) fn on_face(&self, role: &Role, p: &[Rational; 3]) -> bool {
        let (t, p) = self.project(p);
        if !self.through && (t < self.parameter_bounds[0] || t > self.parameter_bounds[1]) {
            return false;
        }
        let mut role = role;
        while let Role::Copy { of, .. } = role {
            role = &of.role;
        }
        if let Some(corner) = &self.corner {
            if corner.location(&p) == Location::Outside {
                return false;
            }
            if matches!(role, Role::Blend { .. }) {
                return circle_value(&p, &corner.centre, &corner.squared).is_zero();
            }
        }
        if let (Some(sector), Role::Side { segment, .. }) = (&self.sector, role) {
            return sector.on_face(segment, &p);
        }
        match role {
            Role::Start | Role::End => self.profile_location(&p) != Location::Outside,
            Role::Side { segment, .. } => self
                .loops
                .iter()
                .flatten()
                .filter(|e| e.ids.contains(segment))
                .any(|e| {
                    if let Some((centre, radius, _)) = &e.circle {
                        circle_value(&p, centre, radius).is_zero()
                            && self
                                .sector
                                .as_ref()
                                .is_none_or(|sector| sector.angle_contains(&sub(&p, centre)))
                    } else {
                        on_segment(&p, e)
                    }
                }),
            _ => false,
        }
    }
    /// Exclude a whole cell using the contact sector, including its arc bounds.
    /// Bounds are outward interval projections, never witness samples.
    pub(crate) fn excludes(&self, cover: &[Interval; 3]) -> bool {
        let Some(corner) = &self.corner else {
            return false;
        };
        let local: [Interval; 3] = std::array::from_fn(|i| {
            (0..3).fold(Interval::from_rational(&self.translation[i]), |s, j| {
                s + Interval::from_rational(&self.matrix[i][j]) * cover[j]
            })
        });
        let t = local[2] / Interval::from_rational(&self.vector.z);
        if !self.through && (t.hi() < 0.0 || t.lo() > 1.0) {
            return true;
        }
        let p: [Interval; 2] = std::array::from_fn(|i| {
            local[i]
                - t * Interval::from_rational(if i == 0 {
                    &self.vector.x
                } else {
                    &self.vector.y
                })
                - Interval::from_rational(&corner.origin[i])
        });
        corner.axes.iter().any(|axis| {
            let coordinate =
                p[0] * Interval::from_rational(&axis[0]) + p[1] * Interval::from_rational(&axis[1]);
            coordinate.hi() < 0.0 || coordinate.lo() > Interval::from_rational(&corner.squared).hi()
        })
    }
    pub(crate) fn ray(
        &self,
        constraints: &[geometry::Constraint],
        p: &[Rational; 3],
        direction: &Vector3,
    ) -> Option<Location> {
        let d = direction.rational_components()?;
        let mut epsilon = q(1);
        for constraint in constraints {
            let coefficients = geometry::native_ray_coefficients(constraint, p, direction)?;
            if let Some((i, c)) = coefficients.iter().enumerate().find(|(_, c)| !c.is_zero()) {
                let tail = coefficients[i + 1..].iter().fold(q(0), |a, c| a + c.abs());
                epsilon = epsilon.min(c.abs() / (c.abs() + tail + q(1)));
            }
        }
        Some(self.location(&std::array::from_fn(|i| &p[i] + &(&epsilon * d[i]))))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use catallax_numbers::{Decider, Real};
    use catallax_shapes::{Point3, RPoint3};

    #[test]
    fn a_nonradial_chord_profile_keeps_the_generic_fallback() {
        let mut source: serde_json::Value = serde_json::from_str(include_str!(
            "../../../catallax-features/tests/fixtures/chamfer-arcs/chord.json"
        ))
        .unwrap();
        source["steps"].as_array_mut().unwrap().pop();
        source["result"] = "base".into();
        let plan =
            crate::Plan::compile(&catallax_recipe::Recipe::parse(&source.to_string()).unwrap())
                .unwrap();
        for leaf in &plan.leaves {
            if let geometry::Leaf::Feature(f) = leaf {
                assert!(f.display_classifier().is_none());
            }
        }
    }

    #[test]
    fn sector_cache_matches_generic_arc_winding_faces_and_germs() {
        for (name, text) in [
            (
                "convex",
                include_str!("../../../catallax-features/tests/fixtures/chamfer-arcs/convex.json"),
            ),
            (
                "major",
                include_str!("../../../catallax-features/tests/fixtures/chamfer-arcs/major.json"),
            ),
            (
                "major-seam",
                include_str!(
                    "../../../catallax-features/tests/fixtures/chamfer-arcs/major-seam.json"
                ),
            ),
        ] {
            for reverse in [false, true] {
                for reflected in [false, true] {
                    let mut source: serde_json::Value = serde_json::from_str(text).unwrap();
                    source["steps"].as_array_mut().unwrap().pop();
                    source["result"] = "base".into();
                    if reverse {
                        let edges = &source["steps"][0]["loops"][0];
                        let start = edges[0]["to"].clone();
                        let end = edges[1]["to"].clone();
                        let centre = edges[1]["centre"].clone();
                        source["steps"][0]["loops"][0] = serde_json::json!([
                            {"kind":"line","id":"radial-a","to":end},
                            {"kind":"arc","id":"wall","centre":centre,"to":start,"turn":"cw"},
                            {"kind":"line","id":"radial-b","to":centre}]);
                    }
                    let map = if reflected {
                        source["steps"].as_array_mut().unwrap().push(serde_json::json!({"id":"image","op":"mirror","body":"base","plane":{"point":["0","0","0"],"normal":["1","0","0"]},"keep_original":false}));
                        source["result"] = "image".into();
                        Isometry::reflection(
                            &RPoint3::new(q(0), q(0), q(0)),
                            &RVector3::new(q(1), q(0), q(0)),
                        )
                        .unwrap()
                    } else {
                        Isometry::identity()
                    };
                    let plan = crate::Plan::compile(
                        &catallax_recipe::Recipe::parse(&source.to_string()).unwrap(),
                    )
                    .unwrap();
                    let directions = [[1, 0, 0], [0, 0, 1], [1, 1, 1]].map(|v| {
                        map.apply_vector(&RVector3::new(q(v[0]), q(v[1]), q(v[2])).to_real())
                    });
                    let mut dec = Decider::default();
                    for (index, leaf) in plan.leaves.iter().enumerate() {
                        if plan.leaf_registers[index].is_none() {
                            continue;
                        }
                        let geometry::Leaf::Feature(feature) = leaf else {
                            continue;
                        };
                        let native = feature.display_classifier().expect("radial sector cache");
                        for xy in [
                            ["0", "0"],
                            ["8", "0"],
                            ["0", "8"],
                            ["0", "-8"],
                            ["-8", "0"],
                            ["6", "0"],
                            ["0", "6"],
                            ["3", "4"],
                            ["-3", "4"],
                            ["-3", "-4"],
                            ["3", "-4"],
                            ["9", "0"],
                            ["24/5", "32/5"],
                        ] {
                            for z in [-1, 0, 6, 12, 13] {
                                let point = map.apply_point(&Point3::new(
                                    Real::from(Rational::parse(xy[0]).unwrap()),
                                    Real::from(Rational::parse(xy[1]).unwrap()),
                                    Real::from_i64(z),
                                ));
                                let point = [
                                    point.x.as_rational().unwrap(),
                                    point.y.as_rational().unwrap(),
                                    point.z.as_rational().unwrap(),
                                ];
                                assert_eq!(
                                    native.location(&point),
                                    leaf.location(&point, &mut dec).unwrap(),
                                    "{name} reverse={reverse} reflected={reflected}"
                                );
                                for (face_index, face) in plan
                                    .faces
                                    .iter()
                                    .enumerate()
                                    .filter(|(_, f)| f.leaf == index)
                                {
                                    let role = plan.definition_role(face_index);
                                    assert_eq!(
                                        Some(native.on_face(role, &point)),
                                        leaf.on_face(&face.origin, role, &point, &mut dec).unwrap()
                                    );
                                }
                                for direction in &directions {
                                    assert_eq!(
                                        native.ray(leaf.constraints(), &point, direction).unwrap(),
                                        leaf.ray(&point, direction, &mut dec).unwrap()
                                    );
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn sector_projection_is_rebuilt_for_edits_and_rigid_placement() {
        use catallax_shapes::{Pose, Rotation};
        for radius in [8, 10] {
            for height in [12_i64, -12] {
                let mut source: serde_json::Value = serde_json::from_str(include_str!(
                    "../../../catallax-features/tests/fixtures/chamfer-arcs/convex.json"
                ))
                .unwrap();
                source["steps"].as_array_mut().unwrap().pop();
                source["result"] = "base".into();
                source["steps"][0]["loops"][0][0]["to"] =
                    serde_json::json!([radius.to_string(), "0"]);
                source["steps"][0]["loops"][0][1]["to"] =
                    serde_json::json!(["0", radius.to_string()]);
                source["steps"][1]["distance"] = height.to_string().into();
                if radius == 10 {
                    source["steps"][1]["direction"] = serde_json::json!(["3", "0", "4"]);
                }
                let original = crate::Plan::compile(
                    &catallax_recipe::Recipe::parse(&source.to_string()).unwrap(),
                )
                .unwrap();
                for quaternion in [[1, 0, 0, 0], [1, 2, 3, 4]] {
                    let [w, x, y, z] = quaternion.map(q);
                    let pose = Pose::from_rational(
                        Rotation::from_quaternion(w, x, y, z).unwrap(),
                        RVector3::new(q(7), q(-3), q(11)),
                    );
                    let image = original
                        .transformed(&pose, &mut Decider::default())
                        .unwrap();
                    for leaf in &image.leaves {
                        let geometry::Leaf::Feature(feature) = leaf else {
                            continue;
                        };
                        let native = feature.display_classifier().unwrap();
                        for [x, y, z] in [
                            [3, 4, 0],
                            [3, 4, height / 2],
                            [3, 4, height],
                            [3, 4, height + height.signum()],
                            [9, 0, height / 2],
                            [6, 8, height / 2],
                        ] {
                            let point = pose
                                .apply_rational_point(&RPoint3::new(q(x), q(y), q(z)))
                                .unwrap();
                            let point = [point.x, point.y, point.z];
                            assert_eq!(
                                native.location(&point),
                                leaf.location(&point, &mut Decider::default()).unwrap()
                            );
                        }
                    }
                }
            }
        }
    }

    #[test]
    fn native_trim_membership_and_germs_match_exact_feature_definitions() {
        let points = [
            [0, 0, 0],
            [60, 40, 5],
            [8, 12, 0],
            [10, 12, 5],
            [10, 12, 2],
            [30, 25, 5],
            [34, 25, 10],
            [38, 25, 15],
            [30, 25, 16],
            [30, 0, 18],
            [30, 5, 18],
            [30, 2, 18],
            [30, 2, 21],
            [30, 2, 15],
        ];
        let directions = [
            RVector3::new(q(1), q(0), q(0)).to_real(),
            RVector3::new(q(0), q(0), q(1)).to_real(),
            RVector3::new(q(1), q(1), q(1)).to_real(),
        ];
        for reflected in [false, true] {
            let mut source: serde_json::Value =
                serde_json::from_str(include_str!("../../../../parts/P1-bracket.json")).unwrap();
            let map = if reflected {
                let body = source["result"].clone();
                source["steps"]
                    .as_array_mut()
                    .unwrap()
                    .push(serde_json::json!({
                        "id":"display-image", "op":"mirror", "body":body,
                        "plane":{"point":["1","2","3"],"normal":["1","2","2"]},
                        "keep_original":false
                    }));
                source["result"] = "display-image".into();
                Isometry::reflection(
                    &RPoint3::new(q(1), q(2), q(3)),
                    &RVector3::new(q(1), q(2), q(2)),
                )
                .unwrap()
            } else {
                Isometry::identity()
            };
            let recipe = catallax_recipe::Recipe::parse(&source.to_string()).unwrap();
            let plan = crate::Plan::compile(&recipe).unwrap();
            let directions = directions.each_ref().map(|d| map.apply_vector(d));
            let mut dec = Decider::default();
            for (index, leaf) in plan.leaves.iter().enumerate() {
                if plan.leaf_registers[index].is_none() {
                    continue;
                }
                let geometry::Leaf::Feature(feature) = leaf else {
                    continue;
                };
                let Some(native) = feature.display_classifier() else {
                    continue;
                };
                for point in points {
                    let point = point.map(q);
                    let point = map.apply_point(&Point3::new(
                        point[0].clone().into(),
                        point[1].clone().into(),
                        point[2].clone().into(),
                    ));
                    let point = [
                        point.x.as_rational().unwrap(),
                        point.y.as_rational().unwrap(),
                        point.z.as_rational().unwrap(),
                    ];
                    assert_eq!(
                        native.location(&point),
                        leaf.location(&point, &mut dec).unwrap()
                    );
                    for (face_index, _) in plan
                        .faces
                        .iter()
                        .enumerate()
                        .filter(|(_, face)| face.leaf == index)
                    {
                        let role = plan.definition_role(face_index);
                        assert_eq!(
                            Some(native.on_face(role, &point)),
                            leaf.on_face(&plan.faces[face_index].origin, role, &point, &mut dec)
                                .unwrap(),
                        );
                    }
                    for direction in &directions {
                        assert_eq!(
                            native.ray(leaf.constraints(), &point, direction).unwrap(),
                            leaf.ray(&point, direction, &mut dec).unwrap(),
                        );
                    }
                }
            }
        }
    }
    fn independent_profile_ray(
        feature: &crate::features::Feature,
        p: &[Rational; 3],
        direction: &Vector3,
        dec: &mut catallax_numbers::Decider,
    ) -> Location {
        // Reference the original finite, sign-cell-preserving profile query,
        // independently of the new direct polynomial germ evaluator.
        let d = direction.rational_components().unwrap();
        let mut epsilon = q(1);
        for constraint in &feature.constraints {
            let coefficients = geometry::native_ray_coefficients(constraint, p, direction).unwrap();
            if let Some((i, c)) = coefficients.iter().enumerate().find(|(_, c)| !c.is_zero()) {
                let tail = coefficients[i + 1..]
                    .iter()
                    .fold(q(0), |sum, c| sum + c.abs());
                epsilon = epsilon.min(c.abs() / (c.abs() + tail + q(1)));
            }
        }
        let point = catallax_shapes::RPoint3::new(
            &p[0] + &(&epsilon * d[0]),
            &p[1] + &(&epsilon * d[1]),
            &p[2] + &(&epsilon * d[2]),
        );
        let mut reference = feature.clone();
        reference.fillet = None;
        reference.classify(&point.to_real(), dec).unwrap()
    }

    #[test]
    fn native_fillet_regions_match_independent_definitions_and_ray_germs() {
        let recipe =
            catallax_recipe::Recipe::parse(include_str!("../../../../parts/P3-filleted-boss.json"))
                .unwrap();
        let plan = crate::Plan::compile(&recipe).unwrap();
        let mut samples = Vec::new();
        plan.visit_candidates(2, |p, id| {
            samples.push((p, id));
            true
        })
        .unwrap();
        let mut dec = Decider::default();
        let mut checked = 0;
        for (p, id) in samples {
            let face = &plan.faces[id.index()];
            let geometry::Leaf::Feature(feature) = &plan.leaves[face.leaf] else {
                continue;
            };
            if feature.fillet.is_none() {
                continue;
            }
            let native = feature
                .fillet
                .as_deref()
                .and_then(Classifier::fillet)
                .or_else(|| feature.display_classifier())
                .unwrap();
            let xyz = [p.x.clone(), p.y.clone(), p.z.clone()];
            // Disable only the query cache trigger in the reference clone.
            // Its retained profile, sweep and frame remain unchanged.
            let mut reference = feature.clone();
            reference.fillet = None;
            assert_eq!(
                native.location(&xyz),
                reference.classify(&p.to_real(), &mut dec).unwrap()
            );
            let inequalities = native.inequalities().unwrap();
            let values: Vec<_> = inequalities.iter().map(|p| p.eval(&xyz)).collect();
            let enclosed = if values.iter().any(|p| p > &q(0)) {
                Location::Outside
            } else if values.iter().any(|p| p.is_zero()) {
                Location::Boundary
            } else {
                Location::Inside
            };
            assert_eq!(enclosed, native.location(&xyz));
            assert_eq!(
                Some(native.on_face(plan.definition_role(id.index()), &xyz)),
                Some(
                    reference
                        .on_face_point(
                            &Role::Side {
                                segment: catallax_recipe::ElementId::new("blend").unwrap(),
                                path: None,
                            },
                            &p.to_real(),
                            &mut dec,
                        )
                        .unwrap()
                )
            );
            let normal = crate::normal(face, &xyz);
            for normal in [
                normal.clone(),
                normal.scaled(&catallax_numbers::Real::from_i64(-1)),
            ] {
                assert_eq!(
                    native
                        .ray(plan.leaves[face.leaf].constraints(), &xyz, &normal)
                        .unwrap(),
                    independent_profile_ray(feature, &xyz, &normal, &mut dec)
                );
            }
            checked += 1;
        }
        assert!(checked >= 20);
    }
    #[test]
    fn quarter_sections_have_independent_inside_outside_and_contact_controls() {
        let recipe =
            catallax_recipe::Recipe::parse(include_str!("../../../../parts/P3-filleted-boss.json"))
                .unwrap();
        let plan = crate::Plan::compile(&recipe).unwrap();
        let ring = plan
            .leaves
            .iter()
            .find_map(|leaf| {
                let geometry::Leaf::Feature(f) = leaf else {
                    return None;
                };
                f.fillet
                    .as_deref()
                    .and_then(Classifier::fillet)
                    .map(|c| (f, c))
            })
            .unwrap();
        let mut dec = Decider::default();
        for (p, expected) in [
            (["28", "15", "10"], Location::Boundary),
            (["142/5", "15", "51/5"], Location::Inside),
            (["144/5", "15", "52/5"], Location::Boundary),
            (["147/5", "15", "51/5"], Location::Outside),
            (["31", "15", "101/10"], Location::Outside),
            (["28", "15", "13"], Location::Outside),
            (["27", "15", "21/2"], Location::Outside),
        ] {
            let p = p.map(|x| Rational::parse(x).unwrap());
            assert_eq!(ring.1.location(&p), expected);
            let point = catallax_shapes::RPoint3::new(p[0].clone(), p[1].clone(), p[2].clone());
            assert_eq!(
                ring.0.classify(&point.to_real(), &mut dec).unwrap(),
                expected
            );
        }
        let corner =
            plan.faces
                .iter()
                .find_map(|face| {
                    let Role::Blend { edge } = &face.origin.role else {
                        return None;
                    };
                    if !edge.between().iter().any(
                        |o| matches!(&o.role,Role::Side{segment,..} if segment.as_str()=="x_min"),
                    ) {
                        return None;
                    }
                    let geometry::Leaf::Feature(f) = &plan.leaves[face.leaf] else {
                        return None;
                    };
                    Some((f, f.display_classifier().unwrap()))
                })
                .unwrap();
        for (p, expected) in [
            (["1/10", "15", "99/10"], Location::Inside),
            (["3/5", "15", "49/5"], Location::Outside),
            (["2/5", "15", "49/5"], Location::Boundary),
            (["0", "15", "9"], Location::Boundary),
            (["0", "15", "899/100"], Location::Outside),
        ] {
            let p = p.map(|x| Rational::parse(x).unwrap());
            assert_eq!(corner.1.location(&p), expected);
            let point = catallax_shapes::RPoint3::new(p[0].clone(), p[1].clone(), p[2].clone());
            assert_eq!(
                corner.0.classify(&point.to_real(), &mut dec).unwrap(),
                expected
            );
        }
    }
}

//! Rational and exact-real points and vectors.

use catallax_numbers::{Rational, Real};

/// A point with exact rational coordinates.
#[derive(Clone, PartialEq, Eq, Hash, Debug)]
pub struct RPoint3 {
    /// x coordinate.
    pub x: Rational,
    /// y coordinate.
    pub y: Rational,
    /// z coordinate.
    pub z: Rational,
}

/// A vector with exact rational components. Directions are not normalised.
#[derive(Clone, PartialEq, Eq, Hash, Debug)]
pub struct RVector3 {
    /// x component.
    pub x: Rational,
    /// y component.
    pub y: Rational,
    /// z component.
    pub z: Rational,
}

/// A point with exact real coordinates.
#[derive(Clone, Debug)]
pub struct Point3 {
    /// x coordinate.
    pub x: Real,
    /// y coordinate.
    pub y: Real,
    /// z coordinate.
    pub z: Real,
}

/// A vector with exact real components.
#[derive(Clone, Debug)]
pub struct Vector3 {
    /// x component.
    pub x: Real,
    /// y component.
    pub y: Real,
    /// z component.
    pub z: Real,
}

impl RPoint3 {
    /// The point `(x, y, z)`.
    pub fn new(x: Rational, y: Rational, z: Rational) -> Self {
        Self { x, y, z }
    }

    /// The origin.
    pub fn origin() -> Self {
        Self::new(Rational::zero(), Rational::zero(), Rational::zero())
    }

    /// `self - other`.
    pub fn minus(&self, other: &RPoint3) -> RVector3 {
        RVector3::new(&self.x - &other.x, &self.y - &other.y, &self.z - &other.z)
    }

    /// `self + v`.
    pub fn plus(&self, v: &RVector3) -> RPoint3 {
        Self::new(&self.x + &v.x, &self.y + &v.y, &self.z + &v.z)
    }

    /// The same point with exact real coordinates.
    pub fn to_real(&self) -> Point3 {
        Point3::new(
            self.x.clone().into(),
            self.y.clone().into(),
            self.z.clone().into(),
        )
    }
}

impl RVector3 {
    /// The vector `(x, y, z)`.
    pub fn new(x: Rational, y: Rational, z: Rational) -> Self {
        Self { x, y, z }
    }

    /// The zero vector.
    pub fn zero() -> Self {
        Self::new(Rational::zero(), Rational::zero(), Rational::zero())
    }

    /// Whether all components are zero.
    pub fn is_zero(&self) -> bool {
        self.x.is_zero() && self.y.is_zero() && self.z.is_zero()
    }

    /// `self + other`.
    pub fn plus(&self, other: &RVector3) -> RVector3 {
        Self::new(&self.x + &other.x, &self.y + &other.y, &self.z + &other.z)
    }

    /// `self - other`.
    pub fn minus(&self, other: &RVector3) -> RVector3 {
        Self::new(&self.x - &other.x, &self.y - &other.y, &self.z - &other.z)
    }

    /// `c * self`.
    pub fn scaled(&self, c: &Rational) -> RVector3 {
        Self::new(&self.x * c, &self.y * c, &self.z * c)
    }

    /// Dot product.
    pub fn dot(&self, other: &RVector3) -> Rational {
        (&self.x * &other.x) + (&self.y * &other.y) + (&self.z * &other.z)
    }

    /// Cross product.
    pub fn cross(&self, other: &RVector3) -> RVector3 {
        Self::new(
            &self.y * &other.z - &self.z * &other.y,
            &self.z * &other.x - &self.x * &other.z,
            &self.x * &other.y - &self.y * &other.x,
        )
    }

    /// Squared length.
    pub fn norm_squared(&self) -> Rational {
        self.dot(self)
    }

    /// The same vector with exact real components.
    pub fn to_real(&self) -> Vector3 {
        Vector3::new(
            self.x.clone().into(),
            self.y.clone().into(),
            self.z.clone().into(),
        )
    }
}

impl Point3 {
    /// The point `(x, y, z)`.
    pub fn new(x: Real, y: Real, z: Real) -> Self {
        Self { x, y, z }
    }

    /// The exact origin.
    pub fn origin() -> Self {
        let zero = Real::zero();
        Self::new(zero.clone(), zero.clone(), zero)
    }

    /// Borrowed cached rational components; never refines or rounds.
    pub fn rational_components(&self) -> Option<[&Rational; 3]> {
        Some([
            self.x.rational_ref()?,
            self.y.rational_ref()?,
            self.z.rational_ref()?,
        ])
    }

    /// Cheap exact rational conversion, without refinement or rounding.
    pub fn as_rational(&self) -> Option<RPoint3> {
        let [x, y, z] = self.rational_components()?;
        Some(RPoint3::new(x.clone(), y.clone(), z.clone()))
    }

    /// `self - other`.
    pub fn minus(&self, other: &Point3) -> Vector3 {
        if let (Some([x, y, z]), Some([a, b, c])) =
            (self.rational_components(), other.rational_components())
        {
            return Vector3::new((x - a).into(), (y - b).into(), (z - c).into());
        }
        Vector3::new(&self.x - &other.x, &self.y - &other.y, &self.z - &other.z)
    }

    /// `self + v`.
    pub fn plus(&self, v: &Vector3) -> Point3 {
        if let (Some([x, y, z]), Some([a, b, c])) =
            (self.rational_components(), v.rational_components())
        {
            return Self::new((x + a).into(), (y + b).into(), (z + c).into());
        }
        Self::new(&self.x + &v.x, &self.y + &v.y, &self.z + &v.z)
    }
}

impl Vector3 {
    /// The vector `(x, y, z)`.
    pub fn new(x: Real, y: Real, z: Real) -> Self {
        Self { x, y, z }
    }

    /// The exact zero vector.
    pub fn zero() -> Self {
        let zero = Real::zero();
        Self::new(zero.clone(), zero.clone(), zero)
    }

    /// Borrowed cached rational components; never refines or rounds.
    pub fn rational_components(&self) -> Option<[&Rational; 3]> {
        Some([
            self.x.rational_ref()?,
            self.y.rational_ref()?,
            self.z.rational_ref()?,
        ])
    }

    /// Cheap exact rational conversion, without refinement or rounding.
    pub fn as_rational(&self) -> Option<RVector3> {
        let [x, y, z] = self.rational_components()?;
        Some(RVector3::new(x.clone(), y.clone(), z.clone()))
    }

    /// `self + other`.
    pub fn plus(&self, other: &Vector3) -> Vector3 {
        if let (Some([x, y, z]), Some([a, b, c])) =
            (self.rational_components(), other.rational_components())
        {
            return Self::new((x + a).into(), (y + b).into(), (z + c).into());
        }
        Self::new(&self.x + &other.x, &self.y + &other.y, &self.z + &other.z)
    }

    /// `self - other`.
    pub fn minus(&self, other: &Vector3) -> Vector3 {
        if let (Some([x, y, z]), Some([a, b, c])) =
            (self.rational_components(), other.rational_components())
        {
            return Self::new((x - a).into(), (y - b).into(), (z - c).into());
        }
        Self::new(&self.x - &other.x, &self.y - &other.y, &self.z - &other.z)
    }

    /// `c * self`.
    pub fn scaled(&self, c: &Real) -> Vector3 {
        if let (Some([x, y, z]), Some(c)) = (self.rational_components(), c.rational_ref()) {
            return Self::new((x * c).into(), (y * c).into(), (z * c).into());
        }
        Self::new(&self.x * c, &self.y * c, &self.z * c)
    }

    /// Dot product.
    pub fn dot(&self, other: &Vector3) -> Real {
        if let (Some([x, y, z]), Some([a, b, c])) =
            (self.rational_components(), other.rational_components())
        {
            return (x * a + y * b + z * c).into();
        }
        (&self.x * &other.x) + (&self.y * &other.y) + (&self.z * &other.z)
    }

    /// Cross product.
    pub fn cross(&self, other: &Vector3) -> Vector3 {
        if let (Some([x, y, z]), Some([a, b, c])) =
            (self.rational_components(), other.rational_components())
        {
            return Self::new(
                (y * c - z * b).into(),
                (z * a - x * c).into(),
                (x * b - y * a).into(),
            );
        }
        Self::new(
            &self.y * &other.z - &self.z * &other.y,
            &self.z * &other.x - &self.x * &other.z,
            &self.x * &other.y - &self.y * &other.x,
        )
    }

    /// Squared length.
    pub fn norm_squared(&self) -> Real {
        if let Some([x, y, z]) = self.rational_components() {
            return (x * x + y * y + z * z).into();
        }
        self.x.square() + self.y.square() + self.z.square()
    }
}

impl From<RPoint3> for Point3 {
    fn from(p: RPoint3) -> Self {
        Self::new(p.x.into(), p.y.into(), p.z.into())
    }
}

impl From<RVector3> for Vector3 {
    fn from(v: RVector3) -> Self {
        Self::new(v.x.into(), v.y.into(), v.z.into())
    }
}

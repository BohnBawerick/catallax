//! Exact numbers and the geometric values built from them.

use std::fmt;

use catallax_numbers::{ParseNumberError, Rational};
use serde::de::{self, Visitor};
use serde::{Deserialize, Deserializer, Serialize, Serializer};

/// An exact number.
///
/// It is read in the rational grammar of piece 1, which the exchange contract also uses
/// (`docs/JUDGE-INTERFACE.md`, ND2): a decimal with an optional exponent, such as `"1.1"`
/// or `"-2.5e-3"`, or a fraction such as `"1/3"`. It is always a JSON string, never a
/// JSON number, so no value passes through a double. It is written in canonical form:
/// a plain decimal when the value has a terminating one, such as `"0.0025"`, and
/// otherwise a fraction in lowest terms, such as `"1/3"`.
#[derive(Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Debug)]
pub struct Number(Box<Rational>);

impl Number {
    /// Wraps an exact value.
    pub fn new(value: Rational) -> Self {
        Number(Box::new(value))
    }

    /// Reads a number in the grammar of piece 1.
    pub fn parse(text: &str) -> Result<Self, ParseNumberError> {
        Rational::parse(text).map(Number::new)
    }

    /// The exact value.
    pub fn value(&self) -> &Rational {
        &self.0
    }

    /// The canonical text: a plain decimal when one is exact, else `p/q` in lowest terms.
    pub fn canonical(&self) -> String {
        self.0
            .to_decimal_exact()
            .unwrap_or_else(|| self.0.to_string())
    }
}

impl From<i64> for Number {
    fn from(value: i64) -> Self {
        Number::new(Rational::from_i64(value))
    }
}

impl From<Rational> for Number {
    fn from(value: Rational) -> Self {
        Number::new(value)
    }
}

impl fmt::Display for Number {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.canonical())
    }
}

impl Serialize for Number {
    fn serialize<S: Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        serializer.serialize_str(&self.canonical())
    }
}

impl<'de> Deserialize<'de> for Number {
    fn deserialize<D: Deserializer<'de>>(deserializer: D) -> Result<Self, D::Error> {
        deserializer.deserialize_str(NumberVisitor)
    }
}

struct NumberVisitor;

impl Visitor<'_> for NumberVisitor {
    type Value = Number;

    fn expecting(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str("an exact number written as a string, such as \"1.5\" or \"1/3\"")
    }

    fn visit_str<E: de::Error>(self, text: &str) -> Result<Number, E> {
        Number::parse(text)
            .map_err(|e| E::custom(format_args!("{text:?} is not an exact number: {e}")))
    }
}

/// A point in the coordinates of a sketch plane, in millimetres.
#[derive(Clone, PartialEq, Eq, Hash, Debug, Serialize, Deserialize)]
#[serde(transparent)]
pub struct Point2(pub [Number; 2]);

/// A point in model space, in millimetres.
#[derive(Clone, PartialEq, Eq, Hash, Debug, Serialize, Deserialize)]
#[serde(transparent)]
pub struct Point3(pub [Number; 3]);

/// A vector in model space. Where only its direction matters, its length is free.
#[derive(Clone, PartialEq, Eq, Hash, Debug, Serialize, Deserialize)]
#[serde(transparent)]
pub struct Vector3(pub [Number; 3]);

/// A rigid placement, as in the exchange contract (ND1): the quaternion `[w, x, y, z]`,
/// not all zero and of any length, rotates a local point, then `translation` is added.
#[derive(Clone, PartialEq, Eq, Hash, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Pose {
    /// Rotation quaternion `[w, x, y, z]`.
    pub rotation: [Number; 4],
    /// Translation added after the rotation.
    pub translation: Vector3,
}

impl Pose {
    /// No rotation and no translation.
    pub fn identity() -> Self {
        let n = Number::from;
        Pose {
            rotation: [n(1), n(0), n(0), n(0)],
            translation: Vector3([n(0), n(0), n(0)]),
        }
    }
}

/// A line in model space: a point on it and a direction.
#[derive(Clone, PartialEq, Eq, Hash, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Axis {
    /// A point on the line.
    pub point: Point3,
    /// Its direction, not zero.
    pub direction: Vector3,
}

/// A plane in model space: a point on it and a normal.
#[derive(Clone, PartialEq, Eq, Hash, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Plane {
    /// A point on the plane.
    pub point: Point3,
    /// Its normal, not zero.
    pub normal: Vector3,
}

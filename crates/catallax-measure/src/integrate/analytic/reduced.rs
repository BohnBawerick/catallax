//! Preserve exact rational antiderivatives behind algebraic construction DAGs.

use catallax_numbers::{Decider, Rational, Real, Refusal};
use catallax_shapes::{Circle, Curve, Ellipse, Line, Point3, RPoint3, RVector3, Vector3};

use super::{Analytic, Node, Program};

fn value(value: &Real, dec: &mut Decider) -> Result<Option<Rational>, Refusal> {
    match value.as_rational() {
        Some(value) => Ok(Some(value)),
        None => dec.rational_value(value),
    }
}

fn coordinates(values: [&Real; 3], dec: &mut Decider) -> Result<Option<[Rational; 3]>, Refusal> {
    let Some(x) = value(values[0], dec)? else {
        return Ok(None);
    };
    let Some(y) = value(values[1], dec)? else {
        return Ok(None);
    };
    let Some(z) = value(values[2], dec)? else {
        return Ok(None);
    };
    Ok(Some([x, y, z]))
}

fn point(point: &Point3, dec: &mut Decider) -> Result<Option<RPoint3>, Refusal> {
    Ok(coordinates([&point.x, &point.y, &point.z], dec)?.map(|[x, y, z]| RPoint3::new(x, y, z)))
}

fn vector(vector: &Vector3, dec: &mut Decider) -> Result<Option<RVector3>, Refusal> {
    Ok(
        coordinates([&vector.x, &vector.y, &vector.z], dec)?
            .map(|[x, y, z]| RVector3::new(x, y, z)),
    )
}

fn curve(curve: &Curve, dec: &mut Decider) -> Result<Option<Curve>, Refusal> {
    Ok(Some(match curve {
        Curve::Line(line) => {
            let Some(origin) = point(&line.origin, dec)? else {
                return Ok(None);
            };
            let Some(dir) = vector(&line.dir, dec)? else {
                return Ok(None);
            };
            Curve::Line(Line {
                origin: origin.into(),
                dir: dir.into(),
            })
        }
        Curve::Circle(circle) => {
            let Some(center) = point(&circle.center, dec)? else {
                return Ok(None);
            };
            let Some(normal) = vector(&circle.normal, dec)? else {
                return Ok(None);
            };
            let Some(ref_dir) = vector(&circle.ref_dir, dec)? else {
                return Ok(None);
            };
            let Some(radius_sq) = value(&circle.radius_sq, dec)? else {
                return Ok(None);
            };
            Curve::Circle(Circle {
                center: center.into(),
                normal: normal.into(),
                ref_dir: ref_dir.into(),
                radius_sq: radius_sq.into(),
            })
        }
        Curve::Ellipse(ellipse) => {
            let Some(center) = point(&ellipse.center, dec)? else {
                return Ok(None);
            };
            let Some(normal) = vector(&ellipse.normal, dec)? else {
                return Ok(None);
            };
            let Some(major_dir) = vector(&ellipse.major_dir, dec)? else {
                return Ok(None);
            };
            let Some(major_sq) = value(&ellipse.major_sq, dec)? else {
                return Ok(None);
            };
            let Some(minor_sq) = value(&ellipse.minor_sq, dec)? else {
                return Ok(None);
            };
            Curve::Ellipse(Ellipse {
                center: center.into(),
                normal: normal.into(),
                major_dir: major_dir.into(),
                major_sq: major_sq.into(),
                minor_sq: minor_sq.into(),
            })
        }
        _ => return Ok(None),
    }))
}

impl Program {
    pub(crate) fn analytic_reduced(
        &self,
        source: &Curve,
        from: &Point3,
        to: &Point3,
        dec: &mut Decider,
    ) -> Result<Option<Analytic>, Refusal> {
        let rational_from = point(from, dec)?;
        let rational_to = point(to, dec)?;
        let Some(curve) = curve(source, dec)? else {
            return Ok(None);
        };
        let mut constants = vec![None; self.exact.len()];
        for node in &self.nodes {
            if let Node::Exact(i) = node
                && constants[*i].is_none()
            {
                let Some(reduced) = value(&self.exact[*i], dec)? else {
                    return Ok(None);
                };
                constants[*i] = Some(reduced);
            }
        }
        if let (Some(from), Some(to)) = (rational_from, rational_to) {
            Ok(self.analytic_rational(&curve, &from, &to, &constants))
        } else {
            self.algebraic_endpoints(&curve, from, to, &constants, dec)
        }
    }
}

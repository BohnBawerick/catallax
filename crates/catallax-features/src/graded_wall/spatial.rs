//! Restrict a grade to each source face before displacing its boundary.

use super::{axis_coordinate, clearance, displaced, physical_faces, unit};
use crate::{
    grade::{Linear, RadiusRamp, distance_projection},
    polyhedron::{ConvexCell, PlaneFace},
    profile::sign,
};
use catallax_numbers::{Decider, Rational, Real, Refusal, Sign};
use catallax_recipe::{Origin, Role, StepId, ThickenSide};
use catallax_shapes::{Plane, Point3, Vector3};

pub(super) enum Spatial {
    Linear(Linear),
    Radius {
        origin: [Rational; 3],
        direction: Option<[Rational; 3]>,
        ramp: Box<RadiusRamp>,
    },
}

enum Piece {
    Constant(Rational),
    Affine {
        intercept: Rational,
        slope: Rational,
    },
}

impl Piece {
    fn at(&self, coordinate: &Real) -> (Real, Real) {
        match self {
            Self::Constant(value) => (value.clone().into(), Real::zero()),
            Self::Affine { intercept, slope } => {
                let slope = Real::from(slope.clone());
                (
                    &Real::from(intercept.clone()) + &(&slope * coordinate),
                    slope,
                )
            }
        }
    }
}

struct Band {
    value_at_origin: Real,
    gradient: Vector3,
    maximum: Real,
    cuts: Vec<PlaneFace>,
}

impl Band {
    fn constant(value: Real) -> Self {
        Self {
            maximum: value.clone(),
            value_at_origin: value,
            gradient: Vector3::zero(),
            cuts: Vec::new(),
        }
    }
}

fn point(coords: &[Rational; 3]) -> Point3 {
    let [x, y, z] = coords.each_ref().map(|q| Real::from(q.clone()));
    Point3::new(x, y, z)
}

fn vector(coords: &[Rational; 3]) -> Vector3 {
    let [x, y, z] = coords.each_ref().map(|q| Real::from(q.clone()));
    Vector3::new(x, y, z)
}

struct Face<'a> {
    plane: &'a Plane,
    normal: Vector3,
    vertices: Vec<&'a Point3>,
}

impl<'a> Face<'a> {
    fn new(source: &'a ConvexCell, index: usize, dec: &mut Decider) -> Result<Self, Refusal> {
        let plane = &source.faces[index].plane;
        Ok(Self {
            plane,
            normal: unit(&plane.normal, dec)?,
            vertices: source.polygons[index]
                .iter()
                .map(|&vertex| &source.vertices()[vertex])
                .collect(),
        })
    }

    fn bands(
        &self,
        origin_value: &Real,
        tangent: &Vector3,
        cuts: &[Rational],
        pieces: &[Piece],
        dec: &mut Decider,
    ) -> Result<Vec<Band>, Refusal> {
        let span = tangent.norm_squared();
        if sign(&span, dec)? == Sign::Zero {
            let mut piece = &pieces[0];
            for (cut, next) in cuts.iter().zip(&pieces[1..]) {
                if sign(&(origin_value - &Real::from(cut.clone())), dec)? != Sign::Positive {
                    break;
                }
                piece = next;
            }
            return Ok(vec![Band::constant(piece.at(origin_value).0)]);
        }
        let coordinate = |p: &Point3| origin_value + &p.minus(&self.plane.origin).dot(tangent);
        let first = coordinate(self.vertices[0]);
        let (mut lo, mut hi) = (first.clone(), first);
        for p in &self.vertices[1..] {
            let t = coordinate(p);
            if sign(&(&t - &lo), dec)? == Sign::Negative {
                lo = t.clone();
            }
            if sign(&(&t - &hi), dec)? == Sign::Positive {
                hi = t;
            }
        }
        let cut = |at: &Real, direction: i64| PlaneFace {
            plane: Plane {
                origin: self
                    .plane
                    .origin
                    .plus(&tangent.scaled(&(&(at - origin_value) / &span))),
                normal: tangent.scaled(&Real::from_i64(direction)),
                ref_dir: self.normal.clone(),
            },
            origins: Vec::new(),
        };
        let mut bands = Vec::new();
        for (i, piece) in pieces.iter().enumerate() {
            if i > 0 && cuts.get(i) == Some(&cuts[i - 1]) {
                continue;
            }
            let mut limits = [lo.clone(), hi.clone()];
            let mut planes = Vec::new();
            if i > 0 {
                let lower = Real::from(cuts[i - 1].clone());
                if sign(&(&hi - &lower), dec)? != Sign::Positive {
                    continue;
                }
                if sign(&(&lo - &lower), dec)? == Sign::Negative {
                    planes.push(cut(&lower, -1));
                    limits[0] = lower;
                }
            }
            if let Some(upper) = cuts.get(i) {
                let upper = Real::from(upper.clone());
                if sign(&(&lo - &upper), dec)? != Sign::Negative {
                    continue;
                }
                if sign(&(&hi - &upper), dec)? == Sign::Positive {
                    planes.push(cut(&upper, 1));
                    limits[1] = upper;
                }
            }
            let (value, slope) = piece.at(origin_value);
            let extreme = usize::from(sign(&slope, dec)? == Sign::Positive);
            bands.push(Band {
                value_at_origin: value,
                gradient: tangent.scaled(&slope),
                maximum: piece.at(&limits[extreme]).0,
                cuts: planes,
            });
        }
        Ok(bands)
    }
}

// A supporting plane through c proves that every projected point has norm at
// least |c|. A polygon's nearest point is in its plane, on an edge, or at a vertex.
fn outside_radius(
    points: &[Vector3],
    radius_sq: &Real,
    dec: &mut Decider,
) -> Result<bool, Refusal> {
    let supports = |c: &Vector3, dec: &mut Decider| -> Result<bool, Refusal> {
        let length_sq = c.norm_squared();
        if sign(&(&length_sq - radius_sq), dec)? == Sign::Negative {
            return Ok(false);
        }
        for p in points {
            if sign(&(&c.dot(p) - &length_sq), dec)? == Sign::Negative {
                return Ok(false);
            }
        }
        Ok(true)
    };

    let a = &points[0];
    for i in 1..points.len() - 1 {
        let normal = points[i].minus(a).cross(&points[i + 1].minus(a));
        let length_sq = normal.norm_squared();
        if sign(&length_sq, dec)? != Sign::Zero {
            let foot = normal.scaled(&(a.dot(&normal) / length_sq));
            if supports(&foot, dec)? {
                return Ok(true);
            }
            break;
        }
    }
    for (a, b) in points.iter().zip(points.iter().cycle().skip(1)) {
        if supports(a, dec)? {
            return Ok(true);
        }
        let edge = b.minus(a);
        let length_sq = edge.norm_squared();
        if sign(&length_sq, dec)? != Sign::Zero {
            let foot = a.minus(&edge.scaled(&(a.dot(&edge) / length_sq)));
            if supports(&foot, dec)? {
                return Ok(true);
            }
        }
    }
    Ok(false)
}

impl Spatial {
    fn bands(&self, face: &Face<'_>, dec: &mut Decider) -> Result<Vec<Band>, Refusal> {
        match self {
            Self::Linear(law) => {
                let direction = vector(&law.direction);
                let gradient =
                    direction.scaled(&(&Real::one() / &Real::from(law.length_sq.clone())));
                let tangent = gradient.minus(&face.normal.scaled(&gradient.dot(&face.normal)));
                let value = face.plane.origin.minus(&point(&law.from)).dot(&gradient);
                face.bands(
                    &value,
                    &tangent,
                    &[Rational::zero(), Rational::one()],
                    &[
                        Piece::Constant(law.start.clone()),
                        Piece::Affine {
                            intercept: law.start.clone(),
                            slope: &law.end - &law.start,
                        },
                        Piece::Constant(law.end.clone()),
                    ],
                    dec,
                )
            }
            Self::Radius {
                origin,
                direction,
                ramp,
            } => {
                let p = point(origin);
                if let Some(axis) = direction {
                    let axis = vector(axis);
                    let offset = face.plane.origin.minus(&p);
                    if sign(&face.normal.dot(&axis), dec)? == Sign::Zero
                        && sign(&face.normal.dot(&offset), dec)? == Sign::Zero
                    {
                        let tangent = unit(&face.normal.cross(&axis), dec)?;
                        let [r0, r1] = &ramp.radii;
                        let [start, end] = &ramp.values;
                        return face.bands(
                            &offset.dot(&tangent),
                            &tangent,
                            &[-r1, -r0, r0.clone(), r1.clone()],
                            &[
                                Piece::Constant(end.clone()),
                                Piece::Affine {
                                    intercept: ramp.intercept.clone(),
                                    slope: -&ramp.slope,
                                },
                                Piece::Constant(start.clone()),
                                Piece::Affine {
                                    intercept: ramp.intercept.clone(),
                                    slope: ramp.slope.clone(),
                                },
                                Piece::Constant(end.clone()),
                            ],
                            dec,
                        );
                    }
                }
                let projection = distance_projection(direction.as_ref())
                    .each_ref()
                    .map(vector);
                let points: Vec<_> = face
                    .vertices
                    .iter()
                    .map(|vertex| {
                        let delta = vertex.minus(&p);
                        let [x, y, z] = projection.each_ref().map(|row| row.dot(&delta));
                        Vector3::new(x, y, z)
                    })
                    .collect();
                let r0_sq = Real::from(ramp.squared[0].clone());
                let mut inside = true;
                for p in &points {
                    if sign(&(&p.norm_squared() - &r0_sq), dec)? == Sign::Positive {
                        inside = false;
                        break;
                    }
                }
                let value = if inside {
                    &ramp.values[0]
                } else if outside_radius(&points, &Real::from(ramp.squared[1].clone()), dec)? {
                    &ramp.values[1]
                } else {
                    return Err(Refusal::new(
                        "distance.grade-curved-wall-unsupported",
                        "the distance ramp on a kept or selected face needs a curved offset surface",
                    ));
                };
                Ok(vec![Band::constant(value.clone().into())])
            }
        }
    }

    pub(super) fn thicken(
        &self,
        source: &ConvexCell,
        index: usize,
        side: ThickenSide,
        step: &StepId,
        dec: &mut Decider,
    ) -> Result<Vec<ConvexCell>, Refusal> {
        let face = Face::new(source, index, dec)?;
        let bands = self.bands(&face, dec)?;
        let faces = source.sweep_faces(index, &bands[0].maximum, side, step, dec)?;
        let rise = if side == ThickenSide::Front {
            face.normal
        } else {
            face.normal.scaled(&Real::from_i64(-1))
        };
        let mut cells = Vec::new();
        for band in bands {
            let mut faces = faces.clone();
            faces[1].plane = displaced(face.plane, &rise, &band.value_at_origin, &band.gradient);
            faces.extend(band.cuts);
            let cell = ConvexCell::new(faces, dec)?;
            if !cell.is_empty() {
                cells.push(physical_faces(cell));
            }
        }
        Ok(cells)
    }

    pub(super) fn cavity(
        &self,
        source: &ConvexCell,
        open: &[bool],
        walls: &[[usize; 2]; 3],
        step: &StepId,
        dec: &mut Decider,
    ) -> Result<Vec<ConvexCell>, Refusal> {
        let mut fixed = source.faces.clone();
        let mut split = Vec::new();
        let mut maxima = vec![Real::zero(); fixed.len()];
        for (index, original) in fixed.iter_mut().enumerate() {
            let origins = std::mem::take(&mut original.origins);
            if open[index] {
                continue;
            }
            let face = Face::new(source, index, dec)?;
            let mut choices = Vec::new();
            for band in self.bands(&face, dec)? {
                if sign(&(&band.maximum - &maxima[index]), dec)? == Sign::Positive {
                    maxima[index] = band.maximum;
                }
                let plane = displaced(
                    face.plane,
                    &face.normal,
                    &(-band.value_at_origin),
                    &band.gradient.scaled(&Real::from_i64(-1)),
                );
                let origins = origins
                    .iter()
                    .map(|of| Origin {
                        step: step.clone(),
                        role: Role::Inner {
                            of: Box::new(of.clone()),
                        },
                    })
                    .collect();
                let mut faces = band.cuts;
                faces.push(PlaneFace { plane, origins });
                choices.push(faces);
            }
            split.push(choices);
        }
        // Opposing maxima leave a rectangular region strictly inside every wall.
        for (axis, &[lo, hi]) in walls.iter().enumerate() {
            let gap = axis_coordinate(&source.faces[hi].plane.origin, axis)
                - axis_coordinate(&source.faces[lo].plane.origin, axis);
            if sign(&(&gap - &(&maxima[lo] + &maxima[hi])), dec)? != Sign::Positive {
                return Err(clearance());
            }
        }
        let mut variable = Vec::new();
        for mut choices in split {
            if choices.len() == 1 {
                fixed.extend(choices.pop().expect("one unsplit wall"));
            } else {
                variable.push(choices);
            }
        }
        let mut cells = vec![ConvexCell::new(fixed, dec)?];
        for choices in variable {
            let mut next = Vec::new();
            for cell in cells {
                let faces: Vec<_> = cell
                    .faces
                    .into_iter()
                    .zip(cell.polygons)
                    .filter_map(|(face, polygon)| (!polygon.is_empty()).then_some(face))
                    .collect();
                for choice in &choices {
                    let cell = ConvexCell::new([faces.clone(), choice.clone()].concat(), dec)?;
                    if !cell.is_empty() {
                        next.push(cell);
                    }
                }
            }
            cells = next;
        }
        Ok(cells.into_iter().map(physical_faces).collect())
    }
}

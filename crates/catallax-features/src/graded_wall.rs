//! Exact normal displacements for a clamped thickness in part coordinates.

use catallax_numbers::{Decider, Rational, Real, Refusal, Sign};
use catallax_recipe::{FaceRef, Grade, Origin, Role, StepId, ThickenSide};
use catallax_shapes::{Plane, Point3, Vector3};

use crate::{
    grade::{self, RadiusRamp},
    polyhedron::{ConvexCell, PlaneFace},
    profile::sign,
};

mod spatial;
use spatial::Spatial;

/// A positive thickness law admitted by the graded distance operations.
pub enum Thickness {
    /// The scalar distance operation covers a constant law.
    Constant(Rational),
    /// A nonconstant law evaluated at the perpendicular foot on each source face.
    Graded(WallGrade),
}

/// Nonconstant wall thickness in part coordinates.
pub struct WallGrade {
    law: Law,
}

enum Law {
    Axis(Box<Linear>),
    Spatial(Box<Spatial>),
}

/// Positive, clamped linear thickness along one part axis.
struct Linear {
    axis: usize,
    anchors: [Rational; 2],
    values: [Rational; 2],
}

/// Refusal for a source outside the exact planar coverage.
pub fn unsupported_source(detail: impl Into<String>) -> Refusal {
    Refusal::new("distance.grade-source-unsupported", detail)
}

/// True for the explicit, unfiltered all-faces opening.
pub fn opens_all(open: &[FaceRef]) -> bool {
    !open.is_empty()
        && open
            .iter()
            .all(|r| r.face.is_none() && r.touching.is_empty() && r.near.is_none())
}

impl Thickness {
    /// Validates the grade and its positive endpoint thicknesses.
    pub fn new(grade: &Grade) -> Result<Self, Refusal> {
        let positive = |value: &Rational| {
            if value.sign() != Sign::Positive {
                Err(Refusal::new(
                    "distance.non-positive-thickness",
                    "every graded wall thickness must be positive",
                ))
            } else {
                Ok(value.clone())
            }
        };
        let law = match grade {
            Grade::Constant(value) => return Ok(Self::Constant(positive(value.value())?)),
            Grade::Linear {
                from,
                to,
                start,
                end,
            } => {
                let mut values = [positive(start.value())?, positive(end.value())?];
                if from == to {
                    return Err(Refusal::new(
                        "distance.grade-coincident-anchors",
                        "a linear thickness needs distinct anchors",
                    ));
                }
                if values[0] == values[1] {
                    return Ok(Self::Constant(values[0].clone()));
                }
                let axes: Vec<_> = (0..3).filter(|&i| from.0[i] != to.0[i]).collect();
                let [axis] = *axes.as_slice() else {
                    let coordinates =
                        |p: &catallax_recipe::Point3| p.0.each_ref().map(|n| n.value().clone());
                    let ramp = grade::Linear::new(
                        coordinates(from),
                        coordinates(to),
                        values[0].clone(),
                        values[1].clone(),
                    );
                    return Ok(Self::graded(Law::Spatial(Box::new(Spatial::Linear(ramp)))));
                };
                let mut anchors = [from.0[axis].value().clone(), to.0[axis].value().clone()];
                if anchors[0] > anchors[1] {
                    anchors.swap(0, 1);
                    values.swap(0, 1);
                }
                Law::Axis(Box::new(Linear {
                    axis,
                    anchors,
                    values,
                }))
            }
            Grade::Point {
                r0, r1, start, end, ..
            }
            | Grade::Radial {
                r0, r1, start, end, ..
            } => {
                let values = [positive(start.value())?, positive(end.value())?];
                crate::lattice::validate_distance_radii(grade, "distance.grade-invalid-radius")?;
                if values[0] == values[1] {
                    return Ok(Self::Constant(values[0].clone()));
                }
                let ramp = RadiusRamp::new([r0.value().clone(), r1.value().clone()], values);
                let (origin, direction) = match grade {
                    Grade::Point { centre, .. } => (&centre.0, None),
                    Grade::Radial {
                        origin, direction, ..
                    } => (&origin.0, Some(&direction.0)),
                    _ => unreachable!("point or radial grade"),
                };
                let coordinates =
                    |p: &[catallax_recipe::Number; 3]| p.each_ref().map(|n| n.value().clone());
                Law::Spatial(Box::new(Spatial::Radius {
                    origin: coordinates(origin),
                    direction: direction.map(coordinates),
                    ramp: Box::new(ramp),
                }))
            }
            Grade::Distance { .. } => {
                return Err(Refusal::new(
                    "distance.grade-function-unsupported",
                    "wall thickness needs a grade defined in part coordinates",
                ));
            }
        };
        Ok(Self::graded(law))
    }

    fn graded(law: Law) -> Self {
        Self::Graded(WallGrade { law })
    }
}

impl WallGrade {
    /// Upper bound on displacement, including the clamps.
    pub fn maximum(&self) -> &Rational {
        let (start, end) = match &self.law {
            Law::Axis(grade) => (&grade.values[0], &grade.values[1]),
            Law::Spatial(grade) => match grade.as_ref() {
                Spatial::Linear(grade) => (&grade.start, &grade.end),
                Spatial::Radius { ramp, .. } => (&ramp.values[0], &ramp.values[1]),
            },
        };
        std::cmp::max(start, end)
    }

    /// Constructs the exact front or back sweep of a selected planar face.
    pub fn thicken(
        &self,
        source: &ConvexCell,
        refs: &[FaceRef],
        side: ThickenSide,
        step: &StepId,
        dec: &mut Decider,
    ) -> Result<Vec<ConvexCell>, Refusal> {
        let index = select_face(source, refs, dec)?;
        match &self.law {
            Law::Axis(grade) => grade.thicken(source, index, side, step, dec),
            Law::Spatial(grade) => grade.thicken(source, index, side, step, dec),
        }
    }

    /// Constructs the exact cavity in an axis-aligned rectangular source.
    pub fn cavity(
        &self,
        source: &ConvexCell,
        open: &[bool],
        step: &StepId,
        dec: &mut Decider,
    ) -> Result<Vec<ConvexCell>, Refusal> {
        let walls = box_walls(source, dec)?;
        match &self.law {
            Law::Axis(grade) => grade.cavity(source, open, &walls, step, dec),
            Law::Spatial(grade) => grade.cavity(source, open, &walls, step, dec),
        }
    }
}

fn basis(axis: usize) -> Vector3 {
    let [x, y, z] = std::array::from_fn(|i| Real::from_i64(i64::from(i == axis)));
    Vector3::new(x, y, z)
}

fn axis_coordinate(point: &Point3, axis: usize) -> &Real {
    [&point.x, &point.y, &point.z][axis]
}

fn unit(vector: &Vector3, dec: &mut Decider) -> Result<Vector3, Refusal> {
    let length = vector.norm_squared().sqrt();
    let length = dec.rational_value(&length)?.map_or(length, Real::from);
    Ok(vector.scaled(&(&Real::one() / &length)))
}

fn physical_faces(mut cell: ConvexCell) -> ConvexCell {
    for (face, polygon) in cell.faces.iter_mut().zip(&cell.polygons) {
        if polygon.is_empty() {
            face.origins.clear();
        }
    }
    cell
}

/// A halfspace with a new support and the source plane's in-plane direction.
fn displaced(source: &Plane, rise: &Vector3, value: &Real, gradient: &Vector3) -> Plane {
    let normal = rise.minus(gradient);
    let ref_dir = source
        .ref_dir
        .minus(&normal.scaled(&(source.ref_dir.dot(&normal) / normal.norm_squared())));
    Plane {
        origin: source.origin.plus(&rise.scaled(value)),
        normal,
        ref_dir,
    }
}

/// The three affine pieces in increasing order along the grade axis.
#[derive(Clone, Copy)]
enum Band {
    Low,
    Ramp,
    High,
}

impl Linear {
    fn law(&self, band: Band, coordinate: &Real) -> (Real, Real) {
        let slope = match band {
            Band::Low => return (self.values[0].clone().into(), Real::zero()),
            Band::High => return (self.values[1].clone().into(), Real::zero()),
            Band::Ramp => (&self.values[1] - &self.values[0])
                .checked_div(&(&self.anchors[1] - &self.anchors[0]))
                .expect("distinct ordered anchors"),
        };
        let slope = Real::from(slope);
        let value = &Real::from(self.values[0].clone())
            + &(&slope * &(coordinate - &Real::from(self.anchors[0].clone())));
        (value, slope)
    }

    fn sample(&self, coordinate: &Real, dec: &mut Decider) -> Result<Real, Refusal> {
        let band = if sign(&(coordinate - &Real::from(self.anchors[0].clone())), dec)?
            != Sign::Positive
        {
            Band::Low
        } else if sign(&(coordinate - &Real::from(self.anchors[1].clone())), dec)? != Sign::Negative
        {
            Band::High
        } else {
            Band::Ramp
        };
        Ok(self.law(band, coordinate).0)
    }

    fn cells(
        &self,
        origin: &Point3,
        tangent: &Vector3,
        ref_dir: &Vector3,
        planes: impl Fn(Band, &mut Decider) -> Result<Vec<PlaneFace>, Refusal>,
        dec: &mut Decider,
    ) -> Result<Vec<ConvexCell>, Refusal> {
        let length_sq = tangent.norm_squared();
        let coordinate = axis_coordinate(origin, self.axis);
        let cut =
            |i: usize, direction: i64| PlaneFace {
                plane: Plane {
                    origin: origin.plus(&tangent.scaled(
                        &(&(&Real::from(self.anchors[i].clone()) - coordinate) / &length_sq),
                    )),
                    normal: tangent.scaled(&Real::from_i64(direction)),
                    ref_dir: ref_dir.clone(),
                },
                origins: Vec::new(),
            };
        let mut result = Vec::new();
        for band in [Band::Low, Band::Ramp, Band::High] {
            let mut planes = planes(band, dec)?;
            match band {
                Band::Low => planes.push(cut(0, 1)),
                Band::Ramp => planes.extend([cut(0, -1), cut(1, 1)]),
                Band::High => planes.push(cut(1, -1)),
            }
            let cell = ConvexCell::new(planes, dec)?;
            if !cell.is_empty() {
                result.push(physical_faces(cell));
            }
        }
        Ok(result)
    }

    /// A single convex planar face's swept solid, with the grade evaluated at
    /// the perpendicular projection onto the source face.
    fn thicken(
        &self,
        source: &ConvexCell,
        index: usize,
        side: ThickenSide,
        step: &StepId,
        dec: &mut Decider,
    ) -> Result<Vec<ConvexCell>, Refusal> {
        let base = &source.faces[index].plane;
        let axis = basis(self.axis);
        let outward = unit(&base.normal, dec)?;
        let tangent = axis.minus(&outward.scaled(&axis.dot(&outward)));
        let thickness = self.sample(axis_coordinate(&base.origin, self.axis), dec)?;
        if sign(&tangent.norm_squared(), dec)? == Sign::Zero {
            return Ok(vec![physical_faces(
                source.face_sweep(index, &thickness, side, step, dec)?,
            )]);
        }
        let rise = outward.scaled(&Real::from_i64(if side == ThickenSide::Front {
            1
        } else {
            -1
        }));
        let faces = source.sweep_faces(index, &thickness, side, step, dec)?;
        self.cells(
            &base.origin,
            &tangent,
            &rise,
            |band, _| {
                let (value, slope) = self.law(band, axis_coordinate(&base.origin, self.axis));
                let mut faces = faces.clone();
                faces[1].plane = displaced(base, &rise, &value, &tangent.scaled(&slope));
                Ok(faces)
            },
            dec,
        )
    }

    /// The box cavity. A kept wall tests the clearance along its inward normal
    /// at the perpendicular foot on that source face. Open walls bound the
    /// cavity at the original box planes.
    fn cavity(
        &self,
        source: &ConvexCell,
        open: &[bool],
        walls: &[[usize; 2]; 3],
        step: &StepId,
        dec: &mut Decider,
    ) -> Result<Vec<ConvexCell>, Refusal> {
        let coordinate = |axis: usize, side: usize| {
            axis_coordinate(&source.faces[walls[axis][side]].plane.origin, axis)
        };
        let mut limits = [
            coordinate(self.axis, 0).clone(),
            coordinate(self.axis, 1).clone(),
        ];
        for (side, limit) in limits.iter_mut().enumerate() {
            if !open[walls[self.axis][side]] {
                let value = self.sample(limit, dec)?;
                *limit = if side == 0 {
                    &*limit + &value
                } else {
                    &*limit - &value
                };
            }
        }
        if sign(&(&limits[1] - &limits[0]), dec)? != Sign::Positive {
            return Err(clearance());
        }
        let end_values = [self.sample(&limits[0], dec)?, self.sample(&limits[1], dec)?];
        let max = if sign(&(&end_values[0] - &end_values[1]), dec)? == Sign::Positive {
            &end_values[0]
        } else {
            &end_values[1]
        };
        for (axis, pair) in walls
            .iter()
            .enumerate()
            .filter(|(axis, _)| *axis != self.axis)
        {
            let kept = pair.iter().filter(|&&wall| !open[wall]).count();
            let gap = coordinate(axis, 1) - coordinate(axis, 0);
            if sign(&(&gap - &(max * &Real::from_i64(kept as i64))), dec)? != Sign::Positive {
                return Err(clearance());
            }
        }
        let axis = basis(self.axis);
        self.cells(
            &Point3::new(Real::zero(), Real::zero(), Real::zero()),
            &axis,
            &basis((self.axis + 1) % 3),
            |band, dec| {
                let mut faces = source.faces.clone();
                for (i, face) in faces.iter_mut().enumerate() {
                    if open[i] {
                        face.origins.clear();
                        continue;
                    }
                    let coord = axis_coordinate(&face.plane.origin, self.axis);
                    let normal = unit(&face.plane.normal, dec)?;
                    let (value, slope) = if sign(&normal.dot(&axis), dec)? == Sign::Zero {
                        self.law(band, coord)
                    } else {
                        (self.sample(coord, dec)?, Real::zero())
                    };
                    // A cavity retains the source plane's outward sense. Distance
                    // measured inward adds the grade's gradient to this normal.
                    face.plane =
                        displaced(&face.plane, &normal, &(-value), &axis.scaled(&(-slope)));
                    for of in &mut face.origins {
                        *of = Origin {
                            step: step.clone(),
                            role: Role::Inner {
                                of: Box::new(of.clone()),
                            },
                        };
                    }
                }
                Ok(faces)
            },
            dec,
        )
    }
}

fn select_face(source: &ConvexCell, refs: &[FaceRef], dec: &mut Decider) -> Result<usize, Refusal> {
    let selected: Vec<_> = source
        .select(refs, dec)?
        .into_iter()
        .enumerate()
        .filter_map(|(i, selected)| selected.then_some(i))
        .collect();
    let [index] = *selected.as_slice() else {
        return Err(if selected.is_empty() {
            Refusal::new(
                "reference.lost-face",
                "thicken needs a surviving selected face",
            )
        } else {
            Refusal::new(
                "distance.grade-selection-unsupported",
                "graded thicken needs exactly one selected face",
            )
        });
    };
    Ok(index)
}

fn box_walls(source: &ConvexCell, dec: &mut Decider) -> Result<[[usize; 2]; 3], Refusal> {
    let mut walls = [[None; 2]; 3];
    if source.faces.len() != 6 || source.vertices().len() != 8 {
        return Err(unsupported_source(
            "graded shell needs an axis-aligned rectangular solid",
        ));
    }
    for (i, face) in source.faces.iter().enumerate() {
        let n = &face.plane.normal;
        let mut components = Vec::new();
        for (axis, component) in [&n.x, &n.y, &n.z].into_iter().enumerate() {
            let direction = sign(component, dec)?;
            if direction != Sign::Zero {
                components.push((axis, usize::from(direction == Sign::Positive)));
            }
        }
        let [(axis, side)] = *components.as_slice() else {
            return Err(unsupported_source(
                "graded shell needs an axis-aligned rectangular solid",
            ));
        };
        if walls[axis][side].replace(i).is_some() {
            return Err(unsupported_source(
                "graded shell needs one source face on each side of the box",
            ));
        }
    }
    Ok(walls.map(|pair| pair.map(|wall| wall.expect("six distinct box walls"))))
}

fn clearance() -> Refusal {
    Refusal::new(
        "distance.grade-clearance",
        "graded shell needs positive cavity width between opposing walls throughout the cavity",
    )
}

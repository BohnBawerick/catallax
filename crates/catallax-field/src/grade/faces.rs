use super::GradeDef;
use crate::{Instruction, Plan, geometry::Leaf, mass};
use catallax_numbers::{Decider, Rational, Refusal};
use catallax_recipe::Grade;
use catallax_shapes::Surface;

fn unsupported() -> Refusal {
    Refusal::new(
        "field.grade-face-coverage",
        "distance grades currently require a block or a rectangular extrusion aligned with the part axes",
    )
}

impl Plan {
    pub(crate) fn compile_face_grade(
        &self,
        region: usize,
        region_faces: &[usize],
        grade: &Grade,
        dec: &mut Decider,
    ) -> Result<Option<GradeDef>, Refusal> {
        let Grade::Distance {
            faces,
            near,
            far,
            span,
        } = grade
        else {
            return Ok(None);
        };
        if faces.is_empty() {
            return Ok(None);
        }
        let Instruction::Leaf(leaf) = self.instructions[region] else {
            return Err(unsupported());
        };
        let bounds = match &self.leaves[leaf] {
            Leaf::Primitive(shape) => mass::axis_block(shape),
            Leaf::Feature(feature) => mass::prism(feature, dec)?,
            _ => None,
        }
        .ok_or_else(unsupported)?;
        if bounds.iter().any(|[lo, hi]| lo >= hi) {
            return Err(unsupported());
        }
        let mut origins: [_; 6] = std::array::from_fn(|_| Vec::new());
        for &id in region_faces {
            dec.check_deadline()?;
            let face = &self.faces[id];
            if face.leaf != leaf {
                return Err(unsupported());
            }
            let Surface::Plane(plane) = &face.surface else {
                return Err(unsupported());
            };
            let position = plane.origin.as_rational().ok_or_else(unsupported)?;
            let normal = plane.normal.as_rational().ok_or_else(unsupported)?;
            let position = [position.x, position.y, position.z];
            let normal = [normal.x, normal.y, normal.z];
            let mut axes = (0..3).filter(|&axis| !normal[axis].is_zero());
            let (Some(axis), None) = (axes.next(), axes.next()) else {
                return Err(unsupported());
            };
            let side = bounds[axis]
                .iter()
                .position(|bound| bound == &position[axis])
                .ok_or_else(unsupported)?;
            origins[2 * axis + side].push(&face.origin);
        }
        if origins.iter().any(Vec::is_empty) {
            return Err(unsupported());
        }
        let mut selected = [false; 6];
        for reference in faces {
            dec.check_deadline()?;
            let candidates = (0..6).filter(|&id| {
                origins[id].iter().any(|origin| reference.matches(origin))
                    && reference.touching.iter().all(|pattern| {
                        (0..6).filter(|&other| id / 2 != other / 2).any(|other| {
                            origins[other].iter().any(|origin| pattern.matches(origin))
                        })
                    })
            });
            let mut candidates = candidates.peekable();
            if candidates.peek().is_none() {
                return Err(Refusal::new(
                    "field.grade-face-missing",
                    "a distance grade reference selects no face of the region",
                ));
            }
            if let Some(near) = &reference.near {
                let point = near.0.each_ref().map(|coordinate| coordinate.value());
                let squared_distance = |id: usize| {
                    (0..3).fold(Rational::zero(), |sum, axis| {
                        let [lo, hi] = &bounds[axis];
                        let coordinate = if axis == id / 2 {
                            &bounds[axis][id % 2]
                        } else {
                            point[axis].clamp(lo, hi)
                        };
                        let difference = point[axis] - coordinate;
                        sum + &difference * &difference
                    })
                };
                let mut ranked: Vec<_> = candidates.map(|id| (squared_distance(id), id)).collect();
                ranked.sort_unstable();
                if ranked.get(1).is_some_and(|next| next.0 == ranked[0].0) {
                    return Err(Refusal::new(
                        "field.grade-face-ambiguous",
                        "a distance grade's nearest face must be strictly nearer than the other matches",
                    ));
                }
                selected[ranked[0].1] = true;
            } else {
                for id in candidates {
                    selected[id] = true;
                }
            }
        }
        Ok(Some(GradeDef::from_faces(
            (0..6)
                .filter(|&id| selected[id])
                .map(|id| (id / 2, id % 2 == 1, bounds[id / 2][id % 2].clone())),
            near.value(),
            far.value(),
            span.value(),
        )))
    }
}

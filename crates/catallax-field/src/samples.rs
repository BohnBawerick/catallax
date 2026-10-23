//! Query-only sampling support. These bounds and candidates never build geometry.
use crate::{
    FaceId, Instruction, Plan,
    draft::{Cell, Region},
    geometry::Leaf,
};
use catallax_numbers::{Decider, Rational, Real, Refusal, Sign};
use catallax_recipe::{BooleanKind, Origin};
use catallax_shapes::{Isometry, Point3, RPoint3, RVector3, Vector3};

mod cell;

pub(crate) type Bounds = [[Rational; 2]; 3];

pub(crate) fn transform(bounds: &Bounds, pose: &Isometry) -> Result<Bounds, Refusal> {
    let mut result: Option<Bounds> = None;
    for corner in 0..8 {
        let [x, y, z] = std::array::from_fn(|i| bounds[i][(corner >> i) & 1].clone());
        let p = pose
            .apply_rational_point(&RPoint3::new(x, y, z))
            .ok_or_else(placement_error)?;
        let coordinates = [p.x, p.y, p.z];
        if let Some(b) = &mut result {
            for i in 0..3 {
                b[i][0] = b[i][0].clone().min(coordinates[i].clone());
                b[i][1] = b[i][1].clone().max(coordinates[i].clone());
            }
        } else {
            result = Some(coordinates.map(|x| [x.clone(), x]));
        }
    }
    Ok(result.unwrap())
}
pub(crate) fn placement_error() -> Refusal {
    Refusal::new(
        "field.sample-algebraic-placement",
        "the rational candidate API needs a native rational placement",
    )
}
pub(crate) fn fraction(i: u32, n: u32) -> Rational {
    Rational::from_ratio(i64::from(2 * i + 1), i64::from(2 * n)).unwrap()
}
pub(crate) fn circle(n: u32) -> Vec<[Rational; 2]> {
    let mut points = Vec::with_capacity((4 * n) as usize);
    for quarter in 0..4 {
        for i in 0..n {
            let t = fraction(i, n);
            let squared = &t * &t;
            let denominator = &Rational::one() + &squared;
            let x = &(Rational::one() - squared) / &denominator;
            let y = Rational::from_i64(2) * t / denominator;
            points.push(match quarter {
                0 => [x, y],
                1 => [-y, x],
                2 => [-x, -y],
                _ => [y, -x],
            });
        }
    }
    points
}
impl Plan {
    /// Visits rational candidates from current primitive, sweep and draft definitions.
    /// Each point carries its exact registry handle: one origin can have pieces
    /// on different supports. No historical primitive is sampled in place of a
    /// replacement definition. Use `face_query` to establish actual membership.
    /// `divisions` is in 1..=64. Returning false from `visit` stops cooperatively
    /// and returns `Ok(false)`, never a claim that all candidates were visited.
    /// Through-all walls are sampled over a conservative bound on the field
    /// result, not an arbitrary finite replacement for the unbounded tool.
    /// An algebraic placement without a native rational frame returns
    /// `field.sample-algebraic-placement`; no candidate coordinates are rounded.
    pub fn visit_candidates(
        &self,
        divisions: u32,
        mut visit: impl FnMut(RPoint3, FaceId) -> bool,
    ) -> Result<bool, Refusal> {
        if !(1..=64).contains(&divisions) {
            return Err(Refusal::new(
                "field.sample-options",
                "sampling divisions must be in 1..=64",
            ));
        }
        let has_derived = self.faces.iter().any(|face| {
            face.in_result
                && self.leaf_registers[face.leaf].is_some()
                && matches!(self.leaves[face.leaf], Leaf::Feature(_) | Leaf::Draft(_) | Leaf::Cell(_))
        });
        if self.is_proven_empty() {
            return Ok(true);
        }
        let mut bounds: Vec<Option<Bounds>> = Vec::with_capacity(if has_derived {
            self.instructions.len()
        } else {
            0
        });
        if has_derived {
            for instruction in &self.instructions {
                let b = match *instruction {
                    Instruction::Empty => Some(std::array::from_fn(|_| {
                        [Rational::zero(), Rational::zero()]
                    })),
                    Instruction::Leaf(i) => match &self.leaves[i] {
                        Leaf::Primitive(p) => Some(p.sample_bounds()?),
                        Leaf::Feature(f) => f.sample_bounds()?,
                        Leaf::Cell(c) => c.sample_bounds()?,
                        Leaf::Draft(region) => Some(region.sample_bounds()?),
                    },
                    Instruction::Boolean { kind, a, b } => match kind {
                        BooleanKind::Subtract => bounds[a].clone(),
                        BooleanKind::Union => match (&bounds[a], &bounds[b]) {
                            (Some(a), Some(b)) => Some(std::array::from_fn(|i| {
                                [
                                    a[i][0].clone().min(b[i][0].clone()),
                                    a[i][1].clone().max(b[i][1].clone()),
                                ]
                            })),
                            _ => None,
                        },
                        BooleanKind::Intersect => match (&bounds[a], &bounds[b]) {
                            (Some(a), Some(b)) => Some(std::array::from_fn(|i| {
                                [
                                    a[i][0].clone().max(b[i][0].clone()),
                                    a[i][1].clone().min(b[i][1].clone()),
                                ]
                            })),
                            (Some(a), None) | (None, Some(a)) => Some(a.clone()),
                            (None, None) => None,
                        },
                    },
                };
                bounds.push(b);
            }
        }
        let result = if has_derived {
            Some(bounds[self.root].as_ref().ok_or_else(|| {
                Refusal::new(
                    "field.sample-unbounded",
                    "result has no finite feature sampling bound",
                )
            })?)
        } else {
            None
        };
        let angular = std::cell::OnceCell::new();
        let radial_grid = std::cell::OnceCell::new();
        for (index, face) in self.faces.iter().enumerate() {
            if !face.in_result || self.leaf_registers[face.leaf].is_none() {
                continue;
            }
            let complete = match &self.leaves[face.leaf] {
                Leaf::Feature(feature) => feature.sample_candidates(
                    self.definition_role(index),
                    divisions,
                    result.expect("derived candidates have a finite sampling bound"),
                    &mut |p| visit(p, FaceId(index)),
                )?,
                Leaf::Draft(region) => {
                    region.sample_candidates(&face.origin, divisions, &mut |p| {
                        visit(p, FaceId(index))
                    })?
                }
                Leaf::Primitive(primitive) => primitive.sample_candidates(
                    self.definition_role(index),
                    divisions,
                    &angular,
                    &mut |p| visit(p, FaceId(index)),
                )?,
                Leaf::Cell(cell) => {
                    let Some(capped) = cell.capped_chart() else {
                        continue;
                    };
                    let support = cell.support_index(&face.polynomial).ok_or_else(|| {
                        Refusal::new(
                            "field.sample-support",
                            "native cell face has no defining support",
                        )
                    })?;
                    capped_candidates(
                        capped,
                        support,
                        divisions,
                        &angular,
                        &radial_grid,
                        &mut |p| visit(p, FaceId(index)),
                    )?
                }
            };
            if !complete {
                return Ok(false);
            }
        }
        // Keep the existing feature traversal and cooperative cancellation order.
        // Cell candidates are a separate, definition-native pass. They are not
        // boundary certificates; callers must use the same face queries as for
        // all other candidate sources.
        cell::visit(self, divisions, result, &mut visit)
    }
}

fn native_scalar(value: &Real, dec: &mut Decider) -> Result<Rational, Refusal> {
    if let Some(value) = value.rational_ref() {
        return Ok(value.clone());
    }
    dec.rational_value(value)?.ok_or_else(placement_error)
}

fn native_point(point: &Point3, dec: &mut Decider) -> Result<RPoint3, Refusal> {
    Ok(RPoint3::new(
        native_scalar(&point.x, dec)?,
        native_scalar(&point.y, dec)?,
        native_scalar(&point.z, dec)?,
    ))
}

fn capped_candidates(
    cell: &catallax_features::draft::capped::CappedCylinder,
    support: usize,
    divisions: u32,
    angular: &std::cell::OnceCell<Vec<[Rational; 2]>>,
    radial_grid: &std::cell::OnceCell<Vec<Rational>>,
    visit: &mut impl FnMut(RPoint3) -> bool,
) -> Result<bool, Refusal> {
    let mut dec = Decider::default();
    if support == 0 {
        for direction in angular.get_or_init(|| circle(divisions)) {
            let q = direction
                .each_ref()
                .map(|x| &cell.radius * &Real::from(x.clone()));
            let heights = cell.heights_at(&q);
            let span = &heights[1] - &heights[0];
            if crate::geometry::sign(&span, &mut dec)? != Sign::Positive {
                continue;
            }
            for i in 0..divisions {
                let height = &heights[0] + &(&span * &Real::from(fraction(i, divisions)));
                if !visit(native_point(&cell.point(&q, &height), &mut dec)?) {
                    return Ok(false);
                }
            }
        }
    } else {
        let grid = radial_grid.get_or_init(|| {
            (0..divisions)
                .map(|i| Rational::from_i64(2) * fraction(i, divisions) - Rational::one())
                .collect()
        });
        for x in grid {
            for y in grid {
                if &(x * x) + &(y * y) >= Rational::one() {
                    continue;
                }
                let q = [x, y].map(|x| &cell.radius * &Real::from(x.clone()));
                let heights = cell.heights_at(&q);
                if crate::geometry::sign(&(&heights[1] - &heights[0]), &mut dec)? != Sign::Positive
                {
                    continue;
                }
                if !visit(native_point(
                    &cell.point(&q, &heights[support - 1]),
                    &mut dec,
                )?) {
                    return Ok(false);
                }
            }
        }
    }
    Ok(true)
}

fn native_vector(vector: &Vector3, dec: &mut Decider) -> Result<RVector3, Refusal> {
    Ok(RVector3::new(
        native_scalar(&vector.x, dec)?,
        native_scalar(&vector.y, dec)?,
        native_scalar(&vector.z, dec)?,
    ))
}

impl Region {
    fn sample_candidates(
        &self,
        binding: &Origin,
        divisions: u32,
        visit: &mut impl FnMut(RPoint3) -> bool,
    ) -> Result<bool, Refusal> {
        let source = self.source_origin(binding).ok_or_else(|| {
            Refusal::new(
                "field.draft-sample-origin",
                "result origin has no cached draft definition",
            )
        })?;
        let mut dec = Decider::default();
        match &self.cell {
            Cell::Planar(cell) => {
                for (support, cycle) in cell.supports.iter().zip(&cell.cycles) {
                    if cycle.is_empty() || !support.origins.contains(source) {
                        continue;
                    }
                    let vertices = cycle
                        .iter()
                        .map(|i| native_point(&cell.vertices[*i], &mut dec))
                        .collect::<Result<Vec<_>, _>>()?;
                    let denominator = Rational::from_i64(i64::from(2 * divisions + 1));
                    for triangle in 1..vertices.len() - 1 {
                        let u = vertices[triangle].minus(&vertices[0]);
                        let v = vertices[triangle + 1].minus(&vertices[0]);
                        for i in 0..divisions {
                            for j in 0..divisions - i {
                                let a =
                                    Rational::from_i64(i64::from(2 * i + 1)) / denominator.clone();
                                let b =
                                    Rational::from_i64(i64::from(2 * j + 1)) / denominator.clone();
                                if !visit(vertices[0].plus(&u.scaled(&a)).plus(&v.scaled(&b))) {
                                    return Ok(false);
                                }
                            }
                        }
                    }
                    return Ok(true);
                }
            }
            Cell::Polyhedral(cell) => {
                let mut matched = false;
                for face in &cell.faces {
                    if !face.origins.contains(source) {
                        continue;
                    }
                    matched = true;
                    for cycle in &face.tiles {
                        let vertices = cycle
                            .iter()
                            .map(|i| native_point(&cell.vertices[*i], &mut dec))
                            .collect::<Result<Vec<_>, _>>()?;
                        let denominator = Rational::from_i64(i64::from(2 * divisions + 1));
                        for triangle in 1..vertices.len() - 1 {
                            let u = vertices[triangle].minus(&vertices[0]);
                            let v = vertices[triangle + 1].minus(&vertices[0]);
                            for i in 0..divisions {
                                for j in 0..divisions - i {
                                    let a = Rational::from_i64(i64::from(2 * i + 1))
                                        / denominator.clone();
                                    let b = Rational::from_i64(i64::from(2 * j + 1))
                                        / denominator.clone();
                                    if !visit(vertices[0].plus(&u.scaled(&a)).plus(&v.scaled(&b))) {
                                        return Ok(false);
                                    }
                                }
                            }
                        }
                    }
                }
                if matched {
                    return Ok(true);
                }
            }
            Cell::Round(cell) => {
                let axis = native_vector(&cell.axis, &mut dec)?;
                let reference = cell
                    .reference
                    .scaled(&(&Real::one() / &cell.reference.norm_squared().sqrt()));
                let u = native_vector(&reference, &mut dec)?;
                let v = native_vector(&cell.axis.cross(&reference), &mut dec)?;
                let origin = native_point(&cell.origin, &mut dec)?;
                let heights = [
                    native_scalar(&cell.heights[0], &mut dec)?,
                    native_scalar(&cell.heights[1], &mut dec)?,
                ];
                let radii = [
                    native_scalar(&cell.radii[0], &mut dec)?,
                    native_scalar(&cell.radii[1], &mut dec)?,
                ];
                let angular = circle(divisions);
                if cell.origins.contains(source) {
                    for i in 0..divisions {
                        let t = fraction(i, divisions);
                        let height = &heights[0] + &(&t * &(&heights[1] - &heights[0]));
                        let radius = &radii[0] + &(&t * &(&radii[1] - &radii[0]));
                        let centre = origin.plus(&axis.scaled(&height));
                        for [c, s] in &angular {
                            if !visit(
                                centre
                                    .plus(&u.scaled(&(&radius * c)))
                                    .plus(&v.scaled(&(&radius * s))),
                            ) {
                                return Ok(false);
                            }
                        }
                    }
                    return Ok(true);
                }
                for i in 0..2 {
                    let Some((_, origins)) = &cell.caps[i] else {
                        continue;
                    };
                    if !origins.contains(source) {
                        continue;
                    }
                    let centre = origin.plus(&axis.scaled(&heights[i]));
                    for ring in 0..divisions {
                        let radius = &radii[i] * &fraction(ring, divisions);
                        for [c, s] in &angular {
                            if !visit(
                                centre
                                    .plus(&u.scaled(&(&radius * c)))
                                    .plus(&v.scaled(&(&radius * s))),
                            ) {
                                return Ok(false);
                            }
                        }
                    }
                    return Ok(true);
                }
            }
        }
        Err(Refusal::new(
            "field.draft-sample-role",
            "result origin has no analytic draft patch",
        ))
    }
}

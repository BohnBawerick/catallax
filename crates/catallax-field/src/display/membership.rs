//! Rational trim predicates for display witnesses. Unsupported profiles keep the
//! shared exact-real path; these caches never replace the feature definition.
use crate::{Location, geometry, profile};
use catallax_numbers::{Rational, Sign};
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

pub(crate) struct Classifier<'a> {
    matrix: [[Rational; 3]; 3],
    translation: [Rational; 3],
    vector: &'a RVector3,
    through: bool,
    loops: &'a [Vec<profile::Edge>],
    sector: Option<Sector<'a>>,
    sector_projection: Option<SectorProjection>,
    parameter_bounds: [Rational; 2],
}
impl<'a> Classifier<'a> {
    pub(crate) fn new(
        inverse: &Isometry,
        vector: &'a RVector3,
        through: bool,
        loops: &'a [Vec<profile::Edge>],
    ) -> Option<Self> {
        let sector = if let [edges] = loops {
            Sector::new(edges)
        } else {
            None
        };
        if sector.is_none()
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
}

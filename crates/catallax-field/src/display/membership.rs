//! Rational trim predicates for display witnesses. Unsupported profiles keep the
//! shared exact-real path; these caches never replace the feature definition.
use crate::{Location, geometry, profile};
use catallax_numbers::Rational;
use catallax_recipe::{Role, Turn};
use catallax_shapes::{Isometry, RVector3, Vector3};

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

struct Sector<'a> {
    centre: &'a P,
    squared: &'a Rational,
    start: P,
    end: P,
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
        let (start, end) = if *turn == Turn::Ccw {
            (&arc.start, &arc.end)
        } else {
            (&arc.end, &arc.start)
        };
        Some(Self {
            centre,
            squared,
            start: sub(start, centre),
            end: sub(end, centre),
        })
    }
    fn angle_contains(&self, v: &P) -> bool {
        let a = cross(&self.start, v);
        let b = cross(v, &self.end);
        if cross(&self.start, &self.end) < q(0) {
            a >= q(0) || b >= q(0)
        } else {
            a >= q(0) && b >= q(0)
        }
    }
    fn location(&self, p: &P) -> Location {
        let value = circle_value(p, self.centre, self.squared);
        let v = sub(p, self.centre);
        if value > q(0) || !self.angle_contains(&v) {
            return Location::Outside;
        }
        let on_ray =
            |ray: &P| cross(ray, &v).is_zero() && (&ray[0] * &v[0] + &ray[1] * &v[1]) >= q(0);
        if value.is_zero() || on_ray(&self.start) || on_ray(&self.end) {
            Location::Boundary
        } else {
            Location::Inside
        }
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
        Some(Self {
            matrix,
            translation: [
                inverse.translation.x.as_rational()?,
                inverse.translation.y.as_rational()?,
                inverse.translation.z.as_rational()?,
            ],
            vector,
            through,
            loops,
            sector,
        })
    }
    fn project(&self, p: &[Rational; 3]) -> (Rational, P) {
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
        if !self.through && (t < q(0) || t > q(1)) {
            return Location::Outside;
        }
        let location = self.profile_location(&p);
        if location == Location::Inside && !self.through && (t == q(0) || t == q(1)) {
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
        if !self.through && (t < q(0) || t > q(1)) {
            return false;
        }
        let mut role = role;
        while let Role::Copy { of, .. } = role {
            role = &of.role;
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

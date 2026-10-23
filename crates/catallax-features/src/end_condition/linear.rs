use super::*;

#[derive(Clone, Copy, Debug)]
pub(super) struct Wall {
    pub plane: usize,
    pub negative: bool,
}
#[derive(Clone, Debug)]
enum Node {
    Convex(Vec<Wall>),
    Sweep {
        profile: Profile,
        vector: RVector3,
        range: [Option<Real>; 2],
    },
    Boolean {
        kind: BooleanKind,
        a: usize,
        b: usize,
    },
}
#[derive(Clone, Debug)]
struct Source {
    plane: usize,
    owner: usize,
    origins: Vec<Origin>,
}
#[derive(Clone, Debug)]
pub(super) struct Scene {
    pub planes: Vec<Plane>,
    nodes: Vec<Node>,
    sources: Vec<Source>,
    pub root: usize,
    pub bound: Rational,
}
#[derive(Clone, Debug)]
pub(super) struct PolyFace {
    pub wall: Wall,
    pub vertices: Vec<usize>,
}
#[derive(Clone, Debug)]
struct Poly {
    faces: Vec<PolyFace>,
}
#[derive(Clone, Debug)]
pub(super) struct Patch {
    pub plane: usize,
    pub negative: bool,
    pub vertices: Vec<usize>,
    pub origins: Vec<Origin>,
    pub group: usize,
}
#[derive(Clone, Debug)]
pub(super) struct Boundary {
    pub points: Vec<Point3>,
    pub patches: Vec<Patch>,
    pub groups: Vec<Vec<usize>>,
    pub touching: Vec<BTreeSet<usize>>,
    cells: Vec<Vec<Wall>>,
}
fn value(plane: &Plane, p: &Point3) -> Real {
    p.minus(&plane.origin).dot(&plane.normal)
}
fn origin(step: &StepId, role: Role) -> Origin {
    Origin {
        step: step.clone(),
        role,
    }
}
fn side(segment: &catallax_recipe::ElementId) -> Role {
    Role::Side {
        segment: segment.clone(),
        path: None,
    }
}
fn average(points: &[Point3], ids: &[usize]) -> Point3 {
    let d = Real::from_i64(ids.len() as i64);
    let v = ids.iter().fold(Vector3::zero(), |a, i| {
        a.plus(&points[*i].minus(&Point3::origin()))
    });
    Point3::origin().plus(&v.scaled(&(&r(1) / &d)))
}
impl Scene {
    pub fn compile(recipe: &Recipe, target: &StepId, dec: &mut Decider) -> Result<Self, Refusal> {
        let mut scene = Self {
            planes: vec![],
            nodes: vec![],
            sources: vec![],
            root: 0,
            bound: recipe_bound(recipe, target, dec)?,
        };
        let mut known = BTreeMap::new();
        scene.root = scene.step(recipe, target, &mut known, dec)?;
        Ok(scene)
    }
    fn add_plane(&mut self, plane: Plane, dec: &mut Decider) -> Result<Wall, Refusal> {
        let mut coeff = [
            plane.normal.x.clone(),
            plane.normal.y.clone(),
            plane.normal.z.clone(),
            -&plane.origin.minus(&Point3::origin()).dot(&plane.normal),
        ];
        let mut pivot = None;
        for c in &coeff[..3] {
            if sign(c, dec)? != Sign::Zero {
                pivot = Some(c.clone());
                break;
            }
        }
        let pivot =
            pivot.ok_or_else(|| bad("extrude.zero-thickness", "limit plane normal is zero"))?;
        let negative = sign(&pivot, dec)? == Sign::Positive;
        for c in &mut coeff {
            *c = &*c / &pivot;
        }
        for (i, p) in self.planes.iter().enumerate() {
            let old = [
                p.normal.x.clone(),
                p.normal.y.clone(),
                p.normal.z.clone(),
                -&p.origin.minus(&Point3::origin()).dot(&p.normal),
            ];
            let mut equal = true;
            for j in 0..4 {
                if sign(&(&old[j] - &coeff[j]), dec)? != Sign::Zero {
                    equal = false;
                    break;
                }
            }
            if equal {
                return Ok(Wall { plane: i, negative });
            }
        }
        let normal = Vector3::new(coeff[0].clone(), coeff[1].clone(), coeff[2].clone());
        let origin = Point3::origin().plus(&normal.scaled(&(-&coeff[3] / normal.norm_squared())));
        let basis = if sign(&(&normal.x.square() + &normal.y.square()), dec)? == Sign::Zero {
            Vector3::new(r(1), r(0), r(0))
        } else {
            Vector3::new(r(0), r(0), r(1))
        };
        let reference = normal.cross(&basis);
        let plane = self.planes.len();
        self.planes.push(Plane {
            origin,
            normal,
            ref_dir: reference,
        });
        Ok(Wall { plane, negative })
    }
    fn node(&mut self, node: Node) -> usize {
        let id = self.nodes.len();
        self.nodes.push(node);
        id
    }
    fn step(
        &mut self,
        recipe: &Recipe,
        id: &StepId,
        known: &mut BTreeMap<StepId, usize>,
        dec: &mut Decider,
    ) -> Result<usize, Refusal> {
        if let Some(node) = known.get(id) {
            return Ok(*node);
        }
        let step = recipe
            .step(id)
            .ok_or_else(|| bad("recipe.input-failed", "limit input is missing"))?;
        let result = match &step.op {
            Operation::Block { size, pose } => {
                if size.iter().any(|n| n.value().sign() != Sign::Positive) {
                    return Err(bad(
                        "topology.non-positive-size",
                        "block sizes must be positive",
                    ));
                }
                let pose = crate::profile::placement(pose, dec)?;
                let mut walls = Vec::with_capacity(6);
                let mut roles = Vec::with_capacity(6);
                for axis in 0..3 {
                    for upper in [false, true] {
                        let mut q = [Rational::zero(), Rational::zero(), Rational::zero()];
                        if upper {
                            q[axis] = size[axis].value().clone();
                        }
                        let [x, y, z] = q;
                        let origin = pose.apply_point(&RPoint3::new(x, y, z).to_real());
                        let [x, y, z] = std::array::from_fn(|i| {
                            Rational::from_i64(if i == axis {
                                if upper { 1 } else { -1 }
                            } else {
                                0
                            })
                        });
                        let normal = pose.apply_vector(&RVector3::new(x, y, z).to_real());
                        let [x, y, z] = std::array::from_fn(|i| {
                            Rational::from_i64(i64::from(i == (axis + 1) % 3))
                        });
                        let reference = pose.apply_vector(&RVector3::new(x, y, z).to_real());
                        walls.push(self.add_plane(
                            Plane {
                                origin,
                                normal,
                                ref_dir: reference,
                            },
                            dec,
                        )?);
                        roles.push(if axis == 2 {
                            if upper { Role::End } else { Role::Start }
                        } else {
                            Role::Side {
                                segment: catallax_recipe::ElementId::new(match (axis, upper) {
                                    (0, false) => "x_min",
                                    (0, true) => "x_max",
                                    (1, false) => "y_min",
                                    _ => "y_max",
                                })
                                .unwrap(),
                                path: None,
                            }
                        });
                    }
                }
                let node = self.node(Node::Convex(walls.clone()));
                for (wall, role) in walls.into_iter().zip(roles) {
                    self.sources.push(Source {
                        plane: wall.plane,
                        owner: node,
                        origins: vec![origin(id, role)],
                    });
                }
                node
            }
            Operation::Extrude {
                profile,
                distance,
                direction,
            } => {
                if distance.value().is_zero() {
                    return Err(bad(
                        "extrude.zero-distance",
                        "extrusion distance must not be zero",
                    ));
                }
                let profile = read_profile(recipe, profile, dec)?;
                let normal = crate::profile::rational_vector(
                    &profile.pose,
                    &RVector3::new(Rational::zero(), Rational::zero(), Rational::one()),
                );
                let direction = direction.as_ref().map(vector).unwrap_or(normal);
                if direction.is_zero() {
                    return Err(bad(
                        "extrude.zero-direction",
                        "extrusion direction must not be zero",
                    ));
                }
                let length = crate::extrude::sqrt(&direction.norm_squared()).ok_or_else(|| {
                    bad(
                        "features.unsupported-algebraic-placement",
                        "planar limit sweep needs its exact normalised frame",
                    )
                })?;
                self.sweep(
                    profile,
                    direction.scaled(&(distance.value() / &length)),
                    [Some(r(0)), Some(r(1))],
                    [Some(Role::Start), Some(Role::End)],
                    id,
                    dec,
                )?
            }
            Operation::ExtrudeThrough {
                target,
                profile,
                direction,
                side: which,
            } => {
                let a = self.step(recipe, target, known, dec)?;
                let profile = read_profile(recipe, profile, dec)?;
                let normal = crate::profile::rational_vector(
                    &profile.pose,
                    &RVector3::new(Rational::zero(), Rational::zero(), Rational::one()),
                );
                let direction = direction.as_ref().map(vector).unwrap_or(normal);
                let (range, caps) = match which.unwrap_or(ThroughSide::Both) {
                    ThroughSide::Both => ([None, None], [None, None]),
                    ThroughSide::Forward => ([Some(r(0)), None], [Some(Role::Start), None]),
                    ThroughSide::Backward => ([None, Some(r(0))], [None, Some(Role::Start)]),
                };
                let b = self.sweep(profile, direction, range, caps, id, dec)?;
                self.node(Node::Boolean {
                    kind: BooleanKind::Subtract,
                    a,
                    b,
                })
            }
            Operation::Boolean { kind, target, tool } => {
                let a = self.step(recipe, target, known, dec)?;
                let b = self.step(recipe, tool, known, dec)?;
                self.node(Node::Boolean { kind: *kind, a, b })
            }
            Operation::ExtrudeTo {
                target,
                profile,
                direction,
                limit,
            } => self.import_cap(
                super::Capped::compile(
                    recipe,
                    target,
                    profile,
                    direction.as_ref(),
                    Some(limit),
                    dec,
                )?,
                id,
                dec,
            )?,
            Operation::ExtrudeNext {
                target,
                profile,
                direction,
            } => self.import_cap(
                super::Capped::compile(recipe, target, profile, direction.as_ref(), None, dec)?,
                id,
                dec,
            )?,
            _ => {
                return Err(bad(
                    "features.unsupported-limit-surface",
                    "first-hit target requires a non-planar boundary carrier",
                ));
            }
        };
        known.insert(id.clone(), result);
        Ok(result)
    }
    fn sweep(
        &mut self,
        profile: Profile,
        vector: RVector3,
        range: [Option<Real>; 2],
        caps: [Option<Role>; 2],
        step: &StepId,
        dec: &mut Decider,
    ) -> Result<usize, Refusal> {
        if vector.is_zero() {
            return Err(bad(
                "extrude.zero-direction",
                "extrusion direction must not be zero",
            ));
        }
        let normal = crate::profile::rational_vector(
            &profile.pose,
            &RVector3::new(Rational::zero(), Rational::zero(), Rational::one()),
        );
        if vector.dot(&normal).is_zero() {
            return Err(bad(
                "extrude.zero-thickness",
                "extrusion direction lies in the sketch plane",
            ));
        }
        if profile
            .loops
            .iter()
            .flat_map(|b| &b.segments)
            .any(|s| s.circle.is_some())
        {
            return Err(bad(
                "features.unsupported-limit-surface",
                "a curved target wall needs a non-planar first-hit carrier",
            ));
        }
        let node = self.node(Node::Sweep {
            profile: profile.clone(),
            vector: vector.clone(),
            range: range.clone(),
        });
        for edge in profile.loops.iter().flat_map(|b| &b.segments) {
            let a = profile.world(&edge.start, &zero_vec());
            let b = profile.world(&edge.end, &zero_vec());
            let plane = Plane {
                origin: a.clone().into(),
                normal: b.minus(&a).cross(&vector).into(),
                ref_dir: vector.clone().into(),
            };
            let wall = self.add_plane(plane, dec)?;
            self.sources.push(Source {
                plane: wall.plane,
                owner: node,
                origins: edge.ids.iter().map(|id| origin(step, side(id))).collect(),
            });
        }
        let normal = profile.pose.apply_vector(
            &RVector3::new(Rational::zero(), Rational::zero(), Rational::one()).to_real(),
        );
        let reference = profile.pose.apply_vector(
            &RVector3::new(Rational::one(), Rational::zero(), Rational::zero()).to_real(),
        );
        for i in 0..2 {
            if let Some(depth) = &range[i] {
                let plane = Plane {
                    origin: profile
                        .pose
                        .apply_point(&Point3::origin())
                        .plus(&vector.to_real().scaled(depth)),
                    normal: normal.clone(),
                    ref_dir: reference.clone(),
                };
                let wall = self.add_plane(plane, dec)?;
                if let Some(role) = &caps[i] {
                    self.sources.push(Source {
                        plane: wall.plane,
                        owner: node,
                        origins: vec![origin(step, role.clone())],
                    });
                }
            }
        }
        Ok(node)
    }
    fn import_cap(
        &mut self,
        cap: super::Capped,
        step: &StepId,
        dec: &mut Decider,
    ) -> Result<usize, Refusal> {
        let profile = cap.profile.clone();
        if profile
            .loops
            .iter()
            .flat_map(|b| &b.segments)
            .any(|s| s.circle.is_some())
        {
            return Err(bad(
                "features.unsupported-limit-surface",
                "a circular capped target has a curved wall",
            ));
        }
        let mut joined = None;
        let mut source_planes = Vec::new();
        for cell in &cap.cells {
            let mut walls = Vec::new();
            for (plane, role) in cap.cell_planes(cell, dec)? {
                let wall = self.add_plane(plane, dec)?;
                walls.push(wall);
                if let Some(role) = role {
                    source_planes.push((wall.plane, role));
                }
            }
            let b = self.node(Node::Convex(walls));
            joined = Some(if let Some(a) = joined {
                self.node(Node::Boolean {
                    kind: BooleanKind::Union,
                    a,
                    b,
                })
            } else {
                b
            });
        }
        let tool = self.sweep(
            profile,
            cap.direction.clone(),
            [None, None],
            [None, None],
            step,
            dec,
        )?;
        let root = if let Some(envelope) = joined {
            self.node(Node::Boolean {
                kind: BooleanKind::Intersect,
                a: envelope,
                b: tool,
            })
        } else {
            self.node(Node::Boolean {
                kind: BooleanKind::Subtract,
                a: tool,
                b: tool,
            })
        };
        for (plane, role) in source_planes {
            self.sources.push(Source {
                plane,
                owner: root,
                origins: vec![origin(step, role)],
            });
        }
        Ok(root)
    }
    fn evaluate(
        &self,
        node: usize,
        p: &Point3,
        germ: Option<&Vector3>,
        dec: &mut Decider,
    ) -> Result<bool, Refusal> {
        match &self.nodes[node] {
            Node::Convex(walls) => {
                for w in walls {
                    let plane = &self.planes[w.plane];
                    let mut s = sign(&value(plane, p), dec)?;
                    if s == Sign::Zero {
                        if let Some(d) = germ {
                            s = sign(&plane.normal.dot(d), dec)?;
                        }
                    }
                    if (w.negative && s == Sign::Positive) || (!w.negative && s == Sign::Negative) {
                        return Ok(false);
                    }
                }
                Ok(true)
            }
            Node::Sweep {
                profile,
                vector,
                range,
            } => {
                let mut p = p.clone();
                if let Some(d) = germ {
                    let mut epsilon = Real::one();
                    for plane in &self.planes {
                        let v = value(plane, &p);
                        if sign(&v, dec)? == Sign::Zero {
                            continue;
                        }
                        let v = if sign(&v, dec)? == Sign::Negative {
                            -v
                        } else {
                            v
                        };
                        let slope = plane.normal.dot(d);
                        let slope = if sign(&slope, dec)? == Sign::Negative {
                            -slope
                        } else {
                            slope
                        };
                        let bound = &v / &(&slope + &r(1));
                        if sign(&(&bound - &epsilon), dec)? == Sign::Negative {
                            epsilon = bound;
                        }
                    }
                    p = p.plus(&d.scaled(&(&epsilon / &r(2))));
                }
                let local = profile.pose.inverse().apply_point(&p);
                let v = crate::profile::rational_vector(&profile.pose.inverse(), vector);
                let t = &local.z / &Real::from(v.z.clone());
                if let Some(lo) = &range[0] {
                    if sign(&(&t - lo), dec)? == Sign::Negative {
                        return Ok(false);
                    }
                }
                if let Some(hi) = &range[1] {
                    if sign(&(&t - hi), dec)? == Sign::Positive {
                        return Ok(false);
                    }
                }
                let q = [
                    &local.x - &(&t * &Real::from(v.x)),
                    &local.y - &(&t * &Real::from(v.y)),
                ];
                Ok(profile_location(profile, &q, dec)? != RegionLocation::Outside)
            }
            Node::Boolean { kind, a, b } => {
                let a = self.evaluate(*a, p, germ, dec)?;
                Ok(match kind {
                    BooleanKind::Union => a || self.evaluate(*b, p, germ, dec)?,
                    BooleanKind::Intersect => a && self.evaluate(*b, p, germ, dec)?,
                    BooleanKind::Subtract => a && !self.evaluate(*b, p, germ, dec)?,
                })
            }
        }
    }
}
pub(super) fn read_profile(
    recipe: &Recipe,
    id: &StepId,
    dec: &mut Decider,
) -> Result<Profile, Refusal> {
    let Some(step) = recipe.step(id) else {
        return Err(bad("recipe.input-failed", "profile is missing"));
    };
    let Operation::Sketch { plane, loops } = &step.op else {
        return Err(bad("recipe.input-failed", "profile is not a sketch"));
    };
    Profile::new(plane, loops, dec)
}
pub(super) fn recipe_bound(
    recipe: &Recipe,
    id: &StepId,
    dec: &mut Decider,
) -> Result<Rational, Refusal> {
    let step = recipe
        .step(id)
        .ok_or_else(|| bad("recipe.input-failed", "bound input is missing"))?;
    let stock = |local: Rational, pose: &catallax_recipe::Pose| {
        local
            + pose
                .translation
                .0
                .iter()
                .fold(Rational::zero(), |a, n| a + n.value().abs())
    };
    Ok(match &step.op {
        Operation::Block { size, pose } => stock(
            size.iter()
                .fold(Rational::zero(), |a, n| a + n.value().abs()),
            pose,
        ),
        Operation::Cylinder {
            radius,
            height,
            pose,
        } => stock(radius.value() + height.value(), pose),
        Operation::Sphere { radius, pose } => stock(radius.value().clone(), pose),
        Operation::Cone {
            base_radius,
            top_radius,
            height,
            pose,
        } => stock(
            base_radius.value().clone().max(top_radius.value().clone()) + height.value().clone(),
            pose,
        ),
        Operation::Torus { major, minor, pose } => stock(major.value() + minor.value(), pose),
        Operation::Extrude {
            profile, distance, ..
        } => {
            crate::evaluate::profile_radius(&read_profile(recipe, profile, dec)?)
                + distance.value().abs()
        }
        Operation::Revolve { profile, axis, .. } => {
            crate::evaluate::profile_radius(&read_profile(recipe, profile, dec)?)
                + axis.point.0.iter().fold(Rational::zero(), |a, n| {
                    a + Rational::from_i64(2) * n.value().abs()
                })
        }
        Operation::ExtrudeThrough { target, .. } | Operation::Hole { target, .. } => {
            recipe_bound(recipe, target, dec)?
        }
        Operation::ExtrudeTo {
            target, profile, ..
        }
        | Operation::ExtrudeNext {
            target, profile, ..
        } => recipe_bound(recipe, target, dec)?.max(crate::evaluate::profile_radius(
            &read_profile(recipe, profile, dec)?,
        )),
        Operation::Boolean { kind, target, tool } => {
            let a = recipe_bound(recipe, target, dec)?;
            match kind {
                BooleanKind::Subtract => a,
                BooleanKind::Union => a.max(recipe_bound(recipe, tool, dec)?),
                BooleanKind::Intersect => a.min(recipe_bound(recipe, tool, dec)?),
            }
        }
        _ => {
            return Err(bad(
                "features.unsupported-limit-surface",
                "target bound requires an unfinished operation",
            ));
        }
    })
}

fn polygon_normal(points: &[Point3], ids: &[usize]) -> Vector3 {
    let a = &points[ids[0]];
    (1..ids.len() - 1).fold(Vector3::zero(), |s, i| {
        s.plus(&points[ids[i]].minus(a).cross(&points[ids[i + 1]].minus(a)))
    })
}
fn clip3(
    poly: &Poly,
    wall: Wall,
    planes: &[Plane],
    points: &mut Vec<Point3>,
    values: &mut Vec<Real>,
    signs: &mut Vec<Sign>,
    cache: &mut BTreeMap<(usize, usize), usize>,
    dec: &mut Decider,
) -> Result<Poly, Refusal> {
    let keep = |s| {
        if wall.negative {
            s != Sign::Positive
        } else {
            s != Sign::Negative
        }
    };
    let mut faces = Vec::with_capacity(poly.faces.len() + 1);
    let mut section = BTreeSet::new();
    for face in &poly.faces {
        let mut ids = Vec::with_capacity(face.vertices.len() + 1);
        for i in 0..face.vertices.len() {
            let a = face.vertices[i];
            let b = face.vertices[(i + 1) % face.vertices.len()];
            if keep(signs[a]) {
                ids.push(a);
            }
            if (signs[a] == Sign::Negative && signs[b] == Sign::Positive)
                || (signs[a] == Sign::Positive && signs[b] == Sign::Negative)
            {
                let key = (a.min(b), a.max(b));
                let id = if let Some(id) = cache.get(&key) {
                    *id
                } else {
                    let (a, b) = key;
                    let t = &values[a] / &(&values[a] - &values[b]);
                    let p = points[a].plus(&points[b].minus(&points[a]).scaled(&t));
                    let id = points.len();
                    points.push(p);
                    values.push(r(0));
                    signs.push(Sign::Zero);
                    cache.insert(key, id);
                    id
                };
                ids.push(id);
            }
        }
        ids.dedup();
        if ids.len() > 1 && ids.first() == ids.last() {
            ids.pop();
        }
        if ids.len() < 3 || sign(&polygon_normal(points, &ids).norm_squared(), dec)? == Sign::Zero {
            continue;
        }
        for i in 0..ids.len() {
            let a = ids[i];
            let b = ids[(i + 1) % ids.len()];
            if signs[a] == Sign::Zero && signs[b] == Sign::Zero {
                section.insert((a.min(b), a.max(b)));
            }
        }
        faces.push(PolyFace {
            wall: face.wall,
            vertices: ids,
        });
    }
    let mut adjacent: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
    for (a, b) in section {
        adjacent.entry(a).or_default().push(b);
        adjacent.entry(b).or_default().push(a);
    }
    if !adjacent.is_empty() {
        if adjacent.values().any(|v| v.len() != 2) {
            return Err(bad(
                "features.limit-arrangement",
                "a plane section did not form a simple cycle",
            ));
        }
        let start = *adjacent.keys().next().unwrap();
        let mut ids = vec![start];
        let mut previous = None;
        let mut current = start;
        loop {
            let next = adjacent[&current]
                .iter()
                .copied()
                .find(|v| Some(*v) != previous)
                .unwrap();
            if next == start {
                break;
            }
            if ids.contains(&next) {
                return Err(bad(
                    "features.limit-arrangement",
                    "a plane section cycle repeated a vertex",
                ));
            }
            ids.push(next);
            previous = Some(current);
            current = next;
        }
        let outward = planes[wall.plane]
            .normal
            .scaled(&r(if wall.negative { 1 } else { -1 }));
        if sign(&polygon_normal(points, &ids).dot(&outward), dec)? == Sign::Negative {
            ids.reverse();
        }
        faces.push(PolyFace {
            wall,
            vertices: ids,
        });
    }
    Ok(Poly { faces })
}
fn root(parents: &mut [usize], mut i: usize) -> usize {
    while parents[i] != i {
        parents[i] = parents[parents[i]];
        i = parents[i];
    }
    i
}
fn merge(parents: &mut [usize], a: usize, b: usize) {
    let a = root(parents, a);
    let b = root(parents, b);
    parents[a.max(b)] = a.min(b);
}

impl Scene {
    pub fn boundary(&mut self, dec: &mut Decider) -> Result<Boundary, Refusal> {
        let count = self.planes.len();
        let bound = &self.bound + &Rational::one();
        let mut points: Vec<_> = (0..8)
            .map(|i| {
                let [x, y, z] = std::array::from_fn(|axis| {
                    Real::from(if i & (4 >> axis) == 0 {
                        -&bound
                    } else {
                        bound.clone()
                    })
                });
                Point3::new(x, y, z)
            })
            .collect();
        let cycles = [
            [0, 1, 3, 2],
            [4, 6, 7, 5],
            [0, 4, 5, 1],
            [2, 3, 7, 6],
            [0, 2, 6, 4],
            [1, 5, 7, 3],
        ];
        let mut faces = Vec::with_capacity(6);
        for (i, cycle) in cycles.into_iter().enumerate() {
            let axis = i / 2;
            let mut p = [r(0), r(0), r(0)];
            p[axis] = Real::from(if i % 2 == 0 { -&bound } else { bound.clone() });
            let [x, y, z] = p;
            let [nx, ny, nz] = std::array::from_fn(|j| {
                r(if j == axis {
                    if i % 2 == 0 { -1 } else { 1 }
                } else {
                    0
                })
            });
            let [rx, ry, rz] = std::array::from_fn(|j| r(i64::from(j == (axis + 1) % 3)));
            let wall = self.add_plane(
                Plane {
                    origin: Point3::new(x, y, z),
                    normal: Vector3::new(nx, ny, nz),
                    ref_dir: Vector3::new(rx, ry, rz),
                },
                dec,
            )?;
            faces.push(PolyFace {
                wall,
                vertices: cycle.to_vec(),
            });
        }
        let mut cells = vec![Poly { faces }];
        for plane in 0..count {
            let mut values: Vec<_> = points
                .iter()
                .map(|p| value(&self.planes[plane], p))
                .collect();
            let mut signs = values
                .iter()
                .map(|v| sign(v, dec))
                .collect::<Result<Vec<_>, _>>()?;
            let mut cache = BTreeMap::new();
            let mut split = Vec::new();
            for cell in cells {
                let mut negative = false;
                let mut positive = false;
                for id in cell.faces.iter().flat_map(|f| &f.vertices) {
                    negative |= signs[*id] == Sign::Negative;
                    positive |= signs[*id] == Sign::Positive;
                }
                if negative && positive {
                    for negative in [true, false] {
                        split.push(clip3(
                            &cell,
                            Wall { plane, negative },
                            &self.planes,
                            &mut points,
                            &mut values,
                            &mut signs,
                            &mut cache,
                            dec,
                        )?);
                    }
                } else {
                    split.push(cell);
                }
            }
            cells = split;
        }
        let mut selected = Vec::new();
        let mut boundary: BTreeMap<(usize, Vec<usize>), PolyFace> = BTreeMap::new();
        for cell in cells {
            let ids: Vec<_> = cell
                .faces
                .iter()
                .flat_map(|f| f.vertices.iter().copied())
                .collect::<BTreeSet<_>>()
                .into_iter()
                .collect();
            if !self.evaluate(self.root, &average(&points, &ids), None, dec)? {
                continue;
            }
            selected.push(cell.faces.iter().map(|f| f.wall).collect());
            for face in cell.faces {
                let mut ids = face.vertices.clone();
                ids.sort_unstable();
                let key = (face.wall.plane, ids);
                if boundary.remove(&key).is_none() {
                    boundary.insert(key, face);
                }
            }
        }
        let mut patches = Vec::with_capacity(boundary.len());
        for (_, face) in boundary {
            let mid = average(&points, &face.vertices);
            let normal = &self.planes[face.wall.plane].normal;
            let mut origins = Vec::new();
            for source in &self.sources {
                if source.plane == face.wall.plane
                    && self.evaluate(source.owner, &mid, Some(normal), dec)?
                        != self.evaluate(source.owner, &mid, Some(&normal.scaled(&r(-1))), dec)?
                {
                    origins.extend(source.origins.iter().cloned());
                }
            }
            origins.sort();
            origins.dedup();
            if origins.is_empty() {
                return Err(bad(
                    "features.limit-arrangement",
                    "a target boundary patch has no defining face origin",
                ));
            }
            patches.push(Patch {
                plane: face.wall.plane,
                negative: face.wall.negative,
                vertices: face.vertices,
                origins,
                group: 0,
            });
        }
        let mut parents: Vec<_> = (0..patches.len()).collect();
        let mut edges = BTreeMap::new();
        for (i, patch) in patches.iter().enumerate() {
            for j in 0..patch.vertices.len() {
                let a = patch.vertices[j];
                let b = patch.vertices[(j + 1) % patch.vertices.len()];
                let key = (patch.plane, patch.negative, a.min(b), a.max(b));
                if let Some(other) = edges.insert(key, i) {
                    merge(&mut parents, i, other);
                }
            }
        }
        let mut group_ids = BTreeMap::new();
        let mut groups: Vec<Vec<usize>> = Vec::new();
        for (i, patch) in patches.iter_mut().enumerate() {
            let p = root(&mut parents, i);
            let next = groups.len();
            let group = *group_ids.entry(p).or_insert(next);
            if group == groups.len() {
                groups.push(vec![]);
            }
            groups[group].push(i);
            patch.group = group;
        }
        for group in &groups {
            let mut origins: Vec<_> = group
                .iter()
                .flat_map(|i| patches[*i].origins.iter().cloned())
                .collect();
            origins.sort();
            origins.dedup();
            for i in group {
                patches[*i].origins = origins.clone();
            }
        }
        let mut edge_groups: BTreeMap<(usize, usize), BTreeSet<usize>> = BTreeMap::new();
        for patch in &patches {
            for j in 0..patch.vertices.len() {
                let a = patch.vertices[j];
                let b = patch.vertices[(j + 1) % patch.vertices.len()];
                edge_groups
                    .entry((a.min(b), a.max(b)))
                    .or_default()
                    .insert(patch.group);
            }
        }
        let mut touching = vec![BTreeSet::new(); groups.len()];
        for groups in edge_groups.values() {
            for a in groups {
                for b in groups {
                    if a != b {
                        touching[*a].insert(*b);
                    }
                }
            }
        }
        Ok(Boundary {
            points,
            patches,
            groups,
            touching,
            cells: selected,
        })
    }
}
fn contains_patch(
    scene: &Scene,
    boundary: &Boundary,
    patch: &Patch,
    p: &Point3,
    dec: &mut Decider,
) -> Result<bool, Refusal> {
    let plane = &scene.planes[patch.plane];
    if sign(&value(plane, p), dec)? != Sign::Zero {
        return Ok(false);
    }
    let normal = plane.normal.scaled(&r(if patch.negative { 1 } else { -1 }));
    for i in 0..patch.vertices.len() {
        let a = &boundary.points[patch.vertices[i]];
        let b = &boundary.points[patch.vertices[(i + 1) % patch.vertices.len()]];
        if sign(&b.minus(a).cross(&p.minus(a)).dot(&normal), dec)? == Sign::Negative {
            return Ok(false);
        }
    }
    Ok(true)
}
impl Boundary {
    pub fn location(
        &self,
        scene: &Scene,
        p: &Point3,
        dec: &mut Decider,
    ) -> Result<RegionLocation, Refusal> {
        for patch in &self.patches {
            if contains_patch(scene, self, patch, p, dec)? {
                return Ok(RegionLocation::Boundary);
            }
        }
        let signs = scene
            .planes
            .iter()
            .map(|plane| sign(&value(plane, p), dec))
            .collect::<Result<Vec<_>, _>>()?;
        for cell in &self.cells {
            if cell.iter().all(|w| {
                if w.negative {
                    signs[w.plane] != Sign::Positive
                } else {
                    signs[w.plane] != Sign::Negative
                }
            }) {
                return Ok(RegionLocation::Inside);
            }
        }
        Ok(RegionLocation::Outside)
    }
    pub fn on_selected(
        &self,
        scene: &Scene,
        selected: &BTreeSet<usize>,
        p: &Point3,
        dec: &mut Decider,
    ) -> Result<bool, Refusal> {
        for patch in &self.patches {
            if selected.contains(&patch.group) && contains_patch(scene, self, patch, p, dec)? {
                return Ok(true);
            }
        }
        Ok(false)
    }
    pub fn resolve(
        &self,
        scene: &Scene,
        reference: &FaceRef,
        dec: &mut Decider,
    ) -> Result<BTreeSet<usize>, Refusal> {
        let matches = |group: usize, pattern: &catallax_recipe::OriginPattern| {
            self.patches[self.groups[group][0]]
                .origins
                .iter()
                .any(|o| pattern.matches(o))
        };
        let mut chosen: Vec<_> = (0..self.groups.len())
            .filter(|g| {
                matches(*g, &reference.face)
                    && reference.touching.iter().all(|pattern| {
                        self.touching[*g]
                            .iter()
                            .any(|other| matches(*other, pattern))
                    })
            })
            .collect();
        if chosen.is_empty() {
            return Err(bad(
                "recipe.missing-reference",
                "the extrusion limit chooses no surviving target face",
            ));
        }
        if let Some(near) = &reference.near {
            let near = point(near);
            let mut best = None;
            let mut tied = false;
            for group in chosen {
                let mut distance = None;
                for i in &self.groups[group] {
                    let patch = &self.patches[*i];
                    let plane = &scene.planes[patch.plane];
                    let v = value(plane, &near);
                    let norm = plane.normal.norm_squared();
                    let projected = near.plus(&plane.normal.scaled(&(&(-&v) / &norm)));
                    let mut candidate = None;
                    if contains_patch(scene, self, patch, &projected, dec)? {
                        candidate = Some(&v.square() / &norm);
                    }
                    for j in 0..patch.vertices.len() {
                        let a = &self.points[patch.vertices[j]];
                        let b = &self.points[patch.vertices[(j + 1) % patch.vertices.len()]];
                        let edge = b.minus(a);
                        let t = &near.minus(a).dot(&edge) / &edge.norm_squared();
                        let t = if sign(&t, dec)? == Sign::Negative {
                            r(0)
                        } else if sign(&(&t - &r(1)), dec)? == Sign::Positive {
                            r(1)
                        } else {
                            t
                        };
                        let d = near.minus(&a.plus(&edge.scaled(&t))).norm_squared();
                        if candidate.as_ref().is_none()
                            || sign(&(&d - candidate.as_ref().unwrap()), dec)? == Sign::Negative
                        {
                            candidate = Some(d);
                        }
                    }
                    let candidate = candidate.unwrap();
                    if distance.as_ref().is_none()
                        || sign(&(&candidate - distance.as_ref().unwrap()), dec)? == Sign::Negative
                    {
                        distance = Some(candidate);
                    }
                }
                let distance = distance.unwrap();
                if let Some((_, old)) = &best {
                    match sign(&(&distance - old), dec)? {
                        Sign::Negative => {
                            best = Some((group, distance));
                            tied = false;
                        }
                        Sign::Zero => tied = true,
                        Sign::Positive => (),
                    }
                } else {
                    best = Some((group, distance));
                }
            }
            if tied {
                return Err(bad(
                    "recipe.ambiguous-reference",
                    "the stored near point ties between extrusion limit faces",
                ));
            }
            chosen = vec![best.unwrap().0];
        }
        Ok(chosen.into_iter().collect())
    }
}

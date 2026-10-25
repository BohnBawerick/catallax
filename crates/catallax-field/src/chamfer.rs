//! Setback edge regions compiled from recipe supports, without B-rep topology.

use catallax_features::chamfer::{
    ChamferDefinition,
    arc::{BoundedArc, ordered_points},
};
use catallax_intersect::{Component, CurveSurface, SurfacePair};
use catallax_numbers::{Decider, Rational, Real, Refusal, Sign};
use catallax_recipe::{BooleanKind, EdgeName, EdgeSelection, ElementId, Origin, Role, StepId};
use catallax_shapes::{Curve, GeomStore, Point3, Surface, Vector3};

use crate::{Instruction, Location, Plan, features::Feature, geometry, profile};

#[derive(Clone, Debug)]
enum Domain {
    Whole,
    Line([Rational; 2]),
    Arc(BoundedArc),
}

pub(crate) struct Setbacks {
    pub distances: [Real; 2],
    pub ordered: bool,
}

fn bad(detail: &str) -> Refusal {
    Refusal::new("field.chamfer-contact", detail)
}
fn q(x: &Real) -> Result<Rational, Refusal> {
    x.as_rational()
        .ok_or_else(|| bad("the section needs algebraic profile coordinates"))
}
fn square_root(x: &Real) -> Real {
    x.as_rational()
        .and_then(|r| profile::sqrt(&r))
        .map(Real::from)
        .unwrap_or_else(|| x.sqrt())
}
fn unit(v: &Vector3) -> Vector3 {
    v.scaled(&(&Real::one() / &square_root(&v.norm_squared())))
}
fn coordinates(p: &Point3) -> Result<[Rational; 3], Refusal> {
    Ok([q(&p.x)?, q(&p.y)?, q(&p.z)?])
}
fn number(x: Rational) -> catallax_recipe::Number {
    x.into()
}
fn point(p: &Point3) -> Result<catallax_recipe::Point3, Refusal> {
    Ok(catallax_recipe::Point3(coordinates(p)?.map(number)))
}
fn vector(v: &Vector3) -> Result<catallax_recipe::Vector3, Refusal> {
    Ok(catallax_recipe::Vector3(
        [q(&v.x)?, q(&v.y)?, q(&v.z)?].map(number),
    ))
}
fn pose(
    origin: &Point3,
    x: &Vector3,
    y: &Vector3,
    z: &Vector3,
    dec: &mut Decider,
) -> Result<catallax_recipe::Pose, Refusal> {
    let m = [
        [x.x.clone(), y.x.clone(), z.x.clone()],
        [x.y.clone(), y.y.clone(), z.y.clone()],
        [x.z.clone(), y.z.clone(), z.z.clone()],
    ];
    let trace = &(&m[0][0] + &m[1][1]) + &m[2][2];
    let one = Real::one();
    let quaternion = if geometry::sign(&(&one + &trace), dec)? == Sign::Positive {
        [
            &one + &trace,
            &m[2][1] - &m[1][2],
            &m[0][2] - &m[2][0],
            &m[1][0] - &m[0][1],
        ]
    } else {
        let mut diagonal = None;
        for (i, row) in m.iter().enumerate() {
            if geometry::sign(&(&one + &row[i]), dec)? == Sign::Positive {
                diagonal = Some(i);
                break;
            }
        }
        let i = diagonal.ok_or_else(|| bad("section frame is not a proper rotation"))?;
        let j = (i + 1) % 3;
        let k = (i + 2) % 3;
        let mut quat = std::array::from_fn(|_| Real::zero());
        quat[0] = &m[k][j] - &m[j][k];
        quat[i + 1] = &(&(&one + &m[i][i]) - &m[j][j]) - &m[k][k];
        quat[j + 1] = &m[j][i] + &m[i][j];
        quat[k + 1] = &m[k][i] + &m[i][k];
        quat
    };
    let [a, b, c, d] = quaternion;
    Ok(catallax_recipe::Pose {
        rotation: [
            number(q(&a)?),
            number(q(&b)?),
            number(q(&c)?),
            number(q(&d)?),
        ],
        translation: catallax_recipe::Vector3(coordinates(origin)?.map(number)),
    })
}

impl Plan {
    fn chamfer_ray(
        &self,
        root: usize,
        p: &[Rational; 3],
        d: &Vector3,
        dec: &mut Decider,
    ) -> Result<Location, Refusal> {
        let mut locations = Vec::with_capacity(root + 1);
        for instruction in &self.instructions[..=root] {
            let location = match *instruction {
                Instruction::Empty => Location::Outside,
                Instruction::Leaf(i) => self.leaves[i].ray(p, d, dec)?,
                Instruction::Boolean { kind, a, b } => {
                    crate::combine(kind, locations[a], locations[b])
                }
            };
            locations.push(location);
        }
        Ok(locations[root])
    }
    fn chamfer_material(
        &self,
        root: usize,
        p: &Point3,
        faces: [usize; 2],
        dec: &mut Decider,
    ) -> Result<Option<(bool, [bool; 2])>, Refusal> {
        let xyz = coordinates(p)?;
        for i in faces {
            let face = &self.faces[i];
            if self.leaves[face.leaf].location(&xyz, dec)? != Location::Boundary
                || self.leaves[face.leaf].on_face(
                    &face.origin,
                    self.definition_role(i),
                    &xyz,
                    dec,
                )? == Some(false)
            {
                return Ok(None);
            }
        }
        let gradients = faces.map(|i| crate::normal(&self.faces[i], &xyz));
        let [a, b] = gradients;
        let dot = a.dot(&b);
        let da = a.minus(&b.scaled(&(&dot / &b.norm_squared())));
        let db = b.minus(&a.scaled(&(&dot / &a.norm_squared())));
        if geometry::sign(&da.norm_squared(), dec)? == Sign::Zero {
            return Ok(None);
        }
        let mut inside = Vec::new();
        let mut outside = Vec::new();
        for signs in [[-1, -1], [-1, 1], [1, -1], [1, 1]] {
            let direction = da
                .scaled(&Real::from_i64(signs[0]))
                .plus(&db.scaled(&Real::from_i64(signs[1])));
            match self.chamfer_ray(root, &xyz, &direction, dec)? {
                Location::Inside => inside.push(signs),
                Location::Outside => outside.push(signs),
                _ => return Ok(None),
            }
        }
        Ok(match (inside.as_slice(), outside.as_slice()) {
            ([signs], outside) if outside.len() == 3 => Some((true, signs.map(|s| s < 0))),
            (inside, [signs]) if inside.len() == 3 => Some((false, signs.map(|s| s > 0))),
            _ => None,
        })
    }
    fn chamfer_domain(
        &self,
        definition: &ChamferDefinition,
        faces: [usize; 2],
        interval: &Domain,
        source_faces: &[usize],
        dec: &mut Decider,
    ) -> Result<Option<([catallax_shapes::Plane; 2], bool)>, Refusal> {
        let mut end_planes = [false; 2];
        let mut arc_planes = [None, None];
        for &i in source_faces {
            let face = &self.faces[i];
            if faces
                .iter()
                .any(|j| face.support.same_definition(&self.faces[*j].support))
            {
                continue;
            }
            if let (Curve::Line(edge), Surface::Plane(p), Domain::Line(interval)) =
                (&definition.edge, &face.surface, interval)
                && geometry::sign(&p.normal.cross(&edge.dir).norm_squared(), dec)? == Sign::Zero
            {
                let t =
                    q(&(p.origin.minus(&edge.origin).dot(&p.normal) / edge.dir.dot(&p.normal)))?;
                if let Some(end) = interval.iter().position(|value| *value == t) {
                    end_planes[end] = true;
                    continue;
                }
            }
            if let (Domain::Arc(arc), Surface::Plane(plane)) = (interval, &face.surface) {
                let mut matched = false;
                for (end, point) in arc.ends.iter().enumerate() {
                    if geometry::sign(&point.minus(&plane.origin).dot(&plane.normal), dec)?
                        == Sign::Zero
                    {
                        for contact in &definition.contacts {
                            let Curve::Circle(contact) = contact else {
                                unreachable!("circular contact")
                            };
                            arc.contact_at(contact, end, plane, dec)?;
                        }
                        end_planes[end] = true;
                        arc_planes[end] = Some(plane.clone());
                        matched = true;
                    }
                }
                if matched {
                    continue;
                }
            }
            if definition.avoids_support(&face.surface, &vertices(definition, interval), dec)? {
                continue;
            }
            return Err(Refusal::new(
                "chamfer.unsupported-overflow",
                "the field affected region may reach another support and needs trimmed continuation",
            ));
        }
        if !matches!(interval, Domain::Whole) && end_planes != [true, true] {
            return Err(bad(
                "finite line strips need proven perpendicular end planes",
            ));
        }
        if let Domain::Arc(arc) = interval {
            return arc
                .trim(
                    definition,
                    &arc_planes.map(|p| p.expect("proven arc end plane")),
                    dec,
                )
                .map(Some);
        }
        Ok(None)
    }

    pub(crate) fn chamfer(
        &mut self,
        root: usize,
        source_faces: &[usize],
        selection: &EdgeSelection,
        setbacks: &Setbacks,
        step: &StepId,
        dec: &mut Decider,
    ) -> Result<usize, Refusal> {
        let distances = &setbacks.distances;
        let unequal = setbacks.ordered;
        for distance in distances {
            if geometry::sign(distance, dec)? != Sign::Positive {
                return Err(Refusal::new(
                    "chamfer.non-positive-distance",
                    "chamfer distances must be positive",
                ));
            }
        }
        if matches!(selection, EdgeSelection::Refs(refs) if refs.is_empty())
            || (matches!(selection, EdgeSelection::All) && source_faces.len() < 2)
        {
            return Ok(root);
        }
        let mut pairs = Vec::new();
        match selection {
            EdgeSelection::All => {
                if unequal {
                    return Err(Refusal::new(
                        "chamfer.ambiguous-order",
                        "two setbacks need ordered references",
                    ));
                }
                for (position, &a) in source_faces.iter().enumerate() {
                    for &b in &source_faces[position + 1..] {
                        pairs.push(([a, b], distances.clone(), false, false));
                    }
                }
            }
            EdgeSelection::Refs(refs) => {
                for reference in refs {
                    if reference.near.is_some() || !reference.touching.is_empty() {
                        return Err(bad(
                            "filtered edge references need field-native endpoint and distance proofs",
                        ));
                    }
                    let mut found = false;
                    for (position, &a) in source_faces.iter().enumerate() {
                        for &b in &source_faces[position + 1..] {
                            let forward = reference.between[0].matches(&self.faces[a].origin)
                                && reference.between[1].matches(&self.faces[b].origin);
                            let backward = reference.between[1].matches(&self.faces[a].origin)
                                && reference.between[0].matches(&self.faces[b].origin);
                            if !forward && !backward {
                                continue;
                            }
                            let ordered = if forward {
                                distances.clone()
                            } else {
                                [distances[1].clone(), distances[0].clone()]
                            };
                            if let Some((_, previous, ambiguous, conflict)) =
                                pairs.iter_mut().find(|(faces, _, _, _)| *faces == [a, b])
                            {
                                *ambiguous |= unequal && forward && backward;
                                for (old, new) in previous.iter().zip(&ordered) {
                                    *conflict |= geometry::sign(&(old - new), dec)? != Sign::Zero;
                                }
                            } else {
                                pairs.push((
                                    [a, b],
                                    ordered,
                                    unequal && forward && backward,
                                    false,
                                ));
                            }
                            found = true;
                        }
                    }
                    if !found {
                        return Err(Refusal::new(
                            "chamfer.lost-edge",
                            "no source support resolves the edge reference",
                        ));
                    }
                }
            }
        }
        let mut bound = Rational::one();
        for leaf in &self.leaves {
            if let Some(r) = leaf.display_radius(dec)? {
                bound = bound + r;
            }
        }
        let mut definitions = Vec::new();
        for (faces, distances, ambiguous, conflict) in pairs {
            let supports = faces.map(|i| self.faces[i].surface.clone());
            let supported = matches!(
                (&supports[0], &supports[1]),
                (Surface::Plane(_), Surface::Plane(_) | Surface::Cylinder(_))
                    | (Surface::Cylinder(_), Surface::Plane(_))
            );
            let mut geom = GeomStore::new();
            let a = geom.add_surface(supports[0].clone());
            let b = geom.add_surface(supports[1].clone());
            let SurfacePair::Meets(components) =
                catallax_intersect::intersect_surfaces(&mut geom, a, b, dec)?.value
            else {
                continue;
            };
            for component in components {
                let Component::Curve { curve, .. } = component else {
                    continue;
                };
                let edge = geom.curve(curve).clone();
                let samples = match &edge {
                    Curve::Circle(c) => {
                        let mut stations = Vec::new();
                        for &other in source_faces {
                            if faces.iter().any(|i| {
                                self.faces[*i]
                                    .support
                                    .same_definition(&self.faces[other].support)
                            }) {
                                continue;
                            }
                            let surface = geom.add_surface(self.faces[other].surface.clone());
                            match catallax_intersect::intersect_curve_surface(
                                &mut geom, curve, surface, dec,
                            )?
                            .value
                            {
                                CurveSurface::Disjoint => (),
                                CurveSurface::Meets(hits) => {
                                    stations
                                        .extend(hits.into_iter().map(|h| geom.position(h.point)));
                                }
                                _ => {
                                    return Err(bad(
                                        "a circular run coincides with another support",
                                    ));
                                }
                            }
                        }
                        let stations = ordered_points(c, stations, dec)?;
                        if stations.is_empty() {
                            vec![(
                                c.center
                                    .plus(&unit(&c.ref_dir).scaled(&square_root(&c.radius_sq))),
                                Domain::Whole,
                            )]
                        } else {
                            if stations.len() < 2 {
                                return Err(bad("a tangent circular end needs continuation"));
                            }
                            let mut samples = Vec::new();
                            for i in 0..stations.len() {
                                let arc = BoundedArc::new(
                                    c.clone(),
                                    [
                                        stations[i].clone(),
                                        stations[(i + 1) % stations.len()].clone(),
                                    ],
                                    dec,
                                )?;
                                samples.push((arc.rational_sample(dec)?, Domain::Arc(arc)));
                            }
                            samples
                        }
                    }
                    Curve::Line(line) => {
                        let length = q(&square_root(&line.dir.norm_squared()))?;
                        let origin_bound = q(&line.origin.x)?.abs()
                            + q(&line.origin.y)?.abs()
                            + q(&line.origin.z)?.abs();
                        let extent = (&bound + &origin_bound) / length;
                        let mut stations = vec![-&extent, extent];
                        for &index in source_faces {
                            let face = &self.faces[index];
                            let surface = geom.add_surface(face.surface.clone());
                            if let CurveSurface::Meets(hits) =
                                catallax_intersect::intersect_curve_surface(
                                    &mut geom, curve, surface, dec,
                                )?
                                .value
                            {
                                for hit in hits {
                                    let p = geom.position(hit.point);
                                    stations.push(q(&(p.minus(&line.origin).dot(&line.dir)
                                        / line.dir.norm_squared()))?);
                                }
                            }
                        }
                        stations.sort();
                        stations.dedup();
                        stations
                            .windows(2)
                            .map(|pair| {
                                (
                                    line.origin.plus(&line.dir.scaled(&Real::from(
                                        (&pair[0] + &pair[1]) / Rational::from_i64(2),
                                    ))),
                                    Domain::Line([pair[0].clone(), pair[1].clone()]),
                                )
                            })
                            .collect()
                    }
                    _ => {
                        return Err(bad(
                            "selected support intersections need a certified line or circle domain",
                        ));
                    }
                };
                let first_definition = definitions.len();
                for (sample, interval) in samples {
                    let Some((convex, senses)) =
                        self.chamfer_material(root, &sample, faces, dec)?
                    else {
                        continue;
                    };
                    if ambiguous {
                        return Err(Refusal::new(
                            "chamfer.ambiguous-order",
                            "an input edge matches its ordered reference both ways",
                        ));
                    }
                    if conflict {
                        return Err(Refusal::new(
                            "chamfer.conflicting-distances",
                            "overlapping references assign different setbacks",
                        ));
                    }
                    if !supported {
                        return Err(bad(
                            "a selected sharp input edge requires an unsupported bevel carrier",
                        ));
                    }
                    if definitions.len() > first_definition {
                        let (previous, previous_interval): &mut (ChamferDefinition, Domain) =
                            definitions.last_mut().unwrap();
                        if let (Domain::Line(previous_interval), Domain::Line(current)) =
                            (&mut *previous_interval, &interval)
                            && previous_interval[1] == current[0]
                            && previous.convex == convex
                            && previous.same_sense == senses
                        {
                            previous_interval[1] = current[1].clone();
                            continue;
                        }
                        if let (Domain::Arc(previous_arc), Domain::Arc(current)) =
                            (&mut *previous_interval, &interval)
                            && geometry::sign(
                                &previous_arc.ends[1].minus(&current.ends[0]).norm_squared(),
                                dec,
                            )? == Sign::Zero
                            && previous.convex == convex
                            && previous.same_sense == senses
                        {
                            if geometry::sign(
                                &previous_arc.ends[0].minus(&current.ends[1]).norm_squared(),
                                dec,
                            )? == Sign::Zero
                            {
                                *previous_interval = Domain::Whole;
                            } else {
                                previous_arc.ends[1] = current.ends[1].clone();
                            }
                            continue;
                        }
                    }
                    let edge_name = EdgeName::new(
                        self.faces[faces[0]].origin.clone(),
                        self.faces[faces[1]].origin.clone(),
                    )
                    .map_err(|e| bad(&e.to_string()))?;
                    let definition = ChamferDefinition::new(
                        Origin {
                            step: step.clone(),
                            role: Role::Bevel { edge: edge_name },
                        },
                        supports.clone(),
                        senses,
                        edge.clone(),
                        distances.clone(),
                        convex,
                        dec,
                    )?;
                    definitions.push((definition, interval));
                }
                // The first and last cells are adjacent on a circle. Joining
                // them preserves a single named run across the carrier seam.
                if definitions.len() > first_definition + 1 {
                    let (first, first_domain) = &definitions[first_definition];
                    let (last, last_domain) = definitions.last().unwrap();
                    if let (Domain::Arc(a), Domain::Arc(b)) = (first_domain, last_domain)
                        && first.convex == last.convex
                        && first.same_sense == last.same_sense
                        && geometry::sign(&b.ends[1].minus(&a.ends[0]).norm_squared(), dec)?
                            == Sign::Zero
                    {
                        let start = b.ends[0].clone();
                        definitions.pop();
                        let Domain::Arc(arc) = &mut definitions[first_definition].1 else {
                            unreachable!("joined circular cells");
                        };
                        arc.ends[0] = start;
                    }
                }
            }
        }
        if definitions.is_empty()
            && !matches!(selection, EdgeSelection::Refs(refs) if refs.is_empty())
        {
            return Err(Refusal::new(
                "chamfer.lost-edge",
                "no selected support intersection is a sharp input boundary",
            ));
        }
        if let EdgeSelection::Refs(refs) = selection {
            for reference in refs {
                if !definitions.iter().any(|(definition, _)| {
                    let Role::Bevel { edge } = &definition.origin.role else {
                        return false;
                    };
                    let [a, b] = edge.between();
                    (reference.between[0].matches(a) && reference.between[1].matches(b))
                        || (reference.between[1].matches(a) && reference.between[0].matches(b))
                }) {
                    return Err(Refusal::new(
                        "chamfer.lost-edge",
                        "an edge reference has no surviving sharp input component",
                    ));
                }
            }
        }
        let mut endpoints: Vec<(Point3, usize)> = Vec::new();
        for (i, (definition, interval)) in definitions.iter().enumerate() {
            if definitions[..i]
                .iter()
                .any(|(other, _)| other.origin == definition.origin)
            {
                return Err(bad(
                    "split named contact runs need a multi-domain face-selection proof",
                ));
            }
            if let (Curve::Line(edge), Domain::Line(interval)) = (&definition.edge, interval) {
                for t in interval {
                    let point = edge.origin.plus(&edge.dir.scaled(&Real::from(t.clone())));
                    let mut found = false;
                    for (other, count) in &mut endpoints {
                        if geometry::sign(&point.minus(other).norm_squared(), dec)? == Sign::Zero {
                            *count += 1;
                            if *count > 2 {
                                return Err(bad(
                                    "three selected runs at a vertex need a corner patch",
                                ));
                            }
                            found = true;
                            break;
                        }
                    }
                    if !found {
                        endpoints.push((point, 1));
                    }
                }
            }
        }
        for (i, (definition, interval)) in definitions.iter().enumerate() {
            for (other, other_interval) in &definitions[i + 1..] {
                let shared_end =
                    if let (Curve::Line(a), Domain::Line(ai), Curve::Line(b), Domain::Line(bi)) =
                        (&definition.edge, interval, &other.edge, other_interval)
                    {
                        let mut shared = false;
                        for at in ai {
                            for bt in bi {
                                let a = a.origin.plus(&a.dir.scaled(&Real::from(at.clone())));
                                let b = b.origin.plus(&b.dir.scaled(&Real::from(bt.clone())));
                                shared |=
                                    geometry::sign(&a.minus(&b).norm_squared(), dec)? == Sign::Zero;
                            }
                        }
                        shared
                    } else {
                        false
                    };
                if shared_end {
                    continue;
                }
                if !definition.avoids_region(
                    other,
                    &vertices(definition, interval),
                    &vertices(other, other_interval),
                    dec,
                )? {
                    return Err(Refusal::new(
                        "chamfer.unsupported-overflow",
                        "nonincident affected regions need a simultaneous overlap proof",
                    ));
                }
            }
        }
        definitions.sort_by_key(|(definition, _)| !definition.convex);
        let mut result = root;
        for (definition, interval) in definitions {
            let faces = [
                self.faces
                    .iter()
                    .position(|f| f.surface.same_definition(&definition.supports[0]))
                    .expect("source support"),
                self.faces
                    .iter()
                    .position(|f| f.surface.same_definition(&definition.supports[1]))
                    .expect("source support"),
            ];
            let mut angular =
                self.chamfer_domain(&definition, faces, &interval, source_faces, dec)?;
            let convex = definition.convex;
            let sweep = finite_sweep(&interval, angular.as_ref(), dec)?;
            if sweep.is_some() {
                angular = None;
            }
            let (feature, faces) = region(definition, interval, sweep, step, dec)?;
            let mut node = self.leaf(geometry::Leaf::Feature(feature), faces, dec)?;
            if let Some((planes, major)) = angular {
                let mut gates = Vec::new();
                for plane in planes {
                    gates.push(self.leaf(
                        geometry::Leaf::Cell(geometry::Cell::new(
                            vec![Surface::Plane(plane)],
                            None,
                        )?),
                        Vec::new(),
                        dec,
                    )?);
                }
                let gate = self.boolean(
                    if major {
                        BooleanKind::Union
                    } else {
                        BooleanKind::Intersect
                    },
                    gates[0],
                    gates[1],
                    dec,
                )?;
                node = self.boolean(BooleanKind::Intersect, node, gate, dec)?;
            }
            result = self.boolean(
                if convex {
                    BooleanKind::Subtract
                } else {
                    BooleanKind::Union
                },
                result,
                node,
                dec,
            )?;
        }
        Ok(result)
    }
}

fn vertices(definition: &ChamferDefinition, interval: &Domain) -> Vec<Point3> {
    let (Curve::Line(edge), [Curve::Line(a), Curve::Line(b)], Domain::Line(interval)) =
        (&definition.edge, &definition.contacts, interval)
    else {
        return Vec::new();
    };
    [edge, a, b]
        .into_iter()
        .flat_map(|line| {
            interval
                .iter()
                .map(|t| line.origin.plus(&edge.dir.scaled(&Real::from(t.clone()))))
        })
        .collect()
}

// Reuse the field's own finite angular sweep when both end planes are radial
// and its existing exact angular classifier covers this domain. Other proven
// end cuts retain the general analytic halfspace construction above.
fn finite_sweep(
    domain: &Domain,
    angular: Option<&([catallax_shapes::Plane; 2], bool)>,
    dec: &mut Decider,
) -> Result<Option<(Vector3, Rational)>, Refusal> {
    let (Domain::Arc(arc), Some((planes, _))) = (domain, angular) else {
        return Ok(None);
    };
    for plane in planes {
        if geometry::sign(
            &arc.circle.center.minus(&plane.origin).dot(&plane.normal),
            dec,
        )? != Sign::Zero
        {
            return Ok(None);
        }
    }
    let [a, b] = arc.ends.each_ref().map(|p| p.minus(&arc.circle.center));
    let dot = geometry::sign(&a.dot(&b), dec)?;
    let cross = geometry::sign(&a.cross(&b).dot(&arc.circle.normal), dec)?;
    let angle = match (dot, cross) {
        (Sign::Zero, Sign::Positive) => 90,
        (Sign::Negative, Sign::Zero) => 180,
        (Sign::Zero, Sign::Negative) => 270,
        _ => return Ok(None),
    };
    Ok(Some((unit(&a), Rational::from_i64(angle))))
}

fn region(
    definition: ChamferDefinition,
    interval: Domain,
    finite: Option<(Vector3, Rational)>,
    step: &StepId,
    dec: &mut Decider,
) -> Result<(Feature, Vec<crate::Face>), Refusal> {
    let (origin, x, y, z, e, a, b, sweep) = match (&definition.edge, &definition.contacts) {
        (Curve::Line(edge), [Curve::Line(a), Curve::Line(b)]) => {
            let z = unit(&edge.dir);
            let x = unit(&a.origin.minus(&edge.origin));
            let y = z.cross(&x);
            let Domain::Line([start, _]) = &interval else {
                return Err(bad("a line bevel needs a finite input interval"));
            };
            let shift = edge.dir.scaled(&Real::from(start.clone()));
            let e = edge.origin.plus(&shift);
            (
                e.clone(),
                x,
                y,
                z,
                e,
                a.origin.plus(&shift),
                b.origin.plus(&shift),
                None,
            )
        }
        (Curve::Circle(edge), [Curve::Circle(a), Curve::Circle(b)]) => {
            let x = finite
                .as_ref()
                .map_or_else(|| unit(&edge.ref_dir), |(x, _)| x.clone());
            let axis = unit(&edge.normal);
            let z = x.cross(&axis);
            let e = edge.center.plus(&x.scaled(&square_root(&edge.radius_sq)));
            let a = a.center.plus(&x.scaled(&square_root(&a.radius_sq)));
            let b = b.center.plus(&x.scaled(&square_root(&b.radius_sq)));
            (
                edge.center.clone(),
                x,
                axis.clone(),
                z,
                e,
                a,
                b,
                Some((edge.center.clone(), axis)),
            )
        }
        _ => return Err(bad("unknown elementary bevel section")),
    };
    let local = |p: &Point3| -> Result<catallax_recipe::Point2, Refusal> {
        let d = p.minus(&origin);
        Ok(catallax_recipe::Point2(
            [q(&d.dot(&x))?, q(&d.dot(&y))?].map(number),
        ))
    };
    let id = |s| ElementId::new(s).expect("fixed section id");
    let loops = vec![catallax_recipe::Loop(vec![
        catallax_recipe::Segment::Line {
            id: id("a"),
            to: local(&a)?,
        },
        catallax_recipe::Segment::Line {
            id: id("bevel"),
            to: local(&b)?,
        },
        catallax_recipe::Segment::Line {
            id: id("b"),
            to: local(&e)?,
        },
    ])];
    let pose = pose(&origin, &x, &y, &z, dec)?;
    let profile = profile::Profile::new(&pose, &loops, dec)?;
    let (mut feature, mut faces) = if let Some((origin, axis)) = sweep {
        Feature::revolve(
            &profile,
            step,
            &catallax_recipe::Axis {
                point: point(&origin)?,
                direction: vector(&axis)?,
            },
            &finite.map_or_else(|| Rational::from_i64(360), |(_, angle)| angle),
            dec,
        )?
    } else {
        let Domain::Line([start, end]) = &interval else {
            unreachable!("line interval")
        };
        let Curve::Line(edge) = &definition.edge else {
            unreachable!("line interval");
        };
        let distance = (end - start) * q(&square_root(&edge.dir.norm_squared()))?;
        Feature::extrude(
            &profile,
            step,
            Some(&distance),
            Some(&vector(&z)?),
            false,
            dec,
        )?
    }
    .ok_or_else(|| bad("empty chamfer section"))?;
    faces.retain(
        |f| matches!(&f.origin.role, Role::Side { segment, .. } if segment.as_str() == "bevel"),
    );
    for face in &mut faces {
        face.origin = definition.origin.clone();
    }
    feature.chamfer = Some(Box::new(definition));
    Ok((feature, faces))
}

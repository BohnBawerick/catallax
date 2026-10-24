//! Rolling-ball edge regions compiled from recipe supports, without B-rep topology.

use catallax_features::fillet::FilletDefinition;
use catallax_intersect::{Component, CurveSurface, SurfacePair};
use catallax_numbers::{Decider, Rational, Real, Refusal, Sign};
use catallax_recipe::{
    BooleanKind, EdgeName, EdgeSelection, ElementId, Origin, Role, StepId, Turn,
};
use catallax_shapes::{Curve, GeomStore, Point3, Surface, Vector3};

use crate::{Instruction, Location, Plan, features::Feature, geometry, profile};

fn bad(detail: &str) -> Refusal {
    Refusal::new("field.fillet-contact", detail)
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
        let i = (0..3)
            .find(|i| geometry::sign(&(&one + &m[*i][*i]), dec).is_ok_and(|s| s == Sign::Positive))
            .ok_or_else(|| bad("section frame is not a proper rotation"))?;
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
    fn fillet_ray(
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
    fn material(
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
            match self.fillet_ray(root, &xyz, &direction, dec)? {
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
    fn fillet_domain(
        &self,
        definition: &FilletDefinition,
        faces: [usize; 2],
        interval: &Option<[Rational; 2]>,
        source_faces: &[usize],
        geom: &mut GeomStore,
        dec: &mut Decider,
    ) -> Result<(), Refusal> {
        let mut end_planes = [false; 2];
        for &i in source_faces {
            let face = &self.faces[i];
            if faces
                .iter()
                .any(|j| face.support.same_definition(&self.faces[*j].support))
            {
                continue;
            }
            if let (Curve::Line(edge), Surface::Plane(p), Some(interval)) =
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
            if definition.avoids_support(&face.surface, dec)? {
                continue;
            }
            let a = geom.add_surface(definition.surface.clone());
            let b = geom.add_surface(self.faces[i].surface.clone());
            if !matches!(
                catallax_intersect::intersect_surfaces(geom, a, b, dec)?.value,
                SurfacePair::Disjoint
            ) {
                return Err(Refusal::new(
                    "fillet.unsupported-overflow",
                    "the field contact region reaches another support and needs trimmed continuation",
                ));
            }
        }
        if interval.is_some() && end_planes != [true, true] {
            return Err(bad(
                "finite line strips need proven perpendicular end planes",
            ));
        }
        Ok(())
    }

    pub(crate) fn fillet(
        &mut self,
        root: usize,
        source_faces: &[usize],
        selection: &EdgeSelection,
        radius: &Real,
        step: &StepId,
        dec: &mut Decider,
    ) -> Result<usize, Refusal> {
        if geometry::sign(radius, dec)? != Sign::Positive {
            return Err(Refusal::new(
                "fillet.non-positive-radius",
                "fillet radius must be positive",
            ));
        }
        let mut pairs = Vec::new();
        match selection {
            EdgeSelection::All => {
                for (position, &a) in source_faces.iter().enumerate() {
                    for &b in &source_faces[position + 1..] {
                        pairs.push([a, b]);
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
                            if (reference.between[0].matches(&self.faces[a].origin)
                                && reference.between[1].matches(&self.faces[b].origin))
                                || (reference.between[1].matches(&self.faces[a].origin)
                                    && reference.between[0].matches(&self.faces[b].origin))
                            {
                                if !pairs.contains(&[a, b]) {
                                    pairs.push([a, b]);
                                }
                                found = true;
                            }
                        }
                    }
                    if !found {
                        return Err(Refusal::new(
                            "fillet.lost-edge",
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
        for faces in pairs {
            let supports = faces.map(|i| self.faces[i].surface.clone());
            if !matches!(
                (&supports[0], &supports[1]),
                (Surface::Plane(_), Surface::Plane(_) | Surface::Cylinder(_))
                    | (Surface::Cylinder(_), Surface::Plane(_))
            ) {
                if matches!(selection, EdgeSelection::All) {
                    continue;
                }
                return Err(bad(
                    "elementary fillets need planes or a perpendicular plane and cylinder",
                ));
            }
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
                    Curve::Circle(c) => vec![(
                        c.center
                            .plus(&unit(&c.ref_dir).scaled(&square_root(&c.radius_sq))),
                        None,
                    )],
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
                                    Some([pair[0].clone(), pair[1].clone()]),
                                )
                            })
                            .collect()
                    }
                    _ => continue,
                };
                let first_definition = definitions.len();
                for (sample, interval) in samples {
                    let Some((convex, senses)) = self.material(root, &sample, faces, dec)? else {
                        continue;
                    };
                    if definitions.len() > first_definition {
                        let (previous, previous_interval): &mut (
                            FilletDefinition,
                            Option<[Rational; 2]>,
                        ) = definitions.last_mut().unwrap();
                        if let (Some(previous_interval), Some(current)) =
                            (previous_interval, &interval)
                            && previous_interval[1] == current[0]
                            && previous.convex == convex
                            && previous.same_sense == senses
                        {
                            previous_interval[1] = current[1].clone();
                            continue;
                        }
                    }
                    let edge_name = EdgeName::new(
                        self.faces[faces[0]].origin.clone(),
                        self.faces[faces[1]].origin.clone(),
                    )
                    .map_err(|e| bad(&e.to_string()))?;
                    let definition = FilletDefinition::new(
                        Origin {
                            step: step.clone(),
                            role: Role::Blend { edge: edge_name },
                        },
                        supports.clone(),
                        senses,
                        edge.clone(),
                        radius.clone(),
                        convex,
                        dec,
                    )?;
                    definitions.push((definition, interval));
                }
            }
        }
        if definitions.is_empty()
            && !matches!(selection, EdgeSelection::Refs(refs) if refs.is_empty())
        {
            return Err(Refusal::new(
                "fillet.lost-edge",
                "no selected support intersection is a sharp input boundary",
            ));
        }
        if let EdgeSelection::Refs(refs) = selection {
            for reference in refs {
                if !definitions.iter().any(|(definition, _)| {
                    let Role::Blend { edge } = &definition.origin.role else {
                        return false;
                    };
                    let [a, b] = edge.between();
                    (reference.between[0].matches(a) && reference.between[1].matches(b))
                        || (reference.between[1].matches(a) && reference.between[0].matches(b))
                }) {
                    return Err(Refusal::new(
                        "fillet.lost-edge",
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
            if let (Curve::Line(edge), Some(interval)) = (&definition.edge, interval) {
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
            self.fillet_domain(
                &definition,
                faces,
                &interval,
                source_faces,
                &mut GeomStore::new(),
                dec,
            )?;
            let convex = definition.convex;
            let (feature, faces) = region(definition, interval, step, dec)?;
            let node = self.leaf(geometry::Leaf::Feature(feature), faces, dec)?;
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

fn region(
    definition: FilletDefinition,
    interval: Option<[Rational; 2]>,
    step: &StepId,
    dec: &mut Decider,
) -> Result<(Feature, Vec<crate::Face>), Refusal> {
    let (origin, x, y, z, e, a, b, centre, sweep) =
        match (&definition.edge, &definition.surface, &definition.contacts) {
            (Curve::Line(edge), Surface::Cylinder(c), [Curve::Line(a), Curve::Line(b)]) => {
                let z = unit(&edge.dir);
                let x = unit(&a.origin.minus(&edge.origin));
                let y = z.cross(&x);
                let [start, _] = interval
                    .as_ref()
                    .ok_or_else(|| bad("a line strip needs a finite input interval"))?;
                let shift = edge.dir.scaled(&Real::from(start.clone()));
                (
                    edge.origin.plus(&shift),
                    x,
                    y,
                    z,
                    edge.origin.plus(&shift),
                    a.origin.plus(&shift),
                    b.origin.plus(&shift),
                    c.origin.plus(&shift),
                    None,
                )
            }
            (Curve::Circle(edge), Surface::Torus(t), [Curve::Circle(a), Curve::Circle(b)]) => {
                let x = unit(&edge.ref_dir);
                let axis = unit(&t.axis);
                let z = x.cross(&axis);
                let e = edge.center.plus(&x.scaled(&square_root(&edge.radius_sq)));
                let a = a.center.plus(&x.scaled(&square_root(&a.radius_sq)));
                let b = b.center.plus(&x.scaled(&square_root(&b.radius_sq)));
                let centre = t.center.plus(&x.scaled(&t.major));
                (
                    edge.center.clone(),
                    x,
                    axis.clone(),
                    z,
                    e,
                    a,
                    b,
                    centre,
                    Some((edge.center.clone(), axis)),
                )
            }
            _ => return Err(bad("unknown elementary contact section")),
        };
    let local = |p: &Point3| -> Result<[Rational; 2], Refusal> {
        let d = p.minus(&origin);
        Ok([q(&d.dot(&x))?, q(&d.dot(&y))?])
    };
    let e = local(&e)?;
    let a = local(&a)?;
    let b = local(&b)?;
    let centre = local(&centre)?;
    let va = [&a[0] - &centre[0], &a[1] - &centre[1]];
    let vb = [&b[0] - &centre[0], &b[1] - &centre[1]];
    let turn = if (&va[0] * &vb[1] - &va[1] * &vb[0]).sign() == Sign::Positive {
        Turn::Ccw
    } else {
        Turn::Cw
    };
    let id = |s| ElementId::new(s).expect("fixed section id");
    let p2 = |p: [Rational; 2]| catallax_recipe::Point2(p.map(number));
    let loops = vec![catallax_recipe::Loop(vec![
        catallax_recipe::Segment::Line {
            id: id("a"),
            to: p2(a),
        },
        catallax_recipe::Segment::Arc {
            id: id("blend"),
            centre: p2(centre),
            to: p2(b),
            turn,
        },
        catallax_recipe::Segment::Line {
            id: id("b"),
            to: p2(e),
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
            &Rational::from_i64(360),
            dec,
        )?
    } else {
        let [start, end] = interval.as_ref().expect("line strip interval");
        let Curve::Line(edge) = &definition.edge else {
            unreachable!("line strip");
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
    .ok_or_else(|| bad("empty rolling-ball section"))?;
    faces.retain(|f| matches!(&f.origin.role,Role::Side {segment,..} if segment.as_str()=="blend"));
    for face in &mut faces {
        face.origin = definition.origin.clone();
    }
    feature.fillet = Some(Box::new(definition));
    Ok((feature, faces))
}

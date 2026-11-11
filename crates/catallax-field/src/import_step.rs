//! Shared STEP admission followed by retained field definitions.

use crate::{Plan, draft, geometry};
use catallax_features::ImportFiles;
use catallax_numbers::{Decider, Rational, Real, Refusal, Sign};
use catallax_recipe::{BooleanKind, Operation, Origin, Recipe, Step};
use catallax_shapes::{Curve, Plane, Point3, Surface};
use catallax_topology::{BodyId, FaceId, Loop, Model, Sense};
use std::collections::{BTreeMap, BTreeSet};

fn unsupported(detail: impl Into<String>) -> Refusal {
    Refusal::new("field.import-domain", detail)
}
fn source_error(error: catallax_import::ImportRefusal) -> Refusal {
    use catallax_import::ImportRefusal;
    match error {
        ImportRefusal::Unsupported { entity, type_name } => Refusal::new(
            "field.import-unsupported",
            format!(
                "{entity}: {type_name}; the shared analytic STEP importer does not build this definition yet"
            ),
        ),
        error => Refusal::new("field.import-source", format!("{error:?}")),
    }
}
fn outward(model: &Model, id: FaceId) -> Result<Plane, Refusal> {
    let face = model.face(id);
    let Surface::Plane(plane) = model.geom().surface(face.surface) else {
        return Err(unsupported("expected a planar import face"));
    };
    let mut plane = plane.clone();
    if face.sense == Sense::Reversed {
        plane.normal = plane.normal.scaled(&Real::from_i64(-1));
    }
    Ok(plane)
}
fn rational(value: &Real, what: &str) -> Result<Rational, Refusal> {
    value
        .as_rational()
        .ok_or_else(|| unsupported(format!("{what} is not rational after placement")))
}

fn rational_sqrt(q: &Rational) -> Option<Rational> {
    let root = |n: &catallax_numbers::Integer| {
        let mut x = catallax_numbers::Integer::one().mul_pow2(n.bit_length().div_ceil(2));
        let two = catallax_numbers::Integer::from_i64(2);
        loop {
            let (quotient, _) = n.div_rem_floor(&x)?;
            let (y, _) = (&x + &quotient).div_rem_floor(&two)?;
            if y >= x {
                return (&x * &x == *n).then_some(x);
            }
            x = y;
        }
    };
    if q.sign() != Sign::Positive {
        return None;
    }
    Rational::new(root(q.numerator())?, root(q.denominator())?)
}

/// A complete imported sphere or torus is the field primitive of the same size
/// and placement. Spherical and toroidal faces with trims are refused.
fn closed_primitive(model: &Model, body: BodyId) -> Result<Option<Operation>, Refusal> {
    use catallax_recipe::{Number, Pose, Vector3};
    let faces = model.faces_of_body(body);
    let closed = |id: &FaceId| {
        matches!(
            model.geom().surface(model.face(*id).surface),
            Surface::Sphere(_) | Surface::Torus(_)
        )
    };
    if !faces.iter().any(closed) {
        return Ok(None);
    }
    let [id] = faces.as_slice() else {
        return Err(unsupported(format!(
            "{} faces: spherical and toroidal faces with trims need a finite curved-domain field; only a complete sphere or torus builds",
            faces.len()
        )));
    };
    let face = model.face(*id);
    if !face.loops.is_empty() {
        return Err(unsupported(
            "a trimmed sphere or torus alone does not enclose material",
        ));
    }
    if face.sense != Sense::Same {
        return Err(unsupported(
            "an inward complete sphere or torus is not a filled solid",
        ));
    }
    let n = |q: Rational| Number::new(q);
    let translation = |p: &Point3| -> Result<Vector3, Refusal> {
        Ok(Vector3([
            n(rational(&p.x, "centre")?),
            n(rational(&p.y, "centre")?),
            n(rational(&p.z, "centre")?),
        ]))
    };
    Ok(Some(match model.geom().surface(face.surface) {
        Surface::Sphere(sphere) => Operation::Sphere {
            radius: n(rational(&sphere.radius, "sphere radius")?),
            pose: Pose {
                rotation: Pose::identity().rotation,
                translation: translation(&sphere.center)?,
            },
        },
        Surface::Torus(torus) => {
            let axis = [&torus.axis.x, &torus.axis.y, &torus.axis.z]
                .map(|c| rational(c, "torus axis"))
                .into_iter()
                .collect::<Result<Vec<_>, _>>()?;
            let length = rational_sqrt(&axis.iter().fold(Rational::zero(), |sum, c| sum + c * c))
                .ok_or_else(|| {
                unsupported(
                    "torus axis length is irrational; no exact rational rotation carries z onto it",
                )
            })?;
            let [x, y, z] = [0, 1, 2].map(|i| &axis[i] / &length);
            // The rotation of z onto the unit axis: quaternion (1 + z.a, z x a).
            let rotation = if z == -Rational::one() {
                [
                    Rational::zero(),
                    Rational::one(),
                    Rational::zero(),
                    Rational::zero(),
                ]
            } else {
                [&Rational::one() + &z, -y, x, Rational::zero()]
            };
            Operation::Torus {
                major: n(rational(&torus.major, "torus major radius")?),
                minor: n(rational(&torus.minor, "torus minor radius")?),
                pose: Pose {
                    rotation: rotation.map(n),
                    translation: translation(&torus.center)?,
                },
            }
        }
        _ => unreachable!("checked spherical or toroidal"),
    }))
}

fn planar(
    model: &Model,
    body: BodyId,
) -> Result<catallax_features::draft::arrangement::Skin, Refusal> {
    use catallax_features::draft::arrangement::{Patch, Skin};
    let mut vertices = Vec::new();
    let mut ids = BTreeMap::new();
    let mut patches = Vec::new();
    for id in model.faces_of_body(body) {
        let face = model.face(id);
        let mut loops = Vec::new();
        for lp in &face.loops {
            let Loop::Coedges(cs) = model.loop_(*lp) else {
                return Err(unsupported("planar import needs finite polygon loops"));
            };
            let mut cycle = Vec::new();
            for c in cs {
                let c = model.coedge(*c);
                let edge = model.edge(c.edge);
                if !matches!(model.geom().curve(edge.curve), Curve::Line(_)) {
                    return Err(unsupported(
                        "planar import with curved boundary needs a finite curved-domain field",
                    ));
                }
                let vertex = if c.sense == Sense::Same {
                    edge.start
                } else {
                    edge.end
                }
                .ok_or_else(|| unsupported("polygon edge has no finite endpoint"))?;
                cycle.push(*ids.entry(vertex).or_insert_with(|| {
                    let index = vertices.len();
                    vertices.push(model.geom().position(model.vertex(vertex).point));
                    index
                }));
            }
            if face.sense == Sense::Reversed {
                cycle.reverse();
            }
            loops.push(cycle);
        }
        patches.push(Patch {
            plane: outward(model, id)?,
            loops,
            origins: model.face_origins(id).to_vec(),
        });
    }
    Ok(Skin { vertices, patches })
}

pub(crate) fn on_named_patch(
    skin: &catallax_features::draft::arrangement::Skin,
    origin: &Origin,
    point: &Point3,
    dec: &mut Decider,
) -> Result<bool, Refusal> {
    for patch in skin.patches.iter().filter(|p| p.origins.contains(origin)) {
        if geometry::sign(
            &point.minus(&patch.plane.origin).dot(&patch.plane.normal),
            dec,
        )? != Sign::Zero
        {
            continue;
        }
        for cycle in &patch.loops {
            for i in 0..cycle.len() {
                let a = point.minus(&skin.vertices[cycle[i]]);
                let b = point.minus(&skin.vertices[cycle[(i + 1) % cycle.len()]]);
                if geometry::sign(&a.cross(&b).norm_squared(), dec)? == Sign::Zero
                    && geometry::sign(&a.dot(&b), dec)? != Sign::Positive
                {
                    return Ok(true);
                }
            }
        }
        if skin.contains_patch(patch, point, dec)? {
            return Ok(true);
        }
    }
    Ok(false)
}

fn trim_planes(
    skin: &catallax_features::draft::arrangement::Skin,
    dec: &mut Decider,
) -> Result<Vec<Plane>, Refusal> {
    let mut planes = Vec::new();
    for patch in &skin.patches {
        for cycle in &patch.loops {
            for i in 0..cycle.len() {
                let a = &skin.vertices[cycle[i]];
                let edge = skin.vertices[cycle[(i + 1) % cycle.len()]].minus(a);
                let normal = edge.cross(&patch.plane.normal);
                if geometry::sign(&normal.norm_squared(), dec)? == Sign::Zero {
                    return Err(unsupported("planar source trim has a degenerate edge"));
                }
                planes.push(Plane {
                    origin: a.clone(),
                    normal,
                    ref_dir: edge,
                });
            }
        }
    }
    Ok(planes)
}
fn ring(
    model: &Model,
    id: FaceId,
) -> Result<(catallax_topology::EdgeId, catallax_shapes::Circle), Refusal> {
    let [lp] = model.face(id).loops.as_slice() else {
        return Err(unsupported(
            "round import needs one complete disk at each end",
        ));
    };
    let Loop::Coedges(cs) = model.loop_(*lp) else {
        return Err(unsupported("round cap is not a disk"));
    };
    let [c] = cs.as_slice() else {
        return Err(unsupported("round cap has additional trims"));
    };
    let id = model.coedge(*c).edge;
    let edge = model.edge(id);
    let Curve::Circle(circle) = model.geom().curve(edge.curve) else {
        return Err(unsupported("round cap is not bounded by a circle"));
    };
    if edge.start.is_some() || edge.end.is_some() {
        return Err(unsupported("round import needs a complete circular rim"));
    }
    Ok((id, circle.clone()))
}
fn round(
    model: &Model,
    body: BodyId,
    dec: &mut Decider,
) -> Result<catallax_features::draft::round::RoundCell, Refusal> {
    use catallax_features::draft::round::RoundCell;
    let fs = model.faces_of_body(body);
    if !(2..=3).contains(&fs.len()) {
        return Err(unsupported(format!(
            "{} faces: curved import currently needs one wall and two axial disks, or a cone wall, its apex and one disk",
            fs.len()
        )));
    }
    let mut wall = None;
    let mut caps = Vec::new();
    for id in fs {
        match model.geom().surface(model.face(id).surface) {
            Surface::Plane(_) => caps.push(id),
            Surface::Cylinder(_) | Surface::Cone(_) if wall.is_none() => wall = Some(id),
            _ => {
                return Err(unsupported(
                    "curved import needs one cylindrical or conical wall",
                ));
            }
        }
    }
    let Some(wall) = wall else {
        return Err(unsupported("curved import needs one wall"));
    };
    let pointed = caps.len() == 1;
    if pointed
        && !matches!(
            model.geom().surface(model.face(wall).surface),
            Surface::Cone(_)
        )
    {
        return Err(unsupported("curved import needs two caps"));
    }
    if model.face(wall).sense != Sense::Same {
        return Err(unsupported("inward round wall is not a filled axial cell"));
    }
    let surface = model.geom().surface(model.face(wall).surface).clone();
    let (origin, axis, reference) = match &surface {
        Surface::Cylinder(c) => (c.origin.clone(), c.axis.clone(), c.ref_dir.clone()),
        Surface::Cone(c) => (c.apex.clone(), c.axis.clone(), c.ref_dir.clone()),
        _ => unreachable!(),
    };
    let axis = axis.scaled(&(&Real::one() / &axis.norm_squared().sqrt()));
    let mut ends = Vec::new();
    for id in caps.iter().copied() {
        let plane = outward(model, id)?;
        let (edge, circle) = ring(model, id)?;
        if geometry::sign(&plane.normal.cross(&axis).norm_squared(), dec)? != Sign::Zero
            || geometry::sign(
                &circle.center.minus(&origin).cross(&axis).norm_squared(),
                dec,
            )? != Sign::Zero
        {
            return Err(unsupported(
                "round import needs coaxial rims and perpendicular caps",
            ));
        }
        let h = circle.center.minus(&origin).dot(&axis);
        ends.push((
            h,
            circle.radius_sq.sqrt(),
            plane,
            model.face_origins(id).to_vec(),
            edge,
        ));
    }
    let mut apex = None;
    if pointed {
        // A pointed cone closes at its apex, the cone's origin: height and
        // radius zero, with no cap.
        let wall_loops = &model.face(wall).loops;
        let vertex = wall_loops.iter().find_map(|lp| match model.loop_(*lp) {
            Loop::Vertex(v) => Some(*v),
            _ => None,
        });
        let Some(vertex) = vertex else {
            return Err(unsupported("pointed round wall has no apex vertex"));
        };
        let at = model.geom().position(model.vertex(vertex).point);
        if geometry::sign(&at.minus(&origin).norm_squared(), dec)? != Sign::Zero {
            return Err(unsupported("pointed round wall does not end at its apex"));
        }
        if geometry::sign(&ends[0].0, dec)? != Sign::Positive {
            return Err(unsupported("pointed cone cap is not on its nappe"));
        }
        apex = Some(vertex);
        let plane = ends[0].2.clone();
        ends.insert(
            0,
            (Real::zero(), Real::zero(), plane, Vec::new(), ends[0].4),
        );
    }
    match geometry::sign(&(&ends[1].0 - &ends[0].0), dec)? {
        Sign::Negative => ends.swap(0, 1),
        Sign::Zero => return Err(unsupported("round import has no positive axial span")),
        Sign::Positive => (),
    }
    for (index, end) in ends.iter().enumerate() {
        if index == 0 && apex.is_some() {
            continue;
        }
        if geometry::sign(&end.2.normal.dot(&axis), dec)?
            != if index == 0 {
                Sign::Negative
            } else {
                Sign::Positive
            }
        {
            return Err(unsupported(
                "round import cap orientation does not enclose material",
            ));
        }
        if matches!(surface, Surface::Cone(_)) && geometry::sign(&end.0, dec)? != Sign::Positive {
            return Err(unsupported("conical import crosses its selected nappe"));
        }
    }
    let wall_loops = &model.face(wall).loops;
    if wall_loops.len() != 2 {
        return Err(unsupported("round wall has additional or missing trims"));
    }
    let mut wall_edges = Vec::new();
    for lp in wall_loops {
        let cs = match model.loop_(*lp) {
            Loop::Coedges(cs) => cs,
            Loop::Vertex(v) if Some(*v) == apex => continue,
            Loop::Vertex(_) => return Err(unsupported("round wall trim is not a ring")),
        };
        let [c] = cs.as_slice() else {
            return Err(unsupported("round wall trim is not one full circle"));
        };
        wall_edges.push(model.coedge(*c).edge);
    }
    if !ends
        .iter()
        .skip(usize::from(apex.is_some()))
        .all(|end| wall_edges.contains(&end.4))
    {
        return Err(unsupported(
            "round wall and caps do not share their complete rims",
        ));
    }
    Ok(RoundCell {
        origin,
        axis,
        surface,
        reference,
        origins: model.face_origins(wall).to_vec(),
        heights: [ends[0].0.clone(), ends[1].0.clone()],
        radii: [ends[0].1.clone(), ends[1].1.clone()],
        caps: [
            apex.is_none()
                .then(|| (ends[0].2.clone(), ends[0].3.clone())),
            Some((ends[1].2.clone(), ends[1].3.clone())),
        ],
    })
}

impl Plan {
    pub(crate) fn import_step(
        &mut self,
        step: &Step,
        files: &ImportFiles<'_>,
        dec: &mut Decider,
    ) -> Result<usize, Refusal> {
        let Operation::ImportStep { file, .. } = &step.op else {
            unreachable!()
        };
        let bytes = files.get(&file.sha256).ok_or_else(|| {
            Refusal::new(
                "field.import-missing-source",
                format!(
                    "{}: missing STEP bytes with sha256 {}",
                    step.id, file.sha256
                ),
            )
        })?;
        let actual = catallax_import::sha256_hex(bytes);
        if actual != file.sha256 {
            return Err(Refusal::new(
                "import.hash-mismatch",
                format!("expected {}, actual {actual}", file.sha256),
            ));
        }
        let source = catallax_import::read_step(bytes).map_err(source_error)?;
        for face in source.solids.iter().flat_map(|solid| &solid.faces) {
            if let catallax_import::Surface::BSpline(_) = &face.surface {
                return Err(Refusal::new(
                    "field.import-spline",
                    format!("{}: spline solid fields wait for M10", face.name),
                ));
            }
        }
        // Admission, saved repairs and placement are shared; field queries retain
        // only the analytic definitions below, never the imported B-rep model.
        let recipe = Recipe::new("field STEP source", vec![step.clone()], step.id.clone());
        let built =
            catallax_features::evaluate_with_files(&recipe, files, dec).map_err(|e| match e {
                catallax_features::BuildError::Geometry { refusal, .. } => refusal,
                error => Refusal::new("field.import-source", error.to_string()),
            })?;
        let model = &built.model;
        let body = built.result;
        if let Some(op) = closed_primitive(model, body)? {
            let (leaf, faces) = geometry::compile(&op, &step.id, dec)?;
            let [face] = faces.as_slice() else {
                unreachable!("a sphere or torus primitive has one face")
            };
            let id = model.faces_of_body(body)[0];
            let faces = model
                .face_origins(id)
                .iter()
                .map(|origin| {
                    let mut named = face.clone();
                    named.origin = origin.clone();
                    named
                })
                .collect();
            return self.leaf(leaf, faces, dec);
        }
        if model.faces_of_body(body).iter().all(|id| {
            matches!(
                model.geom().surface(model.face(*id).surface),
                Surface::Plane(_)
            )
        }) {
            let skin = planar(model, body)?;
            let components = skin
                .components_with_trim_planes(&trim_planes(&skin, dec)?, dec)
                .map_err(|error| {
                    if error.code.starts_with("draft.") {
                        unsupported(format!(
                            "{} source faces: {}: {}",
                            skin.patches.len(),
                            error.code,
                            error.detail
                        ))
                    } else {
                        error
                    }
                })?;
            let expected: BTreeSet<_> = model
                .faces_of_body(body)
                .iter()
                .flat_map(|id| model.face_origins(*id))
                .collect();
            let retained: BTreeSet<_> = components
                .iter()
                .flat_map(|c| &c.faces)
                .flat_map(|f| &f.origins)
                .collect();
            if components.is_empty() || !expected.is_subset(&retained) {
                return Err(unsupported(format!(
                    "retained planar domains preserve {} of {} admitted source face names; no empty or partial import is published",
                    expected.intersection(&retained).count(),
                    expected.len(),
                )));
            }
            let mut result = 0;
            for component in components {
                let mut faces = Vec::new();
                for face in &component.faces {
                    for origin in &face.origins {
                        faces.push(geometry::face(
                            Surface::Plane(component.plane(face)),
                            origin.clone(),
                            dec,
                        )?);
                    }
                }
                let keep_trims = component.faces.iter().any(|f| f.origins.len() > 1);
                let mut region =
                    draft::Region::new(draft::Cell::Polyhedral(Box::new(component)), dec)?;
                if keep_trims {
                    region.imported_patches = Some(Box::new(skin.clone()));
                }
                let node = self.leaf(geometry::Leaf::Draft(region), faces, dec)?;
                result = self.boolean(BooleanKind::Union, result, node, dec)?;
            }
            Ok(result)
        } else {
            let cell = round(model, body, dec)?;
            let mut faces = Vec::new();
            for origin in &cell.origins {
                faces.push(geometry::face(cell.surface.clone(), origin.clone(), dec)?);
            }
            for (plane, origins) in cell.caps.iter().flatten() {
                for origin in origins {
                    faces.push(geometry::face(
                        Surface::Plane(plane.clone()),
                        origin.clone(),
                        dec,
                    )?);
                }
            }
            let region = draft::Region::new(draft::Cell::Round(Box::new(cell)), dec)?;
            self.leaf(geometry::Leaf::Draft(region), faces, dec)
        }
    }
}

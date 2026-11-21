//! Rolling-ball edge regions compiled from recipe supports, without B-rep topology.

use catallax_features::{fillet::FilletDefinition, fillet_corner::CornerDefinition};
use catallax_intersect::{Component, CurveSurface, SurfacePair};
use catallax_numbers::{Decider, Rational, Real, Refusal, Sign};
use catallax_recipe::{
    BooleanKind, EdgeName, EdgeSelection, ElementId, Origin, Role, StepId, Turn, VertexName,
};
use catallax_shapes::{Curve, CurveId, GeomStore, Line, Plane, Point3, Surface, Vector3};

use crate::{Instruction, Location, Plan, features::Feature, geometry, profile};

#[path = "fillet_network.rs"]
mod network;
#[path = "fillet_rim.rs"]
mod rim;
#[path = "fillet_rolling.rs"]
mod rolling;
#[path = "fillet_smooth.rs"]
mod smooth;

/// Edge points with the station interval of each, for line runs.
type Samples = Vec<(Point3, Option<[Rational; 2]>)>;

pub(crate) struct Settings<'a> {
    pub radius: &'a Real,
    pub closure: catallax_recipe::FilletCorner,
}

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
    /// [`Plan::material`] at a point that may be algebraic, such as a point of
    /// a circle of radius `sqrt 21`, which has no rational point. Rational
    /// points take [`Plan::material`] unchanged; algebraic ones use the same
    /// exact membership and ray probes through the field's algebraic queries.
    fn material_point(
        &self,
        root: usize,
        p: &Point3,
        faces: [usize; 2],
        dec: &mut Decider,
    ) -> Result<Option<(bool, [bool; 2])>, Refusal> {
        if coordinates(p).is_ok() {
            return self.material(root, p, faces, dec);
        }
        for i in faces {
            let face = &self.faces[i];
            let leaf = &self.leaves[face.leaf];
            if crate::distance::location(leaf, p, dec)? != Location::Boundary
                || leaf.on_face_algebraic(&face.origin, self.definition_role(i), p, dec)?
                    == Some(false)
            {
                return Ok(None);
            }
        }
        let [a, b] = faces.map(|i| self.gradient_at(i, p));
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
            let mut locations = Vec::with_capacity(root + 1);
            for instruction in &self.instructions[..=root] {
                let location = match *instruction {
                    Instruction::Empty => Location::Outside,
                    Instruction::Leaf(i) => {
                        crate::distance::ray(&self.leaves[i], p, &direction, dec)?
                    }
                    Instruction::Boolean { kind, a, b } => {
                        crate::combine(kind, locations[a], locations[b])
                    }
                };
                locations.push(location);
            }
            match locations[root] {
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
    /// A face's gradient at a point that may be algebraic.
    fn gradient_at(&self, face: usize, p: &Point3) -> Vector3 {
        let values = [p.x.clone(), p.y.clone(), p.z.clone()];
        let gradient = &self.faces[face].gradient;
        Vector3::new(
            gradient[0].eval_real(&values),
            gradient[1].eval_real(&values),
            gradient[2].eval_real(&values),
        )
    }
    fn fillet_domain(
        &self,
        definition: &FilletDefinition,
        faces: [usize; 2],
        interval: &mut Option<[Rational; 2]>,
        source_faces: &[usize],
        geom: &mut GeomStore,
        dec: &mut Decider,
    ) -> Result<(), Refusal> {
        let mut end_planes = [false; 2];
        let mut span = interval.clone();
        let mut others = Vec::with_capacity(source_faces.len());
        for &i in source_faces {
            let face = &self.faces[i];
            if faces
                .iter()
                .any(|j| face.support.same_definition(&self.faces[*j].support))
            {
                continue;
            }
            if let (Curve::Line(edge), Surface::Plane(p), Some(ends)) =
                (&definition.edge, &face.surface, interval.as_ref())
                && geometry::sign(&edge.dir.dot(&p.normal), dec)? != Sign::Zero
            {
                let station =
                    |o: &Point3| q(&(p.origin.minus(o).dot(&p.normal) / edge.dir.dot(&p.normal)));
                let at = station(&edge.origin)?;
                if let Some(end) = ends.iter().position(|value| *value == at) {
                    end_planes[end] = true;
                    // An oblique end plane leans past the vertex at a contact
                    // rail; the region's section spans the edge and both rails.
                    for contact in &definition.contacts {
                        let Curve::Line(rail) = contact else {
                            return Err(bad("a line strip needs straight contact rails"));
                        };
                        let t = station(&rail.origin)?;
                        let span = span.as_mut().expect("line strip interval");
                        span[end] = if end == 0 {
                            t.min(span[0].clone())
                        } else {
                            t.max(span[1].clone())
                        };
                    }
                    continue;
                }
            }
            others.push(i);
        }
        // The region spans its stations, so a support clear of that prism is
        // clear of the region.
        let hull = match (&definition.edge, &span) {
            (Curve::Line(edge), Some(span)) if end_planes == [true, true] => {
                let [first, second] = span
                    .each_ref()
                    .map(|t| edge.origin.plus(&edge.dir.scaled(&Real::from(t.clone()))));
                definition.cut_hull([(&first, &edge.dir), (&second, &edge.dir)], dec)?
            }
            _ => None,
        };
        let ring_hull = if hull.is_none() {
            definition.ring_hull(dec)?
        } else {
            None
        };
        let hull = hull
            .as_ref()
            .map(|points| points.as_slice())
            .or(ring_hull.as_deref());
        for i in others {
            let face = &self.faces[i];
            if definition.avoids_support(&face.surface, dec)? {
                continue;
            }
            if let Some(hull) = hull
                && catallax_features::fillet::hull_clear(hull, &face.surface, dec)?
            {
                continue;
            }
            let a = geom.add_surface(definition.surface.clone());
            let b = geom.add_surface(self.faces[i].surface.clone());
            if !matches!(
                catallax_intersect::intersect_surfaces(geom, a, b, dec)?.value,
                SurfacePair::Disjoint
            ) {
                if finite_round_clear(definition, &self.leaves[face.leaf], dec)? {
                    continue;
                }
                return Err(Refusal::new(
                    "fillet.unsupported-overflow",
                    "the field contact region reaches another support and needs trimmed continuation",
                ));
            }
        }
        if interval.is_some() && end_planes != [true, true] {
            return Err(bad("finite line strips need proven transversal end planes"));
        }
        *interval = span;
        Ok(())
    }

    /// Midpoints of the stations where a support line crosses the source
    /// supports, each with its station interval.
    fn line_samples(
        &self,
        line: &Line,
        curve: CurveId,
        geom: &mut GeomStore,
        source_faces: &[usize],
        bound: &Rational,
        dec: &mut Decider,
    ) -> Result<Samples, Refusal> {
        // A rational lower bound on the length suffices for
        // finite station bracketing. Requiring its square root
        // to be rational excludes oblique planar edges.
        let squared = q(&line.dir.norm_squared())?;
        if squared.sign() != Sign::Positive {
            return Err(bad("a line intersection needs a nonzero direction"));
        }
        let length = squared.min(Rational::one());
        let origin_bound =
            q(&line.origin.x)?.abs() + q(&line.origin.y)?.abs() + q(&line.origin.z)?.abs();
        let extent = (bound + &origin_bound) / length;
        let mut stations = vec![-&extent, extent];
        for &index in source_faces {
            let face = &self.faces[index];
            let surface = geom.add_surface(face.surface.clone());
            if let CurveSurface::Meets(hits) =
                catallax_intersect::intersect_curve_surface(geom, curve, surface, dec)?.value
            {
                for hit in hits {
                    let p = geom.position(hit.point);
                    stations.push(q(
                        &(p.minus(&line.origin).dot(&line.dir) / line.dir.norm_squared())
                    )?);
                }
            }
        }
        stations.sort();
        stations.dedup();
        Ok(stations
            .windows(2)
            .map(|pair| {
                (
                    line.origin.plus(
                        &line
                            .dir
                            .scaled(&Real::from((&pair[0] + &pair[1]) / Rational::from_i64(2))),
                    ),
                    Some([pair[0].clone(), pair[1].clone()]),
                )
            })
            .collect())
    }

    pub(crate) fn fillet(
        &mut self,
        root: usize,
        source_faces: &[usize],
        selection: &EdgeSelection,
        settings: Settings<'_>,
        step: &StepId,
        dec: &mut Decider,
    ) -> Result<usize, Refusal> {
        let Settings { radius, closure } = settings;
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
        let mut runs = Vec::new();
        for faces in pairs {
            let supports = faces.map(|i| self.faces[i].surface.clone());
            let elementary = matches!(
                (&supports[0], &supports[1]),
                (Surface::Plane(_), Surface::Plane(_) | Surface::Cylinder(_))
                    | (Surface::Cylinder(_), Surface::Plane(_))
            );
            let named = matches!(selection, EdgeSelection::Refs(_));
            if !elementary && !rolling::admits(&supports, named, dec)? {
                if matches!(selection, EdgeSelection::All) {
                    continue;
                }
                return Err(bad(
                    "field fillets need a plane with a plane, cylinder, cone or sphere, two cylinders, a cylinder with a sphere, or two coaxial cylinders, cones or spheres",
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
                // The elementary route keeps the lines and circles it sections
                // today; the rolling route takes the rest.
                if !elementary || !matches!(edge, Curve::Line(_) | Curve::Circle(_)) {
                    if let Some(found) = self.rolling_component(
                        root,
                        faces,
                        &supports,
                        &edge,
                        curve,
                        &mut geom,
                        source_faces,
                        &bound,
                        radius,
                        step,
                        dec,
                    )? {
                        runs.extend(found);
                    }
                    continue;
                }
                let samples = match &edge {
                    Curve::Circle(c) => vec![(
                        c.center
                            .plus(&unit(&c.ref_dir).scaled(&square_root(&c.radius_sq))),
                        None,
                    )],
                    Curve::Line(line) => {
                        self.line_samples(line, curve, &mut geom, source_faces, &bound, dec)?
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
                    // The field has no outer-sheet or sphere-strip section;
                    // keep the refusal the elementary reduction gave before
                    // those carriers existed.
                    if matches!(
                        definition.surface,
                        Surface::TorusSheet(_) | Surface::Sphere(_)
                    ) {
                        return Err(Refusal::new(
                            "fillet.unsupported-contact",
                            "this contact requires a horn or spindle torus and its sheet continuation",
                        ));
                    }
                    definitions.push((definition, interval));
                }
            }
        }
        if definitions.is_empty()
            && runs.is_empty()
            && !matches!(selection, EdgeSelection::Refs(refs) if refs.is_empty())
        {
            return Err(Refusal::new(
                "fillet.lost-edge",
                "no selected support intersection is a sharp input boundary",
            ));
        }
        if let EdgeSelection::Refs(refs) = selection {
            for reference in refs {
                let named = |origin: &Origin| {
                    let Role::Blend { edge } = &origin.role else {
                        return false;
                    };
                    let [a, b] = edge.between();
                    (reference.between[0].matches(a) && reference.between[1].matches(b))
                        || (reference.between[1].matches(a) && reference.between[0].matches(b))
                };
                if !definitions
                    .iter()
                    .any(|(definition, _)| named(&definition.origin))
                    && !runs.iter().any(|run| named(run.piece.origin()))
                {
                    return Err(Refusal::new(
                        "fillet.lost-edge",
                        "an edge reference has no surviving sharp input component",
                    ));
                }
            }
        }
        // Straight runs that share one name, which the elementary route
        // refuses, become generator pieces; they register no named faces.
        let mut kept = Vec::with_capacity(definitions.len());
        let shared = |definition: &FilletDefinition,
                      all: &[(FilletDefinition, Option<[Rational; 2]>)]| {
            all.iter()
                .filter(|(other, _)| other.origin == definition.origin)
                .count()
                > 1
        };
        let split = definitions
            .iter()
            .any(|(definition, _)| shared(definition, &definitions))
            && definitions.iter().all(|(definition, interval)| {
                !shared(definition, &definitions)
                    || (matches!(definition.edge, Curve::Line(_)) && interval.is_some())
            });
        if split {
            let all = std::mem::take(&mut definitions);
            for (definition, interval) in &all {
                if !shared(definition, &all) {
                    kept.push((definition.clone(), interval.clone()));
                    continue;
                }
                runs.push(self.converted_run(definition, interval, source_faces, step, dec)?);
            }
            definitions = kept;
        }
        let mut endpoints: Vec<(Point3, Vec<usize>)> = Vec::new();
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
                    for (other, runs) in &mut endpoints {
                        if geometry::sign(&point.minus(other).norm_squared(), dec)? == Sign::Zero {
                            runs.push(i);
                            if runs.len() > 4 {
                                return Err(bad(
                                    "five or more selected runs at a vertex need a general corner patch",
                                ));
                            }
                            found = true;
                            break;
                        }
                    }
                    if !found {
                        endpoints.push((point, vec![i]));
                    }
                }
            }
        }
        let mut result = if let Some(result) =
            self.network_fillet(root, source_faces, &definitions, &endpoints, step, dec)?
        {
            result
        } else if closure == catallax_recipe::FilletCorner::Smooth
            && let Some(result) =
                self.smooth_fillet(root, source_faces, &definitions, &endpoints, step, dec)?
        {
            result
        } else {
            self.elementary_fillets(
                root,
                source_faces,
                definitions,
                &endpoints,
                &mut runs,
                step,
                dec,
            )?
        };
        for run in runs {
            result = self.apply_run(result, run, source_faces, dec)?;
        }
        Ok(result)
    }

    /// The elementary route's strips, rings and trihedral corners. A straight
    /// strip whose section it cannot write exactly, and that meets no other
    /// run, becomes a generator run instead.
    #[allow(clippy::too_many_arguments)]
    fn elementary_fillets(
        &mut self,
        root: usize,
        source_faces: &[usize],
        definitions: Vec<(FilletDefinition, Option<[Rational; 2]>)>,
        endpoints: &[(Point3, Vec<usize>)],
        runs: &mut Vec<rolling::Run>,
        step: &StepId,
        dec: &mut Decider,
    ) -> Result<usize, Refusal> {
        let mut corners = Vec::new();
        for (_, meeting) in endpoints {
            let [a, b, c] = meeting.as_slice() else {
                continue;
            };
            let strips = [*a, *b, *c].map(|i| &definitions[i].0);
            let mut faces: Vec<Origin> = strips
                .iter()
                .flat_map(|strip| &strip.supports)
                .map(|support| {
                    self.faces
                        .iter()
                        .find(|f| f.surface.same_definition(support))
                        .expect("source support")
                        .origin
                        .clone()
                })
                .collect();
            faces.sort();
            faces.dedup();
            let vertex = VertexName::new(faces).map_err(|e| bad(&e.to_string()))?;
            corners.push(CornerDefinition::new(
                Origin {
                    step: step.clone(),
                    role: Role::Corner { vertex },
                },
                strips,
                dec,
            )?);
        }
        let mut result = root;
        for (index, (definition, mut interval)) in definitions.into_iter().enumerate() {
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
            let run = interval.clone();
            let convex = definition.convex;
            let alone = !endpoints
                .iter()
                .any(|(_, meeting)| meeting.len() > 1 && meeting.contains(&index));
            let fallback = (alone && matches!(definition.edge, Curve::Line(_)))
                .then(|| (definition.clone(), run));
            // A concave ring that another support reaches may meet mouths in
            // its plane; the rim route runs only where this one refuses.
            let ring = (!convex && matches!(definition.edge, Curve::Circle(_)))
                .then(|| definition.clone());
            let built = self
                .fillet_domain(
                    &definition,
                    faces,
                    &mut interval,
                    source_faces,
                    &mut GeomStore::new(),
                    dec,
                )
                .and_then(|()| region(definition, interval, step, dec));
            let (feature, faces) = match (built, fallback) {
                (Ok(found), _) => found,
                (Err(refusal), Some((definition, run)))
                    if refusal.code == "field.fillet-contact" =>
                {
                    runs.push(self.converted_run(&definition, &run, source_faces, step, dec)?);
                    continue;
                }
                (Err(refusal), _) if refusal.code == "fillet.unsupported-overflow" => match &ring {
                    Some(ring) => match self.rim_run(root, ring, source_faces, step, dec)? {
                        Some(found) => {
                            runs.push(found);
                            continue;
                        }
                        None => return Err(refusal),
                    },
                    None => return Err(refusal),
                },
                (Err(refusal), _) => return Err(refusal),
            };
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
        for corner in corners {
            let node = self.corner(&corner, step, dec)?;
            result = self.boolean(
                if corner.sector.convex {
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

    /// The material a trihedral corner removes or adds beyond its three strips:
    /// the cone of contact directions at the ball centre, short of all three
    /// supports, less the ball. That cell has the centre, the three feet, the
    /// vertex and one point on each edge as vertices, so it lies within the
    /// centre-to-vertex distance.
    fn corner(
        &mut self,
        corner: &CornerDefinition,
        step: &StepId,
        dec: &mut Decider,
    ) -> Result<usize, Refusal> {
        let radius = number(q(&corner.sector.radius)?);
        let centre = corner.centre();
        let mut planes = Vec::with_capacity(6);
        for (i, (normal, dual)) in corner
            .sector
            .directions
            .iter()
            .zip(corner.sector.dual_basis())
            .enumerate()
        {
            let next = &corner.sector.directions[(i + 1) % 3];
            planes.push(Surface::Plane(Plane {
                origin: centre.plus(&normal.scaled(&corner.sector.radius)),
                normal: normal.clone(),
                ref_dir: normal.cross(next),
            }));
            planes.push(Surface::Plane(Plane {
                origin: centre.clone(),
                normal: dual.scaled(&Real::from_i64(-1)),
                ref_dir: next.clone(),
            }));
        }
        let reach = coordinates(&corner.vertex)?
            .iter()
            .zip(coordinates(centre)?)
            .fold(Rational::zero(), |sum, (v, c)| sum + (v - &c).abs());
        let bound = coordinates(centre)?
            .into_iter()
            .fold(reach, |bound, c| bound + c.abs());
        let cell = geometry::Leaf::Cell(geometry::Cell::new(planes, Some(bound))?);
        let (ball, mut faces) = geometry::compile(
            &catallax_recipe::Operation::Sphere {
                radius,
                pose: catallax_recipe::Pose {
                    rotation: [1, 0, 0, 0].map(|v| number(Rational::from_i64(v))),
                    translation: catallax_recipe::Vector3(coordinates(centre)?.map(number)),
                },
            },
            step,
            dec,
        )?;
        for face in &mut faces {
            face.origin = corner.origin.clone();
        }
        let cell = self.leaf(cell, Vec::new(), dec)?;
        let ball = self.leaf(ball, faces, dec)?;
        self.boolean(BooleanKind::Subtract, cell, ball, dec)
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
                // A planar rail keeps the section frame rational. A rail on a
                // cylinder can sit an irrational distance from the edge.
                let flat = match &definition.supports[0] {
                    Surface::Cylinder(_) => b,
                    _ => a,
                };
                let x = unit(&flat.origin.minus(&edge.origin));
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
    // The support curves between the edge and each contact. An arc centre is
    // the support axis in this section.
    let walls = match definition.supports.each_ref().map(|support| match support {
        Surface::Cylinder(c) if sweep.is_none() => {
            let foot = c.origin.plus(
                &c.axis
                    .scaled(&(origin.minus(&c.origin).dot(&c.axis) / c.axis.norm_squared())),
            );
            local(&foot).map(Some)
        }
        _ => Ok(None),
    }) {
        [a, b] => [a?, b?],
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
    let id = |s: &str| ElementId::new(s).expect("fixed section id");
    let p2 = |p: [Rational; 2]| catallax_recipe::Point2(p.map(number));
    let wall = |name: &'static str,
                centre: &Option<[Rational; 2]>,
                from: &[Rational; 2],
                to: [Rational; 2]| {
        match centre {
            None => catallax_recipe::Segment::Line {
                id: id(name),
                to: p2(to),
            },
            Some(centre) => {
                let start = [&from[0] - &centre[0], &from[1] - &centre[1]];
                let end = [&to[0] - &centre[0], &to[1] - &centre[1]];
                catallax_recipe::Segment::Arc {
                    id: id(name),
                    centre: p2(centre.clone()),
                    to: p2(to),
                    turn: if (&start[0] * &end[1] - &start[1] * &end[0]).sign() == Sign::Positive {
                        Turn::Ccw
                    } else {
                        Turn::Cw
                    },
                }
            }
        }
    };
    let loops = vec![catallax_recipe::Loop(vec![
        wall("a", &walls[0], &e, a.clone()),
        catallax_recipe::Segment::Arc {
            id: id("blend"),
            centre: p2(centre),
            to: p2(b.clone()),
            turn,
        },
        wall("b", &walls[1], &b, e),
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

/// A finite round primitive can miss a concave boss-foot strip even when its
/// infinite cone/cap supports cross the torus. The strip lies in the radial
/// band [body radius, body radius + rolling radius], including its two rails.
/// A parallel round primitive lies within its maximum-radius tube. Strict
/// radial separation of these closed sets excludes every face owned by that
/// leaf, independently of Boolean trimming. Equality does not prove clearance.
/// This does NOT admit a touched mouth: rim-contact continuation still refuses.
fn finite_round_clear(
    definition: &FilletDefinition,
    leaf: &geometry::Leaf,
    dec: &mut Decider,
) -> Result<bool, Refusal> {
    let (false, Curve::Circle(edge), Surface::Torus(torus), geometry::Leaf::Primitive(primitive)) = (
        definition.convex,
        &definition.edge,
        &definition.surface,
        leaf,
    ) else {
        return Ok(false);
    };
    if primitive.kind != crate::metric::Kind::Round
        || !matches!(
            (&definition.supports[0], &definition.supports[1]),
            (Surface::Plane(_), Surface::Cylinder(_)) | (Surface::Cylinder(_), Surface::Plane(_))
        )
    {
        return Ok(false);
    }
    dec.charge_work(1)?;
    let body = square_root(&edge.radius_sq);
    let rail = &body + &definition.radius;
    if geometry::sign(&(&torus.major - &rail), dec)? != Sign::Zero {
        return Ok(false);
    }
    let map = primitive.inverse_placement().inverse();
    let axis = map.apply_vector(&Vector3::new(Real::zero(), Real::zero(), Real::one()));
    if geometry::sign(&axis.cross(&edge.normal).norm_squared(), dec)? != Sign::Zero {
        return Ok(false);
    }
    let origin = map.apply_point(&Point3::new(Real::zero(), Real::zero(), Real::zero()));
    let delta = origin.minus(&edge.center);
    let offset_sq = &delta.cross(&edge.normal).norm_squared() / &edge.normal.norm_squared();
    let radius = Real::from(
        primitive.parameters[0]
            .clone()
            .max(primitive.parameters[1].clone()),
    );
    if geometry::sign(&(&offset_sq - &(&rail + &radius).square()), dec)? == Sign::Positive {
        return Ok(true);
    }
    let gap = &body - &radius;
    Ok(geometry::sign(&gap, dec)? == Sign::Positive
        && geometry::sign(&(&gap.square() - &offset_sq), dec)? == Sign::Positive)
}

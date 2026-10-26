use crate::{
    IgesError, IgesOptions, IgesReport, exact as ex,
    format::{Graph, Value, lengths, text},
};
use catallax_numbers::{Decider, Outcome, Rational, Real, sin_cos_degrees};
use catallax_shapes::{Circle, Curve, Point3, Surface, SurfaceId, Vector3};
use catallax_topology::{BodyId, BodyKind, EdgeId, FaceId, Loop, Model, Sense, ShellId, VertexId};
use std::collections::{BTreeMap, BTreeSet};

struct Exporter<'a> {
    model: &'a Model,
    dec: &'a mut Decider,
    graph: Graph,
    vertices: BTreeMap<VertexId, usize>,
    points: Vec<[Rational; 3]>,
    edges: BTreeMap<EdgeId, (usize, bool)>,
    edge_rows: Vec<Value>,
    surfaces: BTreeMap<SurfaceId, usize>,
    vertex_list: usize,
    edge_list: usize,
    bound: Rational,
    trim_error: Rational,
    parameter_error: Rational,
    added: usize,
    seams: usize,
    faces: usize,
}

pub(crate) fn prepare(
    model: &Model,
    bodies: &[BodyId],
    options: &IgesOptions,
    dec: &mut Decider,
) -> Result<(String, IgesReport), IgesError> {
    let mut seen = BTreeSet::new();
    for body in bodies {
        if !seen.insert(*body) {
            return Err(IgesError::InvalidTopology("duplicate body".into()));
        }
        if !model.body(*body).marks.is_empty() {
            return Err(IgesError::UncertifiedModel);
        }
        let faults = model.check(*body);
        if !faults.is_empty() {
            return Err(IgesError::InvalidTopology(format!("{faults:?}")));
        }
        // Validate kinds before attempting numeric conversion.
        for f in model.faces_of_body(*body) {
            let surface = model.geom().surface(model.face(f).surface);
            if !matches!(
                surface,
                Surface::Plane(_) | Surface::Cylinder(_) | Surface::Cone(_)
            ) {
                return Err(IgesError::UnsupportedSurface(surface.kind()));
            }
        }
        for e in model.edges_of_body(*body) {
            let curve = model.geom().curve(model.edge(e).curve);
            if !matches!(curve, Curve::Line(_) | Curve::Circle(_)) {
                return Err(IgesError::UnsupportedCurve(curve.kind()));
            }
        }
    }
    let mut graph = Graph::default();
    let vertex_list = graph.add(502, 1, vec![]);
    let edge_list = graph.add(504, 1, vec![]);
    let mut e = Exporter {
        model,
        dec,
        graph,
        vertex_list,
        edge_list,
        vertices: BTreeMap::new(),
        points: vec![],
        edges: BTreeMap::new(),
        edge_rows: vec![],
        surfaces: BTreeMap::new(),
        bound: Rational::zero(),
        trim_error: Rational::zero(),
        parameter_error: Rational::zero(),
        added: 0,
        seams: 0,
        faces: 0,
    };
    let mut solids = 0;
    for body in bodies {
        match model.body(*body).kind {
            BodyKind::Empty => continue,
            BodyKind::Solid => {}
            _ => return Err(IgesError::InvalidTopology("not a solid".into())),
        }
        for lump in &model.body(*body).lumps {
            let shells = &model.lump(*lump).shells;
            let outer = e.shell(shells[0])?;
            let mut values = vec![outer.into(), 1usize.into(), (shells.len() - 1).into()];
            for shell in &shells[1..] {
                values.push(e.shell(*shell)?.into());
                values.push(1usize.into());
            }
            let id = e.graph.add(186, 0, values);
            e.graph.entities[(id - 1) / 2].root = true;
            solids += 1;
        }
    }
    let mut vertices = vec![e.points.len().into()];
    for p in &e.points {
        vertices.extend(lengths(p.clone()));
    }
    e.graph.entities[(vertex_list - 1) / 2].values = vertices;
    let edge_count = e.edge_rows.len() / 5;
    let mut edges = vec![edge_count.into()];
    edges.append(&mut e.edge_rows);
    e.graph.entities[(edge_list - 1) / 2].values = edges;
    // All metric parameters share one rational scale. Their denominators divide D,
    // so every q * D * decimal(1/D) is an exactly representable finite decimal.
    let denominator = ex::denominator(e.graph.lengths().map(|q| q / options.units.millimetres()));
    let scale = ex::scale(denominator, options.scale_digits)?;
    let trim_deviation = &e.trim_error + &e.parameter_error;
    let deviation = &(&(&scale - &Rational::one()).abs() * &e.bound) + &(&scale * &trim_deviation);
    if deviation > options.tolerance {
        return Err(IgesError::ToleranceExceeded {
            proven: Box::new(deviation),
            requested: Box::new(options.tolerance.clone()),
        });
    }
    let file = e.graph.file(&scale, options.units)?;
    let report = IgesReport {
        max_deviation: deviation,
        trim_deviation,
        scale,
        coordinate_bound: e.bound,
        units: options.units,
        solids,
        faces: e.faces,
        edges: edge_count,
        added_edges: e.seams,
        vertices: e.points.len(),
        added_vertices: e.added,
    };
    Ok((file, report))
}

impl Exporter<'_> {
    fn include(&mut self, p: &[Rational; 3]) {
        self.bound = self.bound.clone().max(ex::bound(p));
    }
    fn add_vertex(&mut self, p: [Rational; 3]) -> usize {
        self.include(&p);
        self.points.push(p);
        self.points.len()
    }
    fn vertex(&mut self, v: VertexId) -> Result<usize, IgesError> {
        if let Some(i) = self.vertices.get(&v) {
            return Ok(*i);
        }
        let p = ex::point(
            &self.model.geom().position(self.model.vertex(v).point),
            self.dec,
        )?;
        let index = self.add_vertex(p);
        self.vertices.insert(v, index);
        Ok(index)
    }
    fn point(&mut self, p: &Point3) -> Result<usize, IgesError> {
        let p = ex::point(p, self.dec)?;
        let mut values = lengths(p);
        values.push(0usize.into());
        Ok(self.graph.add(116, 0, values))
    }
    fn direction(&mut self, v: &Vector3) -> Result<usize, IgesError> {
        let v = ex::integer_vector(ex::vector(v, self.dec)?);
        if ex::dot(&v, &v).is_zero() {
            return Err(IgesError::UnrepresentableDefinition("zero direction"));
        }
        let values = v
            .iter()
            .map(ex::decimal)
            .collect::<Result<Vec<_>, _>>()?
            .into_iter()
            .map(text)
            .collect();
        Ok(self.graph.add(123, 0, values))
    }
    fn surface(&mut self, id: SurfaceId) -> Result<usize, IgesError> {
        if let Some(i) = self.surfaces.get(&id) {
            return Ok(*i);
        }
        let entity = match self.model.geom().surface(id).clone() {
            Surface::Plane(s) => {
                let p = self.point(&s.origin)?;
                let n = self.direction(&s.normal)?;
                let x = self.direction(&s.ref_dir)?;
                self.graph.add(190, 1, vec![p.into(), n.into(), x.into()])
            }
            Surface::Cylinder(s) => {
                let p = self.point(&s.origin)?;
                let n = self.direction(&s.axis)?;
                let r = ex::rational(&s.radius(), self.dec)?;
                let x = self.direction(&s.ref_dir)?;
                self.graph
                    .add(192, 1, vec![p.into(), n.into(), Value::Length(r), x.into()])
            }
            Surface::Cone(s) => {
                let p = self.point(&s.apex)?;
                let n = self.direction(&s.axis)?;
                let angle = if matches!(
                    self.dec.compare(&s.slope, &Real::one()),
                    Outcome::ProvenCoincident(_)
                ) {
                    Rational::from_i64(45)
                } else {
                    if s.slope.as_rational().is_some() {
                        return Err(IgesError::UnrepresentableDefinition(
                            "rational cone slope other than one is not covered",
                        ));
                    }
                    // A candidate degree angle is accepted only after an exact
                    // trigonometric identity proof. No atan enclosure is a proof.
                    let angle = self.dec.enclose_at(&s.slope, 128)?.atan_checked(self.dec)?;
                    let degrees = (angle
                        * catallax_numbers::Ball::from_rational(&Rational::from_i64(180), 128))
                    .checked_div(&catallax_numbers::Ball::pi_checked(128, self.dec)?)
                    .unwrap();
                    let (a, b) = degrees.mid().to_decimal_bounds(3);
                    let mut found = None;
                    for t in [a, b] {
                        let q = Rational::parse(&t).unwrap();
                        let [sin, cos] = sin_cos_degrees(&q, self.dec)?;
                        if matches!(
                            self.dec.compare(&sin, &(&s.slope * &cos)),
                            Outcome::ProvenCoincident(_)
                        ) {
                            found = Some(q);
                            break;
                        }
                    }
                    found.ok_or(IgesError::UnrepresentableDefinition(
                        "cone semi-angle has no proved finite-decimal degree encoding",
                    ))?
                };
                {
                    let x = self.direction(&s.ref_dir)?;
                    self.graph.add(
                        194,
                        1,
                        vec![
                            p.into(),
                            n.into(),
                            Value::Length(Rational::zero()),
                            text(ex::decimal(&angle)?),
                            x.into(),
                        ],
                    )
                }
            }
            s => return Err(IgesError::UnsupportedSurface(s.kind())),
        };
        self.surfaces.insert(id, entity);
        Ok(entity)
    }
    fn edge(&mut self, id: EdgeId) -> Result<(usize, bool), IgesError> {
        if let Some(x) = self.edges.get(&id) {
            return Ok(*x);
        }
        let edge = self.model.edge(id).clone();
        let curve = self.model.geom().curve(edge.curve).clone();
        let (a, b) = match (edge.start, edge.end) {
            (Some(a), Some(b)) => (self.vertex(a)?, self.vertex(b)?),
            (None, None) => {
                let Curve::Circle(ref c) = curve else {
                    return Err(IgesError::UnrepresentableDefinition("unbounded line trim"));
                };
                let p = ex::point(
                    &Curve::Circle(c.clone()).parameter_origin(self.model.geom()),
                    self.dec,
                )?;
                let i = self.add_vertex(p);
                self.added += 1;
                (i, i)
            }
            _ => return Err(IgesError::InvalidTopology("one-ended edge".into())),
        };
        let (curve, same) = match curve {
            Curve::Line(_) => {
                if a == b {
                    return Err(IgesError::UnrepresentableDefinition("collapsed line"));
                }
                let mut values = lengths(self.points[a - 1].clone());
                values.extend(lengths(self.points[b - 1].clone()));
                (self.graph.add(110, 0, values), true)
            }
            Curve::Circle(c) => (
                self.circle(&c, a, b, edge.sense == Sense::Same)?,
                edge.sense == Sense::Same,
            ),
            c => return Err(IgesError::UnsupportedCurve(c.kind())),
        };
        // IGES EdgeList follows the geometric curve, not the Catallax edge sense.
        let (start, end) = if same { (a, b) } else { (b, a) };
        let index = self.edge_rows.len() / 5 + 1;
        self.edge_rows.extend([
            curve.into(),
            self.vertex_list.into(),
            start.into(),
            self.vertex_list.into(),
            end.into(),
        ]);
        self.edges.insert(id, (index, same));
        Ok((index, same))
    }
    fn circle(&mut self, c: &Circle, a: usize, b: usize, same: bool) -> Result<usize, IgesError> {
        let center = ex::point(&c.center, self.dec)?;
        let mut x = ex::unit(&c.ref_dir, self.dec)?;
        let z = ex::unit(&c.normal, self.dec)?;
        let radius = ex::rational(&c.radius_sq.sqrt(), self.dec)?;
        if radius <= Rational::zero() {
            return Err(IgesError::UnrepresentableDefinition(
                "nonpositive circle radius",
            ));
        }
        if a != b {
            // A type-100 arc's native range starts at its local start angle.
            // Give finite arcs a zero start so [0,T] pcurves have the same
            // range even when the retained circle uses another reference ray.
            let start = if same { a } else { b };
            x = ex::sub(&self.points[start - 1], &center).map(|v| &v / &radius);
            if ex::dot(&x, &x) != Rational::one() || !ex::dot(&x, &z).is_zero() {
                return Err(IgesError::UnrepresentableDefinition(
                    "arc start is not on its circle",
                ));
            }
        }
        let y = ex::cross(&z, &x);
        self.bound = self
            .bound
            .clone()
            .max(ex::bound(&center) + &Rational::from_i64(3) * &radius);
        // Keep the matrix orthonormal. Some independent readers discard scale
        // on circular-arc placements; a nonterminating frame fails closed.
        let mut matrix = Vec::new();
        for i in 0..3 {
            matrix.extend([
                text(ex::decimal(&x[i])?),
                text(ex::decimal(&y[i])?),
                text(ex::decimal(&z[i])?),
                Value::Length(center[i].clone()),
            ]);
        }
        let matrix = self.graph.add(124, 0, matrix);
        let (a, b) = if same { (a, b) } else { (b, a) };
        let local = |v: usize| {
            let d = ex::sub(&self.points[v - 1], &center);
            [ex::dot(&d, &x), ex::dot(&d, &y)]
        };
        let from = local(a);
        let to = local(b);
        let id = self.graph.add(
            100,
            0,
            vec![
                Value::Length(Rational::zero()),
                Value::Length(Rational::zero()),
                Value::Length(Rational::zero()),
                Value::Length(from[0].clone()),
                Value::Length(from[1].clone()),
                Value::Length(to[0].clone()),
                Value::Length(to[1].clone()),
            ],
        );
        self.graph.entities[(id - 1) / 2].transform = matrix;
        Ok(id)
    }
    /// A finite axial band has one rectangular, unwrapped parameter wire.
    /// Recover its branch from the retained directed arcs, never from a fit.
    fn axial_sector(&mut self, f: &catallax_topology::Face) -> Result<usize, IgesError> {
        let bad = || {
            IgesError::UnrepresentableDefinition(
                "axial sector requires two circular arcs and two generators with cardinal endpoints",
            )
        };
        if f.loops.len() != 1 {
            return Err(bad());
        }
        let Loop::Coedges(mut cs) = self.model.loop_(f.loops[0]).clone() else {
            return Err(bad());
        };
        if cs.len() != 4 {
            return Err(bad());
        }
        if f.sense == Sense::Reversed {
            cs.reverse();
        }
        let (origin, axis, reference, radius, slope) = match self.model.geom().surface(f.surface) {
            Surface::Cylinder(c) => (
                ex::point(&c.origin, self.dec)?,
                ex::unit(&c.axis, self.dec)?,
                ex::unit(&c.ref_dir, self.dec)?,
                Some(ex::rational(&c.radius(), self.dec)?),
                None,
            ),
            Surface::Cone(c) => (
                ex::point(&c.apex, self.dec)?,
                ex::unit(&c.axis, self.dec)?,
                ex::unit(&c.ref_dir, self.dec)?,
                None,
                Some(ex::rational(&c.slope, self.dec)?),
            ),
            _ => unreachable!(),
        };
        if !ex::dot(&axis, &reference).is_zero() {
            return Err(bad());
        }
        let transverse = ex::cross(&axis, &reference);
        let mut uses = Vec::new();
        let mut levels = BTreeMap::new();
        for coedge in cs {
            let c = self.model.coedge(coedge).clone();
            let edge = self.model.edge(c.edge).clone();
            let (Some(start), Some(end)) = (edge.start, edge.end) else {
                return Err(bad());
            };
            let (index, same) = self.edge(c.edge)?;
            let orientation = ((c.sense == Sense::Same) == same) == (f.sense == Sense::Same);
            let forward = (c.sense == Sense::Same) == (f.sense == Sense::Same);
            let (a, b) = if forward { (start, end) } else { (end, start) };
            let (a, b) = (self.vertex(a)?, self.vertex(b)?);
            for v in [a, b] {
                if levels.contains_key(&v) {
                    continue;
                }
                let d = ex::sub(&self.points[v - 1], &origin);
                let height = ex::dot(&d, &axis);
                let radial = std::array::from_fn(|i| &d[i] - &(&height * &axis[i]));
                let r = if let Some(slope) = &slope {
                    if height <= Rational::zero() || slope <= &Rational::zero() {
                        return Err(bad()); // Apex and opposite nappe are outside this slice.
                    }
                    &height * slope
                } else {
                    radius.clone().unwrap()
                };
                if r <= Rational::zero() || ex::dot(&radial, &radial) != &r * &r {
                    return Err(bad());
                }
                let x = ex::dot(&radial, &reference);
                let y = ex::dot(&radial, &transverse);
                let u = if x == r && y.is_zero() {
                    0i64
                } else if x.is_zero() && y == r {
                    90
                } else if x == -r.clone() && y.is_zero() {
                    180
                } else if x.is_zero() && y == -r {
                    270
                } else {
                    return Err(bad());
                };
                let level = if slope.is_some() {
                    let root = Real::from(ex::dot(&d, &d)).sqrt();
                    let ball = self.dec.enclose_at(&root, 160)?;
                    let (lower, _) = ball.mid().to_decimal_bounds(25);
                    let q = Rational::parse(&lower).unwrap();
                    let mut error = (&q - &ball.lo()).abs().max((&q - &ball.hi()).abs());
                    match self.dec.compare(&root, &Real::from(q.clone())) {
                        Outcome::ProvenCoincident(_) => error = Rational::zero(),
                        Outcome::Refused(r) => return Err(r.into()),
                        _ => {}
                    }
                    self.trim_error = self.trim_error.clone().max(error);
                    q
                } else {
                    height.clone()
                };
                levels.insert(v, (u, level, height));
            }
            let (ua, va, ha) = &levels[&a];
            let (ub, vb, hb) = &levels[&b];
            let (delta, end) = match self.model.geom().curve(edge.curve) {
                Curve::Circle(circle) => {
                    let center = ex::point(&circle.center, self.dec)?;
                    let normal = ex::vector(&circle.normal, self.dec)?;
                    let expected = std::array::from_fn(|i| &origin[i] + &(ha * &axis[i]));
                    if ha != hb
                        || center != expected
                        || !ex::cross(&normal, &axis).iter().all(Rational::is_zero)
                    {
                        return Err(bad());
                    }
                    let increasing = (ex::dot(&normal, &axis) > Rational::zero()) == orientation;
                    let delta = if increasing {
                        (ub - ua).rem_euclid(360)
                    } else {
                        -((ua - ub).rem_euclid(360))
                    };
                    if delta == 0 {
                        return Err(bad());
                    }
                    let pi = catallax_numbers::Ball::pi_checked(160, self.dec)?;
                    let exact = pi
                        * catallax_numbers::Ball::from_rational(
                            &Rational::from_ratio(delta.abs(), 180).unwrap(),
                            160,
                        );
                    let (t, _) = exact.mid().to_decimal_bounds(17);
                    let t = Rational::parse(&t).unwrap();
                    let error = (&t - &exact.lo()).abs().max((&t - &exact.hi()).abs());
                    self.parameter_error = self.parameter_error.clone().max(&self.bound * &error);
                    (delta, t)
                }
                Curve::Line(_) => {
                    if ua != ub || va == vb || ha == hb {
                        return Err(bad());
                    }
                    // Use the same bounded [0,1] range for the shared 3D
                    // generator and its pcurve, including reversed coedges.
                    let Value::Text(curve) = &self.edge_rows[(index - 1) * 5] else {
                        unreachable!()
                    };
                    let curve = curve.parse::<usize>().unwrap();
                    let entity = &mut self.graph.entities[(curve - 1) / 2];
                    if entity.kind == 110 {
                        let mut segment = vec![
                            text(1),
                            text(1),
                            text(0),
                            text(0),
                            text(1),
                            text(0),
                            text(0),
                            text(0),
                            text(1),
                            text(1),
                            text(1),
                            text(1),
                        ];
                        segment.append(&mut entity.values);
                        segment.extend([text(0), text(1), text(0), text(0), text(0)]);
                        entity.kind = 126;
                        entity.values = segment;
                    }
                    (0, Rational::one())
                }
                _ => return Err(bad()),
            };
            uses.push((index, orientation, a, b, delta, end));
        }
        // Two opposing arcs and two generators must close one finite rectangle.
        if (0..4).any(|i| {
            (uses[i].4 == 0) == (uses[(i + 1) % 4].4 == 0) || uses[i].3 != uses[(i + 1) % 4].2
        }) || uses.iter().map(|x| x.4).sum::<i64>() != 0
        {
            return Err(bad());
        }
        let mut u = levels[&uses[0].2].0;
        let mut values = vec![4usize.into()];
        for (index, orientation, a, b, delta, end) in uses {
            let next = u + delta;
            let pc = self.graph.add(
                126,
                0,
                vec![
                    text(1),
                    text(1),
                    text(1),
                    text(0),
                    text(1),
                    text(0),
                    text(0),
                    text(0),
                    text(ex::decimal(&end)?),
                    text(ex::decimal(&end)?),
                    text(1),
                    text(1),
                    text(u),
                    Value::Length(levels[&a].1.clone()),
                    text(0),
                    text(next),
                    Value::Length(levels[&b].1.clone()),
                    text(0),
                    text(0),
                    text(ex::decimal(&end)?),
                    text(0),
                    text(0),
                    text(1),
                ],
            );
            values.extend([
                0usize.into(),
                self.edge_list.into(),
                index.into(),
                usize::from(orientation).into(),
                1usize.into(),
                1usize.into(),
                pc.into(),
            ]);
            u = next;
        }
        Ok(self.graph.add(508, 1, values))
    }
    fn face(&mut self, id: FaceId) -> Result<usize, IgesError> {
        let f = self.model.face(id).clone();
        let mut loops = Vec::new();
        if f.loops.is_empty() {
            return Err(IgesError::UnrepresentableDefinition("unbounded face"));
        }
        // A periodic band uses one seam twice in opposite senses, joining the
        // two circular boundaries into a single closed parameter-space wire.
        let rings: Option<Vec<_>> = f
            .loops
            .iter()
            .map(|lp| {
                let Loop::Coedges(cs) = self.model.loop_(*lp) else {
                    return None;
                };
                if cs.len() != 1 {
                    return None;
                }
                let c = self.model.coedge(cs[0]).clone();
                let e = self.model.edge(c.edge);
                if e.start.is_some() || e.end.is_some() {
                    return None;
                }
                Some(c)
            })
            .collect();
        let periodic = matches!(
            self.model.geom().surface(f.surface),
            Surface::Cylinder(_) | Surface::Cone(_)
        );
        if let Some(rings) = rings.filter(|r| periodic && r.len() == 2) {
            let mut uses = Vec::new();
            let mut vertices = Vec::new();
            for c in rings {
                let (index, same) = self.edge(c.edge)?;
                let Value::Text(v) = &self.edge_rows[(index - 1) * 5 + 2] else {
                    unreachable!()
                };
                vertices.push(v.parse::<usize>().unwrap());
                let orient = ((c.sense == Sense::Same) == same) == (f.sense == Sense::Same);
                uses.push((index, orient));
            }
            let [a, b] = [vertices[0], vertices[1]];
            let pa = &self.points[a - 1];
            let pb = &self.points[b - 1];
            let residual = match self.model.geom().surface(f.surface) {
                Surface::Cylinder(c) => {
                    ex::cross(&ex::sub(pa, pb), &ex::vector(&c.axis, self.dec)?)
                }
                Surface::Cone(c) => {
                    let p = ex::point(&c.apex, self.dec)?;
                    ex::cross(&ex::sub(pa, &p), &ex::sub(pb, &p))
                }
                _ => unreachable!(),
            };
            if residual.iter().any(|x| !x.is_zero()) {
                return Err(IgesError::UnrepresentableDefinition(
                    "ring vertices do not share a proved generatrix",
                ));
            }
            let mut line = lengths(pa.clone());
            line.extend(lengths(pb.clone()));
            let mut segment = vec![
                text(1),
                text(1),
                text(0),
                text(0),
                text(1),
                text(0),
                text(0),
                text(0),
                text(1),
                text(1),
                text(1),
                text(1),
            ];
            segment.extend(line);
            segment.extend([text(0), text(1), text(0), text(0), text(0)]);
            let curve = self.graph.add(126, 0, segment);
            let seam = self.edge_rows.len() / 5 + 1;
            self.edge_rows.extend([
                curve.into(),
                self.vertex_list.into(),
                a.into(),
                self.vertex_list.into(),
                b.into(),
            ]);
            self.seams += 1;
            // IGES 192/194 use degrees for U. Integer 0 and 360 encode a
            // full period exactly, without a rounded value of pi.
            let mut levels = Vec::new();
            let (origin, axis, reference) = match self.model.geom().surface(f.surface) {
                Surface::Cylinder(c) => (&c.origin, &c.axis, &c.ref_dir),
                Surface::Cone(c) => (&c.apex, &c.axis, &c.ref_dir),
                _ => unreachable!(),
            };
            let origin = ex::point(origin, self.dec)?;
            let axis = ex::unit(axis, self.dec)?;
            let reference = ex::unit(reference, self.dec)?;
            for v in [a, b] {
                let d = ex::sub(&self.points[v - 1], &origin);
                let height = ex::dot(&d, &axis);
                let radial = std::array::from_fn(|i| &d[i] - &(&height * &axis[i]));
                if ex::cross(&radial, &reference).iter().any(|v| !v.is_zero())
                    || ex::dot(&radial, &reference) <= Rational::zero()
                {
                    return Err(IgesError::UnrepresentableDefinition(
                        "periodic seam differs from the retained reference ray",
                    ));
                }
                let level = if matches!(self.model.geom().surface(f.surface), Surface::Cylinder(_))
                {
                    height
                } else {
                    let squared = ex::dot(&d, &d);
                    let root = Real::from(squared).sqrt();
                    let ball = self.dec.enclose_at(&root, 160)?;
                    let (lower, _) = ball.mid().to_decimal_bounds(25);
                    let q = Rational::parse(&lower).unwrap();
                    let mut error = (&q - &ball.lo()).abs().max((&q - &ball.hi()).abs());
                    match self.dec.compare(&root, &Real::from(q.clone())) {
                        Outcome::ProvenCoincident(_) => error = Rational::zero(),
                        Outcome::Refused(refusal) => return Err(refusal.into()),
                        _ => {}
                    }
                    self.trim_error = self.trim_error.clone().max(error);
                    q
                };
                levels.push(level);
            }
            if levels[0] == levels[1] {
                return Err(IgesError::UnrepresentableDefinition(
                    "rounded periodic domain collapsed",
                ));
            }
            let Loop::Coedges(cs) = self.model.loop_(f.loops[0]) else {
                unreachable!()
            };
            let Curve::Circle(first_circle) = self
                .model
                .geom()
                .curve(self.model.edge(self.model.coedge(cs[0]).edge).curve)
            else {
                unreachable!()
            };
            let circle_normal = ex::vector(&first_circle.normal, self.dec)?;
            let agrees = ex::dot(&circle_normal, &axis) > Rational::zero();
            let ua = if uses[0].1 == agrees { 0 } else { 360 };
            let ub = 360 - ua;
            let mut values = vec![4usize.into()];
            for ((index, orientation), (u0, v0, u1, v1)) in [
                (uses[0], (ua, levels[0].clone(), ub, levels[0].clone())),
                ((seam, true), (ub, levels[0].clone(), ub, levels[1].clone())),
                (uses[1], (ub, levels[1].clone(), ua, levels[1].clone())),
                (
                    (seam, false),
                    (ua, levels[1].clone(), ua, levels[0].clone()),
                ),
            ] {
                // Parameter curves follow the LOOP direction (IGES 508).
                // A bounded degree-one curve keeps a positive parameter range
                // when the reader reverses a coedge; an infinite type-110 line
                // reverses its range to [-L,0] and breaks shared SameRange data.
                let end = if index == seam {
                    Rational::one()
                } else {
                    let pi = catallax_numbers::Ball::pi_checked(160, self.dec)?;
                    let tau = &pi + &pi;
                    let (t, _) = tau.mid().to_decimal_bounds(17);
                    let t = Rational::parse(&t).unwrap();
                    let angle_error = (&t - &tau.lo()).abs().max((&t - &tau.hi()).abs());
                    self.parameter_error =
                        self.parameter_error.clone().max(&self.bound * &angle_error);
                    t
                };
                let end = ex::decimal(&end)?;
                let pc = self.graph.add(
                    126,
                    0,
                    vec![
                        text(1),
                        text(1),
                        text(1),
                        text(0),
                        text(1),
                        text(0),
                        text(0),
                        text(0),
                        text(&end),
                        text(&end),
                        text(1),
                        text(1),
                        text(u0),
                        Value::Length(v0),
                        text(0),
                        text(u1),
                        Value::Length(v1),
                        text(0),
                        text(0),
                        text(end),
                        text(0),
                        text(0),
                        text(1),
                    ],
                );
                values.extend([
                    0usize.into(),
                    self.edge_list.into(),
                    index.into(),
                    usize::from(orientation).into(),
                    1usize.into(),
                    1usize.into(),
                    pc.into(),
                ]);
            }
            loops.push(self.graph.add(508, 1, values));
        } else if periodic {
            loops.push(self.axial_sector(&f)?);
        } else {
            for lp in &f.loops {
                let mut values = vec![];
                match self.model.loop_(*lp).clone() {
                    Loop::Vertex(v) => {
                        let v = self.vertex(v)?;
                        values.extend([
                            1usize.into(),
                            1usize.into(),
                            self.vertex_list.into(),
                            v.into(),
                            1usize.into(),
                            0usize.into(),
                        ]);
                    }
                    Loop::Coedges(mut cs) => {
                        if f.sense == Sense::Reversed {
                            cs.reverse();
                        }
                        values.push(cs.len().into());
                        for c in cs {
                            let c = self.model.coedge(c).clone();
                            let (e, same) = self.edge(c.edge)?;
                            let orient =
                                ((c.sense == Sense::Same) == same) == (f.sense == Sense::Same);
                            values.extend([
                                0usize.into(),
                                self.edge_list.into(),
                                e.into(),
                                usize::from(orient).into(),
                                0usize.into(),
                            ]);
                        }
                    }
                }
                loops.push(self.graph.add(508, 1, values));
            }
        }
        let s = self.surface(f.surface)?;
        // Catallax boundaries have no distinguished outer loop; IGES permits
        // flag zero. Their orientation selects the finite face and its holes.
        let mut values = vec![s.into(), loops.len().into(), 0usize.into()];
        values.extend(loops.into_iter().map(Value::from));
        self.faces += 1;
        Ok(self.graph.add(510, 1, values))
    }
    fn shell(&mut self, id: ShellId) -> Result<usize, IgesError> {
        let faces = self.model.shell(id).faces.clone();
        let mut values = vec![faces.len().into()];
        for f in faces {
            values.push(self.face(f)?.into());
            values.push(usize::from(self.model.face(f).sense == Sense::Same).into());
        }
        Ok(self.graph.add(514, 1, values))
    }
}

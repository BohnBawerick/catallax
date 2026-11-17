use std::collections::BTreeMap;

use catallax_numbers::{Decider, Rational};

use crate::geometry_source::{ParsedStep, Topology};
use crate::part21::{Entity, Graph, Param};
use crate::{Curve, Placement, Surface, Vec3};

use super::{Error, TolerantSolid, invalid};
mod uv;

struct Builder<'a> {
    parsed: &'a ParsedStep,
    changes: Graph,
    next: u64,
}
impl Builder<'_> {
    fn add(&mut self, kind: &str, parameters: Vec<Param>) -> Result<u64, Error> {
        self.complex(vec![(kind.into(), parameters)])
    }
    fn complex(&mut self, parts: Vec<(String, Vec<Param>)>) -> Result<u64, Error> {
        let id = self.next;
        self.next = self.next.checked_add(1).ok_or_else(|| {
            invalid(
                "import.tolerant-broken-reference",
                "STEP entity index overflow",
            )
        })?;
        self.changes.insert(id, Entity { parts });
        Ok(id)
    }
    fn set(&mut self, id: u64, kind: &str, parameters: Vec<Param>) -> Result<(), Error> {
        let value = Entity {
            parts: vec![(kind.into(), parameters)],
        };
        if self.changes.get(&id).is_some_and(|prior| *prior != value) {
            return Err(invalid(
                "import.tolerant-shared-repair-conflict",
                format!("independent repairs require different bindings for shared entity #{id}"),
            ));
        }
        self.changes.insert(id, value);
        Ok(())
    }
    fn name(&self, source_name: &str) -> Param {
        let source = source_name
            .strip_prefix('#')
            .and_then(|value| value.parse().ok())
            .and_then(|id| self.parsed.source_graph().get(&id));
        let value = source.and_then(|entity| {
            entity
                .part("REPRESENTATION_ITEM")
                .and_then(|parameters| parameters.first())
                .or_else(|| {
                    entity.parts.iter().find_map(|(_, parameters)| {
                        parameters
                            .first()
                            .filter(|value| matches!(value, Param::Str(_)))
                    })
                })
        });
        value
            .filter(|value| matches!(value, Param::Str(_)))
            .cloned()
            .unwrap_or_else(|| text(source_name))
    }
    fn point(&mut self, position: &Vec3, scale: &Rational) -> Result<u64, Error> {
        let coordinates = position
            .iter()
            .map(|coordinate| decimal(&(coordinate / scale)))
            .collect::<Result<Vec<_>, _>>()?;
        self.add("CARTESIAN_POINT", vec![text(""), Param::List(coordinates)])
    }
    fn direction(&mut self, direction: &Vec3) -> Result<u64, Error> {
        self.add(
            "DIRECTION",
            vec![
                text(""),
                Param::List(direction.iter().map(decimal).collect::<Result<_, _>>()?),
            ],
        )
    }
    fn placement(&mut self, placement: &Placement, scale: &Rational) -> Result<u64, Error> {
        let origin = self.point(&placement.origin, scale)?;
        let axis = self.direction(&placement.axis)?;
        let reference = self.direction(&placement.ref_dir)?;
        self.add(
            "AXIS2_PLACEMENT_3D",
            vec![
                text(""),
                reference_to(origin),
                reference_to(axis),
                reference_to(reference),
            ],
        )
    }
    fn source_curve(&self, name: &str, dec: &mut Decider) -> Result<Option<u64>, Error> {
        let Some(entity) = name.strip_prefix('#').and_then(|value| value.parse().ok()) else {
            return Ok(None);
        };
        Ok(match self.parsed.topology(entity, dec)? {
            Topology::Edge { curve, .. } => curve,
            _ => None,
        })
    }
    fn source_surface(&self, name: &str, dec: &mut Decider) -> Result<Option<u64>, Error> {
        let Some(entity) = name.strip_prefix('#').and_then(|value| value.parse().ok()) else {
            return Ok(None);
        };
        Ok(match self.parsed.topology(entity, dec)? {
            Topology::Face { surface, .. } => surface,
            _ => None,
        })
    }
    fn shell_id(&self, mut id: u64, dec: &mut Decider) -> Result<u64, Error> {
        for _ in 0..32 {
            match self.parsed.topology(id, dec)? {
                Topology::OrientedShell(parent, _) => id = parent,
                Topology::Shell { .. } => return Ok(id),
                _ => {
                    return Err(invalid(
                        "import.tolerant-broken-reference",
                        "source shell reference is not a shell",
                    ));
                }
            }
        }
        Err(invalid(
            "import.tolerant-broken-reference",
            "cyclic source shell reference",
        ))
    }
}

fn text(value: &str) -> Param {
    Param::Str(value.into())
}
fn reference_to(id: u64) -> Param {
    Param::Ref(id)
}
fn boolean(value: bool) -> Param {
    Param::Enum(if value { "T" } else { "F" }.into())
}
fn references(values: impl IntoIterator<Item = u64>) -> Param {
    Param::List(values.into_iter().map(reference_to).collect())
}
fn decimal(value: &Rational) -> Result<Param, Error> {
    let number = value.to_f64_nearest();
    if !number.is_finite() {
        return Err(invalid(
            "import.tolerant-broken-geometry",
            "non-finite tolerant STEP coordinate",
        ));
    }
    Ok(Param::Number(format!("{number:.17e}")))
}

pub(super) struct Rebuilt {
    pub changes: Graph,
    pub body_roots: Vec<Vec<u64>>,
    pub checks: Vec<super::SampledCheck>,
}

pub(super) fn rebuild(
    import: &super::TolerantImport,
    bodies: &[TolerantSolid],
    dec: &mut Decider,
) -> Result<Rebuilt, Error> {
    let parsed = import.parsed.as_ref();
    let scale = &import.source.mm_per_unit;
    let next = parsed
        .entities()
        .last()
        .unwrap_or(0)
        .checked_add(1)
        .ok_or_else(|| {
            invalid(
                "import.tolerant-broken-reference",
                "STEP entity index overflow",
            )
        })?;
    let mut builder = Builder {
        parsed,
        changes: Graph::new(),
        next,
    };
    let mut expanded_roots = BTreeMap::<u64, usize>::new();
    let mut body_roots = Vec::with_capacity(bodies.len());
    let mut uv = uv::Writer::new(import);
    for solid in bodies {
        let body = &solid.body;
        let mut trims = BTreeMap::<&str, Vec<&super::NativeTrim>>::new();
        for trim in import
            .evidence
            .native_trims
            .iter()
            .filter(|trim| trim.body == body.name)
        {
            trims.entry(&trim.edge).or_default().push(trim);
        }
        let root: u64 = body
            .name
            .strip_prefix('#')
            .and_then(|value| value.parse().ok())
            .ok_or_else(|| {
                invalid(
                    "import.tolerant-broken-reference",
                    "native body lacks a source root",
                )
            })?;
        if !matches!(
            parsed.topology(root, dec)?,
            Topology::Solid { .. } | Topology::SurfaceModel(_) | Topology::FaceModel(_)
        ) {
            return Err(invalid(
                "import.tolerant-broken-topology",
                "native body root has no finite shells",
            ));
        }
        let mut active = vec![false; body.edges.len()];
        let mut active_vertices = vec![false; body.vertices.len()];
        for shell in &body.shells {
            for face in &shell.faces {
                for boundary in &face.loops {
                    if let Some(vertex) = boundary.vertex {
                        active_vertices[vertex] = true;
                    }
                    for usage in &boundary.edges {
                        active[usage.edge] = true;
                        for vertex in body.edges[usage.edge].vertices {
                            active_vertices[vertex] = true;
                        }
                    }
                }
            }
        }
        let mut vertices = vec![None; body.vertices.len()];
        for (index, vertex) in body.vertices.iter().enumerate() {
            if !active_vertices[index] {
                continue;
            }
            dec.charge_work(1)?;
            let point = builder.point(&vertex.source.point, scale)?;
            let id = builder.add(
                "VERTEX_POINT",
                vec![builder.name(&vertex.source.name), reference_to(point)],
            )?;
            vertices[index] = Some(id);
        }
        let mut edges = vec![None; body.edges.len()];
        for (index, edge) in body.edges.iter().enumerate() {
            if !active[index] {
                continue;
            }
            dec.charge_work(1)?;
            let endpoints = edge
                .vertices
                .map(|vertex| vertices[vertex].expect("active endpoint"));
            let curve = match builder.source_curve(&edge.name, dec)? {
                Some(curve) => curve,
                None => match &edge.curve {
                    Curve::Line(origin, direction) => {
                        let origin = builder.point(origin, scale)?;
                        let direction = builder.direction(direction)?;
                        let vector = builder.add(
                            "VECTOR",
                            vec![
                                text(""),
                                reference_to(direction),
                                Param::Number("1.".into()),
                            ],
                        )?;
                        builder.add(
                            "LINE",
                            vec![text(""), reference_to(origin), reference_to(vector)],
                        )?
                    }
                    Curve::BSpline(spline) if edge.degenerate => {
                        let point = builder.point(&spline.geometry.controls()[0], scale)?;
                        builder.add(
                            "B_SPLINE_CURVE_WITH_KNOTS",
                            vec![
                                builder.name(&edge.name),
                                Param::Number("1".into()),
                                references([point, point]),
                                Param::Enum("UNSPECIFIED".into()),
                                boolean(true),
                                Param::Enum("U".into()),
                                Param::List(vec![Param::Number("2".into()); 2]),
                                Param::List(vec![
                                    Param::Number("0.".into()),
                                    Param::Number("1.".into()),
                                ]),
                                Param::Enum("UNSPECIFIED".into()),
                            ],
                        )?
                    }
                    _ => {
                        return Err(invalid(
                            "import.tolerant-broken-reference",
                            "native edge has no exportable source carrier",
                        ));
                    }
                },
            };
            let curve = uv.edge_geometry(
                &mut builder,
                edge,
                curve,
                trims.get(edge.name.as_str()).map_or(&[], Vec::as_slice),
                dec,
            )?;
            let parameters = vec![
                builder.name(&edge.name),
                reference_to(endpoints[0]),
                reference_to(endpoints[1]),
                reference_to(curve),
                boolean(edge.curve_same_sense),
            ];
            edges[index] = Some(builder.add("EDGE_CURVE", parameters)?);
        }
        let mut shells = Vec::with_capacity(body.shells.len());
        for shell in &body.shells {
            let mut faces = Vec::with_capacity(shell.faces.len());
            for face in &shell.faces {
                dec.charge_work(1)?;
                let surface = match builder.source_surface(&face.name, dec)? {
                    Some(surface) => surface,
                    None => match &face.surface {
                        Surface::Plane(placement) => {
                            let placement = builder.placement(placement, scale)?;
                            builder.add("PLANE", vec![text(""), reference_to(placement)])?
                        }
                        _ => {
                            return Err(invalid(
                                "import.tolerant-broken-reference",
                                "native face has no exportable source surface",
                            ));
                        }
                    },
                };
                let mut bounds = Vec::with_capacity(face.loops.len());
                for boundary in &face.loops {
                    let loop_id = if let Some(vertex) = boundary.vertex {
                        builder.add(
                            "VERTEX_LOOP",
                            vec![
                                text(""),
                                reference_to(vertices[vertex].expect("active vertex loop")),
                            ],
                        )?
                    } else {
                        let mut uses = Vec::with_capacity(boundary.edges.len());
                        for usage in &boundary.edges {
                            let edge = edges[usage.edge].ok_or_else(|| {
                                invalid(
                                    "import.tolerant-broken-reference",
                                    "native trim references an inactive edge",
                                )
                            })?;
                            uses.push(builder.add(
                                "ORIENTED_EDGE",
                                vec![
                                    text(""),
                                    Param::Derived,
                                    Param::Derived,
                                    reference_to(edge),
                                    boolean(usage.forward),
                                ],
                            )?);
                        }
                        builder.add("EDGE_LOOP", vec![text(""), references(uses)])?
                    };
                    bounds.push(builder.add(
                        if boundary.outer {
                            "FACE_OUTER_BOUND"
                        } else {
                            "FACE_BOUND"
                        },
                        vec![
                            text(""),
                            reference_to(loop_id),
                            boolean(boundary.orientation),
                        ],
                    )?);
                }
                let parameters = vec![
                    builder.name(&face.name),
                    references(bounds),
                    reference_to(surface),
                    boolean(face.same_sense),
                ];
                faces.push(builder.add("ADVANCED_FACE", parameters)?);
            }
            let source: u64 = shell
                .name
                .strip_prefix('#')
                .and_then(|name| name.parse().ok())
                .ok_or_else(|| {
                    invalid(
                        "import.tolerant-broken-reference",
                        "native shell lacks a source identity",
                    )
                })?;
            let source = builder.shell_id(source, dec)?;
            let entity = builder.add(
                "CLOSED_SHELL",
                vec![builder.name(&format!("#{source}")), references(faces)],
            )?;
            shells.push(entity);
        }
        if body.declared_solid && shells.len() > 1 {
            let mut voids = Vec::with_capacity(shells.len() - 1);
            for &shell in &shells[1..] {
                voids.push(builder.add(
                    "ORIENTED_CLOSED_SHELL",
                    vec![text(""), Param::Derived, reference_to(shell), boolean(true)],
                )?);
            }
            builder.set(
                root,
                "BREP_WITH_VOIDS",
                vec![
                    builder.name(&body.name),
                    reference_to(shells[0]),
                    references(voids),
                ],
            )?;
            body_roots.push(vec![root]);
        } else {
            let mut roots = Vec::with_capacity(shells.len());
            for (index, shell) in shells.into_iter().enumerate() {
                let parameters = vec![builder.name(&body.name), reference_to(shell)];
                roots.push(if index == 0 {
                    builder.set(root, "MANIFOLD_SOLID_BREP", parameters)?;
                    root
                } else {
                    builder.add("MANIFOLD_SOLID_BREP", parameters)?
                });
            }
            if roots.len() > 1 {
                expanded_roots.insert(root, body_roots.len());
            }
            body_roots.push(roots);
        }
    }
    // Preserve the original occurrence/assembly graph. Only representation items
    // which held a promoted multi-shell surface set need item-list expansion.
    for (&entity, source) in parsed.source_graph() {
        for (kind, parameters) in &source.parts {
            if !kind.ends_with("REPRESENTATION") || kind == "DEFINITIONAL_REPRESENTATION" {
                continue;
            }
            let Some(Param::List(items)) = parameters.get(1) else {
                continue;
            };
            let mut changed = false;
            let mut replacement = Vec::new();
            for item in items {
                if let Param::Ref(root) = item
                    && let Some(roots) = expanded_roots.get(root)
                {
                    replacement.extend(body_roots[*roots].iter().copied().map(reference_to));
                    changed = true;
                } else {
                    replacement.push(item.clone());
                }
            }
            let promote = matches!(kind.as_str(), "MANIFOLD_SURFACE_SHAPE_REPRESENTATION" | "GEOMETRICALLY_BOUNDED_SURFACE_SHAPE_REPRESENTATION")
                && items.iter().any(|item| matches!(item, Param::Ref(root) if bodies.iter().any(|body| body.name().strip_prefix('#').and_then(|name|name.parse::<u64>().ok()) == Some(*root))));
            if changed || promote {
                let mut value = source.clone();
                for (name, parameters) in &mut value.parts {
                    if name == kind {
                        parameters[1] = Param::List(replacement.clone());
                        if promote {
                            *name = "ADVANCED_BREP_SHAPE_REPRESENTATION".into();
                        }
                    }
                }
                builder.changes.insert(entity, value);
            }
        }
    }
    Ok(Rebuilt {
        changes: builder.changes,
        body_roots,
        checks: uv.into_checks(),
    })
}

pub(super) fn standalone_frame(
    parsed: &ParsedStep,
    changes: &mut Graph,
    name: &str,
    roots: &[u64],
    dec: &mut Decider,
) -> Result<u64, Error> {
    if roots.is_empty() {
        return Err(invalid(
            "import.tolerant-broken-topology",
            "selected body has no reconstructed solid root",
        ));
    }
    let next = parsed
        .entities()
        .chain(changes.keys().copied())
        .max()
        .unwrap_or(0)
        .checked_add(1)
        .ok_or_else(|| {
            invalid(
                "import.tolerant-broken-reference",
                "STEP entity index overflow",
            )
        })?;
    let mut builder = Builder {
        parsed,
        changes: std::mem::take(changes),
        next,
    };
    let mut context = None;
    let mut units = Vec::new();
    for entity in parsed.entities() {
        dec.charge_work(1)?;
        let kinds = parsed.entity_types(entity)?;
        if kinds.contains(&"GEOMETRIC_REPRESENTATION_CONTEXT") {
            let candidate = parsed.representation_context(entity)?;
            if candidate.dimension == Some(3) && !candidate.units.is_empty() {
                context = Some(entity);
                break;
            }
        }
        if kinds.iter().any(|kind| {
            matches!(
                *kind,
                "LENGTH_UNIT" | "PLANE_ANGLE_UNIT" | "SOLID_ANGLE_UNIT"
            )
        }) {
            units.push(entity);
        }
    }
    let context = if let Some(context) = context {
        context
    } else {
        if units.is_empty() {
            let unit = builder.complex(vec![
                ("LENGTH_UNIT".into(), vec![]),
                ("NAMED_UNIT".into(), vec![Param::Derived]),
                (
                    "SI_UNIT".into(),
                    vec![Param::Enum("MILLI".into()), Param::Enum("METRE".into())],
                ),
            ])?;
            units.push(unit);
            let unit = builder.complex(vec![
                ("NAMED_UNIT".into(), vec![Param::Derived]),
                ("PLANE_ANGLE_UNIT".into(), vec![]),
                (
                    "SI_UNIT".into(),
                    vec![Param::Unset, Param::Enum("RADIAN".into())],
                ),
            ])?;
            units.push(unit);
        }
        builder.complex(vec![
            (
                "GEOMETRIC_REPRESENTATION_CONTEXT".into(),
                vec![Param::Number("3".into())],
            ),
            (
                "GLOBAL_UNIT_ASSIGNED_CONTEXT".into(),
                vec![references(units)],
            ),
            ("REPRESENTATION_CONTEXT".into(), vec![text(""), text("")]),
        ])?
    };
    let application = builder.add("APPLICATION_CONTEXT", vec![text("automotive_design")])?;
    let product_context = builder.add(
        "PRODUCT_CONTEXT",
        vec![text(""), reference_to(application), text("mechanical")],
    )?;
    let product = builder.add(
        "PRODUCT",
        vec![
            text(name),
            text(name),
            text("tolerant body-local geometry"),
            references([product_context]),
        ],
    )?;
    let formation = builder.add(
        "PRODUCT_DEFINITION_FORMATION",
        vec![text(""), text(""), reference_to(product)],
    )?;
    let definition_context = builder.add(
        "PRODUCT_DEFINITION_CONTEXT",
        vec![
            text("part definition"),
            reference_to(application),
            text("design"),
        ],
    )?;
    let definition = builder.add(
        "PRODUCT_DEFINITION",
        vec![
            text("design"),
            text(""),
            reference_to(formation),
            reference_to(definition_context),
        ],
    )?;
    let shape = builder.add(
        "PRODUCT_DEFINITION_SHAPE",
        vec![text(""), text(""), reference_to(definition)],
    )?;
    let representation = builder.add(
        "ADVANCED_BREP_SHAPE_REPRESENTATION",
        vec![
            text(name),
            references(roots.iter().copied()),
            reference_to(context),
        ],
    )?;
    let root = builder.add(
        "SHAPE_DEFINITION_REPRESENTATION",
        vec![reference_to(shape), reference_to(representation)],
    )?;
    *changes = builder.changes;
    Ok(root)
}

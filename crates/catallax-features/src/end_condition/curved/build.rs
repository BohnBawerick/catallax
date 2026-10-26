//! Analytic convex columns and the existing regularised Boolean engine.

use super::*;
use catallax_numbers::{Fact, FactKind, FactOrigin};
use catallax_recipe::{ElementId, Operation, Origin, Role, StepId};
use catallax_shapes::{Curve, GeomKey, Line, PointDef, Sphere};
use catallax_topology::{
    Body, BodyId, BodyKind, Coedge, Edge, Face, Loop, Lump, Model, Sense, Shell, Vertex,
};
use std::collections::BTreeMap;

/// A finite analytic intersection with its certified enclosing sphere.
#[derive(Clone, Debug)]
pub struct ConvexCell {
    surfaces: Vec<Surface>,
    roles: Vec<Role>,
    radius: Rational,
}

impl ConvexCell {
    /// Inside is the nonpositive side of every exact support polynomial.
    pub fn surfaces(&self) -> &[Surface] {
        &self.surfaces
    }
    /// Construction-face contributions, in support order.
    pub fn roles(&self) -> &[Role] {
        &self.roles
    }
    /// Certified enclosing sphere about world zero.
    pub fn radius(&self) -> &Rational {
        &self.radius
    }
}

/// Independent analytic cells for the same regularised cap definition.
#[derive(Clone, Debug)]
pub struct CapConstruction {
    /// The finite exact profile column, including its artificial upper cap.
    pub column: ConvexCell,
    /// The closed trimmed stock target, renamed by this extrusion.
    pub target: ConvexCell,
    /// A finite negative-midpoint halfspace, needed for an entry cap.
    pub clipper: Option<ConvexCell>,
}

fn prism_cell(
    base: &Point3,
    vectors: &[Vector3; 3],
    roles: Vec<Role>,
    dec: &mut Decider,
) -> Result<ConvexCell, Refusal> {
    let orientation = sign(&vectors[0].cross(&vectors[1]).dot(&vectors[2]), dec)?;
    if orientation == Sign::Zero {
        return Err(coverage("construction prism is degenerate"));
    }
    let radius = prism_radius(base, vectors, dec)?;
    let mut surfaces = Vec::new();
    for i in 0..3 {
        let mut n = vectors[(i + 1) % 3].cross(&vectors[(i + 2) % 3]);
        if orientation == Sign::Negative {
            n = n.scaled(&Real::from_i64(-1));
        }
        for upper in [false, true] {
            let p = if upper {
                base.plus(&vectors[i])
            } else {
                base.clone()
            };
            let normal = if upper {
                n.clone()
            } else {
                n.scaled(&Real::from_i64(-1))
            };
            surfaces.push(Surface::Plane(Plane::new(
                p,
                normal,
                vectors[(i + 1) % 3].clone(),
                dec,
            )?));
        }
    }
    Ok(ConvexCell {
        surfaces,
        roles,
        radius,
    })
}

fn prism_radius(
    base: &Point3,
    vectors: &[Vector3; 3],
    dec: &mut Decider,
) -> Result<Rational, Refusal> {
    let mut radius = Rational::zero();
    for corner in 0..8 {
        let mut p = base.clone();
        for (i, v) in vectors.iter().enumerate() {
            if corner & (1 << i) != 0 {
                p = p.plus(v);
            }
        }
        let mut upper = Rational::zero();
        for x in [&p.x, &p.y, &p.z] {
            let b = projected::enclose(x, dec)?;
            upper = upper + b.lo().abs().max(b.hi().abs());
        }
        radius = radius.max(upper);
    }
    Ok(radius)
}

fn clipping_prism(
    cap: &CurvedCapped,
    dec: &mut Decider,
) -> Result<(Point3, [Vector3; 3]), Refusal> {
    let n = &cap.midpoint.normal;
    let nn = n.norm_squared();
    let center = cap.midpoint.origin.minus(&Point3::origin()).dot(n);
    let p = Point3::origin().plus(&n.scaled(&(center / nn)));
    let u = projected::perpendicular(n, dec)?;
    let v = n.cross(&u);
    let bound =
        cap.scene.radius.clone() + frame_radius(&cap.frame, &cap.bounds, dec)? + Rational::one();
    let scaled = [u, v, n.clone()]
        .map(|v| -> Result<Vector3, Refusal> {
            let length = projected::enclose(&v.norm_squared().sqrt(), dec)?.lo();
            if length.sign() != Sign::Positive {
                return Err(coverage("construction frame norm is not separated"));
            }
            Ok(v.scaled(&Real::from(&bound / &length)))
        })
        .into_iter()
        .collect::<Result<Vec<_>, _>>()?;
    let base = p
        .plus(&scaled[0].scaled(&Real::from_i64(-1)))
        .plus(&scaled[1].scaled(&Real::from_i64(-1)))
        .plus(&scaled[2].scaled(&Real::from_i64(-2)));
    Ok((
        base,
        std::array::from_fn(|i| scaled[i].scaled(&Real::from_i64(2))),
    ))
}

impl CurvedCapped {
    fn polygon_sphere_source(
        &self,
        model: &mut Model,
        dec: &mut Decider,
    ) -> Result<BodyId, Refusal> {
        let Surface::Sphere(sphere) = &self.scene.surface else {
            unreachable!("a stock sphere operation has a spherical support")
        };
        let w = &self.frame.direction;
        let pole = projected::perpendicular(w, dec)?.cross(w);
        let seam = if self.branch == Branch::Entry {
            w.clone()
        } else {
            w.scaled(&Real::from_i64(-1))
        };
        // At the certified root, grad(F).w is -sqrt(D) for entry and
        // +sqrt(D) for exit. Hence every retained cap normal has a strictly
        // negative component along `seam`. The poles (perpendicular to w)
        // and positive seam meridian are outside the entire cap. Only the
        // sphere's chart changes: its implicit polynomial, centre and radius
        // are identical to the independently compiled analytic definition.
        let surface = model.geom_mut().add_surface(Surface::Sphere(Sphere::new(
            sphere.center.clone(),
            sphere.radius.clone(),
            pole,
            seam,
            dec,
        )?));
        let face = model.add_face(Face {
            surface,
            sense: Sense::Same,
            loops: Vec::new(),
        });
        let shell = model.add_shell(Shell { faces: vec![face] });
        let lump = model.add_lump(Lump {
            shells: vec![shell],
        });
        Ok(model.add_body(Body {
            kind: BodyKind::Solid,
            lumps: vec![lump],
            marks: Default::default(),
        }))
    }

    /// Exposes analytic construction cells, not B-rep faces or containment.
    /// Entry is `(column intersect clipper) subtract target`; exit is
    /// `column intersect target`. The certificate excludes every artificial cap.
    pub fn construction(
        &self,
        sides: &[ElementId],
        dec: &mut Decider,
    ) -> Result<CapConstruction, Refusal> {
        self.check_side_ids(sides)?;
        let base = self.frame.point(&[
            self.bounds[0][0].clone().into(),
            self.bounds[1][0].clone().into(),
        ]);
        let vectors = [
            self.frame
                .u
                .scaled(&Real::from(&self.bounds[0][1] - &self.bounds[0][0])),
            self.frame
                .v
                .scaled(&Real::from(&self.bounds[1][1] - &self.bounds[1][0])),
            self.frame
                .direction
                .scaled(&Real::from(self.extent.clone())),
        ];
        let mut roles: Vec<_> = sides
            .iter()
            .map(|id| Role::Side {
                segment: id.clone(),
                path: None,
            })
            .collect();
        roles.extend([Role::Start, Role::End]);
        let column = if let Some(polygon) = &self.polygon {
            ConvexCell {
                surfaces: polygon.surfaces(&self.frame, &self.extent, dec)?,
                roles,
                // The whole enclosing parallelepiped also bounds the polygon.
                radius: prism_radius(&base, &vectors, dec)?,
            }
        } else {
            prism_cell(&base, &vectors, roles, dec)?
        };
        let mut surfaces = vec![self.scene.surface.clone()];
        surfaces.extend(self.scene.trims.iter().cloned().map(Surface::Plane));
        let target = ConvexCell {
            roles: vec![Role::End; surfaces.len()],
            surfaces,
            radius: self.scene.radius.clone(),
        };
        let clipper = if self.branch == Branch::Entry {
            let (base, vectors) = clipping_prism(self, dec)?;
            Some(prism_cell(&base, &vectors, vec![Role::End; 6], dec)?)
        } else {
            None
        };
        Ok(CapConstruction {
            column,
            target,
            clipper,
        })
    }
}

fn named(model: &mut Model, source: BodyId, origin: Origin) -> BodyId {
    let old = model.body(source).clone();
    let lumps = old
        .lumps
        .into_iter()
        .map(|id| {
            let old = model.lump(id).clone();
            let shells = old
                .shells
                .into_iter()
                .map(|id| {
                    let old = model.shell(id).clone();
                    let faces = old
                        .faces
                        .into_iter()
                        .map(|id| {
                            model.add_named_face(model.face(id).clone(), vec![origin.clone()])
                        })
                        .collect();
                    model.add_shell(Shell { faces })
                })
                .collect();
            model.add_lump(Lump { shells })
        })
        .collect();
    model.add_body(Body { lumps, ..old })
}

/// Constructs a parallelepiped from three independent edge vectors. All
/// incidences recorded here follow from its shared exact affine definition.
fn prism(
    model: &mut Model,
    base: Point3,
    vectors: [Vector3; 3],
    names: [Vec<Origin>; 6],
    dec: &mut Decider,
) -> Result<BodyId, Refusal> {
    let determinant = vectors[0].cross(&vectors[1]).dot(&vectors[2]);
    let orientation = sign(&determinant, dec)?;
    if orientation == Sign::Zero {
        return Err(coverage("construction prism is degenerate"));
    }
    let points: Vec<_> = (0..8)
        .map(|i| {
            let mut p = base.clone();
            for (j, v) in vectors.iter().enumerate() {
                if i & (1 << j) != 0 {
                    p = p.plus(v);
                }
            }
            p
        })
        .collect();
    let vertices: Vec<_> = points
        .iter()
        .map(|p| {
            let point = model.geom_mut().add_point(PointDef::Explicit(p.clone()));
            model.add_vertex(Vertex { point })
        })
        .collect();
    let mut edges = BTreeMap::new();
    let mut faces = Vec::new();
    // Outward cycles for a positive determinant, in x-/x+/y-/y+/z-/z+ order.
    for (face_index, mut ids) in [
        [0, 4, 6, 2],
        [1, 3, 7, 5],
        [0, 1, 5, 4],
        [2, 6, 7, 3],
        [0, 2, 3, 1],
        [4, 5, 7, 6],
    ]
    .into_iter()
    .enumerate()
    {
        if orientation == Sign::Negative {
            ids.reverse();
        }
        let a = &points[ids[0]];
        let u = points[ids[1]].minus(a);
        let v = points[ids[2]].minus(a);
        let plane = Plane::new(a.clone(), u.cross(&v), u, dec)?;
        let surface = model.geom_mut().add_surface(Surface::Plane(plane));
        let mut coedges = Vec::new();
        for j in 0..4 {
            let (a, b) = (ids[j], ids[(j + 1) % 4]);
            let key = (a.min(b), a.max(b));
            let edge = *edges.entry(key).or_insert_with(|| {
                let curve = model.geom_mut().add_curve(Curve::Line(Line {
                    origin: points[key.0].clone(),
                    dir: points[key.1].minus(&points[key.0]),
                }));
                for i in [key.0, key.1] {
                    let point = model.vertex(vertices[i]).point;
                    model.geom_mut().record_fact(Fact {
                        kind: FactKind::Incident,
                        a: GeomKey::Point(point),
                        b: GeomKey::Curve(curve),
                        at: None,
                        origin: FactOrigin::Construction("features.curved-prism"),
                    });
                }
                model.add_edge(Edge {
                    curve,
                    sense: Sense::Same,
                    start: Some(vertices[key.0]),
                    end: Some(vertices[key.1]),
                })
            });
            let curve = model.edge(edge).curve;
            model.geom_mut().record_fact(Fact {
                kind: FactKind::Incident,
                a: GeomKey::Curve(curve),
                b: GeomKey::Surface(surface),
                at: None,
                origin: FactOrigin::Construction("features.curved-prism"),
            });
            coedges.push(model.add_coedge(Coedge {
                edge,
                sense: if a < b { Sense::Same } else { Sense::Reversed },
            }));
        }
        let boundary = model.add_loop(Loop::Coedges(coedges));
        faces.push(model.add_named_face(
            Face {
                surface,
                sense: Sense::Same,
                loops: vec![boundary],
            },
            names[face_index].clone(),
        ));
    }
    let shell = model.add_shell(Shell { faces });
    let lump = model.add_lump(Lump {
        shells: vec![shell],
    });
    Ok(model.add_body(Body {
        kind: BodyKind::Solid,
        lumps: vec![lump],
        marks: Default::default(),
    }))
}

impl CurvedCapped {
    fn check_side_ids(&self, sides: &[ElementId]) -> Result<(), Refusal> {
        if sides.len() != self.sides.len() {
            return Err(coverage(
                "construction side ids do not match the exact profile",
            ));
        }
        Ok(())
    }

    /// Builds the exact capped polygon. Side ids follow `side_constraints` and
    /// come from the authored sketch, never subdivision indices. Rectangle
    /// callers can continue to supply their four-element array.
    pub fn build(
        &self,
        model: &mut Model,
        step: &StepId,
        side_ids: impl AsRef<[ElementId]>,
        dec: &mut Decider,
    ) -> Result<BodyId, Refusal> {
        dec.check_deadline()?;
        let side_ids = side_ids.as_ref();
        self.check_side_ids(side_ids)?;
        let origin = |role| Origin {
            step: step.clone(),
            role,
        };
        let tool = if let Some(polygon) = &self.polygon {
            let mut names: Vec<_> = side_ids
                .iter()
                .map(|id| {
                    vec![origin(Role::Side {
                        segment: id.clone(),
                        path: None,
                    })]
                })
                .collect();
            names.extend([vec![origin(Role::Start)], vec![origin(Role::End)]]);
            polygon.build(model, &self.frame, &self.extent, &names, dec)?
        } else {
            let names = std::array::from_fn(|i| {
                vec![origin(if i < 4 {
                    Role::Side {
                        segment: side_ids[i].clone(),
                        path: None,
                    }
                } else if i == 4 {
                    Role::Start
                } else {
                    Role::End
                })]
            });
            let base = self.frame.point(&[
                self.bounds[0][0].clone().into(),
                self.bounds[1][0].clone().into(),
            ]);
            let vectors = [
                self.frame
                    .u
                    .scaled(&Real::from(&self.bounds[0][1] - &self.bounds[0][0])),
                self.frame
                    .v
                    .scaled(&Real::from(&self.bounds[1][1] - &self.bounds[1][0])),
                self.frame
                    .direction
                    .scaled(&Real::from(self.extent.clone())),
            ];
            prism(model, base, vectors, names, dec)?
        };
        let source = match &self.scene.operation {
            Operation::Sphere { .. } if self.polygon.is_some() => {
                self.polygon_sphere_source(model, dec)?
            }
            Operation::Sphere { radius, pose } => {
                catallax_topology::sphere(model, radius.value().clone(), &scene::placement(pose)?)?
            }
            Operation::Cylinder {
                radius,
                height,
                pose,
            } => catallax_topology::cylinder(
                model,
                radius.value().clone(),
                height.value().clone(),
                &scene::placement(pose)?,
            )?,
            Operation::Cone {
                base_radius,
                top_radius,
                height,
                pose,
            } => {
                let pose = scene::placement(pose)?;
                if base_radius == top_radius {
                    catallax_topology::cylinder(
                        model,
                        base_radius.value().clone(),
                        height.value().clone(),
                        &pose,
                    )?
                } else {
                    catallax_topology::cone(
                        model,
                        base_radius.value().clone(),
                        top_radius.value().clone(),
                        height.value().clone(),
                        &pose,
                    )?
                }
            }
            _ => unreachable!("scene constructor only accepts stock quadrics"),
        };
        let source = named(model, source, origin(Role::End));
        let result = if self.branch == Branch::Exit {
            catallax_boolean::boolean(
                model,
                tool,
                source,
                catallax_boolean::BoolOp::Intersect,
                dec,
            )
            .map_err(|e| e.refusal)?
            .body
        } else {
            // Construct a finite halfspace aligned with the midpoint plane.
            // The plane may be oblique to the authored sweep direction.
            let (base, vectors) = clipping_prism(self, dec)?;
            let clipper = prism(
                model,
                base,
                vectors,
                std::array::from_fn(|_| vec![origin(Role::End)]),
                dec,
            )?;
            let slab = catallax_boolean::boolean(
                model,
                tool,
                clipper,
                catallax_boolean::BoolOp::Intersect,
                dec,
            )
            .map_err(|e| e.refusal)?
            .body;
            catallax_boolean::boolean(model, slab, source, catallax_boolean::BoolOp::Subtract, dec)
                .map_err(|e| e.refusal)?
                .body
        };
        if !model.check(result).is_empty() {
            return Err(Refusal::new(
                "features.invalid-topology",
                "curved cap failed topology validation",
            ));
        }
        Ok(result)
    }
}

pub(super) fn frame_radius(
    frame: &Frame,
    bounds: &[[Rational; 2]; 2],
    dec: &mut Decider,
) -> Result<Rational, Refusal> {
    let mut radius = Rational::zero();
    for u in &bounds[0] {
        for v in &bounds[1] {
            let p = frame.point(&[u.clone().into(), v.clone().into()]);
            let mut upper = Rational::zero();
            for x in [&p.x, &p.y, &p.z] {
                let b = projected::enclose(x, dec)?;
                upper = upper + b.lo().abs().max(b.hi().abs());
            }
            radius = radius.max(upper);
        }
    }
    Ok(radius)
}

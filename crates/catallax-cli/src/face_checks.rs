//! Opt-in face-definition checks for the recipe verification adapter.

use std::collections::BTreeSet;
use std::time::Instant;

use catallax_field::measurement::{CentrePath, FaceDefinition, MeasurementOptions, measure_face};
use catallax_numbers::{Decider, Rational, Real, Refusal};
use catallax_recipe::{Number, OriginPattern, Recipe};
use catallax_shapes::{
    Circle, Cone, Curve, Cylinder, Ellipse, EllipticCylinder, Line, Plane, Point3, Sphere, Surface,
    Torus, Vector3,
};
use serde::Deserialize;
use serde_json::{Value, json};

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
pub(crate) struct FaceChecks {
    pub(crate) format: String,
    pub(crate) recipe: Value,
    definitions: Vec<Attachment>,
    #[serde(default)]
    options: Controls,
    edit: Option<Edit>,
}
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct Attachment {
    face: OriginPattern,
    definition: Definition,
}
#[derive(Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
enum Definition {
    Fillet {
        path: PathInput,
        radius: Number,
    },
    Offset {
        original: SurfaceInput,
        distance: Number,
    },
}
#[derive(Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
enum PathInput {
    Line {
        origin: [Number; 3],
        direction: [Number; 3],
    },
    Circle {
        center: [Number; 3],
        normal: [Number; 3],
        radius: Number,
        reference: [Number; 3],
    },
    Ellipse {
        center: [Number; 3],
        normal: [Number; 3],
        major_direction: [Number; 3],
        major: Number,
        minor: Number,
    },
    /// Uses a branch and exact arc endpoints already owned by the B-rep.
    Edge { between: [OriginPattern; 2] },
    /// Component of two exact outward-offset supports, in intersection order.
    Intersection {
        surfaces: Box<[SurfaceInput; 2]>,
        component: usize,
        from: Option<[Number; 3]>,
        to: Option<[Number; 3]>,
    },
}
#[derive(Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case", deny_unknown_fields)]
enum SurfaceInput {
    Plane {
        origin: [Number; 3],
        normal: [Number; 3],
        reference: [Number; 3],
    },
    Cylinder {
        origin: [Number; 3],
        axis: [Number; 3],
        radius: Number,
        reference: [Number; 3],
    },
    Cone {
        apex: [Number; 3],
        axis: [Number; 3],
        slope: Number,
        reference: [Number; 3],
    },
    Sphere {
        center: [Number; 3],
        axis: [Number; 3],
        radius: Number,
        reference: [Number; 3],
    },
    Torus {
        center: [Number; 3],
        axis: [Number; 3],
        major: Number,
        minor: Number,
        reference: [Number; 3],
    },
    EllipticCylinder {
        origin: [Number; 3],
        axis: [Number; 3],
        major_direction: [Number; 3],
        major: Number,
        minor: Number,
    },
}
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct Controls {
    #[serde(default = "width")]
    width: Number,
    #[serde(default = "cells")]
    max_cells: usize,
    #[serde(default = "segments")]
    path_segments: usize,
}
fn width() -> Number {
    Number::parse("0.0000001").unwrap()
}
fn cells() -> usize {
    4096
}
fn segments() -> usize {
    64
}
impl Default for Controls {
    fn default() -> Self {
        Self {
            width: width(),
            max_cells: cells(),
            path_segments: segments(),
        }
    }
}
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct Edit {
    step: String,
    parameter: String,
    value: Number,
}
fn point(p: &[Number; 3]) -> Point3 {
    Point3::new(real(&p[0]), real(&p[1]), real(&p[2]))
}
fn vector(v: &[Number; 3]) -> Vector3 {
    Vector3::new(real(&v[0]), real(&v[1]), real(&v[2]))
}
fn real(x: &Number) -> Real {
    Real::from_rational(x.value().clone())
}
fn refusal(detail: impl Into<String>) -> Refusal {
    Refusal::new("cli.face-definition", detail)
}
impl SurfaceInput {
    fn build(&self, dec: &mut Decider) -> Result<Surface, Refusal> {
        Ok(match self {
            Self::Plane {
                origin,
                normal,
                reference,
            } => Surface::Plane(Plane::new(
                point(origin),
                vector(normal),
                vector(reference),
                dec,
            )?),
            Self::Cylinder {
                origin,
                axis,
                radius,
                reference,
            } => Surface::Cylinder(Cylinder::new(
                point(origin),
                vector(axis),
                real(radius),
                vector(reference),
                dec,
            )?),
            Self::Cone {
                apex,
                axis,
                slope,
                reference,
            } => Surface::Cone(Cone::new(
                point(apex),
                vector(axis),
                real(slope),
                vector(reference),
                dec,
            )?),
            Self::Sphere {
                center,
                axis,
                radius,
                reference,
            } => Surface::Sphere(Sphere::new(
                point(center),
                real(radius),
                vector(axis),
                vector(reference),
                dec,
            )?),
            Self::Torus {
                center,
                axis,
                major,
                minor,
                reference,
            } => Surface::Torus(Torus::new(
                point(center),
                vector(axis),
                real(major),
                real(minor),
                vector(reference),
                dec,
            )?),
            Self::EllipticCylinder {
                origin,
                axis,
                major_direction,
                major,
                minor,
            } => {
                if minor.value() <= &Rational::zero() || major.value() <= minor.value() {
                    return Err(refusal(
                        "original elliptic cylinder requires 0 < minor < major",
                    ));
                }
                Surface::EllipticCylinder(EllipticCylinder::new(
                    point(origin),
                    vector(axis),
                    vector(major_direction),
                    real(major).square(),
                    real(minor).square(),
                    dec,
                )?)
            }
        })
    }
}
impl PathInput {
    fn build(
        &self,
        built: &mut catallax_features::Built,
        dec: &mut Decider,
    ) -> Result<CentrePath, Refusal> {
        Ok(match self {
            Self::Line { origin, direction } => {
                CentrePath::Line(Line::new(point(origin), vector(direction), dec)?)
            }
            Self::Circle {
                center,
                normal,
                radius,
                reference,
            } => {
                if radius.value() <= &Rational::zero() {
                    return Err(refusal("path radius must be positive"));
                }
                CentrePath::Circle(Circle::new(
                    point(center),
                    vector(normal),
                    real(radius).square(),
                    vector(reference),
                    dec,
                )?)
            }
            Self::Ellipse {
                center,
                normal,
                major_direction,
                major,
                minor,
            } => {
                if minor.value() <= &Rational::zero() || major.value() <= minor.value() {
                    return Err(refusal("path ellipse requires 0 < minor < major"));
                }
                let ellipse = Ellipse::new(
                    point(center),
                    vector(normal),
                    vector(major_direction),
                    real(major).square(),
                    real(minor).square(),
                    dec,
                )?;
                let curve = built.model.geom_mut().add_curve(Curve::Ellipse(ellipse));
                CentrePath::Arc {
                    curve,
                    from: None,
                    to: None,
                }
            }
            Self::Intersection {
                surfaces,
                component,
                from,
                to,
            } => {
                let a = surfaces[0].build(dec)?;
                let b = surfaces[1].build(dec)?;
                let a = built.model.geom_mut().add_surface(a);
                let b = built.model.geom_mut().add_surface(b);
                let intersection =
                    catallax_intersect::intersect_surfaces(built.model.geom_mut(), a, b, dec)?;
                if !intersection.marks.is_empty() {
                    return Err(refusal("centre-path intersection was not exact"));
                }
                let catallax_intersect::SurfacePair::Meets(components) = intersection.value else {
                    return Err(refusal("centre-path supports do not meet in curves"));
                };
                let curve = components
                    .into_iter()
                    .filter_map(|c| match c {
                        catallax_intersect::Component::Curve { curve, .. } => Some(curve),
                        _ => None,
                    })
                    .nth(*component)
                    .ok_or_else(|| refusal("centre-path component does not exist"))?;
                if from.is_none() && to.is_none() {
                    match built.model.geom().curve(curve) {
                        Curve::Line(line) => return Ok(CentrePath::Line(line.clone())),
                        Curve::Circle(circle) => return Ok(CentrePath::Circle(circle.clone())),
                        _ => {}
                    }
                }
                let mut endpoint = |value: &Option<[Number; 3]>| -> Result<_, Refusal> {
                    value
                        .as_ref()
                        .map(|value| {
                            let id = built
                                .model
                                .geom_mut()
                                .add_point(catallax_shapes::PointDef::Explicit(point(value)));
                            match built.model.geom().point_on_curve(curve, id, dec) {
                                catallax_numbers::Outcome::ProvenCoincident(proof)
                                    if proof.proof
                                        != catallax_numbers::ZeroProof::SeparationBound =>
                                {
                                    Ok(id)
                                }
                                catallax_numbers::Outcome::Refused(error) => Err(error),
                                _ => Err(refusal(
                                    "centre-path endpoint is not proven on the selected curve",
                                )),
                            }
                        })
                        .transpose()
                };
                CentrePath::Arc {
                    curve,
                    from: endpoint(from)?,
                    to: endpoint(to)?,
                }
            }
            Self::Edge { between } => {
                let faces = built.model.faces_of_body(built.result);
                let mut matches = Vec::new();
                for edge in built.model.edges_of_body(built.result) {
                    let mut found = [false; 2];
                    for face in &faces {
                        let incident = built.model.face(*face).loops.iter().any(|lp| {
                            match built.model.loop_(*lp) {
                                catallax_topology::Loop::Coedges(cs) => {
                                    cs.iter().any(|c| built.model.coedge(*c).edge == edge)
                                }
                                _ => false,
                            }
                        });
                        if incident {
                            for i in 0..2 {
                                found[i] |= built
                                    .model
                                    .face_origins(*face)
                                    .iter()
                                    .any(|o| between[i].matches(o));
                            }
                        }
                    }
                    if found == [true, true] {
                        matches.push(edge);
                    }
                }
                let [edge] = matches.as_slice() else {
                    return Err(refusal(
                        "centre path edge reference must resolve to exactly one edge",
                    ));
                };
                let edge = built.model.edge(*edge);
                CentrePath::Arc {
                    curve: edge.curve,
                    from: edge.start.map(|v| built.model.vertex(v).point),
                    to: edge.end.map(|v| built.model.vertex(v).point),
                }
            }
        })
    }
}
impl Definition {
    fn build(
        &self,
        built: &mut catallax_features::Built,
        dec: &mut Decider,
    ) -> Result<FaceDefinition, Refusal> {
        Ok(match self {
            Self::Fillet { path, radius } => FaceDefinition::Fillet {
                path: path.build(built, dec)?,
                radius: real(radius),
            },
            Self::Offset { original, distance } => FaceDefinition::Offset {
                original: original.build(dec)?,
                distance: real(distance),
            },
        })
    }
}

impl FaceChecks {
    pub(crate) fn parse(problem: &Value) -> Result<Self, Refusal> {
        let checks: Self =
            serde_json::from_value(problem.clone()).map_err(|e| refusal(e.to_string()))?;
        if checks.format != "catallax-face-check/1" || checks.definitions.is_empty() {
            return Err(refusal(
                "expected catallax-face-check/1 and at least one face definition",
            ));
        }
        Ok(checks)
    }

    pub(crate) fn report(
        &self,
        recipe: &Recipe,
        built: &mut catallax_features::Built,
        started: Instant,
        dec: &mut Decider,
    ) -> Result<Value, Refusal> {
        // Mesh first. This is the first complete display mesh from this actual
        // B-rep build, not the time to answer or the definition proof time.
        let mesh_options = catallax_mesh::MeshOptions::default();
        let mesh = catallax_mesh::mesh_body(&built.model, built.result, &mesh_options, dec)?;
        let first_mesh_ms = started.elapsed().as_secs_f64() * 1000.0;
        let options = MeasurementOptions {
            width: self.options.width.value().clone(),
            max_cells: self.options.max_cells,
            path_segments: self.options.path_segments,
        };
        let faces = built.model.faces_of_body(built.result);
        let mut seen = BTreeSet::new();
        let mut measurements = Vec::new();
        for attachment in &self.definitions {
            let selected: Vec<_> = faces
                .iter()
                .copied()
                .filter(|face| {
                    built
                        .model
                        .face_origins(*face)
                        .iter()
                        .any(|o| attachment.face.matches(o))
                })
                .collect();
            if selected.is_empty() {
                return Err(refusal(
                    "defined face reference resolved to no result faces",
                ));
            }
            let definition = attachment.definition.build(built, dec)?;
            for face in selected {
                if !seen.insert(face) {
                    return Err(refusal("a result face has more than one definition"));
                }
                let span = catallax_numbers::OperationSpan::begin(
                    "field",
                    &recipe.result,
                    "face_measure",
                    dec,
                );
                let measured =
                    measure_face(&mut built.model, face, Some(&definition), &options, dec);
                span.finish(dec, measured.as_ref().err(), None);
                let measured = measured?;
                measurements.push(json!({"face":face.index(),"origins":built.model.face_origins(face),
                    "kind":match definition {FaceDefinition::Fillet {..}=>"fillet",FaceDefinition::Offset {..}=>"offset"},
                    "deviation_mm":{"lower":measured.lower.to_string(),"upper":measured.upper.to_string()},
                    "local_tolerance_mm":measured.upper.to_string(),"converged":measured.converged,
                    "coverage":"whole-face supporting-chart cover; exact trimmed-face lower witnesses",
                    "cells":measured.cells,"witnesses":measured.witnesses}));
            }
        }
        let edit = self
            .edit
            .as_ref()
            .map(|edit| self.edit_report(recipe, edit, dec))
            .transpose()?;
        Ok(
            json!({"format":"catallax-face-measurements/1","faces":measurements,
            "display":{"first_mesh_ms":first_mesh_ms,"vertices":mesh.value.positions.len(),"triangles":mesh.value.triangles.len(),"chord_tolerance_mm":mesh_options.chord_tolerance.to_string()},
            "edit":edit}),
        )
    }

    fn edit_report(
        &self,
        recipe: &Recipe,
        edit: &Edit,
        dec: &mut Decider,
    ) -> Result<Value, Refusal> {
        // Edit the declared scalar in a copy, then parse and rebuild normally.
        // No internal parameter setter or reused face handles bypasses naming.
        let mut value: Value =
            serde_json::from_str(&recipe.to_text()).map_err(|e| refusal(e.to_string()))?;
        let steps = value["steps"].as_array_mut().unwrap();
        let step = steps
            .iter_mut()
            .find(|s| s["id"] == edit.step)
            .ok_or_else(|| refusal("edited step does not exist"))?;
        let parameter = step
            .get_mut(&edit.parameter)
            .filter(|v| v.is_string())
            .ok_or_else(|| refusal("edit must name an existing scalar parameter"))?;
        *parameter = json!(edit.value);
        let recipe = Recipe::parse(&value.to_string()).map_err(|e| refusal(e.to_string()))?;
        let edited =
            catallax_features::evaluate(&recipe, dec).map_err(|e| refusal(e.to_string()))?;
        let faces = edited.model.faces_of_body(edited.result);
        let resolved: Vec<_> = self.definitions.iter().map(|a| {
            let matches: Vec<_> = faces.iter().filter(|f| edited.model.face_origins(**f).iter().any(|o| a.face.matches(o))).map(|f| json!({"face":f.index(),"origins":edited.model.face_origins(*f)})).collect();
            json!({"reference":a.face,"status":if matches.is_empty() {"missing"} else {"resolved"},"faces":matches})
        }).collect();
        Ok(
            json!({"step":edit.step,"parameter":edit.parameter,"value":edit.value,"resolved":resolved}),
        )
    }
}

//! RWStepGeom field layouts, including complex instance components.
use super::*;
use crate::part21::Param;
use catallax_numbers::Rational;

fn arg(a: &[Param], i: usize, id: EntityId) -> Result<&Param, ImportRefusal> {
    a.get(i)
        .ok_or_else(|| malformed(id, format!("missing argument {i}")))
}
fn list(a: &[Param], i: usize, id: EntityId) -> Result<&[Param], ImportRefusal> {
    match arg(a, i, id)? {
        Param::List(v) => Ok(v),
        _ => Err(malformed(id, format!("argument {i} must be a list"))),
    }
}
pub(super) fn reference(a: &[Param], i: usize, id: EntityId) -> Result<EntityId, ImportRefusal> {
    match arg(a, i, id)? {
        Param::Ref(r) => Ok(*r),
        _ => Err(malformed(
            id,
            format!("argument {i} must reference an entity"),
        )),
    }
}
fn optional_ref(a: &[Param], i: usize, id: EntityId) -> Result<Option<EntityId>, ImportRefusal> {
    if matches!(arg(a, i, id)?, Param::Unset) {
        Ok(None)
    } else {
        reference(a, i, id).map(Some)
    }
}
pub(super) fn references(
    a: &[Param],
    i: usize,
    id: EntityId,
) -> Result<Vec<EntityId>, ImportRefusal> {
    list(a, i, id)?
        .iter()
        .map(|p| match p {
            Param::Ref(r) => Ok(*r),
            _ => Err(malformed(id, "expected entity references")),
        })
        .collect()
}
fn number(p: &Param, id: EntityId) -> Result<Rational, ImportRefusal> {
    let Param::Number(s) = p else {
        return Err(malformed(id, "expected a number"));
    };
    let mut s = s.replace(".E", ".0E").replace(".e", ".0e");
    if s.ends_with('.') {
        s.push('0');
    }
    Rational::parse(&s).map_err(|_| malformed(id, "unreadable source number"))
}
fn real(a: &[Param], i: usize, id: EntityId) -> Result<Rational, ImportRefusal> {
    number(arg(a, i, id)?, id)
}
fn optional_real(a: &[Param], i: usize, id: EntityId) -> Result<Option<Rational>, ImportRefusal> {
    if matches!(arg(a, i, id)?, Param::Unset) {
        Ok(None)
    } else {
        real(a, i, id).map(Some)
    }
}
pub(super) fn integer(a: &[Param], i: usize, id: EntityId) -> Result<usize, ImportRefusal> {
    let Param::Number(s) = arg(a, i, id)? else {
        return Err(malformed(id, "expected a nonnegative integer"));
    };
    s.parse()
        .map_err(|_| malformed(id, "expected a nonnegative integer"))
}
fn enumeration(a: &[Param], i: usize, id: EntityId) -> Result<String, ImportRefusal> {
    match arg(a, i, id)? {
        Param::Enum(s) => Ok(s.clone()),
        _ => Err(malformed(id, "expected an enumeration")),
    }
}
fn logical(a: &[Param], i: usize, id: EntityId) -> Result<Option<bool>, ImportRefusal> {
    match enumeration(a, i, id)?.as_str() {
        "T" => Ok(Some(true)),
        "F" => Ok(Some(false)),
        "U" => Ok(None),
        _ => Err(malformed(id, "expected a STEP logical")),
    }
}
fn boolean(a: &[Param], i: usize, id: EntityId) -> Result<bool, ImportRefusal> {
    logical(a, i, id)?.ok_or_else(|| malformed(id, "expected a STEP boolean"))
}
fn reals(a: &[Param], i: usize, id: EntityId) -> Result<Vec<Rational>, ImportRefusal> {
    list(a, i, id)?.iter().map(|p| number(p, id)).collect()
}
fn reference_grid(
    a: &[Param],
    i: usize,
    id: EntityId,
) -> Result<Vec<Vec<EntityId>>, ImportRefusal> {
    list(a, i, id)?
        .iter()
        .map(|p| references(std::slice::from_ref(p), 0, id))
        .collect()
}
fn trim(a: &[Param], i: usize, id: EntityId) -> Result<Vec<TrimSelect>, ImportRefusal> {
    list(a, i, id)?
        .iter()
        .map(|p| match p {
            Param::Ref(r) => Ok(TrimSelect::Point(*r)),
            Param::Typed(kind, values) if kind == "PARAMETER_VALUE" && values.len() == 1 => {
                number(&values[0], id).map(TrimSelect::Parameter)
            }
            _ => Err(malformed(
                id,
                "expected a Cartesian point or PARAMETER_VALUE trim",
            )),
        })
        .collect()
}
fn component<'a>(e: &'a Entity, kind: &str, id: EntityId) -> Result<&'a [Param], ImportRefusal> {
    e.part(kind)
        .ok_or_else(|| malformed(id, format!("missing {kind} component")))
}
fn placement(
    e: &Entity,
    a: &[Param],
    parent: &str,
    id: EntityId,
) -> Result<EntityId, ImportRefusal> {
    if e.parts.len() == 1 {
        reference(a, 1, id)
    } else {
        reference(component(e, parent, id)?, 0, id)
    }
}

pub(super) fn geometry(
    source: &ParsedStep,
    id: EntityId,
    dec: &mut Decider,
) -> Result<Geometry, ImportRefusal> {
    let e = source.entity(id)?;
    let mut pending: Vec<&Param> = e.parts.iter().flat_map(|(_, a)| a).collect();
    while let Some(p) = pending.pop() {
        let work = match p {
            Param::Number(s) | Param::Str(s) | Param::Enum(s) => s.len() as u64 + 1,
            _ => 1,
        };
        dec.charge_work(work)
            .map_err(|_| ImportRefusal::ResourceLimit)?;
        match p {
            Param::List(v) | Param::Typed(_, v) => pending.extend(v),
            _ => {}
        }
    }
    if let Some(spline) = spline(e, id)? {
        return Ok(spline);
    }
    let kinds = [
        "CARTESIAN_POINT",
        "DIRECTION",
        "VECTOR",
        "AXIS1_PLACEMENT",
        "AXIS2_PLACEMENT_2D",
        "AXIS2_PLACEMENT_3D",
        "CARTESIAN_TRANSFORMATION_OPERATOR_2D",
        "CARTESIAN_TRANSFORMATION_OPERATOR_3D",
        "CARTESIAN_TRANSFORMATION_OPERATOR",
        "LINE",
        "CIRCLE",
        "ELLIPSE",
        "HYPERBOLA",
        "PARABOLA",
        "POLYLINE",
        "TRIMMED_CURVE",
        "CURVE_REPLICA",
        "OFFSET_CURVE_3D",
        "SEAM_CURVE",
        "INTERSECTION_CURVE",
        "SURFACE_CURVE",
        "PCURVE",
        "OUTER_BOUNDARY_CURVE",
        "BOUNDARY_CURVE",
        "COMPOSITE_CURVE",
        "COMPOSITE_CURVE_SEGMENT",
        "PLANE",
        "CYLINDRICAL_SURFACE",
        "CONICAL_SURFACE",
        "SPHERICAL_SURFACE",
        "TOROIDAL_SURFACE",
        "SURFACE_OF_LINEAR_EXTRUSION",
        "SURFACE_OF_REVOLUTION",
        "RECTANGULAR_TRIMMED_SURFACE",
        "SURFACE_REPLICA",
        "OFFSET_SURFACE",
        "CURVE_BOUNDED_SURFACE",
        "RECTANGULAR_COMPOSITE_SURFACE",
        "SURFACE_PATCH",
        "DEFINITIONAL_REPRESENTATION",
    ];
    let Some((kind, a)) = kinds
        .iter()
        .find_map(|kind| e.part(kind).map(|a| (*kind, a)))
    else {
        return Err(ImportRefusal::Unsupported {
            entity: format!("#{id}"),
            type_name: e.type_name().into(),
        });
    };
    declaration(source, e, kind, a, id)
}

fn declaration(
    source: &ParsedStep,
    e: &Entity,
    kind: &str,
    a: &[Param],
    id: EntityId,
) -> Result<Geometry, ImportRefusal> {
    let named = usize::from(e.parts.len() == 1);
    let conic = || placement(e, a, "CONIC", id);
    let elementary = || placement(e, a, "ELEMENTARY_SURFACE", id);
    let sweep = || placement(e, a, "SWEPT_SURFACE", id);
    Ok(match kind {
        "CARTESIAN_POINT" | "DIRECTION" => {
            let values = reals(a, named, id)?;
            if !(usize::from(kind == "DIRECTION") + 1..=3).contains(&values.len()) {
                return Err(malformed(
                    id,
                    "invalid Cartesian coordinate or direction dimension",
                ));
            }
            if kind == "DIRECTION" {
                Geometry::Direction(values)
            } else {
                Geometry::Point(values)
            }
        }
        "VECTOR" => Geometry::Vector(reference(a, named, id)?, real(a, named + 1, id)?),
        "AXIS1_PLACEMENT" => Geometry::Axis1Placement(
            placement(e, a, "PLACEMENT", id)?,
            optional_ref(a, if named == 1 { 2 } else { 0 }, id)?,
        ),
        "AXIS2_PLACEMENT_2D" => Geometry::Axis2Placement2d(
            placement(e, a, "PLACEMENT", id)?,
            optional_ref(a, if named == 1 { 2 } else { 0 }, id)?,
        ),
        "AXIS2_PLACEMENT_3D" => {
            let offset = if named == 1 { 2 } else { 0 };
            Geometry::Axis2Placement3d(
                placement(e, a, "PLACEMENT", id)?,
                optional_ref(a, offset, id)?,
                optional_ref(a, offset + 1, id)?,
            )
        }
        "CARTESIAN_TRANSFORMATION_OPERATOR"
        | "CARTESIAN_TRANSFORMATION_OPERATOR_2D"
        | "CARTESIAN_TRANSFORMATION_OPERATOR_3D" => {
            Geometry::Transformation(transformation(e, kind, a, id)?)
        }
        "LINE" => Geometry::Curve(Curve::Line(
            reference(a, named, id)?,
            reference(a, named + 1, id)?,
        )),
        "CIRCLE" => Geometry::Curve(Curve::Circle(
            conic()?,
            real(a, if named == 1 { 2 } else { 0 }, id)?,
        )),
        "ELLIPSE" | "HYPERBOLA" => {
            let offset = if named == 1 { 2 } else { 0 };
            let values = (conic()?, real(a, offset, id)?, real(a, offset + 1, id)?);
            Geometry::Curve(if kind == "ELLIPSE" {
                Curve::Ellipse(values.0, values.1, values.2)
            } else {
                Curve::Hyperbola(values.0, values.1, values.2)
            })
        }
        "PARABOLA" => Geometry::Curve(Curve::Parabola(
            conic()?,
            real(a, if named == 1 { 2 } else { 0 }, id)?,
        )),
        "POLYLINE" => Geometry::Curve(Curve::Polyline(references(a, named, id)?)),
        "TRIMMED_CURVE" => {
            let preference = match enumeration(a, named + 4, id)?.as_str() {
                "CARTESIAN" => TrimmingPreference::Cartesian,
                "PARAMETER" => TrimmingPreference::Parameter,
                "UNSPECIFIED" => TrimmingPreference::Unspecified,
                _ => return Err(malformed(id, "unknown trimming preference")),
            };
            Geometry::Curve(Curve::Trimmed(
                reference(a, named, id)?,
                [trim(a, named + 1, id)?, trim(a, named + 2, id)?],
                boolean(a, named + 3, id)?,
                preference,
            ))
        }
        "CURVE_REPLICA" => Geometry::Curve(Curve::Replica(
            reference(a, named, id)?,
            reference(a, named + 1, id)?,
        )),
        "OFFSET_CURVE_3D" => Geometry::Curve(Curve::Offset3d(
            reference(a, named, id)?,
            real(a, named + 1, id)?,
            logical(a, named + 2, id)?,
            reference(a, named + 3, id)?,
        )),
        "SURFACE_CURVE" | "SEAM_CURVE" | "INTERSECTION_CURVE" => {
            let a = if named == 0 {
                component(e, "SURFACE_CURVE", id)?
            } else {
                a
            };
            Geometry::Curve(Curve::SurfaceCurve(
                reference(a, named, id)?,
                references(a, named + 1, id)?,
                enumeration(a, named + 2, id)?,
                kind == "SEAM_CURVE" || e.part("SEAM_CURVE").is_some(),
            ))
        }
        "PCURVE" => Geometry::Curve(Curve::Pcurve(
            reference(a, named, id)?,
            reference(a, named + 1, id)?,
        )),
        "COMPOSITE_CURVE" | "BOUNDARY_CURVE" | "OUTER_BOUNDARY_CURVE" => {
            let a = if named == 0 {
                component(e, "COMPOSITE_CURVE", id)?
            } else {
                a
            };
            Geometry::Curve(Curve::Composite(
                references(a, named, id)?,
                logical(a, named + 1, id)?,
            ))
        }
        "COMPOSITE_CURVE_SEGMENT" => Geometry::CompositeSegment(
            enumeration(a, 0, id)?,
            boolean(a, 1, id)?,
            reference(a, 2, id)?,
        ),
        "PLANE" => Geometry::Surface(Surface::Plane(elementary()?)),
        "CYLINDRICAL_SURFACE" | "SPHERICAL_SURFACE" => {
            let values = (elementary()?, real(a, if named == 1 { 2 } else { 0 }, id)?);
            Geometry::Surface(if kind == "CYLINDRICAL_SURFACE" {
                Surface::Cylinder(values.0, values.1)
            } else {
                Surface::Sphere(values.0, values.1)
            })
        }
        "CONICAL_SURFACE" | "TOROIDAL_SURFACE" => {
            let offset = if named == 1 { 2 } else { 0 };
            let values = (
                elementary()?,
                real(a, offset, id)?,
                real(a, offset + 1, id)?,
            );
            Geometry::Surface(if kind == "CONICAL_SURFACE" {
                Surface::Cone(values.0, values.1, values.2)
            } else {
                Surface::Torus(values.0, values.1, values.2)
            })
        }
        "SURFACE_OF_LINEAR_EXTRUSION" | "SURFACE_OF_REVOLUTION" => {
            let values = (sweep()?, reference(a, if named == 1 { 2 } else { 0 }, id)?);
            Geometry::Surface(if kind == "SURFACE_OF_REVOLUTION" {
                Surface::Revolution(values.0, values.1)
            } else {
                Surface::Extrusion(values.0, values.1)
            })
        }
        "RECTANGULAR_TRIMMED_SURFACE" => Geometry::Surface(Surface::RectangularTrim(
            reference(a, named, id)?,
            [
                real(a, named + 1, id)?,
                real(a, named + 2, id)?,
                real(a, named + 3, id)?,
                real(a, named + 4, id)?,
            ],
            [boolean(a, named + 5, id)?, boolean(a, named + 6, id)?],
        )),
        "SURFACE_REPLICA" => Geometry::Surface(Surface::Replica(
            reference(a, named, id)?,
            reference(a, named + 1, id)?,
        )),
        "OFFSET_SURFACE" => Geometry::Surface(Surface::Offset(
            reference(a, named, id)?,
            real(a, named + 1, id)?,
            logical(a, named + 2, id)?,
        )),
        "CURVE_BOUNDED_SURFACE" => Geometry::Surface(Surface::CurveBounded(
            reference(a, named, id)?,
            references(a, named + 1, id)?,
            boolean(a, named + 2, id)?,
        )),
        "RECTANGULAR_COMPOSITE_SURFACE" => {
            Geometry::Surface(Surface::RectangularComposite(reference_grid(a, named, id)?))
        }
        "SURFACE_PATCH" => Geometry::SurfacePatch(
            reference(a, 0, id)?,
            [enumeration(a, 1, id)?, enumeration(a, 2, id)?],
            [boolean(a, 3, id)?, boolean(a, 4, id)?],
        ),
        "DEFINITIONAL_REPRESENTATION" => {
            let r = source.representation(id)?;
            Geometry::DefinitionalRepresentation(r.items, r.context)
        }
        _ => {
            return Err(ImportRefusal::Unsupported {
                entity: format!("#{id}"),
                type_name: kind.into(),
            });
        }
    })
}

fn transformation(
    e: &Entity,
    kind: &str,
    a: &[Param],
    id: EntityId,
) -> Result<Transformation, ImportRefusal> {
    let dimension = if kind.ends_with("_3D") { 3 } else { 2 };
    if e.parts.len() == 1 {
        let modern = if dimension == 3 { 6 } else { 5 };
        let offset = match a.len() {
            n if n == modern => 1,
            n if n == modern + 2 => 3,
            _ => {
                return Err(malformed(
                    id,
                    "unrecognized Cartesian transformation field layout",
                ));
            }
        };
        Ok(Transformation {
            dimension: (kind != "CARTESIAN_TRANSFORMATION_OPERATOR").then_some(dimension),
            origin: reference(a, offset + 2, id)?,
            axes: [
                optional_ref(a, offset, id)?,
                optional_ref(a, offset + 1, id)?,
                if dimension == 3 {
                    optional_ref(a, offset + 4, id)?
                } else {
                    None
                },
            ],
            scales: [optional_real(a, offset + 3, id)?, None, None],
        })
    } else {
        let base = component(e, "CARTESIAN_TRANSFORMATION_OPERATOR", id)?;
        Ok(Transformation {
            dimension: (kind != "CARTESIAN_TRANSFORMATION_OPERATOR").then_some(dimension),
            origin: reference(base, 2, id)?,
            axes: [
                optional_ref(base, 0, id)?,
                optional_ref(base, 1, id)?,
                if dimension == 3 {
                    optional_ref(a, 0, id)?
                } else {
                    None
                },
            ],
            scales: [optional_real(base, 3, id)?, None, None],
        })
    }
}

fn integers(a: &[Param], i: usize, id: EntityId) -> Result<Vec<usize>, ImportRefusal> {
    list(a, i, id)?
        .iter()
        .map(|p| integer(std::slice::from_ref(p), 0, id))
        .collect()
}
fn real_grid(a: &[Param], i: usize, id: EntityId) -> Result<Vec<Vec<Rational>>, ImportRefusal> {
    list(a, i, id)?
        .iter()
        .map(|p| reals(std::slice::from_ref(p), 0, id))
        .collect()
}
fn implicit_knots(e: &Entity, curve: bool, id: EntityId) -> Result<KnotData, ImportRefusal> {
    let suffix = if curve { "CURVE" } else { "SURFACE" };
    for (prefix, value) in [
        ("BEZIER", KnotData::Bezier),
        ("UNIFORM", KnotData::Uniform),
        ("QUASI_UNIFORM", KnotData::QuasiUniform),
    ] {
        if e.part(&format!("{prefix}_{suffix}")).is_some() {
            return Ok(value);
        }
    }
    Err(malformed(id, "B-spline has no concrete knot encoding"))
}
fn spline(e: &Entity, id: EntityId) -> Result<Option<Geometry>, ImportRefusal> {
    let simple = e.parts.len() == 1;
    let curve_names = [
        "B_SPLINE_CURVE_WITH_KNOTS",
        "BEZIER_CURVE",
        "UNIFORM_CURVE",
        "QUASI_UNIFORM_CURVE",
    ];
    let surface_names = [
        "B_SPLINE_SURFACE_WITH_KNOTS",
        "BEZIER_SURFACE",
        "UNIFORM_SURFACE",
        "QUASI_UNIFORM_SURFACE",
    ];
    let curve_base = e.part("B_SPLINE_CURVE").or_else(|| {
        if simple {
            curve_names.iter().find_map(|name| e.part(name))
        } else {
            None
        }
    });
    if let Some(a) = curve_base {
        let n = usize::from(simple);
        let knots = if let Some(k) = e.part("B_SPLINE_CURVE_WITH_KNOTS") {
            let offset = if simple { 6 } else { 0 };
            KnotData::Explicit(
                integers(k, offset, id)?,
                reals(k, offset + 1, id)?,
                enumeration(k, offset + 2, id)?,
            )
        } else {
            implicit_knots(e, true, id)?
        };
        let weights = e
            .part("RATIONAL_B_SPLINE_CURVE")
            .map(|a| reals(a, 0, id))
            .transpose()?;
        return Ok(Some(Geometry::Curve(Curve::BSpline(Box::new(
            SplineCurve {
                degree: integer(a, n, id)?,
                controls: references(a, n + 1, id)?,
                form: enumeration(a, n + 2, id)?,
                closed: logical(a, n + 3, id)?,
                self_intersect: logical(a, n + 4, id)?,
                weights,
                knots,
            },
        )))));
    }
    let surface_base = e.part("B_SPLINE_SURFACE").or_else(|| {
        if simple {
            surface_names.iter().find_map(|name| e.part(name))
        } else {
            None
        }
    });
    if let Some(a) = surface_base {
        let n = usize::from(simple);
        let knots = if let Some(k) = e.part("B_SPLINE_SURFACE_WITH_KNOTS") {
            let offset = if simple { 8 } else { 0 };
            let spec = enumeration(k, offset + 4, id)?;
            [
                KnotData::Explicit(
                    integers(k, offset, id)?,
                    reals(k, offset + 2, id)?,
                    spec.clone(),
                ),
                KnotData::Explicit(
                    integers(k, offset + 1, id)?,
                    reals(k, offset + 3, id)?,
                    spec,
                ),
            ]
        } else {
            let knots = implicit_knots(e, false, id)?;
            [knots.clone(), knots]
        };
        let weights = e
            .part("RATIONAL_B_SPLINE_SURFACE")
            .map(|a| real_grid(a, 0, id))
            .transpose()?;
        return Ok(Some(Geometry::Surface(Surface::BSpline(Box::new(
            SplineSurface {
                degree: [integer(a, n, id)?, integer(a, n + 1, id)?],
                controls: reference_grid(a, n + 2, id)?,
                form: enumeration(a, n + 3, id)?,
                closed: [logical(a, n + 4, id)?, logical(a, n + 5, id)?],
                self_intersect: logical(a, n + 6, id)?,
                weights,
                knots,
            },
        )))));
    }
    Ok(None)
}

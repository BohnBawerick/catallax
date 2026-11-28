//! Checks of datum targets and position controls against the part's exact
//! geometry, made before they are written and available to a caller that maps
//! imported items onto an imported model.

use std::collections::{BTreeMap, BTreeSet};

use catallax_intersect::{CurveForm, curve_form};
use catallax_numbers::{Decider, Real, Refusal};
use catallax_shapes::distance::{PlanarBoundary, PlanarEdge};
use catallax_shapes::{Circle, Plane, Point3, Surface, Vector3};
use catallax_topology::{BodyId, BodyKind, FaceId, Loop, Model, Sense};

use super::measure::{
    Shape, face_edges, feature_shape, is_positive, is_zero, is_zero_vector, oriented_plane,
};
use super::target::{self, DatumTarget};
use super::{Feature, ItemKind, Pmi, PositionControl, PositionFrame};
use crate::graph::{Field, Graph, Id};

fn refuse(item: &str, what: impl std::fmt::Display) -> Refusal {
    Refusal::new("step.pmi-invalid", format!("{item}: {what}"))
}

/// The complete trim of a planar face from its loops, or `None` when the face
/// is not planar or an edge is neither a line nor a circle.
pub(super) fn face_boundary(
    model: &Model,
    face: FaceId,
    dec: &mut Decider,
) -> Result<Option<PlanarBoundary>, Refusal> {
    let Surface::Plane(plane) = model.geom().surface(model.face(face).surface) else {
        return Ok(None);
    };
    let mut vertices = Vec::new();
    let mut indices = BTreeMap::new();
    let mut edges = Vec::new();
    for boundary in &model.face(face).loops {
        let Loop::Coedges(coedges) = model.loop_(*boundary) else {
            continue;
        };
        for coedge in coedges {
            let edge = model.edge(model.coedge(*coedge).edge);
            let ends = match (edge.start, edge.end) {
                (None, None) => None,
                (Some(a), Some(b)) => {
                    let mut pair = [a, b].map(|vertex| {
                        *indices.entry(vertex).or_insert_with(|| {
                            vertices.push(model.geom().position(model.vertex(vertex).point));
                            vertices.len() - 1
                        })
                    });
                    if edge.sense == Sense::Reversed {
                        pair.swap(0, 1);
                    }
                    Some(pair)
                }
                _ => return Ok(None),
            };
            match curve_form(model.geom(), edge.curve, dec)? {
                CurveForm::Line { .. } => {
                    let Some(pair) = ends else {
                        return Ok(None);
                    };
                    edges.push(PlanarEdge::Segment(pair));
                }
                CurveForm::Circle {
                    center,
                    normal,
                    radius_sq,
                    ref_dir,
                } => edges.push(PlanarEdge::Arc {
                    circle: Circle {
                        center,
                        normal,
                        radius_sq,
                        ref_dir,
                    },
                    // A circle closed at one vertex is the whole circle.
                    ends: ends.filter(|[a, b]| a != b),
                }),
                _ => return Ok(None),
            }
        }
    }
    Ok((!edges.is_empty()).then(|| PlanarBoundary {
        plane: plane.clone(),
        vertices,
        edges,
    }))
}

/// A planar face and its normal pointing out of the material.
pub(super) fn outward_plane(model: &Model, face: FaceId) -> Option<(&Plane, Vector3)> {
    let face = model.face(face);
    let Surface::Plane(plane) = model.geom().surface(face.surface) else {
        return None;
    };
    let normal = if face.sense == Sense::Same {
        plane.normal.clone()
    } else {
        plane.normal.scaled(&Real::from_i64(-1))
    };
    Some((plane, normal))
}

fn check_target(
    model: &Model,
    target: &DatumTarget<FaceId>,
    item: &str,
    dec: &mut Decider,
) -> Result<(), Refusal> {
    let boundary = face_boundary(model, target.face, dec)?.ok_or_else(|| {
        refuse(
            item,
            "a datum target needs a planar face bounded by lines and circles",
        )
    })?;
    match target::off_face(&boundary, target, dec)? {
        Some(reason) => Err(refuse(item, reason)),
        None => Ok(()),
    }
}

fn fields<'a>(graph: &'a Graph, field: &Field, name: &str) -> Option<&'a [Field]> {
    let Field::Ref(id) = field else { return None };
    let (kind, fields) = graph.record(*id);
    (kind == name).then_some(fields)
}

fn written_real(value: &str) -> Real {
    let value = value.replace(".E", ".0E");
    let value = if value.ends_with('.') {
        format!("{value}0")
    } else {
        value
    };
    Real::from_rational(crate::number::formatted_value(&value))
}

fn values(fields: &[Field]) -> Option<[Real; 3]> {
    let values = match fields.get(1)? {
        Field::Text(values) => values
            .strip_prefix('(')?
            .strip_suffix(')')?
            .split(',')
            .collect::<Vec<_>>(),
        Field::List(values) => values
            .iter()
            .map(|v| match v {
                Field::Text(text) => Some(text.as_ref()),
                _ => None,
            })
            .collect::<Option<Vec<_>>>()?,
        _ => return None,
    };
    let values: [&str; 3] = values.try_into().ok()?;
    Some(values.map(written_real))
}

fn written_plane(graph: &Graph, placement: &Field) -> Option<Plane> {
    let frame = fields(graph, placement, "AXIS2_PLACEMENT_3D")?;
    let [x, y, z] = values(fields(graph, frame.get(1)?, "CARTESIAN_POINT")?)?;
    let origin = Point3::new(x, y, z);
    let [x, y, z] = values(fields(graph, frame.get(2)?, "DIRECTION")?)?;
    let normal = Vector3::new(x, y, z);
    let [x, y, z] = values(fields(graph, frame.get(3)?, "DIRECTION")?)?;
    let reference = Vector3::new(x, y, z);
    let ref_dir = normal.cross(&reference.cross(&normal));
    Some(Plane {
        origin,
        normal,
        ref_dir,
    })
}

fn written_face(graph: &Graph, face: Id) -> Option<(&[Field], &Field, bool)> {
    let ("ADVANCED_FACE", [_, Field::List(bounds), surface, Field::Text(sense)]) =
        graph.record(face)
    else {
        return None;
    };
    Some((bounds, surface, sense.as_ref() == ".T."))
}

fn written_edges(graph: &Graph, bounds: &[Field]) -> Option<BTreeSet<Id>> {
    let mut edges = BTreeSet::new();
    for bound in bounds {
        let cycle = fields(graph, bound, "FACE_BOUND")
            .or_else(|| fields(graph, bound, "FACE_OUTER_BOUND"))?
            .get(1)?;
        if fields(graph, cycle, "VERTEX_LOOP").is_some() {
            continue;
        }
        let Field::List(oriented) = fields(graph, cycle, "EDGE_LOOP")?.get(1)? else {
            return None;
        };
        for edge in oriented {
            let Field::Ref(edge) = fields(graph, edge, "ORIENTED_EDGE")?.get(3)? else {
                return None;
            };
            edges.insert(*edge);
        }
    }
    Some(edges)
}

fn written_point(graph: &Graph, vertex: &Field) -> Option<Point3> {
    let point = fields(graph, vertex, "VERTEX_POINT")?.get(1)?;
    let [x, y, z] = values(fields(graph, point, "CARTESIAN_POINT")?)?;
    Some(Point3::new(x, y, z))
}

fn written_curve<'a>(graph: &'a Graph, curve: &'a Field) -> &'a Field {
    fields(graph, curve, "SEAM_CURVE")
        .or_else(|| fields(graph, curve, "SURFACE_CURVE"))
        .and_then(|fields| fields.get(1))
        .unwrap_or(curve)
}

fn written_boundary(graph: &Graph, face: Id) -> Option<PlanarBoundary> {
    let (bounds, surface, _) = written_face(graph, face)?;
    let plane = written_plane(graph, fields(graph, surface, "PLANE")?.get(1)?)?;
    let mut vertices = Vec::new();
    let mut edges = Vec::new();
    for id in written_edges(graph, bounds)? {
        let (
            "EDGE_CURVE",
            [
                _,
                Field::Ref(start),
                Field::Ref(end),
                curve,
                Field::Text(same),
            ],
        ) = graph.record(id)
        else {
            return None;
        };
        let ends = if start == end {
            None
        } else {
            let mut ends = [0, 0];
            for (slot, vertex) in ends.iter_mut().zip([*start, *end]) {
                *slot = vertices.len();
                vertices.push(written_point(graph, &Field::Ref(vertex))?);
            }
            if same.as_ref() == ".F." {
                ends.swap(0, 1);
            }
            Some(ends)
        };
        let curve = written_curve(graph, curve);
        if fields(graph, curve, "LINE").is_some() {
            edges.push(PlanarEdge::Segment(ends?));
        } else {
            let [_, placement, Field::Text(radius)] = fields(graph, curve, "CIRCLE")? else {
                return None;
            };
            let circle = written_plane(graph, placement)?;
            edges.push(PlanarEdge::Arc {
                circle: Circle {
                    center: circle.origin,
                    normal: circle.normal,
                    radius_sq: written_real(radius).square(),
                    ref_dir: circle.ref_dir,
                },
                ends,
            });
        }
    }
    (!edges.is_empty()).then_some(PlanarBoundary {
        plane,
        vertices,
        edges,
    })
}

pub(super) fn reversed_faces(graph: &Graph) -> BTreeSet<Id> {
    graph
        .records()
        .filter_map(|(kind, fields)| {
            let ("ORIENTED_CLOSED_SHELL", [_, _, shell, Field::Text(sense)]) = (kind, fields)
            else {
                return None;
            };
            if sense.as_ref() != ".F." {
                return None;
            }
            let Field::List(faces) = self::fields(graph, shell, "CLOSED_SHELL")?.get(1)? else {
                return None;
            };
            Some(faces)
        })
        .flatten()
        .filter_map(|face| match face {
            Field::Ref(id) => Some(*id),
            _ => None,
        })
        .collect()
}

pub(super) fn check_written_target(
    graph: &Graph,
    faces: &[Id],
    reversed: &BTreeSet<Id>,
    target: &DatumTarget<FaceId>,
    normal: &Vector3,
    item: &str,
    dec: &mut Decider,
) -> Result<(), Refusal> {
    let [face] = faces else {
        return Err(refuse(item, "the datum target's face was split on export"));
    };
    let no_trim = || refuse(item, "the datum target's written trim cannot be checked");
    let boundary = written_boundary(graph, *face).ok_or_else(no_trim)?;
    let (_, _, same) = written_face(graph, *face).ok_or_else(no_trim)?;
    let same = same != reversed.contains(face);
    let plane = &boundary.plane;
    if is_zero_vector(&plane.ref_dir, dec)?
        || !is_zero_vector(&normal.cross(&plane.normal), dec)?
        || is_positive(&normal.dot(&plane.normal), dec)? != same
    {
        return Err(refuse(
            item,
            "the datum target's direction differs from the written face",
        ));
    }
    for vertex in &boundary.vertices {
        if !is_zero(&vertex.minus(&plane.origin).dot(&plane.normal), dec)? {
            return Err(refuse(
                item,
                "rounding changed the datum target's written trim plane",
            ));
        }
    }
    for edge in &boundary.edges {
        if let PlanarEdge::Arc { circle, ends } = edge {
            if is_zero_vector(&circle.ref_dir, dec)?
                || !is_zero_vector(&circle.normal.cross(&plane.normal), dec)?
                || !is_zero(&circle.center.minus(&plane.origin).dot(&plane.normal), dec)?
            {
                return Err(refuse(
                    item,
                    "rounding changed the datum target's written circular trim",
                ));
            }
            if let Some(ends) = ends {
                for end in ends {
                    let d = boundary.vertices[*end].minus(&circle.center).norm_squared();
                    if !is_zero(&(&d - &circle.radius_sq), dec)? {
                        return Err(refuse(
                            item,
                            "rounding changed the datum target's written circular trim",
                        ));
                    }
                }
            }
        }
    }
    if let Some(reason) = target::off_face(&boundary, target, dec)? {
        return Err(refuse(item, format!("on the written face, {reason}")));
    }
    Ok(())
}

/// The stretch of a projected zone's axis the zone covers: from `at` on the
/// entry face, along `direction`, for the zone's height. A face blocks the
/// zone only where it meets this stretch; a face wholly behind the start or
/// wholly beyond the end leaves it free.
#[doc(hidden)]
pub struct ZoneSpan {
    at: Point3,
    direction: Vector3,
    /// The height in lengths of `direction`.
    reach: Real,
}

/// Where the axis past a zone's end meets one face of the part.
#[doc(hidden)]
#[derive(Clone, Debug)]
pub enum Crossing {
    /// It passes into the material, or out of it, at `t` lengths of the
    /// zone's direction from the start.
    At {
        /// The crossing's place along the axis.
        t: Real,
        /// Whether the axis leaves the material there.
        exiting: bool,
    },
    /// It may meet the face from `t` on, at a place this check cannot settle.
    From(Real),
}

/// Which end of a [`ZoneSpan`] a whole curve lies past.
#[doc(hidden)]
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum Side {
    /// At or behind the start, against the zone's direction.
    Behind,
    /// At or past the end of the zone's height.
    Beyond,
}

/// Where a projected zone's axis meets a plane.
#[doc(hidden)]
pub enum PlaneMeeting {
    /// The forward axis misses the plane.
    Miss,
    /// The axis lies in the plane.
    Coplanar,
    /// The axis meets the plane within the projected span.
    Within(Point3),
    /// The axis meets the plane at or beyond the span's end.
    Beyond(Real, Point3),
}

impl ZoneSpan {
    /// The span from `at` along `direction` for `height` in the same
    /// length unit as `at`.
    pub fn new(at: Point3, direction: Vector3, height: &Real) -> Self {
        let reach = height / &direction.dot(&direction).sqrt();
        Self {
            at,
            direction,
            reach,
        }
    }

    /// The zone's height in lengths of its direction.
    pub fn reach(&self) -> Real {
        self.reach.clone()
    }

    fn end(&self) -> Point3 {
        self.at.plus(&self.direction.scaled(&self.reach))
    }

    /// Where the forward axis meets a plane relative to the zone's height.
    pub fn plane_meeting(&self, plane: &Plane, dec: &mut Decider) -> Result<PlaneMeeting, Refusal> {
        let height = plane.origin.minus(&self.at).dot(&plane.normal);
        let along = self.direction.dot(&plane.normal);
        if is_zero(&along, dec)? {
            return Ok(if is_zero(&height, dec)? {
                PlaneMeeting::Coplanar
            } else {
                PlaneMeeting::Miss
            });
        }
        let t = height / along;
        if is_positive(&(-t.clone()), dec)? {
            return Ok(PlaneMeeting::Miss);
        }
        let point = self.at.plus(&self.direction.scaled(&t));
        Ok(if is_positive(&(&self.reach - &t), dec)? {
            PlaneMeeting::Within(point)
        } else {
            PlaneMeeting::Beyond(t, point)
        })
    }

    /// Whether the span, from its start up to but not at its end, stays off
    /// a cylinder's surface. A span that touches it counts as meeting it.
    pub fn misses_cylinder(
        &self,
        cylinder: &catallax_shapes::Cylinder,
        dec: &mut Decider,
    ) -> Result<bool, Refusal> {
        // |u + t v|^2 = r^2 |axis|^2 on the surface: a t^2 + 2 b t + c = 0.
        let u = self.at.minus(&cylinder.origin).cross(&cylinder.axis);
        let v = self.direction.cross(&cylinder.axis);
        let a = v.dot(&v);
        let b = u.dot(&v);
        let c = u.dot(&u) - &cylinder.radius_sq * &cylinder.axis.dot(&cylinder.axis);
        if is_zero(&a, dec)? {
            return Ok(!is_zero(&c, dec)?);
        }
        if is_positive(&(&a * &c - &b * &b), dec)? {
            return Ok(true);
        }
        if is_zero(&c, dec)? {
            return Ok(false);
        }
        let slope = &a * &self.reach + b.clone() + b.clone();
        let at_end = &slope * &self.reach + c.clone();
        if is_positive(&(-c.clone()), dec)? {
            // The start is inside: the span leaves at the larger root.
            return Ok(!is_positive(&at_end, dec)?);
        }
        if is_positive(&b, dec)? {
            // Both roots are behind the start.
            return Ok(true);
        }
        // Both roots are ahead: clear when the first is at or past the end.
        Ok(!is_positive(&(-at_end), dec)? && !is_positive(&(&a * &self.reach + b), dec)?)
    }

    /// Whether the span, start and end included, stays clear of a ball.
    pub fn misses_ball(
        &self,
        centre: &Point3,
        radius_sq: &Real,
        dec: &mut Decider,
    ) -> Result<bool, Refusal> {
        let offset = centre.minus(&self.at);
        let along = offset.dot(&self.direction);
        let length_sq = self.direction.dot(&self.direction);
        let nearest = if !is_positive(&along, dec)? {
            self.at.clone()
        } else if !is_positive(&(&self.reach * &length_sq - along.clone()), dec)? {
            self.end()
        } else {
            self.at.plus(&self.direction.scaled(&(along / length_sq)))
        };
        let gap = centre.minus(&nearest);
        is_positive(&(gap.dot(&gap) - radius_sq.clone()), dec)
    }

    /// The first place at or past the span's end where the axis may meet a
    /// cylinder's surface.
    pub fn cylinder_beyond(
        &self,
        cylinder: &catallax_shapes::Cylinder,
        dec: &mut Decider,
    ) -> Result<Option<Real>, Refusal> {
        let u = self.at.minus(&cylinder.origin).cross(&cylinder.axis);
        let v = self.direction.cross(&cylinder.axis);
        let a = v.dot(&v);
        let b = u.dot(&v);
        let c = u.dot(&u) - &cylinder.radius_sq * &cylinder.axis.dot(&cylinder.axis);
        if is_zero(&a, dec)? {
            return Ok(is_zero(&c, dec)?.then(|| self.reach.clone()));
        }
        let discriminant = &b * &b - &a * &c;
        if is_positive(&(-discriminant.clone()), dec)? {
            return Ok(None);
        }
        let root = discriminant.sqrt();
        let last = (-b.clone() + root.clone()) / a.clone();
        if is_positive(&(&self.reach - &last), dec)? {
            return Ok(None);
        }
        let first = (-b - root) / a;
        Ok(Some(if is_positive(&(&self.reach - &first), dec)? {
            self.reach.clone()
        } else {
            first
        }))
    }

    /// The first place at or past the span's end where the axis may meet a
    /// ball.
    pub fn ball_beyond(
        &self,
        centre: &Point3,
        radius_sq: &Real,
        dec: &mut Decider,
    ) -> Result<Option<Real>, Refusal> {
        let offset = centre.minus(&self.at);
        let along = offset.dot(&self.direction);
        let length_sq = self.direction.dot(&self.direction);
        let spare = radius_sq * &length_sq - (&offset.dot(&offset) * &length_sq - &along * &along);
        if is_positive(&(-spare.clone()), dec)? {
            return Ok(None);
        }
        let root = spare.sqrt();
        let last = (along.clone() + root.clone()) / length_sq.clone();
        if is_positive(&(&self.reach - &last), dec)? {
            return Ok(None);
        }
        let first = (along - root) / length_sq;
        Ok(Some(if is_positive(&(&self.reach - &first), dec)? {
            self.reach.clone()
        } else {
            first
        }))
    }

    /// Whether the zone starts outside the part, from where the axis past the
    /// span's end first crosses a face: entering the material there means it
    /// was outside, leaving means it started inside. `None` when a face this
    /// check cannot place may come first, or two crossings tie.
    pub fn starts_free(crossings: &[Crossing], dec: &mut Decider) -> Result<Option<bool>, Refusal> {
        let mut first: Option<(&Real, bool)> = None;
        let mut tied = false;
        for crossing in crossings {
            if let Crossing::At { t, exiting } = crossing {
                match first {
                    Some((best, _)) if is_zero(&(t - best), dec)? => tied = true,
                    Some((best, _)) if !is_positive(&(t - best), dec)? => {
                        first = Some((t, *exiting));
                        tied = false;
                    }
                    Some(_) => {}
                    None => first = Some((t, *exiting)),
                }
            }
        }
        for crossing in crossings {
            if let Crossing::From(from) = crossing {
                let later = match first {
                    Some((best, _)) => is_positive(&(from - best), dec)?,
                    None => false,
                };
                if !later {
                    return Ok(None);
                }
            }
        }
        if tied {
            return Ok(None);
        }
        Ok(Some(first.is_none_or(|(_, exiting)| !exiting)))
    }

    /// The side shared by every trim edge.
    pub fn trim_side(sides: &[Option<Side>]) -> Option<Side> {
        let first = sides.first().copied().flatten()?;
        sides
            .iter()
            .all(|side| *side == Some(first))
            .then_some(first)
    }

    /// Whether the axial coordinate can have extrema only on a cylinder's trim.
    pub fn cylinder_extrema_on_trim(
        &self,
        cylinder: &catallax_shapes::Cylinder,
        dec: &mut Decider,
    ) -> Result<bool, Refusal> {
        Ok(!is_zero(&cylinder.axis.dot(&self.direction), dec)?)
    }

    /// The side containing a cone's trim and apex.
    pub fn cone_trim_side(
        &self,
        sides: &[Option<Side>],
        apex: &Point3,
        dec: &mut Decider,
    ) -> Result<Option<Side>, Refusal> {
        let Some(side) = Self::trim_side(sides) else {
            return Ok(None);
        };
        Ok((self.point_side(apex, dec)? == Some(side)).then_some(side))
    }

    /// The end a point lies past, or `None` when it is within the span.
    pub fn point_side(&self, point: &Point3, dec: &mut Decider) -> Result<Option<Side>, Refusal> {
        let along = point.minus(&self.at).dot(&self.direction);
        if !is_positive(&along, dec)? {
            return Ok(Some(Side::Behind));
        }
        let past = along - &self.reach * &self.direction.dot(&self.direction);
        Ok((!is_positive(&(-past), dec)?).then_some(Side::Beyond))
    }

    /// The end a whole circle lies past, or `None` when part of it is within
    /// the span.
    pub fn circle_side(
        &self,
        centre: &Point3,
        axis: &Vector3,
        radius_sq: &Real,
        dec: &mut Decider,
    ) -> Result<Option<Side>, Refusal> {
        if circle_behind(&self.at, &self.direction, centre, axis, radius_sq, dec)? {
            return Ok(Some(Side::Behind));
        }
        let back = self.direction.scaled(&Real::from_i64(-1));
        Ok(
            circle_behind(&self.end(), &back, centre, axis, radius_sq, dec)?
                .then_some(Side::Beyond),
        )
    }
}

/// Whether a whole circle lies at or behind the plane through `at` square to
/// `direction`.
fn circle_behind(
    at: &Point3,
    direction: &Vector3,
    centre: &Point3,
    axis: &Vector3,
    radius_sq: &Real,
    dec: &mut Decider,
) -> Result<bool, Refusal> {
    let height = at.minus(centre).dot(direction);
    let radial = axis.cross(direction);
    let enough = &height * &height * axis.dot(axis) - radius_sq * &radial.dot(&radial);
    Ok(!is_positive(&(-height), dec)? && !is_positive(&(-enough), dec)?)
}

pub(super) fn check_written_projection(
    graph: &Graph,
    faces: &BTreeMap<FaceId, Vec<Id>>,
    reversed: &BTreeSet<Id>,
    control: &PositionControl<FaceId>,
    zone: &super::ProjectedZone<FaceId>,
    item: &str,
    dec: &mut Decider,
) -> Result<(), Refusal> {
    let invalid = || {
        refuse(
            item,
            "the projected zone cannot be proved on the written geometry",
        )
    };
    let entry_ids: BTreeSet<_> = zone
        .entry_faces
        .iter()
        .flat_map(|face| faces[face].iter().copied())
        .collect();
    let shells: Vec<_> = graph
        .records()
        .filter_map(|(kind, fields)| match (kind, fields) {
            ("CLOSED_SHELL", [_, Field::List(members)]) => Some(members),
            _ => None,
        })
        .collect();
    let direction = target::vector(&zone.direction);
    for (feature, entry) in control.features.iter().zip(&zone.entry_faces) {
        let [entry] = faces[entry].as_slice() else {
            return Err(invalid());
        };
        let (bounds, _, same) = written_face(graph, *entry).ok_or_else(invalid)?;
        let entry_edges = written_edges(graph, bounds).ok_or_else(invalid)?;
        let boundary = written_boundary(graph, *entry).ok_or_else(invalid)?;
        let plane = &boundary.plane;
        let normal = if same != reversed.contains(entry) {
            plane.normal.clone()
        } else {
            plane.normal.scaled(&Real::from_i64(-1))
        };
        if !is_zero_vector(&direction.cross(&normal), dec)?
            || !is_positive(&direction.dot(&normal), dec)?
        {
            return Err(invalid());
        }
        let mut first: Option<(Plane, Real)> = None;
        let mut trim = BTreeSet::new();
        for &face in feature.iter().flat_map(|face| &faces[face]) {
            let (bounds, surface, _) = written_face(graph, face).ok_or_else(invalid)?;
            let [_, frame, Field::Text(radius)] =
                fields(graph, surface, "CYLINDRICAL_SURFACE").ok_or_else(invalid)?
            else {
                return Err(invalid());
            };
            let axis = written_plane(graph, frame).ok_or_else(invalid)?;
            let radius = written_real(radius);
            if !is_zero_vector(&axis.normal.cross(&normal), dec)? {
                return Err(invalid());
            }
            if let Some((first, first_radius)) = &first {
                if !is_zero_vector(&axis.origin.minus(&first.origin).cross(&normal), dec)?
                    || !is_zero(&(&radius - first_radius), dec)?
                {
                    return Err(invalid());
                }
            } else {
                first = Some((axis, radius));
            }
            trim.extend(written_edges(graph, bounds).ok_or_else(invalid)?);
        }
        if trim.is_disjoint(&entry_edges) {
            return Err(invalid());
        }
        for &other in entry_ids.iter().filter(|id| *id != entry) {
            let (bounds, _, _) = written_face(graph, other).ok_or_else(invalid)?;
            if !trim.is_disjoint(&written_edges(graph, bounds).ok_or_else(invalid)?) {
                return Err(invalid());
            }
        }
        let mut inward = false;
        for id in trim {
            let ("EDGE_CURVE", [_, a, b, curve, _]) = graph.record(id) else {
                return Err(invalid());
            };
            let curve = written_curve(graph, curve);
            let points = if fields(graph, curve, "LINE").is_some() {
                vec![
                    written_point(graph, a).ok_or_else(invalid)?,
                    written_point(graph, b).ok_or_else(invalid)?,
                ]
            } else {
                let circle = fields(graph, curve, "CIRCLE")
                    .and_then(|fields| written_plane(graph, fields.get(1)?))
                    .ok_or_else(invalid)?;
                if !is_zero_vector(&circle.normal.cross(&normal), dec)? {
                    return Err(invalid());
                }
                vec![circle.origin]
            };
            for point in points {
                let distance = point.minus(&plane.origin).dot(&normal);
                if is_positive(&distance, dec)? {
                    return Err(invalid());
                }
                inward |= is_positive(&(-distance), dec)?;
            }
        }
        if !inward {
            return Err(invalid());
        }
        let (axis, _) = first.ok_or_else(invalid)?;
        let offset = plane.origin.minus(&axis.origin).dot(&normal) / axis.normal.dot(&normal);
        let at = axis.origin.plus(&axis.normal.scaled(&offset));
        let span = ZoneSpan::new(
            at.clone(),
            normal.clone(),
            &Real::from_rational(zone.height_mm.clone()),
        );
        let owners = shells.iter().filter(|members| {
            members
                .iter()
                .any(|f| matches!(f, Field::Ref(face) if face == entry))
        });
        if owners.count() != 1 {
            return Err(invalid());
        }
        let mut crossings = Vec::new();
        for member in shells.iter().flat_map(|shell| shell.iter()) {
            let Field::Ref(face) = member else {
                return Err(invalid());
            };
            let (bounds, surface, same) = written_face(graph, *face).ok_or_else(invalid)?;
            if let Some(plane) = fields(graph, surface, "PLANE")
                .and_then(|fields| written_plane(graph, fields.get(1)?))
            {
                match span.plane_meeting(&plane, dec)? {
                    PlaneMeeting::Miss => continue,
                    PlaneMeeting::Coplanar => {}
                    PlaneMeeting::Within(hit) => {
                        if written_boundary(graph, *face)
                            .ok_or_else(invalid)?
                            .contains(&hit, dec)?
                        {
                            return Err(invalid());
                        }
                        continue;
                    }
                    PlaneMeeting::Beyond(t, hit) => {
                        if written_boundary(graph, *face)
                            .ok_or_else(invalid)?
                            .contains(&hit, dec)?
                        {
                            let along = plane.normal.dot(&normal);
                            let outward = if same != reversed.contains(face) {
                                along
                            } else {
                                -along
                            };
                            crossings.push(Crossing::At {
                                t,
                                exiting: is_positive(&outward, dec)?,
                            });
                        }
                        continue;
                    }
                }
            }
            let edges = written_edges(graph, bounds).ok_or_else(invalid)?;
            let mut sides = Vec::new();
            for id in edges {
                let ("EDGE_CURVE", [_, a, b, curve, _]) = graph.record(id) else {
                    return Err(invalid());
                };
                let curve = written_curve(graph, curve);
                if fields(graph, curve, "LINE").is_some() {
                    for end in [a, b] {
                        let end = written_point(graph, end).ok_or_else(invalid)?;
                        sides.push(span.point_side(&end, dec)?);
                    }
                } else if let Some([_, frame, Field::Text(radius)]) = fields(graph, curve, "CIRCLE")
                {
                    let circle = written_plane(graph, frame).ok_or_else(invalid)?;
                    let radius = written_real(radius);
                    sides.push(span.circle_side(
                        &circle.origin,
                        &circle.normal,
                        &(&radius * &radius),
                        dec,
                    )?);
                } else {
                    sides.push(None);
                }
            }
            let cylinder = match fields(graph, surface, "CYLINDRICAL_SURFACE") {
                Some([_, frame, Field::Text(radius)]) => {
                    let plane = written_plane(graph, frame).ok_or_else(invalid)?;
                    let radius = written_real(radius);
                    Some(catallax_shapes::Cylinder {
                        origin: plane.origin,
                        axis: plane.normal,
                        ref_dir: plane.ref_dir,
                        radius_sq: &radius * &radius,
                    })
                }
                _ => None,
            };
            let ball = match (
                fields(graph, surface, "SPHERICAL_SURFACE"),
                fields(graph, surface, "TOROIDAL_SURFACE"),
            ) {
                (Some([_, frame, Field::Text(radius)]), _) => Some((frame, written_real(radius))),
                (_, Some([_, frame, Field::Text(major), Field::Text(minor)])) => {
                    Some((frame, written_real(major) + written_real(minor)))
                }
                _ => None,
            };
            let ball = match ball {
                Some((frame, radius)) => Some((
                    written_plane(graph, frame).ok_or_else(invalid)?.origin,
                    &radius * &radius,
                )),
                None => None,
            };
            let trim = if fields(graph, surface, "PLANE").is_some()
                || cylinder
                    .as_ref()
                    .map(|cylinder| span.cylinder_extrema_on_trim(cylinder, dec))
                    .transpose()?
                    .unwrap_or(false)
            {
                ZoneSpan::trim_side(&sides)
            } else if let Some([_, frame, Field::Text(radius), _]) =
                fields(graph, surface, "CONICAL_SURFACE")
                && is_zero(&written_real(radius), dec)?
            {
                let apex = written_plane(graph, frame).ok_or_else(invalid)?.origin;
                span.cone_trim_side(&sides, &apex, dec)?
            } else {
                None
            };
            if trim == Some(Side::Behind) {
                continue;
            }
            let clear = trim == Some(Side::Beyond)
                || match (&cylinder, &ball) {
                    (Some(cylinder), _) => span.misses_cylinder(cylinder, dec)?,
                    (_, Some((center, radius_sq))) => span.misses_ball(center, radius_sq, dec)?,
                    _ => false,
                };
            if !clear {
                return Err(invalid());
            }
            let from = match (&cylinder, &ball) {
                (Some(cylinder), _) => span.cylinder_beyond(cylinder, dec)?,
                (_, Some((center, radius_sq))) => span.ball_beyond(center, radius_sq, dec)?,
                _ => Some(span.reach()),
            };
            crossings.extend(from.map(Crossing::From));
        }
        if ZoneSpan::starts_free(&crossings, dec)? != Some(true) {
            return Err(invalid());
        }
    }
    Ok(())
}

/// Whether a datum's feature is a feature of size: one cylinder, or two
/// opposed parallel planes.
fn feature_of_size(
    model: &Model,
    feature: &Feature<FaceId>,
    item: &str,
    dec: &mut Decider,
) -> Result<bool, Refusal> {
    if feature.iter().all(|face| {
        matches!(
            model.geom().surface(model.face(*face).surface),
            Surface::Cylinder(_)
        )
    }) {
        return Ok(matches!(
            feature_shape(model, feature, item, dec),
            Ok(Shape::Axis { .. })
        ));
    }
    let [a, b] = feature.as_slice() else {
        return Ok(false);
    };
    let (Ok(a), Ok(b)) = (
        oriented_plane(model, &[*a], item, dec),
        oriented_plane(model, &[*b], item, dec),
    ) else {
        return Ok(false);
    };
    Ok(is_zero_vector(&a.normal.cross(&b.normal), dec)?
        && is_positive(&(-a.normal.dot(&b.normal)), dec)?
        && !is_zero(&b.point.minus(&a.point).dot(&a.normal), dec)?)
}

fn check_open_entry(
    model: &Model,
    body: BodyId,
    span: &ZoneSpan,
    item: &str,
    dec: &mut Decider,
) -> Result<(), Refusal> {
    if model.body(body).kind != BodyKind::Solid {
        return Err(refuse(item, "the projected zone entry has no owning shell"));
    }
    let blocked = || refuse(item, "a projected zone's axis meets a closed end");
    let unknown = || refuse(item, "the projected zone's exterior start cannot be proved");
    let mut crossings = Vec::new();
    for face in model.faces_of_body(body) {
        let surface = model.geom().surface(model.face(face).surface);
        if let Surface::Plane(plane) = surface {
            match span.plane_meeting(plane, dec)? {
                PlaneMeeting::Miss => continue,
                PlaneMeeting::Coplanar => {}
                PlaneMeeting::Within(hit) => {
                    if face_boundary(model, face, dec)?
                        .ok_or_else(unknown)?
                        .contains(&hit, dec)?
                    {
                        return Err(blocked());
                    }
                    continue;
                }
                PlaneMeeting::Beyond(t, hit) => {
                    if face_boundary(model, face, dec)?
                        .ok_or_else(unknown)?
                        .contains(&hit, dec)?
                    {
                        let (_, outward) = outward_plane(model, face).ok_or_else(unknown)?;
                        crossings.push(Crossing::At {
                            t,
                            exiting: is_positive(&outward.dot(&span.direction), dec)?,
                        });
                    }
                    continue;
                }
            }
        }
        let mut sides = Vec::new();
        for id in face_edges(model, face) {
            let edge = model.edge(id);
            match curve_form(model.geom(), edge.curve, dec)? {
                CurveForm::Line { .. } => {
                    for end in [edge.start, edge.end] {
                        let end = end.ok_or_else(unknown)?;
                        let end = model.geom().position(model.vertex(end).point);
                        sides.push(span.point_side(&end, dec)?);
                    }
                }
                CurveForm::Circle {
                    center,
                    normal: axis,
                    radius_sq,
                    ..
                } => sides.push(span.circle_side(&center, &axis, &radius_sq, dec)?),
                _ => sides.push(None),
            }
        }
        let trim = match surface {
            Surface::Plane(_) => ZoneSpan::trim_side(&sides),
            Surface::Cylinder(cylinder) if span.cylinder_extrema_on_trim(cylinder, dec)? => {
                ZoneSpan::trim_side(&sides)
            }
            Surface::Cone(cone) => span.cone_trim_side(&sides, &cone.apex, dec)?,
            _ => None,
        };
        if trim == Some(Side::Behind) {
            continue;
        }
        let ball = |center: &Point3, radius: Real| (center.clone(), &radius * &radius);
        let bound = match surface {
            Surface::Sphere(sphere) => Some(ball(&sphere.center, sphere.radius.clone())),
            Surface::Torus(torus) => Some(ball(
                &torus.center,
                torus.major.clone() + torus.minor.clone(),
            )),
            _ => None,
        };
        let clear = trim == Some(Side::Beyond)
            || match (surface, &bound) {
                (Surface::Cylinder(cylinder), _) => span.misses_cylinder(cylinder, dec)?,
                (_, Some((center, radius_sq))) => span.misses_ball(center, radius_sq, dec)?,
                _ => false,
            };
        if !clear {
            return Err(unknown());
        }
        let from = match (surface, &bound) {
            (Surface::Cylinder(cylinder), _) => span.cylinder_beyond(cylinder, dec)?,
            (_, Some((center, radius_sq))) => span.ball_beyond(center, radius_sq, dec)?,
            _ => Some(span.reach()),
        };
        crossings.extend(from.map(Crossing::From));
    }
    match ZoneSpan::starts_free(&crossings, dec)? {
        Some(true) => Ok(()),
        Some(false) => Err(refuse(
            item,
            "the projected zone's axis starts inside the part",
        )),
        None => Err(unknown()),
    }
}

fn check_inward_trim(
    model: &Model,
    edges: &BTreeSet<catallax_topology::EdgeId>,
    origin: &Point3,
    normal: &Vector3,
    item: &str,
    dec: &mut Decider,
) -> Result<(), Refusal> {
    let unknown = || refuse(item, "the cylinder's axial trim cannot be checked");
    let mut inward = false;
    let mut check = |point: &Point3, dec: &mut Decider| -> Result<(), Refusal> {
        let distance = point.minus(origin).dot(normal);
        if is_positive(&distance, dec)? {
            return Err(refuse(
                item,
                "the projected zone points back into its controlled cylinder",
            ));
        }
        inward |= is_positive(&(-distance), dec)?;
        Ok(())
    };
    for edge in edges {
        let edge = model.edge(*edge);
        match curve_form(model.geom(), edge.curve, dec)? {
            CurveForm::Line { .. } => {
                for vertex in [edge.start, edge.end] {
                    let vertex = vertex.ok_or_else(unknown)?;
                    check(&model.geom().position(model.vertex(vertex).point), dec)?;
                }
            }
            CurveForm::Circle {
                center,
                normal: axis,
                ..
            } if is_zero_vector(&axis.cross(normal), dec)? => {
                check(&center, dec)?;
            }
            _ => return Err(unknown()),
        }
    }
    if !inward {
        return Err(unknown());
    }
    Ok(())
}

fn check_position(
    model: &Model,
    control: &PositionControl<FaceId>,
    datums: &BTreeMap<&str, &Feature<FaceId>>,
    targeted: &BTreeSet<&str>,
    item: &str,
    body: BodyId,
    dec: &mut Decider,
) -> Result<(), Refusal> {
    let mut axes = Vec::with_capacity(control.features.len());
    for feature in &control.features {
        match feature_shape(model, feature, item, dec) {
            Ok(Shape::Axis {
                point, direction, ..
            }) => axes.push((point, direction)),
            Ok(Shape::Plane { .. }) | Err(_) => {
                return Err(refuse(
                    item,
                    "a position control's features are each one cylinder",
                ));
            }
        }
    }
    for segment in control.frame.segments() {
        for datum in &segment.datums {
            if datum.material.is_none() {
                continue;
            }
            let letter = datum.letter.as_str();
            if targeted.contains(letter) {
                return Err(refuse(
                    item,
                    format!(
                        "datum {letter} is set up by datum targets and has no material boundary"
                    ),
                ));
            }
            let eligible = match datums.get(letter) {
                Some(feature) => feature_of_size(model, feature, item, dec)?,
                None => false,
            };
            if !eligible {
                return Err(refuse(
                    item,
                    format!("datum {letter} is not a feature of size and has no material boundary"),
                ));
            }
        }
    }
    if let PositionFrame::Single {
        projected_zone: Some(zone),
        ..
    } = &control.frame
    {
        let direction = target::vector(&zone.direction);
        if zone.entry_faces.len() != control.features.len() {
            return Err(refuse(
                item,
                "a projected zone needs one entry face per feature",
            ));
        }
        let entries: BTreeMap<_, BTreeSet<_>> = zone
            .entry_faces
            .iter()
            .map(|entry| (*entry, face_edges(model, *entry).collect()))
            .collect();
        for ((feature, (point, axis)), entry) in
            control.features.iter().zip(&axes).zip(&zone.entry_faces)
        {
            let (plane, normal) = outward_plane(model, *entry)
                .ok_or_else(|| refuse(item, "a projected zone's entry face is planar"))?;
            if !is_zero_vector(&axis.cross(&normal), dec)? {
                return Err(refuse(
                    item,
                    "a controlled axis is not square to its entry face",
                ));
            }
            if !is_zero_vector(&direction.cross(&normal), dec)?
                || !is_positive(&direction.dot(&normal), dec)?
            {
                return Err(refuse(
                    item,
                    "the projected zone does not point out of the part from its entry face",
                ));
            }
            let edges: BTreeSet<_> = feature
                .iter()
                .flat_map(|face| face_edges(model, *face))
                .collect();
            let mut meeting = entries.iter().filter(|(_, own)| !own.is_disjoint(&edges));
            if meeting.next().map(|(face, _)| face) != Some(entry) || meeting.next().is_some() {
                return Err(refuse(
                    item,
                    "each controlled feature needs exactly one projected zone entry face it meets",
                ));
            }
            let offset = plane.origin.minus(point).dot(&normal) / axis.dot(&normal);
            check_inward_trim(model, &edges, &plane.origin, &normal, item, dec)?;
            let span = ZoneSpan::new(
                point.plus(&axis.scaled(&offset)),
                normal.clone(),
                &Real::from_rational(zone.height_mm.clone()),
            );
            check_open_entry(model, body, &span, item, dec)?;
        }
    }
    Ok(())
}

/// Checks the datum targets and position controls of `pmi` against `body`:
/// each target lies on its trimmed planar face, each controlled feature is a
/// cylinder, a material boundary is only on a datum feature of size that no
/// target sets up, and a projected zone starts on a face that meets its
/// feature square to its axis and points out of the part. Refuses
/// `step.pmi-invalid` with the item's name. Other items are not checked here.
pub fn check_on_model(
    model: &Model,
    body: BodyId,
    pmi: &Pmi<FaceId>,
    dec: &mut Decider,
) -> Result<(), Refusal> {
    let faces: BTreeSet<_> = model.faces_of_body(body).into_iter().collect();
    let mut datums = BTreeMap::new();
    let mut targeted = BTreeSet::new();
    for item in &pmi.items {
        if item
            .kind
            .faces()
            .into_iter()
            .any(|face| !faces.contains(face))
        {
            return Err(refuse(
                &item.id,
                "a face it is attached to is not in the selected body",
            ));
        }
        match &item.kind {
            ItemKind::Datum(datum) => {
                datums.insert(datum.letter.as_str(), &datum.feature);
            }
            ItemKind::DatumTarget(target) => {
                targeted.insert(target.datum.as_str());
            }
            _ => {}
        }
    }
    for item in &pmi.items {
        match &item.kind {
            ItemKind::DatumTarget(target) => check_target(model, target, &item.id, dec)?,
            ItemKind::Position(control) => {
                check_position(model, control, &datums, &targeted, &item.id, body, dec)?;
            }
            _ => {}
        }
    }
    Ok(())
}

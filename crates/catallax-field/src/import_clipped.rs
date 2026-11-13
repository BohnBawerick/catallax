//! Spheres with disjoint planar cuts and convex meridian sectors of ring tori.

use crate::{Face, geometry, import_step};
use catallax_numbers::{Decider, Rational, Refusal, Sign};
use catallax_shapes::{Curve, Plane, RPoint3, RVector3, Surface};
use catallax_topology::{BodyId, EdgeId, FaceId, Loop, Model, Sense};
use std::collections::BTreeSet;

fn q(n: i64) -> Rational {
    Rational::from_i64(n)
}
fn fail(detail: impl Into<String>) -> Refusal {
    import_step::unsupported(detail)
}
fn point(p: &catallax_shapes::Point3) -> Result<RPoint3, Refusal> {
    p.as_rational()
        .ok_or_else(|| fail("clip point is not rational after placement"))
}
fn unit(v: &catallax_shapes::Vector3) -> Result<RVector3, Refusal> {
    let v = v
        .as_rational()
        .ok_or_else(|| fail("clip direction is not rational"))?;
    let squared = v.norm_squared();
    let length = import_step::rational_sqrt(&squared).ok_or_else(|| {
        fail(format!(
            "clip direction squared length {squared} has no positive rational square root"
        ))
    })?;
    Ok(v.scaled(&(&q(1) / &length)))
}

struct Disk {
    id: FaceId,
    edge: EdgeId,
    plane: Plane,
    normal: RVector3,
    centre: RPoint3,
    radius_sq: Rational,
}

fn disk(model: &Model, id: FaceId) -> Result<Disk, Refusal> {
    let plane = import_step::outward(model, id)?;
    let normal = unit(&plane.normal)?;
    let [lp] = model.face(id).loops.as_slice() else {
        return Err(fail(format!(
            "face {id:?}: expected one disk boundary, found {}",
            model.face(id).loops.len()
        )));
    };
    let Loop::Coedges(cs) = model.loop_(*lp) else {
        return Err(fail("a clipping disk has a vertex boundary"));
    };
    let [c] = cs.as_slice() else {
        return Err(fail(format!(
            "a clipping disk has {} boundary arcs; expected one whole circle",
            cs.len()
        )));
    };
    let c = model.coedge(*c);
    let edge = model.edge(c.edge);
    let Curve::Circle(circle) = model.geom().curve(edge.curve) else {
        return Err(fail("a clipping disk boundary is not circular"));
    };
    if edge.start.is_some() || edge.end.is_some() {
        return Err(fail(
            "a clipping disk boundary is split; expected one whole circle",
        ));
    }
    let centre = point(&circle.center)?;
    let circle_normal = unit(&circle.normal)?;
    let alignment = circle_normal.dot(&normal);
    let travel = if c.sense == edge.sense { q(1) } else { q(-1) };
    let offset = centre.minus(&point(&plane.origin)?).dot(&normal);
    if !circle_normal.cross(&normal).is_zero() || alignment * travel != q(1) || offset != q(0) {
        return Err(fail(format!(
            "disk {id:?}: boundary orientation or plane incidence failed; plane offset {offset}"
        )));
    }
    let radius_sq = import_step::rational(&circle.radius_sq, "clip circle radius squared")?;
    if radius_sq.sign() != Sign::Positive {
        return Err(fail(format!(
            "disk radius squared {radius_sq} must be positive"
        )));
    }
    Ok(Disk {
        id,
        edge: c.edge,
        plane,
        normal,
        centre,
        radius_sq,
    })
}

/// Every curved boundary must meet exactly one disk, with opposite travel.
fn boundaries(model: &Model, curved: FaceId, disks: &[Disk]) -> Result<(), Refusal> {
    let mut edges = BTreeSet::new();
    for lp in &model.face(curved).loops {
        let Loop::Coedges(cs) = model.loop_(*lp) else {
            return Err(fail("clipped curved face has a vertex boundary"));
        };
        let [c] = cs.as_slice() else {
            return Err(fail(format!(
                "curved boundary has {} arcs; expected one whole circle",
                cs.len()
            )));
        };
        let c = model.coedge(*c);
        let Some(disk) = disks.iter().find(|d| d.edge == c.edge) else {
            return Err(fail(format!(
                "curved boundary {:?} has no closing disk",
                c.edge
            )));
        };
        let Loop::Coedges(dc) = model.loop_(model.face(disk.id).loops[0]) else {
            unreachable!()
        };
        if c.sense == model.coedge(dc[0]).sense || !edges.insert(c.edge) {
            return Err(fail(
                "clipped face boundaries repeat or travel with their closing disks",
            ));
        }
    }
    if edges.len() != disks.len() {
        return Err(fail(format!(
            "{} curved boundaries close {} disks",
            edges.len(),
            disks.len()
        )));
    }
    Ok(())
}

fn sphere_volume(
    centre: &RPoint3,
    radius: &Rational,
    disks: &[Disk],
) -> Result<geometry::CellVolume, Refusal> {
    let mut factor = q(4) * radius.clone() * radius.clone() * radius.clone() / q(3);
    for disk in disks {
        let delta = disk.centre.minus(centre);
        let h = delta.dot(&disk.normal);
        let residual = &h * &h + disk.radius_sq.clone() - radius * radius;
        if !delta.cross(&disk.normal).is_zero() || residual != q(0) || h <= q(0) || h >= *radius {
            return Err(fail(format!(
                "sphere disk {:?}: height {h}, radius squared {}, sphere radius {radius}, incidence residual {residual}; cuts must remove minor caps",
                disk.id, disk.radius_sq
            )));
        }
        for other in disks.iter().filter(|d| d.id != disk.id) {
            let gap = other.centre.minus(&disk.centre).dot(&other.normal);
            let cosine = disk.normal.dot(&other.normal);
            let margin = &gap * &gap - &disk.radius_sq * &(q(1) - &cosine * &cosine);
            if gap <= q(0) || margin <= q(0) {
                return Err(fail(format!(
                    "sphere disks {:?} and {:?}: separation gap {gap}, squared margin {margin}; crossing or touching cap boundaries need arc domains",
                    disk.id, other.id
                )));
            }
        }
        let cap = radius - &h;
        factor = factor - &cap * &cap * (&(&q(3) * radius) - &cap) / q(3);
    }
    Ok(geometry::CellVolume {
        factor,
        pi_power: 1,
        angle: None,
    })
}

fn torus_volume(
    torus: &catallax_shapes::Torus,
    disks: &[Disk],
) -> Result<geometry::CellVolume, Refusal> {
    let centre = point(&torus.center)?;
    let axis = unit(&torus.axis)?;
    let major = import_step::rational(&torus.major, "torus major radius")?;
    let minor = import_step::rational(&torus.minor, "torus minor radius")?;
    if minor <= q(0) || major <= minor || disks.len() != 2 {
        return Err(fail(format!(
            "torus major {major}, minor {minor}, {} disks; expected a ring torus and two meridian disks",
            disks.len()
        )));
    }
    let mut radial = Vec::new();
    for disk in disks {
        let v = disk.centre.minus(&centre);
        let axial = v.dot(&axis);
        let plane_offset = v.dot(&disk.normal);
        let radial_sq = v.norm_squared();
        if axial != q(0)
            || plane_offset != q(0)
            || disk.normal.dot(&axis) != q(0)
            || radial_sq != &major * &major
            || disk.radius_sq != &minor * &minor
        {
            return Err(fail(format!(
                "torus disk {:?}: axial offset {axial}, meridian offset {plane_offset}, centre radius squared {radial_sq}, disk radius squared {}; expected {major} and {minor} radii",
                disk.id, disk.radius_sq
            )));
        }
        radial.push(v.scaled(&(&q(1) / &major)));
    }
    let a = disks[0].normal.dot(&radial[1]);
    let b = disks[1].normal.dot(&radial[0]);
    if a >= q(0) || b >= q(0) {
        return Err(fail(format!(
            "torus sector cross-plane signs {a}, {b}; expected a convex sector strictly between 0 and pi with distinct cap planes"
        )));
    }
    let cosine = radial[0].dot(&radial[1]);
    let sine = radial[0].cross(&radial[1]).dot(&axis).abs();
    Ok(geometry::CellVolume {
        factor: major * (&minor * &minor),
        pi_power: 1,
        angle: Some([sine, cosine]),
    })
}

pub(crate) fn cell(
    model: &Model,
    body: BodyId,
    dec: &mut Decider,
) -> Result<(geometry::Leaf, Vec<Face>), Refusal> {
    let ids = model.faces_of_body(body);
    let curved: Vec<_> = ids
        .iter()
        .copied()
        .filter(|id| {
            matches!(
                model.geom().surface(model.face(*id).surface),
                Surface::Sphere(_) | Surface::Torus(_)
            )
        })
        .collect();
    let [curved] = curved.as_slice() else {
        return Err(fail(format!(
            "{} analytic faces, {} spherical or toroidal carriers; clipped import needs exactly one curved carrier and closing disks",
            ids.len(),
            curved.len()
        )));
    };
    if model.face(*curved).sense != Sense::Same {
        return Err(fail("clipped sphere or torus faces must point outward"));
    }
    let mut disks = Vec::new();
    for id in ids.iter().copied().filter(|id| id != curved) {
        if !matches!(
            model.geom().surface(model.face(id).surface),
            Surface::Plane(_)
        ) {
            return Err(fail(format!(
                "{} faces: clipped analytic import needs planar closing disks",
                ids.len()
            )));
        }
        disks.push(disk(model, id)?);
    }
    if disks.len() < 2 {
        return Err(fail(format!(
            "{} disks; two-axis clipped import needs at least two",
            disks.len()
        )));
    }
    boundaries(model, *curved, &disks)?;
    let surface = model.geom().surface(model.face(*curved).surface);
    let (centre, extent, volume) = match surface {
        Surface::Sphere(sphere) => {
            let centre = point(&sphere.center)?;
            let radius = import_step::rational(&sphere.radius, "sphere radius")?;
            let volume = sphere_volume(&centre, &radius, &disks)?;
            (centre, radius, volume)
        }
        Surface::Torus(torus) => (
            point(&torus.center)?,
            import_step::rational(&(&torus.major + &torus.minor), "torus outer radius")?,
            torus_volume(torus, &disks)?,
        ),
        _ => unreachable!(),
    };
    let radius = extent + centre.x.abs() + centre.y.abs() + centre.z.abs();
    let mut surfaces = vec![surface.clone()];
    surfaces.extend(disks.iter().map(|d| Surface::Plane(d.plane.clone())));
    let mut faces = Vec::new();
    for (id, surface) in std::iter::once(*curved)
        .chain(disks.iter().map(|d| d.id))
        .zip(&surfaces)
    {
        for origin in model.face_origins(id) {
            faces.push(geometry::face(surface.clone(), origin.clone(), dec)?);
        }
    }
    let mut cell = geometry::Cell::new(surfaces, Some(radius))?;
    cell.set_volume(volume);
    Ok((geometry::Leaf::Cell(cell), faces))
}

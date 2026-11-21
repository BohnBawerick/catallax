//! A plane with a cylinder, cone or sphere in any pose (4.2a step 2).
//!
//! The field derives the ball path from the two supports and the radius alone,
//! as step 1 did for lines and perpendicular circles: offset both supports by
//! `r` toward the ball and intersect the offsets. A cylinder or cone nappe of
//! slope `k` offset by `r` is the coaxial surface whose radius at the same
//! axial station grows by `b r sqrt(1 + k^2)`, with `b` = +1 for a ball outside
//! the support and -1 inside. A sphere offset is the concentric sphere of
//! radius `R + b r`. The offset plane meets the offset surface in a conic: a
//! circle when the plane is perpendicular to the axis (always, for a sphere),
//! otherwise an ellipse. Nothing here reads a B-rep fillet module or a fillet
//! definition; the B-rep's `rolling_pipe` builds its ellipse its own way.
//!
//! With the unit plane normal `n` toward the ball, the plane `n.x = lambda`,
//! the unit axis `a`, and `u1`, `u2` perpendicular to `a` with `n.u2 = 0`, the
//! station at angle `theta` lies on the meridian `m = cos u1 + sin u2`. With
//! `D = n.a + k cos (n.u1)`, the offset plane fixes its axial coordinate:
//! `c = o + s a + (G / D) m` and `s = (N - rho' cos (n.u1)) / D`, where `rho'`
//! is the offset radius at the origin `o`, `N = lambda + r - n.o` and
//! `G = rho' n.a + k N` are constants. `(n.a)^2 (1 + k^2) > k^2` keeps `D`
//! away from zero, so the conic is closed, and `G D > 0` keeps its radius
//! `G / D` positive. About the pivot `p = o + (N / n.a) a`, where the axis
//! meets the offset plane, `s - N / n.a = -(G / D) cos (n.u1) / n.a`, so
//! `c = p + (G / D) (m - a cos (n.u1) / n.a)`: one divided scalar, which keeps
//! interval enclosures tight however far the apex lies.
//!
//! Circles have the closed forms of step 1: the nearest station is the
//! meridian through the point, as in OpenCASCADE's `Extrema_ExtPElC`, and the
//! separation from the balls is `(rho - C)^2 + h^2 - r^2`. For `C > r` this has
//! the sign of the ring-torus polynomial; for a spindle (`0 < C <= r`) the
//! square root of `rho^2` is kept. A ball exactly as wide as a cylinder it
//! rolls inside, on a plane normal to the axis, has `C = 0` (4.2a step 5):
//! every ball is one ball on the axis, the separation is `|x - c|^2 - r^2`
//! and the face is a sphere, as OpenCASCADE's `ChFiKPart_ComputeData_FilPlnCyl`
//! makes it. A torus cut by a plane normal to its axis (4.2a step 6) is a
//! sphere whose meridian circle sits at the core radius `R` off the axis:
//! its offset is the torus of tube radius `a + b r`, the offset plane meets
//! it in the circle `rho = R +- sqrt((a + b r)^2 - h0^2)` on the edge's side
//! of the core circle, or in one circle where the offset plane touches the
//! offset tube, and the normal line of the meridian wall passes through the
//! core circle's point `(R, h = 0)`. Ellipses use the whole-interval panel cover
//! of step 1. OpenCASCADE reaches the nearest point of an ellipse through
//! `math_TrigonometricFunctionRoots`, a quartic in the half angle; a sampled or
//! root-polished nearest point gives no guaranteed range, the cover does.
//!
//! The fan is the edge's Voronoi zone in the solid of ball centres: the ball
//! side of both supports, the points whose projection along `n` onto the offset
//! plane falls beyond the offset conic (outside the plane offset's face), and,
//! in the point's own meridian, the edge side of the normal line of the offset
//! surface through the station (outside its face). For a sphere that normal
//! line passes through the centre. Two bands keep far zones of these
//! unbounded walls out, and are redundant inside the region: within the
//! normal section of a station the region lies in the quadrilateral of the
//! edge point, the two feet and the centre, so it lies within `2r` of the plane
//! (a linear function, at most `r (1 - n.N) <= 2r` at the far foot) and, for a
//! ball outside the support, within `2r` of the support (a convex distance,
//! whose maximum is at a vertex, at most `r + r`). A ball inside the support
//! needs no support band: the inside of a cylinder, a cone nappe or a sphere
//! within `2r` of the cutting plane is already bounded.

use super::elementary::{
    Flat, balls, length, local, point_balls, position, positive, reach, shift, vector,
};
use super::path::{Cover, Frame, Path};
use super::{BITS, gap, interval, q};
use crate::geometry::{enclose, sign};
use catallax_numbers::{Ball, Decider, Rational, Real, Refusal, Sign};
use catallax_shapes::{Isometry, Point3, Vector3};

/// A support of revolution for a [`Revolved`] path.
#[derive(Clone, Debug)]
pub enum Revolution {
    /// A circular cylinder of positive `radius` about the axis through `point`.
    Cylinder {
        /// Any point of the axis.
        point: Point3,
        /// Nonzero axis direction.
        axis: Vector3,
        /// Positive radius in mm.
        radius: Real,
    },
    /// The cone nappe from `apex` in the direction `axis`, whose radius grows by
    /// the positive `slope` per unit of axial distance.
    Cone {
        /// Apex.
        apex: Point3,
        /// Nonzero direction into the nappe.
        axis: Vector3,
        /// Positive radius gained per unit axial distance.
        slope: Real,
    },
    /// A sphere of positive `radius`.
    Sphere {
        /// Centre.
        centre: Point3,
        /// Positive radius in mm.
        radius: Real,
    },
    /// A ring torus about the axis through `centre` along `axis`, whose tube
    /// of radius `minor` circles the axis at `major` (4.2a step 6). Only a
    /// plane normal to the axis is taken. `edge` is a point of the edge: its
    /// side of the tube's core circle chooses the centre circle.
    Torus {
        /// Centre.
        centre: Point3,
        /// Nonzero axis direction.
        axis: Vector3,
        /// Radius of the core circle in mm.
        major: Real,
        /// Radius of the tube in mm, less than `major`.
        minor: Real,
        /// A point of the edge.
        edge: Point3,
    },
}

#[derive(Clone, Debug)]
enum Shape {
    /// Cylinder (`slope` zero) or cone: radius `rim + slope s` at axial
    /// distance `s` from `origin` along the unit `axis`.
    Linear {
        origin: Point3,
        axis: Vector3,
        rim: Real,
        slope: Real,
        cone: bool,
        // sqrt(1 + k^2), the radial growth of a unit normal offset.
        root: Real,
        // rho' = rim + b r root.
        offset: Real,
        // N = lambda + r - n.o and G = rho' n.a + k N.
        start: Real,
        gain: Real,
        // The same constants for the edge conic, without the offsets.
        edge_start: Real,
        edge_gain: Real,
        // n.a and n.u1.
        cross: Real,
        tilt: Real,
    },
    Sphere {
        centre: Point3,
        radius: Real,
        // R' = R + b r, the offset plane's height h0 above the centre, and
        // the edge circle's height and radius.
        offset: Real,
        height: Real,
        edge: Real,
        edge_major: Real,
    },
    Torus {
        centre: Point3,
        // The core circle's radius R and the tube's radius a.
        major: Real,
        minor: Real,
        // a' = a + b r, the offset plane's height h0 above the centre, the
        // edge plane's height and the edge circle's radius.
        offset: Real,
        height: Real,
        edge: Real,
        edge_major: Real,
        // +1 or -1: the sign of h0 (rho - R) - (C - R) h at the edge, the
        // meridian wall's edge side.
        facing: Real,
    },
}

#[derive(Clone, Debug)]
struct Circle {
    centre: Point3,
    axis: Vector3,
    major: Real,
    major_sq: Real,
    // C > r: the ring-torus polynomial has the separation's sign.
    ring: bool,
    // C = 0: the centre stays on the axis and every ball is one ball.
    axial: bool,
}

/// The rolling-ball path between a plane and a cylinder, cone, sphere or
/// torus.
#[derive(Clone, Debug)]
pub struct Revolved {
    normal: Vector3,
    point: Point3,
    unit: Vector3,
    level: Real,
    radius: Real,
    outside: bool,
    shape: Shape,
    basis: [Vector3; 2],
    circle: Option<Circle>,
    sector: Option<Sector>,
}

/// The part of a circle path between two meridian half-planes, at most half
/// a turn (4.2a step 4). `ends` are the unit radial directions at the
/// sector's start and end, in turning order about the circle's axis;
/// `normals` point from each end's half-plane into the sector. A meridian
/// direction `m` lies in the sector when `normals[k] . m >= 0` for both.
#[derive(Clone, Debug)]
struct Sector {
    centre: Point3,
    ends: [Vector3; 2],
    normals: [Vector3; 2],
}

pub(super) fn root(x: &Real) -> Real {
    x.as_rational()
        .and_then(|s| crate::profile::sqrt(&s))
        .map(Real::from)
        .unwrap_or_else(|| x.sqrt())
}
fn scale(v: &Vector3) -> Vector3 {
    v.scaled(&(&Real::one() / &length(v)))
}
fn side(outside: bool) -> Real {
    Real::from_i64(if outside { 1 } else { -1 })
}
fn magnitude(x: &Real, dec: &mut Decider) -> Result<Real, Refusal> {
    Ok(if sign(x, dec)? == Sign::Negative {
        -x
    } else {
        x.clone()
    })
}
/// An orthonormal pair perpendicular to the unit `axis`.
pub(crate) fn basis_of(axis: &Vector3, dec: &mut Decider) -> Result<[Vector3; 2], Refusal> {
    for k in [[1, 0, 0], [0, 1, 0]] {
        let k = Vector3::new(k[0].into(), k[1].into(), k[2].into());
        let e = axis.cross(&k);
        if positive(&e.norm_squared(), dec)? {
            let e = scale(&e);
            return Ok([axis.cross(&e), e]);
        }
    }
    Err(gap("no basis perpendicular to the axis"))
}
fn sum(a: &[Ball; 3], b: &[Ball; 3]) -> [Ball; 3] {
    std::array::from_fn(|i| &a[i] + &b[i])
}
fn scaled(t: &Ball, a: &[Ball; 3]) -> [Ball; 3] {
    std::array::from_fn(|i| t * &a[i])
}
fn divide(a: &Ball, b: &Ball) -> Result<Ball, Refusal> {
    a.checked_div(b)
        .ok_or_else(|| gap("the conic's denominator enclosure contains zero"))
}
pub(crate) fn hull_points(points: &[[Ball; 3]]) -> [Ball; 3] {
    std::array::from_fn(|i| {
        points
            .iter()
            .map(|p| p[i].clone())
            .reduce(|a, b| a.hull(&b))
            .expect("a nonempty point list")
    })
}
/// The cube around the disc whose diameter is `a b`. A minor circular arc
/// from `a` to `b` sees `a b` at an obtuse angle, so it lies in that disc.
pub(super) fn chord_cube(a: &[Ball; 3], b: &[Ball; 3]) -> Result<[[Ball; 3]; 2], Refusal> {
    let half = interval(&(Rational::one() / q(2)), &(Rational::one() / q(2)));
    let mid: [Ball; 3] = std::array::from_fn(|i| &(&a[i] + &b[i]) * &half);
    let d: [Ball; 3] = std::array::from_fn(|i| &a[i] - &b[i]);
    let reach = (&(&d[0].square() + &d[1].square()) + &d[2].square())
        .sqrt()
        .ok_or_else(|| gap("negative chord enclosure"))?;
    let reach = &reach * &half;
    Ok([
        std::array::from_fn(|i| &mid[i] - &reach),
        std::array::from_fn(|i| &mid[i] + &reach),
    ])
}

impl Revolved {
    /// Derives the ball path from the plane, whose `normal` points to the
    /// ball, the support of revolution, the ball's side of it and the radius.
    pub fn new(
        plane: Flat,
        support: Revolution,
        outside: bool,
        radius: Real,
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        dec.charge_work(1)?;
        if !positive(&radius, dec)? {
            return Err(reach("the rolling radius must be positive"));
        }
        if !positive(&plane.normal.norm_squared(), dec)? {
            return Err(reach("the plane needs a nonzero normal"));
        }
        let unit = scale(&plane.normal);
        let level = unit.dot(&position(&plane.point));
        let b = side(outside);
        let (shape, basis, circle) = match support {
            Revolution::Sphere {
                centre,
                radius: big,
            } => {
                if !positive(&big, dec)? {
                    return Err(reach("the sphere radius must be positive"));
                }
                let offset = &big + &(&b * &radius);
                if !positive(&offset, dec)? {
                    return Err(reach(
                        "the ball is no smaller than the sphere it rolls inside",
                    ));
                }
                let lift = unit.dot(&position(&centre));
                let height = &(&level + &radius) - &lift;
                let major_sq = &offset.square() - &height.square();
                if !positive(&major_sq, dec)? {
                    return Err(reach("the offset plane misses the offset sphere"));
                }
                let edge = &level - &lift;
                let edge_sq = &big.square() - &edge.square();
                if !positive(&edge_sq, dec)? {
                    return Err(reach("the plane does not cut the sphere"));
                }
                let basis = basis_of(&unit, dec)?;
                let ring = positive(&(&major_sq - &radius.square()), dec)?;
                let circle = Circle {
                    centre: centre.plus(&unit.scaled(&height)),
                    axis: unit.clone(),
                    major: root(&major_sq),
                    major_sq,
                    ring,
                    axial: false,
                };
                (
                    Shape::Sphere {
                        centre,
                        radius: big,
                        offset,
                        height,
                        edge,
                        edge_major: root(&edge_sq),
                    },
                    basis,
                    Some(circle),
                )
            }
            Revolution::Torus {
                centre,
                axis,
                major: big,
                minor,
                edge: point,
            } => {
                // The offset torus has tube radius a' = a + b r about the same
                // core circle. A plane normal to the axis meets it where
                // (rho - R)^2 + h0^2 = a'^2: two circles, or one where the
                // offset plane touches the offset tube's top (4.2a step 6).
                if !positive(&minor, dec)? || !positive(&(&big - &minor), dec)? {
                    return Err(reach("the torus needs a tube inside its core circle"));
                }
                if sign(&axis.cross(&unit).norm_squared(), dec)? != Sign::Zero {
                    return Err(reach(
                        "a plane not normal to the torus axis meets its offset in no circle",
                    ));
                }
                let offset = &minor + &(&b * &radius);
                if !positive(&offset, dec)? {
                    return Err(reach(
                        "the ball is no smaller than the tube it rolls inside",
                    ));
                }
                let lift = unit.dot(&position(&centre));
                let height = &(&level + &radius) - &lift;
                let spread_sq = &offset.square() - &height.square();
                if sign(&spread_sq, dec)? == Sign::Negative {
                    return Err(reach("the offset plane misses the offset torus"));
                }
                let edge = &level - &lift;
                let edge_sq = &minor.square() - &edge.square();
                if !positive(&edge_sq, dec)? {
                    return Err(reach("the plane does not cut the torus"));
                }
                // The edge's side of the core circle picks the centre circle.
                let v = point.minus(&centre);
                let w = v.minus(&unit.scaled(&unit.dot(&v)));
                let branch = match sign(&(&w.norm_squared() - &big.square()), dec)? {
                    Sign::Positive => Real::one(),
                    Sign::Negative => Real::from_i64(-1),
                    Sign::Zero => return Err(reach("the edge lies over the core circle")),
                };
                let spread = root(&spread_sq);
                let reach_edge = root(&edge_sq);
                let major = &big + &(&branch * &spread);
                if !positive(&major, dec)? {
                    return Err(reach("the centre circle reaches the torus axis"));
                }
                let facing = match sign(
                    &(&branch * &(&(&height * &reach_edge) - &(&spread * &edge))),
                    dec,
                )? {
                    Sign::Positive => Real::one(),
                    Sign::Negative => Real::from_i64(-1),
                    Sign::Zero => {
                        return Err(reach("the edge lies on the offset tube's normal line"));
                    }
                };
                let basis = basis_of(&unit, dec)?;
                let ring = positive(&(&major - &radius), dec)?;
                let circle = Circle {
                    centre: centre.plus(&unit.scaled(&height)),
                    axis: unit.clone(),
                    major_sq: major.square(),
                    major,
                    ring,
                    axial: false,
                };
                (
                    Shape::Torus {
                        centre,
                        edge_major: &big + &(&branch * &reach_edge),
                        major: big,
                        minor,
                        offset,
                        height,
                        edge,
                        facing,
                    },
                    basis,
                    Some(circle),
                )
            }
            Revolution::Cylinder { .. } | Revolution::Cone { .. } => {
                let (origin, axis, rim, slope, cone) = match support {
                    Revolution::Cylinder {
                        point,
                        axis,
                        radius,
                    } => (point, axis, radius, Real::zero(), false),
                    Revolution::Cone { apex, axis, slope } => {
                        (apex, axis, Real::zero(), slope, true)
                    }
                    Revolution::Sphere { .. } | Revolution::Torus { .. } => {
                        unreachable!("matched above")
                    }
                };
                if !positive(&axis.norm_squared(), dec)? {
                    return Err(reach("the support needs a nonzero axis"));
                }
                if (!cone && !positive(&rim, dec)?) || (cone && !positive(&slope, dec)?) {
                    return Err(reach("a cylinder radius and a cone slope must be positive"));
                }
                let axis = scale(&axis);
                let growth = &Real::one() + &slope.square();
                let root_growth = if cone { root(&growth) } else { Real::one() };
                let offset = &rim + &(&(&b * &radius) * &root_growth);
                // A ball exactly as wide as the cylinder it rolls inside keeps
                // its centre on the axis (4.2a step 5): the offset cylinder is
                // the axis itself, so the path is one point and the face a
                // sphere. RECIPE 4.5.1: "A round exactly as wide as the ball is
                // not a cliff"; OpenCASCADE's `ChFiKPart_ComputeData_FilPlnCyl`
                // makes the same sphere when `cylrad - Radius` is zero. Only a
                // plane normal to the axis is taken, where every meridian's
                // section is the same quarter disc.
                let normal_axis = axis.cross(&unit);
                let on_axis = if cone {
                    false
                } else {
                    match sign(&offset, dec)? {
                        Sign::Positive => false,
                        Sign::Zero if sign(&normal_axis.norm_squared(), dec)? == Sign::Zero => true,
                        _ => {
                            return Err(reach(
                                "the ball is no smaller than the cylinder it rolls inside",
                            ));
                        }
                    }
                };
                let cross = unit.dot(&axis);
                let facing = sign(&cross, dec)?;
                if facing == Sign::Zero {
                    return Err(reach(
                        "a plane parallel to the axis meets the offset in generators, not a conic",
                    ));
                }
                if cone && !positive(&(&(&cross.square() * &growth) - &slope.square()), dec)? {
                    return Err(reach(
                        "the offset plane meets the offset cone in an open conic",
                    ));
                }
                let lift = unit.dot(&position(&origin));
                let start = &(&level + &radius) - &lift;
                let gain = &(&offset * &cross) + &(&slope * &start);
                let edge_start = &level - &lift;
                let edge_gain = &(&rim * &cross) + &(&slope * &edge_start);
                // On the axis the gain is zero: the centre circle has radius 0.
                if (sign(&gain, dec)? != facing && !on_axis) || sign(&edge_gain, dec)? != facing {
                    return Err(reach(
                        "the edge or the offset conic leaves the support's nappe",
                    ));
                }
                let (basis, tilt, circle) = if sign(&normal_axis.norm_squared(), dec)? == Sign::Zero
                {
                    let major = if on_axis {
                        Real::zero()
                    } else {
                        &gain / &cross
                    };
                    // The foot on a cone lies `b r / sqrt(1 + k^2)` nearer the
                    // axis than the centre; past the axis it is on the other
                    // nappe, beyond the apex, and no ball touches the face.
                    if cone
                        && !positive(
                            &(&magnitude(&major, dec)? - &(&(&b * &radius) / &root_growth)),
                            dec,
                        )?
                    {
                        return Err(reach("the ball's foot on the cone lies past its apex"));
                    }
                    let ring = positive(&(&major - &radius), dec)?;
                    let circle = Circle {
                        centre: origin.plus(&axis.scaled(&(&start / &cross))),
                        axis: axis.clone(),
                        major_sq: major.square(),
                        major,
                        ring,
                        axial: on_axis,
                    };
                    (basis_of(&axis, dec)?, Real::zero(), Some(circle))
                } else {
                    let u2 = scale(&normal_axis);
                    let u1 = u2.cross(&axis);
                    let tilt = unit.dot(&u1);
                    // The ellipse's smallest radius |G| / max |D| must exceed r,
                    // so every ball stays off the axis and the meridian of a
                    // region point is defined.
                    let widest = &magnitude(&cross, dec)? + &(&slope * &magnitude(&tilt, dec)?);
                    if !positive(&(&magnitude(&gain, dec)? - &(&radius * &widest)), dec)? {
                        return Err(reach(
                            "the offset ellipse comes within the ball radius of its axis",
                        ));
                    }
                    ([u1, u2], tilt, None)
                };
                (
                    Shape::Linear {
                        origin,
                        axis,
                        rim,
                        slope,
                        cone,
                        root: root_growth,
                        offset,
                        start,
                        gain,
                        edge_start,
                        edge_gain,
                        cross,
                        tilt,
                    },
                    basis,
                    circle,
                )
            }
        };
        Ok(Self {
            normal: plane.normal,
            point: plane.point,
            unit,
            level,
            radius,
            outside,
            shape,
            basis,
            circle,
            sector: None,
        })
    }

    /// Restricts a circle path to the sector from the radial direction
    /// `start` to `end`, turning positively about `turn` (any vector along
    /// the axis), at most half a turn. Both directions must be unit and
    /// perpendicular to the axis. The sector's half-planes are junctions:
    /// each meets a tangent neighbour piece there.
    pub fn with_sector(
        mut self,
        start: Vector3,
        end: Vector3,
        turn: &Vector3,
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        let circle = self
            .circle
            .as_ref()
            .ok_or_else(|| reach("only a circle path takes a sector"))?;
        for direction in [&start, &end] {
            if sign(&(&direction.norm_squared() - &Real::one()), dec)? != Sign::Zero
                || sign(&direction.dot(&circle.axis), dec)? != Sign::Zero
            {
                return Err(reach(
                    "a sector end must be a unit direction perpendicular to the axis",
                ));
            }
        }
        // The sector turns from start to end through less than half a turn,
        // or exactly half when the ends are opposite.
        let wound = sign(&turn.dot(&start.cross(&end)), dec)?;
        let opposite = sign(&start.dot(&end), dec)? == Sign::Negative;
        if wound == Sign::Negative || (wound == Sign::Zero && !opposite) {
            return Err(reach("a ring sector over half a turn has no convex fan"));
        }
        let normals = [turn.cross(&start), end.cross(turn)];
        self.sector = Some(Sector {
            centre: circle.centre.clone(),
            ends: [start, end],
            normals,
        });
        Ok(self)
    }

    /// The sector's two junction half-planes, each as a point and a normal
    /// into the sector, when the path has one.
    pub fn sector_walls(&self) -> Option<[(Point3, Vector3); 2]> {
        let sector = self.sector.as_ref()?;
        Some(
            sector
                .normals
                .clone()
                .map(|normal| (sector.centre.clone(), normal)),
        )
    }

    /// Whether the path is a circle with closed forms.
    pub fn is_circle(&self) -> bool {
        self.circle.is_some()
    }

    /// The exact surface that holds the face, for the face registry: the ring
    /// torus about a centre circle wider than the ball, or the ball itself
    /// when the centre stays on the axis. Its equation has the separation's
    /// sign. A horn or spindle torus has no global equation, and an ellipse
    /// path's pipe is no elementary surface, so both have none.
    pub fn face(&self) -> Option<catallax_shapes::Surface> {
        let circle = self.circle.as_ref()?;
        if circle.ring {
            Some(catallax_shapes::Surface::Torus(catallax_shapes::Torus {
                center: circle.centre.clone(),
                axis: circle.axis.clone(),
                major: circle.major.clone(),
                minor: self.radius.clone(),
                ref_dir: self.basis[0].clone(),
            }))
        } else if circle.axial {
            Some(catallax_shapes::Surface::Sphere(catallax_shapes::Sphere {
                center: circle.centre.clone(),
                radius: self.radius.clone(),
                axis: circle.axis.clone(),
                ref_dir: self.basis[0].clone(),
            }))
        } else {
            None
        }
    }

    /// Ball radius in mm.
    pub fn radius(&self) -> &Real {
        &self.radius
    }

    /// Applies an exact rigid motion, reflections included. Dot products of
    /// the stored vectors and every scalar are invariant.
    pub fn transform(&mut self, pose: &Isometry) {
        let v = |v: &mut Vector3| *v = pose.apply_vector(v);
        let p = |p: &mut Point3| *p = pose.apply_point(p);
        v(&mut self.normal);
        p(&mut self.point);
        v(&mut self.unit);
        self.level = self.unit.dot(&position(&self.point));
        for e in &mut self.basis {
            v(e);
        }
        match &mut self.shape {
            Shape::Linear { origin, axis, .. } => {
                p(origin);
                v(axis);
            }
            Shape::Sphere { centre, .. } | Shape::Torus { centre, .. } => p(centre),
        }
        if let Some(circle) = &mut self.circle {
            p(&mut circle.centre);
            v(&mut circle.axis);
        }
        if let Some(sector) = &mut self.sector {
            p(&mut sector.centre);
            for e in sector.ends.iter_mut().chain(&mut sector.normals) {
                v(e);
            }
        }
    }

    fn frame_at(&self, direction: &Vector3, dec: &mut Decider) -> Result<Frame, Refusal> {
        let circle = self.circle.as_ref().expect("closed forms need a circle");
        let centre = circle.centre.plus(&direction.scaled(&circle.major));
        let plane = centre.plus(&self.unit.scaled(&-&self.radius));
        let b = side(self.outside);
        let wall = match &self.shape {
            Shape::Linear {
                axis, slope, root, ..
            } => centre.plus(
                &direction
                    .minus(&axis.scaled(slope))
                    .scaled(&-&(&(&b * &self.radius) / root)),
            ),
            Shape::Sphere {
                centre: o, offset, ..
            } => centre.plus(&centre.minus(o).scaled(&-&(&(&b * &self.radius) / offset))),
            // The foot lies on the line from the core circle's point in the
            // station's meridian through the station.
            Shape::Torus {
                centre: o,
                major: big,
                offset,
                ..
            } => {
                let core = o.plus(&direction.scaled(big));
                centre.plus(
                    &centre
                        .minus(&core)
                        .scaled(&-&(&(&b * &self.radius) / offset)),
                )
            }
        };
        // A torus's feet are opposite where the offset plane touches the
        // offset tube's top; the face arc is the half toward the edge.
        let side = match &self.shape {
            Shape::Torus {
                centre: o,
                edge,
                edge_major,
                ..
            } => {
                let point = o
                    .plus(&self.unit.scaled(edge))
                    .plus(&direction.scaled(edge_major));
                Some(balls(&point.minus(&centre), dec)?)
            }
            _ => None,
        };
        Ok(Frame {
            centre: point_balls(&centre, dec)?,
            tangent: balls(&circle.axis.cross(direction), dec)?,
            feet: [point_balls(&plane, dec)?, point_balls(&wall, dec)?],
            side,
        })
    }

    /// Enclosure of the station in the meridian of a point off the axis, an
    /// upper bound of the distance to an ellipse path. `None` for circles,
    /// whose nearest station is exact, and on the axis.
    pub(super) fn meridian_station(
        &self,
        local: &[Real; 3],
        dec: &mut Decider,
    ) -> Result<Option<[Ball; 3]>, Refusal> {
        let Shape::Linear {
            origin,
            axis,
            slope,
            offset,
            start,
            gain,
            cross,
            ..
        } = &self.shape
        else {
            return Ok(None);
        };
        if self.circle.is_some() {
            return Ok(None);
        }
        let v = vector(local).minus(&position(origin));
        let s = axis.dot(&v);
        let w = v.minus(&axis.scaled(&s));
        let rho2 = w.norm_squared();
        if sign(&rho2, dec)? == Sign::Zero {
            return Ok(None);
        }
        let rho = root(&rho2);
        let nr = &self.unit.dot(&w) / &rho;
        let d = cross + &(slope * &nr);
        let station = origin
            .plus(&axis.scaled(&(&(start - &(offset * &nr)) / &d)))
            .plus(&w.scaled(&(&(gain / &d) / &rho)));
        let precision = crate::metric::PRECISION;
        Ok(Some([
            enclose(&station.x, precision, dec)?,
            enclose(&station.y, precision, dec)?,
            enclose(&station.z, precision, dec)?,
        ]))
    }

    /// An axis-aligned box holding the region: the hull, over a cover of the
    /// conic by 64 panels of `pi / 32`, of the edge, both rails and the
    /// centre path, and, for a sphere, of the disc on each chord from the edge
    /// to the sphere's rail. The region's section in each meridian lies in the
    /// hull of its edge point, rail points and station, plus the segments
    /// along the plane normal between the centre path and the plane rail, and
    /// a sphere's arc in that disc. The panel width only sets how loose the
    /// box is.
    pub fn hull(&self, dec: &mut Decider) -> Result<[Ball; 3], Refusal> {
        let cover = self.cover(dec)?;
        let mut points = Vec::new();
        for u in self.panels(&cover, dec)? {
            points.extend(Self::panel_points(&cover, &u)?);
        }
        Ok(hull_points(&points))
    }

    /// The panels of [`panels`] that may hold a station: all of them, or for
    /// a sector those whose meridians are not proven outside one of its
    /// half-planes.
    fn panels(&self, cover: &RevolvedCover, dec: &mut Decider) -> Result<Vec<Ball>, Refusal> {
        let Some(sector) = &self.sector else {
            return Ok(panels().collect());
        };
        let normals = [
            balls(&sector.normals[0], dec)?,
            balls(&sector.normals[1], dec)?,
        ];
        let mut kept = Vec::new();
        for u in panels() {
            let (_, _, direction) = cover.angle(&u);
            let outside = normals.iter().any(|n| {
                (&(&(&n[0] * &direction[0]) + &(&n[1] * &direction[1])) + &(&n[2] * &direction[2]))
                    .hi()
                    < Rational::zero()
            });
            if !outside {
                kept.push(u);
            }
        }
        Ok(kept)
    }

    /// The box of the region over each of the 64 panels of [`Revolved::hull`].
    /// A box round a whole ring holds the axis, so a support near the ring's
    /// inside can reach it while staying clear of every panel's box.
    pub fn boxes(&self, dec: &mut Decider) -> Result<Vec<[Ball; 3]>, Refusal> {
        let cover = self.cover(dec)?;
        self.panels(&cover, dec)?
            .iter()
            .map(|u| Ok(hull_points(&Self::panel_points(&cover, u)?)))
            .collect()
    }

    fn panel_points(cover: &RevolvedCover, u: &Ball) -> Result<Vec<[Ball; 3]>, Refusal> {
        let mut points = Vec::new();
        let frame = cover.frame(u)?;
        let edge = cover.edge(u)?;
        if let Kind::Round { .. } = cover.kind {
            points.extend(chord_cube(&edge, &frame.feet[1])?);
        }
        points.push(edge);
        points.push(frame.centre);
        points.extend(frame.feet);
        Ok(points)
    }

    /// Tentative binary64 value for display: negative inside the region. The
    /// square root is the caller's.
    pub(crate) fn fast(&self, p: [f64; 3], sqrt: fn(f64) -> f64) -> f64 {
        let f = |x: &Real| x.to_f64();
        let fv = |v: &Vector3| [v.x.to_f64(), v.y.to_f64(), v.z.to_f64()];
        let fp = |p: &Point3| [p.x.to_f64(), p.y.to_f64(), p.z.to_f64()];
        let dot = |a: [f64; 3], b: [f64; 3]| a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
        let sub = |a: [f64; 3], b: [f64; 3]| [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
        let add =
            |a: [f64; 3], t: f64, b: [f64; 3]| [a[0] + t * b[0], a[1] + t * b[1], a[2] + t * b[2]];
        let norm = |a: [f64; 3]| sqrt(dot(a, a));
        let n = fv(&self.unit);
        let r = f(&self.radius);
        let b = if self.outside { 1.0 } else { -1.0 };
        let height = dot(n, p) - f(&self.level);
        let u1 = fv(&self.basis[0]);
        let (support, station, walls) = match &self.shape {
            Shape::Linear {
                origin,
                axis,
                rim,
                slope,
                root,
                offset,
                start,
                gain,
                cross,
                ..
            } => {
                let (o, a, k) = (fp(origin), fv(axis), f(slope));
                let v = sub(p, o);
                let s = dot(a, v);
                let w = add(v, -s, a);
                let rho = norm(w);
                let dir = if rho > 0.0 {
                    [w[0] / rho, w[1] / rho, w[2] / rho]
                } else {
                    u1
                };
                let nr = dot(n, dir);
                let d = f(cross) + k * nr;
                let rc = f(gain) / d;
                let sc = (f(start) - f(offset) * nr) / d;
                let station = add(add(o, sc, a), rc, dir);
                let support = b * (rho - f(rim) - k * s) / f(root);
                let lifted = add(p, -(height - r), n);
                let vp = sub(lifted, o);
                let sp = dot(a, vp);
                let face = norm(add(vp, -sp, a)) - (f(offset) + k * sp);
                let along = (s - sc) + k * (rho - rc);
                let walls = [-b * face, -along * d.signum()];
                (support, station, walls)
            }
            Shape::Sphere {
                centre,
                radius,
                offset,
                height: h0,
                ..
            } => {
                let o = fp(centre);
                let v = sub(p, o);
                let distance = norm(v);
                let h = dot(n, v);
                let circle = self.circle.as_ref().expect("a sphere path is a circle");
                let c = fp(&circle.centre);
                let w = add(sub(p, c), -dot(n, sub(p, c)), n);
                let rho = norm(w);
                let dir = if rho > 0.0 {
                    [w[0] / rho, w[1] / rho, w[2] / rho]
                } else {
                    u1
                };
                let station = add(c, f(&circle.major), dir);
                let support = b * (distance - f(radius));
                let lifted = add(p, -(height - r), n);
                let face = norm(sub(lifted, o)) - f(offset);
                let walls = [-b * face, f(h0) * distance - f(offset) * h];
                (support, station, walls)
            }
            Shape::Torus {
                centre,
                major: big,
                minor,
                offset,
                height: h0,
                facing,
                ..
            } => {
                let o = fp(centre);
                let h = dot(n, sub(p, o));
                let circle = self.circle.as_ref().expect("a torus path is a circle");
                let c = fp(&circle.centre);
                let w = add(sub(p, c), -dot(n, sub(p, c)), n);
                let rho = norm(w);
                let dir = if rho > 0.0 {
                    [w[0] / rho, w[1] / rho, w[2] / rho]
                } else {
                    u1
                };
                let station = add(c, f(&circle.major), dir);
                let spread = rho - f(big);
                let support = b * (sqrt(spread * spread + h * h) - f(minor));
                let face = sqrt(spread * spread + f(h0) * f(h0)) - f(offset);
                let meridian = f(facing) * (f(h0) * spread - (f(&circle.major) - f(big)) * h);
                (support, station, [-b * face, meridian])
            }
        };
        let m = norm(sub(p, station));
        let in_sector = self.sector.as_ref().is_none_or(|sector| {
            let v = sub(p, fp(&sector.centre));
            sector.normals.iter().all(|n| dot(fv(n), v) >= 0.0)
        });
        if in_sector && walls.iter().all(|w| *w >= 0.0) {
            (-height).max(-support).max(r - m)
        } else {
            (m - r).abs().max(f64::MIN_POSITIVE)
        }
    }
}

/// A cover of `[-1, 1]` by 64 panels of width `1 / 32`, for the boxes of a
/// closed path.
pub(super) fn panels() -> impl Iterator<Item = Ball> {
    const PANELS: i64 = 64;
    (0..PANELS).map(|i| {
        let lo = Rational::from_i64(2 * i - PANELS) / Rational::from_i64(PANELS);
        let hi = Rational::from_i64(2 * i + 2 - PANELS) / Rational::from_i64(PANELS);
        interval(&lo, &hi)
    })
}

/// Interval cover of a [`Revolved`] path; `u` in `[-1, 1]` is the angle `pi u`
/// from the first basis direction.
#[derive(Clone, Debug)]
pub struct RevolvedCover {
    unit: [Ball; 3],
    radius: Ball,
    pi: Ball,
    basis: [[Ball; 3]; 2],
    // b r / sqrt(1 + k^2) for a cylinder or cone, b r / R' for a sphere.
    foot: Ball,
    kind: Kind,
}
// One cover is built per query and never stored in bulk, so the size gap
// between the variants costs nothing worth a box.
#[allow(clippy::large_enum_variant)]
#[derive(Clone, Debug)]
enum Kind {
    Linear {
        axis: [Ball; 3],
        slope: Ball,
        cross: Ball,
        tilt: Ball,
        // Where the axis meets the offset plane and the plane.
        pivot: [Ball; 3],
        edge_pivot: [Ball; 3],
        gain: Ball,
        edge_gain: Ball,
    },
    Round {
        centre: [Ball; 3],
        circle: [Ball; 3],
        major: Ball,
        edge: Ball,
        edge_major: Ball,
        // A torus's core radius; `None` for a sphere.
        core: Option<Ball>,
    },
}
/// Cosine, sine and meridian `cos e1 + sin e2` over the panel of angles
/// `pi u`, by the mean-value form at its midpoint: with `|h| <= delta`,
/// `cos(t + h) = cos t - h sin t` and `sin(t + h) = sin t + h cos t`, each
/// within `delta^2 / 2`. A Taylor series over the whole panel loses much more
/// near `pi`.
pub(super) fn angle(pi: &Ball, basis: &[[Ball; 3]; 2], u: &Ball) -> (Ball, Ball, [Ball; 3]) {
    let centre = pi * &Ball::from_rational(&u.mid(), BITS);
    let (c, s) = (centre.cos(), centre.sin());
    let delta = (pi * &Ball::from_rational(&u.rad(), BITS)).hi();
    let h = interval(&-&delta, &delta);
    let rest = &delta * &delta / q(2);
    let rest = interval(&-&rest, &rest);
    let cos = &(&c - &(&s * &h)) + &rest;
    let sin = &(&s + &(&c * &h)) + &rest;
    let direction = sum(&scaled(&cos, &basis[0]), &scaled(&sin, &basis[1]));
    (cos, sin, direction)
}

impl RevolvedCover {
    fn angle(&self, u: &Ball) -> (Ball, Ball, [Ball; 3]) {
        angle(&self.pi, &self.basis, u)
    }
    /// `p + (G / D) (m - a cos (n.u1) / n.a)` for a pivot and gain.
    fn conic(&self, u: &Ball, pivot: &[Ball; 3], gain: &Ball) -> Result<[Ball; 3], Refusal> {
        let Kind::Linear {
            axis,
            slope,
            cross,
            tilt,
            ..
        } = &self.kind
        else {
            unreachable!("a conic of a cylinder or cone")
        };
        let (cos, _, direction) = self.angle(u);
        let nr = &cos * tilt;
        let rho = divide(gain, &(cross + &(slope * &nr)))?;
        let lean = divide(&nr, cross)?;
        let spoke: [Ball; 3] = std::array::from_fn(|i| &direction[i] - &(&lean * &axis[i]));
        Ok(sum(pivot, &scaled(&rho, &spoke)))
    }
    /// Edge point in the meridian of the panel.
    fn edge(&self, u: &Ball) -> Result<[Ball; 3], Refusal> {
        let (_, _, direction) = self.angle(u);
        Ok(match &self.kind {
            Kind::Linear {
                edge_pivot,
                edge_gain,
                ..
            } => self.conic(u, edge_pivot, edge_gain)?,
            Kind::Round {
                centre,
                edge,
                edge_major,
                ..
            } => sum(
                &sum(centre, &scaled(edge, &self.unit)),
                &scaled(edge_major, &direction),
            ),
        })
    }
}
impl Cover for RevolvedCover {
    fn radius(&self) -> &Ball {
        &self.radius
    }
    fn centre(&self, u: &Ball) -> Result<[Ball; 3], Refusal> {
        Ok(match &self.kind {
            Kind::Linear { pivot, gain, .. } => self.conic(u, pivot, gain)?,
            Kind::Round { circle, major, .. } => {
                let (_, _, direction) = self.angle(u);
                sum(circle, &scaled(major, &direction))
            }
        })
    }
    fn frame(&self, u: &Ball) -> Result<Frame, Refusal> {
        let (cos, sin, direction) = self.angle(u);
        let turn = sum(
            &scaled(&-&sin, &self.basis[0]),
            &scaled(&cos, &self.basis[1]),
        );
        let centre = self.centre(u)?;
        let plane = shift(&centre, &-&self.radius, &self.unit);
        let (tangent, wall) = match &self.kind {
            Kind::Linear {
                axis,
                slope,
                cross,
                tilt,
                ..
            } => {
                // c' times D^2 / G: -(n.t) (a + k m) + D t, with t the
                // meridian's derivative and n.t = -sin (n.u1).
                let d = cross + &(slope * &(&cos * tilt));
                let along = &sin * tilt;
                let tangent = sum(
                    &scaled(&along, &sum(axis, &scaled(slope, &direction))),
                    &scaled(&d, &turn),
                );
                let normal: [Ball; 3] =
                    std::array::from_fn(|i| &direction[i] - &(slope * &axis[i]));
                (tangent, shift(&centre, &-&self.foot, &normal))
            }
            Kind::Round {
                centre: o, core, ..
            } => {
                let radial: [Ball; 3] = match core {
                    None => std::array::from_fn(|i| &centre[i] - &o[i]),
                    Some(core) => {
                        let point = sum(o, &scaled(core, &direction));
                        std::array::from_fn(|i| &centre[i] - &point[i])
                    }
                };
                (turn, shift(&centre, &-&self.foot, &radial))
            }
        };
        let side = match &self.kind {
            Kind::Round { core: Some(_), .. } => {
                let edge = self.edge(u)?;
                Some(std::array::from_fn(|i| &edge[i] - &centre[i]))
            }
            _ => None,
        };
        Ok(Frame {
            centre,
            tangent,
            feet: [plane, wall],
            side,
        })
    }
}

impl Path for Revolved {
    type Cover = RevolvedCover;
    fn local(&self, point: &Point3) -> [Real; 3] {
        local(point)
    }
    fn walls(&self, local: &[Real; 3], dec: &mut Decider) -> Result<Vec<Real>, Refusal> {
        let x = vector(local);
        let two = Real::from_i64(2);
        let b = side(self.outside);
        let r = &self.radius;
        let height = &self.unit.dot(&x) - &self.level;
        let mut walls = vec![
            self.normal.dot(&x.minus(&position(&self.point))),
            &(&two * r) - &height,
        ];
        // Projection along the plane normal onto the offset plane.
        let lifted = x.minus(&self.unit.scaled(&(&height - r)));
        match &self.shape {
            Shape::Linear {
                origin,
                axis,
                rim,
                slope,
                cone,
                root: growth,
                offset,
                start,
                gain,
                cross,
                ..
            } => {
                let v = x.minus(&position(origin));
                let s = axis.dot(&v);
                let w = v.minus(&axis.scaled(&s));
                let rho2 = w.norm_squared();
                // The region keeps a positive distance from the axis (see
                // `Revolved::new`), so an axis point is outside it. With the
                // centre on the axis, the fan's axis points lie between the
                // plane and the centre, inside the one closed ball; only the
                // plane contact is on the region's boundary, where `-height`
                // is zero and the separation decides it.
                if sign(&rho2, dec)? == Sign::Zero {
                    let axial = self.circle.as_ref().is_some_and(|c| c.axial);
                    walls.push(if axial { -&height } else { Real::from_i64(-1) });
                    return Ok(walls);
                }
                let rho = root(&rho2);
                let vp = lifted.minus(&position(origin));
                let sp = axis.dot(&vp);
                let rho2p = vp.minus(&axis.scaled(&sp)).norm_squared();
                if *cone {
                    let gap = &rho - &(rim + &(slope * &s));
                    walls.push(&b * &gap);
                    if self.outside {
                        walls.push(&(&(&two * r) * growth) - &gap);
                    }
                    walls.push(-&(&b * &(&root(&rho2p) - &(offset + &(slope * &sp)))));
                } else {
                    walls.push(&b * &(&rho2 - &rim.square()));
                    if self.outside {
                        walls.push(&(rim + &(&two * r)).square() - &rho2);
                    }
                    walls.push(-&(&b * &(&rho2p - &offset.square())));
                }
                // The meridian's normal line through the station: with
                // L = (s - s_c) + k (rho - rho_c), L D rho = rho P1 + P2 and
                // the edge side is where -sign(D) L > 0, so -(rho P1 + P2) > 0.
                let nw = self.unit.dot(&w);
                let p1 = &(&(&s * cross) + &(&slope.square() * &nw)) - &(start + &(slope * gain));
                let p2 = &(&(&(&s * slope) * &nw) + &(&(slope * cross) * &rho2)) + &(offset * &nw);
                walls.push(-&(&(&rho * &p1) + &p2));
            }
            Shape::Sphere {
                centre,
                radius: big,
                offset,
                height: h0,
                ..
            } => {
                let v = x.minus(&position(centre));
                let d2 = v.norm_squared();
                walls.push(&b * &(&d2 - &big.square()));
                if self.outside {
                    walls.push(&(big + &(&two * r)).square() - &d2);
                }
                let vp = lifted.minus(&position(centre));
                walls.push(-&(&b * &(&vp.norm_squared() - &offset.square())));
                // The offset sphere's normals through the circle meet at the
                // centre; the edge side has h / |v| < h0 / R'.
                walls.push(&(h0 * &root(&d2)) - &(offset * &self.unit.dot(&v)));
            }
            Shape::Torus {
                centre,
                major: big,
                minor,
                offset,
                height: h0,
                facing,
                ..
            } => {
                // In the meridian half-plane the torus is the circle of radius
                // a about (R, 0), as a sphere's is about (0, 0).
                let v = x.minus(&position(centre));
                let h = self.unit.dot(&v);
                let spread = &root(&(&v.norm_squared() - &h.square())) - big;
                let tube = &spread.square() + &h.square();
                walls.push(&b * &(&tube - &minor.square()));
                if self.outside {
                    walls.push(&(minor + &(&two * r)).square() - &tube);
                }
                walls.push(-&(&b * &(&(&spread.square() + &h0.square()) - &offset.square())));
                // The offset tube's normal line through the station passes
                // through the core circle's point (R, 0) of the meridian.
                let circle = self.circle.as_ref().expect("a torus path is a circle");
                walls.push(facing * &(&(h0 * &spread) - &(&(&circle.major - big) * &h)));
            }
        }
        Ok(walls)
    }
    fn walls_are_faces(&self) -> bool {
        true
    }
    fn junctions(&self, local: &[Real; 3], _: &mut Decider) -> Result<Vec<Real>, Refusal> {
        let Some(sector) = &self.sector else {
            return Ok(Vec::new());
        };
        let v = vector(local).minus(&position(&sector.centre));
        Ok(sector.normals.iter().map(|n| n.dot(&v)).collect())
    }
    fn uniform(&self) -> bool {
        self.circle.is_some()
    }
    fn separation(&self, local: &[Real; 3]) -> Option<Real> {
        let circle = self.circle.as_ref()?;
        let v = vector(local).minus(&position(&circle.centre));
        let h = circle.axis.dot(&v);
        let d2 = v.norm_squared();
        let rho2 = &d2 - &h.square();
        let s = &(&d2 + &circle.major_sq) - &self.radius.square();
        Some(if circle.ring {
            &s.square() - &(&(&Real::from_i64(4) * &circle.major_sq) * &rho2)
        } else if circle.axial {
            // One ball: `|x - c|^2 - r^2`.
            s
        } else {
            &s - &(&(&Real::from_i64(2) * &circle.major) * &root(&rho2))
        })
    }
    fn cover(&self, dec: &mut Decider) -> Result<RevolvedCover, Refusal> {
        let ball = |x: &Real, dec: &mut Decider| enclose(x, BITS, dec);
        let b = side(self.outside);
        let (foot, kind) = match &self.shape {
            Shape::Linear {
                origin,
                axis,
                slope,
                root: growth,
                start,
                gain,
                edge_start,
                edge_gain,
                cross,
                tilt,
                ..
            } => (
                ball(&(&(&b * &self.radius) / growth), dec)?,
                Kind::Linear {
                    axis: balls(axis, dec)?,
                    slope: ball(slope, dec)?,
                    cross: ball(cross, dec)?,
                    tilt: ball(tilt, dec)?,
                    pivot: point_balls(&origin.plus(&axis.scaled(&(start / cross))), dec)?,
                    edge_pivot: point_balls(
                        &origin.plus(&axis.scaled(&(edge_start / cross))),
                        dec,
                    )?,
                    gain: ball(gain, dec)?,
                    edge_gain: ball(edge_gain, dec)?,
                },
            ),
            Shape::Sphere {
                centre,
                offset,
                edge,
                edge_major,
                ..
            } => {
                let circle = self.circle.as_ref().expect("a sphere path is a circle");
                (
                    ball(&(&(&b * &self.radius) / offset), dec)?,
                    Kind::Round {
                        centre: point_balls(centre, dec)?,
                        circle: point_balls(&circle.centre, dec)?,
                        major: ball(&circle.major, dec)?,
                        edge: ball(edge, dec)?,
                        edge_major: ball(edge_major, dec)?,
                        core: None,
                    },
                )
            }
            Shape::Torus {
                centre,
                major: big,
                offset,
                edge,
                edge_major,
                ..
            } => {
                let circle = self.circle.as_ref().expect("a torus path is a circle");
                (
                    ball(&(&(&b * &self.radius) / offset), dec)?,
                    Kind::Round {
                        centre: point_balls(centre, dec)?,
                        circle: point_balls(&circle.centre, dec)?,
                        major: ball(&circle.major, dec)?,
                        edge: ball(edge, dec)?,
                        edge_major: ball(edge_major, dec)?,
                        core: Some(ball(big, dec)?),
                    },
                )
            }
        };
        Ok(RevolvedCover {
            unit: balls(&self.unit, dec)?,
            radius: ball(&self.radius, dec)?,
            pi: Ball::pi_checked(BITS, dec)?,
            basis: [balls(&self.basis[0], dec)?, balls(&self.basis[1], dec)?],
            foot,
            kind,
        })
    }
    fn nearest(&self, local: &[Real; 3], dec: &mut Decider) -> Result<Option<Frame>, Refusal> {
        let Some(circle) = &self.circle else {
            return Ok(None);
        };
        let v = vector(local).minus(&position(&circle.centre));
        let w = v.minus(&circle.axis.scaled(&circle.axis.dot(&v)));
        let rho2 = w.norm_squared();
        // On the axis every meridian is nearest; take the first basis direction.
        let on_axis = sign(&rho2, dec)? == Sign::Zero;
        let direction = match &self.sector {
            // The distance from a point to a station's centre grows with the
            // angle between their meridians, up to half a turn, so outside
            // the sector the nearer end is nearest: the one whose radial
            // direction is closer to the point's.
            Some(sector) => {
                let mut inside = !on_axis;
                for n in &sector.normals {
                    inside = inside && sign(&n.dot(&w), dec)? != Sign::Negative;
                }
                if inside {
                    w.scaled(&(&Real::one() / &root(&rho2)))
                } else if on_axis
                    || sign(&w.dot(&sector.ends[0].minus(&sector.ends[1])), dec)? != Sign::Negative
                {
                    sector.ends[0].clone()
                } else {
                    sector.ends[1].clone()
                }
            }
            None if on_axis => self.basis[0].clone(),
            None => w.scaled(&(&Real::one() / &root(&rho2))),
        };
        self.frame_at(&direction, dec).map(Some)
    }
}

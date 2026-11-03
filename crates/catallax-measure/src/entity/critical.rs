//! Critical pairs between the exact cores of two strata.
//!
//! A core is the set a carrier is a constant-distance offset of: a sphere's
//! centre, a cylinder's axis, a torus's core circle, or the carrier itself for
//! planes, cones, edges and vertices. Every pair of points at which the core
//! distance is stationary, including the pairs where the cores meet, is
//! returned. One-parameter sets of such pairs are returned as families.

use super::{root, sign, unit, zero};
use catallax_numbers::{Decider, Real, Refusal, Sign};
use catallax_shapes::{Point3, Vector3};

/// An exact core. Directions are unit vectors; `reference` is perpendicular.
#[derive(Clone, Debug)]
pub(super) enum Core {
    Point(Point3),
    Line {
        origin: Point3,
        dir: Vector3,
    },
    Plane {
        origin: Point3,
        normal: Vector3,
    },
    Circle {
        centre: Point3,
        normal: Vector3,
        radius: Real,
        reference: Vector3,
    },
    /// The positive nappe. The apex is a singular point, never a smooth critical point.
    Cone {
        apex: Point3,
        axis: Vector3,
        slope: Real,
    },
}

impl Core {
    fn rank(&self) -> u8 {
        match self {
            Core::Point(_) => 0,
            Core::Line { .. } => 1,
            Core::Plane { .. } => 2,
            Core::Circle { .. } => 3,
            Core::Cone { .. } => 4,
        }
    }
    /// Whether this core is a cone, whose distance has no offset lower bound.
    pub(super) fn is_cone(&self) -> bool {
        matches!(self, Core::Cone { .. })
    }
}

/// Unit directions from the first core point toward the second along which
/// both carriers' offset points lie.
#[derive(Clone, Debug)]
pub(super) enum Dirs {
    /// The pair is apart; the direction is `(b - a) / dist`.
    Along,
    /// The cores meet; each listed unit vector is a common normal.
    Fixed(Vec<Vector3>),
    /// The cores meet and every unit direction is a common normal.
    All,
}

/// One stationary pair of core points.
#[derive(Clone, Debug)]
pub(super) struct CorePair {
    pub(super) a: Point3,
    pub(super) b: Point3,
    /// Exact nonnegative core distance `|b - a|`.
    pub(super) dist: Real,
    pub(super) dirs: Dirs,
}

/// A one- or two-parameter set of stationary pairs with one core distance:
/// `pair` translated along `slide` and turned about the line through `turn.0`
/// along `turn.1`.
#[derive(Clone, Debug)]
pub(super) struct Family {
    pub(super) pair: CorePair,
    pub(super) slide: Option<Vector3>,
    pub(super) turn: Option<(Point3, Vector3)>,
}

/// Every stationary pair of two cores.
#[derive(Clone, Debug, Default)]
pub(super) struct Critical {
    pub(super) pairs: Vec<CorePair>,
    pub(super) families: Vec<Family>,
}

/// Normal spaces of a core at one of its points.
#[derive(Clone, Debug)]
enum Normals {
    All,
    /// All vectors perpendicular to this unit vector.
    Perp(Vector3),
    /// Multiples of this unit vector.
    Span(Vector3),
}

fn parallel(a: &Vector3, b: &Vector3, dec: &mut Decider) -> Result<bool, Refusal> {
    zero(&a.cross(b).norm_squared(), dec)
}

/// The unit common normals of two meeting cores.
enum Meet {
    Fixed(Vec<Vector3>),
    Ring(Vector3),
    All,
}

fn meet(a: Normals, b: Normals, dec: &mut Decider) -> Result<Meet, Refusal> {
    let pm = |v: Vector3| Meet::Fixed(vec![v.clone(), v.scaled(&Real::from_i64(-1))]);
    Ok(match (a, b) {
        (Normals::All, Normals::All) => Meet::All,
        (Normals::All, Normals::Perp(v)) | (Normals::Perp(v), Normals::All) => Meet::Ring(v),
        (Normals::All, Normals::Span(v)) | (Normals::Span(v), Normals::All) => pm(v),
        (Normals::Perp(v), Normals::Perp(w)) => {
            if parallel(&v, &w, dec)? {
                Meet::Ring(v)
            } else {
                pm(unit(&v.cross(&w)))
            }
        }
        (Normals::Perp(v), Normals::Span(w)) | (Normals::Span(w), Normals::Perp(v)) => {
            if zero(&v.dot(&w), dec)? {
                pm(w)
            } else {
                Meet::Fixed(Vec::new())
            }
        }
        (Normals::Span(v), Normals::Span(w)) => {
            if parallel(&v, &w, dec)? {
                pm(v)
            } else {
                Meet::Fixed(Vec::new())
            }
        }
    })
}

/// A unit vector perpendicular to the unit vector `v`, chosen exactly.
pub(super) fn perpendicular(v: &Vector3, dec: &mut Decider) -> Result<Vector3, Refusal> {
    let one = Real::one;
    let zero_r = Real::zero;
    for e in [
        Vector3::new(one(), zero_r(), zero_r()),
        Vector3::new(zero_r(), one(), zero_r()),
    ] {
        let c = v.cross(&e);
        if !zero(&c.norm_squared(), dec)? {
            return Ok(unit(&c));
        }
    }
    unreachable!("a unit vector is parallel to at most one coordinate axis")
}

struct Out {
    critical: Critical,
}

impl Out {
    /// Adds a pair of core points. Meeting points take their common normals;
    /// a ring of normals becomes a turned family about that normal axis.
    fn pair(
        &mut self,
        a: Point3,
        b: Point3,
        normals: impl FnOnce(&mut Decider) -> Result<(Normals, Normals), Refusal>,
        dec: &mut Decider,
    ) -> Result<(), Refusal> {
        let squared = b.minus(&a).norm_squared();
        if !zero(&squared, dec)? {
            let dist = root(&squared);
            self.critical.pairs.push(CorePair {
                a,
                b,
                dist,
                dirs: Dirs::Along,
            });
            return Ok(());
        }
        let (na, nb) = normals(dec)?;
        match meet(na, nb, dec)? {
            Meet::Fixed(list) => self.critical.pairs.push(CorePair {
                a,
                b,
                dist: Real::zero(),
                dirs: Dirs::Fixed(list),
            }),
            Meet::All => self.critical.pairs.push(CorePair {
                a,
                b,
                dist: Real::zero(),
                dirs: Dirs::All,
            }),
            Meet::Ring(axis) => {
                let e = perpendicular(&axis, dec)?;
                self.critical.families.push(Family {
                    pair: CorePair {
                        a: a.clone(),
                        b,
                        dist: Real::zero(),
                        dirs: Dirs::Fixed(vec![e]),
                    },
                    slide: None,
                    turn: Some((a, axis)),
                });
            }
        }
        Ok(())
    }

    /// Adds a family of apart pairs, or of meeting pairs with fixed normals.
    fn family(
        &mut self,
        pair: (Point3, Point3),
        meeting: Vec<Vector3>,
        slide: Option<Vector3>,
        turn: Option<(Point3, Vector3)>,
        dec: &mut Decider,
    ) -> Result<(), Refusal> {
        let squared = pair.1.minus(&pair.0).norm_squared();
        let (dist, dirs) = if zero(&squared, dec)? {
            (Real::zero(), Dirs::Fixed(meeting))
        } else {
            (root(&squared), Dirs::Along)
        };
        self.critical.families.push(Family {
            pair: CorePair {
                a: pair.0,
                b: pair.1,
                dist,
                dirs,
            },
            slide,
            turn,
        });
        Ok(())
    }
}

fn pm(v: &Vector3) -> Vec<Vector3> {
    vec![v.clone(), v.scaled(&Real::from_i64(-1))]
}

/// Every stationary pair of `x` and `y`, or `None` when no closed form is
/// implemented for this pair of cores.
pub(super) fn critical(x: &Core, y: &Core, dec: &mut Decider) -> Result<Option<Critical>, Refusal> {
    if x.rank() > y.rank() {
        return Ok(critical(y, x, dec)?.map(|c| Critical {
            pairs: c.pairs.into_iter().map(swap).collect(),
            families: c
                .families
                .into_iter()
                .map(|f| Family {
                    pair: swap(f.pair),
                    ..f
                })
                .collect(),
        }));
    }
    let mut out = Out {
        critical: Critical::default(),
    };
    let supported = ordered(x, y, &mut out, dec)?;
    Ok(supported.then_some(out.critical))
}

fn swap(p: CorePair) -> CorePair {
    let dirs = match p.dirs {
        Dirs::Fixed(list) => Dirs::Fixed(
            list.into_iter()
                .map(|v| v.scaled(&Real::from_i64(-1)))
                .collect(),
        ),
        other => other,
    };
    CorePair {
        a: p.b,
        b: p.a,
        dist: p.dist,
        dirs,
    }
}

fn foot_on_line(p: &Point3, origin: &Point3, dir: &Vector3) -> Point3 {
    origin.plus(&dir.scaled(&p.minus(origin).dot(dir)))
}

fn positive(value: &Real, dec: &mut Decider) -> Result<bool, Refusal> {
    Ok(sign(value, dec)? == Sign::Positive)
}

#[allow(clippy::too_many_lines)]
fn ordered(x: &Core, y: &Core, out: &mut Out, dec: &mut Decider) -> Result<bool, Refusal> {
    match (x, y) {
        (Core::Point(a), Core::Point(b)) => {
            out.pair(
                a.clone(),
                b.clone(),
                |_| Ok((Normals::All, Normals::All)),
                dec,
            )?;
        }
        (Core::Point(a), Core::Line { origin, dir }) => {
            let f = foot_on_line(a, origin, dir);
            out.pair(
                a.clone(),
                f,
                |_| Ok((Normals::All, Normals::Perp(dir.clone()))),
                dec,
            )?;
        }
        (Core::Point(a), Core::Plane { origin, normal }) => {
            let h = a.minus(origin).dot(normal);
            let f = a.plus(&normal.scaled(&-h));
            out.pair(
                a.clone(),
                f,
                |_| Ok((Normals::All, Normals::Span(normal.clone()))),
                dec,
            )?;
        }
        (
            Core::Point(a),
            Core::Circle {
                centre,
                normal,
                radius,
                reference,
            },
        ) => {
            let w = a.minus(centre);
            let planar = w.minus(&normal.scaled(&w.dot(normal)));
            if zero(&planar.norm_squared(), dec)? {
                let x = centre.plus(&reference.scaled(radius));
                out.family(
                    (a.clone(), x),
                    Vec::new(),
                    None,
                    Some((centre.clone(), normal.clone())),
                    dec,
                )?;
            } else {
                let e = unit(&planar);
                for s in [1, -1] {
                    let x = centre.plus(&e.scaled(&(radius * &Real::from_i64(s))));
                    let t = normal.cross(&e);
                    out.pair(a.clone(), x, |_| Ok((Normals::All, Normals::Perp(t))), dec)?;
                }
            }
        }
        (Core::Point(a), Core::Cone { apex, axis, slope }) => {
            let w = a.minus(apex);
            let planar = w.minus(&axis.scaled(&w.dot(axis)));
            let on_axis = zero(&planar.norm_squared(), dec)?;
            let radials = if on_axis {
                vec![perpendicular(axis, dec)?]
            } else {
                let e = unit(&planar);
                vec![e.clone(), e.scaled(&Real::from_i64(-1))]
            };
            for e in radials {
                let g = axis.plus(&e.scaled(slope));
                let t = w.dot(&g) / g.norm_squared();
                if !positive(&t, dec)? {
                    continue;
                }
                let foot = apex.plus(&g.scaled(&t));
                let normal = unit(&e.minus(&axis.scaled(slope)));
                if on_axis {
                    out.family(
                        (a.clone(), foot),
                        Vec::new(),
                        None,
                        Some((apex.clone(), axis.clone())),
                        dec,
                    )?;
                } else {
                    out.pair(
                        a.clone(),
                        foot,
                        |_| Ok((Normals::All, Normals::Span(normal))),
                        dec,
                    )?;
                }
            }
        }
        (
            Core::Line {
                origin: o1,
                dir: d1,
            },
            Core::Line {
                origin: o2,
                dir: d2,
            },
        ) => {
            let c = d1.cross(d2);
            if zero(&c.norm_squared(), dec)? {
                let f = foot_on_line(o1, o2, d2);
                let squared = f.minus(o1).norm_squared();
                if zero(&squared, dec)? {
                    // Coincident axes: every axial position and normal direction.
                    let e = perpendicular(d1, dec)?;
                    out.critical.families.push(Family {
                        pair: CorePair {
                            a: o1.clone(),
                            b: o1.clone(),
                            dist: Real::zero(),
                            dirs: Dirs::Fixed(vec![e]),
                        },
                        slide: Some(d1.clone()),
                        turn: Some((o1.clone(), d1.clone())),
                    });
                } else {
                    out.family((o1.clone(), f), Vec::new(), Some(d1.clone()), None, dec)?;
                }
            } else {
                let w = o1.minus(o2);
                let b = d1.dot(d2);
                let den = Real::one() - b.square();
                let (dw1, dw2) = (d1.dot(&w), d2.dot(&w));
                let s = (&(&b * &dw2) - &dw1) / den.clone();
                let t = (&dw2 - &(&b * &dw1)) / den;
                let k1 = o1.plus(&d1.scaled(&s));
                let k2 = o2.plus(&d2.scaled(&t));
                let n = unit(&c);
                out.pair(
                    k1,
                    k2,
                    |_| Ok((Normals::Span(n.clone()), Normals::Span(n))),
                    dec,
                )?;
            }
        }
        (Core::Line { origin, dir }, Core::Plane { origin: o, normal }) => {
            let dn = dir.dot(normal);
            let h = origin.minus(o).dot(normal);
            if zero(&dn, dec)? {
                let f = origin.plus(&normal.scaled(&-h));
                out.family(
                    (origin.clone(), f),
                    pm(normal),
                    Some(dir.clone()),
                    None,
                    dec,
                )?;
            } else {
                let k = origin.plus(&dir.scaled(&(-h / dn)));
                out.pair(
                    k.clone(),
                    k,
                    |_| Ok((Normals::Perp(dir.clone()), Normals::Span(normal.clone()))),
                    dec,
                )?;
            }
        }
        (
            Core::Line { origin, dir },
            Core::Circle {
                centre,
                normal,
                radius,
                reference,
            },
        ) => {
            if parallel(dir, normal, dec)? {
                let f = foot_on_line(centre, origin, dir);
                let planar = f.minus(centre);
                if zero(&planar.norm_squared(), dec)? {
                    let x = centre.plus(&reference.scaled(radius));
                    out.family(
                        (f, x),
                        Vec::new(),
                        None,
                        Some((centre.clone(), normal.clone())),
                        dec,
                    )?;
                } else {
                    let e = unit(&planar);
                    for s in [1, -1] {
                        let x = centre.plus(&e.scaled(&(radius * &Real::from_i64(s))));
                        let t = normal.cross(&e);
                        out.pair(
                            f.clone(),
                            x,
                            |_| Ok((Normals::Perp(dir.clone()), Normals::Perp(t))),
                            dec,
                        )?;
                    }
                }
            } else if zero(&dir.dot(normal), dec)? && zero(&origin.minus(centre).dot(normal), dec)?
            {
                // The line lies in the circle's plane.
                let k = foot_on_line(centre, origin, dir);
                let planar = k.minus(centre);
                let h2 = planar.norm_squared();
                let e = if zero(&h2, dec)? {
                    normal.cross(dir)
                } else {
                    unit(&planar)
                };
                for s in [1, -1] {
                    let x = centre.plus(&e.scaled(&(radius * &Real::from_i64(s))));
                    let t = normal.cross(&e);
                    out.pair(
                        k.clone(),
                        x,
                        |_| Ok((Normals::Perp(dir.clone()), Normals::Perp(t))),
                        dec,
                    )?;
                }
                let m2 = radius.square() - h2;
                if positive(&m2, dec)? {
                    let m = root(&m2);
                    for s in [1, -1] {
                        let x = k.plus(&dir.scaled(&(&m * &Real::from_i64(s))));
                        let t = normal.cross(&x.minus(centre));
                        out.pair(
                            x.clone(),
                            x,
                            |_| Ok((Normals::Perp(dir.clone()), Normals::Perp(unit(&t)))),
                            dec,
                        )?;
                    }
                }
            } else {
                return Ok(false);
            }
        }
        (
            Core::Plane {
                origin: o1,
                normal: n1,
            },
            Core::Plane {
                origin: o2,
                normal: n2,
            },
        ) => {
            let c = n1.cross(n2);
            if zero(&c.norm_squared(), dec)? {
                // Parallel planes form a two-parameter family.
                return Ok(false);
            }
            let b = n1.dot(n2);
            let (c1, c2) = (
                o1.minus(&Point3::origin()).dot(n1),
                o2.minus(&Point3::origin()).dot(n2),
            );
            let den = Real::one() - b.square();
            let alpha = (&c1 - &(&b * &c2)) / den.clone();
            let beta = (&c2 - &(&b * &c1)) / den;
            let k = Point3::origin()
                .plus(&n1.scaled(&alpha))
                .plus(&n2.scaled(&beta));
            out.family((k.clone(), k), Vec::new(), Some(unit(&c)), None, dec)?;
        }
        (
            Core::Plane { origin, normal: n },
            Core::Circle {
                centre,
                normal,
                radius,
                reference,
            },
        ) => {
            if parallel(n, normal, dec)? {
                let x = centre.plus(&reference.scaled(radius));
                let h = x.minus(origin).dot(n);
                let f = x.plus(&n.scaled(&-h));
                out.family(
                    (f, x),
                    pm(n),
                    None,
                    Some((centre.clone(), normal.clone())),
                    dec,
                )?;
            } else {
                let w = n.minus(&normal.scaled(&n.dot(normal)));
                let e = unit(&w);
                for s in [1, -1] {
                    let x = centre.plus(&e.scaled(&(radius * &Real::from_i64(s))));
                    let h = x.minus(origin).dot(n);
                    let f = x.plus(&n.scaled(&-h));
                    let t = normal.cross(&e);
                    out.pair(
                        f,
                        x,
                        |_| Ok((Normals::Span(n.clone()), Normals::Perp(t))),
                        dec,
                    )?;
                }
                // Crossing points: height h_c + radius |w| cos(phi) = 0.
                let hc = centre.minus(origin).dot(n);
                let wn = root(&w.norm_squared());
                let cos = -hc / (radius * &wn);
                let sin2 = Real::one() - cos.square();
                if positive(&sin2, dec)? {
                    let sin = root(&sin2);
                    let e2 = normal.cross(&e);
                    for s in [1, -1] {
                        let x = centre.plus(
                            &e.scaled(&(radius * &cos))
                                .plus(&e2.scaled(&(&(radius * &sin) * &Real::from_i64(s)))),
                        );
                        let t = unit(&normal.cross(&x.minus(centre)));
                        out.pair(
                            x.clone(),
                            x,
                            |_| Ok((Normals::Span(n.clone()), Normals::Perp(t))),
                            dec,
                        )?;
                    }
                }
            }
        }
        (Core::Plane { origin, normal: n }, Core::Cone { apex, axis, slope }) => {
            // A cone normal at azimuth e is e - slope*axis. It is parallel to n
            // only along one generator, which then runs parallel to the plane.
            let nz = n.dot(axis);
            let planar = n.minus(&axis.scaled(&nz));
            let p2 = planar.norm_squared();
            if zero(&p2, dec)? || !zero(&(nz.square() - slope.square() * p2), dec)? {
                return Ok(true);
            }
            let e = match sign(&nz, dec)? {
                Sign::Positive => unit(&planar).scaled(&Real::from_i64(-1)),
                _ => unit(&planar),
            };
            let g = unit(&axis.plus(&e.scaled(slope)));
            let h = apex.minus(origin).dot(n);
            let f = apex.plus(&n.scaled(&-h));
            out.family((f, apex.clone()), pm(n), Some(g), None, dec)?;
        }
        (
            Core::Circle {
                centre: c1,
                normal: n1,
                radius: r1,
                reference: ref1,
            },
            Core::Circle {
                centre: c2,
                normal: n2,
                radius: r2,
                ..
            },
        ) => {
            if !parallel(n1, n2, dec)? {
                return Ok(false);
            }
            let d = c2.minus(c1);
            let z = d.dot(n1);
            let planar = d.minus(&n1.scaled(&z));
            let p2 = planar.norm_squared();
            if zero(&p2, dec)? {
                for s in [1, -1] {
                    let a = c1.plus(&ref1.scaled(r1));
                    let b = c2.plus(&ref1.scaled(&(r2 * &Real::from_i64(s))));
                    let squared = b.minus(&a).norm_squared();
                    if zero(&squared, dec)? {
                        // The same circle: every point meets itself, with a ring
                        // of normals, so offset carriers coincide there.
                        out.critical.families.push(Family {
                            pair: CorePair {
                                a,
                                b,
                                dist: Real::zero(),
                                dirs: Dirs::All,
                            },
                            slide: None,
                            turn: Some((c1.clone(), n1.clone())),
                        });
                        continue;
                    }
                    out.family(
                        (a, b),
                        Vec::new(),
                        None,
                        Some((c1.clone(), n1.clone())),
                        dec,
                    )?;
                }
            } else {
                let e = unit(&planar);
                let t = n1.cross(&e);
                for s1 in [1, -1] {
                    for s2 in [1, -1] {
                        let a = c1.plus(&e.scaled(&(r1 * &Real::from_i64(s1))));
                        let b = c2.plus(&e.scaled(&(r2 * &Real::from_i64(s2))));
                        let t = t.clone();
                        out.pair(
                            a,
                            b,
                            |_| Ok((Normals::Perp(t.clone()), Normals::Perp(t))),
                            dec,
                        )?;
                    }
                }
                // Crossings of the projected circles are stationary at height |z|.
                let dp = root(&p2);
                let along = (&p2 + &(r1.square() - r2.square())) / (&dp * &Real::from_i64(2));
                let h2 = r1.square() - along.square();
                if positive(&h2, dec)? {
                    let m = root(&h2);
                    let e2 = n1.cross(&e);
                    for s in [1, -1] {
                        let a = c1
                            .plus(&e.scaled(&along))
                            .plus(&e2.scaled(&(&m * &Real::from_i64(s))));
                        let b = a.plus(&n1.scaled(&z));
                        let t1 = unit(&n1.cross(&a.minus(c1)));
                        let t2 = unit(&n1.cross(&b.minus(c2)));
                        out.pair(a, b, |_| Ok((Normals::Perp(t1), Normals::Perp(t2))), dec)?;
                    }
                }
            }
        }
        _ => return Ok(false),
    }
    Ok(true)
}

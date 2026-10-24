//! Convex conical cells with complete planar and coaxial circular trims.

use catallax_intersect::{Component, CurveForm, SurfacePair};
use catallax_numbers::{Decider, Real, Refusal, Sign};
use catallax_recipe::{FaceRef, Origin, Role, StepId};
use catallax_shapes::{
    Circle, Cone, Curve, Cylinder, Ellipse, GeomStore, Plane, Point3, Surface, Vector3,
};

use crate::profile::sign;

/// One named support of a bounded convex conical cell.
#[derive(Clone, Debug)]
pub struct ConeFace {
    /// Outward support, retaining the cone's positive nappe.
    pub surface: Surface,
    /// Physical generators of this face.
    pub origins: Vec<Origin>,
}

/// A complete shared circular or elliptic boundary.
#[derive(Clone, Debug)]
pub struct ConeBoundary {
    /// Native analytic meeting curve, oriented along the cone axis.
    pub curve: Curve,
    /// Faces below and above this section, respectively.
    pub faces: [usize; 2],
}

/// A convex cone, optionally limited by a coaxial cylinder and opposing planes.
/// Every retained section is closed and disjoint from the other sections.
#[derive(Clone, Debug)]
pub struct ConeSkin {
    /// Active physical faces only.
    pub faces: Vec<ConeFace>,
    /// Complete physical meeting curves.
    pub boundaries: Vec<ConeBoundary>,
    /// Singular cone tip and its incident face, when retained.
    pub apex: Option<(Point3, usize)>,
    /// Unit axial direction, towards increasing cone radius.
    pub axis: Vector3,
}

fn bad(detail: &str) -> Refusal {
    Refusal::new("distance.cone-domain", detail)
}
pub(crate) fn reduced(value: Real, dec: &mut Decider) -> Result<Real, Refusal> {
    if value.rational_ref().is_some() {
        return Ok(value);
    }
    Ok(dec.rational_value(&value)?.map(Real::from).unwrap_or(value))
}
fn unit(v: &Vector3, dec: &mut Decider) -> Result<Vector3, Refusal> {
    let length = reduced(v.norm_squared().sqrt(), dec)?;
    let v = v.scaled(&(&Real::one() / &length));
    Ok(Vector3::new(
        reduced(v.x, dec)?,
        reduced(v.y, dec)?,
        reduced(v.z, dec)?,
    ))
}
fn extrema(
    curve: &Curve,
    normal: &Vector3,
    origin: &Point3,
    dec: &mut Decider,
) -> Result<[Real; 2], Refusal> {
    let (center, u, v, a2, b2) = match curve {
        Curve::Circle(c) => {
            let u = unit(&c.ref_dir, dec)?;
            let v = unit(&c.normal, dec)?.cross(&u);
            (&c.center, u, v, &c.radius_sq, &c.radius_sq)
        }
        Curve::Ellipse(e) => {
            let u = unit(&e.major_dir, dec)?;
            let v = unit(&e.normal, dec)?.cross(&u);
            (&e.center, u, v, &e.major_sq, &e.minor_sq)
        }
        _ => unreachable!("closed conical sections are conics"),
    };
    let middle = center.minus(origin).dot(normal);
    let radius = reduced(
        (&(a2 * &normal.dot(&u).square()) + &(b2 * &normal.dot(&v).square())).sqrt(),
        dec,
    )?;
    Ok([
        reduced(&middle - &radius, dec)?,
        reduced(&middle + &radius, dec)?,
    ])
}
fn section(a: &Surface, b: &Surface, axis: &Vector3, dec: &mut Decider) -> Result<Curve, Refusal> {
    let mut geom = GeomStore::new();
    let a = geom.add_surface(a.clone());
    let b = geom.add_surface(b.clone());
    let SurfacePair::Meets(components) =
        catallax_intersect::intersect_surfaces(&mut geom, a, b, dec)?.value
    else {
        return Err(bad("the retained cap needs one closed conical section"));
    };
    let [
        Component::Curve {
            curve, singular, ..
        },
    ] = components.as_slice()
    else {
        return Err(bad(
            "split or singular conical sections need a trim arrangement",
        ));
    };
    if !singular.is_empty() {
        return Err(bad("singular conical cap sections need a trim arrangement"));
    }
    let form = catallax_intersect::curve_form(&geom, *curve, dec)?;
    Ok(match form {
        CurveForm::Circle {
            center,
            mut normal,
            radius_sq,
            ref_dir,
        } => {
            if sign(&normal.dot(axis), dec)? == Sign::Negative {
                normal = normal.scaled(&Real::from_i64(-1));
            }
            Curve::Circle(Circle {
                center,
                normal,
                radius_sq,
                ref_dir,
            })
        }
        CurveForm::Ellipse {
            center,
            mut normal,
            major_dir,
            major_sq,
            minor_sq,
        } => {
            if sign(&normal.dot(axis), dec)? == Sign::Negative {
                normal = normal.scaled(&Real::from_i64(-1));
            }
            Curve::Ellipse(Ellipse {
                center,
                normal,
                major_dir,
                major_sq,
                minor_sq,
            })
        }
        _ => return Err(bad("open cone-plane conics need a finite trim arrangement")),
    })
}

impl ConeSkin {
    /// Recognises outward convex supports. The caller separately proves that its
    /// trim inventory consists of complete sections, not selected curve arcs.
    pub fn new(
        surfaces: &[(Surface, Vec<Origin>)],
        dec: &mut Decider,
    ) -> Result<Option<Self>, Refusal> {
        let mut cone = None;
        let mut cylinder = None;
        let mut caps: [Option<(Plane, Vec<Origin>)>; 2] = [None, None];
        for (surface, origins) in surfaces {
            match surface {
                Surface::Cone(c) if cone.is_none() => cone = Some((c.clone(), origins.clone())),
                Surface::Cylinder(c) if cylinder.is_none() => {
                    cylinder = Some((c.clone(), origins.clone()))
                }
                Surface::Plane(_) => (),
                _ => return Ok(None),
            }
        }
        let Some((cone, names)) = cone else {
            return Ok(None);
        };
        let axis = unit(&cone.axis, dec)?;
        if let Some((c, _)) = &cylinder {
            if sign(&c.axis.cross(&axis).norm_squared(), dec)? != Sign::Zero
                || sign(&c.origin.minus(&cone.apex).cross(&axis).norm_squared(), dec)? != Sign::Zero
            {
                return Err(bad(
                    "noncoaxial cone-cylinder joins need a traced trim arrangement",
                ));
            }
        }
        let mut ignored = Vec::new();
        for (surface, origins) in surfaces {
            let Surface::Plane(p) = surface else { continue };
            let rate = p.normal.dot(&axis);
            let slot = match sign(&rate, dec)? {
                Sign::Negative => 0,
                Sign::Positive => 1,
                Sign::Zero => {
                    ignored.push(p.clone());
                    continue;
                }
            };
            let lateral = p.normal.minus(&axis.scaled(&rate));
            if sign(
                &(&rate.square() - &(&cone.slope.square() * &lateral.norm_squared())),
                dec,
            )? != Sign::Positive
            {
                ignored.push(p.clone());
                continue;
            }
            let apex_side = sign(&cone.apex.minus(&p.origin).dot(&p.normal), dec)?;
            if slot == 0 && apex_side != Sign::Positive {
                continue;
            }
            if slot == 1 && apex_side != Sign::Negative {
                return Err(Refusal::new(
                    "distance.empty-cone",
                    "a closing constraint consumes the cone nappe",
                ));
            }
            if let Some((previous, names)) = &mut caps[slot] {
                if sign(&previous.normal.cross(&p.normal).norm_squared(), dec)? != Sign::Zero {
                    let new = section(
                        &Surface::Cone(cone.clone()),
                        &Surface::Plane(p.clone()),
                        &axis,
                        dec,
                    )?;
                    if sign(
                        &extrema(&new, &previous.normal, &previous.origin, dec)?[1],
                        dec,
                    )? == Sign::Negative
                    {
                        *previous = p.clone();
                        *names = origins.clone();
                        continue;
                    }
                    let old = section(
                        &Surface::Cone(cone.clone()),
                        &Surface::Plane(previous.clone()),
                        &axis,
                        dec,
                    )?;
                    if sign(&extrema(&old, &p.normal, &p.origin, dec)?[1], dec)? == Sign::Negative {
                        continue;
                    }
                    return Err(bad(
                        "crossing nonparallel caps need a planar trim arrangement",
                    ));
                }
                match sign(&p.origin.minus(&previous.origin).dot(&previous.normal), dec)? {
                    Sign::Negative => {
                        *previous = p.clone();
                        *names = origins.clone();
                    }
                    Sign::Zero => {
                        names.extend(origins.iter().cloned());
                        names.sort();
                        names.dedup();
                    }
                    Sign::Positive => (),
                }
            } else {
                caps[slot] = Some((p.clone(), origins.clone()));
            }
        }
        if caps[1].is_none() {
            return Err(bad(
                "a conical cell needs a closing cap towards increasing radius",
            ));
        }
        let skin = Self::construct(cone, names, cylinder, caps, dec)?;
        for p in ignored {
            for b in &skin.boundaries {
                if sign(&extrema(&b.curve, &p.normal, &p.origin, dec)?[1], dec)? == Sign::Positive {
                    return Err(bad(
                        "an active open-conic plane needs a finite trim arrangement",
                    ));
                }
            }
            if let Some((a, _)) = &skin.apex
                && sign(&a.minus(&p.origin).dot(&p.normal), dec)? == Sign::Positive
            {
                return Err(bad("an axial plane clips the retained cone tip"));
            }
        }
        Ok(Some(skin))
    }

    fn construct(
        cone: Cone,
        names: Vec<Origin>,
        cylinder: Option<(Cylinder, Vec<Origin>)>,
        mut caps: [Option<(Plane, Vec<Origin>)>; 2],
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        let axis = unit(&cone.axis, dec)?;
        let mut faces = vec![ConeFace {
            surface: Surface::Cone(cone.clone()),
            origins: names,
        }];
        if let Some((c, names)) = &cylinder {
            faces.push(ConeFace {
                surface: Surface::Cylinder(c.clone()),
                origins: names.clone(),
            });
        }
        let join = cylinder.as_ref().map(|(c, _)| &c.radius() / &cone.slope);
        let cone_surface = Surface::Cone(cone.clone());
        let mut rings: [Option<(Curve, usize, usize)>; 2] = [None, None];
        for i in 0..2 {
            let Some((plane, origins)) = caps[i].take() else {
                continue;
            };
            let apex_side = sign(&cone.apex.minus(&plane.origin).dot(&plane.normal), dec)?;
            if i == 0 && apex_side != Sign::Positive {
                continue;
            }
            if i == 1 && apex_side != Sign::Negative {
                return Err(Refusal::new(
                    "distance.empty-cone",
                    "the closing plane consumes the cone nappe",
                ));
            }
            let conic = section(&cone_surface, &Surface::Plane(plane.clone()), &axis, dec)?;
            let wall = if let Some(join) = &join {
                let range = extrema(&conic, &axis, &cone.apex, dec)?;
                let low = sign(&(&range[0] - join), dec)?;
                let high = sign(&(&range[1] - join), dec)?;
                match (low, high) {
                    (_, Sign::Negative) => 0,
                    (Sign::Positive, _) => 1,
                    (Sign::Zero, Sign::Zero) => usize::from(i == 0),
                    _ => {
                        return Err(bad(
                            "a cap crossing or touching the cone-cylinder circle needs split conic trims",
                        ));
                    }
                }
            } else {
                0
            };
            let curve = if wall == 0 {
                conic
            } else {
                section(
                    &faces[wall].surface,
                    &Surface::Plane(plane.clone()),
                    &axis,
                    dec,
                )?
            };
            let cap = faces.len();
            faces.push(ConeFace {
                surface: Surface::Plane(plane),
                origins,
            });
            rings[i] = Some((curve, wall, cap));
        }
        if let (Some((lower, _, lower_id)), Some((upper, _, upper_id))) = (&rings[0], &rings[1]) {
            let Surface::Plane(top) = &faces[*upper_id].surface else {
                unreachable!()
            };
            let Surface::Plane(bottom) = &faces[*lower_id].surface else {
                unreachable!()
            };
            let upper_in_bottom = extrema(upper, &bottom.normal, &bottom.origin, dec)?;
            if sign(&upper_in_bottom[0], dec)? != Sign::Negative {
                return Err(Refusal::new(
                    "distance.empty-cone",
                    "opposing caps consume the conical cell",
                ));
            }
            if sign(&extrema(lower, &top.normal, &top.origin, dec)?[1], dec)? != Sign::Negative
                || sign(&upper_in_bottom[1], dec)? != Sign::Negative
            {
                return Err(bad("intersecting cap conics need a split trim arrangement"));
            }
        }
        let mut boundaries = Vec::new();
        if let Some((curve, wall, cap)) = &rings[0] {
            boundaries.push(ConeBoundary {
                curve: curve.clone(),
                faces: [*cap, *wall],
            });
        }
        let upper_wall = rings[1].as_ref().expect("closing cap").1;
        let lower_wall = rings[0].as_ref().map_or(0, |r| r.1);
        if lower_wall == 0 && upper_wall == 1 {
            boundaries.push(ConeBoundary {
                curve: section(&faces[0].surface, &faces[1].surface, &axis, dec)?,
                faces: [0, 1],
            });
        } else if lower_wall == 1 && upper_wall == 0 {
            return Err(bad("opposing caps reverse the cone-cylinder transition"));
        }
        let (curve, wall, cap) = rings[1].take().expect("closing cap");
        boundaries.push(ConeBoundary {
            curve,
            faces: [wall, cap],
        });
        let mut apex = rings[0].is_none().then_some((cone.apex, 0));
        let active: Vec<_> = (0..faces.len())
            .filter(|i| {
                boundaries.iter().any(|b| b.faces.contains(i))
                    || apex.as_ref().is_some_and(|a| a.1 == *i)
            })
            .collect();
        for b in &mut boundaries {
            b.faces = b.faces.map(|i| {
                active
                    .iter()
                    .position(|j| *j == i)
                    .expect("active boundary")
            });
        }
        if let Some((_, face)) = &mut apex {
            *face = active.iter().position(|i| i == face).expect("active apex");
        }
        Ok(Self {
            faces: active.into_iter().map(|i| faces[i].clone()).collect(),
            boundaries,
            apex: apex.take(),
            axis,
        })
    }

    /// Resolves active faces and their exact full-section adjacency.
    pub fn select(&self, refs: &[FaceRef]) -> Result<Vec<bool>, Refusal> {
        let mut selected = vec![false; self.faces.len()];
        for reference in refs {
            let candidates: Vec<_> = self
                .faces
                .iter()
                .enumerate()
                .filter(|(i, f)| {
                    f.origins.iter().any(|o| reference.face.matches(o))
                        && reference.touching.iter().all(|p| {
                            self.boundaries
                                .iter()
                                .filter(|b| b.faces.contains(i))
                                .any(|b| {
                                    b.faces.iter().any(|j| {
                                        j != i
                                            && self.faces[*j].origins.iter().any(|o| p.matches(o))
                                    })
                                })
                        })
                })
                .map(|(i, _)| i)
                .collect();
            if candidates.is_empty() {
                return Err(Refusal::new(
                    "reference.lost-face",
                    format!("no conical face resolves {reference:?}"),
                ));
            }
            if reference.near.is_some() && candidates.len() > 1 {
                return Err(bad(
                    "nearest conical trim selection needs a certified face-distance comparison",
                ));
            }
            for i in candidates {
                selected[i] = true;
            }
        }
        Ok(selected)
    }

    /// Intersects inward parallel convex supports. The result can consume a cap,
    /// a cylindrical band or the whole cell without extending a second nappe.
    pub fn eroded(&self, distance: &Real, dec: &mut Decider) -> Result<Option<Self>, Refusal> {
        let surfaces = self
            .faces
            .iter()
            .map(|f| {
                Ok(f.surface
                    .parallel_support(&(-distance), dec)?
                    .map(|s| (s, f.origins.clone())))
            })
            .collect::<Result<Option<Vec<_>>, Refusal>>()?;
        let Some(surfaces) = surfaces else {
            return Ok(None);
        };
        match Self::new(&surfaces, dec) {
            Err(e) if e.code == "distance.empty-cone" => Ok(None),
            other => other,
        }
    }

    /// Proves kept-face clearance witnesses for a closed cell or a perpendicular
    /// opening towards increasing radius. Kept cone projections decrease axial
    /// height; cap and cylinder projections cannot cross that open halfspace.
    pub fn cavity(
        &self,
        open: &[bool],
        thickness: &Real,
        step: &StepId,
        dec: &mut Decider,
    ) -> Result<Option<Self>, Refusal> {
        if open.len() != self.faces.len() {
            return Err(bad("selection and conical face inventory differ"));
        }
        for (f, opened) in self.faces.iter().zip(open) {
            if !opened {
                continue;
            }
            let Surface::Plane(p) = &f.surface else {
                return Err(bad(
                    "opening a conical wall needs its trimmed edge-distance envelope",
                ));
            };
            if sign(&p.normal.cross(&self.axis).norm_squared(), dec)? != Sign::Zero
                || sign(&p.normal.dot(&self.axis), dec)? != Sign::Positive
            {
                return Err(bad(
                    "this conical opening needs its circular-edge distance envelope",
                ));
            }
        }
        let surfaces = self
            .faces
            .iter()
            .zip(open)
            .map(|(f, opened)| {
                if *opened {
                    Ok(Some((f.surface.clone(), f.origins.clone())))
                } else {
                    let origins = f
                        .origins
                        .iter()
                        .map(|of| Origin {
                            step: step.clone(),
                            role: Role::Inner {
                                of: Box::new(of.clone()),
                            },
                        })
                        .collect();
                    Ok(f.surface
                        .parallel_support(&(-thickness), dec)?
                        .map(|s| (s, origins)))
                }
            })
            .collect::<Result<Option<Vec<_>>, Refusal>>()?;
        let Some(surfaces) = surfaces else {
            return Ok(None);
        };
        match Self::new(&surfaces, dec) {
            Err(e) if e.code == "distance.empty-cone" => Ok(None),
            other => other,
        }
    }

    /// Convex defining supports plus an unregistered halfspace excluding the
    /// algebraic cone's other nappe. That halfspace is not a physical face.
    pub fn constraints(&self) -> Vec<Surface> {
        let mut surfaces: Vec<_> = self.faces.iter().map(|f| f.surface.clone()).collect();
        if let Some(ConeFace {
            surface: Surface::Cone(c),
            ..
        }) = self
            .faces
            .iter()
            .find(|f| matches!(f.surface, Surface::Cone(_)))
        {
            surfaces.push(Surface::Plane(Plane {
                origin: c.apex.clone(),
                normal: self.axis.scaled(&Real::from_i64(-1)),
                ref_dir: c.ref_dir.clone(),
            }));
        }
        surfaces
    }

    /// Exact coordinate extrema of the finite sections and retained apex.
    /// This is a query, not work performed while constructing a part.
    pub fn bounds(&self, dec: &mut Decider) -> Result<[[Real; 2]; 3], Refusal> {
        let mut result = [
            [Real::zero(), Real::zero()],
            [Real::zero(), Real::zero()],
            [Real::zero(), Real::zero()],
        ];
        for (i, bounds) in result.iter_mut().enumerate() {
            let mut n = [Real::zero(), Real::zero(), Real::zero()];
            n[i] = Real::one();
            let n = Vector3::new(n[0].clone(), n[1].clone(), n[2].clone());
            let mut initialized = false;
            for b in &self.boundaries {
                let range = extrema(&b.curve, &n, &Point3::origin(), dec)?;
                if !initialized {
                    *bounds = range;
                    initialized = true;
                } else {
                    if sign(&(&range[0] - &bounds[0]), dec)? == Sign::Negative {
                        bounds[0] = range[0].clone();
                    }
                    if sign(&(&range[1] - &bounds[1]), dec)? == Sign::Positive {
                        bounds[1] = range[1].clone();
                    }
                }
            }
            if let Some((p, _)) = &self.apex {
                let value = [&p.x, &p.y, &p.z][i];
                if sign(&(value - &bounds[0]), dec)? == Sign::Negative {
                    bounds[0] = value.clone();
                }
                if sign(&(value - &bounds[1]), dec)? == Sign::Positive {
                    bounds[1] = value.clone();
                }
            }
        }
        Ok(result)
    }
}

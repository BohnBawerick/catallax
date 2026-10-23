//! Exact meridians and named tools for recipe holes.

use catallax_numbers::{
    Decider, Outcome, Rational, Real, Refusal, Sign, StrictSign, sin_cos_degrees,
};
use catallax_recipe::{
    HoleDepth, HoleHead, HolePart, HolePlacement, HoleTip, Origin, Role, StepId,
};
use catallax_shapes::{Circle, Cone, Curve, Cylinder, Plane, Point3, PointDef, Surface, Vector3};
use catallax_topology::{Body, BodyId, BodyKind, Edge, Lump, Model, Sense, Shell, Vertex};

/// One positive-height cylindrical or conical band of a drill tool.
#[derive(Clone, Debug)]
pub struct HoleBand {
    /// Axial depths, increasing in the drilling direction.
    pub depth: [Real; 2],
    /// Radii at the two depths.
    pub radius: [Real; 2],
    /// The named lateral surface.
    pub part: HolePart,
    /// Real planar ledges at the upper end, not construction caps.
    pub end: Vec<HolePart>,
}

/// Exact finite tool whose construction caps are strictly outside its target.
/// The bound is a certified sphere about world zero supplied by the evaluator.
#[derive(Clone, Debug)]
pub struct HoleTool {
    /// The authored placement point.
    pub point: Point3,
    /// Unit drilling direction, including algebraic normalisation.
    pub axis: Vector3,
    /// Deterministic perpendicular seam direction.
    pub reference: Vector3,
    /// Bands in increasing depth order, before regularised subtraction.
    pub bands: Vec<HoleBand>,
}

fn sign(value: &Real, dec: &mut Decider) -> Result<Sign, Refusal> {
    if let Some(q) = value.rational_ref() {
        return Ok(q.sign());
    }
    match dec.sign(value) {
        Outcome::Proven(StrictSign::Positive) => Ok(Sign::Positive),
        Outcome::Proven(StrictSign::Negative) => Ok(Sign::Negative),
        Outcome::ProvenCoincident(_) => Ok(Sign::Zero),
        Outcome::NotExact(n) => Err(Refusal::new(n.code, n.detail)),
        Outcome::Refused(r) => Err(r),
    }
}
fn positive(q: &Rational, code: &'static str) -> Result<(), Refusal> {
    if q.sign() == Sign::Positive {
        Ok(())
    } else {
        Err(Refusal::new(code, "hole dimension must be positive"))
    }
}
fn angle_range(angle: &Rational) -> Result<(), Refusal> {
    if *angle <= Rational::zero() || *angle >= Rational::from_i64(180) {
        Err(Refusal::new(
            "hole.invalid-angle",
            "full hole angle must lie strictly between 0 and 180 degrees",
        ))
    } else {
        Ok(())
    }
}

/// Checks dimensions even when a hole has no placements or its target is empty.
pub fn validate(
    diameter: &Rational,
    depth: &HoleDepth,
    tip: &HoleTip,
    head: Option<&HoleHead>,
) -> Result<(), Refusal> {
    positive(diameter, "hole.non-positive-diameter")?;
    if let HoleDepth::Blind(h) = depth {
        positive(h.value(), "hole.non-positive-depth")?;
    }
    if let HoleTip::Cone(a) = tip {
        angle_range(a.value())?;
    }
    if let Some(head) = head {
        let d = match head {
            HoleHead::Counterbore { diameter, depth } => {
                positive(depth.value(), "hole.non-positive-head-depth")?;
                diameter
            }
            HoleHead::Countersink { diameter, angle } => {
                angle_range(angle.value())?;
                diameter
            }
        };
        if d.value() <= diameter {
            return Err(Refusal::new(
                "hole.head-too-small",
                "head diameter must exceed the main diameter",
            ));
        }
    }
    Ok(())
}

fn slope(angle: &Rational, dec: &mut Decider) -> Result<Real, Refusal> {
    angle_range(angle)?;
    let [sine, cosine] = sin_cos_degrees(angle, dec)?;
    Ok(sine / (Real::one() + cosine))
}

impl HoleTool {
    /// Compiles the recipe definition without constructing topology.
    pub fn new(
        at: &HolePlacement,
        diameter: &Rational,
        depth: &HoleDepth,
        tip: &HoleTip,
        head: Option<&HoleHead>,
        bound: &Rational,
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        validate(diameter, depth, tip, head)?;
        let tip_slope = match (depth, tip) {
            (HoleDepth::Blind(_), HoleTip::Cone(a)) => Some(slope(a.value(), dec)?),
            _ => None,
        };
        let half = Rational::from_ratio(1, 2).expect("fixed denominator");
        let r = Real::from(diameter * &half);
        let (entry, head_depth, head_slope) = match head {
            None => (r.clone(), Real::zero(), None),
            Some(HoleHead::Counterbore {
                diameter: d,
                depth: h,
            }) => {
                if d.value() <= diameter {
                    return Err(Refusal::new(
                        "hole.head-too-small",
                        "head diameter must exceed the main diameter",
                    ));
                }
                positive(h.value(), "hole.non-positive-head-depth")?;
                (
                    Real::from(d.value() * &half),
                    Real::from(h.value().clone()),
                    None,
                )
            }
            Some(HoleHead::Countersink { diameter: d, angle }) => {
                if d.value() <= diameter {
                    return Err(Refusal::new(
                        "hole.head-too-small",
                        "head diameter must exceed the main diameter",
                    ));
                }
                let k = slope(angle.value(), dec)?;
                let entry = Real::from(d.value() * &half);
                let h = &(&entry - &r) / &k;
                (entry, h, Some(k))
            }
        };
        let [x, y, z] = at
            .direction
            .0
            .each_ref()
            .map(|n| Real::from(n.value().clone()));
        let direction = Vector3::new(x, y, z);
        if direction
            .norm_squared()
            .rational_ref()
            .is_some_and(Rational::is_zero)
        {
            return Err(Refusal::new(
                "hole.zero-direction",
                "drilling direction must not be zero",
            ));
        }
        let length = direction.norm_squared().sqrt();
        // Preserve rational Pythagorean frames as rational expressions.
        let length = direction
            .norm_squared()
            .rational_ref()
            .and_then(crate::extrude::sqrt)
            .map(Real::from)
            .unwrap_or(length);
        let axis = direction.scaled(&(&Real::one() / &length));
        let basis = if at.direction.0[0].value().is_zero() && at.direction.0[1].value().is_zero() {
            Vector3::new(Real::one(), Real::zero(), Real::zero())
        } else {
            Vector3::new(Real::zero(), Real::zero(), Real::one())
        };
        let reference = axis.cross(&basis);
        let [x, y, z] = at.point.0.each_ref().map(|n| Real::from(n.value().clone()));
        let point = Point3::new(x, y, z);
        let extent = bound.clone()
            + at.point
                .0
                .iter()
                .fold(Rational::one(), |s, n| s + n.value().abs());
        let h = match depth {
            HoleDepth::Through => Real::from(extent.clone()),
            HoleDepth::Blind(h) => Real::from(h.value().clone()),
        };
        let mut bands = vec![HoleBand {
            depth: [Real::from(-extent), Real::zero()],
            radius: [entry.clone(), entry.clone()],
            part: HolePart::Entry,
            end: vec![],
        }];
        let mut cursor = Real::zero();
        if head.is_some() {
            let complete = sign(&(&head_depth - &h), dec)? != Sign::Positive;
            let end = if complete {
                head_depth.clone()
            } else {
                h.clone()
            };
            let end_radius = if let Some(k) = &head_slope {
                &entry - &(&end * k)
            } else {
                entry.clone()
            };
            bands.push(HoleBand {
                depth: [Real::zero(), end.clone()],
                radius: [entry, end_radius],
                part: if head_slope.is_some() {
                    HolePart::Countersink
                } else {
                    HolePart::CounterboreWall
                },
                end: if complete && head_slope.is_none() {
                    vec![HolePart::CounterboreFloor]
                } else {
                    vec![]
                },
            });
            cursor = end;
        }
        if sign(&(&h - &cursor), dec)? == Sign::Positive {
            bands.push(HoleBand {
                depth: [cursor, h.clone()],
                radius: [r.clone(), r.clone()],
                part: HolePart::Wall,
                end: vec![],
            });
        }
        if let HoleDepth::Blind(_) = depth {
            let last = bands.last_mut().expect("positive depth");
            if let Some(k) = tip_slope {
                if sign(&(&last.radius[1] - &r), dec)? == Sign::Positive {
                    last.end.push(HolePart::Floor);
                }
                bands.push(HoleBand {
                    depth: [h.clone(), &h + &(&r / &k)],
                    radius: [r, Real::zero()],
                    part: HolePart::Tip,
                    end: vec![],
                });
            } else {
                last.end.push(HolePart::Floor);
            }
        }
        Ok(Self {
            point,
            axis,
            reference,
            bands,
        })
    }

    /// Proves disjoint tools from their enclosing parallel cylinders.
    pub fn disjoint(&self, other: &Self, dec: &mut Decider) -> Result<bool, Refusal> {
        if sign(&self.axis.cross(&other.axis).norm_squared(), dec)? != Sign::Zero {
            return Ok(false);
        }
        let distance = self
            .point
            .minus(&other.point)
            .cross(&self.axis)
            .norm_squared();
        let radii = &self.bands[0].radius[0] + &other.bands[0].radius[0];
        Ok(sign(&(&distance - &radii.square()), dec)? == Sign::Positive)
    }

    /// A model-space axis point at an exact axial depth.
    pub fn at_depth(&self, depth: &Real) -> Point3 {
        self.point.plus(&self.axis.scaled(depth))
    }

    /// The lateral support of a band, with the cone's retained nappe.
    pub fn surface(&self, band: &HoleBand) -> Surface {
        let difference = &band.radius[0] - &band.radius[1];
        if difference.rational_ref().is_some_and(Rational::is_zero)
            || band.radius[0].same_definition(&band.radius[1])
        {
            Surface::Cylinder(Cylinder {
                origin: self.at_depth(&band.depth[0]),
                axis: self.axis.clone(),
                radius: band.radius[0].clone(),
                ref_dir: self.reference.clone(),
            })
        } else {
            let slope = &difference / &(&band.depth[1] - &band.depth[0]);
            Surface::Cone(Cone {
                apex: self.at_depth(&(&band.depth[0] + &(&band.radius[0] / &slope))),
                axis: self.axis.scaled(&Real::from_i64(-1)),
                slope,
                ref_dir: self.reference.clone(),
            })
        }
    }

    /// A cap support, outward from the tool.
    pub fn cap(&self, depth: &Real, upper: bool) -> Surface {
        Surface::Plane(Plane {
            origin: self.at_depth(depth),
            normal: self
                .axis
                .scaled(&Real::from_i64(if upper { 1 } else { -1 })),
            ref_dir: self.reference.clone(),
        })
    }

    /// Builds the connected meridian directly, without intermediate tool unions.
    pub fn build(
        &self,
        model: &mut Model,
        step: &StepId,
        at: &catallax_recipe::ElementId,
    ) -> BodyId {
        let origins = |parts: &[HolePart]| {
            parts
                .iter()
                .map(|part| Origin {
                    step: step.clone(),
                    role: Role::Hole {
                        at: at.clone(),
                        part: *part,
                    },
                })
                .collect()
        };
        let ring = |model: &mut Model, depth: &Real, radius: &Real| {
            let curve = model.geom_mut().add_curve(Curve::Circle(Circle {
                center: self.at_depth(depth),
                normal: self.axis.clone(),
                radius_sq: radius.square(),
                ref_dir: self.reference.clone(),
            }));
            model.add_edge(Edge {
                curve,
                sense: Sense::Same,
                start: None,
                end: None,
            })
        };
        let first = &self.bands[0];
        let mut lower = ring(model, &first.depth[0], &first.radius[0]);
        let mut faces = vec![crate::extrude::face(
            model,
            self.cap(&first.depth[0], false),
            Sense::Same,
            vec![vec![(lower, Sense::Reversed)]],
            None,
            vec![],
        )];
        for (i, band) in self.bands.iter().enumerate() {
            let apex = band.radius[1].rational_ref().is_some_and(Rational::is_zero);
            let upper = if apex {
                None
            } else {
                Some(ring(model, &band.depth[1], &band.radius[1]))
            };
            let mut paths = vec![vec![(lower, Sense::Same)]];
            let vertex = if let Some(edge) = upper {
                paths.push(vec![(edge, Sense::Reversed)]);
                None
            } else {
                let point = model
                    .geom_mut()
                    .add_point(PointDef::Explicit(self.at_depth(&band.depth[1])));
                Some(model.add_vertex(Vertex { point }))
            };
            faces.push(crate::extrude::face(
                model,
                self.surface(band),
                Sense::Same,
                paths,
                vertex,
                origins(&[band.part]),
            ));
            if let Some(upper) = upper {
                if let Some(next) = self.bands.get(i + 1) {
                    if next.radius[0].same_definition(&band.radius[1])
                        || (&next.radius[0] - &band.radius[1])
                            .rational_ref()
                            .is_some_and(Rational::is_zero)
                    {
                        lower = upper;
                    } else {
                        lower = ring(model, &next.depth[0], &next.radius[0]);
                        faces.push(crate::extrude::face(
                            model,
                            self.cap(&band.depth[1], true),
                            Sense::Same,
                            vec![vec![(upper, Sense::Same)], vec![(lower, Sense::Reversed)]],
                            None,
                            origins(&band.end),
                        ));
                    }
                } else {
                    faces.push(crate::extrude::face(
                        model,
                        self.cap(&band.depth[1], true),
                        Sense::Same,
                        vec![vec![(upper, Sense::Same)]],
                        None,
                        origins(&band.end),
                    ));
                }
            }
        }
        let shell = model.add_shell(Shell { faces });
        let lump = model.add_lump(Lump {
            shells: vec![shell],
        });
        model.add_body(Body {
            kind: BodyKind::Solid,
            lumps: vec![lump],
            marks: Default::default(),
        })
    }
}

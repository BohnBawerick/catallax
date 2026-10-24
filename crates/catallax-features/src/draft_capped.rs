//! Convex cylindrical material between two affine-height planes.
//! The carrier, cap definitions and source origins are retained exactly.

use super::{reduced, sign, unit};
use catallax_numbers::{Decider, Rational, Real, Refusal, Sign};
use catallax_recipe::Origin;
use catallax_shapes::{Cylinder, Isometry, Plane, Point3, Surface, Vector3};

#[derive(Clone, Debug)]
struct Height {
    offset: Real,
    linear: [Real; 2],
}

fn dot(a: &[Real; 2], b: &[Real; 2]) -> Real {
    &a[0] * &b[0] + &a[1] * &b[1]
}

/// Circle-domain topology after intersecting the two cap height laws.
#[derive(Clone, Debug)]
pub enum FootprintBoundary {
    /// The cap ellipses are disjoint closed curves.
    Closed,
    /// A positive-length chord joins these counter-clockwise arc endpoints.
    Chord([[Real; 2]; 2]),
    /// Both closed cap curves meet at this single circle point.
    Tangent([Real; 2]),
}

/// A full cylindrical side clipped by two opposing, possibly oblique planes.
#[derive(Clone, Debug)]
pub struct CappedCylinder {
    /// Unchanged cylindrical carrier.
    pub cylinder: Cylinder,
    /// Physical side origins, not construction-tool names.
    pub origins: Vec<Origin>,
    /// Lower and upper outward cap definitions and their origins.
    pub caps: [(Plane, Vec<Origin>); 2],
    /// Unit source axis and radial chart, independent of the draft pull.
    pub axis: Vector3,
    /// First unit radial chart direction.
    pub u: Vector3,
    /// Second unit radial chart direction.
    pub v: Vector3,
    /// Exact source radius.
    pub radius: Real,
    heights: [Height; 2],
    thickness: Height,
    extent: [Real; 2],
}

impl CappedCylinder {
    /// Constructs the bounded full-dimensional intersection. No fitted cap or
    /// restriction to the source's former height interval is introduced.
    pub fn new(
        cylinder: Cylinder,
        origins: Vec<Origin>,
        caps: [(Plane, Vec<Origin>); 2],
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        let axis = unit(&cylinder.axis, dec)?;
        let u = unit(&cylinder.ref_dir, dec)?;
        let v = axis.cross(&u);
        let radius = reduced(cylinder.radius_sq.sqrt(), dec)?;
        let mut ordered = [None, None];
        let mut heights = [None, None];
        for (plane, names) in caps {
            let rate = plane.normal.dot(&axis);
            let index = match sign(&rate, dec)? {
                Sign::Negative => 0,
                Sign::Positive => 1,
                Sign::Zero => {
                    return Err(Refusal::new(
                        "draft.unsupported-round-cap",
                        "vertical cylindrical cuts need a different bounded cell decomposition",
                    ));
                }
            };
            if ordered[index].is_some() {
                return Err(Refusal::new(
                    "draft.nonclosing-boundary",
                    "cylindrical material needs opposing axial caps",
                ));
            }
            heights[index] = Some(Height {
                offset: &plane.origin.minus(&cylinder.origin).dot(&plane.normal) / &rate,
                linear: [
                    &(-plane.normal.dot(&u)) / &rate,
                    &(-plane.normal.dot(&v)) / &rate,
                ],
            });
            ordered[index] = Some((plane, names));
        }
        let [Some(lower), Some(upper)] = heights else {
            return Err(Refusal::new(
                "draft.nonclosing-boundary",
                "cylindrical material needs opposing axial caps",
            ));
        };
        let thickness = Height {
            offset: &upper.offset - &lower.offset,
            linear: std::array::from_fn(|i| &upper.linear[i] - &lower.linear[i]),
        };
        let thickness_norm = reduced(dot(&thickness.linear, &thickness.linear).sqrt(), dec)?;
        if sign(&(&thickness.offset + &(&radius * &thickness_norm)), dec)? != Sign::Positive {
            return Err(Refusal::new(
                "draft.reversed-solid",
                "replacement caps leave no positive cylindrical interior",
            ));
        }
        let lower_norm = reduced(dot(&lower.linear, &lower.linear).sqrt(), dec)?;
        let upper_norm = reduced(dot(&upper.linear, &upper.linear).sqrt(), dec)?;
        let extent = [
            &lower.offset - &(&radius * &lower_norm),
            &upper.offset + &(&radius * &upper_norm),
        ];
        Ok(Self {
            cylinder,
            origins,
            caps: ordered.map(|p| p.expect("opposing caps checked")),
            axis,
            u,
            v,
            radius,
            heights: [lower, upper],
            thickness,
            extent,
        })
    }

    /// Conservative exact axial bounds, including all extended neighbours.
    pub fn extent(&self) -> &[Real; 2] {
        &self.extent
    }

    /// Maps the native radial coordinates and axial height into world space.
    pub fn point(&self, q: &[Real; 2], h: &Real) -> Point3 {
        self.cylinder
            .origin
            .plus(&self.u.scaled(&q[0]))
            .plus(&self.v.scaled(&q[1]))
            .plus(&self.axis.scaled(h))
    }

    /// Exact heights above a radial point, before rejecting an empty fibre.
    pub fn heights_at(&self, q: &[Real; 2]) -> [Real; 2] {
        self.heights
            .each_ref()
            .map(|h| &h.offset + &dot(&h.linear, q))
    }

    /// Exact cap conic, positively oriented in the local radial chart.
    pub fn cap_curve(
        &self,
        index: usize,
        dec: &mut Decider,
    ) -> Result<catallax_shapes::Curve, Refusal> {
        let height = &self.heights[index];
        let beta2 = dot(&height.linear, &height.linear);
        let center = self.cylinder.origin.plus(&self.axis.scaled(&height.offset));
        let first = self.u.plus(&self.axis.scaled(&height.linear[0]));
        let second = self.v.plus(&self.axis.scaled(&height.linear[1]));
        let normal = first.cross(&second);
        if sign(&beta2, dec)? == Sign::Zero {
            return Ok(catallax_shapes::Curve::Circle(catallax_shapes::Circle {
                center,
                normal,
                radius_sq: self.cylinder.radius_sq.clone(),
                ref_dir: self.u.clone(),
            }));
        }
        let tilt = self
            .u
            .scaled(&height.linear[0])
            .plus(&self.v.scaled(&height.linear[1]));
        Ok(catallax_shapes::Curve::Ellipse(catallax_shapes::Ellipse {
            center,
            normal,
            major_dir: tilt.plus(&self.axis.scaled(&beta2)),
            major_sq: &self.cylinder.radius_sq * &(&Real::one() + &beta2),
            minor_sq: self.cylinder.radius_sq.clone(),
        }))
    }

    /// Intersections of the cap equality chord with the exact source circle.
    pub fn footprint_boundary(&self, dec: &mut Decider) -> Result<FootprintBoundary, Refusal> {
        let t = &self.thickness.linear;
        let norm2 = dot(t, t);
        if sign(&norm2, dec)? == Sign::Zero {
            return Ok(FootprintBoundary::Closed);
        }
        let norm = reduced(norm2.sqrt(), dec)?;
        let minimum = &self.thickness.offset - &(&self.radius * &norm);
        let side = sign(&minimum, dec)?;
        if side == Sign::Positive {
            return Ok(FootprintBoundary::Closed);
        }
        let centre = t
            .each_ref()
            .map(|x| x * &(&(-&self.thickness.offset) / &norm2));
        if side == Sign::Zero {
            return Ok(FootprintBoundary::Tangent(centre));
        }
        let chord2 = &self.cylinder.radius_sq - &dot(&centre, &centre);
        let scale = reduced((&chord2 / &norm2).sqrt(), dec)?;
        let tangent = [-&t[1], t[0].clone()];
        Ok(FootprintBoundary::Chord([
            std::array::from_fn(|i| &centre[i] - &(&tangent[i] * &scale)),
            std::array::from_fn(|i| &centre[i] + &(&tangent[i] * &scale)),
        ]))
    }

    /// Native defining supports, in side/lower/upper order.
    pub fn surfaces(&self) -> [Surface; 3] {
        [
            Surface::Cylinder(self.cylinder.clone()),
            Surface::Plane(self.caps[0].0.clone()),
            Surface::Plane(self.caps[1].0.clone()),
        ]
    }

    /// Strict retained-region membership; boundaries cannot choose a component.
    pub fn contains_interior(&self, point: &Point3, dec: &mut Decider) -> Result<bool, Refusal> {
        if sign(
            &(&point
                .minus(&self.cylinder.origin)
                .cross(&self.axis)
                .norm_squared()
                - &self.cylinder.radius_sq),
            dec,
        )? != Sign::Negative
        {
            return Ok(false);
        }
        for (plane, _) in &self.caps {
            if sign(&point.minus(&plane.origin).dot(&plane.normal), dec)? != Sign::Negative {
                return Ok(false);
            }
        }
        Ok(true)
    }

    /// A strict algebraic witness, including a footprint not containing the axis.
    pub fn interior_witness(&self, dec: &mut Decider) -> Result<Point3, Refusal> {
        let q = if sign(&self.thickness.offset, dec)? == Sign::Positive {
            [Real::zero(), Real::zero()]
        } else {
            let norm = dot(&self.thickness.linear, &self.thickness.linear).sqrt();
            let scale = &(&Real::one() - &(&self.thickness.offset / &(&self.radius * &norm)))
                / &Real::from_i64(2);
            self.thickness
                .linear
                .each_ref()
                .map(|x| x * &(&(&scale * &self.radius) / &norm))
        };
        let h = self.heights_at(&q);
        Ok(self.point(&q, &((&h[0] + &h[1]) / Real::from_i64(2))))
    }

    /// Rigid copies preserve the local affine-height law, including reflections.
    pub fn transform(&mut self, pose: &Isometry) {
        let Surface::Cylinder(cylinder) = pose.surface(&Surface::Cylinder(self.cylinder.clone()))
        else {
            unreachable!("an isometry preserves a cylinder")
        };
        self.cylinder = cylinder;
        for (plane, _) in &mut self.caps {
            let Surface::Plane(transformed) = pose.surface(&Surface::Plane(plane.clone())) else {
                unreachable!("an isometry preserves a plane")
            };
            *plane = transformed;
        }
        self.axis = pose.apply_vector(&self.axis);
        self.u = pose.apply_vector(&self.u);
        self.v = pose.apply_vector(&self.v);
    }

    // Maximises a linear functional over the closed footprint: a disk clipped
    // by the affine nonnegative-thickness halfplane. If the disk maximiser is
    // excluded, the maximum is at an endpoint of the clipping chord.
    fn maximum(&self, linear: &[Real; 2], dec: &mut Decider) -> Result<Real, Refusal> {
        let norm2 = dot(linear, linear);
        if sign(&norm2, dec)? == Sign::Zero {
            return Ok(Real::zero());
        }
        let norm = norm2.sqrt();
        let q = linear.each_ref().map(|x| x * &(&self.radius / &norm));
        if sign(
            &(&self.thickness.offset + &dot(&self.thickness.linear, &q)),
            dec,
        )? != Sign::Negative
        {
            return Ok(&self.radius * &norm);
        }
        let t = &self.thickness.linear;
        let t2 = dot(t, t);
        // A valid interior and an excluded disk maximum imply a proper chord.
        let centre = t
            .each_ref()
            .map(|x| x * &(&(-&self.thickness.offset) / &t2));
        let chord2 = &self.cylinder.radius_sq - &dot(&centre, &centre);
        let transverse = &linear[0] * &t[1] - &linear[1] * &t[0];
        Ok(dot(linear, &centre) + (&(transverse.square() * chord2) / &t2).sqrt())
    }

    /// Tests two-dimensional neutral contact, including a coincident cap but
    /// excluding a tangent generator, rim or vertex.
    pub fn anchored(&self, neutral: &Plane, dec: &mut Decider) -> Result<bool, Refusal> {
        let rate = neutral.normal.dot(&self.axis);
        let radial = [neutral.normal.dot(&self.u), neutral.normal.dot(&self.v)];
        let offset = self
            .cylinder
            .origin
            .minus(&neutral.origin)
            .dot(&neutral.normal);
        let lower = usize::from(sign(&rate, dec)? == Sign::Negative);
        let mut extrema = [Real::zero(), Real::zero()];
        for (i, height) in [&self.heights[lower], &self.heights[1 - lower]]
            .into_iter()
            .enumerate()
        {
            let linear = std::array::from_fn(|j| &radial[j] + &(&rate * &height.linear[j]));
            let centre = &offset + &(&rate * &height.offset);
            extrema[i] = if i == 0 {
                let negative = linear.each_ref().map(|x| -x);
                centre - self.maximum(&negative, dec)?
            } else {
                centre + self.maximum(&linear, dec)?
            };
        }
        if sign(&extrema[0], dec)? == Sign::Negative && sign(&extrema[1], dec)? == Sign::Positive {
            return Ok(true);
        }
        for (plane, _) in &self.caps {
            if sign(&plane.normal.cross(&neutral.normal).norm_squared(), dec)? == Sign::Zero
                && sign(
                    &plane.origin.minus(&neutral.origin).dot(&neutral.normal),
                    dec,
                )? == Sign::Zero
            {
                return Ok(true);
            }
        }
        Ok(false)
    }

    /// Defining world-space enclosure for field queries and construction tools.
    pub fn bound_radius(&self, dec: &mut Decider) -> Result<Rational, Refusal> {
        let mut bound = Rational::zero();
        for height in &self.extent {
            let centre = self.cylinder.origin.plus(&self.axis.scaled(height));
            let mut axial_bound = Rational::zero();
            for x in [&centre.x, &centre.y, &centre.z] {
                let interval = dec.enclose_at(x, 128)?;
                axial_bound = axial_bound + interval.lo().abs().max(interval.hi().abs());
            }
            bound = bound.max(axial_bound);
        }
        Ok(bound + Rational::from_i64(3) * dec.enclose_at(&self.radius, 128)?.hi())
    }
}

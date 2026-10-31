//! Early preview refusal for a capsule carrier buried in empty space.
//!
//! Trim constraints are complete supports, so a capsule's unused half circle
//! still cuts every planar cell it crosses. Its sign changes in each such cell,
//! which keeps the cell unstable even where the capsule's membership is constant.
//! Where the result has no boundary either, nothing emits or covers that cell and
//! refinement only spends the allowance. This module proves that case in exact
//! arithmetic at the antipode of each capsule arc, over an open cube around it:
//!
//! - the carrier has a regular zero at the antipode, on a cap chart's domain;
//! - the convex capsule contains the cube, so the arc itself is inactive there;
//! - every operand is constant on both open sides of the cap, and the result
//!   agrees across them, so the cube holds no boundary;
//! - the cap's symbolic side states still depend on the capsule operand;
//! - the cube's half width exceeds twice the tolerance. A vertex lies within one
//!   tolerance of the boundary and covers cells within one tolerance of itself.
//!
//! A failed step proves nothing and leaves the traversal untouched.

use super::{Builder, Patch, failure, membership, preview, q};
use crate::{Instruction, Location, Plan, geometry::Leaf, profile};
use catallax_numbers::{Decider, Interval, Rational, Refusal};
use catallax_recipe::Turn;
use catallax_shapes::{Isometry, RPoint3, RVector3};

type P = [Rational; 2];
type Local = ([[Rational; 3]; 3], [Rational; 3]);

/// A swept stadium: two equal outward semicircles joined by tangent segments.
pub(crate) struct Capsule<'a> {
    frame: Local,
    placement: Isometry,
    vector: &'a RVector3,
    through: bool,
    centres: [P; 2],
    squared: Rational,
    // The point of each arc's circle diametrically opposite its active middle.
    antipodes: [P; 2],
}

// Middle of a semicircle, a quarter turn from its start about the centre.
fn middle(arc: &profile::Edge) -> Option<P> {
    let (centre, _, turn) = arc.circle.as_ref()?;
    let [x, y] = profile::sub(&arc.start, centre);
    let (x, y) = if *turn == Turn::Ccw { (-y, x) } else { (y, -x) };
    Some([&centre[0] + &x, &centre[1] + &y])
}

impl<'a> Capsule<'a> {
    pub(crate) fn new(
        inverse: &Isometry,
        vector: &'a RVector3,
        through: bool,
        loops: &'a [Vec<profile::Edge>],
    ) -> Option<Self> {
        let [edges] = loops else {
            return None;
        };
        let [a, b, c, d] = edges.as_slice() else {
            return None;
        };
        let (arcs, lines) = if a.circle.is_some() {
            ([a, c], [b, d])
        } else {
            ([b, d], [c, a])
        };
        let (first, squared, turn) = arcs[0].circle.as_ref()?;
        let (second, other, other_turn) = arcs[1].circle.as_ref()?;
        let axis = profile::sub(second, first);
        let semicircle = |arc: &profile::Edge, centre: &P| {
            !arc.full && (0..2).all(|i| &arc.start[i] + &arc.end[i] == centre[i].mul_pow2(1))
        };
        if squared <= &q(0)
            || squared != other
            || turn != other_turn
            || lines.iter().any(|line| line.circle.is_some())
            || !semicircle(arcs[0], first)
            || !semicircle(arcs[1], second)
            || profile::sub(&lines[0].end, &lines[0].start) != axis
            || profile::sub(&lines[1].start, &lines[1].end) != axis
            || !profile::dot(&profile::sub(&arcs[0].end, first), &axis).is_zero()
        {
            return None;
        }
        let middles = [middle(arcs[0])?, middle(arcs[1])?];
        let outward = [
            profile::dot(&profile::sub(&middles[0], first), &axis),
            profile::dot(&profile::sub(&middles[1], second), &axis),
        ];
        if outward[0] >= q(0) || outward[1] <= q(0) {
            return None;
        }
        let [left, right] = middles;
        let opposite =
            |centre: &P, middle: &P| std::array::from_fn(|i| &centre[i].mul_pow2(1) - &middle[i]);
        Some(Self {
            frame: membership::rational_frame(inverse)?,
            placement: inverse.inverse(),
            vector,
            through,
            antipodes: [opposite(first, &left), opposite(second, &right)],
            centres: [first.clone(), second.clone()],
            squared: squared.clone(),
        })
    }

    // Sweep parameter and profile point, as in the feature definition.
    fn project(&self, p: &[Rational; 3]) -> (Rational, P) {
        let (matrix, translation) = &self.frame;
        let local: [Rational; 3] = std::array::from_fn(|i| {
            (0..3).fold(translation[i].clone(), |s, j| s + &matrix[i][j] * &p[j])
        });
        let t = &local[2] / &self.vector.z;
        (
            t.clone(),
            [
                &local[0] - &(&t * &self.vector.x),
                &local[1] - &(&t * &self.vector.y),
            ],
        )
    }

    fn world(&self, p: &P, t: &Rational) -> Option<[Rational; 3]> {
        let local = RPoint3::new(
            &p[0] + &(t * &self.vector.x),
            &p[1] + &(t * &self.vector.y),
            t * &self.vector.z,
        );
        let world = self.placement.apply_rational_point(&local)?;
        Some([world.x, world.y, world.z])
    }

    // Strictly inside the stadium: nearer its axis segment than the radius.
    fn contains(&self, p: &P) -> bool {
        let [a, b] = &self.centres;
        let axis = profile::sub(b, a);
        let offset = profile::sub(p, a);
        let along = profile::dot(&offset, &axis) / profile::dot(&axis, &axis);
        let along = along.max(q(0)).min(q(1));
        let nearest: P = std::array::from_fn(|i| &a[i] + &(&along * &axis[i]));
        let d = profile::sub(p, &nearest);
        profile::dot(&d, &d) < self.squared
    }

    // Whether membership is constant on each open side of the cap within the
    // cube: the convex stadium holds the cube's projection, and the only end
    // plane crossing the cube is a coincident cap through its centre.
    fn constant_sides(
        &self,
        point: &[Rational; 3],
        corners: &[[Rational; 3]],
        capped: bool,
        dec: &mut Decider,
    ) -> Result<bool, Refusal> {
        let (centre, _) = self.project(point);
        let mut range = [centre.clone(), centre.clone()];
        for corner in corners {
            dec.check_deadline()?;
            let (t, p) = self.project(corner);
            if !self.contains(&p) {
                return Ok(false);
            }
            range = [range[0].clone().min(t.clone()), range[1].clone().max(t)];
        }
        Ok(self.through
            || [q(0), q(1)]
                .iter()
                .all(|end| *end <= range[0] || *end >= range[1] || (capped && *end == centre)))
    }
}

// Evaluates the Boolean tape over one location per operand.
fn evaluate(plan: &Plan, operand: impl Fn(usize) -> Location) -> Location {
    let mut tape = Vec::with_capacity(plan.instructions.len());
    for instruction in &plan.instructions {
        let value = match *instruction {
            Instruction::Empty => Location::Outside,
            Instruction::Leaf(leaf) => operand(leaf),
            Instruction::Boolean { kind, a, b } => crate::combine(kind, tape[a], tape[b]),
        };
        tape.push(value);
    }
    tape[plan.root]
}

// An operand's location at the cap point and on the two open sides of the cap.
#[derive(Clone, Copy)]
struct Germs {
    cap: Location,
    sides: [Location; 2],
}
impl Germs {
    fn constant(location: Location) -> Self {
        Self {
            cap: location,
            sides: [location; 2],
        }
    }
}

struct Neighbourhood {
    point: [Rational; 3],
    cube: [Interval; 3],
    corners: Vec<[Rational; 3]>,
    sides: [[Rational; 3]; 2],
    // Sweep parameter of the cap point along the capsule under test.
    parameter: Rational,
}

impl Builder<'_> {
    // An exact refusal leaves the point unclassified; only allowances surface.
    fn exact_location(
        &self,
        leaf: usize,
        point: &[Rational; 3],
        dec: &mut Decider,
    ) -> Result<Option<Location>, Refusal> {
        match self.leaf_location(leaf, point, dec) {
            Ok(Location::Indeterminate) => Ok(None),
            Ok(location) => Ok(Some(location)),
            Err(refusal) if preview::is_allowance(&refusal) => Err(refusal),
            Err(_) => Ok(None),
        }
    }

    fn germs(
        &self,
        patch: &Patch,
        leaf: usize,
        capsules: &[Option<Capsule<'_>>],
        place: &Neighbourhood,
        dec: &mut Decider,
    ) -> Result<Option<Germs>, Refusal> {
        let Neighbourhood {
            point,
            cube,
            corners,
            sides,
            ..
        } = place;
        if self.source_boxes[leaf].as_ref().is_some_and(|bounds| {
            (0..3).any(|i| cube[i].hi() < bounds[i].lo() || cube[i].lo() > bounds[i].hi())
        }) {
            return Ok(Some(Germs::constant(Location::Outside)));
        }
        let constant = match &capsules[leaf] {
            Some(capsule) => {
                let capped = patch.coincident[leaf].iter().any(|&same| same);
                capsule.constant_sides(point, corners, capped, dec)?
            }
            None => {
                let mut constant = true;
                for (j, constraint) in self.constraints[leaf].iter().enumerate() {
                    dec.check_deadline()?;
                    // Only an affine coincident support vanishes on the cap alone.
                    constant &= if patch.coincident[leaf][j] {
                        constraint.linear.is_some()
                    } else {
                        constraint.range(cube).sign().is_some()
                    };
                }
                constant
            }
        };
        if !constant {
            return Ok(None);
        }
        let Some(cap) = self.exact_location(leaf, point, dec)? else {
            return Ok(None);
        };
        let mut locations = [Location::Outside; 2];
        for (location, side) in locations.iter_mut().zip(sides) {
            *location = match self.exact_location(leaf, side, dec)? {
                Some(location @ (Location::Inside | Location::Outside)) => location,
                _ => return Ok(None),
            };
        }
        Ok(Some(Germs {
            cap,
            sides: locations,
        }))
    }

    fn all_germs(
        &self,
        patch: &Patch,
        capsules: &[Option<Capsule<'_>>],
        place: &Neighbourhood,
        dec: &mut Decider,
    ) -> Result<Option<Vec<Germs>>, Refusal> {
        let mut germs = Vec::with_capacity(self.plan.leaves.len());
        for leaf in 0..self.plan.leaves.len() {
            dec.check_deadline()?;
            germs.push(if self.plan.leaf_registers[leaf].is_none() {
                Germs::constant(Location::Indeterminate)
            } else if let Some(germ) = self.germs(patch, leaf, capsules, place, dec)? {
                germ
            } else {
                return Ok(None);
            });
        }
        Ok(Some(germs))
    }

    // Lifts the capsule point to this cap and builds the cube around it.
    fn neighbourhood(
        &self,
        patch: &Patch,
        capsule: &Capsule<'_>,
        antipode: &P,
    ) -> Option<Neighbourhood> {
        let face = &self.plan.faces[patch.faces[0]];
        let support = face.polynomial.rational()?;
        let start = capsule.world(antipode, &q(0))?;
        let end = capsule.world(antipode, &q(1))?;
        let height = support.eval(&start);
        let rise = support.eval(&end) - height.clone();
        if rise.is_zero() {
            return None;
        }
        let parameter = -&height / rise;
        if !capsule.through && (parameter < q(0) || parameter > q(1)) {
            return None;
        }
        let point = capsule.world(antipode, &parameter)?;
        let half = &self.options.tolerance * &q(3);
        let sweep: [Rational; 3] = std::array::from_fn(|i| &end[i] - &start[i]);
        // A move of half the cube along the sweep stays inside it on each side.
        let step = &half / &sweep.iter().fold(q(0), |s, x| s + x.abs()).mul_pow2(1);
        let sides =
            [step.clone(), -&step].map(|s| std::array::from_fn(|i| &point[i] + &(&s * &sweep[i])));
        let cube = std::array::from_fn(|i| {
            let lo = Interval::from_rational(&(&point[i] - &half)).lo();
            let hi = Interval::from_rational(&(&point[i] + &half)).hi();
            Interval::new(lo, hi).expect("ordered cube bounds")
        });
        let corners = (0..8)
            .map(|corner| {
                std::array::from_fn(|i| {
                    if corner & (1 << i) == 0 {
                        &point[i] - &half
                    } else {
                        &point[i] + &half
                    }
                })
            })
            .collect();
        Some(Neighbourhood {
            point,
            cube,
            corners,
            sides,
            parameter,
        })
    }

    // The traversal reaches the cap point when the chart's closed domain holds
    // it and some source face's box cannot skip the cells around it.
    fn reaches(&self, patch: &Patch, point: &[Rational; 3]) -> bool {
        patch.chart.plane_contains(point)
            && patch.faces.iter().any(|&face| {
                self.source_boxes[self.plan.faces[face].leaf]
                    .as_ref()
                    .is_none_or(|bounds| {
                        (0..3).all(|i| {
                            Rational::from_f64_exact(bounds[i].lo()).is_none_or(|lo| lo <= point[i])
                                && Rational::from_f64_exact(bounds[i].hi())
                                    .is_none_or(|hi| point[i] <= hi)
                        })
                    })
            })
    }

    // Index of a capsule carrier with a regular zero at the cap point. The carrier
    // is constant along the sweep, so its restriction to this transverse cap
    // matches the radial line. A quadratic vanishing at this point with opposite
    // signs either side has a nonzero linear term there.
    fn carrier(
        &self,
        patch: &Patch,
        leaf: usize,
        capsule: &Capsule<'_>,
        centre: &P,
        place: &Neighbourhood,
        dec: &mut Decider,
    ) -> Result<Option<usize>, Refusal> {
        let Some(centre) = capsule.world(centre, &place.parameter) else {
            return Ok(None);
        };
        let radial: [Rational; 3] =
            std::array::from_fn(|i| (&place.point[i] - &centre[i]).mul_pow2(-1));
        let [outer, inner] = [q(1), q(-1)]
            .map(|s| std::array::from_fn::<_, 3, _>(|i| &place.point[i] + &(&s * &radial[i])));
        for (j, constraint) in self.constraints[leaf].iter().enumerate() {
            dec.check_deadline()?;
            if !patch.coincident[leaf][j]
                && constraint.polynomial.as_ref().is_some_and(|p| {
                    p.total_degree() == Some(2)
                        && p.eval(&place.point).is_zero()
                        && p.eval(&outer) * p.eval(&inner) < q(0)
                })
            {
                return Ok(Some(j));
            }
        }
        Ok(None)
    }

    // Whether the cap's side states stay unresolved for a free capsule operand,
    // while the operand's actual germs put the whole cube on one side of the result.
    fn cap_undecided(&self, patch: &Patch, leaf: usize, germs: &[Germs]) -> bool {
        let plan = self.plan;
        let actual = [0, 1].map(|side| evaluate(plan, |i| germs[i].sides[side]));
        let free = [0, 1].map(|side| {
            [Location::Inside, Location::Outside].map(|choice| {
                evaluate(plan, |i| {
                    if i == leaf {
                        choice
                    } else {
                        germs[i].sides[side]
                    }
                })
            })
        });
        let varies = free.iter().any(|[a, b]| a != b);
        let cap = evaluate(plan, |i| {
            if i == leaf {
                Location::Indeterminate
            } else {
                germs[i].cap
            }
        });
        actual[0] == actual[1]
            && matches!(actual[0], Location::Inside | Location::Outside)
            && varies
            && !(patch.shared_sides[leaf] && free[0] == free[1])
            && !matches!(cap, Location::Inside | Location::Outside)
    }

    fn obstruction(
        &self,
        patch: &Patch,
        leaf: usize,
        capsules: &[Option<Capsule<'_>>],
        dec: &mut Decider,
    ) -> Result<Option<Refusal>, Refusal> {
        let plan = self.plan;
        let Some(capsule) = &capsules[leaf] else {
            return Ok(None);
        };
        for (antipode, centre) in capsule.antipodes.iter().zip(&capsule.centres) {
            dec.check_deadline()?;
            let Some(place) = self.neighbourhood(patch, capsule, antipode) else {
                continue;
            };
            if !self.reaches(patch, &place.point) {
                continue;
            }
            let Some(carrier) = self.carrier(patch, leaf, capsule, centre, &place, dec)? else {
                continue;
            };
            let Some(germs) = self.all_germs(patch, capsules, &place, dec)? else {
                continue;
            };
            if !self.cap_undecided(patch, leaf, &germs) {
                continue;
            }
            let polynomial = &plan.leaves[leaf].constraints()[carrier].polynomial;
            let cap = patch
                .faces
                .iter()
                .find(|&&face| plan.faces[face].leaf == leaf)
                .unwrap_or(&patch.faces[0]);
            let cap = &plan.faces[*cap].origin;
            let source = plan
                .faces
                .iter()
                .find(|face| face.leaf == leaf && face.polynomial.same_definition(polynomial))
                .map_or(cap, |face| &face.origin);
            return Ok(Some(failure(
                "field.display-trim-obstruction",
                format!(
                    "unused half of capsule carrier {:?} crosses cap support {:?} at ({}, {}, {}); \
                     the result has no boundary within {}, more than twice the tolerance {}, \
                     so no cover of that carrier can resolve",
                    source,
                    cap,
                    place.point[0],
                    place.point[1],
                    place.point[2],
                    &self.options.tolerance * &q(3),
                    self.options.tolerance,
                ),
            )));
        }
        Ok(None)
    }
}

/// Refuses when a capsule carrier provably cannot obtain the required cover on
/// a planar support. Costs are bounded by supports, operands and carriers, and
/// charge the active allowance. Plans without a capsule profile do no work here.
pub(super) fn check(
    builder: &Builder<'_>,
    patches: &[Patch],
    dec: &mut Decider,
) -> Result<(), Refusal> {
    let capsules: Vec<_> = builder
        .plan
        .leaves
        .iter()
        .map(|leaf| match leaf {
            Leaf::Feature(feature) => feature.display_capsule(),
            _ => None,
        })
        .collect();
    if capsules.iter().all(Option::is_none) {
        return Ok(());
    }
    for patch in patches
        .iter()
        .filter(|p| p.native.is_none() && p.chart.is_plane())
    {
        for leaf in 0..capsules.len() {
            if builder.plan.leaf_registers[leaf].is_some()
                && let Some(refusal) = builder.obstruction(patch, leaf, &capsules, dec)?
            {
                return Err(refusal);
            }
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use crate::{
        Plan,
        display::{PreviewOptions, preview_plan},
    };
    use catallax_numbers::Decider;
    use catallax_recipe::Recipe;
    use serde_json::{Value, json};

    fn extrude(id: &str, loops: Value) -> [Value; 2] {
        let profile = format!("{id}-profile");
        [
            json!({"id":profile, "op":"sketch", "loops":loops,
                "plane":{"rotation":["-1","-1","0","0"],"translation":["0","0","0"]}}),
            json!({"id":id, "op":"extrude", "profile":profile, "distance":"10",
                "direction":["0","1","0"]}),
        ]
    }
    fn plate() -> [Value; 2] {
        extrude(
            "plate",
            json!([[
                {"kind":"line","id":"a","to":["120","0"]},
                {"kind":"line","id":"b","to":["120","60"]},
                {"kind":"line","id":"c","to":["0","60"]},
                {"kind":"line","id":"d","to":["0","0"]}
            ]]),
        )
    }
    fn slot() -> [Value; 2] {
        extrude(
            "slot",
            json!([[
                {"kind":"arc","id":"left","centre":["46","35"],"to":["46","29"],"turn":"ccw"},
                {"kind":"line","id":"low","to":["74","29"]},
                {"kind":"arc","id":"right","centre":["74","35"],"to":["74","41"],"turn":"ccw"},
                {"kind":"line","id":"high","to":["46","41"]}
            ]]),
        )
    }
    fn boolean(id: &str, kind: &str, target: &str, tool: &str) -> Value {
        json!({"id":id, "op":"boolean", "kind":kind, "target":target, "tool":tool})
    }
    fn preview(mut steps: Vec<Value>, tolerance: &str) -> Result<usize, &'static str> {
        let result = steps.last().unwrap()["id"].clone();
        let source = json!({"format":"catallax-recipe/4", "name":"t", "steps":steps.split_off(0),
            "result":result});
        let plan = Plan::compile(&Recipe::parse(&source.to_string()).unwrap()).unwrap();
        let options = PreviewOptions {
            tolerance: tolerance.parse().unwrap(),
            ..PreviewOptions::default()
        };
        preview_plan(&plan, &options, &mut Decider::default())
            .map(|mesh| mesh.unresolved_cells)
            .map_err(|refusal| refusal.code)
    }

    #[test]
    fn only_a_buried_carrier_refuses_a_capsule_profile() {
        let cut = [
            plate().to_vec(),
            slot().to_vec(),
            vec![boolean("cut", "subtract", "plate", "slot")],
        ]
        .concat();
        assert_eq!(
            preview(cut.clone(), "1/2"),
            Err("field.display-trim-obstruction")
        );
        // Beyond half the capsule width no cube fits inside it, so nothing is proved.
        assert_eq!(preview(slot().to_vec(), "3"), Ok(0));
        // The carrier stays, but the cap is real material here.
        let kept = [plate(), slot()].concat();
        let kept = [kept, vec![boolean("kept", "intersect", "plate", "slot")]].concat();
        assert_eq!(preview(kept, "1/2"), Ok(0));
        // The cut is filled again, so the operand no longer decides the cap.
        let mut fill = extrude("fill", plate()[0]["loops"].clone()).to_vec();
        fill.push(boolean("whole", "union", "cut", "fill"));
        assert_eq!(preview([cut, fill].concat(), "1/2"), Ok(0));
    }
}

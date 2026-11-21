//! Fillets at a mixed vertex (4.2a step 4, part 3): the corpus row
//! `l-mixed-all`.
//!
//! RECIPE section 4.5.1 settles a vertex where one selected concave edge and
//! two selected convex edges meet: the concave blend forms first and the
//! convex rounds roll over it. The concave blend leaves an arc against the
//! face `P` the two convex edges share, tangent to both convex edges, and the
//! convex rounds continue around that arc. The patch that rolls over the
//! concave blend is the vertex's junction sector, named `corner.vertex`.
//!
//! The field builds that patch from its own pieces. The concave strip between
//! the faces `F1` and `F2` keeps its region. Each convex strip between `P`
//! and `F1` or `F2` stops where the concave strip's foot on that face crosses
//! `P`: from there its ball rolls on `P` and the concave blend's cylinder,
//! whose axis is the concave strip's centre line and whose radius is the ball
//! radius. With `P` normal to that line, this is a ring piece (`Revolved`)
//! about the centre line, cut to the sector between the meridians through the
//! two feet. Those meridian half-planes are junctions: at each, the ring's
//! ball is the convex strip's ball, because the concave centre line lies `r`
//! from `F1` on the air side, so the ring's centre, `2r` from the axis through
//! the foot, lies `r` inside `F1` and `r` from `P`.
//!
//! Why the ring's fan holds only material once the concave strip is added.
//! The ring needs `P` normal to the concave line, and every other support
//! clear of its box except `F1` and `F2` (checked when the run is applied).
//! Then inside the box only `P`, `F1` and `F2` can bound the body, so the body
//! there is a union of the octants they cut at the vertex. Eight rays from the
//! vertex, one into each octant, prove it is the concave corner: material on
//! `P`'s ball side and on the material side of `F1` or `F2`. The fan lies on
//! `P`'s ball side and outside the concave cylinder, and the sector between
//! the meridians through the feet is the wedge between the concave strip's
//! rail planes. So a fan point on the air side of both `F1` and `F2` lies in
//! the concave strip's region (the run must reach `r` past `P`), and every
//! other fan point is material. The ring then removes its fan less the balls,
//! as a convex piece does.
//!
//! How the reference kernels do this. OpenCASCADE fillets all the edges at
//! once and its result for `l-mixed-all` equals the concave edge filleted
//! first and every other edge after (raw volumes within 8e-10 mm^3, recorded
//! with RECIPE 4.5.1). SMLib's `SmFilletExecutive::ProcessSmoothEdgeSequences`
//! joins smoothly meeting filleted edges into one sequence, as part 2 notes.
//! The field follows the recipe's order and cuts the convex sequence into
//! pieces with closed forms: strips and the ring sector, so nothing is walked.

use catallax_features::fillet::FilletDefinition;
use catallax_numbers::{Decider, Rational, Real, Refusal, Sign};
use catallax_recipe::{Origin, Role, StepId, VertexName};
use catallax_shapes::{Curve, Cylinder, Line, Point3, Surface, Vector3};

use super::rolling::Run;
use super::{bad, coordinates, q};
use crate::rolling_region::generator::Generator;
use crate::rolling_region::piece::{Kind, Piece};
use crate::rolling_region::revolved::{Revolution, Revolved};
use crate::{Location, Plan, geometry};

/// A mixed vertex the field builds.
pub(super) struct Mixed {
    /// The endpoint at the vertex.
    pub(super) at: usize,
    /// Each convex run, its end at the vertex and the station of the
    /// junction that replaces that end.
    pub(super) junctions: [(usize, usize, Rational); 2],
    /// The ring sector that rolls over the concave blend.
    pub(super) run: Run,
}

type Definitions = [(FilletDefinition, Option<[Rational; 2]>)];
/// A straight run's line, interval and the end at a vertex.
type RunEnd<'a> = (&'a Line, &'a [Rational; 2], usize);

fn same(a: &Point3, b: &Point3, dec: &mut Decider) -> Result<bool, Refusal> {
    Ok(geometry::sign(&a.minus(b).norm_squared(), dec)? == Sign::Zero)
}

/// The end of a straight run at `vertex`, with its line and interval.
fn end_at<'a>(
    (definition, interval): &'a (FilletDefinition, Option<[Rational; 2]>),
    vertex: &Point3,
    dec: &mut Decider,
) -> Result<Option<RunEnd<'a>>, Refusal> {
    let (Curve::Line(line), Some(interval)) = (&definition.edge, interval) else {
        return Ok(None);
    };
    for (k, t) in interval.iter().enumerate() {
        if same(
            &line.origin.plus(&line.dir.scaled(&Real::from(t.clone()))),
            vertex,
            dec,
        )? {
            return Ok(Some((line, interval, k)));
        }
    }
    Ok(None)
}

impl Plan {
    /// The mixed vertices among `endpoints`: one concave and two convex
    /// straight runs between planes, the convex runs sharing a plane `P`
    /// normal to the concave line and each sharing its other plane with the
    /// concave run. Other meetings are left to the trihedral corners.
    pub(super) fn mixed_vertices(
        &self,
        root: usize,
        definitions: &Definitions,
        endpoints: &[(Point3, Vec<usize>)],
        step: &StepId,
        dec: &mut Decider,
    ) -> Result<Vec<Mixed>, Refusal> {
        let mut found = Vec::new();
        for (at, (vertex, meeting)) in endpoints.iter().enumerate() {
            if meeting.len() != 3 {
                continue;
            }
            let (concave, convex): (Vec<usize>, Vec<usize>) =
                meeting.iter().partition(|&&i| !definitions[i].0.convex);
            let ([concave], [first, second]) = (concave.as_slice(), convex.as_slice()) else {
                continue;
            };
            let (concave, first, second) = (*concave, *first, *second);
            if let Some(mixed) = self.mixed_vertex(
                root,
                at,
                vertex,
                concave,
                [first, second],
                definitions,
                step,
                dec,
            )? {
                found.push(mixed);
            }
        }
        Ok(found)
    }

    /// The source face of a support.
    fn face_of(&self, support: &Surface) -> usize {
        self.faces
            .iter()
            .position(|f| f.surface.same_definition(support))
            .expect("source support")
    }

    #[allow(clippy::too_many_arguments)]
    fn mixed_vertex(
        &self,
        root: usize,
        at: usize,
        vertex: &Point3,
        concave: usize,
        convex: [usize; 2],
        definitions: &Definitions,
        step: &StepId,
        dec: &mut Decider,
    ) -> Result<Option<Mixed>, Refusal> {
        let c = &definitions[concave].0;
        let planes = [c, &definitions[convex[0]].0, &definitions[convex[1]].0]
            .iter()
            .all(|d| d.supports.iter().all(|s| matches!(s, Surface::Plane(_))));
        let Some((c_line, c_interval, c_end)) = end_at(&definitions[concave], vertex, dec)? else {
            return Ok(None);
        };
        if !planes {
            return Ok(None);
        }
        // Each convex run shares one support with the concave run; the other
        // is the face `P` both convex runs share.
        let mut shared = [0; 2];
        let mut end_face: Option<&Surface> = None;
        for (slot, &index) in shared.iter_mut().zip(&convex) {
            let supports = &definitions[index].0.supports;
            let Some((k, j)) = (0..2)
                .flat_map(|k| (0..2).map(move |j| (k, j)))
                .find(|&(k, j)| supports[k].same_definition(&c.supports[j]))
            else {
                return Ok(None);
            };
            *slot = j;
            let other = &supports[1 - k];
            match end_face {
                None => end_face = Some(other),
                Some(p) if p.same_definition(other) => {}
                Some(_) => return Ok(None),
            }
        }
        let Some(end_face) = end_face else {
            return Ok(None);
        };
        let Surface::Plane(p) = end_face else {
            return Ok(None);
        };
        if shared[0] == shared[1]
            || geometry::sign(&p.normal.cross(&c_line.dir).norm_squared(), dec)? != Sign::Zero
        {
            return Ok(None);
        }
        let r = &c.radius;
        let [f1, f2] = c.supports.each_ref().map(|s| self.face_of(s));
        let fp = self.face_of(end_face);

        // The concave run, derived by the field: its centre line and feet.
        let midpoint = |line: &Line, interval: &[Rational; 2]| {
            line.origin.plus(&line.dir.scaled(&Real::from(
                (&interval[0] + &interval[1]) / Rational::from_i64(2),
            )))
        };
        let c_sample = midpoint(c_line, c_interval);
        let c_ball = self.ball_direction([f1, f2], &c_sample, false, c.same_sense)?;
        let walls = Self::generator_walls(&c.supports, &c_ball, &c_sample, dec)?;
        let concave_path = Generator::new(
            walls,
            c_line.origin.clone(),
            c_line.dir.clone(),
            c_interval.clone().map(Real::from),
            r.clone(),
            dec,
        )?;
        let at_vertex = Real::from(c_interval[c_end].clone());
        let [edge, foot_1, foot_2] = concave_path.section_at(&at_vertex);
        let axis_point = concave_path.centre_at(&at_vertex);
        if !same(&edge, vertex, dec)? {
            return Err(bad("a mixed vertex is not on its concave run"));
        }
        // The ring's fan reaches `r` along the concave run.
        let length = &c_interval[1] - &c_interval[0];
        let reach = &(&Real::from(&length * &length) * &c_line.dir.norm_squared()) - &r.square();
        if geometry::sign(&reach, dec)? == Sign::Negative {
            return Err(bad(
                "a concave run shorter than the ball ends at a mixed vertex",
            ));
        }

        // Each convex run now ends where the concave foot on its face
        // crosses `P`.
        let feet = [&foot_1, &foot_2];
        let mut junctions = Vec::with_capacity(2);
        let mut convex_ball = None;
        for (&index, &j) in convex.iter().zip(&shared) {
            let Some((line, interval, k)) = end_at(&definitions[index], vertex, dec)? else {
                return Ok(None);
            };
            let foot = feet[j];
            if geometry::sign(
                &foot.minus(&line.origin).cross(&line.dir).norm_squared(),
                dec,
            )? != Sign::Zero
            {
                return Err(bad(
                    "a concave foot at a mixed vertex leaves the convex edge",
                ));
            }
            let station =
                q(&(&foot.minus(&line.origin).dot(&line.dir) / &line.dir.norm_squared()))?;
            if station <= interval[0] || station >= interval[1] {
                return Err(bad(
                    "a convex run shorter than the concave blend ends at a mixed vertex",
                ));
            }
            if convex_ball.is_none() {
                let definition = &definitions[index].0;
                let faces = definition.supports.each_ref().map(|s| self.face_of(s));
                convex_ball = Some(self.ball_direction(
                    faces,
                    &midpoint(line, interval),
                    true,
                    definition.same_sense,
                )?);
            }
            junctions.push((index, k, station));
        }
        let convex_ball = convex_ball.expect("two convex runs");

        self.concave_corner(root, vertex, p, &c.supports, &convex_ball, &c_ball, dec)?;

        // The ring sector between the meridians through the two feet.
        let flat = Self::flat(end_face, &convex_ball, dec)?.expect("a plane");
        let inverse = &Real::one() / r;
        let mut ends = feet.map(|foot| foot.minus(&axis_point).scaled(&inverse));
        if geometry::sign(&c_line.dir.dot(&ends[0].cross(&ends[1])), dec)? == Sign::Negative {
            ends.swap(0, 1);
        }
        let [start, end] = ends;
        let ring = Revolved::new(
            flat,
            Revolution::Cylinder {
                point: axis_point.clone(),
                axis: c_line.dir.clone(),
                radius: r.clone(),
            },
            true,
            r.clone(),
            dec,
        )?
        .with_sector(start.clone(), end, &c_line.dir, dec)?;
        let cylinder = Surface::Cylinder(Cylinder {
            origin: axis_point,
            axis: c_line.dir.clone(),
            radius_sq: r.square(),
            ref_dir: start,
        });
        let mut names: Vec<Origin> = [fp, f1, f2]
            .iter()
            .map(|&i| self.faces[i].origin.clone())
            .collect();
        names.sort();
        names.dedup();
        let vertex_name = VertexName::new(names).map_err(|e| bad(&e.to_string()))?;
        let piece = Piece::new(
            Kind::Revolved(ring),
            Origin {
                step: step.clone(),
                role: Role::Corner {
                    vertex: vertex_name,
                },
            },
            true,
            [end_face, &cylinder],
            dec,
        )?;
        let [first, second]: [(usize, usize, Rational); 2] =
            junctions.try_into().map_err(|_| bad("two convex runs"))?;
        Ok(Some(Mixed {
            at,
            junctions: [first, second],
            run: Run::new(piece, [fp, fp], vec![f1, f2]),
        }))
    }

    /// Eight rays from the vertex, one into each octant that `P`, `F1` and
    /// `F2` cut, find material exactly on `P`'s ball side and on the
    /// material side of `F1` or `F2`: the concave corner of `F1` and `F2`
    /// cut by `P`.
    #[allow(clippy::too_many_arguments)]
    fn concave_corner(
        &self,
        root: usize,
        vertex: &Point3,
        p: &catallax_shapes::Plane,
        concave_supports: &[Surface; 2],
        convex_ball: &Vector3,
        concave_ball: &Vector3,
        dec: &mut Decider,
    ) -> Result<(), Refusal> {
        let corner = || bad("the body at a mixed vertex is not a concave corner cut by one face");
        let [Surface::Plane(a), Surface::Plane(b)] = concave_supports else {
            return Err(corner());
        };
        let xyz = coordinates(vertex)?;
        let normals = [&p.normal, &a.normal, &b.normal];
        if geometry::sign(&normals[0].dot(&normals[1].cross(normals[2])), dec)? == Sign::Zero {
            return Err(corner());
        }
        let duals = [
            normals[1].cross(normals[2]),
            normals[2].cross(normals[0]),
            normals[0].cross(normals[1]),
        ];
        // The material side of `P` holds the convex ball; that of `F1` and
        // `F2` lies away from the concave ball.
        let toward = [
            geometry::sign(&normals[0].dot(convex_ball), dec)?,
            geometry::sign(&normals[1].dot(concave_ball), dec)?,
            geometry::sign(&normals[2].dot(concave_ball), dec)?,
        ];
        for signs in 0..8u8 {
            let mut direction = Vector3::zero();
            for (k, dual) in duals.iter().enumerate() {
                let s = if signs & (1 << k) == 0 { 1 } else { -1 };
                direction = direction.plus(&dual.scaled(&Real::from_i64(s)));
            }
            let mut side = [Sign::Zero; 3];
            for (k, normal) in normals.iter().enumerate() {
                side[k] = geometry::sign(&normal.dot(&direction), dec)?;
            }
            let material = side[0] == toward[0] && (side[1] != toward[1] || side[2] != toward[2]);
            let found = match self.fillet_ray(root, &xyz, &direction, dec)? {
                Location::Inside => true,
                Location::Outside => false,
                _ => return Err(corner()),
            };
            if found != material {
                return Err(corner());
            }
        }
        Ok(())
    }
}

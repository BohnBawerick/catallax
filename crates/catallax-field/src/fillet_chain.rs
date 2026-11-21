//! Fillet runs that end at a tangent neighbour or run past their rails
//! (4.2a step 4, part 2): the slot file.
//!
//! Two kinds of run end are taken here, and only where the other routes
//! refuse, so every earlier answer is unchanged.
//!
//! A tangent junction: two selected runs share one face and meet at a point
//! where their other faces are tangent, with the same material side, and
//! leave it in opposite directions. Their balls there are the same ball, so
//! their sections in the plane normal to both runs at that point are the
//! same, and that plane is a junction of the two pieces: a straight piece
//! ends there with an end plane normal to its line, a circle piece with a
//! meridian half-plane (a sector of at most half a turn). Recipe version 4
//! continues a selected run across each such junction, as far as the faces
//! stay tangent; version 3 rounds only the runs it names.
//!
//! A run past its rails (recipe version 3): a straight convex run whose ends
//! meet cylindrical end walls runs on past them and is cut by them. The
//! field's piece is the strip over the whole line, past every leaf, with the
//! end walls exempt from its clearance. Each end wall's face must lie
//! beyond its run end, and the ball at each run end must lie inside the
//! wall's solid on the side that holds the ball's centre: a ball that pokes
//! through the end wall does not fit there.
//!
//! A face that meets a piece's box near a junction is cleared by its leaf:
//! an extrude's side face lies in its profile segment swept along the
//! extrusion, so it is shown to lie on the far side of the junction, where
//! it can meet the piece's region only in the junction plane.
//!
//! How the reference kernels do this. OpenCASCADE's
//! `ChFi3d_Builder::PerformElement` (V7_8_0) collects "all mutually tangent
//! edges" into one spine, requiring the faces opposite two adjacent edges to
//! be tangent, and walks the blend along that spine. SMLib's
//! `SmFilletExecutive::ProcessSmoothEdgeSequences` collects smoothly meeting
//! filleted edges (`CollectTangentEdgeSequence`) and joins their curves into
//! one, and adjusts the tangent corners between them (`AdjustTangentCorner`).
//! The field keeps the same chain rule but cuts the chain at each tangent
//! junction into pieces with their own closed forms, a straight strip and a
//! ring sector, which meet in one section, so no walk is needed.

use catallax_numbers::{Decider, Rational, Real, Refusal, Sign};
use catallax_recipe::{EdgeSelection, StepId};
use catallax_shapes::{Circle, Curve, GeomStore, Line, Point3, Surface, Vector3};

use super::rolling::{Run, oriented, revolution};
use super::{coordinates, square_root, unit};
use crate::rolling_region::generator::Generator;
use crate::rolling_region::piece::{Kind, Piece};
use crate::rolling_region::revolved::Revolved;
use crate::{Location, Plan, geometry};

/// Where a span runs.
#[derive(Clone, Debug)]
enum Track {
    /// The stations `t0 < t1` of a line `o + t d`.
    Line { line: Line, stations: [Rational; 2] },
    /// The arc of a circle from the unit radial direction `ends[0]` to
    /// `ends[1]`, turning positively about the circle's normal.
    Arc { circle: Circle, ends: [Vector3; 2] },
}

/// One sharp edge run of a selected support pair.
#[derive(Clone, Debug)]
struct Span {
    faces: [usize; 2],
    convex: bool,
    senses: [bool; 2],
    // An edge point inside the run, for the ball's side.
    sample: Point3,
    track: Track,
}

impl Span {
    fn end(&self, k: usize) -> Point3 {
        match &self.track {
            Track::Line { line, stations } => line
                .origin
                .plus(&line.dir.scaled(&Real::from(stations[k].clone()))),
            Track::Arc { circle, ends } => circle
                .center
                .plus(&ends[k].scaled(&square_root(&circle.radius_sq))),
        }
    }
    /// The run's tangent at end `k`, pointing out of the run.
    fn outward(&self, k: usize) -> Vector3 {
        let forward = match &self.track {
            Track::Line { line, .. } => line.dir.clone(),
            Track::Arc { circle, ends } => circle.normal.cross(&ends[k]),
        };
        if k == 0 {
            forward.scaled(&Real::from_i64(-1))
        } else {
            forward
        }
    }
}

/// One arc between neighbouring cuts: its material, its first and last
/// cut, and its midpoint.
type ArcProbe = (Option<(bool, [bool; 2])>, usize, usize, Point3);

fn same(a: &Point3, b: &Point3, dec: &mut Decider) -> Result<bool, Refusal> {
    Ok(geometry::sign(&a.minus(b).norm_squared(), dec)? == Sign::Zero)
}
fn parallel(a: &Vector3, b: &Vector3, dec: &mut Decider) -> Result<bool, Refusal> {
    Ok(geometry::sign(&a.cross(b).norm_squared(), dec)? == Sign::Zero)
}
fn unsupported(detail: &str) -> Refusal {
    Refusal::new("field.fillet-contact", detail)
}

/// Cyclic order of radial vectors about a circle: the half-turn each starts
/// in, from its first frame axis, then the turn between them.
fn before(a: &[Real; 2], b: &[Real; 2], dec: &mut Decider) -> Result<bool, Refusal> {
    let half = |v: &[Real; 2], dec: &mut Decider| -> Result<u8, Refusal> {
        Ok(match geometry::sign(&v[1], dec)? {
            Sign::Positive => 0,
            Sign::Negative => 1,
            Sign::Zero => u8::from(geometry::sign(&v[0], dec)? == Sign::Negative),
        })
    };
    let (ha, hb) = (half(a, dec)?, half(b, dec)?);
    if ha != hb {
        return Ok(ha < hb);
    }
    let turn = &(&a[0] * &b[1]) - &(&a[1] * &b[0]);
    Ok(geometry::sign(&turn, dec)? == Sign::Positive)
}

/// Whether the chain route may take a fillet the other routes refused: a
/// refusal about the fillet's shape, not about exhausted work.
pub(super) fn retried(refusal: &Refusal) -> bool {
    matches!(
        refusal.code,
        "fillet.unsupported-overflow" | "fillet.unsupported-contact" | "field.fillet-contact"
    )
}

impl Plan {
    /// Builds the selected fillet from junction and run-past pieces, when
    /// every selected run ends at a tangent junction, or, in version 3, every
    /// selected run is a straight convex run past its rails. `None` when the
    /// selection is of neither kind; the caller then keeps its refusal.
    #[allow(clippy::too_many_arguments)]
    pub(super) fn chain_fillet(
        &mut self,
        root: usize,
        source_faces: &[usize],
        selection: &EdgeSelection,
        radius: &Real,
        roll_over: bool,
        step: &StepId,
        dec: &mut Decider,
    ) -> Result<Option<usize>, Refusal> {
        let EdgeSelection::Refs(refs) = selection else {
            return Ok(None);
        };
        // A body an earlier fillet or chamfer rounded is step 5's: rolling
        // over earlier rounds is not a tangent junction of this route.
        for &i in source_faces {
            match &self.leaves[self.faces[i].leaf] {
                geometry::Leaf::Feature(feature)
                    if feature.fillet.is_none() && feature.chamfer.is_none() => {}
                geometry::Leaf::Primitive(_) => {}
                _ => return Ok(None),
            }
        }
        let mut pairs = Vec::new();
        for reference in refs {
            if reference.near.is_some() || !reference.touching.is_empty() {
                return Ok(None);
            }
            for (position, &a) in source_faces.iter().enumerate() {
                for &b in &source_faces[position + 1..] {
                    if ((reference.between[0].matches(&self.faces[a].origin)
                        && reference.between[1].matches(&self.faces[b].origin))
                        || (reference.between[1].matches(&self.faces[a].origin)
                            && reference.between[0].matches(&self.faces[b].origin)))
                        && !pairs.contains(&[a, b])
                    {
                        pairs.push([a, b]);
                    }
                }
            }
        }
        let mut bound = Rational::one();
        for leaf in &self.leaves {
            if let Some(r) = leaf.display_radius(dec)? {
                bound = bound + r;
            }
        }
        let mut spans: Vec<Span> = Vec::new();
        for faces in pairs {
            let Some(found) = self.spans(root, faces, source_faces, &bound, dec)? else {
                return Ok(None);
            };
            spans.extend(found);
        }
        if spans.is_empty() {
            return Ok(None);
        }
        if roll_over && !self.continue_tangent(root, &mut spans, source_faces, &bound, dec)? {
            return Ok(None);
        }
        // The partner of each run end at a tangent junction.
        let mut links: Vec<[Option<(usize, usize)>; 2]> = vec![[None, None]; spans.len()];
        for i in 0..spans.len() {
            for k in 0..2 {
                let p = spans[i].end(k);
                for j in i + 1..spans.len() {
                    for l in 0..2 {
                        if !same(&p, &spans[j].end(l), dec)? {
                            continue;
                        }
                        if links[i][k].is_some() || links[j][l].is_some() {
                            return Err(unsupported(
                                "three or more selected runs meet at one point",
                            ));
                        }
                        links[i][k] = Some((j, l));
                        links[j][l] = Some((i, k));
                    }
                }
            }
        }
        let chained = links.iter().all(|ends| ends.iter().all(Option::is_some));
        let past = !roll_over
            && links.iter().all(|ends| ends.iter().all(Option::is_none))
            && spans
                .iter()
                .all(|span| span.convex && matches!(span.track, Track::Line { .. }));
        if !chained && !past {
            return Ok(None);
        }
        let mut runs = Vec::with_capacity(spans.len());
        for (i, span) in spans.iter().enumerate() {
            if chained {
                for (k, link) in links[i].iter().enumerate() {
                    let (j, l) = link.expect("a chained end");
                    if j > i {
                        self.check_junction(span, k, &spans[j], l, dec)?;
                    }
                }
                runs.push(self.junction_run(span, radius, step, dec)?);
            } else {
                runs.push(self.past_run(span, source_faces, &bound, radius, step, dec)?);
            }
        }
        let mut result = root;
        for run in runs {
            result = self.apply_run(result, run, source_faces, dec)?;
        }
        Ok(Some(result))
    }

    /// The sharp runs of one support pair: lines and circles only. `None`
    /// for any other intersection.
    fn spans(
        &self,
        root: usize,
        faces: [usize; 2],
        source_faces: &[usize],
        bound: &Rational,
        dec: &mut Decider,
    ) -> Result<Option<Vec<Span>>, Refusal> {
        let supports = faces.map(|i| self.faces[i].surface.clone());
        let mut geom = GeomStore::new();
        let a = geom.add_surface(supports[0].clone());
        let b = geom.add_surface(supports[1].clone());
        let catallax_intersect::SurfacePair::Meets(components) =
            catallax_intersect::intersect_surfaces(&mut geom, a, b, dec)?.value
        else {
            return Ok(Some(Vec::new()));
        };
        let mut spans = Vec::new();
        for component in components {
            let catallax_intersect::Component::Curve { curve, .. } = component else {
                continue;
            };
            match geom.curve(curve).clone() {
                Curve::Line(line) => {
                    let samples =
                        self.line_samples(&line, curve, &mut geom, source_faces, bound, dec)?;
                    let mut found: Vec<Span> = Vec::new();
                    for (sample, interval) in samples {
                        let interval = interval.expect("a line sample has its stations");
                        let Some((convex, senses)) = self.material(root, &sample, faces, dec)?
                        else {
                            continue;
                        };
                        if let Some(last) = found.last_mut()
                            && let Track::Line { stations, .. } = &mut last.track
                            && stations[1] == interval[0]
                            && last.convex == convex
                            && last.senses == senses
                        {
                            stations[1] = interval[1].clone();
                            continue;
                        }
                        found.push(Span {
                            faces,
                            convex,
                            senses,
                            sample,
                            track: Track::Line {
                                line: line.clone(),
                                stations: interval,
                            },
                        });
                    }
                    spans.extend(found);
                }
                Curve::Circle(circle) => {
                    let Some(found) =
                        self.arcs(root, faces, &circle, curve, &mut geom, source_faces, dec)?
                    else {
                        return Ok(None);
                    };
                    spans.extend(found);
                }
                _ => return Ok(None),
            }
        }
        Ok(Some(spans))
    }

    /// The sharp arcs of a circle where two supports meet: the circle cut at
    /// every point where another source support crosses or touches it, each
    /// arc probed for its material at its midpoint, and neighbouring arcs of
    /// the same material joined. `None` when the whole circle is one run,
    /// which has no ends.
    #[allow(clippy::too_many_arguments)]
    fn arcs(
        &self,
        root: usize,
        faces: [usize; 2],
        circle: &Circle,
        curve: catallax_shapes::CurveId,
        geom: &mut GeomStore,
        source_faces: &[usize],
        dec: &mut Decider,
    ) -> Result<Option<Vec<Span>>, Refusal> {
        let radius = square_root(&circle.radius_sq);
        let first = circle.ref_dir.clone();
        let second = circle.normal.cross(&circle.ref_dir);
        let mut cuts: Vec<(Vector3, [Real; 2])> = Vec::new();
        for &index in source_faces {
            let face = &self.faces[index];
            if faces.contains(&index)
                || faces
                    .iter()
                    .any(|j| face.support.same_definition(&self.faces[*j].support))
            {
                continue;
            }
            let surface = geom.add_surface(face.surface.clone());
            let catallax_intersect::CurveSurface::Meets(hits) =
                catallax_intersect::intersect_curve_surface(geom, curve, surface, dec)?.value
            else {
                continue;
            };
            for hit in hits {
                let v = geom.position(hit.point).minus(&circle.center);
                let mut known = false;
                for (other, _) in &cuts {
                    if geometry::sign(&v.minus(other).norm_squared(), dec)? == Sign::Zero {
                        known = true;
                        break;
                    }
                }
                if !known {
                    let at = [v.dot(&first), v.dot(&second)];
                    cuts.push((v, at));
                }
            }
        }
        if cuts.is_empty() {
            return Ok(None);
        }
        // Insertion sort by angle; the decider settles each comparison.
        for i in 1..cuts.len() {
            let mut j = i;
            while j > 0 && before(&cuts[j].1, &cuts[j - 1].1, dec)? {
                cuts.swap(j, j - 1);
                j -= 1;
            }
        }
        let mut found: Vec<ArcProbe> = Vec::new();
        for i in 0..cuts.len() {
            let j = (i + 1) % cuts.len();
            let (a, b) = (&cuts[i].0, &cuts[j].0);
            let turn = geometry::sign(&circle.normal.dot(&a.cross(b)), dec)?;
            let middle = match turn {
                Sign::Positive => a.plus(b),
                Sign::Negative => a.plus(b).scaled(&Real::from_i64(-1)),
                // One cut, or two opposite ones.
                Sign::Zero if cuts.len() == 1 => a.scaled(&Real::from_i64(-1)),
                Sign::Zero => circle.normal.cross(a),
            };
            let sample = circle.center.plus(&unit(&middle).scaled(&radius));
            let material = self.material_point(root, &sample, faces, dec)?;
            found.push((material, i, j, sample));
        }
        // Join neighbouring arcs of one material, round the circle.
        let n = found.len();
        let Some(start) = (0..n).find(|&i| found[i].0 != found[(i + n - 1) % n].0) else {
            return Ok(if found[0].0.is_some() {
                None
            } else {
                Some(Vec::new())
            });
        };
        let mut spans = Vec::new();
        let mut i = start;
        loop {
            let (material, from, _, sample) = found[i].clone();
            let mut to = found[i].2;
            let mut k = (i + 1) % n;
            while k != start && found[k].0 == material {
                to = found[k].2;
                k = (k + 1) % n;
            }
            if let Some((convex, senses)) = material {
                let ends = [from, to].map(|c| cuts[c].0.scaled(&(&Real::one() / &radius)));
                spans.push(Span {
                    faces,
                    convex,
                    senses,
                    sample,
                    track: Track::Arc {
                        circle: circle.clone(),
                        ends,
                    },
                });
            }
            i = k;
            if i == start {
                break;
            }
        }
        Ok(Some(spans))
    }

    /// Whether face `i` passes through `p`: its support does and its leaf
    /// does not exclude the point from the face.
    fn through(&self, i: usize, p: &Point3, dec: &mut Decider) -> Result<bool, Refusal> {
        let face = &self.faces[i];
        let leaf = &self.leaves[face.leaf];
        if let Ok(xyz) = coordinates(p) {
            return Ok(leaf.location(&xyz, dec)? == Location::Boundary
                && leaf.on_face(&face.origin, self.definition_role(i), &xyz, dec)? != Some(false));
        }
        Ok(
            crate::distance::location(leaf, p, dec)? == Location::Boundary
                && leaf.on_face_algebraic(&face.origin, self.definition_role(i), p, dec)?
                    != Some(false),
        )
    }

    /// Version 4 continues each run across every end where the face it
    /// leaves is tangent to another face through that end: the run between
    /// the face it keeps and that face joins the selection, as long as the
    /// faces stay tangent. `false` when a continued pair meets in a curve
    /// other than a line or a circle.
    fn continue_tangent(
        &self,
        root: usize,
        spans: &mut Vec<Span>,
        source_faces: &[usize],
        bound: &Rational,
        dec: &mut Decider,
    ) -> Result<bool, Refusal> {
        let mut pending: Vec<usize> = (0..spans.len()).collect();
        while let Some(i) = pending.pop() {
            for k in 0..2 {
                let p = spans[i].end(k);
                let mut linked = false;
                for (j, other) in spans.iter().enumerate() {
                    if j != i && (same(&p, &other.end(0), dec)? || same(&p, &other.end(1), dec)?) {
                        linked = true;
                        break;
                    }
                }
                if linked {
                    continue;
                }
                let [f0, f1] = spans[i].faces;
                for (stay, leave) in [(f0, f1), (f1, f0)] {
                    let gradient = self.gradient_at(leave, &p);
                    for &next in source_faces {
                        let support = &self.faces[next].support;
                        if [stay, leave]
                            .iter()
                            .any(|&f| f == next || support.same_definition(&self.faces[f].support))
                            || !self.through(next, &p, dec)?
                            || !parallel(&gradient, &self.gradient_at(next, &p), dec)?
                        {
                            continue;
                        }
                        let pair = if source_faces.iter().position(|&f| f == stay)
                            < source_faces.iter().position(|&f| f == next)
                        {
                            [stay, next]
                        } else {
                            [next, stay]
                        };
                        let Some(found) = self.spans(root, pair, source_faces, bound, dec)? else {
                            return Ok(false);
                        };
                        for span in found {
                            if !(same(&p, &span.end(0), dec)? || same(&p, &span.end(1), dec)?) {
                                continue;
                            }
                            let mut known = false;
                            for old in spans.iter() {
                                if old.faces == span.faces
                                    && same(&old.end(0), &span.end(0), dec)?
                                    && same(&old.end(1), &span.end(1), dec)?
                                {
                                    known = true;
                                    break;
                                }
                            }
                            if !known {
                                spans.push(span);
                                pending.push(spans.len() - 1);
                            }
                        }
                    }
                }
            }
        }
        Ok(true)
    }

    /// Two runs meet in one ball at a tangent junction: they share one face
    /// with the same material side, are both convex or both concave, their
    /// other faces are tangent there with the same material side, and they
    /// leave the point in opposite directions.
    fn check_junction(
        &self,
        a: &Span,
        k: usize,
        b: &Span,
        l: usize,
        dec: &mut Decider,
    ) -> Result<(), Refusal> {
        let not_tangent = || unsupported("runs that meet at a point need a tangent junction");
        let p = a.end(k);
        let shared: Vec<(usize, usize)> = (0..2)
            .flat_map(|x| (0..2).map(move |y| (x, y)))
            .filter(|&(x, y)| a.faces[x] == b.faces[y])
            .collect();
        let [(x, y)] = shared[..] else {
            return Err(not_tangent());
        };
        if a.convex != b.convex || a.senses[x] != b.senses[y] {
            return Err(not_tangent());
        }
        let (ga, gb) = (
            self.gradient_at(a.faces[1 - x], &p),
            self.gradient_at(b.faces[1 - y], &p),
        );
        if !parallel(&ga, &gb, dec)? {
            return Err(not_tangent());
        }
        let aligned = geometry::sign(&ga.dot(&gb), dec)? == Sign::Positive;
        if aligned != (a.senses[1 - x] == b.senses[1 - y]) {
            return Err(not_tangent());
        }
        let (oa, ob) = (a.outward(k), b.outward(l));
        if !parallel(&oa, &ob, dec)? || geometry::sign(&oa.dot(&ob), dec)? != Sign::Negative {
            return Err(not_tangent());
        }
        Ok(())
    }

    /// The piece of a run between two tangent junctions: a straight piece
    /// whose end planes are junctions, or a ring sector.
    fn junction_run(
        &self,
        span: &Span,
        radius: &Real,
        step: &StepId,
        dec: &mut Decider,
    ) -> Result<Run, Refusal> {
        let supports = span.faces.map(|i| self.faces[i].surface.clone());
        let ball = self.ball_direction(span.faces, &span.sample, span.convex, span.senses)?;
        let junctions: Vec<(Point3, Vector3)> = (0..2)
            .map(|k| (span.end(k), span.outward(k).scaled(&Real::from_i64(-1))))
            .collect();
        let kind = match &span.track {
            Track::Line { line, stations } => {
                let walls = Self::generator_walls(&supports, &ball, &span.sample, dec)?;
                Kind::Generator(
                    Generator::new(
                        walls,
                        line.origin.clone(),
                        line.dir.clone(),
                        stations.clone().map(Real::from),
                        radius.clone(),
                        dec,
                    )?
                    .with_junctions([true, true]),
                )
            }
            Track::Arc { circle, ends } => {
                let (plane, other) = match (
                    Self::flat(&supports[0], &ball, dec)?,
                    Self::flat(&supports[1], &ball, dec)?,
                ) {
                    (Some(plane), None) => (plane, &supports[1]),
                    (None, Some(plane)) => (plane, &supports[0]),
                    _ => {
                        return Err(unsupported(
                            "a ring sector needs a plane and a cylinder, cone or sphere",
                        ));
                    }
                };
                let (support, outside) =
                    revolution(other, &span.sample, &ball).ok_or_else(|| {
                        unsupported("a ring sector needs a plane and a cylinder, cone or sphere")
                    })?;
                let outside =
                    oriented(&outside, "the ball lies in a support's tangent plane", dec)?;
                Kind::Revolved(
                    Revolved::new(plane, support, outside, radius.clone(), dec)?.with_sector(
                        ends[0].clone(),
                        ends[1].clone(),
                        &circle.normal,
                        dec,
                    )?,
                )
            }
        };
        let piece = Piece::new(
            kind,
            self.blend_origin(span.faces, step)?,
            span.convex,
            [&supports[0], &supports[1]],
            dec,
        )?;
        Ok(Run::new(piece, span.faces, Vec::new()).with_junctions(junctions))
    }

    /// A straight convex run past its rails (version 3): the strip over the
    /// whole line, cut by the cylindrical end wall through each run end.
    fn past_run(
        &self,
        span: &Span,
        source_faces: &[usize],
        bound: &Rational,
        radius: &Real,
        step: &StepId,
        dec: &mut Decider,
    ) -> Result<Run, Refusal> {
        let Track::Line { line, stations } = &span.track else {
            unreachable!("a run past its rails is straight");
        };
        let supports = span.faces.map(|i| self.faces[i].surface.clone());
        let ball = self.ball_direction(span.faces, &span.sample, span.convex, span.senses)?;
        let walls = Self::generator_walls(&supports, &ball, &span.sample, dec)?;
        // Past every leaf, as `line_samples` bounds the line.
        let squared = super::q(&line.dir.norm_squared())?;
        let origin = coordinates(&line.origin)?;
        let reach = (bound + &origin.iter().fold(Rational::zero(), |sum, c| sum + c.abs()))
            / squared.min(Rational::one());
        if stations[0] <= -&reach || stations[1] >= reach {
            return Err(unsupported("a run past its rails reaches past the part"));
        }
        let path = Generator::new(
            walls,
            line.origin.clone(),
            line.dir.clone(),
            [-&reach, reach.clone()].map(Real::from),
            radius.clone(),
            dec,
        )?;
        let mut ends = Vec::new();
        for k in 0..2 {
            let p = span.end(k);
            let outward = span.outward(k);
            let level = outward.dot(&p.minus(&Point3::origin()));
            let centre = path.centre_at(&Real::from(stations[k].clone()));
            let mut found = false;
            for &i in source_faces {
                if span.faces.contains(&i) {
                    continue;
                }
                let face = &self.faces[i];
                let values = [p.x.clone(), p.y.clone(), p.z.clone()];
                if geometry::sign(&face.polynomial.eval_real(&values), dec)? != Sign::Zero {
                    continue;
                }
                let Surface::Cylinder(wall) = &face.surface else {
                    return Err(unsupported(
                        "a run past its rails needs cylindrical end walls",
                    ));
                };
                let geometry::Leaf::Feature(feature) = &self.leaves[face.leaf] else {
                    return Err(unsupported(
                        "a run past its rails needs end walls of a swept profile",
                    ));
                };
                if !feature.face_above(self.definition_role(i), &outward, &level, dec)? {
                    return Err(Refusal::new(
                        "fillet.unsupported-overflow",
                        "an end wall reaches back over the run",
                    ));
                }
                // The ball at the run end stays on the side of the wall
                // that holds its centre: inside, `rho + r <= R`; outside,
                // `rho - r >= R`.
                let w = super::rolling::radial(&centre.minus(&wall.origin), &wall.axis);
                let rho2 = w.norm_squared();
                let wall_radius = square_root(&wall.radius_sq);
                let fits = match geometry::sign(&(&rho2 - &wall.radius_sq), dec)? {
                    Sign::Negative => {
                        let room = &wall_radius - radius;
                        geometry::sign(&room, dec)? != Sign::Negative
                            && geometry::sign(&(&room.square() - &rho2), dec)? != Sign::Negative
                    }
                    Sign::Positive => {
                        geometry::sign(&(&rho2 - &(&wall_radius + radius).square()), dec)?
                            != Sign::Negative
                    }
                    Sign::Zero => false,
                };
                if !fits {
                    return Err(Refusal::new(
                        "fillet.unsupported-overflow",
                        "the ball at a run end crosses its end wall",
                    ));
                }
                ends.push(i);
                found = true;
            }
            if !found {
                return Err(unsupported(
                    "a run past its rails needs an end wall at each end",
                ));
            }
        }
        let piece = Piece::new(
            Kind::Generator(path),
            self.blend_origin(span.faces, step)?,
            span.convex,
            [&supports[0], &supports[1]],
            dec,
        )?;
        Ok(Run::new(piece, span.faces, ends))
    }
}

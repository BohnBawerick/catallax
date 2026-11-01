//! Exact plane subdivision, oriented patch winding and retained components.
//!
//! Input adjacency supplies a closed oriented carrier-patch chain, not the
//! output topology. Subdivision, winding and component extraction replace it.

use super::{DraftDefinition, reduced_point, sign};
use catallax_numbers::{Decider, Rational, Real, Refusal, Sign, StrictSign};
use catallax_recipe::Origin;
use catallax_shapes::{Plane, Point3, Surface, Vector3};
use std::collections::{BTreeMap, BTreeSet, VecDeque};

// Carrier and canonical vertex set identify the incident cell/face pairs.
type FacetOwners = BTreeMap<(usize, Vec<usize>), Vec<(usize, usize)>>;
type Classifier = fn(&Point3, &mut Decider) -> Result<bool, Refusal>;
type OriginResolver = fn(&Point3, &Vector3, &mut Decider) -> Result<Vec<Origin>, Refusal>;

/// One oriented input patch; holes have the opposite loop winding.
#[derive(Clone, Debug)]
pub struct Patch {
    /// Outward input carrier.
    pub plane: Plane,
    /// Boundary loops indexing the skin's exact points.
    pub loops: Vec<Vec<usize>>,
    /// Complete source-face origin set.
    pub origins: Vec<Origin>,
}

/// Boundary definitions collected independently by each kernel.
#[derive(Clone, Debug)]
pub struct Skin {
    /// Source positions, shared by incident patches.
    pub vertices: Vec<Point3>,
    /// Oriented physical source faces, not construction-history supports.
    pub patches: Vec<Patch>,
}

/// A connected maximal planar boundary patch.
#[derive(Clone, Debug)]
pub struct PolyFace {
    /// Index of its carrier in `Polyhedron::planes`.
    pub support: usize,
    /// Whether the stored carrier normal points out of material.
    pub same_sense: bool,
    /// Outer loop followed by holes, oriented outwards.
    pub loops: Vec<Vec<usize>>,
    /// Exact convex subdivision retained for independent sampling and trims.
    pub tiles: Vec<Vec<usize>>,
    /// Union of origins on the merged physical patch.
    pub origins: Vec<Origin>,
}

/// One full-dimensional convex term, with a certified strict interior point.
#[derive(Clone, Debug)]
pub struct ConvexTerm {
    /// Carrier indices and outward senses of its bounding halfspaces.
    pub sides: Vec<(usize, bool)>,
    /// A positive average of all its vertices, retained by definition.
    pub witness: Point3,
}

/// A point on a face's carrier, read against that face's convex tiles.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum Foot {
    /// Strictly inside a tile, so on the face.
    Inside,
    /// Strictly outside every tile, so off the closed face.
    Outside,
    /// Unresolved against these tile edges, each a segment of the face.
    Near(Vec<(usize, usize)>),
}

/// One connected regular material component of a plane arrangement.
#[derive(Clone, Debug)]
pub struct Polyhedron {
    /// Exact arrangement carriers; construction-box planes never bound a face.
    pub planes: Vec<Plane>,
    /// Boundary vertices, with no construction-box or unused points.
    pub vertices: Vec<Point3>,
    /// Outward, maximally merged physical faces.
    pub faces: Vec<PolyFace>,
    /// Convex closed cells as oriented carrier intersections. Their union is
    /// the component; internal subdivision planes are not boundary geometry.
    pub terms: Vec<ConvexTerm>,
}

fn value(plane: &Plane, point: &Point3) -> Real {
    point.minus(&plane.origin).dot(&plane.normal)
}
fn failure(code: &'static str, detail: &'static str) -> Refusal {
    Refusal::new(code, detail)
}
fn centre(vertices: &[Point3], cycle: &[usize]) -> Point3 {
    let sum = cycle.iter().fold(Vector3::zero(), |s, i| {
        s.plus(&vertices[*i].minus(&Point3::origin()))
    });
    Point3::origin().plus(&sum.scaled(&(&Real::one() / &Real::from_i64(cycle.len() as i64))))
}
fn area(vertices: &[Point3], cycle: &[usize], normal: &Vector3) -> Real {
    let origin = &vertices[cycle[0]];
    let mut sum = Real::zero();
    for i in 1..cycle.len() - 1 {
        sum = sum
            + vertices[cycle[i]]
                .minus(origin)
                .cross(&vertices[cycle[i + 1]].minus(origin))
                .dot(normal);
    }
    sum
}
fn triple(a: &Plane, b: &Plane, c: &Plane, dec: &mut Decider) -> Result<Option<Point3>, Refusal> {
    let bc = b.normal.cross(&c.normal);
    let determinant = a.normal.dot(&bc);
    if sign(&determinant, dec)? == Sign::Zero {
        return Ok(None);
    }
    let offset = |p: &Plane| p.origin.minus(&Point3::origin()).dot(&p.normal);
    let v = bc
        .scaled(&offset(a))
        .plus(&c.normal.cross(&a.normal).scaled(&offset(b)))
        .plus(&a.normal.cross(&b.normal).scaled(&offset(c)))
        .scaled(&(&Real::one() / &determinant));
    reduced_point(Point3::origin().plus(&v), dec).map(Some)
}

/// Signed winding of a possibly self-intersecting planar chain. Every chain
/// edge is on another arrangement carrier, so a facet interior misses it.
fn patch_winding(
    patch: &Patch,
    vertices: &[Point3],
    point: &Point3,
    dec: &mut Decider,
) -> Result<i32, Refusal> {
    let u = &patch.plane.ref_dir;
    let v = patch.plane.normal.cross(u);
    let coordinate = |i: usize| {
        let d = vertices[i].minus(point);
        [d.dot(u), d.dot(&v)]
    };
    let mut winding = 0;
    for cycle in &patch.loops {
        for i in 0..cycle.len() {
            let a = coordinate(cycle[i]);
            let b = coordinate(cycle[(i + 1) % cycle.len()]);
            let (ay, by) = (sign(&a[1], dec)?, sign(&b[1], dec)?);
            let cross = &a[0] * &b[1] - &a[1] * &b[0];
            let turn = sign(&cross, dec)?;
            if turn == Sign::Zero && sign(&(&a[0] * &b[0] + &a[1] * &b[1]), dec)? != Sign::Positive
            {
                return Err(failure(
                    "draft.unsplit-patch-boundary",
                    "a patch boundary has no separating arrangement carrier",
                ));
            }
            if ay != Sign::Positive && by == Sign::Positive && turn == Sign::Positive {
                winding += 1;
            }
            if ay == Sign::Positive && by != Sign::Positive && turn == Sign::Negative {
                winding -= 1;
            }
        }
    }
    Ok(winding)
}

#[derive(Clone)]
struct Facet {
    support: usize,
    same_sense: bool,
    cycle: Vec<usize>,
}
struct Cell {
    faces: Vec<Facet>,
}
struct Interface {
    facet: Facet,
    from: usize,
    to: usize,
    delta: i32,
    origins: Vec<Origin>,
}

fn loops_from_tiles(tiles: &[Vec<usize>]) -> Result<Vec<Vec<usize>>, Refusal> {
    let mut edges = BTreeSet::new();
    for cycle in tiles {
        for i in 0..cycle.len() {
            let edge = (cycle[i], cycle[(i + 1) % cycle.len()]);
            if !edges.remove(&(edge.1, edge.0)) && !edges.insert(edge) {
                return Err(failure(
                    "draft.contradictory-orientation",
                    "same-sense boundary edges overlap",
                ));
            }
        }
    }
    let mut next = BTreeMap::new();
    for (a, b) in edges {
        if next.insert(a, b).is_some() {
            return Err(failure(
                "draft.self-touching-shell",
                "a planar boundary has a branching wire",
            ));
        }
    }
    let mut loops = Vec::new();
    while let Some((&start, _)) = next.first_key_value() {
        let mut cycle = vec![start];
        let mut at = start;
        loop {
            let Some(to) = next.remove(&at) else {
                return Err(failure(
                    "draft.nonclosing-boundary",
                    "a planar boundary wire is open",
                ));
            };
            if to == start {
                break;
            }
            cycle.push(to);
            at = to;
        }
        if cycle.len() < 3 {
            return Err(failure(
                "draft.nonclosing-boundary",
                "a boundary wire has no area",
            ));
        }
        loops.push(cycle);
    }
    Ok(loops)
}

struct Clipper<'a> {
    support: usize,
    planes: &'a [Plane],
    vertices: &'a mut Vec<Point3>,
    signs: Vec<Sign>,
    edges: BTreeMap<(usize, usize), usize>,
    dec: &'a mut Decider,
}
impl Clipper<'_> {
    fn intersection(&mut self, cell: &Cell, a: usize, b: usize) -> Result<usize, Refusal> {
        let key = (a.min(b), a.max(b));
        if let Some(i) = self.edges.get(&key) {
            return Ok(*i);
        }
        // A crossed edge joins two facets of its convex cell. Their line meets this
        // plane once, and meeting it directly keeps every vertex three planes deep.
        // A chain of interpolations hides the exact zeros of later predicates.
        let mut line = cell
            .faces
            .iter()
            .filter(|face| face.cycle.contains(&a) && face.cycle.contains(&b))
            .map(|face| face.support);
        let (Some(first), Some(second)) = (line.next(), line.next()) else {
            return Err(failure(
                "draft.nonclosing-boundary",
                "a convex cell edge lost its carriers",
            ));
        };
        let support = self.support;
        let Some(point) = triple(
            &self.planes[first],
            &self.planes[second],
            &self.planes[support],
            self.dec,
        )?
        else {
            return Err(failure(
                "draft.nonclosing-boundary",
                "a crossed cell edge is parallel to its clipping plane",
            ));
        };
        let id = self.vertices.len();
        self.vertices.push(point);
        // This is zero by the line/plane construction, not an enclosure guess.
        self.signs.push(Sign::Zero);
        self.edges.insert(key, id);
        Ok(id)
    }
    fn clip(&mut self, cell: &Cell, negative: bool) -> Result<Cell, Refusal> {
        let strict = if negative {
            Sign::Negative
        } else {
            Sign::Positive
        };
        let mut faces = Vec::new();
        let mut cap_edges = Vec::new();
        for face in &cell.faces {
            if !face.cycle.iter().any(|i| self.signs[*i] == strict) {
                continue;
            }
            let mut cycle = Vec::with_capacity(face.cycle.len() + 1);
            for i in 0..face.cycle.len() {
                let (a, b) = (face.cycle[i], face.cycle[(i + 1) % face.cycle.len()]);
                let (sa, sb) = (self.signs[a], self.signs[b]);
                let a_kept = sa == strict || sa == Sign::Zero;
                let b_kept = sb == strict || sb == Sign::Zero;
                if a_kept {
                    cycle.push(a);
                }
                if sa != Sign::Zero && sb != Sign::Zero && a_kept != b_kept {
                    cycle.push(self.intersection(cell, a, b)?);
                }
            }
            for i in 0..cycle.len() {
                let (a, b) = (cycle[i], cycle[(i + 1) % cycle.len()]);
                if self.signs[a] == Sign::Zero && self.signs[b] == Sign::Zero {
                    cap_edges.push((b, a));
                }
            }
            faces.push(Facet {
                support: face.support,
                same_sense: face.same_sense,
                cycle,
            });
        }
        let mut next = BTreeMap::new();
        for (a, b) in cap_edges {
            if next.insert(a, b).is_some() {
                return Err(failure(
                    "draft.nonclosing-boundary",
                    "convex clipping cap branches",
                ));
            }
        }
        let Some((&start, _)) = next.first_key_value() else {
            return Err(failure(
                "draft.nonclosing-boundary",
                "convex clipping cap is empty",
            ));
        };
        let mut cap = vec![start];
        let mut at = start;
        loop {
            let Some(to) = next.remove(&at) else {
                return Err(failure(
                    "draft.nonclosing-boundary",
                    "convex clipping cap is open",
                ));
            };
            if to == start {
                break;
            }
            cap.push(to);
            at = to;
        }
        if !next.is_empty() || cap.len() < 3 {
            return Err(failure(
                "draft.nonclosing-boundary",
                "convex clipping cap is disconnected",
            ));
        }
        faces.push(Facet {
            support: self.support,
            same_sense: negative,
            cycle: cap,
        });
        Ok(Cell { faces })
    }
}

fn root(parent: &mut [usize], mut i: usize) -> usize {
    while parent[i] != i {
        parent[i] = parent[parent[i]];
        i = parent[i];
    }
    i
}
fn join(parent: &mut [usize], a: usize, b: usize) {
    let (a, b) = (root(parent, a), root(parent, b));
    parent[a.max(b)] = a.min(b);
}

fn cell_witness(vertices: &[Point3], cell: &Cell, dec: &mut Decider) -> Result<Point3, Refusal> {
    let (sum, count) = cell
        .faces
        .iter()
        .flat_map(|f| &f.cycle)
        .fold((Vector3::zero(), 0_usize), |(sum, n), v| {
            (sum.plus(&vertices[*v].minus(&Point3::origin())), n + 1)
        });
    reduced_point(
        Point3::origin().plus(&sum.scaled(&(&Real::one() / &Real::from_i64(count as i64)))),
        dec,
    )
}

impl Skin {
    /// Exact membership of a carrier-interior point in a source patch.
    /// Arrangement facet interiors cannot lie on the patch's trim edges.
    pub fn contains_patch(
        &self,
        patch: &Patch,
        point: &Point3,
        dec: &mut Decider,
    ) -> Result<bool, Refusal> {
        Ok(sign(&value(&patch.plane, point), dec)? == Sign::Zero
            && patch_winding(patch, &self.vertices, point, dec)? != 0)
    }

    /// Whether every boundary vertex satisfies every material halfspace.
    pub fn is_convex(&self, dec: &mut Decider) -> Result<bool, Refusal> {
        for point in &self.vertices {
            for patch in &self.patches {
                if sign(&value(&patch.plane, point), dec)? == Sign::Positive {
                    return Ok(false);
                }
            }
        }
        Ok(true)
    }

    /// Resolve physical patches before expanding any coactive aliases.
    /// Near points compare exact distances to closed trims, not carriers.
    pub fn resolve_faces(
        &self,
        references: &[catallax_recipe::FaceRef],
        dec: &mut Decider,
    ) -> Result<BTreeSet<usize>, Refusal> {
        let mut neighbours = Vec::new();
        if references
            .iter()
            .any(|reference| !reference.touching.is_empty())
        {
            neighbours.resize_with(self.patches.len(), Vec::new);
            let mut edges = BTreeMap::<_, Vec<usize>>::new();
            for (index, patch) in self.patches.iter().enumerate() {
                for cycle in &patch.loops {
                    for i in 0..cycle.len() {
                        let (a, b) = (cycle[i], cycle[(i + 1) % cycle.len()]);
                        edges.entry((a.min(b), a.max(b))).or_default().push(index);
                    }
                }
            }
            for uses in edges.values() {
                for &a in uses {
                    for &b in uses {
                        if a != b && !neighbours[a].contains(&b) {
                            neighbours[a].push(b);
                        }
                    }
                }
            }
        }
        let mut selected = BTreeSet::new();
        for reference in references {
            let matches: Vec<_> = self
                .patches
                .iter()
                .enumerate()
                .filter_map(|(i, patch)| {
                    (patch
                        .origins
                        .iter()
                        .any(|origin| reference.face.matches(origin))
                        && reference.touching.iter().all(|pattern| {
                            neighbours[i].iter().any(|j| {
                                self.patches[*j]
                                    .origins
                                    .iter()
                                    .any(|origin| pattern.matches(origin))
                            })
                        }))
                    .then_some(i)
                })
                .collect();
            if matches.is_empty() {
                return Err(Refusal::new(
                    "recipe.missing-reference",
                    format!("no input face resolves {reference:?}"),
                ));
            }
            if let Some(near) = &reference.near
                && matches.len() > 1
            {
                let [x, y, z] = near.0.each_ref().map(|n| Real::from(n.value().clone()));
                let point = Point3::new(x, y, z);
                let mut nearest: Option<(usize, Real)> = None;
                let mut tied = false;
                for index in matches {
                    let distance =
                        self.patch_distance_squared(&self.patches[index], &point, dec)?;
                    let order = if let Some((_, previous)) = &nearest {
                        sign(&(&distance - previous), dec)?
                    } else {
                        Sign::Negative
                    };
                    match order {
                        Sign::Negative => {
                            nearest = Some((index, distance));
                            tied = false;
                        }
                        Sign::Zero => tied = true,
                        Sign::Positive => (),
                    }
                }
                if tied {
                    return Err(Refusal::new(
                        "recipe.ambiguous-reference",
                        "near point has equal minimum distance to several trimmed faces",
                    ));
                }
                selected.insert(nearest.expect("nonempty face matches").0);
            } else {
                selected.extend(matches);
            }
        }
        Ok(selected)
    }

    fn patch_distance_squared(
        &self,
        patch: &Patch,
        point: &Point3,
        dec: &mut Decider,
    ) -> Result<Real, Refusal> {
        let norm = patch.plane.normal.norm_squared();
        if sign(&norm, dec)? != Sign::Positive {
            return Err(failure(
                "draft.zero-face-normal",
                "face carrier has a zero normal",
            ));
        }
        let offset = value(&patch.plane, point);
        let plane_distance = &offset.square() / &norm;
        let projected = point.plus(&patch.plane.normal.scaled(&(-offset / norm)));
        let mut edge_distance: Option<Real> = None;
        for cycle in &patch.loops {
            for i in 0..cycle.len() {
                let a = &self.vertices[cycle[i]];
                let edge = self.vertices[cycle[(i + 1) % cycle.len()]].minus(a);
                let delta = point.minus(a);
                let norm = edge.norm_squared();
                let mut parameter = if sign(&norm, dec)? == Sign::Zero {
                    Real::zero()
                } else {
                    delta.dot(&edge) / norm
                };
                if sign(&parameter, dec)? == Sign::Negative {
                    parameter = Real::zero();
                } else if sign(&(&parameter - &Real::one()), dec)? == Sign::Positive {
                    parameter = Real::one();
                }
                let distance = delta.minus(&edge.scaled(&parameter)).norm_squared();
                let replace = match &edge_distance {
                    None => true,
                    Some(previous) => sign(&(&distance - previous), dec)? == Sign::Negative,
                };
                if replace {
                    edge_distance = Some(distance);
                }
            }
        }
        let edge_distance = edge_distance.ok_or_else(|| {
            failure(
                "draft.unsupported-arrangement",
                "face has no finite trimmed boundary",
            )
        })?;
        if sign(&(&edge_distance - &plane_distance), dec)? == Sign::Zero
            || self.contains_patch(patch, &projected, dec)?
        {
            Ok(plane_distance)
        } else {
            Ok(edge_distance)
        }
    }

    /// Rebuilds all candidate cells, discards swallowed patches and separates
    /// connected interiors. Neutral-plane or explicit retention is applied by
    /// the caller to these components, never to individual convex pieces.
    pub fn draft(
        &self,
        selected: &BTreeSet<usize>,
        definition: &DraftDefinition,
        dec: &mut Decider,
    ) -> Result<Vec<Polyhedron>, Refusal> {
        self.arrange::<Classifier, OriginResolver>(Some((selected, definition)), None, dec)
    }

    /// Extracts the physical boundary of a bounded regular solid whose carriers
    /// occur in this skin. The classifier labels full-dimensional cells; the
    /// resolver supplies surviving input origins on each outward boundary tile.
    /// Construction-box faces are never allowed to become output faces.
    pub fn partition<C, R>(
        &self,
        classify: C,
        resolve: R,
        dec: &mut Decider,
    ) -> Result<Vec<Polyhedron>, Refusal>
    where
        C: FnMut(&Point3, &mut Decider) -> Result<bool, Refusal>,
        R: FnMut(&Point3, &Vector3, &mut Decider) -> Result<Vec<Origin>, Refusal>,
    {
        self.arrange(None, Some((classify, resolve)), dec)
    }

    fn arrange<C, R>(
        &self,
        replacement: Option<(&BTreeSet<usize>, &DraftDefinition)>,
        mut classified: Option<(C, R)>,
        dec: &mut Decider,
    ) -> Result<Vec<Polyhedron>, Refusal>
    where
        C: FnMut(&Point3, &mut Decider) -> Result<bool, Refusal>,
        R: FnMut(&Point3, &Vector3, &mut Decider) -> Result<Vec<Origin>, Refusal>,
    {
        let empty = BTreeSet::new();
        let selected = replacement.map_or(&empty, |(selected, _)| selected);
        let mut patches = self.patches.clone();
        for (i, patch) in patches.iter_mut().enumerate() {
            if sign(
                &area(&self.vertices, &patch.loops[0], &patch.plane.normal),
                dec,
            )? == Sign::Negative
            {
                for cycle in &mut patch.loops {
                    cycle.reverse();
                }
            }
            if selected.contains(&i) {
                let Surface::Plane(plane) = replacement
                    .expect("a selection is supplied only with its draft law")
                    .1
                    .support(&Surface::Plane(patch.plane.clone()), true, dec)?
                    .surface
                else {
                    unreachable!("plane support remains planar")
                };
                patch.plane = plane;
            }
        }
        let mut incident = vec![BTreeSet::new(); self.vertices.len()];
        for (f, patch) in patches.iter().enumerate() {
            for i in patch.loops.iter().flatten() {
                incident[*i].insert(f);
            }
        }
        // Arrangement subdivisions may leave degree-two points on a physical
        // edge. They are not vertex fans and must not constrain its new line.
        let mut redundant = BTreeSet::new();
        for (i, faces) in incident.iter().enumerate() {
            if faces.len() != 2 || !faces.iter().any(|f| selected.contains(f)) {
                continue;
            }
            let mut straight = true;
            for face in faces {
                for cycle in &patches[*face].loops {
                    if let Some(j) = cycle.iter().position(|v| *v == i) {
                        let a = &self.vertices[cycle[(j + cycle.len() - 1) % cycle.len()]];
                        let b = &self.vertices[i];
                        let c = &self.vertices[cycle[(j + 1) % cycle.len()]];
                        if sign(&b.minus(a).cross(&c.minus(b)).norm_squared(), dec)? != Sign::Zero {
                            straight = false;
                            break;
                        }
                    }
                }
            }
            if straight {
                redundant.insert(i);
            }
        }
        if !redundant.is_empty() {
            for patch in &mut patches {
                for cycle in &mut patch.loops {
                    cycle.retain(|v| !redundant.contains(v));
                }
            }
            for i in redundant {
                incident[i].clear();
            }
        }
        let mut moved = self.vertices.clone();
        for (i, faces) in incident.iter().enumerate() {
            if !faces.iter().any(|f| selected.contains(f)) {
                continue;
            }
            let faces: Vec<_> = faces.iter().copied().collect();
            let mut replacement = None;
            'triple: for a in 0..faces.len() {
                for b in a + 1..faces.len() {
                    for c in b + 1..faces.len() {
                        if let Some(point) = triple(
                            &patches[faces[a]].plane,
                            &patches[faces[b]].plane,
                            &patches[faces[c]].plane,
                            dec,
                        )? {
                            replacement = Some((point, [a, b, c]));
                            break 'triple;
                        }
                    }
                }
            }
            let Some((point, defining)) = replacement else {
                return Err(failure(
                    "draft.unsupported-vertex-arrangement",
                    "incident moved carriers need an unbounded or resolved vertex fan",
                ));
            };
            // Its three defining carriers contain it by construction. Only a
            // further incident carrier can disagree, as in a non-simple vertex.
            for (index, f) in faces.into_iter().enumerate() {
                if !defining.contains(&index)
                    && sign(&value(&patches[f].plane, &point), dec)? != Sign::Zero
                {
                    return Err(failure(
                        "draft.unsupported-nonsimple-vertex",
                        "a non-simple moved vertex needs a resolved carrier fan",
                    ));
                }
            }
            moved[i] = point;
        }
        // A strict construction enclosure of the finite patch chain. No chain
        // crosses its boundary; every box interface therefore has winding zero.
        let mut bound = Rational::zero();
        for point in &moved {
            for x in [&point.x, &point.y, &point.z] {
                let b = dec.enclose_at(x, 128)?;
                bound = bound.max(b.lo().abs().max(b.hi().abs()));
            }
        }
        bound = bound + Rational::one();
        let mut vertices: Vec<_> = (0..8)
            .map(|i| {
                let [x, y, z] = [0, 1, 2].map(|axis| {
                    Real::from(if i & (1 << axis) == 0 {
                        -&bound
                    } else {
                        bound.clone()
                    })
                });
                Point3::new(x, y, z)
            })
            .collect();
        let mut planes = Vec::new();
        for axis in 0..3 {
            for direction in [-1, 1] {
                let mut coordinate = [Real::zero(), Real::zero(), Real::zero()];
                coordinate[axis] = Real::from(&bound * &Rational::from_i64(direction));
                let mut normal = [Real::zero(), Real::zero(), Real::zero()];
                normal[axis] = Real::from_i64(direction);
                let mut reference = [Real::zero(), Real::zero(), Real::zero()];
                reference[(axis + 1) % 3] = Real::one();
                planes.push(Plane {
                    origin: Point3::new(
                        coordinate[0].clone(),
                        coordinate[1].clone(),
                        coordinate[2].clone(),
                    ),
                    normal: Vector3::new(normal[0].clone(), normal[1].clone(), normal[2].clone()),
                    ref_dir: Vector3::new(
                        reference[0].clone(),
                        reference[1].clone(),
                        reference[2].clone(),
                    ),
                });
            }
        }
        let cycles = [
            [0, 4, 6, 2],
            [1, 3, 7, 5],
            [0, 1, 5, 4],
            [2, 6, 7, 3],
            [0, 2, 3, 1],
            [4, 5, 7, 6],
        ];
        let mut cells = vec![Cell {
            faces: cycles
                .into_iter()
                .enumerate()
                .map(|(support, c)| Facet {
                    support,
                    same_sense: true,
                    cycle: c.to_vec(),
                })
                .collect(),
        }];
        let mut carriers: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
        for (patch_id, patch) in patches.iter().enumerate() {
            let mut known = None;
            for (i, plane) in planes.iter().enumerate() {
                if sign(&plane.normal.cross(&patch.plane.normal).norm_squared(), dec)? == Sign::Zero
                    && sign(&value(plane, &patch.plane.origin), dec)? == Sign::Zero
                {
                    known = Some(i);
                    break;
                }
            }
            if let Some(i) = known {
                carriers.entry(i).or_default().push(patch_id);
                continue;
            }
            let support = planes.len();
            planes.push(patch.plane.clone());
            carriers.entry(support).or_default().push(patch_id);
            let signs = vertices
                .iter()
                .map(|p| sign(&value(&patch.plane, p), dec))
                .collect::<Result<Vec<_>, _>>()?;
            let mut clipper = Clipper {
                support,
                planes: &planes,
                vertices: &mut vertices,
                signs,
                edges: BTreeMap::new(),
                dec,
            };
            let mut next = Vec::new();
            for cell in cells {
                let mut positive = false;
                let mut negative = false;
                for i in cell.faces.iter().flat_map(|f| &f.cycle) {
                    positive |= clipper.signs[*i] == Sign::Positive;
                    negative |= clipper.signs[*i] == Sign::Negative;
                }
                if positive && negative {
                    next.push(clipper.clip(&cell, true)?);
                    next.push(clipper.clip(&cell, false)?);
                } else {
                    next.push(cell);
                }
            }
            cells = next;
        }
        let exterior = cells.len();
        let mut matches = FacetOwners::new();
        for (i, cell) in cells.iter().enumerate() {
            for (j, face) in cell.faces.iter().enumerate() {
                let mut cycle = face.cycle.clone();
                cycle.sort_unstable();
                matches
                    .entry((face.support, cycle))
                    .or_default()
                    .push((i, j));
            }
        }
        let mut interfaces = Vec::new();
        let mut adjacent = vec![Vec::new(); exterior + 1];
        for ((support, _), owners) in matches {
            if owners.len() > 2 || (owners.len() == 1 && support >= 6) {
                return Err(failure(
                    "draft.nonclosing-boundary",
                    "plane subdivision has unmatched facets",
                ));
            }
            let (from, face) = owners[0];
            let facet = cells[from].faces[face].clone();
            let to = owners.get(1).map_or(exterior, |x| x.0);
            let q = centre(&vertices, &facet.cycle);
            let mut positive = BTreeSet::new();
            let mut negative = BTreeSet::new();
            if let Some(raw) = carriers.get(&support).filter(|_| classified.is_none()) {
                for i in raw {
                    let patch = &patches[*i];
                    let winding = patch_winding(patch, &moved, &q, dec)?;
                    if winding == 0 {
                        continue;
                    }
                    let aligned = sign(&patch.plane.normal.dot(&planes[support].normal), dec)?
                        == Sign::Positive;
                    let outward = (winding > 0) == (aligned == facet.same_sense);
                    if outward {
                        positive.extend(patch.origins.iter().cloned());
                    } else {
                        negative.extend(patch.origins.iter().cloned());
                    }
                }
            }
            let (delta, origins) = if !positive.is_empty() && negative.is_empty() {
                (1, positive.into_iter().collect())
            } else if positive.is_empty() && !negative.is_empty() {
                (-1, negative.into_iter().collect())
            } else {
                (0, Vec::new())
            };
            let id = interfaces.len();
            interfaces.push(Interface {
                facet,
                from,
                to,
                delta,
                origins,
            });
            adjacent[from].push(id);
            adjacent[to].push(id);
        }
        let mut winding = vec![None; exterior + 1];
        winding[exterior] = Some(0_i32);
        let mut witnesses = if classified.is_some() {
            vec![None; exterior]
        } else {
            Vec::new()
        };
        if let Some((classify, resolve)) = &mut classified {
            for (i, cell) in cells.iter().enumerate() {
                let witness = cell_witness(&vertices, cell, dec)?;
                let inside = classify(&witness, dec)?;
                winding[i] = Some(i32::from(inside));
                if inside {
                    witnesses[i] = Some(witness);
                }
            }
            for interface in &mut interfaces {
                let from = winding[interface.from] == Some(1);
                let to = winding[interface.to] == Some(1);
                if from == to {
                    continue;
                }
                let facet = &interface.facet;
                let normal = planes[facet.support].normal.scaled(&Real::from_i64(
                    if facet.same_sense == from { 1 } else { -1 },
                ));
                interface.origins = resolve(&centre(&vertices, &facet.cycle), &normal, dec)?;
                interface.delta = i32::from(from) - i32::from(to);
            }
        } else {
            let mut queue = VecDeque::from([exterior]);
            while let Some(at) = queue.pop_front() {
                for i in &adjacent[at] {
                    let face = &interfaces[*i];
                    let (to, delta) = if at == face.from {
                        (face.to, -face.delta)
                    } else {
                        (face.from, face.delta)
                    };
                    let next = winding[at].expect("queued label") + delta;
                    if let Some(previous) = winding[to] {
                        if previous != next {
                            return Err(failure(
                                "draft.contradictory-orientation",
                                "candidate patch winding is path dependent",
                            ));
                        }
                    } else {
                        winding[to] = Some(next);
                        queue.push_back(to);
                    }
                }
            }
        }
        if winding
            .iter()
            .any(|w| w.is_none() || w.is_some_and(|v| v > 1))
        {
            return Err(failure(
                "draft.contradictory-orientation",
                "candidate cells do not have a single material winding",
            ));
        }
        let mut visited = vec![false; exterior];
        let mut result = Vec::new();
        for first in 0..exterior {
            if visited[first] || winding[first] != Some(1) {
                continue;
            }
            let mut component = BTreeSet::new();
            let mut queue = VecDeque::from([first]);
            visited[first] = true;
            while let Some(at) = queue.pop_front() {
                component.insert(at);
                for i in &adjacent[at] {
                    let f = &interfaces[*i];
                    let to = if at == f.from { f.to } else { f.from };
                    if to != exterior && !visited[to] && winding[to] == Some(1) {
                        visited[to] = true;
                        queue.push_back(to);
                    }
                }
            }
            let mut boundary = Vec::<(Facet, Vec<Origin>)>::new();
            for interface in &interfaces {
                let (a, b) = (
                    component.contains(&interface.from),
                    component.contains(&interface.to),
                );
                if a == b {
                    continue;
                }
                let mut facet = interface.facet.clone();
                if b {
                    facet.same_sense = !facet.same_sense;
                    facet.cycle.reverse();
                }
                if facet.support < 6 || interface.origins.is_empty() {
                    return Err(failure(
                        "draft.nonclosing-boundary",
                        "material requires an undeclared boundary patch",
                    ));
                }
                boundary.push((facet, interface.origins.clone()));
            }
            let mut parent: Vec<_> = (0..boundary.len()).collect();
            let mut edges = BTreeMap::new();
            for (i, (face, _)) in boundary.iter().enumerate() {
                for j in 0..face.cycle.len() {
                    let (a, b) = (face.cycle[j], face.cycle[(j + 1) % face.cycle.len()]);
                    let key = (face.support, face.same_sense, a.min(b), a.max(b));
                    if let Some(previous) = edges.insert(key, i) {
                        join(&mut parent, i, previous);
                    }
                }
            }
            let mut groups: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
            for i in 0..boundary.len() {
                groups.entry(root(&mut parent, i)).or_default().push(i);
            }
            let used: BTreeSet<_> = boundary
                .iter()
                .flat_map(|(f, _)| f.cycle.iter().copied())
                .collect();
            let remap: BTreeMap<_, _> = used
                .iter()
                .enumerate()
                .map(|(new, old)| (*old, new))
                .collect();
            let points: Vec<_> = used.iter().map(|i| vertices[*i].clone()).collect();
            let mut faces = Vec::new();
            let mut incidence: BTreeMap<(usize, usize), (usize, i32)> = BTreeMap::new();
            for group in groups.into_values() {
                let carrier = &boundary[group[0]].0;
                let mut origins = BTreeSet::new();
                let mut tiles = Vec::new();
                for i in group {
                    origins.extend(boundary[i].1.iter().cloned());
                    tiles.push(boundary[i].0.cycle.iter().map(|v| remap[v]).collect());
                }
                let mut loops = loops_from_tiles(&tiles)?;
                let normal = if carrier.same_sense {
                    planes[carrier.support].normal.clone()
                } else {
                    planes[carrier.support].normal.scaled(&Real::from_i64(-1))
                };
                let mut outer = None;
                for (i, cycle) in loops.iter().enumerate() {
                    if sign(&area(&points, cycle, &normal), dec)? == Sign::Positive
                        && outer.replace(i).is_some()
                    {
                        return Err(failure(
                            "draft.self-touching-shell",
                            "a merged patch has disconnected outer wires",
                        ));
                    }
                }
                let Some(outer) = outer else {
                    return Err(failure(
                        "draft.nonclosing-boundary",
                        "a physical face has no outer wire",
                    ));
                };
                loops.swap(0, outer);
                for cycle in &loops {
                    for i in 0..cycle.len() {
                        let (a, b) = (cycle[i], cycle[(i + 1) % cycle.len()]);
                        let entry = incidence.entry((a.min(b), a.max(b))).or_default();
                        entry.0 += 1;
                        entry.1 += if a < b { 1 } else { -1 };
                    }
                }
                faces.push(PolyFace {
                    support: carrier.support,
                    same_sense: carrier.same_sense,
                    loops,
                    tiles,
                    origins: origins.into_iter().collect(),
                });
            }
            if incidence
                .values()
                .any(|(count, balance)| *count != 2 || *balance != 0)
            {
                return Err(failure(
                    "draft.self-touching-shell",
                    "candidate material has a non-manifold boundary edge",
                ));
            }
            let terms = component
                .into_iter()
                .map(|i| -> Result<ConvexTerm, Refusal> {
                    let witness = if let Some(point) = witnesses.get_mut(i).and_then(Option::take) {
                        point
                    } else {
                        cell_witness(&vertices, &cells[i], dec)?
                    };
                    Ok(ConvexTerm {
                        sides: cells[i]
                            .faces
                            .iter()
                            .map(|f| (f.support, f.same_sense))
                            .collect(),
                        witness,
                    })
                })
                .collect::<Result<Vec<_>, _>>()?;
            result.push(Polyhedron {
                planes: planes.clone(),
                vertices: points,
                faces,
                terms,
            });
        }
        if result.is_empty() && replacement.is_some() {
            return Err(failure(
                "draft.reversed-solid",
                "replaced boundary has no full-dimensional material cell",
            ));
        }
        Ok(result)
    }
}

impl Polyhedron {
    /// Outward carrier for a physical face.
    pub fn plane(&self, face: &PolyFace) -> Plane {
        let mut plane = self.planes[face.support].clone();
        if !face.same_sense {
            plane.normal = plane.normal.scaled(&Real::from_i64(-1));
        }
        plane
    }
    /// Reconstructs source patches without carrying internal subdivision edges
    /// into a subsequent draft's adjacency seed.
    pub fn skin(&self) -> Skin {
        Skin {
            vertices: self.vertices.clone(),
            patches: self
                .faces
                .iter()
                .map(|f| Patch {
                    plane: self.plane(f),
                    loops: f.loops.clone(),
                    origins: f.origins.clone(),
                })
                .collect(),
        }
    }
    /// Closed membership at a point, or at its positive infinitesimal ray.
    pub fn contains(
        &self,
        p: &Point3,
        direction: Option<&Vector3>,
        dec: &mut Decider,
    ) -> Result<bool, Refusal> {
        for term in &self.terms {
            let mut inside = true;
            for (i, same) in &term.sides {
                let plane = &self.planes[*i];
                let mut s = sign(&value(plane, p), dec)?;
                if s == Sign::Zero
                    && let Some(d) = direction
                {
                    s = sign(&plane.normal.dot(d), dec)?;
                }
                if (*same && s == Sign::Positive) || (!*same && s == Sign::Negative) {
                    inside = false;
                    break;
                }
            }
            if inside {
                return Ok(true);
            }
        }
        Ok(false)
    }
    /// Closed trimmed-face membership, including a tangent infinitesimal ray.
    pub fn on_face(
        &self,
        face: &PolyFace,
        p: &Point3,
        direction: Option<&Vector3>,
        dec: &mut Decider,
    ) -> Result<bool, Refusal> {
        let plane = &self.planes[face.support];
        if sign(&value(plane, p), dec)? != Sign::Zero {
            return Ok(false);
        }
        if let Some(d) = direction
            && sign(&plane.normal.dot(d), dec)? != Sign::Zero
        {
            return Ok(false);
        }
        let normal = if face.same_sense {
            plane.normal.clone()
        } else {
            plane.normal.scaled(&Real::from_i64(-1))
        };
        for tile in &face.tiles {
            let mut inside = true;
            for i in 0..tile.len() {
                let a = &self.vertices[tile[i]];
                let edge = self.vertices[tile[(i + 1) % tile.len()]].minus(a);
                let mut s = sign(&edge.cross(&p.minus(a)).dot(&normal), dec)?;
                if s == Sign::Zero
                    && let Some(d) = direction
                {
                    s = sign(&edge.cross(d).dot(&normal), dec)?;
                }
                if s == Sign::Negative {
                    inside = false;
                    break;
                }
            }
            if inside {
                return Ok(true);
            }
        }
        Ok(false)
    }
    /// Reads a point on a face's carrier against its trim from enclosures alone.
    /// A point on a tile edge costs no exact zero, because the caller only needs
    /// that edge, which lies in the face, and not the answer.
    pub fn foot(
        &self,
        face: &PolyFace,
        p: &Point3,
        precision: u32,
        dec: &mut Decider,
    ) -> Result<Foot, Refusal> {
        let plane = &self.planes[face.support];
        let normal = if face.same_sense {
            plane.normal.clone()
        } else {
            plane.normal.scaled(&Real::from_i64(-1))
        };
        let mut near = Vec::new();
        'tiles: for tile in &face.tiles {
            let mut open = Vec::new();
            for i in 0..tile.len() {
                let (a, b) = (tile[i], tile[(i + 1) % tile.len()]);
                let edge = self.vertices[b].minus(&self.vertices[a]);
                let side = edge.cross(&p.minus(&self.vertices[a])).dot(&normal);
                match dec.enclose_budgeted_at(&side, precision)?.sign() {
                    Some(StrictSign::Negative) => continue 'tiles,
                    Some(StrictSign::Positive) => (),
                    None => open.push((a, b)),
                }
            }
            if open.is_empty() {
                return Ok(Foot::Inside);
            }
            near.append(&mut open);
        }
        Ok(if near.is_empty() {
            Foot::Outside
        } else {
            Foot::Near(near)
        })
    }
    /// Boundary membership excludes internal convex-decomposition seams.
    pub fn on_boundary(
        &self,
        p: &Point3,
        direction: Option<&Vector3>,
        dec: &mut Decider,
    ) -> Result<bool, Refusal> {
        for f in &self.faces {
            if self.on_face(f, p, direction, dec)? {
                return Ok(true);
            }
        }
        Ok(false)
    }
    /// Strict retained-region picks never accept a boundary point.
    pub fn contains_interior(&self, p: &Point3, dec: &mut Decider) -> Result<bool, Refusal> {
        Ok(self.contains(p, None, dec)? && !self.on_boundary(p, None, dec)?)
    }
    /// Positive-area neutral contact, not the convex hull of isolated contacts.
    pub fn anchored(&self, neutral: &Plane, dec: &mut Decider) -> Result<bool, Refusal> {
        // A face on the neutral plane already places its vertices there.
        let mut contained = BTreeSet::new();
        for face in &self.faces {
            let plane = &self.planes[face.support];
            if sign(&plane.normal.cross(&neutral.normal).norm_squared(), dec)? == Sign::Zero
                && sign(&value(neutral, &plane.origin), dec)? == Sign::Zero
            {
                contained.extend(face.tiles.iter().flatten().copied());
            }
        }
        let signs = self
            .vertices
            .iter()
            .enumerate()
            .map(|(i, p)| {
                if contained.contains(&i) {
                    Ok(Sign::Zero)
                } else {
                    sign(&value(neutral, p), dec)
                }
            })
            .collect::<Result<Vec<_>, _>>()?;
        if signs.contains(&Sign::Positive) && signs.contains(&Sign::Negative) {
            return Ok(true);
        }
        Ok(self
            .faces
            .iter()
            .flat_map(|f| &f.tiles)
            .any(|tile| tile.iter().all(|i| signs[*i] == Sign::Zero)))
    }
}

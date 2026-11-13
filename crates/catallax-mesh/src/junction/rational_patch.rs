//! Exact boundary-station propagation before any rounded mesh is published.

use std::{cmp::Ordering, collections::BTreeSet};

use catallax_numbers::{Outcome, StrictSign, ZeroProof};
use catallax_shapes::{BSplineCurve, RationalPatch, SplineBoundary};
use catallax_topology::RationalPatchFace;

use super::*;

pub(super) struct Proof {
    pub normal_boxes: Option<Vec<[Ball; 3]>>,
    pub embedding: Embedding,
}

pub(super) enum Embedding {
    Patch(Box<catallax_shapes::RationalGridEmbedding>),
    PositionOnly,
}

impl Embedding {
    pub(super) fn accepts_rounding(&self, error: &Rational) -> bool {
        match self {
            Self::Patch(proof) => proof.accepts_rounding(error),
            Self::PositionOnly => false,
        }
    }
}

pub(super) struct Grid {
    pub grid: EnvelopeGrid,
    pub parameters: [Vec<Real>; 2],
    pub normal_boxes: Option<Vec<[Ball; 3]>>,
    pub embedding: Embedding,
}

pub(super) enum Support {
    Patch(RationalPatch),
    Measured(catallax_shapes::BSplineSurface),
    Fillet(Box<super::spline_fillet::Mapping>),
}

pub(super) struct Definition {
    pub face: FaceId,
    pub support: Support,
    pub rows: [(EdgeId, BoundaryRow, Sense, bool); 4],
}

impl Definition {
    fn face(&self) -> FaceId {
        self.face
    }
    fn rows(&self) -> [(EdgeId, BoundaryRow, Sense, bool); 4] {
        self.rows
    }
    fn boundary(&self, row: BoundaryRow) -> BSplineCurve {
        match &self.support {
            Support::Patch(patch) => patch.boundary(axis_and_side(row).1),
            Support::Fillet(mapping) => mapping.boundary(row),
            Support::Measured(surface) => BSplineCurve::boundary(surface, axis_and_side(row).1),
        }
    }
    fn grid(
        &self,
        tolerance: &Rational,
        options: &MeshOptions,
        required: [&[Real]; 2],
        measured: bool,
        dec: &mut Decider,
    ) -> Result<Grid, Refusal> {
        if std::env::var_os("CX_F04_PROFILE").is_some() {
            eprintln!(
                "F04 face={:?} measured={measured} start={:?}",
                self.face,
                catallax_numbers::current_operation_work()
            );
        }
        if measured {
            return match &self.support {
                Support::Patch(patch) => {
                    measured_grid(patch.carrier(), tolerance, options, required, dec)
                }
                Support::Measured(surface) => {
                    measured_grid(surface, tolerance, options, required, dec)
                }
                Support::Fillet(mapping) => {
                    mapping.position_grid(tolerance, options.max_vertices, required, dec)
                }
            };
        }
        match &self.support {
            Support::Measured(_) => Err(bad("a positional surface has no certified embedding")),
            Support::Patch(patch) => {
                let grid = patch.display_grid(
                    tolerance,
                    options.angular_deviation.as_ref(),
                    options.max_edge_length.as_ref(),
                    options.max_vertices,
                    required,
                    dec,
                )?;
                Ok(Grid {
                    grid: grid.grid,
                    parameters: grid.parameters,
                    normal_boxes: grid.normal_boxes,
                    embedding: Embedding::Patch(Box::new(grid.embedding)),
                })
            }
            Support::Fillet(_) => Err(bad("a positional fillet has no certified embedding")),
        }
    }
}

fn measured_grid(
    surface: &catallax_shapes::BSplineSurface,
    tolerance: &Rational,
    options: &MeshOptions,
    required: [&[Real]; 2],
    dec: &mut Decider,
) -> Result<Grid, Refusal> {
    let evaluator = super::position::PositionSurface::new(surface, dec)?;
    let axes = std::array::from_fn(|a| {
        let domain = surface.domain()[a];
        let mut knots = surface.knots()[a]
            .iter()
            .filter(|k| *k >= domain[0] && *k <= domain[1])
            .cloned()
            .collect::<Vec<_>>();
        knots.sort();
        knots.dedup();
        knots.into_iter().map(Real::from).collect()
    });
    evaluator.grid(
        axes,
        required,
        super::position::Controls {
            tolerance,
            max_vertices: options.max_vertices,
            displacement: &Rational::zero(),
            cartesian: true,
        },
        |at| Ok(at.clone()),
        dec,
    )
}

struct Patch {
    checked: Definition,
    required: [Vec<Real>; 2],
    grid: Option<Grid>,
    dirty: bool,
}

struct Row {
    patch: usize,
    row: BoundaryRow,
    axis: usize,
    forward: bool,
    curve: BSplineCurve,
}

fn bad(detail: &str) -> Refusal {
    Refusal::new("mesh.rational-patch-unproven", detail)
}

fn compare(a: &Real, b: &Real, dec: &mut Decider) -> Result<Ordering, Refusal> {
    match dec.sign(&(a - b)) {
        Outcome::Proven(StrictSign::Negative) => Ok(Ordering::Less),
        Outcome::Proven(StrictSign::Positive) => Ok(Ordering::Greater),
        Outcome::ProvenCoincident(proof) if proof.proof != ZeroProof::SeparationBound => {
            Ok(Ordering::Equal)
        }
        Outcome::ProvenCoincident(_) => Err(bad(
            "shared stations need exact equality, not bounded coincidence",
        )),
        Outcome::NotExact(error) => Err(Refusal::from(error)),
        Outcome::Refused(error) => Err(error),
    }
}

pub(super) fn insert(
    values: &mut Vec<Real>,
    value: Real,
    dec: &mut Decider,
) -> Result<bool, Refusal> {
    let (mut lo, mut hi) = (0, values.len());
    while lo < hi {
        let mid = lo + (hi - lo) / 2;
        match compare(&value, &values[mid], dec)? {
            Ordering::Less => hi = mid,
            Ordering::Equal => return Ok(false),
            Ordering::Greater => lo = mid + 1,
        }
    }
    values.insert(lo, value);
    Ok(true)
}

fn additions(current: &[Real], previous: &[Real], dec: &mut Decider) -> Result<Vec<Real>, Refusal> {
    let mut old = 0;
    let mut added = Vec::new();
    for value in current {
        while old < previous.len() && compare(&previous[old], value, dec)? == Ordering::Less {
            old += 1;
        }
        if old == previous.len() || compare(&previous[old], value, dec)? != Ordering::Equal {
            added.push(value.clone());
        }
    }
    Ok(added)
}

fn axis_and_side(row: BoundaryRow) -> (usize, SplineBoundary) {
    match row {
        BoundaryRow::Rail => (0, SplineBoundary::VMin),
        BoundaryRow::Mouth => (0, SplineBoundary::VMax),
        BoundaryRow::Negative => (1, SplineBoundary::UMin),
        BoundaryRow::Positive => (1, SplineBoundary::UMax),
    }
}

struct Shared {
    rows: Vec<Row>,
    stations: usize,
}

fn reserve(count: &mut usize, extra: usize, maximum: usize) -> Result<(), Refusal> {
    *count = count
        .checked_add(extra)
        .filter(|n| *n <= maximum)
        .ok_or_else(|| {
            Refusal::new(
                "mesh.vertex-limit",
                "shared rational patch grids exceed max_vertices",
            )
        })?;
    Ok(())
}

pub(super) fn prepare(
    model: &Model,
    faces: &[FaceId],
    options: &MeshOptions,
    dec: &mut Decider,
) -> Result<Junction, Refusal> {
    let mut definitions = Vec::with_capacity(faces.len());
    for &face in faces {
        let checked = RationalPatchFace::of_face(model, face, dec)?
            .ok_or_else(|| bad("a generic patch mesh requires its checked finite domain"))?;
        definitions.push(Definition {
            face,
            support: Support::Patch(checked.support().clone()),
            rows: checked.rows(),
        });
    }
    prepare_definitions(model, definitions, options, false, dec)
}

pub(super) fn prepare_definitions(
    model: &Model,
    definitions: Vec<Definition>,
    options: &MeshOptions,
    measured: bool,
    dec: &mut Decider,
) -> Result<Junction, Refusal> {
    let mut patches = Vec::with_capacity(definitions.len());
    let mut edges: BTreeMap<EdgeId, Shared> = BTreeMap::new();
    let mut vertices = BTreeSet::new();
    for checked in definitions {
        let index = patches.len();
        for (edge, row, _, forward) in checked.rows() {
            let (axis, _) = axis_and_side(row);
            let group = edges.entry(edge).or_insert_with(|| Shared {
                rows: Vec::with_capacity(2),
                stations: 2,
            });
            if group.rows.len() == 2 {
                return Err(bad("a rational mesh edge has more than two face uses"));
            }
            group.rows.push(Row {
                patch: index,
                row,
                axis,
                forward,
                curve: checked.boundary(row),
            });
            vertices.extend(model.edge(edge).start);
            vertices.extend(model.edge(edge).end);
        }
        patches.push(Patch {
            checked,
            required: [Vec::new(), Vec::new()],
            grid: None,
            dirty: true,
        });
    }
    let mut count = 0;
    reserve(&mut count, vertices.len(), options.max_vertices)?;
    let tolerance = if measured {
        &options.chord_tolerance * &Rational::from_ratio(63, 64).unwrap()
    } else {
        options.chord_tolerance.mul_pow2(-1)
    };
    // Parameters need not agree on a shared geometric row: even a straight
    // row can have a nonlinear rational speed. Transfer exact physical points
    // through each neighbour's complete inverse, rather than merging indices.
    while patches.iter().any(|patch| patch.dirty) {
        for index in 0..patches.len() {
            if !patches[index].dirty {
                continue;
            }
            dec.check_deadline()?;
            let patch = &patches[index];
            let grid = patch.checked.grid(
                &tolerance,
                options,
                patch.required.each_ref().map(Vec::as_slice),
                measured,
                dec,
            )?;
            let previous = patch.grid.as_ref();
            let added = [
                additions(
                    &grid.parameters[0][1..grid.parameters[0].len() - 1],
                    previous.map_or(&[], |g| &g.parameters[0]),
                    dec,
                )?,
                additions(
                    &grid.parameters[1][1..grid.parameters[1].len() - 1],
                    previous.map_or(&[], |g| &g.parameters[1]),
                    dec,
                )?,
            ];
            let own_rows = patch.checked.rows();
            let previous_interior =
                previous.map_or(0, |g| (g.grid.columns - 2) * (g.grid.rows - 2));
            reserve(
                &mut count,
                (grid.grid.columns - 2) * (grid.grid.rows - 2) - previous_interior,
                options.max_vertices,
            )?;
            for &(edge, row, _, _) in &own_rows {
                let group = edges.get_mut(&edge).expect("checked row registered");
                let stations = grid.parameters[axis_and_side(row).0].len();
                if stations > group.stations {
                    reserve(&mut count, stations - group.stations, options.max_vertices)?;
                    group.stations = stations;
                }
            }
            patches[index].required = grid.parameters.clone();
            patches[index].grid = Some(grid);
            patches[index].dirty = false;
            for (edge, _, _, _) in own_rows {
                let rows = &edges[&edge].rows;
                let source = rows
                    .iter()
                    .find(|row| row.patch == index)
                    .expect("checked row registered");
                for target in rows.iter().filter(|row| row.patch != index) {
                    for parameter in &added[source.axis] {
                        let other = if matches!(
                            model.geom().curve(model.edge(edge).curve),
                            catallax_shapes::Curve::BSpline(_)
                        ) {
                            // Both checked rows have an exact normalized
                            // identity to this same finite B-spline edge.
                            // Reuse that proof instead of solving each point's
                            // neighbouring inverse again.
                            let [lo, hi] = source.curve.domain();
                            let mut fraction = &(parameter - &Real::from_rational(lo.clone()))
                                / &Real::from_rational(hi - lo);
                            if source.forward != target.forward {
                                fraction = &Real::one() - &fraction;
                            }
                            let [lo, hi] = target.curve.domain();
                            &Real::from_rational(lo.clone())
                                + &(&fraction * &Real::from_rational(hi - lo))
                        } else {
                            let point = source.curve.point_at(parameter, dec)?;
                            target.curve.unique_parameter_of(&point, dec)?
                                .ok_or_else(|| bad("a sewn boundary point has no unique exact neighbouring preimage"))?
                                .value().clone()
                        };
                        if insert(&mut patches[target.patch].required[target.axis], other, dec)? {
                            patches[target.patch].dirty = true;
                        }
                    }
                }
            }
        }
    }
    let mut seeds = BTreeMap::new();
    for (&edge, group) in &edges {
        let uses = &group.rows;
        let source = &uses[0];
        let grid = patches[source.patch]
            .grid
            .as_ref()
            .expect("all grids accepted");
        let count = grid.parameters[source.axis].len();
        for other in &uses[1..] {
            if patches[other.patch]
                .grid
                .as_ref()
                .expect("all grids accepted")
                .parameters[other.axis]
                .len()
                != count
            {
                return Err(bad(
                    "shared rational rows retain different certified station counts",
                ));
            }
        }
        let exact = match &patches[source.patch].checked.support {
            Support::Fillet(mapping) => {
                mapping.exact_row(source.row, &grid.parameters[source.axis], dec)?
            }
            Support::Patch(_) | Support::Measured(_) => grid.parameters[source.axis]
                .iter()
                .map(|parameter| {
                    source
                        .curve
                        .point_at(parameter, dec)
                        .map(|point| [point.x, point.y, point.z])
                })
                .collect::<Result<Vec<_>, Refusal>>()?,
        };
        let mut points = row_points(&grid.grid, source.row);
        let mut error = grid.grid.chord.clone();
        if matches!(patches[source.patch].checked.support, Support::Fillet(_))
            && matches!(source.row, BoundaryRow::Negative | BoundaryRow::Positive)
        {
            let projected = exact
                .iter()
                .map(|p| {
                    let [x, y, z] = p.each_ref().map(|v| dec.enclose_at(v, BITS));
                    Ok([x?, y?, z?])
                })
                .collect::<Result<Vec<_>, Refusal>>()?;
            let displacement = points
                .iter()
                .zip(&projected)
                .map(|(a, b)| {
                    (0..3).fold(Rational::zero(), |sum, i| {
                        sum + (&a[i] - &b[i]).lo().abs().max((&a[i] - &b[i]).hi().abs())
                    })
                })
                .max()
                .unwrap_or_else(Rational::zero);
            error = error + displacement;
            points = projected;
        }
        seeds.insert(
            edge,
            Seed {
                points,
                error,
                forward: source.forward,
                exact,
            },
        );
    }
    let mut envelopes = BTreeMap::new();
    for patch in patches {
        let rows = patch.checked.rows();
        let along_u = (rows[0].2 == Sense::Same) == rows[0].3;
        let grid = patch.grid.expect("all grids accepted");
        envelopes.insert(
            patch.checked.face(),
            Envelope {
                grid: grid.grid,
                rows: rows
                    .into_iter()
                    .map(|(edge, row, _, forward)| (edge, row, forward))
                    .collect(),
                along_u,
                ends: Ends::Meridians,
                rational: Some(Proof {
                    normal_boxes: grid.normal_boxes,
                    embedding: grid.embedding,
                }),
            },
        );
    }
    Ok(Junction { envelopes, seeds })
}

/// atan(|t × n| / (t · n)) is no larger than the ratio when the dot is
/// positive. All n in the complete cell box are covered, including rounding of t.
pub(super) fn angular(
    normal: &[Ball; 3],
    field: &[Ball; 3],
    forward: bool,
) -> Result<Rational, Refusal> {
    let normal = normal
        .each_ref()
        .map(|x| if forward { x.clone() } else { -x });
    let dot = (0..3).fold(Ball::from_rational(&Rational::zero(), BITS), |sum, i| {
        &sum + &(&normal[i] * &field[i])
    });
    if dot.lo() <= Rational::zero() {
        return Err(bad(
            "rounded facet normal orientation is not proved over its entire cell",
        ));
    }
    let cross = [
        &normal[1] * &field[2] - &normal[2] * &field[1],
        &normal[2] * &field[0] - &normal[0] * &field[2],
        &normal[0] * &field[1] - &normal[1] * &field[0],
    ];
    let numerator = cross.iter().fold(Rational::zero(), |sum, component| {
        sum + component.lo().abs().max(component.hi().abs())
    });
    Ok(&numerator / &dot.lo())
}

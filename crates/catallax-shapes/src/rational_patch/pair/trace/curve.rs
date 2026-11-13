//! Point queries on a component retained by the complete cell-pair search.

use super::*;

mod germ;

type Location = (usize, Real, [Real; 4]);

/// One complete regular component with both source covers and its exact graph
/// pieces. Only a successful complete patch-pair meeting can construct it.
#[derive(Clone, Debug)]
pub struct PatchPairCurve {
    pieces: Arc<[PairTracePiece]>,
}

impl PatchPairIntersection {
    /// Retain a component for exact point queries and subsequent tracing.
    pub fn curve_branch(
        &self,
        branch: usize,
        dec: &mut Decider,
    ) -> Result<PatchPairCurve, Refusal> {
        Ok(PatchPairCurve {
            pieces: self.trace_branch(branch, dec)?.into(),
        })
    }
}

impl PatchPairCurve {
    /// Exact graph pieces in first-normal cross second-normal traversal order.
    pub fn pieces(&self) -> &[PairTracePiece] {
        &self.pieces
    }

    /// Whether the complete endpoint graph is a cycle.
    pub fn is_closed(&self) -> bool {
        self.pieces[0].definition.branch.closed
    }

    /// Endpoints of an open component, or its seam point twice when closed.
    pub fn endpoints(&self) -> [Point3; 2] {
        [
            endpoint(&self.pieces[0], 0),
            endpoint(self.pieces.last().expect("a certified branch has arcs"), 1),
        ]
    }

    /// Exact join witnesses in retained traversal order.
    pub fn witnesses(&self) -> impl Iterator<Item = Point3> + '_ {
        self.pieces.iter().map(|piece| endpoint(piece, 0))
    }

    /// Whether the point belongs to this finite connected component.
    pub fn contains(&self, point: &Point3, dec: &mut Decider) -> Result<bool, Refusal> {
        Ok(self.locate(point, dec)?.is_some())
    }

    /// Open-component traversal order. Positive means `a` precedes `b`.
    pub fn compare(&self, a: &Point3, b: &Point3, dec: &mut Decider) -> Result<Sign, Refusal> {
        if self.is_closed() {
            return Err(failure("linear order requires an open pair component"));
        }
        let (i, a, _) = self
            .locate(a, dec)?
            .ok_or_else(|| failure("order point is outside the pair component"))?;
        let (j, b, _) = self
            .locate(b, dec)?
            .ok_or_else(|| failure("order point is outside the pair component"))?;
        match i.cmp(&j) {
            std::cmp::Ordering::Less => Ok(Sign::Positive),
            std::cmp::Ordering::Greater => Ok(Sign::Negative),
            std::cmp::Ordering::Equal => exact_sign(&(b - a), dec),
        }
    }

    /// Positive means `a,b,c` follow the closed component's traversal.
    pub fn cyclic_sign(&self, points: &[Point3; 3], dec: &mut Decider) -> Result<Sign, Refusal> {
        if !self.is_closed() {
            return Err(failure("cyclic order requires a closed pair component"));
        }
        let mut at = std::array::from_fn::<_, 3, _>(|_| Real::zero());
        for (point, position) in points.iter().zip(&mut at) {
            let (index, fraction, _) = self
                .locate(point, dec)?
                .ok_or_else(|| failure("cyclic-order point is outside the pair component"))?;
            *position = Real::from_i64(index as i64) + fraction;
        }
        exact_sign(
            &(&(&at[1] - &at[0]) * &(&(&at[2] - &at[0]) * &(&at[2] - &at[1]))),
            dec,
        )
    }

    /// Arriving and leaving tangents, retaining both sides of a chart join.
    pub fn tangents(&self, point: &Point3, dec: &mut Decider) -> Result<[Vector3; 2], Refusal> {
        let mut arriving = None;
        let mut leaving = None;
        for piece in self.pieces.iter() {
            for end in 0..2 {
                if same_point(&endpoint(piece, end), point, dec)? {
                    let tangent = tangent(piece, piece.arc().ends[end].parameters(), dec)?;
                    if end == 0 {
                        leaving = Some(tangent);
                    } else {
                        arriving = Some(tangent);
                    }
                }
            }
        }
        match (arriving, leaving) {
            (Some(a), Some(b)) => Ok([a, b]),
            (Some(t), None) | (None, Some(t)) => Ok([t.clone(), t]),
            (None, None) => {
                let (index, _, at) = self
                    .locate(point, dec)?
                    .ok_or_else(|| failure("tangent point is outside the pair component"))?;
                let value = tangent(&self.pieces[index], &at, dec)?;
                Ok([value.clone(), value])
            }
        }
    }

    fn locate(&self, point: &Point3, dec: &mut Decider) -> Result<Option<Location>, Refusal> {
        for (index, piece) in self.pieces.iter().enumerate() {
            for end in 0..2 {
                if same_point(&endpoint(piece, end), point, dec)? {
                    return Ok(Some((
                        index,
                        Real::from_i64(end as i64),
                        piece.arc().ends[end].parameters().clone(),
                    )));
                }
            }
        }
        let mut preimages = BTreeMap::new();
        for (index, piece) in self.pieces.iter().enumerate() {
            let arc = piece.arc();
            let Some(first) = preimage(
                &mut preimages,
                (0, arc.first_chart),
                &piece.definition.first,
                point,
                dec,
            )?
            else {
                continue;
            };
            let Some(second) = preimage(
                &mut preimages,
                (1, arc.second_chart),
                &piece.definition.second,
                point,
                dec,
            )?
            else {
                continue;
            };
            let [u, v] = first;
            let [s, t] = second;
            let at = [u, v, s, t];
            if !contains(&arc.first_cell, &[at[0].clone(), at[1].clone()], dec)?
                || !contains(&arc.second_cell, &[at[2].clone(), at[3].clone()], dec)?
            {
                continue;
            }
            let values = graph_values(&at);
            let constant = piece.parameter.eval_real(&values);
            let slope = piece.parameter.partial(4).eval_real(&values);
            if exact_sign(&slope, dec)? == Sign::Zero {
                return Err(failure(
                    "the retained graph parameter lost its nonzero slope",
                ));
            }
            let fraction = -constant / slope;
            if exact_sign(&fraction, dec)? == Sign::Negative
                || exact_sign(&(&Real::one() - &fraction), dec)? == Sign::Negative
            {
                return Err(failure(
                    "a meeting in the certified cell is outside its graph arc",
                ));
            }
            return Ok(Some((index, fraction, at)));
        }
        Ok(None)
    }
}

fn endpoint(piece: &PairTracePiece, end: usize) -> Point3 {
    let at = piece.arc().ends[end].parameters();
    piece.first().point(&[at[0].clone(), at[1].clone()])
}

fn same_point(a: &Point3, b: &Point3, dec: &mut Decider) -> Result<bool, Refusal> {
    Ok(exact_sign(&a.minus(b).norm_squared(), dec)? == Sign::Zero)
}

fn tangent(piece: &PairTracePiece, at: &[Real; 4], dec: &mut Decider) -> Result<Vector3, Refusal> {
    let arc = piece.arc();
    let [a, b] = piece.first().derivatives(&[at[0].clone(), at[1].clone()]);
    let [c, d] = piece
        .definition
        .second
        .span(arc.second_chart, arc.second_span)
        .derivatives(&[at[2].clone(), at[3].clone()]);
    let first_normal = a.cross(&b);
    let second_normal = c.cross(&d);
    let velocity = [
        -b.dot(&second_normal),
        a.dot(&second_normal),
        first_normal.dot(&d),
        -first_normal.dot(&c),
    ];
    let values = graph_values(at);
    let mut directional = Real::zero();
    for (column, speed) in velocity.iter().enumerate() {
        directional = directional + piece.parameter_jacobian[column].eval_real(&values) * speed;
    }
    if exact_sign(&directional, dec)? == Sign::Zero {
        return Err(failure(
            "the regular pair tangent lost its strict graph coordinate",
        ));
    }
    let slope = piece.parameter.partial(4).eval_real(&values);
    Ok(first_normal
        .cross(&second_normal)
        .scaled(&(-slope / directional)))
}

type Preimages = BTreeMap<(usize, usize), Option<[Real; 2]>>;

fn preimage(
    cache: &mut Preimages,
    key: (usize, usize),
    source: &RetainedSource,
    point: &Point3,
    dec: &mut Decider,
) -> Result<Option<[Real; 2]>, Refusal> {
    let at = match cache.entry(key) {
        std::collections::btree_map::Entry::Vacant(entry) => {
            entry.insert(source.patch(key.1).parameters_at(point, dec)?)
        }
        std::collections::btree_map::Entry::Occupied(entry) => entry.into_mut(),
    };
    Ok(at.clone())
}

fn graph_values(at: &[Real; 4]) -> [Real; 5] {
    [
        at[0].clone(),
        at[1].clone(),
        at[2].clone(),
        at[3].clone(),
        Real::zero(),
    ]
}

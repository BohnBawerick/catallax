//! Finite natural arcs retain their exact parameter endpoints and graph pin.

use super::*;

impl PatchPairCurve {
    /// Trace a finite natural arc. Equal endpoints select the whole closed
    /// component; an open component requires two ordered endpoints.
    pub fn trace_arc(
        &self,
        from: Option<&Point3>,
        to: Option<&Point3>,
        dec: &mut Decider,
    ) -> Result<Vec<PairTracePiece>, Refusal> {
        let endpoint = |point: &Point3, dec: &mut Decider| -> Result<(Real, [Real; 4]), Refusal> {
            let (index, fraction, parameters) = self
                .locate(point, dec)?
                .ok_or_else(|| failure("trace endpoint is outside its retained pair component"))?;
            Ok((Real::from_i64(index as i64) + fraction, parameters))
        };
        let count = Real::from_i64(self.pieces.len() as i64);
        let (start, mut end, start_at, end_at) = match (from, to) {
            (None, None) if self.is_closed() => (
                Real::zero(),
                count.clone(),
                self.pieces[0].arc().ends[0].parameters().clone(),
                self.pieces
                    .last()
                    .expect("a certified branch has arcs")
                    .arc()
                    .ends[1]
                    .parameters()
                    .clone(),
            ),
            (Some(from), Some(to)) => {
                let (start, start_at) = endpoint(from, dec)?;
                let (end, end_at) = endpoint(to, dec)?;
                (start, end, start_at, end_at)
            }
            _ => return Err(failure("an open pair trace needs both finite endpoints")),
        };
        if exact_sign(&(&end - &start), dec)? != Sign::Positive {
            if !self.is_closed() {
                return Err(failure(
                    "open pair trace endpoints oppose its natural direction",
                ));
            }
            end = &end + &count;
        }
        let mut pieces = Vec::new();
        for turn in 0..=usize::from(self.is_closed()) {
            for (index, piece) in self.pieces.iter().enumerate() {
                dec.check_deadline()?;
                let low = Real::from_i64((turn * self.pieces.len() + index) as i64);
                let high = &low + &Real::one();
                let start_inside = exact_sign(&(&start - &low), dec)? == Sign::Positive;
                let end_inside = exact_sign(&(&end - &high), dec)? == Sign::Negative;
                let a = if start_inside {
                    start.clone()
                } else {
                    low.clone()
                };
                let b = if end_inside {
                    end.clone()
                } else {
                    high.clone()
                };
                if exact_sign(&(&b - &a), dec)? != Sign::Positive {
                    continue;
                }
                if !start_inside && !end_inside {
                    pieces.push(piece.clone());
                    continue;
                }
                let first = if start_inside {
                    start_at.clone()
                } else {
                    piece.arc().ends[0].parameters().clone()
                };
                let last = if end_inside {
                    end_at.clone()
                } else {
                    piece.arc().ends[1].parameters().clone()
                };
                let from = &a - &low;
                let to = &b - &low;
                let affine = &RealPoly::constant(5, from.clone())
                    + &RealPoly::var(5, 4).scale(&(&to - &from));
                let parameter = piece.parameter.substitute(4, &affine);
                let parameter_jacobian = std::array::from_fn(|column| parameter.partial(column));
                pieces.push(PairTracePiece {
                    parameter,
                    parameter_jacobian,
                    clipped_ends: Some(Arc::new([first, last])),
                    ..piece.clone()
                });
            }
        }
        if pieces.is_empty() {
            return Err(failure("pair trace clipping found no exact owned arc"));
        }
        Ok(pieces)
    }
}

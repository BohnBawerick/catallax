//! Complete control-hull extrema with actual subdivision corners as witnesses.

use super::*;

type Box3 = [[Rational; 2]; 3];

struct Net {
    count: [usize; 2],
    controls: Vec<[Ball; 4]>,
    bits: u32,
    depth: usize,
}

impl Net {
    fn halves(self, dec: &Decider) -> Result<[Self; 2], Refusal> {
        let axis = self.depth % 2;
        let count = self.count;
        let degree = count[axis] - 1;
        let half = Ball::from_rational(&Rational::from_ratio(1, 2).expect("nonzero"), self.bits);
        let mut left = self.controls;
        let mut right = left.clone();
        let mut work = Vec::with_capacity(count[axis]);
        for fixed in 0..count[1 - axis] {
            let at = |i| {
                if axis == 0 {
                    i * count[1] + fixed
                } else {
                    fixed * count[1] + i
                }
            };
            work.clear();
            work.extend((0..count[axis]).map(|i| left[at(i)].clone()));
            for level in 1..=degree {
                dec.check_deadline()?;
                for i in 0..=degree - level {
                    let (before, after) = work.split_at_mut(i + 1);
                    for (value, next) in before[i].iter_mut().zip(&after[0]) {
                        *value = &(&*value + next) * &half;
                    }
                }
                left[at(level)] = work[0].clone();
                right[at(degree - level)] = work[degree - level].clone();
            }
        }
        Ok([left, right].map(|controls| Self {
            count,
            controls,
            bits: self.bits,
            depth: self.depth + 1,
        }))
    }
}

struct Cell {
    net: Net,
    bounds: Box3,
    low_witness: [Rational; 3],
    high_witness: [Rational; 3],
    rounding: Rational,
}

fn cartesian(h: &[Ball; 4]) -> Option<[Ball; 3]> {
    if h[3].lo() <= Rational::zero() {
        return None;
    }
    Some([
        h[0].checked_div(&h[3])?,
        h[1].checked_div(&h[3])?,
        h[2].checked_div(&h[3])?,
    ])
}

impl Cell {
    // None means insufficient positive-weight enclosure precision, not a pole.
    fn new(net: Net, dec: &Decider) -> Result<Option<Self>, Refusal> {
        let Some(first) = cartesian(&net.controls[0]) else {
            return Ok(None);
        };
        let mut bounds = first.each_ref().map(|x| [x.lo(), x.hi()]);
        let mut low_witness = first.each_ref().map(Ball::hi);
        let mut high_witness = first.each_ref().map(Ball::lo);
        let mut rounding = first
            .iter()
            .map(|x| x.hi() - x.lo())
            .max()
            .expect("three coordinates");
        for (index, h) in net.controls.iter().enumerate().skip(1) {
            dec.check_deadline()?;
            let Some(point) = cartesian(h) else {
                return Ok(None);
            };
            let corner = (index / net.count[1] == 0 || index / net.count[1] + 1 == net.count[0])
                && (index % net.count[1] == 0 || index % net.count[1] + 1 == net.count[1]);
            for axis in 0..3 {
                let lo = point[axis].lo();
                let hi = point[axis].hi();
                rounding = rounding.max(&hi - &lo);
                if lo < bounds[axis][0] {
                    bounds[axis][0] = lo.clone();
                }
                if hi > bounds[axis][1] {
                    bounds[axis][1] = hi.clone();
                }
                if corner {
                    if hi < low_witness[axis] {
                        low_witness[axis] = hi;
                    }
                    if lo > high_witness[axis] {
                        high_witness[axis] = lo;
                    }
                }
            }
        }
        Ok(Some(Self {
            net,
            bounds,
            low_witness,
            high_witness,
            rounding,
        }))
    }

    fn complete(&self, low: &[Rational; 3], high: &[Rational; 3], width: &Rational) -> bool {
        (0..3).all(|i| {
            &low[i] - &self.bounds[i][0] <= *width && &self.bounds[i][1] - &high[i] <= *width
        })
    }
}

fn witnesses(cell: &Cell, low: &mut [Rational; 3], high: &mut [Rational; 3]) {
    for i in 0..3 {
        if cell.low_witness[i] < low[i] {
            low[i] = cell.low_witness[i].clone();
        }
        if cell.high_witness[i] > high[i] {
            high[i] = cell.high_witness[i].clone();
        }
    }
}

impl RationalPatch {
    /// Smallest and largest world coordinates, each as a proved [lo,hi] pair
    /// no wider than `width`. Positive homogeneous weights bound complete
    /// subpatches by their control hulls; only actual corners set witnesses.
    pub fn ranges(
        &self,
        width: &Rational,
        dec: &mut Decider,
    ) -> Result<[[[Rational; 2]; 2]; 3], Refusal> {
        if width <= &Rational::zero() {
            return Err(Refusal::new(
                "measure.invalid-accuracy",
                "rational patch range width must be positive",
            ));
        }
        let maximum = dec.policy().max_precision_bits;
        if maximum == 0 {
            return Err(Refusal::new(
                "numbers.precision-limit",
                "no range enclosure precision",
            ));
        }
        let wanted = width
            .denominator()
            .bit_length()
            .saturating_sub(width.numerator().bit_length())
            .saturating_add(80)
            .max(128);
        let mut bits = wanted.min(u64::from(maximum)) as u32;
        let count = self.carrier().degrees().map(|degree| degree + 1);
        let rounding_goal = width.mul_pow2(-12);
        let mut open = loop {
            let mut cells = Vec::new();
            let mut enough = true;
            for (_, controls) in self.bezier_spans(dec)? {
                let controls = controls
                    .iter()
                    .map(|h| {
                        let [x, y, z, w] = h.each_ref().map(|x| dec.enclose_budgeted_at(x, bits));
                        Ok([x?, y?, z?, w?])
                    })
                    .collect::<Result<_, Refusal>>()?;
                let Some(cell) = Cell::new(
                    Net {
                        count,
                        controls,
                        bits,
                        depth: 0,
                    },
                    dec,
                )?
                else {
                    enough = false;
                    break;
                };
                if cell.rounding > rounding_goal {
                    enough = false;
                    break;
                }
                cells.push(cell);
            }
            if enough {
                break cells;
            }
            if bits == maximum {
                return Err(Refusal::new(
                    "numbers.precision-limit",
                    "range control rounding exceeds the requested width",
                ));
            }
            bits = bits.saturating_mul(2).min(maximum);
        };
        let first = open
            .first()
            .ok_or_else(|| unproven("a rational patch has no finite spans"))?;
        let mut low = first.low_witness.clone();
        let mut high = first.high_witness.clone();
        for cell in &open {
            witnesses(cell, &mut low, &mut high);
        }
        let mut accepted = None::<Box3>;
        let mut created = open.len();
        while !open.is_empty() {
            let mut next = Vec::new();
            for cell in open {
                dec.check_deadline()?;
                if cell.complete(&low, &high, width) {
                    if let Some(bounds) = &mut accepted {
                        for (bound, cell_bound) in bounds.iter_mut().zip(&cell.bounds) {
                            if cell_bound[0] < bound[0] {
                                bound[0] = cell_bound[0].clone();
                            }
                            if cell_bound[1] > bound[1] {
                                bound[1] = cell_bound[1].clone();
                            }
                        }
                    } else {
                        accepted = Some(cell.bounds);
                    }
                    continue;
                }
                created += 2;
                if created > MAX_PROOF_CELLS {
                    return Err(unproven(
                        "rational patch extrema exceeded their finite hull-subdivision allowance",
                    ));
                }
                for net in cell.net.halves(dec)? {
                    let child = Cell::new(net, dec)?.ok_or_else(|| {
                        Refusal::new(
                            "numbers.precision-limit",
                            "subdivided positive weight has no separated enclosure",
                        )
                    })?;
                    witnesses(&child, &mut low, &mut high);
                    next.push(child);
                }
            }
            open = next;
        }
        let bounds = accepted.ok_or_else(|| unproven("no complete extrema cover was retained"))?;
        Ok(std::array::from_fn(|i| {
            [
                [bounds[i][0].clone(), low[i].clone()],
                [high[i].clone(), bounds[i][1].clone()],
            ]
        }))
    }
}

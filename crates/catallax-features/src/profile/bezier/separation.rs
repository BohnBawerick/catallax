use super::*;

fn bounds(points: &[P], axis: &P) -> (Rational, Rational) {
    let values: Vec<_> = points.iter().map(|p| dot(p, axis)).collect();
    (
        values.iter().min().unwrap().clone(),
        values.iter().max().unwrap().clone(),
    )
}

fn on_line(point: &P, line: &[P]) -> bool {
    let chord = sub(line.last().unwrap(), &line[0]);
    let delta = sub(point, &line[0]);
    cross(&chord, &delta).is_zero()
        && dot(&delta, &chord).sign() != Sign::Negative
        && dot(&sub(point, line.last().unwrap()), &chord).sign() != Sign::Positive
}

// A nonconstant Bernstein row attains its hull's supporting value only at
// endpoints. A flat row is a line and needs a finite segment check instead.
fn contact(a: &[P], b: &[P], axis: &P, cut: &Rational, allowed: &[P]) -> bool {
    let flat_a = a.iter().all(|p| dot(p, axis) == *cut);
    let flat_b = b.iter().all(|p| dot(p, axis) == *cut);
    if flat_a && flat_b {
        return false;
    }
    let endpoints = |p: &[P]| [p[0].clone(), p.last().unwrap().clone()];
    if flat_a || flat_b {
        let (flat, curved) = if flat_a { (a, b) } else { (b, a) };
        return endpoints(curved)
            .iter()
            .filter(|p| dot(p, axis) == *cut && on_line(p, flat))
            .all(|p| allowed.contains(p));
    }
    endpoints(a)
        .iter()
        .filter(|p| dot(p, axis) == *cut)
        .all(|p| !endpoints(b).contains(p) || allowed.contains(p))
}

fn separated(a: &[P], b: &[P], allowed: &[P]) -> bool {
    let mut axes = vec![
        [Rational::one(), Rational::zero()],
        [Rational::zero(), Rational::one()],
    ];
    for points in [a, b] {
        for i in 0..points.len() {
            for j in i + 1..points.len() {
                let d = sub(&points[j], &points[i]);
                if !d.iter().all(Rational::is_zero) {
                    axes.push([-d[1].clone(), d[0].clone()]);
                }
            }
        }
    }
    for axis in axes {
        let (a_lo, a_hi) = bounds(a, &axis);
        let (b_lo, b_hi) = bounds(b, &axis);
        if a_hi < b_lo || b_hi < a_lo {
            return true;
        }
        if a_hi == b_lo && contact(a, b, &axis, &a_hi, allowed)
            || b_hi == a_lo && contact(a, b, &axis, &a_lo, allowed)
        {
            return true;
        }
    }
    false
}

fn split(points: &[P]) -> (Vec<P>, Vec<P>) {
    let mut row = points.to_vec();
    let mut left = vec![row[0].clone()];
    let mut right = vec![row.last().unwrap().clone()];
    for count in (1..row.len()).rev() {
        for i in 0..count {
            row[i] = std::array::from_fn(|j| (&row[i][j] + &row[i + 1][j]) / Rational::from_i64(2));
        }
        left.push(row[0].clone());
        right.push(row[count - 1].clone());
    }
    right.reverse();
    (left, right)
}

fn extent(points: &[P]) -> Rational {
    let (x_lo, x_hi) = bounds(points, &[Rational::one(), Rational::zero()]);
    let (y_lo, y_hi) = bounds(points, &[Rational::zero(), Rational::one()]);
    x_hi - x_lo + y_hi - y_lo
}

pub(super) fn prove(a: &[P], b: &[P], allowed: &[P], dec: &mut Decider) -> Result<(), Refusal> {
    let mut pending = vec![(a.to_vec(), b.to_vec(), 0)];
    let mut cells = 0;
    while let Some((a, b, depth)) = pending.pop() {
        dec.check_deadline()?;
        if separated(&a, &b, allowed) {
            continue;
        }
        if a.len() == 2 && b.len() == 2 {
            return Err(bad(
                "sketch.bezier-separation",
                "profile segments cross or overlap",
            ));
        }
        cells += 1;
        if cells > 4096 || depth >= 48 {
            return Err(bad(
                "sketch.bezier-proof",
                "control hull subdivision did not prove disjoint profile segments",
            ));
        }
        if b.len() == 2 || a.len() != 2 && extent(&a) >= extent(&b) {
            let (left, right) = split(&a);
            pending.push((right, b.clone(), depth + 1));
            pending.push((left, b, depth + 1));
        } else {
            let (left, right) = split(&b);
            pending.push((a.clone(), right, depth + 1));
            pending.push((a, left, depth + 1));
        }
    }
    Ok(())
}

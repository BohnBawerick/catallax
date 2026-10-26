//! Exact ray parity checks the producer's outer-shell and cavity ownership.
use catallax_numbers::{Decider, Rational, Refusal};
type Point = [Rational; 3];
fn sub(a: &Point, b: &Point) -> Point {
    std::array::from_fn(|i| &a[i] - &b[i])
}
fn cross(a: &Point, b: &Point) -> Point {
    std::array::from_fn(|i| {
        &(&a[(i + 1) % 3] * &b[(i + 2) % 3]) - &(&a[(i + 2) % 3] * &b[(i + 1) % 3])
    })
}
fn dot(a: &Point, b: &Point) -> Rational {
    (0..3).fold(Rational::zero(), |sum, i| &sum + &(&a[i] * &b[i]))
}
fn inside(
    point: &Point,
    shell: &[[usize; 3]],
    positions: &[Point],
    dec: &Decider,
) -> Result<bool, Refusal> {
    // Every exceptional ray direction is a root of a nonzero degree-at-most-two
    // polynomial in k. Retry exact edge/vertex hits, never count them twice.
    for k in 1..i64::MAX {
        dec.check_deadline()?;
        let k = Rational::from_i64(k);
        let direction = [Rational::one(), k.clone(), &k * &k];
        let mut hits = 0;
        let mut exceptional = false;
        for t in shell {
            dec.check_deadline()?;
            let [a, b, c] = t.map(|i| &positions[i]);
            let ab = sub(b, a);
            let ac = sub(c, a);
            let h = cross(&direction, &ac);
            let mut determinant = dot(&ab, &h);
            if determinant.is_zero() {
                exceptional = true;
                break;
            }
            let s = sub(point, a);
            let mut u = dot(&s, &h);
            let q = cross(&s, &ab);
            let mut v = dot(&direction, &q);
            let mut distance = dot(&ac, &q);
            if determinant < Rational::zero() {
                determinant = -determinant;
                u = -u;
                v = -v;
                distance = -distance;
            }
            if distance < Rational::zero()
                || u < Rational::zero()
                || v < Rational::zero()
                || &u + &v > determinant
            {
                continue;
            }
            if distance.is_zero() {
                return Err(Refusal::new(
                    "print.brep-orientation",
                    "encoded shells meet during nesting classification",
                ));
            }
            if u.is_zero() || v.is_zero() || &u + &v == determinant {
                exceptional = true;
                break;
            }
            hits += 1;
        }
        if !exceptional {
            return Ok(hits % 2 == 1);
        }
    }
    Err(Refusal::new(
        "print.brep-orientation",
        "no nonexceptional exact nesting ray within the work allowance",
    ))
}
pub(super) fn check(
    positions: &[Point],
    shells: &[Vec<[usize; 3]>],
    owners: &[(usize, bool)],
    dec: &Decider,
) -> Result<(), Refusal> {
    for (i, &(owner, outer)) in owners.iter().enumerate() {
        dec.check_deadline()?;
        let first = shells[i]
            .first()
            .ok_or_else(|| Refusal::new("print.brep-orientation", "empty represented shell"))?;
        let point = &positions[first[0]];
        if !outer {
            let parent = owners
                .iter()
                .position(|&(lump, is_outer)| lump == owner && is_outer)
                .expect("model lump has an outer shell");
            if !inside(point, &shells[parent], positions, dec)? {
                return Err(Refusal::new(
                    "print.brep-orientation",
                    "encoded cavity is outside its material outer shell",
                ));
            }
            for (j, &(lump, is_outer)) in owners.iter().enumerate() {
                if j != i
                    && lump == owner
                    && !is_outer
                    && inside(point, &shells[j], positions, dec)?
                {
                    return Err(Refusal::new(
                        "print.brep-orientation",
                        "encoded cavities of one lump overlap by nesting",
                    ));
                }
            }
        } else {
            for (j, &(other, is_outer)) in owners.iter().enumerate() {
                if other == owner || !is_outer || !inside(point, &shells[j], positions, dec)? {
                    continue;
                }
                let mut in_void = false;
                for (k, &(lump, outer)) in owners.iter().enumerate() {
                    if lump == other && !outer && inside(point, &shells[k], positions, dec)? {
                        in_void = true;
                        break;
                    }
                }
                if !in_void {
                    return Err(Refusal::new(
                        "print.brep-orientation",
                        "encoded material components overlap by nesting",
                    ));
                }
            }
        }
    }
    Ok(())
}

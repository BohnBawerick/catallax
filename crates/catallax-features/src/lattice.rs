//! Lattice cell templates and parameter checks shared by both kernels.
//!
//! RECIPE.md section 4.15 defines a lattice as the beams `P(a (k + n_i))` to
//! `P(a (k + n_j))` over every integer vector `k`. The preset templates below are
//! reduced: a beam that lies entirely in a cell's upper face is the same segment
//! as a beam of the neighbouring cell's lower face, so it is listed once. The
//! reduced template generates exactly the infinite beam set of the full one.
//! Custom cells are kept as written, because nothing proves their faces repeat.

use catallax_numbers::{Rational, Refusal};
use catallax_recipe::{Grade, LatticeCell};

/// One beam template: two node positions in unit-cell coordinates.
#[derive(Clone, PartialEq, Eq, Debug)]
pub struct BeamTemplate {
    /// First node, in unit-cell coordinates.
    pub start: [Rational; 3],
    /// Second node, in unit-cell coordinates.
    pub end: [Rational; 3],
}

/// The cell template after validation.
#[derive(Clone, PartialEq, Eq, Debug)]
pub struct CellTemplate {
    /// Distinct beams of one cell. Empty for a custom cell with no beams.
    pub beams: Vec<BeamTemplate>,
    /// How far, in cells, a beam reaches outside its own unit cube.
    pub overhang: Rational,
    /// Whether the template is a preset whose covering radius is below one cell.
    pub preset: bool,
}

fn q(n: i64) -> Rational {
    Rational::from_i64(n)
}
fn half() -> Rational {
    Rational::from_ratio(1, 2).expect("nonzero denominator")
}
fn corner(bits: usize) -> [Rational; 3] {
    std::array::from_fn(|axis| q(((bits >> axis) & 1) as i64))
}
fn face_centre(axis: usize, upper: bool) -> [Rational; 3] {
    std::array::from_fn(|i| if i == axis { q(upper as i64) } else { half() })
}
fn beam(start: [Rational; 3], end: [Rational; 3]) -> BeamTemplate {
    BeamTemplate { start, end }
}

/// The cell's beams. Presets are the reduced templates; custom cells are checked.
pub fn cell_template(cell: &LatticeCell) -> Result<CellTemplate, Refusal> {
    let mut beams = Vec::new();
    match cell {
        LatticeCell::Cubic => {
            for axis in 0..3 {
                beams.push(beam(corner(0), corner(1 << axis)));
            }
        }
        LatticeCell::BodyCentred => {
            let centre = [half(), half(), half()];
            for bits in 0..8 {
                beams.push(beam(centre.clone(), corner(bits)));
            }
        }
        LatticeCell::FaceCentred | LatticeCell::Octet => {
            for axis in 0..3 {
                let centre = face_centre(axis, false);
                for bits in 0..8 {
                    if (bits >> axis) & 1 == 0 {
                        beams.push(beam(centre.clone(), corner(bits)));
                    }
                }
            }
            if *cell == LatticeCell::Octet {
                for a in 0..3 {
                    for b in (a + 1)..3 {
                        for upper_a in [false, true] {
                            for upper_b in [false, true] {
                                beams.push(beam(face_centre(a, upper_a), face_centre(b, upper_b)));
                            }
                        }
                    }
                }
            }
        }
        LatticeCell::Custom {
            nodes,
            beams: pairs,
        } => {
            let nodes: Vec<[Rational; 3]> = nodes
                .iter()
                .map(|p| p.0.each_ref().map(|n| n.value().clone()))
                .collect();
            for [i, j] in pairs {
                let (Some(start), Some(end)) = (nodes.get(*i as usize), nodes.get(*j as usize))
                else {
                    return Err(Refusal::new(
                        "lattice.node-index",
                        format!("beam [{i}, {j}] names a node the cell does not have"),
                    ));
                };
                if start == end {
                    return Err(Refusal::new(
                        "lattice.impossible",
                        format!("beam [{i}, {j}] joins coincident nodes"),
                    ));
                }
                beams.push(beam(start.clone(), end.clone()));
            }
            beams.sort_by(|a, b| (&a.start, &a.end).cmp(&(&b.start, &b.end)));
            beams.dedup();
            let mut overhang = q(0);
            for b in &beams {
                for p in [&b.start, &b.end] {
                    for x in p {
                        overhang = overhang.max(-x).max(x - &q(1));
                    }
                }
            }
            return Ok(CellTemplate {
                beams,
                overhang,
                preset: false,
            });
        }
    }
    Ok(CellTemplate {
        beams,
        overhang: q(0),
        preset: true,
    })
}

/// The smallest and largest radius a grade can take, or the unsupported grade.
pub fn radius_range(grade: &Grade) -> Result<[Rational; 2], Refusal> {
    match grade {
        Grade::Constant(r) => Ok([r.value().clone(), r.value().clone()]),
        Grade::Linear {
            from,
            to,
            start,
            end,
        } => {
            if from == to {
                return Err(Refusal::new(
                    "lattice.impossible",
                    "a linear grade needs two distinct points",
                ));
            }
            let (a, b) = (start.value().clone(), end.value().clone());
            Ok([a.clone().min(b.clone()), a.max(b)])
        }
        Grade::Distance { .. } => Err(Refusal::new(
            "lattice.unsupported-grade",
            "distance grades need trimmed-face distance; unfinished coverage",
        )),
    }
}

/// Checks the impossible cases of RECIPE.md section 4.15 and returns the cell.
/// The radius check is conservative: a linear grade must stay positive over its
/// whole clamped range, not only where the region reaches.
pub fn validate(
    cell: &LatticeCell,
    size: &Rational,
    grade: &Grade,
) -> Result<CellTemplate, Refusal> {
    if size.sign() != catallax_numbers::Sign::Positive {
        return Err(Refusal::new(
            "lattice.impossible",
            "cell size must be positive",
        ));
    }
    let [low, _] = radius_range(grade)?;
    if low.sign() != catallax_numbers::Sign::Positive {
        return Err(Refusal::new(
            "lattice.impossible",
            "beam radius must be positive everywhere",
        ));
    }
    cell_template(cell)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn preset_templates_have_the_reduced_beam_counts() {
        let counts = [
            (LatticeCell::Cubic, 3),
            (LatticeCell::BodyCentred, 8),
            (LatticeCell::FaceCentred, 12),
            (LatticeCell::Octet, 24),
        ];
        for (cell, expected) in counts {
            let template = cell_template(&cell).unwrap();
            assert_eq!(template.beams.len(), expected, "{cell:?}");
            assert!(template.preset);
            let mut distinct = template.beams.clone();
            distinct.sort_by(|a, b| (&a.start, &a.end).cmp(&(&b.start, &b.end)));
            distinct.dedup();
            assert_eq!(distinct.len(), expected, "{cell:?} repeats a beam");
        }
    }

    #[test]
    fn full_cell_expansion_matches_the_reduced_template() {
        // Expand each reduced template over the 27 neighbouring cells and check
        // that every beam of RECIPE.md's full template appears exactly once.
        fn key(a: &[Rational; 3], b: &[Rational; 3]) -> ([Rational; 3], [Rational; 3]) {
            if a <= b {
                (a.clone(), b.clone())
            } else {
                (b.clone(), a.clone())
            }
        }
        let full_cubic: Vec<([Rational; 3], [Rational; 3])> = (0..8)
            .flat_map(|bits| (0..3).map(move |axis| (bits, axis)))
            .filter(|(bits, axis)| (bits >> axis) & 1 == 0)
            .map(|(bits, axis)| key(&corner(bits), &corner(bits | (1 << axis))))
            .collect();
        assert_eq!(full_cubic.len(), 12);
        let template = cell_template(&LatticeCell::Cubic).unwrap();
        let mut expanded = std::collections::BTreeSet::new();
        for kx in -1..=1 {
            for ky in -1..=1 {
                for kz in -1..=1 {
                    for b in &template.beams {
                        let shift =
                            |p: &[Rational; 3]| [&p[0] + &q(kx), &p[1] + &q(ky), &p[2] + &q(kz)];
                        assert!(expanded.insert(key(&shift(&b.start), &shift(&b.end))));
                    }
                }
            }
        }
        for beam in &full_cubic {
            assert!(expanded.contains(beam), "{beam:?} missing");
        }
    }

    #[test]
    fn impossible_lattices_refuse() {
        let grade = Grade::Constant(1.into());
        assert_eq!(
            validate(&LatticeCell::Cubic, &q(0), &grade)
                .unwrap_err()
                .code,
            "lattice.impossible"
        );
        assert_eq!(
            validate(&LatticeCell::Cubic, &q(2), &Grade::Constant(0.into()))
                .unwrap_err()
                .code,
            "lattice.impossible"
        );
        let custom = LatticeCell::Custom {
            nodes: vec![
                catallax_recipe::Point3([0.into(), 0.into(), 0.into()]),
                catallax_recipe::Point3([0.into(), 0.into(), 0.into()]),
            ],
            beams: vec![[0, 1]],
        };
        assert_eq!(
            validate(&custom, &q(2), &grade).unwrap_err().code,
            "lattice.impossible"
        );
        let custom = LatticeCell::Custom {
            nodes: vec![catallax_recipe::Point3([0.into(), 0.into(), 0.into()])],
            beams: vec![[0, 3]],
        };
        assert_eq!(
            validate(&custom, &q(2), &grade).unwrap_err().code,
            "lattice.node-index"
        );
    }
}

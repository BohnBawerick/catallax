use super::{
    Decider, Domain, Grid, Rational, Refusal, Value, add, less, mul, restricted_net, restriction,
    sub,
};

mod hull;

// Saturation and nonhierarchical rectangles use direct exact restriction.
const MAX_CACHED_CELLS: usize = 64;
type Cuts = [[Rational; 2]; 2];

#[derive(Clone)]
struct Split {
    axis: usize,
    middle: Rational,
    children: [usize; 2],
}

struct Cell {
    parent: Option<usize>,
    // Weight, du.xyz, dv.xyz; each component is populated only when requested.
    grids: [Option<Grid>; 7],
    split: Option<Split>,
}

pub(super) struct Prepared {
    here: Option<usize>,
    cuts: Option<Cuts>,
}

#[derive(Default)]
pub(super) struct RestrictionCache {
    root: Option<Split>,
    cells: Vec<Cell>,
    monotone: [Option<[bool; 2]>; 7],
}

impl RestrictionCache {
    fn division(&self, here: Option<usize>) -> Option<&Split> {
        match here {
            None => self.root.as_ref(),
            Some(index) => self.cells[index].split.as_ref(),
        }
    }

    fn side(
        cell: &Domain,
        axis: usize,
        middle: &Rational,
        dec: &mut Decider,
    ) -> Result<Option<usize>, Refusal> {
        if !less(middle, &cell[axis][1], dec)? {
            Ok(Some(0))
        } else if !less(&cell[axis][0], middle, dec)? {
            Ok(Some(1))
        } else {
            Ok(None)
        }
    }

    pub(super) fn prepare(
        &mut self,
        domain: &Domain,
        cell: &Domain,
        original_cuts: Cuts,
        dec: &mut Decider,
    ) -> Result<Prepared, Refusal> {
        let mut bounds = domain.clone();
        let mut here = None;
        let half = Rational::from_ratio(1, 2).expect("nonzero denominator");
        loop {
            let mut same = [false; 2];
            for axis in 0..2 {
                same[axis] = <Rational as Value>::same(&bounds[axis][0], &cell[axis][0], dec)?
                    && <Rational as Value>::same(&bounds[axis][1], &cell[axis][1], dec)?;
            }
            if same == [true; 2] {
                return Ok(Prepared { here, cuts: None });
            }
            let next = if let Some(next) = self.division(here) {
                Self::side(cell, next.axis, &next.middle, dec)?.map(|side| (next.clone(), side))
            } else if self.cells.len() + 2 <= MAX_CACHED_CELLS {
                let axis = usize::from(same[0]);
                let sum = add(&bounds[axis][0], &bounds[axis][1], dec)?;
                let middle = mul(&sum, &half, dec)?;
                if let Some(side) = Self::side(cell, axis, &middle, dec)? {
                    let children = [self.cells.len(), self.cells.len() + 1];
                    for _ in 0..2 {
                        self.cells.push(Cell {
                            parent: here,
                            grids: std::array::from_fn(|_| None),
                            split: None,
                        });
                    }
                    let next = Split {
                        axis,
                        middle,
                        children,
                    };
                    match here {
                        None => self.root = Some(next.clone()),
                        Some(index) => self.cells[index].split = Some(next.clone()),
                    }
                    Some((next, side))
                } else {
                    None
                }
            } else {
                None
            };
            if let Some((next, side)) = next {
                bounds[next.axis][1 - side] = next.middle;
                here = Some(next.children[side]);
                continue;
            }
            let cuts = if here.is_none() {
                original_cuts
            } else {
                let widths = [
                    sub(&bounds[0][1], &bounds[0][0], dec)?,
                    sub(&bounds[1][1], &bounds[1][0], dec)?,
                ];
                restriction(&bounds, &widths, cell, dec)?
            };
            return Ok(Prepared {
                here,
                cuts: Some(cuts),
            });
        }
    }

    fn fill(
        &mut self,
        root: &Grid,
        component: usize,
        index: usize,
        dec: &mut Decider,
    ) -> Result<(), Refusal> {
        if self.cells[index].grids[component].is_some() {
            return Ok(());
        }
        // Each division creates two cells, bounding the ancestor stack.
        let mut pending = [0; MAX_CACHED_CELLS / 2];
        let mut count = 0;
        let mut ancestor = Some(index);
        while let Some(index) = ancestor {
            if self.cells[index].grids[component].is_some() {
                break;
            }
            pending[count] = index;
            count += 1;
            ancestor = self.cells[index].parent;
        }
        let half = Rational::from_ratio(1, 2).expect("nonzero denominator");
        for &index in pending[..count].iter().rev() {
            let parent = self.cells[index].parent;
            let division = self
                .division(parent)
                .expect("cached child has a parent split");
            let axis = division.axis;
            let children = division.children;
            let grid = match parent {
                None => root,
                Some(index) => self.cells[index].grids[component]
                    .as_ref()
                    .expect("prepared parent grid"),
            };
            let parts = grid.split(axis, &half, dec)?;
            for (child, grid) in children.into_iter().zip(parts) {
                self.cells[child].grids[component] = Some(grid);
            }
        }
        Ok(())
    }

    pub(super) fn range(
        &mut self,
        root: &Grid,
        component: usize,
        prepared: &Prepared,
        dec: &mut Decider,
    ) -> Result<[Rational; 2], Refusal> {
        if root.deg == [0, 0] {
            return hull::range(root, [true; 2], dec);
        }
        let monotone = if let Some(monotone) = self.monotone[component] {
            monotone
        } else {
            let monotone = hull::monotone(root, dec)?;
            self.monotone[component] = Some(monotone);
            monotone
        };
        let grid = match prepared.here {
            None => root,
            Some(index) => {
                self.fill(root, component, index, dec)?;
                self.cells[index].grids[component]
                    .as_ref()
                    .expect("prepared component grid")
            }
        };
        match &prepared.cuts {
            None => hull::range(grid, monotone, dec),
            Some(cuts) => {
                let net = restricted_net(grid, cuts, dec)?;
                hull::range(net.as_ref(), monotone, dec)
            }
        }
    }
}

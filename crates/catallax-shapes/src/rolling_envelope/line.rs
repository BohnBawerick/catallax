//! Complete transverse crossings of a line with the represented carrier.
//!
//! Over the chart, two planes through the line give a continuous piecewise
//! polynomial map G whose zeros are the line's hits. Bernstein hulls discard
//! each cell that keeps G from zero, so every hit stays in a surviving cell.
//! Surviving cells join into rectangles. Where a preconditioner Y gives
//! |I - Y J| < 1 over a rectangle's interval Jacobian, G is injective there:
//! the mean of J along any segment, knots included, stays in that hull. Such a
//! rectangle owns at most one hit, and its Jacobian keeps one sign, so the line
//! crosses the carrier there. A rectangle with no such Y refuses.
//!
//! The mean-value image of a rectangle then decides its hit. Clear of every
//! knot and row, the hit gets a certified root of the span's own equations. On
//! a knot or row it is an exact point of that section, or the image narrows off
//! it. Nothing is accepted from a nearby enclosure alone.

use super::*;
use catallax_numbers::{BoxRoots, MPoly, PolySystem, RealPoly, RootBox};

// Proposal bounds, not tolerances: running out refuses and proves nothing.
const MAX_DEPTH: usize = 32;
const MAX_CELLS: usize = 256;
const MAX_NARROWING: usize = 16;

type Rect = [[Rational; 2]; 2];

fn unresolved_line(detail: &str) -> Refusal {
    Refusal::new("intersect.parametric-surface", detail)
}

fn meets_rect(a: &Rect, c: &Rect) -> bool {
    (0..2).all(|axis| a[axis][0] <= c[axis][1] && c[axis][0] <= a[axis][1])
}
fn common(a: &Rect, c: &Rect) -> Option<Rect> {
    meets_rect(a, c).then(|| {
        std::array::from_fn(|axis| {
            [
                a[axis][0].clone().max(c[axis][0].clone()),
                a[axis][1].clone().min(c[axis][1].clone()),
            ]
        })
    })
}
fn join_rect(a: &Rect, c: &Rect) -> Rect {
    std::array::from_fn(|axis| {
        [
            a[axis][0].clone().min(c[axis][0].clone()),
            a[axis][1].clone().max(c[axis][1].clone()),
        ]
    })
}
fn width(rect: &Rect) -> Rational {
    (&rect[0][1] - &rect[0][0]).max(&rect[1][1] - &rect[1][0])
}

impl Patch {
    // The sub-patch over a rectangle inside this cell, by two De Casteljau cuts
    // along each axis. A degenerate rectangle gives the net of one point.
    fn clip(&self, rect: &Rect, bits: u32) -> Self {
        let mut net = self.net.clone();
        for (axis, [from, to]) in rect.iter().enumerate() {
            let [lo, hi] = &self.cell[axis];
            let cut = |net: &[[Ball; 4]], ratio: Rational| {
                let right = b(&ratio, bits);
                let left = b(&(Rational::one() - ratio), bits);
                divide_net(net, [DEGREE, 2], axis, |a: &[Ball; 4], c: &[Ball; 4]| {
                    std::array::from_fn(|i| &(&a[i] * &left) + &(&c[i] * &right))
                })
            };
            if from > lo {
                let [_, right] = cut(&net, (from - lo) / (hi - lo));
                net = right;
            }
            if to < hi && hi > from {
                let [left, _] = cut(&net, (to - from) / (hi - from));
                net = left;
            }
        }
        Self {
            cell: rect.clone(),
            net,
        }
    }
}

// What a contact may lie on: the knot section of a node, or a long boundary row.
#[derive(Clone, Copy, PartialEq, Debug)]
enum Section {
    Knot(usize),
    Mouth,
    Rail,
}

// A rectangle owning at most one hit: `cells` holds any hit, and `domain`, which
// holds `cells` with a margin, is where the mean-value image `image` was proved.
struct Piece {
    cells: Rect,
    domain: Rect,
    image: Rect,
}

struct Contact {
    parameter: Real,
    local: [Ball; 3],
}

struct Query<'a> {
    envelope: &'a RollingEnvelope,
    bits: u32,
    origin: &'a Point3,
    direction: &'a Vector3,
    // The same line in the definition frame, exactly and enclosed.
    exact: [[Real; 3]; 2],
    enclosed: [[Ball; 3]; 2],
    // Two independent exact normals of the line, enclosed, with their offsets.
    normals: [[Real; 3]; 2],
    enclosed_normals: [[Ball; 3]; 2],
    offsets: [Ball; 2],
    spans: Vec<Patch>,
}

impl RollingEnvelope {
    /// Every point where a whole line meets the closed finite carrier, in line
    /// order. Each is an exact transverse crossing and appears once, including
    /// at a knot and on a boundary row. A tangency, or a contact this cannot
    /// separate, refuses. It never reports a hit or a miss from an enclosure.
    pub fn line_hits(
        &self,
        origin: &Point3,
        direction: &Vector3,
        dec: &mut Decider,
    ) -> Result<Vec<Point3>, Refusal> {
        dec.check_deadline()?;
        math::validate_point(origin, dec)?;
        math::validate_vector(direction, dec)?;
        let query = Query::new(self, origin, direction, dec)?;
        let pieces = query.pieces(dec)?;
        let mut hits: Vec<Real> = Vec::with_capacity(pieces.len());
        for index in 0..pieces.len() {
            if let Some(parameter) = query.decide(&pieces, index, dec)? {
                let at = hits
                    .iter()
                    .try_fold(0, |at, hit| -> Result<usize, Refusal> {
                        Ok(if math::sign(&(&parameter - hit), dec)? == Sign::Positive {
                            at + 1
                        } else {
                            at
                        })
                    })?;
                hits.insert(at, parameter);
            }
        }
        Ok(hits.iter().map(|hit| query.world(hit)).collect())
    }
}

impl<'a> Query<'a> {
    fn new(
        envelope: &'a RollingEnvelope,
        origin: &'a Point3,
        direction: &'a Vector3,
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        let bits = envelope.query_bits(dec)?;
        let def = &envelope.data.definition;
        let frame = LocalFrame::budgeted(def, bits, dec)?;
        let exact = [
            def.local(origin),
            [
                direction.dot(&def.toward),
                direction.dot(&def.axis),
                direction.dot(&def.across),
            ],
        ];
        let along = math::budgeted_vector(direction, bits, dec)?;
        let enclosed = [
            frame.local(&math::budgeted_point(origin, bits, dec)?),
            frame.axes.each_ref().map(|axis| dot(&along, axis)),
        ];
        // Crossing with the axis the direction leans on least keeps both normals long.
        let least = (0..3)
            .min_by_key(|&i| magnitude(&enclosed[1][i]))
            .expect("three axes");
        let unit: [Real; 3] = std::array::from_fn(|i| {
            if i == least {
                Real::one()
            } else {
                Real::zero()
            }
        });
        let first = cross(&exact[1], &unit);
        let second = cross(&exact[1], &first);
        let normals = [first, second];
        let mut enclosed_normals = Vec::with_capacity(2);
        let mut offsets = Vec::with_capacity(2);
        for normal in &normals {
            let normal: [Ball; 3] = [
                dec.enclose_budgeted_at(&normal[0], bits)?,
                dec.enclose_budgeted_at(&normal[1], bits)?,
                dec.enclose_budgeted_at(&normal[2], bits)?,
            ];
            offsets.push(dot(&normal, &enclosed[0]));
            enclosed_normals.push(normal);
        }
        Ok(Self {
            envelope,
            bits,
            origin,
            direction,
            exact,
            enclosed,
            normals,
            enclosed_normals: enclosed_normals.try_into().expect("two normals"),
            offsets: offsets.try_into().expect("two normals"),
            spans: Patch::at(&envelope.data.carrier, &envelope.data.net, bits),
        })
    }

    fn world(&self, parameter: &Real) -> Point3 {
        self.origin.plus(&self.direction.scaled(parameter))
    }

    // G on the controls of a patch. Positive weights keep each hit a zero of G.
    fn equations(&self, patch: &Patch) -> Vec<[Ball; 2]> {
        patch
            .net
            .iter()
            .map(|h| {
                std::array::from_fn(|k| {
                    let along = dot(
                        &[h[0].clone(), h[1].clone(), h[2].clone()],
                        &self.enclosed_normals[k],
                    );
                    &along - &(&self.offsets[k] * &h[3])
                })
            })
            .collect()
    }

    // Rectangles that together hold every hit, each with its mean-value image.
    fn pieces(&self, dec: &mut Decider) -> Result<Vec<Piece>, Refusal> {
        let mut cells = self.spans.clone();
        for depth in 0..MAX_DEPTH {
            let mut kept = Vec::with_capacity(cells.len());
            for cell in cells {
                dec.check_deadline()?;
                let values = self.equations(&cell);
                let apart = (0..2).any(|k| {
                    values.iter().all(|value| value[k].lo() > Rational::zero())
                        || values.iter().all(|value| value[k].hi() < Rational::zero())
                });
                if !apart {
                    kept.push(cell);
                }
            }
            cells = kept;
            if cells.is_empty() {
                return Ok(Vec::new());
            }
            if cells.len() > MAX_CELLS {
                break;
            }
            if depth >= 2
                && let Some(pieces) = self.resolved(&cells, dec)?
            {
                return Ok(pieces);
            }
            cells = cells.into_iter().flat_map(Patch::split).collect();
        }
        Err(unresolved_line(
            "a line contact on this finite spline is not a separated transverse crossing",
        ))
    }

    fn resolved(&self, cells: &[Patch], dec: &mut Decider) -> Result<Option<Vec<Piece>>, Refusal> {
        let chart = [
            [Rational::from_i64(-1), Rational::one()],
            [Rational::zero(), Rational::one()],
        ];
        // Each cell's margin is its own size, so a hit sits inside its rectangle.
        let mut parts: Vec<(Rect, Rect)> = cells
            .iter()
            .map(|cell| {
                let wider = std::array::from_fn(|axis| {
                    let [lo, hi] = &cell.cell[axis];
                    let margin = hi - lo;
                    [
                        (lo - &margin).max(chart[axis][0].clone()),
                        (hi + &margin).min(chart[axis][1].clone()),
                    ]
                });
                (cell.cell.clone(), wider)
            })
            .collect();
        let mut index = 0;
        while index < parts.len() {
            let mut other = index + 1;
            let mut merged = false;
            while other < parts.len() {
                if meets_rect(&parts[index].1, &parts[other].1) {
                    let (cells, domain) = parts.swap_remove(other);
                    parts[index] = (
                        join_rect(&parts[index].0, &cells),
                        join_rect(&parts[index].1, &domain),
                    );
                    merged = true;
                } else {
                    other += 1;
                }
            }
            if !merged {
                index += 1;
            }
        }
        let mut pieces = Vec::with_capacity(parts.len());
        for (cells, domain) in parts {
            dec.check_deadline()?;
            let Some(image) = self.image(&domain) else {
                return Ok(None);
            };
            pieces.push(Piece {
                cells,
                domain,
                image,
            });
        }
        Ok(Some(pieces))
    }

    // m - Y G(m) + (I - Y J)(R - m), or `None` without a contracting Y on R.
    // J's hull covers every span of R, so the bound holds across knots.
    fn image(&self, rect: &Rect) -> Option<Rect> {
        let bits = self.bits;
        let mut hull: [[Option<Ball>; 2]; 2] = [[None, None], [None, None]];
        let mut centre = None;
        let middle: [Rational; 2] =
            std::array::from_fn(|axis| (&rect[axis][0] + &rect[axis][1]).mul_pow2(-1));
        for span in &self.spans {
            let Some(part) = common(&span.cell, rect) else {
                continue;
            };
            if span.cell[0][0] <= middle[0] && middle[0] <= span.cell[0][1] && centre.is_none() {
                let point =
                    std::array::from_fn(|axis| [middle[axis].clone(), middle[axis].clone()]);
                centre = Some(self.equations(&span.clip(&point, bits))[0].clone());
            }
            if part[0][0] == part[0][1] || part[1][0] == part[1][1] {
                continue;
            }
            let values = self.equations(&span.clip(&part, bits));
            let steps = [
                b(&(q(DEGREE as i64) / (&part[0][1] - &part[0][0])), bits),
                b(&(q(2) / (&part[1][1] - &part[1][0])), bits),
            ];
            let at = |i: usize, j: usize| &values[i * 3 + j];
            for i in 0..=DEGREE {
                for j in 0..3 {
                    for (k, row) in hull.iter_mut().enumerate() {
                        let mut include = |slope: Ball, axis: usize| {
                            let entry = &mut row[axis];
                            *entry = Some(match entry.take() {
                                Some(wide) => wide.hull(&slope),
                                None => slope,
                            });
                        };
                        if i < DEGREE {
                            include(&(&at(i + 1, j)[k] - &at(i, j)[k]) * &steps[0], 0);
                        }
                        if j < 2 {
                            include(&(&at(i, j + 1)[k] - &at(i, j)[k]) * &steps[1], 1);
                        }
                    }
                }
            }
        }
        let [[Some(a), Some(c)], [Some(d), Some(e)]] = hull else {
            return None;
        };
        let jacobian = [[a, c], [d, e]];
        let value: [Ball; 2] = centre?;
        let mid = [
            [jacobian[0][0].mid(), jacobian[0][1].mid()],
            [jacobian[1][0].mid(), jacobian[1][1].mid()],
        ];
        let determinant = &mid[0][0] * &mid[1][1] - &mid[0][1] * &mid[1][0];
        let inverse = determinant.recip()?;
        let y = [
            [
                b(&(&mid[1][1] * &inverse), bits),
                b(&-(&mid[0][1] * &inverse), bits),
            ],
            [
                b(&-(&mid[1][0] * &inverse), bits),
                b(&(&mid[0][0] * &inverse), bits),
            ],
        ];
        let spread: [Ball; 2] = std::array::from_fn(|axis| {
            b(&rect[axis][0], bits).hull(&b(&rect[axis][1], bits)) - b(&middle[axis], bits)
        });
        let mut image = rect.clone();
        for row in 0..2 {
            let mut total = Rational::zero();
            let mut value_at =
                b(&middle[row], bits) - (&(&y[row][0] * &value[0]) + &(&y[row][1] * &value[1]));
            for column in 0..2 {
                let identity = if row == column { 1 } else { 0 };
                let residual = b(&q(identity), bits)
                    - (&(&y[row][0] * &jacobian[0][column]) + &(&y[row][1] * &jacobian[1][column]));
                total = total + magnitude(&residual);
                value_at = value_at + &residual * &spread[column];
            }
            if total >= Rational::one() {
                return None;
            }
            image[row] = [value_at.lo(), value_at.hi()];
        }
        Some(image)
    }

    fn touched(&self, image: &Rect) -> Vec<Section> {
        let mut touched: Vec<Section> = self
            .envelope
            .data
            .nodes
            .iter()
            .enumerate()
            .filter(|(_, node)| image[0][0] <= node.s && node.s <= image[0][1])
            .map(|(index, _)| Section::Knot(index))
            .collect();
        if image[1][0] <= Rational::zero() {
            touched.push(Section::Rail);
        }
        if image[1][1] >= Rational::one() {
            touched.push(Section::Mouth);
        }
        touched
    }

    // The line parameter of the hit of one piece, if it has one.
    fn decide(
        &self,
        pieces: &[Piece],
        index: usize,
        dec: &mut Decider,
    ) -> Result<Option<Real>, Refusal> {
        let piece = &pieces[index];
        let mut domain = piece.domain.clone();
        let mut image = piece.image.clone();
        let mut tried = Vec::new();
        for _ in 0..MAX_NARROWING {
            dec.check_deadline()?;
            let Some(core) = common(&image, &piece.cells) else {
                return Ok(None);
            };
            let touched = self.touched(&core);
            if touched.is_empty()
                && width(&core) < Rational::one().mul_pow2(-i64::from(self.bits / 4))
            {
                return self.interior(piece, &core, dec);
            }
            for section in touched {
                if tried.contains(&section) {
                    continue;
                }
                tried.push(section);
                if let Some(contact) = self.owned_contact(section, pieces, index, dec)? {
                    return Ok(Some(contact));
                }
            }
            let Some(next) = common(&image, &domain) else {
                return Ok(None);
            };
            let Some(narrower) = self.image(&next) else {
                break;
            };
            if next == domain && narrower == image {
                break;
            }
            domain = next;
            image = narrower;
        }
        Err(unresolved_line(
            "a line contact on this finite spline is not separated from a chart knot or boundary",
        ))
    }

    // The contact on a section that belongs to this piece. Every hit lies in some
    // piece, so a contact clear of all the others is the piece's unique hit.
    fn owned_contact(
        &self,
        section: Section,
        pieces: &[Piece],
        index: usize,
        dec: &mut Decider,
    ) -> Result<Option<Real>, Refusal> {
        for parameter in self.contacts(section, dec)? {
            let contact = self.contact(parameter, dec)?;
            let mut owners = Vec::new();
            for (other, piece) in pieces.iter().enumerate() {
                if self.may_hold(&piece.domain, &contact)? {
                    owners.push(other);
                }
            }
            match owners.as_slice() {
                [owner] if *owner == index => return Ok(Some(contact.parameter)),
                [_] => {}
                [] => {
                    return Err(unresolved_line(
                        "an exact line contact has no chart cell on this finite spline",
                    ));
                }
                _ if owners.contains(&index) => {
                    return Err(unresolved_line(
                        "line contacts on this finite spline are not separated",
                    ));
                }
                _ => {}
            }
        }
        Ok(None)
    }

    fn contact(&self, parameter: Real, dec: &mut Decider) -> Result<Contact, Refusal> {
        let along = dec.enclose_budgeted_at(&parameter, self.bits)?;
        Ok(Contact {
            parameter,
            local: std::array::from_fn(|i| &self.enclosed[0][i] + &(&along * &self.enclosed[1][i])),
        })
    }

    fn may_hold(&self, rect: &Rect, contact: &Contact) -> Result<bool, Refusal> {
        for span in &self.spans {
            let Some(part) = common(&span.cell, rect) else {
                continue;
            };
            let hull = span.clip(&part, self.bits).hull()?;
            if (0..3).all(|i| meets(&hull[i], &contact.local[i])) {
                return Ok(true);
            }
        }
        Ok(false)
    }

    // Exact line parameters where the line meets one section of the chart.
    fn contacts(&self, section: Section, dec: &mut Decider) -> Result<Vec<Real>, Refusal> {
        match section {
            Section::Knot(index) => self.knot_contacts(index, dec),
            Section::Mouth => {
                let [origin, direction] = &self.exact;
                if math::sign(&direction[1], dec)? == Sign::Zero {
                    if math::sign(&origin[1], dec)? != Sign::Zero {
                        return Ok(Vec::new());
                    }
                    return Err(unresolved_line(
                        "a line in the mouth plane needs a circle contact solver",
                    ));
                }
                self.on_boundary(vec![-(&origin[1] / &direction[1])], dec)
            }
            Section::Rail => {
                let [o, d] = &self.exact;
                let body = &self.envelope.data.definition.body;
                let quadratic = RealPoly::from_terms(
                    1,
                    vec![
                        (vec![0], &(&o[0].square() + &o[2].square()) - &body.square()),
                        (
                            vec![1],
                            &Real::from_i64(2) * &(&(&o[0] * &d[0]) + &(&o[2] * &d[2])),
                        ),
                        (vec![2], &d[0].square() + &d[2].square()),
                    ],
                );
                if quadratic.is_zero(dec)? {
                    return Err(unresolved_line(
                        "a line on the rail cylinder needs a tangency solver",
                    ));
                }
                self.on_boundary(quadratic.real_roots(dec)?, dec)
            }
        }
    }

    fn on_boundary(&self, candidates: Vec<Real>, dec: &mut Decider) -> Result<Vec<Real>, Refusal> {
        let mut kept = Vec::with_capacity(candidates.len());
        for candidate in candidates {
            if self.envelope.on_boundary(&self.world(&candidate), dec)? {
                kept.push(candidate);
            }
        }
        Ok(kept)
    }

    // The section at a knot is one rational quadratic in v with exact controls,
    // so its contacts are the common roots of the two equations of the line.
    fn knot_contacts(&self, index: usize, dec: &mut Decider) -> Result<Vec<Real>, Refusal> {
        let data = &self.envelope.data;
        let (pair, column) = if index + 1 < data.nodes.len() {
            (&data.nodes[index..index + 2], 0)
        } else {
            (&data.nodes[index - 1..index + 1], DEGREE)
        };
        let t = data.definition.transition();
        let section: [[Real; 4]; 3] = span_rows(&data.definition, pair)
            .map(|row| std::array::from_fn(|k| row[k].at(column, &t.t)));
        let [origin, direction] = &self.exact;
        let along = |normal: &[Real; 3], row: &[Real; 4]| -> Real {
            let offset = dot_real(normal, origin);
            &dot_real(normal, &[row[0].clone(), row[1].clone(), row[2].clone()])
                - &(&offset * &row[3])
        };
        let equation = |normal: &[Real; 3]| {
            let [c0, c1, c2] = [
                along(normal, &section[0]),
                along(normal, &section[1]),
                along(normal, &section[2]),
            ];
            let two = Real::from_i64(2);
            let linear = &two * &(&c1 - &c0);
            let quadratic = &(&c0 - &(&two * &c1)) + &c2;
            RealPoly::from_terms(
                1,
                vec![(vec![0], c0), (vec![1], linear), (vec![2], quadratic)],
            )
        };
        let shared = equation(&self.normals[0]).gcd(&equation(&self.normals[1]), dec)?;
        if shared.is_zero(dec)? {
            return Err(unresolved_line(
                "a line along a knot section needs an overlap solver",
            ));
        }
        let mut parameters = Vec::new();
        for v in shared.real_roots(dec)? {
            if math::sign(&v, dec)? == Sign::Negative
                || math::sign(&(&v - &Real::one()), dec)? == Sign::Positive
            {
                continue;
            }
            let rest = &Real::one() - &v;
            let basis = [
                rest.square(),
                &Real::from_i64(2) * &(&rest * &v),
                v.square(),
            ];
            let weighted = |k: usize| {
                (0..3).fold(Real::zero(), |sum, j| {
                    plus(&sum, &times(&section[j][k], &basis[j]))
                })
            };
            let weight = weighted(3);
            let relative: [Real; 3] =
                std::array::from_fn(|k| &(&weighted(k) / &weight) - &origin[k]);
            parameters.push(&dot_real(&relative, direction) / &dot_real(direction, direction));
        }
        Ok(parameters)
    }

    // A hit clear of every knot and row: one root of the span's exact equations.
    fn interior(
        &self,
        piece: &Piece,
        core: &Rect,
        dec: &mut Decider,
    ) -> Result<Option<Real>, Refusal> {
        let data = &self.envelope.data;
        let index = data
            .nodes
            .windows(2)
            .position(|pair| pair[0].s <= core[0][0] && core[0][1] <= pair[1].s)
            .expect("a core clear of every knot lies in one span");
        let span = [&data.nodes[index].s, &data.nodes[index + 1].s];
        // Room for the certificate, stopping short of the rows, knots and domain.
        let limits = [
            [
                span[0].clone().max(piece.domain[0][0].clone()),
                span[1].clone().min(piece.domain[0][1].clone()),
            ],
            [
                Rational::zero().max(piece.domain[1][0].clone()),
                Rational::one().min(piece.domain[1][1].clone()),
            ],
        ];
        let room: [Rational; 2] = std::array::from_fn(|axis| {
            let reach = (&core[axis][0] - &limits[axis][0]).min(&limits[axis][1] - &core[axis][1]);
            (reach / q(2))
                .min((&core[axis][1] - &core[axis][0]) * q(4) + Rational::one().mul_pow2(-24))
        });
        if room.iter().any(|room| room.sign() != Sign::Positive) {
            return Err(unresolved_line(
                "a line contact leaves no certificate room on this finite spline",
            ));
        }
        let box_: Rect = std::array::from_fn(|axis| {
            [&core[axis][0] - &room[axis], &core[axis][1] + &room[axis]]
        });
        let Some((system, candidate)) = self.span_system(index, &box_, dec)? else {
            return Err(unresolved_line(
                "this finite spline crossing has no rational equations",
            ));
        };
        match dec.certify_root_at(&system, &candidate, self.bits) {
            Outcome::Proven(BoxRoots::One(root)) => Ok(Some(root.coordinate(0))),
            Outcome::Proven(BoxRoots::None) => Ok(None),
            Outcome::NotExact(n) => Err(Refusal::new(n.code, n.detail)),
            Outcome::ProvenCoincident(_) => {
                Err(unresolved_line("a root certificate reported coincidence"))
            }
            Outcome::Refused(refusal) => Err(refusal),
        }
    }

    // The exact equations of one span in (line parameter, span parameter, v,
    // t, end weight). The carrier's only irrationals are the transition tangent
    // t and the shared end weight, so they enter as variables fixed by their
    // squares, and each root coordinate stays a leaf of a rational system.
    fn span_system(
        &self,
        index: usize,
        rect: &Rect,
        dec: &mut Decider,
    ) -> Result<Option<(Arc<PolySystem>, RootBox)>, Refusal> {
        let data = &self.envelope.data;
        let def = &data.definition;
        let t = def.transition();
        let rational = |value: &Real| value.rational_ref().cloned();
        let vector = |value: &[Real; 3]| -> Option<[Rational; 3]> {
            Some([
                rational(&value[0])?,
                rational(&value[1])?,
                rational(&value[2])?,
            ])
        };
        let (Some(origin), Some(direction), Some(first), Some(second)) = (
            vector(&self.exact[0]),
            vector(&self.exact[1]),
            vector(&self.normals[0]),
            vector(&self.normals[1]),
        ) else {
            return Ok(None);
        };
        let Some(tangent_square) = rational(&t.square) else {
            return Ok(None);
        };
        let pair = &data.nodes[index..index + 2];
        let ends = pair.iter().any(|node| node.s.abs() == Rational::one());
        let with_end = |value: i64| -> Vec<Node> {
            pair.iter()
                .map(|node| {
                    let mut node = node.clone();
                    if node.s.abs() == Rational::one() {
                        node.values[WEIGHT] = Real::from_i64(value);
                    }
                    node
                })
                .collect()
        };
        let nvars = if ends { 5 } else { 4 };
        let (Some(base), Some(end)) = (
            span_terms(def, &with_end(0), nvars),
            span_terms(def, &with_end(1), nvars),
        ) else {
            return Ok(None);
        };
        let homogeneous: [MPoly; 4] = std::array::from_fn(|k| {
            if ends {
                let weight = MPoly::var(nvars, 4);
                &base[k] + &(&(&end[k] - &base[k]) * &weight)
            } else {
                base[k].clone()
            }
        });
        let combine = |normal: &[Rational; 3]| {
            let along = (0..3).fold(MPoly::zero(nvars), |sum, k| {
                &sum + &homogeneous[k].scale(&normal[k])
            });
            let offset = rational_dot(normal, &origin);
            &along - &homogeneous[3].scale(&offset)
        };
        let lambda = MPoly::var(nvars, 0);
        let length = rational_dot(&direction, &direction);
        let square = |variable: usize, value: &Rational| {
            let x = MPoly::var(nvars, variable);
            &(&x * &x) - &MPoly::constant(nvars, value.clone())
        };
        let mut equations = vec![
            combine(&first),
            combine(&second),
            &(&lambda * &homogeneous[3]).scale(&length) - &combine(&direction),
            square(3, &tangent_square),
        ];
        let tangent = dec.enclose_budgeted_at(&t.t, self.bits)?;
        let slack = Rational::one().mul_pow2(-i64::from(self.bits / 2));
        let mut lo = vec![Rational::zero(); nvars];
        let mut hi = vec![Rational::zero(); nvars];
        lo[3] = tangent.lo() - slack.clone();
        hi[3] = tangent.hi() + slack.clone();
        if ends {
            let t_mouth = rational(&t.mouth);
            let Some(mouth) = t_mouth else {
                return Ok(None);
            };
            let value = (Rational::one() + tangent_square.clone() * mouth.clone() * mouth)
                / (q(2) * (Rational::one() + tangent_square.clone()));
            equations.push(square(4, &value));
            let weight = pair
                .iter()
                .find(|node| node.s.abs() == Rational::one())
                .expect("an end span has an end node");
            let weight = dec.enclose_budgeted_at(&weight.values[WEIGHT], self.bits)?;
            lo[4] = weight.lo() - slack.clone();
            hi[4] = weight.hi() + slack.clone();
        }
        let width = &data.nodes[index + 1].s - &data.nodes[index].s;
        lo[1] = (&rect[0][0] - &data.nodes[index].s) / width.clone();
        hi[1] = (&rect[0][1] - &data.nodes[index].s) / width;
        lo[2] = rect[1][0].clone();
        hi[2] = rect[1][1].clone();
        // Along the line, the rectangle's whole image bounds the parameter.
        let hull = self.spans[index].clip(rect, self.bits).hull()?;
        let relative: [Ball; 3] = std::array::from_fn(|i| &hull[i] - &self.enclosed[0][i]);
        let along = dot(&relative, &self.enclosed[1]);
        let scale = dot(&self.enclosed[1], &self.enclosed[1]);
        let Some(parameter) = along.checked_div(&scale) else {
            return Ok(None);
        };
        lo[0] = parameter.lo() - slack.clone();
        hi[0] = parameter.hi() + slack;
        let Some(system) = PolySystem::new(equations) else {
            return Ok(None);
        };
        Ok(RootBox::new(lo, hi).map(|candidate| (Arc::new(system), candidate)))
    }
}

fn cross(a: &[Real; 3], c: &[Real; 3]) -> [Real; 3] {
    std::array::from_fn(|i| {
        let (j, k) = ((i + 1) % 3, (i + 2) % 3);
        &times(&a[j], &c[k]) - &times(&a[k], &c[j])
    })
}
fn dot_real(a: &[Real; 3], c: &[Real; 3]) -> Real {
    (0..3).fold(Real::zero(), |sum, i| plus(&sum, &times(&a[i], &c[i])))
}
fn rational_dot(a: &[Rational; 3], c: &[Rational; 3]) -> Rational {
    (0..3).fold(Rational::zero(), |sum, i| sum + &a[i] * &c[i])
}

// One span's homogeneous local coordinates as polynomials in the span
// parameter (variable 1), v (2) and t (3), or `None` off rational node data.
fn span_terms(def: &RollingEnvelopeDefinition, pair: &[Node], nvars: usize) -> Option<[MPoly; 4]> {
    let rows = span_rows(def, pair);
    let bernstein = |variable: usize, degree: usize, i: usize| -> MPoly {
        let x = MPoly::var(nvars, variable);
        let rest = &MPoly::constant(nvars, Rational::one()) - &x;
        let mut value = MPoly::constant(nvars, binomial(degree, i));
        for _ in 0..i {
            value = &value * &x;
        }
        for _ in i..degree {
            value = &value * &rest;
        }
        value
    };
    let along: Vec<MPoly> = (0..=DEGREE).map(|i| bernstein(1, DEGREE, i)).collect();
    let across: Vec<MPoly> = (0..3).map(|j| bernstein(2, 2, j)).collect();
    let tangent = MPoly::var(nvars, 3);
    let mut components = Vec::with_capacity(4);
    for k in 0..4 {
        let mut sum = MPoly::zero(nvars);
        for (j, row) in rows.iter().enumerate() {
            let basis = &across[j];
            let mut line = MPoly::zero(nvars);
            for (i, along) in along.iter().enumerate() {
                let even = row[k].even[i].rational_ref()?;
                let odd = row[k].odd[i].rational_ref()?;
                let coefficient = &MPoly::constant(nvars, even.clone()) + &tangent.scale(odd);
                line = &line + &(along * &coefficient);
            }
            sum = &sum + &(&line * basis);
        }
        components.push(sum);
    }
    components.try_into().ok()
}

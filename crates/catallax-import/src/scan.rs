//! Open scans as reference data, RECIPE.md section 4.19.1.
//!
//! [`ScanReference::read`] reads the same binary STL, ASCII STL and Wavefront
//! OBJ bytes as closed mesh import, but asks nothing of them beyond being
//! readable: boundary edges (holes), inconsistent orientations, crossings,
//! degenerate triangles and isolated samples are all kept. The value is the
//! placed union of the file's closed triangles, segments and points, a
//! reference surface with no inside or outside. Nothing is closed, joined by
//! distance or inferred; the report counts what the file holds, holes included.
//!
//! Coordinates stay the file's exact rationals under exact scale and exact
//! rational placement. Distance queries for deviation reports run on the nearest
//! binary64 copies and are unproved.

use std::collections::BTreeMap;

use catallax_numbers::{Decider, Rational, Refusal, Sign};
use catallax_recipe::{Operation, Pose};

use crate::mesh::bvh::Bvh;
use crate::mesh::{fast_box_distance, fast_triangle_distance_sq, outward_box, parse};

/// The stated triangle limit for a scan, the same as closed mesh import's.
pub const MAX_SCAN_TRIANGLES: usize = crate::mesh::MAX_TRIANGLES;

/// The stated limit on vertices an OBJ scan may state, so that a file of
/// isolated samples is bounded too. A scan builds no B-rep, so its limit is
/// memory for exact coordinates, not face count.
pub const MAX_SCAN_VERTICES: usize = 1_000_000;

/// One element of the reference surface, by the file's own indices.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ScanItem {
    /// A triangle with positive area.
    Triangle([u32; 3]),
    /// A triangle of zero area whose corners are not all equal: its longest
    /// side, which contains the other corner.
    Segment([u32; 2]),
    /// A triangle whose corners are equal, or an isolated sample.
    Point(u32),
}

/// What the file holds, holes and all. Counts, never repairs.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct ScanReport {
    /// Distinct points, welded by exact equality only.
    pub vertices: usize,
    /// File triangles, degenerate ones included.
    pub triangles: usize,
    /// Triangles of zero area with distinct corners, kept as segments.
    pub degenerate_segments: usize,
    /// Triangles with three equal corners, kept as points.
    pub degenerate_points: usize,
    /// OBJ vertices no triangle uses, kept as samples.
    pub isolated_points: usize,
    /// Edges of positive-area triangles used by exactly one of them: the
    /// borders of holes and of the scanned region.
    pub boundary_edges: usize,
    /// Edges used by two triangles in the same direction.
    pub inconsistent_edges: usize,
    /// Edges used by three or more triangles.
    pub nonmanifold_edges: usize,
}

/// A placed scan: reference data, never a solid.
#[derive(Clone)]
pub struct ScanReference {
    sha256: String,
    points: Vec<[Rational; 3]>,
    triangles: Vec<[u32; 3]>,
    report: ScanReport,
    surface: FastSurface,
}

/// Triangles, segments and points in binary64 with a bounding-volume tree,
/// for unproved nearest-distance queries such as sampled deviation reports.
#[derive(Clone, Debug)]
pub struct FastSurface {
    points: Vec<[f64; 3]>,
    items: Vec<ScanItem>,
    tree: Bvh,
}

impl FastSurface {
    /// Indexes `items`, whose vertices are indices into `points`.
    pub fn new(points: Vec<[f64; 3]>, items: Vec<ScanItem>) -> Self {
        Self::new_with_check(points, items, |_| Ok::<(), std::convert::Infallible>(()))
            .unwrap_or_else(|never| match never {})
    }

    /// Builds the index while charging tree construction and checking the deadline.
    pub fn new_checked(
        points: Vec<[f64; 3]>,
        items: Vec<ScanItem>,
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        Self::new_with_check(points, items, |units| dec.charge_work(units))
    }

    fn new_with_check<E>(
        points: Vec<[f64; 3]>,
        items: Vec<ScanItem>,
        mut check: impl FnMut(u64) -> Result<(), E>,
    ) -> Result<Self, E> {
        check(items.len() as u64)?;
        let boxes = items
            .iter()
            .map(|item| {
                check(1)?;
                let at = |v: u32| &points[v as usize];
                Ok(match *item {
                    ScanItem::Triangle([a, b, c]) => outward_box([at(a), at(b), at(c)]),
                    ScanItem::Segment([a, b]) => outward_box([at(a), at(b), at(b)]),
                    ScanItem::Point(a) => outward_box([at(a), at(a), at(a)]),
                })
            })
            .collect::<Result<_, E>>()?;
        Ok(Self {
            tree: Bvh::new_checked(boxes, check)?,
            points,
            items,
        })
    }

    /// The vertices.
    pub fn points(&self) -> &[[f64; 3]] {
        &self.points
    }

    /// The elements.
    pub fn items(&self) -> &[ScanItem] {
        &self.items
    }

    /// The nearest element to `p` and its distance, or `None` when there are
    /// no elements. Of equally near elements the first one met is kept; the
    /// search order depends only on the input, so the answer is reproducible.
    pub fn nearest(&self, p: [f64; 3]) -> Option<(f64, usize)> {
        self.nearest_checked(p, |_| Ok::<(), std::convert::Infallible>(()))
            .unwrap_or_else(|never| match never {})
    }

    fn nearest_checked<E>(
        &self,
        p: [f64; 3],
        check: impl FnMut(u64) -> Result<(), E>,
    ) -> Result<Option<(f64, usize)>, E> {
        let mut best: Option<(f64, usize)> = None;
        self.tree.nearest_checked(
            |b| fast_box_distance(b, &p),
            |i| {
                let at = |v: u32| self.points[v as usize];
                let d = match self.items[i as usize] {
                    ScanItem::Triangle([a, b, c]) => {
                        fast_triangle_distance_sq(p, at(a), at(b), at(c))
                    }
                    ScanItem::Segment([a, b]) => fast_segment_distance_sq(p, at(a), at(b)),
                    ScanItem::Point(a) => fast_segment_distance_sq(p, at(a), at(a)),
                };
                if best.is_none_or(|(e, _)| d < e) {
                    best = Some((d, i as usize));
                }
                Some(d)
            },
            check,
        )?;
        Ok(best.map(|(d, i)| (d.sqrt(), i)))
    }

    /// Distance from `p` to the nearest element; infinite when there is none.
    pub fn distance(&self, p: [f64; 3]) -> f64 {
        self.nearest(p).map_or(f64::INFINITY, |(d, _)| d)
    }

    /// Unproved distance with work charged before tree nodes and candidates
    /// are examined. Refuses on an exhausted allowance or expired deadline.
    pub fn distance_checked(&self, p: [f64; 3], dec: &mut Decider) -> Result<f64, Refusal> {
        self.nearest_checked(p, |units| dec.charge_work(units))
            .map(|nearest| nearest.map_or(f64::INFINITY, |(distance, _)| distance))
    }
}

impl std::fmt::Debug for ScanReference {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("ScanReference")
            .field("sha256", &self.sha256)
            .field("report", &self.report)
            .finish()
    }
}

/// The rational orthogonal matrix of a quaternion of any nonzero length.
fn rotation(pose: &Pose) -> Result<[[Rational; 3]; 3], Refusal> {
    let [w, x, y, z] = pose.rotation.clone().map(|n| n.value().clone());
    let norm = &(&(&w * &w) + &(&x * &x)) + &(&(&y * &y) + &(&z * &z));
    if norm.sign() == Sign::Zero {
        return Err(Refusal::new(
            "shapes.zero-quaternion",
            "the scan pose's rotation quaternion is zero",
        ));
    }
    let two = Rational::from(2i64);
    let entry = |v: Rational| v.checked_div(&norm).expect("nonzero norm");
    let (ww, xx, yy, zz) = (&w * &w, &x * &x, &y * &y, &z * &z);
    let m = [
        [
            &(&ww + &xx) - &(&yy + &zz),
            &two * &(&(&x * &y) - &(&w * &z)),
            &two * &(&(&x * &z) + &(&w * &y)),
        ],
        [
            &two * &(&(&x * &y) + &(&w * &z)),
            &(&ww - &xx) + &(&yy - &zz),
            &two * &(&(&y * &z) - &(&w * &x)),
        ],
        [
            &two * &(&(&x * &z) - &(&w * &y)),
            &two * &(&(&y * &z) + &(&w * &x)),
            &(&ww - &xx) - &(&yy - &zz),
        ],
    ];
    Ok(m.map(|row| row.map(entry)))
}

fn sub(a: &[Rational; 3], b: &[Rational; 3]) -> [Rational; 3] {
    std::array::from_fn(|i| &a[i] - &b[i])
}

fn fast_segment_distance_sq(p: [f64; 3], a: [f64; 3], b: [f64; 3]) -> f64 {
    let ab: [f64; 3] = std::array::from_fn(|i| b[i] - a[i]);
    let ap: [f64; 3] = std::array::from_fn(|i| p[i] - a[i]);
    let length = ab.iter().map(|x| x * x).sum::<f64>();
    let t = if length > 0.0 {
        (ab.iter().zip(&ap).map(|(x, y)| x * y).sum::<f64>() / length).clamp(0.0, 1.0)
    } else {
        0.0
    };
    (0..3).map(|i| (ap[i] - ab[i] * t).powi(2)).sum()
}

impl ScanReference {
    /// Reads and places a scan under the stated limits. Work is charged to
    /// `dec`'s operation allowance and checked against its deadline.
    pub fn read(
        bytes: &[u8],
        scale: &Rational,
        pose: &Pose,
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        Self::read_with_limits(
            bytes,
            scale,
            pose,
            MAX_SCAN_TRIANGLES,
            MAX_SCAN_VERTICES,
            dec,
        )
    }

    /// [`ScanReference::read`] under smaller limits. Limits above the stated
    /// ones are lowered to them. A file over a limit refuses with
    /// `import.scan-too-large`.
    pub fn read_with_limits(
        bytes: &[u8],
        scale: &Rational,
        pose: &Pose,
        triangle_limit: usize,
        vertex_limit: usize,
        dec: &mut Decider,
    ) -> Result<Self, Refusal> {
        if scale.sign() != Sign::Positive {
            return Err(Refusal::new(
                "import.scan-scale",
                format!("scale {scale} is not positive; a scan import needs scale > 0"),
            ));
        }
        let matrix = rotation(pose)?;
        dec.charge_work(bytes.len() as u64 / 16 + 1)?;
        let soup = parse::read_with(
            bytes,
            triangle_limit.min(MAX_SCAN_TRIANGLES),
            vertex_limit.min(MAX_SCAN_VERTICES),
        )
        .map_err(|r| match r.code {
            "import.mesh-too-large" => Refusal::new("import.scan-too-large", r.detail),
            _ => r,
        })?;
        let corners = soup.points.len();
        let shift = pose.translation.0.clone().map(|n| n.value().clone());
        let mut points = Vec::with_capacity(corners + soup.isolated.len());
        for chunk in soup
            .points
            .into_iter()
            .chain(soup.isolated)
            .collect::<Vec<_>>()
            .chunks(1024)
        {
            dec.charge_work(chunk.len() as u64)?;
            for p in chunk {
                let scaled: [Rational; 3] = std::array::from_fn(|i| &p[i] * scale);
                points.push(std::array::from_fn(|r| {
                    (0..3).fold(shift[r].clone(), |s, c| &s + &(&matrix[r][c] * &scaled[c]))
                }));
            }
        }
        let approx: Vec<[f64; 3]> = points
            .iter()
            .map(|p| p.clone().map(|x: Rational| x.to_f64_nearest()))
            .collect();
        let mut report = ScanReport {
            vertices: points.len(),
            triangles: soup.triangles.len(),
            isolated_points: points.len() - corners,
            ..ScanReport::default()
        };
        let mut items = Vec::with_capacity(soup.triangles.len() + report.isolated_points);
        // Undirected edge -> (uses a->b with a < b, uses b->a).
        let mut edges: BTreeMap<(u32, u32), (u32, u32)> = BTreeMap::new();
        for chunk in soup.triangles.chunks(1024) {
            dec.charge_work(chunk.len() as u64 * 8)?;
            for &t in chunk {
                let [a, b, c] = t.map(|v| &points[v as usize]);
                let (u, v) = (sub(b, a), sub(c, a));
                let normal = [
                    &(&u[1] * &v[2]) - &(&u[2] * &v[1]),
                    &(&u[2] * &v[0]) - &(&u[0] * &v[2]),
                    &(&u[0] * &v[1]) - &(&u[1] * &v[0]),
                ];
                if normal.iter().any(|x| x.sign() != Sign::Zero) {
                    items.push(ScanItem::Triangle(t));
                    for j in 0..3 {
                        let (p, q) = (t[j], t[(j + 1) % 3]);
                        let entry = edges.entry((p.min(q), p.max(q))).or_default();
                        if p < q {
                            entry.0 += 1;
                        } else {
                            entry.1 += 1;
                        }
                    }
                } else if t[0] == t[1] && t[1] == t[2] {
                    report.degenerate_points += 1;
                    items.push(ScanItem::Point(t[0]));
                } else {
                    report.degenerate_segments += 1;
                    let length = |p: u32, q: u32| {
                        let d = sub(&points[p as usize], &points[q as usize]);
                        &(&(&d[0] * &d[0]) + &(&d[1] * &d[1])) + &(&d[2] * &d[2])
                    };
                    let sides = [(t[0], t[1]), (t[1], t[2]), (t[2], t[0])];
                    let mut longest = sides[0];
                    for side in &sides[1..] {
                        if length(side.0, side.1) > length(longest.0, longest.1) {
                            longest = *side;
                        }
                    }
                    items.push(ScanItem::Segment([longest.0, longest.1]));
                }
            }
        }
        for (forward, backward) in edges.into_values() {
            match forward + backward {
                1 => report.boundary_edges += 1,
                2 if forward != 1 => report.inconsistent_edges += 1,
                2 => {}
                _ => report.nonmanifold_edges += 1,
            }
        }
        items.extend((corners..points.len()).map(|v| ScanItem::Point(v as u32)));
        Ok(Self {
            sha256: crate::sha256_hex(bytes),
            points,
            triangles: soup.triangles,
            report,
            surface: FastSurface::new(approx, items),
        })
    }

    /// Reads the scan a recipe's `scan_import` step names from that file's
    /// bytes. Bytes whose SHA-256 differs from the step's are a recipe error,
    /// refused with `import.scan-hash` before anything is read.
    pub fn from_step(op: &Operation, bytes: &[u8], dec: &mut Decider) -> Result<Self, Refusal> {
        let Operation::ScanImport { file, pose, scale } = op else {
            return Err(Refusal::new(
                "import.scan-step",
                format!("a {} step is not a scan_import step", op.name()),
            ));
        };
        let actual = crate::sha256_hex(bytes);
        if actual != file.sha256 {
            return Err(Refusal::new(
                "import.scan-hash",
                format!(
                    "scan source {} sha256 mismatch: expected {}, actual {actual}",
                    file.name, file.sha256
                ),
            ));
        }
        Self::read(bytes, scale.value(), pose, dec)
    }

    /// SHA-256 of the source bytes, the scan's identity.
    pub fn sha256(&self) -> &str {
        &self.sha256
    }

    /// What the file holds.
    pub fn report(&self) -> &ScanReport {
        &self.report
    }

    /// Whether the reference holds no point at all.
    pub fn is_empty(&self) -> bool {
        self.surface.items.is_empty()
    }

    /// Placed exact points in millimetres: triangle corners in first-use
    /// order, then isolated samples in file order.
    pub fn points(&self) -> &[[Rational; 3]] {
        &self.points
    }

    /// Nearest binary64 copies of [`ScanReference::points`].
    pub fn approx_points(&self) -> &[[f64; 3]] {
        &self.surface.points
    }

    /// The file's triangles in file order, degenerate ones included.
    pub fn triangles(&self) -> &[[u32; 3]] {
        &self.triangles
    }

    /// The reference surface's elements: triangles in file order, with
    /// degenerate ones as segments or points, then isolated samples.
    pub fn items(&self) -> &[ScanItem] {
        &self.surface.items
    }

    /// The reference surface for unproved binary64 distance queries.
    pub fn surface(&self) -> &FastSurface {
        &self.surface
    }

    /// Unproved distance in millimetres from `p` to the reference surface,
    /// computed on binary64 copies. Infinite for an empty reference, as
    /// RECIPE.md section 4.19.1 defines distance to an empty set.
    pub fn fast_distance(&self, p: [f64; 3]) -> f64 {
        self.surface.distance(p)
    }
}

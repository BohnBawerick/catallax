//! The graphic callout of a surface texture, as AP242 polyline presentation.
//!
//! OpenCASCADE and several other readers have no class for the semantic
//! surface texture, so the writer also draws it: a leader to a point on the
//! face, the ISO 1302 symbol for the process condition, and the requirement,
//! such as `Ra 1.6`, in stroked text. The CAx-IF recommended practices for PMI
//! (version 4.1, sections 8.1 and 9.1) give the form: polylines in one plane
//! parallel to an annotation plane whose x axis is the reading direction.
//!
//! Every callout is drawn in the front plane: text runs along +x, upwards is
//! +z, and the reader looks along +y. The drawing is presentation only, so its
//! points are decimals; no exact decision is made here.

use catallax_numbers::{Integer, Rational};

use super::{Camera, Process, SurfaceTexture};
use crate::number::integer_sqrt;

/// A point of the drawing in file units, millimetres.
pub(crate) type Point = [Rational; 3];

/// The drawing of one callout.
pub(crate) struct Drawing {
    /// Each polyline of the leader, the symbol and the text.
    pub(crate) polylines: Vec<Vec<Point>>,
    /// The requirement as the callout writes it, such as `Ra 1.6`.
    pub(crate) text: String,
    junction: Point,
}

fn q(numerator: i64, denominator: i64) -> Rational {
    Rational::from_ratio(numerator, denominator).expect("a nonzero denominator")
}

/// tan 30 degrees to five places: the ISO 1302 legs stand at 60 degrees.
fn tan30() -> Rational {
    q(57_735, 100_000)
}

/// One grid step of the stroke font, so a capital is 3 mm high (ISO 3098).
fn unit() -> Rational {
    q(1, 2)
}

/// The stroke font: polylines on a grid 4 wide and 6 high, with the x-height
/// at 4 and descenders to -2. Polylines are separated by `;`.
fn glyph(c: char) -> Option<&'static str> {
    Some(match c {
        'A' => "0,0 2,6 4,0;1,3 3,3",
        'B' => "0,0 0,6 3,6 4,5 4,4 3,3 0,3;3,3 4,2 4,1 3,0 0,0",
        'C' => "4,1 3,0 1,0 0,1 0,5 1,6 3,6 4,5",
        'D' => "0,0 0,6 3,6 4,5 4,1 3,0 0,0",
        'E' => "4,0 0,0 0,6 4,6;0,3 3,3",
        'F' => "0,0 0,6 4,6;0,3 3,3",
        'G' => "4,5 3,6 1,6 0,5 0,1 1,0 3,0 4,1 4,3 2,3",
        'H' => "0,0 0,6;4,0 4,6;0,3 4,3",
        'I' => "1,0 3,0;2,0 2,6;1,6 3,6",
        'J' => "0,1 1,0 2,0 3,1 3,6;2,6 4,6",
        'K' => "0,0 0,6;4,6 0,2;1,3 4,0",
        'L' => "0,6 0,0 4,0",
        'M' => "0,0 0,6 2,3 4,6 4,0",
        'N' => "0,0 0,6 4,0 4,6",
        'O' => "1,0 0,1 0,5 1,6 3,6 4,5 4,1 3,0 1,0",
        // A slashed zero, so it differs from O.
        '0' => "1,0 0,1 0,5 1,6 3,6 4,5 4,1 3,0 1,0;1,1 3,5",
        'P' => "0,0 0,6 3,6 4,5 4,4 3,3 0,3",
        'Q' => "1,0 0,1 0,5 1,6 3,6 4,5 4,1 3,0 1,0;2,2 4,0",
        'R' => "0,0 0,6 3,6 4,5 4,4 3,3 0,3;2,3 4,0",
        'S' => "4,5 3,6 1,6 0,5 0,4 1,3 3,3 4,2 4,1 3,0 1,0 0,1",
        'T' => "0,6 4,6;2,6 2,0",
        'U' => "0,6 0,1 1,0 3,0 4,1 4,6",
        'V' => "0,6 2,0 4,6",
        'W' => "0,6 1,0 2,3 3,0 4,6",
        'X' => "0,0 4,6;0,6 4,0",
        'Y' => "0,6 2,3 4,6;2,3 2,0",
        'Z' => "0,6 4,6 0,0 4,0",
        'a' => "4,4 4,0;4,3 3,4 1,4 0,3 0,1 1,0 3,0 4,1",
        'b' => "0,6 0,0;0,1 1,0 3,0 4,1 4,3 3,4 1,4 0,3",
        'c' => "4,3 3,4 1,4 0,3 0,1 1,0 3,0 4,1",
        'd' => "4,6 4,0;4,3 3,4 1,4 0,3 0,1 1,0 3,0 4,1",
        'e' => "0,2 4,2 4,3 3,4 1,4 0,3 0,1 1,0 3,0 4,1",
        'f' => "1,0 1,5 2,6 3,6;0,4 3,4",
        'g' => "4,4 4,-1 3,-2 1,-2 0,-1;4,3 3,4 1,4 0,3 0,1 1,0 3,0 4,1",
        'h' => "0,6 0,0;0,3 1,4 3,4 4,3 4,0",
        'i' => "2,0 2,4;2,5 2,6",
        'j' => "3,4 3,-1 2,-2 1,-2;3,5 3,6",
        'k' => "0,6 0,0;3,4 0,1;1,2 3,0",
        'l' => "1,6 2,6 2,1 3,0",
        'm' => "0,0 0,4;0,3 1,4 2,3 2,0;2,3 3,4 4,3 4,0",
        'n' => "0,0 0,4;0,3 1,4 3,4 4,3 4,0",
        'o' => "1,0 0,1 0,3 1,4 3,4 4,3 4,1 3,0 1,0",
        'p' => "0,4 0,-2;0,3 1,4 3,4 4,3 4,1 3,0 1,0 0,1",
        'q' => "4,4 4,-2;4,3 3,4 1,4 0,3 0,1 1,0 3,0 4,1",
        'r' => "0,0 0,4;0,3 1,4 3,4",
        's' => "4,3 3,4 1,4 0,3 1,2 3,2 4,1 3,0 1,0 0,1",
        't' => "1,6 1,1 2,0 3,0;0,4 3,4",
        'u' => "0,4 0,1 1,0 3,0 4,1;4,4 4,0",
        'v' => "0,4 2,0 4,4",
        'w' => "0,4 1,0 2,2 3,0 4,4",
        'x' => "0,0 4,4;0,4 4,0",
        'y' => "0,4 2,0;4,4 1,-2",
        'z' => "0,4 4,4 0,0 4,0",
        '1' => "1,5 2,6 2,0;1,0 3,0",
        '2' => "0,5 1,6 3,6 4,5 4,4 0,0 4,0",
        '3' => "0,5 1,6 3,6 4,5 4,4 3,3 4,2 4,1 3,0 1,0 0,1;1,3 3,3",
        '4' => "3,0 3,6 0,2 4,2",
        '5' => "4,6 0,6 0,3 3,3 4,2 4,1 3,0 1,0 0,1",
        '6' => "4,5 3,6 1,6 0,5 0,1 1,0 3,0 4,1 4,2 3,3 0,3",
        '7' => "0,6 4,6 1,0",
        '8' => "1,3 0,4 0,5 1,6 3,6 4,5 4,4 3,3 1,3 0,2 0,1 1,0 3,0 4,1 4,2 3,3",
        '9' => "0,1 1,0 3,0 4,1 4,5 3,6 1,6 0,5 0,4 1,3 4,3",
        '.' => "1,0 1,1",
        '\u{2300}' => "1,1 0,2 0,4 1,5 3,5 4,4 4,2 3,1 1,1;0,0 4,6",
        '(' => "2,-1 1,1 1,4 2,6",
        ')' => "2,-1 3,1 3,4 2,6",
        '|' => "2,-1 2,6",
        '+' => "0,3 4,3;2,1 2,5",
        '-' => "0,3 4,3",
        ' ' => "",
        _ => return None,
    })
}

/// How far the pen moves after a character, in grid steps.
fn advance(c: char) -> i64 {
    match c {
        '.' => 3,
        ' ' => 4,
        _ => 6,
    }
}

/// The polylines of `text` with its baseline starting at `(u, v)` of the
/// plane, and the end of the last character. Every character is one the
/// font draws: the writer validated the parameter as letters and digits and
/// writes the value as a decimal.
fn stroke(text: &str, u: &Rational, v: &Rational) -> (Vec<Vec<[Rational; 2]>>, Rational) {
    let mut lines = Vec::new();
    let mut pen = u.clone();
    for c in text.chars() {
        let strokes = glyph(c).expect("the callout draws letters, digits and decimals");
        for line in strokes.split(';').filter(|line| !line.is_empty()) {
            let points = line
                .split(' ')
                .map(|point| {
                    let (x, y) = point.split_once(',').expect("a glyph point is x,y");
                    let at = |text: &str| {
                        Rational::from_i64(text.parse::<i64>().expect("a glyph coordinate"))
                            * unit()
                    };
                    [&pen + &at(x), v + &at(y)]
                })
                .collect();
            lines.push(points);
        }
        pen = &pen + &(Rational::from_i64(advance(c)) * unit());
    }
    // The last character's trailing space is not part of the text.
    let end = &pen - &(Rational::from_i64(2) * unit());
    (lines, end)
}

/// The process condition in the PMI Unicode string, CAx-IF recommended
/// practice version 1.0, section 7.4.
fn process_code(process: Process) -> &'static str {
    match process {
        Process::Any => "APA",
        Process::Removal => "MRR",
        Process::NoRemoval => "NMR",
    }
}

/// The equivalent Unicode string of a surface texture: the keyword with its
/// process and `N` (not all around), and the requirement in zone 4, as in
/// `SRFMRRN\w\w\w\wRa 1.6` (PMI Unicode string practice, sections 7.1 to
/// 7.5).
pub(crate) fn equivalent_string(process: Process, requirement: &str) -> String {
    format!("SRF{}N\\w\\w\\w\\w{requirement}", process_code(process))
}

/// Draws the callout of `texture` with the requirement `value` written as
/// a decimal. `anchor` is a point on the texture's face, where the leader ends;
/// `slot` vertically separates callouts in their stable item order.
pub(crate) fn draw<F>(
    texture: &SurfaceTexture<F>,
    value: &str,
    anchor: Option<&Point>,
    slot: usize,
) -> Drawing {
    let zero = Rational::zero();
    let text = format!("{} {value}", texture.parameter);
    // Plane coordinates (u along x, v along z) of the drawing, from the
    // symbol's point.
    let mut flat: Vec<Vec<[Rational; 2]>> = Vec::new();
    let origin = leader(anchor, slot, &mut flat);
    // ISO 1302 for a 3 mm text height: legs H1 = 4.2 and H2 = 9 high.
    let short = q(21, 5);
    let long = Rational::from_i64(9);
    let gap = q(3, 5);
    let left = [&zero - &(&short * &tan30()), short.clone()];
    let right = [&long * &tan30(), long.clone()];
    let start = &right[0] + &gap;
    let baseline = &(&long - &gap) - &(Rational::from_i64(6) * unit());
    let (letters, end) = stroke(&text, &start, &baseline);
    let line_end = [&end + &gap, long.clone()];
    flat.push(vec![
        left.clone(),
        [zero.clone(), zero.clone()],
        right.clone(),
        line_end,
    ]);
    match texture.process {
        Process::Any => {}
        // A bar across the short leg closes the symbol.
        Process::Removal => flat.push(vec![left, [&short * &tan30(), short.clone()]]),
        // A circle in the V, touching both legs, its top at H1.
        Process::NoRemoval => {
            let centre = &short * &q(2, 3);
            let radius = &short * &q(1, 3);
            let (c, s) = (q(86_603, 100_000), q(1, 2));
            let one = Rational::one();
            let turn = [
                (one.clone(), zero.clone()),
                (c.clone(), s.clone()),
                (s.clone(), c.clone()),
                (zero.clone(), one.clone()),
            ];
            let mut circle = Vec::new();
            for quarter in 0..4 {
                for (index, (x, y)) in turn.iter().enumerate() {
                    if (quarter > 0 && index == 0) || (quarter == 3 && index == turn.len() - 1) {
                        continue;
                    }
                    let (x, y) = match quarter {
                        0 => (x.clone(), y.clone()),
                        1 => (&zero - y, x.clone()),
                        2 => (&zero - x, &zero - y),
                        _ => (y.clone(), &zero - x),
                    };
                    circle.push([&x * &radius, &centre + &(&y * &radius)]);
                }
            }
            circle.push(circle[0].clone());
            flat.push(circle);
        }
    }
    flat.extend(letters);
    Drawing {
        polylines: place(&origin, flat),
        text,
        junction: origin,
    }
}

/// The leader of a callout in slot `slot`, up and to the right of `anchor`
/// with an open arrow on the face, added to `flat` in plane coordinates from
/// the returned start of the callout's symbol or text.
fn leader(anchor: Option<&Point>, slot: usize, flat: &mut Vec<Vec<[Rational; 2]>>) -> Point {
    let zero = Rational::zero();
    let offset = Rational::from_i64(i64::try_from(slot).expect("callout slot fits i64") * 14);
    match anchor {
        Some(anchor) => {
            let reach = Rational::from_i64(7);
            let symbol = [
                &anchor[0] + &reach,
                anchor[1].clone(),
                &(&anchor[2] + &reach) + &offset,
            ];
            let back = |u: Rational, v: Rational| {
                [&zero - &reach + u, &zero - &reach - offset.clone() + v]
            };
            let tip = back(zero.clone(), zero.clone());
            // Arrow wings 2.5 mm long at 15 degrees either side of the leader.
            let wings = [
                back(q(216_506, 100_000), q(5, 4)),
                back(q(5, 4), q(216_506, 100_000)),
            ];
            flat.push(vec![tip.clone(), [zero.clone(), zero.clone()]]);
            flat.push(vec![wings[0].clone(), tip, wings[1].clone()]);
            symbol
        }
        None => [zero.clone(), zero, offset],
    }
}

fn place(origin: &Point, flat: Vec<Vec<[Rational; 2]>>) -> Vec<Vec<Point>> {
    flat.into_iter()
        .map(|line| {
            line.into_iter()
                .map(|[u, v]| [&origin[0] + &u, origin[1].clone(), &origin[2] + &v])
                .collect()
        })
        .collect()
}

/// Whether the stroke font draws every character of `text`.
pub(crate) fn drawable(text: &str) -> bool {
    text.chars().all(|c| glyph(c).is_some())
}

pub(crate) fn label(line: &str) -> Drawing {
    let zero = Rational::zero();
    let gap = Rational::one();
    let height = Rational::from_i64(5);
    let (letters, end) = stroke(line, &gap, &gap);
    let right = &end + &gap;
    let mut flat = vec![vec![
        [zero.clone(), zero.clone()],
        [right.clone(), zero.clone()],
        [right, height.clone()],
        [zero.clone(), height],
        [zero.clone(), zero.clone()],
    ]];
    flat.extend(letters);
    Drawing {
        polylines: place(&[zero.clone(), zero.clone(), zero], flat),
        text: String::new(),
        junction: std::array::from_fn(|_| Rational::zero()),
    }
}

/// A coordinate or measure of the drawing to six decimal places, a
/// micrometre's thousandth.
pub(crate) fn presented(value: &Rational) -> Rational {
    let scale = Integer::from_i64(1_000_000);
    let shifted = &(value * &Rational::from_integer(scale.clone())) + &q(1, 2);
    Rational::new(shifted.floor(), scale).expect("a nonzero denominator")
}

/// The square root of a nonnegative rational to nine decimal places, below.
fn root(value: &Rational) -> Rational {
    let scale = Integer::from_i64(10).pow(9);
    let squared = Rational::from_integer(&scale * &scale);
    let whole = (value * &squared).floor();
    Rational::new(integer_sqrt(&whole), scale).expect("a nonzero denominator")
}

fn unit_axis(axis: &Point) -> Point {
    let axis = super::decimal_direction(axis);
    let largest = axis
        .iter()
        .map(Rational::abs)
        .max()
        .expect("three coordinates");
    let scaled = axis.each_ref().map(|value| value / &largest);
    let squared = dot(&scaled, &scaled);
    let scale = crate::number::rational(&(&Rational::one() / &(&largest * &root(&squared))), 17);
    axis.each_ref().map(|value| value * &scale.value)
}

fn dot(a: &Point, b: &Point) -> Rational {
    a.iter()
        .zip(b)
        .fold(Rational::zero(), |sum, (a, b)| &sum + &(a * b))
}

pub(crate) fn cross(a: &Point, b: &Point) -> Point {
    std::array::from_fn(|i| {
        let j = (i + 1) % 3;
        let k = (i + 2) % 3;
        &(&a[j] * &b[k]) - &(&a[k] * &b[j])
    })
}

impl Drawing {
    pub(crate) fn in_view(
        &mut self,
        anchor: &Point,
        camera: &Camera,
        slot: usize,
        count: usize,
    ) -> Result<[Point; 2], &'static str> {
        let across = cross(&camera.direction, &camera.up);
        let right = unit_axis(&across);
        let up = unit_axis(&camera.up);
        let offset = std::array::from_fn(|i| &anchor[i] - &camera.eye[i]);
        let depth = dot(&offset, &unit_axis(&camera.direction));
        if depth <= Rational::zero() {
            return Err("the callout is behind the camera");
        }
        let scale = if camera.perspective {
            &depth / &camera.distance
        } else {
            Rational::one()
        };
        let [width, height] = &camera.window;
        let rows = Rational::from_i64(i64::try_from(count).expect("view items fit i64"));
        let cell = &(height * &q(9, 10)) / &rows;
        let row = Rational::from_i64(i64::try_from(slot + 1).expect("view slot fits i64"));
        let left = width * &q(-9, 20);
        let bottom = &(height * &q(9, 20)) - &(&row * &cell);
        let bounds: [(Rational, Rational); 2] = [0, 2].map(|axis| {
            let mut values = self.polylines.iter().flatten().map(|point| &point[axis]);
            let first = values.next().expect("a drawing has strokes");
            let (lo, hi) = values.fold((first, first), |(lo, hi), v| (lo.min(v), hi.max(v)));
            (lo.clone(), hi - lo)
        });
        let fit = Rational::one()
            .min(&(width * &q(9, 10)) / &bounds[0].1)
            .min(&(&cell * &q(4, 5)) / &bounds[1].1);
        let at = |x: &Rational, y: &Rational| {
            let x = &(&(&left + &(&(x - &bounds[0].0) * &fit)) * &scale) - &dot(&offset, &right);
            let y = &(&(&bottom + &(&(y - &bounds[1].0) * &fit)) * &scale) - &dot(&offset, &up);
            let x = crate::number::rational(&x, 17).value;
            let y = crate::number::rational(&y, 17).value;
            std::array::from_fn(|i| &anchor[i] + &(&(&x * &right[i]) + &(&y * &up[i])))
        };
        for point in self.polylines.iter_mut().flatten() {
            *point = at(&point[0], &point[2]);
        }
        if self
            .polylines
            .iter()
            .any(|line| line.windows(2).all(|pair| pair[0] == pair[1]))
        {
            return Err("a callout stroke collapses in the saved view");
        }
        let end = at(&self.junction[0], &self.junction[2]);
        if *anchor != end {
            self.polylines.insert(0, vec![anchor.clone(), end]);
        }
        Ok([across, camera.direction.clone()])
    }

    /// The total length of the polylines and their centroid weighted by
    /// length, the `polyline curve length` and `polyline centre point`
    /// validation properties (recommended practices, section 10.2.3.1).
    pub(crate) fn length_and_centre(&self) -> (Rational, Point) {
        let mut length = Rational::zero();
        let mut moment = [Rational::zero(), Rational::zero(), Rational::zero()];
        for line in &self.polylines {
            for pair in line.windows(2) {
                let [a, b] = pair else {
                    unreachable!("windows of two")
                };
                let squared = (0..3).fold(Rational::zero(), |sum, i| {
                    let d = &b[i] - &a[i];
                    &sum + &(&d * &d)
                });
                let piece = root(&squared);
                for i in 0..3 {
                    let middle = (&a[i] + &b[i]).mul_pow2(-1);
                    moment[i] = &moment[i] + &(&middle * &piece);
                }
                length = &length + &piece;
            }
        }
        let centre = moment.map(|m| if length.is_zero() { m } else { &m / &length });
        (length, centre)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn texture() -> SurfaceTexture<()> {
        SurfaceTexture {
            parameter: "Ra".to_owned(),
            value: Rational::from_i64(16),
            process: Process::NoRemoval,
            feature: Vec::new(),
        }
    }

    #[test]
    fn no_removal_circle_has_no_repeated_adjacent_points() {
        let drawing = draw(&texture(), "1.6", None, 0);
        let circle = drawing
            .polylines
            .iter()
            .find(|line| line.first() == line.last())
            .unwrap();
        assert_eq!(circle.first(), circle.last());
        assert!(
            circle.windows(2).all(|pair| pair[0] != pair[1]),
            "circle contains a zero-length segment"
        );
    }

    #[test]
    fn anchored_callout_slots_separate_same_face_textures() {
        let anchor = [
            Rational::from_i64(3),
            Rational::zero(),
            Rational::from_i64(5),
        ];
        let first = draw(&texture(), "1.6", Some(&anchor), 0);
        let mut second_texture = texture();
        second_texture.parameter = "Rz".to_owned();
        let second = draw(&second_texture, "6.3", Some(&anchor), 1);
        let second_again = draw(&second_texture, "6.3", Some(&anchor), 1);

        assert_eq!(first.polylines[0][0], anchor);
        assert_eq!(second.polylines[0][0], anchor);
        assert_ne!(first.polylines[2], second.polylines[2]);
        assert_eq!(second.polylines, second_again.polylines);
    }
}

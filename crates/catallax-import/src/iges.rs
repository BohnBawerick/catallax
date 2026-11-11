//! Recognition of IGES 5.3 manifold solid B-rep objects as exact source data.
//!
//! An IGES file becomes the same [`StepFile`] that a STEP file does, so the
//! builder, its exact incidence checks and its refusals are shared. Entity
//! names are directory-entry sequence numbers, such as `"D17"`. Decimal text
//! becomes an exact rational, lengths are converted to millimetres exactly, and
//! anything outside the subset is refused with its entity and IGES type.
//!
//! Read: 186 manifold solid, 514 closed shell, 510 face, 508 loop, 504 edge list,
//! 502 vertex list; 190 plane, 192 cylinder, 194 cone, 196 sphere, 198 torus and
//! 128 rational B-spline surfaces; 110 line, 100 circular arc and 126 rational
//! B-spline curves; 116 points, 123 directions and 124 orthonormal transforms.
//! Parameter-space curves of loops are not read: the model-space curve decides,
//! as it does for STEP.

use std::collections::BTreeMap;

use catallax_numbers::{Decider, Integer, Rational};
use catallax_shapes::spline_tools::{RationalSplineCurve, RationalSplineSurface};

use crate::{
    Curve, Face, FaceLoop, ImportRefusal, OrientedEdge, Placement, SemiAngle, Solid, StepFile,
    StepSplineCurve, StepSplineSurface, Surface, Vec3, Vertex,
};

/// One parameter of an entity or of the global section.
#[derive(Clone, Debug)]
enum Token {
    Empty,
    Text(String),
    Number(String),
}

struct Entry {
    kind: u32,
    form: u32,
    transform: usize,
    params: Vec<Token>,
}

/// Longest chain of transformation matrices followed before refusing.
const MAX_TRANSFORM_DEPTH: usize = 8;

/// Reads and recognizes an IGES 5.3 file.
pub fn read_iges(bytes: &[u8]) -> Result<StepFile, ImportRefusal> {
    read_iges_with_decider(bytes, &Decider::default())
}

/// Reads exact IGES source geometry using caller-owned cancellation and work
/// policy. Spline data do not acquire incidence or solid proofs through reading.
pub fn read_iges_with_decider(bytes: &[u8], dec: &Decider) -> Result<StepFile, ImportRefusal> {
    let text = std::str::from_utf8(bytes).map_err(|_| malformed("IGES text is not ASCII"))?;
    let mut sections: BTreeMap<char, Vec<&str>> = BTreeMap::new();
    for line in text.lines() {
        let line = line.trim_end_matches('\r');
        if line.trim().is_empty() {
            continue;
        }
        let section = line
            .chars()
            .nth(72)
            .ok_or_else(|| malformed("IGES record shorter than 73 columns"))?;
        if !"SGDPTC".contains(section) {
            return Err(malformed(format!(
                "unknown IGES section letter `{section}`"
            )));
        }
        if section == 'C' {
            return Err(malformed("compressed IGES files are not read"));
        }
        sections.entry(section).or_default().push(line);
    }
    let field = |line: &str, i: usize| line.get(i * 8..i * 8 + 8).unwrap_or("").trim().to_string();
    let global: String = sections
        .get(&'G')
        .ok_or_else(|| malformed("IGES file has no global section"))?
        .iter()
        .map(|l| &l[..72])
        .collect();
    let (pd, rd) = delimiters(&global);
    let global = tokens(&global, pd, rd)?;
    let scale = units(&global)?;
    let directory = sections.get(&'D').map_or(&[][..], Vec::as_slice);
    if !directory.len().is_multiple_of(2) {
        return Err(malformed("IGES directory entries come in pairs of records"));
    }
    let parameters = sections.get(&'P').map_or(&[][..], Vec::as_slice);
    let mut entries = BTreeMap::new();
    for (index, pair) in directory.chunks(2).enumerate() {
        let int = |line: &str, i: usize| -> Result<usize, ImportRefusal> {
            let text = field(line, i);
            if text.is_empty() {
                return Ok(0);
            }
            text.parse()
                .map_err(|_| malformed(format!("unreadable IGES directory field `{text}`")))
        };
        let sequence = 2 * index + 1;
        let kind = int(pair[0], 0)?;
        let start = int(pair[0], 1)?;
        let transform = int(pair[0], 6)?;
        let lines = int(pair[1], 3)?;
        let form = int(pair[1], 4)?;
        if start == 0 || start + lines - 1 > parameters.len() || lines == 0 {
            return Err(malformed(format!(
                "D{sequence} points outside the parameter section"
            )));
        }
        let data: String = parameters[start - 1..start - 1 + lines]
            .iter()
            .map(|l| &l[..64.min(l.len())])
            .collect();
        let mut params = tokens(&data, pd, rd)?;
        if params.is_empty() || number_text(&params[0]) != Some(kind.to_string().as_str()) {
            return Err(malformed(format!(
                "D{sequence} parameter data does not start with its type {kind}"
            )));
        }
        params.remove(0);
        entries.insert(
            sequence,
            Entry {
                kind: u32::try_from(kind).map_err(|_| malformed("IGES type out of range"))?,
                form: u32::try_from(form).map_err(|_| malformed("IGES form out of range"))?,
                transform,
                params,
            },
        );
    }
    Reader {
        entries: &entries,
        scale,
        dec,
    }
    .file()
}

fn malformed(detail: impl Into<String>) -> ImportRefusal {
    ImportRefusal::Malformed {
        detail: detail.into(),
    }
}

fn unsupported(sequence: usize, what: impl Into<String>) -> ImportRefusal {
    ImportRefusal::Unsupported {
        entity: format!("D{sequence}"),
        type_name: what.into(),
    }
}

fn type_name(kind: u32) -> &'static str {
    match kind {
        100 => "CIRCULAR ARC",
        102 => "COMPOSITE CURVE",
        104 => "CONIC ARC",
        106 => "COPIOUS DATA",
        108 => "PLANE",
        110 => "LINE",
        112 => "PARAMETRIC SPLINE CURVE",
        114 => "PARAMETRIC SPLINE SURFACE",
        116 => "POINT",
        118 => "RULED SURFACE",
        120 => "SURFACE OF REVOLUTION",
        122 => "TABULATED CYLINDER",
        123 => "DIRECTION",
        124 => "TRANSFORMATION MATRIX",
        126 => "RATIONAL B-SPLINE CURVE",
        128 => "RATIONAL B-SPLINE SURFACE",
        130 => "OFFSET CURVE",
        140 => "OFFSET SURFACE",
        142 => "CURVE ON A PARAMETRIC SURFACE",
        143 => "BOUNDED SURFACE",
        144 => "TRIMMED SURFACE",
        186 => "MANIFOLD SOLID B-REP OBJECT",
        190 => "PLANE SURFACE",
        192 => "RIGHT CIRCULAR CYLINDRICAL SURFACE",
        194 => "RIGHT CIRCULAR CONICAL SURFACE",
        196 => "SPHERICAL SURFACE",
        198 => "TOROIDAL SURFACE",
        502 => "VERTEX",
        504 => "EDGE",
        508 => "LOOP",
        510 => "FACE",
        514 => "SHELL",
        _ => "ENTITY",
    }
}

fn kind_text(kind: u32, form: u32) -> String {
    format!("IGES {kind} form {form} ({})", type_name(kind))
}

/// The parameter and record delimiters: the first two global parameters,
/// each a one-character Hollerith string, defaulting to `,` and `;`.
fn delimiters(global: &str) -> (char, char) {
    // A one-character Hollerith string `1Hx` at the start of `text`.
    let hollerith = |text: &str| text.strip_prefix("1H").and_then(|t| t.chars().next());
    let (pd, rest) = match hollerith(global) {
        Some(pd) => (pd, &global[2 + pd.len_utf8()..]),
        None => (',', global),
    };
    let rest = rest.strip_prefix(pd).unwrap_or(rest);
    (pd, hollerith(rest).unwrap_or(';'))
}

/// Splits parameters at `pd` up to the first `rd`, reading Hollerith strings.
fn tokens(text: &str, pd: char, rd: char) -> Result<Vec<Token>, ImportRefusal> {
    let chars: Vec<char> = text.chars().collect();
    let mut out = Vec::new();
    let mut i = 0;
    loop {
        while i < chars.len() && chars[i] == ' ' {
            i += 1;
        }
        let start = i;
        while i < chars.len() && chars[i].is_ascii_digit() {
            i += 1;
        }
        let token = if i > start && i < chars.len() && chars[i] == 'H' {
            let count: usize = chars[start..i]
                .iter()
                .collect::<String>()
                .parse()
                .map_err(|_| malformed("unreadable Hollerith count"))?;
            let body = chars
                .get(i + 1..i + 1 + count)
                .ok_or_else(|| malformed("Hollerith string runs past its record"))?;
            i += 1 + count;
            Token::Text(body.iter().collect())
        } else {
            while i < chars.len() && chars[i] != pd && chars[i] != rd {
                i += 1;
            }
            let value: String = chars[start..i].iter().collect();
            let value = value.trim();
            if value.is_empty() {
                Token::Empty
            } else {
                Token::Number(value.to_string())
            }
        };
        while i < chars.len() && chars[i] == ' ' {
            i += 1;
        }
        out.push(token);
        match chars.get(i) {
            Some(c) if *c == pd => i += 1,
            Some(c) if *c == rd => return Ok(out),
            Some(c) => return Err(malformed(format!("unexpected `{c}` in IGES parameters"))),
            None => return Err(malformed("IGES parameters have no record delimiter")),
        }
    }
}

fn number_text(token: &Token) -> Option<&str> {
    match token {
        Token::Number(text) => Some(text),
        _ => None,
    }
}

/// An exact decimal: IGES reals may use a `D` exponent and omit digits around
/// the decimal point.
fn real(token: &Token) -> Result<Rational, ImportRefusal> {
    let text = number_text(token).ok_or_else(|| malformed("expected an IGES number"))?;
    let mut fixed = text.replace(['D', 'd'], "E");
    let (sign, body) = match fixed.strip_prefix(['+', '-']) {
        Some(body) => (&text[..1], body.to_string()),
        None => ("", fixed.clone()),
    };
    let body = if body.starts_with('.') {
        format!("0{body}")
    } else {
        body
    };
    let body = body.replace(".E", ".0E");
    let body = if body.ends_with('.') {
        format!("{body}0")
    } else {
        body
    };
    fixed = format!("{}{body}", if sign == "-" { "-" } else { "" });
    Rational::parse(&fixed).map_err(|_| malformed(format!("unreadable IGES number `{text}`")))
}

fn integer(token: &Token) -> Result<usize, ImportRefusal> {
    let text = number_text(token).ok_or_else(|| malformed("expected an IGES integer"))?;
    text.parse()
        .map_err(|_| malformed(format!("unreadable IGES integer `{text}`")))
}

/// Millimetres per model unit from global parameters 13 to 15.
fn units(global: &[Token]) -> Result<Rational, ImportRefusal> {
    let get = |i: usize| global.get(i).cloned().unwrap_or(Token::Empty);
    match get(12) {
        Token::Empty => {}
        token if real(&token)? == Rational::one() => {}
        _ => {
            return Err(ImportRefusal::UnknownUnit {
                found: "model space scale other than 1".into(),
            });
        }
    }
    let flag = match get(13) {
        Token::Empty => 1,
        token => integer(&token)?,
    };
    let name = match get(14) {
        Token::Text(name) => name.to_ascii_uppercase(),
        _ => String::new(),
    };
    let ratio = |n: i64, d: i64| Rational::from_ratio(n, d).expect("nonzero unit denominator");
    let by_name = |name: &str| match name {
        "IN" | "INCH" => Some(ratio(254, 10)),
        "MM" => Some(Rational::one()),
        "FT" => Some(ratio(3048, 10)),
        "MI" => Some(ratio(1_609_344, 1)),
        "M" => Some(ratio(1000, 1)),
        "KM" => Some(ratio(1_000_000, 1)),
        "MIL" => Some(ratio(254, 10_000)),
        "UM" => Some(ratio(1, 1000)),
        "CM" => Some(ratio(10, 1)),
        "UIN" => Some(ratio(254, 10_000_000)),
        _ => None,
    };
    let unit = match flag {
        1 => Some(ratio(254, 10)),
        2 => Some(Rational::one()),
        3 => by_name(&name),
        4 => Some(ratio(3048, 10)),
        5 => Some(ratio(1_609_344, 1)),
        6 => Some(ratio(1000, 1)),
        7 => Some(ratio(1_000_000, 1)),
        8 => Some(ratio(254, 10_000)),
        9 => Some(ratio(1, 1000)),
        10 => Some(ratio(10, 1)),
        11 => Some(ratio(254, 10_000_000)),
        _ => None,
    };
    unit.ok_or_else(|| ImportRefusal::UnknownUnit {
        found: format!("IGES unit flag {flag} `{name}`"),
    })
}

/// An exact rigid motion `x -> R x + t` with orthonormal `R`.
#[derive(Clone)]
struct Motion {
    rows: [[Rational; 3]; 3],
    shift: Vec3,
}

impl Motion {
    fn identity() -> Self {
        let unit = |i: usize, j: usize| {
            if i == j {
                Rational::one()
            } else {
                Rational::zero()
            }
        };
        Self {
            rows: std::array::from_fn(|i| std::array::from_fn(|j| unit(i, j))),
            shift: std::array::from_fn(|_| Rational::zero()),
        }
    }

    fn rotate(&self, v: &Vec3) -> Vec3 {
        std::array::from_fn(|i| {
            (0..3).fold(Rational::zero(), |sum, j| sum + &self.rows[i][j] * &v[j])
        })
    }

    fn point(&self, p: &Vec3) -> Vec3 {
        let r = self.rotate(p);
        std::array::from_fn(|i| &r[i] + &self.shift[i])
    }

    /// `self` after `inner`.
    fn after(&self, inner: &Self) -> Self {
        let rows = std::array::from_fn(|i| {
            std::array::from_fn(|j| {
                (0..3).fold(Rational::zero(), |sum, k| {
                    sum + &self.rows[i][k] * &inner.rows[k][j]
                })
            })
        });
        Self {
            rows,
            shift: self.point(&inner.shift),
        }
    }
}

fn dot(a: &Vec3, b: &Vec3) -> Rational {
    (0..3).fold(Rational::zero(), |sum, i| sum + &a[i] * &b[i])
}

fn cross(a: &Vec3, b: &Vec3) -> Vec3 {
    [
        &a[1] * &b[2] - &a[2] * &b[1],
        &a[2] * &b[0] - &a[0] * &b[2],
        &a[0] * &b[1] - &a[1] * &b[0],
    ]
}

/// A placement with its reference direction projected off the axis, as STEP's
/// `first_proj_axis` does. Without a reference the STEP default is used.
fn placement(origin: Vec3, axis: Vec3, ref_dir: Option<Vec3>) -> Result<Placement, ImportRefusal> {
    if axis.iter().all(Rational::is_zero) {
        return Err(malformed("zero IGES axis direction"));
    }
    let ref_dir = ref_dir.unwrap_or_else(|| {
        if axis[1].is_zero() && axis[2].is_zero() {
            [Rational::zero(), Rational::one(), Rational::zero()]
        } else {
            [Rational::one(), Rational::zero(), Rational::zero()]
        }
    });
    let along = &dot(&axis, &ref_dir) / &dot(&axis, &axis);
    let ref_dir: Vec3 = std::array::from_fn(|i| &ref_dir[i] - &(&along * &axis[i]));
    if ref_dir.iter().all(Rational::is_zero) {
        return Err(malformed("IGES axis and reference direction are parallel"));
    }
    Ok(Placement {
        origin,
        axis,
        ref_dir,
    })
}

/// The floor of the square root of a nonnegative integer, by Newton's method
/// from above.
fn isqrt(n: &Integer) -> Integer {
    if n.is_zero() {
        return Integer::zero();
    }
    let two = Integer::from_i64(2);
    let mut x = Integer::one().mul_pow2(n.bit_length().div_ceil(2));
    loop {
        let (q, _) = n.div_rem_floor(&x).expect("nonzero divisor");
        let (y, _) = (&x + &q).div_rem_floor(&two).expect("nonzero divisor");
        if y >= x {
            return x;
        }
        x = y;
    }
}

/// The exact square root of a nonnegative rational, if it has one.
fn exact_sqrt(q: &Rational) -> Option<Rational> {
    if q.sign() == catallax_numbers::Sign::Negative {
        return None;
    }
    let (n, d) = (isqrt(q.numerator()), isqrt(q.denominator()));
    (&(&n * &n) == q.numerator() && &(&d * &d) == q.denominator())
        .then(|| Rational::new(n, d))
        .flatten()
}

fn source_error(error: catallax_numbers::Refusal) -> ImportRefusal {
    if error.code.starts_with("numbers.") {
        ImportRefusal::ResourceLimit
    } else {
        malformed(error.detail)
    }
}

struct Reader<'a> {
    entries: &'a BTreeMap<usize, Entry>,
    scale: Rational,
    dec: &'a Decider,
}

impl Reader<'_> {
    fn entry(&self, sequence: usize) -> Result<&Entry, ImportRefusal> {
        self.dec
            .check_deadline()
            .map_err(|_| ImportRefusal::ResourceLimit)?;
        self.entries
            .get(&sequence)
            .ok_or_else(|| malformed(format!("reference to missing directory entry D{sequence}")))
    }

    /// The parameters of an entry of one of `kinds`, refused otherwise.
    fn params(&self, sequence: usize, kinds: &[u32]) -> Result<&[Token], ImportRefusal> {
        let e = self.entry(sequence)?;
        if !kinds.contains(&e.kind) {
            return Err(unsupported(sequence, kind_text(e.kind, e.form)));
        }
        Ok(&e.params)
    }

    fn param<'p>(&self, params: &'p [Token], i: usize) -> Result<&'p Token, ImportRefusal> {
        params
            .get(i)
            .ok_or_else(|| malformed(format!("IGES entity has no parameter {}", i + 1)))
    }

    fn int(&self, params: &[Token], i: usize) -> Result<usize, ImportRefusal> {
        integer(self.param(params, i)?)
    }

    fn real(&self, params: &[Token], i: usize) -> Result<Rational, ImportRefusal> {
        real(self.param(params, i)?)
    }

    fn length(&self, params: &[Token], i: usize) -> Result<Rational, ImportRefusal> {
        Ok(&self.real(params, i)? * &self.scale)
    }

    fn positive(&self, params: &[Token], i: usize) -> Result<Rational, ImportRefusal> {
        let v = self.length(params, i)?;
        if v.sign() != catallax_numbers::Sign::Positive {
            return Err(malformed("IGES radius must be positive"));
        }
        Ok(v)
    }

    fn triple(&self, params: &[Token], at: usize) -> Result<Vec3, ImportRefusal> {
        Ok([
            self.real(params, at)?,
            self.real(params, at + 1)?,
            self.real(params, at + 2)?,
        ])
    }

    /// The motion an entry's transformation matrix pointer gives, in millimetres.
    fn motion(&self, sequence: usize) -> Result<Motion, ImportRefusal> {
        let mut motion = Motion::identity();
        let mut next = self.entry(sequence)?.transform;
        for _ in 0..MAX_TRANSFORM_DEPTH {
            if next == 0 {
                return Ok(motion);
            }
            let e = self.entry(next)?;
            if e.kind != 124 {
                return Err(unsupported(next, kind_text(e.kind, e.form)));
            }
            let p = &e.params;
            let rows: [[Rational; 3]; 3] = [
                [self.real(p, 0)?, self.real(p, 1)?, self.real(p, 2)?],
                [self.real(p, 4)?, self.real(p, 5)?, self.real(p, 6)?],
                [self.real(p, 8)?, self.real(p, 9)?, self.real(p, 10)?],
            ];
            let shift = [self.length(p, 3)?, self.length(p, 7)?, self.length(p, 11)?];
            // Only a rigid motion keeps circles, cylinders and their radii exact.
            let orthonormal = (0..3).all(|i| {
                (0..3).all(|j| {
                    let column = |k: usize| -> Vec3 { std::array::from_fn(|r| rows[r][k].clone()) };
                    dot(&column(i), &column(j))
                        == if i == j {
                            Rational::one()
                        } else {
                            Rational::zero()
                        }
                })
            });
            if !orthonormal {
                return Err(unsupported(
                    next,
                    "IGES 124 (TRANSFORMATION MATRIX that is not exactly orthonormal)",
                ));
            }
            motion = Motion { rows, shift }.after(&motion);
            next = e.transform;
        }
        Err(unsupported(
            sequence,
            "IGES 124 (TRANSFORMATION MATRIX chain longer than eight)",
        ))
    }

    fn untransformed(&self, sequence: usize) -> Result<(), ImportRefusal> {
        let e = self.entry(sequence)?;
        if e.transform != 0 {
            return Err(unsupported(
                sequence,
                format!("{} with a transformation matrix", kind_text(e.kind, e.form)),
            ));
        }
        Ok(())
    }

    fn point(&self, sequence: usize) -> Result<Vec3, ImportRefusal> {
        let p = self.params(sequence, &[116])?;
        let local = [self.length(p, 0)?, self.length(p, 1)?, self.length(p, 2)?];
        Ok(self.motion(sequence)?.point(&local))
    }

    fn direction(&self, sequence: usize) -> Result<Vec3, ImportRefusal> {
        let p = self.params(sequence, &[123])?;
        let v = self.triple(p, 0)?;
        if v.iter().all(Rational::is_zero) {
            return Err(malformed(format!("zero IGES direction D{sequence}")));
        }
        Ok(self.motion(sequence)?.rotate(&v))
    }

    fn optional_direction(&self, p: &[Token], i: usize) -> Result<Option<Vec3>, ImportRefusal> {
        match p.get(i) {
            None | Some(Token::Empty) => Ok(None),
            Some(token) => match integer(token)? {
                0 => Ok(None),
                sequence => self.direction(sequence).map(Some),
            },
        }
    }

    fn file(&self) -> Result<StepFile, ImportRefusal> {
        let mut solids = Vec::new();
        for (&sequence, e) in self.entries {
            if e.kind == 186 {
                solids.push(self.solid(sequence)?);
            }
        }
        if solids.is_empty() {
            // A surface model is not a solid; name its first bounded surface.
            if let Some((&sequence, e)) = self
                .entries
                .iter()
                .find(|(_, e)| matches!(e.kind, 143 | 144 | 514))
            {
                return Err(unsupported(sequence, kind_text(e.kind, e.form)));
            }
            return Err(malformed(
                "no MANIFOLD SOLID B-REP OBJECT (type 186) in the file",
            ));
        }
        Ok(StepFile {
            mm_per_unit: self.scale.clone(),
            solids,
        })
    }

    fn solid(&self, sequence: usize) -> Result<Solid, ImportRefusal> {
        self.untransformed(sequence)?;
        let p = self.params(sequence, &[186])?;
        let shell = self.int(p, 0)?;
        let agrees = self.int(p, 1)? == 1;
        if self.int(p, 2)? != 0 {
            return Err(unsupported(
                sequence,
                "IGES 186 (MANIFOLD SOLID B-REP OBJECT with voids)",
            ));
        }
        self.untransformed(shell)?;
        let e = self.entry(shell)?;
        if e.kind != 514 || e.form != 1 {
            return Err(unsupported(shell, kind_text(e.kind, e.form)));
        }
        let s = &e.params;
        let count = self.int(s, 0)?;
        let mut faces = Vec::with_capacity(count);
        for i in 0..count {
            let face = self.int(s, 1 + 2 * i)?;
            let same = (self.int(s, 2 + 2 * i)? == 1) == agrees;
            faces.push(self.face(face, same)?);
        }
        Ok(Solid {
            name: format!("D{sequence}"),
            faces,
        })
    }

    fn face(&self, sequence: usize, same_sense: bool) -> Result<Face, ImportRefusal> {
        self.untransformed(sequence)?;
        let p = self.params(sequence, &[510])?;
        let surface = self.surface(self.int(p, 0)?)?;
        let count = self.int(p, 1)?;
        let outer = self.int(p, 2)? == 1;
        let mut loops = Vec::with_capacity(count);
        for i in 0..count {
            // IGES loops run with the face on their left about the surface
            // normal; a face against its surface reverses them.
            let mut bound = self.loop_(self.int(p, 3 + i)?)?;
            bound.outer = outer && i == 0;
            bound.orientation = same_sense;
            loops.push(bound);
        }
        Ok(Face {
            name: format!("D{sequence}"),
            surface,
            same_sense,
            loops,
        })
    }

    fn loop_(&self, sequence: usize) -> Result<FaceLoop, ImportRefusal> {
        self.untransformed(sequence)?;
        let p = self.params(sequence, &[508])?;
        let count = self.int(p, 0)?;
        let mut at = 1;
        let mut edges = Vec::with_capacity(count);
        let mut vertex = None;
        for _ in 0..count {
            let kind = self.int(p, at)?;
            let list = self.int(p, at + 1)?;
            let index = self.int(p, at + 2)?;
            let forward = self.int(p, at + 3)? == 1;
            let curves = self.int(p, at + 4)?;
            at += 5 + 2 * curves;
            match kind {
                0 => edges.push(self.edge(list, index, forward)?),
                1 if count == 1 => vertex = Some(self.vertex(list, index)?),
                // A vertex among edges is a degenerate edge, such as a sphere's
                // pole between two uses of its seam. It bounds nothing; the
                // chain check still requires the edges to meet there.
                1 => {
                    self.vertex(list, index)?;
                }
                _ => return Err(malformed(format!("D{sequence} has loop use type {kind}"))),
            }
        }
        Ok(FaceLoop {
            outer: false,
            orientation: true,
            edges,
            vertex,
        })
    }

    fn vertex(&self, list: usize, index: usize) -> Result<Vertex, ImportRefusal> {
        self.untransformed(list)?;
        let e = self.entry(list)?;
        if e.kind != 502 || e.form != 1 {
            return Err(unsupported(list, kind_text(e.kind, e.form)));
        }
        let count = self.int(&e.params, 0)?;
        if index == 0 || index > count {
            return Err(malformed(format!("D{list} has no vertex {index}")));
        }
        let at = 1 + 3 * (index - 1);
        Ok(Vertex {
            name: format!("D{list}:{index}"),
            point: [
                self.length(&e.params, at)?,
                self.length(&e.params, at + 1)?,
                self.length(&e.params, at + 2)?,
            ],
        })
    }

    fn edge(
        &self,
        list: usize,
        index: usize,
        forward: bool,
    ) -> Result<OrientedEdge, ImportRefusal> {
        self.untransformed(list)?;
        let e = self.entry(list)?;
        if e.kind != 504 || e.form != 1 {
            return Err(unsupported(list, kind_text(e.kind, e.form)));
        }
        let count = self.int(&e.params, 0)?;
        if index == 0 || index > count {
            return Err(malformed(format!("D{list} has no edge {index}")));
        }
        let row = 1 + 5 * (index - 1);
        let start = self.vertex(self.int(&e.params, row + 1)?, self.int(&e.params, row + 2)?)?;
        let end = self.vertex(self.int(&e.params, row + 3)?, self.int(&e.params, row + 4)?)?;
        let (curve, curve_same_sense) = self.curve(self.int(&e.params, row)?, &start, &end)?;
        Ok(OrientedEdge {
            edge: format!("D{list}:{index}"),
            start,
            end,
            curve_same_sense,
            curve,
            forward,
        })
    }

    /// The model-space curve of an edge, and whether the edge runs along it.
    fn curve(
        &self,
        sequence: usize,
        start: &Vertex,
        end: &Vertex,
    ) -> Result<(Curve, bool), ImportRefusal> {
        let e = self.entry(sequence)?;
        let motion = self.motion(sequence)?;
        let p = &e.params;
        let line = |a: Vec3, b: Vec3| -> Result<(Curve, bool), ImportRefusal> {
            let direction: Vec3 = std::array::from_fn(|i| &b[i] - &a[i]);
            let travel: Vec3 = std::array::from_fn(|i| &end.point[i] - &start.point[i]);
            if direction.iter().all(Rational::is_zero) {
                return Err(malformed(format!("D{sequence} is a line of zero length")));
            }
            let along = dot(&direction, &travel);
            Ok((
                Curve::Line(a, direction),
                along.sign() != catallax_numbers::Sign::Negative,
            ))
        };
        match e.kind {
            110 => {
                let a = motion.point(&[self.length(p, 0)?, self.length(p, 1)?, self.length(p, 2)?]);
                let b = motion.point(&[self.length(p, 3)?, self.length(p, 4)?, self.length(p, 5)?]);
                line(a, b)
            }
            100 => {
                let z = self.length(p, 0)?;
                let centre = [self.length(p, 1)?, self.length(p, 2)?, z];
                let dx = &self.length(p, 3)? - &centre[0];
                let dy = &self.length(p, 4)? - &centre[1];
                let radius = exact_sqrt(&(&dx * &dx + &dy * &dy)).ok_or_else(|| {
                    unsupported(
                        sequence,
                        "IGES 100 (CIRCULAR ARC whose radius is not rational)",
                    )
                })?;
                if radius.is_zero() {
                    return Err(malformed(format!("D{sequence} is an arc of zero radius")));
                }
                let x = motion.rotate(&[Rational::one(), Rational::zero(), Rational::zero()]);
                let y = motion.rotate(&[Rational::zero(), Rational::one(), Rational::zero()]);
                // A reflection reverses the arc's turning sense with its axis.
                let placement = placement(motion.point(&centre), cross(&x, &y), Some(x))?;
                Ok((Curve::Circle(placement, radius), true))
            }
            126 => {
                let spline = self.spline_curve(sequence, &motion)?;
                let geometry = &spline.geometry;
                // A straight two-point segment, as seams and generators are written.
                if geometry.degree() == 1 && geometry.controls().len() == 2 {
                    let [a, b] = [&geometry.controls()[0], &geometry.controls()[1]];
                    if geometry.weights()[0] == geometry.weights()[1] {
                        return line(a.clone(), b.clone());
                    }
                }
                Ok((Curve::BSpline(spline), true))
            }
            kind => Err(unsupported(sequence, kind_text(kind, e.form))),
        }
    }

    fn spline_curve(
        &self,
        sequence: usize,
        motion: &Motion,
    ) -> Result<StepSplineCurve, ImportRefusal> {
        let p = &self.entry(sequence)?.params;
        let upper = self.int(p, 0)?;
        let degree = self.int(p, 1)?;
        let closed = self.int(p, 3)? == 1;
        let count = upper + 1;
        let knot_count = upper + degree + 2;
        let mut at = 6;
        let knots = (0..knot_count)
            .map(|i| self.real(p, at + i))
            .collect::<Result<Vec<_>, _>>()?;
        at += knot_count;
        let weights = (0..count)
            .map(|i| self.real(p, at + i))
            .collect::<Result<Vec<_>, _>>()?;
        at += count;
        let controls = (0..count)
            .map(|i| {
                Ok(motion.point(&[
                    self.length(p, at + 3 * i)?,
                    self.length(p, at + 3 * i + 1)?,
                    self.length(p, at + 3 * i + 2)?,
                ]))
            })
            .collect::<Result<Vec<_>, ImportRefusal>>()?;
        let geometry = RationalSplineCurve::new(degree, knots, controls, weights, self.dec)
            .map_err(source_error)?
            .into_clamped(self.dec)
            .map_err(source_error)?;
        Ok(StepSplineCurve {
            geometry,
            closed: Some(closed),
            self_intersect: None,
        })
    }

    fn surface(&self, sequence: usize) -> Result<Surface, ImportRefusal> {
        let e = self.entry(sequence)?;
        let motion = self.motion(sequence)?;
        let p = &e.params;
        let located = |at: usize| -> Result<Vec3, ImportRefusal> {
            Ok(motion.point(&self.point(self.int(p, at)?)?))
        };
        let oriented = |at: usize| -> Result<Vec3, ImportRefusal> {
            Ok(motion.rotate(&self.direction(self.int(p, at)?)?))
        };
        let reference = |at: usize| -> Result<Option<Vec3>, ImportRefusal> {
            Ok(self.optional_direction(p, at)?.map(|v| motion.rotate(&v)))
        };
        Ok(match (e.kind, e.form) {
            (190, 0 | 1) => Surface::Plane(placement(located(0)?, oriented(1)?, reference(2)?)?),
            (192, 0 | 1) => Surface::Cylinder(
                placement(located(0)?, oriented(1)?, reference(3)?)?,
                self.positive(p, 2)?,
            ),
            (194, 0 | 1) => {
                let place = placement(located(0)?, oriented(1)?, reference(4)?)?;
                let radius = self.length(p, 2)?;
                if radius.sign() == catallax_numbers::Sign::Negative {
                    return Err(malformed("IGES cone radius must not be negative"));
                }
                let angle = self.real(p, 3)?;
                if angle == Rational::from_i64(45) {
                    Surface::Cone45(place, radius)
                } else {
                    Surface::Cone(
                        place,
                        radius,
                        SemiAngle {
                            value: angle,
                            degrees: true,
                        },
                    )
                }
            }
            (196, 0 | 1) => {
                let axis = match self.optional_direction(p, 2)? {
                    Some(axis) => motion.rotate(&axis),
                    None => motion.rotate(&[Rational::zero(), Rational::zero(), Rational::one()]),
                };
                Surface::Sphere(
                    placement(located(0)?, axis, reference(3)?)?,
                    self.positive(p, 1)?,
                )
            }
            (198, 0 | 1) => {
                let (major, minor) = (self.positive(p, 2)?, self.positive(p, 3)?);
                if minor >= major {
                    return Err(unsupported(
                        sequence,
                        "IGES 198 (TOROIDAL SURFACE whose minor radius is not below its major radius)",
                    ));
                }
                Surface::Torus(
                    placement(located(0)?, oriented(1)?, reference(4)?)?,
                    major,
                    minor,
                )
            }
            (128, _) => Surface::BSpline(self.spline_surface(sequence, &motion)?),
            (kind, form) => return Err(unsupported(sequence, kind_text(kind, form))),
        })
    }

    fn spline_surface(
        &self,
        sequence: usize,
        motion: &Motion,
    ) -> Result<StepSplineSurface, ImportRefusal> {
        let p = &self.entry(sequence)?.params;
        let upper = [self.int(p, 0)?, self.int(p, 1)?];
        let degrees = [self.int(p, 2)?, self.int(p, 3)?];
        let closed = [Some(self.int(p, 4)? == 1), Some(self.int(p, 5)? == 1)];
        let count = [upper[0] + 1, upper[1] + 1];
        let mut at = 9;
        let mut knots: [Vec<Rational>; 2] = Default::default();
        for d in 0..2 {
            let n = upper[d] + degrees[d] + 2;
            knots[d] = (0..n)
                .map(|i| self.real(p, at + i))
                .collect::<Result<_, _>>()?;
            at += n;
        }
        let total = count[0] * count[1];
        // IGES lists weights and controls with the first index varying fastest;
        // the spline tools take rows of constant u.
        let index = |u: usize, v: usize| u + count[0] * v;
        let mut weights = Vec::with_capacity(total);
        let mut controls = Vec::with_capacity(total);
        for u in 0..count[0] {
            for v in 0..count[1] {
                weights.push(self.real(p, at + index(u, v))?);
                let c = at + total + 3 * index(u, v);
                controls.push(motion.point(&[
                    self.length(p, c)?,
                    self.length(p, c + 1)?,
                    self.length(p, c + 2)?,
                ]));
            }
        }
        let [u_knots, v_knots] = knots;
        let geometry = RationalSplineSurface::new(
            degrees,
            [u_knots, v_knots],
            count,
            controls,
            weights,
            self.dec,
        )
        .map_err(source_error)?
        .into_clamped(self.dec)
        .map_err(source_error)?;
        Ok(StepSplineSurface {
            geometry,
            closed,
            self_intersect: None,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn hollerith_and_delimiters() {
        let (pd, rd) = delimiters("1H,,1H;,4HPart,2HMM;");
        assert_eq!((pd, rd), (',', ';'));
        let t = tokens("1H,,1H;,4HP,rt,,3.5D0;", pd, rd).unwrap();
        assert!(matches!(&t[2], Token::Text(s) if s == "P,rt"));
        assert!(matches!(t[3], Token::Empty));
        assert_eq!(real(&t[4]).unwrap(), Rational::from_ratio(7, 2).unwrap());
        let (pd, rd) = delimiters("1H//1H#/");
        assert_eq!((pd, rd), ('/', '#'));
    }

    #[test]
    fn reals_read_exactly() {
        for (text, n, d) in [
            ("1.", 1, 1),
            (".5", 1, 2),
            ("-.25", -1, 4),
            ("+2.5D-1", 1, 4),
            ("1.E1", 10, 1),
            ("3", 3, 1),
        ] {
            assert_eq!(
                real(&Token::Number(text.into())).unwrap(),
                Rational::from_ratio(n, d).unwrap(),
                "{text}"
            );
        }
    }

    #[test]
    fn units_convert_exactly() {
        let global = |flag: &str, name: &str| {
            let mut t = vec![Token::Empty; 12];
            t.push(Token::Number("1.0".into()));
            t.push(Token::Number(flag.into()));
            t.push(Token::Text(name.into()));
            t
        };
        assert_eq!(
            units(&global("1", "IN")).unwrap(),
            Rational::from_ratio(127, 5).unwrap()
        );
        assert_eq!(units(&global("2", "MM")).unwrap(), Rational::one());
        assert_eq!(units(&global("3", "CM")).unwrap(), Rational::from_i64(10));
        assert!(matches!(
            units(&global("3", "FURLONG")),
            Err(ImportRefusal::UnknownUnit { .. })
        ));
        let mut scaled = global("2", "MM");
        scaled[12] = Token::Number("2.0".into());
        assert!(matches!(
            units(&scaled),
            Err(ImportRefusal::UnknownUnit { .. })
        ));
    }
}

//! Binary STL, ASCII STL and Wavefront OBJ triangle readers.
//!
//! Every coordinate is read as the exact rational the file states: a binary32
//! value in binary STL, a decimal in ASCII STL and OBJ. Vertices are welded by
//! exact coordinate equality, never by distance.

use std::collections::HashMap;

use catallax_numbers::{Rational, Refusal};

/// Exact vertices in first-use order and triangles as vertex indices.
pub(crate) struct Soup {
    pub points: Vec<[Rational; 3]>,
    pub triangles: Vec<[u32; 3]>,
    /// OBJ vertices that no triangle uses, in file order, each once and none
    /// equal to a triangle corner. A closed mesh ignores them; a scan keeps
    /// them as isolated samples.
    pub isolated: Vec<[Rational; 3]>,
}

fn malformed(detail: impl Into<String>) -> Refusal {
    Refusal::new("import.mesh-malformed", detail)
}

fn too_large(limit: usize, count: Option<usize>) -> Refusal {
    Refusal::new(
        "import.mesh-too-large",
        match count {
            Some(n) => format!("the file holds {n} triangles; the stated limit is {limit}"),
            None => format!("the file holds more than {limit} triangles, the stated limit"),
        },
    )
}

/// Reads a file of any accepted format. Formats are told apart by content.
pub(crate) fn read(bytes: &[u8], limit: usize) -> Result<Soup, Refusal> {
    read_with(bytes, limit, usize::MAX)
}

/// [`read`] that also refuses an OBJ file stating more than `vertex_limit`
/// vertices, so that a file of isolated samples is bounded too.
pub(crate) fn read_with(bytes: &[u8], limit: usize, vertex_limit: usize) -> Result<Soup, Refusal> {
    if bytes.len() >= 84 {
        let count = u32::from_le_bytes(bytes[80..84].try_into().expect("four bytes")) as usize;
        if count
            .checked_mul(50)
            .and_then(|n| n.checked_add(84))
            .is_some_and(|n| n == bytes.len())
        {
            return binary_stl(bytes, count, limit);
        }
    }
    let text = std::str::from_utf8(bytes)
        .map_err(|_| malformed("not binary STL, and not UTF-8 text for ASCII STL or OBJ"))?;
    match text.split_ascii_whitespace().next() {
        Some("solid") => ascii_stl(text, limit),
        _ => obj(text, limit, vertex_limit),
    }
}

#[derive(Default)]
struct Welder<K> {
    index: HashMap<K, u32>,
    points: Vec<[Rational; 3]>,
}

impl<K: std::hash::Hash + Eq> Welder<K> {
    fn add(&mut self, key: K, point: impl FnOnce() -> [Rational; 3]) -> u32 {
        let next = self.points.len() as u32;
        *self.index.entry(key).or_insert_with(|| {
            self.points.push(point());
            next
        })
    }
}

fn binary_stl(bytes: &[u8], count: usize, limit: usize) -> Result<Soup, Refusal> {
    if count > limit {
        return Err(too_large(limit, Some(count)));
    }
    let mut welder: Welder<[u32; 3]> = Welder::default();
    let mut triangles = Vec::with_capacity(count);
    for t in 0..count {
        let record = &bytes[84 + 50 * t..84 + 50 * (t + 1)];
        let mut corner = [0u32; 3];
        for (v, slot) in corner.iter_mut().enumerate() {
            let mut bits = [0u32; 3];
            let mut values = [0f32; 3];
            for axis in 0..3 {
                let at = 12 + 12 * v + 4 * axis;
                let value = f32::from_le_bytes(record[at..at + 4].try_into().expect("four bytes"));
                if !value.is_finite() {
                    return Err(malformed(format!(
                        "triangle {t} has a coordinate that is not a finite number"
                    )));
                }
                // Positive and negative zero are the same point; adding
                // positive zero turns negative zero into positive zero.
                let value = value + 0.0;
                values[axis] = value;
                bits[axis] = value.to_bits();
            }
            *slot = welder.add(bits, || {
                values.map(|x| Rational::from_f64_exact(f64::from(x)).expect("finite binary32"))
            });
        }
        triangles.push(corner);
    }
    Ok(Soup {
        points: welder.points,
        triangles,
        isolated: Vec::new(),
    })
}

/// Parses one decimal coordinate exactly. Exponents beyond the binary64 range
/// are refused so that no single number can exhaust memory.
fn decimal(token: &str) -> Result<Rational, Refusal> {
    let bad = || malformed(format!("{token:?} is not a decimal number"));
    if token.len() > 64 {
        return Err(bad());
    }
    let (mantissa, exponent) = match token.find(['e', 'E']) {
        Some(at) => (&token[..at], Some(&token[at + 1..])),
        None => (token, None),
    };
    if let Some(e) = exponent {
        let value: i64 = e.parse().map_err(|_| bad())?;
        if value.abs() > 400 {
            return Err(bad());
        }
    }
    let (sign, body) = match mantissa.strip_prefix(['+', '-']) {
        Some(rest) => (&mantissa[..1], rest),
        None => ("", mantissa),
    };
    let mut body = body.to_string();
    if body.starts_with('.') {
        body.insert(0, '0');
    }
    if body.ends_with('.') {
        body.push('0');
    }
    let normalised = match exponent {
        Some(e) => format!("{sign}{body}e{e}"),
        None => format!("{sign}{body}"),
    };
    Rational::parse(&normalised).map_err(|_| bad())
}

fn ascii_stl(text: &str, limit: usize) -> Result<Soup, Refusal> {
    let mut tokens = text.split_ascii_whitespace().peekable();
    let mut welder: Welder<[Rational; 3]> = Welder::default();
    let mut triangles = Vec::new();
    let expect = |tokens: &mut std::iter::Peekable<std::str::SplitAsciiWhitespace<'_>>,
                  word: &str|
     -> Result<(), Refusal> {
        match tokens.next() {
            Some(t) if t == word => Ok(()),
            Some(t) => Err(malformed(format!(
                "ASCII STL: expected {word:?}, found {t:?}"
            ))),
            None => Err(malformed(format!(
                "ASCII STL: expected {word:?}, found the end"
            ))),
        }
    };
    while let Some(token) = tokens.next() {
        if token != "solid" {
            return Err(malformed(format!(
                "ASCII STL: expected \"solid\", found {token:?}"
            )));
        }
        // The solid's name runs up to its first facet or its end.
        while tokens
            .peek()
            .is_some_and(|t| *t != "facet" && *t != "endsolid")
        {
            tokens.next();
        }
        loop {
            match tokens.next() {
                Some("facet") => {}
                Some("endsolid") => break,
                Some(t) => {
                    return Err(malformed(format!(
                        "ASCII STL: expected \"facet\" or \"endsolid\", found {t:?}"
                    )));
                }
                None => return Err(malformed("ASCII STL: the file ends inside a solid")),
            }
            expect(&mut tokens, "normal")?;
            for _ in 0..3 {
                // The stated normal is redundant with the vertex order and is not read.
                tokens
                    .next()
                    .ok_or_else(|| malformed("ASCII STL: the file ends inside a facet normal"))?;
            }
            expect(&mut tokens, "outer")?;
            expect(&mut tokens, "loop")?;
            let mut corner = [0u32; 3];
            for slot in &mut corner {
                expect(&mut tokens, "vertex")?;
                let mut point: [Rational; 3] = Default::default();
                for value in &mut point {
                    let token = tokens
                        .next()
                        .ok_or_else(|| malformed("ASCII STL: the file ends inside a vertex"))?;
                    *value = decimal(token)?;
                }
                *slot = welder.add(point.clone(), || point);
            }
            expect(&mut tokens, "endloop")?;
            expect(&mut tokens, "endfacet")?;
            if triangles.len() == limit {
                return Err(too_large(limit, None));
            }
            triangles.push(corner);
        }
        // The end line may repeat the name; another solid may follow.
        while tokens.peek().is_some_and(|t| *t != "solid") {
            tokens.next();
        }
    }
    Ok(Soup {
        points: welder.points,
        triangles,
        isolated: Vec::new(),
    })
}

fn obj(text: &str, limit: usize, vertex_limit: usize) -> Result<Soup, Refusal> {
    let mut stated: Vec<[Rational; 3]> = Vec::new();
    let mut welded: Vec<Option<u32>> = Vec::new();
    let mut welder: Welder<[Rational; 3]> = Welder::default();
    let mut triangles = Vec::new();
    for (number, line) in text.lines().enumerate() {
        let line = line.split('#').next().unwrap_or("");
        let mut words = line.split_ascii_whitespace();
        let Some(keyword) = words.next() else {
            continue;
        };
        let at = |detail: String| malformed(format!("OBJ line {}: {detail}", number + 1));
        match keyword {
            "v" => {
                let values = words.map(decimal).collect::<Result<Vec<_>, _>>()?;
                let point = match values.as_slice() {
                    [x, y, z] => [x.clone(), y.clone(), z.clone()],
                    [x, y, z, w] if w.sign() != catallax_numbers::Sign::Zero => [
                        x.checked_div(w).expect("nonzero"),
                        y.checked_div(w).expect("nonzero"),
                        z.checked_div(w).expect("nonzero"),
                    ],
                    _ => {
                        return Err(at(
                            "a vertex needs three coordinates and an optional nonzero weight"
                                .into(),
                        ));
                    }
                };
                if stated.len() == vertex_limit {
                    return Err(Refusal::new(
                        "import.mesh-too-large",
                        format!(
                            "the file holds more than {vertex_limit} vertices, the stated limit"
                        ),
                    ));
                }
                stated.push(point);
                welded.push(None);
            }
            "f" => {
                let corners = words
                    .map(|word| {
                        let index: i64 = word
                            .split('/')
                            .next()
                            .and_then(|i| i.parse().ok())
                            .ok_or_else(|| at(format!("{word:?} is not a vertex reference")))?;
                        let resolved = if index > 0 {
                            index - 1
                        } else {
                            stated.len() as i64 + index
                        };
                        usize::try_from(resolved)
                            .ok()
                            .filter(|&i| index != 0 && i < stated.len())
                            .ok_or_else(|| at(format!("vertex reference {index} is out of range")))
                    })
                    .collect::<Result<Vec<_>, _>>()?;
                let [a, b, c] = corners[..] else {
                    return Err(at(format!(
                        "a face has {} corners; only triangles are read, never split",
                        corners.len()
                    )));
                };
                let mut corner = [0u32; 3];
                for (slot, i) in corner.iter_mut().zip([a, b, c]) {
                    *slot = *welded[i]
                        .get_or_insert_with(|| welder.add(stated[i].clone(), || stated[i].clone()));
                }
                if triangles.len() == limit {
                    return Err(too_large(limit, None));
                }
                triangles.push(corner);
            }
            "vn" | "vt" | "vp" | "g" | "o" | "s" | "usemtl" | "mtllib" => {}
            other => {
                return Err(at(format!(
                    "{other:?} is not read; a closed mesh file holds vertices and triangles"
                )));
            }
        }
    }
    if triangles.is_empty() && stated.is_empty() {
        return Err(malformed(
            "the file holds no vertices or triangles in a known format",
        ));
    }
    // Unused vertices are welded after every triangle corner, so one equal
    // to a corner is already a point of the surface and is not repeated.
    let corners = welder.points.len();
    for (point, used) in stated.iter().zip(&welded) {
        if used.is_none() {
            welder.add(point.clone(), || point.clone());
        }
    }
    let isolated = welder.points.split_off(corners);
    Ok(Soup {
        points: welder.points,
        triangles,
        isolated,
    })
}

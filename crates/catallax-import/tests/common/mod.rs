//! Public synthetic STEP fixtures: axis-aligned boxes written as Part 21 text.
#![allow(dead_code)]

use catallax_numbers::Rational;

#[allow(unused_imports)]
pub use soup::{Caps, Hole, Soup};
mod soup;

pub fn r(s: &str) -> Rational {
    Rational::parse(s).unwrap()
}

pub(crate) fn dec(v: &Rational) -> String {
    let t = v.to_decimal_exact().expect("terminating decimal");
    if t.contains('.') { t } else { format!("{t}.0") }
}

#[derive(Clone, Copy)]
pub enum Unit {
    Millimetre,
    Metre,
    Inch,
}

/// One box. `shift` gives face `index` its own copy of its vertices and edges,
/// moved by the offset, which opens a gap along its four seams.
pub struct BoxFixture {
    pub origin: [Rational; 3],
    pub size: [Rational; 3],
    pub unit: Unit,
    pub shift: Option<(usize, [Rational; 3])>,
    /// A second shifted face.
    pub also: Option<(usize, [Rational; 3])>,
    /// Leaves this face out of the shell, which opens the shell.
    pub omit: Option<usize>,
}

impl BoxFixture {
    pub fn new(size: [&str; 3]) -> Self {
        Self {
            origin: [r("0"), r("0"), r("0")],
            size: size.map(r),
            unit: Unit::Millimetre,
            shift: None,
            also: None,
            omit: None,
        }
    }
}

pub(crate) struct Writer {
    pub(crate) lines: Vec<String>,
}

impl Writer {
    fn add(&mut self, body: String) -> u64 {
        self.lines.push(body);
        self.lines.len() as u64
    }
    fn point(&mut self, p: &[Rational; 3]) -> u64 {
        let n = self.lines.len() + 1;
        self.add(format!(
            "CARTESIAN_POINT('p{n}',({},{},{}))",
            dec(&p[0]),
            dec(&p[1]),
            dec(&p[2])
        ))
    }
    fn dir(&mut self, axis: usize) -> u64 {
        let v: Vec<&str> = (0..3)
            .map(|i| if i == axis { "1." } else { "0." })
            .collect();
        self.add(format!("DIRECTION('',({},{},{}))", v[0], v[1], v[2]))
    }
}

/// Entity names of the solids in the returned text, `#id` per box in order.
pub fn boxes_step(boxes: &[BoxFixture]) -> String {
    let mut w = Writer { lines: Vec::new() };
    let mut solids = Vec::new();
    for b in boxes {
        solids.push(write_box(&mut w, b));
    }
    let mut out = String::from(
        "ISO-10303-21;\nHEADER;\nFILE_DESCRIPTION(('synthetic'),'2;1');\nFILE_NAME('b.step','',(''),(''),'','','');\nFILE_SCHEMA(('AUTOMOTIVE_DESIGN'));\nENDSEC;\nDATA;\n",
    );
    for (i, line) in w.lines.iter().enumerate() {
        out.push_str(&format!("#{}={};\n", i + 1, line));
    }
    let base = w.lines.len() as u64;
    let unit = match boxes[0].unit {
        Unit::Millimetre => "( LENGTH_UNIT() NAMED_UNIT(*) SI_UNIT(.MILLI.,.METRE.) )".to_string(),
        Unit::Metre => "( LENGTH_UNIT() NAMED_UNIT(*) SI_UNIT($,.METRE.) )".to_string(),
        Unit::Inch => format!(
            "( CONVERSION_BASED_UNIT('INCH',#{}) LENGTH_UNIT() NAMED_UNIT(*) )",
            base + 3
        ),
    };
    out.push_str(&format!(
        "#{}=( GEOMETRIC_REPRESENTATION_CONTEXT(3) GLOBAL_UNIT_ASSIGNED_CONTEXT((#{})) REPRESENTATION_CONTEXT('',''));\n",
        base + 1,
        base + 2
    ));
    out.push_str(&format!("#{}={};\n", base + 2, unit));
    if matches!(boxes[0].unit, Unit::Inch) {
        out.push_str(&format!(
            "#{}=LENGTH_MEASURE_WITH_UNIT(LENGTH_MEASURE(25.4),#{});\n",
            base + 3,
            base + 4
        ));
        out.push_str(&format!(
            "#{}=( LENGTH_UNIT() NAMED_UNIT(*) SI_UNIT(.MILLI.,.METRE.) );\n",
            base + 4
        ));
    }
    out.push_str("ENDSEC;\nEND-ISO-10303-21;\n");
    let _ = solids;
    out
}

pub fn box_step(b: BoxFixture) -> String {
    boxes_step(&[b])
}

fn corner(b: &BoxFixture, c: [usize; 3]) -> [Rational; 3] {
    std::array::from_fn(|i| {
        if c[i] == 1 {
            &b.origin[i] + &b.size[i]
        } else {
            b.origin[i].clone()
        }
    })
}

struct Shared {
    vertex: [Option<u64>; 8],
    edge: std::collections::BTreeMap<(usize, usize), u64>,
}

fn corner_index(c: [usize; 3]) -> usize {
    c[0] + 2 * c[1] + 4 * c[2]
}

fn write_box(w: &mut Writer, b: &BoxFixture) -> u64 {
    let mut shared = Shared {
        vertex: [None; 8],
        edge: Default::default(),
    };
    let mut faces = Vec::new();
    for (index, (axis, side)) in (0..3).flat_map(|a| [(a, 0), (a, 1)]).enumerate() {
        if b.omit == Some(index) {
            continue;
        }
        let (bb, cc) = ((axis + 1) % 3, (axis + 2) % 3);
        let mut ring = [[0usize; 3]; 4];
        for (k, (ib, ic)) in [(0, 0), (1, 0), (1, 1), (0, 1)].into_iter().enumerate() {
            ring[k][axis] = side;
            ring[k][bb] = ib;
            ring[k][cc] = ic;
        }
        if side == 0 {
            ring.reverse();
        }
        let offset = [&b.shift, &b.also]
            .into_iter()
            .flatten()
            .find(|(i, _)| *i == index)
            .map(|(_, d)| d.clone());
        let mut own: std::collections::BTreeMap<usize, u64> = Default::default();
        let mut own_edge: std::collections::BTreeMap<(usize, usize), u64> = Default::default();
        let mut vertex = |w: &mut Writer, shared: &mut Shared, c: [usize; 3]| -> u64 {
            let ci = corner_index(c);
            if let Some(d) = &offset {
                *own.entry(ci).or_insert_with(|| {
                    let p: [Rational; 3] = std::array::from_fn(|i| &corner(b, c)[i] + &d[i]);
                    let pt = w.point(&p);
                    w.add(format!("VERTEX_POINT('',#{pt})"))
                })
            } else {
                if shared.vertex[ci].is_none() {
                    let pt = w.point(&corner(b, c));
                    shared.vertex[ci] = Some(w.add(format!("VERTEX_POINT('',#{pt})")));
                }
                shared.vertex[ci].unwrap()
            }
        };
        let mut oriented = Vec::new();
        for k in 0..4 {
            let (from, to) = (ring[k], ring[(k + 1) % 4]);
            let (fi, ti) = (corner_index(from), corner_index(to));
            let (lo, hi, forward) = if fi < ti {
                (from, to, true)
            } else {
                (to, from, false)
            };
            let key = (corner_index(lo), corner_index(hi));
            let existing = if offset.is_some() {
                own_edge.get(&key).copied()
            } else {
                shared.edge.get(&key).copied()
            };
            let edge = match existing {
                Some(e) => e,
                None => {
                    let v0 = vertex(w, &mut shared, lo);
                    let v1 = vertex(w, &mut shared, hi);
                    let diff = (0..3).find(|&i| lo[i] != hi[i]).unwrap();
                    let mut start = corner(b, lo);
                    if let Some(d) = &offset {
                        start = std::array::from_fn(|i| &start[i] + &d[i]);
                    }
                    let pt = w.point(&start);
                    let dir = w.dir(diff);
                    let vec = w.add(format!("VECTOR('',#{dir},{})", dec(&b.size[diff])));
                    let line = w.add(format!("LINE('',#{pt},#{vec})"));
                    let e = w.add(format!("EDGE_CURVE('',#{v0},#{v1},#{line},.T.)"));
                    if offset.is_some() {
                        own_edge.insert(key, e);
                    } else {
                        shared.edge.insert(key, e);
                    }
                    e
                }
            };
            oriented.push(w.add(format!(
                "ORIENTED_EDGE('',*,*,#{edge},.{}.)",
                if forward { "T" } else { "F" }
            )));
        }
        let list = oriented
            .iter()
            .map(|o| format!("#{o}"))
            .collect::<Vec<_>>()
            .join(",");
        let lp = w.add(format!("EDGE_LOOP('',({list}))"));
        let bound = w.add(format!("FACE_OUTER_BOUND('',#{lp},.T.)"));
        let mut at = corner(b, {
            let mut c = [0; 3];
            c[axis] = side;
            c
        });
        if let Some(d) = &offset {
            at = std::array::from_fn(|i| &at[i] + &d[i]);
        }
        let origin = w.point(&at);
        let (n, rd) = (w.dir(axis), w.dir(bb));
        let pl = w.add(format!("AXIS2_PLACEMENT_3D('',#{origin},#{n},#{rd})"));
        let plane = w.add(format!("PLANE('',#{pl})"));
        faces.push(w.add(format!(
            "ADVANCED_FACE('',(#{bound}),#{plane},.{}.)",
            if side == 1 { "T" } else { "F" }
        )));
    }
    let list = faces
        .iter()
        .map(|f| format!("#{f}"))
        .collect::<Vec<_>>()
        .join(",");
    let shell = w.add(format!("CLOSED_SHELL('',({list}))"));
    w.add(format!("MANIFOLD_SOLID_BREP('box',#{shell})"))
}

/// A solid given as planar faces, each a ring of corners counter-clockwise seen
/// from outside. Corners and straight edges are shared by coordinates.
pub type Ring = Vec<[&'static str; 3]>;

pub fn polyhedra_step(solids: &[Vec<Ring>]) -> String {
    let mut w = Writer { lines: Vec::new() };
    let mut names = Vec::new();
    for faces in solids {
        let mut vertex: std::collections::BTreeMap<[Rational; 3], u64> = Default::default();
        let mut edge: std::collections::BTreeMap<([Rational; 3], [Rational; 3]), u64> =
            Default::default();
        let mut face_ids = Vec::new();
        for ring in faces {
            let pts: Vec<[Rational; 3]> = ring.iter().map(|p| p.map(r)).collect();
            let mut oriented = Vec::new();
            for k in 0..pts.len() {
                let (from, to) = (&pts[k], &pts[(k + 1) % pts.len()]);
                let forward = from < to;
                let (lo, hi) = if forward { (from, to) } else { (to, from) };
                let key = (lo.clone(), hi.clone());
                if !edge.contains_key(&key) {
                    let mut ids = [0u64; 2];
                    for (slot, p) in [lo, hi].into_iter().enumerate() {
                        ids[slot] = *vertex.entry(p.clone()).or_insert_with(|| {
                            let pt = w.point(p);
                            w.add(format!("VERTEX_POINT('',#{pt})"))
                        });
                    }
                    let d: [Rational; 3] = std::array::from_fn(|i| &hi[i] - &lo[i]);
                    let origin = w.point(lo);
                    let dir = w.add(format!(
                        "DIRECTION('',({},{},{}))",
                        dec(&d[0]),
                        dec(&d[1]),
                        dec(&d[2])
                    ));
                    let vec = w.add(format!("VECTOR('',#{dir},1.0)"));
                    let line = w.add(format!("LINE('',#{origin},#{vec})"));
                    let e = w.add(format!(
                        "EDGE_CURVE('',#{},#{},#{line},.T.)",
                        ids[0], ids[1]
                    ));
                    edge.insert(key.clone(), e);
                }
                oriented.push(w.add(format!(
                    "ORIENTED_EDGE('',*,*,#{},.{}.)",
                    edge[&key],
                    if forward { "T" } else { "F" }
                )));
            }
            let list = oriented
                .iter()
                .map(|o| format!("#{o}"))
                .collect::<Vec<_>>()
                .join(",");
            let lp = w.add(format!("EDGE_LOOP('',({list}))"));
            let bound = w.add(format!("FACE_OUTER_BOUND('',#{lp},.T.)"));
            // Newell normal of the ring, scaled to integers is not needed: any exact normal.
            let n = {
                let mut n = [r("0"), r("0"), r("0")];
                for k in 0..pts.len() {
                    let (a, b) = (&pts[k], &pts[(k + 1) % pts.len()]);
                    n[0] = &n[0] + &(&(&a[1] - &b[1]) * &(&a[2] + &b[2]));
                    n[1] = &n[1] + &(&(&a[2] - &b[2]) * &(&a[0] + &b[0]));
                    n[2] = &n[2] + &(&(&a[0] - &b[0]) * &(&a[1] + &b[1]));
                }
                n
            };
            let origin = w.point(&pts[0]);
            let axis = w.add(format!(
                "DIRECTION('',({},{},{}))",
                dec(&n[0]),
                dec(&n[1]),
                dec(&n[2])
            ));
            let e0 = {
                let d: [Rational; 3] = std::array::from_fn(|i| &pts[1][i] - &pts[0][i]);
                w.add(format!(
                    "DIRECTION('',({},{},{}))",
                    dec(&d[0]),
                    dec(&d[1]),
                    dec(&d[2])
                ))
            };
            let pl = w.add(format!("AXIS2_PLACEMENT_3D('',#{origin},#{axis},#{e0})"));
            let plane = w.add(format!("PLANE('',#{pl})"));
            face_ids.push(w.add(format!("ADVANCED_FACE('',(#{bound}),#{plane},.T.)")));
        }
        let list = face_ids
            .iter()
            .map(|f| format!("#{f}"))
            .collect::<Vec<_>>()
            .join(",");
        let shell = w.add(format!("CLOSED_SHELL('',({list}))"));
        names.push(w.add(format!("MANIFOLD_SOLID_BREP('p',#{shell})")));
    }
    let mut out = String::from("ISO-10303-21;\nHEADER;\nENDSEC;\nDATA;\n");
    for (i, line) in w.lines.iter().enumerate() {
        out.push_str(&format!("#{}={};\n", i + 1, line));
    }
    let base = w.lines.len() as u64;
    out.push_str(&format!(
        "#{}=( GEOMETRIC_REPRESENTATION_CONTEXT(3) GLOBAL_UNIT_ASSIGNED_CONTEXT((#{})) REPRESENTATION_CONTEXT('',''));\n#{}=( LENGTH_UNIT() NAMED_UNIT(*) SI_UNIT(.MILLI.,.METRE.) );\nENDSEC;\nEND-ISO-10303-21;\n",
        base + 1,
        base + 2,
        base + 2
    ));
    out
}

/// The six faces of an axis-aligned box `[lo, hi]`, outward counter-clockwise.
pub fn box_rings(lo: [&'static str; 3], hi: [&'static str; 3]) -> Vec<Ring> {
    let p =
        |x: usize, y: usize, z: usize| [[lo[0], hi[0]][x], [lo[1], hi[1]][y], [lo[2], hi[2]][z]];
    vec![
        vec![p(0, 0, 0), p(0, 0, 1), p(0, 1, 1), p(0, 1, 0)],
        vec![p(1, 0, 0), p(1, 1, 0), p(1, 1, 1), p(1, 0, 1)],
        vec![p(0, 0, 0), p(1, 0, 0), p(1, 0, 1), p(0, 0, 1)],
        vec![p(0, 1, 0), p(0, 1, 1), p(1, 1, 1), p(1, 1, 0)],
        vec![p(0, 0, 0), p(0, 1, 0), p(1, 1, 0), p(1, 0, 0)],
        vec![p(0, 0, 1), p(1, 0, 1), p(1, 1, 1), p(0, 1, 1)],
    ]
}

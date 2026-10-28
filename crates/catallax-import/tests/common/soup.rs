//! Face soups with curved faces: boxes, cylinders, spheres and tori written into
//! one shell. Faces are written as stated; the soup is not required to be a valid
//! solid, which is what the self-intersection tests need.
#![allow(dead_code)]

use std::collections::BTreeMap;

use catallax_numbers::Rational;

use super::{Writer, dec, r};

type P3 = [Rational; 3];

/// A ring hole: centre, axis index, reference index, radius, traversed forward.
pub type Hole = (P3, usize, usize, Rational, bool);

#[derive(Clone, Copy, PartialEq, Eq)]
pub enum Caps {
    Both,
    TopOnly,
    None,
}

pub struct Soup {
    w: Writer,
    faces: Vec<u64>,
    vertex: BTreeMap<P3, u64>,
    line: BTreeMap<(P3, P3), u64>,
    circle: BTreeMap<String, u64>,
}

impl Soup {
    pub fn new() -> Self {
        Self {
            w: Writer { lines: Vec::new() },
            faces: Vec::new(),
            vertex: BTreeMap::new(),
            line: BTreeMap::new(),
            circle: BTreeMap::new(),
        }
    }

    fn vertex(&mut self, p: &P3) -> u64 {
        if let Some(v) = self.vertex.get(p) {
            return *v;
        }
        let pt = self.w.point(p);
        let v = self.w.add(format!("VERTEX_POINT('',#{pt})"));
        self.vertex.insert(p.clone(), v);
        v
    }

    fn dir(&mut self, d: &P3) -> u64 {
        self.w.add(format!(
            "DIRECTION('',({},{},{}))",
            dec(&d[0]),
            dec(&d[1]),
            dec(&d[2])
        ))
    }

    fn placement(&mut self, origin: &P3, axis: &P3, reference: &P3) -> u64 {
        let o = self.w.point(origin);
        let (a, b) = (self.dir(axis), self.dir(reference));
        self.w.add(format!("AXIS2_PLACEMENT_3D('',#{o},#{a},#{b})"))
    }

    fn unit(axis: usize) -> P3 {
        std::array::from_fn(|i| if i == axis { r("1") } else { r("0") })
    }

    fn at(base: &P3, axis: usize, by: &Rational) -> P3 {
        let mut p = base.clone();
        p[axis] = &p[axis] + by;
        p
    }

    /// A straight shared edge and its orientation as used from `from` to `to`.
    fn segment(&mut self, from: &P3, to: &P3) -> (u64, bool) {
        let forward = from < to;
        let (lo, hi) = if forward { (from, to) } else { (to, from) };
        let key = (lo.clone(), hi.clone());
        if let Some(e) = self.line.get(&key) {
            return (*e, forward);
        }
        let (v0, v1) = (self.vertex(lo), self.vertex(hi));
        let d: P3 = std::array::from_fn(|i| &hi[i] - &lo[i]);
        let origin = self.w.point(lo);
        let dir = self.dir(&d);
        let vec = self.w.add(format!("VECTOR('',#{dir},1.0)"));
        let line = self.w.add(format!("LINE('',#{origin},#{vec})"));
        let e = self
            .w
            .add(format!("EDGE_CURVE('',#{v0},#{v1},#{line},.T.)"));
        self.line.insert(key, e);
        (e, forward)
    }

    /// A circle edge about `axis` through its first vertex at `centre + radius * unit(ref_axis)`.
    fn ring(&mut self, centre: &P3, axis: usize, ref_axis: usize, radius: &Rational) -> u64 {
        let key = format!("{centre:?}{axis}{radius:?}");
        if let Some(e) = self.circle.get(&key) {
            return *e;
        }
        let first = Self::at(centre, ref_axis, radius);
        let v = self.vertex(&first);
        let pl = self.placement(centre, &Self::unit(axis), &Self::unit(ref_axis));
        let c = self.w.add(format!("CIRCLE('',#{pl},{})", dec(radius)));
        let e = self.w.add(format!("EDGE_CURVE('',#{v},#{v},#{c},.T.)"));
        self.circle.insert(key, e);
        e
    }

    fn oriented(&mut self, edge: u64, forward: bool) -> u64 {
        self.w.add(format!(
            "ORIENTED_EDGE('',*,*,#{edge},.{}.)",
            if forward { "T" } else { "F" }
        ))
    }

    /// A bound of one loop.
    fn bound(&mut self, oriented: &[u64], outer: bool) -> u64 {
        let list = oriented
            .iter()
            .map(|o| format!("#{o}"))
            .collect::<Vec<_>>()
            .join(",");
        let lp = self.w.add(format!("EDGE_LOOP('',({list}))"));
        self.w.add(format!(
            "{}('',#{lp},.T.)",
            if outer {
                "FACE_OUTER_BOUND"
            } else {
                "FACE_BOUND"
            }
        ))
    }

    fn finish_face(&mut self, bounds: &[u64], surface: u64, same_sense: bool) {
        let list = bounds
            .iter()
            .map(|b| format!("#{b}"))
            .collect::<Vec<_>>()
            .join(",");
        let f = self.w.add(format!(
            "ADVANCED_FACE('',({list}),#{surface},.{}.)",
            if same_sense { "T" } else { "F" }
        ));
        self.faces.push(f);
    }

    /// A planar face from an outer ring (outward counter-clockwise) and ring loops
    /// for holes, as `(centre, axis index, ref index, radius, forward)`.
    pub fn plane_face(
        &mut self,
        ring: &[P3],
        holes: &[Hole],
        normal_axis: usize,
        outward_positive: bool,
        offset: &P3,
    ) {
        let mut oriented = Vec::new();
        for k in 0..ring.len() {
            let (e, fwd) = self.segment(&ring[k], &ring[(k + 1) % ring.len()]);
            oriented.push(self.oriented(e, fwd));
        }
        let mut bounds = vec![self.bound(&oriented, true)];
        for (c, a, rf, rad, fwd) in holes {
            let e = self.ring(c, *a, *rf, rad);
            let o = self.oriented(e, *fwd);
            bounds.push(self.bound(&[o], false));
        }
        let pl = self.placement(
            offset,
            &Self::unit(normal_axis),
            &Self::unit((normal_axis + 1) % 3),
        );
        let surface = self.w.add(format!("PLANE('',#{pl})"));
        self.finish_face(&bounds, surface, outward_positive);
    }

    /// The six faces of a box. `holes` maps a face index (as in `box_rings`) to ring holes.
    pub fn box_with_holes(
        &mut self,
        lo: [&str; 3],
        hi: [&str; 3],
        holes: &BTreeMap<usize, Vec<Hole>>,
    ) {
        let (lo, hi) = (lo.map(r), hi.map(r));
        let p = |x: usize, y: usize, z: usize| -> P3 {
            [
                [lo[0].clone(), hi[0].clone()][x].clone(),
                [lo[1].clone(), hi[1].clone()][y].clone(),
                [lo[2].clone(), hi[2].clone()][z].clone(),
            ]
        };
        let rings: [Vec<P3>; 6] = [
            vec![p(0, 0, 0), p(0, 0, 1), p(0, 1, 1), p(0, 1, 0)],
            vec![p(1, 0, 0), p(1, 1, 0), p(1, 1, 1), p(1, 0, 1)],
            vec![p(0, 0, 0), p(1, 0, 0), p(1, 0, 1), p(0, 0, 1)],
            vec![p(0, 1, 0), p(0, 1, 1), p(1, 1, 1), p(1, 1, 0)],
            vec![p(0, 0, 0), p(0, 1, 0), p(1, 1, 0), p(1, 0, 0)],
            vec![p(0, 0, 1), p(1, 0, 1), p(1, 1, 1), p(0, 1, 1)],
        ];
        for (index, ring) in rings.iter().enumerate() {
            let (axis, side) = (index / 2, index % 2);
            let offset = ring[0].clone();
            let none = Vec::new();
            let h = holes.get(&index).unwrap_or(&none);
            self.plane_face(ring, h, axis, side == 1, &offset);
        }
    }

    pub fn box_(&mut self, lo: [&str; 3], hi: [&str; 3]) {
        self.box_with_holes(lo, hi, &BTreeMap::new());
    }

    /// A closed cylinder along `axis` from `base` to `base + height`, with both caps.
    /// `inward` flips the wall normal, for a hole wall (use `caps = false` then).
    pub fn cylinder(
        &mut self,
        base: [&str; 3],
        axis: usize,
        height: &str,
        radius: &str,
        caps: Caps,
    ) {
        let base: P3 = base.map(r);
        let (h, rad) = (r(height), r(radius));
        let ref_axis = (axis + 1) % 3;
        let top = Self::at(&base, axis, &h);
        let bottom_ring = self.ring(&base, axis, ref_axis, &rad);
        let top_ring = self.ring(&top, axis, ref_axis, &rad);
        let pl = self.placement(&base, &Self::unit(axis), &Self::unit(ref_axis));
        let wall = self
            .w
            .add(format!("CYLINDRICAL_SURFACE('',#{pl},{})", dec(&rad)));
        let (o1, o2) = (
            self.oriented(bottom_ring, true),
            self.oriented(top_ring, false),
        );
        let (b1, b2) = (self.bound(&[o1], true), self.bound(&[o2], false));
        self.finish_face(&[b1, b2], wall, true);
        {
            for (centre, ring, up) in [(&base, bottom_ring, false), (&top, top_ring, true)] {
                if (!up && caps != Caps::Both) || caps == Caps::None {
                    continue;
                }
                let o = self.oriented(ring, up);
                let b = self.bound(&[o], true);
                let pl = self.placement(centre, &Self::unit(axis), &Self::unit(ref_axis));
                let plane = self.w.add(format!("PLANE('',#{pl})"));
                self.finish_face(&[b], plane, up);
            }
        }
    }

    /// A cylinder wall only, with its normal reversed, as a hole wall.
    pub fn hole_wall(&mut self, base: [&str; 3], axis: usize, height: &str, radius: &str) {
        let base: P3 = base.map(r);
        let (h, rad) = (r(height), r(radius));
        let ref_axis = (axis + 1) % 3;
        let top = Self::at(&base, axis, &h);
        let (bottom_ring, top_ring) = (
            self.ring(&base, axis, ref_axis, &rad),
            self.ring(&top, axis, ref_axis, &rad),
        );
        let pl = self.placement(&base, &Self::unit(axis), &Self::unit(ref_axis));
        let wall = self
            .w
            .add(format!("CYLINDRICAL_SURFACE('',#{pl},{})", dec(&rad)));
        let (o1, o2) = (
            self.oriented(bottom_ring, false),
            self.oriented(top_ring, true),
        );
        let (b1, b2) = (self.bound(&[o1], true), self.bound(&[o2], false));
        self.finish_face(&[b1, b2], wall, false);
    }

    /// A flat elliptical face about `axis`, semi-axes `a` along the next axis and `b`.
    pub fn elliptical_disc(&mut self, centre: [&str; 3], axis: usize, a: &str, b: &str) {
        let c: P3 = centre.map(r);
        let ref_axis = (axis + 1) % 3;
        let first = Self::at(&c, ref_axis, &r(a));
        let v = self.vertex(&first);
        let pl = self.placement(&c, &Self::unit(axis), &Self::unit(ref_axis));
        let ellipse = self.w.add(format!("ELLIPSE('',#{pl},{a},{b})"));
        let e = self
            .w
            .add(format!("EDGE_CURVE('',#{v},#{v},#{ellipse},.T.)"));
        let o = self.oriented(e, true);
        let bound = self.bound(&[o], true);
        let pl = self.placement(&c, &Self::unit(axis), &Self::unit(ref_axis));
        let plane = self.w.add(format!("PLANE('',#{pl})"));
        self.finish_face(&[bound], plane, true);
    }

    /// A plane face with no boundary at all.
    pub fn unbounded_plane(&mut self, origin: [&str; 3], axis: usize) {
        let o: P3 = origin.map(r);
        let pl = self.placement(&o, &Self::unit(axis), &Self::unit((axis + 1) % 3));
        let plane = self.w.add(format!("PLANE('',#{pl})"));
        self.finish_face(&[], plane, true);
    }

    /// A whole sphere face.
    pub fn sphere(&mut self, centre: [&str; 3], radius: &str) {
        let c: P3 = centre.map(r);
        let pl = self.placement(&c, &Self::unit(2), &Self::unit(0));
        let s = self
            .w
            .add(format!("SPHERICAL_SURFACE('',#{pl},{})", radius));
        self.finish_face(&[], s, true);
    }

    /// A whole torus face about `axis`.
    pub fn torus(&mut self, centre: [&str; 3], axis: usize, major: &str, minor: &str) {
        let c: P3 = centre.map(r);
        let pl = self.placement(&c, &Self::unit(axis), &Self::unit((axis + 1) % 3));
        let s = self
            .w
            .add(format!("TOROIDAL_SURFACE('',#{pl},{major},{minor})"));
        self.finish_face(&[], s, true);
    }

    pub fn finish(mut self) -> String {
        let list = self
            .faces
            .iter()
            .map(|f| format!("#{f}"))
            .collect::<Vec<_>>()
            .join(",");
        let shell = self.w.add(format!("CLOSED_SHELL('',({list}))"));
        self.w.add(format!("MANIFOLD_SOLID_BREP('soup',#{shell})"));
        let mut out = String::from("ISO-10303-21;\nHEADER;\nENDSEC;\nDATA;\n");
        for (i, line) in self.w.lines.iter().enumerate() {
            out.push_str(&format!("#{}={};\n", i + 1, line));
        }
        let base = self.w.lines.len() as u64;
        out.push_str(&format!(
            "#{}=( GEOMETRIC_REPRESENTATION_CONTEXT(3) GLOBAL_UNIT_ASSIGNED_CONTEXT((#{})) REPRESENTATION_CONTEXT('',''));\n#{}=( LENGTH_UNIT() NAMED_UNIT(*) SI_UNIT(.MILLI.,.METRE.) );\nENDSEC;\nEND-ISO-10303-21;\n",
            base + 1,
            base + 2,
            base + 2
        ));
        out
    }
}

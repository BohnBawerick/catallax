//! Authored features bounded by complete circles: rims, planar caps and
//! cylindrical walls between two rims, all about axes along z. Each rim has
//! its one vertex at angle zero.
use catallax_numbers::{Ball, Rational};
use catallax_shapes::{Circle, Curve, Cylinder, Line, Plane, RPoint3, RVector3, Surface};
use catallax_topology::{
    Sense,
    numerical::{
        CopyTag, NumericalDefinition, NumericalEdge, NumericalFace, NumericalTrim, NumericalVertex,
    },
};

pub fn q(s: &str) -> Rational {
    Rational::parse(s).unwrap()
}

fn point(p: &[Rational; 3]) -> catallax_shapes::Point3 {
    RPoint3::new(p[0].clone(), p[1].clone(), p[2].clone()).to_real()
}

fn vector(p: [i64; 3]) -> catallax_shapes::Vector3 {
    RVector3::new(p[0].into(), p[1].into(), p[2].into()).to_real()
}

fn period() -> Rational {
    Ball::pi(128).mid() * q("2")
}

/// A complete circle in a plane of constant z.
#[derive(Clone)]
pub struct Rim {
    edge: usize,
    vertex: usize,
    centre: [Rational; 3],
    radius: Rational,
}

impl Rim {
    /// Where the rim starts and ends, in the x and y of its plane.
    fn uv(&self) -> [Rational; 2] {
        [&self.centre[0] + &self.radius, self.centre[1].clone()]
    }
}

/// A definition with no entities yet.
pub fn empty(name: &str) -> NumericalDefinition {
    NumericalDefinition {
        name: name.into(),
        vertices: Vec::new(),
        edges: Vec::new(),
        faces: Vec::new(),
    }
}

/// Adds the circle of `radius` about `centre`, keyed `key`, and its vertex,
/// keyed `key` with `-v` after it.
pub fn rim(d: &mut NumericalDefinition, key: &str, centre: [&str; 3], radius: &str) -> Rim {
    let centre = centre.map(q);
    let radius = q(radius);
    let vertex = d.vertices.len();
    d.vertices.push(NumericalVertex {
        key: format!("{key}-v"),
        position: [&centre[0] + &radius, centre[1].clone(), centre[2].clone()],
    });
    let edge = d.edges.len();
    d.edges.push(NumericalEdge {
        key: key.into(),
        curve: Curve::Circle(Circle {
            center: point(&centre),
            normal: vector([0, 0, 1]),
            radius_sq: (&radius * &radius).into(),
            ref_dir: vector([1, 0, 0]),
        }),
        vertices: [vertex, vertex],
        parameters: [q("0"), period()],
    });
    Rim {
        edge,
        vertex,
        centre,
        radius,
    }
}

fn circle_loop(rim: &Rim, sense: Sense) -> Vec<NumericalTrim> {
    vec![NumericalTrim {
        tag: CopyTag::Independent,
        edge: rim.edge,
        sense,
        uv: [rim.uv(), rim.uv()],
    }]
}

/// Adds `rim` as a hole of the planar face at index `face`, whose
/// parameters are x and y. `up` says whether the face's material normal is
/// plus z.
pub fn hole(d: &mut NumericalDefinition, face: usize, up: bool, rim: &Rim) {
    let sense = if up { Sense::Reversed } else { Sense::Same };
    d.faces[face].loops.push(circle_loop(rim, sense));
}

/// Adds the planar face in the plane of `outer`, bounded by it, with
/// `holes` cut out. `up` says whether its material normal is plus z.
pub fn cap(d: &mut NumericalDefinition, key: &str, up: bool, outer: &Rim, holes: &[&Rim]) {
    let origin = [q("0"), q("0"), outer.centre[2].clone()];
    let face = d.faces.len();
    d.faces.push(NumericalFace {
        key: key.into(),
        surface: Surface::Plane(Plane {
            origin: point(&origin),
            normal: vector([0, 0, 1]),
            ref_dir: vector([1, 0, 0]),
        }),
        sense: if up { Sense::Same } else { Sense::Reversed },
        loops: vec![circle_loop(
            outer,
            if up { Sense::Same } else { Sense::Reversed },
        )],
    });
    for rim in holes {
        hole(d, face, up, rim);
    }
}

/// Adds the cylindrical face between two rims of one radius about one axis,
/// `lower` under `upper`, and its seam, keyed `key` with `-seam` after it.
/// `outward` says whether the material is inside the cylinder (a boss) or
/// outside it (a bore).
pub fn wall(d: &mut NumericalDefinition, key: &str, lower: &Rim, upper: &Rim, outward: bool) {
    let (z0, z1) = (lower.centre[2].clone(), upper.centre[2].clone());
    let seam = d.edges.len();
    d.edges.push(NumericalEdge {
        key: format!("{key}-seam"),
        curve: Curve::Line(Line {
            origin: point(&[lower.uv()[0].clone(), lower.uv()[1].clone(), z0.clone()]),
            dir: vector([0, 0, 1]),
        }),
        vertices: [lower.vertex, upper.vertex],
        parameters: [q("0"), &z1 - &z0],
    });
    let turn = period();
    let zero = q("0");
    let trim =
        |edge: usize, sense: Sense, from: [&Rational; 2], to: [&Rational; 2]| NumericalTrim {
            tag: CopyTag::Independent,
            edge,
            sense,
            uv: [
                [from[0].clone(), from[1].clone()],
                [to[0].clone(), to[1].clone()],
            ],
        };
    let loops = if outward {
        vec![
            trim(lower.edge, Sense::Same, [&zero, &z0], [&turn, &z0]),
            trim(seam, Sense::Same, [&turn, &z0], [&turn, &z1]),
            trim(upper.edge, Sense::Reversed, [&turn, &z1], [&zero, &z1]),
            trim(seam, Sense::Reversed, [&zero, &z1], [&zero, &z0]),
        ]
    } else {
        vec![
            trim(lower.edge, Sense::Reversed, [&turn, &z0], [&zero, &z0]),
            trim(seam, Sense::Same, [&zero, &z0], [&zero, &z1]),
            trim(upper.edge, Sense::Same, [&zero, &z1], [&turn, &z1]),
            trim(seam, Sense::Reversed, [&turn, &z1], [&turn, &z0]),
        ]
    };
    d.faces.push(NumericalFace {
        key: key.into(),
        surface: Surface::Cylinder(Cylinder {
            origin: point(&[
                lower.centre[0].clone(),
                lower.centre[1].clone(),
                zero.clone(),
            ]),
            axis: vector([0, 0, 1]),
            radius_sq: (&lower.radius * &lower.radius).into(),
            ref_dir: vector([1, 0, 0]),
        }),
        sense: if outward {
            Sense::Same
        } else {
            Sense::Reversed
        },
        loops: vec![loops],
    });
}

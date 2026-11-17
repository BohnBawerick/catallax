//! Two closed boxes share one non-affinely parameterized straight edge.
//! Its authored pcurves refer only to the first box's distinct plane carriers.
use std::collections::BTreeMap;

use catallax_import::geometry_source::{AngleUnit, ParsedStep, Topology};
use catallax_import::tolerant::native::NativeSurface;
use catallax_numbers::Decider;

use super::common::{BoxFixture, boxes_step, r};

pub struct Fixture {
    pub source: String,
    pub roots: [u64; 2],
    pub first_surfaces: Vec<u64>,
}

pub fn uses(parsed: &ParsedStep, root: u64, dec: &mut Decider) -> BTreeMap<u64, Vec<u64>> {
    let Topology::Solid { outer, .. } = parsed.topology(root, dec).unwrap() else {
        unreachable!()
    };
    let Topology::Shell { faces, .. } = parsed.topology(outer, dec).unwrap() else {
        unreachable!()
    };
    let mut result = BTreeMap::<u64, Vec<u64>>::new();
    for face in faces {
        let Topology::Face {
            bounds,
            surface: Some(surface),
            ..
        } = parsed.topology(face, dec).unwrap()
        else {
            unreachable!()
        };
        for bound in bounds {
            let Topology::FaceBound { loop_id, .. } = parsed.topology(bound, dec).unwrap() else {
                unreachable!()
            };
            let Topology::EdgeLoop(edges) = parsed.topology(loop_id, dec).unwrap() else {
                unreachable!()
            };
            for edge in edges {
                let Topology::OrientedEdge(edge, _) = parsed.topology(edge, dec).unwrap() else {
                    unreachable!()
                };
                result.entry(edge).or_default().push(surface);
            }
        }
    }
    result
}

struct Append {
    next: u64,
    records: Vec<String>,
}
impl Append {
    fn add(&mut self, value: String) -> u64 {
        let id = self.next;
        self.next += 1;
        self.records.push(format!("#{id}={value};"));
        id
    }
    fn point(&mut self, coordinates: &[f64]) -> u64 {
        self.add(format!(
            "CARTESIAN_POINT('',({}))",
            coordinates
                .iter()
                .map(|value| format!("{value:.17e}"))
                .collect::<Vec<_>>()
                .join(",")
        ))
    }
    fn spline(&mut self, points: [u64; 3]) -> u64 {
        self.add(format!("B_SPLINE_CURVE_WITH_KNOTS('',2,(#{},#{},#{}),.UNSPECIFIED.,.F.,.F.,(3,3),(0.,1.),.UNSPECIFIED.)", points[0], points[1], points[2]))
    }
}

fn references(value: &str, replacements: &BTreeMap<u64, u64>) -> String {
    let mut result = String::new();
    let mut index = 0;
    while index < value.len() {
        if value.as_bytes()[index] == b'#' {
            let mut end = index + 1;
            while end < value.len() && value.as_bytes()[end].is_ascii_digit() {
                end += 1;
            }
            let id = value[index + 1..end].parse::<u64>().unwrap();
            result.push_str(&format!(
                "#{}",
                replacements.get(&id).copied().unwrap_or(id)
            ));
            index = end;
        } else {
            result.push(value.as_bytes()[index] as char);
            index += 1;
        }
    }
    result
}

pub fn boxes() -> Fixture {
    let source = boxes_step(&[
        BoxFixture::new(["2", "3", "4"]),
        BoxFixture::new(["2", "3", "4"]),
    ]);
    let mut dec = Decider::default();
    let parsed = ParsedStep::read(source.as_bytes(), &mut dec).unwrap();
    let roots = parsed
        .entities()
        .filter(|&id| matches!(parsed.topology(id, &mut dec), Ok(Topology::Solid { .. })))
        .collect::<Vec<_>>();
    let a = uses(&parsed, roots[0], &mut dec);
    let b = uses(&parsed, roots[1], &mut dec);
    let (&edge, surfaces) = a.first_key_value().unwrap();
    let (&other, _) = b.first_key_value().unwrap();
    let Topology::Edge {
        vertices,
        curve: Some(curve),
        ..
    } = parsed.topology(edge, &mut dec).unwrap()
    else {
        unreachable!()
    };
    let Topology::Edge {
        vertices: other_vertices,
        ..
    } = parsed.topology(other, &mut dec).unwrap()
    else {
        unreachable!()
    };
    let mut replacements = BTreeMap::from([(other, edge)]);
    for (before, after) in other_vertices.into_iter().zip(vertices) {
        replacements.insert(before, after);
    }
    let points = vertices.map(|vertex| {
        let Topology::Vertex(point) = parsed.topology(vertex, &mut dec).unwrap() else {
            unreachable!()
        };
        parsed
            .cartesian_point(point, &[r("1"), r("1"), r("1")], &mut dec)
            .unwrap()
            .map(|coordinate| coordinate.to_f64_nearest())
    });
    let middle =
        std::array::from_fn::<_, 3, _>(|axis| points[0][axis] * 0.75 + points[1][axis] * 0.25);
    let controls = [points[0], middle, points[1]];
    let mut append = Append {
        next: parsed.entities().last().unwrap() + 1,
        records: Vec::new(),
    };
    let context = append.add("(GEOMETRIC_REPRESENTATION_CONTEXT(2)PARAMETRIC_REPRESENTATION_CONTEXT()REPRESENTATION_CONTEXT('',''))".into());
    let mut pcurves = Vec::new();
    for &surface in surfaces {
        let carrier =
            NativeSurface::from_source(&parsed, surface, &r("1"), AngleUnit::Radians, &mut dec)
                .unwrap();
        let uv = controls.map(|point| {
            carrier
                .project(point, None, &mut dec)
                .unwrap()
                .location
                .parameters
        });
        let points = uv.map(|point| append.point(&point));
        let spline = append.spline(points);
        let representation = append.add(format!(
            "DEFINITIONAL_REPRESENTATION('',(#{spline}),#{context})"
        ));
        pcurves.push(append.add(format!("PCURVE('',#{surface},#{representation})")));
    }
    let points = controls.map(|point| append.point(&point));
    let basis = append.spline(points);
    let wrapper = format!(
        "SURFACE_CURVE('shared nonlinear parameter',#{basis},(#{},#{}),.CURVE_3D.)",
        pcurves[0], pcurves[1]
    );
    let source = source
        .lines()
        .map(|line| {
            let Some((left, right)) = line.split_once('=') else {
                return line.to_owned();
            };
            if left == format!("#{curve}") {
                format!("{left}={wrapper};")
            } else {
                format!("{left}={}", references(right, &replacements))
            }
        })
        .collect::<Vec<_>>()
        .join("\n");
    let source = source.replace(
        "ENDSEC;\nEND-ISO-10303-21;",
        &format!("{}\nENDSEC;\nEND-ISO-10303-21;", append.records.join("\n")),
    );
    Fixture {
        source,
        roots: [roots[0], roots[1]],
        first_surfaces: a.values().flatten().copied().collect(),
    }
}

//! Display translation of the draft's cached processor geometry.
use super::{Cell, Fast, Region};
use crate::unproved_wgsl::{number, unsupported, vector};
use catallax_numbers::Refusal;

fn conjunction(terms: impl Iterator<Item = String>) -> String {
    terms.collect::<Vec<_>>().join(" && ")
}

impl Region {
    pub(crate) fn wgsl(&self) -> Result<String, Refusal> {
        let mut s = String::from("var best=1e30; var inside=false;");
        match &self.fast {
            Fast::Oblique(fast) => return Ok(fast.wgsl()),
            Fast::Planar { vertices, planes } => {
                let halfspaces = |p: &str| {
                    conjunction(
                        planes
                            .iter()
                            .map(|(o, n)| format!("dot({p}-{}, {})<=0.0", vector(o), vector(n))),
                    )
                };
                match &self.cell {
                    Cell::Planar(_) => {
                        s.push_str(&format!("inside={};", halfspaces("p")));
                        for (o, n) in planes {
                            s.push_str(&format!(
                                "{{let n={}; let value=dot(p-{},n); let norm=dot(n,n); \
                                 let foot=p-n*(value/norm); if({}) {{best=min(best,value*value/norm);}}}}",
                                vector(n), vector(o), halfspaces("foot"),
                            ));
                        }
                    }
                    Cell::Polyhedral(cell) => {
                        let terms = cell
                            .terms
                            .iter()
                            .map(|term| {
                                format!(
                                    "({})",
                                    conjunction(term.sides.iter().map(|(i, same)| {
                                        let (o, n) = planes[*i];
                                        format!(
                                            "dot(p-{}, {})*{}<=0.0",
                                            vector(&o),
                                            vector(&n),
                                            if *same { "1.0" } else { "-1.0" }
                                        )
                                    }))
                                )
                            })
                            .collect::<Vec<_>>()
                            .join(" || ");
                        s.push_str(&format!("inside={terms};"));
                        for face in &cell.faces {
                            let (o, n) = planes[face.support];
                            let tiles = face
                                .tiles
                                .iter()
                                .map(|tile| {
                                    format!(
                                        "({})",
                                        conjunction(tile.iter().enumerate().map(|(i, a)| {
                                            let b = tile[(i + 1) % tile.len()];
                                            let edge = std::array::from_fn::<_, 3, _>(|j| {
                                                vertices[b][j] - vertices[*a][j]
                                            });
                                            format!(
                                                "dot(cross({},foot-{}),n)*{}>=0.0",
                                                vector(&edge),
                                                vector(&vertices[*a]),
                                                if face.same_sense { "1.0" } else { "-1.0" }
                                            )
                                        }))
                                    )
                                })
                                .collect::<Vec<_>>()
                                .join(" || ");
                            s.push_str(&format!(
                                "{{let n={}; let value=dot(p-{},n); let norm=dot(n,n); \
                                 let foot=p-n*(value/norm); if({tiles}) {{best=min(best,value*value/norm);}}}}",
                                vector(&n), vector(&o),
                            ));
                        }
                    }
                    _ => return Err(unsupported()),
                }
                for (a, b) in &self.edges {
                    let direction =
                        std::array::from_fn::<_, 3, _>(|i| vertices[*b][i] - vertices[*a][i]);
                    s.push_str(&format!(
                        "{{let a={}; let edge={}; let t=clamp(dot(p-a,edge)/dot(edge,edge),0.0,1.0); \
                         let delta=p-(a+t*edge); best=min(best,dot(delta,delta));}}",
                        vector(&vertices[*a]), vector(&direction),
                    ));
                }
            }
            Fast::Round {
                origin,
                axis,
                heights,
                radii,
            } => {
                let Cell::Round(cell) = &self.cell else {
                    return Err(unsupported());
                };
                s.push_str(&radial_coordinates(origin, axis, heights));
                s.push_str(&wall(*radii));
                for i in 0..2 {
                    if cell.caps[i].is_some() {
                        s.push_str(&format!(
                            "{{let r=max(radius-{},0.0); let z=h-{}; best=min(best,r*r+z*z);}}",
                            number(radii[i]),
                            number(heights[i]),
                        ));
                    }
                }
                s.push_str(&format!(
                    "inside=h>=lo && h<=hi && radius<={}+(h-lo)/dh*({}-{});",
                    number(radii[0]),
                    number(radii[1]),
                    number(radii[0]),
                ));
            }
            Fast::Annular {
                origin,
                axis,
                heights,
                radii,
            } => {
                s.push_str(&radial_coordinates(origin, axis, heights));
                for r in radii {
                    s.push_str(&wall(*r));
                }
                for i in 0..2 {
                    s.push_str(&format!(
                        "{{let r=max(max({}-radius,0.0),radius-{}); let z=h-{}; best=min(best,r*r+z*z);}}",
                        number(radii[1][i]), number(radii[0][i]), number(heights[i]),
                    ));
                }
                s.push_str(&format!(
                    "let t=(h-lo)/dh; let outer={}+t*({}-{}); let inner={}+t*({}-{}); \
                     inside=h>=lo && h<=hi && radius>=inner && radius<=outer;",
                    number(radii[0][0]),
                    number(radii[0][1]),
                    number(radii[0][0]),
                    number(radii[1][0]),
                    number(radii[1][1]),
                    number(radii[1][0]),
                ));
            }
        }
        s.push_str("let distance=sqrt(best); return select(distance,-distance,inside);");
        Ok(s)
    }
}

fn radial_coordinates(origin: &[f64; 3], axis: &[f64; 3], heights: &[f64; 2]) -> String {
    format!(
        "let delta=p-{}; let axis={}; let h=dot(delta,axis); \
         let radial=delta-axis*h; let radius=length(radial); let lo={}; let hi={}; let dh=hi-lo;",
        vector(origin),
        vector(axis),
        number(heights[0]),
        number(heights[1]),
    )
}

fn wall(radii: [f64; 2]) -> String {
    format!(
        "{{let rb={}; let rt={}; let dr=rt-rb; \
         let t=clamp(((radius-rb)*dr+(h-lo)*dh)/(dr*dr+dh*dh),0.0,1.0); \
         let r=radius-(rb+t*dr); let z=h-(lo+t*dh); best=min(best,r*r+z*z);}}",
        number(radii[0]),
        number(radii[1]),
    )
}

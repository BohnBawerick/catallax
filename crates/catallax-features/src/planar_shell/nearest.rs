//! Closed nearest-feature halfspaces intersected with every source material cell.

use catallax_numbers::{Answer, Decider, MarkLog, Refusal, Sign, current_operation_work};
use catallax_shapes::{Plane, Point3};
use catallax_topology::{BodyId, FaceId, Model};

use crate::{draft::arrangement::Skin, polyhedron::{ConvexCell, PlaneFace}, profile::sign};

struct SourceCell {
    geometry: ConvexCell,
    body: Option<BodyId>,
}
pub(super) struct Source {
    cells: Vec<SourceCell>,
}

fn value(plane: &Plane, point: &Point3) -> catallax_numbers::Real {
    point.minus(&plane.origin).dot(&plane.normal)
}
fn additional(cell: &ConvexCell, planes: &[Plane], dec: &mut Decider) -> Result<Option<Vec<PlaneFace>>, Refusal> {
    let mut cuts=Vec::new();
    for plane in planes {
        let mut negative=false;
        let mut positive=false;
        for point in cell.vertices() {
            match sign(&value(plane,point),dec)? {
                Sign::Negative => negative=true,
                Sign::Positive => positive=true,
                Sign::Zero => {}
            }
        }
        // No strict feasible side means no volume. A positive-radius curved
        // carrier can meet the remaining planar contact only below face dimension.
        if !negative { return Ok(None); }
        if positive { cuts.push(PlaneFace {plane:plane.clone(),origins:Vec::new()}); }
    }
    Ok(Some(cuts))
}
fn interior(cell: &ConvexCell, point: &Point3, dec: &mut Decider) -> Result<Option<bool>, Refusal> {
    let mut boundary=false;
    for face in &cell.faces {
        match sign(&value(&face.plane,point),dec)? {
            Sign::Positive => return Ok(Some(false)),
            Sign::Zero => boundary=true,
            Sign::Negative => {}
        }
    }
    Ok((!boundary).then_some(true))
}

impl Source {
    pub(super) fn new(skin: &Skin, dec: &mut Decider) -> Result<Self,Refusal> {
        let before=current_operation_work().used;
        let result:Result<Self,Refusal>=(||{
            let mut cells=Vec::new();
            for component in skin.components(dec)? {
                for term in &component.terms {
                    let mut faces=Vec::with_capacity(term.sides.len());
                    for &(index,same) in &term.sides {
                        let mut plane=component.planes[index].clone();
                        if !same { plane.normal=plane.normal.scaled(&catallax_numbers::Real::from_i64(-1)); }
                        faces.push(PlaneFace {plane,origins:Vec::new()});
                    }
                    let geometry=ConvexCell::new(faces,dec)?;
                    if !geometry.is_empty() { cells.push(SourceCell {geometry,body:None}); }
                }
            }
            if cells.is_empty() {
                return Err(Refusal::new("distance.unsupported-envelope", "the nonempty source has no full-dimensional winding cell"));
            }
            Ok(Self {cells})
        })();
        eprintln!("NEAREST_SOURCE phase=partition work={} cumulative={} cells={} outcome={:?}",current_operation_work().used-before,current_operation_work().used,result.as_ref().map(|source|source.cells.len()).unwrap_or(0),result.as_ref().map(|_|"ok").map_err(|e|e.code));
        result
    }

    pub(super) fn clip(
        &mut self,
        model: &mut Model,
        kind: &str,
        planes: &[Plane],
        make: impl FnOnce(&mut Model,&mut Decider)->Result<FaceId,Refusal>,
        dec: &mut Decider,
    ) -> Result<Answer<Vec<FaceId>>,Refusal> {
        let before=current_operation_work().used;
        let mut make=Some(make);
        let mut patch=None;
        let result:Result<Answer<Vec<FaceId>>,Refusal>=(||{
            let mut retained=Vec::new();
            let mut marks=MarkLog::new();
            for source in &mut self.cells {
                let Some(extra)=additional(&source.geometry,planes,dec)? else {continue;};
                let clipped;
                let (region,body)=if extra.is_empty() {
                    if source.body.is_none() {source.body=Some(source.geometry.build(model));}
                    (&source.geometry,source.body.expect("initialized source cell"))
                } else {
                    let mut faces=source.geometry.faces.clone();
                    faces.extend(extra);
                    clipped=ConvexCell::new(faces,dec)?;
                    if clipped.is_empty() {continue;}
                    let body=clipped.build(model);
                    (&clipped,body)
                };
                let face=if let Some(face)=patch {face} else {
                    let face=make.take().expect("primitive is built once")(model,dec)?;
                    patch=Some(face);
                    face
                };
                let answer=catallax_boolean::trim_patch(model,body,face,|point,dec|interior(region,point,dec),dec)
                    .map_err(|error|error.refusal)?;
                retained.extend(answer.value);
                marks.extend(&answer.marks);
            }
            Ok(Answer::with_marks(retained,marks))
        })();
        eprintln!("NEAREST_DOMAIN kind={kind} work={} cumulative={} generated={} retained={} outcome={:?}",current_operation_work().used-before,current_operation_work().used,patch.is_some(),result.as_ref().map(|answer|answer.value.len()).unwrap_or(0),result.as_ref().map(|_|"ok").map_err(|e|e.code));
        result
    }
}

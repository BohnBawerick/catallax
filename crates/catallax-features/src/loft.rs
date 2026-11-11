//! Definition-native rational lofts with complete embedding certificates.

use crate::Profile;
use catallax_numbers::{
    Decider, Outcome, Rational, Real, RealPoly, Refusal, Sign, StrictSign, ZeroProof,
};
use catallax_recipe::{LoftParameterization, LoftSection, LoftShape, Operation, Recipe, StepId};
use catallax_topology::{BodyId, Model};

mod charts;
mod math;
mod proof;
pub(crate) mod topology;

fn unproven(detail: &str) -> Refusal {
    Refusal::new("loft.admissibility-unproven", detail)
}

fn sign(value: &Real, dec: &mut Decider) -> Result<Sign, Refusal> {
    match dec.sign(value) {
        Outcome::Proven(StrictSign::Positive) => Ok(Sign::Positive),
        Outcome::Proven(StrictSign::Negative) => Ok(Sign::Negative),
        Outcome::ProvenCoincident(proof) if proof.proof != ZeroProof::SeparationBound => {
            Ok(Sign::Zero)
        }
        Outcome::ProvenCoincident(_) => Err(unproven(
            "loft authority requires exact equality, not bounded coincidence",
        )),
        Outcome::NotExact(error) => Err(Refusal::from(error)),
        Outcome::Refused(error) => Err(error),
    }
}

pub(crate) fn build(
    model: &mut Model,
    recipe: &Recipe,
    sections: &[LoftSection],
    shape: LoftShape,
    step: &StepId,
    profiles: &im::OrdMap<StepId, Profile>,
    dec: &mut Decider,
) -> Result<(BodyId, Rational), Refusal> {
    if !matches!(
        recipe.loft_parameterization(),
        Ok(LoftParameterization::RationalCharts)
    ) {
        return Err(Refusal::new(
            "loft.legacy-parameterization",
            "legacy arc-length lofts cannot be reinterpreted as rational charts",
        ));
    }
    if sections.len() < 2 {
        return Err(Refusal::new(
            "loft.sections",
            "a loft needs at least two sections",
        ));
    }
    let validated = sections
        .iter()
        .map(|section| {
            profiles
                .get(&section.profile)
                .ok_or_else(|| unproven("the section has no validated profile"))
        })
        .collect::<Result<Vec<_>, Refusal>>()?;
    if validated.iter().all(|p| p.is_empty()) {
        return Ok((catallax_topology::empty(model), Rational::zero()));
    }
    if validated.iter().any(|p| p.is_empty()) {
        return Err(Refusal::new(
            "loft.incompatible-sections",
            "a nonempty loop cannot correspond to an empty section",
        ));
    }
    let mut authored = Vec::with_capacity(sections.len());
    let mut infos: Vec<topology::LoopInfo> = Vec::new();
    for (station, (section, profile)) in sections.iter().zip(&validated).enumerate() {
        let sketch = recipe
            .steps
            .iter()
            .find(|step| step.id == section.profile)
            .ok_or_else(|| unproven("the retained section definition is absent"))?;
        let Operation::Sketch { loops, .. } = &sketch.op else {
            return Err(unproven("a loft section must retain its authored sketch"));
        };
        let mut boundary_indices = Vec::new();
        let mut section_charts = Vec::new();
        let mut section_infos = Vec::new();
        for start in &section.start {
            let lp = loops
                .iter()
                .find(|lp| lp.0.iter().any(|s| s.id() == start))
                .ok_or_else(|| unproven("a designated start does not belong to a section loop"))?;
            let boundary = profile
                .loops
                .iter()
                .position(|b| b.segments.iter().any(|s| s.ids.contains(start)))
                .ok_or_else(|| unproven("a designated start has no validated loop"))?;
            let (charts, ids, material_forward) = charts::build(profile, lp, start, dec)?;
            boundary_indices.push(boundary);
            section_infos.push(topology::LoopInfo {
                ids,
                depth: profile.loops[boundary].depth,
                parent: profile.loops[boundary].parent,
                material_forward,
            });
            section_charts.push(charts);
        }
        for info in &mut section_infos {
            info.parent = info.parent.map(|parent| {
                boundary_indices
                    .iter()
                    .position(|index| *index == parent)
                    .expect("complete loop correspondence")
            });
        }
        if station == 0 {
            infos = section_infos;
        } else if section_infos.len() != infos.len()
            || section_infos
                .iter()
                .zip(&infos)
                .any(|(a, b)| a.depth != b.depth || a.parent != b.parent)
        {
            return Err(Refusal::new(
                "loft.incompatible-nesting",
                "corresponding loops do not have the same material nesting",
            ));
        }
        authored.push(section_charts);
    }
    let mut common = Vec::with_capacity(infos.len());
    for loop_index in 0..infos.len() {
        let mut cuts: Vec<_> = authored
            .iter()
            .flat_map(|station| {
                station[loop_index]
                    .iter()
                    .flat_map(|chart| chart.bounds.iter().cloned())
            })
            .collect();
        cuts.sort();
        cuts.dedup();
        let mut intervals = Vec::new();
        for cut in cuts.windows(2) {
            let mut inputs = Vec::with_capacity(sections.len());
            for station in &authored {
                let chart = station[loop_index]
                    .iter()
                    .find(|chart| chart.bounds[0] <= cut[0] && chart.bounds[1] >= cut[1])
                    .ok_or_else(|| unproven("the authored charts leave a common parameter gap"))?;
                inputs.push(chart.curve.span_bezier_controls_on(
                    chart.curve.degree(),
                    [&cut[0], &cut[1]],
                    dec,
                )?);
            }
            math::rational_interval(&mut inputs, dec)?;
            intervals.push(math::common([cut[0].clone(), cut[1].clone()], inputs, dec)?);
        }
        common.push(intervals);
    }
    if let Some(vector) = proof::prism(&common, validated[0], shape, dec)? {
        let bits = dec.policy().max_precision_bits.min(128);
        let mut radius = Rational::zero();
        for point in common.iter().flatten().flat_map(|part| &part.sections[0]) {
            let mut extent = Rational::zero();
            for coordinate in &point[..3] {
                let bound = dec.enclose_budgeted_at(&(coordinate / &point[3]), bits)?;
                extent = extent + bound.lo().abs().max(bound.hi().abs());
            }
            radius = radius.max(extent);
        }
        radius = radius + vector.x.abs() + vector.y.abs() + vector.z.abs();
        let body = crate::extrude::extrude_vector(
            model,
            validated[0],
            step,
            &vector,
            [
                Some(catallax_recipe::Role::Start),
                Some(catallax_recipe::Role::End),
            ],
            None,
            dec,
        )?;
        return Ok((body, radius));
    }
    let sources = common
        .iter()
        .map(|parts| math::surface(parts, shape, dec))
        .collect::<Result<Vec<_>, Refusal>>()?;
    let authority = proof::certify(&common, &sources, &validated, shape, dec)?;
    let proof::Proof::Solid {
        positive,
        end_positive,
    } = authority
    else {
        return Ok((catallax_topology::empty(model), Rational::zero()));
    };
    let bits = dec.policy().max_precision_bits.min(128);
    let mut radius = Rational::zero();
    for source in &sources {
        for point in source.controls() {
            let mut extent = Rational::zero();
            for coordinate in [&point.x, &point.y, &point.z] {
                let bound = dec.enclose_budgeted_at(coordinate, bits)?;
                extent = extent + bound.lo().abs().max(bound.hi().abs());
            }
            radius = radius.max(extent);
        }
    }
    let atlases = sources
        .into_iter()
        .map(|source| catallax_shapes::RationalAtlas::partition(source, true, dec))
        .collect::<Result<Vec<_>, Refusal>>()?;
    let body = topology::assemble(
        model,
        atlases,
        &infos,
        &validated,
        [positive, end_positive],
        step,
        dec,
    )?;
    Ok((body, radius))
}

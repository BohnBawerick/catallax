//! Capture the existing scan tests' answers for the Phase 4.4 differential.

use catallax_field::reconstruction::scan::{ScanFit, ScanFitOptions, fit_scan};
use catallax_numbers::{Decider, Refusal, operation_budget_active, with_operation_budget};
use catallax_recipe::{Operation, Recipe};

pub fn fit_recorded(
    step: &Operation,
    bytes: &[u8],
    options: &ScanFitOptions,
    dec: &mut Decider,
) -> Result<ScanFit, Refusal> {
    if std::env::var_os("CATALLAX_IMPORTS_DIFFERENTIAL").is_none() {
        return fit_scan(step, bytes, options, dec);
    }
    let (answer, work) = if operation_budget_active() {
        (fit_scan(step, bytes, options, dec), None)
    } else {
        let (answer, work) = with_operation_budget(dec.policy().max_work, || {
            fit_scan(step, bytes, options, dec)
        });
        (answer, Some(work))
    };
    let snapshot = answer.as_ref().map(|fit| {
        (
            &fit.sha256,
            &fit.samples,
            fit.sigma_mm.to_bits(),
            &fit.surfaces,
            &fit.assignment,
            &fit.unexplained,
            &fit.refusals,
            fit.recipe.as_ref().map(Recipe::to_text),
            fit.part.as_ref().map(|part| {
                (
                    part.recipe.to_text(),
                    part.brep.model.check(part.brep.result),
                    part.field
                        .faces()
                        .iter()
                        .map(|face| &face.origin)
                        .collect::<Vec<_>>(),
                    &part.deviation,
                )
            }),
        )
    });
    let hash =
        catallax_import::sha256_hex(format!("{snapshot:?} {:?} {work:?}", dec.stats()).as_bytes());
    println!(
        "imports-differential scan {} {} {hash}",
        catallax_import::sha256_hex(bytes),
        options.workers
    );
    answer
}

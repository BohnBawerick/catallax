use crate::{Plan, PreviewHardware, UnprovedValue, Workspace};
use catallax_numbers::Refusal;

/// Reusable storage for threaded, unproved processor preview values.
/// Results have input order and carry no proof or geometric error bound.
#[derive(Default, Debug)]
pub struct UnprovedBatch {
    values: Vec<UnprovedValue>,
}

impl Plan {
    /// Evaluates unproved processor preview points on all available field workers.
    /// Shares the pool used by deterministic field batches and meshing.
    pub fn unproved_preview_batch_default<'a>(
        &self,
        points: &[[f64; 3]],
        storage: &'a mut UnprovedBatch,
    ) -> Result<&'a [UnprovedValue], Refusal> {
        self.unproved_preview_batch(points, crate::workers::worker_count(), storage)
    }

    /// Evaluates unproved preview points on at most `workers` processor threads.
    /// The immutable plan is shared; each worker owns its scratch registers.
    /// Values are bit-identical to [`Self::unproved_preview`] on this machine,
    /// independent of worker count. Zero workers is a typed refusal, even for
    /// empty input. Any failure clears all results; the first failing input
    /// input wins over later evaluation failures, regardless of scheduling.
    ///
    /// Like the scalar preview, this does no proof work and creates no numeric
    /// allowance or proof cache. It cannot feed exports or proven measurements.
    pub fn unproved_preview_batch<'a>(
        &self,
        points: &[[f64; 3]],
        workers: usize,
        storage: &'a mut UnprovedBatch,
    ) -> Result<&'a [UnprovedValue], Refusal> {
        storage.values.clear();
        if workers == 0 {
            return Err(Refusal::new(
                "field.preview-workers",
                "unproved preview requires at least one worker",
            ));
        }
        if points.is_empty() {
            return Ok(&storage.values);
        }
        let results = crate::workers::with_workers(workers, || {
            crate::workers::map_init(points, Workspace::default, |workspace, &point| {
                self.tentative::<crate::metric::Processor>(point, workspace)
            })
        })?;
        for result in results {
            match result {
                Ok(value) => storage.values.push(UnprovedValue {
                    value,
                    true_distance_formula: self.true_distance_formula(),
                    hardware: PreviewHardware::Processor,
                }),
                Err(error) => {
                    storage.values.clear();
                    return Err(error);
                }
            }
        }
        Ok(&storage.values)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use catallax_recipe::Recipe;

    #[test]
    fn refusal_discards_previous_and_partial_results() {
        let plan = Plan::compile(
            &Recipe::parse(include_str!("../../../parts/P3-filleted-boss.json")).unwrap(),
        )
        .unwrap();
        let mut storage = UnprovedBatch::default();
        for workers in [1, 2, 4, 8] {
            plan.unproved_preview_batch(&[[1.0; 3]; 17], workers, &mut storage)
                .unwrap();
            // Earlier range failure must win over a later non-finite input.
            let mut points = [[1.0; 3]; 17];
            points[1] = [f64::MAX; 3];
            points[16] = [f64::NAN; 3];
            let scalar = plan
                .unproved_preview(points[1], &mut Workspace::default())
                .unwrap_err();
            let error = plan
                .unproved_preview_batch(&points, workers, &mut storage)
                .unwrap_err();
            assert_eq!(error.code, scalar.code);
            assert_eq!(error.detail, scalar.detail);
            assert!(storage.values.is_empty());
            plan.unproved_preview_batch(&[[1.0; 3]], workers, &mut storage)
                .unwrap();
            assert_eq!(storage.values.len(), 1);
        }
        assert_eq!(
            plan.unproved_preview_batch(&[], 0, &mut storage)
                .unwrap_err()
                .code,
            "field.preview-workers"
        );
        assert!(storage.values.is_empty());
    }
}

use crate::{Plan, PreviewHardware, UnprovedValue, Workspace};
use catallax_numbers::Refusal;

/// Reusable storage for threaded, unproved processor preview values.
/// Results have input order and carry no proof or geometric error bound.
#[derive(Default, Debug)]
pub struct UnprovedBatch {
    workspaces: Vec<Workspace>,
    values: Vec<UnprovedValue>,
}

impl Plan {
    /// Evaluates unproved preview points on at most `workers` processor threads.
    /// The immutable plan is shared; each worker owns its scratch registers.
    /// Values are bit-identical to [`Self::unproved_preview`] on this machine,
    /// independent of worker count. Zero workers is a typed refusal, even for
    /// empty input. Any failure clears all results; the first failing input
    /// chunk wins over later evaluation failures, regardless of scheduling.
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
        let workers = workers.min(points.len());
        storage.workspaces.resize_with(workers, Workspace::default);
        storage.values.resize(
            points.len(),
            UnprovedValue {
                value: 0.0,
                true_distance_formula: self.true_distance_formula(),
                hardware: PreviewHardware::Processor,
            },
        );
        let evaluate = |points: &[[f64; 3]], values: &mut [UnprovedValue], w: &mut Workspace| {
            for (&point, value) in points.iter().zip(values) {
                value.value = self.tentative::<crate::metric::Processor>(point, w)?;
            }
            Ok(())
        };
        let result = if workers == 1 {
            evaluate(points, &mut storage.values, &mut storage.workspaces[0])
        } else {
            let chunk = points.len().div_ceil(workers);
            std::thread::scope(|scope| {
                let mut handles = Vec::with_capacity(workers);
                let mut spawn_error = None;
                for ((points, values), workspace) in points
                    .chunks(chunk)
                    .zip(storage.values.chunks_mut(chunk))
                    .zip(&mut storage.workspaces)
                {
                    match std::thread::Builder::new()
                        .spawn_scoped(scope, move || evaluate(points, values, workspace))
                    {
                        Ok(handle) => handles.push(handle),
                        Err(error) => {
                            spawn_error = Some(Refusal::new(
                                "field.preview-thread",
                                format!("could not start unproved preview worker: {error}"),
                            ));
                            break;
                        }
                    }
                }
                let mut result = Ok(());
                // Join every worker before returning, including after a refusal.
                for handle in handles {
                    let worker = handle.join().unwrap_or_else(|_| {
                        Err(Refusal::new(
                            "field.preview-thread",
                            "unproved preview worker panicked",
                        ))
                    });
                    result = result.and(worker);
                }
                result.and(spawn_error.map_or(Ok(()), Err))
            })
        };
        if let Err(error) = result {
            storage.values.clear();
            return Err(error);
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

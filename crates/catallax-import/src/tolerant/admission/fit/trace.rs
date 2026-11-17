use super::*;

pub(super) mod planar;

pub(super) struct Fit<'a> {
    pub curve: &'a PhysicalCurve,
    pub surface: &'a NativeSurface,
    pub precision: f64,
    pub maximum_residual: f64,
    pub samples: usize,
    pub converged: usize,
}

impl Fit<'_> {
    fn interval(
        &mut self,
        left: &UvSample,
        right: &UvSample,
        depth: u8,
        output: &mut Vec<UvSample>,
        dec: &mut Decider,
    ) -> Result<(), Error> {
        dec.charge_work(1)?;
        if output.len() >= 65_536 {
            return Err(invalid(
                "import.tolerant-trace-resolution",
                "native UV reconstruction exceeds its finite sample allowance",
            ));
        }
        let parameter = left.curve.parameters[0] * 0.5 + right.curve.parameters[0] * 0.5;
        let (location, target) = self.curve.sample_at(parameter, left.curve.chart, dec)?;
        let uv = std::array::from_fn(|axis| {
            left.surface.parameters[axis] * 0.5 + right.surface.parameters[axis] * 0.5
        });
        let chord_location = self.surface.location(uv)?;
        let chord = self.surface.sample(chord_location, dec)?;
        let projection = self
            .surface
            .project(target.point_mm, Some(chord_location), dec)?;
        let uv = self
            .surface
            .align_parameters(projection.location.parameters, Some(uv));
        let surface_location = self.surface.location(uv)?;
        let projected = self.surface.sample(surface_location, dec)?;
        let error = norm(subtract(chord.point_mm, projected.point_mm))
            + chord.uncertainty_mm
            + projected.uncertainty_mm;
        let residual = norm(subtract(chord.point_mm, target.point_mm))
            + chord.uncertainty_mm
            + target.uncertainty_mm;
        self.samples += 1;
        self.converged += usize::from(projection.converged);
        let floor = 64.0 * f64::EPSILON * norm(projected.point_mm).max(1.0);
        if error <= self.precision.max(floor) {
            self.maximum_residual = self.maximum_residual.max(residual);
            output.push(right.clone());
            return Ok(());
        }
        if depth >= 24
            || output.len() >= 65_536
            || parameter <= left.curve.parameters[0].min(right.curve.parameters[0])
            || parameter >= left.curve.parameters[0].max(right.curve.parameters[0])
        {
            return Err(invalid(
                "import.tolerant-trace-resolution",
                &format!(
                    "native UV reconstruction did not reach its stated precision at depth {depth}: chord deviation {error} mm, source residual {residual} mm"
                ),
            ));
        }
        let middle = UvSample {
            fraction: (left.fraction + right.fraction) * 0.5,
            curve: location,
            surface: surface_location,
            residual_mm: projection.residual_mm,
            converged: projection.converged,
        };
        self.interval(left, &middle, depth + 1, output, dec)?;
        self.interval(&middle, right, depth + 1, output, dec)
    }

    pub fn rebuild(&mut self, samples: &mut Vec<UvSample>, dec: &mut Decider) -> Result<(), Error> {
        let mut rebuilt = Vec::with_capacity(samples.len());
        rebuilt.push(samples[0].clone());
        for pair in samples.windows(2) {
            self.interval(&pair[0], &pair[1], 0, &mut rebuilt, dec)?;
        }
        *samples = rebuilt;
        Ok(())
    }
}

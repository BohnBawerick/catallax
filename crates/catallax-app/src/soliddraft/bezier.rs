use super::*;

pub(super) fn convert(input: &[InputSegment], li: usize, step: &str) -> Result<Loop> {
    if input
        .iter()
        .any(|s| !matches!(s, InputSegment::Line { .. } | InputSegment::Bezier { .. }))
    {
        return Err(not_yet(
            step,
            "soliddraft.bezier-mix",
            "Bezier loops currently accept lines and polynomial spans",
        ));
    }
    let mut segments = Vec::with_capacity(input.len());
    for (i, segment) in input.iter().enumerate() {
        let from = match segment {
            InputSegment::Line { from, .. } | InputSegment::Bezier { from, .. } => point(*from)?,
            _ => unreachable!("checked segment family"),
        };
        let previous = end(&input[(i + input.len() - 1) % input.len()])?;
        let to = end(segment)?;
        if !close(&previous, &from) {
            return Err(invalid(Some(step), "open profile loop"));
        }
        if close(&from, &to) {
            return Err(invalid(Some(step), "degenerate profile segment"));
        }
        let id = ElementId::new(format!("l{li}s{i}")).unwrap();
        segments.push(match segment {
            InputSegment::Line { .. } => Segment::Line { id, to },
            InputSegment::Bezier { controls, .. } => {
                if !(1..=2).contains(&controls.len()) {
                    return Err(invalid(
                        Some(step),
                        "a Bezier span needs one or two interior controls",
                    ));
                }
                Segment::Bezier {
                    id,
                    controls: controls.iter().map(|p| point(*p)).collect::<Result<_>>()?,
                    to,
                }
            }
            _ => unreachable!("checked segment family"),
        });
    }
    Ok(Loop(segments))
}

//! Versioned numerical records; v02 has no sampled-carrier declarations and
//! v03 no tolerance report.
use super::{
    BodyId, Model, NumericalBody, NumericalDefinition, NumericalHistory, NumericalNativeTrim,
    NumericalOptions, ToleranceReport,
};
use catallax_numbers::{
    NumericalEvidence, Refusal,
    transfer::{Codec, Reader, Registry, Writer},
};
use std::collections::BTreeSet;

const BODY: &[u8; 8] = b"CXNUMB04";
const HISTORY: &[u8; 8] = b"CXNUMH03";
fn version(
    reader: &mut Reader<'_>,
    current: &[u8; 8],
    previous: &[u8; 8],
    kind: &str,
) -> Result<bool, Refusal> {
    let tag = reader.take(8)?;
    if tag == current {
        Ok(true)
    } else if tag == previous {
        Ok(false)
    } else {
        Err(Refusal::new(
            "handover.version-mismatch",
            format!("unsupported {kind} record layout"),
        ))
    }
}
impl Codec for NumericalBody {
    fn write(&self, w: &mut Writer<'_>) -> Result<(), Refusal> {
        w.raw(BODY)?;
        self.source.write(w)?;
        self.source_native_trims.write(w)?;
        self.source_sampled_carriers.write(w)?;
        self.sampled_carriers.write(w)?;
        self.native_trims.write(w)?;
        self.history.write(w)?;
        self.definition.write(w)?;
        self.options.write(w)?;
        self.model.write(w)?;
        self.body.write(w)?;
        self.evidence.write(w)?;
        self.tolerances.write(w)
    }
    fn read(r: &mut Reader<'_>) -> Result<Self, Refusal> {
        // v04 adds the tolerance report. A v03 body reads without one: its
        // report is empty, and a rebuild makes it again.
        let tag = r.take(8)?;
        let (carriers, report) = match tag {
            t if t == BODY => (true, true),
            t if t == b"CXNUMB03" => (true, false),
            t if t == b"CXNUMB02" => (false, false),
            _ => {
                return Err(Refusal::new(
                    "handover.version-mismatch",
                    "unsupported numerical body record layout",
                ));
            }
        };
        Ok(Self {
            source: Codec::read(r)?,
            source_native_trims: Codec::read(r)?,
            source_sampled_carriers: if carriers {
                Codec::read(r)?
            } else {
                BTreeSet::new()
            },
            sampled_carriers: if carriers {
                Codec::read(r)?
            } else {
                BTreeSet::new()
            },
            native_trims: Codec::read(r)?,
            history: Codec::read(r)?,
            definition: Codec::read(r)?,
            options: Codec::read(r)?,
            model: Codec::read(r)?,
            body: Codec::read(r)?,
            evidence: Codec::read(r)?,
            tolerances: if report {
                Codec::read(r)?
            } else {
                ToleranceReport::default()
            },
            // The rows of the gap rule for imports describe one admission
            // and are not stored.
            import_gaps: Vec::new(),
        })
        .map(|mut body: NumericalBody| {
            // The options skip the vertex cap; the report keeps it.
            body.options.vertex_tolerance_mm = body.tolerances.vertex_cap_mm.clone();
            body
        })
    }
    fn register(r: &mut Registry) {
        r.add::<NumericalDefinition>();
        r.add::<Vec<NumericalNativeTrim>>();
        r.add::<BTreeSet<String>>();
        r.add::<Vec<NumericalHistory>>();
        r.add::<NumericalOptions>();
        r.add::<Model>();
        r.add::<BodyId>();
        r.add::<NumericalEvidence>();
        r.add::<ToleranceReport>();
    }
}
impl Codec for NumericalHistory {
    fn write(&self, w: &mut Writer<'_>) -> Result<(), Refusal> {
        w.raw(HISTORY)?;
        self.definition.write(w)?;
        self.native_trims.write(w)?;
        self.sampled_carriers.write(w)?;
        self.options.write(w)?;
        self.evidence.write(w)
    }
    fn read(r: &mut Reader<'_>) -> Result<Self, Refusal> {
        let carriers = version(r, HISTORY, b"CXNUMH02", "numerical history")?;
        Ok(Self {
            definition: Codec::read(r)?,
            native_trims: Codec::read(r)?,
            sampled_carriers: if carriers {
                Codec::read(r)?
            } else {
                BTreeSet::new()
            },
            options: Codec::read(r)?,
            evidence: Codec::read(r)?,
        })
    }
    fn register(r: &mut Registry) {
        r.add::<NumericalDefinition>();
        r.add::<Vec<NumericalNativeTrim>>();
        r.add::<BTreeSet<String>>();
        r.add::<NumericalOptions>();
        r.add::<NumericalEvidence>();
    }
}

//! Per-record layout tags reject earlier positional numerical payloads.
use super::{
    BodyId, Model, NumericalBody, NumericalDefinition, NumericalHistory, NumericalNativeTrim,
    NumericalOptions,
};
use catallax_numbers::{
    NumericalEvidence, Refusal,
    transfer::{Codec, Reader, Registry, Writer},
};
use std::collections::BTreeSet;

const BODY: &[u8; 8] = b"CXNUMB03";
const HISTORY: &[u8; 8] = b"CXNUMH03";
fn version(reader: &mut Reader<'_>, expected: &[u8; 8], kind: &str) -> Result<(), Refusal> {
    if reader.take(8)? != expected {
        return Err(Refusal::new(
            "handover.version-mismatch",
            format!("unsupported {kind} record layout"),
        ));
    }
    Ok(())
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
        self.evidence.write(w)
    }
    fn read(r: &mut Reader<'_>) -> Result<Self, Refusal> {
        version(r, BODY, "numerical body")?;
        Ok(Self {
            source: Codec::read(r)?,
            source_native_trims: Codec::read(r)?,
            source_sampled_carriers: Codec::read(r)?,
            sampled_carriers: Codec::read(r)?,
            native_trims: Codec::read(r)?,
            history: Codec::read(r)?,
            definition: Codec::read(r)?,
            options: Codec::read(r)?,
            model: Codec::read(r)?,
            body: Codec::read(r)?,
            evidence: Codec::read(r)?,
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
        version(r, HISTORY, "numerical history")?;
        Ok(Self {
            definition: Codec::read(r)?,
            native_trims: Codec::read(r)?,
            sampled_carriers: Codec::read(r)?,
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

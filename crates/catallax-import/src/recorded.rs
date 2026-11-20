//! The tolerance a Catallax numerical export records in its file.
//!
//! Check once (roadmap item 3.0A, step 5; Checkpoint 2 decisions, section
//! 4.1 item 5): "Reopening trusts the tolerance recorded at export." A
//! numerical STEP export stamps its header's `FILE_DESCRIPTION` with
//! `catallax numerical approximate` and a list of `key=value` fields. This
//! module reads the tolerance out of that stamp. Nothing binds the stamp to
//! the geometry under it: the file is trusted as Parasolid trusts the edge
//! and vertex tolerances of an XT file.

use catallax_numbers::{Rational, Sign};

/// The stamp every numerical export begins its description with.
const STAMP: &str = "catallax numerical approximate";

/// Which field of the stamp the tolerance was read from.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum RecordedField {
    /// `recorded_tolerance_mm`: the body's largest entity tolerance plus the
    /// rounding error of the written numbers.
    Recorded,
    /// `source_tolerance_mm`, in a file written before check once. That
    /// export admitted the written values under this limit, so every gap of
    /// the file was measured at or under it.
    AdmittedLimit,
}

/// The tolerance a numerical export recorded.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RecordedTolerance {
    /// In millimetres as a reader takes the file's lengths: a scaled
    /// export's stamp is in the model's millimetres, and its reader takes
    /// every length that many times larger.
    pub tolerance_mm: Rational,
    /// The field it was read from.
    pub field: RecordedField,
}

/// The tolerance recorded in a file's header, or none: a file another
/// system wrote, or a stamp without a positive tolerance. Such a file is
/// measured as every import is (rule T3).
///
/// Only the header is read. Every caller here has charged the file's bytes
/// before it asks.
pub fn recorded_tolerance(bytes: &[u8]) -> Option<RecordedTolerance> {
    let description = crate::part21::header_description(bytes);
    let stamp = description.first()?;
    let mut fields = stamp.split(';');
    if fields.next()? != STAMP {
        return None;
    }
    // The body's name comes first and is free text. The writer's own fields
    // follow it, so the last field of a name is the writer's.
    let fields: Vec<&str> = fields.collect();
    let last = |key: &str| {
        fields
            .iter()
            .rev()
            .find_map(|field| field.strip_prefix(key)?.strip_prefix('='))
    };
    let (text, field) = match last("recorded_tolerance_mm") {
        Some(text) => (text, RecordedField::Recorded),
        None => (last("source_tolerance_mm")?, RecordedField::AdmittedLimit),
    };
    let written = Rational::parse(text).ok()?;
    let scale = match last("scaled_export") {
        Some(factor) => Rational::parse(factor.strip_prefix('x')?).ok()?,
        None => Rational::one(),
    };
    let tolerance_mm = written * scale;
    (tolerance_mm.sign() == Sign::Positive).then_some(RecordedTolerance {
        tolerance_mm,
        field,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn file(description: &str) -> Vec<u8> {
        format!(
            "ISO-10303-21;\nHEADER;\nFILE_DESCRIPTION(({description}),'2;1');\nFILE_NAME('a.step','',(''),(''),'','','');\nFILE_SCHEMA(('X'));\nENDSEC;\nDATA;\n#1=CARTESIAN_POINT('',(0.,0.,0.));\nENDSEC;\nEND-ISO-10303-21;\n"
        )
        .into_bytes()
    }

    fn q(text: &str) -> Rational {
        Rational::parse(text).unwrap()
    }

    #[test]
    fn a_numerical_export_records_its_tolerance() {
        let read = recorded_tolerance(&file(
            "'catallax numerical approximate;source=part;source_tolerance_mm=1/100;rounding_error_mm=2e-15;recorded_tolerance_mm=1.1000000000000003e-7'",
        ))
        .unwrap();
        assert_eq!(read.tolerance_mm, q("1.1000000000000003e-7"));
        assert_eq!(read.field, RecordedField::Recorded);
    }

    #[test]
    fn a_file_written_before_check_once_records_the_limit_it_was_admitted_under() {
        let read = recorded_tolerance(&file(
            "'catallax numerical approximate;source=SolidDraft;source_tolerance_mm=1/10000000;source_checks=237'",
        ))
        .unwrap();
        assert_eq!(read.tolerance_mm, q("1/10000000"));
        assert_eq!(read.field, RecordedField::AdmittedLimit);
    }

    #[test]
    fn a_scaled_export_is_read_at_its_scale() {
        let read = recorded_tolerance(&file(
            "'catallax numerical approximate;source=part;source_tolerance_mm=1/100;recorded_tolerance_mm=2e-6;scaled_export=x1000;written_values=millimetres'",
        ))
        .unwrap();
        assert_eq!(read.tolerance_mm, q("0.002"));
    }

    #[test]
    fn a_name_cannot_stand_in_for_the_writers_field() {
        // The name holds a field of its own; the writer's comes after it.
        let read = recorded_tolerance(&file(
            "'catallax numerical approximate;source=a;recorded_tolerance_mm=5;source_tolerance_mm=1/100;recorded_tolerance_mm=3e-7'",
        ))
        .unwrap();
        assert_eq!(read.tolerance_mm, q("3e-7"));
    }

    #[test]
    fn other_files_record_nothing() {
        for description in [
            "'display'",
            "'catallax export'",
            "'catallax numerical approximate;source=part'",
            "'catallax numerical approximate;source=part;recorded_tolerance_mm=0'",
            "'catallax numerical approximate;source=part;recorded_tolerance_mm=wide'",
            "'not;catallax numerical approximate;recorded_tolerance_mm=1'",
        ] {
            assert_eq!(
                recorded_tolerance(&file(description)),
                None,
                "{description}"
            );
        }
        assert_eq!(
            recorded_tolerance(b"ISO-10303-21;\nHEADER;\nENDSEC;\n"),
            None
        );
        // A string that looks like the start of the data section is a string.
        let read = recorded_tolerance(&file(
            "'catallax numerical approximate;source=DATA;ENDSEC;recorded_tolerance_mm=1e-6'",
        ))
        .unwrap();
        assert_eq!(read.tolerance_mm, q("1e-6"));
    }
}

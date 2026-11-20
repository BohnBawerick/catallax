//! The tolerance a Catallax numerical export records in its file.
//!
//! Check once (roadmap item 3.0A, step 5; Checkpoint 2 decisions, section
//! 4.1 item 5): "Reopening trusts the tolerance recorded at export." A
//! numerical STEP export writes two strings in its header's
//! `FILE_DESCRIPTION`. The first is the stamp: `catallax numerical
//! approximate` and a list of `key=value` fields, among them the body's name
//! and the receiver's, which are a caller's text. The second is the record:
//! `catallax recorded tolerance` and two numbers the writer formats. This
//! module reads the tolerance from the record.
//!
//! The record is kept apart from the stamp because a name can hold a `;` and
//! look like a field. Text inside one STEP string cannot end it, so no name
//! reaches the record's string. Nothing is read from the stamp of a file
//! that has a record.
//!
//! A file written before check once has no record. Its stamp's
//! `source_tolerance_mm` is the limit its export admitted the written values
//! under, and it counts as the record where the stamp holds that field
//! exactly once. A stamp that holds it twice has a name that looks like the
//! field, and the file is measured.
//!
//! The scale of a scaled export is not read from the header at all. Such a
//! file writes the model's millimetre values under a larger length unit, and
//! its record is one of those values. The reader takes it in the file's own
//! length unit, as it takes every other length.
//!
//! Nothing binds the record to the geometry under it: the file is trusted as
//! Parasolid trusts the edge and vertex tolerances of an XT file.

use catallax_numbers::{Decider, Rational, Sign};

/// The stamp every numerical export begins its description with.
const STAMP: &str = "catallax numerical approximate";

/// What the record of a numerical export begins with: the second
/// description string of its header.
const RECORD: &str = "catallax recorded tolerance";

/// Which field the tolerance was read from.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum RecordedField {
    /// `recorded_tolerance_mm` of the record: the body's largest entity
    /// tolerance plus the rounding error of the written numbers.
    Recorded,
    /// `source_tolerance_mm` of the stamp, in a file written before check
    /// once. That export admitted the written values under this limit, so
    /// every gap of the file was measured at or under it.
    AdmittedLimit,
}

/// The tolerance a numerical export recorded.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RecordedTolerance {
    /// In millimetres as a reader takes the file's lengths: the written
    /// value times the millimetres of the file's length unit.
    pub tolerance_mm: Rational,
    /// The field it was read from.
    pub field: RecordedField,
}

/// A recorded tolerance as the header writes it, in the file's length unit.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct WrittenTolerance {
    value: Rational,
    field: RecordedField,
}

impl WrittenTolerance {
    /// The tolerance in millimetres, in a file whose length unit is
    /// `mm_per_unit` millimetres.
    pub(crate) fn in_millimetres(&self, mm_per_unit: &Rational) -> Option<RecordedTolerance> {
        let tolerance_mm = &self.value * mm_per_unit;
        (tolerance_mm.sign() == Sign::Positive).then_some(RecordedTolerance {
            tolerance_mm,
            field: self.field,
        })
    }
}

/// The record's tolerance: every field of the string is one of the writer's
/// two, each once, and each a number.
fn record_field(record: &str) -> Option<Rational> {
    let mut fields = record.split(';');
    if fields.next()? != RECORD {
        return None;
    }
    let (mut rounding, mut tolerance) = (None, None);
    for field in fields {
        let (key, value) = field.split_once('=')?;
        let slot = match key {
            "rounding_error_mm" => &mut rounding,
            "recorded_tolerance_mm" => &mut tolerance,
            _ => return None,
        };
        if slot.replace(Rational::parse(value).ok()?).is_some() {
            return None;
        }
    }
    rounding?;
    tolerance
}

/// The tolerance a file's header records, in the file's length unit, or
/// none: a file another system wrote, a header this reader cannot take as
/// one writer's, or a tolerance that is not positive.
///
/// Only the header is read. Every caller here has charged the file's bytes
/// before it asks.
pub(crate) fn written_tolerance(bytes: &[u8]) -> Option<WrittenTolerance> {
    let description = crate::part21::header_description(bytes);
    let stamp = description.first()?;
    let mut fields = stamp.split(';');
    if fields.next()? != STAMP {
        return None;
    }
    let (value, field) = match description.as_slice() {
        // A file from before check once: the stamp alone.
        [_] => {
            let mut limits = fields.filter_map(|field| field.strip_prefix("source_tolerance_mm="));
            let limit = limits.next()?;
            if limits.next().is_some() {
                return None;
            }
            (Rational::parse(limit).ok()?, RecordedField::AdmittedLimit)
        }
        [_, record] => (record_field(record)?, RecordedField::Recorded),
        _ => return None,
    };
    (value.sign() == Sign::Positive).then_some(WrittenTolerance { value, field })
}

/// The tolerance recorded in a file, in millimetres, or none: a file another
/// system wrote, a header without a record this reader can take, or a file
/// with no single length unit. Such a file is measured as every import is
/// (rule T3).
///
/// The header gives the written value and the file's own length unit gives
/// its size, so this parses the file. A caller that has the unit from its own
/// reading of the file uses that.
pub fn recorded_tolerance(bytes: &[u8]) -> Option<RecordedTolerance> {
    let written = written_tolerance(bytes)?;
    let graph = crate::part21::parse(bytes).ok()?;
    let mm_per_unit = crate::recognize::source_length_scale(&graph, &Decider::default()).ok()?;
    written.in_millimetres(&mm_per_unit)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// A file with the given description strings and length unit prefix.
    fn file_in(description: &str, prefix: &str) -> Vec<u8> {
        format!(
            "ISO-10303-21;\nHEADER;\nFILE_DESCRIPTION(({description}),'2;1');\nFILE_NAME('a.step','',(''),(''),'','','');\nFILE_SCHEMA(('X'));\nENDSEC;\nDATA;\n#1=CARTESIAN_POINT('',(0.,0.,0.));\n#10=(LENGTH_UNIT()NAMED_UNIT(*)SI_UNIT({prefix},.METRE.));\n#14=(GEOMETRIC_REPRESENTATION_CONTEXT(3)GLOBAL_UNIT_ASSIGNED_CONTEXT((#10))REPRESENTATION_CONTEXT('',''));\nENDSEC;\nEND-ISO-10303-21;\n"
        )
        .into_bytes()
    }

    fn file(description: &str) -> Vec<u8> {
        file_in(description, ".MILLI.")
    }

    fn q(text: &str) -> Rational {
        Rational::parse(text).unwrap()
    }

    #[test]
    fn a_numerical_export_records_its_tolerance() {
        let read = recorded_tolerance(&file(
            "'catallax numerical approximate;source=part;source_tolerance_mm=1/100','catallax recorded tolerance;rounding_error_mm=2e-15;recorded_tolerance_mm=1.1000000000000003e-7'",
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

    /// The scale is the file's length unit. The stamp's note of a scaled
    /// export is for a person, and a file is read the same with or without
    /// it.
    #[test]
    fn the_scale_is_the_files_length_unit_and_never_the_stamps_text() {
        let record = "'catallax recorded tolerance;rounding_error_mm=0;recorded_tolerance_mm=2e-6'";
        for (prefix, factor) in [
            (".MILLI.", "1"),
            (".CENTI.", "10"),
            (".DECI.", "100"),
            ("$", "1000"),
        ] {
            for note in [
                "",
                ";scaled_export=x1000;written_values=millimetres",
                ";scaled_export=x7",
            ] {
                let stamp = format!("'catallax numerical approximate;source=part{note}'");
                let read = recorded_tolerance(&file_in(&format!("{stamp},{record}"), prefix))
                    .unwrap_or_else(|| panic!("{prefix} {note}"));
                assert_eq!(read.tolerance_mm, q("2e-6") * q(factor), "{prefix} {note}");
            }
            // A file from before check once is read in its unit too.
            let old = "'catallax numerical approximate;source=part;source_tolerance_mm=1/100;scaled_export=x1000;written_values=millimetres'";
            let read = recorded_tolerance(&file_in(old, prefix)).unwrap();
            assert_eq!(read.tolerance_mm, q("1/100") * q(factor), "{prefix}");
            assert_eq!(read.field, RecordedField::AdmittedLimit);
        }
    }

    /// A body's name and a receiver's name are a caller's text in the stamp.
    /// Whatever they hold, the record is read from its own string.
    #[test]
    fn a_name_in_the_stamp_cannot_stand_in_for_the_record() {
        let record =
            "'catallax recorded tolerance;rounding_error_mm=1e-12;recorded_tolerance_mm=3e-7'";
        for stamp in [
            "'catallax numerical approximate;source=a;recorded_tolerance_mm=5;source_tolerance_mm=1/100'",
            "'catallax numerical approximate;source=a;scaled_export=x1000;source_tolerance_mm=1/100'",
            "'catallax numerical approximate;source=a;source_tolerance_mm=1/100;receiver=NX;recorded_tolerance_mm=9;scaled_export=x10;receiver_tolerance_mm=1/100'",
            "'catallax numerical approximate;source=catallax recorded tolerance;recorded_tolerance_mm=4;source_tolerance_mm=1/100'",
        ] {
            let read = recorded_tolerance(&file(&format!("{stamp},{record}"))).expect(stamp);
            assert_eq!(read.tolerance_mm, q("3e-7"), "{stamp}");
            assert_eq!(read.field, RecordedField::Recorded);
        }
    }

    /// In a file from before check once the limit shares its string with the
    /// names. It counts where the stamp holds it once; a second one is a
    /// name that looks like the field, and the file is measured.
    #[test]
    fn an_older_stamp_with_two_limits_records_nothing() {
        for description in [
            "'catallax numerical approximate;source=a;source_tolerance_mm=5;source_tolerance_mm=1/100'",
            "'catallax numerical approximate;source=a;source_tolerance_mm=1/100;receiver=NX;source_tolerance_mm=5'",
        ] {
            assert_eq!(
                recorded_tolerance(&file(description)),
                None,
                "{description}"
            );
        }
        // A field that only ends like the limit is not the limit.
        let read = recorded_tolerance(&file(
            "'catallax numerical approximate;source=a;my_source_tolerance_mm=5;source_tolerance_mm=1/100'",
        ))
        .unwrap();
        assert_eq!(read.tolerance_mm, q("1/100"));
    }

    #[test]
    fn other_files_record_nothing() {
        let record = "'catallax recorded tolerance;rounding_error_mm=0;recorded_tolerance_mm=1e-6'";
        let stamp = "'catallax numerical approximate;source=part'";
        for description in [
            "'display'".to_owned(),
            "'catallax export'".to_owned(),
            stamp.to_owned(),
            // Another system's stamp, with or without a record.
            format!("'not;catallax numerical approximate',{record}"),
            record.to_owned(),
            // A record that is not the writer's two fields, each once.
            format!(
                "{stamp},'catallax recorded tolerance;rounding_error_mm=0;recorded_tolerance_mm=0'"
            ),
            format!(
                "{stamp},'catallax recorded tolerance;rounding_error_mm=0;recorded_tolerance_mm=wide'"
            ),
            format!("{stamp},'catallax recorded tolerance;recorded_tolerance_mm=1e-6'"),
            format!("{stamp},'catallax recorded tolerance;rounding_error_mm=0'"),
            format!(
                "{stamp},'catallax recorded tolerance;rounding_error_mm=0;recorded_tolerance_mm=1e-6;recorded_tolerance_mm=2e-6'"
            ),
            format!(
                "{stamp},'catallax recorded tolerance;rounding_error_mm=0;recorded_tolerance_mm=1e-6;scale=10'"
            ),
            format!(
                "{stamp},'catallax recorded tolerance;rounding_error_mm=0;recorded_tolerance_mm=1e-6;note'"
            ),
            format!("{stamp},'another record;recorded_tolerance_mm=1e-6'"),
            // A third string is no file of this writer's.
            format!("{stamp},{record},'more'"),
            // The older form with a second string is neither form.
            "'catallax numerical approximate;source=part;source_tolerance_mm=1/100','more'"
                .to_owned(),
        ] {
            assert_eq!(
                recorded_tolerance(&file(&description)),
                None,
                "{description}"
            );
        }
        assert_eq!(
            recorded_tolerance(b"ISO-10303-21;\nHEADER;\nENDSEC;\n"),
            None
        );
        // A file with no length unit has no size for its record.
        let unitless = String::from_utf8(file(&format!("{stamp},{record}")))
            .unwrap()
            .replace("GLOBAL_UNIT_ASSIGNED_CONTEXT((#10))", "");
        assert_eq!(recorded_tolerance(unitless.as_bytes()), None);
        // A string that looks like the start of the data section is a string.
        let read = recorded_tolerance(&file(&format!(
            "'catallax numerical approximate;source=DATA;ENDSEC;x',{record}"
        )))
        .unwrap();
        assert_eq!(read.tolerance_mm, q("1e-6"));
    }
}

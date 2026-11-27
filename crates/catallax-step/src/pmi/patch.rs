//! Narrow edits of a STEP file's PMI records that keep every other record as
//! it was written.
//!
//! [`Records`] splits the data section of a Part 21 file into its records
//! by their exact bytes. A [`Patch`] removes or rewrites some of them and adds
//! new ones under fresh numbers above the file's highest. Everything else,
//! including the header, comments, spacing and line ends, is copied byte for
//! byte.

use std::collections::{BTreeMap, BTreeSet};
use std::ops::Range;

use catallax_numbers::Refusal;

fn refuse(what: impl std::fmt::Display) -> Refusal {
    Refusal::new("step.pmi-patch", what.to_string())
}

/// One record of the data section.
#[derive(Clone, Debug)]
pub struct Record {
    /// The entity number.
    pub id: u64,
    /// The record's bytes from its `#` to its `;` and the line end after it.
    pub span: Range<usize>,
    /// The bytes after `=` and before `;`.
    pub body: Range<usize>,
    /// The entities it refers to, in order, without repeats.
    pub refs: Vec<u64>,
    reference_spans: Vec<(u64, Range<usize>)>,
}

/// The records of a Part 21 file.
pub struct Records<'a> {
    text: &'a str,
    records: Vec<Record>,
    index: BTreeMap<u64, usize>,
    users: BTreeMap<u64, Vec<u64>>,
    /// Where new records go: the start of the data section's `ENDSEC;`.
    end: usize,
    newline: &'static str,
}

/// Skips a string or comment at `at`, returning the byte after it.
fn skip(bytes: &[u8], at: usize) -> Option<usize> {
    match bytes[at] {
        b'\'' => {
            let mut i = at + 1;
            loop {
                match bytes.get(i)? {
                    b'\'' if bytes.get(i + 1) == Some(&b'\'') => i += 2,
                    b'\'' => return Some(i + 1),
                    _ => i += 1,
                }
            }
        }
        b'/' if bytes.get(at + 1) == Some(&b'*') => {
            let mut i = at + 2;
            while bytes.get(i..i + 2)? != b"*/" {
                i += 1;
            }
            Some(i + 2)
        }
        _ => None,
    }
}

fn layout(bytes: &[u8], mut at: usize) -> usize {
    while let Some(byte) = bytes.get(at) {
        if byte.is_ascii_whitespace() {
            at += 1;
        } else if *byte == b'/'
            && let Some(next) = skip(bytes, at)
        {
            at = next;
        } else {
            break;
        }
    }
    at
}

fn token(bytes: &[u8], at: usize, word: &[u8]) -> Option<usize> {
    let tail = bytes.get(at..)?.strip_prefix(word)?;
    let end = layout(bytes, bytes.len() - tail.len());
    (bytes.get(end) == Some(&b';')).then_some(end + 1)
}

impl<'a> Records<'a> {
    /// Splits a file into records. Refuses text that is not UTF-8, a file
    /// with other than one data section, and a record that does not close.
    pub fn parse(bytes: &'a [u8]) -> Result<Self, Refusal> {
        let text = std::str::from_utf8(bytes).map_err(|_| refuse("the file is not UTF-8 text"))?;
        let newline = if text.contains("\r\n") { "\r\n" } else { "\n" };
        let start = token(bytes, layout(bytes, 0), b"ISO-10303-21")
            .and_then(|at| token(bytes, layout(bytes, at), b"HEADER"))
            .ok_or_else(|| refuse("the Part 21 header is missing"))?;
        let mut cursor = start;
        let start = loop {
            cursor = layout(bytes, cursor);
            if cursor >= bytes.len() {
                return Err(refuse("the Part 21 header does not end"));
            }
            if let Some(end) = token(bytes, cursor, b"ENDSEC") {
                break token(bytes, layout(bytes, end), b"DATA")
                    .ok_or_else(|| refuse("unsupported exchange section before DATA"))?;
            }
            loop {
                if cursor >= bytes.len() {
                    return Err(refuse("a header record does not close"));
                }
                if let Some(next) = skip(bytes, cursor) {
                    cursor = next;
                } else if bytes[cursor] == b';' {
                    cursor += 1;
                    break;
                } else {
                    cursor += 1;
                }
            }
        };
        let mut records = Vec::new();
        let mut i = start;
        let end = loop {
            i = layout(bytes, i);
            if i >= bytes.len() {
                return Err(refuse("the data section does not end"));
            }
            if let Some(after_data) = token(bytes, i, b"ENDSEC") {
                let end = token(bytes, layout(bytes, after_data), b"END-ISO-10303-21")
                    .ok_or_else(|| refuse("unsupported exchange section after DATA"))?;
                if layout(bytes, end) != bytes.len() {
                    return Err(refuse(
                        "unsupported exchange section after the closing token",
                    ));
                }
                break i;
            }
            if bytes[i] != b'#' {
                return Err(refuse(format!("unexpected text at byte {i}")));
            }
            let first = i;
            let mut j = i + 1;
            while bytes.get(j).is_some_and(u8::is_ascii_digit) {
                j += 1;
            }
            let id: u64 = text[i + 1..j]
                .parse()
                .map_err(|_| refuse(format!("a record at byte {i} has no number")))?;
            j = layout(bytes, j);
            if bytes.get(j) != Some(&b'=') {
                return Err(refuse(format!("record #{id} has no `=`")));
            }
            let body_start = j + 1;
            let mut refs = Vec::new();
            let mut seen_refs = BTreeSet::new();
            let mut reference_spans = Vec::new();
            let mut k = body_start;
            let body_end = loop {
                let Some(&c) = bytes.get(k) else {
                    return Err(refuse(format!("record #{id} does not close")));
                };
                if let Some(next) = skip(bytes, k) {
                    k = next;
                    continue;
                }
                match c {
                    b';' => break k,
                    b'#' => {
                        let mut m = k + 1;
                        while bytes.get(m).is_some_and(u8::is_ascii_digit) {
                            m += 1;
                        }
                        if let Ok(target) = text[k + 1..m].parse::<u64>() {
                            reference_spans.push((target, k..m));
                            if seen_refs.insert(target) {
                                refs.push(target);
                            }
                        }
                        k = m;
                    }
                    _ => k += 1,
                }
            };
            let mut span_end = body_end + 1;
            if bytes[span_end..].starts_with(b"\r\n") {
                span_end += 2;
            } else if bytes.get(span_end) == Some(&b'\n') {
                span_end += 1;
            }
            records.push(Record {
                id,
                span: first..span_end,
                body: body_start..body_end,
                refs,
                reference_spans,
            });
            i = span_end;
        };
        let mut index = BTreeMap::new();
        let mut users: BTreeMap<u64, Vec<u64>> = BTreeMap::new();
        for (position, record) in records.iter().enumerate() {
            if index.insert(record.id, position).is_some() {
                return Err(refuse(format!("record #{} is numbered twice", record.id)));
            }
            for target in &record.refs {
                users.entry(*target).or_default().push(record.id);
            }
        }
        Ok(Self {
            text,
            records,
            index,
            users,
            end,
            newline,
        })
    }

    /// The record numbered `id`.
    pub fn get(&self, id: u64) -> Option<&Record> {
        self.index.get(&id).map(|position| &self.records[*position])
    }

    /// A record's full text from `#` to and with its line end.
    pub fn text(&self, id: u64) -> Option<&'a str> {
        self.get(id).map(|record| &self.text[record.span.clone()])
    }

    /// The records that refer to `id`.
    pub fn users(&self, id: u64) -> &[u64] {
        self.users.get(&id).map_or(&[][..], Vec::as_slice)
    }

    /// Every record, in file order.
    pub fn all(&self) -> &[Record] {
        &self.records
    }

    /// Substitutes a unique reference token in a record's body.
    pub fn replace_reference(&self, id: u64, old: u64, new: u64) -> Result<String, Refusal> {
        let record = self
            .get(id)
            .ok_or_else(|| refuse(format!("record #{id} is missing")))?;
        let mut matches = record
            .reference_spans
            .iter()
            .filter(|(target, _)| *target == old);
        let (Some((_, span)), None) = (matches.next(), matches.next()) else {
            return Err(refuse(format!(
                "record #{id} does not refer to #{old} exactly once"
            )));
        };
        Ok(format!(
            "{}#{new}{}",
            &self.text[record.body.start..span.start],
            &self.text[span.end..record.body.end],
        ))
    }
}

/// Removals, rewrites and new records for one file.
pub struct Patch {
    removed: BTreeSet<u64>,
    rewritten: BTreeMap<u64, String>,
    added: Vec<(u64, String)>,
    highest: u64,
}

impl Patch {
    /// An empty patch of `records`; new records are numbered above its highest.
    pub fn new(records: &Records<'_>) -> Self {
        let highest = records
            .records
            .iter()
            .map(|record| record.id)
            .max()
            .unwrap_or(0);
        Self {
            removed: BTreeSet::new(),
            rewritten: BTreeMap::new(),
            added: Vec::new(),
            highest,
        }
    }

    /// Removes a record.
    pub fn remove(&mut self, id: u64) {
        self.rewritten.remove(&id);
        self.removed.insert(id);
    }

    /// Gives a record a new body under its own number.
    pub fn rewrite(&mut self, id: u64, body: String) {
        self.rewritten.insert(id, body);
    }

    /// Adds a record under a fresh number and returns the number.
    pub fn add(&mut self, body: String) -> Result<u64, Refusal> {
        self.highest = self
            .highest
            .checked_add(1)
            .ok_or_else(|| refuse("no fresh record number is available"))?;
        self.added.push((self.highest, body));
        Ok(self.highest)
    }

    /// The records this patch removes or rewrites.
    pub fn touched(&self) -> BTreeSet<u64> {
        self.removed
            .iter()
            .chain(self.rewritten.keys())
            .copied()
            .collect()
    }

    /// The file with the patch applied. Every byte outside the removed and
    /// rewritten records is copied; new records go just before the data
    /// section's `ENDSEC;`.
    pub fn apply(&self, records: &Records<'_>) -> Vec<u8> {
        let text = records.text;
        let mut out = String::with_capacity(text.len());
        let mut at = 0;
        for record in &records.records {
            if self.removed.contains(&record.id) {
                out.push_str(&text[at..record.span.start]);
                at = record.span.end;
            } else if let Some(body) = self.rewritten.get(&record.id) {
                out.push_str(&text[at..record.body.start]);
                out.push_str(body);
                at = record.body.end;
            }
        }
        out.push_str(&text[at..records.end]);
        for (id, body) in &self.added {
            out.push_str(&format!("#{id}={body};{}", records.newline));
        }
        out.push_str(&text[records.end..]);
        out.into_bytes()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const FILE: &str = "ISO-10303-21;\r\nHEADER;\r\nFILE_NAME('DATA; #9=X;','',(''),(''),'','','');\r\nENDSEC;\r\nDATA;\r\n#1=A('a;b',#2);\r\n/* #7=GHOST(); */\r\n#2=(B(#3)C());\r\n#3=D('it''s',#1,#2);\r\nENDSEC;\r\nEND-ISO-10303-21;\r\n";

    #[test]
    fn records_split_by_their_exact_bytes() {
        let records = Records::parse(FILE.as_bytes()).unwrap();
        let ids: Vec<u64> = records.all().iter().map(|record| record.id).collect();
        assert_eq!(ids, [1, 2, 3]);
        assert_eq!(records.get(3).unwrap().refs, [1, 2]);
        assert_eq!(records.users(2), [1, 3]);
        assert_eq!(records.text(1), Some("#1=A('a;b',#2);\r\n"));
        assert_eq!(Patch::new(&records).apply(&records), FILE.as_bytes());
    }

    #[test]
    fn a_patch_changes_only_its_records() {
        let records = Records::parse(FILE.as_bytes()).unwrap();
        let mut patch = Patch::new(&records);
        let fresh = patch.add("E()".into()).unwrap();
        assert_eq!(fresh, 4);
        patch.rewrite(3, format!("D('x',#1,#{fresh})"));
        patch.remove(2);
        let out = String::from_utf8(patch.apply(&records)).unwrap();
        assert_eq!(
            out,
            FILE.replace("#2=(B(#3)C());\r\n", "")
                .replace("#3=D('it''s',#1,#2);", "#3=D('x',#1,#4);")
                .replace("ENDSEC;\r\nEND", "#4=E();\r\nENDSEC;\r\nEND")
        );
        assert!(Records::parse(b"DATA;\n#1=A(;\n").is_err());
    }
}

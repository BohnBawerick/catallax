//! Narrow edits of a STEP file's PMI records that keep every other record as
//! it was written.
//!
//! [`Records`] splits the data section of a Part 21 file into its records
//! by their exact bytes. A [`Patch`] removes or rewrites some of them and adds
//! new ones under fresh numbers above the file's highest. Everything else,
//! including the header, comments, spacing and line ends, is copied byte for
//! byte. Which records an edit may touch is decided by [`Records::owned`]: the
//! records only one PMI item uses.

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
    /// The entity type, or the first partial type of a complex instance.
    pub kind: String,
    /// The entities it refers to, in order, without repeats.
    pub refs: Vec<u64>,
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

impl<'a> Records<'a> {
    /// Splits a file into records. Refuses text that is not UTF-8, a file
    /// with other than one data section, and a record that does not close.
    pub fn parse(bytes: &'a [u8]) -> Result<Self, Refusal> {
        let text = std::str::from_utf8(bytes).map_err(|_| refuse("the file is not UTF-8 text"))?;
        let newline = if text.contains("\r\n") { "\r\n" } else { "\n" };
        let mut i = 0;
        let mut data = None;
        let mut sections = 0;
        // Find the data section outside strings and comments.
        while i < bytes.len() {
            if let Some(next) = skip(bytes, i) {
                i = next;
                continue;
            }
            if bytes[i..].starts_with(b"DATA") {
                let mut j = i + 4;
                while bytes.get(j).is_some_and(u8::is_ascii_whitespace) {
                    j += 1;
                }
                let preceded =
                    i == 0 || !(bytes[i - 1].is_ascii_alphanumeric() || bytes[i - 1] == b'_');
                if preceded && bytes.get(j) == Some(&b';') {
                    sections += 1;
                    data.get_or_insert(j + 1);
                    i = j + 1;
                    continue;
                }
            }
            i += 1;
        }
        let start = match (data, sections) {
            (Some(start), 1) => start,
            _ => return Err(refuse("the file does not have exactly one data section")),
        };
        let mut records = Vec::new();
        let mut i = start;
        let end = loop {
            while bytes.get(i).is_some_and(u8::is_ascii_whitespace) {
                i += 1;
            }
            if i >= bytes.len() {
                return Err(refuse("the data section does not end"));
            }
            if let Some(next) = skip(bytes, i).filter(|_| bytes[i] == b'/') {
                i = next;
                continue;
            }
            if bytes[i..].starts_with(b"ENDSEC") {
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
            while bytes.get(j).is_some_and(u8::is_ascii_whitespace) {
                j += 1;
            }
            if bytes.get(j) != Some(&b'=') {
                return Err(refuse(format!("record #{id} has no `=`")));
            }
            let body_start = j + 1;
            let mut refs = Vec::new();
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
                        if let Ok(target) = text[k + 1..m].parse::<u64>()
                            && !refs.contains(&target)
                        {
                            refs.push(target);
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
            let body = text[body_start..body_end].trim_start();
            let body = body.strip_prefix('(').unwrap_or(body).trim_start();
            let kind: String = body
                .chars()
                .take_while(|c| c.is_ascii_alphanumeric() || *c == '_')
                .collect();
            records.push(Record {
                id,
                span: first..span_end,
                body: body_start..body_end,
                kind,
                refs,
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

    /// A record's bytes between `=` and `;`.
    pub fn body(&self, id: u64) -> Option<&'a str> {
        self.get(id).map(|record| &self.text[record.body.clone()])
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

    /// The records one PMI item owns, starting from its root entity: the
    /// records only those records use, and the records of the kinds that
    /// attach to an item (its name, zone, limits, usages, properties and
    /// relationships) together with what only they use. The part, its
    /// geometry, contexts, units, datums and the shared presentation are never
    /// owned.
    pub fn owned(&self, root: u64) -> BTreeSet<u64> {
        let mut owned = BTreeSet::from([root]);
        loop {
            let mut grown = Vec::new();
            for id in &owned {
                let Some(record) = self.get(*id) else {
                    continue;
                };
                for target in &record.refs {
                    if owned.contains(target) {
                        continue;
                    }
                    let Some(child) = self.get(*target) else {
                        continue;
                    };
                    if !shared(&child.kind)
                        && self.users(*target).iter().all(|user| owned.contains(user))
                    {
                        grown.push(*target);
                    }
                }
                for user in self.users(*id) {
                    if !owned.contains(user)
                        && self.get(*user).is_some_and(|record| attaches(&record.kind))
                    {
                        grown.push(*user);
                    }
                }
            }
            if grown.is_empty() {
                return owned;
            }
            owned.extend(grown);
        }
    }

    /// The records outside `set` that refer to a record in it.
    pub fn outside_users(&self, set: &BTreeSet<u64>) -> Vec<u64> {
        let mut found = BTreeSet::new();
        for id in set {
            for user in self.users(*id) {
                if !set.contains(user) {
                    found.insert(*user);
                }
            }
        }
        found.into_iter().collect()
    }
}

/// Kinds a PMI item shares with the rest of the file and never owns.
fn shared(kind: &str) -> bool {
    matches!(
        kind,
        "PRODUCT_DEFINITION_SHAPE"
            | "PRODUCT_DEFINITION"
            | "PRODUCT"
            | "ADVANCED_FACE"
            | "FACE_SURFACE"
            | "DATUM"
            | "DATUM_FEATURE"
            | "DRAUGHTING_MODEL"
            | "ANNOTATION_PLANE"
            | "PRESENTATION_STYLE_ASSIGNMENT"
            | "APPLICATION_CONTEXT"
            | "DOCUMENT"
            | "DIMENSIONAL_EXPONENTS"
    ) || kind.ends_with("_CONTEXT")
        || kind.ends_with("UNIT")
        || kind.ends_with("SHAPE_REPRESENTATION")
}

/// Kinds that attach to a PMI item and go with it.
fn attaches(kind: &str) -> bool {
    matches!(
        kind,
        "ID_ATTRIBUTE"
            | "TOLERANCE_ZONE"
            | "PROJECTED_ZONE_DEFINITION"
            | "GEOMETRIC_TOLERANCE_RELATIONSHIP"
            | "PLUS_MINUS_TOLERANCE"
            | "DIMENSIONAL_CHARACTERISTIC_REPRESENTATION"
            | "PROPERTY_DEFINITION"
            | "PROPERTY_DEFINITION_REPRESENTATION"
            | "SHAPE_DEFINITION_REPRESENTATION"
            | "GEOMETRIC_ITEM_SPECIFIC_USAGE"
            | "ITEM_IDENTIFIED_REPRESENTATION_USAGE"
            | "SHAPE_ASPECT_RELATIONSHIP"
            | "SHAPE_ASPECT_DERIVING_RELATIONSHIP"
            | "FEATURE_FOR_DATUM_TARGET_RELATIONSHIP"
            | "GENERAL_PROPERTY_ASSOCIATION"
            | "CHARACTERIZED_ITEM_WITHIN_REPRESENTATION"
            | "DRAUGHTING_MODEL_ITEM_ASSOCIATION"
    )
}

/// Removals, rewrites and new records for one file.
pub struct Patch {
    removed: BTreeSet<u64>,
    rewritten: BTreeMap<u64, String>,
    added: Vec<(u64, String)>,
    next: u64,
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
            next: highest + 1,
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
    pub fn add(&mut self, body: String) -> u64 {
        let id = self.next;
        self.next += 1;
        self.added.push((id, body));
        id
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
        assert_eq!(records.get(2).unwrap().kind, "B");
        assert_eq!(records.get(3).unwrap().refs, [1, 2]);
        assert_eq!(records.users(2), [1, 3]);
        assert_eq!(records.text(1), Some("#1=A('a;b',#2);\r\n"));
        assert_eq!(Patch::new(&records).apply(&records), FILE.as_bytes());
    }

    #[test]
    fn a_patch_changes_only_its_records() {
        let records = Records::parse(FILE.as_bytes()).unwrap();
        let mut patch = Patch::new(&records);
        let fresh = patch.add("E()".into());
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

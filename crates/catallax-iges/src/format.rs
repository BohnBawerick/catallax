//! Fixed-column IGES S/G/D/P/T records, with stable order and metadata.
use crate::{IgesError, Units, exact};
use catallax_numbers::Rational;

pub(crate) enum Value {
    Length(Rational),
    Text(String),
}
impl From<usize> for Value {
    fn from(x: usize) -> Self {
        Self::Text(x.to_string())
    }
}
pub(crate) struct Entity {
    pub kind: usize,
    pub form: usize,
    pub transform: usize,
    pub root: bool,
    pub values: Vec<Value>,
}
#[derive(Default)]
pub(crate) struct Graph {
    pub entities: Vec<Entity>,
}
impl Graph {
    pub fn add(&mut self, kind: usize, form: usize, values: Vec<Value>) -> usize {
        let id = self.entities.len() * 2 + 1;
        self.entities.push(Entity {
            kind,
            form,
            transform: 0,
            root: false,
            values,
        });
        id
    }
    pub fn lengths(&self) -> impl Iterator<Item = Rational> + '_ {
        self.entities
            .iter()
            .flat_map(|e| &e.values)
            .filter_map(|v| match v {
                Value::Length(q) => Some(q.clone()),
                _ => None,
            })
    }
    pub fn file(&self, scale: &Rational, units: Units) -> Result<String, IgesError> {
        let mut file = String::new();
        record(&mut file, "Catallax exact-surface B-rep", 'S', 1);
        let (flag, name) = match units {
            Units::Millimetres => (2, "2HMM"),
        };
        // IGES 5.3 global parameters 1..26. No host paths or wall-clock metadata.
        let global = format!(
            "1H,,1H;,8HCatallax,8Hgeometry,8HCatallax,8HCatallax,32,38,6,308,15,8HCatallax,1.0,{flag},{name},1,1.0,15H19700101.000000,1.0E-9,1.0E99,0H,0H,11,0,15H19700101.000000,0H;"
        );
        let gcount = chunks(&global, 72)?.len();
        for (i, c) in chunks(&global, 72)?.into_iter().enumerate() {
            record(&mut file, &c, 'G', i + 1);
        }
        let mut pdata = Vec::new();
        let mut pline = 1;
        let unit = units.millimetres();
        for (i, e) in self.entities.iter().enumerate() {
            let id = 2 * i + 1;
            let mut fields = vec![e.kind.to_string()];
            for value in &e.values {
                fields.push(match value {
                    Value::Length(q) => exact::decimal(&(&(q * scale) / &unit))?,
                    Value::Text(text) => text.clone(),
                });
            }
            let parameters = fields.join(",") + ";";
            let lines = chunks(&parameters, 64)?;
            // Status: visible, physically dependent except MSBO roots, geometry,
            // hierarchy top-down. Entity sequence values occupy columns 74..80.
            let status = if e.root { "00000000" } else { "00010000" };
            let d1 = format!(
                "{:>8}{:>8}{:>8}{:>8}{:>8}{:>8}{:>8}{:>8}{status}",
                e.kind, pline, 0, 0, 0, 0, e.transform, 0
            );
            let d2 = format!(
                "{:>8}{:>8}{:>8}{:>8}{:>8}{:>8}{:>8}{:>8}{:>8}",
                e.kind,
                0,
                0,
                lines.len(),
                e.form,
                "",
                "",
                "",
                0
            );
            record(&mut file, &d1, 'D', id);
            record(&mut file, &d2, 'D', id + 1);
            for line in lines {
                pdata.push(format!("{line:<64} {id:>7}P{pline:>7}\n"));
                pline += 1;
            }
        }
        for line in pdata {
            file.push_str(&line);
        }
        record(
            &mut file,
            &format!(
                "S{:>7}G{:>7}D{:>7}P{:>7}",
                1,
                gcount,
                2 * self.entities.len(),
                pline - 1
            ),
            'T',
            1,
        );
        Ok(file)
    }
}
fn record(file: &mut String, body: &str, section: char, sequence: usize) {
    file.push_str(&format!("{body:<72}{section}{sequence:>7}\n"));
}
fn chunks(mut s: &str, width: usize) -> Result<Vec<String>, IgesError> {
    let mut out = Vec::new();
    while s.len() > width {
        let split = s[..width]
            .rfind([',', ';'])
            .ok_or(IgesError::UnrepresentableDefinition(
                "numeric field exceeds IGES record width",
            ))?
            + 1;
        out.push(s[..split].to_owned());
        s = &s[split..];
    }
    if !s.is_empty() {
        out.push(s.to_owned());
    }
    Ok(out)
}
pub(crate) fn lengths(p: [Rational; 3]) -> Vec<Value> {
    p.into_iter().map(Value::Length).collect()
}
pub(crate) fn text(s: impl ToString) -> Value {
    Value::Text(s.to_string())
}

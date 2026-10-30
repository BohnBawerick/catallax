//! Entity references stay symbolic until the required depth-first numbering pass.

use std::borrow::Cow;
use std::fmt::Write as _;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(crate) struct Id(usize);

#[derive(Clone)]
pub(crate) enum Field {
    Text(Cow<'static, str>),
    Ref(Id),
    List(Vec<Field>),
    /// One named partial entity of an external-mapping complex instance.
    Part(&'static str, Vec<Field>),
}

impl From<Id> for Field {
    fn from(id: Id) -> Self {
        Self::Ref(id)
    }
}

pub(crate) fn text(value: impl Into<Cow<'static, str>>) -> Field {
    Field::Text(value.into())
}

pub(crate) fn refs(ids: impl IntoIterator<Item = Id>) -> Field {
    Field::List(ids.into_iter().map(Field::Ref).collect())
}

struct Entity {
    name: &'static str,
    fields: Vec<Field>,
}

#[derive(Default)]
pub(crate) struct Graph {
    entities: Vec<Entity>,
}

impl Graph {
    pub(crate) fn add(&mut self, name: &'static str, fields: Vec<Field>) -> Id {
        let id = Id(self.entities.len());
        self.entities.push(Entity { name, fields });
        id
    }

    /// A complex instance. STEP requires its parts in alphabetical name order.
    pub(crate) fn add_complex(&mut self, parts: Vec<(&'static str, Vec<Field>)>) -> Id {
        debug_assert!(parts.windows(2).all(|pair| pair[0].0 < pair[1].0));
        self.add(
            "",
            parts
                .into_iter()
                .map(|(name, fields)| Field::Part(name, fields))
                .collect(),
        )
    }

    pub(crate) fn duplicate(&mut self, id: Id) -> Id {
        let entity = &self.entities[id.0];
        self.add(entity.name, entity.fields.clone())
    }

    pub(crate) fn replace(&mut self, id: Id, name: &'static str, fields: Vec<Field>) {
        self.entities[id.0] = Entity { name, fields };
    }

    pub(crate) fn write(&self, roots: &[Id], header: impl FnOnce(&[usize]) -> String) -> String {
        let mut numbers = vec![0; self.entities.len()];
        let mut order = Vec::with_capacity(self.entities.len());
        // Iterative preorder avoids a call stack proportional to a spline's size.
        let mut pending: Vec<_> = roots.iter().rev().copied().collect();
        while let Some(id) = pending.pop() {
            if numbers[id.0] != 0 {
                continue;
            }
            numbers[id.0] = order.len() + 21;
            order.push(id);
            for field in self.entities[id.0].fields.iter().rev() {
                field.references_reversed(&mut pending);
            }
        }
        let roots = roots.iter().map(|id| numbers[id.0]).collect::<Vec<_>>();
        let mut output = header(&roots);
        for id in order {
            let entity = &self.entities[id.0];
            write!(output, "#{}={}(", numbers[id.0], entity.name).unwrap();
            write_fields(&entity.fields, &numbers, &mut output);
            output.push_str(");\n");
        }
        output
    }
}

impl Field {
    fn references_reversed(&self, children: &mut Vec<Id>) {
        match self {
            Self::Ref(id) => children.push(*id),
            Self::List(fields) | Self::Part(_, fields) => {
                for field in fields.iter().rev() {
                    field.references_reversed(children);
                }
            }
            Self::Text(_) => {}
        }
    }

    fn write(&self, numbers: &[usize], out: &mut String) {
        match self {
            Self::Text(text) => out.push_str(text),
            Self::Ref(id) => write!(out, "#{}", numbers[id.0]).unwrap(),
            Self::List(fields) => {
                out.push('(');
                write_fields(fields, numbers, out);
                out.push(')');
            }
            Self::Part(name, fields) => {
                out.push_str(name);
                out.push('(');
                write_fields(fields, numbers, out);
                out.push(')');
            }
        }
    }
}

fn write_fields(fields: &[Field], numbers: &[usize], out: &mut String) {
    for (i, field) in fields.iter().enumerate() {
        if i != 0 && !matches!(field, Field::Part(..)) {
            out.push(',');
        }
        field.write(numbers, out);
    }
}

pub(crate) fn boolean(value: bool) -> Field {
    text(if value { ".T." } else { ".F." })
}

pub(crate) fn string(value: &str) -> String {
    let mut out = String::from("'");
    for c in value.chars() {
        match c {
            '\'' => out.push_str("''"),
            '\\' => out.push_str("\\\\"),
            ' '..='~' => out.push(c),
            _ => {
                out.push_str("\\X2\\");
                for unit in c.encode_utf16(&mut [0; 2]) {
                    write!(out, "{unit:04X}").unwrap();
                }
                out.push_str("\\X0\\");
            }
        }
    }
    out.push('\'');
    out
}

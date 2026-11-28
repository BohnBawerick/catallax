//! Entity references stay symbolic until the required depth-first numbering pass.

use std::borrow::Cow;
use std::collections::{BTreeMap, BTreeSet};
use std::fmt::Write as _;

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub(crate) struct Id(usize);

#[derive(Clone)]
pub(crate) enum Field {
    Text(Cow<'static, str>),
    Ref(Id),
    QuotedNumber(Id),
    List(Vec<Field>),
    /// One named partial entity of an external-mapping complex instance.
    Part(&'static str, Vec<Field>),
    /// A typed parameter, such as `SET_REPRESENTATION_ITEM((#1,#2))`.
    Typed(&'static str, Vec<Field>),
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

#[derive(Clone, Copy)]
pub(crate) struct Offset(usize);

impl Offset {
    pub(crate) fn apply(self, id: Id) -> Id {
        Id(self.0 + id.0)
    }
}

enum Placement {
    Entity(Id),
    Gap,
    Reserve(Id),
    Emit(Id),
}

/// The Part 21 numbering and emission order. Some legacy slots are empty,
/// and subassembly definitions are emitted later than their numbers.
#[derive(Default)]
pub(crate) struct Layout(Vec<Placement>);

impl Layout {
    pub(crate) fn push(&mut self, id: Id) {
        self.0.push(Placement::Entity(id));
    }

    pub(crate) fn extend(&mut self, ids: impl IntoIterator<Item = Id>) {
        self.0.extend(ids.into_iter().map(Placement::Entity));
    }

    pub(crate) fn prefix(&mut self, ids: &[Option<Id>]) {
        self.0.extend(
            ids.iter()
                .map(|id| id.map_or(Placement::Gap, Placement::Entity)),
        );
    }

    pub(crate) fn reserve(&mut self, ids: &[Id]) {
        self.0.extend(ids.iter().copied().map(Placement::Reserve));
    }

    pub(crate) fn emit(&mut self, ids: &[Id]) {
        self.0.extend(ids.iter().copied().map(Placement::Emit));
    }

    pub(crate) fn append(&mut self, local: Self, offset: Offset) {
        self.0.extend(local.0.into_iter().map(|entry| match entry {
            Placement::Entity(id) => Placement::Entity(offset.apply(id)),
            Placement::Gap => Placement::Gap,
            Placement::Reserve(id) => Placement::Reserve(offset.apply(id)),
            Placement::Emit(id) => Placement::Emit(offset.apply(id)),
        }));
    }

    /// Span of a part's layout, with no deferred entities.
    pub(crate) fn part_span(&self) -> usize {
        self.0.len()
    }
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

    pub(crate) fn record(&self, id: Id) -> (&'static str, &[Field]) {
        let entity = &self.entities[id.0];
        (entity.name, &entity.fields)
    }

    pub(crate) fn records(&self) -> impl Iterator<Item = (&'static str, &[Field])> {
        self.entities
            .iter()
            .map(|entity| (entity.name, entity.fields.as_slice()))
    }

    pub(crate) fn duplicate(&mut self, id: Id) -> Id {
        let entity = &self.entities[id.0];
        self.add(entity.name, entity.fields.clone())
    }

    pub(crate) fn replace(&mut self, id: Id, name: &'static str, fields: Vec<Field>) {
        self.entities[id.0] = Entity { name, fields };
    }

    pub(crate) fn extend_list(
        &mut self,
        id: Id,
        field: usize,
        items: impl IntoIterator<Item = Id>,
    ) {
        let Field::List(list) = &mut self.entities[id.0].fields[field] else {
            unreachable!("expected a representation item list");
        };
        list.extend(items.into_iter().map(Field::Ref));
    }

    pub(crate) fn append(&mut self, mut local: Self) -> Offset {
        let offset = Offset(self.entities.len());
        for entity in &mut local.entities {
            for field in &mut entity.fields {
                field.offset(offset);
            }
        }
        self.entities.append(&mut local.entities);
        offset
    }

    pub(crate) fn write(&self, roots: &[Id], header: impl FnOnce(&[usize]) -> String) -> String {
        let mut layout = Layout::default();
        layout.prefix(&[None; 20]);
        self.depth_first(&mut layout, &[], roots.iter().copied());
        let numbers = self.number(&layout);
        let roots = roots.iter().map(|id| numbers[id.0]).collect::<Vec<_>>();
        let mut output = header(&roots);
        self.write_data(&layout, &numbers, &mut output);
        output
    }

    /// Extend a layout with the reachable geometry and annotations in preorder,
    /// without following references back into the standing prefix.
    pub(crate) fn depth_first(
        &self,
        layout: &mut Layout,
        prefix: &[Option<Id>],
        roots: impl DoubleEndedIterator<Item = Id>,
    ) {
        let mut seen = vec![false; self.entities.len()];
        for id in prefix.iter().flatten() {
            seen[id.0] = true;
        }
        let mut pending: Vec<_> = roots.rev().collect();
        while let Some(id) = pending.pop() {
            if std::mem::replace(&mut seen[id.0], true) {
                continue;
            }
            layout.push(id);
            for field in self.entities[id.0].fields.iter().rev() {
                field.references_reversed(&mut pending);
            }
        }
    }

    fn number(&self, layout: &Layout) -> Vec<usize> {
        let mut numbers = vec![0; self.entities.len()];
        let mut next = 1;
        for entry in &layout.0 {
            match *entry {
                Placement::Entity(id) | Placement::Reserve(id) => numbers[id.0] = next,
                Placement::Gap => {}
                Placement::Emit(_) => continue,
            }
            next += 1;
        }
        numbers
    }

    pub(crate) fn unique_draughting_model_names(&mut self) {
        let (mut names, cameras): (Vec<_>, Vec<_>) = self
            .entities
            .iter_mut()
            .filter_map(|entity| {
                if entity.name == "DRAUGHTING_MODEL"
                    && let [Field::Text(name), Field::List(members), _] = &mut entity.fields[..]
                {
                    let camera = match members.last() {
                        Some(Field::Ref(id)) => Some(*id),
                        _ => None,
                    };
                    Some((name, camera))
                } else {
                    None
                }
            })
            .unzip();
        unique_quoted_names(&mut names);
        let names: Vec<_> = names.into_iter().map(|name| (*name).clone()).collect();
        for (name, camera) in names.into_iter().zip(cameras) {
            if let Some(camera) = camera {
                let entity = &mut self.entities[camera.0];
                if entity.name == "CAMERA_MODEL_D3" {
                    entity.fields[0] = Field::Text(name);
                }
            }
        }
    }

    pub(crate) fn write_layout(&self, layout: &Layout, mut output: String) -> String {
        self.write_data(layout, &self.number(layout), &mut output);
        output
    }

    fn write_data(&self, layout: &Layout, numbers: &[usize], output: &mut String) {
        for entry in &layout.0 {
            let id = match *entry {
                Placement::Entity(id) | Placement::Emit(id) => id,
                Placement::Gap | Placement::Reserve(_) => continue,
            };
            let entity = &self.entities[id.0];
            write!(output, "#{}={}(", numbers[id.0], entity.name).unwrap();
            write_fields(&entity.fields, numbers, output);
            output.push_str(");\n");
        }
    }
}

impl Field {
    fn offset(&mut self, offset: Offset) {
        match self {
            Self::Ref(id) | Self::QuotedNumber(id) => *id = offset.apply(*id),
            Self::List(fields) | Self::Part(_, fields) | Self::Typed(_, fields) => {
                for field in fields {
                    field.offset(offset);
                }
            }
            Self::Text(_) => {}
        }
    }

    fn references_reversed(&self, children: &mut Vec<Id>) {
        match self {
            Self::Ref(id) => children.push(*id),
            Self::List(fields) | Self::Part(_, fields) | Self::Typed(_, fields) => {
                for field in fields.iter().rev() {
                    field.references_reversed(children);
                }
            }
            Self::Text(_) | Self::QuotedNumber(_) => {}
        }
    }

    fn write(&self, numbers: &[usize], out: &mut String) {
        match self {
            Self::Text(text) => out.push_str(text),
            Self::Ref(id) => write!(out, "#{}", numbers[id.0]).unwrap(),
            Self::QuotedNumber(id) => write!(out, "'{}'", numbers[id.0]).unwrap(),
            Self::List(fields) => {
                out.push('(');
                write_fields(fields, numbers, out);
                out.push(')');
            }
            Self::Part(name, fields) | Self::Typed(name, fields) => {
                out.push_str(name);
                out.push('(');
                write_fields(fields, numbers, out);
                out.push(')');
            }
        }
    }
}

fn unique_quoted_names(names: &mut [&mut Cow<'static, str>]) -> usize {
    let bases = names.iter().map(|name| {
        let name = name
            .strip_prefix('\'')
            .and_then(|name| name.strip_suffix('\''))
            .expect("a quoted model name");
        (name, if name.is_empty() { "PMI" } else { name })
    });
    let (allocated, tried) = crate::pmi::unique_view_names(bases);
    for (name, allocated) in names.iter_mut().zip(allocated) {
        **name = Cow::Owned(format!("'{allocated}'"));
    }
    tried
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
            '\u{0}'..='\u{1F}' | '\u{7F}' => write!(out, "\\X\\{:02X}", u32::from(c)).unwrap(),
            '\u{80}'..='\u{FFFF}' => write!(out, "\\X2\\{:04X}\\X0\\", u32::from(c)).unwrap(),
            _ => write!(out, "\\X4\\{:08X}\\X0\\", u32::from(c)).unwrap(),
        }
    }
    out.push('\'');
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn quoted(names: &[&str]) -> Vec<Cow<'static, str>> {
        names
            .iter()
            .map(|name| Cow::Owned(format!("'{name}'")))
            .collect()
    }

    #[test]
    fn repeated_model_names_take_free_suffixes() {
        let mut names = quoted(&["", "Top", "", "Top", "Top (2)", "", "PMI (3)"]);
        let mut refs: Vec<_> = names.iter_mut().collect();
        unique_quoted_names(&mut refs);
        let got: Vec<&str> = names.iter().map(|name| name.as_ref()).collect();
        assert_eq!(
            got,
            [
                "''",
                "'Top'",
                "'PMI (2)'",
                "'Top (3)'",
                "'Top (2)'",
                "'PMI (4)'",
                "'PMI (3)'"
            ]
        );
    }

    #[test]
    fn many_blank_models_are_named_in_linear_work() {
        // Every finished assembly leaf may write a blank global model.
        let count = 4000;
        let mut names = quoted(&vec![""; count]);
        let mut refs: Vec<_> = names.iter_mut().collect();
        let tried = unique_quoted_names(&mut refs);
        // One candidate per repeat; restarting at 2 each time tried about
        // count * count / 2.
        assert_eq!(tried, count - 1);
        let distinct: BTreeSet<_> = names.iter().collect();
        assert_eq!(distinct.len(), count);
        assert_eq!(names[count - 1], format!("'PMI ({count})'"));
    }
}

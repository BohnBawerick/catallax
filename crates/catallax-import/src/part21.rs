//! A small ISO 10303-21 reader: entity instances as typed parameter trees.

use std::collections::BTreeMap;

use crate::ImportRefusal;

/// One STEP parameter.
#[derive(Clone, PartialEq, Eq, Debug)]
pub(crate) enum Param {
    Ref(u64),
    /// Integer or real, kept as written text.
    Number(String),
    Str(String),
    Enum(String),
    List(Vec<Param>),
    /// A typed value such as `LENGTH_MEASURE(1.)`.
    Typed(String, Vec<Param>),
    Unset,
    Derived,
}

/// One entity instance. A complex instance has several parts, in file order.
#[derive(Clone, Debug)]
pub(crate) struct Entity {
    pub parts: Vec<(String, Vec<Param>)>,
}

impl Entity {
    pub fn part(&self, name: &str) -> Option<&[Param]> {
        self.parts
            .iter()
            .find(|(n, _)| n == name)
            .map(|(_, p)| p.as_slice())
    }
    /// The type name of a simple instance, or the first part of a complex one.
    pub fn type_name(&self) -> &str {
        &self.parts[0].0
    }
}

pub(crate) type Graph = BTreeMap<u64, Entity>;

fn malformed(detail: impl Into<String>) -> ImportRefusal {
    ImportRefusal::Malformed {
        detail: detail.into(),
    }
}

/// Parses the DATA section of a Part 21 file.
pub(crate) fn parse(bytes: &[u8]) -> Result<Graph, ImportRefusal> {
    let text = std::str::from_utf8(bytes).map_err(|_| malformed("file is not UTF-8 text"))?;
    let text = strip_comments(text)?;
    let upper_start = text
        .find("DATA;")
        .ok_or_else(|| malformed("no DATA section"))?;
    let body = &text[upper_start + 5..];
    let end = body
        .find("ENDSEC;")
        .ok_or_else(|| malformed("DATA section is not closed"))?;
    let mut graph = Graph::new();
    for statement in split_statements(&body[..end])? {
        let statement = statement.trim();
        if statement.is_empty() {
            continue;
        }
        let rest = statement
            .strip_prefix('#')
            .ok_or_else(|| malformed(format!("statement without an entity name: {statement}")))?;
        let (id, value) = rest
            .split_once('=')
            .ok_or_else(|| malformed(format!("no '=' in #{rest}")))?;
        let id: u64 = id
            .trim()
            .parse()
            .map_err(|_| malformed(format!("bad entity name #{id}")))?;
        let mut cursor = Cursor::new(value.trim());
        let parts = cursor.entity_parts()?;
        cursor.skip_space();
        if !cursor.done() {
            return Err(malformed(format!("trailing text in #{id}")));
        }
        if graph.insert(id, Entity { parts }).is_some() {
            return Err(malformed(format!("entity #{id} defined twice")));
        }
    }
    Ok(graph)
}

fn strip_comments(text: &str) -> Result<String, ImportRefusal> {
    let mut out = String::with_capacity(text.len());
    let (mut quoted, mut rest) = (false, text);
    while let Some(c) = rest.chars().next() {
        if !quoted && rest.starts_with("/*") {
            let end = rest
                .find("*/")
                .ok_or_else(|| malformed("unterminated comment"))?;
            rest = &rest[end + 2..];
            out.push(' ');
            continue;
        }
        if c == '\'' {
            quoted = !quoted;
        }
        out.push(c);
        rest = &rest[c.len_utf8()..];
    }
    Ok(out)
}

fn split_statements(text: &str) -> Result<Vec<&str>, ImportRefusal> {
    let (mut quoted, mut start, mut out) = (false, 0, Vec::new());
    for (i, c) in text.char_indices() {
        match c {
            '\'' => quoted = !quoted,
            ';' if !quoted => {
                out.push(&text[start..i]);
                start = i + 1;
            }
            _ => {}
        }
    }
    if quoted {
        return Err(malformed("unterminated string"));
    }
    if !text[start..].trim().is_empty() {
        return Err(malformed("statement without terminating ';'"));
    }
    Ok(out)
}

struct Cursor<'a> {
    s: &'a str,
    i: usize,
}

impl<'a> Cursor<'a> {
    fn new(s: &'a str) -> Self {
        Self { s, i: 0 }
    }
    fn done(&self) -> bool {
        self.i >= self.s.len()
    }
    fn peek(&self) -> Option<char> {
        self.s[self.i..].chars().next()
    }
    fn skip_space(&mut self) {
        while self.peek().is_some_and(char::is_whitespace) {
            self.i += 1;
        }
    }
    fn expect(&mut self, c: char) -> Result<(), ImportRefusal> {
        self.skip_space();
        if self.peek() == Some(c) {
            self.i += 1;
            Ok(())
        } else {
            Err(malformed(format!("expected '{c}' in `{}`", self.s)))
        }
    }
    fn keyword(&mut self) -> String {
        self.skip_space();
        let start = self.i;
        while self
            .peek()
            .is_some_and(|c| c.is_ascii_alphanumeric() || c == '_')
        {
            self.i += 1;
        }
        self.s[start..self.i].to_ascii_uppercase()
    }
    /// `TYPE(args)` or `( TYPE(args) TYPE(args) ... )`.
    fn entity_parts(&mut self) -> Result<Vec<(String, Vec<Param>)>, ImportRefusal> {
        self.skip_space();
        if self.peek() == Some('(') {
            self.i += 1;
            let mut parts = Vec::new();
            loop {
                self.skip_space();
                if self.peek() == Some(')') {
                    self.i += 1;
                    break;
                }
                let name = self.keyword();
                if name.is_empty() {
                    return Err(malformed(format!("bad complex entity `{}`", self.s)));
                }
                parts.push((name, self.args()?));
            }
            if parts.is_empty() {
                return Err(malformed("empty complex entity"));
            }
            Ok(parts)
        } else {
            let name = self.keyword();
            if name.is_empty() {
                return Err(malformed(format!("missing entity type in `{}`", self.s)));
            }
            Ok(vec![(name, self.args()?)])
        }
    }
    fn args(&mut self) -> Result<Vec<Param>, ImportRefusal> {
        self.expect('(')?;
        self.list_tail()
    }
    /// Parses parameters up to and including the closing ')'.
    fn list_tail(&mut self) -> Result<Vec<Param>, ImportRefusal> {
        let mut out = Vec::new();
        self.skip_space();
        if self.peek() == Some(')') {
            self.i += 1;
            return Ok(out);
        }
        loop {
            out.push(self.param()?);
            self.skip_space();
            match self.peek() {
                Some(',') => self.i += 1,
                Some(')') => {
                    self.i += 1;
                    return Ok(out);
                }
                _ => return Err(malformed(format!("bad parameter list in `{}`", self.s))),
            }
        }
    }
    fn param(&mut self) -> Result<Param, ImportRefusal> {
        self.skip_space();
        match self.peek() {
            Some('#') => {
                self.i += 1;
                let start = self.i;
                while self.peek().is_some_and(|c| c.is_ascii_digit()) {
                    self.i += 1;
                }
                self.s[start..self.i]
                    .parse()
                    .map(Param::Ref)
                    .map_err(|_| malformed("bad entity reference"))
            }
            Some('$') => {
                self.i += 1;
                Ok(Param::Unset)
            }
            Some('*') => {
                self.i += 1;
                Ok(Param::Derived)
            }
            Some('\'') => {
                self.i += 1;
                let mut out = String::new();
                loop {
                    match self.peek() {
                        None => return Err(malformed("unterminated string")),
                        Some('\'') => {
                            self.i += 1;
                            if self.peek() == Some('\'') {
                                self.i += 1;
                                out.push('\'');
                            } else {
                                return Ok(Param::Str(out));
                            }
                        }
                        Some(c) => {
                            self.i += c.len_utf8();
                            out.push(c);
                        }
                    }
                }
            }
            Some('.') => {
                self.i += 1;
                let start = self.i;
                while self.peek().is_some_and(|c| c != '.') {
                    self.i += 1;
                }
                let word = self.s[start..self.i].to_ascii_uppercase();
                self.expect('.')?;
                Ok(Param::Enum(word))
            }
            Some('(') => {
                self.i += 1;
                Ok(Param::List(self.list_tail()?))
            }
            Some(c) if c.is_ascii_digit() || c == '-' || c == '+' => {
                let start = self.i;
                self.i += 1;
                while self
                    .peek()
                    .is_some_and(|c| c.is_ascii_digit() || matches!(c, '.' | 'e' | 'E' | '-' | '+'))
                {
                    self.i += 1;
                }
                Ok(Param::Number(self.s[start..self.i].to_string()))
            }
            Some(c) if c.is_ascii_alphabetic() => {
                let name = self.keyword();
                Ok(Param::Typed(name, self.args()?))
            }
            _ => Err(malformed(format!("unreadable parameter in `{}`", self.s))),
        }
    }
}

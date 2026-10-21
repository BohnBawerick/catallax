//! The canonical JSON layout.
//!
//! Two-space indentation. Each object member, and each element of an array of objects
//! or arrays, goes on its own line. An array of strings or numbers stays on one line,
//! as `["1", "0"]`. `": "` separates a key from its value. Empty arrays and objects are
//! `[]` and `{}`.

use std::io;

use serde_json::ser::Formatter;

enum Frame {
    Object { members: bool },
    Array { layout: Layout },
}

#[derive(Clone, Copy, PartialEq, Eq)]
enum Layout {
    /// No element seen yet.
    Unknown,
    /// Elements are strings or numbers, on one line.
    Inline,
    /// Elements are objects or arrays, one per line.
    Lines,
}

/// A formatter for the canonical layout.
pub(crate) struct Canonical {
    indent: usize,
    stack: Vec<Frame>,
}

impl Canonical {
    pub(crate) fn new() -> Self {
        Canonical {
            indent: 0,
            stack: Vec::new(),
        }
    }

    fn newline<W: ?Sized + io::Write>(&self, w: &mut W) -> io::Result<()> {
        w.write_all(b"\n")?;
        for _ in 0..self.indent {
            w.write_all(b"  ")?;
        }
        Ok(())
    }

    /// Called before every value. The first element of an array sets the array's layout.
    fn value<W: ?Sized + io::Write>(&mut self, w: &mut W, container: bool) -> io::Result<()> {
        if let Some(Frame::Array { layout }) = self.stack.last_mut()
            && *layout == Layout::Unknown
        {
            if container {
                *layout = Layout::Lines;
                self.indent += 1;
                self.newline(w)?;
            } else {
                *layout = Layout::Inline;
            }
        }
        Ok(())
    }
}

macro_rules! scalar {
    ($($method:ident: $ty:ty),* $(,)?) => {
        $(
            fn $method<W: ?Sized + io::Write>(&mut self, w: &mut W, value: $ty) -> io::Result<()> {
                self.value(w, false)?;
                serde_json::ser::CompactFormatter.$method(w, value)
            }
        )*
    };
}

impl Formatter for Canonical {
    scalar!(
        write_bool: bool,
        write_i8: i8,
        write_i16: i16,
        write_i32: i32,
        write_i64: i64,
        write_i128: i128,
        write_u8: u8,
        write_u16: u16,
        write_u32: u32,
        write_u64: u64,
        write_u128: u128,
        write_f32: f32,
        write_f64: f64,
        write_number_str: &str,
    );

    fn write_null<W: ?Sized + io::Write>(&mut self, w: &mut W) -> io::Result<()> {
        self.value(w, false)?;
        w.write_all(b"null")
    }

    fn begin_string<W: ?Sized + io::Write>(&mut self, w: &mut W) -> io::Result<()> {
        self.value(w, false)?;
        w.write_all(b"\"")
    }

    fn begin_array<W: ?Sized + io::Write>(&mut self, w: &mut W) -> io::Result<()> {
        self.value(w, true)?;
        self.stack.push(Frame::Array {
            layout: Layout::Unknown,
        });
        w.write_all(b"[")
    }

    fn end_array<W: ?Sized + io::Write>(&mut self, w: &mut W) -> io::Result<()> {
        if let Some(Frame::Array {
            layout: Layout::Lines,
        }) = self.stack.pop()
        {
            self.indent -= 1;
            self.newline(w)?;
        }
        w.write_all(b"]")
    }

    fn begin_array_value<W: ?Sized + io::Write>(
        &mut self,
        w: &mut W,
        first: bool,
    ) -> io::Result<()> {
        if first {
            return Ok(());
        }
        match self.stack.last() {
            Some(Frame::Array {
                layout: Layout::Lines,
            }) => {
                w.write_all(b",")?;
                self.newline(w)
            }
            _ => w.write_all(b", "),
        }
    }

    fn end_array_value<W: ?Sized + io::Write>(&mut self, _w: &mut W) -> io::Result<()> {
        Ok(())
    }

    fn begin_object<W: ?Sized + io::Write>(&mut self, w: &mut W) -> io::Result<()> {
        self.value(w, true)?;
        self.stack.push(Frame::Object { members: false });
        self.indent += 1;
        w.write_all(b"{")
    }

    fn end_object<W: ?Sized + io::Write>(&mut self, w: &mut W) -> io::Result<()> {
        self.indent -= 1;
        if let Some(Frame::Object { members: true }) = self.stack.pop() {
            self.newline(w)?;
        }
        w.write_all(b"}")
    }

    fn begin_object_key<W: ?Sized + io::Write>(
        &mut self,
        w: &mut W,
        first: bool,
    ) -> io::Result<()> {
        if let Some(Frame::Object { members }) = self.stack.last_mut() {
            *members = true;
        }
        if !first {
            w.write_all(b",")?;
        }
        self.newline(w)
    }

    fn begin_object_value<W: ?Sized + io::Write>(&mut self, w: &mut W) -> io::Result<()> {
        w.write_all(b": ")
    }

    fn end_object_value<W: ?Sized + io::Write>(&mut self, _w: &mut W) -> io::Result<()> {
        Ok(())
    }
}

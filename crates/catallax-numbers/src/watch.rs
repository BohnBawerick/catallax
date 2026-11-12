//! A list of what work on a copy of a store read.

use std::fmt;
use std::sync::{Arc, Mutex, PoisonError};

/// What a store has been asked since [`Watch::begin`], while a list is kept.
///
/// A store that only grows can still change its answers: a fact gets
/// recorded, a table entry gains a member. Work that ran on a copy is good
/// for the store the copy came from only while that store still answers what
/// the work asked as the copy did. A store keeps one of these, writes each
/// such question into it, and checks the list against the original afterwards.
///
/// Copies of a watched store share the list, since work on a copy of a
/// watched copy is the same work. Writing the list spends no decision work,
/// and nothing reads it before the work is over.
pub struct Watch<R> {
    list: Option<Arc<Mutex<Vec<R>>>>,
}

impl<R> Default for Watch<R> {
    fn default() -> Self {
        Self { list: None }
    }
}

impl<R> Clone for Watch<R> {
    fn clone(&self) -> Self {
        Self {
            list: self.list.clone(),
        }
    }
}

impl<R> fmt::Debug for Watch<R> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(if self.list.is_some() {
            "Watch(on)"
        } else {
            "Watch(off)"
        })
    }
}

impl<R: Ord> Watch<R> {
    /// Writes a read into the list, while one is kept.
    pub fn note(&self, read: impl FnOnce() -> R) {
        if let Some(list) = &self.list {
            list.lock()
                .unwrap_or_else(PoisonError::into_inner)
                .push(read());
        }
    }

    /// Starts a list.
    pub fn begin(&mut self) {
        self.list = Some(Arc::default());
    }

    /// Ends the list and returns it, each read once. Panics when none was
    /// begun.
    pub fn end(&mut self) -> Vec<R> {
        let list = self.list.take().expect("a watch was begun");
        let mut reads = std::mem::take(&mut *list.lock().unwrap_or_else(PoisonError::into_inner));
        reads.sort_unstable();
        reads.dedup();
        reads
    }
}

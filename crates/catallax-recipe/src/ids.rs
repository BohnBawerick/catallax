//! Step ids and element ids.

use std::fmt;

use serde::{Deserialize, Serialize};

/// Why a text is not an id.
#[derive(Clone, PartialEq, Eq, Debug)]
pub struct IdError {
    /// The rejected text.
    pub text: String,
}

impl fmt::Display for IdError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{:?} is not an id: use 1 to 64 characters from a-z 0-9 _ -",
            self.text
        )
    }
}

impl std::error::Error for IdError {}

/// The id grammar of the exchange contract's step names: 1 to 64 characters from
/// `a-z 0-9 _ -`.
fn is_id(text: &str) -> bool {
    (1..=64).contains(&text.len())
        && text
            .bytes()
            .all(|c| c.is_ascii_lowercase() || c.is_ascii_digit() || c == b'_' || c == b'-')
}

macro_rules! id_type {
    ($(#[$doc:meta])* $name:ident) => {
        $(#[$doc])*
        #[derive(Clone, PartialEq, Eq, PartialOrd, Ord, Hash, Debug, Serialize, Deserialize)]
        #[serde(try_from = "String", into = "String")]
        pub struct $name(String);

        impl $name {
            /// Checks `text` against the id grammar.
            pub fn new(text: impl Into<String>) -> Result<Self, IdError> {
                let text = text.into();
                if is_id(&text) {
                    Ok($name(text))
                } else {
                    Err(IdError { text })
                }
            }

            /// The id text.
            pub fn as_str(&self) -> &str {
                &self.0
            }
        }

        impl TryFrom<String> for $name {
            type Error = IdError;
            fn try_from(text: String) -> Result<Self, IdError> {
                $name::new(text)
            }
        }

        impl From<$name> for String {
            fn from(id: $name) -> String {
                id.0
            }
        }

        impl fmt::Display for $name {
            fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
                f.write_str(&self.0)
            }
        }
    };
}

id_type!(
    /// The id of a step. Later steps name its value by it, and every name of a face
    /// the step made starts with it. It never changes when the step's parameters do.
    StepId
);

id_type!(
    /// The id of an element inside a step: a sketch segment, a path segment or a hole
    /// placement. Unique within its step. Faces made from the element are named by it.
    ElementId
);

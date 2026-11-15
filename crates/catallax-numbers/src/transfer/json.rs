use super::{Codec, Decider, Refusal, Writer, invalid};
use serde::{Serialize, Serializer, ser};
use std::{
    cell::{Cell, RefCell},
    io,
};

thread_local! {
    static NUMBER_METER: Cell<bool> = const { Cell::new(false) };
}

struct NumberMeter(bool);
impl NumberMeter {
    fn enter() -> Self {
        Self(NUMBER_METER.replace(true))
    }
}
impl Drop for NumberMeter {
    fn drop(&mut self) {
        NUMBER_METER.set(self.0);
    }
}

// Recipe numbers allocate their canonical text before entering serialize_str.
// Charge their existing bit lengths before that conversion, only during transfer.
pub fn charge_json_number(value: &crate::Rational) -> Result<(), Refusal> {
    if !NUMBER_METER.get() {
        return Ok(());
    }
    let units = value.n.0.bits().saturating_add(value.d.0.bits());
    if crate::operation::charge(units) {
        return Err(Refusal::new(
            "numbers.operation-budget",
            "number formatting exhausted the transfer allowance",
        ));
    }
    Ok(())
}

struct Budget<'a> {
    decider: RefCell<&'a mut Decider>,
    refusal: RefCell<Option<Refusal>>,
}
impl Budget<'_> {
    fn charge(&self, units: usize) -> Result<(), Refusal> {
        self.decider
            .borrow_mut()
            .charge_work(units as u64)
            .inspect_err(|error| {
                self.refusal.replace(Some(error.clone()));
            })
    }
    fn serde<E: ser::Error>(&self, units: usize) -> Result<(), E> {
        self.charge(units).map_err(E::custom)
    }
}

struct Output<'a, 'b> {
    bytes: &'a mut Vec<u8>,
    budget: &'a Budget<'b>,
}
impl io::Write for Output<'_, '_> {
    fn write(&mut self, bytes: &[u8]) -> io::Result<usize> {
        self.budget.charge(bytes.len()).map_err(io::Error::other)?;
        if bytes.len() > self.bytes.capacity().saturating_sub(self.bytes.len()) {
            self.budget
                .charge(self.bytes.len())
                .map_err(io::Error::other)?;
        }
        self.bytes
            .try_reserve(bytes.len())
            .map_err(io::Error::other)?;
        self.budget
            .decider
            .borrow()
            .check_deadline()
            .map_err(io::Error::other)?;
        self.bytes.extend_from_slice(bytes);
        Ok(bytes.len())
    }
    fn flush(&mut self) -> io::Result<()> {
        Ok(())
    }
}

impl Writer<'_> {
    pub fn json<T: Serialize + ?Sized>(&mut self, value: &T) -> Result<(), Refusal> {
        self.charge(1)?;
        let mut bytes = Vec::new();
        let (result, refusal) = {
            let _number_meter = NumberMeter::enter();
            let budget = Budget {
                decider: RefCell::new(&mut *self.decider),
                refusal: RefCell::new(None),
            };
            let result = {
                let mut output = Output {
                    bytes: &mut bytes,
                    budget: &budget,
                };
                let mut serializer = serde_json::Serializer::new(&mut output);
                value.serialize(Metered {
                    inner: &mut serializer,
                    budget: &budget,
                })
            };
            (result, budget.refusal.into_inner())
        };
        if result.is_err() {
            if let Some(refusal) = refusal {
                return Err(refusal);
            }
            self.decider.charge_work(0)?;
            return Err(invalid("metadata encoding failed"));
        }
        bytes.len().write(self)?;
        self.raw(&bytes)
    }
}

struct Metered<'a, 'b, S> {
    inner: S,
    budget: &'a Budget<'b>,
}
struct Value<'a, 'b, T: ?Sized> {
    inner: &'a T,
    budget: &'a Budget<'b>,
}

struct Text<'a, 'b> {
    value: String,
    budget: &'a Budget<'b>,
}
impl std::fmt::Write for Text<'_, '_> {
    fn write_str(&mut self, text: &str) -> std::fmt::Result {
        self.budget
            .charge(text.len())
            .map_err(|_| std::fmt::Error)?;
        if text.len() > self.value.capacity().saturating_sub(self.value.len()) {
            self.budget
                .charge(self.value.len())
                .map_err(|_| std::fmt::Error)?;
            self.value
                .try_reserve(text.len())
                .map_err(|_| std::fmt::Error)?;
        }
        self.budget
            .decider
            .borrow()
            .check_deadline()
            .map_err(|_| std::fmt::Error)?;
        self.value.push_str(text);
        Ok(())
    }
}
impl<T: Serialize + ?Sized> Serialize for Value<'_, '_, T> {
    fn serialize<S: Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
        self.budget.serde(1)?;
        self.inner.serialize(Metered {
            inner: serializer,
            budget: self.budget,
        })
    }
}

macro_rules! scalar {
    ($($method:ident: $ty:ty),* $(,)?) => {$ (
        fn $method(self, value: $ty) -> Result<Self::Ok, Self::Error> {
            self.budget.serde(1)?;
            self.inner.$method(value)
        }
    )*};
}
impl<'a, 'b, S: Serializer> Serializer for Metered<'a, 'b, S> {
    type Ok = S::Ok;
    type Error = S::Error;
    type SerializeSeq = Metered<'a, 'b, S::SerializeSeq>;
    type SerializeTuple = Metered<'a, 'b, S::SerializeTuple>;
    type SerializeTupleStruct = Metered<'a, 'b, S::SerializeTupleStruct>;
    type SerializeTupleVariant = Metered<'a, 'b, S::SerializeTupleVariant>;
    type SerializeMap = Metered<'a, 'b, S::SerializeMap>;
    type SerializeStruct = Metered<'a, 'b, S::SerializeStruct>;
    type SerializeStructVariant = Metered<'a, 'b, S::SerializeStructVariant>;
    scalar!(serialize_bool: bool, serialize_i8: i8, serialize_i16: i16, serialize_i32: i32,
        serialize_i64: i64, serialize_i128: i128, serialize_u8: u8, serialize_u16: u16,
        serialize_u32: u32, serialize_u64: u64, serialize_u128: u128, serialize_f32: f32,
        serialize_f64: f64, serialize_char: char);
    fn serialize_str(self, value: &str) -> Result<Self::Ok, Self::Error> {
        self.budget.serde(value.len())?;
        self.inner.serialize_str(value)
    }
    fn serialize_bytes(self, value: &[u8]) -> Result<Self::Ok, Self::Error> {
        self.budget.serde(value.len())?;
        self.inner.serialize_bytes(value)
    }
    fn serialize_none(self) -> Result<Self::Ok, Self::Error> {
        self.budget.serde(1)?;
        self.inner.serialize_none()
    }
    fn serialize_some<T: Serialize + ?Sized>(self, value: &T) -> Result<Self::Ok, Self::Error> {
        self.budget.serde(1)?;
        self.inner.serialize_some(&Value {
            inner: value,
            budget: self.budget,
        })
    }
    fn serialize_unit(self) -> Result<Self::Ok, Self::Error> {
        self.budget.serde(1)?;
        self.inner.serialize_unit()
    }
    fn serialize_unit_struct(self, name: &'static str) -> Result<Self::Ok, Self::Error> {
        self.budget.serde(1)?;
        self.inner.serialize_unit_struct(name)
    }
    fn serialize_unit_variant(
        self,
        name: &'static str,
        index: u32,
        variant: &'static str,
    ) -> Result<Self::Ok, Self::Error> {
        self.budget.serde(variant.len())?;
        self.inner.serialize_unit_variant(name, index, variant)
    }
    fn serialize_newtype_struct<T: Serialize + ?Sized>(
        self,
        name: &'static str,
        value: &T,
    ) -> Result<Self::Ok, Self::Error> {
        self.budget.serde(1)?;
        self.inner.serialize_newtype_struct(
            name,
            &Value {
                inner: value,
                budget: self.budget,
            },
        )
    }
    fn serialize_newtype_variant<T: Serialize + ?Sized>(
        self,
        name: &'static str,
        index: u32,
        variant: &'static str,
        value: &T,
    ) -> Result<Self::Ok, Self::Error> {
        self.budget.serde(variant.len())?;
        self.inner.serialize_newtype_variant(
            name,
            index,
            variant,
            &Value {
                inner: value,
                budget: self.budget,
            },
        )
    }
    fn serialize_seq(self, len: Option<usize>) -> Result<Self::SerializeSeq, Self::Error> {
        self.budget.serde(len.unwrap_or(1))?;
        Ok(Metered {
            inner: self.inner.serialize_seq(len)?,
            budget: self.budget,
        })
    }
    fn serialize_tuple(self, len: usize) -> Result<Self::SerializeTuple, Self::Error> {
        self.budget.serde(len)?;
        Ok(Metered {
            inner: self.inner.serialize_tuple(len)?,
            budget: self.budget,
        })
    }
    fn serialize_tuple_struct(
        self,
        name: &'static str,
        len: usize,
    ) -> Result<Self::SerializeTupleStruct, Self::Error> {
        self.budget.serde(len)?;
        Ok(Metered {
            inner: self.inner.serialize_tuple_struct(name, len)?,
            budget: self.budget,
        })
    }
    fn serialize_tuple_variant(
        self,
        name: &'static str,
        index: u32,
        variant: &'static str,
        len: usize,
    ) -> Result<Self::SerializeTupleVariant, Self::Error> {
        self.budget.serde(len.saturating_add(variant.len()))?;
        Ok(Metered {
            inner: self
                .inner
                .serialize_tuple_variant(name, index, variant, len)?,
            budget: self.budget,
        })
    }
    fn serialize_map(self, len: Option<usize>) -> Result<Self::SerializeMap, Self::Error> {
        self.budget.serde(len.unwrap_or(1))?;
        Ok(Metered {
            inner: self.inner.serialize_map(len)?,
            budget: self.budget,
        })
    }
    fn serialize_struct(
        self,
        name: &'static str,
        len: usize,
    ) -> Result<Self::SerializeStruct, Self::Error> {
        self.budget.serde(len)?;
        Ok(Metered {
            inner: self.inner.serialize_struct(name, len)?,
            budget: self.budget,
        })
    }
    fn serialize_struct_variant(
        self,
        name: &'static str,
        index: u32,
        variant: &'static str,
        len: usize,
    ) -> Result<Self::SerializeStructVariant, Self::Error> {
        self.budget.serde(len.saturating_add(variant.len()))?;
        Ok(Metered {
            inner: self
                .inner
                .serialize_struct_variant(name, index, variant, len)?,
            budget: self.budget,
        })
    }
    fn is_human_readable(&self) -> bool {
        self.inner.is_human_readable()
    }
    fn collect_str<T: std::fmt::Display + ?Sized>(
        self,
        value: &T,
    ) -> Result<Self::Ok, Self::Error> {
        self.budget.serde(1)?;
        let mut text = Text {
            value: String::new(),
            budget: self.budget,
        };
        std::fmt::write(&mut text, format_args!("{value}")).map_err(ser::Error::custom)?;
        self.serialize_str(&text.value)
    }
}

macro_rules! sequence {
    ($trait:ident, $method:ident) => {
        impl<S: ser::$trait> ser::$trait for Metered<'_, '_, S> {
            type Ok = S::Ok;
            type Error = S::Error;
            fn $method<T: Serialize + ?Sized>(&mut self, value: &T) -> Result<(), Self::Error> {
                self.budget.serde(1)?;
                self.inner.$method(&Value {
                    inner: value,
                    budget: self.budget,
                })
            }
            fn end(self) -> Result<Self::Ok, Self::Error> {
                self.budget.serde(1)?;
                self.inner.end()
            }
        }
    };
}
sequence!(SerializeSeq, serialize_element);
sequence!(SerializeTuple, serialize_element);
sequence!(SerializeTupleStruct, serialize_field);
sequence!(SerializeTupleVariant, serialize_field);

macro_rules! structure {
    ($trait:ident) => {
        impl<S: ser::$trait> ser::$trait for Metered<'_, '_, S> {
            type Ok = S::Ok;
            type Error = S::Error;
            fn serialize_field<T: Serialize + ?Sized>(
                &mut self,
                key: &'static str,
                value: &T,
            ) -> Result<(), Self::Error> {
                self.budget.serde(key.len())?;
                self.inner.serialize_field(
                    key,
                    &Value {
                        inner: value,
                        budget: self.budget,
                    },
                )
            }
            fn skip_field(&mut self, key: &'static str) -> Result<(), Self::Error> {
                self.budget.serde(1)?;
                self.inner.skip_field(key)
            }
            fn end(self) -> Result<Self::Ok, Self::Error> {
                self.budget.serde(1)?;
                self.inner.end()
            }
        }
    };
}
structure!(SerializeStruct);
structure!(SerializeStructVariant);
impl<S: ser::SerializeMap> ser::SerializeMap for Metered<'_, '_, S> {
    type Ok = S::Ok;
    type Error = S::Error;
    fn serialize_key<T: Serialize + ?Sized>(&mut self, key: &T) -> Result<(), Self::Error> {
        self.budget.serde(1)?;
        self.inner.serialize_key(&Value {
            inner: key,
            budget: self.budget,
        })
    }
    fn serialize_value<T: Serialize + ?Sized>(&mut self, value: &T) -> Result<(), Self::Error> {
        self.budget.serde(1)?;
        self.inner.serialize_value(&Value {
            inner: value,
            budget: self.budget,
        })
    }
    fn end(self) -> Result<Self::Ok, Self::Error> {
        self.budget.serde(1)?;
        self.inner.end()
    }
}

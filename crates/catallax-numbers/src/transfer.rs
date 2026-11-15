//! Internal codecs for body handover. The application module owns the public API.
//!
//! Shared objects are numbered, written once, then restored in dependency order.
//! Expression depth therefore does not become decoder call-stack depth.

#![allow(missing_docs)]

use crate::{Decider, Refusal};
use sha2::{Digest, Sha256};
use std::{
    any::{Any, TypeId, type_name},
    cell::RefCell,
    collections::{BTreeMap, BTreeSet, HashMap, HashSet},
    sync::{Arc, Mutex, OnceLock},
};

const MAGIC: &[u8; 8] = b"CXBODY\0\0";
const VERSION: u32 = 1;

pub fn invalid(detail: &str) -> Refusal {
    Refusal::new("handover.invalid-data", detail)
}

pub trait Codec: Sized + 'static {
    fn write(&self, writer: &mut Writer) -> Result<(), Refusal>;
    fn read(reader: &mut Reader<'_>) -> Result<Self, Refusal>;
    fn register(_registry: &mut Registry) {}
}

/// Retained metadata already has a lossless JSON representation.
#[doc(hidden)]
#[macro_export]
macro_rules! transfer_json {
    ($($ty:ty),* $(,)?) => {$ (
        impl $crate::transfer::Codec for $ty {
            fn write(&self, writer: &mut $crate::transfer::Writer) -> Result<(), $crate::Refusal> {
                let bytes = serde_json::to_vec(self).map_err(|_| $crate::transfer::invalid("metadata encoding failed"))?;
                $crate::transfer::Codec::write(&bytes, writer)
            }
            fn read(reader: &mut $crate::transfer::Reader<'_>) -> Result<Self, $crate::Refusal> {
                let bytes = <Vec<u8> as $crate::transfer::Codec>::read(reader)?;
                serde_json::from_slice(&bytes).map_err(|_| $crate::transfer::invalid("invalid retained metadata"))
            }
        }
    )*};
}

type Factory = fn(&mut Reader<'_>) -> Result<Box<dyn Any>, Refusal>;
type Emit = Box<dyn FnOnce(&mut Writer) -> Result<(), Refusal>>;
type HandleRanges = RefCell<HashMap<TypeId, (Option<usize>, Option<u32>)>>;

#[derive(Default)]
pub struct Registry {
    visited: HashSet<TypeId>,
    factories: HashMap<&'static str, Factory>,
}

impl Registry {
    pub fn add<T: Codec>(&mut self) {
        if self.visited.insert(TypeId::of::<T>()) {
            T::register(self);
        }
    }
}

#[derive(Default)]
pub struct Writer {
    bytes: Vec<u8>,
    keys: HashMap<(TypeId, usize), u32>,
    pending: Vec<Option<(&'static str, Emit)>>,
    dependencies: Vec<Vec<usize>>,
    active: Option<usize>,
}

impl Writer {
    fn shared<T: 'static>(&mut self, pointer: usize, emit: Emit) -> Result<(), Refusal> {
        let key = (TypeId::of::<T>(), pointer);
        let id = match self.keys.get(&key) {
            Some(id) => *id,
            None => {
                let id = u32::try_from(self.pending.len())
                    .map_err(|_| invalid("too many shared objects"))?;
                self.keys.insert(key, id);
                self.pending.push(Some((type_name::<T>(), emit)));
                self.dependencies.push(Vec::new());
                id
            }
        };
        if let Some(active) = self.active {
            self.dependencies[active].push(id as usize);
        }
        id.write(self)
    }

    pub fn raw(&mut self, bytes: &[u8]) {
        self.bytes.extend_from_slice(bytes);
    }
}

pub struct Reader<'a> {
    data: &'a [u8],
    at: usize,
    objects: &'a HashMap<u32, Box<dyn Any>>,
    decider: &'a mut Decider,
    handles: &'a HandleRanges,
}

impl Reader<'_> {
    pub fn handle<T: 'static>(&self, index: u32) {
        let mut handles = self.handles.borrow_mut();
        let entry = handles.entry(TypeId::of::<T>()).or_default();
        entry.1 = Some(entry.1.map_or(index, |old| old.max(index)));
    }

    pub fn arena<T: 'static>(&self, length: usize) -> Result<(), Refusal> {
        let mut handles = self.handles.borrow_mut();
        let entry = handles.entry(TypeId::of::<T>()).or_default();
        if entry.0.replace(length).is_some() {
            return Err(invalid("multiple arenas of the same kind"));
        }
        Ok(())
    }
    pub fn take(&mut self, count: usize) -> Result<&[u8], Refusal> {
        let end = self
            .at
            .checked_add(count)
            .ok_or_else(|| invalid("length overflow"))?;
        let bytes = self
            .data
            .get(self.at..end)
            .ok_or_else(|| invalid("truncated record"))?;
        self.decider.charge_work(count as u64)?;
        self.at = end;
        Ok(bytes)
    }

    pub fn count(&mut self) -> Result<usize, Refusal> {
        let count = usize::read(self)?;
        if count > self.data.len().saturating_sub(self.at) {
            return Err(invalid("record length exceeds input"));
        }
        Ok(count)
    }

    fn finish(&self) -> Result<(), Refusal> {
        if self.at == self.data.len() {
            Ok(())
        } else {
            Err(invalid("trailing record bytes"))
        }
    }
}

pub fn encode<T: Codec>(body: &T, decider: &mut Decider) -> Result<Vec<u8>, Refusal> {
    let mut writer = Writer::default();
    body.write(&mut writer)?;
    let root = std::mem::take(&mut writer.bytes);
    let mut records = Vec::new();
    let mut i = 0;
    while i < writer.pending.len() {
        decider.charge_work(1)?;
        writer.active = Some(i);
        let (tag, emit) = writer.pending[i]
            .take()
            .ok_or_else(|| invalid("missing shared object"))?;
        emit(&mut writer)?;
        records.push((tag, std::mem::take(&mut writer.bytes)));
        i += 1;
    }
    // Iterative postorder preserves graph sharing without recursive traversal.
    let mut order = Vec::new();
    let mut state = vec![0_u8; records.len()];
    for start in 0..records.len() {
        let mut stack = vec![(start, false)];
        while let Some((id, exit)) = stack.pop() {
            if exit {
                state[id] = 2;
                order.push(id);
            } else if state[id] != 2 {
                if state[id] == 1 {
                    return Err(invalid("cyclic shared definitions"));
                }
                state[id] = 1;
                stack.push((id, true));
                for &child in writer.dependencies[id].iter().rev() {
                    if state[child] != 2 {
                        stack.push((child, false));
                    }
                }
            }
        }
    }
    let mut output = Writer::default();
    output.raw(MAGIC);
    VERSION.write(&mut output)?;
    type_name::<T>().to_string().write(&mut output)?;
    let tags: BTreeSet<_> = records.iter().map(|r| r.0).collect();
    let tags: Vec<_> = tags.into_iter().collect();
    tags.len().write(&mut output)?;
    for tag in &tags {
        tag.to_string().write(&mut output)?;
    }
    records.len().write(&mut output)?;
    for id in order {
        (id as u32).write(&mut output)?;
        tags.binary_search(&records[id].0)
            .map_err(|_| invalid("missing record type"))?
            .write(&mut output)?;
        records[id].1.len().write(&mut output)?;
        output.raw(&records[id].1);
    }
    root.len().write(&mut output)?;
    output.raw(&root);
    decider.charge_work(output.bytes.len() as u64)?;
    let digest = Sha256::digest(&output.bytes);
    output.raw(&digest);
    Ok(output.bytes)
}

pub fn decode<T: Codec>(bytes: &[u8], decider: &mut Decider) -> Result<T, Refusal> {
    decider.charge_work(bytes.len() as u64)?;
    if bytes.get(..8) != Some(MAGIC.as_slice()) {
        return Err(invalid("body handover magic is missing"));
    }
    let version = bytes
        .get(8..12)
        .ok_or_else(|| invalid("truncated version"))?;
    if version != VERSION.to_le_bytes() {
        return Err(Refusal::new(
            "handover.version-mismatch",
            "unsupported body handover format version",
        ));
    }
    let end = bytes
        .len()
        .checked_sub(32)
        .filter(|end| *end >= 12)
        .ok_or_else(|| invalid("truncated digest"))?;
    if Sha256::digest(&bytes[..end]).as_slice() != &bytes[end..] {
        return Err(invalid("body handover checksum mismatch"));
    }
    let mut registry = Registry::default();
    registry.add::<T>();
    let mut objects = HashMap::new();
    let handles = RefCell::new(HashMap::new());
    let mut reader = Reader {
        data: &bytes[12..end],
        at: 0,
        objects: &objects,
        decider,
        handles: &handles,
    };
    if String::read(&mut reader)? != type_name::<T>() {
        return Err(invalid("wrong body handover kind"));
    }
    let count = reader.count()?;
    let mut factories = Vec::new();
    for _ in 0..count {
        let tag = String::read(&mut reader)?;
        factories.push(
            *registry
                .factories
                .get(tag.as_str())
                .ok_or_else(|| invalid("unknown shared record type"))?,
        );
    }
    let count = reader.count()?;
    let mut at = reader.at;
    for _ in 0..count {
        let mut reader = Reader {
            data: &bytes[12..end],
            at,
            objects: &objects,
            decider,
            handles: &handles,
        };
        let id = u32::read(&mut reader)?;
        if id as usize >= count || objects.contains_key(&id) {
            return Err(invalid("invalid shared object number"));
        }
        let kind = usize::read(&mut reader)?;
        let factory = factories
            .get(kind)
            .ok_or_else(|| invalid("invalid shared record type"))?;
        let length = reader.count()?;
        let start = reader.at;
        reader.take(length)?;
        at = reader.at;
        let mut record = Reader {
            data: &bytes[12 + start..12 + start + length],
            at: 0,
            objects: &objects,
            decider,
            handles: &handles,
        };
        let value = factory(&mut record)?;
        record.finish()?;
        objects.insert(id, value);
    }
    let mut reader = Reader {
        data: &bytes[12..end],
        at,
        objects: &objects,
        decider,
        handles: &handles,
    };
    let length = reader.count()?;
    let start = reader.at;
    reader.take(length)?;
    reader.finish()?;
    let mut root = Reader {
        data: &bytes[12 + start..12 + start + length],
        at: 0,
        objects: &objects,
        decider,
        handles: &handles,
    };
    let value = T::read(&mut root)?;
    root.finish()?;
    for (length, index) in handles.borrow().values() {
        if let Some(index) = index
            && length.is_none_or(|length| *index as usize >= length)
        {
            return Err(invalid("handle outside its arena"));
        }
    }
    Ok(value)
}

macro_rules! integer {
    ($($ty:ty),*) => {$ (
        impl Codec for $ty {
            fn write(&self, w: &mut Writer) -> Result<(), Refusal> { w.raw(&self.to_le_bytes()); Ok(()) }
            fn read(r: &mut Reader<'_>) -> Result<Self, Refusal> {
                let bytes = r.take(std::mem::size_of::<Self>())?;
                Ok(Self::from_le_bytes(bytes.try_into().map_err(|_| invalid("invalid scalar width"))?))
            }
        }
    )*};
}
integer!(i8, u8, u16, u32, u64, i32, i64, i128, u128);

impl Codec for usize {
    fn write(&self, w: &mut Writer) -> Result<(), Refusal> {
        (*self as u64).write(w)
    }
    fn read(r: &mut Reader<'_>) -> Result<Self, Refusal> {
        Self::try_from(u64::read(r)?).map_err(|_| invalid("index exceeds platform width"))
    }
}
impl Codec for num_bigint::BigInt {
    fn write(&self, w: &mut Writer) -> Result<(), Refusal> {
        let bytes = self.to_signed_bytes_le();
        bytes.write(w)
    }
    fn read(r: &mut Reader<'_>) -> Result<Self, Refusal> {
        Ok(Self::from_signed_bytes_le(&Vec::<u8>::read(r)?))
    }
}
impl Codec for crate::Rational {
    fn write(&self, w: &mut Writer) -> Result<(), Refusal> {
        self.n.write(w)?;
        self.d.write(w)
    }
    fn read(r: &mut Reader<'_>) -> Result<Self, Refusal> {
        let n = crate::Integer::read(r)?;
        let d = crate::Integer::read(r)?;
        if d.sign() != crate::Sign::Positive {
            return Err(invalid("non-positive rational denominator"));
        }
        let value =
            Self::new(n.clone(), d.clone()).ok_or_else(|| invalid("zero rational denominator"))?;
        if value.n != n || value.d != d {
            return Err(invalid("noncanonical rational"));
        }
        Ok(value)
    }
}
impl Codec for f64 {
    fn write(&self, w: &mut Writer) -> Result<(), Refusal> {
        self.to_bits().write(w)
    }
    fn read(r: &mut Reader<'_>) -> Result<Self, Refusal> {
        Ok(Self::from_bits(u64::read(r)?))
    }
}
impl Codec for bool {
    fn write(&self, w: &mut Writer) -> Result<(), Refusal> {
        u8::from(*self).write(w)
    }
    fn read(r: &mut Reader<'_>) -> Result<Self, Refusal> {
        match u8::read(r)? {
            0 => Ok(false),
            1 => Ok(true),
            _ => Err(invalid("invalid boolean")),
        }
    }
}
impl Codec for String {
    fn write(&self, w: &mut Writer) -> Result<(), Refusal> {
        self.len().write(w)?;
        w.raw(self.as_bytes());
        Ok(())
    }
    fn read(r: &mut Reader<'_>) -> Result<Self, Refusal> {
        let n = r.count()?;
        String::from_utf8(r.take(n)?.to_vec()).map_err(|_| invalid("invalid UTF-8"))
    }
}

impl<T: Codec> Codec for Option<T> {
    fn write(&self, w: &mut Writer) -> Result<(), Refusal> {
        self.is_some().write(w)?;
        if let Some(x) = self {
            x.write(w)?;
        }
        Ok(())
    }
    fn read(r: &mut Reader<'_>) -> Result<Self, Refusal> {
        if bool::read(r)? {
            Ok(Some(T::read(r)?))
        } else {
            Ok(None)
        }
    }
    fn register(r: &mut Registry) {
        r.add::<T>();
    }
}
impl<T: Codec> Codec for Vec<T> {
    fn write(&self, w: &mut Writer) -> Result<(), Refusal> {
        self.len().write(w)?;
        for x in self {
            x.write(w)?;
        }
        Ok(())
    }
    fn read(r: &mut Reader<'_>) -> Result<Self, Refusal> {
        let n = r.count()?;
        (0..n).map(|_| T::read(r)).collect()
    }
    fn register(r: &mut Registry) {
        r.add::<T>();
    }
}
impl<T: Codec, const N: usize> Codec for [T; N] {
    fn write(&self, w: &mut Writer) -> Result<(), Refusal> {
        for x in self {
            x.write(w)?;
        }
        Ok(())
    }
    fn read(r: &mut Reader<'_>) -> Result<Self, Refusal> {
        (0..N)
            .map(|_| T::read(r))
            .collect::<Result<Vec<_>, _>>()?
            .try_into()
            .map_err(|_| invalid("invalid array width"))
    }
    fn register(r: &mut Registry) {
        r.add::<T>();
    }
}
impl<T: Codec> Codec for Arc<T> {
    fn write(&self, w: &mut Writer) -> Result<(), Refusal> {
        let value = self.clone();
        w.shared::<Self>(
            Arc::as_ptr(self) as usize,
            Box::new(move |w| value.as_ref().write(w)),
        )
    }
    fn read(r: &mut Reader<'_>) -> Result<Self, Refusal> {
        let id = u32::read(r)?;
        r.objects
            .get(&id)
            .and_then(|x| x.downcast_ref::<Self>())
            .cloned()
            .ok_or_else(|| invalid("unresolved or mistyped shared reference"))
    }
    fn register(r: &mut Registry) {
        r.factories
            .insert(type_name::<Self>(), |r| Ok(Box::new(Arc::new(T::read(r)?))));
        r.add::<T>();
    }
}
impl<T: Codec> Codec for Arc<[T]> {
    fn write(&self, w: &mut Writer) -> Result<(), Refusal> {
        let value = self.clone();
        w.shared::<Self>(
            Arc::as_ptr(self) as *const () as usize,
            Box::new(move |w| {
                value.len().write(w)?;
                for x in value.iter() {
                    x.write(w)?;
                }
                Ok(())
            }),
        )
    }
    fn read(r: &mut Reader<'_>) -> Result<Self, Refusal> {
        let id = u32::read(r)?;
        r.objects
            .get(&id)
            .and_then(|x| x.downcast_ref::<Self>())
            .cloned()
            .ok_or_else(|| invalid("unresolved slice reference"))
    }
    fn register(r: &mut Registry) {
        r.factories.insert(type_name::<Self>(), |r| {
            let value: Arc<[T]> = Vec::<T>::read(r)?.into();
            Ok(Box::new(value))
        });
        r.add::<T>();
    }
}
impl<T: Codec> Codec for Box<T> {
    fn write(&self, w: &mut Writer) -> Result<(), Refusal> {
        self.as_ref().write(w)
    }
    fn read(r: &mut Reader<'_>) -> Result<Self, Refusal> {
        Ok(Box::new(T::read(r)?))
    }
    fn register(r: &mut Registry) {
        r.add::<T>();
    }
}
impl<T: Codec> Codec for OnceLock<T> {
    fn write(&self, w: &mut Writer) -> Result<(), Refusal> {
        let value = self.get();
        value.is_some().write(w)?;
        if let Some(x) = value {
            x.write(w)?;
        }
        Ok(())
    }
    fn read(r: &mut Reader<'_>) -> Result<Self, Refusal> {
        Ok(Option::<T>::read(r)?.map_or_else(Self::new, Self::from))
    }
    fn register(r: &mut Registry) {
        r.add::<T>();
    }
}
impl<T: Codec> Codec for Mutex<T> {
    fn write(&self, w: &mut Writer) -> Result<(), Refusal> {
        self.lock()
            .map_err(|_| invalid("poisoned retained data"))?
            .write(w)
    }
    fn read(r: &mut Reader<'_>) -> Result<Self, Refusal> {
        Ok(Self::new(T::read(r)?))
    }
    fn register(r: &mut Registry) {
        r.add::<T>();
    }
}
impl<A: Codec, B: Codec> Codec for (A, B) {
    fn write(&self, w: &mut Writer) -> Result<(), Refusal> {
        self.0.write(w)?;
        self.1.write(w)
    }
    fn read(r: &mut Reader<'_>) -> Result<Self, Refusal> {
        Ok((A::read(r)?, B::read(r)?))
    }
    fn register(r: &mut Registry) {
        r.add::<A>();
        r.add::<B>();
    }
}
impl<A: Codec, B: Codec, C: Codec> Codec for (A, B, C) {
    fn write(&self, w: &mut Writer) -> Result<(), Refusal> {
        self.0.write(w)?;
        self.1.write(w)?;
        self.2.write(w)
    }
    fn read(r: &mut Reader<'_>) -> Result<Self, Refusal> {
        Ok((A::read(r)?, B::read(r)?, C::read(r)?))
    }
    fn register(r: &mut Registry) {
        r.add::<A>();
        r.add::<B>();
        r.add::<C>();
    }
}
impl<A: Codec, B: Codec, C: Codec, D: Codec> Codec for (A, B, C, D) {
    fn write(&self, w: &mut Writer) -> Result<(), Refusal> {
        self.0.write(w)?;
        self.1.write(w)?;
        self.2.write(w)?;
        self.3.write(w)
    }
    fn read(r: &mut Reader<'_>) -> Result<Self, Refusal> {
        Ok((A::read(r)?, B::read(r)?, C::read(r)?, D::read(r)?))
    }
    fn register(r: &mut Registry) {
        r.add::<A>();
        r.add::<B>();
        r.add::<C>();
        r.add::<D>();
    }
}

impl<K: Codec + Ord, V: Codec> Codec for BTreeMap<K, V> {
    fn write(&self, w: &mut Writer) -> Result<(), Refusal> {
        self.len().write(w)?;
        for (k, v) in self {
            k.write(w)?;
            v.write(w)?;
        }
        Ok(())
    }
    fn read(r: &mut Reader<'_>) -> Result<Self, Refusal> {
        let n = r.count()?;
        let mut map = Self::new();
        for _ in 0..n {
            if map.insert(K::read(r)?, V::read(r)?).is_some() {
                return Err(invalid("duplicate map key"));
            }
        }
        Ok(map)
    }
    fn register(r: &mut Registry) {
        r.add::<K>();
        r.add::<V>();
    }
}
impl<T: Codec + Ord> Codec for BTreeSet<T> {
    fn write(&self, w: &mut Writer) -> Result<(), Refusal> {
        self.len().write(w)?;
        for x in self {
            x.write(w)?;
        }
        Ok(())
    }
    fn read(r: &mut Reader<'_>) -> Result<Self, Refusal> {
        let n = r.count()?;
        let mut set = Self::new();
        for _ in 0..n {
            if !set.insert(T::read(r)?) {
                return Err(invalid("duplicate set entry"));
            }
        }
        Ok(set)
    }
    fn register(r: &mut Registry) {
        r.add::<T>();
    }
}
impl<T: Codec + Clone> Codec for im::Vector<T> {
    fn write(&self, w: &mut Writer) -> Result<(), Refusal> {
        self.len().write(w)?;
        for x in self {
            x.write(w)?;
        }
        Ok(())
    }
    fn read(r: &mut Reader<'_>) -> Result<Self, Refusal> {
        Ok(Vec::<T>::read(r)?.into_iter().collect())
    }
    fn register(r: &mut Registry) {
        r.add::<T>();
    }
}
impl<K: Codec + Clone + Ord, V: Codec + Clone> Codec for im::OrdMap<K, V> {
    fn write(&self, w: &mut Writer) -> Result<(), Refusal> {
        self.len().write(w)?;
        for (k, v) in self {
            k.write(w)?;
            v.write(w)?;
        }
        Ok(())
    }
    fn read(r: &mut Reader<'_>) -> Result<Self, Refusal> {
        Ok(BTreeMap::<K, V>::read(r)?.into_iter().collect())
    }
    fn register(r: &mut Registry) {
        r.add::<K>();
        r.add::<V>();
    }
}
impl<T: Codec + Clone + Ord> Codec for im::OrdSet<T> {
    fn write(&self, w: &mut Writer) -> Result<(), Refusal> {
        self.len().write(w)?;
        for x in self {
            x.write(w)?;
        }
        Ok(())
    }
    fn read(r: &mut Reader<'_>) -> Result<Self, Refusal> {
        Ok(BTreeSet::<T>::read(r)?.into_iter().collect())
    }
    fn register(r: &mut Registry) {
        r.add::<T>();
    }
}

// Existing public mark/fact codes are static strings. Bound the interner so
// repeated imports reuse them and arbitrary input cannot leak unbounded memory.
impl Codec for &'static str {
    fn write(&self, w: &mut Writer) -> Result<(), Refusal> {
        self.to_string().write(w)
    }
    fn read(r: &mut Reader<'_>) -> Result<Self, Refusal> {
        static STRINGS: OnceLock<Mutex<BTreeSet<&'static str>>> = OnceLock::new();
        let text = String::read(r)?;
        let mut strings = STRINGS
            .get_or_init(Default::default)
            .lock()
            .map_err(|_| invalid("code interner poisoned"))?;
        if let Some(text) = strings.get(text.as_str()) {
            return Ok(*text);
        }
        if text.len() + strings.iter().map(|s| s.len()).sum::<usize>() > 65_536 {
            return Err(invalid("retained code vocabulary exceeds limit"));
        }
        let text = Box::leak(text.into_boxed_str());
        strings.insert(text);
        Ok(text)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{Rational, Real};

    #[test]
    fn shared_exact_definitions_keep_their_graph_and_integer_bits() {
        let rational = Rational::parse("123456789012345678901234567890123456789012345678901/700000000000000000000000000000000000000000000000003").unwrap();
        let real = Real::from(rational.clone()).sqrt();
        let bytes = encode(&(real.clone(), real), &mut Decider::default()).unwrap();
        let (first, second): (Real, Real) = decode(&bytes, &mut Decider::default()).unwrap();
        assert!(Arc::ptr_eq(&first.node, &second.node));
        let crate::real::Expr::Sqrt(inner) = &first.node.expr else {
            panic!("expression definition changed")
        };
        assert_eq!(inner.rational_ref().unwrap(), &rational);
        assert_eq!(
            encode(&(first, second), &mut Decider::default()).unwrap(),
            bytes
        );
    }

    #[test]
    fn deep_definitions_use_numbered_records() {
        let mut real = Real::from_i64(2).sqrt();
        for _ in 0..2000 {
            real = &real + &Real::one();
        }
        let bytes = encode(&real, &mut Decider::default()).unwrap();
        let imported: Real = decode(&bytes, &mut Decider::default()).unwrap();
        assert_eq!(encode(&imported, &mut Decider::default()).unwrap(), bytes);
    }

    #[test]
    fn checksummed_bad_lengths_and_tags_refuse_without_panicking() {
        for value in [vec![true], vec![false, true]] {
            let mut bytes = encode(&value, &mut Decider::default()).unwrap();
            let end = bytes.len() - 32;
            // The last root byte is a boolean. Keep the envelope checksum valid.
            bytes[end - 1] = 7;
            let digest = Sha256::digest(&bytes[..end]);
            bytes[end..].copy_from_slice(&digest);
            assert_eq!(
                decode::<Vec<bool>>(&bytes, &mut Decider::default())
                    .unwrap_err()
                    .code,
                "handover.invalid-data"
            );
        }
        let mut bytes = encode(&Vec::<u8>::new(), &mut Decider::default()).unwrap();
        let end = bytes.len() - 32;
        bytes[end - 8..end].copy_from_slice(&u64::MAX.to_le_bytes());
        let digest = Sha256::digest(&bytes[..end]);
        bytes[end..].copy_from_slice(&digest);
        assert!(decode::<Vec<u8>>(&bytes, &mut Decider::default()).is_err());
    }

    #[test]
    fn concurrent_cache_fills_keep_the_tag_and_value_consistent() {
        for _ in 0..16 {
            let cells: Vec<_> = (0..4096)
                .map(|_| Arc::new(OnceLock::<u64>::new()))
                .collect();
            std::thread::scope(|scope| {
                scope.spawn(|| {
                    for (i, cell) in cells.iter().enumerate() {
                        cell.set(i as u64).unwrap();
                        std::thread::yield_now();
                    }
                });
                let bytes = encode(&cells, &mut Decider::default()).unwrap();
                let imported: Vec<Arc<OnceLock<u64>>> =
                    decode(&bytes, &mut Decider::default()).unwrap();
                for (i, cell) in imported.iter().enumerate() {
                    if let Some(value) = cell.get() {
                        assert_eq!(*value, i as u64);
                    }
                }
            });
        }
    }
}

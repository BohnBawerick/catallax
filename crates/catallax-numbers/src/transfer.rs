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
const VERSION: u32 = 2;

mod json;
pub use json::charge_json_number;

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
                writer.json(self)
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
struct RationalCache {
    values: HashMap<Vec<u8>, crate::Rational>,
    bytes: usize,
}

fn reserve<T>(decider: &mut Decider, values: &mut Vec<T>, extra: usize) -> Result<(), Refusal> {
    if extra > values.capacity().saturating_sub(values.len()) {
        decider.charge_work(values.len() as u64)?;
        values
            .try_reserve(extra)
            .map_err(|_| Refusal::new("handover.resource-limit", "transfer allocation failed"))?;
        decider.check_deadline()?;
    }
    Ok(())
}

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

pub struct Writer<'a> {
    decider: &'a mut Decider,
    bytes: Vec<u8>,
    keys: HashMap<(TypeId, usize), u32>,
    pending: Vec<Option<(&'static str, Emit)>>,
    owners: Vec<Box<dyn Any>>,
    dependencies: Vec<Vec<usize>>,
    active: Option<usize>,
}

impl<'a> Writer<'a> {
    pub fn new(decider: &'a mut Decider) -> Self {
        Self {
            decider,
            bytes: Vec::new(),
            keys: HashMap::new(),
            pending: Vec::new(),
            owners: Vec::new(),
            dependencies: Vec::new(),
            active: None,
        }
    }

    pub fn charge(&mut self, units: usize) -> Result<(), Refusal> {
        self.decider.charge_work(units as u64)
    }

    pub fn bytes_written(&self) -> usize {
        self.bytes.len()
    }

    fn shared<T: 'static>(
        &mut self,
        pointer: usize,
        owner: Box<dyn Any>,
        emit: Emit,
    ) -> Result<(), Refusal> {
        self.charge(1)?;
        let key = (TypeId::of::<T>(), pointer);
        let id = match self.keys.get(&key) {
            Some(id) => *id,
            None => {
                let id = u32::try_from(self.pending.len())
                    .map_err(|_| invalid("too many shared objects"))?;
                self.charge(4)?;
                if self.keys.len() == self.keys.capacity() {
                    self.charge(self.keys.len())?;
                    self.keys.try_reserve(1).map_err(|_| {
                        Refusal::new(
                            "handover.resource-limit",
                            "transfer index allocation failed",
                        )
                    })?;
                    self.decider.check_deadline()?;
                }
                reserve(self.decider, &mut self.owners, 1)?;
                reserve(self.decider, &mut self.pending, 1)?;
                reserve(self.decider, &mut self.dependencies, 1)?;
                self.keys.insert(key, id);
                // A mutable cache may replace its Arc while later records are written.
                self.owners.push(owner);
                self.pending.push(Some((type_name::<T>(), emit)));
                self.dependencies.push(Vec::new());
                id
            }
        };
        if let Some(active) = self.active {
            self.charge(1)?;
            reserve(self.decider, &mut self.dependencies[active], 1)?;
            self.dependencies[active].push(id as usize);
        }
        id.write(self)
    }

    pub fn raw(&mut self, bytes: &[u8]) -> Result<(), Refusal> {
        self.charge(bytes.len())?;
        reserve(self.decider, &mut self.bytes, bytes.len())?;
        self.decider.check_deadline()?;
        self.bytes.extend_from_slice(bytes);
        Ok(())
    }
}

pub struct Reader<'a> {
    version: u32,
    data: &'a [u8],
    at: usize,
    objects: &'a [Option<Box<dyn Any>>],
    decider: &'a mut Decider,
    handles: &'a HandleRanges,
    rationals: &'a mut RationalCache,
}

impl<'a> Reader<'a> {
    pub fn version(&self) -> u32 {
        self.version
    }

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
    #[inline]
    pub fn take(&mut self, count: usize) -> Result<&'a [u8], Refusal> {
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

    #[inline]
    pub fn bytes(&mut self) -> Result<&'a [u8], Refusal> {
        let count = self.count()?;
        self.take(count)
    }

    #[inline]
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
    let (root, records, order) = {
        let mut writer = Writer::new(decider);
        body.write(&mut writer)?;
        let root = std::mem::take(&mut writer.bytes);
        let mut records = Vec::new();
        let mut i = 0;
        while i < writer.pending.len() {
            writer.charge(1)?;
            writer.active = Some(i);
            let (tag, emit) = writer.pending[i]
                .take()
                .ok_or_else(|| invalid("missing shared object"))?;
            emit(&mut writer)?;
            writer.charge(1)?;
            reserve(writer.decider, &mut records, 1)?;
            records.push((tag, std::mem::take(&mut writer.bytes)));
            i += 1;
        }
        // Iterative postorder preserves graph sharing without recursive traversal.
        writer.charge(records.len())?;
        let mut order = Vec::new();
        let mut state = vec![0_u8; records.len()];
        writer.charge(records.len())?;
        for start in 0..records.len() {
            if state[start] == 2 {
                continue;
            }
            writer.charge(1)?;
            let mut stack = vec![(start, false)];
            while let Some((id, exit)) = stack.pop() {
                writer.charge(1)?;
                if exit {
                    state[id] = 2;
                    reserve(writer.decider, &mut order, 1)?;
                    order.push(id);
                } else if state[id] != 2 {
                    if state[id] == 1 {
                        return Err(invalid("cyclic shared definitions"));
                    }
                    state[id] = 1;
                    reserve(writer.decider, &mut stack, 1)?;
                    stack.push((id, true));
                    writer.charge(writer.dependencies[id].len())?;
                    for &child in writer.dependencies[id].iter().rev() {
                        if state[child] != 2 {
                            reserve(writer.decider, &mut stack, 1)?;
                            stack.push((child, false));
                        }
                    }
                }
            }
        }
        (root, records, order)
    };
    let mut output = Writer::new(decider);
    output.raw(MAGIC)?;
    VERSION.write(&mut output)?;
    type_name::<T>().write(&mut output)?;
    output.charge(records.len())?;
    let mut tags = BTreeSet::new();
    for (tag, _) in &records {
        output.charge(tag.len())?;
        tags.insert(*tag);
    }
    output.charge(tags.len())?;
    let tags: Vec<_> = tags.into_iter().collect();
    tags.len().write(&mut output)?;
    output.charge(tags.len())?;
    for tag in &tags {
        (*tag).write(&mut output)?;
    }
    records.len().write(&mut output)?;
    output.charge(order.len())?;
    for id in order {
        (id as u32).write(&mut output)?;
        output.charge(tags.len().max(1).ilog2() as usize + 1)?;
        tags.binary_search(&records[id].0)
            .map_err(|_| invalid("missing record type"))?
            .write(&mut output)?;
        records[id].1.len().write(&mut output)?;
        output.raw(&records[id].1)?;
    }
    root.len().write(&mut output)?;
    output.raw(&root)?;
    output.charge(output.bytes.len())?;
    let digest = Sha256::digest(&output.bytes);
    output.raw(&digest)?;
    Ok(output.bytes)
}

pub fn decode<T: Codec>(bytes: &[u8], decider: &mut Decider) -> Result<T, Refusal> {
    decider.charge_work(bytes.len() as u64)?;
    if bytes.get(..8) != Some(MAGIC.as_slice()) {
        return Err(invalid("body handover magic is missing"));
    }
    let version = u32::from_le_bytes(
        bytes
            .get(8..12)
            .ok_or_else(|| invalid("truncated version"))?
            .try_into()
            .map_err(|_| invalid("truncated version"))?,
    );
    if version != 1 && version != VERSION {
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
    let mut objects = Vec::new();
    let handles = RefCell::new(HashMap::new());
    let mut rationals = RationalCache::default();
    let mut reader = Reader {
        version,
        data: &bytes[12..end],
        at: 0,
        objects: &objects,
        decider,
        handles: &handles,
        rationals: &mut rationals,
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
    // Each record has a u32 ID and two u64 fields before its payload.
    if count > reader.data.len().saturating_sub(reader.at) / 20 {
        return Err(invalid("shared object count exceeds input"));
    }
    let mut at = reader.at;
    decider.charge_work(count as u64)?;
    objects
        .try_reserve_exact(count)
        .map_err(|_| Refusal::new("handover.resource-limit", "shared object allocation failed"))?;
    decider.check_deadline()?;
    objects.resize_with(count, || None);
    for _ in 0..count {
        let mut reader = Reader {
            version,
            data: &bytes[12..end],
            at,
            objects: &objects,
            decider,
            handles: &handles,
            rationals: &mut rationals,
        };
        let id = u32::read(&mut reader)?;
        if id as usize >= count || objects[id as usize].is_some() {
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
            version,
            data: &bytes[12 + start..12 + start + length],
            at: 0,
            objects: &objects,
            decider,
            handles: &handles,
            rationals: &mut rationals,
        };
        let value = factory(&mut record)?;
        record.finish()?;
        objects[id as usize] = Some(value);
    }
    let mut reader = Reader {
        version,
        data: &bytes[12..end],
        at,
        objects: &objects,
        decider,
        handles: &handles,
        rationals: &mut rationals,
    };
    let length = reader.count()?;
    let start = reader.at;
    reader.take(length)?;
    reader.finish()?;
    let mut root = Reader {
        version,
        data: &bytes[12 + start..12 + start + length],
        at: 0,
        objects: &objects,
        decider,
        handles: &handles,
        rationals: &mut rationals,
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
            fn write(&self, w: &mut Writer) -> Result<(), Refusal> { w.raw(&self.to_le_bytes())?; Ok(()) }
            #[inline]
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
    #[inline]
    fn read(r: &mut Reader<'_>) -> Result<Self, Refusal> {
        Self::try_from(u64::read(r)?).map_err(|_| invalid("index exceeds platform width"))
    }
}
impl Codec for num_bigint::BigInt {
    fn write(&self, w: &mut Writer) -> Result<(), Refusal> {
        let length = self.bits().saturating_div(8).saturating_add(1);
        w.decider.charge_work(length)?;
        let bytes = self.to_signed_bytes_le();
        bytes.len().write(w)?;
        w.raw(&bytes)
    }
    fn read(r: &mut Reader<'_>) -> Result<Self, Refusal> {
        Ok(Self::from_signed_bytes_le(r.bytes()?))
    }
}
impl Codec for crate::Rational {
    fn write(&self, w: &mut Writer) -> Result<(), Refusal> {
        self.n.write(w)?;
        self.d.write(w)
    }
    fn read(r: &mut Reader<'_>) -> Result<Self, Refusal> {
        let start = r.at;
        let numerator = r.bytes()?;
        let denominator = r.bytes()?;
        let key = &r.data[start..r.at];
        r.decider.charge_work(key.len() as u64)?;
        // Only identical encodings reuse a completed canonicality check.
        if let Some(value) = r.rationals.values.get(key) {
            return Ok(value.clone());
        }
        let n = crate::Integer(num_bigint::BigInt::from_signed_bytes_le(numerator));
        let d = crate::Integer(num_bigint::BigInt::from_signed_bytes_le(denominator));
        if d.sign() != crate::Sign::Positive {
            return Err(invalid("non-positive rational denominator"));
        }
        let value =
            Self::new(n.clone(), d.clone()).ok_or_else(|| invalid("zero rational denominator"))?;
        if value.n != n || value.d != d {
            return Err(invalid("noncanonical rational"));
        }
        let cache = &mut *r.rationals;
        if cache.values.len() == cache.values.capacity() {
            r.decider.charge_work(cache.bytes as u64)?;
            cache.values.try_reserve(1).map_err(|_| {
                Refusal::new(
                    "handover.resource-limit",
                    "rational cache allocation failed",
                )
            })?;
            r.decider.check_deadline()?;
        }
        r.decider.charge_work(key.len() as u64)?;
        cache.bytes += key.len();
        cache.values.insert(key.to_vec(), value.clone());
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
        w.raw(self.as_bytes())?;
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
        w.charge(self.len())?;
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
        w.charge(self.len())?;
        for x in self {
            x.write(w)?;
        }
        Ok(())
    }
    fn read(r: &mut Reader<'_>) -> Result<Self, Refusal> {
        let mut values = std::array::from_fn(|_| None);
        for value in &mut values {
            *value = Some(T::read(r)?);
        }
        Ok(values.map(Option::unwrap))
    }
    fn register(r: &mut Registry) {
        r.add::<T>();
    }
}
impl<T: Codec> Codec for Arc<T> {
    fn write(&self, w: &mut Writer) -> Result<(), Refusal> {
        w.charge(3)?;
        let value = self.clone();
        w.shared::<Self>(
            Arc::as_ptr(self) as usize,
            Box::new(self.clone()),
            Box::new(move |w| value.as_ref().write(w)),
        )
    }
    fn read(r: &mut Reader<'_>) -> Result<Self, Refusal> {
        let id = u32::read(r)?;
        r.objects
            .get(id as usize)
            .and_then(Option::as_ref)
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
        w.charge(3)?;
        let value = self.clone();
        w.shared::<Self>(
            Arc::as_ptr(self) as *const () as usize,
            Box::new(self.clone()),
            Box::new(move |w| {
                w.charge(value.len())?;
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
            .get(id as usize)
            .and_then(Option::as_ref)
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
        w.charge(self.len())?;
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
        w.charge(self.len())?;
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
        w.charge(self.len())?;
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
        w.charge(self.len())?;
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
        w.charge(self.len())?;
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
        self.len().write(w)?;
        w.raw(self.as_bytes())
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
    fn fixed_arrays_decode_empty_values_and_refuse_invalid_elements() {
        let bytes = encode(&[0_u64; 0], &mut Decider::default()).unwrap();
        assert_eq!(
            decode::<[u64; 0]>(&bytes, &mut Decider::default()).unwrap(),
            [0_u64; 0]
        );
        let mut bytes = encode(&[true, false, true], &mut Decider::default()).unwrap();
        let end = bytes.len() - 32;
        bytes[end - 1] = 7;
        let digest = Sha256::digest(&bytes[..end]);
        bytes[end..].copy_from_slice(&digest);
        assert_eq!(
            decode::<[bool; 3]>(&bytes, &mut Decider::default())
                .unwrap_err()
                .code,
            "handover.invalid-data"
        );
    }

    #[test]
    fn numbered_objects_refuse_duplicate_and_out_of_range_ids() {
        fn length(bytes: &[u8], at: &mut usize) -> usize {
            let value = u64::from_le_bytes(bytes[*at..*at + 8].try_into().unwrap()) as usize;
            *at += 8;
            value
        }
        let bytes = encode(
            &vec![Arc::new(1_u64), Arc::new(2_u64)],
            &mut Decider::default(),
        )
        .unwrap();
        let mut at = 12;
        at += 8 + u64::from_le_bytes(bytes[at..at + 8].try_into().unwrap()) as usize;
        for _ in 0..length(&bytes, &mut at) {
            let n = length(&bytes, &mut at);
            at += n;
        }
        assert_eq!(length(&bytes, &mut at), 2);
        let first = at;
        at += 12;
        let n = length(&bytes, &mut at);
        at += n;
        let second = at;
        for (offset, id) in [(first, 2_u32), (second, 0_u32)] {
            let mut corrupt = bytes.clone();
            corrupt[offset..offset + 4].copy_from_slice(&id.to_le_bytes());
            let end = corrupt.len() - 32;
            let digest = Sha256::digest(&corrupt[..end]);
            corrupt[end..].copy_from_slice(&digest);
            assert_eq!(
                decode::<Vec<Arc<u64>>>(&corrupt, &mut Decider::default())
                    .unwrap_err()
                    .code,
                "handover.invalid-data"
            );
        }
    }

    #[test]
    fn repeated_rationals_keep_their_bytes_and_pay_for_every_read() {
        let rational = Rational::parse("123456789012345678901234567890123456789012345678901/700000000000000000000000000000000000000000000000003").unwrap();
        let source = vec![rational; 64];
        let bytes = encode(&source, &mut Decider::default()).unwrap();
        let (restored, work) = crate::with_operation_budget(1_000_000, || {
            decode::<Vec<Rational>>(&bytes, &mut Decider::default())
        });
        let restored = restored.unwrap();
        assert_eq!(restored, source);
        assert_eq!(encode(&restored, &mut Decider::default()).unwrap(), bytes);
        let (refused, spent) = crate::with_operation_budget(work.used - 1, || {
            decode::<Vec<Rational>>(&bytes, &mut Decider::default())
        });
        assert!(refused.is_err());
        assert!(spent.exhausted);
    }

    #[test]
    fn cached_rationals_do_not_admit_noncanonical_or_invalid_encodings() {
        let half = Rational::from_ratio(1, 2).unwrap();
        for (n, d) in [(2, 4), (0, 2), (1, 0), (-1, -2)] {
            let invalid = Rational {
                n: n.into(),
                d: d.into(),
            };
            let bytes = encode(
                &vec![half.clone(), half.clone(), invalid],
                &mut Decider::default(),
            )
            .unwrap();
            assert_eq!(
                decode::<Vec<Rational>>(&bytes, &mut Decider::default())
                    .unwrap_err()
                    .code,
                "handover.invalid-data"
            );
        }
    }

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

    #[test]
    fn registered_objects_stay_alive_until_the_snapshot_finishes() {
        let value = Arc::new(42_u64);
        let weak = Arc::downgrade(&value);
        let mut decider = Decider::default();
        let mut writer = Writer::new(&mut decider);
        value.write(&mut writer).unwrap();
        drop(value);
        let (_, emit) = writer.pending[0].take().unwrap();
        emit(&mut writer).unwrap();
        assert_eq!(weak.upgrade().as_deref(), Some(&42));
        drop(writer);
        assert!(weak.upgrade().is_none());
    }

    #[test]
    fn collection_walks_stop_before_visiting_a_large_payload() {
        use std::sync::atomic::{AtomicUsize, Ordering};
        static VISITED: AtomicUsize = AtomicUsize::new(0);
        struct Probe;
        impl Codec for Probe {
            fn write(&self, w: &mut Writer) -> Result<(), Refusal> {
                VISITED.fetch_add(1, Ordering::Relaxed);
                1_u8.write(w)
            }
            fn read(_: &mut Reader<'_>) -> Result<Self, Refusal> {
                unreachable!()
            }
        }
        let values: Vec<_> = (0..100_000).map(|_| Probe).collect();
        let (result, work) =
            crate::with_operation_budget(16, || encode(&values, &mut Decider::default()));
        assert!(result.is_err());
        assert!(work.exhausted);
        assert_eq!(VISITED.load(Ordering::Relaxed), 0);
    }

    #[test]
    fn json_walks_and_string_scans_refuse_before_materialization() {
        use std::sync::atomic::{AtomicUsize, Ordering};
        static VISITED: AtomicUsize = AtomicUsize::new(0);
        struct Probe;
        impl serde::Serialize for Probe {
            fn serialize<S: serde::Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
                VISITED.fetch_add(1, Ordering::Relaxed);
                serializer.serialize_u8(1)
            }
        }
        let values: Vec<_> = (0..100_000).map(|_| Probe).collect();
        let (result, work) = crate::with_operation_budget(16, || {
            let mut decider = Decider::default();
            let mut writer = Writer::new(&mut decider);
            writer.json(&values)
        });
        assert!(result.is_err());
        assert!(work.exhausted);
        assert_eq!(VISITED.load(Ordering::Relaxed), 0);
        let text = "unescaped-text".repeat(100_000);
        let ((result, written), work) = crate::with_operation_budget(16, || {
            let mut decider = Decider::default();
            let mut writer = Writer::new(&mut decider);
            (writer.json(&text), writer.bytes_written())
        });
        assert!(result.is_err());
        assert!(work.exhausted);
        assert_eq!(written, 0);
    }

    #[test]
    fn integer_conversion_is_charged_before_allocating_encoded_digits() {
        let value: num_bigint::BigInt = num_bigint::BigInt::from(1_u8) << 1_000_000_u32;
        let ((result, written), work) = crate::with_operation_budget(16, || {
            let mut decider = Decider::default();
            let mut writer = Writer::new(&mut decider);
            (value.write(&mut writer), writer.bytes_written())
        });
        assert!(result.is_err());
        assert!(work.exhausted);
        assert_eq!(written, 0);
    }

    #[test]
    fn metered_json_preserves_payload_and_reads_version_one() {
        struct Metadata(Vec<String>);
        impl Codec for Metadata {
            fn write(&self, w: &mut Writer) -> Result<(), Refusal> {
                w.json(&self.0)
            }
            fn read(r: &mut Reader<'_>) -> Result<Self, Refusal> {
                Ok(Self(serde_json::from_slice(&Vec::<u8>::read(r)?).unwrap()))
            }
        }
        let value = Metadata(vec!["line\nquote\"".into(), "é".into()]);
        let json = serde_json::to_vec(&value.0).unwrap();
        let kind = type_name::<Metadata>();
        let mut legacy = MAGIC.to_vec();
        legacy.extend(1_u32.to_le_bytes());
        legacy.extend((kind.len() as u64).to_le_bytes());
        legacy.extend(kind.as_bytes());
        legacy.extend(0_u64.to_le_bytes());
        legacy.extend(0_u64.to_le_bytes());
        legacy.extend((8 + json.len() as u64).to_le_bytes());
        legacy.extend((json.len() as u64).to_le_bytes());
        legacy.extend(json);
        legacy.extend(Sha256::digest(&legacy).as_slice());
        let mut current = legacy[..legacy.len() - 32].to_vec();
        current[8..12].copy_from_slice(&2_u32.to_le_bytes());
        current.extend(Sha256::digest(&current).as_slice());
        assert_eq!(encode(&value, &mut Decider::default()).unwrap(), current);
        let restored: Metadata = decode(&legacy, &mut Decider::default()).unwrap();
        assert_eq!(restored.0, value.0);
    }
}

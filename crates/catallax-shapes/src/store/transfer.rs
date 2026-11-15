use super::*;
use catallax_numbers::transfer::{Codec, Reader, Registry, Writer};

impl<T: 'static> Codec for Handle<T> {
    fn write(&self, writer: &mut Writer) -> Result<(), Refusal> {
        self.index.write(writer)
    }
    fn read(reader: &mut Reader<'_>) -> Result<Self, Refusal> {
        let index = u32::read(reader)?;
        reader.handle::<T>(index);
        Ok(Self {
            index,
            _kind: PhantomData,
        })
    }
}

impl<T: Codec + Clone> Codec for Arena<T> {
    fn write(&self, writer: &mut Writer) -> Result<(), Refusal> {
        self.values.write(writer)
    }
    fn read(reader: &mut Reader<'_>) -> Result<Self, Refusal> {
        let mut arena = Self::new();
        for value in Vector::<T>::read(reader)? {
            arena.push(value);
        }
        reader.arena::<T>(arena.len())?;
        Ok(arena)
    }
    fn register(registry: &mut Registry) {
        registry.add::<T>();
    }
}

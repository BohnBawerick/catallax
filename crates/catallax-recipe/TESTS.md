# Recipe format tests

The tests use `docs/RECIPE.md` and the crate's public interface. They need no kernel and evaluate no geometry.

```sh
cargo test -p catallax-recipe --release
```

`tests/fixtures/every-operation.json` is a recipe with one step of every operation and a name of every role, in canonical form. The tests check that:

- parsing it and writing it back gives the same bytes;
- a copy with keys in reverse order, every number respelled, no whitespace and every edge and vertex name listed backwards also writes back to the same bytes;
- writing is repeatable, and every number has one canonical spelling;
- the reader rejects JSON numbers in geometric places, unknown keys, operations and roles, bad ids, malformed names and other formats;
- the structural check rejects each structural error with its stable code.

The fixture's geometry is meant to read sensibly, but nothing has built it. A change to the canonical form changes the fixture's bytes; that is a format change and needs a new format version.

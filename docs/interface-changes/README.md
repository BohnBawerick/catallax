# Interface-change fragments

Record every change to a public interface here: a type, field, function signature, trait, outcome or refusal code, budget default, or file format that crosses a piece boundary or the judge boundary.

## Add a change

- Public interface changes require review and a migration for affected callers.
- Add one new Markdown fragment per change in the same pull request that makes it, before or with the code. Use `YYYY-MM-DD-short-slug.md`, for example `2026-10-26-curve-query-result.md`. Choose a descriptive slug unique to the change; do not allocate a sequence number.
- Do not edit `docs/INTERFACE-CHANGES.md`, this README or a shared list to register a fragment. Independent pull requests add independent files.
- Entries are never edited or removed. A correction is a new fragment that names and links to the one it corrects.
- Name the affected pieces in the pull request so their callers can migrate together.
- Relative links in fragments resolve from this directory. Use `../ARCHITECTURE.md` to link to the architecture document.

Entry format:

```markdown
## YYYY-MM-DD: short title

Pieces: 1, 2, ... (or "judge interface", "all")
Change: what changed, in one or two sentences.
Why: the reason.
Migration: what a worker on an affected piece must do.
```

## Read the history

The `0001-` through `0173-` filenames preserve the original log order at migration, not date order. Dates and any pull-request references remain in the entries. Read later fragments by their dates; same-day independent changes have no implied ordering. A dependent change should link to its prerequisite fragment.

Release notes are generated from merged pull requests by GitHub, as described in [Releasing](../RELEASING.md). They do not assemble this interface record. No release tooling or generated release output changes with this layout.

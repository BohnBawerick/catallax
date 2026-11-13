# P4 enclosure reuse and remaining budget gap

P4's complete verification request falls from 1,540,488 to **1,366,808 work
units**, saving 173,680. The 1M scoreboard gate remains unmet by 366,808 units.
It still refuses during the rounding Boolean with `numbers.operation-budget`.
The measured 1M stop is 999,993 units. M4 remains open.

The later [remaining-work verdict](parked/m4-p4-scoreboard.md) parks the
unfinished 1M work, inventories the residual costs and states what resuming
must prove. The enclosure-reuse measurements below remain their original record.

The baseline is main `4d51b5c36bdad29632a448585290cfaa627098ba`. It already
contains the [cylindrical probe proofs](https://github.com/BohnBawerick/catallax/pull/481).
Fresh runs reproduce that landing's 1,540,488-unit total and 999,995-unit 1M
stop. This change retains those proofs.

## Reused proof

Shell recipes now enable the existing request-scoped enclosure cache, as fillet
recipes already do. Completed enclosures of the same immutable expression node
at the same precision can be reused through the remainder of the request.
The existing cache pins definitions, bounds its size, clears between requests,
and preserves cancellation and bounded numeric-decision checks. No new cache or
numeric method is added.

Reuse saves 157,581 units in the rounding Boolean and 16,099 in STEP export.
The cache does not replace the bounded coefficient, sign or root proofs. Their
remaining costs are recorded below.

## Complete request profile

The input is the unchanged [canonical P4 recipe](../parts/P4-shelled-tray.json).
[Portable evidence](evidence/m4-p4-enclosure-reuse.json) records the source and
binary hashes, every charged caller, operation counts and production comparisons.
Both development observer builds use release optimisation and the existing
`work-profile` feature. That feature disables helper execution. Production
thread comparisons are separate from these observations.

| Operation | Main | With enclosure reuse |
| --- | ---: | ---: |
| Shell | 62,617 | 62,617 |
| Rounding Boolean | 1,299,354 | 1,141,773 |
| Exact construction | 1,361,971 | 1,204,390 |
| Summary | 32 | 32 |
| STEP export | 178,485 | 162,386 |
| Complete verification | 1,540,488 | 1,366,808 |

The construction row includes the shell and rounding Boolean. The observer
attributes every work unit, with zero unattributed work.

| Proof cost | Main | With enclosure reuse |
| --- | ---: | ---: |
| Probe coefficient enclosures | 205,620 | 205,620 |
| Circular-chart constant enclosures | 180,471 | 154,203 |
| Root isolation | 110,435 | 110,435 |
| Boolean geometric signs | 273,311 | 273,311 |
| Stored point/curve incidence signs | 143,448 | 143,448 |
| Stored point/curve incidence prefilter | 50,692 | 50,692 |
| All remaining calls | 576,511 | 429,099 |
| Total | 1,540,488 | 1,366,808 |

The three previously identified callers account for 496,526 units on main.
The other calls account for **1,043,962**, already exceeding the complete
scoreboard allowance by 43,962 units. Their work cannot be ignored when choosing
the next reduction. Of this change's saving, 26,268 units come from the known
chart caller and **147,412 from other calls**.

After reuse the three known callers still cost 470,258 units. All other calls
cost 896,550. The 366,808-unit gap needs further certified reductions. Keeping
the current summary and export costs would require exact construction to fit
837,582 units, including a rounding Boolean of at most 774,965 if shell work
stays fixed. Current construction alone exceeds 1M by 204,390 units.

## Production checks

Production builds have no profiling features. The baseline is main at the
revision above; the candidate code is `a565f69774044fbe4a9f77b4d4db748540764e4b`.

| Request | Main work | Candidate work | Result |
| --- | ---: | ---: | --- |
| P4, 20M or 100M cap | 1,540,488 | 1,366,808 | Answered, identical manifest and STEP bytes |
| P4 thickness 1.5, 100M cap | 1,617,846 | 1,444,166 | Answered, identical manifest and STEP bytes |
| P4, 1M cap | 999,995 | 999,993 | `numbers.operation-budget`, no STEP |
| P4 thickness 1.5, 1M cap | 999,993 | 999,995 | `numbers.operation-budget`, no STEP |
| P4, 100k cap | 99,984 | 99,984 | `numbers.operation-budget`, no STEP |
| Near-miss control, 1M cap | 6,980 | 6,980 | Answered, identical manifest and STEP bytes |

The answered comparisons exclude build identity and numeric statistics. Each
candidate configuration also compares complete responses with only clocks
removed, including work, numeric statistics and operation logs. All 48 fresh
candidate requests agree at 1, 2 and 4 kernel threads, with 1, 2 and 4 persistent
workers, in normal and Boolean run-ahead modes. The baseline adds eight fresh
requests. Every worker exits zero; no cached artifact is used.

The original P4 STEP SHA-256 remains
`f9a8b9a15a576a4989690f75e1bd9f697e052299b80b33c1bfb4fdc4634c7381`.
The new CLI regression check retains that hash and caps completed-request work
at the measured 1,366,808 units under the existing 100M verification allowance.
This check guards the saving while the separate 1M gate remains open.

## Verdict and reproduction

The 1M build line is unfinished. The retained change reduces repeated arithmetic
while preserving the proof methods. No further certified reduction was
established in this pass. No tolerance, target recipe or work allowance changes.
The 100M request measures completed work under the existing verification cap;
it does not qualify the 1M gate.

Build the observer, set an external `PROOF` directory and reproduce both caps:

```sh
cargo build --release -p catallax-cli --features work-profile \
  --example work_profile -j 2
CATALLAX_THREADS=1 target/release/examples/work_profile \
  parts/P4-shelled-tray.json --out "$PROOF/profile-100m" --max-work 100000000
CATALLAX_THREADS=1 target/release/examples/work_profile \
  parts/P4-shelled-tray.json --out "$PROOF/profile-1m" --max-work 1000000
```

For production validation, build `catallax-verify` without profiling features
and use the [verification driver](VERIFY-INTERFACE.md). Compare complete P4 and
its thickness-1.5 edit at the existing 20M and 100M caps. Separately inspect the
unchanged 1M gate. Fresh runs at 1, 2 and 4 kernel threads and persistent workers
must agree on responses with clocks removed, work, refusals and strict STEP bytes.
SolidDraft interaction, field qualification and the complete M4 exit remain
separate gates. Raw evidence retained privately.

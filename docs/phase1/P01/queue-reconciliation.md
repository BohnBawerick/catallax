# Queue reconciliation: records are not live jobs

Read-only Firstmate mailbox scan, seq 1–1285, followed by an empty terminal page. No replies, acknowledgements, subscriptions or queue mutations. Metadata only is retained in `evidence/phase1/P01/queue-inventory.json`; unrelated or sensitive message bodies are not republished. Repository queue sources read: archived plan §§11.3, 11.4, 13.2 and ROADMAP phase/order/closed-work sections.

## Coordinator-defined historical census

The coordinator's verified earlier census defines **29 historical pending directive/disposition records** as the 32 non-fleet pending records minus three nonrouted records: seq **3, 5 and 410**. The five correction records **897, 899, 901, 908 and 910 remain within the 29**. This corrects the earlier missing-selection/count discrepancy; no mailbox rescan was performed for this correction.

Exact 29 sequences: **62, 64, 65, 146, 248, 296, 513, 743, 820, 821, 822, 823, 828, 831, 897, 899, 901, 908, 910, 968, 969, 984, 988, 1000, 1012, 1023, 1024, 1149, 1155**.

This is a historical pending-record classification, not 29 verified live jobs or 29 definitively stale jobs. Mailbox pending status establishes neither execution state nor technical acceptance. The metadata manifest records the coordinator-defined selection and the three excluded nonrouted records. The table below preserves distinctions for all 32 inspected records; seq 3, 5 and 410 are context outside the 29-record census.

## Classification of available non-fleet pending records

Classification is a dispatch recommendation, not a mailbox-status update, new execution authorisation, or a current technical acceptance verdict. Cross-project records are excluded from this Catallax job; no companion repository was accessed.

| Sequence | Classification | Effect on proposed Catallax dispatch |
|---|---|---|
| 3, 5 | Historical bridge setup/coordination, not product work | Exclude; subsequent mailbox history exists; do not retry setup. |
| 62 | Historical multi-clause external-lane rulings | Exclude from P01; retain dispositions with original owner; no retry/handoff. |
| 64, 65 | Duplicate references to prior rulings | Coalesce with referenced ruling; no independent work item. |
| 146 | External-lane policy disposition | Preserve as history; no Catallax coverage exception inferred. |
| 248 | Operational reconciliation and closures | Preserve explicit “completed is not validated”; no new task or acceptance inferred. |
| 296 | Historical Catallax file-ownership/conflict ruling | Preserve ownership history; do not resurrect batch 7/PR48. Live-owner confirmation needed if relevant again. |
| 410 | Existing Catallax PR146 checkpoint, explicitly no new assignment | Exclude duplicate dispatch; partial test/source report is not final qualification. |
| 513 | Delivery-closed disposition | Exclude duplicate delivery/release-evidence job; no rerun requested. |
| 743 | Repository-location notice | Use current referenced companion URLs in future reports; no remote edit or companion access. |
| 820, 821, 822 | Three distinct M/N/O design evidence attachments for one companion PR review | Group under original review, not three build jobs; distinct scenarios remain distinct evidence. Companion access is outside P01. |
| 823 | Wrapper requesting completion of the above review | Coalesce with 820–822 and original review; not another review/build job. |
| 828 | Historical exact-head Catallax PR241 premerge review | Frozen main log includes merged work; later mailbox seq833 records stand-in review/merge. Remove duplicate old premerge dispatch; no new validation claimed. |
| 831 | Historical exact-head Catallax PR242 premerge review | Same later seq833 disposition; retain query-only/export limitation. Remove duplicate old premerge dispatch, not qualification requirements. |
| 897 | Report count correction | Preserve correction; no work dispatch. |
| 899 | Report timestamp correction | Preserve correction; no work dispatch. |
| 901 | Timing-claim correction | Preserve retraction of unsound one-second verdict; no rerun or speed qualification inferred. |
| 908 | Report/test-method correction outside current work | No load/probe task; no historical instructions revived. |
| 910 | Report timestamp correction | Preserve correction; no work dispatch. |
| 968, 969 | Duplicate external-port resumption directives | Coalesce historical directive, exclude from Catallax P01; do not resume external port. |
| 984 | Findings routed to another owner, explicitly no further action | Retain owner disposition; exclude duplicate fix task. |
| 988 | External merged/accepted disposition | Exclude duplicate action; no fresh validation inferred. |
| 1000 | Sensitive operational incident disposition outside project | Exclude; no details republished, no security/account action. |
| 1012 | Other-project delivery confirmation, explicitly no duplicate nudge | Exclude duplicate communication. |
| 1023, 1024 | Other-project deployment request followed by completion report | Coalesce as request/result history; exclude from Catallax, no deployment action. |
| 1149 | Catallax publication notice plus separate PR364 review request | Preserve publication/review distinction. Historical head is not current; no new review from P01. |
| 1155 | Exact-head correction for PR366 | Retain corrected SHA in historical review context, not a new job or gate retry. |

## Repository queue deduplication

Archived plan §11.3's five queue entries are historical acceptance context. Display of newly built kinds maps to the app-alongside lane; primitive/name equivalence maps to edit/name acceptance; comparison gaps map to the standing resolve lane; error-policy enforcement maps to the current owner/API work, with N01 ownership preserved; cloud offload is capacity rather than product scope. None warrants a duplicate P01 algorithm job.

Archived §11.4 distinguishes Done, Partly done, parked and superseded. Keep those distinctions. ROADMAP explicitly closes or parks selected research, moves timing to phase 5, and parks FEA pending the owner. These are source-level dispositions, not evidence that every historically pending directive is stale. Missing whole-part, export, app or final-head evidence stays missing.

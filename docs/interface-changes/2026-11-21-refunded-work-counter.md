## 2026-11-21: refunded work has its own counter

Pieces: 1 (numbers), judge interface (`catallax-verify` responses).
Change: `catallax_numbers::OperationWork` gains `refunded: u64`. `refund_on_refusal` still takes a refused attempt's charges back out of `used`, and now adds them to `refunded`. Nested refused attempts are each counted once. An `OperationTask` item carries its refunded work to the request when its record is committed, refused at admission or discarded; a discarded item's other charges stay out of both counters. A `catallax-response/1` `work` object gains `refunded`.
Why: Checkpoint 2 decision Q13, predictor step 1 ([section 4.5](../checkpoints/checkpoint-2-decisions.md#45-item-30b-one-model-and-the-predictor)). A refused exact attempt is machine time. Until now its count was erased, so a doomed attempt cost nothing on the record. The predictor's cap tiers are judged against this counter.
Migration: code that builds an `OperationWork` value adds `refunded: 0`. Nothing else changes: `used`, `exhausted`, every refusal and every allowance stay as they were, and a request with no refused attempt reports `refunded` 0.

---
name: blind-lot-subagents
description: Run and compare benchmark-grade order-only multiple-myeloma LOT experiments with Codex subagents instead of direct provider API calls. Use for development model/prompt/retrieval searches, sanitized worker-bundle generation, repeated blind runs, locked final evaluation, provenance validation, model leaderboards, and aggregate findings.
---

# Benchmark blind LOT models with subagents

Keep preparation, validation, selection, and evaluation in Python. Use subagents
only for blinded classification. Never expose reviewer/COTA answers or earlier
results to workers.

## Protocol

Use three phases:

- Screening: target fold 0 for inexpensive elimination of weak model/prompt
  ideas. Use exact provenance and sanitized bundles. Never report screening as
  final performance.
- Development: target folds 0–2; retrieval examples also restricted to folds
  0–2 while excluding each target's fold. Compare prompts, models, reasoning,
  retrieval k, and repeats here.
- Final: target folds 3–4; retrieval examples restricted to folds 0–2. Freeze
  the chosen development condition before running this phase. Never retune or
  select a winner using final results.

Primary model-selection metric: mean generated-total exact accuracy across
repeats. Tie-break with lower MAE, then higher within-one accuracy. Report
abstention, usable coverage, and three-way routing separately.

Require an exact model label for screening, development, and final runs. Use
`account-default` only for smoke tests. Run at least two development repeats and
three final repeats per condition when account limits permit.

## Prepare a sanitized run

Choose a new external bundle path outside this repository:

```bash
python3 src/py/prepare_blind_lot_agent_run.py \
  --folds 0,1,2 --retrieval-training-folds 0,1,2 \
  --run-purpose development --repeat-index 1 \
  --retrieval-k 3 --shards 4 --orchestrator codex \
  --model-label <exact-model> --reasoning-effort <exact-effort> \
  --worker-bundle-dir /tmp/blind-lot-<condition>-r1
```

For a final run, use `--folds 3,4`, `--run-purpose final`, no `--limit`,
and the same `--retrieval-training-folds 0,1,2`.

The preparation command must make zero provider calls. Record both the run
directory and external worker-bundle directory.

## Run workers in a separate chat

Start a new Codex chat with its working directory set to the external bundle.
Invoke:

```text
Use $blind-lot-bundle to process every shard. Wait for all workers and report
the attempt count and agent ID for each shard. Do not inspect anything outside
this bundle.
```

The bundle contains prompt packets, schemas, a pinned worker configuration, and
its own skill, but no evaluation data. Do not run workers from the main
repository for development or final benchmarks.

Workers must include truthful `worker_provenance`. Preserve failed attempts
before retrying; never silently repair an invalid output.

## Import, validate, and evaluate

Back in the repository chat, run:

```bash
python3 src/py/finalize_blind_lot_agent_run.py \
  --run-dir <run-directory> \
  --worker-bundle-dir /tmp/blind-lot-<condition>-r1
```

Stop on packet hash changes, incomplete coverage, duplicate cases, schema
errors, unpermitted retrieval citations, model/effort mismatches, or missing
provenance. Never overwrite a finalized run.

## Compare conditions

Use `--phase screening` to rank broad screening candidates. After every planned
development finalist and repeat is complete:

```bash
python3 src/py/compare_blind_lot_model_runs.py \
  --phase development --expected-repeats <n> \
  --run-dir <run-1> --run-dir <run-2> \
  --output artifacts/benchmarks/comparisons/<development-comparison>
```

Use the selected development condition unchanged for final repeats. Build the
final report with `--phase final`; treat its ranking as descriptive only.

## Report

Read only public aggregate/leaderboard files plus finalization provenance.
Report exact model and effort, repeat count, cohort/folds, retrieval training
folds, prompt and knowledge versions, exact accuracy, MAE, within-one accuracy,
abstention, usable coverage, and three-way coverage/accuracy. Separate
development selection from locked final performance.

State that Python made zero direct provider API requests while Codex subagents
used hosted inference through the signed-in account. Treat `codex-subagents` as
a distinct execution provider, not an API-identical replication.

Do not commit or publish restricted artifacts or worker bundles.

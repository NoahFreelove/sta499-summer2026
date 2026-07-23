---
name: blind-lot-subagents
description: Run and compare benchmark-grade order-only multiple-myeloma LOT experiments with Claude Code subagents instead of direct provider API calls, including sanitized worker bundles, repeated development searches, locked final evaluation, provenance validation, and aggregate leaderboards.
---

# Benchmark blind LOT models with Claude Code subagents

Use fold 0 to screen broad ideas, folds 0–2 for development finalists, and folds
3–4 for locked final evaluation.
Restrict retrieval examples to folds 0–2 in both phases. Select models, prompts,
reasoning, and retrieval k only on development results. Freeze the chosen
condition before final evaluation.

Require exact model labels outside smoke tests. Use at least two development
repeats and three final repeats when limits permit. Select by mean generated-total
exact accuracy, breaking ties with lower MAE and higher within-one accuracy.

Prepare each development run with:

```bash
python3 src/py/prepare_blind_lot_agent_run.py \
  --folds 0,1,2 --retrieval-training-folds 0,1,2 \
  --run-purpose development --repeat-index 1 \
  --retrieval-k 3 --shards 4 --orchestrator claude-code \
  --model-label <exact-model> --reasoning-effort <exact-effort> \
  --worker-bundle-dir /tmp/blind-lot-<condition>-r1
```

For final runs, use folds 3–4, `--run-purpose final`, no patient limit, and the
same retrieval-training folds.

Start a separate Claude Code chat rooted at the external bundle and invoke
`$blind-lot-bundle`. Do not run workers from the main repository. The bundle
contains no evaluation answers and pins the worker model. Preserve failed
attempts and include worker provenance.

Back in the repository, import and finalize:

```bash
python3 src/py/finalize_blind_lot_agent_run.py \
  --run-dir <run-directory> --worker-bundle-dir <bundle-directory>
```

Compare completed conditions with `compare_blind_lot_model_runs.py`, using
`--phase development` for selection and `--phase final` for report-only locked
results. Read only public outputs and finalization provenance when summarizing.

Stop on missing exact model provenance, packet changes, incomplete coverage,
schema errors, retrieval-citation violations, or worker model mismatches. Never
expose evaluation data to workers, overwrite finalized runs, or publish
restricted artifacts.

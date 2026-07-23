# Order-only blind AI LOT benchmark

This experiment evaluates an independent LLM vote alongside the unchanged modified rule algorithm and COTA LOT. It is intentionally limited to the 136 adjudicated COTA patients and does not run the 557 unadjudicated patients. It does not provide a production API, UI, continuous learning, deployment, or vector database.

The primary policy auto-accepts only non-abstained cases where algorithm, blind AI, and COTA totals all agree; every other case routes to human review. Reviewer consensus is ground truth. COTA is an independent vote and is never a training label. Single-source and pairwise policies are secondary exploratory analyses added after inspecting fold-0 outcomes; they do not redefine the original endpoint and must be prespecified before use on untouched folds.

## What “order-only” means

The AI receives only the flattened sequence of normalized drug sets. It does not receive vendor line numbers or dates, COTA LOT, reviewer answers for the target, algorithm output or flags, exclusion-group identifiers, or raw patient identifiers. Retrieved training examples may contain a pseudonymous case ID, the same blind trajectory representation, and reviewer-consensus patient-total LOT. The source workbooks do not support reliable reviewer transition labels independent of COTA boundaries, so none are fabricated; patient-total exact accuracy is primary and transition output is structurally validated explanatory output.

The benchmark cannot apply or infer a 28-day regimen-build window, gap rules, date overlap, progression/relapse timing, conditioning adjacency, post-transplant maintenance timing, CAR-T lymphodepletion timing, dose-only changes, route-only changes, or treatment intent. Timing- or intent-dependent boundaries should be marked `INSUFFICIENT_INFORMATION`, which forces abstention.

The compact runtime knowledge is versioned in `src/py/blind_lot/knowledge/order_only_v1.json`. The stable prompt is separately versioned in `src/py/blind_lot/prompts/order_only_v1.txt`. Neither reproduces `textbook_algo_cota.py`.

## Install and entrypoints

```bash
python3 -m pip install -r requirements.txt
python3 src/py/run_blind_lot_benchmark.py --help
```

Every invocation creates a new `runs/<run-id>` directory; prior experiments are never overwritten. Cache reads are opt-in with `--use-cache`, while validated new responses are always cached for later resumption. The benchmark automatically loads credentials from `.env` in the project root; values explicitly exported in the shell take precedence. Copy `.env.example` to `.env` and fill in `OPENAI_API_KEY`, `ANTHROPIC_API_KEY`, and either `KIMI_API_KEY` or the official `MOONSHOT_API_KEY` name for Kimi. `.env` is Git-ignored, and no key is stored in code or artifacts. OpenAI Responses API storage is disabled. `--temperature` and `--reasoning-effort` are optional and should be omitted when the selected model does not support them.

The real-model providers are selected with `--provider openai`, `--provider claude`, or `--provider kimi`. For example:

```bash
python3 src/py/run_blind_lot_benchmark.py \
  --fold 0 --limit 5 --retrieval-k 0 \
  --provider claude --model <claude-model>

python3 src/py/run_blind_lot_benchmark.py \
  --fold 0 --limit 5 --retrieval-k 0 \
  --provider kimi --model <kimi-model>
```

## Coding-agent subagent handoff

The repository also supports a two-phase, file-based handoff for Codex or Claude
Code subagents. The Python process makes no model-provider API requests in this
mode. The signed-in coding-agent client still uses hosted inference and its own
account limits.

Prepare immutable blind prompt shards:

```bash
python3 src/py/prepare_blind_lot_agent_run.py \
  --fold 0 --limit 5 --retrieval-k 0 --shards 2 \
  --orchestrator codex --model-label account-default
```

Invoke the repository skill in a Codex chat with `$blind-lot-subagents`, or the
same-named skill in Claude Code. After the workers create every declared output
shard, validate and evaluate the frozen predictions:

```bash
python3 src/py/finalize_blind_lot_agent_run.py \
  --run-dir artifacts/benchmarks/runs/<run-id>
```

Prompt packets are hashed, shard coverage must be exact, retrieved-example
citations are allowlisted, and strict response validation occurs before reviewer
or COTA answers are loaded. Coding-agent runs use the distinct provider labels
`codex-subagents` and `claude-code-subagents`; they are not API-identical
replications of direct-provider runs.

### Benchmark-grade model selection

Use fold 0 for inexpensive screening, folds 0–2 for development, and folds 3–4
for locked final evaluation. Screening eliminates weak ideas but is never
reported as final performance.
Restrict retrieval candidates to folds 0–2 in both phases. This permits prompt,
model, reasoning, and retrieval tuning on development cases without using final
outcomes for selection.

Development and final subagent runs require an exact model label and a sanitized
worker bundle outside the repository. The bundle includes only blind packets,
schemas, a pinned worker configuration, and a bundle-local skill. It contains no
reviewer/COTA evaluation records.

Example development preparation:

```bash
python3 src/py/prepare_blind_lot_agent_run.py \
  --folds 0,1,2 --retrieval-training-folds 0,1,2 \
  --run-purpose development --repeat-index 1 \
  --retrieval-k 3 --shards 4 --orchestrator codex \
  --model-label <exact-model> --reasoning-effort <exact-effort> \
  --worker-bundle-dir /tmp/blind-lot-<condition>-r1
```

Start a separate coding-agent chat rooted at the printed bundle directory and
invoke `$blind-lot-bundle`. Then import, validate, and evaluate:

```bash
python3 src/py/finalize_blind_lot_agent_run.py \
  --run-dir <run-directory> \
  --worker-bundle-dir /tmp/blind-lot-<condition>-r1
```

Repeat each development condition at least twice when account limits permit.
Build the development leaderboard with:

```bash
python3 src/py/compare_blind_lot_model_runs.py \
  --phase development --expected-repeats <n> \
  --run-dir <run-1> --run-dir <run-2> \
  --output artifacts/benchmarks/comparisons/<development-comparison>
```

Selection uses mean generated-total exact accuracy, then lower MAE and higher
within-one accuracy as tie-breakers. Freeze the selected condition and run it on
folds 3–4 with `--run-purpose final`, no patient limit, and three repeats when
feasible. A final-phase leaderboard is report-only and must not trigger retuning.

## Reproducible commands

Run these from the repository root.

1. Dry-run leakage validation, with no API calls:

```bash
python3 src/py/run_blind_lot_benchmark.py \
  --dry-run --fold 0 --limit 5 --retrieval-k 0
```

2. Five-patient deterministic mock-provider test:

```bash
python3 src/py/run_blind_lot_benchmark.py \
  --fold 0 --limit 5 --retrieval-k 0 --provider mock
```

3. Five-patient real-model smoke test with no retrieval:

```bash
python3 src/py/run_blind_lot_benchmark.py \
  --fold 0 --limit 5 --retrieval-k 0 \
  --provider openai --model <configured-model> --reasoning-effort <configured-effort>
```

4. Five-patient real-model smoke test with three retrieved examples:

```bash
python3 src/py/run_blind_lot_benchmark.py \
  --fold 0 --limit 5 --retrieval-k 3 \
  --provider openai --model <configured-model> --reasoning-effort <configured-effort>
```

5. One complete held-out fold, only after inspecting smoke outputs:

```bash
python3 src/py/run_blind_lot_benchmark.py \
  --fold 0 --retrieval-k 3 \
  --provider openai --model <configured-model> --reasoning-effort <configured-effort>
```

6. All five folds, only with explicit approval after one fold succeeds:

```bash
python3 src/py/run_blind_lot_benchmark.py \
  --all-folds --retrieval-k 3 \
  --provider openai --model <configured-model> --reasoning-effort <configured-effort>
```

7. Reaggregate a completed run using saved joined records only. This command does not initialize a provider, execute retrieval, issue model requests, or modify predictions. It preserves the source experiment metadata and writes suffixed evaluation artifacts:

```bash
python3 src/py/evaluate_blind_lot_run.py \
  --run-dir artifacts/benchmarks/runs/<run-id>
```

8. Compare completed runs after reaggregation. Folds, model, reasoning effort, prompt and knowledge versions, input hashes, and the evaluation cohort must match; retrieval `k` must differ. The output directory receives `policy_comparison.json`, `policy_comparison.csv`, and `vote_pattern_comparison.csv`:

```bash
python3 src/py/compare_blind_lot_runs.py \
  --run-dir artifacts/benchmarks/runs/<k0-run> \
  --run-dir artifacts/benchmarks/runs/<k3-run> \
  --run-dir artifacts/benchmarks/runs/<k5-run> \
  --output artifacts/benchmarks/comparisons/fold0-k-comparison
```

9. Resume using the shared validated cache (a new run directory is still created):

```bash
python3 src/py/run_blind_lot_benchmark.py \
  --fold 0 --retrieval-k 3 --use-cache \
  --provider openai --model <configured-model> --reasoning-effort <configured-effort>
```

10. Locate results:

```text
artifacts/benchmarks/runs/<run-id>/public/aggregate_evaluation.json
artifacts/benchmarks/runs/<run-id>/restricted/experiment_metadata.json
artifacts/benchmarks/runs/<run-id>/restricted/blind_predictions.jsonl
artifacts/benchmarks/runs/<run-id>/restricted/retrieval_debug.jsonl
artifacts/benchmarks/runs/<run-id>/restricted/joined_evaluation.jsonl
artifacts/benchmarks/runs/<run-id>/restricted/joined_evaluation.reaggregated.jsonl
artifacts/benchmarks/runs/<run-id>/restricted/evaluation_reaggregation.json
artifacts/benchmarks/runs/<run-id>/public/aggregate_evaluation.reaggregated.json
```

Dry runs instead write `public/dry_run_report.json` and restricted prompt previews. Public files contain aggregates only. Restricted predictions contain pseudonymous case IDs, fold, the validated model response, permitted retrieved case IDs, cache information, and provider-attempt counts. The restricted evaluation join is created only after `blind_predictions.jsonl` has been finalized.

## Validation and uncertainty

The response contract rejects extra fields, malformed enums, missing transitions, nonconsecutive transition indices, inconsistent LOT totals, and insufficient-information decisions without abstention. Cached responses pass the current contract again when loaded. Retrieval deterministically combines length, transition-type sequence, aligned regimen-set Jaccard, transplant, CAR-T, bispecific, unknown/investigational, steroid-only, and reappearance features. Targets, their exclusion groups, and the entire held-out fold are ineligible.

Aggregate reports separate generated-total diagnostics from selective AI-vote performance. Generated-total coverage uses all eligible patients as its denominator; its accuracy and error metrics use patients with a valid integer total, including numeric totals emitted with abstention. Selective coverage uses all eligible patients; selective accuracy and error metrics use only valid, non-abstained AI votes. AI pairwise and three-way agreement coverage also use all eligible patients, while reviewer accuracy uses only the corresponding usable agreement subset. `metrics.ai_alone` remains as a deprecated, explicitly labeled alias of `metrics.generated_total`; it must not be interpreted as non-abstained accuracy.

Abstained and invalid outputs always route to human review. Three-way consensus requires a valid non-abstained AI vote matching both the unchanged algorithm and COTA votes. Restricted joined records include the derived prediction status, agreement flags, and deterministic routing reason. Retrieval debug records preserve the exact already-rendered demonstrations and retrieval context with SHA-256 provenance; generating these diagnostics does not alter ranking or model input.

## Routing policy and vote-pattern analysis

Every policy uses all joined, reviewer-adjudicated patients as its coverage denominator. Algorithm-only and COTA-only require a valid non-negative integer from that source. Usable-AI-only requires a valid non-abstained AI vote. Pairwise policies require both named votes to be valid and equal; AI-containing policies additionally require a usable AI vote. The primary three-way policy requires all three valid votes to agree. Accuracy and false-accept rate use accepted patients as their denominator, and are `null` when no patient is accepted. Review proportion uses all eligible patients.

Vote patterns are mutually exclusive and exhaustive. With a usable AI vote and valid algorithm/COTA votes they are: all three agree; algorithm+AI agree; COTA+AI agree; algorithm+COTA agree; or all three differ. Valid numeric AI abstentions are split by whether algorithm and COTA agree. Invalid AI output is separate from abstention. Two additional missing-non-AI-vote categories keep malformed evaluation inputs exhaustive. Each stratum reports source accuracy, exact reviewer matches, unique winners, and majority-versus-dissent attribution where two votes agree.

The conditional sections partition algorithm+AI agreement by COTA agreement, algorithm+COTA agreement by AI agreement/disagreement/abstention/invalidity, and COTA+AI agreement by algorithm agreement. These are retrospective diagnostics of veto behavior, not a production-rule selection procedure.

Percentile 95% confidence intervals use patient-level nonparametric bootstrap resampling of the full eligible cohort with replacement, a fixed seed of `4992026`, and 2,000 replicates. Each metric applies its subset filter after resampling. Undefined empty-subset replicates are omitted and reported. Small consensus subsets must not be overinterpreted.

## Later controlled full-fold retrieval comparison

After selecting one complete held-out fold and freezing model, reasoning effort, prompt version, knowledge version, inputs, bootstrap settings, policy definitions, decision criteria, and the handling of multiple comparisons, run all three conditions with identical arguments except `--retrieval-k`. Do not combine smoke tests with the controlled comparison.

```bash
COMMON="--fold <fold> --provider openai --model <model> --reasoning-effort <effort> --bootstrap-seed 4992026 --bootstrap-replicates 2000"
python3 src/py/run_blind_lot_benchmark.py $COMMON --retrieval-k 0
python3 src/py/run_blind_lot_benchmark.py $COMMON --retrieval-k 3
python3 src/py/run_blind_lot_benchmark.py $COMMON --retrieval-k 5

python3 src/py/compare_blind_lot_runs.py \
  --run-dir artifacts/benchmarks/runs/<k0-run> \
  --run-dir artifacts/benchmarks/runs/<k3-run> \
  --run-dir artifacts/benchmarks/runs/<k5-run> \
  --output artifacts/benchmarks/comparisons/fold0-k-comparison
```

## Remaining limitations

- Order alone cannot resolve timing- or intent-dependent boundaries.
- Reviewer transition labels are unavailable without risking COTA-boundary leakage, so transition explanations are not scored as ground truth.
- Retrieved examples expose reviewer patient totals and may influence decomposition despite explicit prompt controls; cross-validation prevents target and overlap leakage but does not remove all few-shot sensitivity.
- Exact duplicate exclusion reduces direct trajectory leakage, but broader clinical similarity is expected and is the purpose of retrieval.
- Confidence intervals quantify patient resampling uncertainty, not model/provider nondeterminism.
- A five-case smoke test validates plumbing, not clinical performance. Full conclusions require the same frozen configuration across all five held-out folds.

## Frozen pooled out-of-fold k=3 versus k=5 analysis

`configs/pooled_oof_k3_k5.json` is the versioned source of truth for the ten frozen
fold runs. It names every run explicitly; pooling never discovers or selects a newer
run. The following commands only read completed metadata, predictions, retrieval
diagnostics, and joined evaluation records. They do not initialize a provider, read a
model cache, execute retrieval, issue a model request, or modify saved predictions.

```bash
env -u OPENAI_API_KEY python3 src/py/pool_blind_lot_runs.py \
  --manifest configs/pooled_oof_k3_k5.json \
  --runs-root artifacts/benchmarks/runs \
  --output artifacts/benchmarks/pooled/frozen-oof

env -u OPENAI_API_KEY python3 src/py/compare_blind_lot_predictions.py \
  --k3-run artifacts/benchmarks/pooled/frozen-oof/k3-all-folds \
  --k5-run artifacts/benchmarks/pooled/frozen-oof/k5-all-folds \
  --output artifacts/benchmarks/comparisons/frozen-oof-k3-v-k5
```

Pooling validates all expected `(retrieval_k, fold)` pairs, per-fold row counts and
fold labels, cohort identity, frozen configuration, and byte hashes for source
metadata, predictions, retrieval diagnostics, and joined rows. Fold-0 source runs
predate evaluation 1.2.0; their saved predictions are reaggregated through the current
1.2.0 evaluator, and both original and effective versions are preserved in provenance.
The frozen regression assertions stop publication unless the following values recur:

- Algorithm+COTA: 76/136 accepted, 63 correct and 13 incorrect.
- k=3: generated total 45/136 correct; 38/136 usable and 18/38 correct;
  algorithm+AI 14/16 correct; COTA+AI 12/15; three-way 9/10.
- k=5: generated total 48/136 correct; 47/136 usable and 21/47 correct;
  algorithm+AI 16/21 correct; COTA+AI 14/16; three-way 10/11.

Patient-level joined, paired, provenance, and diagnostic files are restricted.
Public JSON and CSV files contain aggregates only and are checked for identifier
markers and absolute local paths. Fold 0 was development-influenced. Pairwise policies
remain exploratory, and small accepted subsets must not be overinterpreted. Strict
three-way agreement had one false acceptance under both retrieval settings; it did not
have perfect precision.

Restricted error inspection uses saved responses and diagnostics only:

```bash
python3 src/py/inspect_blind_lot_errors.py \
  --run-dir artifacts/benchmarks/pooled/frozen-oof/k3-all-folds \
  --compare-run-dir artifacts/benchmarks/pooled/frozen-oof/k5-all-folds \
  --policy three_way_agreement --incorrect-only \
  --output artifacts/benchmarks/diagnostics/three-way-false-accepts
```

The default test suite uses committed synthetic fixtures and needs no ignored
artifacts. Local restricted regression tests are opt-in:

```bash
env -u OPENAI_API_KEY python3 -m unittest discover -s tests -v
RUN_RESTRICTED_REGRESSION_TESTS=1 \
  env -u OPENAI_API_KEY python3 -m unittest discover -s tests -v
```

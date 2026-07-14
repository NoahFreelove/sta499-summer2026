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

Every invocation creates a new `runs/<run-id>` directory; prior experiments are never overwritten. Cache reads are opt-in with `--use-cache`, while validated new responses are always cached for later resumption. OpenAI access uses `OPENAI_API_KEY`; no key is stored in code or artifacts, and Responses API storage is disabled. `--temperature` is optional and should be omitted for models that do not support it.

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

# Agentic COTA Line-of-Therapy Row Misclassification Detection

## Overview

This project uses the OpenAI Agents SDK to identify **individual COTA line-of-therapy rows whose split/merge pattern is likely incorrect**.

The main research objective is:

> For each transition between vendor LoT rows, did COTA correctly start a new line, or should that row have remained part of the previous line?

The system uses:

- a Worker agent;
- an independent Auditor agent;
- historical doctor-adjudicated COTA examples;
- transition-aware split/merge comparison;
- a review-routing policy for suspicious cases.

The goal is to catch as many truly misclassified row transitions as possible while reducing unnecessary physician review.

## Why row-pattern scoring matters

Absolute LoT numbers can overcount errors because one earlier numbering difference may carry forward.

| Vendor LoT | 1 | 2 | 3 | 4 |
|---|---:|---:|---:|---:|
| Doctor LoT | 1 | 2 | 4 | 5 |

A simple number comparison would mark rows 3 and 4 as wrong. But both COTA and the doctor still chose to split at the later transitions. The actual disagreement occurred only at the earlier boundary.

The pipeline therefore converts each patient into:

```text
start
split
split
merge
...
```

Each transition is classified as:

- `split` - start a new line;
- `merge` - keep the row in the previous line;
- `start` - first line for the patient.

This prevents downstream numbering carryover from being counted as repeated mistakes.

## Research goal

The pipeline is evaluated mainly on:

1. How accurately the agents reproduce the doctor's true split/merge pattern.
2. How many truly misclassified transitions are caught by the review process.

The primary metrics are:

- Worker split/merge accuracy;
- Auditor split/merge accuracy;
- total misclassified transitions;
- misclassified transition catch rate;
- misclassified transitions missed through auto-acceptance.

Patient-level results are secondary.

## Latest evaluation

The latest held-out run used:

```text
116 training patients
20 evaluation patients
114 total split/merge transitions
```

### Row-pattern results

| Metric | Result |
|---|---:|
| Worker correct transitions | 98 of 114 |
| Worker split/merge accuracy | 85.96% |
| Auditor correct transitions | 100 of 114 |
| Auditor split/merge accuracy | 87.72% |
| Truly misclassified transitions | 13 |
| Misclassified transitions caught through review | 13 of 13 |
| Misclassified transitions missed through auto-acceptance | 0 |
| Misclassified transition catch rate | 100% |

The Auditor was the stronger row-pattern classifier in this run.

### Secondary patient-level results

| Metric | Result |
|---|---:|
| Worker exact patient-pattern accuracy | 11 of 20, 55% |
| Auditor exact patient-pattern accuracy | 14 of 20, 70% |
| Misclassified patients caught for review | 5 of 5 |
| Correct patients unnecessarily reviewed | 12 of 15 |
| Overall human-review rate | 85% |
| Safe auto-accept rate | 100% |

The current pipeline is highly sensitive to row-level errors, but still conservative because many correct patients are routed to review.

## Pipeline flow

```text
New COTA patient rows
        |
        v
Preprocess and normalize treatment lines
        |
        v
Retrieve similar historical labeled patients
        |
        +----------------------+
        |                      |
        v                      v
Worker Agent             Auditor Agent
        |                      |
        +----------+-----------+
                   |
                   v
Convert outputs into split/merge patterns
                   |
                   v
Compare transition by transition
                   |
        +----------+----------+
        |                     |
        v                     v
Safe agreement          Suspicious row,
auto-accept             low confidence, or
                        disagreement
                              |
                              v
                       Human review
```

## How it works

### 1. Preprocessing

Each patient's COTA rows are converted into a structured sequence containing:

- vendor LoT number;
- regimen;
- start date;
- end date;
- discontinuation reason.

Continuation rows in Excel are collapsed into the preceding LoT row.

### 2. Historical reference retrieval

The old doctor-adjudicated COTA workbook is used as a reference library.

For each new patient, the system retrieves similar historical sequences based on:

- overlapping drugs;
- drug additions and removals;
- regimen reductions and expansions;
- steroid-only changes;
- identical adjacent regimens;
- previous oversplitting patterns.

During evaluation, held-out patients are excluded from the reference library to prevent leakage.

### 3. Worker agent

The Worker independently predicts whether each vendor row should:

- remain a new LoT;
- merge with the previous LoT;
- remain uncertain.

It also returns rationale and confidence.

### 4. Auditor agent

The Auditor receives the same patient and historical examples but does not see the Worker output.

It independently produces its own split/merge pattern.

### 5. Transition-aware comparison

The system compares split/merge decisions rather than raw LoT numbers.

```text
Worker:  start, split, split, merge
Auditor: start, split, split, merge
```

This counts as agreement even if absolute corrected LoT numbers differ because of numbering carryover.

### 6. Review routing

A case is sent to review when:

- Worker and Auditor disagree on a transition;
- either agent explicitly requests review;
- a row is marked uncertain;
- confidence is below the configured threshold.

Generic descriptive uncertainty does not automatically force review unless it affects a split/merge boundary.

## Project structure

```text
Agentic_OpenAI_LoT_Adj/
├── .env
├── requirements.txt
├── config.py
├── schemas.py
├── preprocessing.py
├── reference_retrieval.py
├── lot_rules.py
├── agent_definitions.py
├── decision.py
├── pipeline.py
├── eval_pipeline.py
├── eval_metrics.py
├── data/
│   ├── new_cota_data(1).xlsx
│   └── LoT Adjudication Datasets(6).xlsx
├── output/
└── eval_output/
```

## Main files

- `preprocessing.py` - cleans raw COTA rows.
- `reference_retrieval.py` - builds the labeled reference library and retrieves similar cases.
- `agent_definitions.py` - defines Worker and Auditor agents.
- `decision.py` - converts outputs into transition patterns and routes cases.
- `eval_pipeline.py` - runs held-out evaluation without leakage.
- `eval_metrics.py` - calculates row-level and patient-level metrics.

## Setup

```powershell
python -m venv .venv
Set-ExecutionPolicy -Scope Process -ExecutionPolicy Bypass
.\.venv\Scripts\Activate.ps1
.\.venv\Scripts\python.exe -m pip install -r requirements.txt
```

Create `.env`:

```text
OPENAI_API_KEY=your_key_here
OPENAI_MODEL=your_model_here
```

## Run the evaluation

```powershell
.\.venv\Scripts\python.exe eval_pipeline.py
```

Outputs:

```text
eval_output/
├── split.json
├── evaluation_details.csv
├── evaluation_metrics.json
└── patients/
```

## Important metrics

### `worker_split_merge_accuracy`

Percentage of row transitions where the Worker matched the doctor.

### `auditor_split_merge_accuracy`

Percentage of row transitions where the Auditor matched the doctor.

### `true_cota_misclassified_transitions`

Number of transitions where COTA and the doctor disagreed on split versus merge.

### `misclassified_transitions_caught_via_patient_review`

Number of true row errors contained in patients routed to review.

### `misclassified_transitions_missed_via_autoaccept`

Number of true row errors inside patients that were auto-accepted.

### `misclassified_transition_catch_rate`

```text
caught misclassified transitions / all misclassified transitions
```

## Strategies tested

| Strategy | Row-error result | Review behavior |
|---|---|---|
| Very permissive triage | Missed many row errors | Low review rate |
| Very strict risk gate | Caught all row errors | Reviewed almost every patient |
| Lower confidence threshold | Missed some row errors | Review remained high |
| Hand-weighted risk score | Missed most row errors | Too many unsafe auto-accepts |
| Transition-aware Boolean policy | Caught all known row errors in best runs | Best safety/coverage balance so far |

The weighted-risk version was rejected because agent confidence scores were not calibrated probabilities.

## Best current interpretation

The strongest current result is:

> The Auditor correctly classified 100 of 114 split/merge transitions, or 87.72%.

The review process caught:

```text
13 of 13 truly misclassified transitions
```

with:

```text
0 misclassified transitions missed through auto-acceptance
```

The main remaining weakness is specificity: many correct patients are still sent to review.

## Comparison with Noah's pipeline

Noah's pipeline uses one AI prediction as a third vote alongside the algorithm and COTA. It mainly predicts same line, new line, insufficient information, and total LoT count.

This pipeline differs because it focuses directly on row-level split/merge decisions.

| Noah's pipeline | This pipeline |
|---|---|
| One AI vote | Independent Worker and Auditor |
| Focuses heavily on total LoT count | Focuses on each row transition |
| Uses drug order only | Uses vendor rows, regimens, dates, and available context |
| Measures total-count agreement | Measures boundary correctness |
| AI acts mainly as a conservative third vote | Agents identify suspicious row transitions |

The key advantage is that this pipeline evaluates the exact boundary decision the doctor makes. It can identify which vendor row is likely wrong and does not count numbering carryover as repeated errors.

## Current limitations

- The latest evaluation contains only 20 held-out patients.
- Worker and Auditor may have correlated errors.
- Agent confidence values are not calibrated probabilities.
- Historical retrieval quality has not been evaluated independently.
- Many correct patients are still sent to review.
- Final performance should be measured with grouped cross-validation across all 136 patients.

## Next steps

1. Run grouped five-fold cross-validation across all patients.
2. Aggregate out-of-fold split/merge accuracy.
3. Measure row-error recall and false-negative rate.
4. Analyze which transition types are most often missed.
5. Compare Auditor-only, Worker-only, and combined routing.
6. Reduce unnecessary review without lowering row-error catch rate.
7. Evaluate retrieval quality for oversplitting patterns.

Primary optimization target:

```text
Maintain near-100% misclassified-transition catch rate
while reducing unnecessary patient review.
```

## Current status

```text
Auditor split/merge accuracy: 100/114 = 87.72%
Worker split/merge accuracy: 98/114 = 85.96%
True misclassified transitions caught: 13/13 = 100%
True misclassified transitions missed: 0
```

Further grouped validation is required before making general performance claims.

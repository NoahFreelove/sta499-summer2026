# Agentic COTA Line-of-Therapy Adjudication

## Overview

This project uses the OpenAI Agents SDK to identify COTA patients whose vendor-assigned lines of therapy are likely to be misclassified and should be sent for physician review.

The system is not intended to replace physician adjudication. Its purpose is to reduce manual review by automatically accepting straightforward cases while preserving a high catch rate for likely COTA errors.

The current architecture uses:

- a **Worker agent** to independently assess each patient;
- an **Auditor agent** to independently review the same patient;
- historical doctor-adjudicated COTA cases as retrieved reference examples;
- deterministic split/merge comparison logic;
- a conservative routing policy that sends suspicious cases to human review.

---

## Research Objective

The central question is:

> Can an agentic workflow identify the COTA patients most likely to contain line-of-therapy errors, so that physicians review only a smaller, higher-risk subset?

The main performance goal is not exact reproduction of every doctor label. The main goal is triage:

1. Catch as many truly misclassified COTA patients as possible.
2. Avoid auto-accepting patients whose COTA classifications are wrong.
3. Reduce unnecessary review of patients whose original COTA classification was already correct.

The most important evaluation metrics are therefore:

- **Misclassification catch rate**
- **Missed-error rate**
- **Safe auto-accept rate**
- **Unnecessary review rate among correct patients**
- **Overall human-review rate**

---

## Why Transition-Level Scoring Is Necessary

A direct comparison of absolute LoT numbers can overcount errors because one earlier numbering difference may carry forward into later rows.

Example:

| Vendor LoT | 1 | 2 | 3 | 4 |
|---|---:|---:|---:|---:|
| Doctor LoT | 1 | 2 | 4 | 5 |

A row-level number comparison would count the third and fourth lines as mismatches. However, both the vendor and doctor still decided to split at both later transitions.

The pipeline therefore converts each patient into transition decisions:

```text
start
split
split
merge
...
```

It evaluates whether the vendor and doctor made the same split/merge decision at each boundary. Downstream numbering carryover is not counted as a new error.

---

## Pipeline Architecture

```text
New COTA patient
        |
        v
Preprocessing
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
Compare split/merge decisions
                   |
        +----------+----------+
        |                     |
        v                     v
Safe agreement          Boundary risk,
auto-accept             uncertainty, or
                        disagreement
                              |
                              v
                       Human review
```

---

## How the Pipeline Works

### 1. Preprocess the patient

The new COTA workbook is converted into a structured patient-level representation.

Each patient contains vendor-assigned LoT rows with fields such as:

- patient ID;
- vendor LoT number;
- regimen;
- start date;
- end date;
- discontinuation reason.

Continuation rows in Excel are collapsed into the correct preceding regimen.

Sensitive identifiers should remain pseudonymous.

### 2. Build the historical reference library

The old doctor-adjudicated COTA workbook is loaded from the `Cota` sheet.

For each historical patient, the reference library stores:

- vendor LoT number;
- regimen;
- adjudicated LoT;
- whether the vendor line was kept or merged;
- treatment transition patterns;
- whether an oversplit correction occurred.

The reference library is represented using TF-IDF features and cosine similarity.

### 3. Retrieve similar historical cases

For every new patient, the pipeline retrieves a small set of similar historical patients.

Similarity considers:

- overlapping drugs;
- drug additions;
- drug removals;
- regimen reductions;
- regimen expansions;
- steroid-only changes;
- identical adjacent regimens;
- prior oversplitting patterns.

During evaluation, the held-out patient is excluded from the reference library to prevent leakage.

### 4. Worker agent

The Worker independently reviews the patient and predicts the corrected split/merge sequence.

The Worker returns structured output including:

- corrected vendor-to-adjudicated LoT mapping;
- action for each vendor LoT;
- confidence;
- rationale;
- uncertainty;
- historical references used.

### 5. Auditor agent

The Auditor receives the same patient and reference examples but does not see the Worker output.

It independently assesses the same boundaries and produces its own structured adjudication.

### 6. Deterministic comparison

The Python decision layer compares Worker and Auditor outputs using transition decisions rather than raw LoT numbers.

It checks:

- split versus merge decisions;
- missing vendor LoTs;
- low-confidence lines;
- explicit uncertainty;
- whether either agent requests human review.

### 7. Final routing

The current best-performing policy automatically accepts a patient only when:

- Worker and Auditor agree on all split/merge transitions;
- neither result contains a meaningful boundary-level review trigger;
- neither result explicitly requests human review;
- no line is marked uncertain;
- confidence requirements are met.

Generic descriptive uncertainty does not automatically trigger review unless it affects a treatment boundary.

---

## Project Structure

```text
Agentic_OpenAI_LoT_Adj/
├── .env
├── .gitignore
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
├── validate_setup.py
├── data/
│   ├── new_cota_data(1).xlsx
│   └── LoT Adjudication Datasets(6).xlsx
├── output/
└── eval_output/
```

---

## Main Files

### `config.py`

Stores:

- workbook paths;
- sheet names;
- model name;
- output directories;
- reference count;
- confidence thresholds;
- concurrency settings.

### `schemas.py`

Defines structured Pydantic models for:

- corrected LoT lines;
- adjudication results;
- audit critique;
- pipeline decisions.

### `preprocessing.py`

Cleans new COTA data and converts it into patient-level structured input.

### `reference_retrieval.py`

Loads the old labeled COTA dataset, builds the reference library, and retrieves similar historical patients.

### `lot_rules.py`

Contains the LoT adjudication rules used by the agents.

### `agent_definitions.py`

Defines the Worker, Auditor, and optional Critic agents using the OpenAI Agents SDK.

### `decision.py`

Implements transition-aware Worker/Auditor comparison and the final review-routing policy.

### `pipeline.py`

Runs the production workflow on the new COTA workbook.

### `eval_pipeline.py`

Creates a held-out evaluation split, removes evaluation patients from the reference library, runs the full pipeline, and saves results.

### `eval_metrics.py`

Calculates both:

- adjudication accuracy metrics;
- triage metrics aligned with the research objective.

---

## Setup

### 1. Create a virtual environment

```powershell
python -m venv .venv
```

If PowerShell blocks activation:

```powershell
Set-ExecutionPolicy -Scope Process -ExecutionPolicy Bypass
.\.venv\Scripts\Activate.ps1
```

Activation is optional. The virtual environment can also be used directly:

```powershell
.\.venv\Scripts\python.exe
```

### 2. Install dependencies

```powershell
.\.venv\Scripts\python.exe -m pip install --upgrade pip
.\.venv\Scripts\python.exe -m pip install -r requirements.txt
```

Expected packages include:

```text
openai-agents
pandas
pydantic
python-dotenv
openpyxl
scikit-learn
```

The package is installed as `scikit-learn` but imported as `sklearn`.

### 3. Add the OpenAI API key

Create `.env`:

```text
OPENAI_API_KEY=your_key_here
OPENAI_MODEL=your_model_here
```

Do not commit `.env`.

### 4. Configure workbook paths

In `config.py`, confirm:

```python
NEW_COTA_FILE = BASE_DIR / "data" / "new_cota_data(1).xlsx"
OLD_ADJUDICATED_FILE = (
    BASE_DIR / "data" / "LoT Adjudication Datasets(6).xlsx"
)

NEW_COTA_SHEET = "Sheet1"
OLD_COTA_SHEET = "Cota"
```

Check workbook sheet names:

```powershell
.\.venv\Scripts\python.exe -c "import pandas as pd; print(pd.ExcelFile('data/new_cota_data(1).xlsx').sheet_names)"
```

---

## Running the Production Pipeline

```powershell
.\.venv\Scripts\python.exe pipeline.py
```

Outputs are written under:

```text
output/
├── adjudication_summary.csv
├── human_review_queue.csv
└── <patient_id>/
    ├── preprocessed_input.json
    ├── reference_examples.json
    ├── worker.json
    ├── auditor.json
    ├── comparison.json
    ├── decision.json
    └── selected_result.json
```

`selected_result.json` is only created for auto-accepted patients.

---

## Running the Evaluation

```powershell
.\.venv\Scripts\python.exe eval_pipeline.py
```

The evaluation:

1. splits patients by `cpid`;
2. holds out a patient-level test subset;
3. hides adjudicated labels from the agents;
4. builds the reference library only from training patients;
5. excludes evaluation patients from retrieval;
6. compares predictions against hidden doctor labels.

Outputs:

```text
eval_output/
├── split.json
├── evaluation_details.csv
├── evaluation_metrics.json
└── patients/
```

---

## Evaluation Metrics

### Triage metrics

#### `misclassification_catch_rate`

Of all truly misclassified COTA patients, how many were sent to review?

```text
caught misclassified patients / all misclassified patients
```

#### `missed_error_rate`

Of all truly misclassified patients, how many were incorrectly auto-accepted?

```text
missed misclassified patients / all misclassified patients
```

#### `review_precision`

Of all patients sent to review, how many actually contained a COTA error?

#### `unnecessary_review_rate_among_correct_patients`

Of all patients whose original COTA classification was correct, how many were still sent to review?

#### `safe_autoaccept_rate`

Of all auto-accepted patients, how many were truly safe to accept?

### Secondary accuracy metrics

- Worker exact patient accuracy
- Auditor exact patient accuracy
- Worker LoT classification accuracy
- Auditor LoT classification accuracy
- Worker/Auditor agreement
- Split/merge transition accuracy
- Merge precision, recall, and F1

---

## Strategies Tested

| Strategy | Misclassified patients caught | Correct patients unnecessarily reviewed | Interpretation |
|---|---:|---:|---|
| Initial conservative agreement policy | 100% | 70% | Safe but reviewed many correct cases |
| Permissive triage-focused agents | 40% | 6.7% | Too permissive; missed too many errors |
| Strict risk-gated policy | 100% | 93.3% | Safe but almost no reduction in review |
| Confidence threshold reduced to 0.80 | 80% | 80% | Lower threshold reduced safety without enough review reduction |
| Hand-weighted risk formula | 25% | 20% | Agent confidence was not calibrated; unsafe |
| **Best Boolean transition-aware policy** | **100%** | **60%** | **Best observed safety and coverage trade-off** |

The runs above used different holdout samples, so they are directional comparisons rather than final benchmark results.

---

## Best Result So Far

On the best 20-patient evaluation:

- 5 patients were truly misclassified.
- The pipeline sent all 5 to review.
- 0 misclassified patients were auto-accepted.
- 15 patients were truly correct.
- 6 correct patients were safely auto-accepted.
- 9 correct patients were unnecessarily reviewed.

Key metrics:

```text
Misclassification catch rate: 100%
Missed-error rate: 0%
Safe auto-accept rate: 100%
Correct-patient auto-accept rate: 40%
Overall auto-accept rate: 30%
Overall review rate: 70%
Unnecessary review among correct patients: 60%
```

This was the best observed configuration because it reduced review by 30% without missing a known COTA error in that sample.

---

## Comparison With Noah's Pipeline

Noah's pipeline uses one AI prediction as a third vote alongside:

- the existing algorithm;
- COTA;
- the AI.

It retrieves similar drug-order trajectories and predicts:

- same line;
- new line;
- insufficient information.

It then compares total LoT counts.

Reported results from Noah's pipeline:

| Policy | Accepted | Correct | Accuracy among accepted | Coverage |
|---|---:|---:|---:|---:|
| Algorithm alone | 136 | 89 | 65.4% | 100% |
| Algorithm + COTA agree | 76 | 63 | 82.9% | 55.9% |
| Three-way agreement, k=3 | 10 | 9 | 90.0% | 7.4% |
| Three-way agreement, k=5 | 11 | 10 | 90.9% | 8.1% |

The AI itself had usable standalone accuracy below 50%, so its main value was as a conservative selector rather than an independent expert.

### Key differences

| Noah's pipeline | This pipeline |
|---|---|
| One AI vote | Independent Worker and Auditor |
| Primarily drug-order input | Uses vendor LoTs, regimens, dates, and available context |
| Predicts total LoT count | Predicts split/merge decisions at each boundary |
| Compares three total values | Compares detailed transition decisions |
| Mainly consensus filtering | Specifically detects patients needing physician review |
| About 8% three-way acceptance coverage | Best observed holdout auto-accepted 30% |
| AI often abstains | Worker/Auditor produce structured, auditable reasoning |

This pipeline is better aligned with the operational adjudication task because it focuses on where a boundary is likely wrong rather than only whether total LoT counts agree.

However, a definitive comparison requires grouped cross-validation over all 136 patients.

---

## Current Limitations

### Small evaluation samples

The strongest result currently comes from one 20-patient holdout. Results may vary between splits.

### Shared model dependence

Worker and Auditor may use the same underlying model, rules, and references. Their errors are therefore not fully independent.

### Confidence is not calibrated

The agents' confidence values should not be treated as true probabilities. A hand-weighted risk score performed poorly because the agents were sometimes confidently wrong.

### Historical retrieval quality

Poorly matched historical examples may bias the agents. Retrieval should eventually be validated separately.

### Missing clinical context

The dataset may not include every clinical fact needed for definitive LoT adjudication.

### Generalizability

The current system is tailored to the historical COTA dataset and its labeling conventions.

---

## Recommended Next Steps

1. Run grouped five-fold cross-validation across all 136 patients.
2. Keep similar or duplicate treatment trajectories in the same fold.
3. Report aggregate out-of-fold triage metrics.
4. Analyze why correct patients are still sent to review.
5. Separate line-level and patient-level confidence thresholds.
6. Validate whether retrieved examples are clinically relevant.
7. Compare Worker-only, Auditor-only, and Worker-Auditor policies.
8. Consider training a calibrated logistic regression model on pipeline features only after enough out-of-fold predictions are available.
9. Preserve misclassification catch rate while gradually lowering unnecessary review.

The next optimization target is:

```text
Maintain near-100% misclassification catch rate
while reducing unnecessary review below 60%.
```

---

## Research Interpretation

The current results suggest that the agentic pipeline may be useful as a selective triage layer.

The strongest tested version did not replace physician adjudication. Instead, it:

- caught every known misclassified patient in the holdout;
- safely removed 30% of patients from manual review;
- preserved a complete audit trail for every decision.

The main remaining challenge is improving specificity: safely auto-accepting more already-correct patients without missing true COTA errors.

---

## Privacy and Security

- Use pseudonymous patient IDs.
- Do not send names, addresses, medical record numbers, or direct identifiers.
- Do not commit API keys.
- Keep `.env` in `.gitignore`.
- Store outputs only in approved locations.
- Confirm institutional and contractual approval before processing protected health information through external APIs.

---

## Status

Current best policy:

```text
Transition-aware Worker + Auditor
with deterministic split/merge comparison
and boundary-specific review triggers
```

Current best observed result:

```text
100% misclassification catch rate
0 missed misclassified patients
30% overall auto-accept rate
100% safe auto-accept rate
```

Further validation across all patients is required before deployment or clinical use.

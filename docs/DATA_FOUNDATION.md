# LOT data-preparation and evaluation foundation

This foundation converts the two adjudication workbooks into normalized,
patient-level JSON Lines records. Raw patient identifiers are used only in memory;
derived records and fold manifests contain deterministic source-scoped hashes.
These hashes are pseudonyms, not a claim of irreversible anonymization.

## Outputs

All generated outputs live under `artifacts/` and are intentionally git-ignored
because they are derived from sensitive patient data.

| Output | Purpose |
|---|---|
| `artifacts/normalized/cota_patients.jsonl` | 136 COTA patients, ordered normalized treatment sets, dates, reviewer consensus, and separate COTA LOT |
| `artifacts/normalized/preamble_patients.jsonl` | New-workbook patients in the same normalized representation |
| `artifacts/reports/data_quality_report.json` | Aggregate continuation, malformed-treatment, duplicate, overlap, and fold counts |
| `artifacts/cv/cota_5fold_manifest.csv` | Deterministic patient-level folds for all 136 adjudicated COTA patients |
| `artifacts/reports/baseline_evaluation.json` | Current algorithm metrics and pseudonymous agree-but-wrong cases |

The normalized patient and fold-row contracts are defined by
`schemas/lot_patient.schema.json` and `schemas/cv_manifest.schema.json`.
The manifest's `fold` is the held-out/test fold for that patient; the other four
folds form the corresponding training set. Assignment is deterministic and
stratified by reviewer-consensus LOT.

## Reproduce

From the repository root:

```bash
python3 -m pip install -r requirements.txt
python3 src/py/prepare_lot_data.py
python3 src/py/evaluate_baseline.py
python3 -m unittest discover -s tests -v
```

Preparation and evaluation outputs are deterministically serialized. Re-running
the commands against unchanged workbooks produces byte-identical files.

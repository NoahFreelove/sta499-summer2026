# LOT data-preparation and evaluation foundation

This foundation reads the adjudicated COTA sheet in `LoT Adjudication
Datasets.xlsx` and the separate, unadjudicated `new_cota_data.xlsx`. It does not
contain model calls, retrieval, APIs, continuous learning, Docker, or a user
interface.

Raw patient identifiers are used only in memory. Every generated patient-level
artifact uses deterministic, source-scoped pseudonyms and lives under the
git-ignored `artifacts/restricted/` or `artifacts/blind/` directories.
These keys are operational pseudonyms, not a claim of irreversible anonymization.

## Data contracts

Two deliberately separate record contracts prevent answer leakage:

| Contract | Schema | Contents |
|---|---|---|
| Restricted evaluation record | `schemas/evaluation_record.schema.json` | Pseudonymous patient key, flattened normalized trajectory, source line/group metadata, vendor line dates, COTA LOT, reviewer labels when present, and quality metadata |
| Blind model input | `schemas/blind_model_input.schema.json` | Pseudonymous case key, flattened normalized trajectory, and diagnosis date only |

Blind inputs never contain source line numbers, COTA LOT, reviewer labels,
algorithm output/flags, quality flags, or vendor-derived line start/end dates.

## Vendor-boundary policy

Each reconstructed COTA treatment string may contain several bracketed regimen
groups. Preparation emits each group as a separate chronological treatment event
in source order. This prevents a future model from recovering COTA LOT merely by
counting vendor rows.

The source does not give separate dates for groups within a vendor line. Its only
treatment dates are vendor-line start/end dates, so those dates are retained only
in restricted evaluation records and omitted from blind inputs. The group order
is recoverable, but exact within-line group timing cannot be recovered from these
workbooks.

## Leakage-safe folds and overlap handling

Preparation computes two trajectory signatures from the flattened events:

- exact ordered trajectory;
- ordered trajectory after collapsing consecutive duplicate events.

Records sharing either signature form one exclusion group. Old/new records with
the same raw patient ID are also joined internally before the raw ID is discarded.
Connected exclusion groups are assigned as indivisible units to deterministic
folds. Therefore duplicate trajectories and cross-workbook versions cannot cross
train/test boundaries.

The patient-level exclusion manifest and CV manifest are restricted. The public
overlap report contains aggregate counts only.

## Generated outputs

| Output | Classification |
|---|---|
| `artifacts/restricted/evaluation/cota_adjudicated.jsonl` | Restricted |
| `artifacts/restricted/evaluation/cota_unadjudicated.jsonl` | Restricted |
| `artifacts/restricted/overlap/exclusion_groups.jsonl` | Restricted |
| `artifacts/restricted/cv/cota_adjudicated_5fold_manifest.csv` | Restricted |
| `artifacts/blind/cota_adjudicated.jsonl` | Blind model input |
| `artifacts/blind/cota_unadjudicated.jsonl` | Blind model input |
| `artifacts/public/data_quality_report.json` | Public aggregate |
| `artifacts/public/baseline_evaluation.json` | Public aggregate |

Patient-level baseline debugging rows are not written by default. To create a
clearly restricted debugging artifact, pass `--restricted-cases` to the evaluation
command.

## Reproduce from the repository root

```bash
python3 -m pip install -r requirements.txt
python3 src/py/prepare_lot_data.py
python3 src/py/evaluate_baseline.py
python3 -m unittest discover -s tests -v
```

Explicit paths are supported:

```bash
python3 src/py/prepare_lot_data.py \
  --project-root . \
  --adjudicated-workbook data/LoT\ Adjudication\ Datasets.xlsx \
  --new-workbook data/new_cota_data.xlsx \
  --output-dir artifacts

python3 src/py/evaluate_baseline.py \
  --project-root . \
  --patients artifacts/restricted/evaluation/cota_adjudicated.jsonl \
  --output artifacts/public/baseline_evaluation.json
```

`LOT_PROJECT_ROOT` may be used instead of `--project-root`. With unchanged inputs,
preparation and evaluation serialize byte-identical outputs.

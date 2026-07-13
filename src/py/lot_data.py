"""Data preparation utilities for adjudicated LOT workbooks.

This module deliberately keeps raw patient identifiers in memory only.  Public
records contain deterministic, source-scoped pseudonymous keys.
"""

from __future__ import annotations

import csv
import hashlib
import json
import math
import re
from collections import Counter, OrderedDict, defaultdict
from dataclasses import dataclass, field
from datetime import date, datetime
from pathlib import Path
from typing import Any, Iterable

from openpyxl import load_workbook

from textbook_algo_cota import SCT_TOKENS, normalize_drug


SCHEMA_VERSION = "1.0.0"
KEY_NAMESPACE = "lot-data-foundation-v1"


def is_blank(value: Any) -> bool:
    return value is None or (isinstance(value, str) and not value.strip())


def as_int(value: Any) -> int | None:
    if is_blank(value):
        return None
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return int(number) if math.isfinite(number) and number.is_integer() else None


def iso_date(value: Any) -> str | None:
    if is_blank(value):
        return None
    if isinstance(value, datetime):
        return value.date().isoformat()
    if isinstance(value, date):
        return value.isoformat()
    text = str(value).strip()
    return text or None


def patient_key(source: str, raw_patient_id: Any) -> str:
    """Return a deterministic source-scoped key without exposing the raw ID."""
    payload = f"{KEY_NAMESPACE}|{source}|{str(raw_patient_id).strip()}"
    digest = hashlib.sha256(payload.encode("utf-8")).hexdigest()[:16]
    return f"{source.upper()}-{digest}"


def parse_treatment(text: Any, *, bracketed: bool) -> tuple[frozenset[str], list[str]]:
    """Normalize a treatment string into a drug set and parser quality flags."""
    if is_blank(text):
        return frozenset(), ["missing_treatment"]

    value = str(text).strip()
    flags: list[str] = []
    if value.count("[") != value.count("]"):
        flags.append("unbalanced_brackets")

    if bracketed and "[" in value and "unbalanced_brackets" not in flags:
        groups = re.findall(r"\[([^\]]*)\]", value)
        if not groups:
            groups = [re.sub(r"[\[\]]", "", value)]
    else:
        groups = [re.sub(r"[\[\]]", "", value)]

    drugs: set[str] = set()
    for group in groups:
        for token in group.split(","):
            drug = normalize_drug(token)
            if drug and drug not in SCT_TOKENS:
                drugs.add(drug)

    if not drugs:
        flags.append("empty_drug_set")
    return frozenset(drugs), sorted(set(flags))


@dataclass
class TreatmentEvent:
    order: int
    source_line_number: int | None
    drugs: frozenset[str]
    start_date: str | None = None
    end_date: str | None = None
    line_start_date: str | None = None
    line_end_date: str | None = None
    fragment_count: int = 1
    continuation_row_count: int = 0
    quality_flags: list[str] = field(default_factory=list)

    def public_dict(self) -> dict[str, Any]:
        return {
            "order": self.order,
            "source_line_number": self.source_line_number,
            "drugs": sorted(self.drugs),
            "start_date": self.start_date,
            "end_date": self.end_date,
            "line_start_date": self.line_start_date,
            "line_end_date": self.line_end_date,
            "fragment_count": self.fragment_count,
            "continuation_row_count": self.continuation_row_count,
            "quality_flags": sorted(set(self.quality_flags)),
        }


@dataclass
class PatientRecord:
    source: str
    raw_patient_id: str
    trajectory: list[TreatmentEvent]
    reviewer_a_lot: int | None
    reviewer_b_lot: int | None
    reviewer_consensus_lot: int | None
    reviewer_consensus_status: str
    cota_lot: int | None = None
    preamble_lot: int | None = None
    quality_flags: list[str] = field(default_factory=list)

    @property
    def key(self) -> str:
        return patient_key(self.source, self.raw_patient_id)

    def public_dict(self) -> dict[str, Any]:
        return {
            "schema_version": SCHEMA_VERSION,
            "patient_key": self.key,
            "source": self.source,
            "trajectory": [event.public_dict() for event in self.trajectory],
            "reviewer_lot": {
                "reviewer_a": self.reviewer_a_lot,
                "reviewer_b": self.reviewer_b_lot,
                "consensus": self.reviewer_consensus_lot,
                "consensus_status": self.reviewer_consensus_status,
            },
            "cota_lot": self.cota_lot,
            "preamble_lot": self.preamble_lot,
            "quality_flags": sorted(set(self.quality_flags)),
        }


@dataclass
class LoadStats:
    source: str
    raw_rows: int = 0
    patient_count: int = 0
    continuation_rows: int = 0
    orphan_continuation_rows: int = 0
    rows_without_patient_or_treatment: int = 0
    treatment_events: int = 0

    def public_dict(self) -> dict[str, Any]:
        return dict(vars(self))


def consensus(a: int | None, b: int | None) -> tuple[int | None, str]:
    if a is None or b is None:
        return None, "missing_reviewer"
    if a != b:
        return None, "disagreement"
    return a, "agreed"


def _headers(worksheet: Any) -> tuple[dict[str, int], Iterable[tuple[Any, ...]]]:
    rows = worksheet.iter_rows(values_only=True)
    header = next(rows)
    return {str(value): index for index, value in enumerate(header) if value is not None}, rows


def _first(values: Iterable[Any]) -> Any:
    return next((value for value in values if not is_blank(value)), None)


def load_cota(path: str | Path, sheet: str = "Cota") -> tuple[list[PatientRecord], LoadStats]:
    """Load the COTA sheet, reconstructing blank-ID continuation rows."""
    worksheet = load_workbook(path, read_only=True, data_only=True)[sheet]
    columns, rows = _headers(worksheet)
    required = {
        "cpid", "line_of_therapy_c", "line_of_therapy_name",
        "date_start_line_of_therapy", "date_end_line_of_therapy",
        "Alpesh 1 LoT", "Alberto LOT ",
    }
    missing = required - columns.keys()
    if missing:
        raise ValueError(f"COTA sheet missing columns: {sorted(missing)}")

    stats = LoadStats(source="cota")
    grouped: OrderedDict[tuple[str, int | None], dict[str, Any]] = OrderedDict()
    reviewer_values: dict[str, dict[str, list[int]]] = defaultdict(lambda: {"a": [], "b": []})
    current_patient: str | None = None
    current_line: int | None = None

    for row_number, row in enumerate(rows, start=2):
        stats.raw_rows += 1
        raw_id = row[columns["cpid"]]
        treatment = row[columns["line_of_therapy_name"]]
        continuation = is_blank(raw_id)
        if continuation:
            if is_blank(treatment):
                stats.rows_without_patient_or_treatment += 1
                continue
            stats.continuation_rows += 1
            if current_patient is None:
                stats.orphan_continuation_rows += 1
                continue
            raw_id = current_patient
            line_number = current_line
        else:
            current_patient = str(raw_id).strip()
            line_number = as_int(row[columns["line_of_therapy_c"]])
            current_line = line_number

        raw_id_text = str(raw_id).strip()
        a = as_int(row[columns["Alpesh 1 LoT"]])
        b = as_int(row[columns["Alberto LOT "]])
        if a is not None:
            reviewer_values[raw_id_text]["a"].append(a)
        if b is not None:
            reviewer_values[raw_id_text]["b"].append(b)

        key = (raw_id_text, line_number)
        group = grouped.setdefault(key, {
            "first_row": row_number, "parts": [], "continuations": 0,
            "start_dates": [], "end_dates": [],
        })
        if not is_blank(treatment):
            group["parts"].append(str(treatment))
        if continuation:
            group["continuations"] += 1
        group["start_dates"].append(row[columns["date_start_line_of_therapy"]])
        group["end_dates"].append(row[columns["date_end_line_of_therapy"]])

    by_patient: OrderedDict[str, list[tuple[int | None, dict[str, Any]]]] = OrderedDict()
    for (raw_id, line_number), group in grouped.items():
        by_patient.setdefault(raw_id, []).append((line_number, group))

    records: list[PatientRecord] = []
    for raw_id, entries in by_patient.items():
        entries.sort(key=lambda item: (
            item[0] is None,
            item[0] if item[0] is not None else item[1]["first_row"],
            item[1]["first_row"],
        ))
        trajectory: list[TreatmentEvent] = []
        patient_flags: list[str] = []
        for order, (line_number, group) in enumerate(entries, start=1):
            full_treatment = "".join(group["parts"])
            drugs, flags = parse_treatment(full_treatment, bracketed=True)
            patient_flags.extend(flags)
            trajectory.append(TreatmentEvent(
                order=order,
                source_line_number=line_number,
                drugs=drugs,
                start_date=iso_date(_first(group["start_dates"])),
                end_date=iso_date(_first(group["end_dates"])),
                fragment_count=len(group["parts"]),
                continuation_row_count=group["continuations"],
                quality_flags=flags,
            ))

        values = reviewer_values[raw_id]
        a = max(values["a"], default=None)
        b = max(values["b"], default=None)
        agreed, status = consensus(a, b)
        cota_values = [line for line, _ in entries if line is not None]
        records.append(PatientRecord(
            source="cota",
            raw_patient_id=raw_id,
            trajectory=trajectory,
            reviewer_a_lot=a,
            reviewer_b_lot=b,
            reviewer_consensus_lot=agreed,
            reviewer_consensus_status=status,
            cota_lot=max(cota_values, default=None),
            quality_flags=patient_flags,
        ))

    stats.patient_count = len(records)
    stats.treatment_events = sum(len(record.trajectory) for record in records)
    return records, stats


def load_preamble(
    path: str | Path, sheet: str = "Preamble_treatment"
) -> tuple[list[PatientRecord], LoadStats]:
    """Load the newer PREAMBLE adjudication workbook at treatment-row granularity."""
    worksheet = load_workbook(path, read_only=True, data_only=True)[sheet]
    columns, rows = _headers(worksheet)
    required = {
        "patid", "Treatment_Regimens", "Starting_date", "Ending_date",
        "LINE_START_DATE", "LINE_END_DATE", "JK_line", "AR_Line",
        "Line_number_PREAMBLE",
    }
    missing = required - columns.keys()
    if missing:
        raise ValueError(f"PREAMBLE sheet missing columns: {sorted(missing)}")

    stats = LoadStats(source="preamble")
    events: OrderedDict[str, list[dict[str, Any]]] = OrderedDict()
    reviewer_values: dict[str, dict[str, list[int]]] = defaultdict(lambda: {"a": [], "b": []})
    vendor_values: dict[str, list[int]] = defaultdict(list)

    for row_number, row in enumerate(rows, start=2):
        stats.raw_rows += 1
        raw_id = row[columns["patid"]]
        treatment = row[columns["Treatment_Regimens"]]
        if is_blank(raw_id):
            if not is_blank(treatment):
                stats.orphan_continuation_rows += 1
            else:
                stats.rows_without_patient_or_treatment += 1
            continue
        raw_id_text = str(raw_id).strip()
        events.setdefault(raw_id_text, []).append({
            "row": row_number,
            "treatment": treatment,
            "start": row[columns["Starting_date"]],
            "end": row[columns["Ending_date"]],
            "line_start": row[columns["LINE_START_DATE"]],
            "line_end": row[columns["LINE_END_DATE"]],
            "vendor_line": as_int(row[columns["Line_number_PREAMBLE"]]),
        })
        a = as_int(row[columns["JK_line"]])
        b = as_int(row[columns["AR_Line"]])
        vendor = as_int(row[columns["Line_number_PREAMBLE"]])
        if a is not None:
            reviewer_values[raw_id_text]["a"].append(a)
        if b is not None:
            reviewer_values[raw_id_text]["b"].append(b)
        if vendor is not None:
            vendor_values[raw_id_text].append(vendor)

    records: list[PatientRecord] = []
    for raw_id, raw_events in events.items():
        trajectory: list[TreatmentEvent] = []
        patient_flags: list[str] = []
        for order, item in enumerate(raw_events, start=1):
            drugs, flags = parse_treatment(item["treatment"], bracketed=False)
            patient_flags.extend(flags)
            trajectory.append(TreatmentEvent(
                order=order,
                source_line_number=item["vendor_line"],
                drugs=drugs,
                start_date=iso_date(item["start"]),
                end_date=iso_date(item["end"]),
                line_start_date=iso_date(item["line_start"]),
                line_end_date=iso_date(item["line_end"]),
                quality_flags=flags,
            ))
        values = reviewer_values[raw_id]
        a = max(values["a"], default=None)
        b = max(values["b"], default=None)
        agreed, status = consensus(a, b)
        records.append(PatientRecord(
            source="preamble",
            raw_patient_id=raw_id,
            trajectory=trajectory,
            reviewer_a_lot=a,
            reviewer_b_lot=b,
            reviewer_consensus_lot=agreed,
            reviewer_consensus_status=status,
            preamble_lot=max(vendor_values[raw_id], default=None),
            quality_flags=patient_flags,
        ))

    stats.patient_count = len(records)
    stats.treatment_events = sum(len(record.trajectory) for record in records)
    return records, stats


def trajectory_signature(record: PatientRecord, *, collapse_consecutive: bool = False) -> str:
    sequence: list[tuple[str, ...]] = []
    for event in record.trajectory:
        drugs = tuple(sorted(event.drugs))
        if not drugs:
            continue
        if collapse_consecutive and sequence and drugs == sequence[-1]:
            continue
        sequence.append(drugs)
    return json.dumps(sequence, separators=(",", ":"), ensure_ascii=True)


def duplicate_summary(records: list[PatientRecord]) -> dict[str, int]:
    signatures: dict[str, list[PatientRecord]] = defaultdict(list)
    collapsed_signatures: dict[str, list[PatientRecord]] = defaultdict(list)
    for record in records:
        signatures[trajectory_signature(record)].append(record)
        collapsed_signatures[trajectory_signature(record, collapse_consecutive=True)].append(record)
    groups = [group for group in signatures.values() if len(group) > 1]
    collapsed_groups = [group for group in collapsed_signatures.values() if len(group) > 1]
    consecutive_duplicate_events = 0
    patients_with_consecutive_duplicates = 0
    for record in records:
        duplicate_count = sum(
            bool(current.drugs) and current.drugs == previous.drugs
            for previous, current in zip(record.trajectory, record.trajectory[1:])
        )
        consecutive_duplicate_events += duplicate_count
        patients_with_consecutive_duplicates += duplicate_count > 0
    return {
        "exact_duplicate_trajectory_groups": len(groups),
        "patients_in_exact_duplicate_trajectory_groups": sum(len(group) for group in groups),
        "duplicate_groups_after_collapsing_consecutive_repeats": len(collapsed_groups),
        "patients_in_duplicate_groups_after_collapsing_consecutive_repeats": sum(
            len(group) for group in collapsed_groups
        ),
        "patients_with_consecutive_duplicate_treatments": patients_with_consecutive_duplicates,
        "consecutive_duplicate_treatment_events": consecutive_duplicate_events,
    }


def malformed_summary(records: list[PatientRecord]) -> dict[str, Any]:
    event_counts: Counter[str] = Counter()
    patient_counts: Counter[str] = Counter()
    malformed_events = 0
    for record in records:
        seen: set[str] = set()
        for event in record.trajectory:
            event_counts.update(event.quality_flags)
            seen.update(event.quality_flags)
            malformed_events += bool(event.quality_flags)
        patient_counts.update(seen)
    return {
        "event_counts_by_flag": dict(sorted(event_counts.items())),
        "patient_counts_by_flag": dict(sorted(patient_counts.items())),
        "malformed_event_count": malformed_events,
        "patients_with_malformed_treatments": len({
            record.key for record in records
            if any(event.quality_flags for event in record.trajectory)
        }),
    }


def overlap_summary(old: list[PatientRecord], new: list[PatientRecord]) -> dict[str, int]:
    old_ids = {record.raw_patient_id for record in old}
    new_ids = {record.raw_patient_id for record in new}
    old_exact = Counter(trajectory_signature(record) for record in old)
    new_exact = Counter(trajectory_signature(record) for record in new)
    old_collapsed = Counter(trajectory_signature(record, collapse_consecutive=True) for record in old)
    new_collapsed = Counter(trajectory_signature(record, collapse_consecutive=True) for record in new)
    old_regimens = {tuple(sorted(event.drugs)) for record in old for event in record.trajectory if event.drugs}
    new_regimens = {tuple(sorted(event.drugs)) for record in new for event in record.trajectory if event.drugs}
    return {
        "raw_patient_id_overlap_count": len(old_ids & new_ids),
        "exact_trajectory_signatures_shared": len(old_exact.keys() & new_exact.keys()),
        "exact_trajectory_patient_pair_count": sum(
            old_exact[sig] * new_exact[sig] for sig in old_exact.keys() & new_exact.keys()
        ),
        "collapsed_trajectory_signatures_shared": len(old_collapsed.keys() & new_collapsed.keys()),
        "collapsed_trajectory_patient_pair_count": sum(
            old_collapsed[sig] * new_collapsed[sig]
            for sig in old_collapsed.keys() & new_collapsed.keys()
        ),
        "unique_normalized_regimen_sets_old": len(old_regimens),
        "unique_normalized_regimen_sets_new": len(new_regimens),
        "unique_normalized_regimen_sets_shared": len(old_regimens & new_regimens),
    }


def make_folds(records: list[PatientRecord], n_folds: int = 5) -> list[dict[str, Any]]:
    """Assign patients to stable folds stratified by consensus LOT."""
    if n_folds < 2:
        raise ValueError("n_folds must be at least 2")
    strata: dict[int | None, list[PatientRecord]] = defaultdict(list)
    for record in records:
        strata[record.reviewer_consensus_lot].append(record)

    fold_sizes = [0] * n_folds
    manifest: list[dict[str, Any]] = []
    for lot in sorted(strata, key=lambda value: (value is None, value or 0)):
        ordered = sorted(
            strata[lot],
            key=lambda record: hashlib.sha256(
                f"lot-cv-v1|{record.key}".encode("utf-8")
            ).hexdigest(),
        )
        start_fold = min(range(n_folds), key=lambda fold: (fold_sizes[fold], fold))
        for index, record in enumerate(ordered):
            fold = (start_fold + index) % n_folds
            fold_sizes[fold] += 1
            manifest.append({
                "patient_key": record.key,
                "fold": fold,
                "reviewer_consensus_lot": record.reviewer_consensus_lot,
            })
    return manifest


def write_jsonl(path: str | Path, records: list[PatientRecord]) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for record in sorted(records, key=lambda item: item.key):
            handle.write(json.dumps(record.public_dict(), sort_keys=True, separators=(",", ":")))
            handle.write("\n")


def write_json(path: str | Path, value: Any) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        json.dump(value, handle, indent=2, sort_keys=True)
        handle.write("\n")


def write_fold_manifest(path: str | Path, rows: list[dict[str, Any]]) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["patient_key", "fold", "reviewer_consensus_lot"],
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(sorted(rows, key=lambda row: row["patient_key"]))


def read_public_jsonl(path: str | Path) -> list[dict[str, Any]]:
    with Path(path).open(encoding="utf-8") as handle:
        return [json.loads(line) for line in handle if line.strip()]

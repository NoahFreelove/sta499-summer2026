from __future__ import annotations

from typing import Literal

from pydantic import BaseModel, Field, model_validator


class CorrectedVendorLine(BaseModel):
    vendor_lot: int = Field(ge=1)
    adjudicated_lot: int = Field(ge=1)
    action: Literal["keep", "merge_with_previous", "uncertain"]
    normalized_regimen: list[str] = Field(default_factory=list)
    rationale: str
    applied_rules: list[str] = Field(default_factory=list)
    supporting_reference_ids: list[str] = Field(default_factory=list)
    confidence: float = Field(ge=0.0, le=1.0)


class AdjudicationResult(BaseModel):
    patient_id: str
    corrected_lines: list[CorrectedVendorLine]
    summary_rationale: str
    overall_confidence: float = Field(ge=0.0, le=1.0)
    uncertainty_flags: list[str] = Field(default_factory=list)
    requires_human_review: bool = False

    @model_validator(mode="after")
    def validate_mapping(self):
        vendor_lots = [line.vendor_lot for line in self.corrected_lines]
        if vendor_lots != sorted(vendor_lots) or len(vendor_lots) != len(set(vendor_lots)):
            raise ValueError("vendor_lot values must be unique and sorted.")

        corrected = [line.adjudicated_lot for line in self.corrected_lines]
        if corrected:
            if corrected[0] != 1:
                raise ValueError("The first adjudicated LoT must be 1.")
            for previous, current in zip(corrected, corrected[1:]):
                if current not in {previous, previous + 1}:
                    raise ValueError(
                        "adjudicated_lot may stay the same for a merge or increase by exactly one."
                    )
        return self


class Disagreement(BaseModel):
    category: Literal[
        "patient_id",
        "missing_vendor_line",
        "mapping",
        "action",
        "confidence",
        "other",
    ]
    vendor_lot: int | None = None
    description: str
    worker_value: str | None = None
    auditor_value: str | None = None
    severity: Literal["minor", "major", "critical"]


class AuditCritique(BaseModel):
    patient_id: str
    worker_valid: bool
    auditor_valid: bool
    disagreements: list[Disagreement] = Field(default_factory=list)
    preferred_result: Literal["worker", "auditor", "neither", "equivalent"]
    critique: str
    human_review_required: bool


class PipelineDecision(BaseModel):
    patient_id: str
    status: Literal[
        "accepted_agreement",
        "accepted_worker",
        "accepted_auditor",
        "human_review",
        "processing_error",
    ]
    selected_result: Literal["worker", "auditor", "none"]
    reason: str
    worker_confidence: float | None = None
    auditor_confidence: float | None = None
    major_disagreement: bool = False

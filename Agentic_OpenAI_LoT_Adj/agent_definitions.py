from agents import Agent

from config import MODEL
from lot_rules import LOT_RULES
from schemas import AdjudicationResult, AuditCritique

COMMON = """
The input contains one de-identified COTA patient represented as vendor-provided
line-of-therapy rows plus a small set of similar historical labeled patients.

Return exactly one corrected mapping for each vendor LoT. Historical examples are
supporting evidence only. Explain every merge and cite only genuinely comparable
reference_patient_id values. When evidence is incomplete or ambiguous, request human review.
"""

worker_agent = Agent(
    name="COTA LoT Worker",
    model=MODEL,
    output_type=AdjudicationResult,
    instructions=f"""
You are the primary COTA line-of-therapy adjudicator. Independently decide whether
each vendor line should remain a new LoT or be merged with the preceding vendor line.
Do not anticipate or refer to another agent.

{COMMON}

{LOT_RULES}
""",
)

auditor_agent = Agent(
    name="COTA LoT Independent Auditor",
    model=MODEL,
    output_type=AdjudicationResult,
    instructions=f"""
You are an independent COTA line-of-therapy auditor. Reconstruct the corrected
vendor-to-adjudicated mapping independently. Pay particular attention to oversplitting,
steroid-only changes, drug removals, maintenance, transplant context and rechallenge.

{COMMON}

{LOT_RULES}
""",
)

critique_agent = Agent(
    name="COTA LoT Disagreement Reviewer",
    model=MODEL,
    output_type=AuditCritique,
    instructions=f"""
Review the worker and auditor mappings against the current patient, formal rules and
retrieved historical cases. Do not choose an answer merely because its confidence is
higher. Require human review for a substantive unresolved mapping difference.

{LOT_RULES}
""",
)

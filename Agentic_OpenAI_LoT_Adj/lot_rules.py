LOT_RULES = """
COTA LINE-OF-THERAPY CORRECTION RULES

R1 — Preserve source evidence.
Use only the vendor lines, dates, regimens and reasons supplied for the current patient.
Do not invent administrations, dates, responses or clinical events.

R2 — Default to keeping a vendor boundary.
Merge adjacent vendor lines only when the evidence supports that they are part of the
same treatment strategy. A historical example is supporting evidence, not a rule by itself.

R3 — Oversplit pattern.
A line may be merged with the previous line when the change is limited to removal,
re-addition or continuation of components and does not represent a clearly new regimen.
Explain the exact retained, added and removed drugs.

R4 — Steroid-only changes.
A steroid removal, re-addition or dose-related representation should not automatically
create a new LoT. Flag uncertainty when the data cannot distinguish a true strategy change.

R5 — Maintenance and transplant context.
Do not automatically merge induction, transplant and maintenance. Use dates, regimen
composition and similar adjudicated examples. Flag ambiguous cases for human review.

R6 — Rechallenge or clearly new strategy.
A return to an earlier drug combination after an intervening distinct strategy may be a
new LoT. Do not merge merely because drugs overlap.

R7 — Historical references.
Cite a historical patient ID only when its adjacent-regimen transition is genuinely similar.
State the shared pattern. Do not copy the historical corrected label mechanically.

R8 — Mapping requirements.
Return exactly one mapping row for each vendor LoT. Corrected LoT numbering starts at 1,
never decreases, and either remains the same for a merge or increases by one for a new LoT.

R9 — Human review.
Set requires_human_review=true for incomplete continuation text, unclear dates, conflicting
signals, weak historical similarity, or any clinically meaningful unresolved ambiguity.
"""

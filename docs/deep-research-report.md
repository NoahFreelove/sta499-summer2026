# Multiple Myeloma Line of Therapy Classification Knowledge Base

## Executive synthesis

This report is for **historical treatment-line classification** in multiple myeloma, not treatment recommendation. The literature shows that oncology real-world evidence does **not** have a single universal line-of-therapy definition. Instead, most algorithms require explicit conventions for: a regimen-defining window, rules for drug addition or replacement, rules for treatment gaps, and disease-specific exceptions. Those conventions vary materially across sources, especially for gap thresholds and partial regimen changes. citeturn29search1turn34search0turn30view0turn41search13

The strongest **multiple-myeloma-specific** finding is that frontline therapy for transplant-eligible patients is commonly conceptualized as a **single treatment line** that spans **induction → autologous stem-cell transplant → post-transplant maintenance**, with transplant and high-dose melphalan functioning as part of the frontline treatment workflow rather than independent LOTs. This view appears in myeloma-specific educational guidance, transplant guidance, and real-world MM studies that explicitly analyze stem-cell transplant as occurring **within first-line therapy**. citeturn27view0turn37search1turn19search0turn43search0turn26search5

For deterministic classification, the best-supported hard rules are these: a regimen starts on the first active anti-myeloma treatment date; drugs started within a short regimen-build window belong to the same regimen; dose reduction or route change alone does not create a new line; temporary holds shorter than the chosen discontinuation gap do not create a new line; full non-overlap switches are new lines; CAR-T lymphodepletion and CAR-T infusion belong to the same cellular-therapy line; and bispecific step-up dosing belongs to the same bispecific line. citeturn28search9turn41search7turn29search1turn21search0turn21search1turn22search1turn48search0turn48search2turn48search6

The weakest or most conflict-prone areas are: exact gap length for discontinuation and restart; whether a partial drug drop stays in the same line outside a clear maintenance or toxicity context; whether bridging or holding therapy before CAR-T should count as its own line; and whether steroid-only treatment represents a true anti-myeloma line versus temporary disease-control or palliative bridging. Those scenarios should default to **human review** unless the project adopts an explicit convention. citeturn34search0turn41search13turn47search4turn40search9

### Evidence hierarchy and conflict map

At the top of the hierarchy are formal oncology LOT methodology sources, which consistently describe LOT classification as a rules engine built from regimen definition, drug-change logic, and gap logic, but they also show that thresholds are **conventions**, not natural laws. Hess and colleagues emphasize drug additions, discontinuations, exchanges, and therapy gaps as the core components of a LOT algorithm, while the 2026 scoping review shows that published algorithms use common gap definitions of **60, 90, 120, and 180 days**, depending on context. citeturn30view0turn29search1turn34search0

Myeloma-specific guidance adds disease-specific exceptions that general oncology algorithms would miss. In particular, myeloma sources treat induction, transplant, and maintenance as a linked therapeutic arc in frontline care, and transplant-focused guidelines continue to describe ASCT as **consolidation** after induction, with maintenance after recovery. That matters because a generic “drug-set changed” rule will often misclassify post-transplant maintenance as a new line unless an MM exception is added. citeturn27view0turn37search1turn19search0turn19search3turn32search14

The literature also shows that some clinically validated general-oncology algorithms treat stopping one or more drugs, dose changes, route changes, and later restart of the same drug as **not** advancing the LOT, whereas other frameworks more aggressively advance the line after any substantive regimen change outside the build window. That disagreement is tolerable in solid tumors with simpler episodes, but it becomes a major source of misclassification in myeloma because frontline therapy is intentionally **multipart** and because maintenance and cellular-therapy workflows are role-driven rather than composition-driven. citeturn41search7turn29search1turn36search3

### What should be hard-coded

The most defensible hard rules are those tied either to explicit product labels or to stable MM care structure. CAR-T labeling clearly specifies a preplanned lymphodepletion course with fludarabine plus cyclophosphamide immediately preceding ide-cel or cilta-cel administration, which supports grouping lymphodepletion and CAR-T infusion into one cellular-therapy line. Likewise, teclistamab, elranatamab, and talquetamab labels all describe **step-up dosing** as part of treatment initiation, not as separate therapy lines. citeturn21search0turn21search1turn22search1turn48search0turn48search2turn48search6

The other strong hard rule is the MM transplant workflow exception: if a patient follows the standard frontline transplant path, the classifier should keep induction, mobilization/collection, conditioning high-dose melphalan, ASCT, consolidation if present, and post-transplant maintenance within the same LOT unless there is documented relapse/progression or a clearly non-planned salvage switch. citeturn27view0turn37search1turn19search0turn20search4turn20search5

### What should stay review-based

A deterministic engine should **not** pretend that the literature resolves everything. Steroid-only episodes require review because corticosteroids can be used alone either as temporary disease control, symptom control, or as true antimyeloma treatment when a patient is too unwell for other therapy. Bridging and holding therapy before CAR-T also require review because recent consensus documents define them by timing relative to leukapheresis and lymphodepletion, but they do not establish a universal prior-line counting rule for databases. Partial substitutions—especially toxicity-driven ones—also remain context-sensitive in myeloma. citeturn40search9turn47search4turn25search3turn41search13

## Boundary logic from the literature

### Core oncology LOT mechanics

Formal LOT methodology papers converge on a common structure. First, define the regimen; second, decide what changes advance the line; third, decide what length of treatment gap ends a line. Hess’s framework and subsequent workshop material use a short initial window—commonly **28 days**, sometimes **35 days** in other publications—to assemble the regimen. After that window, adding a new active drug usually advances the line, while a sufficiently long therapy gap followed by restart also advances the line. The 2026 scoping review confirms that published algorithms differ substantially on the discontinuation interval, most commonly using 60, 90, 120, or 180 days. citeturn28search9turn29search1turn34search0turn30view0

That general architecture is appropriate for multiple myeloma, but it must be modified for MM-specific workflows. Myeloma treatment often contains phase transitions that are **planned**, not signs of failure: induction to transplant, transplant to consolidation, consolidation to maintenance, or CAR-T lymphodepletion to infusion. A pure composition-based classifier will overcount LOTs unless it recognizes role and temporal adjacency. The Ailawadhi validation study makes this explicit by noting that MM treatment complexity required identification of smaller treatment components or “mini-LOTs” before assembling the overall line, and the maintenance algorithm performed less well than simple LOT1 identification—evidence that phase handling is a real source of error. citeturn36search3turn35search17

**Practical consequence:** use a **28-day regimen-build window** as a project convention, but label it as a convention rather than a universal truth. Keep a configurable discontinuation gap parameter because the literature does not justify a single immutable number across all use cases. For MM, if the project needs one default, 90 days is the most common general-oncology convention among the reviewed formal methods, but that still remains a convention rather than an MM-validated gold standard. citeturn29search1turn34search0turn30view0

### Progression, relapse, gaps, and rechallenge

IMWG-based relapse language separates **biochemical relapse**, **clinical relapse**, and **refractory disease**. The IMF’s IMWG-aligned relapse page describes relapse as reappearance of signs and symptoms after improvement, notes biochemical relapse as progression by M-protein change without organ-damage symptoms, and defines refractory disease as progression during treatment or within **60 days** of treatment. citeturn49view0

Those clinical categories matter, but they do **not** by themselves solve database line counting. In real-world data, progression is often missing or poorly structured, so subsequent therapy initiation is commonly used as a proxy for progression. Some conceptual frameworks also note that objective progression does not always immediately equal treatment failure or immediate switching. Therefore, for a treatment-record classifier, **documented progression plus treatment change** is strong evidence for a new LOT, whereas **documented progression without treatment change** should usually be a review flag instead of an automatic new-line event. citeturn38search0turn38search9turn38search16turn49view0

For **rechallenge**, the formal LOT literature supports treating restart of the **same regimen after a long gap** as a new line under the chosen gap rule. However, a short planned interruption with later restart of the same drugs belongs to the same line, and MM adds a special case: collection/recovery around transplant can create apparent gaps that are procedural, not true treatment discontinuation. The “same regimen restarted after the break is still the same line” view appears in expert interviews, but only when the break is interpreted as a continuation of the same therapeutic episode rather than a completed line followed by later retreatment. citeturn29search1turn41search13turn27view0

**Practical consequence:**  
A restart should be classified as:

- **SAME_LINE** if the interruption is below the gap threshold or is clearly embedded in a transplant or step-up workflow.
- **NEW_LINE** if the interruption exceeds the project’s discontinuation threshold and the same regimen is reinitiated.
- **HUMAN_REVIEW** if dates are incomplete, or if progression/intent information is missing and the break could represent either planned interruption or true retreatment. citeturn29search1turn34search0turn41search13

### Multiple myeloma transplant and maintenance workflow

This is the most important MM-specific override. The IMF frontline page states that, for transplant-eligible patients, frontline therapy is a **three-step process that comprises one LOT**: induction, ASCT, and maintenance. The MMRF treatment page makes the same point in nearly identical terms, adding that induction is followed by ASCT and maintenance, and that these are collectively one line of therapy. Real-world MM studies also explicitly analyze **SCT in line 1** rather than as an independent line. citeturn27view0turn37search1turn43search0turn26search5

Transplant guidelines are consistent with that interpretation. The ASTCT consensus describes early autologous transplantation as **consolidation** after 4–6 cycles of induction, and transplant-focused reviews define consolidation in MM as limited-duration systemic therapy **after ASCT and before maintenance**. The front-line pathway is therefore a single multi-phase therapeutic strategy, not a chain of unrelated LOTs. citeturn19search0turn19search3turn32search14

High-dose melphalan should be interpreted as **conditioning**, not as a stand-alone salvage chemotherapy line, when it occurs in transplant adjacency. Multiple transplant sources describe high-dose melphalan—typically Mel200—as the standard conditioning regimen before ASCT in MM. In classification terms, that means high-dose melphalan immediately adjacent to ASCT should be normalized as a **conditioning/procedural token** and kept in the same LOT. citeturn20search4turn20search5turn20search15

Post-transplant maintenance should also remain in the same LOT when it follows planned frontline ASCT recovery. The IMF frontline page states that maintenance is recommended after recovery from ASCT and that medications used for induction may also be used for maintenance with changes in dose and schedule; ASTCT recommends lenalidomide maintenance after autologous transplantation; and MM reviews define maintenance as the therapy that follows consolidation or transplant recovery to preserve remission. citeturn27view0turn19search0turn19search9turn35search0

### Cellular therapy, bispecifics, and bridging

CAR-T workflows have unusually clear labeled sequencing. The ABECMA and CARVYKTI labels both specify lymphodepleting chemotherapy with **cyclophosphamide 300 mg/m² plus fludarabine 30 mg/m² daily for 3 days** before CAR-T administration, with the product label treating this as a single treatment course centered on the CAR-T infusion. The IMWG CAR-T consensus likewise addresses patient selection, bridging therapy, lymphodepletion, and response assessment as parts of CAR-T practice rather than separate diseasespecific LOT entities. citeturn21search0turn21search1turn25search0turn25search2turn25search3

Bispecific antibodies also have a clear internal initiation structure. The TECVAYLI, ELREXFIO, and TALVEY labels all define **step-up dosing** followed by the first full treatment dose, after which therapy continues until progression or unacceptable toxicity. That supports treating the step-up phase and first full dose as the same LOT rather than separate lines. citeturn22search1turn48search0turn48search2turn48search6

Bridging before CAR-T is structurally important but methodologically unresolved for LOT counting. The 2026 bridging consensus defines **bridging therapy** as therapy after leukapheresis and before lymphodepleting chemotherapy, and distinguishes it from **pre-bridging/holding therapy** before apheresis. The short-term goals are disease control and preservation of performance status until infusion, not necessarily establishment of a durable new line. That role-based intent argues for representing bridging as a **separate role taxonomy**, but defaulting to **human review** for LOT counting unless the project explicitly chooses to count or exclude bridge episodes. citeturn47search4turn25search3

Steroid-only bridging is especially ambiguous. Cancer.ca notes that dexamethasone or prednisone are often used with other anti-myeloma drugs, but can also be used **alone** when a person is too sick for treatment, to slow myeloma growth and reduce M-protein. Because the same drug pattern can reflect true anti-myeloma intent, temporary bridge, or palliative control, a steroid-only episode should not be deterministically classified without timing and intent context. citeturn40search9

## Machine-readable rule library

The YAML below is a **project-ready rule library**. It intentionally separates hard rules from review-based rules and maps every rule to source IDs listed later in the source manifest. The strongest rules are anchored to MM-specific workflow sources, IMWG/ASTCT guidance, and product labels; the weaker rules are explicitly marked as limited or conflicting. citeturn27view0turn19search0turn21search0turn22search1turn48search0turn47search4

```yaml
mm_lot_start_first_active_regimen:
  title: First active anti-myeloma treatment starts a line
  scope: multiple myeloma, all settings
  required_inputs:
    - canonical_drug_names
    - administration_dates
    - anti_myeloma_activity_flag
  applies_when: >
    The record contains the earliest identifiable active anti-myeloma systemic therapy
    for a patient and no prior active line has already been established.
  default_decision: NEW_LINE
  exceptions:
    - Do not start a line from mobilization, leukapheresis, stem-cell collection, or transplant procedure tokens alone.
    - Do not start a line from supportive-care-only records.
  reason_codes:
    - FIRST_ACTIVE_MM_THERAPY
    - EXCLUDE_PROCEDURE_ONLY
    - EXCLUDE_SUPPORTIVE_ONLY
  evidence_strength: strong
  source_ids:
    - SRC_HESS_2021
    - SRC_IMF_FRONTLINE
    - SRC_MMRF_TREATMENT
  known_failure_modes:
    - Earliest therapy missing from source data
    - Oral therapy start date captured late
    - Supportive dexamethasone miscoded as active treatment
  last_verified: 2026-07-13

mm_lot_regimen_build_up_28d:
  title: Regimen build-up within 28 days stays in the same line
  scope: multiple myeloma, systemic therapy construction
  required_inputs:
    - canonical_drug_names
    - administration_dates
    - current_line_start_date
  applies_when: >
    A new active anti-myeloma drug first appears on or before day 28 from the current
    line start date and there is no explicit documentation that the addition represents
    salvage escalation after failure.
  default_decision: SAME_LINE
  exceptions:
    - If explicit progression, relapse, or treatment-failure documentation states the addition is salvage escalation, send to HUMAN_REVIEW.
    - If the project intentionally uses a different regimen-building window, this rule must be parameterized rather than edited.
  reason_codes:
    - BUILD_UP_WITHIN_28D
    - WINDOW_PARAMETERIZED
    - POSSIBLE_SALVAGE_ESCALATION
  evidence_strength: moderate
  source_ids:
    - SRC_HESS_2021
    - SRC_JMIR_VALIDATION_2021
    - SRC_FLATIRON_LOT_2026
  known_failure_modes:
    - Late charting of oral agents
    - Missing explicit intent note
    - True salvage add-on occurring inside the window
  last_verified: 2026-07-13

mm_lot_add_new_agent_outside_build_up:
  title: New active drug added after the build-up window starts a new line
  scope: multiple myeloma, systemic therapy transitions
  required_inputs:
    - canonical_drug_names
    - administration_dates
    - current_line_start_date
    - regimen_role_context
  applies_when: >
    A new active anti-myeloma drug that was not part of the established regimen first
    appears after day 28 from line start and the event is not an MM-specific planned
    transition such as post-transplant maintenance, CAR-T lymphodepletion, or bispecific
    step-up dosing.
  default_decision: NEW_LINE
  exceptions:
    - If the change is an MM-specific planned workflow covered by another rule, apply that other rule first.
    - If explicit toxicity-driven temporary overlap is documented and no new long-term regimen is intended, send to HUMAN_REVIEW.
  reason_codes:
    - NEW_ACTIVE_DRUG_AFTER_WINDOW
    - EXCLUDE_MM_PLANNED_WORKFLOW
    - POSSIBLE_TOXICITY_OVERLAP
  evidence_strength: moderate
  source_ids:
    - SRC_HESS_2021
    - SRC_AILAWADHI_2024
    - SRC_ASTCT_2022
  known_failure_modes:
    - Maintenance misread as salvage
    - Missing role metadata
    - Investigational agent not recognized as active drug
  last_verified: 2026-07-13

mm_lot_full_nonoverlap_switch:
  title: Full switch to a non-overlapping active regimen is a new line
  scope: multiple myeloma, all settings
  required_inputs:
    - prior_active_drug_set
    - next_active_drug_set
    - administration_dates
  applies_when: >
    The subsequent active regimen shares no continuing core anti-myeloma drug with the
    prior regimen and is not merely conditioning, lymphodepletion, or another procedural
    step.
  default_decision: NEW_LINE
  exceptions:
    - If dates are too sparse to determine overlap, send to HUMAN_REVIEW.
  reason_codes:
    - FULL_SWITCH_NONOVERLAP
    - INSUFFICIENT_OVERLAP_DATA
  evidence_strength: strong
  source_ids:
    - SRC_HESS_2021
    - SRC_IMF_FRONTLINE
    - SRC_ASTCT_2022
  known_failure_modes:
    - Drug normalization failures
    - Missed supportive drugs treated as core drugs
  last_verified: 2026-07-13

mm_lot_explicit_de_escalation_to_maintenance:
  title: Explicit de-escalation to maintenance remains the same line
  scope: multiple myeloma, maintenance transitions
  required_inputs:
    - prior_regimen_role
    - next_regimen_role
    - transplant_context
    - canonical_drug_names
    - administration_dates
  applies_when: >
    The treatment record explicitly identifies a transition to maintenance, or the
    sequence is a standard post-transplant recovery-to-maintenance transition in MM.
  default_decision: SAME_LINE
  exceptions:
    - If maintenance intent is not explicit and transplant context is absent, send to HUMAN_REVIEW.
    - If relapse/progression is documented before maintenance start, classify with the progression rule instead.
  reason_codes:
    - EXPLICIT_MAINTENANCE_ROLE
    - POST_ASCT_MAINTENANCE
    - MAINTENANCE_INTENT_MISSING
    - RELAPSE_BEFORE_MAINTENANCE
  evidence_strength: strong
  source_ids:
    - SRC_IMF_FRONTLINE
    - SRC_MMRF_TREATMENT
    - SRC_ASTCT_2022
    - SRC_ROLE_CONSOLIDATION_MAINTENANCE_2024
  known_failure_modes:
    - Maintenance not labeled in source system
    - Delayed transplant causing temporal ambiguity
  last_verified: 2026-07-13

mm_lot_temporary_component_suppression_reintroduced:
  title: Temporary component suppression with later reintroduction stays in the same line
  scope: multiple myeloma, toxicity management and temporary de-escalation
  required_inputs:
    - drug_level_administration_history
    - administration_dates
    - reintroduction_observed
  applies_when: >
    One or more drugs in a combination regimen are held or suppressed while another
    component continues, and the suppressed component later reappears without an
    intervening clearly new regimen.
  default_decision: SAME_LINE
  exceptions:
    - If a clearly new salvage regimen starts before reintroduction, apply new-line logic.
    - If the component never returns and maintenance or toxicity intent is not documented, send to HUMAN_REVIEW.
  reason_codes:
    - TEMP_COMPONENT_SUPPRESSION
    - REINTRODUCED_WITHOUT_NEW_REGIMEN
    - NONRETURNING_COMPONENT
  evidence_strength: moderate
  source_ids:
    - SRC_HEALTHCORE_WORKSHOP_2022
    - SRC_OHDSI_LOT_2024
  known_failure_modes:
    - Sparse outpatient medication capture
    - Long undocumented off-treatment periods
  last_verified: 2026-07-13

mm_lot_drop_component_no_context:
  title: Partial drug removal without context requires review
  scope: multiple myeloma, partial regimen simplification
  required_inputs:
    - prior_active_drug_set
    - next_active_drug_set
    - regimen_role_context
    - progression_context
  applies_when: >
    One or more active anti-myeloma drugs stop while at least one core anti-myeloma drug
    continues, but there is no explicit maintenance, toxicity, suppression-and-return, or
    transplant context.
  default_decision: HUMAN_REVIEW
  exceptions:
    - If maintenance or planned de-escalation is explicit, apply the maintenance rule.
    - If all active drugs stop and a wholly new set starts, apply the full-switch rule.
  reason_codes:
    - PARTIAL_DROP_NO_CONTEXT
    - POSSIBLE_MAINTENANCE
    - POSSIBLE_TOXICITY_DEESCALATION
  evidence_strength: conflicting
  source_ids:
    - SRC_AILAWADHI_2024
    - SRC_HEALTHCORE_WORKSHOP_2022
    - SRC_FALCHETTO_2024
  known_failure_modes:
    - Maintenance intent absent from structured data
    - Simplification after response mistaken for progression
  last_verified: 2026-07-13

mm_lot_dose_or_route_change_only:
  title: Dose reduction or route change alone does not start a new line
  scope: multiple myeloma, within-drug modifications
  required_inputs:
    - canonical_drug_names
    - dose_fields
    - route_fields
    - administration_dates
  applies_when: >
    The canonical drug identity is unchanged and the observable change is limited to dose,
    schedule intensity, or route of administration.
  default_decision: SAME_LINE
  exceptions:
    - If the apparent route change actually reflects a different drug product with a different canonical active agent, re-normalize first.
    - If additional active drugs are introduced beyond the build-up window, apply the add-new-agent rule.
  reason_codes:
    - DOSE_CHANGE_ONLY
    - ROUTE_CHANGE_ONLY
    - SAME_CANONICAL_ACTIVE_AGENT
  evidence_strength: moderate
  source_ids:
    - SRC_OHDSI_LOT_2024
    - SRC_SAINI_2021
  known_failure_modes:
    - Product-level normalization errors
    - Different formulations stored as separate drug identities
  last_verified: 2026-07-13

mm_lot_short_hold_below_gap_threshold:
  title: Temporary hold shorter than the configured discontinuation gap stays in the same line
  scope: multiple myeloma, interruptions and holidays
  required_inputs:
    - last_administration_date
    - next_administration_date
    - configured_gap_days
    - same_regimen_restart_flag
  applies_when: >
    Treatment resumes with the same regimen after an interruption shorter than the
    project's configured discontinuation gap.
  default_decision: SAME_LINE
  exceptions:
    - If progression or salvage intent is explicitly documented during the gap, send to HUMAN_REVIEW.
    - If the pause is part of transplant recovery, apply MM transplant workflow rules first.
  reason_codes:
    - SHORT_HOLD_BELOW_GAP
    - GAP_PARAMETERIZED
    - SAME_REGIMEN_RESUME
  evidence_strength: moderate
  source_ids:
    - SRC_HESS_2021
    - SRC_ONASANYA_2026
    - SRC_FALCHETTO_2024
  known_failure_modes:
    - Missing administrations create false long gaps
    - Oral fills do not equal actual use
  last_verified: 2026-07-13

mm_lot_long_gap_restart_same_regimen:
  title: Restart of the same regimen after a long gap starts a new line
  scope: multiple myeloma, discontinuation and rechallenge
  required_inputs:
    - last_administration_date
    - next_administration_date
    - configured_gap_days
    - same_regimen_restart_flag
  applies_when: >
    The same active regimen is restarted after an interruption equal to or longer than the
    configured discontinuation gap.
  default_decision: NEW_LINE
  exceptions:
    - If the apparent gap is explained by transplant collection, conditioning, ASCT, or other covered MM workflow, apply those rules first.
    - If dates are unreliable, send to HUMAN_REVIEW.
  reason_codes:
    - LONG_GAP_RESTART_SAME_REGIMEN
    - DISCONTINUATION_THRESHOLD_MET
    - EXCLUDE_TRANSPLANT_WORKFLOW
  evidence_strength: moderate
  source_ids:
    - SRC_HESS_2021
    - SRC_ONASANYA_2026
    - SRC_FALCHETTO_2024
  known_failure_modes:
    - Incomplete dates
    - Uncaptured oral therapy bridging the gap
  last_verified: 2026-07-13

mm_lot_progression_with_regimen_change:
  title: Documented relapse or progression followed by regimen change is a new line
  scope: multiple myeloma, clinically annotated transitions
  required_inputs:
    - progression_or_relapse_flag
    - progression_date
    - regimen_change_flag
    - administration_dates
  applies_when: >
    IMWG-consistent progression, relapse, or treatment failure is documented and a
    substantive anti-myeloma regimen change follows.
  default_decision: NEW_LINE
  exceptions:
    - If the only change is transition to conditioning, lymphodepletion, or planned maintenance, apply the disease-specific workflow rules first.
  reason_codes:
    - PROGRESSION_WITH_CHANGE
    - RELAPSE_WITH_CHANGE
    - REFRACTORY_WITH_CHANGE
  evidence_strength: strong
  source_ids:
    - SRC_IMWG_RESPONSE_2016
    - SRC_IMF_RELAPSE_2025
  known_failure_modes:
    - Progression date present but therapy unchanged
    - Ambiguous relapse wording in notes
  last_verified: 2026-07-13

mm_lot_progression_without_regimen_change:
  title: Documented relapse or progression without regimen change requires review
  scope: multiple myeloma, clinically annotated but treatment-static episodes
  required_inputs:
    - progression_or_relapse_flag
    - progression_date
    - regimen_change_flag
  applies_when: >
    Progression, relapse, or refractory disease is documented, but no substantive anti-myeloma
    regimen change is yet observable in the treatment records.
  default_decision: HUMAN_REVIEW
  exceptions:
    - None
  reason_codes:
    - PROGRESSION_NO_CHANGE
    - RELAPSE_NO_CHANGE
    - NEED_CLINICAL_ADJUDICATION
  evidence_strength: moderate
  source_ids:
    - SRC_IMF_RELAPSE_2025
    - SRC_PROGRESS_PROXY_2023
    - SRC_OXNARD_2012
  known_failure_modes:
    - Delay between decision and administration
    - Progression event abstracted from notes but therapy line not yet switched
  last_verified: 2026-07-13

mm_lot_frontline_transplant_workflow:
  title: Standard frontline MM transplant workflow remains one line
  scope: newly diagnosed transplant-eligible multiple myeloma
  required_inputs:
    - regimen_role_sequence
    - transplant_tokens
    - administration_dates
    - relapse_progression_context
  applies_when: >
    The patient follows a standard frontline sequence of induction with or without stem-cell
    mobilization/collection, conditioning high-dose melphalan, autologous stem-cell transplant,
    consolidation if present, and maintenance, without intervening relapse or salvage switch.
  default_decision: SAME_LINE
  exceptions:
    - If relapse or clear salvage intent is documented before the next phase starts, send to HUMAN_REVIEW.
    - Tandem or salvage transplant outside the frontline context requires HUMAN_REVIEW.
  reason_codes:
    - NDMM_FRONTLINE_TRANSPLANT_WORKFLOW
    - INDUCTION_TO_ASCT_TO_MAINTENANCE
    - TANDEM_OR_SALVAGE_TRANSPLANT_REVIEW
  evidence_strength: strong
  source_ids:
    - SRC_IMF_FRONTLINE
    - SRC_MMRF_TREATMENT
    - SRC_ASTCT_2022
    - SRC_RICHTER_2023
  known_failure_modes:
    - Delayed transplant after long induction
    - Missing transplant date
    - Salvage transplant mistaken for frontline ASCT
  last_verified: 2026-07-13

mm_lot_high_dose_melphalan_conditioning:
  title: High-dose melphalan adjacent to ASCT is conditioning, not a new line
  scope: multiple myeloma, transplant conditioning
  required_inputs:
    - melphalan_dose_or_token
    - transplant_tokens
    - temporal_adjacency_to_asct
  applies_when: >
    High-dose melphalan is recorded in the immediate conditioning window adjacent to planned
    autologous stem-cell transplant.
  default_decision: SAME_LINE
  exceptions:
    - If melphalan is given outside transplant adjacency and appears to be standalone anti-myeloma therapy, send to HUMAN_REVIEW.
  reason_codes:
    - HDM_CONDITIONING
    - MEL200_ADJACENT_TO_ASCT
    - NONADJACENT_MELPHALAN_REVIEW
  evidence_strength: strong
  source_ids:
    - SRC_ASTCT_CONDITIONING_2022
    - SRC_ASCO_OH_2025
    - SRC_ALBERTA_MM_TRANSPLANT
  known_failure_modes:
    - Melphalan dose not stored
    - ASCT not recorded in same source
  last_verified: 2026-07-13

mm_lot_post_transplant_maintenance:
  title: Maintenance started after ASCT recovery remains the same line
  scope: multiple myeloma, post-ASCT maintenance
  required_inputs:
    - transplant_date
    - maintenance_start_date
    - maintenance_role_flag
    - relapse_progression_context
  applies_when: >
    Maintenance therapy begins after ASCT recovery in the absence of an intervening relapse,
    progression event, or clearly new salvage regimen.
  default_decision: SAME_LINE
  exceptions:
    - If maintenance intent cannot be inferred and the start occurs after a very prolonged gap, send to HUMAN_REVIEW.
  reason_codes:
    - POST_ASCT_MAINTENANCE
    - NO_INTERVENING_RELAPSE
    - MAINTENANCE_CONTEXT_UNCERTAIN
  evidence_strength: strong
  source_ids:
    - SRC_IMF_FRONTLINE
    - SRC_ASTCT_2022
    - SRC_ROLE_CONSOLIDATION_MAINTENANCE_2024
  known_failure_modes:
    - Missing transplant date
    - Maintenance begun after undocumented relapse
  last_verified: 2026-07-13

mm_lot_cart_lymphodepletion_and_infusion:
  title: CAR-T lymphodepletion and CAR-T administration are one line
  scope: multiple myeloma, ide-cel and cilta-cel workflows
  required_inputs:
    - cellular_therapy_product
    - lymphodepletion_drugs
    - lymphodepletion_dates
    - cart_infusion_date
  applies_when: >
    Fludarabine/cyclophosphamide lymphodepletion occurs in labeled temporal proximity to
    ide-cel or cilta-cel infusion.
  default_decision: SAME_LINE
  exceptions:
    - If lymphodepletion occurs but the intended CAR-T product is never administered, send to HUMAN_REVIEW.
    - If another active anti-myeloma regimen starts between lymphodepletion and infusion, send to HUMAN_REVIEW.
  reason_codes:
    - CART_LD_PLUS_INFUSION
    - FLU_CY_LYMPHODEPLETION
    - CART_NOT_ADMINISTERED
  evidence_strength: strong
  source_ids:
    - SRC_ABECMA_PI
    - SRC_CARVYKTI_PI
    - SRC_IMWG_CART_2024
  known_failure_modes:
    - CAR-T product delayed or canceled
    - Flu/Cy used for another reason
  last_verified: 2026-07-13

mm_lot_bispecific_step_up_phase:
  title: Bispecific step-up dosing and first full dose are one line
  scope: multiple myeloma, teclistamab, elranatamab, talquetamab
  required_inputs:
    - bispecific_agent
    - step_up_dates
    - first_full_dose_date
  applies_when: >
    The record shows labeled step-up dosing followed by the first treatment dose for an
    approved bispecific antibody in multiple myeloma.
  default_decision: SAME_LINE
  exceptions:
    - If the agent identity is uncertain or the labeled step-up sequence is incomplete, send to HUMAN_REVIEW.
  reason_codes:
    - BISPECIFIC_STEP_UP
    - FIRST_FULL_DOSE
    - AGENT_OR_SEQUENCE_UNCERTAIN
  evidence_strength: strong
  source_ids:
    - SRC_TECVAYLI_PI
    - SRC_ELREXFIO_PI
    - SRC_TALVEY_PI
    - SRC_IMWG_BISPECIFIC_2024
  known_failure_modes:
    - Partial inpatient capture of step-up doses
    - Brand and generic names both present without deduplication
  last_verified: 2026-07-13

mm_lot_bridging_or_holding_before_cart:
  title: Holding or bridging therapy before CAR-T requires review
  scope: multiple myeloma, pre-CAR-T disease-control therapy
  required_inputs:
    - leukapheresis_date
    - lymphodepletion_date
    - bridging_or_holding_dates
    - regimen_role_context
  applies_when: >
    Anti-myeloma therapy is given between progression and leukapheresis as holding therapy,
    or after leukapheresis and before lymphodepletion as bridging therapy.
  default_decision: HUMAN_REVIEW
  exceptions:
    - None
  reason_codes:
    - PRE_CAR_T_HOLDING
    - PRE_CAR_T_BRIDGING
    - LOT_COUNTING_CONVENTION_REQUIRED
  evidence_strength: conflicting
  source_ids:
    - SRC_IMWG_CART_2024
    - SRC_RASCHE_BRIDGING_2026
    - SRC_AILAWADHI_BRIDGING_2026
  known_failure_modes:
    - Leukapheresis date absent
    - Bridging intent only in notes
    - True salvage therapy misidentified as bridge
  last_verified: 2026-07-13

mm_lot_steroid_only_episode:
  title: Steroid-only episode requires review unless intent is explicit
  scope: multiple myeloma, dexamethasone or prednisone monotherapy episodes
  required_inputs:
    - canonical_drug_names
    - regimen_role_context
    - symptom_or_bridge_context
    - dates
  applies_when: >
    The observable treatment consists only of corticosteroid therapy without another active
    anti-myeloma drug.
  default_decision: HUMAN_REVIEW
  exceptions:
    - If explicit documentation states the steroid course is supportive-only or symptom-control-only, do not count it as a line.
    - If explicit documentation states steroid monotherapy is being used as active anti-myeloma therapy, classify as NEW_LINE.
  reason_codes:
    - STEROID_ONLY_AMBIGUOUS
    - SUPPORTIVE_ONLY_EXCLUSION
    - ACTIVE_STEROID_MONOTHERAPY
  evidence_strength: limited
  source_ids:
    - SRC_CANCER_CA_SUPPORTIVE
    - SRC_RASCHE_BRIDGING_2026
  known_failure_modes:
    - Dexamethasone entered without companion oral drug
    - Emergency pulses for cord compression or renal failure
  last_verified: 2026-07-13

mm_lot_investigational_or_unmappable_regimen:
  title: Investigational or unmappable regimen requires review
  scope: multiple myeloma, trial and emerging-therapy records
  required_inputs:
    - raw_regimen_text
    - normalized_drug_match_status
    - trial_identifier
  applies_when: >
    The regimen contains an investigational code, emerging therapy, or unmappable agent
    that cannot be placed confidently into a known anti-myeloma regimen role.
  default_decision: HUMAN_REVIEW
  exceptions:
    - If the investigational regimen can be normalized to a known agent class and role with high confidence, classify using the applicable role rule.
  reason_codes:
    - INVESTIGATIONAL_CODE
    - UNMAPPABLE_AGENT
    - LOW_CONFIDENCE_NORMALIZATION
  evidence_strength: limited
  source_ids:
    - SRC_AILAWADHI_2024
    - SRC_IMWG_BISPECIFIC_2024
    - SRC_IMWG_CART_2024
  known_failure_modes:
    - Trial codenames without active moiety
    - Misspelled brand or generic names
  last_verified: 2026-07-13

mm_lot_rechallenge_after_intervening_line:
  title: Rechallenge after an intervening line is a new line
  scope: multiple myeloma, retreatment episodes
  required_inputs:
    - canonical_regimen_identity
    - prior_line_history
    - intervening_line_count
    - administration_dates
  applies_when: >
    A regimen previously used in an earlier line is restarted after one or more clearly
    intervening lines have already been delivered.
  default_decision: NEW_LINE
  exceptions:
    - If the apparent intervening activity consists only of transplant workflow or bispecific/CAR-T initiation steps that belong to the same overarching line, send to HUMAN_REVIEW.
  reason_codes:
    - RECHALLENGE_AFTER_INTERVENING_LINE
    - PRIOR_REGIMEN_REUSED
    - OVERARCHING_WORKFLOW_AMBIGUITY
  evidence_strength: moderate
  source_ids:
    - SRC_IMF_RELAPSE_2025
    - SRC_HESS_2021
  known_failure_modes:
    - Historical lines incompletely captured
    - Regimen identity changed only by brand/generic spelling
  last_verified: 2026-07-13
```

## Normalization assets and regimen-role taxonomy

The normalization tables below are designed to be copied into a pipeline. They synthesize naming patterns from myeloma treatment overviews, transplant guidance, and FDA labeling for currently relevant MM cellular and bispecific therapies. The role taxonomy that follows explicitly marks which roles are inferable from drug composition alone and which require context such as dates, transplant status, or intent. citeturn27view0turn37search1turn19search0turn21search0turn21search1turn22search1turn48search0turn48search2turn48search6

### Canonical drugs

```csv
canonical_generic_name,common_aliases,drug_class,normalization_notes
lenalidomide,"revlimid|len|lena|r",IMiD,"Normalize brand and single-letter shorthand to lenalidomide"
pomalidomide,"pomalyst|pomalid|pom|p",IMiD,"Normalize shorthand carefully because P may also mean prednisone or protocol placeholders"
thalidomide,"thalomid|thal|t",IMiD,"Common in VTd, DVTd, DT-PACE, VTD-PACE"
bortezomib,"velcade|bor|btz|v",Proteasome inhibitor,"Do not split SQ and IV into separate agents"
carfilzomib,"kyprolis|carfil|k",Proteasome inhibitor,"K-based regimen abbreviations usually indicate carfilzomib"
ixazomib,"ninlaro|ixa|ixa-z",Proteasome inhibitor,"Often oral maintenance or relapse therapy"
daratumumab,"darzalex|darzalex faspro|dara|d",Anti-CD38 monoclonal antibody,"Normalize daratumumab/hyaluronidase-fihj to daratumumab"
isatuximab,"sarclisa|isa|isatuximab-irfc",Anti-CD38 monoclonal antibody,"Common in Isa-Kd and Isa-Pd"
elotuzumab,"empliciti|elo",SLAMF7 monoclonal antibody,"Less common in recent sequencing but still normalize"
melphalan,"alkeran|hdm|mel200|mel140",Alkylator,"If temporally adjacent to ASCT, usually conditioning token rather than standalone regimen anchor"
cyclophosphamide,"cytoxan|cyclo|cy|c",Alkylator,"Appears in CyBorD/VCd, DCEP, PACE, lymphodepletion combinations"
fludarabine,"fludara|flu",Purine analog,"Most important as CAR-T lymphodepletion component"
dexamethasone,"dex|dexa|dxm",Corticosteroid,"Normalize as active MM drug when clearly part of regimen; may also appear in supportive or bridge contexts"
prednisone,"pred|pdn",Corticosteroid,"May be active or supportive; intent-sensitive"
prednisolone,"predsol|prednisol",Corticosteroid,"Intent-sensitive"
selinexor,"xpovio|sel",XPO1 inhibitor,"Commonly paired with dexamethasone-based relapse regimens"
belantamab_mafodotin,"blenrep|belamaf|belantamab",BCMA-directed antibody-drug conjugate,"Normalize historical brand use"
bendamustine,"treanda|bendeka|benda",Alkylator,"Occasional salvage or bridge component"
doxorubicin,"adriamycin|doxo",Anthracycline,"PACE-family regimens"
liposomal_doxorubicin,"doxil|pegylated liposomal doxorubicin|pld",Anthracycline,"Normalize separately from conventional doxorubicin if source distinguishes"
etoposide,"vp-16|etop",Topoisomerase inhibitor,"PACE-family regimens"
cisplatin,"cis|cdDP",Platinum,"PACE-family regimens"
vincristine,"oncovin|vcr",Vinca alkaloid,"Historic VAD-like records; normalize if present"
teclistamab,"tecvayli|teclistamab-cqyv",BCMAxCD3 bispecific antibody,"Step-up phase belongs to same line"
elranatamab,"elrexfio|elranatamab-bcmm",BCMAxCD3 bispecific antibody,"Step-up phase belongs to same line"
talquetamab,"talvey|talquetamab-tgvs",GPRC5DxCD3 bispecific antibody,"Step-up phase belongs to same line"
idecabtagene_vicleucel,"abecma|ide-cel|bb2121",BCMA-directed CAR-T,"Group with preceding labeled lymphodepletion"
ciltacabtagene_autoleucel,"carvykti|cilta-cel|jnj-68284528",BCMA-directed CAR-T,"Group with preceding labeled lymphodepletion"
venetoclax,"venclexta|ven",BCL2 inhibitor,"In MM often biomarker-selected or investigational/off-label context"
mezigdomide,"mezigdomide|mezi",Cereblon E3 ligase modulator,"Investigational/emerging"
iberdomide,"iberdomide|iber",Cereblon E3 ligase modulator,"Investigational/emerging"
linvoseltamab,"lynozyfic|linvoseltamab-gcpt",BCMAxCD3 bispecific antibody,"Normalize if present in recent records"
```

### Common regimen abbreviations

```csv
regimen_abbreviation,alternate_forms,canonical_components,primary_role_pattern,normalization_notes
VRd,"RVd|VRD|RVD","bortezomib|lenalidomide|dexamethasone","induction or relapse","Normalize VRd and RVd to same regimen"
D-VRd,"Dara-VRd|DVRd|D-RVd","daratumumab|bortezomib|lenalidomide|dexamethasone","frontline induction/consolidation","Common transplant-eligible frontline regimen"
KRd,"KRD|K-Rd","carfilzomib|lenalidomide|dexamethasone","induction or relapse","Role depends on line and transplant context"
DRd,"D-Rd|Dara-Rd","daratumumab|lenalidomide|dexamethasone","frontline non-transplant or relapse","Context required"
DKd,"D-Kd|Dara-Kd","daratumumab|carfilzomib|dexamethasone","relapse","Context required"
Isa-Kd,"IsaKd|Isa-Kd","isatuximab|carfilzomib|dexamethasone","relapse","Anti-CD38 plus PI doublet backbone"
DVd,"D-Vd|Dara-Vd","daratumumab|bortezomib|dexamethasone","relapse","Do not confuse with D-VRd"
Pd,"P-d|Pom-d","pomalidomide|dexamethasone","relapse","Very common simplified relapse or maintenance-like salvage"
DPd,"D-Pd|Dara-Pd","daratumumab|pomalidomide|dexamethasone","relapse","Context required"
Isa-Pd,"IsaPd|Isa-Pom-d","isatuximab|pomalidomide|dexamethasone","relapse","Context required"
Kd,"K-d","carfilzomib|dexamethasone","relapse","Sometimes bridge; review intent when near CAR-T"
Rd,"R-d|Len-d","lenalidomide|dexamethasone","frontline or relapse","Context required"
CyBorD,"VCd|VCD|CyBord","cyclophosphamide|bortezomib|dexamethasone","induction or relapse","Normalize VCd and CyBorD together"
VTd,"VTD|V-Td","bortezomib|thalidomide|dexamethasone","induction","More common in older/European-style records"
DVTd,"D-VTd","daratumumab|bortezomib|thalidomide|dexamethasone","frontline induction","Context required"
DVCd,"Dara-CyBorD|Dara-VCd","daratumumab|cyclophosphamide|bortezomib|dexamethasone","frontline induction","Normalize Dara-CyBorD to DVCd components"
DCEP,"","dexamethasone|cyclophosphamide|etoposide|cisplatin","salvage or bridge","Often aggressive cytoreduction/bridge"
DT-PACE,"DTPACE|DTPACE","dexamethasone|thalidomide|cisplatin|doxorubicin|cyclophosphamide|etoposide","salvage or bridge","Context required; may precede transplant or CAR-T"
VTD-PACE,"VTDPACE|VTD-PACE","bortezomib|thalidomide|dexamethasone|cisplatin|doxorubicin|cyclophosphamide|etoposide","salvage or bridge","Context required"
PACE,"","cisplatin|doxorubicin|cyclophosphamide|etoposide","salvage or bridge","Rare standalone shorthand; often with modifier"
D-R_maintenance,"DR maintenance|D-R","daratumumab|lenalidomide","maintenance","Usually post-transplant or response-adapted maintenance"
R_maintenance,"len maintenance|R","lenalidomide","maintenance","Keep same line after frontline ASCT when planned"
Tecvayli,"","teclistamab","bispecific immunotherapy","Step-up doses and full doses are one line"
Elrexfio,"","elranatamab","bispecific immunotherapy","Step-up doses and full doses are one line"
Talvey,"","talquetamab","bispecific immunotherapy","Step-up doses and full doses are one line"
ABECMA,"ide-cel","idecabtagene_vicleucel","cellular therapy","Group with preceding labeled Flu/Cy lymphodepletion"
CARVYKTI,"cilta-cel","ciltacabtagene_autoleucel","cellular therapy","Group with preceding labeled Flu/Cy lymphodepletion"
```

### Procedure, transplant, CAR-T, and lymphodepletion tokens

```csv
token_type,canonical_token,aliases,classification_use
non_drug_procedure,stem_cell_mobilization,"mobilization|g-csf mobilization|plerixafor mobilization","procedural role; does not start new LOT by itself"
non_drug_procedure,leukapheresis,"apheresis|collection|pbsc collection|stem-cell collection","procedural role; does not start new LOT by itself"
transplant,autologous_stem_cell_transplant,"ASCT|auto-SCT|auto-HCT|PBSCT|auto-PBSCT|stem-cell rescue","transplant role; part of MM frontline workflow unless salvage/tandem context requires review"
transplant,allogeneic_stem_cell_transplant,"allo-SCT|allo-HCT","rare in MM; usually review context"
conditioning,high_dose_melphalan,"HDM|Mel200|Mel140|conditioning melphalan","conditioning role; same LOT when adjacent to ASCT"
consolidation,consolidation,"post-transplant consolidation|consolidation therapy","same LOT in standard MM transplant workflow"
maintenance,maintenance,"maintenance therapy|post-ASCT maintenance|continuous therapy","same LOT when explicit planned maintenance"
bridging,holding_therapy,"pre-apheresis therapy|holding","requires context; default HUMAN_REVIEW"
bridging,bridging_therapy,"bridge|BT|post-apheresis bridge","requires context; default HUMAN_REVIEW"
lymphodepletion,flu_cy_lymphodepletion,"lymphodepletion|LD chemo|Flu/Cy|fludarabine/cyclophosphamide","same LOT as CAR-T when label-consistent"
cellular_therapy,car_t_infusion,"CAR-T|CAR T infusion|ide-cel infusion|cilta-cel infusion","same LOT with preceding lymphodepletion"
investigational,trial_regimen,"study drug|protocol therapy|investigational regimen","default HUMAN_REVIEW unless mapped confidently"
supportive_or_local,radiation_only,"RT|radiotherapy|bridging radiation","not systemic LOT by itself unless project defines otherwise"
supportive_or_local,steroid_only_episode,"dex pulse|prednisone pulse|dex only","intent-sensitive; default HUMAN_REVIEW"
```

### JSON-compatible regimen-role taxonomy

Most roles in myeloma require **context**, not just drug composition. The exceptions are conditioning, lymphodepletion, transplant, and cellular therapy, which are often identifiable from token plus temporal structure. Maintenance is sometimes inferable from composition alone only when the pattern is a canonical post-ASCT maintenance regimen and the transplant date is known. citeturn27view0turn19search0turn20search4turn21search0turn48search0

```json
{
  "induction": {
    "description": "Initial anti-myeloma treatment intended to reduce disease burden before plateau, transplant, or ongoing therapy.",
    "inferable_from_composition_alone": "sometimes",
    "context_required": true,
    "typical_examples": ["VRd", "D-VRd", "KRd", "CyBorD", "VTd"],
    "notes": "Requires diagnosis setting and transplant candidacy context."
  },
  "salvage": {
    "description": "A new active anti-myeloma regimen used after relapse, progression, or refractory disease.",
    "inferable_from_composition_alone": "rarely",
    "context_required": true,
    "typical_examples": ["DKd", "Isa-Kd", "DPd", "Pd", "DCEP", "DT-PACE", "VTD-PACE"],
    "notes": "Usually requires prior-line history or explicit relapse/progression context."
  },
  "mobilization": {
    "description": "Stem-cell mobilization and collection workflow before transplant.",
    "inferable_from_composition_alone": "sometimes",
    "context_required": true,
    "typical_examples": ["stem_cell_mobilization", "leukapheresis"],
    "notes": "Procedural; not a separate systemic LOT by default."
  },
  "conditioning": {
    "description": "Preparative therapy immediately preceding transplant.",
    "inferable_from_composition_alone": "often",
    "context_required": true,
    "typical_examples": ["high_dose_melphalan"],
    "notes": "In MM, high-dose melphalan adjacent to ASCT should usually be classified as conditioning."
  },
  "transplant": {
    "description": "Autologous or allogeneic stem-cell transplant procedure.",
    "inferable_from_composition_alone": "yes",
    "context_required": false,
    "typical_examples": ["autologous_stem_cell_transplant", "allogeneic_stem_cell_transplant"],
    "notes": "In MM frontline workflows, ASCT is part of the same LOT as induction/maintenance."
  },
  "consolidation": {
    "description": "Limited-duration therapy after ASCT and before maintenance to deepen response.",
    "inferable_from_composition_alone": "rarely",
    "context_required": true,
    "typical_examples": ["post-transplant consolidation"],
    "notes": "Requires transplant adjacency and role context."
  },
  "maintenance": {
    "description": "Lower-intensity or continuing therapy intended to preserve remission after induction or transplant.",
    "inferable_from_composition_alone": "sometimes",
    "context_required": true,
    "typical_examples": ["R_maintenance", "D-R_maintenance", "ixazomib maintenance"],
    "notes": "Do not infer from single-agent lenalidomide alone unless timing and role context support maintenance."
  },
  "bridging": {
    "description": "Disease-control therapy used while awaiting CAR-T or another definitive next step.",
    "inferable_from_composition_alone": "no",
    "context_required": true,
    "typical_examples": ["bridging_therapy", "holding_therapy", "radiation_only", "steroid_only_episode"],
    "notes": "Role depends on timing relative to progression, leukapheresis, lymphodepletion, and infusion."
  },
  "lymphodepletion": {
    "description": "Label-specified chemotherapy given immediately before CAR-T infusion.",
    "inferable_from_composition_alone": "often",
    "context_required": true,
    "typical_examples": ["flu_cy_lymphodepletion"],
    "notes": "Should be grouped with CAR-T administration in the same LOT."
  },
  "cellular_therapy": {
    "description": "CAR-T administration or similar cellular therapy event.",
    "inferable_from_composition_alone": "yes",
    "context_required": false,
    "typical_examples": ["ABECMA", "CARVYKTI", "car_t_infusion"],
    "notes": "Represents the principal therapy event; preceding labeled lymphodepletion stays in same line."
  },
  "investigational_unknown": {
    "description": "Regimen with trial code, emerging agent, or ambiguous composition/intent.",
    "inferable_from_composition_alone": "no",
    "context_required": true,
    "typical_examples": ["trial_regimen", "unmappable agent string"],
    "notes": "Default to HUMAN_REVIEW until normalized."
  }
}
```

## Decision matrix and source manifest

The matrix below is intended as a high-level operational summary. Where the literature is disease-specific and strong, the default can be deterministic. Where the literature is indirect, conflicting, or strongly intent-dependent, the default should be review-based. citeturn27view0turn29search1turn21search0turn48search0turn47search4turn40search9

### Decision matrix

| transition_type | default_result | required_context | key_exceptions | source_ids |
|---|---|---|---|---|
| same drug set continues | SAME_LINE | dates, normalized drug names | long gap beyond configured threshold can still create NEW_LINE | SRC_HESS_2021, SRC_ONASANYA_2026 |
| add active drug within 28 days of line start | SAME_LINE | line start date | explicit salvage escalation inside window → HUMAN_REVIEW | SRC_HESS_2021, SRC_JMIR_VALIDATION_2021 |
| add active drug after build window | NEW_LINE | line start date, role context | post-ASCT maintenance, CAR-T lymphodepletion, bispecific step-up are exceptions | SRC_HESS_2021, SRC_ASTCT_2022, SRC_ABECMA_PI, SRC_TECVAYLI_PI |
| drop one or more drugs with explicit maintenance or de-escalation intent | SAME_LINE | maintenance/de-escalation context | if relapse/progression already documented, use progression rule | SRC_IMF_FRONTLINE, SRC_MMRF_TREATMENT, SRC_ASTCT_2022 |
| drop one or more drugs with no role context | HUMAN_REVIEW | intent, progression, toxicity context | temporary suppression with later reintroduction → SAME_LINE | SRC_AILAWADHI_2024, SRC_HEALTHCORE_WORKSHOP_2022 |
| partial substitution | HUMAN_REVIEW | toxicity intent, progression intent, overlap timing | complete non-overlap switch → NEW_LINE | SRC_FALCHETTO_2024, SRC_OHDSI_LOT_2024 |
| full switch to non-overlapping regimen | NEW_LINE | normalized active drug sets | procedural/conditioning sequences are exceptions | SRC_HESS_2021, SRC_ASTCT_2022 |
| dose reduction only | SAME_LINE | canonical drug identity | if a new active drug is introduced, use drug-addition logic | SRC_OHDSI_LOT_2024 |
| route change only | SAME_LINE | canonical drug identity | if formulation actually changes active agent identity, re-normalize first | SRC_OHDSI_LOT_2024 |
| temporary hold below gap threshold | SAME_LINE | configured gap threshold | if explicit relapse/progression drives switch, review | SRC_HESS_2021, SRC_ONASANYA_2026 |
| restart same regimen after long gap | NEW_LINE | configured gap threshold | transplant workflow gap may be exception | SRC_HESS_2021, SRC_FALCHETTO_2024 |
| documented progression/relapse plus regimen change | NEW_LINE | progression flag/date, regimen change | planned maintenance or CAR-T/transplant workflow exceptions | SRC_IMWG_RESPONSE_2016, SRC_IMF_RELAPSE_2025 |
| documented progression without regimen change | HUMAN_REVIEW | progression flag/date | none | SRC_IMF_RELAPSE_2025, SRC_PROGRESS_PROXY_2023 |
| induction → mobilization/collection → conditioning → ASCT → consolidation → maintenance | SAME_LINE | transplant setting, dates, no intervening relapse | tandem or salvage transplant requires review | SRC_IMF_FRONTLINE, SRC_MMRF_TREATMENT, SRC_ASTCT_2022 |
| high-dose melphalan adjacent to ASCT | SAME_LINE | ASCT adjacency | non-adjacent melphalan requires review | SRC_ASTCT_CONDITIONING_2022, SRC_ALBERTA_MM_TRANSPLANT |
| maintenance starting after ASCT recovery | SAME_LINE | transplant date, maintenance role, no relapse before start | very prolonged gap without intent label → HUMAN_REVIEW | SRC_IMF_FRONTLINE, SRC_ASTCT_2022 |
| CAR-T lymphodepletion → CAR-T infusion | SAME_LINE | product identity, Flu/Cy dates, infusion date | CAR-T canceled after lymphodepletion → HUMAN_REVIEW | SRC_ABECMA_PI, SRC_CARVYKTI_PI, SRC_IMWG_CART_2024 |
| bispecific step-up doses → first full dose | SAME_LINE | agent identity, labeled step-up sequence | incomplete or uncertain sequence → HUMAN_REVIEW | SRC_TECVAYLI_PI, SRC_ELREXFIO_PI, SRC_TALVEY_PI |
| bridging or holding therapy before CAR-T | HUMAN_REVIEW | leukapheresis, lymphodepletion, infusion timing | none | SRC_IMWG_CART_2024, SRC_RASCHE_BRIDGING_2026 |
| steroid-only bridge | HUMAN_REVIEW | explicit intent | supportive-only excludes LOT; explicit active monotherapy can be NEW_LINE | SRC_CANCER_CA_SUPPORTIVE |
| investigational or unmappable regimen | HUMAN_REVIEW | trial identifier, mapping confidence | if normalized with high confidence, apply role rule | SRC_AILAWADHI_2024, SRC_IMWG_BISPECIFIC_2024, SRC_IMWG_CART_2024 |
| rechallenge after an intervening line | NEW_LINE | full prior-line history | if “intervening line” is actually part of the same MM workflow, review | SRC_HESS_2021, SRC_IMF_RELAPSE_2025 |

### Source manifest

```csv
source_id,title,organization_or_authors,publication_date,source_type,url,relevant_sections,authority_level,rules_supported
SRC_HESS_2021,Defining Treatment Regimens and Lines of Therapy Using Real-World Data in Oncology,"Hess et al.",2021-02-25,peer_reviewed_methodology,"https://www.tandfonline.com/doi/full/10.2217/fon-2020-1041","regimen definition; additions/removals; 90-day gap rule",high,"mm_lot_start_first_active_regimen|mm_lot_regimen_build_up_28d|mm_lot_add_new_agent_outside_build_up|mm_lot_full_nonoverlap_switch|mm_lot_short_hold_below_gap_threshold|mm_lot_long_gap_restart_same_regimen|mm_lot_rechallenge_after_intervening_line"
SRC_SAINI_2021,Determining lines of therapy in patients with solid cancers: a proposed new systematic and comprehensive framework,"Saini & Twelves",2021-04-12,peer_reviewed_methodology,"https://www.nature.com/articles/s41416-021-01319-8","definitions; minimum dataset; clinical progression discussion",high,"mm_lot_dose_or_route_change_only|mm_lot_start_first_active_regimen"
SRC_FALCHETTO_2024,Concepts of lines of therapy in cancer treatment: findings from an expert interview-based study,"Falchetto et al.",2024-05-02,peer_reviewed_methodology,"https://pmc.ncbi.nlm.nih.gov/articles/PMC11094945/","experts on breaks, same therapy restart, ambiguity",moderate,"mm_lot_drop_component_no_context|mm_lot_short_hold_below_gap_threshold|mm_lot_long_gap_restart_same_regimen"
SRC_ONASANYA_2026,Identifying the Initiation of a New Line of Therapy for Metastatic Lung Breast and Colorectal Cancer in Real-World Data: A Scoping Review,"Onasanya et al.",2026-04-27,peer_reviewed_methodology,"https://doi.org/10.1002/pds.70370","variation in gap definitions; algorithm heterogeneity",high,"mm_lot_short_hold_below_gap_threshold|mm_lot_long_gap_restart_same_regimen"
SRC_JMIR_VALIDATION_2021,An Automated Line-of-Therapy Algorithm for Adults With Metastatic Non–Small Cell Lung Cancer: Validation Study Using Blinded Manual Chart Review,"Meng et al.",2021-10-12,peer_reviewed_methodology,"https://medinform.jmir.org/2021/10/e29017/","28-day line regimen defining window",moderate,"mm_lot_regimen_build_up_28d"
SRC_AILAWADHI_2024,Development and validation of algorithms for identifying lines of therapy in multiple myeloma using real-world data,"Ailawadhi et al.",2024-01-17,peer_reviewed_mm_methodology,"https://pubmed.ncbi.nlm.nih.gov/38231002/","MM line-identification validation; mini-LOT concept; maintenance handling",high,"mm_lot_add_new_agent_outside_build_up|mm_lot_drop_component_no_context|mm_lot_investigational_or_unmappable_regimen"
SRC_BERGER_EHA_2021,Validation of algorithms to identify first-line therapy and use of maintenance therapy in multiple myeloma,"Berger et al.",2021-06-09,conference_abstract_mm_methodology,"https://library.ehaweb.org/eha/2021/eha2021-virtual-congress/324313/ariel.berger.validation.of.algorithms.to.identify.first-line.therapy.html","maintenance and LOT1 validation",moderate,"mm_lot_explicit_de_escalation_to_maintenance|mm_lot_post_transplant_maintenance"
SRC_ASTCT_2022,ASTCT Clinical Practice Recommendations for Transplantation and Cellular Therapies in Multiple Myeloma,"Dhakal et al./ASTCT",2022-06-01,professional_society_guideline,"https://pubmed.ncbi.nlm.nih.gov/35306217/","early autologous transplant as consolidation; maintenance after ASCT; CAR-T sequencing",very_high,"mm_lot_add_new_agent_outside_build_up|mm_lot_explicit_de_escalation_to_maintenance|mm_lot_frontline_transplant_workflow|mm_lot_post_transplant_maintenance"
SRC_ASTCT_CONDITIONING_2022,High-dose conditioning regimens used prior to autologous transplantation for multiple myeloma,"Ali et al./ASTCT journal",2022-06-23,peer_reviewed_review,"https://pubmed.ncbi.nlm.nih.gov/35750284/","high-dose melphalan as standard MM conditioning",high,"mm_lot_high_dose_melphalan_conditioning"
SRC_ASCO_OH_2025,Treatment of Multiple Myeloma: ASCO–Ontario Health Living Guideline,"ASCO and Ontario Health",2025-01-06,professional_society_guideline,"https://www.cancercareontario.ca/sites/ccocancercare/files/guidelines/full/pebc6-a-2023-1v2f.pdf","frontline transplant, high-dose melphalan conditioning, maintenance recommendations",very_high,"mm_lot_high_dose_melphalan_conditioning"
SRC_IMWG_RESPONSE_2016,International Myeloma Working Group consensus criteria for response and minimal residual disease assessment in multiple myeloma,"Kumar et al./IMWG",2016-08-27,professional_consensus,"https://pubmed.ncbi.nlm.nih.gov/27511158/","response, progression, relapse framework",very_high,"mm_lot_progression_with_regimen_change"
SRC_IMF_RELAPSE_2025,What Is Multiple Myeloma Relapse?,"International Myeloma Foundation",2025-06-12,IMWG_aligned_foundation_resource,"https://www.myeloma.org/treatment/relapse-definition","relapse, biochemical relapse, refractory disease, sequencing after relapse",high,"mm_lot_progression_with_regimen_change|mm_lot_progression_without_regimen_change|mm_lot_rechallenge_after_intervening_line"
SRC_IMF_FRONTLINE,First-Line Treatment Options for Multiple Myeloma,"International Myeloma Foundation",2025-07-30,foundation_resource,"https://www.myeloma.org/frontline-treatment-options","one LOT = induction + ASCT + maintenance; post-ASCT maintenance timing",high,"mm_lot_explicit_de_escalation_to_maintenance|mm_lot_frontline_transplant_workflow|mm_lot_post_transplant_maintenance"
SRC_MMRF_TREATMENT,Multiple Myeloma Treatments | Drugs & Treatment by Stage,"Multiple Myeloma Research Foundation",2026-07-13,current_foundation_resource,"https://themmrf.org/diagnosis-and-treatment/treatment-options/","collectively one line; induction-transplant-maintenance sequence",high,"mm_lot_frontline_transplant_workflow|mm_lot_explicit_de_escalation_to_maintenance"
SRC_RICHTER_2023,Real-world multiple myeloma front-line treatment and outcomes associated with stem cell transplant in the first line of therapy,"Richter et al.",2023-11-30,peer_reviewed_mm_real_world,"https://pubmed.ncbi.nlm.nih.gov/38024614/","SCT in first line of therapy",high,"mm_lot_frontline_transplant_workflow"
SRC_ROLE_CONSOLIDATION_MAINTENANCE_2024,Role of Consolidation and Maintenance,"Kumar & Chari",2024-01-01,peer_reviewed_review,"https://www.binasss.sa.cr/abr24/52.pdf","consolidation after ASCT and before maintenance",high,"mm_lot_explicit_de_escalation_to_maintenance|mm_lot_post_transplant_maintenance"
SRC_ALBERTA_MM_TRANSPLANT,Multiple Myeloma: Transplant Eligible,"Alberta Health Services",2023-10-26,regional_guideline,"https://www.albertahealthservices.ca/assets/info/hp/cancer/if-hp-cancer-guide-lyhe013-mmte.pdf","standard conditioning regimen is high-dose melphalan 200 mg/m2",high,"mm_lot_high_dose_melphalan_conditioning"
SRC_ABECMA_PI,Package Insert and Medication Guide - ABECMA,"U.S. FDA",2025-12-04,regulatory_label,"https://www.fda.gov/media/147055/download","lymphodepleting chemotherapy and timing before ide-cel",very_high,"mm_lot_cart_lymphodepletion_and_infusion"
SRC_CARVYKTI_PI,Package Insert and Medication Guide - CARVYKTI,"U.S. FDA",2025-10-01,regulatory_label,"https://www.fda.gov/media/156560/download","lymphodepleting chemotherapy for cilta-cel",very_high,"mm_lot_cart_lymphodepletion_and_infusion"
SRC_IMWG_CART_2024,Consensus guidelines and recommendations for the management and response assessment of CAR T-cell therapy in clinical practice for relapsed and refractory multiple myeloma,"Lin et al./IMWG",2024-06-01,professional_consensus,"https://pubmed.ncbi.nlm.nih.gov/38821074/","bridging therapy, lymphodepletion, response assessment",very_high,"mm_lot_cart_lymphodepletion_and_infusion|mm_lot_bridging_or_holding_before_cart|mm_lot_investigational_or_unmappable_regimen"
SRC_TECVAYLI_PI,TECVAYLI Prescribing Information,"U.S. FDA / Janssen",2022-10-25,regulatory_label,"https://www.accessdata.fda.gov/drugsatfda_docs/label/2022/761291s000lbl.pdf","step-up dosing then full dose until progression",very_high,"mm_lot_bispecific_step_up_phase"
SRC_ELREXFIO_PI,ELREXFIO Prescribing Information,"U.S. FDA",2025-07-16,regulatory_label,"https://www.accessdata.fda.gov/drugsatfda_docs/label/2025/761345s003lbl.pdf","step-up dosing schedule then full-dose treatment",very_high,"mm_lot_bispecific_step_up_phase"
SRC_TALVEY_PI,TALVEY Prescribing Information,"U.S. FDA",2023-08-09,regulatory_label,"https://www.accessdata.fda.gov/drugsatfda_docs/label/2023/761342s000lbl.pdf","step-up phase dosing and first treatment dose",very_high,"mm_lot_bispecific_step_up_phase"
SRC_IMWG_BISPECIFIC_2024,International Myeloma Working Group immunotherapy committee consensus guidelines and recommendations for optimal use of T-cell-engaging bispecific antibodies in multiple myeloma,"Rodriguez-Otero et al./IMWG",2024-05-01,professional_consensus,"https://pubmed.ncbi.nlm.nih.gov/38697166/","bispecific practice guidance",very_high,"mm_lot_bispecific_step_up_phase|mm_lot_investigational_or_unmappable_regimen"
SRC_HEALTHCORE_WORKSHOP_2022,Approaches to Algorithm Development for the Estimation of Lines of Therapy,"ISPOR workshop / HealthCore et al.",2022-10-21,formal_methodology_workshop,"https://www.ispor.org/docs/default-source/euro2022/ispor-eu-2022-workshop-slidesv521oct2022.pdf","drug addition, discontinuation, exchanges, gaps, component suppression",moderate,"mm_lot_temporary_component_suppression_reintroduced|mm_lot_drop_component_no_context"
SRC_OHDSI_LOT_2024,Clinically validated line of therapy algorithm for patients with metastatic NSCLC,"OHDSI poster / Moreira et al.",2024-10-01,validated_algorithm_poster,"https://www.ohdsi.org/wp-content/uploads/2024/10/85-Moreira-Clinically-validated-line-of-therapy-algorithm-for-patients-with-metastatic-Non-Small-Cell-Lung-Cancer-Jie-Yeap.pdf","dose change, route change, stopping one or more drugs while LoT continues",moderate,"mm_lot_temporary_component_suppression_reintroduced|mm_lot_dose_or_route_change_only"
SRC_RASCHE_BRIDGING_2026,Experiences and expert panel consensus on bridging therapy before anti-BCMA CAR-T cell therapy in multiple myeloma,"Rasche et al.",2026-05-01,peer_reviewed_mm_consensus,"https://www.nature.com/articles/s41409-026-02905-1","definition of bridging vs holding; BT goals; timing relative to leukapheresis and lymphodepletion",high,"mm_lot_bridging_or_holding_before_cart|mm_lot_steroid_only_episode"
SRC_AILAWADHI_BRIDGING_2026,Optimizing Selection of Bridging Therapies Prior to CAR-T Therapy for Multiple Myeloma,"Ailawadhi et al.",2026-02-01,peer_reviewed_review,"https://www.clinical-lymphoma-myeloma-leukemia.com/article/S2152-2650(25)00268-X/fulltext","bridging therapy as disease-control before CAR-T",high,"mm_lot_bridging_or_holding_before_cart"
SRC_CANCER_CA_SUPPORTIVE,Supportive therapy for multiple myeloma,"Canadian Cancer Society",2025-04-25,high_quality_foundation_resource,"https://cancer.ca/en/cancer-information/cancer-types/multiple-myeloma/treatment/supportive-therapy","corticosteroids used alone in some patients; intent ambiguity",moderate,"mm_lot_steroid_only_episode"
SRC_PROGRESS_PROXY_2023,Overview of approaches to estimate real-world disease progression in patients with cancer,"Amorrortu et al.",2023-11-10,peer_reviewed_rwd_methods,"https://pmc.ncbi.nlm.nih.gov/articles/PMC10637832/","progression data limitations in RWD",moderate,"mm_lot_progression_without_regimen_change"
SRC_OXNARD_2012,When Progressive Disease Does Not Mean Treatment Failure,"Oxnard et al.",2012-08-10,peer_reviewed_conceptual,"https://pmc.ncbi.nlm.nih.gov/articles/PMC3708548/","progression not always identical to immediate treatment failure",moderate,"mm_lot_progression_without_regimen_change"
SRC_FLATIRON_LOT_2026,Lines of therapy across geographies: Bridging the gap between clinical practice and real-world data,"Flatiron Health",2025-12-01,vendor_specific_methodology,"https://resources.flatiron.com/real-world-evidence/lines-of-therapy-across-geographies-bridging-the-gap-between-clinical-practice-real-world-data","vendor-specific LOT operationalization context",vendor_specific,"mm_lot_regimen_build_up_28d"
```

## Open questions and validation checklist

### Open questions for clinical reviewer adjudication

The literature does **not** resolve the following issues well enough for universal automation, so they should remain on the project’s clinical-review queue. citeturn34search0turn41search13turn47search4turn40search9

- **What discontinuation gap should govern MM rechallenge?** Published LOT methods use materially different thresholds, and none emerged as a universally accepted MM-specific standard. A parameter must be chosen and validated locally. citeturn34search0turn29search1
- **How should partial substitutions be counted in MM?** The literature is stronger for full switches and weaker for toxicity-driven partial swaps. Without explicit intent, many such events remain ambiguous. citeturn41search7turn41search13
- **Should bridging or holding therapy before CAR-T count toward LOT totals, be excluded from LOT totals, or be stored only as role metadata?** Recent MM consensus documents define the role clearly, but do not establish a uniform database-counting convention. citeturn47search4turn25search3
- **How should steroid-only episodes be handled when intent is absent?** Steroids can be active monotherapy, supportive therapy, or bridging. Dates alone are not enough. citeturn40search9
- **How should delayed or salvage ASCT be grouped?** Frontline ASCT is well supported as part of first-line therapy, but tandem ASCT and later salvage ASCT need project-specific conventions and likely reviewer input. citeturn19search0turn27view0
- **How should investigational regimens be represented when only protocol codenames appear?** A normalization layer can map many cases, but some trial records will still require manual curation. citeturn36search3turn25search3
- **Can progression alone create a boundary in the absence of a treatment change?** Published evidence and clinical logic both argue that progression and treatment boundaries are related but not identical. For treatment-sequence classification, progression alone should usually not increment LOT without adjudication. citeturn38search0turn38search16turn49view0

### Validation checklist for adjudicated COTA data

For each hard rule below, the test is deliberately falsifiable: it can succeed or fail against chart-adjudicated truth.

| rule_id | falsifiable_test_on_adjudicated_cota_data | adjudication_fields_required | expected_pass_condition |
|---|---|---|---|
| mm_lot_start_first_active_regimen | Identify earliest active anti-myeloma therapy per patient; compare algorithmic first-line start date to adjudicated LOT1 start date. | adjudicated LOT1 start date, active-vs-supportive flag | ≥90% exact-date match or clinically acceptable ±3 day match threshold |
| mm_lot_regimen_build_up_28d | For episodes where an added drug appears within 28 days of LOT start and no explicit salvage note exists, compare algorithm’s SAME_LINE output to adjudication. | adjudicated regimen membership, explicit salvage flag | High PPV for same-line grouping |
| mm_lot_add_new_agent_outside_build_up | For additions after day 28 that are not maintenance/CAR-T/bispecific exceptions, compare algorithmic NEW_LINE to adjudicated boundary. | adjudicated boundary flag, role context | High PPV for new-line detection |
| mm_lot_full_nonoverlap_switch | Among transitions with zero overlap in normalized active drug sets, compare NEW_LINE output to adjudication. | normalized active drug sets, adjudicated line boundary | Very high specificity for NEW_LINE |
| mm_lot_dose_or_route_change_only | Find records where only dose/schedule/route changed and canonical drug identity remained constant; verify no adjudicated LOT increment. | canonical drug identity, dose/route history, adjudicated line number | False-positive new-line rate near zero |
| mm_lot_short_hold_below_gap_threshold | For same-regimen interruptions shorter than configured gap, verify algorithm preserves line number versus adjudication. | last admin date, restart date, adjudicated line number | High NPV for false discontinuation |
| mm_lot_long_gap_restart_same_regimen | For same-regimen restarts at or beyond configured gap, compare NEW_LINE output to adjudication. | threshold value, restart timing, adjudicated line number | High PPV under chosen project convention |
| mm_lot_progression_with_regimen_change | For adjudicated progression/relapse plus regimen-change events, verify algorithm increments LOT exactly once. | progression date, adjudicated reason for switch, line boundary | High sensitivity and no double-counting |
| mm_lot_frontline_transplant_workflow | In transplant-eligible frontline cases with induction→ASCT→maintenance and no relapse, verify algorithm returns one LOT total across that workflow. | transplant date, conditioning date, maintenance start, adjudicated LOT totals | Near-perfect agreement on single-line grouping |
| mm_lot_high_dose_melphalan_conditioning | Detect melphalan within transplant adjacency window; confirm algorithm does not create a separate LOT relative to adjudication. | melphalan date/dose, ASCT date, adjudicated boundary | False-positive boundary rate near zero |
| mm_lot_post_transplant_maintenance | For post-ASCT maintenance starts without intervening relapse, verify no new LOT is created. | ASCT date, maintenance start, relapse adjudication | High PPV for SAME_LINE |
| mm_lot_cart_lymphodepletion_and_infusion | For ide-cel/cilta-cel cases with labeled Flu/Cy before infusion, verify algorithm groups lymphodepletion and infusion into one LOT. | lymphodepletion dates, infusion date, adjudicated cellular-therapy line | Near-perfect agreement |
| mm_lot_bispecific_step_up_phase | For teclistamab/elranatamab/talquetamab starts with step-up dosing, verify algorithm keeps step-up doses and first full dose in one LOT. | step-up dates, full-dose date, adjudicated line boundary | Near-perfect agreement |
| mm_lot_rechallenge_after_intervening_line | For patients who later reuse a prior regimen after at least one adjudicated intervening line, verify algorithm assigns a new LOT. | full adjudicated line history, normalized regimen identity | High PPV for NEW_LINE |

The best implementation strategy is to freeze the **strong MM rules** as hard logic, expose **gap length** as a parameter, and route the remaining conflict-heavy scenarios—especially partial drops, bridging, steroid-only episodes, delayed transplant, and unmappable investigational regimens—to a review queue. That division aligns with the current evidence base and should minimize false precision in downstream LOT totals. citeturn27view0turn29search1turn47search4turn40search9
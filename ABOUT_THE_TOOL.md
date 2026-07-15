# How the pipeline works
## 1. Prepare the data

Each patient is converted into something like:

Order 1: Drug A + Drug B
Order 2: Drug A + Drug B + Drug C
Order 3: Drug D

Sensitive identifiers are replaced with pseudonymous case IDs.

Patients are divided into five folds. Similar or duplicate patient trajectories are kept together so that one version cannot appear in training while another appears in testing.

## 2. Give the AI a blind case

The AI receives only the drug order.

It does not receive:

reviewer LOT;
COTA LOT;
algorithm LOT;
vendor line numbers;
vendor treatment dates.

This makes it a genuine independent vote, but it also limits performance because many real LOT decisions depend on timing.

## 3. Retrieve examples

For k=3 or k=5, the system finds three or five similar patients from the other four folds.

Those examples show:

drug trajectory → reviewer total LOT

The target patient and related duplicate trajectories are excluded.

## 4. Generate the AI output

The AI predicts whether each transition represents:

same line
new line
insufficient information

It then outputs a total LOT count.

When the evidence is too ambiguous, it can abstain. An abstained numeric answer is retained for diagnostics but does not count as a vote.

## 5. Evaluate and route

The evaluator compares the three values against reviewer consensus.

### Examples:

Algorithm = 4
COTA = 4
AI = 4, usable
→ three-way consensus candidate
Algorithm = 4
COTA = 4
AI abstains
→ human review
Algorithm = 4
COTA = 5
AI = 4
→ disagreement; human review
What we found

Across all 136 patients:

Policy	Accepted	Correct
Algorithm alone	136	89
Algorithm + COTA agree	76	63
Three-way, k=3	10	9
Three-way, k=5	11	10

So three-way agreement increased the point-estimate accuracy from 82.9% to about 90%, but coverage fell from 55.9% to about 8%. It also still produced one incorrect consensus case, so it is not perfectly safe.

The AI itself is not a strong standalone classifier:

k=3 usable AI accuracy: 18/38, or 47.4%
k=5 usable AI accuracy: 21/47, or 44.7%

Its potential value is mainly as a selective confidence filter, not as a replacement for the existing algorithm.

## The most important interpretation

> The AI did not reliably catch errors by explicitly disagreeing.

In fact, when the algorithm and COTA agreed but the AI disagreed:

k=3: algorithm+COTA were correct in 8/8 cases;
k=5: algorithm+COTA were correct in 12/12 cases.

Most of the apparent safety improvement came because the AI abstained frequently, leaving only a small, easier subset.

A good way to phrase that is:

> The AI currently functions more like a conservative selector than an accurate third expert. It identifies a very small subset where all sources agree, but it does not yet provide reliable independent correction of the algorithm or COTA.

## Main limitation

The benchmark is deliberately order-only.

### The AI cannot see:

treatment gaps;
dates;
overlap between therapies;
transplant timing;
maintenance timing;
CAR-T timing;
progression timing.

Those are often essential for determining whether a treatment change represents a new LOT.

So the poor standalone AI performance may partly reflect missing information rather than only poor reasoning.

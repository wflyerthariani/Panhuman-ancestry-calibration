# DMS assay review: score direction and assay quality

**Who this is for:** someone with a molecular biology or genetics background. No coding needed.
**Time needed:** about 5–10 minutes per assay. **Assays to review at this scope (`smoke`):** 23.

## Why this review is needed

The audit compares AVI predictions with deep mutational scanning (DMS) measurements. Each DMS
assay reports a score per variant, but assays disagree on what a *high* score means. In a
growth assay, a high score can mean the variant works (it grows). In a drug-resistance screen,
a high score can mean the variant is broken (it escapes the drug). Before scores can be pooled,
every assay must be oriented so that **a higher score always means more functional damage**
(protocol Step 13). The protocol also requires an explicit judgement of **assay quality**
(Step 12).

A wrong direction silently turns a good predictor into a bad one for that gene, so a person
has to sign off on each assay. The pipeline has gathered the evidence and made a proposal.
Your job is to confirm or correct it.

## The three kinds of evidence on each card

- **A. Calibration.** Some MaveDB score sets have been calibrated, often against ClinVar
  variants, into ranges labelled *normal* and *abnormal* function. If the abnormal range sits
  below the normal range, lower scores mean damage. Only a minority of assays have one.
- **B. Built-in controls.** Nonsense variants (premature stop, e.g. `p.Arg213Ter`) almost always
  destroy function, and synonymous variants (e.g. `p.Arg213=`) almost never do. The pipeline
  compares their scores:
  - **AUC** is the probability that a random nonsense variant scores higher than a random
    synonymous one.
  - **AUC near 1:** high score = damaging. **AUC near 0:** low score = damaging.
  - **AUC near 0.5:** the assay does not separate them, which is a warning sign for quality,
    or the assay measures something where truncation is not "damaging" in the assay's terms.
  - **Caveat:** nonsense variants near the protein's C-terminus may escape nonsense-mediated
    decay and behave like missense variants. Assays that start downstream of the protein may
    show little separation.
- **C. Author text.** Sentences from the MaveDB description that mention the score. They are
  shown verbatim; the pipeline does not interpret them. Please check them against A and B,
  and follow the MaveDB link or the paper when they are unclear.

## What to decide for each assay

Fill in these columns of **`config/dms_curation.tsv`**. It opens in Excel or Google Sheets.
The rows under review have `in_scope = true` and `prescreen_status = candidate`, and they are
listed first.

| Column | Allowed values | Meaning |
|---|---|---|
| `reviewer_decision` | `include` / `exclude` | Should this assay enter the DMS benchmark? |
| `reviewer_direction` | `higher_is_more_damaging` / `lower_is_more_damaging` | The direction of the **raw** MaveDB score. Required when `include`. The pipeline reverses the score if it is `lower_is_more_damaging`. |
| `reviewer_notes` | free text | **Required** if you change the proposed direction, exclude an assay, or have a caveat (e.g. "measures abundance only, not activity"). |
| `reviewer_name` | your name | |
| `reviewed_on` | YYYY-MM-DD | |

Please do **not** edit the other columns. They are regenerated each time the pipeline runs, and
your five columns are kept.

### Questions to ask yourself

1. **Does the assay measure a loss of the protein's normal function?** For example: growth that
   depends on the protein, reporter activity, abundance or stability, binding, or drug response
   through the protein.
   - Assays of gain of function, or of a neomorphic or dominant-negative effect, can still be
     included if "damage" is defined consistently.
   - Say so in the notes. Such assays deserve extra care over direction.
2. **Do A, B and C agree on the direction?** If they conflict, or the proposal is `unknown`,
   decide from the paper and explain in the notes.
3. **Is the quality acceptable?** Warning signs:
   - `weak_control_separation`, or an AUC between about 0.35 and 0.65
   - very few scored missense variants
   - no replicate or error information
   - an assay the authors themselves describe as noisy or preliminary
4. **Is it the right version?** When several score sets come from the same experiment (for
   example nucleotide-level and amino-acid-level scores, or per-condition and combined scores),
   prefer the **amino-acid-level, primary** score. Exclude the others with the note
   "duplicate of <urn>". The protocol pools variants within a gene, so counting the same
   measurement twice would bias results.
5. **Anything that makes this assay unsuitable for comparing populations?** For example: only a
   short domain was tested, or an unusual isoform or construct was used, such as a common
   polymorphism in the backbone (`TP53 (P72R)`). Note it; it matters for Stage 9.

### Flags you may see

| Flag | Meaning |
|---|---|
| `weak_control_separation` | Nonsense vs synonymous AUC is between 0.25 and 0.75 |
| `no_nonsense_synonymous_controls` | Fewer than 10 of either class had a score, or no protein-level HGVS |
| `fewer_than_20_scored_missense` | Too few single missense variants to be useful |
| `genomic_mapping_incomplete` / `_failed` | MaveDB could not map every variant to GRCh38. Stage 6 falls back to protein-level mapping |
| `licence_CC_BY-NC-SA_4.0` | Non-commercial licence (see limitation L-05) |
| `gene_resolution_conflict` | Different identifiers pointed to different genes. Please check the target |

## Proposal summary

- `conflict`: 2
- `higher_is_more_damaging`: 2
- `lower_is_more_damaging`: 12
- `unknown`: 7

## Assays to review

### 1. BRCA1 — Enrich2 nucleotide variant scores for BRCA1 E3

`urn:mavedb:00000003-a-1` · [MaveDB page](https://www.mavedb.org/score-sets/urn:mavedb:00000003-a-1) · PubMed:28784151 (10.1186/s13059-017-1272-5) · licence CC0

- **What was measured:** Nucleotide variant scores for deep mutational scan of the BRCA1 RING domain using autoubiquitination calculated by Enrich2.
- **Assay:** —
- **Variants:** 20724 in MaveDB, 11530 with a numeric score, 1222 single missense
- **Replicates:** replicate columns: SE_PlusE2NewRep3, score_PlusE2NewRep3, SE_PlusE2NewRep4, score_PlusE2NewRep4, SE_PlusE2NewRep5, score_PlusE2NewRep5; error columns: SE, SE_PlusE2NewRep3, SE_PlusE2NewRep4, SE_PlusE2NewRep5, SE_PlusE2Rep3, SE_PlusE2Rep4
- **Genome mapping by MaveDB:** incomplete

| Evidence | What it shows | Implied direction |
|---|---|---|
| A. Calibration | none | — |
| B. Controls | nonsense (n=60, median -2.42) vs synonymous (n=377, median -0.308); AUC 0.16 | lower_is_more_damaging |

C. What the authors say about the score:
   > (no sentence about score meaning found — see the MaveDB page)

**Proposed:** `lower_is_more_damaging` (confidence: medium) · **Flags:** genomic_mapping_incomplete · **Review status:** **not yet reviewed**

### 2. BRCA1 — Enrich2 amino acid variant scores for BRCA1 E3

`urn:mavedb:00000003-a-2` · [MaveDB page](https://www.mavedb.org/score-sets/urn:mavedb:00000003-a-2) · PubMed:28784151 (10.1186/s13059-017-1272-5) · licence CC0

- **What was measured:** Amino acid variant scores for deep mutational scan of the BRCA1 RING domain using autoubiquitination calculated by Enrich2.
- **Assay:** —
- **Variants:** 12316 in MaveDB, 6847 with a numeric score, 3710 single missense
- **Replicates:** replicate columns: SE_PlusE2NewRep3, score_PlusE2NewRep3, SE_PlusE2NewRep4, score_PlusE2NewRep4, SE_PlusE2NewRep5, score_PlusE2NewRep5; error columns: SE, SE_PlusE2NewRep3, SE_PlusE2NewRep4, SE_PlusE2NewRep5, SE_PlusE2Rep3, SE_PlusE2Rep4
- **Genome mapping by MaveDB:** incomplete

| Evidence | What it shows | Implied direction |
|---|---|---|
| A. Calibration | none | — |
| B. Controls | too few controls (nonsense 147, synonymous 1; need 10 each) | — |

C. What the authors say about the score:
   > (no sentence about score meaning found — see the MaveDB page)

**Proposed:** `unknown` (confidence: none) · **Flags:** no_nonsense_synonymous_controls,genomic_mapping_incomplete · **Review status:** **not yet reviewed**

### 3. BRCA1 — Enrich2 nucleotide variant scores for BRCA1 Y2H

`urn:mavedb:00000003-b-1` · [MaveDB page](https://www.mavedb.org/score-sets/urn:mavedb:00000003-b-1) · PubMed:28784151 (10.1186/s13059-017-1272-5) · licence CC0

- **What was measured:** Nucleotide variant scores for deep mutational scan of the BRCA1 RING domain using yeast two-hybrid calculated by Enrich2.
- **Assay:** —
- **Variants:** 20724 in MaveDB, 17165 with a numeric score, 1354 single missense
- **Replicates:** replicate columns: SE_Y2H_1_Rep1, score_Y2H_1_Rep1, SE_Y2H_1_Rep2, score_Y2H_1_Rep2, SE_Y2H_1_Rep3, score_Y2H_1_Rep3; error columns: SE, SE_Y2H_1_Rep1, SE_Y2H_1_Rep2, SE_Y2H_1_Rep3, SE_Y2H_2_Rep1, SE_Y2H_2_Rep2
- **Genome mapping by MaveDB:** incomplete

| Evidence | What it shows | Implied direction |
|---|---|---|
| A. Calibration | none | — |
| B. Controls | nonsense (n=79, median 0.438) vs synonymous (n=396, median -0.0188); AUC 0.65 | — |

C. What the authors say about the score:
   > (no sentence about score meaning found — see the MaveDB page)

**Proposed:** `unknown` (confidence: none) · **Flags:** weak_control_separation,genomic_mapping_incomplete · **Review status:** **not yet reviewed**

### 4. BRCA1 — Enrich2 amino acid variant scores for BRCA1 Y2H

`urn:mavedb:00000003-b-2` · [MaveDB page](https://www.mavedb.org/score-sets/urn:mavedb:00000003-b-2) · PubMed:28784151 (10.1186/s13059-017-1272-5) · licence CC0

- **What was measured:** Amino acid variant scores for deep mutational scan of the BRCA1 RING domain using yeast two-hybrid calculated by Enrich2.
- **Assay:** —
- **Variants:** 12316 in MaveDB, 10402 with a numeric score, 4142 single missense
- **Replicates:** replicate columns: SE_Y2H_1_Rep1, score_Y2H_1_Rep1, SE_Y2H_1_Rep2, score_Y2H_1_Rep2, SE_Y2H_1_Rep3, score_Y2H_1_Rep3; error columns: SE, SE_Y2H_1_Rep1, SE_Y2H_1_Rep2, SE_Y2H_1_Rep3, SE_Y2H_2_Rep1, SE_Y2H_2_Rep2
- **Genome mapping by MaveDB:** incomplete

| Evidence | What it shows | Implied direction |
|---|---|---|
| A. Calibration | none | — |
| B. Controls | too few controls (nonsense 196, synonymous 1; need 10 each) | — |

C. What the authors say about the score:
   > (no sentence about score meaning found — see the MaveDB page)

**Proposed:** `unknown` (confidence: none) · **Flags:** no_nonsense_synonymous_controls,genomic_mapping_incomplete · **Review status:** **not yet reviewed**

### 5. BRCA1 — BRCA1 RING and BRCT domains depletion scores

`urn:mavedb:00000081-a-1` · [MaveDB page](https://www.mavedb.org/score-sets/urn:mavedb:00000081-a-1) · PubMed:30219179 (10.1016/j.ajhg.2018.07.016) · licence CC0

- **What was measured:** Depletion scores for BRCA1 RING and BRCT domain variants in homology-directed DNA repair function assay.
- **Assay:** —
- **Variants:** 1061 in MaveDB, 1061 with a numeric score, 1055 single missense
- **Replicates:** replicate columns: replicates
- **Genome mapping by MaveDB:** incomplete

| Evidence | What it shows | Implied direction |
|---|---|---|
| A. Calibration | none | — |
| B. Controls | too few controls (nonsense 4, synonymous 2; need 10 each) | — |

C. What the authors say about the score:
   > The depletion score reported here is the number of replicates where the variant was depleted relative to the corresponding control siRNA replicate.

**Proposed:** `unknown` (confidence: none) · **Flags:** no_nonsense_synonymous_controls,genomic_mapping_incomplete · **Review status:** **not yet reviewed**

### 6. BRCA1 — BRCA1 RING and BRCT domains mean HDR fluorescence score

`urn:mavedb:00000081-a-2` · [MaveDB page](https://www.mavedb.org/score-sets/urn:mavedb:00000081-a-2) · PubMed:30219179 (10.1016/j.ajhg.2018.07.016) · licence CC0

- **What was measured:** Mean HDR fluorescence score for BRCA1 RING and BRCT domain variants in homology-directed DNA repair functional assay
- **Assay:** —
- **Variants:** 2820 in MaveDB, 2820 with a numeric score, 2749 single missense
- **Replicates:** none reported in score columns
- **Genome mapping by MaveDB:** complete

| Evidence | What it shows | Implied direction |
|---|---|---|
| A. Calibration | none | — |
| B. Controls | too few controls (nonsense 71, synonymous 0; need 10 each) | — |

C. What the authors say about the score:
   > The frequency of each variant was determined by dividing the variant counts by the total variant count ($\sum{c_{+,i}}$) for each sample: $$f_{+,v} = \frac{c_{+,v}}{\sum{c_{+,i}}}$$ Selection log ratio is calculated as: $$r_v = \ln{\frac{f_{+,v}}{f_{-,v}}}$$ The final score is normalized by wildtype
   > This score set describes an alternative set of scores for this dataset that is based on the frequency of fluorescent cells rather than the number of replicates where a variant was depleted (as originally published).

**Proposed:** `unknown` (confidence: none) · **Flags:** no_nonsense_synonymous_controls · **Review status:** **not yet reviewed**

### 7. BRCA1 — BASE_ACMG quantities derived from the RING domain M2H functional assay

`urn:mavedb:00000093-a-1` · [MaveDB page](https://www.mavedb.org/score-sets/urn:mavedb:00000093-a-1) · doi:10.1101/092619v2 · licence CC0

- **What was measured:** BASE_ACMG are the quantitatively defined ACMG strength of evidence categories and ACMG-scaled points. This score set reports the BASE_ACMG values for variants used in the RING domain M2H functional assays across multiple experiments.
- **Assay:** —
- **Variants:** 853 in MaveDB, 853 with a numeric score, 659 single missense
- **Replicates:** none reported in score columns
- **Genome mapping by MaveDB:** incomplete

| Evidence | What it shows | Implied direction |
|---|---|---|
| A. Calibration | none | — |
| B. Controls | too few controls (nonsense 0, synonymous 0; need 10 each) | — |

C. What the authors say about the score:
   > (no sentence about score meaning found — see the MaveDB page)

**Proposed:** `unknown` (confidence: none) · **Flags:** no_nonsense_synonymous_controls,genomic_mapping_incomplete · **Review status:** **not yet reviewed**

### 8. BRCA1 — Scores from growth assay of BRCA1 variants

`urn:mavedb:00001222-a-2` · [MaveDB page](https://www.mavedb.org/score-sets/urn:mavedb:00001222-a-2) · PubMed:35196514 (10.1016/j.ajhg.2022.01.019) · licence CC0

- **What was measured:** Growth assay separating functionally normal and functionally abnormal BRCA1 variants using cisplatin resistance as a selection
- **Assay:** Cell fitness · Loss of function · Immortalized human cells
- **Variants:** 1427 in MaveDB, 1427 with a numeric score, 1359 single missense
- **Replicates:** replicate columns: rep1_score, rep2_score, rep3_score, rep4_score, control_rep1_score, control_rep2_score; error columns: std, var, control_std, control_var
- **Genome mapping by MaveDB:** complete

| Evidence | What it shows | Implied direction |
|---|---|---|
| A. Calibration | abnormal range(s) below normal range(s) — IGVF Coding Variant Focus Group -- Controls: All Variants; Investigator-provided functional classes; IGVF Coding Variant Focus Group -- Controls: Missense Variants Only; ExCALIBR calibration [research use only]; ExCALIBR calibration (ClinVar 2018) [research use only] | lower_is_more_damaging |
| B. Controls | too few controls (nonsense 68, synonymous 0; need 10 each) | — |

C. What the authors say about the score:
   > (no sentence about score meaning found — see the MaveDB page)

**Proposed:** `lower_is_more_damaging` (confidence: medium) · **Flags:** no_nonsense_synonymous_controls · **Review status:** **not yet reviewed**

### 9. BRCA1 — Scores from multiplexed functional assay of BRCA1 variants

`urn:mavedb:00001222-b-2` · [MaveDB page](https://www.mavedb.org/score-sets/urn:mavedb:00001222-b-2) · PubMed:35196514 (10.1016/j.ajhg.2022.01.019) · licence CC0

- **What was measured:** Multiplexed assay of BRCA1 variants measuring homology directed repair activity
- **Assay:** Reporter · Loss of function · Immortalized human cells
- **Variants:** 2271 in MaveDB, 2271 with a numeric score, 2150 single missense
- **Replicates:** replicate columns: rep1_score, rep2_score, rep3_score, rep4_score, control_rep1_score, control_rep2_score; error columns: std, var, control_std, control_var
- **Genome mapping by MaveDB:** complete

| Evidence | What it shows | Implied direction |
|---|---|---|
| A. Calibration | abnormal range(s) below normal range(s) — Investigator-provided functional classes; IGVF Coding Variant Focus Group -- Controls: Missense Variants Only; IGVF Coding Variant Focus Group -- Controls: All Variants; ExCALIBR calibration [research use only]; ExCALIBR calibration (ClinVar 2018) [research use only] | lower_is_more_damaging |
| B. Controls | too few controls (nonsense 121, synonymous 0; need 10 each) | — |

C. What the authors say about the score:
   > (no sentence about score meaning found — see the MaveDB page)

**Proposed:** `lower_is_more_damaging` (confidence: medium) · **Flags:** no_nonsense_synonymous_controls · **Review status:** **not yet reviewed**

### 10. TP53 — p53 variant effect measured by cell growth

`urn:mavedb:00000059-a-1` · [MaveDB page](https://www.mavedb.org/score-sets/urn:mavedb:00000059-a-1) · PubMed:29979965 (10.1016/j.molcel.2018.06.012) · licence CC0

- **What was measured:** The impact of p53 variant is represented by the frequency change of cells harboring that mutation over several time point.
- **Assay:** —
- **Variants:** 9273 in MaveDB, 9273 with a numeric score, 4069 single missense
- **Replicates:** none reported in score columns
- **Genome mapping by MaveDB:** incomplete

| Evidence | What it shows | Implied direction |
|---|---|---|
| A. Calibration | none | — |
| B. Controls | nonsense (n=344, median 0.414) vs synonymous (n=572, median -2.42); AUC 1.00 | higher_is_more_damaging |

C. What the authors say about the score:
   > (no sentence about score meaning found — see the MaveDB page)

**Proposed:** `higher_is_more_damaging` (confidence: medium) · **Flags:** genomic_mapping_incomplete · **Review status:** **not yet reviewed**

### 11. TP53 — Mutated p53 paired with wildtype under nutlin-3

`urn:mavedb:00000068-a-1` · [MaveDB page](https://www.mavedb.org/score-sets/urn:mavedb:00000068-a-1) · PubMed:30224644 (10.1038/s41588-018-0204-y) · licence CC0

- **What was measured:** TP53 saturation mutagenesis screens with TP53 wild-type using nutlin-3 which interrupts a negative regulator of p53
- **Assay:** Cell fitness · Dominant-negative effect · Immortalized human cells
- **Variants:** 8274 in MaveDB, 8258 with a numeric score, 7487 single missense
- **Replicates:** none reported in score columns
- **Genome mapping by MaveDB:** incomplete

| Evidence | What it shows | Implied direction |
|---|---|---|
| A. Calibration | abnormal range(s) above normal range(s) — ExCALIBR calibration [research use only]; ExCALIBR calibration (ClinVar 2018) [research use only]; Investigator-provided functional classes | higher_is_more_damaging |
| B. Controls | nonsense (n=393, median -0.264) vs synonymous (n=378, median 0.0495); AUC 0.22 | lower_is_more_damaging |

C. What the authors say about the score:
   > (no sentence about score meaning found — see the MaveDB page)

**Proposed:** `conflict` (confidence: none) · **Flags:** genomic_mapping_incomplete · **Review status:** **not yet reviewed**

### 12. TP53 — Mutated p53 with nutlin-3

`urn:mavedb:00000068-b-1` · [MaveDB page](https://www.mavedb.org/score-sets/urn:mavedb:00000068-b-1) · PubMed:30224644 (10.1038/s41588-018-0204-y) · licence CC0

- **What was measured:** TP53 saturation mutagenesis screens with TP53 null cell line using nutlin-3 which interrupts a negative regulator of p53
- **Assay:** Cell fitness · Loss of function · Immortalized human cells
- **Variants:** 8274 in MaveDB, 8258 with a numeric score, 7487 single missense
- **Replicates:** none reported in score columns
- **Genome mapping by MaveDB:** incomplete

| Evidence | What it shows | Implied direction |
|---|---|---|
| A. Calibration | abnormal range(s) above normal range(s) — ExCALIBR calibration [research use only]; ExCALIBR calibration (ClinVar 2018) [research use only]; Investigator-provided functional classes | higher_is_more_damaging |
| B. Controls | nonsense (n=393, median 1.16) vs synonymous (n=378, median -0.254); AUC 0.84 | higher_is_more_damaging |

C. What the authors say about the score:
   > (no sentence about score meaning found — see the MaveDB page)

**Proposed:** `higher_is_more_damaging` (confidence: high) · **Flags:** genomic_mapping_incomplete · **Review status:** **not yet reviewed**

### 13. TP53 — Mutated p53 with etoposide

`urn:mavedb:00000068-c-1` · [MaveDB page](https://www.mavedb.org/score-sets/urn:mavedb:00000068-c-1) · PubMed:30224644 (10.1038/s41588-018-0204-y) · licence CC0

- **What was measured:** TP53 saturation mutagenesis screens with TP53 null cell line using etoposide which is a DNA double-strand break-inducing agent.
- **Assay:** Cell fitness · Loss of function · Immortalized human cells
- **Variants:** 8274 in MaveDB, 8258 with a numeric score, 7487 single missense
- **Replicates:** none reported in score columns
- **Genome mapping by MaveDB:** incomplete

| Evidence | What it shows | Implied direction |
|---|---|---|
| A. Calibration | abnormal range(s) below normal range(s) — ExCALIBR calibration [research use only]; ExCALIBR calibration (ClinVar 2018) [research use only]; Investigator-provided functional classes | lower_is_more_damaging |
| B. Controls | nonsense (n=393, median -1.52) vs synonymous (n=378, median 0.365); AUC 0.10 | lower_is_more_damaging |

C. What the authors say about the score:
   > (no sentence about score meaning found — see the MaveDB page)

**Proposed:** `lower_is_more_damaging` (confidence: high) · **Flags:** genomic_mapping_incomplete · **Review status:** **not yet reviewed**

### 14. TP53 — Deep mutational scan of human TP53 in HCT116 colorectal cancer cells covering Exons 5, 6, 7 and 8

`urn:mavedb:00001213-a-1` · [MaveDB page](https://www.mavedb.org/score-sets/urn:mavedb:00001213-a-1) · no publication listed · licence CC0

- **What was measured:** Deep mutational scan of 9225 mutations in human TP53 covering Exons 5, 6, 7 and 8 plus flanking 13nt intronic regions in HCT116 colorectal cancer cells. Introduced where SNVs, DNVs, TNVs, insertions of size 1 and deletions up to size 3.
- **Assay:** —
- **Variants:** 8052 in MaveDB, 8052 with a numeric score, 0 single missense
- **Replicates:** error columns: SE_RFS_enrich2
- **Genome mapping by MaveDB:** complete

| Evidence | What it shows | Implied direction |
|---|---|---|
| A. Calibration | none | — |
| B. Controls | no protein-level HGVS; controls not classifiable | — |

C. What the authors say about the score:
   > (no sentence about score meaning found — see the MaveDB page)

**Proposed:** `unknown` (confidence: none) · **Flags:** no_nonsense_synonymous_controls,fewer_than_20_scored_missense · **Review status:** **not yet reviewed**

### 15. TP53 — Scores from arrayed yeast-based assay of TP53 p21WAF1 promoter

`urn:mavedb:00001234-a-1` · [MaveDB page](https://www.mavedb.org/score-sets/urn:mavedb:00001234-a-1) · PubMed:12826609 (10.1073/pnas.1431692100) · licence CC0

- **What was measured:** Arrayed yeast-based assay of TP53 transcriptional activity of the p21WAF1 promoter
- **Assay:** Reporter · Gain or loss of function · Yeast
- **Variants:** 2314 in MaveDB, 2314 with a numeric score, 2314 single missense
- **Replicates:** none reported in score columns
- **Genome mapping by MaveDB:** complete

| Evidence | What it shows | Implied direction |
|---|---|---|
| A. Calibration | abnormal range(s) below normal range(s) — ExCALIBR calibration [research use only]; ExCALIBR calibration (ClinVar 2018) [research use only]; Investigator-provided functional classes | lower_is_more_damaging |
| B. Controls | too few controls (nonsense 0, synonymous 0; need 10 each) | — |

C. What the authors say about the score:
   > (no sentence about score meaning found — see the MaveDB page)

**Proposed:** `lower_is_more_damaging` (confidence: medium) · **Flags:** no_nonsense_synonymous_controls · **Review status:** **not yet reviewed**

### 16. TP53 — Scores from arrayed yeast-based assay of TP53 MDM2 promoter

`urn:mavedb:00001234-b-1` · [MaveDB page](https://www.mavedb.org/score-sets/urn:mavedb:00001234-b-1) · PubMed:12826609 (10.1073/pnas.1431692100) · licence CC0

- **What was measured:** Arrayed yeast-based assay of TP53 transcriptional activity of the MDM2 promoter
- **Assay:** Reporter · Gain or loss of function · Yeast
- **Variants:** 2314 in MaveDB, 2314 with a numeric score, 2314 single missense
- **Replicates:** none reported in score columns
- **Genome mapping by MaveDB:** complete

| Evidence | What it shows | Implied direction |
|---|---|---|
| A. Calibration | abnormal range(s) below normal range(s) — ExCALIBR calibration [research use only]; ExCALIBR calibration (ClinVar 2018) [research use only]; Investigator-provided functional classes | lower_is_more_damaging |
| B. Controls | too few controls (nonsense 0, synonymous 0; need 10 each) | — |

C. What the authors say about the score:
   > (no sentence about score meaning found — see the MaveDB page)

**Proposed:** `lower_is_more_damaging` (confidence: medium) · **Flags:** no_nonsense_synonymous_controls · **Review status:** **not yet reviewed**

### 17. TP53 — Scores from arrayed yeast-based assay of TP53 BAX promoter

`urn:mavedb:00001234-c-1` · [MaveDB page](https://www.mavedb.org/score-sets/urn:mavedb:00001234-c-1) · PubMed:12826609 (10.1073/pnas.1431692100) · licence CC0

- **What was measured:** Arrayed yeast-based assay of TP53 transcriptional activity of the BAX promoter
- **Assay:** Reporter · Gain or loss of function · Yeast
- **Variants:** 2314 in MaveDB, 2314 with a numeric score, 2314 single missense
- **Replicates:** none reported in score columns
- **Genome mapping by MaveDB:** complete

| Evidence | What it shows | Implied direction |
|---|---|---|
| A. Calibration | abnormal range(s) below normal range(s) — ExCALIBR calibration [research use only]; ExCALIBR calibration (ClinVar 2018) [research use only]; Investigator-provided functional classes | lower_is_more_damaging |
| B. Controls | too few controls (nonsense 0, synonymous 0; need 10 each) | — |

C. What the authors say about the score:
   > (no sentence about score meaning found — see the MaveDB page)

**Proposed:** `lower_is_more_damaging` (confidence: medium) · **Flags:** no_nonsense_synonymous_controls · **Review status:** **not yet reviewed**

### 18. TP53 — Arrayed yeast-based assay of TP53 h1433snWT promoter

`urn:mavedb:00001234-d-1` · [MaveDB page](https://www.mavedb.org/score-sets/urn:mavedb:00001234-d-1) · PubMed:12826609 (10.1073/pnas.1431692100) · licence CC0

- **What was measured:** Arrayed yeast-based assay of TP53 transcriptional activity of the h1433snWT promoter
- **Assay:** Reporter · Gain or loss of function · Yeast
- **Variants:** 2314 in MaveDB, 2314 with a numeric score, 2314 single missense
- **Replicates:** none reported in score columns
- **Genome mapping by MaveDB:** complete

| Evidence | What it shows | Implied direction |
|---|---|---|
| A. Calibration | abnormal range(s) below normal range(s) — ExCALIBR calibration [research use only]; ExCALIBR calibration (ClinVar 2018) [research use only]; Investigator-provided functional classes | lower_is_more_damaging |
| B. Controls | too few controls (nonsense 0, synonymous 0; need 10 each) | — |

C. What the authors say about the score:
   > (no sentence about score meaning found — see the MaveDB page)

**Proposed:** `lower_is_more_damaging` (confidence: medium) · **Flags:** no_nonsense_synonymous_controls · **Review status:** **not yet reviewed**

### 19. TP53 — Scores from arrayed yeast-based assay of TP53 p53AIP1 promoter

`urn:mavedb:00001234-e-1` · [MaveDB page](https://www.mavedb.org/score-sets/urn:mavedb:00001234-e-1) · PubMed:12826609 (10.1073/pnas.1431692100) · licence CC0

- **What was measured:** Arrayed yeast-based assay of TP53 transcriptional activity of the p53AIP1 promoter
- **Assay:** Reporter · Gain or loss of function · Yeast
- **Variants:** 2314 in MaveDB, 2314 with a numeric score, 2314 single missense
- **Replicates:** none reported in score columns
- **Genome mapping by MaveDB:** complete

| Evidence | What it shows | Implied direction |
|---|---|---|
| A. Calibration | abnormal range(s) below normal range(s) — ExCALIBR calibration [research use only]; ExCALIBR calibration (ClinVar 2018) [research use only]; Investigator-provided functional classes | lower_is_more_damaging |
| B. Controls | too few controls (nonsense 0, synonymous 0; need 10 each) | — |

C. What the authors say about the score:
   > (no sentence about score meaning found — see the MaveDB page)

**Proposed:** `lower_is_more_damaging` (confidence: medium) · **Flags:** no_nonsense_synonymous_controls · **Review status:** **not yet reviewed**

### 20. TP53 — Scores from arrayed yeast-based assay of TP53 GADD45 promoter

`urn:mavedb:00001234-f-1` · [MaveDB page](https://www.mavedb.org/score-sets/urn:mavedb:00001234-f-1) · PubMed:12826609 (10.1073/pnas.1431692100) · licence CC0

- **What was measured:** Arrayed yeast-based assay of TP53 transcriptional activity of the GADD45 promoter
- **Assay:** Reporter · Gain or loss of function · Yeast
- **Variants:** 2314 in MaveDB, 2314 with a numeric score, 2314 single missense
- **Replicates:** none reported in score columns
- **Genome mapping by MaveDB:** complete

| Evidence | What it shows | Implied direction |
|---|---|---|
| A. Calibration | abnormal range(s) below normal range(s) — ExCALIBR calibration [research use only]; ExCALIBR calibration (ClinVar 2018) [research use only]; Investigator-provided functional classes | lower_is_more_damaging |
| B. Controls | too few controls (nonsense 0, synonymous 0; need 10 each) | — |

C. What the authors say about the score:
   > (no sentence about score meaning found — see the MaveDB page)

**Proposed:** `lower_is_more_damaging` (confidence: medium) · **Flags:** no_nonsense_synonymous_controls · **Review status:** **not yet reviewed**

### 21. TP53 — Scores from arrayed yeast-based assay of TP53 Noxa promoter

`urn:mavedb:00001234-g-1` · [MaveDB page](https://www.mavedb.org/score-sets/urn:mavedb:00001234-g-1) · PubMed:12826609 (10.1073/pnas.1431692100) · licence CC0

- **What was measured:** Arrayed yeast-based assay of TP53 transcriptional activity of the Noxa promoter
- **Assay:** Reporter · Gain or loss of function · Yeast
- **Variants:** 2314 in MaveDB, 2314 with a numeric score, 2314 single missense
- **Replicates:** none reported in score columns
- **Genome mapping by MaveDB:** complete

| Evidence | What it shows | Implied direction |
|---|---|---|
| A. Calibration | abnormal range(s) below normal range(s) — ExCALIBR calibration [research use only]; ExCALIBR calibration (ClinVar 2018) [research use only]; Investigator-provided functional classes | lower_is_more_damaging |
| B. Controls | too few controls (nonsense 0, synonymous 0; need 10 each) | — |

C. What the authors say about the score:
   > (no sentence about score meaning found — see the MaveDB page)

**Proposed:** `lower_is_more_damaging` (confidence: medium) · **Flags:** no_nonsense_synonymous_controls · **Review status:** **not yet reviewed**

### 22. TP53 — Scores from arrayed yeast-based assay of TP53 p53R2 promoter

`urn:mavedb:00001234-h-1` · [MaveDB page](https://www.mavedb.org/score-sets/urn:mavedb:00001234-h-1) · PubMed:12826609 (10.1073/pnas.1431692100) · licence CC0

- **What was measured:** Arrayed yeast-based assay of TP53 transcriptional activity of the p53R2 promoter
- **Assay:** Reporter · Gain or loss of function · Yeast
- **Variants:** 2314 in MaveDB, 2314 with a numeric score, 2314 single missense
- **Replicates:** none reported in score columns
- **Genome mapping by MaveDB:** complete

| Evidence | What it shows | Implied direction |
|---|---|---|
| A. Calibration | abnormal range(s) below normal range(s) — ExCALIBR calibration [research use only]; ExCALIBR calibration (ClinVar 2018) [research use only]; Investigator-provided functional classes | lower_is_more_damaging |
| B. Controls | too few controls (nonsense 0, synonymous 0; need 10 each) | — |

C. What the authors say about the score:
   > (no sentence about score meaning found — see the MaveDB page)

**Proposed:** `lower_is_more_damaging` (confidence: medium) · **Flags:** no_nonsense_synonymous_controls · **Review status:** **not yet reviewed**

### 23. TP53 — Scores from Nutlin-3a treatment of TP53 variant-expressing K562 reporter cells

`urn:mavedb:00001235-a-1` · [MaveDB page](https://www.mavedb.org/score-sets/urn:mavedb:00001235-a-1) · PubMed:31395785 (10.1126/science.aax3649) · licence CC0

- **What was measured:** Nutlin-3a treatment of TP53 variant-expressing K562 reporter cells was used to enrich for dominant-negative variants
- **Assay:** Reporter · Dominant-negative effect · Immortalized human cells
- **Variants:** 7893 in MaveDB, 7893 with a numeric score, 7466 single missense
- **Replicates:** none reported in score columns
- **Genome mapping by MaveDB:** incomplete

| Evidence | What it shows | Implied direction |
|---|---|---|
| A. Calibration | abnormal range(s) above normal range(s) — ExCALIBR calibration [research use only]; ExCALIBR calibration (ClinVar 2018) [research use only]; Investigator-provided functional classes | higher_is_more_damaging |
| B. Controls | nonsense (n=391, median -0.3) vs synonymous (n=36, median -0.0493); AUC 0.22 | lower_is_more_damaging |

C. What the authors say about the score:
   > (no sentence about score meaning found — see the MaveDB page)

**Proposed:** `conflict` (confidence: none) · **Flags:** genomic_mapping_incomplete · **Review status:** **not yet reviewed**


## When you are done

Tell Arman. The next pipeline stage (DMS harmonization) refuses to run until every in-scope
assay has a `reviewer_decision`, and every included assay has a `reviewer_direction`.

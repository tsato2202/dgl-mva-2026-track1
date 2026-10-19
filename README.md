# DGL-COMP-01 — MVA Hackathon 2026 / Track 1

**STATUS: SYNTHETIC PREFLIGHT ONLY — NO CLINICAL PREDICTIONS, NO REAL SUBMISSIONS.**

This repository contains reproducibility scaffolding for a Track 1 variant-prioritization competition submission. It does **not** contain a patient VCF, patient-derived genotype tables, clinical notes, access tokens, or clinically validated causal hypotheses.

## Evidence completed on 2026-10-19 (Google Colab CPU)

| Stage | Input | Checks completed | Status |
|---|---|---|---|
| E00 | synthetic VCF | 7 records; 5 QC pass; 2 excluded; 12-column CSV | PASS (synthetic only) |
| E01 | synthetic VCF | 12 input records; 13 expanded ALT; QC 7 pass / 1 review / 5 exclude; 1 synthetic trans pair, 1 cis exclusion | PASS (synthetic only) |

E00/E01 are *software checks*, not clinical validation and not official Track 1 submissions. Official submission quota last observed: 0 of 6 used.

## Contents

- `notebooks/`: synthetic-only Colab exercises. The notebooks use invented examples. **Do not upload their generated example CSV to the competition.**
- `scripts/validate_track1_official.py`: CSV syntax/shape preflight, **not** a variant classifier or official evaluator.
- `templates/track1_submission_BLANK.csv`: header-only template; **not valid for submission without supported candidates**.
- `reports/E00_Methods_Report_DRAFT.md`: method disclosure draft to be completed with evidence; not yet ready for submission.
- `docs/GITHUB_SETUP_JA.md`: steps to create a private repository and upload this scaffold.
- `docs/DATA_HANDLING.md`: controlled-data boundaries and deletion reminders.

## Reproducibility

Prerequisites: Python 3.11+ and Google Colab CPU for notebooks (standard-library only). From the repository root:

```bash
python scripts/validate_track1_official.py --self-test
```

Expected: `SELF-TEST PASS (synthetic data only)`.

The controlled case must **only** be processed after the data-use terms for the specific compute service are verified. Use a separate, access-controlled environment with no patient data in this repository. Never publish any clinical free text or genome-scale variant table.

## Competition

Official Space: https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026

Official dataset: https://huggingface.co/datasets/SageBio/mva-hackathon-2026-data

Deadline: 2026-10-24 23:59 UTC (2026-10-25 08:59 JST).

Official Track 1 form requires a GitHub URL, a ranked predictions CSV (up to ten rows, GRCh38 coordinates) and a methods report. All final submissions require manual review and authorization. Public availability of submissions and final repository must follow organizer rules.

**Scientific caveat:** The synthetic notebook `GENE`, `CSQ`, `AF`, `HPO` fields are exercise-specific and not expected to exist in a real VCF. The synthetic ranking number is not a calibrated causal probability. Phased `trans` evidence requires verification in real data and may remain unknown.

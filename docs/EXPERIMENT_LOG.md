# Experiment log — DGL-COMP-01

## E00 — Synthetic baseline / 2026-10-19

- environment: Google Colab, CPU
- 7 input variants; QC pass 5, exclude 2; 5 ranked hypotheses
- official CSV 12 columns: PASS
- outcome: local synthetic only; no official submission
- limitations: artificial annotations and ranking, no real phenotype concordance, no clinically confirmed variant

## E01 — Synthetic edge cases / 2026-10-19

- environment: Google Colab, CPU
- input records 12; multi-allelic expansion 13
- QC pass 7 / review 1 / exclude 5
- one trans compound-heterozygous synthetic pair, one cis pair excluded
- ranked hypotheses 8; official CSV header PASS
- outcome: local synthetic only; no official submission
- limitations: phase and quality field interpretations tested in controlled synthetic examples, no real patient analysis

## E02 — GitHub reproducibility scaffold

- status: package prepared, repository not yet created
- goal: verified publishable methods/code separated from controlled genomic data

## S01 — Official first submission

- status: NOT STARTED, 0/6 slots used last observed
- preconditions: authorized processing environment; evidence-backed patient candidates; tested CSV; methods report; GitHub URL; explicit user approval

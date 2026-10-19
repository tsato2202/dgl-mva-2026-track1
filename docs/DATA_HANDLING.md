# Controlled patient genomic data boundary

**Never store in this repository (private repository included):**
- VCF/BCF/BAM/CRAM/FASTQ/TBI and any copies or slices
- filtered/annotated genome-scale genotype tables
- raw patient clinical narratives and protected identifiers
- token, credentials, environment dumps, code-cell output containing sensitive data
- notebook state holding per-variant genotypes, cache snapshots and model weights fitted on raw case genomes

`gitignore` is a backup barrier, **not authorization** to add patient-derived data. Manual review mandatory.

Only code, synthetic examples, licensed public references, and permitted submission-scale findings may be committed. Keep any actual identifiers and sensitive data in authorized processing storage, not in ChatGPT or a public repository.

## Third-party processing decision

The Competition organizer has not yet answered the participant-specific question about consumer Google AI Pro Colab: https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026/discussions/30

Do not process gated genome data in Google Colab until applicable service terms/permissions have been verified. Google AI Pro entitlement != enterprise data processing agreement.

See organizer's Processor-vs-Recipient guidance: https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026/discussions/2

## Deletion

Follow official deletion window, including cloud runtime disks, caches, intermediate genotype-level tables, raw downloads, model state and user-controlled copies. Keep an internal inventory and send required deletion attestation. Code/method descriptions and permitted candidate-scale results may remain subject to official rules.

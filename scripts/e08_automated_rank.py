#!/usr/bin/env python3
"""E08 local ranking of compound-heterozygote *hypotheses* for MVA Track 1.

Input:  patient_private/E04I/E04I_ALL_CANDIDATES_PRIVATE_DO_NOT_UPLOAD.csv
Output: patient_private/E08/ ONLY. No network; no raw patient details on stdout.

Uses explicit evidence REVIEW weights. Does NOT estimate posterior causal probabilities,
read genome-wide VCF, prove phasing, annotate transcript consequence, or query gnomAD.
"""
from __future__ import annotations
import argparse
import csv
import html
import itertools
import json
import math
import os
import re
import sys
import tempfile
from collections import Counter, defaultdict
from pathlib import Path

SOURCE='E04I_ALL_CANDIDATES_PRIVATE_DO_NOT_UPLOAD.csv'
OUTPUT='E08_TRACK1_AUTO_RANKED_DRAFT_REVIEW_BEFORE_UPLOAD.csv'
AUDIT='E08_PRIVATE_EVIDENCE_AUDIT_DO_NOT_UPLOAD.csv'
PREVIEW='E08_PRIVATE_REVIEW_DO_NOT_UPLOAD.html'
REPORT='E08_METHODS_DRAFT_REVIEW_BEFORE_UPLOAD.md'
SUMMARY='E08_AGGREGATE_ONLY.json'
OFFICIAL=['proband_id','chrom_1','pos_1','ref_1','alt_1','chrom_2','pos_2','ref_2','alt_2','epcr','finding_type','notes']
NEED={'gene','chrom','pos','ref','alt','GT','DP','GQ','FILTER','public_HPO_exact_n',
      'public_disease_max_HPO_n','genecode_location_proxy','clinvar_status'}
EVIDENCE_FIELDS=['candidate_rank','selection','hypothesis_type','gene','chrom_1','pos_1','ref_1','alt_1',
    'chrom_2','pos_2','ref_2','alt_2','evidence_score_NOT_PROBABILITY','HPO_disease_overlap','HPO_gene_overlap',
    'HPO_gene_IDF_proxy','variant1_GT','variant2_GT','variant1_DP','variant2_DP','variant1_GQ','variant2_GQ',
    'variant1_position_proxy','variant2_position_proxy','variant1_ClinVar','variant2_ClinVar',
    'variant1_prior_tier','variant2_prior_tier','variant1_ClinVar_condition','variant2_ClinVar_condition',
    'phasing','population_AF','coding_consequence','reason_flags']
PHASE='UNVALIDATED_NOT_PROVEN_TRANS'
CODE='TRANSCRIPT_SPECIFIC_CONSEQUENCE_NOT_EVALUATED'
AF='GNOMAD_AF_NOT_EVALUATED'
REF_BASE=re.compile(r'^[ACGTN]{1,1000}$')


def check_private(path:Path):
    path=path.resolve()
    if any((q/'.git').exists() for q in (path,*path.parents)):
        raise ValueError('Blocked: private analysis data is inside a git repository')
    return path


def fnum(s,default=0.0):
    try:
        n=float(s)
        return n if math.isfinite(n) else default
    except (ValueError,TypeError): return default


def canon_chrom(s):
    s=str(s).strip()
    if s.lower().startswith('chr'):s=s[3:]
    if s in {str(i) for i in range(1,23)} | {'X','Y'}:return 'chr'+s
    if s in {'M','MT'}:return 'chrM'
    raise ValueError('Unsupported chromosome')


def loc(row):
    pos=int(row['pos'])
    ref=str(row['ref']).upper().strip();alt=str(row['alt']).upper().strip()
    if pos<1 or not REF_BASE.fullmatch(ref) or not REF_BASE.fullmatch(alt) or ref==alt:
        raise ValueError('Invalid variant locus')
    return canon_chrom(row['chrom']),pos,ref,alt


def qualified(r):
    if r.get('GT') not in {'0/1','1/0','0|1','1|0'}:return False
    if r.get('FILTER')!='PASS':return False
    if fnum(r.get('DP'),-1)<10 or fnum(r.get('GQ'),-1)<20:return False
    if r.get('clinvar_status')=='BENIGN_LABEL_CHECK_CONDITION':return False
    return True


def is_prior(r):
    return str(r.get('previous_review_tier','')).startswith(('P1_','P2_'))


def positional(r):
    return {
        'CDS_POSITIONAL_OVERLAP':4.0,
        'NONCDS_EXON_OVERLAP':1.5,
        'WITHIN_2BP_EXON_EDGE_PROXY':1.5,
    }.get(r.get('genecode_location_proxy'),0.0)


def variant_score(r):
    # Neither missense nor LoF is inferred from CDS overlap / indel length.
    points=positional(r)
    if r.get('clinvar_status')=='PATHOGENIC_LABEL_REVIEW_CONDITION':points+=2.0  # condition UNVERIFIED
    if is_prior(r):points+=1.0   # historical review, not independent clinical evidence
    if positional(r)>0 and len(r['ref'])!=len(r['alt']):points+=0.5  # length change only
    return points


def gene_score(a,b):
    disease=min(max(fnum(a.get('public_disease_max_HPO_n')),
                    fnum(b.get('public_disease_max_HPO_n'))),8.0)
    gene=min(max(fnum(a.get('public_HPO_exact_n')),
                 fnum(b.get('public_HPO_exact_n'))),8.0)
    idf=min(max(fnum(a.get('public_HPO_IDF_score')),
                fnum(b.get('public_HPO_IDF_score'))),40.0)
    # Separate overlapping, correlated HPO terms; heuristic capped signals.
    return round(2.0*disease+0.75*gene+0.05*idf,4)


def variant_key(v):return loc(v)


def read_input(root):
    inp=check_private(root/'patient_private'/'E04I'/SOURCE)
    if not inp.is_file():raise ValueError('E04I input missing')
    with inp.open(encoding='utf-8-sig',newline='') as f:
        rd=csv.DictReader(f)
        if not NEED.issubset(set(rd.fieldnames or [])):
            raise ValueError('E04I input schema mismatch')
        rows=list(rd)
    if not rows:raise ValueError('E04I input empty')
    if len(rows)>1000000:raise ValueError('Unexpected input size')
    return rows


def hypothesis(a,b,typ='pair'):
    v1=loc(a);v2=loc(b) if b else None
    if b and (a['gene']!=b['gene'] or v1[0]!=v2[0] or v1==v2 or v1[1]==v2[1]):
        raise ValueError('Invalid same-gene pair')
    # Canonicalize pair allele ordering to make output reproducible.
    if b and v2<v1:
        a,b=b,a
        v1,v2=v2,v1
    s=gene_score(a,b or a)+variant_score(a)
    if b:
        s+=variant_score(b)+0.25*min(positional(a),positional(b))
    else:
        s-=2.0  # single does not encode known two-variant key
    score=round(s,4)
    return {'gene':a['gene'],'a':a,'b':b,'v1':v1,'v2':v2,'type':typ,
            'score':score,'gscore':gene_score(a,b or a),
            'pos_score_a':positional(a),'pos_score_b':positional(b) if b else 0.0,
            'prior_present':is_prior(a) or (b is not None and is_prior(b))}


def hyp_key(h):
    return (h['v1'],h['v2']) if h['v2'] is None else tuple(sorted((h['v1'],h['v2'])))


def sort_hyp(h):
    # HPO, annotation, review status, and location give deterministic tie-breaking.
    return (-h['score'],-h['gscore'],-h['pos_score_a']-h['pos_score_b'],
            h['gene'],h['v1'],h['v2'] or ('',0,'',''))


def candidate_pool(rows,top_per_gene=4):
    by_gene=defaultdict(dict)
    prior_seen=set()
    malformed=0
    eligible=0
    functional_eligible=0
    functional_eligible_genes=set()
    excluded_no_priority_evidence=0
    eligible_genes=set()
    for r in rows:
        try:k=loc(r)
        except (ValueError,KeyError,TypeError):
            malformed+=1
            continue
        if is_prior(r):prior_seen.add((r['gene'],k))
        if not qualified(r):continue
        eligible+=1
        eligible_genes.add(r['gene'])
        if not (positional(r)>0 or is_prior(r) or r.get('clinvar_status')=='PATHOGENIC_LABEL_REVIEW_CONDITION'):
            excluded_no_priority_evidence+=1
            continue
        functional_eligible+=1
        functional_eligible_genes.add(r['gene'])
        best=by_gene[r['gene']].get(k)
        if best is None or variant_score(r)>variant_score(best):by_gene[r['gene']][k]=r
    hyps=[];pairs_by_gene=Counter(); eligible_prior=set()
    for gene,loci in by_gene.items():
        candidates=sorted(loci.values(),key=lambda r:(-variant_score(r),loc(r)))
        # Ensure a previously reviewed P1/P2 allele can compete, but do not force it into selection.
        selected=candidates[:top_per_gene]
        extra=[r for r in candidates if is_prior(r) and loc(r) not in {loc(q) for q in selected}]
        selected+=extra
        eligible_prior.update((gene,loc(r)) for r in candidates if is_prior(r))
        for a,b in itertools.combinations(selected,2):
            if loc(a)[1]==loc(b)[1]:continue
            h=hypothesis(a,b)
            hyps.append(h)
            pairs_by_gene[gene]+=1
    hyps.sort(key=sort_hyp)
    # A known 2-variant answer makes paired hypotheses preferable, but include
    # singles only when there aren't enough distinct pair hypotheses.
    singles=[]
    for gene,loci in by_gene.items():
        for r in loci.values():
            if positional(r)>0 or is_prior(r) or r.get('clinvar_status')=='PATHOGENIC_LABEL_REVIEW_CONDITION':
                singles.append(hypothesis(r,None,'single'))
    singles.sort(key=sort_hyp)
    counts={'input_rows':len(rows),'malformed_rows_excluded':malformed,
            'quality_eligible_allele_rows':eligible,'quality_eligible_genes':len(eligible_genes),
            'annotation_or_prior_evidence_eligible_rows':functional_eligible,
            'annotation_or_prior_evidence_eligible_genes':len(functional_eligible_genes),
            'quality_eligible_but_no_prior_position_or_clinvar_signal_excluded':excluded_no_priority_evidence,
            'pair_hypotheses_scored':len(hyps),'genes_with_pair_hypotheses':len(pairs_by_gene),
            'single_hypotheses_scored':len(singles),'prior_P1_P2_distinct_alleles_seen':len(prior_seen),
            'prior_P1_P2_distinct_alleles_quality_eligible':len(eligible_prior)}
    return hyps,singles,counts,prior_seen


def rank_select(pairs,singles,max_rows=10,max_per_gene=2):
    chosen=[];per_gene=Counter();seen=set()
    for pool in (pairs,singles):
        for h in pool:
            if len(chosen)>=max_rows:break
            k=hyp_key(h)
            if k in seen or per_gene[h['gene']]>=max_per_gene:continue
            seen.add(k);chosen.append(h);per_gene[h['gene']]+=1
        if len(chosen)>=max_rows:break
    if not chosen:raise ValueError('No eligible ranked candidates')
    return chosen


def epcr_placeholder(i,max_rows):
    # Fixed rank-to-number mapping, NOT posterior probability calibration.
    # 10 rows: 0.91,0.82,...,0.10; smaller lists remain strictly ordered.
    return f'{(max_rows-i)/(max_rows+1):.6f}'


def official_row(h,rank,max_rows):
    a=h['v1'];b=h['v2'] or ('','','','')
    notes=(f"{h['gene']}; automated uncalibrated evidence rank {rank}; "
           +('two heterozygous loci; trans not validated; ' if h['v2'] else 'single variant; ')
           +'population AF, transcript consequence and condition concordance not confirmed')
    return {'proband_id':'PROBAND01','chrom_1':a[0],'pos_1':str(a[1]),'ref_1':a[2],'alt_1':a[3],
            'chrom_2':b[0],'pos_2':str(b[1]) if b[1]!='' else '',
            'ref_2':b[2],'alt_2':b[3],
            'epcr':epcr_placeholder(rank-1,max_rows),'finding_type':'primary','notes':notes}


def validate_output(rows,known_gene):
    if not 1<=len(rows)<=10:raise ValueError('Submission row count outside 1-10')
    seen=set();prev=2.0
    for r in rows:
        if r['proband_id']!='PROBAND01':raise ValueError('Unexpected proband')
        x=float(r['epcr'])
        if not 0<x<=1 or x>=prev:raise ValueError('EPCR not strictly descending in (0,1]')
        prev=x
        v1=(r['chrom_1'],int(r['pos_1']),r['ref_1'],r['alt_1'])
        if v1 not in known_gene:raise ValueError('Unknown locus')
        p2=[r[z] for z in ['chrom_2','pos_2','ref_2','alt_2']]
        if any(p2) and not all(p2):raise ValueError('Incomplete second allele')
        if all(p2):
            v2=(r['chrom_2'],int(r['pos_2']),r['ref_2'],r['alt_2'])
            if v2 not in known_gene or v1==v2:raise ValueError('Unknown/reused pair allele')
            k=tuple(sorted((v1,v2)))
        else:k=(v1,)
        if k in seen:raise ValueError('Duplicate hypothesis')
        seen.add(k)
        if r['finding_type']!='primary':raise ValueError('Unsupported finding type')
    return True


def reasons(h):
    a=h['a'];b=h['b']
    parts=['HPO_ASSOCIATION_NOT_DISEASE_PROBABILITY']
    if h['prior_present']:parts.append('PREVIOUS_P1_P2_REVIEW')
    if any(r.get('clinvar_status')=='PATHOGENIC_LABEL_REVIEW_CONDITION' for r in (a,b) if r):
        parts.append('CLINVAR_PATHOGENIC_LABEL_CONDITION_NOT_VERIFIED')
    if any(r and len(r['ref'])!=len(r['alt']) for r in (a,b)):
        parts.append('ALLELE_LENGTH_DIFF_NOT_FRAMESHIFT_PREDICTION')
    parts+=['EXON_CDS_INTERVAL_ONLY','PHASE_NOT_VALIDATED','GNOMAD_AF_NOT_EVALUATED',
            'FUNCTIONAL_EFFECT_NOT_CALLED','MISSING_INHERITANCE_AND_FAMILY_DATA']
    return ';'.join(parts)


def row_audit(h,rank,is_selected):
    a=h['a'];b=h['b'];v1=h['v1'];v2=h['v2'] or ('','','','')
    return {'candidate_rank':str(rank),'selection':'YES' if is_selected else 'NO',
        'hypothesis_type':h['type'],'gene':h['gene'],
        'chrom_1':v1[0],'pos_1':v1[1],'ref_1':v1[2],'alt_1':v1[3],
        'chrom_2':v2[0],'pos_2':v2[1],'ref_2':v2[2],'alt_2':v2[3],
        'evidence_score_NOT_PROBABILITY':f"{h['score']:.3f}",
        'HPO_disease_overlap':str(max(fnum(a.get('public_disease_max_HPO_n')),fnum(b.get('public_disease_max_HPO_n')) if b else 0)),
        'HPO_gene_overlap':str(max(fnum(a.get('public_HPO_exact_n')),fnum(b.get('public_HPO_exact_n')) if b else 0)),
        'HPO_gene_IDF_proxy':str(max(fnum(a.get('public_HPO_IDF_score')),fnum(b.get('public_HPO_IDF_score')) if b else 0)),
        'variant1_GT':a.get('GT',''),'variant2_GT':b.get('GT','') if b else '',
        'variant1_DP':a.get('DP',''),'variant2_DP':b.get('DP','') if b else '',
        'variant1_GQ':a.get('GQ',''),'variant2_GQ':b.get('GQ','') if b else '',
        'variant1_position_proxy':a.get('genecode_location_proxy',''),
        'variant2_position_proxy':b.get('genecode_location_proxy','') if b else '',
        'variant1_ClinVar':a.get('clinvar_status',''),
        'variant2_ClinVar':b.get('clinvar_status','') if b else '',
        'variant1_prior_tier':a.get('previous_review_tier',''),
        'variant2_prior_tier':b.get('previous_review_tier','') if b else '',
        'variant1_ClinVar_condition':a.get('clinvar_CLNDN',''),
        'variant2_ClinVar_condition':b.get('clinvar_CLNDN','') if b else '',
        'phasing':PHASE if b else 'NOT_APPLICABLE','population_AF':AF,
        'coding_consequence':CODE,'reason_flags':reasons(h)}


def draw_html(audited):
    css=('body{font-family:system-ui,"Meiryo",sans-serif;margin:2em;line-height:1.55}'
         'table{border-collapse:collapse;display:block;overflow-x:auto}td,th{border:1px solid #ccc;'
         'padding:5px;vertical-align:top;white-space:nowrap}th{background:#edf2f7;position:sticky;top:0}'
         '.warn{font-weight:bold;color:#8a3800}')
    wanted=['candidate_rank','selection','hypothesis_type','gene','chrom_1','pos_1',
            'chrom_2','pos_2','evidence_score_NOT_PROBABILITY',
            'HPO_disease_overlap','HPO_gene_overlap','variant1_GT','variant2_GT',
            'variant1_position_proxy','variant2_position_proxy','variant1_ClinVar',
            'variant2_ClinVar','phasing']
    text=['<!doctype html><html lang="ja"><head><meta charset="utf-8">',
          "<meta http-equiv=\"Content-Security-Policy\" content=\"default-src 'none'; style-src 'unsafe-inline'\">",
          '<style>'+css+'</style></head><body><h1>E08｜非公開・自動評価順位</h1>',
          '<p class="warn">このHTML・CSVは患者固有の変異情報を含みます。ChatGPT、Colab、GitHubへ送信しないでください。</p>',
          '<p>この順位は限定的な証拠によるヒューリスティックです。EPCRの確率較正は未実施です。'
          '各行の二変異の位相・変異機能・集団頻度・疾患整合性は未検証です。</p>',
          '<p>選択はプログラムによる自動順位です。上位20件のレビュー表示です。詳細監査は同フォルダの監査CSVを参照してください。</p>',
          '<table><thead><tr>']
    text+=['<th>'+html.escape(k)+'</th>' for k in wanted]
    text+=['</tr></thead><tbody>']
    for r in audited[:20]:
        text+=['<tr>']+['<td>'+html.escape(str(r.get(k,'')),quote=True)+'</td>' for k in wanted]+['</tr>']
    text+=['</tbody></table></body></html>']
    return ''.join(text)


def methods(selected,allpairs,stats,max_per_gene):
    n_pair=sum(h['v2'] is not None for h in selected)
    n_single=len(selected)-n_pair
    return f'''# DGL-COMP-01 | Track 1 Automated Evidence Ranking — Methods DRAFT

## Research objective
Exploratory, nonclinical research submission for the Rare Disease, Real Kid MVA Hackathon 2026. The automated answer key is based on a clinically validated two-variant compound-heterozygous diagnosis, but this analysis does NOT assert any proposed pair is in trans, pathogenic, or the true causal pair. **Submission not performed by the program.**

## Data provenance, access and privacy
Restricted challenge GRCh38 VCF and the clinical report were processed only on the participant's Windows 11 device under a protected non-Git `patient_private` folder. The clinical source lists eight curated HPO IDs. Public HPO gene–phenotype associations, GENCODE v48 Basic GTF (GRCh38), and ClinVar GRCh38 VCF were downloaded as public references. Google Colab (Google AI Pro) was used for public reference downloads only. The local E08 script has no network functions and does not transmit any variant, phenotype, or patient-level information. Selected-variant CSV is restricted and must be uploaded ONLY to the authorized official competition submission portal by a human after reviewing the challenge terms and privacy requirements. Never commit that CSV, patient genotype files, audit reports, or clinical records to GitHub.

## Input pipeline and previously completed checks
- Original patient VCF: 5,012,204 records scanned in E03, one sample, GRCh38. Variant QC fields GT, DP, GQ and FILTER read locally. PASS is technical filtering, not clinical pathogenicity.
- Initial 7 MVA-related genes were expanded to 50 additional public HPO-associated genes (57 targeted total). E04-H generated 5,467 candidate allele rows from 56 genes with candidate variants. These 56 genes are NOT an exhaustive genome-wide causal search.
- E04-I annotated candidate positions with GENCODE v48 Basic GTF: 67 CDS-position overlaps, 272 non-CDS exon overlaps, 5,128 other gene-region candidates. Positional overlaps do NOT determine HGVS consequences such as missense, stop-gain, frameshift or splice damage.
- Local ClinVar full-GRCh38 VCF lookup scanned 4,599,326 public rows. Labels in the candidate table: 453 benign-type, 2 conflicting, 3 VUS, 1 pathogenic-type label requiring review of the variant-condition relationship, 5,008 with no exact match. These counts can include labels retained from earlier 7-gene review; no-match does not imply benignity. REF/ALT left-normalization and transcript-specific ClinVar matching were not performed.
- E06/E07 preserved three prior P1/P2 review candidates. All E08 scores were recalculated from the existing E04-I rows, not copied from E07's subjective EPCR or E04-I's review_score.

## Automated candidate construction and ranking (E08 v0.1)
E08 reads `E04I_ALL_CANDIDATES_PRIVATE_DO_NOT_UPLOAD.csv` and deduplicates alternate alleles by gene plus exact CHROM/POS/REF/ALT. To qualify for the *pair* shortlist, both alleles must have GT 0/1 or 1/0 (or a phased but non-trans-proving equivalent), FILTER PASS, DP ≥10, GQ ≥20, and neither must have a ClinVar benign-type label. Each must also have a GENCODE CDS/exon/near-exon-position proxy, a ClinVar pathogenic-type label pending condition review, or an original P1/P2 review marker. This screening can miss intronic and noncoding disease alleles; the exclusion is not a clinical benignity judgment. Genotypes 1/2 are not interpreted as individual heterozygous alleles. Pair alleles must be at distinct positions within one annotated gene. Four top evidence-positionally supported distinct loci per gene are considered, with previously prioritized loci retained for review even if outside these four. All within-gene two-locus combinations are scored. If insufficient paired hypotheses exist, eligible single-variant hypotheses are allowed as fallback.

The REVIEW evidence scoring rule is prespecified and deterministic, not trained and not a probability model:

- **Gene-level HPO**: `2 × min(max same-gene disease-term overlap, 8)` + `0.75 × min(gene-term overlap, 8)` + `0.05 × min(HPO-IDF proxy, 40)`.
- **Per-variant position support**: CDS interval overlap +4, non-CDS exon overlap or ±2bp exon-boundary proxy +1.5, other interval +0.
- **Per-variant supplementary support**: +2 for a ClinVar pathogenic-type **label pending disease/review-level confirmation**; +1 for previously flagged P1/P2; +0.5 for allele-length difference at a positional coding/exon locus (this is NOT a frameshift prediction).
- **Pair balance**: +0.25 × minimum of the two positional-support scores. Quality fields are gates, not pathogenicity bonuses. Scores are not ACMG classifications.
- Sort descending by total evidence score with deterministic tie-breaks. Select at most {max_per_gene} candidate pairs per gene to diversify the 10-row shortlist. Output top {len(selected)} candidates, including {n_pair} two-variant hypotheses and {n_single} single-variant fallback hypotheses. Candidate pairs scored before gene-diversified selection: {len(allpairs)}. No manual EPCR or candidate-order edits are required.

## EPCR interpretation
The official Track 1 field `epcr` is syntactically required to be in (0,1]. **True Estimated Probability of Causal Relationship cannot be reliably estimated from the available unphased, non-calibrated evidence**. To preserve evidence-based *ordering*, E08 places equally spaced **rank-only placeholder numbers** in the EPCR column (formula `(10 - zero_based_rank) / 11` for up to ten rows, not a statistically estimated posterior). Thus these numerical values are **NOT calibrated probabilities**, should not be reported as clinical risk, and may not support scientific claims about the numerical differences between rows. The official scorer sorts by these numbers; rankings are attributable to the separate transparent evidence score. This is an exploratory baseline and the EPCR values should be interpreted exclusively as a submission-interface encoding of ranks.

## Key limitations and falsification priorities
1. Neither trans phasing nor segregation with parental genotypes was validated. Two heterozygous loci in one gene may be on the same chromosome (cis). No clinical diagnosis is established.
2. Transcript-specific variant effects, in silico predictions, actual frameshift status, splice effects, structure variants, CNV, and regulatory variants were not evaluated.
3. gnomAD allele frequencies were NOT measured; the original INFO AF is not a substitute for population allele frequency.
4. Pathogenic-type ClinVar label may concern a different disease; VUS or no exact match does not establish benignity. All variant–disease matches require clinical and expert review.
5. A 57-gene interval search can miss a causal gene, especially variants not included by HPO or outside the initial prioritization. HPO terms themselves may be correlated or nonspecific.
6. There is no independent labeled training dataset, no posterior calibration, no known sensitivity or specificity and no validated clinical interpretation.

## Reproducibility and environment
Python 3.11 standard library, Windows 11, offline execution, all scoring thresholds/weights fixed in `e08_automated_rank.py` (no external APIs). Public references: GENCODE v48 Basic, ClinVar GRCh38, HPO ontology and phenotype association tables. The specific local output may be reproduced using the restricted input under the original challenge terms. An independent verifier should reproduce CSV syntax and review privacy exclusions before submission.

## Generative AI disclosure — verify before submission
OpenAI ChatGPT (GPT-6; ChatGPT Plus) assisted with the analysis scripts, scoring design, methods documentation, and interpretation of **aggregate** outputs; patient-specific VCF/variant tables/individual HPO identifiers and clinical document contents were not supplied in this conversation. Google AI Pro/Colab was used only for public references. User should verify their ChatGPT data-control/training settings and disclose them accurately. **Settings not verified in this report: [CONFIRM BEFORE SUBMISSION].** No generative-AI posterior cause probability is claimed.

## Official track specification and source references
- https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026
- https://huggingface.co/spaces/SageBio/rare-disease-real-kid-mva-hackathon-2026/blob/main/evaluation.py
- https://www.gencodegenes.org/human/release_48.html
- https://www.ncbi.nlm.nih.gov/clinvar/
- https://hpo.jax.org/

*This is a methods **DRAFT**; it requires review before official submission. The patient-level submission CSV and audit HTML must never be shared outside the authorized competition process.*
'''


def atomic_text(path,text):
    path.parent.mkdir(parents=True,exist_ok=True)
    with tempfile.NamedTemporaryFile('w',dir=path.parent,delete=False,encoding='utf-8',newline='',prefix='.tmp_E08_',suffix='.txt') as f:
        tmp=Path(f.name);f.write(text)
    try:os.replace(tmp,path)
    finally:tmp.unlink(missing_ok=True)


def atomic_csv(path,rows,fields):
    path.parent.mkdir(parents=True,exist_ok=True)
    with tempfile.NamedTemporaryFile('w',dir=path.parent,delete=False,encoding='utf-8',newline='',prefix='.tmp_E08_',suffix='.csv') as f:
        tmp=Path(f.name);wr=csv.DictWriter(f,fieldnames=fields,extrasaction='raise')
        wr.writeheader();wr.writerows(rows)
    try:os.replace(tmp,path)
    finally:tmp.unlink(missing_ok=True)


def run(root,overwrite=False,max_per_gene=2,verbose=True):
    root=Path(root).resolve()
    check_private(root/'patient_private')
    if not (1<=max_per_gene<=5):raise ValueError('max-per-gene must be 1..5')
    out=check_private(root/'patient_private'/'E08')
    if out.exists() and not overwrite and any((out/x).exists() for x in (OUTPUT,AUDIT,REPORT)):
        raise ValueError('E08 output already exists; use --overwrite only after backup/review')
    rows=read_input(root)
    pairs,singles,base,prior_seen=candidate_pool(rows)
    selected=rank_select(pairs,singles,max_per_gene=max_per_gene)
    known=set()
    for r in rows:
        try:known.add(loc(r))
        except (KeyError,ValueError,TypeError):pass
    selected_rows=[official_row(x,i+1,10) for i,x in enumerate(selected)]
    validate_output(selected_rows,known)
    selection_keys={hyp_key(x) for x in selected}
    audited=[]
    seen=set()
    for hyp in pairs+singles:
        k=hyp_key(hyp)
        if k in seen:continue
        seen.add(k)
        audited.append(row_audit(hyp,len(audited)+1,k in selection_keys))
    selected_prior=set()
    for hyp in selected:
        for r in (hyp['a'],hyp['b']):
            if r and is_prior(r):selected_prior.add((r['gene'],loc(r)))
    # No incomplete artifacts written if basic validation fails.
    output_summary={'experiment':'E08_AUTOMATED_EVIDENCE_RANKING_v0.1',
        **base,'selected_rows':len(selected),
        'selected_pair_rows':sum(x['v2'] is not None for x in selected),
        'selected_single_rows':sum(x['v2'] is None for x in selected),
        'selected_distinct_genes':len({x['gene'] for x in selected}),
        'prior_P1_P2_alleles_in_selected':len(selected_prior),
        'prior_P1_P2_alleles_not_selected':len(prior_seen-selected_prior),
        'rank_generated_automatically':True,
        'all_selected_alleles_revalidated_locally':True,
        'official_header_and_epcr_syntax_checked':True,
        'evidence_score_is_probability':False,
        'epcr_probability_calibrated':False,
        'epcr_numeric_values_are_rank_placeholders':True,
        'gnomad_AF_evaluated':False,
        'phase_validated':False,
        'transcript_consequence_computed':False,
        'patient_details_printed':False,'network_requests':0,
        'official_submission_performed':False,
        'output_private_only':True}
    out.mkdir(parents=True,exist_ok=True)
    # Private source hash deliberately only written within private methods text.
    atomic_csv(out/OUTPUT,selected_rows,OFFICIAL)
    atomic_csv(out/AUDIT,audited,EVIDENCE_FIELDS)
    private_preview=[row_audit(x,i+1,True) for i,x in enumerate(selected)]
    private_preview += [r for r in audited if r['selection']=='NO']
    atomic_text(out/PREVIEW,draw_html(private_preview))
    atomic_text(out/REPORT,methods(selected,pairs,output_summary,max_per_gene))
    atomic_text(out/SUMMARY,json.dumps(output_summary,ensure_ascii=False,indent=2))
    if verbose:
        print(json.dumps(output_summary,ensure_ascii=False,indent=2))
        print('LOCAL PRIVATE DRAFT ONLY. All clinical details stay in patient_private/E08. NO SUBMISSION.')
    return output_summary



def diagnose(root):
    root=Path(root).resolve()
    check_private(root/'patient_private')
    out=check_private(root/'patient_private'/'E08')
    inp=root/'patient_private'/'E04I'/SOURCE
    stats={'diagnostic':'E08_PREFLIGHT_v0.1','input_exists':inp.is_file(),
           'existing_E08_draft_exists':(out/OUTPUT).is_file(),
           'network_requests':0,'patient_variant_details_printed':False}
    if inp.is_file():
        data=read_input(root)
        pairs,singles,counters,prior=candidate_pool(data)
        stats.update({'candidate_rows':len(data),'pair_hypotheses_scored':len(pairs),
                      'eligible_genes_with_pairs':counters['genes_with_pair_hypotheses'],
                      'prior_P1_P2_distinct_alleles_seen':len(prior),
                      'ready_to_rank':bool(pairs or singles)})
    print(json.dumps(stats,ensure_ascii=False,indent=2))
    return stats

def test_fixture(root):
    dest=root/'patient_private'/'E04I';dest.mkdir(parents=True)
    names=sorted(NEED|{'public_HPO_IDF_score','previous_review_tier','clinvar_CLNDN'})
    cases=[]
    for i in range(1,13):
        # Distinct locus-pairs per gene, with some 3rd lower-quality candidates.
        for j in range(2):
            cases.append({'gene':f'SYNTH_GENE_{i:02d}','chrom':str(i),'pos':str(10000+i*100+j*50),
                'ref':'A','alt':'G' if j==0 else 'T','GT':'0/1',
                'DP':'30','GQ':'70','FILTER':'PASS','public_HPO_exact_n':str(5 if i<4 else 2),
                'public_disease_max_HPO_n':str(4 if i<4 else 2),
                'public_HPO_IDF_score':'9',
                'genecode_location_proxy':'CDS_POSITIONAL_OVERLAP' if i<7 else 'NONCDS_EXON_OVERLAP',
                'clinvar_status':'PATHOGENIC_LABEL_REVIEW_CONDITION' if i==3 and j==0 else 'NO_CLINVAR_EXACT_MATCH_NOT_BENIGN',
                'previous_review_tier':'P1_SYNTH' if i==3 and j==0 else '',
                'clinvar_CLNDN':'SYNTHETIC_OTHER_DISEASE'})
    cases += [dict(cases[0])|{'pos':'3000','GT':'1/2'},
              dict(cases[0])|{'pos':'4000','clinvar_status':'BENIGN_LABEL_CHECK_CONDITION'}]
    with (dest/SOURCE).open('w',newline='',encoding='utf-8') as f:
        wr=csv.DictWriter(f,fieldnames=names);wr.writeheader();wr.writerows(cases)
    return cases


def self_test():
    with tempfile.TemporaryDirectory() as td:
        root=Path(td);cases=test_fixture(root)
        a=run(root,False,2,False)
        assert a['input_rows']==len(cases)
        assert a['selected_rows']==10 and a['selected_pair_rows']==10
        assert a['official_header_and_epcr_syntax_checked']
        assert a['prior_P1_P2_distinct_alleles_seen']==1
        out=root/'patient_private'/'E08'
        with (out/OUTPUT).open(encoding='utf-8',newline='') as f:
            rd=csv.DictReader(f);testrows=list(rd)
        assert rd.fieldnames==OFFICIAL and len(testrows)==10
        assert all(x['chrom_2'] and x['pos_2'] for x in testrows)
        e=[float(x['epcr']) for x in testrows]
        assert all(e[i]>e[i+1] for i in range(len(e)-1))
        old=(out/OUTPUT).read_bytes()
        try:run(root,False,2,False)
        except ValueError:pass
        else:raise AssertionError('Existing file overwritten without authorization')
        b=run(root,True,2,False)
        assert b==a and (out/OUTPUT).read_bytes()==old
        assert not (root/'patient_private'/'E07').exists()
        # Explicitly prohibit putting patient outputs inside Git repository.
        (root/'.git').mkdir()
        try:run(root,True,2,False)
        except ValueError:pass
        else:raise AssertionError('Git privacy guard not triggered')
    print('SELF_TEST_PASS: True (synthetic patient alleles only)')


def main():
    p=argparse.ArgumentParser(description='Offline and private evidence-based E08 ranking for Track1')
    p.add_argument('--root',help='Private MVA project root folder')
    p.add_argument('--self-test',action='store_true')
    p.add_argument('--diagnose',action='store_true',help='Sanitized local-input preflight, writes no patient outputs')
    p.add_argument('--max-per-gene',type=int,default=2,help='Diversification cap 1..5 (default 2)')
    p.add_argument('--overwrite',action='store_true',help='Replace prior E08 output after explicit review')
    args=p.parse_args()
    try:
        if args.self_test:self_test();return
        if not args.root:raise ValueError('--root required')
        if args.diagnose:diagnose(args.root);return
        run(args.root,args.overwrite,args.max_per_gene)
    except Exception as e:
        # Error messages deliberately omit all patient-specific variants and gene names.
        print('ERROR: '+type(e).__name__+' — details suppressed for private data safety.',file=sys.stderr)
        raise SystemExit(1)

if __name__=='__main__':main()

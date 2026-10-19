#!/usr/bin/env python3
"""DGL-COMP-01 Track 1: format-only preflight against official 2026 template.
No scoring, biological validation, or upload is performed. Stdlib only.
"""
import argparse, csv, io, math, re, sys
from pathlib import Path

COLUMNS = [
    "proband_id", "chrom_1", "pos_1", "ref_1", "alt_1", "chrom_2", "pos_2",
    "ref_2", "alt_2", "epcr", "finding_type", "notes"
]
CHROM_REGEX = re.compile(r"^chr(?:[1-9]|1[0-9]|2[0-2]|X|Y|M)$")


def check_csv(f):
    errors, warnings = [], []
    reader = csv.reader(f)
    try:
        header = next(reader)
    except StopIteration:
        return 0, ["CSV is empty"], warnings
    if header != COLUMNS:
        errors.append("Header does not match official 12-column template/order")
        return 0, errors, warnings
    rows = list(reader)
    if not rows:
        errors.append("No candidate rows; the header-only file cannot be submitted")
    if len(rows) > 10:
        errors.append("More than 10 candidate rows")
    last_epcr = None
    unique_rows = set()
    proband_ids = set()
    for num, vals in enumerate(rows, start=2):
        if len(vals) != 12:
            errors.append(f"Line {num}: expected 12 fields; got {len(vals)}")
            continue
        row = dict(zip(COLUMNS, (v.strip() for v in vals)))
        if not row['proband_id']:
            errors.append(f"Line {num}: proband_id is empty")
        proband_ids.add(row['proband_id'])
        # Explicitly reject the two mock candidates shipped in the official template.
        if (row['proband_id'] == 'PROBAND01' and
            (row['chrom_1'], row['pos_1'], row['ref_1'], row['alt_1']) in
            {('chr1','100000','A','T'), ('chr7','300000','C','T')}):
            errors.append(f"Line {num}: resembles an official illustrative/example row; replace with actual evidence-backed candidate")
        for idx in (1,):
            for field in (f'chrom_{idx}',f'pos_{idx}',f'ref_{idx}',f'alt_{idx}'):
                if not row[field]: errors.append(f"Line {num}: {field} is empty")
        pair = [row[c] for c in ('chrom_2','pos_2','ref_2','alt_2')]
        if any(pair) and not all(pair):
            errors.append(f"Line {num}: second variant must have all 4 fields or none")
        for idx in (1,2):
            if idx == 2 and not any(pair): continue
            chrom, pos = row[f'chrom_{idx}'], row[f'pos_{idx}']
            if chrom and not CHROM_REGEX.fullmatch(chrom):
                warnings.append(f"Line {num}: unusual GRCh38 chromosome label: {chrom}")
            try:
                if int(pos) < 1: raise ValueError()
            except ValueError:
                errors.append(f"Line {num}: pos_{idx} must be a positive integer")
            if row[f'ref_{idx}'].upper() == row[f'alt_{idx}'].upper() and row[f'ref_{idx}']:
                warnings.append(f"Line {num}: REF and ALT are identical for variant {idx}")
        try:
            epcr = float(row['epcr'])
            if not math.isfinite(epcr) or not (0 <= epcr <= 1): raise ValueError()
            if last_epcr is not None and epcr > last_epcr:
                errors.append(f"Line {num}: EPCR is not sorted in descending order")
            last_epcr = epcr
        except ValueError:
            errors.append(f"Line {num}: EPCR must be a finite number between 0 and 1")
        if row['finding_type'] not in ('primary','secondary'):
            errors.append(f"Line {num}: finding_type must be primary or secondary")
        key = tuple(row[c] for c in COLUMNS if c not in ('epcr','notes'))
        if key in unique_rows: warnings.append(f"Line {num}: duplicate candidate")
        unique_rows.add(key)
    if len(proband_ids) > 1:
        warnings.append("More than one proband_id (this competition concerns a single case)")
    return len(rows),errors,warnings


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('csv_path',nargs='?',help='Your candidate predictions CSV')
    p.add_argument('--self-test',action='store_true',help='No patient data; tests synthetic examples')
    args=p.parse_args()
    if args.self_test:
        sample=','.join(COLUMNS)+'\n'+('SYNTHETIC_CASE,chr1,100000,A,T,,,,,0.9,primary,synthetic-test-only\n')
        n,e,w=check_csv(io.StringIO(sample))
        assert n==1 and not e, (n,e,w)
        assert check_csv(io.StringIO(','.join(COLUMNS)+'\n'))[1]
        print('SELF-TEST PASS (synthetic data only)')
        return 0
    if not args.csv_path:
        p.error('CSV path is required unless --self-test is given')
    path=Path(args.csv_path)
    with path.open('r',newline='',encoding='utf-8-sig') as f:
        n,e,w=check_csv(f)
    print('DGL-COMP-01 Track 1 preflight')
    print('Candidate count:',n)
    print('Header/shape/value checks:', 'PASS' if not e else 'FAIL')
    for msg in e: print('ERROR:',msg)
    for msg in w: print('WARNING:',msg)
    print('Official evaluator and biological validity: NOT TESTED')
    print('No network activity. Candidate data were not printed or uploaded.')
    return 1 if e else 0


if __name__=='__main__':
    sys.exit(main())

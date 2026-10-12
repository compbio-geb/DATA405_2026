# Build script for the DATA405 Lab 1 data package (see data.qmd for provenance).
# Inputs: BIOL419 Lab 4 folder (chr22.fa, chr22.gtf, ENCFF729YAX_subset100k.fastq, STAR BAM, RSEM genes.results)
# and airway_scaledcounts.csv / airway_metadata.csv (bioboot/bimm143_S18 on GitHub).
# Requires: biopython, pysam, pandas. Run from the folder containing those inputs; writes to ../site/data/.

import re, random, collections, csv, sys
from Bio import SeqIO
from Bio.Seq import Seq
import pysam, pandas as pd
out = sys.argv[1] if len(sys.argv) > 1 else "../site/data/"
chr22 = str(next(SeqIO.parse('chr22.fa', 'fasta')).seq).upper()
def attr(a, k):
    m = re.search(k + ' "([^"]+)"', a); return m.group(1) if m else ''
genes = []
for l in open('chr22.gtf'):
    if l.startswith('#'): continue
    r = l.rstrip('\n').split('\t')
    if r[2] != 'gene': continue
    genes.append(dict(gene_id=attr(r[8], 'gene_id').split('.')[0], gene_name=attr(r[8], 'gene_name'), gene_type=attr(r[8], 'gene_type'), start=int(r[3]), end=int(r[4]), strand=r[6]))
with open(out + 'chr22_genes.tsv', 'w') as f:
    w = csv.writer(f, delimiter='\t'); w.writerow(['gene_id', 'gene_name', 'gene_type', 'start', 'end', 'strand', 'length_bp'])
    for g in genes: w.writerow([g['gene_id'], g['gene_name'], g['gene_type'], g['start'], g['end'], g['strand'], g['end'] - g['start'] + 1])
gid2 = {g['gene_id']: g for g in genes}
rows = list(csv.DictReader(open('sample_rsem.genes.results'), delimiter='\t'))
with open(out + 'chr22_expression.tsv', 'w') as f:
    w = csv.writer(f, delimiter='\t'); w.writerow(['gene_id', 'gene_name', 'gene_type', 'length', 'expected_count', 'TPM'])
    for r in rows:
        gid = r['gene_id'].split('.')[0]; g = gid2.get(gid, {})
        w.writerow([gid, g.get('gene_name', ''), g.get('gene_type', ''), r['length'], r['expected_count'], r['TPM']])
sel = [l.rstrip('\n').split('\t') for l in open('chr22.gtf') if 'gene_name "COMT";' in l]
tx = [r for r in sel if r[2] == 'transcript' and 'Ensembl_canonical' in r[8]][0]
tid = attr(tx[8], 'transcript_id'); gs, ge = int(tx[3]), int(tx[4])
cds = sorted([r for r in sel if r[2] in ('CDS', 'stop_codon') and attr(r[8], 'transcript_id') == tid], key=lambda r: int(r[3]))
exons = sorted([r for r in sel if r[2] == 'exon' and attr(r[8], 'transcript_id') == tid], key=lambda r: int(r[3]))
cdsseq = ''.join(chr22[int(r[3]) - 1:int(r[4])] for r in cds); prot = str(Seq(cdsseq).translate())
assert prot[157] == 'V' and cdsseq[471:474] == 'GTG' and prot.count('*') == 1 and prot.endswith('*')
gene_seq = chr22[gs - 1:ge]
wrap = lambda s: '\n'.join(s[i:i + 70] for i in range(0, len(s), 70))
open(out + 'COMT_gene.fa', 'w').write(f'>COMT chr22:{gs}-{ge} + strand GRCh38 (GENCODE canonical transcript {tid})\n' + wrap(gene_seq) + '\n')
open(out + 'COMT_cds.fa', 'w').write(f'>COMT_CDS {tid} MB-COMT coding sequence incl. stop, {len(cdsseq)} nt, GRCh38\n' + wrap(cdsseq) + '\n')
open(out + 'COMT_protein.fa', 'w').write(f'>COMT_protein MB-COMT, {len(prot) - 1} aa\n' + wrap(prot[:-1]) + '\n')
with open(out + 'COMT_exons.tsv', 'w') as f:
    w = csv.writer(f, delimiter='\t'); w.writerow(['transcript', 'feature', 'chrom_start', 'chrom_end', 'gene_offset_start_0based', 'gene_offset_end_exclusive'])
    for r in exons: w.writerow([tid, 'exon', r[3], r[4], int(r[3]) - gs, int(r[4]) - gs + 1])
    for r in cds: w.writerow([tid, r[2], r[3], r[4], int(r[3]) - gs, int(r[4]) - gs + 1])
random.seed(405)
recs = list(SeqIO.parse('ENCFF729YAX_subset100k.fastq', 'fastq'))
SeqIO.write(random.sample(recs, 5000), out + 'reads_5k.fastq', 'fastq')
b = pysam.AlignmentFile('Aligned.sortedByCoord.out.bam')
cov = collections.Counter(); alt = collections.defaultdict(collections.Counter)
for r in b.fetch(until_eof=True):
    if r.is_unmapped or r.mapping_quality < 20: continue
    for qp, rp in r.get_aligned_pairs(matches_only=True):
        if r.query_qualities[qp] < 20: continue
        cov[rp] += 1; base = r.query_sequence[qp]
        if base != chr22[rp]: alt[rp][base] += 1
pc = [(g['start'], g['end'], g['gene_name']) for g in genes if g['gene_type'] == 'protein_coding']
def gene_at(p):
    h = [n for s, e, n in pc if s <= p <= e]; return h[0] if h else ''
with open(out + 'chr22_mismatches.tsv', 'w') as f:
    w = csv.writer(f, delimiter='\t'); w.writerow(['chrom', 'pos', 'ref', 'alt', 'alt_count', 'depth', 'alt_fraction', 'gene_name'])
    for p in sorted(alt):
        if cov[p] < 5: continue
        for a, c in alt[p].items():
            if c >= 2: w.writerow(['chr22', p + 1, chr22[p], a, c, cov[p], round(c / cov[p], 3), gene_at(p + 1)])
a = pd.read_csv('airway_scaledcounts.csv'); a = a[a.iloc[:, 1:].sum(axis=1) >= 10]
a.rename(columns={'ensgene': 'gene_id'}).to_csv(out + 'airway_counts.tsv', sep='\t', index=False, float_format='%.0f')
pd.read_csv('airway_metadata.csv').to_csv(out + 'airway_metadata.tsv', sep='\t', index=False)
print('done')

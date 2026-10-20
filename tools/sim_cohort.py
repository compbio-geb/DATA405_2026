"""Simulate the Lab 3 GWAS/PRS cohort (seeded). Writes cohort.vcf.gz, pheno.txt, discovery.txt, target.txt, truth.tsv."""
import numpy as np, gzip, sys
rng = np.random.default_rng(405)
out = sys.argv[1] if len(sys.argv) > 1 else "data"
N1, N2 = 1750, 750            # two ancestry groups
N = N1 + N2
NCHR, PER = 22, 500
BLOCK = 10
FST = 0.03
group = np.r_[np.zeros(N1, int), np.ones(N2, int)]
perm = rng.permutation(N); group = group[perm]   # shuffle order
ids = [f"P{i+1:04d}" for i in range(N)]

chrlen = np.array([248,242,198,190,181,171,159,145,138,134,135,133,114,107,102,90,83,80,59,64,47,51]) * 1_000_000  # GRCh38, Mb
G = np.zeros((NCHR * PER, N), np.int8)
chrom = np.repeat(np.arange(1, NCHR + 1), PER)
pos = np.zeros(NCHR * PER, int)
p_anc = np.zeros(NCHR * PER); p_grp = np.zeros((NCHR * PER, 2))
for c in range(NCHR):
    starts = np.sort(rng.choice(np.arange(1_000_000, chrlen[c] - 1_000_000, 20_000), PER // BLOCK, replace=False))
    pp = []
    for s in starts:
        pp += list(s + np.sort(rng.choice(np.arange(0, 15_000, 50), BLOCK, replace=False)))
    pos[c * PER:(c + 1) * PER] = pp
# allele frequencies: ancestral U-shaped-ish, group via Balding-Nichols
M = NCHR * PER
pa = rng.beta(0.6, 0.6, M) * 0.9 + 0.05
rare = rng.random(M) < 0.04
pa[rare] = rng.uniform(0.0005, 0.005, rare.sum())
a = pa * (1 - FST) / FST; b = (1 - pa) * (1 - FST) / FST
for k in range(2):
    p_grp[:, k] = np.where(rare, pa, rng.beta(a, b))
p_grp = np.clip(p_grp, 0.0, 1.0)
RHO = 0.85
for blk in range(M // BLOCK):
    sl = slice(blk * BLOCK, (blk + 1) * BLOCK)
    for k in range(2):
        idx = np.where(group == k)[0]
        for hap in range(2):
            h = (rng.random(len(idx)) < p_grp[blk * BLOCK, k]).astype(np.int8)
            G[blk * BLOCK, idx] += h
            for j in range(1, BLOCK):
                m = blk * BLOCK + j
                copy = rng.random(len(idx)) < RHO
                # keep frequency right: copy, else draw fresh
                fresh = (rng.random(len(idx)) < p_grp[m, k]).astype(np.int8)
                h = np.where(copy, h, fresh)
                G[m, idx] += h
# because copying distorts freq, that's fine (freq ~ mix of block-leader and own)

# HWE failures (genotyping-error SNPs): 25 common SNPs lose most heterozygotes
common = np.where((G.mean(1) / 2 > 0.2) & (G.mean(1) / 2 < 0.8))[0]
hwe_bad = rng.choice(common, 25, replace=False)
for m in hwe_bad:
    het = np.where(G[m] == 1)[0]
    flip = het[rng.random(len(het)) < 0.6]
    G[m, flip] = rng.choice([0, 2], len(flip)).astype(np.int8)

# phenotype
freq = G.mean(1) / 2
ok = (freq > 0.05) & (freq < 0.95)
ok[hwe_bad] = False
big_chr = [2, 6, 11, 16, 19]
big = []
for c in big_chr:
    cand = np.where((chrom == c) & ok & (np.abs(p_grp[:, 0] - p_grp[:, 1]) < 0.1))[0]
    big.append(rng.choice(cand))
big = np.array(big)
small = rng.choice(np.setdiff1d(np.where(ok)[0], big), 400, replace=False)
std = lambda g: (g - g.mean()) / g.std()
beta = np.zeros(M)
beta[big] = np.sqrt(0.018) * rng.choice([-1, 1], len(big))
beta[small] = rng.normal(0, np.sqrt(0.30 / len(small)), len(small))
causal = np.r_[big, small]
gen = sum(beta[m] * std(G[m].astype(float)) for m in causal)
y = gen + 1.0 * group + rng.normal(0, np.sqrt(1 - gen.var()), N)
ldl = 130 + 30 * (y - y.mean()) / y.std()

# missingness
Gm = G.astype(np.int16)
bad_snp = rng.choice(np.setdiff1d(np.arange(M), causal), int(0.015 * M), replace=False)
for m in bad_snp:
    Gm[m, rng.random(N) < rng.uniform(0.05, 0.2)] = -1
bad_ppl = rng.choice(N, 15, replace=False)
for i in bad_ppl:
    Gm[rng.random(M) < 0.08, i] = -1
Gm[rng.random(Gm.shape) < 0.002] = -1   # background

bases = np.array(list("ACGT"))
ref = rng.integers(0, 4, M); alt = (ref + rng.integers(1, 4, M)) % 4
idsnp = [f"{chrom[m]}:{pos[m]}:{bases[ref[m]]}:{bases[alt[m]]}" for m in range(M)]
gt = np.array(["./.", "0/0", "0/1", "1/1"])
with gzip.open(f"{out}/cohort.vcf.gz", "wt", compresslevel=6) as f:
    f.write("##fileformat=VCFv4.2\n##source=DATA405 Lab 3 simulation (not real people)\n")
    for c in range(NCHR):
        f.write(f"##contig=<ID={c+1},length={chrlen[c]}>\n")
    f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
    f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(ids) + "\n")
    for m in range(M):
        f.write(f"{chrom[m]}\t{pos[m]}\t{idsnp[m]}\t{bases[ref[m]]}\t{bases[alt[m]]}\t.\tPASS\t.\tGT\t" + "\t".join(gt[Gm[m] + 1]) + "\n")
with open(f"{out}/pheno.txt", "w") as f:
    f.write("#IID\tLDL\n")
    for i in range(N): f.write(f"{ids[i]}\t{ldl[i]:.1f}\n")
order = rng.permutation(N)
disc, targ = sorted(order[:2000]), sorted(order[2000:])
open(f"{out}/discovery.txt", "w").write("#IID\n" + "".join(ids[i] + "\n" for i in disc))
open(f"{out}/target.txt", "w").write("#IID\n" + "".join(ids[i] + "\n" for i in targ))
# truth (private, for the key)
with open(f"{out}/../truth.tsv", "w") as f:
    f.write("ID\tCHR\tPOS\tbeta_std\tkind\tfreq_g0\tfreq_g1\n")
    for m in causal:
        f.write(f"{idsnp[m]}\t{chrom[m]}\t{pos[m]}\t{beta[m]:.4f}\t{'big' if m in big else 'small'}\t{p_grp[m,0]:.3f}\t{p_grp[m,1]:.3f}\n")
    for m in hwe_bad: f.write(f"{idsnp[m]}\t{chrom[m]}\t{pos[m]}\t0\thwe_bad\t\t\n")
with open(f"{out}/../truth_people.tsv", "w") as f:
    f.write("IID\tgroup\tbad_missing\n")
    for i in range(N): f.write(f"{ids[i]}\t{group[i]}\t{int(i in set(bad_ppl))}\n")
print("SNPs", M, "people", N, "big", [idsnp[m] for m in big], "gen var", round(gen.var(), 3))

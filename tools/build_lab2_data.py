"""Build the Lab 2 data files (DATA405, Lecture 12-13).

rs2814778_by_population.tsv : REAL 1000 Genomes phase 3 genotype counts (Ensembl REST, fetched Oct 6, 2026)
gwas_region.tsv / gwas_snps.tsv : SIMULATED genotypes for 15 SNPs in one region + a trait
stratification.tsv : SIMULATED two-group cohort with no genetic effect
"""
import numpy as np, pandas as pd
out = "data/lab2/"

# --- real counts: rs2814778 (ACKR1 / Duffy), ref T, alt C; genotype counts TT, CT, CC
pops = pd.DataFrame({
    "population": ["AFR", "AMR", "EAS", "EUR", "SAS"],
    "description": ["African", "Admixed American", "East Asian", "European", "South Asian"],
    "TT": [4, 299, 504, 497, 489],
    "CT": [40, 42, 0, 6, 0],
    "CC": [617, 6, 0, 0, 0],
})
assert pops[["TT", "CT", "CC"]].sum().tolist() == [1793, 88, 623]
pops.to_csv(out + "rs2814778_by_population.tsv", sep="\t", index=False)

# --- simulated GWAS region: 15 SNPs in three LD blocks; SNP 8 is causal
rng = np.random.default_rng(2026)
N = 2000
nsnp = 15
blocks = [range(0, 5), range(5, 10), range(10, 15)]
anchor_freq = [0.35, 0.30, 0.40]
copy_prob = [0.95, 0.97, 0.94, 0.90, 0.80,   # how tightly each SNP follows its block anchor
             0.85, 0.93, 1.00, 0.96, 0.75,
             0.92, 0.88, 0.97, 0.83, 0.90]
def haplotypes(n):
    H = np.zeros((n, nsnp), dtype=int)
    for b, idx in enumerate(blocks):
        anchor = rng.random(n) < anchor_freq[b]
        for j in idx:
            follow = rng.random(n) < copy_prob[j]
            other = rng.random(n) < anchor_freq[b]
            H[:, j] = np.where(follow, anchor, other)
    return H
G = haplotypes(N) + haplotypes(N)                    # two copies per person
y = 0.30 * G[:, 7] + rng.normal(0, 1, N)
pos = 1_000_000 + np.array([0, 12, 25, 31, 48, 70, 82, 95, 101, 118, 140, 152, 166, 171, 190]) * 1000
snps = pd.DataFrame({"snp": [f"snp{j+1}" for j in range(nsnp)], "position": pos})
snps.to_csv(out + "gwas_snps.tsv", sep="\t", index=False)
df = pd.DataFrame(G, columns=snps["snp"])
df.insert(0, "person", [f"P{i+1:04d}" for i in range(N)])
df["trait"] = np.round(y, 3)
df.to_csv(out + "gwas_region.tsv", sep="\t", index=False)

# --- simulated stratification: two groups, the variant has NO effect on the trait
rng2 = np.random.default_rng(7)
n = 500
grp = np.repeat([1, 2], n)
freq = np.where(grp == 1, 0.2, 0.7)
g = rng2.binomial(2, freq)
trait = 0.8 * (grp == 2) + rng2.normal(0, 1, 2 * n)
pd.DataFrame({"person": [f"S{i+1:04d}" for i in range(2 * n)], "group": grp, "genotype": g,
              "trait": np.round(trait, 3)}).to_csv(out + "stratification.tsv", sep="\t", index=False)
print("done")

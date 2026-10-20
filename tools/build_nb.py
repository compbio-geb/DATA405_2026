import nbformat as nbf, sys
BASE = "https://compbio-geb.github.io/DATA405_2026/data/lab3"
PLINK = "https://s3.amazonaws.com/plink2-assets/alpha7/plink2_linux_x86_64_20261015.zip"
md = nbf.v4.new_markdown_cell; code = nbf.v4.new_code_cell
def py():
    c = []
    c.append(md("""# Lab 3 · Parts B and C — GWAS and polygenic risk score with PLINK
**DATA405/BIOL480 · Genomic Data Science · Fall 2026**

This notebook runs **PLINK 2**, the program most GWAS studies use, on a cohort of 2,500 people.
The people and their LDL-cholesterol values are **simulated** (no real person is in this file), but the files, commands and outputs are exactly what a real study uses.

**How to use it:** run the cells from top to bottom with ▶ (or Shift+Enter). Lines starting with `!` are terminal commands — those are the PLINK commands you would type on a server. Read the handout section for each step before running it, and write down the numbers marked **✏️ Record** for your report."""))
    c.append(md("## Step 0 · Install PLINK 2 and download the data\nPLINK is one small program; we download it and the four data files (about 9 MB)."))
    c.append(code(f"""import os
if not os.path.exists("plink2"):
    !wget -q {PLINK} -O plink2.zip && unzip -o -q plink2.zip && chmod +x plink2
for f in ["cohort.vcf.gz", "pheno.txt", "discovery.txt", "target.txt"]:
    if not os.path.exists(f):
        !wget -q {BASE}/{{f}}
!./plink2 --version
!ls -lh cohort.vcf.gz pheno.txt discovery.txt target.txt"""))
    c.append(md("## Step 1 · Look at the input\nThe cohort arrives as a **VCF** — the same format that Part A (Galaxy) produced, but with 2,500 people as columns instead of 3. We print the first 5 variants and the first 6 people."))
    c.append(code("""!zcat cohort.vcf.gz | grep -v "^##" | head -6 | cut -f1-15
!echo; head -5 pheno.txt; echo ...; wc -l pheno.txt discovery.txt target.txt"""))
    c.append(md("## Step 2 · Convert the VCF to PLINK's own format\nVCF is text and slow. PLINK converts it once into a compact binary file set: `.pgen` (genotypes), `.pvar` (variant list) and `.psam` (sample list)."))
    c.append(code("""!./plink2 --vcf cohort.vcf.gz --make-pgen --out cohort --silent
!grep -E "loaded from|written" cohort.log
!head -3 cohort.pvar; head -3 cohort.psam"""))
    c.append(md("## Step 3 · Quality control (QC)\nFirst we **measure** each problem, then remove it. ✏️ Record how many people and variants fail each check."))
    c.append(code("""!./plink2 --pfile cohort --missing --out miss --silent     # missing calls per person (.smiss) and per variant (.vmiss)
!./plink2 --pfile cohort --freq    --out freq --silent     # allele frequency of each variant (.afreq)
!./plink2 --pfile cohort --hardy   --out hardy --silent    # Hardy-Weinberg test of each variant (.hardy)
import pandas as pd
smiss = pd.read_csv("miss.smiss", sep="\\t"); vmiss = pd.read_csv("miss.vmiss", sep="\\t")
afreq = pd.read_csv("freq.afreq", sep="\\t"); hardy = pd.read_csv("hardy.hardy", sep="\\t")
maf = afreq["ALT_FREQS"].where(afreq["ALT_FREQS"] <= 0.5, 1 - afreq["ALT_FREQS"])
print("people with > 2% missing calls  :", (smiss["F_MISS"] > 0.02).sum())
print("variants with > 2% missing calls:", (vmiss["F_MISS"] > 0.02).sum())
print("variants with MAF < 1%          :", (maf < 0.01).sum())
print("variants with HWE p < 1e-6      :", (hardy["P"] < 1e-6).sum())"""))
    c.append(md("Now apply all four filters in one command. PLINK's log tells you what it removed."))
    c.append(code("""!./plink2 --pfile cohort --mind 0.02 --geno 0.02 --maf 0.01 --hwe 1e-6 --make-pgen --out qc --silent
!grep -E "removed|remaining" qc.log"""))
    c.append(md("## Step 4 · Who is in the cohort? Principal component analysis (PCA)\nPCA summarises each person's genotypes as a few numbers (PC1, PC2, …). People with similar ancestry land close together. We first keep a set of roughly independent variants (**LD pruning**), because blocks of linked variants would otherwise dominate the PCs."))
    c.append(code("""!./plink2 --pfile qc --indep-pairwise 50 5 0.2 --out prune --silent
!grep "removed" prune.log
!./plink2 --pfile qc --extract prune.prune.in --pca 10 --out pca --silent
import matplotlib.pyplot as plt
pcs = pd.read_csv("pca.eigenvec", sep="\\t")
plt.figure(figsize=(5, 4)); plt.scatter(pcs.PC1, pcs.PC2, s=4, alpha=.5)
plt.xlabel("PC1"); plt.ylabel("PC2"); plt.title("Each dot is one person"); plt.show()"""))
    c.append(md("✏️ Record: how many groups do you see in PC1?"))
    c.append(md("## Step 5 · GWAS, first attempt: no ancestry correction\nFor each of ~10,700 variants PLINK fits the same regression you did by hand in Lab 2: LDL = a + β·(number of ALT alleles). We use the 2,000 **discovery** people only (the other 500 are saved for the PRS in Part C)."))
    c.append(code("""!./plink2 --pfile qc --keep discovery.txt --pheno pheno.txt --glm allow-no-covars --out gwas_nopc --silent
!grep -E "remaining after|Results" gwas_nopc.log
import numpy as np
from scipy import stats
def manhattan_qq(file, title):
    d = pd.read_csv(file, sep="\\t"); d = d[d.TEST == "ADD"].dropna(subset=["P"])
    d["x"] = np.arange(len(d)); lp = -np.log10(d.P)
    lam = np.median(stats.chi2.isf(d.P, 1)) / stats.chi2.ppf(0.5, 1)
    fig, ax = plt.subplots(1, 2, figsize=(12, 3.5), gridspec_kw={"width_ratios": [3, 1]})
    ax[0].scatter(d.x, lp, s=4, c=np.where(d["#CHROM"] % 2 == 0, "#1f4e79", "#7fa7d1"))
    ax[0].axhline(-np.log10(5e-8), color="red", lw=1); ax[0].set_ylabel("-log10(p)")
    mids = d.groupby("#CHROM").x.median(); ax[0].set_xticks(mids.values, mids.index, fontsize=7)
    ax[0].set_xlabel("chromosome"); ax[0].set_title(title)
    e = -np.log10((np.arange(1, len(d) + 1) - 0.5) / len(d)); ax[1].scatter(e, np.sort(lp)[::-1], s=4)
    ax[1].plot([0, e.max()], [0, e.max()], "r-", lw=1); ax[1].set_xlabel("expected -log10(p)")
    ax[1].set_title(f"QQ plot, lambda = {lam:.2f}"); plt.tight_layout(); plt.show()
    print("variants with p < 5e-8:", (d.P < 5e-8).sum(), "   lambda:", round(lam, 2))
    return d
d0 = manhattan_qq("gwas_nopc.LDL.glm.linear", "GWAS without PCs")"""))
    c.append(md("✏️ Record the number of significant variants and lambda (λ). λ ≈ 1 means the p-values behave as expected; λ much larger than 1 means *something* is inflating every test."))
    c.append(md("## Step 6 · GWAS done properly: add the PCs as covariates\nSame command, plus `--covar` with the first 5 PCs. Now each variant is compared among people of similar ancestry."))
    c.append(code("""!./plink2 --pfile qc --keep discovery.txt --pheno pheno.txt --covar pca.eigenvec --covar-name PC1-PC5 --glm hide-covar --out gwas_pc --silent
d1 = manhattan_qq("gwas_pc.LDL.glm.linear", "GWAS with 5 PCs")"""))
    c.append(md("A significant **region** usually holds several variants in LD with each other (Lecture 13: *a hit is a region, not a single SNP*). **Clumping** keeps the best variant per region and lists the rest under it."))
    c.append(code("""!./plink2 --pfile qc --keep discovery.txt --clump gwas_pc.LDL.glm.linear --clump-p1 5e-8 --clump-r2 0.1 --clump-kb 250 --out loci --silent
!grep "clumps formed" loci.log
loci = pd.read_csv("loci.clumps", sep="\\t")
loci[["#CHROM", "POS", "ID", "P", "TOTAL"]]"""))
    c.append(md("✏️ Record: how many independent regions (loci) pass 5 × 10⁻⁸, and the lead variant, β and p of the strongest one."))
    c.append(code("""d1.sort_values("P").head(10)[["#CHROM", "POS", "ID", "A1", "A1_FREQ", "BETA", "SE", "P"]]"""))
    c.append(md("## Part C · Polygenic risk score (PRS)\nA PRS adds up a person's ALT alleles, each weighted by its GWAS β:  PRS = Σ βᵢ · (ALT count)ᵢ.\nWe (1) clump the **discovery** GWAS so we keep one variant per region, (2) score the 500 **target** people, who were *not* in the GWAS, using different p-value cut-offs, and (3) ask how much of the LDL differences the score explains."))
    c.append(code("""!./plink2 --pfile qc --keep discovery.txt --clump gwas_pc.LDL.glm.linear --clump-p1 1 --clump-r2 0.1 --clump-kb 250 --out prs_clump --silent
cl = pd.read_csv("prs_clump.clumps", sep="\\t")
d1.set_index("ID").loc[cl.ID, ["A1", "BETA"]].to_csv("weights.txt", sep="\\t")      # ID, effect allele, beta
cl[["ID", "P"]].to_csv("snp_p.txt", sep="\\t", index=False, header=False)
cuts = {"5e-8": 5e-8, "1e-5": 1e-5, "0.001": 1e-3, "0.01": 0.01, "0.05": 0.05, "0.5": 0.5, "1": 1}
open("ranges.txt", "w").write("".join(f"p{k} 0 {v}\\n" for k, v in cuts.items()))
!./plink2 --pfile qc --keep target.txt --score weights.txt 1 2 3 header --q-score-range ranges.txt snp_p.txt --out prs --silent
print(len(cl), "clumped variants available for the score")"""))
    c.append(code("""ph = pd.read_csv("pheno.txt", sep="\\t"); P5 = [f"PC{i}" for i in range(1, 6)]
def r2(y, X):
    X = np.c_[np.ones(len(y)), X]; b = np.linalg.lstsq(X, y, rcond=None)[0]; return 1 - (y - X @ b).var() / y.var()
rows = []
for k in cuts:
    s = pd.read_csv(f"prs.p{k}.sscore", sep="\\t").merge(ph, on="#IID").merge(pcs, on="#IID")
    n_snp = int(s.ALLELE_CT.iloc[0] // 2)
    gain = r2(s.LDL.values, s[P5 + ["SCORE1_AVG"]].values) - r2(s.LDL.values, s[P5].values)
    rows.append({"p cut-off": k, "variants in score": n_snp, "R2 explained by PRS (beyond PCs)": round(gain, 3)})
res = pd.DataFrame(rows); print(res.to_string(index=False))
plt.figure(figsize=(6, 3)); plt.bar(res["p cut-off"], res.iloc[:, 2]); plt.xlabel("p-value cut-off"); plt.ylabel("R² in target"); plt.show()"""))
    c.append(md("✏️ Record the best cut-off and its R². Then look at one person: the top 10% vs the bottom 10% of PRS."))
    c.append(code("""best = res.iloc[res.iloc[:, 2].idxmax()]["p cut-off"]
s = pd.read_csv(f"prs.p{best}.sscore", sep="\\t").merge(ph, on="#IID")
s["decile"] = pd.qcut(s.SCORE1_AVG, 10, labels=range(1, 11))
print(s.groupby("decile", observed=True).LDL.mean().round(1))"""))
    c.append(md("**You are done with the notebook.** Answer the report questions in the handout using your recorded numbers. Export one Manhattan plot (right-click → save image) for the report."))
    return c
nb = nbf.v4.new_notebook(); nb.cells = py()
nb.metadata = {"kernelspec": {"name": "python3", "display_name": "Python 3"}, "colab": {"provenance": []}}
nbf.write(nb, sys.argv[1])

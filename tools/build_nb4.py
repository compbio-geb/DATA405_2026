import nbformat as nbf, sys
BASE = "https://compbio-geb.github.io/DATA405_2026/data/lab4"
md = nbf.v4.new_markdown_cell; code = nbf.v4.new_code_cell
c = []
c.append(md("""# Lab 4 · RNA-seq data analysis — four quantifications of the same reads
**DATA405/BIOL480 · Genomic Data Science · Fall 2026**

**Runtime → Run all** runs everything (about 20 minutes). Read the handout sections while the tools run, and write down the numbers marked **✏️**.

What the notebook does:
1. installs STAR, RSEM, Cufflinks and Salmon (Ubuntu packages, ~2 min);
2. downloads 2,000,000 yeast RNA-seq reads (Nagalakshmi et al. 2008 — the reads behind Qian et al. 2010) and the yeast reference;
3. quantifies every gene four ways: STAR unique-read counts (no EM), RSEM, Cufflinks, Salmon;
4. joins the results with the duplicate-gene table and asks whether the choice of tool changes the paper's conclusion.

Lines starting with `!` are shell commands — exactly what you would type on a Linux server."""))
c.append(md("## Step 1 · Install the tools (≈ 2 min)"))
c.append(code("""%%bash
set -e
apt-get -qq update > /dev/null
apt-get -qq install -y rna-star rsem salmon cufflinks samtools > /dev/null 2>&1
echo "STAR      $(STAR --version)"
echo "RSEM      $(rsem-calculate-expression --version | head -1)"
echo "Salmon    $(salmon --version 2>&1)"
echo "Cufflinks $(cufflinks 2>&1 | head -1)"
echo "samtools  $(samtools --version | head -1)\""""))
c.append(md("## Step 2 · Download the reads and the reference (≈ 1 min)\nReads: 2,000,000 single-end 33-bp reads from SRR002051 (the poly(A) library of Nagalakshmi et al. 2008). Reference: *S. cerevisiae* R64 genome, gene annotation (GTF) and transcript sequences (Ensembl Fungi 63)."))
c.append(code(f"""%%bash
set -e
for f in yeast_reads_2M.fastq.gz yeast_genome_R64.fa yeast_genes_R64.gtf yeast_transcripts_R64.fa yeast_paralogs.tsv; do
  [ -f $f ] || wget -q {BASE}/$f
done
ls -lh yeast_*
echo; echo "first read:"; zcat yeast_reads_2M.fastq.gz | head -4
echo; echo "genes in the GTF: $(awk '$3==\\"gene\\"' yeast_genes_R64.gtf | wc -l)\""""))
c.append(md("""## Step 3 · STAR: map the reads and count unique reads per gene (≈ 10 min)
STAR builds an index of the genome (with the annotated splice junctions), aligns every read, and with `--quantMode TranscriptomeSAM GeneCounts` also writes (a) **counts of uniquely mapped reads per gene** — our no-EM baseline, the same rule Qian et al. used — and (b) the alignments in transcript coordinates, which RSEM needs."""))
c.append(code("""%%bash
set -e
mkdir -p star_index
STAR --runMode genomeGenerate --genomeDir star_index --genomeFastaFiles yeast_genome_R64.fa \\
     --sjdbGTFfile yeast_genes_R64.gtf --sjdbOverhang 32 --genomeSAindexNbases 10 --runThreadN 2 > star_index.log
STAR --genomeDir star_index --readFilesIn yeast_reads_2M.fastq.gz --readFilesCommand zcat --runThreadN 2 \\
     --outSAMtype BAM Unsorted --quantMode TranscriptomeSAM GeneCounts --outFileNamePrefix star_ > star_run.log
grep -E "input reads|Uniquely mapped reads %|multiple loci|too many loci|unmapped: too short" star_Log.final.out
echo; head -8 star_ReadsPerGene.out.tab"""))
c.append(md("✏️ Record: % uniquely mapped, % multi-mapped, % unmapped. The first four rows of the counts table are not genes: `N_unmapped`, `N_multimapping`, `N_noFeature` (mapped outside any gene), `N_ambiguous`. Column 2 = unstranded counts."))
c.append(md("## Step 4 · RSEM: EM on the transcript alignments (≈ 3 min)\nRSEM prepares its own transcript set from genome + GTF, then runs EM over all reads at once (`expected_count` is fractional because shared reads are split). This STAR → RSEM pair is the ENCODE / TCGA standard."))
c.append(code("""%%bash
set -e
mkdir -p rsem_ref
rsem-prepare-reference --gtf yeast_genes_R64.gtf yeast_genome_R64.fa rsem_ref/yeast > rsem_prep.log 2>&1
rsem-calculate-expression --alignments --no-bam-output -p 2 star_Aligned.toTranscriptome.out.bam rsem_ref/yeast rsem > rsem_run.log 2>&1
head -3 rsem.genes.results"""))
c.append(md("## Step 5 · Cufflinks: EM on the genome alignments (≈ 3 min)\nCufflinks needs a position-sorted BAM. By default it splits each multi-mapping read evenly between its positions (a read at 2 places counts ½ each) and then runs EM across isoforms. (Its `-u` 'multi-read correct' option re-weights by unique-read coverage; for near-identical pairs with almost no unique reads it collapses the estimate — try it later as an extra.)"))
c.append(code("""%%bash
set -e
samtools sort -@ 2 -o star_sorted.bam star_Aligned.out.bam 2> /dev/null
cufflinks -q -p 2 -G yeast_genes_R64.gtf -o cufflinks star_sorted.bam > cufflinks.log 2>&1
cut -f1,5,8,10 cufflinks/genes.fpkm_tracking | head -3"""))
c.append(md("## Step 6 · Salmon: EM without alignment (≈ 1 min)\nSalmon indexes the transcript sequences by k-mers (k = 19 here because the reads are only 33 bases), finds the transcripts each read could come from, and runs EM. `-g` adds up transcripts into genes."))
c.append(code("""%%bash
set -e
salmon index -t yeast_transcripts_R64.fa -i salmon_index -k 19 > salmon_index.log 2>&1
salmon quant -i salmon_index -l U -r yeast_reads_2M.fastq.gz -g yeast_genes_R64.gtf -p 2 -o salmon > salmon_run.log 2>&1
grep "Mapping rate" salmon/logs/salmon_quant.log
head -3 salmon/quant.genes.sf"""))
c.append(md("""## Part B · Compare the four answers
Join the four tables with the duplicate-gene table (every yeast gene's closest paralog from Ensembl, with the protein identity between them) into one table `d`, one row per gene."""))
c.append(code("""import pandas as pd, numpy as np
import matplotlib.pyplot as plt

star = pd.read_csv("star_ReadsPerGene.out.tab", sep="\\t", header=None, names=["gene", "count", "fwd", "rev"])
special = star[star.gene.str.startswith("N_")]
star = star[~star.gene.str.startswith("N_")][["gene", "count"]].rename(columns={"count": "STAR"})
rsem = (pd.read_csv("rsem.genes.results", sep="\\t")[["gene_id", "expected_count", "TPM"]]
          .rename(columns={"gene_id": "gene", "expected_count": "RSEM", "TPM": "RSEM_TPM"}))
cuff = (pd.read_csv("cufflinks/genes.fpkm_tracking", sep="\\t")[["tracking_id", "FPKM"]]
          .rename(columns={"tracking_id": "gene", "FPKM": "Cufflinks_FPKM"}))
salmon = (pd.read_csv("salmon/quant.genes.sf", sep="\\t")[["Name", "NumReads", "TPM"]]
            .rename(columns={"Name": "gene", "NumReads": "Salmon", "TPM": "Salmon_TPM"}))
par = pd.read_csv("yeast_paralogs.tsv", sep="\\t")

d = (star.merge(rsem, on="gene").merge(cuff, on="gene", how="left").merge(salmon, on="gene", how="left")
         .merge(par[["gene", "name", "closest_paralog", "min_pct_id"]], on="gene", how="left"))
d["identity"] = pd.cut(d.min_pct_id.fillna(0), [-1, 30, 70, 95, 101],
                       labels=["no close paralog", "30-70 %", "70-95 %", "95-100 %"])
print(special.to_string(index=False))
print(len(d), "genes;", d.min_pct_id.notna().sum(), "have a paralog")"""))
c.append(md("### Exercise 1 · How many reads did each tool use?"))
c.append(code("""print("STAR, unique reads in genes:", int(d.STAR.sum()))
print("RSEM, expected counts      :", int(d.RSEM.sum()))
print("Salmon, reads              :", int(d.Salmon.sum()))
print("multi-mapping reads STAR threw away:", int(special.loc[special.gene == "N_multimapping", "count"].iloc[0]))"""))
c.append(md("### Exercise 2 · One duplicate pair: TEF1 and TEF2\n*TEF1* (YPR080W) and *TEF2* (YBR118W): proteins 100 % identical. Then *TDH3*/*TDH2* (96 %) and *ADH1*/*ADH2* (93 %). ✏️ Record the four values for TEF1, TEF2 and TDH3."))
c.append(code("""cols = ["gene", "name", "closest_paralog", "min_pct_id", "STAR", "RSEM", "Cufflinks_FPKM", "Salmon"]
pairs = ["YPR080W", "YBR118W", "YGR192C", "YJR009C", "YOL086C", "YMR303C"]
d[d.gene.isin(pairs)][cols].round(1)"""))
c.append(md("### Exercise 3 · From one pair to all genes\n**gain** = RSEM ÷ STAR unique count, for genes with ≥ 20 RSEM reads, by identity to the closest paralog. ✏️ Record the 95–100 % row. Save the box plot for the report (right-click → save image)."))
c.append(code("""e = d[d.RSEM >= 20].copy()
e["gain"] = e.RSEM / e.STAR.clip(lower=1)
summary = e.groupby("identity", observed=True).gain.agg(
    n="size", median="median", q75=lambda x: x.quantile(0.75), share_over_1_5=lambda x: (x > 1.5).mean())
print(summary.round(2))
fig, ax = plt.subplots(figsize=(6, 4))
e.boxplot(column="gain", by="identity", ax=ax, showfliers=False)
ax.set_yscale("log"); ax.set_ylabel("RSEM / STAR unique"); ax.set_xlabel("protein identity to closest paralog")
ax.set_title(""); fig.suptitle(""); plt.show()"""))
c.append(md("### Exercise 4 · The paper's test, with and without EM\nFor each duplicate pair (≥ 40 RSEM reads together): **lost** = 1 − (STAR₁ + STAR₂) / (RSEM₁ + RSEM₂), the share of the pair's expression that disappears under unique-only counting. Qian et al.'s excess of reduced pairs was 31 %. ✏️ Record the share of pairs losing more than 31 % in each bin."))
c.append(code("""p = d.dropna(subset=["closest_paralog"]).merge(
        d[["gene", "STAR", "RSEM"]].rename(columns={"gene": "closest_paralog", "STAR": "STAR2", "RSEM": "RSEM2"}),
        on="closest_paralog")
p = p[p.gene < p.closest_paralog]                       # each pair once
p = p[(p.RSEM + p.RSEM2) >= 40]
p["lost"] = 1 - (p.STAR + p.STAR2) / (p.RSEM + p.RSEM2)
p["identity"] = pd.cut(p.min_pct_id, [0, 70, 95, 101], labels=["<70 %", "70-95 %", "95-100 %"])
print(p.groupby("identity", observed=True).lost.agg(n="size", median="median", share_over_31=lambda x: (x > 0.31).mean()).round(2))
p.sort_values("lost", ascending=False)[["gene", "name", "closest_paralog", "min_pct_id", "STAR", "STAR2", "RSEM", "RSEM2", "lost"]].head(8).round(2)"""))
c.append(md("### Exercise 5 · Do the EM tools agree?"))
c.append(code("""fig, axes = plt.subplots(1, 3, figsize=(15, 4.3))
close = e.identity == "95-100 %"
for ax, x, y, xl, yl in [(axes[0], "STAR", "RSEM", "STAR unique reads", "RSEM expected count"),
                         (axes[1], "Salmon", "RSEM", "Salmon reads", "RSEM expected count"),
                         (axes[2], "Cufflinks_FPKM", "RSEM_TPM", "Cufflinks FPKM", "RSEM TPM")]:
    ax.scatter(e[x][~close] + 1, e[y][~close] + 1, s=4, alpha=.3, color="grey")
    ax.scatter(e[x][close] + 1, e[y][close] + 1, s=8, color="#C00000", label="paralog 95-100 % identical")
    ax.plot([1, 3e4], [1, 3e4], "k--", lw=1); ax.set_xscale("log"); ax.set_yscale("log")
    r = np.corrcoef(np.log10(e[x].fillna(0) + 1), np.log10(e[y] + 1))[0, 1]
    ax.set_xlabel(xl); ax.set_ylabel(yl); ax.set_title(f"r = {r:.3f} (log scale)"); ax.legend(loc="upper left", fontsize=8)
plt.tight_layout(); plt.show()
d.to_csv("lab4_four_quantifications.tsv", sep="\\t", index=False)
print("saved lab4_four_quantifications.tsv — download it from the Files panel if you want to keep it")"""))
c.append(md("**Done.** Answer Q1–Q7 in the handout with your recorded numbers; include the box plot (Ex 3) and one scatter (Ex 5).\n\n*Extra (optional):* rerun Step 5 with `-u` added to the cufflinks command and look at TEF1 again."))
nb = nbf.v4.new_notebook(); nb.cells = c
nb.metadata = {"kernelspec": {"name": "python3", "display_name": "Python 3"}, "colab": {"provenance": []}}
nbf.write(nb, sys.argv[1])

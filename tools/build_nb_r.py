import nbformat as nbf, sys
BASE = "https://compbio-geb.github.io/DATA405_2026/data/lab3"
PLINK = "https://s3.amazonaws.com/plink2-assets/alpha7/plink2_linux_x86_64_20261015.zip"
md = nbf.v4.new_markdown_cell; code = nbf.v4.new_code_cell
c = []
c.append(md("""# Lab 3 · Parts B and C — GWAS and polygenic risk score with PLINK (R version)
**DATA405/BIOL480 · Genomic Data Science · Fall 2026**

**First:** in Colab choose **Runtime → Change runtime type → R**, then Save.

This notebook runs **PLINK 2**, the program most GWAS studies use, on a cohort of 2,500 people. The people and their LDL-cholesterol values are **simulated** (no real person is in this file), but the files, commands and outputs are exactly what a real study uses.

Run the cells from top to bottom. `sh("...")` runs a terminal command — those are the PLINK commands you would type on a server. Write down the numbers marked **✏️ Record** for your report."""))
c.append(md("## Step 0 · Install PLINK 2 and download the data"))
c.append(code(f"""sh <- function(cmd) cat(system(cmd, intern = TRUE), sep = "\\n")
if (!file.exists("plink2")) {{
  download.file("{PLINK}", "plink2.zip", quiet = TRUE); unzip("plink2.zip"); Sys.chmod("plink2", "755")
}}
for (f in c("cohort.vcf.gz", "pheno.txt", "discovery.txt", "target.txt"))
  if (!file.exists(f)) download.file(paste0("{BASE}/", f), f, quiet = TRUE)
sh("./plink2 --version")
sh("ls -lh cohort.vcf.gz pheno.txt discovery.txt target.txt")"""))
c.append(md("## Step 1 · Look at the input\nThe cohort arrives as a **VCF** — the same format Part A (Galaxy) produced, but with 2,500 people as columns instead of 3."))
c.append(code("""sh("zcat cohort.vcf.gz | grep -v '^##' | head -6 | cut -f1-15")
head(read.delim("pheno.txt"), 4)"""))
c.append(md("## Step 2 · Convert the VCF to PLINK's own format\n`.pgen` (genotypes), `.pvar` (variant list), `.psam` (sample list)."))
c.append(code("""sh("./plink2 --vcf cohort.vcf.gz --make-pgen --out cohort --silent")
sh("grep -E 'loaded from|written' cohort.log")
sh("head -3 cohort.pvar; head -3 cohort.psam")"""))
c.append(md("## Step 3 · Quality control (QC)\nMeasure each problem first, then remove it. ✏️ Record how many fail each check."))
c.append(code("""sh("./plink2 --pfile cohort --missing --out miss --silent")   # missing calls per person / variant
sh("./plink2 --pfile cohort --freq    --out freq --silent")   # allele frequencies
sh("./plink2 --pfile cohort --hardy   --out hardy --silent")  # Hardy-Weinberg test
smiss <- read.delim("miss.smiss"); vmiss <- read.delim("miss.vmiss")
af <- read.delim("freq.afreq")$ALT_FREQS; maf <- pmin(af, 1 - af); hardy <- read.delim("hardy.hardy")
cat("people with > 2% missing calls  :", sum(smiss$F_MISS > 0.02), "\\n")
cat("variants with > 2% missing calls:", sum(vmiss$F_MISS > 0.02), "\\n")
cat("variants with MAF < 1%          :", sum(maf < 0.01), "\\n")
cat("variants with HWE p < 1e-6      :", sum(hardy$P < 1e-6), "\\n")"""))
c.append(code("""sh("./plink2 --pfile cohort --mind 0.02 --geno 0.02 --maf 0.01 --hwe 1e-6 --make-pgen --out qc --silent")
sh("grep -E 'removed|remaining' qc.log")"""))
c.append(md("## Step 4 · Who is in the cohort? PCA\nLD pruning first, then PCA on the pruned variants."))
c.append(code("""sh("./plink2 --pfile qc --indep-pairwise 50 5 0.2 --out prune --silent")
sh("grep removed prune.log")
sh("./plink2 --pfile qc --extract prune.prune.in --pca 10 --out pca --silent")
pcs <- read.delim("pca.eigenvec", check.names = FALSE)
plot(pcs$PC1, pcs$PC2, pch = 16, cex = .4, xlab = "PC1", ylab = "PC2", main = "Each dot is one person")"""))
c.append(md("✏️ Record: how many groups do you see in PC1?"))
c.append(md("## Step 5 · GWAS, first attempt: no ancestry correction\nFor each variant PLINK fits LDL = a + β·(ALT count), in the 2,000 **discovery** people."))
c.append(code("""sh("./plink2 --pfile qc --keep discovery.txt --pheno pheno.txt --glm allow-no-covars --out gwas_nopc --silent")
sh("grep -E 'remaining after|Results' gwas_nopc.log")
manhattan_qq <- function(file, title) {
  d <- read.delim(file, check.names = FALSE); d <- d[d$TEST == "ADD" & !is.na(d$P), ]
  lp <- -log10(d$P); lam <- median(qchisq(d$P, 1, lower.tail = FALSE)) / qchisq(0.5, 1)
  layout(matrix(c(1, 1, 1, 2), 1)); par(mar = c(4, 4, 2, 1))
  plot(seq_along(lp), lp, pch = 16, cex = .4, col = ifelse(d[["#CHROM"]] %% 2 == 0, "#1f4e79", "#7fa7d1"),
       xaxt = "n", xlab = "chromosome", ylab = "-log10(p)", main = title)
  mids <- tapply(seq_along(lp), d[["#CHROM"]], median); axis(1, mids, names(mids), cex.axis = .6)
  abline(h = -log10(5e-8), col = "red")
  e <- -log10((seq_along(lp) - 0.5) / length(lp))
  plot(e, sort(lp, decreasing = TRUE), pch = 16, cex = .4, xlab = "expected -log10(p)", ylab = "observed",
       main = sprintf("QQ, lambda = %.2f", lam)); abline(0, 1, col = "red"); layout(1)
  cat("variants with p < 5e-8:", sum(d$P < 5e-8), "   lambda:", round(lam, 2), "\\n"); invisible(d)
}
d0 <- manhattan_qq("gwas_nopc.LDL.glm.linear", "GWAS without PCs")"""))
c.append(md("✏️ Record the number of significant variants and λ (lambda)."))
c.append(md("## Step 6 · GWAS done properly: add the PCs as covariates"))
c.append(code("""sh("./plink2 --pfile qc --keep discovery.txt --pheno pheno.txt --covar pca.eigenvec --covar-name PC1-PC5 --glm hide-covar --out gwas_pc --silent")
d1 <- manhattan_qq("gwas_pc.LDL.glm.linear", "GWAS with 5 PCs")"""))
c.append(md("**Clumping** keeps the best variant per region (a hit is a region, not a single SNP)."))
c.append(code("""sh("./plink2 --pfile qc --keep discovery.txt --clump gwas_pc.LDL.glm.linear --clump-p1 5e-8 --clump-r2 0.1 --clump-kb 250 --out loci --silent")
loci <- read.delim("loci.clumps", check.names = FALSE); loci[, c("#CHROM", "POS", "ID", "P", "TOTAL")]
head(d1[order(d1$P), c("#CHROM", "POS", "ID", "A1", "A1_FREQ", "BETA", "SE", "P")], 10)"""))
c.append(md("✏️ Record the number of loci and the lead variant, β and p of the strongest one."))
c.append(md("## Part C · Polygenic risk score (PRS)\nPRS = Σ βᵢ · (ALT count)ᵢ. Clump the discovery GWAS, score the 500 **target** people (not in the GWAS) at several p-value cut-offs, and measure how much of LDL the score explains."))
c.append(code("""sh("./plink2 --pfile qc --keep discovery.txt --clump gwas_pc.LDL.glm.linear --clump-p1 1 --clump-r2 0.1 --clump-kb 250 --out prs_clump --silent")
cl <- read.delim("prs_clump.clumps", check.names = FALSE)
w <- d1[match(cl$ID, d1$ID), c("ID", "A1", "BETA")]
write.table(w, "weights.txt", sep = "\\t", quote = FALSE, row.names = FALSE)
write.table(cl[, c("ID", "P")], "snp_p.txt", sep = "\\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
cuts <- c("5e-8" = 5e-8, "1e-5" = 1e-5, "0.001" = 1e-3, "0.01" = 0.01, "0.05" = 0.05, "0.5" = 0.5, "1" = 1)
writeLines(paste0("p", names(cuts), " 0 ", cuts), "ranges.txt")
sh("./plink2 --pfile qc --keep target.txt --score weights.txt 1 2 3 header --q-score-range ranges.txt snp_p.txt --out prs --silent")
cat(nrow(cl), "clumped variants available for the score\\n")"""))
c.append(code("""ph <- read.delim("pheno.txt", check.names = FALSE); P5 <- paste0("PC", 1:5)
res <- data.frame()
for (k in names(cuts)) {
  s <- read.delim(paste0("prs.p", k, ".sscore"), check.names = FALSE)
  s <- merge(merge(s, ph, by = "#IID"), pcs, by = "#IID")
  base <- summary(lm(reformulate(P5, "LDL"), s))$r.squared
  full <- summary(lm(reformulate(c(P5, "SCORE1_AVG"), "LDL"), s))$r.squared
  res <- rbind(res, data.frame(cutoff = k, variants = s$ALLELE_CT[1] / 2, R2_beyond_PCs = round(full - base, 3)))
}
print(res, row.names = FALSE)
barplot(res$R2_beyond_PCs, names.arg = res$cutoff, xlab = "p-value cut-off", ylab = "R² in target")"""))
c.append(md("✏️ Record the best cut-off and its R². Then compare the top and bottom 10% of PRS."))
c.append(code("""best <- res$cutoff[which.max(res$R2_beyond_PCs)]
s <- merge(read.delim(paste0("prs.p", best, ".sscore"), check.names = FALSE), ph, by = "#IID")
s$decile <- cut(rank(s$SCORE1_AVG), 10, labels = 1:10)
round(tapply(s$LDL, s$decile, mean), 1)"""))
c.append(md("**You are done with the notebook.** Answer the report questions in the handout. Save one Manhattan plot (right-click → save image) for the report."))
nb = nbf.v4.new_notebook(); nb.cells = c
nb.metadata = {"kernelspec": {"name": "ir", "display_name": "R"}, "language_info": {"name": "R"}, "colab": {"provenance": []}}
nbf.write(nb, sys.argv[1])

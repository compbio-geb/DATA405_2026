# Lab 4 cross-validation: local pipeline vs Galaxy (Oct 27, 2026)
Input: yeast_reads_2M.fastq.gz (SRR002051 subsample, seed 405), Ensembl Fungi 63 R64 genome/GTF, RSEM-derived transcripts.
| Quantity | Local (Ubuntu 24.04 apt / GitHub binaries) | usegalaxy.org (tool-shed wrappers) |
|---|---|---|
| STAR 2.7.11b uniquely mapped | 1,025,300 (51.27 %) | 1,025,300 (51.27 %) |
| STAR multi-mapped | 135,672 (6.78 %) | 135,672 (6.78 %) |
| STAR unmapped too short | 785,502 (39.28 %) | 785,502 (39.28 %) |
| STAR unique count TEF1 / TEF2 / TDH3 / ADH1 | 768 / 859 / 9,111 / 8,178 | 768 / 859 / 9,111 / 8,178 |
| Salmon 1.10 (k 19, -l U, -g GTF) NumReads TEF1/TEF2/TDH3/ADH1 | 5,584.6 / 4,614.4 / 11,719.2 / 6,129.4 | 5,584.6 / 4,614.4 / 11,719.2 / 6,129.44 |
| Cufflinks 2.2.1 default FPKM TEF1/TEF2/TDH3/ADH1 | 4,956.5 / 4,930.8 / 13,749.7 / 8,178.8 | 4,956.52 / 4,930.74 / 13,749.7 / 8,178.77 |
| Cufflinks 2.2.1 with -u (multi-read correct) TEF1/TEF2 | 250.4 / 222.5 | 250.444 / 222.453 |
| RSEM 1.3.x expected_count TEF1/TEF2 | 7,006.5 / 6,217.5 | (RSEM only on usegalaxy.eu; uploads pending) |
Galaxy history (anonymous, usegalaxy.org): bbd44e69cb8906b5c9d2bf58a978ae76, datasets 20–36.

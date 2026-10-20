import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
BLUE="#0070C0"; GREY="#7F7F7F"; ORANGE="#E97132"; GREEN="#00A651"; PURPLE="#7030A0"
# ---------- pipeline ----------
fig, ax = plt.subplots(figsize=(11, 4.4)); ax.set_xlim(0, 11); ax.set_ylim(0, 4.4); ax.axis("off")
def box(x, y, w, h, title, sub, fc, ec):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.02,rounding_size=0.08", fc=fc, ec=ec, lw=1.5))
    ax.text(x + w/2, y + h*0.64, title, ha="center", va="center", fontsize=9.6, weight="bold")
    ax.text(x + w/2, y + h*0.28, sub, ha="center", va="center", fontsize=8.3, color="#333333")
def arrow(x1, y1, x2, y2):
    ax.annotate("", (x2, y2), (x1, y1), arrowprops=dict(arrowstyle="-|>", color="#444444", lw=1.3))
ax.text(0.05, 4.2, "ONE PERSON: reads → variants", fontsize=10, weight="bold", color=GREY)
row1 = [("Sequence", "FASTQ reads\n(Illumina)", "done"), ("Quality check", "FastQC", "done"), ("Map reads", "BWA-MEM → BAM", "done"),
        ("Duplicates", "mark them\n(Picard / GATK)", "done"), ("Call variants", "FreeBayes\n(or GATK)", "A"), ("Filter", "QUAL, depth", "A"), ("Annotate", "SnpEff", "A")]
w, h, gap, y = 1.44, 1.05, 0.12, 2.95
for i, (t, s, k) in enumerate(row1):
    x = 0.05 + i * (w + gap)
    fc, ec = ("#F2F2F2", "#BFBFBF") if k == "done" else ("#DCEBF7", BLUE)
    box(x, y, w, h, t, s, fc, ec)
    if i: arrow(x - gap + 0.01, y + h/2, x - 0.01, y + h/2)
ax.text(0.05 + 2*(w+gap) - 0.3, y - 0.22, "already done for you: the BAM file is your input", fontsize=8.5, color=GREY, ha="center")
ax.text(0.05 + 5*(w+gap) + w/2, y - 0.22, "Part A · Galaxy", fontsize=9.5, color=BLUE, weight="bold", ha="center")
ax.text(0.05, 1.95, "MANY PEOPLE: variants → association → score", fontsize=10, weight="bold", color=GREY)
row2 = [("Cohort VCF", "2,500 people", "in"), ("QC", "missing, MAF, HWE", "B"), ("PCA", "ancestry", "B"), ("GWAS", "PLINK --glm", "B"), ("Clump", "one SNP per region", "C"), ("PRS", "PLINK --score", "C")]
w2, y2 = 1.62, 0.55
for i, (t, s, k) in enumerate(row2):
    x = 0.05 + i * (w2 + 0.2)
    fc, ec = {"in": ("#F2F2F2", "#BFBFBF"), "B": ("#FBE5D6", ORANGE), "C": ("#EADCF4", PURPLE)}[k]
    box(x, y2, w2, h, t, s, fc, ec)
    if i: arrow(x - 0.19, y2 + h/2, x - 0.01, y2 + h/2)
ax.text(0.05 + 2*(w2+0.2) + w2/2, y2 - 0.25, "Part B · PLINK in Colab", fontsize=9.5, color=ORANGE, weight="bold", ha="center")
ax.text(0.05 + 4.5*(w2+0.2) + w2/2 - 0.1, y2 - 0.25, "Part C · PLINK in Colab", fontsize=9.5, color=PURPLE, weight="bold", ha="center")
plt.savefig("/root/lab3/figs/pipeline.png", dpi=200, bbox_inches="tight"); plt.close()

# ---------- pileup ----------
ref = "ACGTTGCATGCCAGTAAGCT"
het, err, hom = 5, 11, 17
reads = []
import random; random.seed(3)
alt = {het: "T", err: "A", hom: "C"}
for r in range(10):
    s = list(ref)
    if r % 2 == 0: s[het] = "T"
    if r == 6: s[err] = "A"
    s[hom] = "C"
    start = random.randint(0, 3); end = random.randint(len(ref) - 2, len(ref))
    reads.append((start, end, s))
col = {"A": GREEN, "C": BLUE, "G": "#C55A11", "T": "#C00000"}
fig, ax = plt.subplots(figsize=(9, 5.2)); ax.axis("off"); ax.set_xlim(-3.2, len(ref) + 0.3); ax.set_ylim(-12.4, 2.6)
for j in (het, err, hom):
    ax.add_patch(plt.Rectangle((j - 0.45, -10.6), 0.9, 11.6, fc="#FFF2CC", ec="none"))
ax.text(-3.1, 0, "reference", fontsize=10, va="center", weight="bold")
for j, b in enumerate(ref): ax.text(j, 0, b, ha="center", va="center", fontsize=13, family="monospace", weight="bold", color=col[b])
for i, (a, b_, s) in enumerate(reads):
    yy = -1.0 - i
    ax.text(-3.1, yy, f"read {i+1}", fontsize=9, va="center", color=GREY)
    for j in range(a, b_):
        if s[j] == ref[j]: ax.text(j, yy, "·", ha="center", va="center", fontsize=13, family="monospace", color="#A6A6A6")
        else: ax.text(j, yy, s[j], ha="center", va="center", fontsize=13, family="monospace", weight="bold", color=col[s[j]])
lab = {het: "5 of 10 reads say T\n→ heterozygous (0/1)", err: "1 of 10 reads says A\n→ probably an error", hom: "10 of 10 reads say C\n→ homozygous ALT (1/1)"}
for j, t in lab.items(): ax.text(j, -11.6, t, ha="center", va="center", fontsize=8.6, color="#333333")
ax.text(len(ref)/2 - 1, 2.0, "A pileup: reads stacked under the reference ( · = same letter as the reference)", ha="center", fontsize=9.6, weight="bold")
plt.savefig("/root/lab3/figs/pileup.png", dpi=200, bbox_inches="tight"); plt.close()
print("ok")

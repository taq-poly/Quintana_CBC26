#!/usr/bin/env python3
"""Per-candidate metrics from split competitive mapping (adults).

microbes: per genome (contig prefix 'tag|'): assigned read pairs, ppm of trimmed pairs, covered bp,
          breadth, mean depth, expected breadth for that depth under even coverage (1-exp(-depth)),
          and breadth/expected (<<1 = reads piled on a few loci: conserved genes/rRNA or expression).
dog/human: assigned pairs, ppm, reads on mtDNA, and number of distinct loci (reads merged within 1 kb).
"""
import csv, math, collections, subprocess, os
B = "/workspace/toxocara_screen"
S = {"SRR26866479": "F_T1", "SRR26866491": "F_T2", "SRR26866490": "F_T3",
     "SRR26866487": "M_T10", "SRR26866489": "M_T4", "SRR26866488": "M_T6"}
MITO = {"dog": "NC_002008.4", "human": "NC_012920.1"}
out = f"{B}/results/summary"; os.makedirs(out, exist_ok=True)
mic_rows, mam_rows = [], []
for r, lab in S.items():
    d = f"{B}/results/{r}/compete"
    trim = int(list(csv.reader(open(f"{B}/results/{r}/counts.tsv"), delimiter="\t"))[1][2])
    # microbes per genome
    g = collections.defaultdict(lambda: [0, 0, 0, 0.0])  # reads, len, covbases, depth*len
    for row in csv.DictReader(open(f"{d}/assigned_microbes.cov.tsv"), delimiter="\t"):
        tag = row["#rname"].split("|")[0]; ln = int(row["endpos"])
        g[tag][0] += int(row["numreads"]); g[tag][1] += ln; g[tag][2] += int(row["covbases"]); g[tag][3] += float(row["meandepth"]) * ln
    for tag, (n, ln, cb, dl) in g.items():
        if n == 0: continue
        depth = dl / ln; br = cb / ln; exp = 1 - math.exp(-depth)
        mic_rows.append([lab, tag, n // 2, f"{1e6*(n/2)/trim:.2f}", cb, f"{100*br:.3f}", f"{depth:.4f}",
                         f"{100*exp:.3f}", f"{br/exp:.2f}" if exp > 0 else "NA"])
    # mammals
    for t in ("dog", "human"):
        bam = f"{d}/assigned_{t}.bam"
        reads = subprocess.run(["samtools", "view", "-c", bam], capture_output=True, text=True).stdout.strip()
        mt = subprocess.run(["samtools", "view", "-c", bam, MITO[t]], capture_output=True, text=True).stdout.strip()
        loci = 0; last = (None, -10**9)
        for l in subprocess.run(["samtools", "view", bam], capture_output=True, text=True).stdout.splitlines():
            f = l.split("\t", 4); c, p = f[2], int(f[3])
            if c != last[0] or p - last[1] > 1000: loci += 1
            last = (c, p)
        n = int(reads) // 2
        mam_rows.append([lab, t, n, f"{1e6*n/trim:.2f}", int(mt) // 2, loci])
with open(f"{out}/adult_competitive_microbes.tsv", "w") as o:
    o.write("sample\tgenome\tread_pairs\tppm_trimmed\tcovered_bp\tbreadth_pct\tmean_depth\texpected_breadth_pct\tbreadth_obs_over_exp\n")
    for r in mic_rows: o.write("\t".join(map(str, r)) + "\n")
with open(f"{out}/adult_competitive_mammals.tsv", "w") as o:
    o.write("sample\ttarget\tread_pairs\tppm_trimmed\tmt_read_pairs\tdistinct_loci_1kb\n")
    for r in mam_rows: o.write("\t".join(map(str, r)) + "\n")
print("written")

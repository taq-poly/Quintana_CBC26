#!/usr/bin/env python3
"""Cross-sample summary of pass-1 screening (Kraken2 genus counts + sendsketch hits) for the adults."""
import csv, collections, os
B = "/workspace/toxocara_screen"
S = {"SRR26866479": "F_T1", "SRR26866491": "F_T2", "SRR26866490": "F_T3",
     "SRR26866487": "M_T10", "SRR26866489": "M_T4", "SRR26866488": "M_T6"}
out = os.path.join(B, "results", "summary")
os.makedirs(out, exist_ok=True)

counts, genus, species = {}, collections.defaultdict(dict), collections.defaultdict(dict)
for r, lab in S.items():
    d = os.path.join(B, "results", r)
    row = list(csv.reader(open(f"{d}/counts.tsv"), delimiter="\t"))[1]
    hm = int(open(f"{d}/half_mapped_reads.txt").read())
    counts[lab] = dict(trim=int(row[2]), nonworm=int(row[3]), half=hm)
    for l in open(f"{d}/kraken.report"):
        f = l.rstrip("\n").split("\t")
        clade, rank, name = int(f[1]), f[3], f[5].strip()
        if rank == "G": genus[name][lab] = clade
        if rank == "S": species[name][lab] = clade

labs = list(S.values())
with open(f"{out}/adult_counts.tsv", "w") as o:
    o.write("sample\ttrim_pairs\tnonworm_pairs\tnonworm_pct\thalf_mapped_reads\n")
    for lab in labs:
        c = counts[lab]
        o.write(f"{lab}\t{c['trim']}\t{c['nonworm']}\t{100*c['nonworm']/c['trim']:.2f}\t{c['half']}\n")

def table(dct, fn, minsum=50):
    rows = sorted(((n, v) for n, v in dct.items() if sum(v.values()) >= minsum),
                  key=lambda x: -sum(x[1].values()))
    with open(f"{out}/{fn}", "w") as o:
        o.write("taxon\t" + "\t".join(labs) + "\tn_samples\tsum\tppm_of_trimmed_pairs_mean\n")
        for n, v in rows:
            ppm = sum(1e6 * v.get(l, 0) / counts[l]["trim"] for l in labs) / len(labs)
            o.write(f"{n}\t" + "\t".join(str(v.get(l, 0)) for l in labs)
                    + f"\t{sum(1 for l in labs if v.get(l,0)>=5)}\t{sum(v.values())}\t{ppm:.2f}\n")
table(genus, "adult_kraken_genus.tsv")
table(species, "adult_kraken_species.tsv", 30)

# sendsketch nt: union of hits, WKID per sample
ss = collections.defaultdict(dict)
for r, lab in S.items():
    p = os.path.join(B, "results", r, "sendsketch_nt.txt")
    for l in open(p):
        if l.startswith("#") or not l.strip(): continue
        f = l.rstrip("\n").split("\t")
        ss[f[1]][lab] = (float(f[2]), float(f[10]))  # ANI, WKID
with open(f"{out}/adult_sendsketch_nt.tsv", "w") as o:
    o.write("ref\tn_samples\tmax_ANI\t" + "\t".join(f"{l}_WKID" for l in labs) + "\n")
    for n, v in sorted(ss.items(), key=lambda x: -len(x[1])):
        o.write(f"{n}\t{len(v)}\t{max(a for a, _ in v.values()):.1f}\t"
                + "\t".join(f"{v[l][1]:.3f}" if l in v else "-" for l in labs) + "\n")
print("written to", out)

#!/usr/bin/env python3
"""Bleed-through test for viruses: per-site base frequencies from samtools mpileup (MAPQ>=20, BQ>=20).
For each virus, sites where well-covered samples have different major alleles are 'diagnostic';
low-level samples are then scored by which sample's alleles their reads carry (or a distinct allele).
"""
import subprocess, collections, sys, re
B = "/workspace/toxocara_screen"
S = {"SRR26866479": "F_T1", "SRR26866491": "F_T2", "SRR26866490": "F_T3",
     "SRR26866487": "M_T10", "SRR26866489": "M_T4", "SRR26866488": "M_T6"}
VIR = {"OZ077912.1": "TcLA_seg1", "OZ077907.1": "RV1", "OZ077897.1": "RV2",
       "CAXMVU010000001.1": "RV3_c1", "CAXMVU010000002.1": "RV3_c2",
       "CAXMVJ010000001.1": "RV4_c1", "CAXMVJ010000002.1": "RV4_c2"}
def pileup(bam, ctg):
    out = subprocess.run(["samtools", "mpileup", "-q", "20", "-Q", "20", "-d", "0", "-B", "-r", ctg,
                          "-f", f"{B}/refs/viruses/all_tcanis_viruses.fa", bam],
                         capture_output=True, text=True).stdout
    res = {}
    for l in out.splitlines():
        f = l.split("\t"); pos, ref, bases = int(f[1]), f[2].upper(), f[4]
        bases = re.sub(r"\^.|\$", "", bases); bases = re.sub(r"[+-](\d+)", lambda m: "#" * 0, bases)
        # remove indel sequences properly
        b2, i = [], 0
        raw = f[4]; raw = re.sub(r"\^.", "", raw).replace("$", "")
        while i < len(raw):
            c = raw[i]
            if c in "+-":
                j = i + 1; n = ""
                while raw[j].isdigit(): n += raw[j]; j += 1
                i = j + int(n); continue
            b2.append(ref if c in ".," else c.upper()); i += 1
        cnt = collections.Counter(x for x in b2 if x in "ACGT")
        if cnt: res[pos] = cnt
    return res
for ctg, vn in VIR.items():
    P = {lab: pileup(f"{B}/results/viruses/{r}.tcv.bam", ctg) for r, lab in S.items()}
    depth = {lab: sum(sum(c.values()) for c in P[lab].values()) for lab in P}
    # sources: samples with >=10x mean depth at >=80% of positions they cover
    major = {}
    for lab in P:
        major[lab] = {p: c.most_common(1)[0][0] for p, c in P[lab].items() if sum(c.values()) >= 10 and c.most_common(1)[0][1] / sum(c.values()) >= 0.8}
    srcs = [lab for lab in P if len(major[lab]) >= 100]
    if len(srcs) < 1: continue
    print(f"\n== {vn}: well-covered samples: {srcs}")
    # pairwise differences among sources
    for i in range(len(srcs)):
        for j in range(i + 1, len(srcs)):
            a, b = srcs[i], srcs[j]
            sh = set(major[a]) & set(major[b]); d = sum(1 for p in sh if major[a][p] != major[b][p])
            print(f"   {a} vs {b}: {d} differing sites of {len(sh)} shared ({100*d/max(1,len(sh)):.2f}%)")
    # score every sample (incl. low-level) against each source's alleles at sites where sources differ
    for lab in P:
        n = sum(sum(c.values()) for c in P[lab].values())
        if n == 0: continue
        score = collections.Counter(); tot = 0
        diag = set.intersection(*[set(major[s]) for s in srcs]) if len(srcs) > 1 else set(major[srcs[0]])
        if len(srcs) > 1: diag = {p for p in diag if len({major[s][p] for s in srcs}) > 1}
        for p in diag:
            if p not in P[lab]: continue
            for base, k in P[lab][p].items():
                hit = [s for s in srcs if major[s][p] == base and s != lab]
                score[hit[0] if len(hit) == 1 else ("shared" if hit else "other")] += k; tot += k
        if len(srcs) == 1:  # single source: count identity to it at all sites
            agree = sum(P[lab][p].get(major[srcs[0]][p], 0) for p in major[srcs[0]] if p in P[lab])
            tot1 = sum(sum(P[lab][p].values()) for p in major[srcs[0]] if p in P[lab])
            print(f"   {lab:6} bases={n:>9}  identity to {srcs[0]} major alleles: {100*agree/max(1,tot1):.2f}% over {tot1} bases")
        else:
            print(f"   {lab:6} bases={n:>9}  at {len(diag)} diagnostic sites: " + ", ".join(f"{k}={v}" for k, v in score.most_common()) + f" (of {tot})")

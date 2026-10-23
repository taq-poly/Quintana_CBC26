#!/bin/bash
# Screen T. canis reference contigs for contaminants before using it to subtract worm reads.
# A contaminant contig left in the reference would silently absorb contaminant reads as "worm".
# Flags: (a) contigs aligning to dog/human/PhiX/E. coli over >=50% of length (minimap2 asm10),
#        (b) contigs Kraken2 assigns to a non-eukaryote / human with >=30% of k-mers.
# Output: refs/tcanis.clean.fa, refs/tcanis.flagged.tsv; then HISAT2 index.
set -euo pipefail
export PATH=/workspace/toxocara_screen/env/bin:$PATH
B=/workspace/toxocara_screen; R=$B/refs; W=$R/refqc; mkdir -p $W $R/idx
cd $W
seqkit seq -m 1 $R/tcanis.fna.gz > tcanis.fa
seqkit fx2tab -n -i -l tcanis.fa > contig_len.tsv

for g in dog human phix ecoli_k12; do
  [ -s $g.paf ] || minimap2 -x asm10 -t 32 --secondary=no $R/$g.fna.gz tcanis.fa > $g.paf 2> /dev/null
done

# C. perfringens added after Kraken2 flagged Clostridium contigs (asm20: allows strain divergence)
[ -s cperf.paf ] || minimap2 -x asm20 -c -t 32 $R/cperf.fna.gz tcanis.fa > cperf.paf 2> /dev/null

[ -s kraken_contigs.out ] || kraken2 --db $B/db/k2_pluspf16 --threads 32 --confidence 0.1 \
  --report kraken_contigs.report --output kraken_contigs.out tcanis.fa 2> /dev/null

python3 - <<'PY'
import collections
L={l.split()[0]:int(l.split()[1]) for l in open("contig_len.tsv")}
cov=collections.defaultdict(lambda: collections.defaultdict(list))
for g in ["dog","human","phix","ecoli_k12","cperf"]:
    for l in open(f"{g}.paf"):
        f=l.split("\t")
        if int(f[11])>=20: cov[f[0]][g].append((int(f[2]),int(f[3])))
def frac(iv,n):
    iv.sort(); tot=0; s=e=-1
    for a,b in iv:
        if a>e: tot+=e-s if e>s else 0; s,e=a,b
        else: e=max(e,b)
    tot+=e-s if e>s else 0
    return tot/n
out=open("../tcanis.flagged.tsv","w"); out.write("contig\tlength\treason\n")
flag=set()
for c,d in cov.items():
    for g,iv in d.items():
        fr=frac(iv,L[c])
        if fr>=0.5: flag.add(c); out.write(f"{c}\t{L[c]}\tminimap2_{g}_{fr:.2f}\n")
# kraken: taxid col 3, kmer string col 5 'taxid:count ...'.
# The 16 GB capped DB keeps only a subset of k-mers, so genuinely bacterial contigs hit at ~10%;
# flag any contig assigned under Bacteria/Archaea/Viruses (or human) with >=5% of k-mers classified.
parent={}
for l in open("/workspace/toxocara_screen/db/k2_pluspf16/nodes.dmp"):
    f=l.split("\t|\t"); parent[f[0]]=f[1]
def lineage(t):
    s=set()
    while t in parent and t not in s: s.add(t); t=parent[t] if parent[t]!=t else None
    return s
for l in open("kraken_contigs.out"):
    f=l.rstrip("\n").split("\t")
    if f[0]!="C": continue
    hit=tot=0
    for t in f[4].split():
        if t=="|:|": continue
        tx,n=t.rsplit(":",1); n=int(n); tot+=n
        if tx not in ("0","A"): hit+=n
    lin=lineage(f[2])
    if tot and hit/tot>=0.05 and lin & {"2","2157","10239","9606"}:
        flag.add(f[1]); out.write(f"{f[1]}\t{L.get(f[1],0)}\tkraken_taxid{f[2]}_{hit/tot:.2f}\n")
out.close()
open("flag_ids.txt","w").write("\n".join(sorted(flag))+"\n")
print(len(flag),"contigs flagged,",sum(L[c] for c in flag if c in L),"bp of",sum(L.values()))
PY
seqkit grep -v -f flag_ids.txt tcanis.fa > $R/tcanis.clean.fa
[ -s $R/idx/tcanis_hisat2.1.ht2 ] || hisat2-build -p 32 $R/tcanis.clean.fa $R/idx/tcanis_hisat2 > $W/hisat2_build.log 2>&1
echo REFQC DONE

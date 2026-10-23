#!/usr/bin/env python3
"""Assign each read pair to the best of {worm, microbes, dog, human}.

Score per mate = aligned bases (M/=/X) - 2*NM, so scores from sr and splice:sr presets are
comparable; pair score = sum over mates (0 for an unmapped mate). A pair goes to the best target
only if it beats the runner-up by >= MARGIN; otherwise it is 'ambiguous'. MAPQ is not used across
targets (it only reflects uniqueness within one target).
Writes assigned_<target>.names and compete_summary.tsv.
"""
import subprocess, sys, re, collections
S = sys.argv[1]
T = ["worm", "microbes", "dog", "human"]
MARGIN = 10
cig = re.compile(r"(\d+)([MIDNSHP=X])")
score = collections.defaultdict(lambda: [0] * len(T))
for i, t in enumerate(T):
    p = subprocess.Popen(["samtools", "view", f"{t}.bam"], stdout=subprocess.PIPE, text=True)
    for l in p.stdout:
        f = l.split("\t", 12)
        aln = sum(int(n) for n, op in cig.findall(f[5]) if op in "M=X")
        nm = next((int(x[5:]) for x in f[11].split("\t") if x.startswith("NM:i:")), 0)
        score[f[0]][i] += max(0, aln - 2 * nm)
    p.wait()
names = {t: [] for t in T}
cnt = collections.Counter()
for q, s in score.items():
    order = sorted(range(len(T)), key=lambda i: -s[i])
    best, second = order[0], order[1]
    if s[best] < 40:
        cnt["unassigned_low_score"] += 1
    elif s[best] - s[second] < MARGIN:
        cnt["ambiguous_" + "/".join(sorted([T[best], T[second]]))] += 1
    else:
        names[T[best]].append(q); cnt[T[best]] += 1
for t in T:
    open(f"assigned_{t}.names", "w").write("\n".join(names[t]) + "\n")
with open("compete_summary.tsv", "w") as o:
    o.write("category\tread_pairs\n")
    for k, v in cnt.most_common(): o.write(f"{k}\t{v}\n")
print(S, dict(cnt.most_common()))

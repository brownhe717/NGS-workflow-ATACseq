#summarize hybrid mapping

import pysam
import sys
from collections import defaultdict

bam_path = sys.argv[1]
out_path = sys.argv[2]

parents = {
    "Dmel": "Dmel_",
    "Dsim": "Dsim_"
}

read_parents = defaultdict(set)
unmapped_reads = set()

bam = pysam.AlignmentFile(bam_path, "rb")

for read in bam.fetch(until_eof=True):
    qname = read.query_name

    if read.is_unmapped:
        unmapped_reads.add(qname)
        continue

    rname = bam.get_reference_name(read.reference_id)

    for parent, prefix in parents.items():
        if rname.startswith(prefix):
            read_parents[qname].add(parent)

bam.close()

# Classification
counts = {
    "Dmel_only": 0,
    "Dsim_only": 0,
    "Both": 0,
    "Unmapped": 0
}

for read, parents_seen in read_parents.items():
    if parents_seen == {"Dmel"}:
        counts["Dmel_only"] += 1
    elif parents_seen == {"Dsim"}:
        counts["Dsim_only"] += 1
    elif parents_seen == {"Dmel", "Dsim"}:
        counts["Both"] += 1

# Reads that never appeared as mapped
counts["Unmapped"] = len(unmapped_reads - set(read_parents.keys()))

total = sum(counts.values())

# Write output
with open(out_path, "w") as out:
    out.write("category\tcount\tpercent\n")
    for k, v in counts.items():
        pct = (v / total * 100) if total > 0 else 0
        out.write(f"{k}\t{v}\t{pct:.2f}\n")

    out.write(f"\nTotal_reads\t{total}\n")

#!/usr/bin/env python
import sys

data = dict()
with open(sys.argv[1]) as f:
    for line in f:
        if line.startswith("#"):
            continue
        row = line.strip("\n").split("\t")
        strand = row[6]
        gid = None
        for s in row[8].split(";"):
            if s.startswith("gene_id"):
                gid = s[9:-1]
        assert gid is not None
        data[gid] = strand
        
with open(sys.argv[2], "w+") as fw:
    for gid, strand in data.items():
        fw.write("%s,%s\n" % (gid, strand))
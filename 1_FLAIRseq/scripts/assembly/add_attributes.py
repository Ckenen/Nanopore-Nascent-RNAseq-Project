#!/usr/bin/env python
import sys
from collections import defaultdict
import pandas as pd
from pyBioInfo.IO.File import GtfFile

def main():
    infile, outfile = sys.argv[1:]
    
    with GtfFile(infile) as f:
        records = [x for x in f]
        
    gid2records = defaultdict(list)
    for r in records:
        gid = r.attributes["gene_id"]
        gid2records[gid].append(r)
        
    tid2records = defaultdict(list)
    for r in records:
        if r.feature == "gene":
            continue
        tid = r.attributes["transcript_id"]
        tid2records[tid].append(r)
        
    rows = []
    gid2gname = dict()
    gid2gtype = dict()
    for r in records:
        if r.feature == "gene":
            gid = r.attributes["gene_id"]
            gname = r.attributes["gene_name"]
            gtype = r.attributes["gene_type"]
            chrom = r.chrom
            start = r.start
            end = r.end
            strand = r.strand
            row = [gid, gname, gtype, chrom, start, end, strand]
            rows.append(row)
            gid2gname[gid] = gname
            gid2gtype[gid] = gtype
    m = pd.DataFrame(rows, columns=[
        "GeneID", "GeneName", "GeneType", 
        "Chrom", "Start", "End", "Strand"])
    m.to_csv(outfile + ".gene_info.tsv", sep="\t", index=False)
    
    rows = []
    for r in records:
        if r.feature == "transcript":
            gid = r.attributes["gene_id"]
            gname = gid2gname[gid]
            gtype = gid2gtype[gid]
            tid = r.attributes["transcript_id"]
            tname = r.attributes.get("transcript_name")
            if tname is None:
                tname = "%s.novel.%s" % (gname, tid.split(".")[-1])
            ttype = r.attributes.get("transcript_type")
            if ttype is None:
                ttype = "Unknown"
            for r2 in tid2records[tid]:
                r2.attributes["gene_name"] = gname
                r2.attributes["gene_type"] = gtype
                r2.attributes["transcript_name"] = tname
                r2.attributes["transcript_type"] = ttype
            chrom = r.chrom
            start = r.start
            end = r.end
            strand = r.strand
            row = [tid, tname, ttype, gid, gname, gtype, chrom, start, end, strand]
            rows.append(row)
    m2 = pd.DataFrame(rows, columns=[
        "TranscriptID", "TranscriptName", "TranscriptType", 
        "GeneID", "GeneName", "GeneType", 
        "Chrom", "Start", "End", "Strand"])
    m2.to_csv(outfile + ".transcript_info.tsv", sep="\t", index=False)
    
    with open(outfile, "w+") as fw:
        for r in records:
            fw.write(r.format() + "\n")
            
        
if __name__ == "__main__":
    main()
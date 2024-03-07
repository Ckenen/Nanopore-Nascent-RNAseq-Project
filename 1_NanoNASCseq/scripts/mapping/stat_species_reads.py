#!/usr/bin/env python
import sys
import os


def main():
    infile, outfile = sys.argv[1:]
    
    hg, hg_mito, mm, mm_mito = 0, 0, 0, 0
    with open(infile) as f:
        for line in f:
            chrom, count = line.strip("\n").split("\t")
            count = int(count)
            if chrom.startswith("hg_"):
                hg += count
                if chrom.endswith("_chrM"):
                    hg_mito += count
            elif chrom.startswith("mm_"):
                mm += count
                if chrom.endswith("_chrM"):
                    mm_mito += count
    
    total = hg + mm
    if total > 0:
        r1 = hg / total
        r2 = mm / total
        r = max(r1, r2)
    else:
        r = 0
    with open(outfile, "w+") as fw:
        fw.write("File\tHuman\tHuman_Mito\tMouse\tMouse_Mito\tSpecificity\n")
        fw.write("\t".join(map(str, [
            os.path.basename(infile), 
            hg, hg_mito, mm, mm_mito, r])) + "\n")
    
    
if __name__ == "__main__":
    main()
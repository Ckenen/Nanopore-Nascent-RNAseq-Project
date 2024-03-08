#!/usr/bin/env python
import sys
import pysam
import multiprocessing as mp
from pyBioInfo.IO.File import BamFile
from pyBioInfo.Utils import ShiftLoader, BlockTools


def worker(f_bed, f_bam, chrom):
    rows = []
    with pysam.TabixFile(f_bed) as bed, BamFile(f_bam) as bam:
        loader = ShiftLoader(bam.fetch(chrom))
        for line in bed.fetch(chrom):
            row = line.strip("\n").split("\t")[:6]
            start = int(row[1])
            end = int(row[2])
            key = (start, end)
            names = []
            for align in loader.fetch(chrom=chrom, start=start, end=end):
                if key in set(BlockTools.gaps(align.blocks)):
                    names.append(align.name)
            names = set(names) # Fragment count
            row[4] = len(names)
            rows.append(row)
    return rows
            

def main():
    f_bed, f_bam, threads, outfile = sys.argv[1:]
    threads = int(threads)
    
    results = []
    pool = mp.Pool(threads)
    with pysam.TabixFile(f_bed) as f:
        for chrom in f.contigs:
            args = (f_bed, f_bam, chrom)
            results.append(pool.apply_async(worker, args))
    pool.close()
    pool.join()
    
    with open(outfile, "w+") as fw:
        for r in results:
            for row in r.get():
                fw.write("\t".join(map(str, row)) + "\n")
            

if __name__ == "__main__":
    main()
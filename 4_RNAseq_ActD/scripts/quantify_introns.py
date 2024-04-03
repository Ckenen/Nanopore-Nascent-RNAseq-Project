#!/usr/bin/env python
import sys
import multiprocessing as mp
from pyBioInfo.IO.File import BedFileRandom, BamFileRandom
from pyBioInfo.Utils import ShiftLoader, BlockTools


def worker(f_bed, f_bam, chrom):
    transcripts = []
    with BedFileRandom(f_bed) as f:
        for t in f.fetch(chrom):
            transcripts.append(t)

    introns = []
    for t in transcripts:
        for start, end in BlockTools.gaps(t.blocks):
            introns.append((chrom, start, end, t.strand))
    introns = list(sorted(set(introns)))

    rows = []
    with BamFileRandom(f_bam) as f:
        loader = ShiftLoader(f.fetch(chrom))
        for chrom, start, end, strand in introns:
            key = (start, end)
            names = []
            for align in loader.fetch(chrom=chrom, start=start, end=end):
                if key in set(BlockTools.gaps(align.blocks)):
                    names.append(align.name)
            names = set(names) # Fragment count
            row = [chrom, start, end, "Intron", len(names), strand]
            rows.append(row)
    return rows
            

def main():
    f_bed, f_bam, threads, outfile = sys.argv[1:]
    threads = int(threads)
    
    results = []
    pool = mp.Pool(threads)
    with BedFileRandom(f_bed) as f:
        chroms = list(sorted(f.handle.contigs))
    for chrom in chroms:
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
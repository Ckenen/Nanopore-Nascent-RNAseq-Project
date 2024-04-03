#!/usr/bin/env python
from collections import defaultdict
import optparse
import multiprocessing as mp
from pyBioInfo.IO.File import BamFileRandom, FastaFileRandom
from pyBioInfo.Utils import BlockTools


def worker(bamfile, fafile, chrom, direction):
    counter = defaultdict(int)
    with BamFileRandom(bamfile) as f:
        for a in f.fetch(chrom):
            s = a.segment
            introns = list(BlockTools.gaps(a.blocks))
            strand = a.strand
            if s.is_paired and s.is_read2:
                strand = "-" if strand == "+" else "+"
            if direction == "R":
                strand = "-" if strand == "+" else "+"
            for x, y in introns:
                counter[(chrom, x, y, strand)] += 1

    motifs = dict()
    with FastaFileRandom(fafile) as f:
        for k in counter:
            chrom, start, end, strand = k
            seq1 = f.fetch(chrom, start, start + 2)
            seq2 = f.fetch(chrom, end - 2, end)
            if strand == "-":
                s1 = FastaFileRandom.reverse_complement(seq1)
                s2 = FastaFileRandom.reverse_complement(seq2)
                seq1, seq2 = s2, s1
            motif = "%s-%s" % (seq1, seq2)
            motifs[k] = motif

    return counter, motifs

            
def main():
    parser = optparse.OptionParser(usage="%prog [options] <input.bam> <genome.fa> <out.bed>")
    parser.add_option("-t", "--threads", dest="threads", type="int", default=0)
    parser.add_option("-d", "--direction", dest="direction", default="F")
    options, args = parser.parse_args()
    bamfile, fafile, outfile = args
    threads = options.threads
    direction = options.direction

    results = []
    pool = mp.Pool(threads)
    with BamFileRandom(bamfile) as f:
        chroms = list(sorted(f.references.keys()))
    for chrom in chroms:
        results.append(pool.apply_async(worker, (bamfile, fafile, chrom, direction)))
    pool.close()
    pool.join()

    with open(outfile, "w+") as fw:
        for r in results:
            counter, motifs = r.get()
            for k in sorted(counter.keys()):
                chrom, start, end, strand = k
                count = counter[k]
                motif = motifs[k]
                row = [chrom, start, end, motif, count, strand]
                fw.write("\t".join(map(str, row)) + "\n")


if __name__ == "__main__":
    main()
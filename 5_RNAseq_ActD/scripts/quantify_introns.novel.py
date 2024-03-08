#!/usr/bin/env python
import sys
from collections import defaultdict
from Bio.Seq import Seq
import pysam
from pyBioInfo.Utils import SegmentTools


def get_introns(s):
    items = SegmentTools.parse_cigar(s)
    for item in items:
        if item[0] == "N":
            x, y = item[3]
            yield x, y
            
            
def main():
    f_bam, f_fasta, direction, f_out = sys.argv[1:]
    
    counter = defaultdict(int)
    with pysam.AlignmentFile(f_bam) as f:
        for chrom in f.references:
            for s in f.fetch(chrom):
                introns = list(get_introns(s))
                strand = "+" if s.is_forward else "-"
                if s.is_paired and s.is_read2:
                    strand = "-" if strand == "+" else "+"
                if direction == "R":
                    strand = "-" if strand == "+" else "+"
                for x, y in introns:
                    counter[(chrom, x, y, strand)] += 1
    
    motifs = dict()
    with pysam.FastaFile(f_fasta) as f:
        for k in counter:
            chrom, start, end, strand = k
            seq1 = f.fetch(chrom, start, start + 2)
            seq2 = f.fetch(chrom, end - 2, end)
            if strand == "+":
                motif = "%s-%s" % (seq1, seq2)
            else:
                motif = "%s-%s" % (str(Seq(seq2).reverse_complement()), str(Seq(seq1).reverse_complement()))
            motifs[k] = motif
            
    with open(f_out, "w+") as fw:
        for k in sorted(counter.keys()):
            chrom, start, end, strand = k
            row = [chrom, start, end, motifs[k], counter[k], strand]
            fw.write("\t".join(map(str, row)) + "\n")
               
            
if __name__ == "__main__":
    main()
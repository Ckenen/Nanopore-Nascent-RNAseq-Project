#!/usr/bin/env python
import sys
import os
import shutil
from collections import defaultdict
from Bio.Seq import Seq
import pysam
from spoa import poa
import subprocess
import multiprocessing as mp


def reverse_complement(seq):
    return str(Seq(seq).reverse_complement())


def get_1_read_accuracy(segment):
    nm = segment.get_tag("NM")
    length = segment.infer_query_length()
    accuracy = 1 - nm / length
    return nm, length, accuracy


def get_2_read_accuracy(segment1, segment2):
    if segment1.get_tag("de") < segment2.get_tag("de"):
        segment = segment1
    else:
        segment = segment2
    return get_1_read_accuracy(segment)


def make_polished_sequence(umi_name, tmpdir):
    fa1 = os.path.join(tmpdir, "%s.raw.fasta" % umi_name)
    fa2 = os.path.join(tmpdir, "%s.consensus.fasta" % umi_name)
    fa3 = os.path.join(tmpdir, "%s.polished.fasta" % umi_name)
    sam = os.path.join(tmpdir, "%s.aligned.sam" % umi_name)
    seqs = []
    with pysam.FastxFile(fa1) as f:
        for record in f:
            seqs.append(record.sequence)

    consensus, msa = poa(seqs)
    with open(fa2, "w+") as fw:
        fw.write(">%s\n" % umi_name)
        fw.write("%s\n" % consensus)

    cmd = "minimap2 -ax map-ont -o %s %s %s" % (sam, fa2, fa1)
    subprocess.check_call(cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    cmd = "racon %s %s %s > %s" % (fa1, sam, fa2, fa3)
    subprocess.check_call(cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    os.remove(fa1)
    os.remove(fa2)
    os.remove(sam)

    
def main():
    f_bam, f_mmi, f_bed, threads, outdir = sys.argv[1:]
    threads = int(threads)

    if not os.path.exists(outdir):
        os.mkdir(outdir)

    tmpdir = os.path.join(outdir, "tmp")
    if not os.path.exists(tmpdir):
        os.mkdir(tmpdir)
        
    rows = []
    umi_sizes = dict()
    with pysam.AlignmentFile(f_bam) as f:
        for chrom in f.references:
            segments = defaultdict(list)
            for s in f.fetch(chrom):
                segments[s.get_tag("CN")].append(s)
            for umi_name, ss in segments.items():
                umi_sizes[umi_name] = len(ss)
                if len(ss) == 1:
                    nm, length, accuracy = get_1_read_accuracy(ss[0])
                    rows.append([umi_name, len(ss), nm, length, accuracy])
                elif len(ss) == 2:
                    nm, length, accuracy = get_2_read_accuracy(ss[0], ss[1])
                    rows.append([umi_name, len(ss), nm, length, accuracy])
                else:
                    seqs = []
                    for s in ss:
                        seq = s.query_sequence
                        if s.is_reverse:
                            seq = reverse_complement(seq)
                        seqs.append(seq)

                    with open(os.path.join(tmpdir, "%s.raw.fasta" % umi_name), "w+") as fw:
                        for s, seq in zip(ss, seqs):
                            fw.write(">%s\n" % s.query_name)
                            fw.write("%s\n" % seq)
            # break

    pool = mp.Pool(threads)
    for umi_name, umi_size in umi_sizes.items():
        if umi_size >= 3:
            pool.apply_async(make_polished_sequence, (umi_name, tmpdir))
    pool.close()
    pool.join()

    polished_fa = os.path.join(outdir, "polished.fasta")
    with open(polished_fa, "w+") as fw:
        for umi_name, umi_size in umi_sizes.items():
            if umi_size >= 3:
                fa3 = os.path.join(tmpdir, "%s.polished.fasta" % umi_name)
                with open(fa3) as f:
                    for line in f:
                        fw.write(line)

    polished_sam = os.path.join(outdir, "aligned.sam")
    cmd = "minimap2 -ax splice -u f -Y --MD -o %s --junc-bed %s -t %d %s %s" % (polished_sam, f_bed, threads, f_mmi, polished_fa)
    subprocess.check_call(cmd, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

    with pysam.AlignmentFile(polished_sam) as f:
        for s in f:
            if s.is_unmapped or s.is_supplementary or s.is_secondary or s.mapping_quality < 60:
                continue
            umi_name = s.query_name
            umi_size = umi_sizes[umi_name]
            nm, length, accuracy = get_1_read_accuracy(s)
            rows.append([umi_name, umi_size, nm, length, accuracy])

    f_tsv = os.path.join(outdir, "summary.tsv")
    with open(f_tsv, "w+") as fw:
        fw.write("UMI\tSize\tNM\tLength\tAccuracy\n")
        for row in rows:
            fw.write("\t".join(map(str, row)) + "\n")

    if os.path.exists(tmpdir):
        shutil.rmtree(tmpdir)


if __name__ == "__main__":
    main()
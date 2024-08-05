#!/usr/bin/env python
import sys
from collections import defaultdict, Counter
from functools import cmp_to_key
import edlib
import pysam
from pyBioInfo.IO.File import Alignment
from pyBioInfo.Utils import BlockTools


# Only available for NanoNASC-seq 

MAX_THREE_CLIP = 20 # 5
MAX_POLYA_DIFF = 20 # 10
MAX_CAP_DIFF = 20
MAX_LENGTH_DIFF = 20
MAX_SPLICING_SITE_DIFF = 20
MAX_UMI_EDIT_DISTANCE = 2
DOWNSAMPLE = True
DOWNSAMPLE_SIZE = 10

def get_clip_count(segment):
    c1 = 0
    c2 = 0
    cigars = segment.cigartuples
    flag, count = cigars[0]
    if flag == pysam.CSOFT_CLIP or flag == pysam.CHARD_CLIP:
        c1 = count
    flag, count = cigars[-1]
    if flag == pysam.CSOFT_CLIP or flag == pysam.CHARD_CLIP:
        c2 = count
    if segment.is_reverse:
        c1, c2 = c2, c1
    return c1, c2


def load_alignments(f):
    for segment in f:
        c1, c2 = get_clip_count(segment)
        if c2 > MAX_THREE_CLIP:
            continue
        obj = Alignment(segment)
        obj.polya = obj.end if obj.strand == "+" else obj.start
        obj.cap = obj.start if obj.strand == "+" else obj.end
        obj.umi = segment.get_tag("UM")
        obj.introns = BlockTools.gaps(obj.blocks)
        yield obj


def load_bundles(alignments):
    array = None
    end = None
    for a in alignments:
        if array is None:
            array = [a]
            end = a.end
        else:
            if a.chrom == array[0].chrom and a.start < end:
                array.append(a)
                end = max(end, a.end)
            else:
                yield array
                array = [a]
                end = a.end
    if array is not None:
        yield array
        

def cluster_by_strandness(alignments):
    clusters = []
    tmp1 = []
    tmp2 = []
    for a in alignments:
        if a.strand == "+":
            tmp1.append(a)
        else:
            tmp2.append(a)
    if len(tmp1) > 0:
        clusters.append(tmp1)
    if len(tmp2) > 0:
        clusters.append(tmp2)
    return clusters   


def cluster_by_polya(clusters):
    new_clusters = []
    for alignments in clusters:
        alignments = list(sorted(alignments, key=lambda a: [a.polya, a.name]))
        tmp = None
        for a in alignments:
            if tmp is None:
                tmp = [a]
            elif a.polya - tmp[-1].polya <= MAX_POLYA_DIFF:
                tmp.append(a)
            else:
                new_clusters.append(tmp)
                tmp = [a]
        if tmp is not None:
            new_clusters.append(tmp)
    return new_clusters
                
            
def cluster_by_cap(clusters):
    new_clusters = []
    for alignments in clusters:
        alignments = list(sorted(alignments, key=lambda a: [a.cap, a.name]))
        tmp = None
        for a in alignments:
            if tmp is None:
                tmp = [a]
            elif a.cap - tmp[-1].cap <= MAX_CAP_DIFF:
                tmp.append(a)
            else:
                new_clusters.append(tmp)
                tmp = [a]
        if tmp is not None:
            new_clusters.append(tmp)
    return new_clusters


def cluster_by_length(clusters):
    new_clusters = []
    for alignments in clusters:
        alignments = list(sorted(alignments, key=lambda a: [len(a), a.name]))
        tmp = None
        for a in alignments:
            if tmp is None:
                tmp = [a]
            elif len(a) - len(tmp[-1]) <= MAX_LENGTH_DIFF:
                tmp.append(a)
            else:
                new_clusters.append(tmp)
                tmp = [a]
        if tmp is not None:
            new_clusters.append(tmp)
    return new_clusters


def cluster_by_introns(clusters):
    new_clusters = []
    for alignments0 in clusters:
        data = defaultdict(list)
        for a in alignments0:
            data[len(a.blocks)].append(a)
        for exon_count, alignments in data.items():
            if exon_count == 1: # mono-exon
                new_clusters.append(alignments)
            else: # multi-exons
                graph = defaultdict(list)
                for i in range(len(alignments) - 1):
                    for j in range(i + 1, len(alignments)):
                        a1, a2 = alignments[i], alignments[j]
                        identical = True
                        for intron1, intron2 in zip(a1.introns, a2.introns):
                            if abs(intron1[0] - intron2[0]) > MAX_SPLICING_SITE_DIFF \
                                or abs(intron1[1] - intron2[1]) > MAX_SPLICING_SITE_DIFF:
                                identical = False
                                break
                        if identical:
                            graph[id(a1)].append(a2)
                            graph[id(a2)].append(a1)      
                id_to_alignment = {id(a): a for a in alignments}         
                set_masked = set()
                for a in alignments:
                    a_id = id(a)
                    if a_id in set_masked:
                        continue
                    current_cluster = [] # this cluster
                    set_search = set() # umis will be search
                    set_search.add(a_id)
                    while len(set_search) > 0:
                        a1_id = set_search.pop()
                        a1 = id_to_alignment[a1_id]
                        current_cluster.append(a1)
                        set_masked.add(a1_id)
                        for a2 in graph[a1_id]:
                            a2_id = id(a2)
                            if a2_id in set_masked:
                                continue
                            set_search.add(a2_id)
                    new_clusters.append(current_cluster)
    
    return new_clusters         


def edlib_align(seq1, seq2):
    ed1 = edlib.align(seq1, seq2, mode="HW")["editDistance"]
    ed2 = edlib.align(seq2, seq1, mode="HW")["editDistance"]
    ed = min(ed1, ed2)
    return ed
    

def cluster_by_umi_counts(counter):
    items = list(sorted(counter.items(), key=lambda item: (item[1], item[0]), reverse=True))
    umis = [item[0] for item in items]
    graph = defaultdict(list)
    for i in range(len(umis) - 1):
        for j in range(i + 1, len(umis)):
            umi1 = umis[i]
            umi2 = umis[j]
            ed = edlib_align(umi1, umi2)
            if ed <= MAX_UMI_EDIT_DISTANCE:
                count1, count2 = counter[umi1], counter[umi2]
                if True:
                    # umi-tools directional network rule
                    if counter[umi1] >= counter[umi2] * 2 - 1:
                        graph[umi1].append(umi2)
                    if counter[umi2] >= counter[umi1] * 2 - 1:
                        graph[umi2].append(umi1)
                elif False:
                    # custom directional network rule
                    if count1 >= count2:
                        graph[umi1].append(umi2)
                    if count2 >= count1:
                        graph[umi2].append(umi1)
                else:
                    # bidirectional network
                    graph[umi1].append(umi2)
                    graph[umi2].append(umi1)

    set1 = set() # searched umis
    umi_clusters = []
    for umi in umis:
        if umi in set1:
            continue
        set2 = set() # this umi cluster
        anchor_umi = umi
        # anchor_umi_count = counter[umi]
        set3 = set() # umis will be search
        set3.add(umi)
        while len(set3) > 0:
            umi1 = set3.pop()
            set1.add(umi1)
            set2.add(umi1)
            for umi2 in graph[umi1]:
                if umi2 in set1:
                    continue
                if edlib_align(anchor_umi, umi2) > 4:
                    continue
                set3.add(umi2) 
        umi_clusters.append(set2)
    return umi_clusters
    

def cluster_by_umi(clusters):
    new_clusters = []
    for alignments in clusters:
        data = defaultdict(list) # UMI: [alignments]
        for a in alignments:
            data[a.umi].append(a)
            
        counter = {k: len(v) for k, v in data.items()}
        for cluster in cluster_by_umi_counts(counter): 
            objs = []
            for umi in cluster:
                for obj in data[umi]:
                    objs.append(obj)
            new_clusters.append(objs)
    return new_clusters


def get_representative_umi(alignments):
    counter = Counter([a.umi for a in alignments])
    items = list(sorted(counter.items(), key=lambda item: [item[1], item[0]], reverse=True))
    max_count = items[0][1]
    representative_umis = []
    for umi, umi_count in items:
        if umi_count == max_count:
            representative_umis.append(umi)
        else:
            break
    representative_umis.sort()
    representative_umi = representative_umis[0]
    return representative_umi


def merge_identical_representative_umi_clusters(clusters):
    new_clusters = []
    data = defaultdict(list)
    for cluster in clusters:
        representative_umi = get_representative_umi(cluster)
        data[representative_umi].extend(cluster)
    for k, v in data.items():
        new_clusters.append(v)
    return new_clusters


def compare_best_alignment(a1, a2):
    de1 = a1.segment.get_tag("de")
    de2 = a2.segment.get_tag("de")
    if de1 < de2:
        return -1
    elif de1 == de2:
        length1 = len(a1)
        length2 = len(a2)
        if length1 < length2:
            return 1
        elif length1 == length2:
            if a1.name < a2.name:
                return -1
            elif a1.name == a2.name:
                assert False
                return 0
            else:
                return 1
        else:
            return -1
    else:
        return 1
    

def downsample_alignments(alignments):
    alignments1 = [] # keep
    alignments2 = [] # ignore
    data = defaultdict(list)
    for a in alignments:
        data[a.umi].append(a)
    counter = Counter([a.umi for a in alignments])
    items = list(sorted(counter.items(), key=lambda item: (item[1], item[0]), reverse=True))
    # most abundant UMI
    # lowest de tag value
    for umi, umi_count in items:
        for a in sorted(data[umi], key=cmp_to_key(compare_best_alignment)):
            if len(alignments1) < DOWNSAMPLE_SIZE:
                alignments1.append(a)
            else:
                alignments2.append(a)
    return alignments1, alignments2


def modified_duplicate_flag(alignments):
    # select the lowest de value alignment of the most abundant UMI as the uniq alignment
    data = defaultdict(list)
    for a in alignments:
        data[a.umi].append(a)
    max_count = max([len(v) for v in data.values()])
    array = []
    for k, v in data.items():
        if len(v) == max_count:
            array.extend(v)
    for a in sorted(array, key=cmp_to_key(compare_best_alignment)):
        a.segment.flag = a.segment.flag & 0b101111111111
        break


def main():
    infile, outfile, statfile = sys.argv[1:]
    cell = infile.split("/")[-1][:-4]
        
    print("Infile: %s" % infile)
    print("Outfile: %s" % outfile)
    print("Statfile: %s" % statfile)
    print("Cell: %s" % cell)

    with pysam.AlignmentFile(infile) as h_in_bam, \
        pysam.AlignmentFile(outfile, "wb", h_in_bam) as h_out_bam, \
        open(statfile, "w+") as h_stat:
        h_stat.write("Chrom\tStart\tEnd\tStrand\tUMI\tDownSampleSize\tDownSampleUMIs\tAllSize\tAllUMIs\n")

        umi_id = 0
        for bundle_alignments in load_bundles(load_alignments(h_in_bam)):
            for stranded_alignments in cluster_by_strandness(bundle_alignments):
                for sub_alignments in load_bundles(stranded_alignments):
                    clusters = [sub_alignments]

                    cluster_number = None
                    repeat = 0
                    while True:
                        clusters = cluster_by_polya(clusters)
                        clusters = cluster_by_cap(clusters)
                        clusters = cluster_by_length(clusters)
                        clusters = cluster_by_introns(clusters)
                        clusters = cluster_by_umi(clusters)
                        if repeat == 0:
                            cluster_number = len(clusters)
                            repeat += 1
                            continue
                        if len(clusters) == cluster_number:
                            break
                        else:
                            assert len(clusters) > cluster_number
                            cluster_number = len(clusters)
                            repeat += 1
                            continue

                    clusters = merge_identical_representative_umi_clusters(clusters)
                
                    for alignments in clusters:

                        alignments1, alignments2 = downsample_alignments(alignments)

                        chrom = alignments1[0].chrom
                        start = min([a.start for a in alignments1])
                        end = max([a.end for a in alignments1])
                        strand = alignments1[0].strand
                        cn = "UMI.%d" % umi_id
                        cs = len(alignments)

                        for a in alignments:
                            a.segment.flag = a.segment.flag | 0b010000000000
                            a.segment.set_tag("CN", cn)
                            a.segment.set_tag("CS", cs)
                        modified_duplicate_flag(alignments1) # select the uniq alignment
                        for a in alignments1:
                            a.keep = True
                        for a in alignments2:
                            a.keep = False
                        cs1 = len(alignments1)

                        umi_counter = Counter([a.umi for a in alignments])
                        items = list(sorted(umi_counter.items(), key=lambda item: (item[1], item[0]), reverse=True))
                        s = ";".join(["%s:%d" % (umi, count) for umi, count in items])

                        umi_counter1 = Counter([a.umi for a in alignments1])
                        items1 = list(sorted(umi_counter1.items(), key=lambda item: (item[1], item[0]), reverse=True))
                        s1 = ";".join(["%s:%d" % (umi, count) for umi, count in items1])

                        line = "\t".join(map(str, [chrom, start, end, strand, cn, cs1, s1, cs, s]))
                        h_stat.write(line + "\n")
                        umi_id += 1
                
            for a in bundle_alignments:
                if a.keep:
                    h_out_bam.write(a.segment)

    print("Finished!")


if __name__ == '__main__':
    main()

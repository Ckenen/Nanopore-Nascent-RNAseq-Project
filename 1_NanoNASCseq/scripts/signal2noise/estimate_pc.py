#!/usr/bin/env python
import sys
import pickle
import optparse
from collections import defaultdict
import numpy as np
import pandas as pd
from scipy.stats import binom
import pysam


def createMkn(Akn, p_e): #Left out from Akn
    M = np.zeros(Akn.shape)
    for n in range(Akn.shape[1]):
        for k in range(Akn.shape[0]):
            Ekn = np.sum(Akn[(k+1):, n]) * binom.pmf(k, n, p_e)
            if Ekn > 0.01 * Akn[k, n]:
                M[k, n] = 1
    return M


def EstepAkn(Akn, Mkn, p_c): # Alters Akn in place - modify initial step
    for k in range(Akn.shape[0]):
        for n in range(Akn.shape[1]):
            if Mkn[k, n] == 1:
                num = 0
                denom = 0
                for kp in range(Mkn.shape[0]):
                    if Mkn[kp, n] == 1:
                        num = num + binom.pmf(k, n, p_c) * Akn[kp, n]
                        denom = denom + binom.pmf(kp, n, p_c)
                Akn[k, n] = num / denom
    return Akn


def MstepP_c(Akn):
    num = 0
    denom = 0
    for k in range(Akn.shape[0]):
        for n in range(Akn.shape[1]):
            num = num + k * Akn[k, n]
            denom = denom + n * Akn[k, n]
    p_c = num / denom
    return p_c


def run_on_long_read_mode(model_tsv, ratio_tsv, event_tsv):
    model = pd.read_csv(model_tsv, sep="\t")
    assert len(model) > 0
    ratios = pd.read_csv(ratio_tsv, sep="\t", index_col=0)
    p_e = 0
    for t, k, r, w in model.values:
        p_e += (ratios["Ratio"][t] * k * w)

    dat = pd.read_csv(event_tsv, sep="\t", header=0, index_col=0)
    dat = dat[dat["Size"] >= 2]
    if len(dat) > 0:
        max_n = max(dat["T"])
        max_k = max(dat["T-C"])
        Akn = np.zeros((max_k + 1, max_n + 1))
        for n, k in dat[["T", "T-C"]].values:
            Akn[k, n] += 1
        Akn0 = Akn
        l = 0
        r = 1
        p_c0 = (l + r) / 2
        p_c = p_c0
        Mkn = createMkn(Akn0, p_e)
        Akn = Akn0
        while r - l >= 10e-8:
            Akn = EstepAkn(Akn, Mkn, p_c)
            p_c_old = p_c
            p_c = MstepP_c(Akn)
            if p_c < p_c_old:
                r = p_c
            else:
                l = p_c
    else:
        p_c = np.nan

    return p_e, p_c


def run_on_short_read_mode(model_tsv, ratio_tsv, f_bam):
    model = pd.read_csv(model_tsv, sep="\t")
    assert len(model) > 0
    dat = pd.read_csv(ratio_tsv, sep="\t", index_col=0)
    p_e = 0
    for t, k, r, w in model.values:
        p_e += (dat["Ratio.NoSNP"][t] * k * w)

    mapper = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}
    data = defaultdict(int)
    with pysam.AlignmentFile(f_bam) as f:
        for s in f:
            strand = s.get_tag("ST")
            base_counter = defaultdict(int)
            for x in s.get_tag("RC").split(";"):
                if x.strip() != "":
                    k, v = x.split(",")
                    if strand == "-":
                        k = mapper[k]
                    base_counter[k] = int(v)
            event_counter = defaultdict(int)
            for x in s.get_tag("CE").split(";"):
                if x.strip() != "":
                    pos, ref, alt, qua, dis = x.split(",")
                    if strand == "-":
                        ref, alt = mapper[ref], mapper[alt]
                    event_counter["%s-%s" % (ref, alt)] += 1
            n = base_counter["T"]
            k = event_counter["T-C"]
            data[(n, k)] += 1
            
    max_n = max([n for n, k in data.keys()])
    max_k = max([k for n, k in data.keys()])
    
    Akn = np.zeros((max_k + 1, max_n + 1))
    for (n, k), c in data.items():
        Akn[k, n] = c
    Akn0 = Akn

    l = 0
    r = 1
    p_c0 = (l + r) / 2

    p_c = p_c0
    Mkn = createMkn(Akn0, p_e)
    Akn = Akn0

    while r - l >= 10e-8:
        # print(r-l)
        Akn = EstepAkn(Akn, Mkn, p_c)
        p_c_old = p_c
        p_c = MstepP_c(Akn)
        if p_c < p_c_old:
            r = p_c
        else:
            l = p_c

    return p_e, p_c


def main():
    parser = optparse.OptionParser("%prog [options] <model.tsv> <ratio.tsv>")
    parser.add_option("-m", "--mode", dest="mode", help="Mode: long or short.")
    group = optparse.OptionGroup(parser, "Long-reads:")
    group.add_option("-e", "--event-table", dest="events")
    parser.add_option_group(group)
    group = optparse.OptionGroup(parser, "Short-reads:")
    group.add_option("-b", "--bam", dest="bam")
    parser.add_option_group(group)
    options, args = parser.parse_args()
    model_tsv, ratio_tsv = args

    mode = options.mode
    if mode == "long":
        p_e, p_c = run_on_long_read_mode(model_tsv, ratio_tsv, options.bam)
    elif mode == "short":
        p_e, p_c = run_on_short_read_mode(model_tsv, ratio_tsv, options.bam)
    else:
        raise RuntimeError()
    print("Pe\tPc\tPc/Pe")
    print(p_e, p_c, p_c/p_e, sep="\t")


# def main():
#     ratio_tsv, model_tsv, event_tsv = sys.argv[1:]
    
#     ratios = pd.read_csv(ratio_tsv, sep="\t", index_col=0)
#     model = pd.read_csv(model_tsv, sep="\t")
#     p_e = 0
#     for t, k, r, w in model.values:
#         p_e += (ratios["Ratio"][t] * k * w)

#     dat = pd.read_csv(event_tsv, sep="\t", header=0, index_col=0)
#     dat = dat[dat["Size"] >= 2]
#     if len(dat) > 0:
#         max_n = max(dat["T"])
#         max_k = max(dat["T-C"])
        
#         Akn = np.zeros((max_k + 1, max_n + 1))
#         for n, k in dat[["T", "T-C"]].values:
#             Akn[k, n] += 1
#         Akn0 = Akn

#         l = 0
#         r = 1
#         p_c0 = (l + r) / 2

#         p_c = p_c0
#         Mkn = createMkn(Akn0, p_e)
#         Akn = Akn0

#         while r - l >= 10e-8:
#             # print(r-l)
#             Akn = EstepAkn(Akn, Mkn, p_c)
#             p_c_old = p_c
#             p_c = MstepP_c(Akn)
#             if p_c < p_c_old:
#                 r = p_c
#             else:
#                 l = p_c
#     else:
#         p_c = np.nan

#     print("Pe\tPc\tPc/Pe")
#     print(p_e, p_c, p_c/p_e, sep="\t")


if __name__ == "__main__":
    main()
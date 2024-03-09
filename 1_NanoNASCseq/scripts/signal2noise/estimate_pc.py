#!/usr/bin/env python
import sys
import pickle
import numpy as np
import pandas as pd
from scipy.stats import binom


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


def main():
    ratio_tsv, model_tsv, event_tsv = sys.argv[1:]
    
    ratios = pd.read_csv(ratio_tsv, sep="\t", index_col=0)
    model = pd.read_csv(model_tsv, sep="\t")
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
            # print(r-l)
            Akn = EstepAkn(Akn, Mkn, p_c)
            p_c_old = p_c
            p_c = MstepP_c(Akn)
            if p_c < p_c_old:
                r = p_c
            else:
                l = p_c
    else:
        p_c = np.nan

    print("Pe\tPc\tPc/Pe")
    print(p_e, p_c, p_c/p_e, sep="\t")


if __name__ == "__main__":
    main()
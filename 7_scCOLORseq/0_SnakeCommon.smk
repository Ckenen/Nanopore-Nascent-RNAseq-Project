#!/usr/bin/env runsnakemake
configfile: "config.yaml"
srr_list = config["samples"]
import glob
import pandas as pd
dat = pd.read_csv("data/SraRunTable.scCOLORseq.csv")
srr_batch_list = []
srr_to_batch = dict()
for srr in srr_list:
    srr_to_batch[srr] = list()
    d = "results/prepare/split/%s" % srr
    if os.path.exists(d):
        for path in sorted(glob.glob(d + "/*.fastq.gz")):
            batch = path.split("/")[-1][:-9]
            srr_batch_list.append("%s/%s" % (srr, batch))
            srr_to_batch[srr].append(batch)
# print("SRR batch:", len(srr_batch_list))
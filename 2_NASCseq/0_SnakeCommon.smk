#!/usr/bin/env runsnakemake
configfile: "config.yaml"
import pandas as pd
RUNS = config["RUNS"]
THREADS = config["THREADS"]
DAT = pd.read_excel(config["TABLE"])
DAT.index = DAT["Cell"]
DAT = DAT[DAT["Species"] == config["SPECIES"]]
DAT["RunCell"] = ["%s/%s" % (run, cell) for run, cell in DAT[["Run", "Cell"]].values]
DAT1 = DAT[DAT["Run"].isin(RUNS)]
RUN_CELLS = DAT1["RunCell"].values
RUN_CELLS_PE = DAT1[DAT1["Layout"] == "PE"]["RunCell"].values
RUN_CELLS_SE = DAT1[DAT1["Layout"] == "SE"]["RunCell"].values

def get_layout(cell):
    return DAT[DAT["Cell"] == cell]["Layout"].values[0]


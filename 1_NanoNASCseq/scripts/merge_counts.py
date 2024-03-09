#!/usr/bin/env python
import sys
from collections import defaultdict
import pandas as pd

            
def main():
    filelist, outfile = sys.argv[1:]
    
    data = dict()
    header = None
    for path in open(filelist):
        path = path.strip()
        m = pd.read_csv(path, sep="\t", header=0)
        header = m.columns
        for gid, total, nascent, total_alleles, nascent_alleles in m.values:
            if gid not in data:
                data[gid] = {
                    "Total": 0, 
                    "Nascent": 0, 
                    "Total.Alleles": defaultdict(int), 
                    "Nascent.Alleles": defaultdict(int)}
            data[gid]["Total"] += total
            data[gid]["Nascent"] += nascent
            if isinstance(total_alleles, str):
                for item in total_alleles.split(";"):
                    allele, count = item.split(":")
                    data[gid]["Total.Alleles"][allele] += int(count)
            if isinstance(nascent_alleles, str):
                for item in nascent_alleles.split(";"):
                    allele, count = item.split(":")
                    data[gid]["Nascent.Alleles"][allele] += int(count)
    
    with open(outfile, "w+") as fw:
        fw.write("\t".join(header) + "\n")
        # fw.write("GeneID\tTotal\tNascent\tTotal.Alleles\tNascent.Alleles\n")
        for gid in sorted(data.keys()):
            total = data[gid]["Total"]
            nascent = data[gid]["Nascent"]
            total_alleles = data[gid]["Total.Alleles"]
            nascent_alleles = data[gid]["Nascent.Alleles"]
            s1 = ";".join(["%s:%d" % (k, v) for k, v in sorted(total_alleles.items())])
            s2 = ";".join(["%s:%d" % (k, v) for k, v in sorted(nascent_alleles.items())])
            line = "\t".join(map(str, [gid, total, nascent, s1, s2]))
            fw.write(line + "\n")
        

if __name__ == "__main__":
    main()
import argparse
import re
from collections import defaultdict

parser = argparse.ArgumentParser()
parser.add_argument("-g", "--gene_file")
parser.add_argument("-s", "--snp_file")
parser.add_argument("-o", "--out_file")

arg = parser.parse_args()

gene_dict = defaultdict(bool)
g_f = open(arg.gene_file, "r")
g_l = g_f.readline()
while g_l:
    gene_dict[re.sub("\n", "", g_l)] = True
    g_l = g_f.readline()
g_f.close()

s_written = defaultdict(bool)
s_f = open(arg.snp_file, "r")
s_l = s_f.readline()
s_l = s_f.readline()
while s_l:
    r_l = re.sub("\n", "", s_l).split("\t")
    snps = r_l[1].split(",")
    for s in snps:
        s_written[s] = True
    s_l = s_f.readline()
s_f.close()

o_f = open(arg.out_file, "w+")
for k in s_written.keys():
    o_f.write(k + "\n")

o_f.close()



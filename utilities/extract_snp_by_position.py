import argparse
import re
from collections import defaultdict
import gzip

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--inputfile")
parser.add_argument("-s", "--snpfile")
parser.add_argument("-o", "--outputfile")
args = parser.parse_args()

snp_dict = defaultdict(bool)
s_f = open(args.snpfile, "r")
s_l = s_f.readline()
while s_l:
    snp_dict[":".join(re.sub("\n", "", s_l).split(":")[:2])] = True
    s_l = s_f.readline()
s_f.close()

v_f = gzip.open(args.inputfile, "rb")
header = True
o_f = open(args.outputfile, "w+")
for v_l in v_f:
    if header:
        x = re.search(b'^##', v_l)
        o_f.write(v_l.decode("utf-8"))
        if x:
            continue
        else:
            header = False
            continue
    v_l_d = v_l.decode("utf-8")
    r_v_l = re.sub("\n", "", v_l_d).split("\t")
    snp_key = f'{r_v_l[0]}:{r_v_l[1]}'
    if snp_dict[snp_key]:
        o_f.write(v_l_d)
v_f.close()
o_f.close()

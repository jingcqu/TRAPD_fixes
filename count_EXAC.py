import operator
import re
import sys
import gzip
import bisect
import argparse
from collections import defaultdict


parser = argparse.ArgumentParser()
parser.add_argument("-v", "--vcf")
parser.add_argument("-g", "--gene_file")
parser.add_argument("-p", "--population")
parser.add_argument("-o", "--output")
args = parser.parse_args()
vcffile=gzip.open(args.vcf, "rb")

population = args.population

g_f = open(args.gene_file, "r")
snp_key_check = defaultdict(bool)
snp_info_dict = dict()
gene_lists = defaultdict(list)
g_l = g_f.readline()
g_l = g_f.readline()
while g_l:
    r_g_l = re.sub("\n", "", g_l).split("\t")
    snps = r_g_l[-1].split(",")
    for s in snps:
        snp_key_check[s] = True
        snp_info_dict[s] = [r_g_l[0], 0, 0, 0, 0] # [gene, ac, AC_HET, AC_HOM, an]
        gene_lists[r_g_l[0]].append(s)
    g_l = g_f.readline()
g_f.close()
header = True
for v_l in vcffile:
    if header:
        x = re.search(b'^##', v_l)
        if x:
            continue
        else:
            header = False
            continue
    v_l_d = v_l.decode("utf-8")
    r_v_l = re.sub("\n", "", v_l_d).split("\t")
    variants = r_v_l[4].split(",")
    if len(variants) < 2:
        snp_key = f'{r_v_l[0]}:{r_v_l[1]}:{r_v_l[3]}:{r_v_l[4]}'
        if snp_key_check[snp_key]:
            info_field = r_v_l[7].split(";")
            snp_info_dict[snp_key]
            for info in info_field:
                r_info = info.split("=")
                if population == "ALL":
                    if r_info[0] == "AC":
                        snp_info_dict[snp_key][1] = int(r_info[1])
                    elif r_info[0] == "AN":
                        snp_info_dict[snp_key][4] = int(r_info[1])
                    elif r_info[0] == "AC_Het":
                        snp_info_dict[snp_key][2] = int(r_info[1])
                    elif r_info[0] == "AC_Hom":
                        snp_info_dict[snp_key][3] = int(r_info[1])
                else:
                    if r_info[0] == ("AC_" + population):
                        snp_info_dict[snp_key][1] = int(r_info[1])
                    elif r_info[0] == ("AN_" + population):
                        snp_info_dict[snp_key][4] = int(r_info[1])
                    elif r_info[0] == ("Het_" + population):
                        snp_info_dict[snp_key][2] = int(r_info[1])
                    elif r_info[0] == ("Hom_" + population):
                        snp_info_dict[snp_key][3] = int(r_info[1])
    else:
        for i in range(len(variants)):
            snp_key = f'{r_v_l[0]}:{r_v_l[1]}:{r_v_l[3]}:{variants[i]}'
            if snp_key_check[snp_key]:
                info_field = r_v_l[7].split(";")
                snp_info_dict[snp_key]
                for info in info_field:
                    r_info = info.split("=")
                    if population == "ALL":
                        if r_info[0] == "AC":
                            snp_info_dict[snp_key][1] = int(r_info[1].split(",")[i])
                        elif r_info[0] == "AN":
                            snp_info_dict[snp_key][4] = int(r_info[1])
                        elif r_info[0] == "AC_Het":
                            snp_info_dict[snp_key][2] = int(r_info[1].split(",")[i])
                        elif r_info[0] == "AC_Hom":
                            snp_info_dict[snp_key][3] = int(r_info[1].split(",")[i])
                    else:
                        if r_info[0] == ("AC_" + population):
                            snp_info_dict[snp_key][1] = int(r_info[1].split(",")[i])
                        elif r_info[0] == ("AN_" + population):
                            snp_info_dict[snp_key][4] = int(r_info[1])
                        elif r_info[0] == ("Het_" + population):
                            snp_info_dict[snp_key][2] = int(r_info[1].split(",")[i])
                        elif r_info[0] == ("Hom_" + population):
                            snp_info_dict[snp_key][3] = int(r_info[1].split(",")[i])
vcffile.close()

out_f = open(args.output, "w+")
out_f.write("GENE,AC_Het_sum,AC_Het_avg,AC_Hom_sum,AC_Hom_avg,AC_sum,AC_avg,AN_sum,AN_avg,snp,AC_Het,AC_Hom,AC,AN\n")
for k in gene_lists.keys():
    info = [0, 0, 0, 0, 0, 0, 0, 0] #AC_HET_sum AC_HET_avg AC_HOM_sum AC_HOM_avg AC_sum AC_avg AN_sum AN_avg 
    num_snps = 0
    for s in gene_lists[k]: #following columns are snpname, AC, AC_het, AC_HOM, an
        num_snps += 1
        info[0] += snp_info_dict[s][2]
        info[2] += snp_info_dict[s][3]
        info[4] += snp_info_dict[s][1]
        info[6] += snp_info_dict[s][4]
        info.append(s)
        info.append(snp_info_dict[s][2])
        info.append(snp_info_dict[s][3])
        info.append(snp_info_dict[s][1])
        info.append(snp_info_dict[s][4])
    for i in range(1, 8, 2):
        info[i] = info[i-1] / num_snps
    for i in range(len(info)):
        info[i] = str(info[i])
    out_f.write(f'{k},{",".join(info)}\n')
out_f.close()

        

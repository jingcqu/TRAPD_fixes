import argparse
from collections import defaultdict
import re

"""
THREE filters:

Read depth per sample (DP > 10)

NFE MAF <0.03

Predicted deleterious in at least 1 of a number of gene prediction softwares:
'SIFT_pred="D" || Polyphen2_HDIV_pred="D" || Polyphen2_HDIV_pred="P" || Polyphen2_HVAR_pred="D" || Polyphen2_HVAR_pred="P" || LRT_pred="D" || MutationTaster_pred="A" || MutationTaster_pred="D" || MutationAssessor_pred="H" || MutationAssessor_pred="M" || FATHMM_pred="D" || PROVEAN_pred="D" || MetaSVM_pred="D" || MetaLR_pred="D"
"""


parser = argparse.ArgumentParser()
parser.add_argument("-v", "--vcffile", help="vcf file")
parser.add_argument("-ov", "--outvcffile", help="out vcf file")
parser.add_argument("-os", "--outsnpfile", help="out snp file")
parser.add_argument('-c', "--conditions", help="le <=, ge >=, l <, g >, e. Sample: -c Polyphen2_HVAR_pred:e:D,Polyphen2_HVAR_pred:e:P,ExAC_NFE:l:0.03")
args = parser.parse_args()

vcf_file = args.vcffile
out_vcf_file = args.outvcffile
out_snp_file = args.outsnpfile

bad_lines = "^##"

vcf_file = open(vcf_file, "r")
out_vcf_file = open(out_vcf_file, "w+")
vcf_line = vcf_file.readline()
header = True
while header:
    x = re.search(bad_lines, vcf_line)
    if x:
        out_vcf_file.write(vcf_line)
        vcf_line = vcf_file.readline()
    else:
        out_vcf_file.write(vcf_line)
        vcf_line = vcf_file.readline()
        header = False



"""
snplist
header
GENE chr:pos:Ref:Alt

vcf
GENE: Gene.refGene
chr: index 0 
pos: index 1
ref: index 3
alt index 4

"""
good_rs = defaultdict(list)

while vcf_line:
    rawline = vcf_line.split("\t")
    info = rawline[7]
    info = info.split(";")
    info_values = defaultdict(str)
    for field in info:
        broke = field.split("=")
        if len(broke) > 1:
            info_values[broke[0]] = broke[1]
    #print(info)
    ##test
    #exit(1)
    if info_values["ExAC_NFE"] != ".":
        if float(info_values["ExAC_NFE"]) < 0.03:
            if info_values["SIFT_pred"] == "D" or info_values["Polyphen2_HDIV_pred"] == "D" or info_values["Polyphen2_HDIV_pred"]=="P" or info_values["Polyphen2_HVAR_pred"] == "D" or info_values["Polyphen2_HVAR_pred"] == "P" or info_values["LRT_pred"] == "D" or info_values["MutationTaster_pred"] == "A" or info_values["MutationTaster_pred"] == "D" or info_values["MutationAssessor_pred"] == "H" or info_values["MutationAssessor_pred"] == "M" or info_values["FATHMM_pred"] == "D" or info_values["PROVEAN_pred"] == "D" or info_values["MetaSVM_pred"] == "D" or info_values["MetaLR_pred"] == "D":
                out_vcf_file.write(vcf_line)
                variants = rawline[4].split(",")
                good_rs[info_values["Gene.refGene"]].append(rawline[0] + ":" + rawline[1] + ":" + rawline[3] + ":" + variants[0])
    vcf_line = vcf_file.readline()
vcf_file.close()
out_vcf_file.close()


out_snp_file = open(out_snp_file, "w+")
out_snp_file.write("#GENE\tSNPS\n")

for k in good_rs.keys():
    sep = ","
    out_snp_file.write(f'{k}\t{sep.join(good_rs[k])}\n')
out_snp_file.close()

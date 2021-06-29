
#!/usr/bin/env Rscript
library("argparse")
suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser()
parser$add_argument("--casefile", action="store")
parser$add_argument("--casesize", action="store", type="integer")
parser$add_argument("--controlfile", action="store")
parser$add_argument("--outfile", action="store")

args <- parser$parse_args()

casefile<-args$casefile
controlfile<-args$controlfile
casesize<-args$casesize
outfile<-args$outfile

case.dat<-read.delim(casefile, header=T, stringsAsFactors=F, sep="\t")
names(case.dat)[1]<-"GENE"
control.dat<-read.delim(controlfile, header=T, stringsAsFactors=F, sep=",")
names(control.dat)[1]<-"GENE"

dat<-merge(case.dat, control.dat, by="GENE", all.x=T, all.y=T)
dat[is.na(dat)]<-0

dat$DOM_case<-0
dat$DOM_case1<-0
dat$DOM_control<-0
dat$DOM_control1<-0
dat$REC_case<-0
dat$REC_case1<-0
dat$REC_control<-0
dat$REC_control1<-0
dat$RealREC_case<-0
dat$RealREC_case1<-0
dat$RealREC_control<-0
dat$RealREC_control1<-0
dat$P_DOM<-0
dat$P_REC<-0
dat$P_REC2<-0


for(i in 1:nrow(dat)){
  
  #Dominant model
  case_count<-dat[i,]$CASE_COUNT_HET+dat[i,]$CASE_COUNT_HOM+dat[i,]$CASE_COUNT_CH
  control_count<-dat[i,]$AC_Het_sum+dat[i,]$AC_Hom_sum
  control_size<-dat[i,]$AN_avg/2
  control_size<-round(control_size,0)
  casesize<-casesize  
  
  dat[i,]$DOM_case<-case_count
  dat[i,]$DOM_case1<-(casesize-case_count)
  dat[i,]$DOM_control<-control_count
  dat[i,]$DOM_control1<-(control_size-control_count)
  if(case_count>casesize){
    case_count<-casesize
  }else if(case_count<0){
    case_count<-0
  }
  if (((casesize-case_count) > 0) && ((control_size-control_count) > 0)){
    mat<-cbind(c(case_count, (casesize-case_count)), c(control_count, (control_size-control_count)))
    dat[i,]$P_DOM<-fisher.test(mat, alternative="greater")$p.value
  }
  else{
    dat[i,]$P_DOM<-1
  }
  
  #Recessive model
  case_count_rec<-dat[i,]$CASE_COUNT_HOM
  control_count_rec<-dat[i,]$AC_Hom_sum
  dat[i,]$REC_case<-case_count_rec
  dat[i,]$REC_case1<-(casesize-case_count_rec)
  dat[i,]$REC_control<-control_count_rec
  dat[i,]$REC_control1<-(control_size-control_count_rec)
  
  if(control_count_rec<0){ control_count_rec<-0}
  if(case_count_rec>casesize){case_count_rec<-casesize}
  if (((casesize-case_count_rec) > 0) && ((control_size-control_count_rec) > 0)){
    mat_rec<-cbind(c(case_count_rec, (casesize-case_count_rec)), c(control_count_rec, (control_size-control_count_rec)))
    dat[i,]$P_REC<-fisher.test(mat_rec, alternative="greater")$p.value
  }
  else{
    dat[i,]$P_REC<-1
  }
  #Recessive model
  case_count_rec2<-dat[i,]$CASE_COUNT_CH+dat[i,]$CASE_COUNT_HOM
  if (control_size > 0){
    control_count_rec2<-dat[i,]$AC_Hom_sum+control_size*(((dat[i,]$AC_Het_sum)/(control_size))^2)
    control_count_rec2<-round(control_count_rec2,0)
    if(control_count_rec2<0){ control_count_rec2<-0}
    if(case_count_rec2>casesize){case_count_rec2<-casesize}
    if (((casesize-case_count_rec2) > 0)&&((control_size-control_count_rec2) > 0)){
      
      dat[i,]$RealREC_case<-case_count_rec2
      dat[i,]$RealREC_case1<-casesize-case_count_rec2
      dat[i,]$RealREC_control<-control_count_rec2
      dat[i,]$RealREC_control1<-control_size-control_count_rec2
      mat_rec<-cbind(c(case_count_rec2, (casesize-case_count_rec2)), c(control_count_rec2, (control_size-control_count_rec2)))
      dat[i,]$P_REC2<-fisher.test(mat_rec, alternative="greater")$p.value
    }
    else{
      dat[i,]$P_REC2<-1
    }
  }
  else{
    dat[i,]$P_REC2<-1
  }
 
}

write.table(x=dat,file=outfile, sep="\t", quote=F, row.names=F, col.names=T)

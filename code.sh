# MiXeR

cd ~
cd shared-genetic-architecture
for i in {1..21}
do
python pleioFDR/python_convert-master/sumstats.py csv --sumstats standard-summary-data-from-R/trait${i} --out MiXeR/data/trait${i}.csv --force --auto
python pleioFDR/python_convert-master/sumstats.py zscore --sumstats MiXeR/data/trait${i}.csv | \
python pleioFDR/python_convert-master/sumstats.py qc --exclude-ranges 6:26000000-34000000 --out MiXeR/data/trait${i}_noMHC.csv --force
gzip MiXeR/data/trait${i}_noMHC.csv  
done

cd MiXeR
for i in {1..22}
do
    python mixer/precimed/mixer.py ld \
        --lib mixer/src/build/lib/libbgmg.so \
        --bfile 1000G_EUR_Phase3_plink/1000G.EUR.QC.${i} \
        --out 1000G_EUR_Phase3_plink/1000G.EUR.QC.${i}.ld \
        --r2min 0.05 --ldscore-r2min 0.05 --ld-window-kb 30000
done

for i in {1..20}
do
python mixer/precimed/mixer.py snps \
   --lib mixer/src/build/lib/libbgmg.so \
   --bim-file 1000G_EUR_Phase3_plink/1000G.EUR.QC.@.bim \
   --ld-file 1000G_EUR_Phase3_plink/1000G.EUR.QC.@.ld \
   --out 1000G_EUR_Phase3_plink/1000G.EUR.QC.prune_maf0p05_rand2M_r2p8.rep${i}.snps \
   --maf 0.05 --subset 2000000 --r2 0.8 --seed $i
done

cd shared-genetic-architecture/MiXeR
for j in {13..15}
do
for i in {1..20}
do
    python mixer/precimed/mixer.py fit1 \
        --trait1-file data/trait${j}_noMHC.csv.gz \
        --out data/trait${j}_noMHC.fit.rep${i} \
        --extract 1000G_EUR_Phase3_plink/1000G.EUR.QC.prune_maf0p05_rand2M_r2p8.rep${i}.snps \
        --bim-file 1000G_EUR_Phase3_plink/1000G.EUR.QC.@.bim \
        --ld-file 1000G_EUR_Phase3_plink/1000G.EUR.QC.@.ld \
        --lib  mixer/src/build/lib/libbgmg.so
    python mixer/precimed/mixer.py test1 \
        --trait1-file data/trait${j}_noMHC.csv.gz \
        --load-params-file data/trait${j}_noMHC.fit.rep${i}.json \
        --out data/trait${j}_noMHC.test.rep${i} \
        --bim-file 1000G_EUR_Phase3_plink/1000G.EUR.QC.@.bim \
        --ld-file 1000G_EUR_Phase3_plink/1000G.EUR.QC.@.ld \
        --lib  mixer/src/build/lib/libbgmg.so
done
done

python mixer/precimed/mixer_figures.py combine --json data/trait1_noMHC.fit.rep@.json --out trait1.fit
python mixer/precimed/mixer_figures.py one --json trait1.fit.json --out trait1 --statistic mean std

cd data
ls *.json > json.txt

cd shared-genetic-architecture/MiXeR
for i in {1..20}
do
python mixer/precimed/mixer.py fit2 \
      --trait1-file data/trait5_noMHC.csv.gz \
      --trait2-file data/trait17_noMHC.csv.gz \
      --trait1-params-file data/trait5_noMHC.fit.rep${i}.json \
      --trait2-params-file data/trait17_noMHC.fit.rep${i}.json \
      --out data/trait5_noMHC_vs_trait17_noMHC.fit.rep${i} \
      --extract 1000G_EUR_Phase3_plink/1000G.EUR.QC.prune_maf0p05_rand2M_r2p8.rep${i}.snps \
      --bim-file 1000G_EUR_Phase3_plink/1000G.EUR.QC.@.bim \
      --ld-file 1000G_EUR_Phase3_plink/1000G.EUR.QC.@.ld \
      --lib  mixer/src/build/lib/libbgmg.so
done

# LDSC
source activate
conda activate py2.7

for i in {1..21}
do
  munge_sumstats.py \
  --sumstats shared-genetic-architecture/standard-summary-data-from-R/no-MHC/trait${i} \
  --chunksize 500000 \
  --merge-alleles shared-genetic-architecture/LDSC/lab_11.13.19/w_hm3.snplist \
  --out shared-genetic-architecture/LDSC/trait${i}
done
# heritability
cd shared-genetic-architecture/LDSC
for i in {1..21}
do
  ldsc.py \
  --h2 trait${i}.sumstats.gz \
  --ref-ld-chr lab_11.13.19/ \
  --w-ld-chr lab_11.13.19/ \
  --out trait${i}_h2
done
# genetic correlation
for j in {1..20}
do
    for ((i=j+1;i<=21;i++))
    do
      ldsc.py \
      --rg trait${j}.sumstats.gz,trait${i}.sumstats.gz \
      --ref-ld-chr lab_11.13.19/ --w-ld-chr lab_11.13.19/ --out trait_${j}_${i}_rg
done
done

# heatmap
setwd("/gpfs/share/home/1610306225/shared-genetic-architecture/data/result")
library(pheatmap)
library(data.table)
rg <- fread("Summary of Genetic Correlation Results.txt")
head(rg);str(rg)
rg_reshape <- data.frame(matrix(NA, nrow=21, ncol=21))
name <- c("CAD","HF","AF","CM","AS","AIS","LAS","CES","SVS","ICH",
  "BMI","SBP","T2D","TG","LDL-C","HDL-C","IMT","ASI","baPWV","bfPWV","cfPWV")
raw_name <- c("trait1.sumstats.gz","trait2.sumstats.gz","trait3.sumstats.gz","trait4.sumstats.gz","trait5.sumstats.gz",
  "trait6.sumstats.gz","trait7.sumstats.gz","trait8.sumstats.gz","trait9.sumstats.gz","trait10.sumstats.gz",
  "trait11.sumstats.gz","trait12.sumstats.gz","trait13.sumstats.gz","trait14.sumstats.gz","trait15.sumstats.gz",
  "trait16.sumstats.gz","trait17.sumstats.gz","trait18.sumstats.gz","trait19.sumstats.gz","trait20.sumstats.gz",
  "trait21.sumstats.gz")
rownames(rg_reshape) <- name
colnames(rg_reshape) <- name
head(rg_reshape)
for(i in 1:21){
  for(j in 1:21){
    if(i==j){rg_reshape[i,j] <- 1}
    else if(i>j){rg_reshape[i,j] <- rg$rg[rg$p1==raw_name[j]&rg$p2==raw_name[i]]}
    else{rg_reshape[i,j] <- rg$rg[rg$p1==raw_name[i]&rg$p2==raw_name[j]]}
  }
}
rg_reshape <- round(rg_reshape,2)
rg_reshape <- as.matrix(rg_reshape)
rg_reshape[rg_reshape>1] <- 1
library(corrplot)
library(dplyr)
col2 <- colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582","#FDDBC7",
                           "#FFFFFF", "#D1E5F0", "#92C5DE","#4393C3", "#2166AC",
                           "#053061"))(100) %>% rev()
corrplot(rg_reshape, is.corr = FALSE,
         method="circle", type='upper', outline = FALSE,
         tl.pos = 'd', tl.col = "black", tl.offset = 1, tl.cex = 1,
         cl.pos = 'r',
         addCoef.col = T,
         order = "original",
         number.font = 6, col = col2)



# HDL

head /gpfs/share/home/1610306225/shared-genetic-architecture/standard-summary-data-from-R/trait1
for i in {1..21}
do
Rscript /gpfs/share/home/1610306225/shared-genetic-architecture/HDL/HDL-R/HDL.data.wrangling.R \
gwas.file=/gpfs/share/home/1610306225/shared-genetic-architecture/standard-summary-data-from-R/no-MHC/trait${i} \
LD.path=/gpfs/share/home/1610306225/shared-genetic-architecture/HDL/UKB_imputed_SVD_eigen99_extraction \
SNP=rsID A1=A1 A2=A2 N=N b=BETA se=SE \
output.file=/gpfs/share/home/1610306225/shared-genetic-architecture/HDL/standard-summary-data-from-R-HDL/trait${i} \
log.file=/gpfs/share/home/1610306225/shared-genetic-architecture/HDL/standard-summary-data-from-R-HDL/trait${i}
done

# heritability
for i in {1..21}
do
Rscript /gpfs/share/home/1610306225/shared-genetic-architecture/HDL/HDL-R/HDL.run.R \
gwas.df=/gpfs/share/home/1610306225/shared-genetic-architecture/HDL/standard-summary-data-from-R-HDL/trait${i}.hdl.rds \
LD.path=/gpfs/share/home/1610306225/shared-genetic-architecture/HDL/UKB_imputed_SVD_eigen99_extraction \
output.file=/gpfs/share/home/1610306225/shared-genetic-architecture/HDL/HDL-output/trait${i}-h2.Rout
done

# genetic correlation
for j in {1..21}
do
    for ((i=j+1;i<=21;i++))
    do
    Rscript /gpfs/share/home/1610306225/shared-genetic-architecture/HDL/HDL-R/HDL.run.R \
    gwas1.df=/gpfs/share/home/1610306225/shared-genetic-architecture/HDL/standard-summary-data-from-R-HDL/trait${j}.hdl.rds \
    gwas2.df=/gpfs/share/home/1610306225/shared-genetic-architecture/HDL/standard-summary-data-from-R-HDL/trait${i}.hdl.rds \
    LD.path=/gpfs/share/home/1610306225/shared-genetic-architecture/HDL/UKB_imputed_SVD_eigen99_extraction \
    output.file=/gpfs/share/home/1610306225/shared-genetic-architecture/HDL/HDL-output/trait${j}_trait${i}-rg.Rout
    done
done

# pleioFDR

cd /gpfs/share/home/1610306225/shared-genetic-architecture
for i in {11..21}
do
    python pleioFDR/python_convert-master/sumstats.py csv --auto --sumstats standard-summary-data-from-R/trait${i} --out pleioFDR/standard-summary-data-from-R-pleiofdr/trait${i}.csv --force
    python pleioFDR/python_convert-master/sumstats.py zscore --sumstats pleioFDR/standard-summary-data-from-R-pleiofdr/trait${i}.csv --out pleioFDR/standard-summary-data-from-R-pleiofdr/trait${i}_z.csv 
    python pleioFDR/python_convert-master/sumstats.py mat --sumstats pleioFDR/standard-summary-data-from-R-pleiofdr/trait${i}_z.csv --ref pleioFDR/9545380.ref --out pleioFDR/standard-summary-data-from-R-pleiofdr/trait${i}.mat
done

cd /gpfs/share/home/1610306225/shared-genetic-architecture/pleioFDR/pleiofdr-master
cp config_default.txt config.txt
module load matlab/R2019a

traitname_array=(NA CAD HF AF MI AS AIS LAS CES SVS PAD BMI SBP T2D TG LDL HDL IMT ASI baPWV bfPWV cfPWV)
for i in {1..10}
do
    for j in {11..21}
    do
        sed -i '15c traitfile1='trait${i}.mat'' config.txt 
        sed -i '17c traitfiles={"'trait${j}.mat'"}' config.txt 
        sed -i '16c traitname1='${traitname_array[i]}'' config.txt 
        sed -i '18c traitnames={"'${traitname_array[j]}'"}' config.txt 
        sed -i '21c outputdir=/gpfs/share/home/1610306225/shared-genetic-architecture/pleioFDR/output/'${traitname_array[i]}'_'${traitname_array[j]}'_cond' config.txt   
        matlab -nodisplay -nosplash < runme.m
    done
done

traitname_array=(NA CAD HF AF MI PAD AS AIS LAS CES SVS)
for i in {1..5}
do
    for ((j=i+1;i<=10;i++))
    do
        if [ j>=5 ];then
        sed -i '15c traitfile1='trait${i}.mat'' config.txt 
        sed -i '17c traitfiles={"'trait${j}.mat'"}' config.txt 
        sed -i '16c traitname1='${traitname_array[i]}'' config.txt 
        sed -i '18c traitnames={"'${traitname_array[j]}'"}' config.txt 
        sed -i '21c outputdir=/gpfs/share/home/1610306225/shared-genetic-architecture/pleioFDR/output/'${traitname_array[i]}'_'${traitname_array[j]}'_cond' config.txt   
        matlab -nodisplay -nosplash < runme.m
        else 
        echo "no need"
    done
done

# manhattan plots
library(CMplot)
library(data.table)
setwd("/gpfs/share/home/1610306225/shared-genetic-architecture/pleioFDR/out_beta")

traitname <- c("CAD","HF","AF","MI","AS","AIS","LAS","CES",
"SVS","PAD","BMI","SBP","T2D","TG","LDL","HDL","IMT","ASI","baPWV","bfPWV","cfPWV")

for (u in 1:10){ #0.01
for(v in 11:17){
    dat <- fread(paste0("beta_",traitname[u],"_",traitname[v],".csv"))
    figdat <- dat[,c("SNP","CHR","BP","FDR")]
    colnames(figdat) <- c("SNP","Chromosome","Position","FDR")
    SNPs <- c(
        dat$SNP[dat$direction==1&dat$FDR<0.01&dat$lead==0],
        dat$SNP[dat$direction==0&dat$FDR<0.01&dat$lead==0],
        dat$SNP[dat$direction==1&dat$lead==1],
        dat$SNP[dat$direction==0&dat$lead==1]
    )
    color <- c(
    rep("#00b0eb",length(dat$SNP[dat$direction==1&dat$FDR<0.01&dat$lead==0])),
    rep("#ffd401",length(dat$SNP[dat$direction==0&dat$FDR<0.01&dat$lead==0])),
    rep("black",length(dat$SNP[dat$direction==1&dat$lead==1])),
    rep("red",length(dat$SNP[dat$direction==0&dat$lead==1]))
    )
    pch <- c(
    rep(16,length(dat$SNP[dat$direction==1&dat$FDR<0.01&dat$lead==0])),
    rep(16,length(dat$SNP[dat$direction==0&dat$FDR<0.01&dat$lead==0])),
    rep(18,length(dat$SNP[dat$direction==1&dat$lead==1])),
    rep(18,length(dat$SNP[dat$direction==0&dat$lead==1]))
    )
    cex<- c(
    rep(0.8,length(dat$SNP[dat$direction==1&dat$FDR<0.01&dat$lead==0])),
    rep(0.8,length(dat$SNP[dat$direction==0&dat$FDR<0.01&dat$lead==0])),
    rep(1.2,length(dat$SNP[dat$direction==1&dat$lead==1])),
    rep(1.2,length(dat$SNP[dat$direction==0&dat$lead==1]))
    )

    CMplot(figdat, type="p", plot.type="m", 
        threshold=1e-2, # 注意修改
        threshold.lty=2, col="grey90",
        threshold.lwd=2, cex=0.4,
        threshold.col="green", 
        amplify=TRUE, bin.size=1e6,
        #chr.den.col=c("darkgreen", "yellow", "red"), 
        signal.col="green",
        signal.cex=0.6, file="jpg", dpi=300, file.output=TRUE, verbose=TRUE, LOG10 = T,file.name=paste0(traitname[u],"_",traitname[v]),
        highlight=SNPs,  highlight.col = color, highlight.pch = pch, highlight.cex =cex)
    rm(dat,figdat)
}}

# s-LDSC

cd /gpfs/share/home/1610306225/shared-genetic-architecture

wc -l pleioFDR/output_cond/CAD_BMI_cond/CAD_BMI_condfdr_0.01_all.csv

traitname_array=(NA CAD HF AF MI AS AIS LAS CES SVS PAD BMI SBP T2D TG LDL HDL IMT baPWV)
for i in {1..10}
do
    for j in {11..18}
    do
        cd /gpfs/share/home/1610306225/shared-genetic-architecture/pleioFDR/output_cond/"${traitname_array[i]}"_"${traitname_array[j]}"_cond
        wc -l ./"${traitname_array[i]}"_"${traitname_array[j]}"_condfdr_0.01_all.csv
    done
done

R
library(data.table)
setDTthreads(threads = 0)
getDTthreads(verbose = getOption("datatable.verbose"))
setwd("/gpfs/share/home/1610306225/shared-genetic-architecture")
rm(list=ls())

for (j in 1:22){
	traitname <- c(
		"CAD_BMI","CAD_SBP","CAD_T2D","CAD_TG","CAD_LDL","CAD_HDL","CAD_IMT","CAD_baPWV")
		# "HF_BMI","HF_SBP","HF_T2D","HF_TG","HF_LDL","HF_HDL","HF_IMT","HF_baPWV")
		# "AF_BMI","AF_SBP","AF_T2D","AF_TG","AF_LDL","AF_HDL","AF_IMT","AF_baPWV")
		# "MI_BMI","MI_SBP","MI_T2D","MI_TG","MI_LDL","MI_HDL","MI_IMT","MI_baPWV")
		# "AS_BMI","AS_SBP","AS_T2D","AS_TG","AS_LDL","AS_HDL","AS_IMT","AS_baPWV")
		# "AIS_BMI","AIS_SBP","AIS_T2D","AIS_TG","AIS_LDL","AIS_HDL","AIS_IMT","AIS_baPWV")
		# "LAS_BMI","LAS_SBP","LAS_T2D","LAS_TG","LAS_LDL","LAS_HDL","LAS_IMT","LAS_baPWV")
		# "CES_BMI","CES_SBP","CES_T2D","CES_TG","CES_LDL","CES_HDL","CES_IMT","CES_baPWV")
		# "SVS_BMI","SVS_SBP","SVS_T2D","SVS_TG","SVS_LDL","SVS_HDL","SVS_IMT","SVS_baPWV")
		# "PAD_BMI","PAD_SBP","PAD_T2D","PAD_TG","PAD_LDL","PAD_HDL","PAD_IMT","PAD_baPWV")
	snp <- fread(paste0("MiXeR/1000G_EUR_Phase3_plink/1000G.EUR.QC.",j,".bim"))
	snp <- snp[,4]
	for (i in 1:8){
		annot <- fread(paste0("./pleioFDR/output_cond/",traitname[i],"_cond/",traitname[i],"_condfdr_0.01_all.csv"))
		annot <- annot[annot$chrnum==j,c(5,8)]
		snp_annot <- merge(snp,annot,by.x="V4",by.y="chrpos",all.x=T)
		snp_annot <- data.frame(snp_annot)
		snp_annot[!is.na(snp_annot[,i+1]),i+1] <- 1
		snp_annot[is.na(snp_annot[,i+1]),i+1] <- 0
		snp <- snp_annot
		}
	snp_annot[,1] <- 1
	fwrite(snp_annot,paste0("./s-LDSC/LD-annot/cad.",j,".annot"),sep = "\t",quote = F,row.names = F,col.names = T)
	rm(list=ls())
}


source activate
conda activate py2.7
cd /gpfs/share/home/1610306225/shared-genetic-architecture/s-LDSC

for i in {1..22}
do
ldsc.py --l2 --bfile 1000G_EUR_Phase3_plink/1000G.EUR.QC.${i} --ld-wind-cm 1 --annot LD-annot/cad.${i}.annot --thin-annot --out LD-annot/cad.${i} --print-snps hapmap3_snps/hm.${i}.snp
done

ldsc.py --h2 /gpfs/share/home/1610306225/shared-genetic-architecture/LDSC/trait1.sumstats.gz --ref-ld-chr LD-annot/cad. --w-ld-chr example/weights_hm3_no_hla/weights. --overlap-annot --not-M-5-50 --out output/CAD
# ldsc.py --h2 /gpfs/share/home/1610306225/shared-genetic-architecture/LDSC/trait2.sumstats.gz --ref-ld-chr LD-annot/hf. --w-ld-chr example/weights_hm3_no_hla/weights. --overlap-annot --not-M-5-50 --out output/HF
# ldsc.py --h2 /gpfs/share/home/1610306225/shared-genetic-architecture/LDSC/trait3.sumstats.gz --ref-ld-chr LD-annot/af. --w-ld-chr example/weights_hm3_no_hla/weights. --overlap-annot --not-M-5-50 --out output/AF
# ldsc.py --h2 /gpfs/share/home/1610306225/shared-genetic-architecture/LDSC/trait4.sumstats.gz --ref-ld-chr LD-annot/mi. --w-ld-chr example/weights_hm3_no_hla/weights. --overlap-annot --not-M-5-50 --out output/MI
# ldsc.py --h2 /gpfs/share/home/1610306225/shared-genetic-architecture/LDSC/trait5.sumstats.gz --ref-ld-chr LD-annot/as. --w-ld-chr example/weights_hm3_no_hla/weights. --overlap-annot --not-M-5-50 --out output/AS
# ldsc.py --h2 /gpfs/share/home/1610306225/shared-genetic-architecture/LDSC/trait6.sumstats.gz --ref-ld-chr LD-annot/ais. --w-ld-chr example/weights_hm3_no_hla/weights. --overlap-annot --not-M-5-50 --out output/AIS
# ldsc.py --h2 /gpfs/share/home/1610306225/shared-genetic-architecture/LDSC/trait7.sumstats.gz --ref-ld-chr LD-annot/las. --w-ld-chr example/weights_hm3_no_hla/weights. --overlap-annot --not-M-5-50 --out output/LAS
# ldsc.py --h2 /gpfs/share/home/1610306225/shared-genetic-architecture/LDSC/trait8.sumstats.gz --ref-ld-chr LD-annot/ces. --w-ld-chr example/weights_hm3_no_hla/weights. --overlap-annot --not-M-5-50 --out output/CES
# ldsc.py --h2 /gpfs/share/home/1610306225/shared-genetic-architecture/LDSC/trait9.sumstats.gz --ref-ld-chr LD-annot/svs. --w-ld-chr example/weights_hm3_no_hla/weights. --overlap-annot --not-M-5-50 --out output/SVS
# ldsc.py --h2 /gpfs/share/home/1610306225/shared-genetic-architecture/LDSC/trait10.sumstats.gz --ref-ld-chr LD-annot/pad. --w-ld-chr example/weights_hm3_no_hla/weights. --overlap-annot --not-M-5-50 --out output/PAD


R
library(data.table)
setDTthreads(threads = 0)
getDTthreads(verbose = getOption("datatable.verbose"))
setwd("/gpfs/share/home/1610306225/shared-genetic-architecture")
rm(list=ls())

tonewsum<-function(x,y){
	traitname <- c("CAD","HF","AF","MI","AS","AIS","LAS","CES","SVS","PAD",
                "BMI","SBP","T2D","TG","LDL","HDL","IMT","ASI","baPWV","bfPWV","cfPWV")
	dat1 <- fread(paste0("/gpfs/share/home/1610306225/shared-genetic-architecture/standard-summary-data-from-R/trait",x))
	dat2 <- fread(paste0("/gpfs/share/home/1610306225/shared-genetic-architecture/standard-summary-data-from-R/trait",y))
	dat2 <- dat2[,c("CHR","BP")]
	test <- merge(dat1,dat2,by=c("CHR","BP"))
	dat3 <- fread(paste0("/gpfs/share/home/1610306225/shared-genetic-architecture/pleioFDR/output_cond/",traitname[x],"_",traitname[y],"_cond/result.mat.csv"))
	dat4 <- dat3[is.na(dat3$FDR)==F,]
	dat5 <- dat4[,c(1,3)]
	test <- merge(test,dat5,by=c("CHR","BP")) 
	fwrite(test,paste0("/gpfs/share/home/1610306225/shared-genetic-architecture/s-LDSC/trait-summary/",traitname[x],"_",traitname[y]),sep = "\t",quote = F,row.names = F,col.names = T)
	rm(list=ls())
}	

for (j in 1:10){
for(i in 11:17){
    tonewsum(j,i)
}
}

source activate
conda activate py2.7

traitname_array=(NA CAD HF AF MI AS AIS LAS CES SVS PAD BMI SBP T2D TG LDL HDL IMT baPWV)
for i in {1..10}
do
    for j in {11..17}
    do
		munge_sumstats.py \
		--sumstats shared-genetic-architecture/s-LDSC/trait-summary/"${traitname_array[i]}"_"${traitname_array[j]}" \
		--chunksize 500000 \
		--merge-alleles shared-genetic-architecture/LDSC/lab_11.13.19/w_hm3.snplist \
		--out shared-genetic-architecture/s-LDSC/trait-summary/"${traitname_array[i]}"_"${traitname_array[j]}"
    done
done

traitname_small=(NA cad hf af mi as ais las ces svs pad)
cd /gpfs/share/home/1610306225/shared-genetic-architecture/s-LDSC
for i in {1..10}
do
    for j in {11..17}
    do
		ldsc.py \
		--h2 /gpfs/share/home/1610306225/shared-genetic-architecture/s-LDSC/trait-summary/"${traitname_array[i]}"_"${traitname_array[j]}".sumstats.gz \
		--ref-ld-chr /gpfs/share/home/1610306225/shared-genetic-architecture/s-LDSC/LD-annot/"${traitname_small[i]}". \
		--w-ld-chr /gpfs/share/home/1610306225/shared-genetic-architecture/s-LDSC/example/weights_hm3_no_hla/weights. \
		--overlap-annot --not-M-5-50 \
		--out /gpfs/share/home/1610306225/shared-genetic-architecture/s-LDSC/output/"${traitname_array[i]}"_"${traitname_array[j]}"
    done
done


# COLOC-SUISE

source activate R4.2
R

library(coloc)
library(dplyr)
library(data.table)
library(LDlinkR)
library(Rfast)
library(ieugwasr)
traitname <- c("CAD","HF","AF","MI","AS","AIS")

for (i in 1){      
        for (j in 5){
            block <- fread(paste0("/gpfs/share/home/1610306225/shared-genetic-architecture/pleioFDR/clump/",traitname[i],"_",traitname[j],"_conj_result.clump.loci.csv"))
            block$nsnps_raw <- NA 
            block$bpadd <- NA; block$nsnps <- NA
            block$lower_coverage_s1 <- NA; block$lower_coverage_s2 <- NA
            block$s1_num <- NA; block$s2_num <-NA
            block$PP.H4 <- NA 
            block$sig.pp4.num <- NA
            block$cred.snp <- NA

            trait1 <- fread(paste0("/gpfs/share/home/1610306225/shared-genetic-architecture/standard-summary-data-from-R/trait",i,"_1000G_allele"))
            s1=mean(trait1$Ncases,na.rm=T)/mean(trait1$N,na.rm=T)
            N1=ceiling(mean(trait1$N))
            trait1$rs_id <- paste(trait1$rsID,trait1$A1,trait1$A2,sep="_")
            trait1$rs_id_rev <- paste(trait1$rsID,trait1$A2,trait1$A1,sep="_")
            trait1 <- trait1[,c("rsID","CHR","BP","BETA","SE","P","MAF","rs_id","rs_id_rev")]
            trait1$MAF[trait1$MAF<=0] <- 0.0001
            trait1$MAF[trait1$MAF>=1] <- 0.99
            trait1$varbeta <- trait1$SE^2
            colnames(trait1)[c(4,5,6)] <- c("beta","sdY","pval_nominal") 
            
            trait2 <- fread(paste0("/gpfs/share/home/1610306225/shared-genetic-architecture/standard-summary-data-from-R/trait",j,"_1000G_allele"))
            N2=ceiling(mean(trait2$N))
            trait2$rs_id <- paste(trait2$rsID,trait2$A1,trait2$A2,sep="_")
            trait2$rs_id_rev <- paste(trait2$rsID,trait2$A2,trait2$A1,sep="_")    
            trait2 <- trait2[,c("rsID","CHR","BP","BETA","SE","P","rs_id","rs_id_rev")]
            trait2$varbeta <- trait2$SE^2
            colnames(trait2)[c(4,5,6)] <- c("beta","sdY","pval_nominal") 
            trait1 <- trait1[!duplicated(trait1$rsID),]; trait2 <- trait2[!duplicated(trait2$rsID),]
            head(trait1); head(trait2)

            for (l in 1:nrow(block)){
                subset1 <- trait1[trait1$CHR==block$CHR[l] & trait1$BP<=block$MaxBP[l] & trait1$BP>=block$MinBP[l],]
                subset2 <- trait2[trait2$CHR==block$CHR[l] & trait2$BP<=block$MaxBP[l] & trait2$BP>=block$MinBP[l],]
                subset2 <- subset2[,-c("CHR","BP")]
                subset1 <- subset1[complete.cases(subset1),]; subset2 <- subset2[complete.cases(subset2),]
                input <- merge(subset1, subset2, by="rsID", all=FALSE, suffixes=c("_trait","_rf"))
                head(input)
                block$nsnps_raw[l] <- nrow(input)

                if(nrow(input)<100){
                    block$MaxBP[l] <- block$MaxBP[l]+100000
                    block$MinBP[l] <- block$MinBP[l]-100000
                    block$bpadd[l] <- "T"
                    subset1 <- trait1[trait1$CHR==block$CHR[l] & trait1$BP<=block$MaxBP[l] & trait1$BP>=block$MinBP[l],]
                    subset2 <- trait2[trait2$CHR==block$CHR[l] & trait2$BP<=block$MaxBP[l] & trait2$BP>=block$MinBP[l],]
                    subset2 <- subset2[,-c("CHR","BP")]
                    subset1 <- subset1[complete.cases(subset1),]; subset2 <- subset2[complete.cases(subset2),]
                    input <- merge(subset1, subset2, by="rsID", all=FALSE, suffixes=c("_trait","_rf"))
                }

                ld <- ld_matrix(
                input$rsID, 
                plink_bin = "/gpfs/share/home/1610306225/software/miniconda3/envs/R4.2/lib/R/library/genetics.binaRies/bin/plink",
            
                bfile = "/gpfs/share/home/1610306225/shared-genetic-architecture/coloc_susie/ld_ref/EUR/1000G.EUR.QC")
                
                input$beta_trait[which(rownames(ld) %in% input$rs_id_rev_trait)] <- -input$beta_trait[which(rownames(ld) %in% input$rs_id_rev_trait)]
                input$beta_rf[which(rownames(ld) %in% input$rs_id_rev_trait)] <- -input$beta_rf[which(rownames(ld) %in% input$rs_id_rev_trait)]
                ld_rs <- ld_matrix(
                input$rsID, with_alleles = FALSE,
                plink_bin = "/gpfs/share/home/1610306225/software/miniconda3/envs/R4.2/lib/R/library/genetics.binaRies/bin/plink",
                
                bfile = "/gpfs/share/home/1610306225/shared-genetic-architecture/coloc_susie/ld_ref/EUR/1000G.EUR.QC")
                
                input <- input[input$rsID %in% colnames(ld_rs),]
                block$nsnps[l] <- nrow(input)

                D1 <- list(snp=input$rsID,pvalues=input$pval_nominal_trait,beta=input$beta_trait,varbeta=input$varbeta_trait, sdY= input$sdY_trait, N=N1, s=s1, type="cc", LD=ld_rs) # check_dataset(D1,req="LD")
                D2 <- list(snp=input$rsID,pvalues=input$pval_nominal_rf,beta=input$beta_rf,varbeta=input$varbeta_rf, sdY= input$sdY_rf, N=N2, type="quant", LD=ld_rs) # check_dataset(D2,req="LD")
     
                
                temp1 <- try(runsusie(D1), silent = T)
                if(class(temp1) == "try-error"){S1 <- "bushoulian"}else{S1 <- temp1}
                
                if(class(S1)=="character"){l=l+1}
                else{
                    if(is.null(summary(S1)$cs)==T){
                        temp1 <- try(runsusie(D1,coverage=0.1), silent = T)
                        if(class(temp1) == "try-error"){S1 <- "bushoulian"}else{S1 <- temp1}
                        block$lower_coverage_s1[l] <- T
                        }
                    if(class(S1)=="character"){l=l+1}
                    else if(is.null(summary(S1)$cs)==T){l=l+1}
                    else {
                        temp2 <- try(runsusie(D2), silent = T)
                        if(class(temp2) == "try-error"){S2 <- "bushoulian"}else{S2 <- temp2}
                        if(class(S2)=="character"){l=l+1}
                        else{
                            if(is.null(summary(S2)$cs)==T){
                            temp2 <- try(runsusie(D2,coverage=0.1), silent = T)
                            if(class(temp2) == "try-error"){S2 <- "bushoulian"}else{S2 <- temp2}
                            block$lower_coverage_s2[l] <- T
                            }
                            if(class(S2)=="character"){l=l+1}
                            else if(is.null(summary(S2)$cs)==T){l=l+1}
                            else{
                                block$s1_num[l] <- nrow(summary(S1)$cs)
                                block$s2_num[l] <- nrow(summary(S2)$cs)
                                
                                res=coloc.susie(S1,S2)                     
                                block$PP.H4[l] <- max(res$summary$PP.H4.abf)
                                if(max(res$summary$PP.H4.abf)>=0.8){
                                    num <- which(res$summary$PP.H4.abf>=0.8)
                                    block$sig.pp4.num[l] <- length(num)
                                    snp <- c()
                                    snp.PPA <- data.frame(res$results)
                                    for(n in num){
                                        o <- order(snp.PPA[,n+1],decreasing=TRUE)
                                        cs <- cumsum(snp.PPA[,n+1][o])
                                        w <- which(cs > 0.95)[1]
                                        snp <- c(snp,snp.PPA[o,][1:w,]$snp)}
                                        cred.snp <- c()
                                        for(k in 1:length(snp)){cred.snp <- paste(cred.snp,snp[k])}     
                                        block$cred.snp[l] <- cred.snp
                                    }
                            }
                        }   
                    }
                }             
                fwrite(block,paste0("/gpfs/share/home/1610306225/shared-genetic-architecture/coloc_susie/",traitname[i],"_",traitname[j],"_conj_clump_loci_coloc.csv"),sep = ",",quote = F,row.names = F,col.names = T)
            }
    }                       
}


# MOLOC
source activate R4.2
R

library(moloc)
options(scipen = 1, digits = 2)
library(data.table)
setDTthreads(threads = 0)
getDTthreads(verbose = getOption("datatable.verbose"))
rm(list=ls())
traitname <- c("CAD","HF","AF","MI","AS","AIS","LAS","CES",
                "SVS","PAD","BMI","SBP","T2D","TG","LDL","HDL","IMT")
block <- fread("/gpfs/share/home/1610306225/shared-genetic-architecture/moloc/LD_merge_locus_rf.csv")
block <- block[-1,]
colnames(block)[6] <- "traits_combination" 
output <- data.frame(matrix(NA,1,13))
colnames(output)[1:6] <- colnames(block)[1:6]
colnames(output)[7:13] <- c("prior","sumbf","logBF_locus","PPA","coloc_ppas","best.snp.coloc","nsnp")

for (i in 1:9){
    for (j in 10){ 
        for (k in 13){
            trait1 <- fread(paste0("/gpfs/share/home/1610306225/shared-genetic-architecture/standard-summary-data-from-R/trait",i,"_1000G"))
            trait1 <- trait1[,c("rsID","CHR","BP","BETA","SE","MAF","N","Ncases")]
            colnames(trait1)[c(1,3)] <- c("SNP","POS") 
            trait2 <- fread(paste0("/gpfs/share/home/1610306225/shared-genetic-architecture/standard-summary-data-from-R/trait",j,"_1000G"))
            trait2 <- trait2[,c("rsID","CHR","BP","BETA","SE","MAF","N","Ncases")]
            colnames(trait2)[c(1,3)] <- c("SNP","POS") 
            trait3 <- fread(paste0("/gpfs/share/home/1610306225/shared-genetic-architecture/standard-summary-data-from-R/trait",k,"_1000G"))
            # trait3 <- trait3[,c("rsID","CHR","BP","BETA","SE","MAF","N")] 
            trait3 <- trait3[,c("rsID","CHR","BP","BETA","SE","MAF","N","Ncases")] #T2D
            colnames(trait3)[c(1,3)] <- c("SNP","POS")
            trait1 <- trait1[!duplicated(trait1$SNP),]; trait2 <- trait2[!duplicated(trait2$SNP),]; trait3 <- trait3[!duplicated(trait3$SNP),]
            block_subset <- block[block$traits_combination==paste0(traitname[j],"_",traitname[i],"_",traitname[k]),]
            if(nrow(block_subset)==0){
                output <- output
            }else{
                op <- data.frame(matrix(NA,1,7))
                colnames(op) <- c("prior","sumbf","logBF_locus","PPA","coloc_ppas","best.snp.coloc","nsnp")
                for (l in 1:nrow(block_subset)){
                    subset1 <- trait1[trait1$CHR==block_subset$chr[l] & trait1$POS<=block_subset$locus_max[l] & trait1$POS>=block_subset$locus_min[l],]
                    subset2 <- trait2[trait2$CHR==block_subset$chr[l] & trait2$POS<=block_subset$locus_max[l] & trait2$POS>=block_subset$locus_min[l],]
                    subset3 <- trait3[trait3$CHR==block_subset$chr[l] & trait3$POS<=block_subset$locus_max[l] & trait3$POS>=block_subset$locus_min[l],]
                    subset1 <- subset1[complete.cases(subset1),]; subset2 <- subset2[complete.cases(subset2),]; subset3 <- subset3[complete.cases(subset3),]
                    if(nrow(subset1)<3|nrow(subset2)<3|nrow(subset3)<3){
                        output <- output
                    }else{
                        dat <- list(subset1,subset2,subset3)
                        
                        moloc <- moloc_test(dat, prior_var=c(0.01, 0.1, 0.5), priors=c(1e-04, 1e-06, 1e-07))

                        Posteriors <- moloc[[1]][rownames(moloc[[1]])=="abc",]
                        SNP <- moloc[[3]][rownames(moloc[[3]])=="abc",]
                        rownames(Posteriors) <- ""; rownames(SNP) <- ""
                        nsnp <- data.frame(moloc[[2]])
                        colnames(nsnp) <- "nsnp"
                        add <- cbind(Posteriors,SNP,nsnp)
                        op <- rbind(op,add)
                    }
                }
                op <- op[-1,]
                output_subset <- cbind(block_subset,op)
                output <- rbind(output,output_subset)
            }
        }
    }
}
fwrite(output,"/gpfs/share/home/1610306225/shared-genetic-architecture/moloc/moloc.xls",sep = "\t",quote = F,row.names = F,col.names = T)


# OPERA
cd /gpfs/share/home/1610306225/shared-genetic-architecture/OPERA/QTL/sQTL
Adipose_Subcutaneous Adipose_Visceral_Omentum Adrenal_Gland Artery_Aorta Artery_Coronary Artery_Tibial 
Brain_Amygdala Brain_Anterior_cingulate_cortex_BA24 Brain_Caudate_basal_ganglia Brain_Cerebellar_Hemisphere 
Brain_Cerebellum Brain_Cortex Brain_Frontal_Cortex_BA9 Brain_Hippocampus Brain_Hypothalamus 
Brain_Nucleus_accumbens_basal_ganglia Brain_Putamen_basal_ganglia Brain_Spinal_cord_cervical_c-1 
Brain_Substantia_nigra Breast_Mammary_Tissue 
Cells_Cultured_fibroblasts Cells_EBV-transformed_lymphocytes Colon_Sigmoid Colon_Transverse 
Esophagus_Gastroesophageal_Junction Esophagus_Mucosa Esophagus_Muscularis 
Heart_Atrial_Appendage Heart_Left_Ventricle 
Kidney_Cortex Liver Lung 
Minor_Salivary_Gland Muscle_Skeletal Nerve_Tibial Ovary Pancreas Pituitary Prostate 
Skin_Not_Sun_Exposed_Suprapubic Skin_Sun_Exposed_Lower_leg Small_Intestine_Terminal_Ileum Spleen Stomach 
Testis Thyroid Uterus Vagina Whole_Blood 
# mqTL
cd /gpfs/share/home/1610306225/shared-genetic-architecture/OPERA/QTL/mQTL
wget https://yanglab.westlake.edu.cn/data/SMR/LBC_BSGS_meta.tar.gz

# GWAS
source activate R4.2
R
library(data.table)
setDTthreads(threads = 0)
getDTthreads(verbose = getOption("datatable.verbose"))
traitname <- c("CAD","HF","AF","MI","AS","AIS","LAS","CES","SVS","PAD",
                "BMI","SBP","T2D","TG","LDL","HDL","IMT")

for(i in 1:17){
    dat <- fread(paste0("/gpfs/share/home/1610306225/shared-genetic-architecture/standard-summary-data-from-R/trait",i))
    # head(dat)
    # SNP A1 A2 freq b se p N 
    # Note: 1) For a case-control study, the effect size should be log(odds ratio) with its corresponding standard error. 
    dat <- dat[,c("rsID","A1","A2","EAF","BETA","SE","P","N")]
    colnames(dat) <- c("SNP","A1","A2","freq","b","se","p","N")
    fwrite(dat,paste0("/gpfs/share/home/1610306225/shared-genetic-architecture/OPERA/GWAS/",traitname[i],".ma"),sep = "\t",quote = F,row.names = F,col.names = T)
}

cd /gpfs/share/home/1610306225/shared-genetic-architecture/OPERA/QTL/eQTL
unzip '*.zip'
cd /gpfs/share/home/1610306225/shared-genetic-architecture/OPERA/QTL/sQTL
unzip '*.zip'
cd /gpfs/share/home/1610306225/shared-genetic-architecture/OPERA/QTL/mQTL
tar –zxvf /gpfs/share/home/1610306225/shared-genetic-architecture/OPERA/QTL/mQTL/LBC_BSGS_meta.tar.gz

cd /gpfs/share/home/1610306225/shared-genetic-architecture/OPERA

traitname_array=(NA CAD HF AF MI AS AIS LAS CES SVS PAD BMI SBP T2D TG LDL HDL IMT ASI baPWV bfPWV cfPWV)
for i in {1..10}
do
    for j in {11..21}
    do
    python python_convert-master/fdrmat2csv.py --mat output_cond/"${traitname_array[i]}"_"${traitname_array[j]}"_cond/result.mat --ref 9545380.ref
    done
done

cd /gpfs/share/home/1610306225/shared-genetic-architecture/OPERA/QTL/sQTL
tissue_name=(Adipose_Subcutaneous Adipose_Visceral_Omentum Adrenal_Gland Artery_Aorta Artery_Coronary Artery_Tibial Brain_Amygdala Brain_Anterior_cingulate_cortex_BA24 Brain_Caudate_basal_ganglia Brain_Cerebellar_Hemisphere Brain_Cerebellum Brain_Cortex Brain_Frontal_Cortex_BA9 Brain_Hippocampus Brain_Hypothalamus Brain_Nucleus_accumbens_basal_ganglia Brain_Putamen_basal_ganglia Brain_Spinal_cord_cervical_c-1 Brain_Substantia_nigra Breast_Mammary_Tissue Cells_Cultured_fibroblasts Cells_EBV-transformed_lymphocytes Colon_Sigmoid Colon_Transverse Esophagus_Gastroesophageal_Junction Esophagus_Mucosa Esophagus_Muscularis Heart_Atrial_Appendage Heart_Left_Ventricle Kidney_Cortex Liver Lung Minor_Salivary_Gland Muscle_Skeletal Nerve_Tibial Ovary Pancreas Pituitary Prostate Skin_Not_Sun_Exposed_Suprapubic Skin_Sun_Exposed_Lower_leg Small_Intestine_Terminal_Ileum Spleen Stomach Testis Thyroid Uterus Vagina Whole_Blood)
for i in {1..49}
do
    touch "${tissue_name[i]}"_sqtl-chr-merge.list
    for j in {1..22}
    do
        echo "/gpfs/share/home/1610306225/shared-genetic-architecture/OPERA/QTL/sQTL/sQTL_"${tissue_name[i]}"/chr${j} >> "${tissue_name[i]}"_sqtl-chr-merge.list
    done
    /gpfs/share/home/1610306225/shared-genetic-architecture/OPERA/smr-1.3.1-linux-x86_64/smr-1.3.1 --besd-flist /gpfs/share/home/1610306225/shared-genetic-architecture/OPERA/sqtl-chr-merge.list --make-besd-dense --out /gpfs/share/home/1610306225/shared-genetic-architecture/OPERA/QTL/sQTL/sQTL_"${tissue_name[i]}"/
done

# mQTL
touch mqtl-chr-merge.list
for i in {1..22}
do
echo "/gpfs/share/home/1610306225/shared-genetic-architecture/OPERA/QTL/mQTL/LBC_BSGS_meta/bl_mqtl_chr${i}" >> mqtl-chr-merge.list
done
cat mqtl-chr-merge.list
/gpfs/share/home/1610306225/shared-genetic-architecture/OPERA/smr-1.3.1-linux-x86_64/smr-1.3.1 --besd-flist /gpfs/share/home/1610306225/shared-genetic-architecture/OPERA/mqtl-chr-merge.list --make-besd-dense --out /gpfs/share/home/1610306225/shared-genetic-architecture/OPERA/QTL/mQTL/LBC_BSGS_meta/bl_mqtl

# Run OPERA for stage 1 analysis
cd /gpfs/share/home/1610306225/shared-genetic-architecture/OPERA
./opera_Linux/opera_Linux --besd-flist /gpfs/share/home/1610306225/shared-genetic-architecture/OPERA/qtl-list --gwas-summary /gpfs/share/home/1610306225/shared-genetic-architecture/OPERA/GWAS/CAD.ma --bfile /gpfs/share/home/1610306225/shared-genetic-architecture/1000G_v5a/binary/1000g503eur --estimate-pi --out OUTPUT-1/CAD

# Run OPERA for stage 2 analysis and heterogeneity analysis
opera --besd-flist /gpfs/share/home/1610306225/shared-genetic-architecture/OPERA/qtl-list --extract-gwas-loci myloci --chr 7 --gwas-summary mygwas.ma --bfile mydata --prior-pi-file myopera.pi --prior-var-file myopera.var --out myopera_chr7


# PIPE
# load packages
library(PIPE)
library(tidyverse)
library(patchwork)

# specify built-in data placeholder
placeholder <- "http://www.comptransmed.pro/bigdata_pipe"

# create the output folder 'pipe_showcase'
outdir <- 'pipe_showcase'
if(!dir.exists(outdir)) dir.create(outdir)

data.file <- file.path(placeholder, "GWAS.cross.SNP.txt")
input_data <- read_delim(data.file, delim='\t') %>% filter(pubmedid==26974007) %>% distinct(snp_id_current,pvalue)

input_data <- data.frame(matrix(NA,21,2))
input_data[,1] <- c('rs17035646',
'rs880315',
'rs4151702',
'rs4135240',
'rs2107595',
'rs2891168',
'rs7859727',
'rs635634',
'rs2519093',
'rs532436',
'rs360153',
'rs3184504',
'rs1535791',
'rs9534454',
'rs4942561',
'rs4932373',
'rs35346340',
'rs76774446',
'rs73015024',
'rs55997232',
'rs112374545')
input_data[,2] <- 5e-10
colnames(input_data) <- c('snp_id_current','pvalue')

# LD.customised
LD.customised <- oRDS('GWAS_LD.EUR', placeholder=placeholder)

# All common SNPs (primarily sourced from dbSNP)
options(timeout=120) # please increase the timeout if failed to download
GR.SNP <- oRDS("dbSNP_Common", placeholder=placeholder)

# All known genes (primarily sourced from UCSC)
GR.Gene <- oRDS("UCSC_knownGene", placeholder=placeholder)

# Functional interaction network (primarily sourced from STRING)
# restricted to high-quality interactions ("experiments" and "databases")
network.customised <- oRDS("STRING_fin", placeholder=placeholder)

# include.RGB
include.RGB <- c("PCHiC_PMID27863249_Activated_total_CD4_T_cells","PCHiC_PMID27863249_Macrophages_M0","PCHiC_PMID27863249_Macrophages_M1","PCHiC_PMID27863249_Macrophages_M2","PCHiC_PMID27863249_Megakaryocytes","PCHiC_PMID27863249_Monocytes","PCHiC_PMID27863249_Naive_B_cells","PCHiC_PMID27863249_Naive_CD4_T_cells","PCHiC_PMID27863249_Naive_CD8_T_cells","PCHiC_PMID27863249_Neutrophils","PCHiC_PMID27863249_Nonactivated_total_CD4_T_cells","PCHiC_PMID27863249_Total_B_cells","PCHiC_PMID27863249_Total_CD4_T_cells","PCHiC_PMID27863249_Total_CD8_T_cells","PCHiC_PMID31501517_DorsolateralPrefrontalCortex","PCHiC_PMID31501517_H1NeuralProgenitorCell","PCHiC_PMID31501517_Hippocampus","PCHiC_PMID31501517_LCL", "PCHiC_PMID31367015_astrocytes","PCHiC_PMID31367015_excitatory","PCHiC_PMID31367015_hippocampal","PCHiC_PMID31367015_motor", "PCHiC_PMID25938943_CD34","PCHiC_PMID25938943_GM12878")

# include.QTL
include.QTL <- c("eQTL_eQTLGen", "eQTL_eQTLCatalogue_Alasoo_2018_macrophage_IFNg","eQTL_eQTLCatalogue_Alasoo_2018_macrophage_IFNgSalmonella","eQTL_eQTLCatalogue_Alasoo_2018_macrophage_Salmonella","eQTL_eQTLCatalogue_Alasoo_2018_macrophage_naive","eQTL_eQTLCatalogue_BLUEPRINT_Tcell","eQTL_eQTLCatalogue_BLUEPRINT_monocyte","eQTL_eQTLCatalogue_BLUEPRINT_neutrophil","eQTL_eQTLCatalogue_BrainSeq_brain","eQTL_eQTLCatalogue_CEDAR_microarray_Bcell_CD19","eQTL_eQTLCatalogue_CEDAR_microarray_Tcell_CD4","eQTL_eQTLCatalogue_CEDAR_microarray_Tcell_CD8","eQTL_eQTLCatalogue_CEDAR_microarray_monocyte_CD14","eQTL_eQTLCatalogue_CEDAR_microarray_neutrophil_CD15","eQTL_eQTLCatalogue_CEDAR_microarray_platelet","eQTL_eQTLCatalogue_Fairfax_2012_microarray_Bcell_CD19","eQTL_eQTLCatalogue_Fairfax_2014_microarray_monocyte_IFN24","eQTL_eQTLCatalogue_Fairfax_2014_microarray_monocyte_LPS2","eQTL_eQTLCatalogue_Fairfax_2014_microarray_monocyte_LPS24","eQTL_eQTLCatalogue_Fairfax_2014_microarray_monocyte_naive","eQTL_eQTLCatalogue_GENCORD_LCL","eQTL_eQTLCatalogue_GENCORD_Tcell","eQTL_eQTLCatalogue_GEUVADIS_LCL","eQTL_eQTLCatalogue_GTEx_V8_Brain_Amygdala","eQTL_eQTLCatalogue_GTEx_V8_Brain_Anterior_cingulate_cortex_BA24","eQTL_eQTLCatalogue_GTEx_V8_Brain_Caudate_basal_ganglia","eQTL_eQTLCatalogue_GTEx_V8_Brain_Cerebellar_Hemisphere","eQTL_eQTLCatalogue_GTEx_V8_Brain_Cerebellum","eQTL_eQTLCatalogue_GTEx_V8_Brain_Cortex","eQTL_eQTLCatalogue_GTEx_V8_Brain_Frontal_Cortex_BA9","eQTL_eQTLCatalogue_GTEx_V8_Brain_Hippocampus","eQTL_eQTLCatalogue_GTEx_V8_Brain_Hypothalamus","eQTL_eQTLCatalogue_GTEx_V8_Brain_Nucleus_accumbens_basal_ganglia","eQTL_eQTLCatalogue_GTEx_V8_Brain_Putamen_basal_ganglia","eQTL_eQTLCatalogue_GTEx_V8_Brain_Spinal_cord_cervical_c1","eQTL_eQTLCatalogue_GTEx_V8_Brain_Substantia_nigra","eQTL_eQTLCatalogue_GTEx_V8_Nerve_Tibial","eQTL_eQTLCatalogue_GTEx_V8_Pituitary","eQTL_eQTLCatalogue_GTEx_V8_Whole_Blood","eQTL_eQTLCatalogue_Kasela_2017_microarray_Tcell_CD4","eQTL_eQTLCatalogue_Kasela_2017_microarray_Tcell_CD8","eQTL_eQTLCatalogue_Lepik_2017_blood","eQTL_eQTLCatalogue_Naranbhai_2015_microarray_neutrophil_CD16","eQTL_eQTLCatalogue_Nedelec_2016_macrophage_Listeria","eQTL_eQTLCatalogue_Nedelec_2016_macrophage_Salmonella","eQTL_eQTLCatalogue_Nedelec_2016_macrophage_naive","eQTL_eQTLCatalogue_Quach_2016_monocyte_IAV","eQTL_eQTLCatalogue_Quach_2016_monocyte_LPS","eQTL_eQTLCatalogue_Quach_2016_monocyte_Pam3CSK4","eQTL_eQTLCatalogue_Quach_2016_monocyte_R848","eQTL_eQTLCatalogue_Quach_2016_monocyte_naive","eQTL_eQTLCatalogue_ROSMAP_brain_naive","eQTL_eQTLCatalogue_Schmiedel_2018_Bcell_naive","eQTL_eQTLCatalogue_Schmiedel_2018_CD4_Tcell_antiCD3CD28","eQTL_eQTLCatalogue_Schmiedel_2018_CD4_Tcell_naive","eQTL_eQTLCatalogue_Schmiedel_2018_CD8_Tcell_antiCD3CD28","eQTL_eQTLCatalogue_Schmiedel_2018_CD8_Tcell_naive","eQTL_eQTLCatalogue_Schmiedel_2018_NKcell_naive","eQTL_eQTLCatalogue_Schmiedel_2018_Tfh_memory","eQTL_eQTLCatalogue_Schmiedel_2018_Th117_memory","eQTL_eQTLCatalogue_Schmiedel_2018_Th17_memory","eQTL_eQTLCatalogue_Schmiedel_2018_Th1_memory","eQTL_eQTLCatalogue_Schmiedel_2018_Th2_memory","eQTL_eQTLCatalogue_Schmiedel_2018_Treg_memory","eQTL_eQTLCatalogue_Schmiedel_2018_Treg_naive","eQTL_eQTLCatalogue_Schmiedel_2018_monocyte_CD16_naive","eQTL_eQTLCatalogue_Schmiedel_2018_monocyte_naive","eQTL_eQTLCatalogue_Schwartzentruber_2018_sensory_neuron")

# prepare predictors
ls_pNode_genomic <- oPierSNPsAdv(data=input_data, LD.customised=LD.customised, significance.threshold=5e-8, GR.SNP=GR.SNP, GR.Gene=GR.Gene, distance.max=20000, decay.kernel="constant", include.QTL=include.QTL, include.RGB=include.RGB, network.customised=network.customised, placeholder=placeholder, verbose.details=F)

# do prioritisation（the initial round）
ls_pNode <- Filter(Negate(is.null), ls_pNode_genomic)
dTarget0 <- oPierMatrix(ls_pNode, displayBy="pvalue", aggregateBy="fishers", rangeMax=10, GR.Gene=GR.Gene, placeholder=placeholder)

# GSP (phase 2 or above drug targets for inflammatory disorders)
ig.EF.ChEMBL <- oRDS('ig.EF.ChEMBL', placeholder=placeholder)
codes_name <- c("psoriasis","ankylosing spondylitis","Crohn's disease","ulcerative colitis","sclerosing cholangitis")
codes_name <- c("cardiovascular disease","coronary artery disease","myocardial infarction",
                "heart failure","atrial fibrillation","Ischemic stroke","stroke",
                "hypertension","genetic hypertension",                
                "Disorder of lipid metabolism","hyperlipidemia",
                "Genetic obesity","obesity",
                "diabetes mellitus",
                "internal carotid artery stenosis")
ind <- match(codes_name, V(ig.EF.ChEMBL)$name)
GSP <- base::Reduce(union, V(ig.EF.ChEMBL)$GSP[ind])

# GSN
druggable_targets <- do.call(rbind, V(ig.EF.ChEMBL)$phase) %>% pull(target) %>% unique()
g1 <- oRDS('org.Hs.string_medium', placeholder=placeholder)
g2 <- oRDS('org.Hs.PCommons_UN', placeholder=placeholder)
g <- oCombineNet(list(g1,g2), combineBy='union', attrBy="intersect")
sGS <- oGSsimulator(GSP=GSP, population=druggable_targets, network.customised=g, neighbor.order=1, verbose=F)
GSN <- sGS$GSN
message(sprintf("GSP (%d), GSN (%d)", length(GSP), length(GSN)), appendLF=T)

# predictor matrix
df_predictor <- oPierMatrix(ls_pNode, displayBy="score", aggregateBy="none", GR.Gene=GR.Gene, placeholder=placeholder)

# informative predictors
sClass <- oClassifyRF(df_predictor, GSP=GSP, GSN=GSN, nfold=3, nrepeat=10, mtry=NULL, ntree=1000, cv.aggregateBy="fishers")
df_imp <- sClass$importance %>% as_tibble(rownames='predictor')
df_perf <- sClass$performance %>% as_tibble(rownames='predictor')
imp_perf <- df_imp %>% inner_join(df_perf, by='predictor') %>% mutate(group=str_replace_all(predictor,'_.*','')) %>% arrange(group,-MeanDecreaseAccuracy)
xintercept <- imp_perf %>% filter(group=='nGene') %>% pull(MeanDecreaseAccuracy)
yintercept <- imp_perf %>% filter(group=='nGene') %>% pull(MeanDecreaseGini)
zintercept <- imp_perf %>% filter(group=='nGene') %>% pull(auroc)
imp_perf2 <- imp_perf %>% filter(MeanDecreaseAccuracy>=0, MeanDecreaseGini>=0, direction=='+') %>% mutate(x=MeanDecreaseAccuracy>=xintercept, y=MeanDecreaseGini>=yintercept, z=auroc>=zintercept) %>% arrange(-MeanDecreaseAccuracy)

# predictors selected based the importance no worse than nGene
vec_predictor_selected <- imp_perf2 %>% filter(x|y|z) %>% pull(predictor)
vec_predictor_selected

# do prioritisation（based on informative predictors）
dTarget <- oPierMatrix(ls_pNode[vec_predictor_selected], displayBy="pvalue", aggregateBy="fishers", rangeMax=10, GR.Gene=GR.Gene, placeholder=placeholder, verbose=F)

# write into a file 'CRM_IND_priority.txt.gz' under the folder 'pipe_showcase'
dTarget$priority %>% select(name,rating,rank,seed,nGene,eGene,cGene,description) %>% write_delim(file.path(outdir,'CRM_IND_priority.txt.gz'), delim='\t')

# save into a file 'dTarget_IND.CRM.RData' under the folder 'pipe_showcase'
save(list='dTarget', file=file.path(outdir,'dTarget_IND.CRM.RData'))

# df_priority
df_priority <- dTarget$priority %>% column_to_rownames("name") %>% top_frac(1,rating) %>% transmute(priority=rating, rank=rank)

write.table(V(ig.EF.ChEMBL)$name, "name.csv", quote = F, row.names = F)

# prepare SET
codes_name <- c("cardiovascular disease","coronary artery disease","myocardial infarction",
"heart failure","atrial fibrillation","Ischemic stroke","stroke",
"hypertension","genetic hypertension",                
"Disorder of lipid metabolism","hyperlipidemia",
"Genetic obesity","obesity",
"diabetes mellitus",
"internal carotid artery stenosis")
ind <- match(codes_name, V(ig.EF.ChEMBL)$name)
codes <- c('CVD','CAD','MI','HF','AF','AS','AIS','HYP','gHYP','LM','HLIP',
           'GOS','OS','DM','IMT')

ls_gsp <- V(ig.EF.ChEMBL)$GSP[ind]
names(ls_gsp) <- codes

info <- ls_gsp %>% enframe() %>% transmute(id=name, member=value)

info <- info %>% mutate(name=codes_name, namespace='IND', distance=NA) %>% mutate(n=map_dbl(member,length)) %>% select(name,id,namespace,distance,member,n)

set <- list(info=info, stamp=as.Date(Sys.time()))
class(set) <- "SET"

# do leading prioritisation analysis (LPA)
set.seed(1020)

str(df_priority)
df_cvd <- data.frame(matrix(NA,12,2))
df_cvd[,1] <- 10
df_cvd[,2] <- 1
rownames(df_cvd) <- c('CASZ1', 'CDKN1A', 'TWIST1', 'CDKN2B', 'ABO', 
'SWAP70', 'SH2B3', 'LRCH1', 'FES', 'GOSR2', 'RPRML', 'LDLR')
colnames(df_cvd) <- c('priority','rank')

df <- df_priority[100:120, ]

eLPA <- oPierGSEA(df_priority, customised.genesets=set, size.range=c(10,5000), nperm=500000)

df_summary <- eLPA$df_summary %>% mutate(flag=ifelse(adjp<0.05,'Y','N')) %>% mutate(frac=nLead/nAnno) %>% filter(frac<=1) %>% arrange(tolower(setID))
df_leading <- eLPA$leading[df_summary %>% pull(setID)] %>% enframe() %>% transmute(setID=name,members=value) %>% mutate(rank_members=map_chr(members, function(x) str_c(names(x),' (',x,')', collapse=', ')))
df_summary_leading <- df_summary %>% inner_join(df_leading, by='setID')

# save into a file 'CRM_IND_lpa.txt' under the folder 'pipe_showcase'
df_summary_leading %>% select(-members) %>% write_delim(file.path(outdir,'CRM_IND_lpa.txt'), delim='\t')


pict <- df_summary_leading %>% mutate(pict=log10(nes * frac / adjp)) %>% select(setID, nes, adjp, frac, nAnno, nLead, pict) %>% arrange(tolower(setID))
gp_pict <- oPolarBar(pict %>% transmute(name=setID,value=pict), colormap='sci_locuszoom', size.name=9, size.value=4)

# saved into a file 'CRM_IND_pict.pdf' under the folder 'pipe_showcase'
ggsave(file.path(outdir,'CRM_IND_pict.pdf'), gp_pict, width=6.5, height=6.5)

leading.query <- df_priority %>% as_tibble(rownames='target') %>% top_frac(1,priority) %>% pull(target)
gp_leading <- oGSEAdotplot(eLPA, top='ulcerative colitis', x.scale="normal", leading=T, leading.edge.only=T, leading.query.only=T, leading.query=leading.query, colormap="jet.top", zlim=c(0,10), peak.color='black', leading.size=2, clab='Credit\nscore', leading.label.direction="left") + theme(plot.title=element_text(hjust=0.5,size=10), plot.subtitle=element_text(hjust=0.5,size=8))

# saved into a file 'CRM_IND_leading_uc.pdf' under the folder 'pipe_showcase'
ggsave(file.path(outdir,'CRM_IND_leading_uc.pdf'), gp_leading, width=6, height=8)

# Member genes for each disease defined as clinical proof-of-concept targets recovered at the leading prioritisation
data <- df_summary_leading %>% transmute(name=setID, members=rank_members)

# remove (rank)
data0 <- data %>% separate_rows(members, sep=", ") %>% mutate(members=str_replace_all(members, ' \\(.*', '')) %>% nest(data=-name) %>% mutate(members=map_chr(data, function(x) str_c(unlist(x), collapse=', '))) %>% select(-data)

# bigraph
big <- data %>% separate_rows(members, sep=", ") %>% mutate(flag=1) %>% select(members,name,flag) %>% pivot_wider(names_from=members, values_from=flag, values_fill=list(flag=0)) %>% column_to_rownames("name") %>% oBicreate() %>% oBigraph()
big <- big %>% oLayout("graphlayouts.layout_with_stress")
igraph::degree(big) -> V(big)$degree
if(1){
  V(big)$rating <- str_replace_all(V(big)$name, '.*\\(|\\).*', '') %>% as.numeric()
  ifelse(V(big)$type=='xnode', V(big)$name, "") -> V(big)$label
  ifelse(V(big)$type=='ynode' & V(big)$rating<=nrow(df_priority) * 0.02 | V(big)$type=='xnode', V(big)$name, "") -> V(big)$label
}
gp_big <- big %>% oGGnetwork(node.label="label", node.label.size=2, node.label.force=0.02, node.xcoord="xcoord", node.ycoord="ycoord", colormap='steelblue-steelblue', node.shape="type", node.size="degree", node.size.range=c(0.5,6),edge.color="grey80", edge.color.alpha=0.3, edge.arrow.gap=0, edge.curve=0.05)
## saved into a file 'CRM_IND_big.pdf' under the folder 'pipe_showcase'
ggsave(file.path(outdir,'CRM_IND_big.pdf'), gp_big, width=6, height=6)

# identify the minimum spanning tree (mst)
tig <- data %>% oTIG(method='empirical', empirical.cutoff=0.8)
tig <- tig %>% oLayout("gplot.layout.fruchtermanreingold")
gp_mst <- oGGnetwork(tig, node.label='name', node.label.size=2.5, node.label.color='steelblue', node.label.force=1, node.xcoord='xcoord', node.ycoord='ycoord', colormap="steelblue-steelblue", edge.size='weight', node.shape=18, node.size='n_members', node.size.title="Num of\ngenes", slim=c(0,25), node.size.range=c(2,5), edge.color='gray80', edge.color.alpha=0.3, edge.curve=0,edge.arrow.gap=0.01)
## saved into a file 'CRM_IND_mst.pdf' under the folder 'pipe_showcase'
ggsave(file.path(outdir,'CRM_IND_mst.pdf'), gp_mst, width=5, height=5)

library(ComplexHeatmap)

# df_priority_leading
data <- df_summary_leading %>% transmute(name=setID,members) %>% mutate(target=map(members, function(x) names(x))) %>% mutate(rating=map(members, function(x) x)) %>% select(name,target,rating) %>% unnest(c(target,rating))
df_x <- data %>% pivot_wider(names_from=name, values_from=rating)
df_y <- df_priority %>% as_tibble(rownames='target')
df_priority_leading <- df_y %>% inner_join(df_x, by='target') %>% arrange(-priority)

# write into a file 'CRM_IND_priority_leading.txt' under the folder 'pipe_showcase'
df_priority_leading %>% write_delim(file.path(outdir,'CRM_IND_priority_leading.txt'), delim='\t')

# ht_main
D <- df_priority_leading %>% mutate(rname=str_c(target,' (',rank,')')) %>% select(rname, priority) %>% column_to_rownames('rname')
col_fun <- circlize::colorRamp2(seq(0,6,length=64), oColormap('spectral')(64))
#vec_group <- df_priority_leading %>% pull(rank)
#vec_group <- as.numeric(cut(df_priority_leading$rank, breaks=seq(0,10000,100)))
vec_group <- as.numeric(cut(df_priority_leading$rank, breaks=ceiling(nrow(df_priority) * seq(0,0.5,0.05))))
ht_main <- Heatmap(D, name="Credit\nscore", col=col_fun, border=T, 
                   row_split=vec_group, cluster_rows=F, show_row_dend=T, row_title=NULL, show_row_names=T, row_names_gp=gpar(fontsize=5), row_names_rot=0, row_names_side="left", clustering_distance_rows="euclidean", clustering_method_rows="average",
                   cluster_columns=F, column_order=1:ncol(D), show_column_names=F, column_names_side="top", column_names_gp=gpar(fontsize=8), column_names_rot=0
)

# row_anno: ls_ha_mark
colormap <- "sci_locuszoom"
ncolor <- ncol(df_priority_leading) - 3
ls_ha_mark <- lapply(1:ncolor, function(j){
  node.color <- oColormap(colormap)(ncolor)[j]
  vec_at <- df_priority_leading %>% pull(j+3)
  vec_at[!is.na(vec_at)] <- 1
  vec_at[is.na(vec_at)] <- 0
  at <- which(vec_at==1)
  ha_mark_at <- rowAnnotation(
    D=anno_simple(vec_at, col=c("0"="grey95", "1"=node.color), width=unit(1.5,"mm")),
    mark=anno_mark(at=at, labels=df_priority_leading$target[at], labels_gp=gpar(fontsize=5,col="black"), lines_gp=gpar(col=node.color))
  )
  # D1, D2, D3, ...
  ha_mark_at@anno_list$D@label <- str_c('D',j)
  ha_mark_at
})
names(ls_ha_mark) <- colnames(df_priority_leading)[4:(3+ncolor)]
names(ls_ha_mark)
#ht <- ht_main + ls_ha_mark[[1]] + ls_ha_mark[[2]] + ls_ha_mark[[3]] + ls_ha_mark[[4]] + ls_ha_mark[[5]]
#draw(ht, heatmap_legend_side="left", annotation_legend_side="left")

# row_anno: ha_mark_approved
codes_name <- c("ankylosing spondylitis","Crohn's disease","psoriasis","sclerosing cholangitis","ulcerative colitis")
ChEMBL <- oRDS("ChEMBL_v33", placeholder=placeholder)
ChEMBL_subset <- ChEMBL %>% filter(efo_term %in% codes_name)
### df_target_label
df_target_label <- ChEMBL_subset %>% filter(phase==4,target_number<=5) %>% select(Symbol,pref_name_drug,efo_term) %>% nest(data=pref_name_drug) %>% mutate(drugs=map_chr(data, function(x) {
  y <- x %>% pull(pref_name_drug)
  if(length(y)>1){
    y <- c(y %>% head(1), '...')
    str_c(y,collapse=" | ")
  }else{
    str_c(y %>% head(1),collapse=" | ")
  }
})) %>% mutate(label=str_c(efo_term,' (',drugs,')')) %>% select(Symbol,label) %>% group_by(Symbol) %>% summarise(label=str_c(label,collapse='\n    ')) %>% ungroup() %>% mutate(label=str_c(Symbol,':\n    ',label))
### df_label
df_label <- df_priority_leading %>% select(target) %>% left_join(df_target_label, by=c('target'='Symbol'))
vec_at <- df_label$label
at <- which(!is.na(vec_at))
ha_mark_approved <- rowAnnotation(
  A=anno_simple(0+!is.na(vec_at), col=c("0"="grey90", "1"="black"), width=unit(3,"mm")),
  mark=anno_mark(at=at, labels=df_label$label[at], labels_gp=gpar(fontsize=5,col='black'), lines_gp=gpar(col='black'))
)
#ht <- ht_main + ha_mark_approved
#draw(ht, heatmap_legend_side="left", annotation_legend_side="left")

# do visualisation
## saved into a file 'CRM_IND_heatmap_therapeutics.pdf' under the folder 'pipe_showcase'
ht <- ht_main + ls_ha_mark[[1]] + ls_ha_mark[[2]] + ls_ha_mark[[3]] + ls_ha_mark[[4]] + ls_ha_mark[[5]] + ha_mark_approved
pdf(file.path(outdir, 'CRM_IND_heatmap_therapeutics.pdf'), width=6.7, height=9.5)
draw(ht, heatmap_legend_side="left", annotation_legend_side="left")
dev.off()
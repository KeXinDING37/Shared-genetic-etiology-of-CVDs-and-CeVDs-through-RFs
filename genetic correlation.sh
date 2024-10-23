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
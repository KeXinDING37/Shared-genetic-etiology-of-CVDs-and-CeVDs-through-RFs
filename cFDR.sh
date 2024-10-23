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
        threshold=1e-2, 
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

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

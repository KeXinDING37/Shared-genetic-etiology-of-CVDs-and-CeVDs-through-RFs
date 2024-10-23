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
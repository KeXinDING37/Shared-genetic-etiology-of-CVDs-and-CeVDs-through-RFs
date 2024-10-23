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


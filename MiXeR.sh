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

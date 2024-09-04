
#cd /home/sunhongyu/liran/clusIBD/challengingDNA

##筛选家系中的无关个体
/opt/software/vcftools/bin/vcftools \
--vcf /home/sunhongyu/liran/HBfam_ASA/HBfam_ASA_renamed.vcf \
--maf 0.05 \
--min-alleles 2 \
--max-alleles 2 \
--hwe 0.000001 \
--max-missing 0.5 \
--thin 2000 \
--keep /home/sunhongyu/liran/HBfam_ASA/HBfam_unrelated.txt \
--kept-sites \
--out HBfam_ASA_renamed_unrelated

plink --bfile /home/sunhongyu/liran/HBfam_ASA/degratedDNA/data \
--chr 1-22 \
--recode vcf-iid \
-out HBfam_degratedDNA_thin

##筛选多态性的vcf
/opt/software/vcftools/bin/vcftools \
--vcf HBfam_degratedDNA_thin.vcf \
--positions HBfam_ASA_renamed_unrelated.kept.sites \
--recode \
--out HBfam_degratedDNA
##重命名
bcftools reheader -s /home/sunhongyu/liran/HBfam_ASA/degratedDNA/old_newname.txt HBfam_degratedDNA.recode.vcf -o HBfam_degratedDNA.vcf
plink --vcf HBfam_degratedDNA.vcf  --double-id  --out HBfam_degratedDNA


##分析
num_cores=10
##RUN IBIS
ibis -t $num_cores  -b HBfam_degratedDNA -ibd2 -f ./results/ibis_HBdd
##RUN TRUFFLE
~/software/truffle/truffle --vcf HBfam_degratedDNA.vcf --cpu $num_cores --segments --out ./results/truffle_HBdd
##RUN IBDseq
for j in {1..22}
do
java -Xmx2000m -jar ~/software/IBDseq/ibdseq.r1206.jar \
gt=HBfam_degratedDNA.vcf \
out=./results/ibdseq_HBdd_$j \
nthreads=$num_cores \
chrom=$j

done
find ./results -type f -printf '%p\n' | grep ibdseq_HBdd_ | grep ibd$ |xargs cat >./results/ibdseq_HBdd.ibd

##Run clusIBD
#clusIBD -f HBfam_degratedDNA -o ./results/clusIBD_HBdd
python /home/sunhongyu/liuzhentang/pycharm/clusIBD/IBDpy_test.py -f HBfam_degratedDNA -o ./results/clusIBD_HBdd



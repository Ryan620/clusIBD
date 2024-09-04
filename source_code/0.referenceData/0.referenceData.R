
##cd /home/sunhongyu/zangyu/Pedsim_simulation/simulated_pedigrees/Error/LiRan/ER659k/1kg.Chinese.ER659K
##从千人基因组数据，筛选中国汉族和随机选择的659K SNP数据

for i in {1..22}; do /opt/software/vcftools/bin/vcftools \
--gzvcf ~/public_data/SNP_STR_phased/1kg.snp.str.chr${i}.vcf.gz  \
--keep ~/liran/simulations/1000G_Chinese_samples.txt \
--positions  ~/zangyu/Pedsim_simulation/simulated_pedigrees/SNP_filter/659K/snp_list3.txt \
--maf 0.05 \
--min-alleles 2 \
--max-alleles 2 \
--recode \
--out /home/sunhongyu/zangyu/Pedsim_simulation/simulated_pedigrees/Error/LiRan/ER659k/1kg.Chinese.ER659K/1kg.chinese.er659k.chr${i};done

for i in {1..22}
do bgzip 1kg.chinese.er659k.chr${i}.recode.vcf
bcftools index 1kg.chinese.er659k.chr${i}.recode.vcf.gz
done

#ls -l | grep ".gz$" | bcftools concat -Ov -o 1kg.chinese.er659K.vcf 
bcftools concat  1kg.chinese.er659k.chr1.recode.vcf.gz \
1kg.chinese.er659k.chr2.recode.vcf.gz \
1kg.chinese.er659k.chr3.recode.vcf.gz \
1kg.chinese.er659k.chr4.recode.vcf.gz \
1kg.chinese.er659k.chr5.recode.vcf.gz \
1kg.chinese.er659k.chr6.recode.vcf.gz \
1kg.chinese.er659k.chr7.recode.vcf.gz \
1kg.chinese.er659k.chr8.recode.vcf.gz \
1kg.chinese.er659k.chr9.recode.vcf.gz \
1kg.chinese.er659k.chr10.recode.vcf.gz \
1kg.chinese.er659k.chr11.recode.vcf.gz \
1kg.chinese.er659k.chr12.recode.vcf.gz \
1kg.chinese.er659k.chr13.recode.vcf.gz \
1kg.chinese.er659k.chr14.recode.vcf.gz \
1kg.chinese.er659k.chr15.recode.vcf.gz \
1kg.chinese.er659k.chr16.recode.vcf.gz \
1kg.chinese.er659k.chr17.recode.vcf.gz \
1kg.chinese.er659k.chr18.recode.vcf.gz \
1kg.chinese.er659k.chr19.recode.vcf.gz \
1kg.chinese.er659k.chr20.recode.vcf.gz \
1kg.chinese.er659k.chr21.recode.vcf.gz \
1kg.chinese.er659k.chr22.recode.vcf.gz \
-Ov -o 1kg.chinese.er659K.vcf

##提出FORMAT中GT之外的信息
bcftools annotate --remove ^FORMAT/GT 1kg.chinese.er659K.vcf -Ov -o 1kg.chinese.er659K.clean.vcf

/opt/software/vcftools/bin/vcftools \
--vcf 1kg.chinese.er659K.clean.vcf \
--thin 2000 \
--recode \
--out 1kg.chinese.400K.clean

##compress the vcf file
bgzip 1kg.chinese.400K.clean.recode.vcf 

#plink --vcf 1kg.chinese.er659K.clean.vcf  --indep-pairwise 2000 100 0.5  --out unlinked_sites  ##18W SNPs


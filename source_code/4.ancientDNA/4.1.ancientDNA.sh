
#awk -F' ' '{print $2 "\t" $4}' /home/sunhongyu/zangyu/Sample/AADNA_ENA/1240K/1240K/v54.1.p1_1240K_public.snp >1240k.snps.txt

#cd /home/sunhongyu/liran/clusIBD/ancientDNA

for i in {1..22}
do
echo "chr ${i} start"
 /opt/software/vcftools/bin/vcftools \
--gzvcf ~/public_data/SNP_STR_phased/1kg.snp.str.chr${i}.vcf.gz  \
--keep samples_from_England.txt \
--positions 1240k.snps.txt \
--maf 0.05 \
--thin 2000 \
--min-alleles 2 \
--max-alleles 2 \
--chr 1 \
--chr 2 \
--chr 3 \
--chr 4 \
--chr 5 \
--chr 6 \
--chr 7 \
--chr 8 \
--chr 9 \
--chr 10 \
--chr 11 \
--chr 12 \
--chr 13 \
--chr 14 \
--chr 15 \
--chr 16 \
--chr 17 \
--chr 18 \
--chr 19 \
--chr 20 \
--chr 21 \
--chr 22 \
--recode \
--out 1kg.England.chr${i}
bcftools view 1kg.England.chr${i}.recode.vcf -Oz -o 1kg.England.chr${i}.recode.vcf.gz
bcftools index 1kg.England.chr${i}.recode.vcf.gz
done
##merge VCF files
ls | grep 1kg | grep gz$ >1kg.England.list
bcftools concat -f 1kg.England.list -Ov -o 1kg.England.vcf

####output the SNP sites
/opt/software/vcftools/bin/vcftools \
--vcf 1kg.England.vcf \
--kept-sites \
--out 1kg.England

##skip the first line
awk -F'\t' 'NR > 1 {print "chr" $1 "\t" $2}' 1kg.England.kept.sites >1kg.England.chr.kept.sites
##delete unused files
for i in {1..22}
do
rm 1kg.England.chr${i}.recode.vcf
rm 1kg.England.chr${i}.log 
rm 1kg.England.chr${i}.recode.vcf.gz
rm 1kg.England.chr${i}.recode.vcf.gz.csi
done

#######################################################################
#######################################################################
#古DNA
cd /home/sunhongyu/zangyu/Sample/AADNA_ENA/ENA_PRJEB46958_All_Samples/Rawnochr
###change the chromosome names from 1-22 to chr1-chr22 for each bam file
ls | grep bam$ |  cat | while read id
do 
sh ./add_chr_to_bam.sh $id /home/sunhongyu/zangyu/Sample/AADNA_ENA/ENA_PRJEB46958_All_Samples/Rawchr/chr_$id
echo $id
done

##call SNP
cd /home/sunhongyu/zangyu/Sample/AADNA_ENA/ENA_PRJEB46958_All_Samples/Rawchr

ls | grep bam$ > ancient.sample.list

samtools mpileup \
-uf ~/public_data/ref_genomes/hg19.fa \
-l /home/sunhongyu/liran/clusIBD/ancientDNA/1kg.England.chr.kept.sites \
-b ancient.sample.list | bcftools call -Ou -mv > /home/sunhongyu/liran/clusIBD/ancientDNA/var.raw.bcf

#cd /home/sunhongyu/liran/clusIBD/ancientDNA
#
bcftools view var.raw.bcf | bcftools filter -e 'QUAL<20 || DP<3' > ancientDNA_DP3.vcf

plink --vcf ancientDNA_DP3.vcf --chr 1-22 --geno 0.5  --mind 0.5  --double-id --recode vcf-iid --out ancientDNA_filter
bcftools query -l ancientDNA_filter.vcf > oldnames.txt

####rename the ancient DNA sanmples so that they are consistent with those in  Fowler et al.  Nature 2022
#R
#oldnames <- read.table("oldnames.txt",header = F,sep = "\t",stringsAsFactors = F)
#newnames <- read.table("sample_names.txt",header = T,sep = "\t",stringsAsFactors = F)
#oldnames$new_name <- as.character( oldnames$V1)
#for (i in 1:nrow(newnames)){
#  x <- as.character( newnames$names2[i])
#  y <- grep(x= as.character(oldnames$V1),pattern=x,value = F)
  
#  if (length(y)>0)oldnames$new_name[y] <- as.character(newnames$name1[i])
#}

#write.table(x=oldnames,file = "old2newname.txt",col.names = F,row.names = F,quote = F,sep = "\t")
#quit("no")


bcftools reheader -s old2newname.txt ancientDNA_filter.vcf -o ancientDNA_filter_renamed.vcf

#plink --vcf ancientDNA_filter_renamed.vcf  --double-id --make-bed --out ancientDNA_filter_renamed

####导出保留的位点
/opt/software/vcftools/bin/vcftools \
--vcf ancientDNA_filter_renamed.vcf \
--kept-sites \
--out ancientDNA

/opt/software/vcftools/bin/vcftools \
--vcf 1kg.England.vcf \
--positions ancientDNA.kept.sites \
--recode \
--out England_present_thin
bcftools annotate --remove ^FORMAT/GT England_present_thin.recode.vcf -Ov -o England_present.vcf
###去除phase
plink --vcf England_present.vcf  --double-id --recode vcf-iid --out England_present2


###We merge vcf files with R, because we fail to do so using bcftools
R
library(readr)
x <- read_delim(file = "England_present2.vcf",delim = "\t",skip = 27)
y <- read_delim(file = "ancientDNA_filter_renamed.vcf",delim = "\t",skip = 27)
new_x <- merge(x=y,y=x[,c(1,2,10:ncol(x))],by=c("#CHROM","POS"),all.x=T,sort=F)
new_x$"#CHROM" <- paste0("chr",new_x$"#CHROM")
which(new_x$POS[2:nrow(new_x)]<new_x$POS[1:(nrow(new_x)-1)])
write.table(new_x,file = "ancientDNA_merged.vcf",col.names = T,row.names = F,quote = F,sep = "\t")
quit("no")


plink --vcf ancientDNA_merged.vcf  --double-id --make-bed --out ancientDNA_merged


##分析IBD
num_cores=10
ibis -t $num_cores  -b ancientDNA_merged -ibd2 -f ./results/ibis_ancientDNA
~/software/truffle/truffle --vcf ancientDNA_merged.vcf --cpu $num_cores --segments --out ./results/truffle_ancientDNA

for j in {1..22}
do
java -Xmx2000m -jar ~/software/IBDseq/ibdseq.r1206.jar \
gt=ancientDNA_merged.vcf \
out=./results/ibdseq_ancientDNA_$j \
nthreads=$num_cores \
chrom=$j

done
find ./results -type f -printf '%p\n' | grep ibdseq_ancientDNA_ | grep ibd$ |xargs cat >./results/ibdseq_ancientDNA.ibd

#python /home/sunhongyu/liran/clusIBD/bin/clusIBD2.py -f ancientDNA_merged -n 300  -c 10  -o ./results/clusIBD_ancientDNA
python /home/sunhongyu/liuzhentang/pycharm/clusIBD/IBDpy_test.py -f ancientDNA_merged  -c 10 -o ./results/clusIBD_ancientDNA_default
python /home/sunhongyu/liuzhentang/pycharm/clusIBD/IBDpy_test.py -f ancientDNA_merged  -c 10 -n 300 -o ./results/clusIBD_ancientDNA_n300
python /home/sunhongyu/liuzhentang/pycharm/clusIBD/IBDpy_test.py -f ancientDNA_merged  -c 10 -n 500 -o ./results/clusIBD_ancientDNA_n500

python /home/sunhongyu/liuzhentang/pycharm/clusIBD/IBDpy_test.py -f ancientDNA_merged  -s 100 -c 10 -n 500 -o ./results/clusIBD_ancientDNA_n500





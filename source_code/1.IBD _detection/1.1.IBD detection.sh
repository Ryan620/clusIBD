

#uncompress the reference data
gunzip 1kg.chinese.400K.clean.recode.vcf.gz

##generate IBD segments of 10-30 Mb and with genotype error of 0,0.005,0.01,0.05,0.1,and 0.2
Rscript generate_fake_IBD.R


num_cores=20 
IBDlength=30
##transform related files to plink format files
for i in {0,0.005,0.01,0.05,0.1,0.2}
do
  echo "start"
  plink --vcf fakeibd_error_${i}_IBD_${IBDlength}cM.vcf  --double-id  --out fakeibd_error_${i}_IBD_${IBDlength}cM
done
########################
##detect IBD segments

##Run IBIS
 for i in {0,0.005,0.01,0.05,0.1,0.2}
 do
  ibis -t $num_cores  -b fakeibd_error_${i}_IBD_${IBDlength}cM -f ./results/ibis_error_${i}_${IBDlength}cM
 done

##truffle
for i in {0,0.005,0.01,0.05,0.1,0.2}
do
~/software/truffle/truffle --vcf fakeibd_error_${i}_IBD_${IBDlength}cM.vcf --cpu $num_cores --segments --out ./results/truffle_error_${i}_${IBDlength}cM
done

#clusIBD
for i in {0,0.005,0.01,0.05,0.1,0.2}
do
python /home/sunhongyu/liuzhentang/pycharm/clusIBD/IBDpy_test.py -f fakeibd_error_${i}_IBD_${IBDlength}cM  -c $num_cores  -o ./results/clusIBD_error_${i}_${IBDlength}cM
done

##ibdseq
for i in {0,0.005,0.01,0.05,0.1,0.2}
do
  for j in {1..22}
    do
     java -Xmx2000m -jar ~/software/IBDseq/ibdseq.r1206.jar \
     gt=fakeibd_error_${i}_IBD_${IBDlength}cM.vcf \
     out=./results/ibdseq__error_${i}_${IBDlength}cM_$j \
     nthreads=$num_cores \
     chrom=$j

   done
  find ./results -type f -printf '%p\n' | grep ibdseq__error_${i}_${IBDlength}cM_ | grep ibd$ |xargs cat >./results/ibdseq_error_${i}_${IBDlength}cM.ibd
  
done


#clusIBD
#for i in {0,0.005,0.01,0.05,0.1,0.2}
#do
#python /home/sunhongyu/liran/clusIBD/bin/clusIBD.py -f fakeibd_error_${i}_IBD_${IBDlength}cM -n 200  -c $num_cores  -o ./results/clusIBD_error_${i}_${IBDlength}cM
#done


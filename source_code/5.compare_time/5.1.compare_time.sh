




#!/bin/bash
#cd /home/sunhongyu/liran/clusIBD/comparetime

##
for i in {100,200,300,400,500}
do
vcftools --vcf /home/sunhongyu/liran/clusIBD/kinship/error_0.vcf --max-indv $i --recode --out random_${i}
plink --vcf random_${i}.recode.vcf  --double-id  --out random_${i} 
done


#cd /home/sunhongyu/zangyu/Pedsim_simulation/simulated_pedigrees/Error/LiRan/ER659k/vcf_100_500/results
num_cores=5

##Run IBIS
echo "compare start!"  | tee  compare_time.log
for i in {100,200,300,400,500}
do
    startTime=`date +"%Y-%m-%d %H:%M:%S"`
    
    ibis -b random_${i}  -ibd2 -t $num_cores -f ibis_$i
    
    endTime=`date +"%Y-%m-%d %H:%M:%S"`
    st=`date -d  "$startTime" +%s`
    et=`date -d  "$endTime" +%s`
    sumTime=$(($et-$st))
    echo "IBIS, sampleN==${i},Total time is : $sumTime second."  | tee -a  compare_time.log
done

##Run TRUFFLE
for i in {100,200,300,400,500}
do
    startTime=`date +"%Y-%m-%d %H:%M:%S"`
    
    ~/software/truffle/truffle --vcf random_${i}.recode.vcf --cpu $num_cores --segments --out truffle_$i
    
    endTime=`date +"%Y-%m-%d %H:%M:%S"`
    st=`date -d  "$startTime" +%s`
    et=`date -d  "$endTime" +%s`
    sumTime=$(($et-$st))
    echo "TRUFFLE, sampleN==${i},Total time is : $sumTime second."  | tee -a  compare_time.log
done




##Run clusIBD
for i in {100,200,300,400,500}
do
startTime=`date +"%Y-%m-%d %H:%M:%S"`
     #clusIBD -f /home/sunhongyu/zangyu/Pedsim_simulation/simulated_pedigrees/Error/LiRan/ER659k/vcf_100_500/random_${i}  -n 200  -c $num_cores  -o clusIBD__${i}
    /home/sunhongyu/liuzhentang/pycharm/clusIBD/dist/clusIBD/clusIBD -f random_${i} -c $num_cores  -o clusIBD__${i}
    endTime=`date +"%Y-%m-%d %H:%M:%S"`
    st=`date -d  "$startTime" +%s`
    et=`date -d  "$endTime" +%s`
    sumTime=$(($et-$st))
    echo "clusIBD, sampleN==${i},Total time is : $sumTime second."  | tee -a  compare_time.log
done


##IBDseq
for i in {100,200,300,400,500}
  do
  startTime=`date +"%Y-%m-%d %H:%M:%S"`

      for j in {1..22}
      do
      java -Xmx5000m -jar ~/software/IBDseq/ibdseq.r1206.jar gt=random_${i}.recode.vcf out=IBDseq_${i}_${j} nthreads=$num_cores chrom=$j
      done
  endTime=`date +"%Y-%m-%d %H:%M:%S"`
  st=`date -d  "$startTime" +%s`
  et=`date -d  "$endTime" +%s`
  sumTime=$(($et-$st))
  echo "IBDseq, sampleN==${i},Total time is : $sumTime second."  | tee -a  compare_time.log

done

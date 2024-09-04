
#cd /home/sunhongyu/liran/clusIBD/kinship
num_cores=15 
IBDlength=30



for i in {0,0.005,0.01,0.05,0.1,0.2}
do
  echo "start detect IBD segments with four tools." | tee   compare_time.log
   plink --vcf error_${i}.vcf  --double-id  --out error_${i}
done

startTime=`date +"%Y-%m-%d %H:%M:%S"`

##Run IBIS
 for i in {0,0.005,0.01,0.05,0.1,0.2}
 do
  ibis -t $num_cores -ibd2   -b error_${i} -f ./results/ibis_error_${i}
 done
endTime=`date +"%Y-%m-%d %H:%M:%S"`
st=`date -d  "$startTime" +%s`
et=`date -d  "$endTime" +%s`
sumTime=$(($et-$st))
echo "IBIS, Total time is : $sumTime second." | tee -a  compare_time.log

##truffle
startTime=`date +"%Y-%m-%d %H:%M:%S"`
for i in {0,0.005,0.01,0.05,0.1,0.2}
do
~/software/truffle/truffle --vcf error_${i}.vcf  --cpu $num_cores --segments --out ./results/truffle_error_${i}
done
endTime=`date +"%Y-%m-%d %H:%M:%S"`
st=`date -d  "$startTime" +%s`
et=`date -d  "$endTime" +%s`
sumTime=$(($et-$st))
echo "TRUFFLE, Total time is : $sumTime second."  | tee -a  compare_time.log

#clusIBD
startTime=`date +"%Y-%m-%d %H:%M:%S"`

for i in {0,0.005,0.01,0.05,0.1,0.2}
do
python /home/sunhongyu/liuzhentang/pycharm/clusIBD/IBDpy_test.py -f error_${i} -p target_samplepairs.txt  -c $num_cores  -o ./results/clusIBD_error_${i}
done
endTime=`date +"%Y-%m-%d %H:%M:%S"`
st=`date -d  "$startTime" +%s`
et=`date -d  "$endTime" +%s`
sumTime=$(($et-$st))
echo "clusIBD, Total time is : $sumTime second."  | tee -a  compare_time.log

##ibdseq
startTime=`date +"%Y-%m-%d %H:%M:%S"`
for i in {0,0.005,0.01,0.05,0.1,0.2}
do
  for j in {1..22}
    do
     java -Xmx2000m -jar ~/software/IBDseq/ibdseq.r1206.jar \
     gt=error_${i}.vcf \
     out=./results/ibdseq_error_${i}_$j \
     nthreads=$num_cores \
     chrom=$j

   done
  find ./results -type f -printf '%p\n' | grep ibdseq_error_${i}_ | grep ibd$ |xargs cat >./results/ibdseq_error_${i}.ibd
  
done
endTime=`date +"%Y-%m-%d %H:%M:%S"`
st=`date -d  "$startTime" +%s`
et=`date -d  "$endTime" +%s`
sumTime=$(($et-$st))
echo "IBDseq, Total time is : $sumTime second."  | tee -a  compare_time.log


#clusIBD
#for i in {0,0.005,0.01,0.05,0.1,0.2}
#do
#python /home/sunhongyu/liran/clusIBD/bin/clusIBD.py -f fakeibd_error_${i}_IBD_${IBDlength}cM -n 200  -c $num_cores  -o ./results/clusIBD_error_${i}_${IBDlength}cM
#done


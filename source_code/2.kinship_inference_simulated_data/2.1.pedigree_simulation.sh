   
##
gunzip ../0.referenceData/1kg.chinese.400K.clean.recode.vcf.gz


##simulation of family data
##error0#0
for i in {7799..7808}
 do /home/sunhongyu/software/ped-sim-master/ped-sim \
 -d full_half_2C.def \
 -m sexavg.simmap \
 -i 1kg.chinese.400K.clean.recode.vcf \
 --err_rate 0 \
 --seed ${i} \
 -o error0_sim${i} \
 --pois \
 --keep_phase
done

##error2#0.005
for i in {7799..7808}
 do /home/sunhongyu/software/ped-sim-master/ped-sim \
 -d full_half_2C.def \
 -m sexavg.simmap \
 -i 1kg.chinese.400K.clean.recode.vcf \
 --err_rate 5.0e-03 \
 --seed ${i} \
 -o error2_sim${i} \
 --pois \
 --keep_phase
done


##error3#0.01
for i in {7799..7808}
 do /home/sunhongyu/software/ped-sim-master/ped-sim \
 -d full_half_2C.def \
 -m sexavg.simmap \
 -i 1kg.chinese.400K.clean.recode.vcf \
 --err_rate 1.0e-02 \
 --seed ${i} \
 -o error3_sim${i} \
 --pois \
 --keep_phase
done

#ฒน#error4#0.05
for i in {7799..7808}
 do /home/sunhongyu/software/ped-sim-master/ped-sim \
 -d full_half_2C.def \
 -m sexavg.simmap \
 -i 1kg.chinese.400K.clean.recode.vcf \
 --err_rate 5.0e-02 \
 --seed ${i} \
 -o error4_sim${i} \
 --pois \
 --keep_phase
done

#ฒน#error5#0.10
for i in {7799..7808}
 do /home/sunhongyu/software/ped-sim-master/ped-sim \
 -d full_half_2C.def \
 -m sexavg.simmap \
 -i 1kg.chinese.400K.clean.recode.vcf \
 --err_rate 1.0e-01 \
 --seed ${i} \
 -o error5_sim${i} \
 --pois \
 --keep_phase
done

#ฒน#error5#0.20
for i in {7799..7808}
 do /home/sunhongyu/software/ped-sim-master/ped-sim \
 -d full_half_2C.def \
 -m sexavg.simmap \
 -i 1kg.chinese.400K.clean.recode.vcf \
 --err_rate 2.0e-01 \
 --seed ${i} \
 -o error6_sim${i} \
 --pois 
 #--keep_phase
done



##Genotype error rate:1.0e-02;Missingness rate:1.0e-03
#
find . -name "*sim7799*" -exec sed -i "s/full_half_2C/sim7799/g" {} \;
find . -name "*sim7800*" -exec sed -i "s/full_half_2C/sim7800/g" {} \;
find . -name "*sim7801*" -exec sed -i "s/full_half_2C/sim7801/g" {} \;
find . -name "*sim7802*" -exec sed -i "s/full_half_2C/sim7802/g" {} \;
find . -name "*sim7803*" -exec sed -i "s/full_half_2C/sim7803/g" {} \;
find . -name "*sim7804*" -exec sed -i "s/full_half_2C/sim7804/g" {} \;
find . -name "*sim7805*" -exec sed -i "s/full_half_2C/sim7805/g" {} \;
find . -name "*sim7806*" -exec sed -i "s/full_half_2C/sim7806/g" {} \;
find . -name "*sim7807*" -exec sed -i "s/full_half_2C/sim7807/g" {} \;
find . -name "*sim7808*" -exec sed -i "s/full_half_2C/sim7808/g" {} \;

for i in {7799..7808}; do bgzip error0_sim${i}.vcf; done
for i in {7799..7808}; do bcftools index error0_sim${i}.vcf.gz --threads 10 -f; done
bcftools merge *.vcf.gz -O v -o error0.vcf


 find . -name "*.vcf" -exec sed -i "s/g1-b1-i1/A/g" {} \;
 find . -name "*.vcf" -exec sed -i "s/g1-b1-s1/B/g" {} \;
 find . -name "*.vcf" -exec sed -i "s/g2-b1-s1/C/g" {} \;
 find . -name "*.vcf" -exec sed -i "s/g2-b1-i1/D/g" {} \;
 find . -name "*.vcf" -exec sed -i "s/g2-b2-s1/E/g" {} \;
 find . -name "*.vcf" -exec sed -i "s/g2-b2-i1/F/g" {} \;
 find . -name "*.vcf" -exec sed -i "s/g2-b2-s2/G/g" {} \;
 find . -name "*.vcf" -exec sed -i "s/g3-b1-s1/H/g" {} \;
 find . -name "*.vcf" -exec sed -i "s/g3-b1-i1/I/g" {} \;
 find . -name "*.vcf" -exec sed -i "s/g3-b2-s1/J/g" {} \;
 find . -name "*.vcf" -exec sed -i "s/g3-b2-i1/K/g" {} \;
 find . -name "*.vcf" -exec sed -i "s/g3-b3-i1/L/g" {} \; 
 find . -name "*.vcf" -exec sed -i "s/g3-b3-s1/M/g" {} \;
 find . -name "*.vcf" -exec sed -i "s/g4-b1-s1/N/g" {} \;
 find . -name "*.vcf" -exec sed -i "s/g4-b1-i1/O/g" {} \;
 find . -name "*.vcf" -exec sed -i "s/g4-b2-s1/P/g" {} \;
 find . -name "*.vcf" -exec sed -i "s/g4-b2-i1/Q/g" {} \;
 find . -name "*.vcf" -exec sed -i "s/g4-b3-i1/R/g" {} \;
 find . -name "*.vcf" -exec sed -i "s/g4-b3-s1/S/g" {} \;
 find . -name "*.vcf" -exec sed -i "s/g5-b1-i1/T/g" {} \;
 find . -name "*.vcf" -exec sed -i "s/g5-b2-i1/U/g" {} \;
 find . -name "*.vcf" -exec sed -i "s/g5-b3-i1/V/g" {} \;




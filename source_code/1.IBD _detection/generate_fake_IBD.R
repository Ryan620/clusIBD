library(readr)

generate_fakeIBD <- function(ref_vcf,error_rate=0,start_id=10,IBDlength=10,pair_size=100,IBD_type=1,phased=F,target_chr=NULL){
    
    generate_error_geno <- function(g,error_rate=0.01){
      idx <- which(runif(length(g)) < error_rate)
      g0 <- g[idx]
      gx <- g0
      idx1 <- which(g0 =="0|0")
      idx2 <- which(g0 =="0|1")
      idx3 <- which(g0 =="1|0")
      idx4 <- which(g0 =="1|1")
      gx[idx1] <- "0|1"
      gx[idx1[runif(length(idx1))<0.5]] <- "1|0"
      gx[idx2] <- "0|0"
      gx[idx2[runif(length(idx2))<0.5]] <- "1|1"
      gx[idx3] <- "0|0"
      gx[idx3[runif(length(idx3))<0.5]] <- "1|1"
      gx[idx4] <- "0|1"
      gx[idx4[runif(length(idx4))<0.5]] <- "1|0"
      g[idx] <- gx
      return(g)
    }
    
    geno <- ref_vcf[,start_id:ncol(ref_vcf)]
    samplenames <- colnames(geno)
    chr_pos <- ref_vcf[,1:2]
    colnames(chr_pos) <- c("chr","pos")
    # ##
    
    ##填充遗传距离  
    chrs <- unique(chr_pos$chr)
    if (is.null(target_chr)) target_chr <- chrs
   
    chr_pos$pos_cM <- chr_pos$pos/1000000 
    
    #两两随机组合
    print("generate reference haplotype")
    sample_num <- ncol(geno)
    ref_haplotypes <- matrix("",nrow = nrow(geno),ncol = 2*sample_num)
    
    for (i in 1:ncol(geno)){
      g <- as.character(geno[[i]])
      #g_mat <- matrix(unlist(strsplit(g,split = "|",fixed=T)),ncol = 2,byrow = T)
      g_mat <- cbind(substr(g,1,1),substr(g,3,3))
      ref_haplotypes[,c(2*i-1,2*i)] <- g_mat
    }
    print("generation of reference haplotype finished")
    
    sample_pairs <- t(combn(1:sample_num,2))
    if (nrow(sample_pairs) > pair_size) sample_pairs <- sample_pairs[sample(1:nrow(sample_pairs),pair_size),]
    #
    print("generate IBD")
    real_IBD <- matrix(NA,nrow =nrow(sample_pairs) *22,ncol = 6 )
    colnames(real_IBD) <- c("sample1","sample2","chr","start_real","end_real","length_real")
    real_IBD <- as.data.frame(real_IBD)
    out_genotypes <- matrix("",nrow = nrow(geno),ncol = nrow(sample_pairs)*2)
    colnames(out_genotypes) <- paste("pairs",rep(1: nrow(sample_pairs),each=2 ),1:2,sep = "_")
    
    for (i in 1:nrow(sample_pairs)){
      id1 <- sample_pairs[i,1]
      id2 <- sample_pairs[i,2]
      hap1 <- ref_haplotypes[,c(2*id1-1,2*id1)]
      hap2 <- ref_haplotypes[,c(2*id2-1,2*id2)]
      #real_IBD[(i*22-21):(i*22),1] <- samplenames[id1]
      #real_IBD[(i*22-21):(i*22),2] <- samplenames[id2]
      real_IBD[(i*22-21):(i*22),1] <- paste("pairs",i,"1",sep = "_")
      real_IBD[(i*22-21):(i*22),2] <- paste("pairs",i,"2",sep = "_")
      ##每条染色体随机产生IBD
      #new_geno <- c()
      for (j in chrs){
        if (!(j %in% target_chr)) next
        idx_chr <- chr_pos$chr==j
        h1 <- hap1[idx_chr,]    #个体1的两个单倍型
        h2 <- hap2[idx_chr,]    #个体2的两个单倍型
        sub_pos <- chr_pos[chr_pos$chr==j,]
        idx2 <- NA
        cycle_num=1
        while (is.na(idx2)){
          idx1 <- sample(1:nrow(h1),1)
          idx2 <- which(sub_pos$pos_cM > sub_pos$pos_cM[idx1]+IBDlength)[1]
          
          cycle_num <- cycle_num +1 
          if (cycle_num > 100000) break
        }
         ##如果IBD结束位置始终大于染色体长度，则取整个染色体
        if (is.na(idx2)){idx1=1;idx2=nrow(h1)}
        col_index <- sample(1:2,2,replace = T) 
        #将个体的替换目标区域，替换为个体2
        h1[idx1:idx2,col_index[1]] <- h2[idx1:idx2,col_index[2]]
        #new_geno <- rbind(new_geno,cbind(h1,h2))
        hap1[idx_chr,] <- h1
        hap2[idx_chr,] <- h2
        
        real_IBD[i*22-22+j,"chr"] <- j
        real_IBD[i*22-22+j,"start_real"] <- sub_pos$pos_cM[idx1]
        real_IBD[i*22-22+j,"end_real"] <- sub_pos$pos_cM[idx2]
        real_IBD[i*22-22+j,"length_real"] <- sub_pos$pos_cM[idx2]-sub_pos$pos_cM[idx1]
      }
      
      g1 <- paste(hap1[,1],hap1[,2],sep = "|")
      g2 <- paste(hap2[,1],hap2[,2],sep = "|")
      g1 <- generate_error_geno(g = g1,error_rate = error_rate )
      g2 <- generate_error_geno(g = g2,error_rate = error_rate )
      
      if (!phased){
        g1 <- gsub(g1,pattern = "|",replacement = "/",fixed = T)
        g2 <- gsub(g2,pattern = "|",replacement = "/",fixed = T)
        g1[g1 == "1/0"] <- "0/1"
        g2[g2 == "1/0"] <- "0/1"
      }
      out_genotypes[,2*i-1] <- g1
      out_genotypes[,2*i] <- g2
    }
    out_genotypes <- cbind(ref_vcf[,1:(start_id-1)],out_genotypes)
    out_genotypes <- as.data.frame(out_genotypes)
    print("generaton of IBD finished")
    return(c(IBD=list(real_IBD),genotypes=list(out_genotypes)))
  }

print("load reference panel .......................")
system.time( { ref_vcf <- read_delim(file = paste0("/home/sunhongyu/zangyu/Pedsim_simulation/simulated_pedigrees/Error/LiRan/ER659k/1kg.Chinese.ER659K/1kg.chinese.400K.clean.recode.vcf"),delim = "\t",skip = 35)})
print("load reference panel succeed.")

##不同错误率下，不同IBD片段长度的准确性
  t1 <- Sys.time()
  error_rates <- c(0.1,0.2)
  names(error_rates) <- paste0("error",0:(length(error_rates)-1))
  IBD_length <- seq(5,40,by=5)

  for (i in error_rates) {
    for (j in IBD_length){
      set.seed(111)
      simu_IBDs <- generate_fakeIBD(ref_vcf = ref_vcf,
                                error_rate = i,
                                start_id = 10,
                                pair_size = 100,
                                IBDlength = j,
                                target_chr = 1:5)
      simu_pairs <- matrix(colnames(simu_IBDs[[2]])[4:ncol(simu_IBDs[[2]])],ncol=2,byrow = TRUE)
      simu_geno <- simu_IBDs[[2]]
      write.table(x=simu_geno,file = paste0("fakeibd_error_",i,"_IBD_",j,"cM.vcf"),sep = "\t",row.names = F,col.names = T,quote = F)
      simu_IBDs <- as.data.frame(simu_IBDs[[1]])
      write.table(x=simu_IBDs,file = paste0("realIBD_errors_",i,"_IBD_",j,"cM.txt"),sep = "\t",row.names = F,col.names = T,quote = F)
      print(paste0("error= ",i, ", IBD=",j," Mb done"))
    }
  }



t2 <- Sys.time()

print(paste0("Total time is ",t2-t1))
#f/opa/software/vcftools/bin/vcftools --vcf /home/sunhongyu/zangyu/Pedsim_simulation/simulated_pedigrees/Error/LiRan/ER659k/1kg.Chinese.ER659K/1kg.chinese.er659K.clean.vcf --thin 2000 --recode --out 1kg.chinese.400K.clean




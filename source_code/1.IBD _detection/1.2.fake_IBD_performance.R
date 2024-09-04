
library(readr)
#library(ggplot2)
#library(reshape2)

setwd("/home/sunhongyu/liran/clusIBD/fake_IBD")

#real_ibd_data <- read_delim(file = "~/liran/clusIBD/fake_IBD/realIBD_errors_0.05_IBD_20cM.txt",delim="\t")
#sample_pairs <- data.frame(sample1 = paste("pairs",1:100,1,sep = "_"),sample2 = paste("pairs",1:100,2,sep = "_"))
#real_ibd_data <- merge(x=sample_pairs,real_ibd_data,all.x  = T)
#real_ibd_data <- real_ibd_data[real_ibd_data$IBD_length == 20,c("sample1","sample2","error_rate","chr","IBD_length","pos_start_real","pos_end_real")]


estimate_IBDacc <- function(x){
  
  overlap_start <- apply(cbind(x$start_detect,x$start_real), MARGIN = 1, FUN = function(x)max(x)) 
  overlap_end <- apply(cbind(x$end_detect,x$end_real), MARGIN = 1, FUN = function(x)min(x)) 
  overlap_length <- overlap_end-overlap_start
  overlap_length[overlap_length < 0] <- 0
  overlap_length[is.na(overlap_length) ] <- 0
  ##检测出的IBD片段超过0.5为真实的IBD即为正确
  acc <- length(which(overlap_length/x$length_detect > 0.5))/length(which(!is.na(x$length_detect)))
  overlap_prop<- overlap_length/x$length_detect
  #检测出的IBD index
  idx <- which(!is.na(x$length_detect))
  len.acc <- sum(overlap_prop[idx],na.rm = T)/length(idx)
  
  #真实的IBD index
  recall_rate <- length(which(overlap_length/x$length_real > 0.5))/length(which(!is.na(x$length_real)))
  overlap_prop<- overlap_length/x$length_real
  idx <- which(!is.na(x$length_real))
  power_rate <- sum(overlap_prop[idx],na.rm = T)/length(idx)
  
  return(c("accuracy"=acc,
           "len.accuracy"=len.acc,
           "recall" =recall_rate,
           "power"= power_rate
  ))
}


error_rates <- c(0,0.005,0.01,0.05,0.1,0.2)
IBDlengths <- c(10,15,20,25,30)
min_length <- 5
max_chr <- 6
fake_IBD_all <- c()

##读取IBIS的结果
for (IBDlength in IBDlengths){
ibis_acc <- c()
for (i in error_rates){
  
  ibis_fake_ibd <- read_delim(file = paste0("/home/sunhongyu/liran/clusIBD/fake_IBD/results/ibis_error_",i,"_",IBDlength,"cM.seg"),delim="\t",col_names = F)
  sample1 <- sapply(strsplit(as.character(ibis_fake_ibd$X1),split = ":",fixed = T),FUN = function(x)x[1])
  sample2 <- sapply(strsplit(as.character(ibis_fake_ibd$X2),split = ":",fixed = T),FUN = function(x)x[1])
  ibis_fake_ibd <- data.frame(sample1=sample1,
                              sample2=sample2,
                              chr=ibis_fake_ibd$X3,
                              start_detect = ibis_fake_ibd$X4,
                              end_detect = ibis_fake_ibd$X5)
  ibis_fake_ibd$start_detect <- ibis_fake_ibd$start_detect/1000000
  ibis_fake_ibd$end_detect <- ibis_fake_ibd$end_detect/1000000
  ibis_fake_ibd$length_detect <- ibis_fake_ibd$end_detect - ibis_fake_ibd$start_detect
  ##1-5号染色体大于7cM的片段
  ibis_fake_ibd <- ibis_fake_ibd[ibis_fake_ibd$length_detect > min_length & ibis_fake_ibd$chr < max_chr, ]
  #两个样本来自同一对
  s1 <- substr(as.character(ibis_fake_ibd$sample1),1,nchar(as.character(ibis_fake_ibd$sample1))-2)
  s2 <- substr(as.character(ibis_fake_ibd$sample2),1,nchar(as.character(ibis_fake_ibd$sample2))-2)
  ibis_fake_ibd <- ibis_fake_ibd[s1==s2,]
  
  read_IBD <- read.table(file = paste0("realIBD_errors_",i,"_IBD_",IBDlength,"cM.txt"),sep = "\t",header = T)
  read_IBD <- read_IBD[!is.na(read_IBD$chr),]
  
  x <- merge(read_IBD,y=ibis_fake_ibd,all.x=T,all.y = T)
  ibis_acc <- rbind(ibis_acc,estimate_IBDacc(x))
}


##读取truffle的结果
error_rates <- c(0,0.005,0.01,0.05,0.1,0.2)
truffle_acc <- c()
for (i in error_rates){
  filename <- paste0("/home/sunhongyu/liran/clusIBD/fake_IBD/results/truffle_error_",i,"_",IBDlength,"cM.segments")
  truffle_fake_ibd <- read_table(filename)
  truffle_fake_ibd$sample1 <-truffle_fake_ibd$ID1
  truffle_fake_ibd$sample2 <- truffle_fake_ibd$ID2
  truffle_fake_ibd$chr <- truffle_fake_ibd$CHROM
  truffle_fake_ibd$start_detect <- truffle_fake_ibd$POS
  truffle_fake_ibd$end_detect <- truffle_fake_ibd$start_detect + truffle_fake_ibd$LENGTH
  truffle_fake_ibd$length_detect <- truffle_fake_ibd$LENGTH
  ##1-5号染色体大于7cM的片段
  truffle_fake_ibd <- truffle_fake_ibd[truffle_fake_ibd$length_detect > min_length & truffle_fake_ibd$chr < max_chr, ]
  #两个样本来自同一对
  s1 <- substr(as.character(truffle_fake_ibd$sample1),1,nchar(as.character(truffle_fake_ibd$sample1))-2)
  s2 <- substr(as.character(truffle_fake_ibd$sample2),1,nchar(as.character(truffle_fake_ibd$sample2))-2)
  truffle_fake_ibd <- truffle_fake_ibd[s1==s2,]
  
  read_IBD <- read.table(file = paste0("realIBD_errors_",i,"_IBD_",IBDlength,"cM.txt"),sep = "\t",header = T)
  read_IBD <- read_IBD[!is.na(read_IBD$chr),]
  
  x <- merge(read_IBD,y=truffle_fake_ibd,all.x=T,all.y = T)
  truffle_acc <- rbind(truffle_acc,estimate_IBDacc(x))
  
}

##IBDseq
error_rates <- c(0,0.005,0.01,0.05,0.1,0.2)
ibdseq_acc <- c()
for (i in error_rates){
  filename <- paste0("/home/sunhongyu/liran/clusIBD/fake_IBD/results/ibdseq_error_",i,"_",IBDlength,"cM.ibd")

  ibdseq_fake_ibd <- read_delim(filename,col_names = F,delim = "\t")
  ibdseq_fake_ibd$sample1 <-ibdseq_fake_ibd$X1
  ibdseq_fake_ibd$sample2 <- ibdseq_fake_ibd$X3
  ibdseq_fake_ibd$chr <- ibdseq_fake_ibd$X5
  ibdseq_fake_ibd$start_detect <- ibdseq_fake_ibd$X6/1000000
  ibdseq_fake_ibd$end_detect <- ibdseq_fake_ibd$X7/1000000
  ibdseq_fake_ibd$length_detect <- ibdseq_fake_ibd$end_detect - ibdseq_fake_ibd$start_detect
  ##1-5号染色体大于7cM的片段
  ibdseq_fake_ibd <- ibdseq_fake_ibd[ibdseq_fake_ibd$length_detect > min_length & ibdseq_fake_ibd$chr < max_chr, ]
  #两个样本来自同一对
  s1 <- substr(as.character(ibdseq_fake_ibd$sample1),1,nchar(as.character(ibdseq_fake_ibd$sample1))-2)
  s2 <- substr(as.character(ibdseq_fake_ibd$sample2),1,nchar(as.character(ibdseq_fake_ibd$sample2))-2)
  ibdseq_fake_ibd <- ibdseq_fake_ibd[s1==s2,]
  
  read_IBD <- read.table(file = paste0("realIBD_errors_",i,"_IBD_",IBDlength,"cM.txt"),sep = "\t",header = T)
  read_IBD <- read_IBD[!is.na(read_IBD$chr),]
  
  x <- merge(read_IBD,y=ibdseq_fake_ibd,all.x=T,all.y = T)
  ibdseq_acc <- rbind(ibdseq_acc,estimate_IBDacc(x))
}



clusIBD_acc <- c()
for (i in error_rates){
  
  clusIBD_fake_ibd <- read_delim(file = paste0("/home/sunhongyu/liran/clusIBD/fake_IBD/results/clusIBD_error_",i,"_",IBDlength,"cM.IBD.details"),delim="\t",col_names = F)
  clusIBD_fake_ibd <- data.frame(sample1=clusIBD_fake_ibd$X1,
                              sample2=clusIBD_fake_ibd$X2,
                              chr=clusIBD_fake_ibd$X3,
                              start_detect = clusIBD_fake_ibd$X4,
                              end_detect = clusIBD_fake_ibd$X5)
 
  clusIBD_fake_ibd$length_detect <- clusIBD_fake_ibd$end_detect - clusIBD_fake_ibd$start_detect
  ##1-5号染色体大于7cM的片段
  clusIBD_fake_ibd <- clusIBD_fake_ibd[clusIBD_fake_ibd$length_detect > min_length & clusIBD_fake_ibd$chr < max_chr, ]
  #两个样本来自同一对
  s1 <- substr(as.character(clusIBD_fake_ibd$sample1),1,nchar(as.character(clusIBD_fake_ibd$sample1))-2)
  s2 <- substr(as.character(clusIBD_fake_ibd$sample2),1,nchar(as.character(clusIBD_fake_ibd$sample2))-2)
  clusIBD_fake_ibd <- clusIBD_fake_ibd[s1==s2,]
  
  read_IBD <- read.table(file = paste0("realIBD_errors_",i,"_IBD_",IBDlength,"cM.txt"),sep = "\t",header = T)
  read_IBD <- read_IBD[!is.na(read_IBD$chr),]
  
  x <- merge(read_IBD,y=clusIBD_fake_ibd,all.x=T,all.y = T)
  clusIBD_acc <- rbind(clusIBD_acc,estimate_IBDacc(x))
}

mydata <- rbind(
  cbind(data.frame(method=rep("IBIS",length(error_rates)),
             rates=error_rates),ibis_acc),
  
  cbind(data.frame(method=rep("TRUFFLE",length(error_rates)),
                   rates=error_rates),truffle_acc),
  cbind(data.frame(method=rep("IBDseq",length(error_rates)),
                   rates=error_rates),ibdseq_acc),
  
  cbind(data.frame(method=rep("clusIBD",length(error_rates)),
                   rates=error_rates),clusIBD_acc)
)
mydata <- cbind(IBDlength = rep(IBDlength,nrow(mydata)),mydata)
fake_IBD_all <- rbind(fake_IBD_all,mydata)
}


write.table(fake_IBD_all,file = "fake_IBD_powers_for_different_methods.txt",sep = "\t",quote = F,row.names = F)



##本地电脑
mycol2 = c( "#E64B35B2", "#4DBBD5B2" ,"#00A087B2", "#3C5488B2", "#F39B7FB2",
            "#8491B4B2", "#91D1C2B2", "#DC0000B2", "#7E6148B2")
mycol1 <- c("#BD6263","#8EA325","#A9D179","#84CAC0","#F5AE6B","#BCB8D3","#4387B5")
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
               "#D55E00", "#CC79A7")
mycols <- c("#E64B35B2", "#4DBBD5B2" ,"#E69F00", "#3C5488B2", "#F39B7FB2","#999999")
mydata <- read.table(file = "A:/课题/kinship/errorIBD/data/fake_IBD/fake_IBD_powers_for_different_methods.txt",sep = "\t",header = T)
mydata[mydata$method=="clusIBD",]
x <- melt(mydata,id.vars = c("method","rates","IBDlength"),
          measure.vars = c("accuracy" ,"len.accuracy","recall"  ,"power" ),
          variable.name = "parameters",
          value.name = "values")

x$Error <- factor(paste0(x$rates * 100,"%"),levels = c("0%","0.5%","1%","5%","10%","20%"))

x$parameters <- factor(x$parameters,levels =c("recall"  ,"power","accuracy" ,"len.accuracy") )
x$Error2 <- as.numeric(x$Error)
x$values[ is.na(x$values) & x$parameters %in% c("recall"  ,"power") ] <- 0
x$method2 <- factor(x$method ,levels = c("clusIBD","IBIS","TRUFFLE","IBDseq"))

g_power <- ggplot(data = x,mapping = aes(x=IBDlength,y=values,color=Error))+
  scale_colour_manual(values = mycols)+
  geom_point(alpha=0.9)+
  geom_line(alpha=0.9)+
  facet_grid(method2~parameters)+
  ylim(c(0,1))+
  labs(x="IBD length (Mb)",y="",color="Error rates")

#g_power
#ggsave(filename =  paste0("A:/课题/kinship/errorIBD/results/Figure 1 powr",Sys.Date(),".jpeg"),g_power,height=5,width =6)


x$IBDlength <- factor(paste(x$IBDlength,"Mb"))

g_power_compare <- ggplot(x,mapping = aes(x=Error2,y=values,color=method))+
  geom_point(alpha=0.7,size=0.8)+
  scale_colour_manual(values = mycols)+
  geom_line()+
  lims(y=c(0,1))+
  facet_grid(parameters~IBDlength)+
  scale_x_continuous(breaks = 1:6,labels = c("0%","0.5%","1%","5%","10%","20%"))+
  labs(x="Error rates",y="",color="Methods")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))



g_power_compare
ggsave(filename =  paste0("A:/课题/kinship/errorIBD/results/Figure 2. comparison of the perforamce ",Sys.Date(),".jpeg"),g_power_compare,height=5,width =7)
ggsave(filename =  paste0("A:/课题/kinship/errorIBD/results/Figure 2. comparison of the perforamce.jpeg"),g_power_compare,height=5,width =7)







#####IBD片段的边界
error_rates <- c(0,0.005,0.01,0.05,0.1,0.2)
IBDlength=15
clusIBD_edge <- c()
for (IBDlength in c(10,15,20,25,30)){
  for (i in error_rates){
    
    clusIBD_fake_ibd <- read_delim(file = paste0("A:/课题/kinship/errorIBD/data/fake_IBD/clusIBD_error_",i,"_",IBDlength,"cM.IBD.details"),delim="\t",col_names = F)
    clusIBD_fake_ibd <- data.frame(sample1=clusIBD_fake_ibd$X1,
                                   sample2=clusIBD_fake_ibd$X2,
                                   chr=clusIBD_fake_ibd$X3,
                                   start_detect = clusIBD_fake_ibd$X4,
                                   end_detect = clusIBD_fake_ibd$X5)
    
    clusIBD_fake_ibd$length_detect <- clusIBD_fake_ibd$end_detect - clusIBD_fake_ibd$start_detect
    ##1-5号染色体大于5cM的片段
    clusIBD_fake_ibd <- clusIBD_fake_ibd[clusIBD_fake_ibd$length_detect > 5 & clusIBD_fake_ibd$chr < 6, ]
    #两个样本来自同一对
    s1 <- substr(as.character(clusIBD_fake_ibd$sample1),1,nchar(as.character(clusIBD_fake_ibd$sample1))-2)
    s2 <- substr(as.character(clusIBD_fake_ibd$sample2),1,nchar(as.character(clusIBD_fake_ibd$sample2))-2)
    clusIBD_fake_ibd <- clusIBD_fake_ibd[s1==s2,]
    
    read_IBD <- read.table(file = paste0("A:/课题/kinship/errorIBD/data/fake_IBD/realIBD_errors_",i,"_IBD_",IBDlength,"cM.txt"),sep = "\t",header = T)
    read_IBD <- read_IBD[!is.na(read_IBD$chr),]
    
    x <- merge(read_IBD,y=clusIBD_fake_ibd,all.x=T,all.y = T)
    x$rates <- rep(i,nrow(x))
    x$IBDlength <- rep(IBDlength,nrow(x))
    clusIBD_edge <- rbind(clusIBD_edge,x)
  }
  
}
clusIBD_edge$start_dif <- clusIBD_edge$start_detect-clusIBD_edge$start_real
clusIBD_edge$end_dif <- clusIBD_edge$end_detect-clusIBD_edge$end_real
clusIBD_edge$rates2 <- factor(clusIBD_edge$rates,levels =c(0,0.005,0.01,0.05,0.1,0.2),labels =  c("0%","0.5%","1%","5%","10%","20%"))
clusIBD_edge$IBDlength2 <- factor(clusIBD_edge$IBDlength,levels = unique(clusIBD_edge$IBDlength),labels = paste0(unique(clusIBD_edge$IBDlength),"Mb"))

mycol2 = c( "#E64B35B2", "#4DBBD5B2" ,"#00A087B2", "#3C5488B2", "#F39B7FB2",
            "#8491B4B2", "#91D1C2B2", "#DC0000B2", "#7E6148B2")
g_edge <- ggplot(clusIBD_edge)+
  geom_density(mapping = aes(x=start_dif,fill="Start"),alpha=0.5)+
  geom_density(mapping = aes(x=end_dif,fill="End"),alpha=0.5)+
  scale_fill_manual(values = mycol2[1:2])+
  #geom_abline(xintercept=0, color="red", linetype="dashed",size=0.5)+
  geom_vline(xintercept = 0,color="red", linetype="dashed",size=0.5)+
  facet_grid(rates2~IBDlength2)+
  labs(x="Estimated - Actual (Mb)",fill="")+
  theme_bw()+
  xlim(c(-10,10))
g_edge
ggsave(filename = paste0("A:/课题/kinship/errorIBD/results/Supplementary Figure 4 edges of clusIBD ",Sys.Date(),".jpeg"),
       g_edge,height=7,width = 8)

##the median difference between estiamted and real edges
aggregate(start_dif~rates+IBDlength,clusIBD_edge,FUN = median)
aggregate(end_dif~rates+IBDlength,clusIBD_edge,FUN = median)




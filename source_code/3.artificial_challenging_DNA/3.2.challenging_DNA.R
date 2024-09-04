
library(readr)
library(ggplot2)
library(reshape2)
library(ggpubr)


HBfam_ddASA <- "A:/课题/kinship/errorIBD/data/dilution-degradationASA/HBfam_degratedDNA.vcf"
HBfam_ddASA_data <-  read_delim(file = HBfam_ddASA,delim = "\t",skip = 27)

###估算分型错误率
estimate_genotype_acc <- function(g1,g2,allelesep="/"){
  idx1 <- g1 != paste0(".",allelesep,".")
  idx2 <- g2 != paste0(".",allelesep,".")
  idx <- idx1 & idx2
  g1 <- factor( g1[idx],levels = paste(c(0,0,1),c(0,1,1),sep=allelesep))
  g2 <- factor( g2[idx],levels = paste(c(0,0,1),c(0,1,1),sep=allelesep))
  dif_idx <- which(g1 != g2)
  x <- table(data.frame(g1=g1[dif_idx],g2=g2[dif_idx]))
  return(c(error_rate=list(length(dif_idx)/length(g1)),
           dropin=list(sum(x[1,2],x[3,2])/length(g1)),
           dropout=list(sum(x[2,1],x[2,3])/length(g1)),
           switch_error=list(sum(x[1,3],x[3,1])/length(g1)),
           effictiveSNPnum = list(length(g1)),
           details=list(x)))
}
fam_names <- c("A01","A03","A20","D01","D03","D10")
types <- c("10ng","1ng","0.5ng","0.1ng","1500bp","800bp","400bp","150bp")
acc_outs <- data.frame(ref=rep(paste0(fam_names,"_100ng"),each=length(types)),
                       samples=paste(rep(fam_names,each=length(types)),rep(types,6),sep = "_"),
                       types=rep(types,6))
acc_outs$types <-  factor(acc_outs$types,levels = unique(types)[order(as.numeric(substr(unique(types),1,nchar(unique(types))-2)))])
acc_outs$error_rates=rep(0,nrow(acc_outs))
acc_outs$dropin=rep(0,nrow(acc_outs))
acc_outs$dropout=rep(0,nrow(acc_outs))
acc_outs$switch_error=rep(0,nrow(acc_outs))

for (i in 1:nrow(acc_outs)) {
  x <- estimate_genotype_acc(g1=HBfam_ddASA_data[,acc_outs$ref[i]],g2=HBfam_ddASA_data[,acc_outs$samples[i]])
  acc_outs[i,"error_rates"] <- x[[1]]
  acc_outs[i,"dropin"] <- x[[2]]
  acc_outs[i,"dropout"] <- x[[3]]
  acc_outs[i,"switch_error"] <- x[[4]]
}
x <- melt(acc_outs,id.vars = c("ref","samples","types"),
          measure.vars = c("error_rates" ,"dropin","dropout"  ,"switch_error" ),
          variable.name = "Error_type",
          value.name = "Values")

aggregate(Values~types+Error_type,x,FUN = "mean")
aggregate(Values~types+Error_type,x,FUN = "sd")

x$Error_type <- factor(x$Error_type,
                       levels = c("dropin","dropout"  ,"switch_error","error_rates" ),
                       labels  = c("dropin","dropout"  ,"switch error","all" ))
mycol2 = c( "#E64B35B2", "#4DBBD5B2" ,"#00A087B2", "#3C5488B2", "#F39B7FB2",
            "#8491B4B2", "#91D1C2B2", "#DC0000B2", "#7E6148B2")
acc_data <- x
g_genotype_acc <- ggplot(data = acc_data,mapping = aes(x=Error_type,y=Values,color=Error_type))+
  geom_violin( scale = "width")+
  geom_point()+
  scale_color_manual(values =mycol2[1:4] )+
  #geom_boxplot(alpha=0.5)+
  
  facet_wrap(vars(types),nrow = 2)+
  labs(color="Error types",x="",y="Error rates")+
  theme_bw()+
  theme(axis.text.x = element_blank(),legend.position = "right")


############
##不同软件的比较
min_length <- 5
IBISouts <- read_delim(file ="A:/课题/kinship/errorIBD/data/dilution-degradationASA/ibis_HBdd.seg",delim = "\t",col_names   = F)
IBISouts$sample1 <- sapply(as.character(IBISouts$X1),FUN = function(x)strsplit(x,split=":",fixed = T)[[1]][1])
IBISouts$sample2 <- sapply(as.character(IBISouts$X2),FUN = function(x)strsplit(x,split=":",fixed = T)[[1]][1])
IBISouts$chr <- IBISouts$X3
IBISouts$start_detect <- IBISouts$X4
IBISouts$end_detect <- IBISouts$X5
IBISouts$IBD_type <- IBISouts$X6
IBISouts$length_detect <- (IBISouts$X5 - IBISouts$X4)/1000000
##1-5号染色体大于7cM的片段
IBISouts <- IBISouts[IBISouts$length_detect > min_length, ]
IBISouts$length2 <- IBISouts$length_detect * c(1,2)[(as.character( IBISouts$X6) == "IBD2")+1]
IBISouts <-IBISouts[, c("sample1","sample2","chr","start_detect","end_detect","IBD_type","length2")]
IBISouts <- aggregate(length2 ~ sample1 + sample2,data = IBISouts,FUN = "sum")

truffle_out <-  read_table(file ="A:/课题/kinship/errorIBD/data/dilution-degradationASA/truffle_HBdd.segments")
truffle_out$sample1 <- truffle_out$ID1
truffle_out$sample2 <- truffle_out$ID2
truffle_out$chr <- truffle_out$CHROM
truffle_out$start_detect <- truffle_out$POS
truffle_out$end_detect <- truffle_out$POS+truffle_out$LENGTH
truffle_out$IBD_type <- truffle_out$TYPE
truffle_out <- truffle_out[truffle_out$LENGTH > min_length, ]
truffle_out$length2 <- truffle_out$LENGTH
truffle_out <-truffle_out[, c("sample1","sample2","chr","start_detect","end_detect","IBD_type","length2")]
truffle_out <- aggregate(length2~sample1+sample2,truffle_out,FUN = sum)



ibdseq_d <-  read_delim(file ="A:/课题/kinship/errorIBD/data/dilution-degradationASA/ibdseq_HBdd.ibd",delim = "\t",skip = 0,trim_ws=T,col_names = F)
ibdseq_d$sample1 <-ibdseq_d$X1
ibdseq_d$sample2 <- ibdseq_d$X3
ibdseq_d$chr <- ibdseq_d$X5
ibdseq_d$start_detect <- ibdseq_d$X6
ibdseq_d$end_detect <- ibdseq_d$X7
ibdseq_d$length2 <- (ibdseq_d$end_detect - ibdseq_d$start_detect)/1000000
##1-5号染色体大于7cM的片段
ibdseq_d <- ibdseq_d[ibdseq_d$length2 > min_length, ]
ibdseq_d <-ibdseq_d[, c("sample1","sample2","chr","start_detect","end_detect","length2")]

ibdseq_out <- aggregate(length2 ~ sample1 + sample2,data = ibdseq_d,FUN = "sum")


clusIBD_out <-  read_delim("A:/课题/kinship/errorIBD/data/dilution-degradationASA/clusIBD_HBdd.IBD.details",delim = "\t",col_names = F)
#clusIBD_out <-  read_delim("A:/课题/kinship/errorIBD/data/dilution-degradationASA/HBfam_degratedDNA_20240729_115241.IBD.details",delim = "\t",col_names = F)

clusIBD_out$sample1 <-clusIBD_out$X1
clusIBD_out$sample2 <- clusIBD_out$X2
clusIBD_out$chr <- clusIBD_out$X3
clusIBD_out$start_detect <- clusIBD_out$X4
clusIBD_out$end_detect <- clusIBD_out$X5
clusIBD_out$IBD_type <- clusIBD_out$X7
clusIBD_out$length2 <- clusIBD_out$X6
##1-5号染色体大于7cM的片段
clusIBD_out <- clusIBD_out[clusIBD_out$length2 > min_length, ]
clusIBD_out <-clusIBD_out[, c("sample1","sample2","chr","start_detect","end_detect","IBD_type","length2")]
clusIBD_out <- aggregate(length2 ~ sample1 + sample2,data = clusIBD_out,FUN = "sum")

############################################################################################################################################
####################################################################################################################################
genome_length <- 2878.22

all_samples <- colnames(HBfam_ddASA_data)[10:ncol(HBfam_ddASA_data)]
all_samples_pairs <- data.frame( t(combn(all_samples,m=2)))
colnames(all_samples_pairs) <- c("sample1","sample2")
ddDNA_data1 <- merge(x=all_samples_pairs,y=IBISouts,by=c("sample1","sample2"),all.x=T)
ddDNA_data2 <- merge(x=all_samples_pairs,y=truffle_out,by=c("sample1","sample2"),all.x=T)
ddDNA_data3 <- merge(x=all_samples_pairs,y=ibdseq_out,by=c("sample1","sample2"),all.x=T)
ddDNA_data4 <- merge(x=all_samples_pairs,y=clusIBD_out,by=c("sample1","sample2"),all.x=T)

ddDNA_data1 <- merge(x=ddDNA_data1,y=ddDNA_data2,by=c("sample1","sample2"),suffixes=c(".IBIS",".TRUFFLE"),all=T)
ddDNA_data2 <- merge(x=ddDNA_data3,y=ddDNA_data4,by=c("sample1","sample2"),suffixes=c(".IBDseq",".clusIBD"),all = T)
ddDNA_data <- merge(ddDNA_data1,ddDNA_data2,all = T)
ddDNA_data <- melt(ddDNA_data,variable.name = "method", value.name = "IBDlength") 
ddDNA_data$Methods <- substr(ddDNA_data$method,9,nchar(as.character(ddDNA_data$method)))
ddDNA_data$Methods <- factor(ddDNA_data$Methods,levels = c("IBIS","TRUFFLE","IBDseq","clusIBD"))

ddDNA_data$id1 <- substr(ddDNA_data$sample1,1,3)
ddDNA_data$id2 <- substr(ddDNA_data$sample2,1,3)
ddDNA_data$type1 <- substr(ddDNA_data$sample1,5,nchar(as.character(ddDNA_data$sample1)))
ddDNA_data$type2 <- substr(ddDNA_data$sample2,5,nchar(as.character(ddDNA_data$sample2)))

ddDNA_data$phi <- ddDNA_data$IBDlength/(4*genome_length)
###未检测到表示为0
ddDNA_data$phi[is.na(ddDNA_data$phi)] <- 0
ddDNA_data$pairs <- apply(ddDNA_data, MARGIN = 1,FUN = function(x) paste(sort(x[6:7]),collapse = "-"))

rel_degrees <- c("A01-A01"=0,
                 "A03-A03"=0,
                 "A20-A20"=0,    
                 "A01-A03"=1,
                 "A01-A20"=2,
                 "A03-A20"=1,
                 "A01-D01"=3,
                 "D01-D01"=0,
                 "D03-D03"=0,
                 "D10-D10"=0,
                 "D01-D03"=1,
                 "D01-D10"=2,
                 "D03-D10"=1,
                 "A01-D03"=4,
                 "A01-D10"=5,
                 "A03-D01"=4,
                 "A03-D03"=5,
                 "A03-D10"=6,
                 "A20-D01"=5,
                 "A20-D03"=6,
                 "A20-D10"=7)
ddDNA_data$degrees <- factor( rel_degrees[ddDNA_data$pairs],levels = 0:7)


kin_thresholds <- numeric(9);kin_thresholds[9] <- -1
for (i in 1:8)kin_thresholds[i] <- 1/(2^((2*i+1)/2))
ddDNA_data$predicted_degrees <- factor( sapply(ddDNA_data$phi,FUN = function(x,y=kin_thresholds){which(x>y)[1]-1}),levels = 0:8)

ddDNA_data[ddDNA_data$type1 == "0.1ng" & ddDNA_data$type2 == "150bp",]

alltypes <- unique(ddDNA_data$type1)
in_types <- c("0.1ng","150bp","400bp")
#in_types <- alltypes[!(alltypes %in% c("0.1ng","150bp","400bp"))]

mydata <- ddDNA_data[(ddDNA_data$type1 %in% in_types) & (ddDNA_data$type2 %in% in_types) ,]

outs <- c()
for (m in c("IBIS","TRUFFLE","IBDseq","clusIBD")){
  x <- mydata[mydata$Methods == m,c("degrees","predicted_degrees")]
  x <- table(x)
  x <- melt(x)
  outs <- rbind(outs,data.frame(Methods=rep(m,nrow(x)),x))
}
outs$value2 <- outs$value
outs$value2[outs$degrees > outs$predicted_degrees] <- ""
outs$degrees <- factor(outs$degrees,levels = 0:7,labels = c("Identical","1st","2nd","3rd","4th","5th","6th","7th"))
outs$predicted_degrees <- factor(outs$predicted_degrees,levels = 0:8,labels = c("Identical","1st","2nd","3rd","4th","5th","6th","7th","unrelated"))
outs$prop <- 100 * outs$value/ (aggregate(value~degrees+Methods,outs,"sum")$value)[1:8]
outs$Methods <- factor(outs$Methods,levels = c("clusIBD","IBIS","TRUFFLE","IBDseq"))

mycol2 = c( "#E64B35B2", "#4DBBD5B2" ,"#00A087B2", "#3C5488B2", "#F39B7FB2",
            "#8491B4B2", "#91D1C2B2", "#DC0000B2", "#7E6148B2")

aggregate(Values~types +Error_type,acc_data,mean)
g_genotype_acc <- ggplot(data = acc_data,mapping = aes(x=Error_type,y=Values,color=Error_type))+
  geom_violin( scale = "width",show.legend = FALSE)+
  geom_point()+
  scale_color_manual(values =mycol2[1:4] )+
  #geom_boxplot(alpha=0.5)+
  facet_wrap(vars(types),nrow = 2)+
  labs(color="Error types",x="",y="Error rates")+
  theme_bw()+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "right",
        plot.margin = margin(t=0.3, b=0.1, l=1.5, r=0.25, "cm"),
        title = element_text(size = 10))

g_kinship <- ggplot(outs,mapping = aes(x=predicted_degrees  ,y=degrees))+
  geom_raster(mapping = aes(fill=prop))+
  facet_wrap(vars(Methods),nrow = 2,drop=FALSE)+
  #scale_fill_gradient2(low="#003366", high="#990033", mid="white")+
  scale_fill_gradient2(low="#003366", high="#DC0000B2", mid="gray99",midpoint = 30)+
  geom_text(mapping = aes(label=value2),size=2,vjust=0.5)+
  #geom_text(mapping = aes(label=prop2),size=3,vjust=1)+
  labs(x="Predicted degrees of relatedness",y="True degrees of relatedness",fill="Percentage(%)")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90),
        plot.margin = margin(t=0.1, b=0.3, l=1, r=0.4, "cm"),
        title = element_text(size = 10))

g <- ggarrange(g_genotype_acc,g_kinship,nrow = 2,
               align = c("v"),
               heights = c(1,1.5),
               font.label = list(size = 10),
               labels = c("A","B"))

ggsave(filename =  paste0("A:/课题/kinship/errorIBD/results/Fig 5 kinship inference for degraded DNA from HB family",Sys.Date(),".jpeg"),g,height=6,width =6)

write.csv(acc_data,file = "A:/课题/kinship/errorIBD/results/Supplementary Table 1.csv",row.names = F)

##clusIBD
#sum(c(18,19,11,7,2,0,0,0))/sum((aggregate(value~degrees+Methods,outs,"sum")$value)[1:8])

g_genotype_acc+guides(color = guide_legend(override.aes = list(shape = 16, size = 3, colour = mycol2[1:4] )))  # 注意：colour = NA 可能不会去除边框，因为点默认没有边框
  

sum(c(18,12,10,2,2,0,0,0))/sum((aggregate(value~degrees+Methods,outs,"sum")$value)[1:8])
c(18,12,10,2,2,0,0,0)/(aggregate(value~degrees+Methods,outs,"sum")$value)[1:8]

#IBIS
sum(c(1,6,3,2,2))/sum((aggregate(value~degrees+Methods,outs,"sum")$value)[1:8])
#IBIS
sum(c(18,6,3))/sum((aggregate(value~degrees+Methods,outs,"sum")$value)[1:8])


##clusIBD
c(18,19,11,7,2,0,0,0)/(aggregate(value~degrees+Methods,outs,"sum")$value)[1:8]
#IBIS
sum(c(1,6,3,2,2))/sum((aggregate(value~degrees+Methods,outs,"sum")$value)[1:8])
#IBIS
sum(c(18,6,3))/sum((aggregate(value~degrees+Methods,outs,"sum")$value)[1:8])



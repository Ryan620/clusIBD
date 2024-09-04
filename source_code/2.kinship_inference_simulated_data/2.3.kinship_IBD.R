library(readr)


###############################################################################################################################
###############################################################################################################################
error_rates <- c(0,0.001,0.005,0.01,0.05,0.1,0.2)
names(error_rates) <- paste0("error",0:6)
##模拟家系只需要D F O Q T U S
sim_pairs <- c("O-T","D-F","D-O","D-T","D-U","O-Q","O-U","T-U","D-S","F-S")
sim_degree <- c(1,1,2,3,4,5,6,7,0,0)

pair_num <- length(sim_pairs)
##所有样本对
sim_pairs <- matrix(unlist(strsplit(sim_pairs,split = "-")),ncol = 2,byrow = T)
sample_data <- data.frame(simID=rep(paste0("sim",7799:7808),each=18*pair_num),
                          famID= rep(rep(paste0("def",1:18),each=pair_num),10),
                          id1 = rep(sim_pairs[,1],180),
                          id2 = rep(sim_pairs[,2],180),
                          degree=rep(sim_degree,180)
)
sample_data$sim.fam=paste(sample_data$simID,sample_data$famID,sep = ".")
sample_data$ID1 <- paste(sample_data$sim.fam,sample_data$id1,sample_data$sim.fam,sample_data$id1,sep = "_")
sample_data$ID2 <- paste(sample_data$sim.fam,sample_data$id2,sample_data$sim.fam,sample_data$id2,sep = "_")
sample_data$sample_pair2 <- paste(sample_data$sim.fam,sample_data$id1,sample_data$id2,sep = "_")
sample_data$sample_pair3 <- paste(sample_data$ID1,sample_data$ID2,sep = "-")

write.table(sample_data[,c("ID1","ID2")],file = "target_samplepairs.txt",col.names = F,row.names = F,quote = F,sep = "\t")



y <- read_delim("error_0.005.bim",delim="\t",col_names = F)
y <- aggregate(X4~X1,y,max)
genome_length <- sum(as.numeric(y$X4),na.rm=T)/1000000

error_rates <- c(0,0.005,0.01,0.05,0.1,0.2)
min_length <- 5
#genome_length <- 2878.22
##IBIS结果
IBIS_data <- c()
for (i in 1:length(error_rates)){
  IBISouts <- read_delim(file = paste0("./results/ibis_error_",error_rates[i],".seg"),delim = "\t",col_names   = F)
  #IBISouts$Individual1 <- gsub(IBISouts$Individual1,pattern = ":",fixed = T,replacement ="_" )
  #IBISouts$Individual2 <- gsub(IBISouts$Individual2,pattern = ":",fixed = T,replacement ="_" )
  
  IBISouts$sample1 <- sapply(as.character(IBISouts$X1),FUN = function(x)strsplit(x,split=":",fixed = T)[[1]][1])
  IBISouts$sample2 <- sapply(as.character(IBISouts$X2),FUN = function(x)strsplit(x,split=":",fixed = T)[[1]][1])
  IBISouts$chr <- IBISouts$X3
  IBISouts$start_detect <- IBISouts$X4
  IBISouts$end_detect <- IBISouts$X5
  IBISouts$length_detect <- (IBISouts$X5 - IBISouts$X4)/1000000
  ##1-5号染色体大于7cM的片段
  IBISouts <- IBISouts[IBISouts$length_detect > min_length, ]
  IBISouts$length_detect2 <- IBISouts$length_detect * c(1,2)[(as.character( IBISouts$X6) == "IBD2")+1]
  #IBISouts$error_rate <- rep(error_rates[k],nrow(IBISouts))
  
  IBISouts$sample_pair3 <- paste(IBISouts$sample1,IBISouts$sample2,sep = "-")
  IBISouts <- aggregate(length_detect2~sample_pair3,IBISouts,FUN = sum)
  x <- merge(x=sample_data,IBISouts,by="sample_pair3",all.x=TRUE,sort = F)
  x$Kinship_Coefficient <- x$length_detect2/(genome_length*4)
  #IBD_ibis <- x[,c(2,3,4,5,8,9,6,13)]
  IBD_ibis <- x[,c("simID","famID","id1" ,"id2" ,"degree","Kinship_Coefficient" )]
  colnames(IBD_ibis) <- c("simID","famID","id1" ,"id2" ,"degree","phi_IBIS" )
  IBD_ibis$error_type <- rep(error_rates[i],nrow(IBD_ibis))
  IBIS_data <- rbind(IBIS_data,IBD_ibis)
  print(i)
}

print("finished for IBIS...............................................")

################################################################
#truffle
TRUFFLE_data <- c()
for (i in 1:length(error_rates)){
  truffle_out <-  read_table(file = paste0("./results/truffle_error_",error_rates[i],".segments"))
  truffle_out$sample_pair3 <- paste(truffle_out$ID1,truffle_out$ID2,sep = "-")
  truffle_out <- truffle_out[truffle_out$LENGTH > min_length, ]
  #truffle_out$lengths <- truffle_out$LENGTH * c(1,2)[(as.character( truffle_out$TYPE) == "IBD2")+1]
  truffle_out$lengths <- truffle_out$LENGTH
  truffle_out <- aggregate(lengths~sample_pair3,truffle_out,FUN = sum)
  truffle_out <- merge(x=sample_data[,c(1:5,10)],y=truffle_out,by="sample_pair3",all.x=TRUE,sort = F)
  
  truffle_out$phi_TRUFFLE <- truffle_out$lengths/(genome_length*4)
  truffle_out$error_type <- rep(error_rates[i],nrow(truffle_out))
  TRUFFLE_data <- rbind(TRUFFLE_data,truffle_out)
  print(i)
}
print("finished for TRUFFLE...............................................")

##IBDseq
IBDseq_data <- c()
for (i in 1:length(error_rates)){
  ibdseq_d <-  read_delim(file =paste0("./results/ibdseq_error_",error_rates[i],".ibd"),delim = "\t",skip = 0,trim_ws=T,col_names = F)
  ibdseq_d$ID1 <-ibdseq_d$X1
  ibdseq_d$ID2 <- ibdseq_d$X3
  ibdseq_d$sample_pair3 <- paste(ibdseq_d$ID1,ibdseq_d$ID2,sep = "-")
  ibdseq_d$chr <- ibdseq_d$X5
  ibdseq_d$start_detect <- ibdseq_d$X6
  ibdseq_d$end_detect <- ibdseq_d$X7

  ibdseq_d$length_detect <- (ibdseq_d$end_detect - ibdseq_d$start_detect)/1000000
  ##1-5号染色体大于7cM的片段
  ibdseq_d <- ibdseq_d[ibdseq_d$length_detect > min_length, ]
  
  ibdseq_d <- aggregate(length_detect ~ ID1 + ID2,data = ibdseq_d,FUN = "sum")
  
  x <- merge(x=sample_data,y=ibdseq_d,all.x=TRUE,sort = F)
  idx <- grepl(x$ID1,pattern = "_D_") &  grepl(x$ID2,pattern = "_F_")
  x$length_IBD <- x$length_detect
  x$length_IBD[idx] <- x$length_IBD[idx]*4/3
  x$phi_IBDseq <- x$length_IBD/(genome_length*4)
  x$error_type <- rep(error_rates[i],nrow(x))
  IBDseq_data <- rbind(IBDseq_data,x)
  print(paste(i,"finished"))
}

print("finished for IBDSEQ...............................................")

#clusIBD

clusIBD_data <- c()
for (i in 1:length(error_rates)){
  clusIBD_out <-  read_delim(file = paste0("./results/clusIBD_error_",error_rates[i],".IBD.details"),delim = "\t",col_names = F)
  clusIBD_out$sample_pair3 <- paste(clusIBD_out$X1,clusIBD_out$X2,sep = "-")
  clusIBD_out <- clusIBD_out[clusIBD_out$X6 > min_length, ]
  #clusIBD_out$lengths <- clusIBD_out$X6 * c(1,2)[(as.character( clusIBD_out$X7) == "IBD2")+1]
  clusIBD_out$lengths <- clusIBD_out$X6 
  clusIBD_out <- aggregate(lengths~sample_pair3,clusIBD_out,FUN = sum)
  clusIBD_out <- merge(x=sample_data[,c(1:5,10)],y=clusIBD_out,by="sample_pair3",all.x=TRUE,sort = F)
  
  clusIBD_out$phi_clusIBD <- clusIBD_out$lengths/(genome_length*4)
  clusIBD_out$error_type <- rep(error_rates[i],nrow(clusIBD_out))
  clusIBD_data <- rbind(clusIBD_data,clusIBD_out)
  print(i)
}

print("finished for clusIBD...............................................")


###合并不同软件的数据
x  <- merge(x=IBIS_data[,c("simID", "famID", "id1", "id2", "degree","error_type","phi_IBIS")],
            y=TRUFFLE_data[,c("simID", "famID", "id1", "id2", "degree","error_type","phi_TRUFFLE")])
x  <- merge(x=x,y=IBDseq_data[,c("simID", "famID", "id1", "id2", "degree","error_type","phi_IBDseq")])
y  <- merge(x=x,y=clusIBD_data[,c("simID", "famID", "id1", "id2", "degree","error_type","phi_clusIBD")])

kinship_data <- cbind(IBIS_data[,c("simID", "famID", "id1", "id2", "degree","error_type")],
                     phi_IBIS=IBIS_data$phi_IBIS,
                     phi_TRUFFLE=TRUFFLE_data$phi_TRUFFLE,
                     phi_IBDseq=IBDseq_data$phi_IBDseq,
                     phi_clusIBD=clusIBD_data$phi_clusIBD)

write.table(y,file = "kinship results of four methods.txt",col.names = T,row.names = F,quote = F,sep = "\t")




################################################################################################
###本地电脑
mycols = c( "#CC79A7", "#F39B7FB2", "#00A087B2", "#3C5488B2","#4DBBD5B2" ,
            "#8491B4B2", "#91D1C2B2", "#7E6148B2", "#DC0000B2")[1:9]
kinship_data  <- read.table(file = "A:/课题/kinship/errorIBD/data/kinship/kinship results of four methods.txt",header  = T,sep = "\t")
kinship_data$type <-  paste("degree",kinship_data$degree)
kinship_data$type[kinship_data$id1 == "O" & kinship_data$id2 == "T"] <- "parent-child"
kinship_data$type[kinship_data$id1 == "D" & kinship_data$id2 == "F"] <- "full siblings"
kinship_data$type[kinship_data$id1 == "D" & kinship_data$id2 == "S"] <- "unrelated"
kinship_data$type[kinship_data$id1 == "F" & kinship_data$id2 == "S"] <- "unrelated"
kinship_data$type <- factor(kinship_data$type,levels = c("parent-child","full siblings",paste("degree",2:7),"unrelated"))

##计算相关系数
r_values <- matrix(NA,nrow = 6,ncol = 3)
for (i in 1:6){
  for (j in 1:3){
    subdata <- kinship_data[kinship_data$error_type == c(0,0.005,0.01,0.05,0.1,0.2)[i],]
    r <- cor(x=subdata[[6+j]],y=subdata[[10]],use = "pairwise.complete.obs")
    r_values[i,j] <- r^2
  }
}
colnames(r_values) <- c("IBIS","TRUFFLE","IBDseq")
r_values <- cbind(error_type=c(0,0.005,0.01,0.05,0.1,0.2),r_values)
rlabels <- melt(data.frame(r_values),id.vars = "error_type",variable.name = "Methods",value.name = "Rvalue")
rlabels$Rvalue2 <- paste0("r2=", round(rlabels$Rvalue,3))

x <- melt(kinship_data,id.vars = c("error_type","type","phi_clusIBD"),
          measure.vars = c("phi_IBDseq" ,"phi_IBIS","phi_TRUFFLE"   ),
          variable.name = "Methods",
          value.name = "phi")
x$phi[is.na(x$phi)] <- 0
x$Methods <- gsub(x$Methods,pattern = "phi_",replacement = "")
x$Methods <- factor(x$Methods,levels = c("IBIS","TRUFFLE","IBDseq","KING"))

g <- ggplot(data = x,mapping = aes(x=phi_clusIBD,y=phi))+
  geom_point(alpha=0.3,aes(color=type))+
  scale_colour_manual(values = mycols)+
  #geom_smooth(method = "lm",level=0.99,color="black",se=T,size=0.8)+
  geom_abline(slope = 1,size=1,colour="black",linetype='dotted')+
  lims(y=c(0.02,0.22),y=c(-0.01,0.25))+
  labs(x="kinship coefficients (clusIBD)",y="kinship coefficients",color="")+
  #geom_text(data=dx,mapping = aes(x=x,y=y,label=labels),inherit.aes=F)+
  geom_text(data = rlabels,mapping = aes(x=0.1,y=0.24,label=Rvalue2),inherit.aes = F,size=2)+
  facet_grid(Methods~error_type)+
  ##调整图例形状的大小
  guides(color = guide_legend(override.aes = list(size = 4)))

ggsave(filename = paste0("A:/课题/kinship/errorIBD/results/Supplementary Figure 4 kinship coefficients",Sys.Date(),".jpeg"),width = 10,height = 6,g)  


##准确率
x <- melt(kinship_data,id.vars = c("error_type","degree","type"),
          measure.vars = c("phi_IBDseq" ,"phi_IBIS","phi_TRUFFLE" ,"phi_clusIBD"  ),
          variable.name = "Methods",
          value.name = "phi")
x$phi[is.na(x$phi)] <- 0
kin_thresholds <- numeric(9);kin_thresholds[9] <- -1
for (i in 1:8)kin_thresholds[i] <- 1/(2^((2*i+1)/2))

x$predicted_degree <- sapply(x$phi,FUN = function(x,y=kin_thresholds){which(x>y)[1]-1})
x$predicted_degree[x$predicted_degree==8] <- 0
x$degree <- factor(x$degree)
x$predicted_degree <- factor(x$predicted_degree)

acc_data <- c()
acc_data2 <- c()  ##区分1级的亲子和全同胞
etp <- as.character(unique(x$error_type))
mtp <- as.character(unique(x$Methods))
for (i in 1:length(etp) ){
  for (j in 1:length(mtp)){
    sub_mydata <- x[x$error_type==etp[i] & x$Methods==mtp[j],]
    y <- table(sub_mydata[,c("degree","predicted_degree")])
    y2 <- table(sub_mydata[,c("type","predicted_degree")])
    acc <- diag(y)/rowSums(y)
    acc2 <- c(y2[1,2],y2[2,2],y2[3,3],y2[4,4],y2[5,5],y2[6,6],y2[7,7],y2[8,8],y2[9,1])/rowSums(y2)
    acc_data <- rbind(acc_data,
                      data.frame(error_group=rep(etp[i],length(acc)),
                                 relatedness=names(acc),
                                 Methods=rep(mtp[j],length(acc)),
                                 Accuracy=acc))
    acc_data2 <- rbind(acc_data2,
                       data.frame(error_group=rep(etp[i],length(acc2)),
                                  relatedness=names(acc2),
                                  Methods=rep(mtp[j],length(acc2)),
                                  Accuracy=acc2))
  }
}

acc_data$relatedness <- factor(acc_data$relatedness,levels = c(1:7,0), labels  = c("1st","2nd","3rd","4th","5th","6th","7th","unrelated"))
acc_data$Methods <- gsub(acc_data$Methods,pattern = "phi_",replacement = "")
acc_data2$relatedness <- factor(acc_data2$relatedness,levels = names(acc2), labels  = c("parent-child","full siblings","2nd","3rd","4th","5th","6th","7th","unrelated"))
acc_data2$Methods <- gsub(acc_data2$Methods,pattern = "phi_",replacement = "")
acc_data2$error_group <- factor(acc_data2$error_group,levels = unique(acc_data2$error_group), labels  = c("0%","0.5%","1%","5%","10%","20%"))


g <- ggplot(acc_data2,aes(x= Methods,y=Accuracy,fill=Methods))+
  geom_bar(stat = "identity")+
  scale_fill_manual(values = mycols[1:4])+
  geom_text(aes(label=round(Accuracy,2)),size=2,angle=90,vjust="top",nudge_y =0.1)+
  facet_grid(error_group ~ relatedness)+
  lims(y=c(0,1.2))+
  labs(y="Recall rates",x="")+
  theme_bw()+
  theme(
    #panel.background = element_rect(fill = 'white', color = 'white'), panel.border = element_blank(),
    #text=element_text(size = 5),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    #axis.text.x = element_text(angle = 90))
    axis.text.x = element_blank(), axis.ticks.x = element_blank())

ggsave(filename = paste0("A:/课题/kinship/errorIBD/results/Figure 3. comparison of accurarcies for four methods ",Sys.Date(),".jpeg"),g,height=6,width = 8)
#ggsave(filename = paste0("A:/课题/kinship/errorIBD/results/Figure 3. comparison of accurarcies for four methods.jpeg"),g,height=6,width = 8)



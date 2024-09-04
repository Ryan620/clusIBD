library(ggplot2)
library(readr)
library(gridExtra)
#####
clusIBD_show <- function(s1,s2,IBDfile){
  clusIBD_out <-  read_delim(IBDfile,delim = "\t",col_names = F)
  clusIBD_out$sample1 <-clusIBD_out$X1
  clusIBD_out$sample2 <- clusIBD_out$X2
  clusIBD_out$chr <- clusIBD_out$X3
  clusIBD_out$start_detect <- clusIBD_out$X4
  clusIBD_out$end_detect <- clusIBD_out$X5
  clusIBD_out$IBD_type <- clusIBD_out$X7
  clusIBD_out$length2 <- clusIBD_out$X6
  clusIBD_out$IBD_type <- factor(clusIBD_out$IBD_type,levels = c("IBD1","IBD2"))
  chr_lengths <- c(249168108,243042288,197794954,190912862,180690007,
                   170919684,159119221,146276852,141027784,135449605,
                   134942403,133728845,115072771,107257542,102394005,
                   90140205,81046645,78013879,59045905,62900718,48075537,51203313)/1000000
  bg_data_chrs<-data.frame(chr=factor(1:22,levels=1:22),start_pos=rep(1,22),end_pos=chr_lengths)
  s1 <- as.character(s1);s2 <- as.character(s2)
  subdata <- clusIBD_out[clusIBD_out$sample1 == s1 & clusIBD_out$sample2 == s2,]
  subdata$IBD_type <- factor(subdata$IBD_type,levels = c("IBD1","IBD2"))
  ##交换s1和s2
  if (nrow(subdata) == 0){
    x=s1
    s1=s2
    s2=x
    subdata <- clusIBD_out[clusIBD_out$sample1 == s1 & clusIBD_out$sample2 == s2,]
  }
  g_segment <- NULL
  if (nrow(subdata)>0){
    total_lengths <- round(sum(subdata$length2),2)
    g_segment <- ggplot() +
      geom_rect(data = bg_data_chrs, aes(xmin = start_pos, xmax = end_pos, ymin = 0, ymax = 1),fill = "gray80")+
      geom_rect(data = subdata, aes(xmin = start_detect, xmax = end_detect, ymin = 0, ymax = 1,fill=IBD_type ))+
      facet_wrap(vars(chr),ncol = 1,strip.position= 'left')+
      labs(x = "Position (Mb)", y = "", fill = "IBD type",
           title = paste0(s1," - ",s2,", Total IBD: ",total_lengths," Mb"))+
      theme(
        legend.position = c(0.8,0.2), 
        #legend.key = element_rect(color = "black"),
        panel.background = element_rect(fill = 'white', color = 'white'), panel.border = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.text.y = element_blank(), axis.ticks.y = element_blank()
      )
  }
  return(g_segment)
}
relationship_distance<-function(pedigreedata,s1,s2){
  result<-c(0,0)
  names(result)<-c("miosisN","degreeN")
  lineages<-NA
  pedigree_s1<-pedigree_finder(pedigreedata,s1)
  pedigree_s2<-pedigree_finder(pedigreedata,s2)
  findIDX<-FALSE
  #########################################################
  if (s2 %in% unlist(pedigree_s1)){						##直系相关 s2 is an ancestor of s1
    for (i in 1:length(pedigree_s1)){
      parentnames<-unlist(pedigree_s1[i])
      if (s2 %in% parentnames) {
        miosisN=i;degreeN=i
        findIDX<-TRUE
        break}
    }
    ##样本之前关系单一路径
    if (findIDX){
      result<-c(miosisN,degreeN)
      lineages<-c(s2,pedigree_chain(pedigreedata=pedigreedata,pedigree_list=pedigree_s1,s=s2),s1)
      
    }
  }
  ###################################################
  if (s1 %in% unlist(pedigree_s2)){						###直系相关
    for (i in 1:length(pedigree_s2)){
      parentnames<-unlist(pedigree_s2[i])
      if (s1 %in% parentnames) {
        miosisN=i;degreeN=i				
        findIDX<-TRUE
        break}
    }
    ##样本之前关系单一路径
    if (findIDX){
      result<-c(miosisN,degreeN)
      lineages<-c(s1,pedigree_chain(pedigreedata=pedigreedata,pedigree_list=pedigree_s2,s=s1),s2)
      
    }
  }
  if (findIDX) return(c(relateness=list(result),lineages=list(lineages)))
  if (length(pedigree_s1)==0 | length(pedigree_s2)==0) return(c(relateness=list(result),lineages=list(lineages)))
  #######################################################
  findIDX=FALSE
  for (i in 1:length(pedigree_s1)){						###旁系相关
    for (j in 1:length(pedigree_s2)){
      parentnames1<-pedigree_s1[[i]]
      parentnames2<-pedigree_s2[[j]]
      coreperson<-parentnames1[parentnames1 %in% parentnames2]
      t<-sum(as.numeric(parentnames1 %in% parentnames2),na.rm=TRUE)
      if (t==1){miosisN<-i+j;degreeN<-i+j;ancestor<-"ONE";findIDX<-TRUE}				#一个共同祖先
      if (t==2){miosisN<-i+j;degreeN<-i+j-1;ancestor<-"TWO";findIDX<-TRUE}			#两个共同祖先
      if (findIDX){
        result<-c(miosisN,degreeN)
        lineage1<-c(pedigree_chain(pedigreedata=pedigreedata,pedigree_list=pedigree_s1,s=coreperson[1]),s1)
        lineage2<-c(pedigree_chain(pedigreedata=pedigreedata,pedigree_list=pedigree_s2,s=coreperson[1]),s2)
        lineages<-c(rev(lineage1),ancestor,lineage2)
        return(c(relateness=list(result),lineages=list(lineages)))
      }
    }
    
  }
  
  c(relateness=list(result),lineages=list(lineages))
}

pedigree_chain<-function(pedigreedata,pedigree_list,s){
  lineages<-NULL
  for (i in 1:length(pedigree_list)) {
    if (s %in% pedigree_list[[i]]) {
      idx=i
      if (idx==1) return(NULL)
      parent<-s
      lineages<-c()
      for (j in (idx-1):1) {
        x<-pedigreedata[which(pedigreedata[,1] %in% pedigree_list[[j]]), 1:4]
        for (k in 1:nrow(x)){
          if (parent %in% x[k,]) {lineages<-c(lineages,x[k,1]);parent<-x[k,1]}
        }
        
      }
    }
  }
  
  return(lineages)
}
##返回样本性别
find_sampleSex<-function(pedigreedata,s){
  if (length(s)==0) return(NA)
  if (length(which(is.na(s)))==length(s)) return(NA)
  sex<-c()
  for (x in s){y<-pedigreedata[as.character(pedigreedata[,1])==x,2];if (length(y)==0)y=x;sex<-c(sex,y)}
  return(sex)
}

###样本的上游系谱（父母，父母的父母。。。）
pedigree_finder<-function(pedigreedata,s){
  pedigree_s<-list()
  n=0
  new_s<-s
  parentnames<-c(NA)
  while(length(parentnames)!=0){						#反复寻找样本的父母
    if (n>50) break		
    parentnames<-c()
    for (s1 in new_s){
      x<-as.character(parent_finder(pedigreedata,s1))
      parentnames<-c(parentnames,x[!is.na(x)])
    }	
    if (length(parentnames)!=0){pedigree_s<-c(pedigree_s,list(parentnames));n=n+1}
    new_s<-parentnames	
  }
  if (n>0)names(pedigree_s)<-paste("generation",1:n,sep="_")
  return(pedigree_s)
}

###样本的父母
parent_finder<-function(pedigreedata,s){
  fam_members<-pedigreedata[,1]
  parentnames<-unlist(pedigreedata[fam_members==s,3:4])
  return(parentnames)
}

kin_sex<- function(pedigreedata,samplepairs){
  relateness_para<-c()
  unique_sex_chian <-c()
  for (i in 1:nrow(samplepairs)){
    s1<-as.character(samplepairs[i,1])
    s2<-as.character(samplepairs[i,2])
    x=relationship_distance(pedigreedata,s1=s1,s2=s2)
    relateness<-x$relateness
    lineages<-x$lineages
    if (sum(relateness)==0) {
      mid_ancestor <- paste(sort(c(find_sampleSex(pedigreedata,s=s1),find_sampleSex(pedigreedata,s=s2))),collapse="")
      relateness_para<-rbind(relateness_para,c(relateness,mid_ancestor));next}
    sexes<-find_sampleSex(pedigreedata,s=lineages)
    
    #半同胞系，确定共有祖先性别
    if ("ONE" %in% lineages){
      idx <- which(lineages=="ONE")
      p1 <- parent_finder(pedigreedata,s=lineages[idx-1])
      p2 <- parent_finder(pedigreedata,s=lineages[idx+1])
      shareAcestor <- p1[p1 %in% p2]
      shareAcestor_sex <- find_sampleSex(pedigreedata,s=shareAcestor)
      sexes[idx] <- paste0(1,shareAcestor_sex)
    }
    #全同胞系，确定共有祖先性别
    if ("TWO" %in% lineages){
      idx <- which(lineages=="TWO")
      sexes[idx] <- "FM"
    }
    mid_ancestor<-NA
    
    mid_ancestor<-paste(sexes,collapse="-")
    ##由于某些旁系关系是对称的，需要去除重复的链,
    
    if (!(mid_ancestor %in% unique_sex_chian) & any(c(grepl(mid_ancestor,pattern="1"),grepl(mid_ancestor,pattern="FM")))) {		
      sex_chain<-unlist(strsplit(mid_ancestor,split="-"))			
      sex_chain<-paste(rev(sex_chain),collapse="-")
      
      if (sex_chain %in% unique_sex_chian)mid_ancestor <-sex_chain
      
      
    }
    unique_sex_chian<-c(unique_sex_chian,mid_ancestor)
    
    #排除
    exchange=FALSE
    if (any(c(grepl(mid_ancestor,pattern="1"),grepl(mid_ancestor,pattern="FM")))){
      sex_chain<-unlist(strsplit(mid_ancestor,split="-"))	
      if (length(sex_chain)>2){
        idx <- max(which(sex_chain=="1F"),which(sex_chain=="1M"),which(sex_chain=="FM"))
        if (length(idx) != 0){
          #两侧，谁的数目少排前，相同数目，F排前
          
          sex_chain1<- sex_chain[1:(idx-1)]
          sex_chain2<- rev(sex_chain[(idx+1):length(sex_chain)])
          if (length(sex_chain1) > length(sex_chain2)) exchange=TRUE
          if (length(sex_chain1) == length(sex_chain2)){
            fm <- c("F"=1,"M"=2)
            for (j in 1:length(sex_chain1)){
              if (fm[sex_chain1[j]]<fm[sex_chain2[j]])break
              if (fm[sex_chain1[j]]>fm[sex_chain2[j]]){exchange=TRUE;break}
            }
          }
        }
      }
    }
    if (mid_ancestor=="NA") stop("")
    mid_ancestor <- ifelse(exchange,paste(rev(sex_chain),collapse="-"),mid_ancestor)
    relateness_para<-rbind(relateness_para,c(relateness,mid_ancestor))
    
  }
  sample_relateness<-data.frame(cbind(samplepairs,relateness_para))
  colnames(sample_relateness) <- c("s1","s2","miosisN", "degreeN","sex_chain")
  sample_relateness
}





############
##不同软件的比较
min_length <- 5
IBISouts <- read_delim(file ="A:/课题/kinship/errorIBD/data/ancientDNA/ibis_ancientDNA.seg",delim = "\t",col_names   = F)
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

truffle_out <-  read_table(file ="A:/课题/kinship/errorIBD/data/ancientDNA/truffle_ancientDNA.segments")
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

clusIBD_out <-  read_delim("A:/课题/kinship/errorIBD/data/ancientDNA/clusIBD_ancientDNA.IBD.details",delim = "\t",col_names = F)
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



ibdseq_d <-  read_delim(file ="A:/课题/kinship/errorIBD/data/ancientDNA/ibdseq_ancientDNA.ibd",delim = "\t",skip = 0,trim_ws=T,col_names = F)
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



##############################################################################################
##本地电脑
clusIBD_out <-  read_delim("A:/课题/kinship/errorIBD/data/ancientDNA/clusIBD_ancientDNA_defaut.IBD.details",delim = "\t",col_names = F)
#clusIBD_out <-  read_delim("A:/课题/kinship/errorIBD/data/ancientDNA/clusIBD_ancientDNA_n200.IBD.details",delim = "\t",col_names = F)
clusIBD_out <-  read_delim("A:/课题/kinship/errorIBD/data/ancientDNA/clusIBD_ancientDNA_n300.IBD.details",delim = "\t",col_names = F)
#clusIBD_out <-  read_delim("A:/课题/kinship/errorIBD/data/ancientDNA/clusIBD_ancientDNA_n400.IBD.details",delim = "\t",col_names = F)
clusIBD_out <-  read_delim("A:/课题/kinship/errorIBD/data/ancientDNA/clusIBD_ancientDNA_n500.IBD.details",delim = "\t",col_names = F)
fileout <- c("A:/课题/kinship/errorIBD/data/ancientDNA/clusIBD_ancientDNA_default.IBD.details",
             "A:/课题/kinship/errorIBD/data/ancientDNA/clusIBD_ancientDNA_n300.IBD.details",
             "A:/课题/kinship/errorIBD/data/ancientDNA/clusIBD_ancientDNA_n500.IBD.details")
figout <- vector("list",6)
for (index in 1:3){
  clusIBD_out <-  read_delim(fileout[index],delim = "\t",col_names = F)
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
  
  
  mycol2 = c( "#E64B35B2", "#4DBBD5B2" ,"#00A087B2", "#3C5488B2", "#F39B7FB2",
              "#8491B4B2", "#91D1C2B2", "#7E6148B2")
  sample_names <- read.table("A:/课题/kinship/errorIBD/data/ancientDNA/ancientDNA_filter_renamed.fam",header = F,sep = " ",stringsAsFactors = F)
  sample_names <- data.frame( t(combn(sample_names$V2,2)))
  colnames(sample_names) <- c("sample1","sample2")
  
  ancientDNA_data1 <- merge(sample_names,IBISouts,all.x = T)
  ancientDNA_data2 <- merge(sample_names,truffle_out,all.x = T)
  ancientDNA_data3 <- merge(sample_names,ibdseq_out,all.x = T)
  ancientDNA_data4 <- merge(sample_names,clusIBD_out,all.x = T)
  ancientDNA_data <- data.frame(sample1=sample_names$sample1,
                                sample2=sample_names$sample2,
                                IBIS=ancientDNA_data1$length2,
                                TRUFFLE=ancientDNA_data2$length2,
                                #IBDseq=ancientDNA_data3$length2,
                                IBDseq=rep(0,nrow(sample_names)),
                                clusIBD=ancientDNA_data4$length2)
  ###引入系谱信息
  
  sample_pairs = sample_names
  pedigreedata <- read.table("A:/课题/kinship/errorIBD/data/ancientDNA/pedigree.txt",header = T,sep = "\t",stringsAsFactors = F)
  sample_pairs_data <- cbind(sample_pairs,miosisN=rep(NA,nrow(sample_pairs)),degreeN=rep(NA,nrow(sample_pairs)))
  for (i in 1:nrow(sample_pairs_data)) sample_pairs_data[i,3:4] <- relationship_distance(pedigreedata,s1=sample_pairs_data[i,1],s2=sample_pairs_data[i,2])[[1]]
  
  
  x <- merge(x=sample_pairs_data,y= ancientDNA_data,all.x = T)
  x <- melt(x,measure.vars = c("IBIS","TRUFFLE","IBDseq"))
  ##提出NC5m
  x <- x[x$sample1 != "NC5m" & x$sample2 != "NC5m", ]
  ##
  x <- x[x$degreeN < 8,]
  x$clusIBD[is.na(x$clusIBD)] <- 0
  x$value[is.na(x$value)] <- 0
  
  x$degreeN <- factor(x$degreeN,  levels = c(1:7,0),labels =  c("1st","2nd","3rd","4th","5th","6th","7th","UN"))
  g0 <- ggplot(data = x,mapping = aes(x=clusIBD,y=value,color=degreeN))+ 
    geom_point(mapping = aes(shape=variable))+
    geom_abline(slope = 1,size=0.7,colour="black",linetype='dotted')+
    theme_bw()+
    scale_color_manual(values =mycol2)+
    labs(y="IBD lengths (Mb)",
         x="IBD lengths (Mb) with clusIBD",
         color="Relatedness",shape="Methods")
  
  
  
  ##kinship inference
  kin_thresholds <- numeric(9);kin_thresholds[9] <- -1
  for (i in 1:8)kin_thresholds[i] <- 1/(2^((2*i+1)/2))
  x$predicted_degrees <- factor( sapply(x$clusIBD/(4*2878.22),FUN = function(x,y=kin_thresholds){which(x>y)[1]-1}),levels = 0:8)
  x$predicted_degrees <- factor(x$predicted_degrees,levels = 1:8,labels = c("1st","2nd","3rd","4th","5th","6th","7th","UN"))
  y <- x[x$variable == "IBIS",]
  ancientout <- table(y[,c("degreeN","predicted_degrees")])
  sumrelationships <- rowSums(ancientout)
  for (i in 1:8){
    for (j in 1:8){
      #if (i > j & ancientout[i,j] == 0) ancientout[i,j] <- NA
    }
  }
  
  
  ancientout <- melt(ancientout)
  ancientout$prop <- ancientout$value/sumrelationships * 100
  
  g1 <- ggplot(ancientout,mapping = aes(x=predicted_degrees  ,y=degreeN))+
    geom_raster(mapping = aes(fill=prop))+
    #facet_wrap(vars(Methods),nrow = 2,drop=FALSE)+
    #scale_fill_gradient2(low="#003366", high="#990033", mid="white")+
    scale_fill_gradient2(low="#003366", high="#DC0000B2", mid="gray99",midpoint = 30)+
    geom_text(mapping = aes(label=value),size=2,vjust=0.5)+
    #geom_text(mapping = aes(label=prop2),size=3,vjust=1)+
    labs(x="by clusIBD",y="by Fowler et al",fill="Percentage(%)")+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 0))
  #theme(axis.text.x = element_text(angle = 90),
  #      plot.margin = margin(t=0.1, b=0.3, l=1, r=0.4, "cm"),
  #      title = element_text(size = 10))
  
  #g10 g11 default   g20 g21 n=300 g30 g21 n =500
  figout[[2*index - 1]] <- g0
  figout[[2*index]] <- g1
}
legend1 <- get_legend(g0)
legend2 <- get_legend(g1)

g_ancient1 <- ggarrange(figout[[1]],figout[[3]],figout[[5]],
                        labels = c("A","C","E","D","E","F"),
                        font.label = list(size = 10),
                        ncol = 1,nrow = 3,
                        common.legend = T,legend = "left")
g_ancient2 <- ggarrange(figout[[2]],figout[[4]],figout[[6]],
                        labels = c("B","D","F","D","E","F"),
                        font.label = list(size = 10),
                        ncol = 1,nrow = 3,
                        common.legend = T,legend = "right")
g_ancient <- ggarrange(g_ancient1,g_ancient2)
#g_ancient <- ggarrange(g10,g20,g30,g11,g21,g31,labels = c("A","B","C","D","E","F"),common.legend = T,legend = "right")

ggsave(filename =  paste0("A:/课题/kinship/errorIBD/results/Figure 5 kinship inference for ancientDNA",Sys.Date(),".jpeg"),g_ancient,height=7,width =8)

1/(nrow(sample_pairs)-236)
sum(c(1,3,2))/(nrow(sample_pairs)-236)
sum(c(5,9,9,2,2,1))/(nrow(sample_pairs)-236)


unrelated_ancient <- c("NC10m","NE3m","NE4m","SC10f","SE4m","SE5m","SE6f","HN1f")

##无关
x[x$sample1 %in% unrelated_ancient & x$variable == "IBIS" & x$clusIBD > 0, ]

IBDfile <- "A:/课题/kinship/errorIBD/data/ancientDNA/clusIBD_ancientDNA_n500.IBD.details"

g001 <- clusIBD_show(s1="NE4m",s2="SC2m",IBDfile = IBDfile );print(g001)
g002 <- clusIBD_show(s1="NE4m",s2="SE1m",IBDfile = IBDfile );print(g002)
g003 <- clusIBD_show(s1="NE4m",s2="SP1m",IBDfile = IBDfile);print(g003)
g004 <- clusIBD_show(s1="NC10m",s2="SP3m", IBDfile = IBDfile);print(g004)
g005 <- clusIBD_show(s1="NC10m",s2="NC7f", IBDfile = IBDfile);print(g005)

g_segment <- marrangeGrob(list(g001,g002,g003,g004,g005),ncol = 1,nrow=1,top = quote(paste("Supplementary Figure", g+4)))
ggsave(filename =  paste0("A:/课题/kinship/errorIBD/results/Supplementary Figure 5-9 IBD segments for ancientDNA","",".pdf"),g_segment,height=7,width =4)





library(ggplot2)

rum_times <-  data.frame(sample_size=seq(100,500,by=100),
                         #clusIBD2=c(81,301,735,1204,2035),
                         #clusIBD=c(67,227,535,968,1467),
                         #clusIBD=c(33,115,267,484,733),
                         clusIBD=c(47,144,329,577,905),
                         IBIS=c(2,5,8,9,12),
                         TRUFFLE=c(4,13,25,44,64),
                         IBDseq=c(418,1127,2211,4069,5771)
)

rum_times <-  data.frame(sample_size=seq(100,500,by=100),
                         clusIBD=c(51,129,279,477,759),
                         IBIS=c(1,3,4,5,7),
                         TRUFFLE=c(3,9,16,27,41),
                         IBDseq=c(247,670,1418,2228,3678)
)

rum_times <- reshape2::melt(rum_times,measure.vars = c("clusIBD","IBIS","TRUFFLE","IBDseq"))

mycol2 = c( "#E64B35B2", "#4DBBD5B2" ,"#00A087B2", "#3C5488B2", "#F39B7FB2",
            "#8491B4B2", "#91D1C2B2", "#7E6148B2")
g_time= ggplot2::ggplot(rum_times,mapping = aes(x=sample_size,y=value,color=variable,fill=variable))+
  #geom_bar(stat = "identity",position = position_dodge())+
  geom_point(size=1.4)+
  geom_line(size=1.0)+
  scale_fill_manual(values =mycol2[c(1,2,5,4)])+
  scale_color_manual(values =mycol2[c(1,2,5,4)])+
  #scale_y_log10()+
  theme_bw()+
  theme(legend.position.inside = c(0.15,0.85),
        legend.position="inside",
        legend.background = element_blank())+
  labs(x="Sampe size",y="Running time (s)",color="",fill="")

ggsave(filename = paste0("A:/¿ÎÌâ/kinship/errorIBD/results/Figure 6 Running times ",Sys.Date(),".jpeg"),g_time,width = 5,height = 5)


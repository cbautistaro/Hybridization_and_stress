#You may need to install the following packages
rm(list=ls())
library(dplyr)
library(ggplot2)
library(tidyr)
library(data.table)
library(magrittr)
library(stringr)
library(UpSetR)
library(readr)
library(cowplot)
library(magrittr)
library(tidyverse)
library(readxl)
library(gridExtra)
library(growthcurver)
library(Cairo)
library(Hmisc)
library(nlme)
library(Matrix)
library(lme4)
library(ggpubr)
library(gtable)
library(grid.raster)



#Put the path to the folder where all the .csv are:
setwd("")



################Figure S2################
con <- read.csv("Fig_S2.csv")

table(con$note_fit, con$check)
con <- con %>% mutate(Final_perc_check=ifelse(check==TRUE, 0, Final_perc))
con <- con %>% mutate(Final_perc_check=ifelse(note_fit=="ERROR", 0, Final_perc_check))

################Figure S2A################
#FIGURE PANEL A 
#Scer in cm
scer <- filter(con, strain== "1Scer")%>%
  group_by(rep,strain) %>%
  ggplot(aes(x=as.numeric(NQO_c), y=Final_perc_check, group=interaction(rep, strain)))+ 
  #geom_errorbar(aes(ymin=len-sd, ymax=len+sd), width=.1, 
  #             position=position_dodge(0.05)) +
  geom_line(color="green4") +
  geom_point(color="green4") +
  theme(plot.title = element_text(size=22))+
  ylim(0,100)+
  labs(title="S. cerevisiae",x="UV mimetic (μM)", y = "Relative growth (%)")+
  theme_classic() + theme(plot.title=element_text(face = "italic")) +
  theme(axis.title.x = element_text(size=10), axis.text.x = element_text(size = 8),
        axis.title.y = element_text(size=10), axis.text.y = element_text(size = 8))+
  theme (plot.title = element_text(size = 8))

#Spar in cm
spar <- filter(con, strain== "2Spar")%>%
  group_by(rep,strain) %>%
  ggplot(aes(x=as.numeric(NQO_c), y=Final_perc_check, group=interaction(rep, strain)))+ 
  #geom_errorbar(aes(ymin=len-sd, ymax=len+sd), width=.1, 
  #             position=position_dodge(0.05)) +
  geom_line(color="dodgerblue1") +
  geom_point(color="dodgerblue1") +
  ylim(0,100)+
  labs(title="S. paradoxus",x="UV mimetic (μM)", y = "Relative growth (%)")+
  theme_classic() + theme(plot.title=element_text(face = "italic")) +
  theme(axis.title.x = element_text(size=10), axis.text.x = element_text(size = 8),
        axis.title.y = element_text(size=10), axis.text.y = element_text(size = 8))+
  theme (plot.title = element_text(size = 8))

#Hy in cm
hy <- filter(con, strain== "3Hybrid")%>%
  group_by(rep,strain) %>%
  ggplot(aes(x=as.numeric(NQO_c), y=Final_perc_check, group=interaction(rep, strain)))+ 
  #geom_errorbar(aes(ymin=len-sd, ymax=len+sd), width=.1, 
  #             position=position_dodge(0.05)) +
  geom_line(color="#FF9999") +
  geom_point(color="#FF9999") +
  ylim(0,100)+
  labs(title="Hybrid",x="UV mimetic (μM)", y = "Relative growth (%)")+
  theme_classic() +  theme(axis.title.x = element_text(size=10), axis.text.x = element_text(size = 8),
                           axis.title.y = element_text(size=10), axis.text.y = element_text(size = 8))+
  theme (plot.title = element_text(size = 8))

figA <- grid.arrange(scer, spar, hy, nrow=3 )

################Figure S2B################
gr <- filter(gr, strain!="White")

#Remove the first p-value
toexpr<-function(x) {
  getfun <- function(x) {
    ifelse(x=="Hybrid", "plain", "italic")
  }
  as.expression(unname(Map(function(f,v) substitute(f(v), list(f=as.name(f), v=as.character(v))), getfun(x), x)))
}
con1  <-filter(con, NQO_c != 0)

box <- con %>%
  ggplot(aes(interaction(as.numeric(NQO_c)), Final_perc_check, col=strain)) + 
  geom_boxplot(outlier.shape = NA)+
  geom_point(position=position_jitterdodge(), size = 0.3) +
  theme_classic()+
  theme(axis.text.x = element_text(hjust = 1)) +
  xlab("UV mimetic (μM)") + ylab("Relative growth (%)") +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.position = "top", legend.text.align = 0,
        strip.background = element_rect(fill = "white", color="white"), 
        strip.text.y = element_blank()) + scale_color_manual(values=c("green4", "dodgerblue1", "#FF9999"), 
                                                             labels = toexpr(c("S. cerevisiae", "S. paradoxus", "Hybrid"))) +
  guides(color = guide_legend(title = "Genotype")) +
  theme(axis.title.x = element_text(size=10), axis.text.x = element_text(size = 8),
        axis.title.y = element_text(size=10), axis.text.y = element_text(size = 8)) +
  stat_compare_means(aes(label = ..p.format..), size=2, data = con %>% filter(NQO_c != 0)) + ylim(0,100) +
  theme(legend.text=element_text(size=7), legend.title=element_text(size=8)) 


box1 <- box + guides(
  color = guide_legend(
    title = "Genotype",
    override.aes = aes(label = "")
  )
)

################Save Figure S2################
figA <- plot_grid(scer,spar,hy, labels = c("A","",""), nrow = 3)
box1 <- plot_grid(box1, labels = c("B"), nrow = 1)

figS2 <- plot_grid(figA, box1, ncol=2 )

ggsave (plot = figS2, filename = "figureS2.jpg", units = "cm", device = "jpg",width = 20, height = 15, dpi = 300)
ggsave (plot = figS2, filename = "figureS2.pdf", units = "cm", device = "pdf",width = 20, height = 15, dpi = 300)


















################FigureS3################
################Get Dataset################
fdata <- read.csv("FigS3A.csv")

fdata <- na.omit(fdata)
fdata$day <- as.numeric(fdata$day)
fdata$gen = (fdata$day)*10

################Figure figuring out if variable OD at the end of each transfer################
fdata1 <- fdata %>% 
  dplyr::filter(strain!="White")  %>%
  group_by(day, strain, NQO, rep, plate, tecan, gen) %>% 
  dplyr::summarise(od = max(od),max_temp = max(time)) 

toexpr<-function(x) {
  getfun <- function(x) {
    ifelse(x=="Hybrid", "plain", "italic")
  }
  as.expression(unname(Map(function(f,v) substitute(f(v), list(f=as.name(f), v=as.character(v))), getfun(x), x)))
}

fdata1 <- filter(fdata1, NQO=="Yes")

#Fig S3A
fdata1 <- filter(fdata1, gen<220)
A <- ggplot(fdata1, aes(interaction(as.numeric(gen)), od, colour=strain)) +
  geom_boxplot() +
  geom_point(position=position_jitterdodge(), size=0.7) +
  scale_shape_manual(values=c(16, 1)) + theme_bw() +
  theme(legend.position = "top", axis.title.x = element_text(size=12), axis.text.x = element_text(size = 10),
        axis.title.y = element_text(size=12), axis.text.y = element_text(size = 10)) + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                                                                             panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Time (generations)") + ylab("OD (600 nm)") + scale_color_manual(values=c("green4", "dodgerblue1", "#FF9999"), 
                                                                        labels = toexpr(c("S. cerevisiae", "S. paradoxus", "Hybrid"))) +
  theme(legend.text.align = 0) + guides(color = guide_legend(title = "Genotype", linetype = guide_legend(title = "UV mimetic (μM)"), override.aes=list(fill=NA)), 
                                        shape = guide_legend(title = "UV mimetic"), override.aes=list(fill=NA)) + 
  theme(legend.text=element_text(size=6), legend.title=element_text(size=7), 
        legend.direction = "horizontal", legend.box = "vertical") +
  theme(legend.spacing.y = unit(-0.2, "cm"),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(0,0,-10,-10))

################Panel FigureS3B################
data <- read.csv("FigS3B.csv")
#I need to filter these data to choose the days:
#24 H CONTROL
control24 <- filter(data, Timing==24)
control24 <- filter(control24 , NQO_c==0)

#Graphic 4nqo 16 uM
nqo48 <- filter(data, Timing==48)
nqo48 <- filter(nqo48 , NQO_c==16)

nqo96 <- filter(data, Timing==96)
nqo96 <- filter(nqo96 , NQO_c==16)


fdata <- rbind(nqo96, nqo48)
fdata$Time <- as.numeric(fdata$Time)
fdata$hour = (fdata$Time)/3600
fdata <- filter(fdata, Strain!="White")
fdata <- fdata %>% mutate(strain = ifelse(Strain=="Parent", "2Spar", ifelse(Strain=="Diploid", "1Scer", ifelse(Strain=="Hybrid", "3Hybrid", FALSE))))

toexpr<-function(x) {
  getfun <- function(x) {
    ifelse(x=="Hybrid", "plain", "italic")
  }
  as.expression(unname(Map(function(f,v) substitute(f(v), list(f=as.name(f), v=as.character(v))), getfun(x), x)))
}

#Fig S3B in cm
B <- ggplot(fdata, aes(x=hour, y=DO, colour=strain, group=interaction(Replica,strain,Timing))) +
  #stat_summary(position="identity" ,geom="line", fun.y = mean, size=0.1) +
  #stat_summary(fun.data = mean_sdl, geom = "errorbar") +
  geom_line(aes(linetype = as.factor(Timing))) +
  xlab("Time (h)") + ylab("OD (600 nm)") +theme_minimal()+
  theme(axis.text.x = element_text(hjust = 1)) +
  xlim(0,30)+
  facet_grid(. ~ strain) + 
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),  legend.position = "top",
                     axis.line = element_line(colour = "black"), 
                     strip.background = element_rect(fill = "white", color="white"), 
                     strip.text.x = element_blank()) + scale_color_manual(values=c("green4", "dodgerblue1", "#FF9999"), 
                                                                          labels = toexpr(c("S. cerevisiae", "S. paradoxus", "Hybrid"))) +
  scale_linetype_manual(values = c(1, 96), labels = c("10","20"))+
  theme(axis.title.x = element_text(size=12), axis.text.x = element_text(size = 10),
        axis.title.y = element_text(size=12), axis.text.y = element_text(size = 10))+
  theme(legend.text.align = 0) + guides(color = guide_legend(title = "Genotype"),
                                        linetype = guide_legend(title = "Generation")) +
  theme(legend.text=element_text(size=6), legend.title=element_text(size=7), 
        legend.direction = "horizontal", legend.box = "vertical") +
  theme(legend.spacing.y = unit(-0.2, "cm"),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(0,0,-4,-4))

################Save Figure################
FIGS3 <- plot_grid(A,B, labels = c("A","B"), nrow = 2)
FIGS3

ggsave (plot = FIGS3, filename = "figureS3.jpg", units = "cm", device = "jpg",width = 20, height = 15, dpi = 300)
ggsave (plot = FIGS3, filename = "figureS3.pdf", units = "cm", device = "pdf",width = 20, height = 15, dpi = 300)




################Fig S4################
################Fig S4A################
fdata <- read.csv("Fig_2.csv")

fdata <- na.omit(fdata)
fdata$day <- as.numeric(fdata$day)

#remove days we do not have all time-points of the curves
fdata <- dplyr::filter(fdata, day!=1)
fdata <- dplyr::filter(fdata, day!=2)
fdata <- dplyr::filter(fdata, day!=11)
fdata <- dplyr::filter(fdata, day!=22)
fdata <- dplyr::filter(fdata, day!=23)

fdata$hour = (fdata$time)/3600


growth_rates <- fdata %>% group_by(strain, NQO, rep, day) %>%
  dplyr::mutate(note_fit = SummarizeGrowth(hour, od)$vals$note,
                rval = SummarizeGrowth(hour, od)$vals$r,
                kval = SummarizeGrowth(hour, od)$vals$k,
                tmid = SummarizeGrowth(hour, od)$vals$t_mid,
                aucexp = SummarizeGrowth(hour, od)$vals$auc_l,
                max_size = max(od)) %>%
  mutate(note_fit = ifelse(note_fit=="", "OK", "ERROR")) %>%
  slice(1) %>%
  ungroup() %>%
  select(strain, NQO, rep, day, time, note_fit,rval, kval,tmid,aucexp, max_size, od, plate)

gr <- dplyr::filter(growth_rates, note_fit=="OK")
gr <- dplyr::filter(growth_rates, strain!="White")
fdataAREA <- gr

################Panel S4A################
#Figure in control conditions during exp evolution
fig0 <- filter(fdataAREA, day == "3")
fig0 <- filter(fig0, NQO == "No")
fig0a <- filter(fdataAREA, day == "21")
fig0a <- filter(fig0a, NQO == "No")

fig0 <- rbind(fig0, fig0a)

my_comparisons2 <- list( c("1Scer.3", "1Scer.21"),c("2Spar.3", "2Spar.21"),c("3Hybrid.3", "3Hybrid.21"))

fig0$strain_day <- paste(fig0$strain,".", fig0$day)

#Fig S4 in cm
fig0$strain_day <- paste(fig0$strain,".", fig0$day)
FigS4 <- fig0 %>% mutate(strain_day = as.character(interaction (fig0$strain, fig0$day))) %>%
  mutate(strain_day=factor(strain_day, levels= c("1Scer.3", "1Scer.21", "2Spar.3", "2Spar.21", "3Hybrid.3", "3Hybrid.21"))) %>%
  ggplot(., aes(x=as.factor(strain_day), y=rval, col=strain)) + 
  geom_boxplot() +
  theme_minimal()+
  xlab("Generation") + ylab("Growth rate (OD/hour)")+
  ggtitle("Growth in YPD during experimental evolution")+
  ylim(0,1.1) +
  scale_x_discrete(breaks=c("1Scer.3","1Scer.21","2Spar.3","2Spar.21","3Hybrid.3","3Hybrid.21"), labels=c("T30","T200","T30","T200","T30","T200"))+
  theme(legend.position = "top",axis.title.x = element_text(size=12), axis.text.x = element_text(size = 10),
        axis.title.y = element_text(size=12), axis.text.y = element_text(size = 10)) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.ticks = element_line(colour="black"))+
  stat_compare_means(aes(group = strain_day), comparisons= my_comparisons2, method = "t.test", label.y = 0.85, paired = TRUE, label="p.format", size = 2.5) +
  scale_color_manual(values=c("green4", "dodgerblue1", "#FF9999"), 
                     labels = toexpr(c("S. cerevisiae", "S. paradoxus", "Hybrid"))) +
  theme(legend.text.align = 0) + guides(color = guide_legend(title = "Genotype"))+
  geom_jitter(position=position_jitter(0.2), size=0.4)+
  theme(plot.title = element_text(size = 9, hjust = 0.5))+
  theme(legend.text=element_text(size=7), legend.title=element_text(size=8))

################Fig S4B################
fdata <- read.csv("Fig_3.csv")

#Remove NA data
fdata <- na.omit(fdata)

fdata$hour = (fdata$time)/3600

growth_rates <- fdata %>% group_by(strain, NQO, rep, day, gen) %>%
  mutate(note_fit = SummarizeGrowth(hour, od)$vals$note,
         rval = SummarizeGrowth(hour, od)$vals$r,
         kval = SummarizeGrowth(hour, od)$vals$k,
         tmid = SummarizeGrowth(hour, od)$vals$t_mid,
         aucexp = SummarizeGrowth(hour, od)$vals$auc_l,
         max_size = max(od)) %>%
  mutate(note_fit = ifelse(note_fit=="", "OK", "ERROR")) %>%
  slice(1) %>%
  ungroup() %>%
  select(strain, NQO, rep, day, time, note_fit,rval, kval,tmid,aucexp, max_size, od, gen)


gr <- filter(growth_rates, note_fit=="OK")

#To know the replicates
Scer <- filter(gr, strain=="1Scer")
Scer <- filter(Scer, NQO=="No")
Scer <- filter(Scer, gen=="0")

Sp <- filter(gr, strain=="2Spar")
Sp <- filter(Sp, NQO=="No")
Sp <- filter(Sp, gen=="0")

H <- filter(gr, strain=="3Hybrid")
H <- filter(H, NQO=="No")
H <- filter(H, gen=="0")

################Panel S4B###############
gr_1 <- filter(gr, rep!="31")
gr_1 <- filter(gr_1, rep!="32")

my_comparisons1 <- list( c("0.2Spar", "200.2Spar"),
                         c("0.3Hybrid", "200.3Hybrid"),
                         c("0.1Scer", "200.1Scer"))

toexpr<-function(x) {
  getfun <- function(x) {
    ifelse(x=="Hybrid", "plain", "italic")
  }
  as.expression(unname(Map(function(f,v) substitute(f(v), list(f=as.name(f), v=as.character(v))), getfun(x), x)))
}

gr_1 <- filter(gr_1 , NQO=="No")
FigS4B <- gr_1 %>% ggplot(aes(interaction(gen, strain),rval, col=strain)) + 
  geom_boxplot() + theme_classic()+
  ylim(0,1.1)+
  scale_x_discrete(breaks=c("0.1Scer", "200.1Scer","0.2Spar", "200.2Spar","0.3Hybrid", "200.3Hybrid"),labels=c("T0","T200","T0","T200","T0","T200"))+
  ggtitle("Growth in YPD after experimental evolution")+
  stat_compare_means(comparisons = my_comparisons1, method="t.test", label.y = c(0.9,0.9,0.9), size=2.5, paired = TRUE) +
  theme(legend.position = "top", panel.border = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Generation") + ylab("Growth rate (OD/hour)") + scale_color_manual(values=c("green4", "dodgerblue1", "#FF9999"), 
                                                                          labels = toexpr(c("S. cerevisiae", "S. paradoxus", "Hybrid"))) +
  theme(legend.text.align = 0) + guides(color = guide_legend(title = "Genotype")) +
  theme(axis.title.x = element_text(size=12), axis.text.x = element_text(size = 10),
        axis.title.y = element_text(size=12), axis.text.y = element_text(size = 10))+
  geom_jitter(position=position_jitter(0.2), size=0.4) +
  theme(legend.text=element_text(size=12))+
  theme(plot.title = element_text(size = 9, hjust = 0.5))+
  theme(legend.text=element_text(size=7), legend.title=element_text(size=8))


################Save Figure################
figS4 <- plot_grid(FigS4,FigS4B, labels = c("A","B"), nrow = 1)
figS4
ggsave (plot = figS4, filename = "figureS4.jpg", units = "cm", device = "jpg",width = 20, height = 10, dpi = 300)
ggsave (plot = figS4, filename = "figureS4.pdf", units = "cm", device = "pdf",width = 20, height = 10, dpi = 300)
















################Figure S5################
con <- read.csv("FigS5.csv")

#Fig Scer in cm
scer <- filter(con, strain== "1Scer")%>%
  group_by(rep,evol,strain,curve) %>%
  ggplot(aes(x=as.numeric(NQO_c), y=Final_perc_check, shape=evol, group=interaction(rep, evol,strain,curve), linetype=evol)) + 
  #geom_errorbar(aes(ymin=len-sd, ymax=len+sd), width=.1, 
  #             position=position_dodge(0.05)) +
  geom_line(color="green4") +
  geom_point(color="green4") +
  ylim(0,100)+
  scale_shape_manual(values=c(20, 1)) +
  scale_linetype_manual(values=c("solid", "dashed")) +
  labs(title="S. cerevisiae",x="UV mimetic (μM)", y = "Relative growth (%)")+
  theme_classic() + theme(plot.title=element_text(face = "italic")) +
  theme(axis.title.x = element_text(size=14), axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size=14), axis.text.y = element_text(size = 12)) +
  guides(linetype = guide_legend(title = "Evolved in:"), shape = FALSE)

#Fig Spar in cm
spar <- filter(con, strain== "2Spar")%>%
  group_by(rep,evol,strain,curve) %>%
  ggplot(aes(x=as.numeric(NQO_c), y=Final_perc_check, shape=evol, group=interaction(rep, evol,strain,curve), linetype=evol))+ 
  #geom_errorbar(aes(ymin=len-sd, ymax=len+sd), width=.1, 
  #             position=position_dodge(0.05)) +
  geom_line(color="dodgerblue1") +
  geom_point(color="dodgerblue1") +
  ylim(0,100)+
  scale_shape_manual(values=c(20, 1)) +
  scale_linetype_manual(values=c("solid", "dashed")) +
  labs(title="S. paradoxus",x="UV mimetic (μM)", y = "Relative growth (%)")+
  theme_classic() + theme(plot.title=element_text(face = "italic")) +
  theme(axis.title.x = element_text(size=14), axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size=14), axis.text.y = element_text(size = 12)) +
  guides(linetype = guide_legend(title = "Evolved in:"), shape = FALSE)

#Hybrid in cm
hy <- filter(con, strain== "3Hybrid")%>%
  group_by(rep,evol,strain,curve) %>%
  ggplot(aes(x=as.numeric(NQO_c), y=Final_perc_check, shape=evol, group=interaction(rep, evol,strain,curve), linetype=evol))+ 
  #geom_errorbar(aes(ymin=len-sd, ymax=len+sd), width=.1, 
  #             position=position_dodge(0.05)) +
  geom_line(color="#FF9999") +
  geom_point(color="#FF9999") +
  ylim(0,100)+
  scale_shape_manual(values=c(20, 1)) +
  scale_linetype_manual(values=c("solid", "dashed")) +
  labs(title="Hybrid",x="UV mimetic (μM)", y = "Relative growth (%)")+
  theme_classic()  +
  theme(axis.title.x = element_text(size=14), axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size=14), axis.text.y = element_text(size = 12)) +
  guides(linetype = guide_legend(title = "Evolved in:"), shape = FALSE)


##############Save figure S5
figS5 <- plot_grid(scer,spar,hy, labels = c("A","B","C"), nrow = 3)

ggsave (plot = figS5, filename = "figureS5.jpg", units = "cm", device = "jpg",width = 15, height = 20, dpi = 300)
ggsave (plot = figS5, filename = "figureS5.pdf", units = "cm", device = "pdf",width = 15, height = 20, dpi = 300)

  
  




################Figure S6################
Ploidy_2 <- read.csv("Fig_4.csv")

Ploidy_2$gen <- as.factor(Ploidy_2$gen)
Ploidy_2$rep <- as.factor(Ploidy_2$rep)  
Ploidy_2 <- dplyr::filter(Ploidy_2, rep!= "31")
Ploidy_2 <- dplyr::filter(Ploidy_2, rep!= "32")
Ploidy_2 <- dplyr::filter(Ploidy_2, Ploidy!= "NA")

Ploidy_2$gen = factor(Ploidy_2$gen, levels=c("200_control", "0", "200_UV"))
pd <- position_dodge(0.4)

toexpr<-function(x) {
  getfun <- function(x) {
    ifelse(x=="Hybrid", "plain", "italic")
  }
  as.expression(unname(Map(function(f,v) substitute(f(v), list(f=as.name(f), v=as.character(v))), getfun(x), x)))
}

Ploidy_2 <- Ploidy_2 %>% mutate(colusd= Ploidy + rnorm( n = nrow(Ploidy_2), mean = 0, sd = 0.05))

################Panel S6A################
#Control
NOP <- filter(Ploidy_2, NQO=="No")
no <-ggplot(NOP, aes(as.factor(gen), as.numeric(colusd), color=strain, group=interaction(strain, rep))) + 
  geom_line() +
  geom_point() +
  scale_shape_manual(values=c(16, 8)) + theme_bw() +
  facet_grid(. ~ strain) + theme(plot.title = element_text(hjust = 0.5, face = 'bold'),
                                 panel.background = element_rect(fill = "white"),
                                 axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16),
                                 axis.line = element_line(size = 0.5),
                                 axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 14),
                                 strip.text.x = element_text(size = 15, face = 'bold'),
                                 strip.background = element_rect(fill = 'white')) +
  xlab("Time (generations)") + ylab("Ploidy (n)") + 
  scale_x_discrete(limit = c("0", "200_control"),
                   label= c("Ancestors", "Evolved in control"), 
                   breaks=c("0", "200_control")) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"), 
        strip.background = element_rect(fill = "white", color="white"), 
        strip.text.y = element_blank()) + 
  scale_color_manual(values=c("green4", "dodgerblue1", "#FF9999"),
                     labels = toexpr(c("S. cerevisiae", "S. paradoxus", "Hybrid"))) +
  theme(legend.text.align = 0) + 
  guides(color = guide_legend(title = "Genotype", override.aes=list(fill=NA)),
         shape = FALSE) +
  theme(strip.text.x = element_text(size=0)) +
  scale_y_discrete(limit = c(2, 3),
                   breaks=c(2, 3),
                   label = c("2", "3"))


################Panel S6B################
YES <- filter(Ploidy_2, gen=="200_UV")
YES1 <- filter(Ploidy_2, gen=="0")
YES <- rbind(YES,YES1)
yes <- ggplot(YES, aes(as.factor(gen), as.numeric(colusd), color=strain, group=interaction(strain, rep))) + 
  geom_line() +
  geom_point() +
  scale_shape_manual(values=c(16, 8)) + theme_bw() +
  facet_grid(. ~ strain) + theme(plot.title = element_text(hjust = 0.5, face = 'bold'),
                                 panel.background = element_rect(fill = "white"),
                                 axis.title.x = element_text(size = 16), axis.title.y = element_text(size = 16),
                                 axis.line = element_line(size = 0.5),
                                 axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 14),
                                 strip.text.x = element_text(size = 15, face = 'bold'),
                                 strip.background = element_rect(fill = 'white')) +
  xlab("Time (generations)") + ylab("Ploidy (n)") + 
  scale_x_discrete(limit = c("0", "200_UV"),
                   label= c("Ancestors", "Evolved in UV mimetic"), 
                   breaks=c("0", "200_UV")) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"), 
        strip.background = element_rect(fill = "white", color="white"), 
        strip.text.y = element_blank()) + 
  scale_color_manual(values=c("green4", "dodgerblue1", "#FF9999"),
                     labels = toexpr(c("S. cerevisiae", "S. paradoxus", "Hybrid"))) +
  theme(legend.text.align = 0) + 
  guides(color = guide_legend(title = "Genotype", override.aes=list(fill=NA)),
         shape = FALSE) +
  theme(strip.text.x = element_text(size=0)) +
  scale_y_discrete(limit = c(2, 3),
                   breaks=c(2, 3),
                   label = c("2", "3"))


################Save S6################
figS6 <- plot_grid(no, yes, labels = c("A", "B"), nrow = 2)
figS6

ggsave (plot = figS6, filename = "figureS6.jpg", units = "cm", device = "jpg",width = 20, height = 15, dpi = 300)
ggsave (plot = figS6, filename = "figureS6.pdf", units = "cm", device = "pdf",width = 20, height = 15, dpi = 300)

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
library(readxl)
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
library(drc)
library(EDcomp)
library(reshape2)
library(broom)
library(modelr)
library(egg)
library(ggpubr)
library(rstatix)



#Put the path to the folder where all the .csv are:
setwd("")


################Figure S2################
#Number of generations
Generations <- read.csv("Fig_S2.csv")
#Generations Figure
toexpr<-function(x) {
  getfun <- function(x) {
    ifelse(x=="Hybrid", "plain", "italic")
  }
  as.expression(unname(Map(function(f,v) substitute(f(v), list(f=as.name(f), v=as.character(v))), getfun(x), x)))
}

stat.test <- aov(Generations ~ Specie, data = Generations) %>%
  tukey_hsd()

stat_pvalue_manual(
  stat.test, label = "p.adj", 
  y.position = c(0.7, 0.72, 0.74),
  size=2.5
)

#in cm
FigS2 <- Generations %>% ggplot(., aes(x=as.factor(Specie), y=Generations, col=Specie)) + 
  geom_boxplot(outlier.shape = NA) + ylim(4.5,5.5) +
  theme_minimal()+
  xlab("Genotype") + ylab("Generations")+
  theme(axis.title.x = element_text(size=12), axis.text.x = element_text(size = 8),
        axis.title.y = element_text(size=12), axis.text.y = element_text(size = 10)) +
  theme(legend.position="none",panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.ticks = element_line(colour="black"))+
  stat_compare_means(method = "anova", label.y = 5.3, aes(label = ..p.format..), size=2.5) +
  stat_pvalue_manual(
    stat.test, label = "p.adj", 
    y.position = c(5.13, 5.17, 5.21),
    size=2.5
  ) +
  scale_color_manual(values=c("green4", "dodgerblue1", "#FF9999"), 
                     labels = toexpr(c("S. cerevisiae", "S. paradoxus", "Hybrid"))) +
  theme(legend.text.align = 0) + guides(color = guide_legend(title = "Genotype"))+
  geom_jitter(position=position_jitter(0.2), size=0.6) +
  theme(legend.text=element_text(size=12)) +
  scale_x_discrete("Genotype", labels=expression(italic(S.cerevisiae), italic(S.paradoxus), Hybrid)) +
  theme(legend.text=element_text(size=7), legend.title=element_text(size=8))

FigS2 <- FigS2 + guides(
  color = guide_legend(
    title = "Genotype",
    override.aes = aes(label = "")
  )
)


ggsave (plot = FigS2, filename = "FigS2.jpg", units = "cm", device = "jpg",width = 10, height = 10, dpi = 300)
ggsave (plot = FigS2, filename = "FigS2.pdf", units = "cm", device = "pdf",width = 10, height = 10, dpi = 300)



################Figure S3################
################Get Dataset################
fdata <- read.csv("Fig_S3A.csv")

fdata <- na.omit(fdata)
fdata$day <- as.numeric(fdata$day)
fdata$gen = (fdata$day)*5
fdata <- filter(fdata, day<22)

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

A <- ggplot(fdata1, aes(interaction(as.numeric(gen)), od, colour=strain)) +
  geom_boxplot() +
  geom_point(position=position_jitterdodge(), size=0.7) +
  scale_shape_manual(values=c(16, 1)) + theme_bw() +
  theme(legend.position = "top", axis.title.x = element_text(size=12), axis.text.x = element_text(size = 10),
        axis.title.y = element_text(size=12), axis.text.y = element_text(size = 10)) + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                                                                             panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Time (generations)") + ylab("OD (595 nm)") + scale_color_manual(values=c("green4", "dodgerblue1", "#FF9999"), 
                                                                        labels = toexpr(c("S. cerevisiae", "S. paradoxus", "Hybrid"))) +
  theme(legend.text.align = 0) + guides(color = guide_legend(title = "Genotype", linetype = guide_legend(title = "UV mimetic (μM)"), override.aes=list(fill=NA)), 
                                        shape = guide_legend(title = "UV mimetic"), override.aes=list(fill=NA)) + 
  theme(legend.text=element_text(size=6), legend.title=element_text(size=7), 
        legend.direction = "horizontal", legend.box = "vertical") +
  theme(legend.spacing.y = unit(-0.2, "cm"),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(0,0,-10,-10))


################Panel FigureS3B################
data <- read.csv("Fig_S3B.csv")
#Need to filter these data to choose the days:
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
D <- ggplot(fdata, aes(x=hour, y=DO, colour=strain, group=interaction(Replica,strain,Timing))) +
  geom_line(aes(linetype = as.factor(Timing))) +
  xlab("Time (h)") + ylab("OD (595 nm)") +theme_minimal()+
  theme(axis.text.x = element_text(hjust = 1)) +
  xlim(0,30)+
  facet_grid(. ~ strain) + 
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),  legend.position = "top",
                     axis.line = element_line(colour = "black"), 
                     strip.background = element_rect(fill = "white", color="white"), 
                     strip.text.x = element_blank()) + scale_color_manual(values=c("green4", "dodgerblue1", "#FF9999"), 
                                                                          labels = toexpr(c("S. cerevisiae", "S. paradoxus", "Hybrid"))) +
  scale_linetype_manual(values = c(1, 96), labels = c("6","12"))+
  theme(axis.title.x = element_text(size=12), axis.text.x = element_text(size = 10),
        axis.title.y = element_text(size=12), axis.text.y = element_text(size = 10))+
  theme(legend.text.align = 0) + guides(color = guide_legend(title = "Genotype"),
                                        linetype = guide_legend(title = "Generation")) +
  theme(legend.text=element_text(size=6), legend.title=element_text(size=7), 
        legend.direction = "horizontal", legend.box = "vertical") +
  theme(legend.spacing.y = unit(-0.2, "cm"),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(0,0,-4,-4))

################Save Figure S3################
FIGS3 <- plot_grid(A,D, labels = c("A","B"), nrow = 2)


ggsave (plot = FIGS3, filename = "FigS3.jpg", units = "cm", device = "jpg",width = 20, height = 15, dpi = 300)
ggsave (plot = FIGS3, filename = "FigS3.pdf", units = "cm", device = "pdf",width = 20, height = 15, dpi = 300)


################Figure S4################
con <- read.csv("Fig_S4.csv")

table(con$note_fit, con$check)
con <- con %>% mutate(Final_perc_check=ifelse(check==TRUE, 0, Final_perc))
con <- con %>% mutate(Final_perc_check=ifelse(note_fit=="ERROR", 0, Final_perc_check))

################Figure S4A################
#FIGURE PANEL A 
#Scer in cm
scer <- filter(con, strain== "1Scer")%>%
  group_by(rep,strain) %>%
  ggplot(aes(x=as.numeric(NQO_c), y=Final_perc_check, group=interaction(rep, strain)))+ 
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
  geom_line(color="#FF9999") +
  geom_point(color="#FF9999") +
  ylim(0,100)+
  labs(title="Hybrid",x="UV mimetic (μM)", y = "Relative growth (%)")+
  theme_classic() +  theme(axis.title.x = element_text(size=10), axis.text.x = element_text(size = 8),
                           axis.title.y = element_text(size=10), axis.text.y = element_text(size = 8))+
  theme (plot.title = element_text(size = 8))

figA <- grid.arrange(scer, spar, hy, nrow=3 )

################Figure S4B################
con <- filter(con, strain!="White")

#Remove the first p-value
toexpr<-function(x) {
  getfun <- function(x) {
    ifelse(x=="Hybrid", "plain", "italic")
  }
  as.expression(unname(Map(function(f,v) substitute(f(v), list(f=as.name(f), v=as.character(v))), getfun(x), x)))
}

box <- con %>%
  ggplot(aes(interaction(as.numeric(NQO_c)), Final_perc_check, col=strain)) + 
  geom_boxplot(outlier.shape = NA)+
  geom_point(position=position_jitterdodge(), size = 0.3) +
  theme_classic() +
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


################Save Figure S4################
figA <- plot_grid(scer,spar,hy, labels = c("A","",""), nrow = 3)
box1 <- plot_grid(box1, labels = c("B"), nrow = 1)

figS4 <- plot_grid(figA, box1, ncol=2 )

ggsave (plot = figS4, filename = "FigS4.jpg", units = "cm", device = "jpg",width = 20, height = 15, dpi = 300)
ggsave (plot = figS4, filename = "FigS4.pdf", units = "cm", device = "pdf",width = 20, height = 15, dpi = 300)




################Figure S5################
################Figure S5################
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
  ungroup() 

gr <- dplyr::filter(growth_rates, note_fit=="OK")
gr <- dplyr::filter(growth_rates, strain!="White")
fdataAREA <- gr



################Panel S5A################
#Figure in control conditions during exp evolution
fig0 <- filter(fdataAREA, day == "3")
fig0 <- filter(fig0, NQO == "No")
fig0a <- filter(fdataAREA, day == "21")
fig0a <- filter(fig0a, NQO == "No")

fig0 <- rbind(fig0, fig0a)

my_comparisons2 <- list( c("1Scer.3", "1Scer.21"),c("2Spar.3", "2Spar.21"),c("3Hybrid.3", "3Hybrid.21"))

fig0$strain_day <- paste(fig0$strain,".", fig0$day)

#Fig S5 in cm
fig0$strain_day <- paste(fig0$strain,".", fig0$day)
FigS5 <- fig0 %>% mutate(strain_day = as.character(interaction (fig0$strain, fig0$day))) %>%
  mutate(strain_day=factor(strain_day, levels= c("1Scer.3", "1Scer.21", "2Spar.3", "2Spar.21", "3Hybrid.3", "3Hybrid.21"))) %>%
  ggplot(., aes(x=as.factor(strain_day), y=rval, col=strain)) + 
  geom_boxplot(outlier.shape = NA) +
  theme_minimal()+
  xlab("Generation") + ylab("Growth rate (OD/hour)")+
  ggtitle("Growth in YPD during experimental evolution")+
  ylim(0,1.1) +
  scale_x_discrete(breaks=c("1Scer.3","1Scer.21","2Spar.3","2Spar.21","3Hybrid.3","3Hybrid.21"), labels=c("T15","T100","T15","T100","T15","T100"))+
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

################Fig S5B################
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
  ungroup() 

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

################Panel S5B###############
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

FigS5B <- gr_1 %>% ggplot(aes(interaction(gen, strain),rval, col=strain)) + 
  geom_boxplot(outlier.shape = NA) + theme_classic()+
  ylim(0,1.1)+
  scale_x_discrete(breaks=c("0.1Scer", "200.1Scer","0.2Spar", "200.2Spar","0.3Hybrid", "200.3Hybrid"),labels=c("Ancestor","T100","Ancestors","T100","Ancestors","T100"))+
  ggtitle("Growth in YPD after experimental evolution in UV mimetic")+
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


################Save Figure S5################
figS5 <- plot_grid(FigS5,FigS5B, labels = c("A","B"), nrow = 1)
figS5
ggsave (plot = figS5, filename = "FigS5.jpg", units = "cm", device = "jpg",width = 20, height = 10, dpi = 300)
ggsave (plot = figS5, filename = "FigS5.pdf", units = "cm", device = "pdf",width = 20, height = 10, dpi = 300)



################Figure S6################
con <- read.csv("Fig_S6.csv")

#Fig Scer in cm
scer <- filter(con, strain== "1Scer")%>%
  group_by(rep,evol,strain,curve) %>%
  ggplot(aes(x=as.numeric(NQO_c), y=Final_perc_check, shape=evol, group=interaction(rep, evol,strain,curve), linetype=evol)) + 
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


##############Save Figure S6##############
figS6 <- plot_grid(scer,spar,hy, labels = c("A","B","C"), nrow = 3)

ggsave (plot = figS6, filename = "FigS6.jpg", units = "cm", device = "jpg",width = 15, height = 20, dpi = 300)
ggsave (plot = figS6, filename = "FigS6.pdf", units = "cm", device = "pdf",width = 15, height = 20, dpi = 300)



################Figure S7################
#Graph correlation coefficients through time by replicate
correlation_graph <- read.csv("Fig_S7.csv")
correlation_graph <- correlation_graph %>% mutate(gen= day*5)

toexpr<-function(x) {
  getfun <- function(x) {
    ifelse(x=="Hybrid", "plain", "italic")
  }
  as.expression(unname(Map(function(f,v) substitute(f(v), list(f=as.name(f), v=as.character(v))), getfun(x), x)))
}


#Fig S7 in cm
FigS7 <- correlation_graph %>%
  ggplot(., aes(gen, corr, colour=strain, shape=NQO, linetype=NQO))  + 
  geom_point() +
  scale_linetype_manual(values=c("solid","dotted")) +
  scale_shape_manual(values=c(16, 1))+
  theme_classic()+
  theme(legend.position = "top",panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  ylim(-1, 1) +
  scale_color_manual(values=c("green4", "dodgerblue1", "#FF9999"), 
                     labels = toexpr(c("S. cerevisiae", "S. paradoxus", "Hybrid"))) +
  theme(legend.text.align = 0) + guides(color = guide_legend(title = "Genotype"), shape = guide_legend(title = "UV mimetic"), linetype=FALSE) +
  stat_smooth(method="lm", se=FALSE) + labs(y="r/K correlation", x = "Time (generations)") +
  theme(axis.title.x = element_text(size=12), axis.text.x = element_text(size = 10),
        axis.title.y = element_text(size=12), axis.text.y = element_text(size = 10)) +
  guides(color = guide_legend(title = "Genotype",
                              override.aes=list(fill=NA)), 
         linetype = guide_legend(title = "UV mimetic", 
                                 override.aes=list(fill=NA)),
         shape = guide_legend(title = "UV mimetic", override.aes=list(colour="black"))) +
  theme(legend.text=element_text(size=5), legend.title=element_text(size=6), 
        legend.direction = "horizontal", legend.box = "vertical") +
  theme(legend.spacing.y = unit(-0.2, "cm"),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(0,0,-4,-4))  


FigS7 <- FigS7 + geom_hline(yintercept=0, linetype="dashed", color = "black")

FigS7

ggsave (plot = FigS7, filename = "FigS7.jpg", units = "cm", device = "jpg",width = 15, height = 10, dpi = 300)
ggsave (plot = FigS7, filename = "FigS7.pdf", units = "cm", device = "pdf",width = 15, height = 10, dpi = 300)


################FigureS8################
#Fitness_cost
################Get Dataset with growth rate################
fdata <- read.csv("Fig_3.csv")

growth_rates <- fdata %>% group_by(strain, NQO, rep, day, gen) %>%
  mutate(note_fit = SummarizeGrowth(hour, od)$vals$note,
         rval = SummarizeGrowth(hour, od)$vals$r,
         kval = SummarizeGrowth(hour, od)$vals$k,
         tmid = SummarizeGrowth(hour, od)$vals$t_mid,
         aucexp = SummarizeGrowth(hour, od)$vals$auc_l,
         max_size = max(od)) %>%
  mutate(note_fit = ifelse(note_fit=="", "OK", "ERROR")) %>%
  slice(1) %>%
  ungroup() 

gr <- filter(growth_rates, note_fit=="OK")


################Need to figure S8################
#Filter gen 200 evolved in NQO
fig200 <- filter(gr, gen == "200")
fig0 <- filter(fig200, NQO == "Yes")
fig0 <- arrange(fig0, rep, strain,day, NQO, gen)

#Filter ancestors
fig0a <- filter(gr, gen == "0")
fig0a <- arrange(fig0a, rep, strain,day)

#Remove two replicas were not evolved
fig0a <- filter(fig0a, rep!=31)
fig0a <- filter(fig0a, rep!=32)
fig0a <- arrange(fig0a, rep, strain,day, NQO, gen)

#Change
fig0$subs = ((fig0$rval) - (fig0a$rval)) 


################Figure S8################
#cost: growth in ypd at t0 versus growth in ypd at t final) 
#Adaptation: 
#subtracting the growth in UV mimetic conditions at T0 from the growth in UV mimetic at T200
#T200 evolved in 4-NQO (in 4-NQO conditions) - T0 ancestors (in 4-NQO conditions)
#Cost values: 
#subtracting the growth in control of strains evolved in UV mimetic (T200) from the growth in control at T0
#T200 evolved in 4-NQO (in control conditions) - T0 ancestors (in control conditions)

#Reload the dataset from EXPEVOL (Figure 2)
fdata <- read.csv("Fig_2.csv")

fdata <- na.omit(fdata)
fdata <- dplyr::filter(fdata, day!=1)
fdata <- dplyr::filter(fdata, day!=2)
fdata <- dplyr::filter(fdata, day!=11)
fdata <- dplyr::filter(fdata, day!=22)
fdata <- dplyr::filter(fdata, day!=23)

fdata$hour = (fdata$time)/3600

growth_rates_adap <- fdata %>% group_by(strain, NQO, rep, day) %>%
  mutate(note_fit = SummarizeGrowth(hour, od)$vals$note,
         rval = SummarizeGrowth(hour, od)$vals$r,
         kval = SummarizeGrowth(hour, od)$vals$k,
         tmid = SummarizeGrowth(hour, od)$vals$t_mid,
         aucexp = SummarizeGrowth(hour, od)$vals$auc_l,
         max_size = max(od)) %>%
  mutate(note_fit = ifelse(note_fit=="", "OK", "ERROR")) %>%
  slice(1) %>%
  ungroup() 

gradap <- dplyr::filter(growth_rates_adap, note_fit=="OK")
gradap <- dplyr::filter(gradap, strain!="White")

fdataAREA <- gradap

#Calculated adap coefficients: 
#Adaptation: 
#subtracting the growth in UV mimetic conditions at T0 from the growth in UV mimetic at T200
#T200 evolved in 4-NQO (in 4-NQO conditions) - T0 ancestors (in 4-NQO conditions)

#%adap= rvalT200 - rvaT0
NQO <- filter(gradap, NQO=="Yes") 
NQO$rep <- NQO$rep
NQO$rep <- as.numeric(NQO$rep)
NQO_3 <- filter(NQO, day=="3")
NQO_21 <- filter(NQO, day=="21")
NQO_3 <- NQO_3 %>% arrange(rep)
NQO_21 <- NQO_21 %>% arrange(rep)

NQO_21$adap = (NQO_21$rval - NQO_3$rval) / NQO_21$rval 

#######(cost: growth in ypd at t0 versus growth in ypd at t final) 
#Cost values: 
#subtracting the growth in control of strains evolved in UV mimetic (T200) from the growth in control at T0
#T200 evolved in 4-NQO (in control conditions) - T0 ancestors (in control conditions)
gen0 <- filter(gr, gen=="0")
gen0 <- filter(gen0, rep!="31")
gen0 <- filter(gen0, rep!="32")

#growth in ypd after 200 gen in NQO
gen200 <- filter(gr, gen=="200")
gen200 <- filter(gen200, NQO=="Yes") ###!= negar algo

gen0$rep <- as.numeric(gen0$rep)
gen0 <- gen0 %>% arrange(rep)
gen200$rep <- as.numeric(gen200$rep)
gen200 <- gen200 %>% arrange(rep)

gen200$cost = (gen0$rval - gen200$rval) / gen0$rval 

gen200  <- gen200 %>% 
  rename(repi = rep)
NQO_21  <- NQO_21 %>% 
  rename(repi = rep)

ALL <- inner_join(x=gen200, y=NQO_21, by=c("strain"="strain","NQO"="NQO","repi"="repi"))

toexpr<-function(x) {
  getfun <- function(x) {
    ifelse(x=="Hybrid", "plain", "italic")
  }
  as.expression(unname(Map(function(f,v) substitute(f(v), list(f=as.name(f), v=as.character(v))), getfun(x), x)))
}


#Fig S8 in cm
FigS8 <- ggscatter(ALL, x = "adap", y = "cost", group = "strain",color ="strain",
                   add = "reg.line",  # Add regressin line,
                   add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                   conf.int = FALSE) + stat_cor(method = "spearman", size=4, label.y = 0.34) + xlab("Cost (%)") + ylab("Adaptation (%)") +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Adaptation rates (OD/hour)") + ylab("Trade-off (OD/hour)") + scale_color_manual(values=c("green4", "dodgerblue1", "#FF9999"), 
                                                                                  labels = toexpr(c("S. cerevisiae", "S. paradoxus", "Hybrid"))) +
  theme(legend.text.align = 0) + guides(color = guide_legend(title = "Genotype")) +
  theme(axis.title.x = element_text(size=14), axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size=14), axis.text.y = element_text(size = 12)) +
  theme(legend.text=element_text(size=9), legend.title=element_text(size=10))+
  theme(legend.spacing.y = unit(-0.2, "cm"),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(0,0,-4,-4))+
  
  scale_x_continuous(breaks = scales::pretty_breaks(n = 5))


################Save Figure S8################
ggsave (plot = FigS8, filename = "FigS8.jpg", units = "cm", device = "jpg",width = 20, height = 15, dpi = 300)
ggsave (plot = FigS8, filename = "FigS8.pdf", units = "cm", device = "pdf",width = 20, height = 15, dpi = 300)



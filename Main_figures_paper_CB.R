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

#Put the path to the folder where all the .csv are
setwd("")


################Figure 1################
################Get Dataset################
fdata <- read.csv("Fig1.csv")
################Figure 1B###############
fdata%<>% mutate(day=as.numeric(day))
fdata%<>% mutate(hour=as.numeric(hour))
fdata$NQO_c <- as.factor(fdata$NQO_c)
fdata0 <- filter(fdata, NQO_c==0)
fdata4 <- filter(fdata, NQO_c==4)
fdata <- rbind(fdata0, fdata4)
fdata <- filter(fdata, strain!="White")
fdata <- filter(fdata, hour < 30)
fdata <- fdata %>% arrange(hour)
fdata <- fdata %>% filter(hour < 25)

toexpr<-function(x) {
  getfun <- function(x) {
    ifelse(x=="Hybrid", "plain", "italic")
  }
  as.expression(unname(Map(function(f,v) substitute(f(v), list(f=as.name(f), v=as.character(v))), getfun(x), x)))
}

f0y4 <- ggplot(fdata, aes(x=hour, y=od, colour=strain, group=interaction(rep,strain,NQO_c))) +
  #stat_summary(position="identity" ,geom="line", fun.y = mean, size=0.1) +
  #stat_summary(fun.data = mean_sdl, geom = "errorbar") +
  geom_line(aes(linetype = NQO_c)) +
  xlab("Time (h)") + ylab("OD (595 nm)") +theme_minimal()+
  theme(axis.text.x = element_text(hjust = 1)) + xlim(0,24) +
  facet_grid(. ~ strain) + 
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),  legend.position = "top",
                     axis.line = element_line(colour = "black"), 
                     strip.background = element_rect(fill = "white", color="white"), 
                     strip.text.x = element_blank()) + scale_color_manual(values=c("green4", "dodgerblue1", "#FF9999"), 
                                                                          labels = toexpr(c("S. cerevisiae", "S. paradoxus", "Hybrid"))) +
  theme(axis.title.x = element_text(size=14), axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size=14), axis.text.y = element_text(size = 12))+
  theme(legend.text.align = 0) + guides(color = guide_legend(title = "Genotype"),
                                        linetype = guide_legend(title = "UV mimetic (μM)")) +
  theme(legend.text=element_text(size=7), legend.title=element_text(size=8))

################Arrange figures################
f0y4 <- plot_grid(f0y4, labels = c("B"), nrow=1 )

ggsave (plot = f0y4, filename = "figure1B.jpg", units = "cm", device = "jpg",width = 15, height = 10, dpi = 300)
ggsave (plot = f0y4, filename = "figure1B.pdf", units = "cm", device = "pdf",width = 15, height = 10, dpi = 300)












################Figure2##############
fdata <- read.csv("Fig_2.csv")
################Get growth rates################
#remove days we do not have all time-points of the curves
fdata <- fdata %>% 
  dplyr::filter(strain!="White", day>2)  %>%
  group_by(day, strain, NQO, rep, plate, tecan, time) %>% 
  dplyr::summarise(od = max(od),max_temp = max(time)) 

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
gr <- dplyr::filter(gr, strain!="White")

fdataAREA <- gr

################Statistical models################
mod2 = lm(rval ~ NQO*day*strain,
          data = fdataAREA)
summary (mod2)


#Remove the outlier to get a god fit in the model
fdataAREA <- filter(fdataAREA, rval < 20)

mod2 = lm(rval ~ NQO*day*strain,
          data = fdataAREA)
summary (mod2)

################Figure 2B################
pd <- position_dodge(0.1) # move them .05 to the left and right

#Confidence intervals
new.fdataAREA <- data.frame(fdataAREA)
pred.pred <- predict(mod2, newdata = new.fdataAREA)
pred.conf <- predict(mod2, newdata = new.fdataAREA, interval = "confidence", level = 0.05)
predict(mod2, newdata = new.fdataAREA, interval = "prediction")
pred.int <- predict(mod2, interval = "prediction")
mydata <- cbind(new.fdataAREA, pred.conf)

#Function to have italic font in legend
toexpr<-function(x) {
  getfun <- function(x) {
    ifelse(x=="Hybrid", "plain", "italic")
  }
  as.expression(unname(Map(function(f,v) substitute(f(v), list(f=as.name(f), v=as.character(v))), getfun(x), x)))
}

#Change column time by generations
mydata <- mydata %>% mutate(gen= day*10)

#Fig 2B in cm
Fig2B <- ggplot(mydata, aes(x=gen, y=rval, colour=strain, shape=NQO, linetype=NQO)) +
  geom_point(position=pd, size=1) +
  ylim(0, 0.8) +
  scale_shape_manual(values=c(16, 1)) + theme_bw() +
  stat_smooth(method = lm, se=TRUE, size = 0.5, level = 0.95) +
  theme(legend.position = "top", axis.title.x = element_text(size=12), axis.text.x = element_text(size = 10),
        axis.title.y = element_text(size=12), axis.text.y = element_text(size = 10)) +  
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black")) +
  xlab("Time (generations)") + ylab("Growth rate (OD/hour)") +
  scale_color_manual(values=c("green4", "dodgerblue1", "#FF9999"),
                     labels = toexpr(c("S. cerevisiae", "S. paradoxus", "Hybrid"))) +
  theme(legend.text.align = 0) + 
  scale_linetype_discrete() +
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

################Get slopes from Figure 2B################
NQO <- filter(mydata, NQO=="Yes")
coef(lmList(rval~gen|strain , data = NQO))

NQO <- filter(mydata, NQO=="No")
coef(lmList(rval~gen|strain , data = NQO))

################Fig 2A################
#We get initial and final growth in UV mimetic conditions
fig0 <- filter(fdataAREA, day == "3")
fig0 <- filter(fig0, NQO == "Yes")
fig0a <- filter(fdataAREA, day == "21")
fig0a <- filter(fig0a, NQO == "Yes")

fig0 <- rbind(fig0, fig0a)

my_comparisons2 <- list( c("1Scer.3", "1Scer.21"),c("2Spar.3", "2Spar.21"),c("3Hybrid.3", "3Hybrid.21"))

fig0$strain_day <- paste(fig0$strain,".", fig0$day)

#Fig2A in cm
Fig2A <- fig0 %>% mutate(strain_day = as.character(interaction (fig0$strain, fig0$day))) %>%
  mutate(strain_day=factor(strain_day, levels= c("1Scer.3", "1Scer.21", "2Spar.3", "2Spar.21", "3Hybrid.3", "3Hybrid.21"))) %>%
  ggplot(., aes(x=as.factor(strain_day), y=rval, col=strain)) + 
  geom_boxplot() +
  theme_minimal()+
  xlab("Generation") + ylab("Growth rate (OD/hour)")+
  ylim(0,0.8) +
  scale_x_discrete(breaks=c("1Scer.3","1Scer.21","2Spar.3","2Spar.21","3Hybrid.3","3Hybrid.21"), labels=c("T0","T200","T0","T200","T0","T200"))+
  theme(legend.position = "top",axis.title.x = element_text(size=12), axis.text.x = element_text(size = 10),
        axis.title.y = element_text(size=12), axis.text.y = element_text(size = 10)) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.ticks = element_line(colour="black"))+
  stat_compare_means(aes(group = strain_day), comparisons= my_comparisons2, method = "t.test", label.y = 0.75, label="p.format", size=2.5) +
  scale_color_manual(values=c("green4", "dodgerblue1", "#FF9999"), 
                     labels = toexpr(c("S. cerevisiae", "S. paradoxus", "Hybrid"))) +
  theme(legend.text.align = 0) + guides(color = guide_legend(title = "Genotype"))+
  geom_jitter(position=position_jitter(0.2)) +
  theme(legend.text=element_text(size=7), legend.title=element_text(size=8))

################Fig 2C################
#graph comparing differences among strains 
fig0 <- filter(fdataAREA, day == "3")
fig0 <- filter(fig0, NQO == "Yes")
fig0 <- arrange(fig0, rep, strain,day)
fig0a <- filter(fdataAREA, day == "21")
fig0a <- filter(fig0a, NQO == "Yes")
fig0a <- arrange(fig0a, rep, strain,day)

fig0$subs = ((fig0a$rval) - (fig0$rval))

my_comparisons2 <- list( c("1Scer.3", "2Spar.3"),c("1Scer.3", "3Hybrid.3"),c("2Spar.3", "3Hybrid.3"))

#in cm
Fig2C <- fig0 %>% mutate(strain_day = as.character(interaction (fig0$strain, fig0$day))) %>%
  mutate(strain_day=factor(strain_day, levels= c("1Scer.3", "1Scer.21", "2Spar.3", "2Spar.21", "3Hybrid.3", "3Hybrid.21"))) %>%
  ggplot(., aes(x=as.factor(strain_day), y=subs, col=strain)) + 
  geom_boxplot() +
  theme_minimal()+
  xlab("Generation") + ylab("Fitness gain (OD/hour)")+
  #scale_x_discrete(breaks=c("1Scer.3","1Scer.21","2Spar.3","2Spar.21","3Hybrid.3","3Hybrid.21"), labels=c("T0","T200","T0","T200","T0","T200"))+
  theme(axis.title.x = element_text(size=12), axis.text.x = element_text(size = 10),
        axis.title.y = element_text(size=12), axis.text.y = element_text(size = 10)) +
  theme(legend.position="none",panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.ticks = element_line(colour="black"))+
  stat_compare_means(method = "anova", label.y = 0.75, aes(label = ..p.format..), size=2.5) +
  stat_compare_means(comparisons= my_comparisons2, size=2.5, label = "p.format")+ 
  scale_color_manual(values=c("green4", "dodgerblue1", "#FF9999"), 
                     labels = toexpr(c("S. cerevisiae", "S. paradoxus", "Hybrid"))) +
  theme(legend.text.align = 0) + guides(color = guide_legend(title = "Genotype"))+
  geom_jitter(position=position_jitter(0.2)) +
  theme(legend.text=element_text(size=12)) +
  scale_x_discrete("Genotype", labels=expression(italic(S.cerevisiae), italic(S.paradoxus), Hybrid)) +
  theme(legend.text=element_text(size=7), legend.title=element_text(size=8))

Fig2C <- Fig2C + guides(
  color = guide_legend(
    title = "Genotype",
    override.aes = aes(label = "")
  )
)

################Fig 2D################
#GRAPH CORRELATION COEFFICIENTS THROUGH TIME by replicate
correlation_graph <- read.csv("Fig_2D.csv")
correlation_graph <- correlation_graph %>% mutate(gen= day*10)

toexpr<-function(x) {
  getfun <- function(x) {
    ifelse(x=="Hybrid", "plain", "italic")
  }
  as.expression(unname(Map(function(f,v) substitute(f(v), list(f=as.name(f), v=as.character(v))), getfun(x), x)))
}

#Fig 2D in cm
Fig2D <- correlation_graph %>%
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


Fig2D <- Fig2D + geom_hline(yintercept=0, linetype="dashed", color = "black")

Fig2D

################Arrange figures################
FIG1 <- plot_grid(Fig2A,Fig2B, labels = c("A","B"), nrow = 1)
FIGS <- plot_grid(Fig2C,Fig2D, labels = c("C","D"), nrow = 1)

fig2 <- plot_grid(FIG1,FIGS, nrow=2)
fig2

#Save
ggsave (plot = fig2, filename = "figure2.jpg", units = "cm", device = "jpg",width = 20, height = 15, dpi = 300)
ggsave (plot = fig2, filename = "figure2.pdf", units = "cm", device = "pdf",width = 20, height = 15, dpi = 300)



################Figure 3################
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
  ungroup() %>%
  select(strain, NQO, rep, day, time, note_fit,rval, kval,tmid,aucexp, max_size, od, gen)

gr <- filter(growth_rates, note_fit=="OK")

################Figure 3A################
my_comparisons <- list( c("0.No.3Hybrid", "200.No.3Hybrid"), c("0.No.3Hybrid", "200.Yes.3Hybrid"), c("200.No.3Hybrid", "200.Yes.3Hybrid"),
                        c("0.No.1Scer", "200.No.1Scer"), c("0.No.1Scer", "200.Yes.1Scer"), c("200.No.1Scer", "200.Yes.1Scer"),
                        c("0.No.2Spar", "200.No.2Spar"), c("0.No.2Spar", "200.Yes.2Spar"), c("200.No.2Spar", "200.Yes.2Spar"))

toexpr<-function(x) {
  getfun <- function(x) {
    ifelse(x=="Hybrid", "plain", "italic")
  }
  as.expression(unname(Map(function(f,v) substitute(f(v), list(f=as.name(f), v=as.character(v))), getfun(x), x)))
}

gr1 <- filter(gr, rep!=31)
gr1 <- filter(gr1, rep!=32)

#Fig3A in cm
Fig3A <- gr1 %>% ggplot(aes(interaction(gen, NQO, strain),rval, col=strain)) + 
  geom_boxplot(outlier.shape = NA) + theme_classic()+ 
  scale_x_discrete(limit =c("0.No.1Scer","200.No.1Scer", "200.Yes.1Scer",
                            "0.No.2Spar","200.No.2Spar", "200.Yes.2Spar",
                            "0.No.3Hybrid","200.No.3Hybrid", "200.Yes.3Hybrid"),
                   label= c("Ancestors","Evolved in control", "Evolved in UV mimetic",
                            "Ancestors","Evolved in control", "Evolved in UV mimetic",
                            "Ancestors","Evolved in control", "Evolved in UV mimetic"))+
  stat_compare_means(comparisons = my_comparisons, method="t.test", label.y = c(0.96,1,1.04,0.96,1,1.04,0.96,1,1.04), 
                     size=2, paired = TRUE) +
  theme(legend.position = "top", axis.text.x = element_text(angle = 30, hjust = 1), panel.border = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Group") + ylab("Growth rate (OD/hour)") + scale_color_manual(values=c("green4", "dodgerblue1", "#FF9999"), 
                                                                     labels = toexpr(c("S. cerevisiae", "S. paradoxus", "Hybrid"))) +
  theme(legend.text.align = 0) + guides(color = guide_legend(title = "Genotype")) +
  theme(legend.text=element_text(size=7), legend.title=element_text(size=8)) +
  theme(axis.title.x = element_text(size=12), axis.text.x = element_text(size = 6),
        axis.title.y = element_text(size=12), axis.text.y = element_text(size = 10), plot.title = element_text(hjust = 0.5, face="bold", size = 24))+
  ylim (0.35, 1.1) +
  theme(legend.text=element_text(size=5), legend.title=element_text(size=6)) +
  geom_jitter(position=position_jitter(), size=0.4)+
  theme(legend.spacing.y = unit(-0.2, "cm"),legend.margin=margin(0,0,0,0),legend.box.margin=margin(0,0,-4,-4))

################Figure 3B################
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
#fig0$subs = ((fig0a$rval) - (fig0$rval)) / (fig0a$rval) 

my_comparisons2 <- list( c("1Scer.1", "2Spar.1"),c("1Scer.1", "3Hybrid.1"),c("2Spar.1", "3Hybrid.1"))

Fig3B <- fig0 %>% mutate(strain_day = as.character(interaction (fig0$strain, fig0$day))) %>%
  ggplot(., aes(x=as.factor(strain_day), y=subs, col=strain)) + 
  geom_boxplot() +
  theme_minimal()+
  xlab("Generation") + ylab("Change in growth rate (OD/hour)")+
  #scale_x_discrete(breaks=c("1Scer.3","1Scer.21","2Spar.3","2Spar.21","3Hybrid.3","3Hybrid.21"), labels=c("T0","T200","T0","T200","T0","T200"))+
  theme(axis.title.x = element_text(size=12), axis.text.x = element_text(size = 10),
        axis.title.y = element_text(size=10), axis.text.y = element_text(size = 10)) +
  theme(legend.position="none",panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.ticks = element_line(colour="black"))+
  stat_compare_means(method = "anova", label.y = 0.25, aes(label = ..p.format..), size=2.5) +
  stat_compare_means(comparisons= my_comparisons2, size=2.5)+ 
  scale_color_manual(values=c("green4", "dodgerblue1", "#FF9999"), 
                     labels = toexpr(c("S. cerevisiae", "S. paradoxus", "Hybrid"))) +
  theme(legend.text.align = 0) + guides(color = guide_legend(title = "Genotype"))+
  geom_jitter(position=position_jitter(0.2)) +
  theme(legend.text=element_text(size=7), legend.title=element_text(size=8))+
  scale_x_discrete("Genotype", labels=expression(italic(S.cerevisiae), italic(S.paradoxus), Hybrid))

Fig3B <- Fig3B + guides(
  color = guide_legend(
    title = "Genotype",
    override.aes = aes(label = "")
  )
)

################Figure 3C################
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
  ungroup() %>%
  select(strain, NQO, rep, day, time, note_fit,rval, kval,tmid,aucexp, max_size, od, plate)

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

#Fig 3C in cm
Fig3C <- ggscatter(ALL, x = "adap", y = "cost", group = "strain",color ="strain",
                   add = "reg.line",  # Add regressin line,
                   add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                   conf.int = TRUE) + stat_cor(method = "spearman", size=2.5) + xlab("Cost (%)") + ylab("Adaptation (%)") +
  theme(axis.text.x = element_text(hjust = 1), panel.border = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Adaptation (OD/hour)") + ylab("Trade-off (OD/hour)") + scale_color_manual(values=c("green4", "dodgerblue1", "#FF9999"), 
                                                                                  labels = toexpr(c("S. cerevisiae", "S. paradoxus", "Hybrid"))) +
  theme(legend.text.align = 0) + guides(color = guide_legend(title = "Genotype")) +
  theme(axis.title.x = element_text(size=12), axis.text.x = element_text(size = 10),
        axis.title.y = element_text(size=12), axis.text.y = element_text(size = 10)) +
  theme(legend.text=element_text(size=7), legend.title=element_text(size=8))+
  theme(legend.spacing.y = unit(-0.2, "cm"),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(0,0,-4,-4))


################Figure 3D################
#Ratio figure
ALL$ratio = (ALL$cost)/(ALL$adap)
#Remove outlier
ALL1 <- filter(ALL, ratio < 5)

my_comparisons2 <- list( c("1Scer.1", "2Spar.1"),c("1Scer.1", "3Hybrid.1"),c("2Spar.1", "3Hybrid.1"))

#Fig 3D in cm
fig <- ALL1 %>% mutate(strain_day = as.character(interaction (ALL1$strain, ALL1$day.x))) %>%
  ggplot(., aes(x=as.factor(strain_day), y=ratio, col=strain)) + 
  geom_boxplot() +
  theme_minimal()+ 
  xlab("Generation") + ylab("Ratio trade-off/adaptation") +
  #scale_x_discrete(breaks=c("1Scer.3","1Scer.21","2Spar.3","2Spar.21","3Hybrid.3","3Hybrid.21"), labels=c("T0","T200","T0","T200","T0","T200"))+
  theme(axis.title.x = element_text(size=12), axis.text.x = element_text(size = 10),
        axis.title.y = element_text(size=10), axis.text.y = element_text(size = 10)) +
  theme(legend.position="none",panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.ticks = element_line(colour="black"))+
  stat_compare_means(method = "anova", label.y = 3, aes(label = ..p.format..), size=2.5) +
  stat_compare_means(comparisons= my_comparisons2, size=2.5)+ 
  scale_color_manual(values=c("green4", "dodgerblue1", "#FF9999"), 
                     labels = toexpr(c("S. cerevisiae", "S. paradoxus", "Hybrid"))) +
  theme(legend.text.align = 0) + guides(color = guide_legend(title = "Genotype"))+
  geom_jitter(position=position_jitter(0.2)) +
  theme(legend.text=element_text(size=12))+
  scale_x_discrete("Genotype", labels=expression(italic(S.cerevisiae), italic(S.paradoxus), Hybrid))+
  theme(legend.text=element_text(size=7), legend.title=element_text(size=8))


Fig3D <- fig + guides(
  color = guide_legend(
    title = "Genotype",
    override.aes = aes(label = "")
  )
)

################Arrange Figures################
fig3a <- plot_grid(Fig3A,Fig3B, labels = c("A","B"), nrow = 1)
fig3b <- plot_grid(Fig3C,Fig3D, labels = c("C","D"), nrow = 1)
fig3 <- plot_grid(fig3a, fig3b, nrow = 2)

ggsave (plot = fig3, filename = "figure3.jpg", units = "cm", device = "jpg",width = 20, height = 15, dpi = 300)
ggsave (plot = fig3, filename = "figure3.pdf", units = "cm", device = "pdf",width = 20, height = 15, dpi = 300)




################Figure 4################
################Fitness data############
fdata <- read.csv("Fig_2.csv")

#Remove NA data
fdata <- na.omit(fdata)
fdata$day <- as.numeric(fdata$day)

################Get growth rates################
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
gr <- dplyr::filter(gr, strain!="White")

data_max <- gr

fdataAREA <- filter(fdataAREA, rval < 20)

#graph comparing if differences among strains are different
fig0 <- filter(fdataAREA, day == "3")
fig0 <- filter(fig0, NQO == "Yes")
fig0 <- arrange(fig0, rep, strain,day)
fig0a <- filter(fdataAREA, day == "21")
fig0a <- filter(fig0a, NQO == "Yes")
fig0a <- arrange(fig0a, rep, strain,day)

fig0$subs = ((fig0a$rval) - (fig0$rval))

my_comparisons2 <- list( c("1Scer.3", "2Spar.3"),c("1Scer.3", "3Hybrid.3"),c("2Spar.3", "3Hybrid.3"))


#Function to have italic font in legend
toexpr<-function(x) {
  getfun <- function(x) {
    ifelse(x=="Hybrid", "plain", "italic")
  }
  as.expression(unname(Map(function(f,v) substitute(f(v), list(f=as.name(f), v=as.character(v))), getfun(x), x)))
}

################Figure 4 and Ploidy Dataset################
Ploidy_2 <- read.csv("Fig_4.csv")

#Prepare data
Ploidy_2$rep <- as.character(Ploidy_2$rep)
fig0$strain_day <- paste(fig0$strain,".", fig0$day)
fig0$rep <- as.character(fig0$rep) 
Ploidy_2$rep <- as.character(Ploidy_2$rep) 
data_max <- inner_join(Ploidy_2, fig0, by = c("strain", "NQO", "rep"))

my_comparisons2 <- list( c("1Scer.3", "2Spar.3"),c("1Scer.3", "3Hybrid.3"),c("2Spar.3", "3Hybrid.3"))
data_max <- dplyr::filter(data_max, rep!= "31")
data_max<- dplyr::filter(data_max, rep!= "32")

#Figure
my_comparisons <- list( c("2_2", "2_3"), c("2_2","3_3"), c("2_3","3_3"),
                        c("3_2","2_2"), c("3_2","3_3"), c("3_2","2_3"))

data_max <- data_max %>%  mutate(strain = gsub(strain, pattern = "1Scer",replacement = "S.cerevisiae")) 
data_max <- data_max %>%  mutate(strain = gsub(strain, pattern = "2Spar",replacement = "S.paradoxus"))
data_max <- data_max %>%  mutate(strain = gsub (strain, pattern = "3Hybrid",replacement = "Hybrid")) 

my_comparisons <- list( c("2_2", "2_3"), c("2_2","3_2"), c("2_2","3_3"),
                        c("2_3","3_2"),c("2_3","3_3"),c("3_2","3_3"))

p200 <- filter(data_max, gen=="200_UV")
p200 <- filter(p200, strain=="Hybrid")
#p1 <- p200 %>% ggplot(., aes(Change,subs)) + 
geom_boxplot(col="#FF9999") + 
  theme_bw() +
  ylim(-0.1,0.5)+
  #ggtitle("Hybrid") +
  stat_compare_means(method="kruskal.test", label = "p.format", label.y=0.5) +
  xlab("Ploidy change") + ylab("Fitness gain (OD/hour)") +
  theme(axis.title.x = element_text(size=14), axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size=14), axis.text.y = element_text(size = 12)) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
  theme(strip.text.y = element_text(size=12))+
  theme(plot.title=element_text(size=12))+
  #stat_compare_means(comparisons = my_comparisons, label = "p.format", label.y = c(0.41, 0.44, 0.47,0.50,0.53,0.57)) +
  labs(title="Hybrid")+ geom_jitter(position=position_jitter(0.2), col="#FF9999")+
  annotate("text",x=c("2_2","2_3","3_3","3_2"), y=-0.1, label=c("(26)","(3)","(1)","(1)"))

#get p valor 0.64
p1 <- p200 %>% ggplot(., aes(Change,subs)) + 
  geom_boxplot(col="#FF9999") + 
  theme_bw() +
  ylim(-0.1,0.5)+
  #ggtitle("Hybrid") +
  #stat_compare_means(method="kruskal.test", label = "p.format", label.y=0.5) +
  xlab("Ploidy change") + ylab("Fitness gain (OD/hour)") +
  theme(axis.title.x = element_text(size=14), axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size=14), axis.text.y = element_text(size = 12)) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
  theme(strip.text.y = element_text(size=12))+
  theme(plot.title=element_text(size=12))+
  #stat_compare_means(comparisons = my_comparisons, label = "p.format", label.y = c(0.41, 0.44, 0.47,0.50,0.53,0.57)) +
  labs(title="Hybrid")+ geom_jitter(position=position_jitter(0.2), col="#FF9999")+
  annotate("text",x=c("2_2","2_3","3_3","3_2", "2_3"), y=c(-0.1,-0.1,-0.1,-0.1, 0.5), label=c("(26)","(3)","(1)","(1)", "0.64"))


my_comparisons <- list( c("2_2", "2_3"), c("2_2","3_3"), c("2_3","3_3"))

p200 <- filter(data_max, gen=="200_UV")
p200 <- filter(p200, strain=="S.cerevisiae")

#witohout 3_2
#p2 <- p200 %>% ggplot(., aes(Change,subs)) + 
geom_boxplot(col="green4") + theme_bw() +
  ylim(0,0.6)+
  #ggtitle("Hybrid") +
  stat_compare_means(method="kruskal.test", label = "p.format", label.y = 0.6) +
  #stat_compare_means(comparisons = my_comparisons, label = "p.format", label.y = c(0.45, 0.48, 0.51)) +
  xlab("Ploidy change") + ylab("Fitness gain (OD/hour)") +
  theme(axis.title.x = element_text(size=14), axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size=14), axis.text.y = element_text(size = 12)) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
  theme(strip.text.y = element_text(face = "italic", size=15)) +
  labs(title="Saccharomyces cerevisiae")+
  scale_x_discrete(breaks = c("2_2","2_3","3_3"),
                   labels=as.character(c("2_2","2_3","3_3"))) +
  theme(plot.title=element_text(face = "italic")) +
  annotate("text",x=c("2_2","2_3","3_3"), 
           y=0.42, label=c("(26)","(2)","(2)")) + geom_jitter(position=position_jitter(0.2), col="green4")

#We get the pvalue = 0.62

p2 <- p200 %>% ggplot(., aes(Change,subs)) + 
  geom_boxplot(col="green4") + theme_bw() +
  ylim(-0.1,0.5)+
  #ggtitle("Hybrid") +
  stat_compare_means(method="kruskal.test", label = "p.format", label.y = 0.6) +
  #stat_compare_means(comparisons = my_comparisons, label = "p.format", label.y = c(0.45, 0.48, 0.51)) +
  xlab("Ploidy change") + ylab("Fitness gain (OD/hour)") +
  theme(axis.title.x = element_text(size=14), axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size=14), axis.text.y = element_text(size = 12)) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
  theme(strip.text.y = element_text(face = "italic", size=12)) +
  labs(title="S. cerevisiae")+
  scale_x_discrete(breaks = c("2_2","2_3","3_3","3_2"),
                   labels=as.character(c("2_2","2_3","3_3","3_2"))) +
  theme(plot.title=element_text(face = "italic", size=12)) +
  annotate("text",x=c("2_2","2_3","3_3","3_2","2_3"), 
           y=c(-0.1,-0.1,-0.1,-0.1, 0.5), label=c("(26)","(2)","(2)","(0)","0.62")) +
  geom_jitter(position=position_jitter(0.2), col="green4")

#Spar
my_comparisons <- list(c("2_2","3_3"))

p200 <- filter(data_max, gen=="200_UV")
p200 <- filter(p200, strain=="S.paradoxus")
#p3 <- p200 %>% ggplot(., aes(Change,subs)) +
geom_boxplot(col="dodgerblue1") +
  theme_bw() +
  ylim(0,0.6)+
  #ggtitle("Hybrid") +
  stat_compare_means(method="wilcox.test", label = "p.format", label.y = 0.6) +
  xlab("Ploidy change") + ylab("Fitness gain (OD/hour)")+
  theme(axis.title.x = element_text(size=14), axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size=14), axis.text.y = element_text(size = 12)) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
  theme(strip.text.y = element_text(face = "italic", size=15)) +
  scale_x_discrete(breaks = c("2_2","3_2"),labels=as.character(c("2_2","3_2"))) +
  labs(title="Saccharomyces paradoxus")+
  theme(plot.title=element_text(face = "italic"),size=12) +
  geom_jitter(position=position_jitter(0.2), col="dodgerblue1") +
  annotate("text",x=c("2_2","3_3"), 
           y=0.5, label=c("(28)","(0)")) 

#Get p-value p = 0.23
p3 <- p200 %>% ggplot(., aes(Change,subs)) +
  geom_boxplot(col="dodgerblue1") +
  theme_bw() +
  ylim(-0.1,0.5)+
  #ggtitle("Hybrid") +
  stat_compare_means(method="wilcox.test", label = "p.format", label.y = 0.6) +
  xlab("Ploidy change") + ylab("Fitness gain (OD/hour)")+
  theme(axis.title.x = element_text(size=14), axis.text.x = element_text(size = 12),
        axis.title.y = element_text(size=14), axis.text.y = element_text(size = 12)) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) + 
  theme(strip.text.y = element_text(face = "italic", size=12)) +
  scale_x_discrete(breaks = c("2_2","2_3","3_3","3_2"),labels=as.character(c("2_2","2_3","3_3","3_2"))) +
  labs(title="S. paradoxus")+
  theme(plot.title=element_text(face = "italic",size=12)) +
  geom_jitter(position=position_jitter(0.2), col="dodgerblue1") +
  annotate("text",x=c("2_2","2_3","3_3","3_2", "2_3"), 
           y=c(-0.1, -0.1, -0.1, -0.1, 0.5), label=c("(28)","(0)","(2)","(0)", "0.23")) 
#Figure
fig4 <- plot_grid(p2,p3,p1, labels = c("","",""), nrow = 1)

################Save Figure 4################
ggsave (plot = fig4, filename = "figure4.jpg", units = "cm", device = "jpg",width = 20, height = 10, dpi = 300)
ggsave (plot = fig4, filename = "figure4.pdf", units = "cm", device = "pdf",width = 20, height = 10, dpi = 300)


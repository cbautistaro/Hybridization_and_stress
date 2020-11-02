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
library(rstatix)
library(readxl)
library(gridExtra)
library(grid)
library(ggpubr)
library(growthcurver)
library(Cairo)
library(Hmisc)
library(nlme)
library(Matrix)
library(lme4)
library(drc)
library(HDPenReg)
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

#Put the path to the folder where all the .csv are
setwd("")



################Figure 1################
################Get Dataset################
fdata <- read.csv("Fig_1.csv")
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


growth_rates1 <- fdata %>% group_by(strain, rep, NQO_c) %>%
  dplyr::mutate(note_fit = SummarizeGrowth(hour, od)$vals$note,
                rval = SummarizeGrowth(hour, od)$vals$r,
                kval = SummarizeGrowth(hour, od)$vals$k,
                tmid = SummarizeGrowth(hour, od)$vals$t_mid,
                aucexp = SummarizeGrowth(hour, od)$vals$auc_l,
                max_size = max(od)) %>%
  mutate(note_fit = ifelse(note_fit=="", "OK", "ERROR")) %>%
  slice(1) %>%
  ungroup() 

gr <- dplyr::filter(growth_rates1, note_fit=="OK")
gr <- dplyr::filter(gr, strain!="White")

fdataAREA <- gr
fdataAREA <- filter(fdataAREA, rval<20)
fdataAREA <- filter(fdataAREA, NQO_c=="4")
my_comparisons2 <- list( c("1Scer", "2Spar"),c("1Scer", "3Hybrid"),c("2Spar", "3Hybrid"))

################Figure 1B################
#in cm
toexpr<-function(x) {
  getfun <- function(x) {
    ifelse(x=="Hybrid", "plain", "italic")
  }
  as.expression(unname(Map(function(f,v) substitute(f(v), list(f=as.name(f), v=as.character(v))), getfun(x), x)))
}

f0y4 <- ggplot(fdata, aes(x=hour, y=od, colour=strain, group=interaction(rep,strain,NQO_c))) +
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
  theme(axis.title.x = element_text(size=14), axis.text.x = element_text(size = 10),
        axis.title.y = element_text(size=14), axis.text.y = element_text(size = 10))+
  theme(legend.text.align = 0) + guides(color = guide_legend(title = "Genotype"),
                                        linetype = guide_legend(title = "UV mimetic (μM)")) +
  theme(legend.text=element_text(size=7), legend.title=element_text(size=8)) +
  theme(legend.text=element_text(size=7), legend.title=element_text(size=8), 
        legend.direction = "horizontal", legend.box = "vertical") +
  theme(legend.spacing.y = unit(-0.2, "cm"),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(0,0,-4,-4))
################Figure 1C################

# Data preparation
df <- fdataAREA
df$strain <- as.factor(df$strain)

# Tukey HSD
stat.test <- aov(rval ~ strain, data = df) %>%
  tukey_hsd()

Fig1C <- fdataAREA %>% 
  ggplot(., aes(x=as.factor(strain), y=rval, col=strain)) + 
  geom_boxplot(outlier.shape = NA) +
  theme_minimal()+
  xlab("Generation") + ylab("Growth rate (OD/hour)")+
  theme(axis.title.x = element_text(size=12), axis.text.x = element_text(size = 10),
        axis.title.y = element_text(size=12), axis.text.y = element_text(size = 10)) +
  theme(legend.position="none",panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.ticks = element_line(colour="black"))+
  stat_compare_means(method = "anova", label.y = 0.75, aes(label = ..p.format..), size=2.5) +
  stat_pvalue_manual(
    stat.test, label = "p.adj", 
    y.position = c(0.7, 0.72, 0.74),
    size=2.5
  )+
  scale_color_manual(values=c("green4", "dodgerblue1", "#FF9999"), 
                     labels = toexpr(c("S. cerevisiae", "S. paradoxus", "Hybrid"))) +
  theme(legend.text.align = 0) + guides(color = guide_legend(title = "Genotype"))+
  geom_jitter(position=position_jitter(0.2), size=0.6) +
  theme(legend.text=element_text(size=12)) +
  scale_x_discrete("Genotype", labels=expression(italic(S.cerevisiae), italic(S.paradoxus), Hybrid)) +
  theme(legend.text=element_text(size=7), legend.title=element_text(size=8))

################Arrange figures################
Fig1 <- plot_grid(f0y4, Fig1C, labels = c("B", "C"), nrow=1 )

ggsave (plot = Fig1, filename = "Fig1.jpg", units = "cm", device = "jpg",width = 20, height = 10, dpi = 300)
ggsave (plot = Fig1, filename = "Fig1.pdf", units = "cm", device = "pdf",width = 20, height = 10, dpi = 300)
ggsave (plot = Fig1, filename = "Fig1.tiff", units = "cm", device = "tiff",width = 20, height = 10, dpi = 300)


################Figure2 for the model##############
#Only for the model
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
  ungroup() 

gr <- dplyr::filter(growth_rates, note_fit=="OK")
gr <- dplyr::filter(gr, strain!="White")

fdataAREA <- gr

################Statistical models################
fdataAREA <- fdataAREA %<>% mutate(strain=ifelse(strain=="1Scer", "Scer",
                                    ifelse(strain=="2Spar", "Spar",
                                           ifelse(strain=="3Hybrid", "Hybrid", NA)))) 
                                  
lm_NQOdaystrain = lm(rval ~ NQO*day*strain,
          data = fdataAREA)
summary (lm_NQOdaystrain)


#Remove the outlier to get a better fit in the model
fdataAREA <- filter(fdataAREA, rval < 20)

lm_NQOdaystrain = lm(rval ~ NQO*day*strain,
          data = fdataAREA)
summary (lm_NQOdaystrain)

################Statistical models without NQO effect################
#Linnear
NQOYES <- filter(fdataAREA, NQO=="Yes")

lm_daystrain = lm(rval ~ day*strain,
          data = NQOYES)
summary (lm_daystrain)


################Data for figure 2################
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
  ungroup() 

gr <- dplyr::filter(growth_rates, note_fit=="OK")
gr <- dplyr::filter(gr, strain!="White")

fdataAREA <- gr

################Figure 2B################

#Function to have italic font in legend
toexpr<-function(x) {
  getfun <- function(x) {
    ifelse(x=="Hybrid", "plain", "italic")
  }
  as.expression(unname(Map(function(f,v) substitute(f(v), list(f=as.name(f), v=as.character(v))), getfun(x), x)))
}

#Change column time by generations
fdataAREA <- fdataAREA %>% mutate(gen= day*5)

pd <- position_dodge(0.1)

#Fig 2B in cm
Fig2B <- ggplot(fdataAREA, aes(x=gen, y=rval, colour=strain, shape=NQO, linetype=NQO)) +
  geom_point(position=pd, size=0.5) +
  geom_line(alpha=0.5, aes(group = interaction(strain,rep,NQO))) +
  ylim(0, 0.8) +
  geom_jitter(size=0.8)+
  scale_shape_manual(values=c(16, 1)) + theme_bw() +
  theme(legend.position = "top", axis.title.x = element_text(size=12), axis.text.x = element_text(size = 10),
        axis.title.y = element_text(size=12), axis.text.y = element_text(size = 10)) +  
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black")) +
  xlab("Time (generations)") + ylab("Growth rate (OD/hour)") +
  scale_color_manual(values=c("green4", "dodgerblue1", "#FF9999"),
                     labels = toexpr(c("S. cerevisiae", "S. paradoxus", "Hybrid"))) +
  theme(legend.text.align = 0) + 
  facet_grid(. ~strain)+
  scale_linetype_discrete() +
  guides(color = guide_legend(title = "Genotype",
                              override.aes=list(fill=NA)), 
         linetype = guide_legend(title = "UV mimetic", 
                                 override.aes=list(fill=NA)),
         shape = guide_legend(title = "UV mimetic", override.aes=list(colour="black"))) +
  theme(legend.text=element_text(size=6), legend.title=element_text(size=7), 
        legend.direction = "horizontal", legend.box = "vertical") +
  theme(legend.spacing.y = unit(-0.2, "cm"),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(0,0,-4,-4)) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),  legend.position = "top",
          axis.line = element_line(colour = "black"), 
          strip.background = element_rect(fill = "white", color="white"), 
          strip.text.x = element_blank()) +
  scale_x_continuous(breaks = c(25,50,75,100))

################Figure 2C1################
#To extract legend
Fig2C1 <- fdataAREA %>% filter(NQO =="Yes") %>% 
  ggplot(., aes(x=gen, y=rval, colour=strain, shape=NQO, linetype=NQO)) +
  ylim(0, 0.8) +
  scale_shape_manual(values=16) + theme_bw() +
  stat_smooth(method = lm, se=TRUE, size = 0.5, level = 0.95) +
  theme(legend.position = "top", axis.title.x = element_text(size=12), axis.text.x = element_text(size = 9),
        axis.title.y = element_text(size=12), axis.text.y = element_text(size = 10)) +  
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black")) +
  xlab("Time (generations)") + ylab("Growth rate (OD/hour)") +
  scale_color_manual(values=c("green4", "dodgerblue1", "#FF9999"),
                     labels = toexpr(c("S. cerevisiae", "S. paradoxus", "Hybrid"))) +
  theme(legend.text.align = 0) + 
  scale_linetype_manual(values=c("solid","solid")) +
  guides(color = guide_legend(title = "Genotype",
                              override.aes=list(fill=NA)), 
         linetype = guide_legend(title = "UV mimetic", 
                                 override.aes=list(fill=NA)),
         shape = guide_legend(title = "UV mimetic", override.aes=list(colour="black"))) +
  theme(legend.text=element_text(size=6), legend.title=element_text(size=7), 
        legend.direction = "horizontal", legend.box = "vertical") +
  theme(legend.spacing.y = unit(-0.2, "cm"),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(0,0,-4,-4)) +
  ggtitle("Linear Model") +
  theme(plot.title = element_text(hjust = 0.5))

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

mylegend<-g_legend(Fig2C1)

#Figure 2C1
Fig2C1 <- fdataAREA %>% filter(NQO =="Yes") %>% 
  ggplot(., aes(x=gen, y=rval, colour=strain, shape=NQO, linetype=NQO)) +
  ylim(0, 0.8) +
  scale_shape_manual(values=16) + theme_bw() +
  stat_smooth(method = lm, se = FALSE) +
  theme(legend.position = "NONE", axis.title.x = element_text(size=10), axis.text.x = element_text(size = 9),
        axis.title.y = element_text(size=12), axis.text.y = element_text(size = 10)) +  
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black")) +
  scale_linetype_manual(values=c("solid","solid")) +
  xlab("Time (generations)") + ylab("Growth rate (OD/hour)") +
  scale_color_manual(values=c("green4", "dodgerblue1", "#FF9999"),
                     labels = toexpr(c("S. cerevisiae", "S. paradoxus", "Hybrid"))) +
  theme(legend.text.align = 0) + 
  guides(color = guide_legend(title = "Genotype",
                              override.aes=list(fill=NA)), 
         linetype = guide_legend(title = "UV mimetic", 
                                 override.aes=list(fill=NA)),
         shape = guide_legend(title = "UV mimetic", override.aes=list(colour="black"))) +
  theme(legend.text=element_text(size=6), legend.title=element_text(size=7), 
        legend.direction = "horizontal", legend.box = "vertical") +
  theme(legend.spacing.y = unit(-0.2, "cm"),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(0,0,-4,-4)) +
  ggtitle("Linear Model") +
  theme(plot.title = element_text(hjust = 0.5 , size = 9)) +
  scale_x_continuous(breaks = c(25,50,75,100))

################Figure 2C2################
#Second panel
fdata <- read.csv("Fig_2.csv")
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
  ungroup() 

growth_rates$gen <- growth_rates$day*5

gr <- dplyr::filter(growth_rates, note_fit=="OK")
gr <- dplyr::filter(gr, strain!="White")

fdataAREA <- gr


nitro1 <- filter(gr, NQO=="Yes") 

nitro1 <- filter(nitro1, rval<20)
df <- nitro1

#view the data in normal and log scale for dose
p1 <- df %>% ggplot() + geom_point(aes(gen, rval, color = strain)) + theme_bw()
p2 <- df %>% ggplot() + geom_point(aes(log((gen)), rval, color = strain)) + theme_bw()
ggarrange(p1, p2)



nitro1$gen <- nitro1$day * 5
df <- nitro1

#Asymtotic
# define drm function to use with map
drm.func <- function(x) {
  drm(rval ~ gen, 
      fct = LL2.2(),#names = c("Slope", "Lower Limit", "Upper Limit", "ED50")), 
      data = x)
}

predict.fun <- function(x) {
  add_predictions(data.frame(gen = seq(15,105)), x)
}

coefs.fun <- function(x) {coef(x) %>% tidy}

df2 <- df %>% group_by(strain) %>% nest() %>%
  mutate(drmod = map(data, drm.func), 
         pred = map(drmod, predict.fun),
         coefs = map(drmod, coefs.fun))

# plot raw data, model and ED50 line
toexpr<-function(x) {
  getfun <- function(x) {
    ifelse(x=="Hybrid", "plain", "italic")
  }
  as.expression(unname(Map(function(f,v) substitute(f(v), list(f=as.name(f), v=as.character(v))), getfun(x), x)))
}

Fig2C2<- df2 %>% unnest(data) %>% 
  ggplot() + 
  geom_line(aes(gen, pred, color = strain), data =
              df2 %>% unnest(pred), linetype = "solid", size=0.8) +
  geom_vline(aes(xintercept = log(x), color = strain),
             data = df2 %>% unnest(coefs) %>% filter(names == "ED50:(Intercept)")) +
  theme_bw() +ylim(0,0.8) +xlim(15,130) +
  theme(legend.position = "NONE", axis.title.x = element_text(size=10), axis.text.x = element_text(size = 9),
        axis.title.y = element_text(size=12), axis.text.y = element_text(size = 10)) +  
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black")) +
  xlab("Time (generations)") + ylab("Growth rate (OD/hour)") +
  scale_color_manual(values=c("green4", "dodgerblue1", "#FF9999"),
                     labels = toexpr(c("S. cerevisiae", "S. paradoxus", "Hybrid"))) +
  theme(legend.text.align = 0) + 
  scale_x_continuous(breaks = c(25,50,75,100)) +
  guides(color = guide_legend(title = "Genotype",
                              override.aes=list(fill=NA)), 
         linetype = guide_legend(title = "UV mimetic", 
                                 override.aes=list(fill=NA)),
         shape = guide_legend(title = "UV mimetic", override.aes=list(colour="black"))) +
  theme(legend.text=element_text(size=4), legend.title=element_text(size=5), 
        legend.direction = "horizontal", legend.box = "vertical") +
  theme(legend.spacing.y = unit(-0.2, "cm"),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(0,0,-4,-4)) +
  ggtitle("Asymptotic Model") +
  theme(plot.title = element_text(hjust = 0.5, size = 9))


#Logarithmic
#Fig 2C3
# define drm function to use with map
drm.func <- function(x) {
  drm(rval ~ gen, 
      fct = DRC.logCurve(),#names = c("Slope", "Lower Limit", "Upper Limit", "ED50")), 
      data = x)
}

predict.fun <- function(x) {
  add_predictions(data.frame(gen = seq(15,105)), x)
}

coefs.fun <- function(x) {coef(x) %>% tidy}

df2 <- df %>% group_by(strain) %>% nest() %>%
  mutate(drmod = map(data, drm.func), 
         pred = map(drmod, predict.fun),
         coefs = map(drmod, coefs.fun))

# plot raw data, model and ED50 line
toexpr<-function(x) {
  getfun <- function(x) {
    ifelse(x=="Hybrid", "plain", "italic")
  }
  as.expression(unname(Map(function(f,v) substitute(f(v), list(f=as.name(f), v=as.character(v))), getfun(x), x)))
}

Fig2C3<- df2 %>% unnest(data) %>% 
  ggplot() + 
  geom_line(aes(gen, pred, color = strain), data =
              df2 %>% unnest(pred), linetype = "solid", size=0.8) +
  geom_vline(aes(xintercept = log(x), color = strain),
             data = df2 %>% unnest(coefs) %>% filter(names == "ED50:(Intercept)")) +
  theme_bw() +ylim(0,0.8) +xlim(15,130) +
  theme(legend.position = "NONE", axis.title.x = element_text(size=10), axis.text.x = element_text(size = 9),
        axis.title.y = element_text(size=12), axis.text.y = element_text(size = 10)) +  
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black")) +
  xlab("Time (generations)") + ylab("Growth rate (OD/hour)") +
  scale_color_manual(values=c("green4", "dodgerblue1", "#FF9999"),
                     labels = toexpr(c("S. cerevisiae", "S. paradoxus", "Hybrid"))) +
  theme(legend.text.align = 0) + 
  scale_x_continuous(breaks = c(25,50,75,100)) +
  guides(color = guide_legend(title = "Genotype",
                              override.aes=list(fill=NA)), 
         linetype = guide_legend(title = "UV mimetic", 
                                 override.aes=list(fill=NA)),
         shape = guide_legend(title = "UV mimetic", override.aes=list(colour="black"))) +
  theme(legend.text=element_text(size=4), legend.title=element_text(size=5), 
        legend.direction = "horizontal", legend.box = "vertical") +
  theme(legend.spacing.y = unit(-0.2, "cm"),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(0,0,-4,-4)) +
  ggtitle("Logarithmic Model") +
  theme(plot.title = element_text(hjust = 0.5, size = 9))



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
  geom_boxplot(outlier.shape = NA) +
  theme_minimal()+
  xlab("Generation") + ylab("Growth rate (OD/hour)")+
  ylim(0,0.8) +
  scale_x_discrete(breaks=c("1Scer.3","1Scer.21","2Spar.3","2Spar.21","3Hybrid.3","3Hybrid.21"), labels=c("T15","T100","T15","T100","T15","T100"))+
  theme(legend.position = "top",axis.title.x = element_text(size=12), axis.text.x = element_text(size = 9),
        axis.title.y = element_text(size=12), axis.text.y = element_text(size = 10)) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.ticks = element_line(colour="black"))+
  stat_compare_means(aes(group = strain_day), comparisons= my_comparisons2, method = "t.test", label.y = 0.75, label="p.format", size=2.5) +
  scale_color_manual(values=c("green4", "dodgerblue1", "#FF9999"), 
                     labels = toexpr(c("S. cerevisiae", "S. paradoxus", "Hybrid"))) +
  theme(legend.text.align = 0) + guides(color = guide_legend(title = "Genotype"))+
  geom_jitter(position=position_jitter(0.2), size=0.6) +
  theme(legend.text=element_text(size=6), legend.title=element_text(size=7))

################Fig 2D################
#graph comparing differences among strains 
fig0 <- filter(fdataAREA, day == "3")
fig0 <- filter(fig0, NQO == "Yes")
fig0 <- arrange(fig0, rep, strain,day)
fig0a <- filter(fdataAREA, day == "21")
fig0a <- filter(fig0a, NQO == "Yes")
fig0a <- arrange(fig0a, rep, strain,day)

fig0$subs = ((fig0a$rval) - (fig0$rval))

fig0 <- fig0 %>% mutate(strain_day = as.character(interaction (fig0$strain, fig0$day))) %>%
  mutate(strain_day=factor(strain_day, levels= c("1Scer.3", "1Scer.21", "2Spar.3", "2Spar.21", "3Hybrid.3", "3Hybrid.21"))) 

stat.test <- aov(subs ~ strain_day, data = fig0) %>%
  tukey_hsd()

#in cm
Fig2D <- ggplot(fig0, aes(x=as.factor(strain_day), y=subs, col=strain)) + 
  geom_boxplot(outlier.shape = NA) +
  theme_minimal()+
  xlab("Generation") + ylab("Fitness gain (OD/hour)")+
  theme(axis.title.x = element_text(size=12), axis.text.x = element_text(size = 8),
        axis.title.y = element_text(size=12), axis.text.y = element_text(size = 10)) +
  theme(legend.position="none",panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.ticks = element_line(colour="black"))+
  stat_compare_means(method = "anova", label.y = 0.65, aes(label = ..p.format..), size=2.5) +
  stat_pvalue_manual(
    stat.test, label = "p.adj", 
   y.position = c(0.5, 0.55, 0.6),
   size=2.5
  )+
  scale_color_manual(values=c("green4", "dodgerblue1", "#FF9999"), 
                     labels = toexpr(c("S. cerevisiae", "S. paradoxus", "Hybrid"))) +
  theme(legend.text.align = 0) + guides(color = guide_legend(title = "Genotype"))+
  geom_jitter(position=position_jitter(0.2), size=0.6) +
  theme(legend.text=element_text(size=12)) +
  scale_x_discrete("Genotype", labels=expression(italic(S.cerevisiae), italic(S.paradoxus), Hybrid)) +
  theme(legend.text=element_text(size=7), legend.title=element_text(size=8))

Fig2D <- Fig2D + guides(
  color = guide_legend(
    title = "Genotype",
    override.aes = aes(label = "")
  )
)

################Arrange figures################
Fig2Cnolegend <- plot_grid(Fig2C1,Fig2C2,Fig2C3, nrow = 1)
Fig2C <- plot_grid(mylegend,Fig2Cnolegend, nrow = 2,
                   rel_heights = c(1, 5))


FIG1 <- plot_grid(Fig2A,Fig2D, labels = c("A","D"), nrow = 1)
FigC <- plot_grid(Fig2C, labels = c("C"), nrow = 1)
FIGS <- plot_grid(Fig2B, labels = c("B"), nrow = 1)

fig2rev <- plot_grid(FIG1, FIGS, FigC, nrow=3)
fig2rev

#Save
ggsave (plot = fig2rev, filename = "Fig2.jpg", units = "cm", device = "jpg",width = 20, height = 20, dpi = 300)
ggsave (plot = fig2rev, filename = "Fig2.pdf", units = "cm", device = "pdf",width = 20, height = 20, dpi = 300)
ggsave (plot = fig2rev, filename = "Fig2.tiff", units = "cm", device = "tiff",width = 20, height = 20, dpi = 300)

################Figure 3################
#Fitness_cost
################Get Dataset with growth rate################
fdata <- read.csv("Fig_3.csv")

growth_rates <- fdata %>% group_by(strain, NQO, rep, day, gen) %>%
  dplyr::mutate(note_fit = SummarizeGrowth(hour, od)$vals$note,
         rval = SummarizeGrowth(hour, od)$vals$r,
         kval = SummarizeGrowth(hour, od)$vals$k,
         tmid = SummarizeGrowth(hour, od)$vals$t_mid,
         aucexp = SummarizeGrowth(hour, od)$vals$auc_l,
         max_size = max(od)) %>%
  mutate(note_fit = ifelse(note_fit=="", "OK", "ERROR")) %>%
  slice(1) %>%
  ungroup() 

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

fig0 <- fig0 %>% mutate(strain_day = as.character(interaction (fig0$strain, fig0$day))) 
  
my_comparisons2 <- list( c("1Scer.1", "2Spar.1"),c("1Scer.1", "3Hybrid.1"),c("2Spar.1", "3Hybrid.1"))

stat.test <- aov(subs ~ strain_day, data = fig0) %>%
  tukey_hsd()


Fig3B <- ggplot(fig0, aes(x=as.factor(strain_day), y=subs, col=strain)) + 
  geom_boxplot(outlier.shape = NA) +
  theme_minimal()+
  xlab("Generation") + ylab("Change in growth rate (OD/hour)")+
  theme(axis.title.x = element_text(size=12), axis.text.x = element_text(size = 10),
        axis.title.y = element_text(size=10), axis.text.y = element_text(size = 10)) +
  theme(legend.position="none",panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.ticks = element_line(colour="black"))+
  stat_compare_means(method = "anova", label.y = 0.2, aes(label = ..p.format..), size=2.5) +
  stat_pvalue_manual(
    stat.test, label = "p.adj", 
    y.position = c(0.12, 0.15, 0.18),
    size=2.5
    ) +
  scale_color_manual(values=c("green4", "dodgerblue1", "#FF9999"), 
                     labels = toexpr(c("S. cerevisiae", "S. paradoxus", "Hybrid"))) +
  theme(legend.text.align = 0) + guides(color = guide_legend(title = "Genotype"))+
  geom_jitter(position=position_jitter(0.2), size=0.6) +
  theme(legend.text=element_text(size=7), legend.title=element_text(size=8))+
  scale_x_discrete("Genotype", labels=expression(italic(S.cerevisiae), italic(S.paradoxus), Hybrid))

Fig3B <- Fig3B + guides(
  color = guide_legend(
    title = "Genotype",
    override.aes = aes(label = "")
  )
)

################Figure 3C################
#cost: growth in ypd at t0 versus growth in ypd at t final 
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
  dplyr::mutate(note_fit = SummarizeGrowth(hour, od)$vals$note,
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
fdata <- read.csv("Fig_3.csv")

growth_rates <- fdata %>% group_by(strain, NQO, rep, day, gen) %>%
  dplyr::mutate(note_fit = SummarizeGrowth(hour, od)$vals$note,
         rval = SummarizeGrowth(hour, od)$vals$r,
         kval = SummarizeGrowth(hour, od)$vals$k,
         tmid = SummarizeGrowth(hour, od)$vals$t_mid,
         aucexp = SummarizeGrowth(hour, od)$vals$auc_l,
         max_size = max(od)) %>%
  mutate(note_fit = ifelse(note_fit=="", "OK", "ERROR")) %>%
  slice(1) %>%
  ungroup() #%>%
  #select(strain, NQO, rep, day, time, note_fit,rval, kval,tmid,aucexp, max_size, od, gen)

gr <- filter(growth_rates, note_fit=="OK")

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

pd <- position_dodge(0.1)

#Fig 3C in cm
Fig3C <- ggplot(ALL, aes(x=adap, y=cost, colour=strain)) +
  geom_point(position=pd, size=0.6) +
  theme_bw()+
  stat_smooth(method = lm, se=FALSE, size = 0.5, level = 0.95) +
  theme(axis.text.x = element_text(hjust = 1), panel.border = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  xlab("Adaptation (OD/hour)") + ylab("Trade-off (OD/hour)") + scale_color_manual(values=c("green4", "dodgerblue1", "#FF9999"), 
                                                                                  labels = toexpr(c("S. cerevisiae", "S. paradoxus", "Hybrid"))) +
  theme(legend.position = "top", axis.title.x = element_text(size=12), axis.text.x = element_text(size = 10),
        axis.title.y = element_text(size=12), axis.text.y = element_text(size = 10)) +  
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black")) +
  scale_color_manual(values=c("green4", "dodgerblue1", "#FF9999"),
                     labels = toexpr(c("S. cerevisiae", "S. paradoxus", "Hybrid"))) +
  theme(legend.text.align = 0) + 
  guides(color = guide_legend(title = "Genotype",
                              override.aes=list(fill=NA))) +
  theme(legend.text=element_text(size=5), legend.title=element_text(size=6), 
        legend.direction = "horizontal", legend.box = "vertical") +
  theme(legend.spacing.y = unit(-0.2, "cm"),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(0,0,-4,-4))


#Model lm(cost adap)
mod3 = lm(cost ~ adap*strain,
          data = ALL)
summary (mod3)


################Figure 3D################
#Ratio figure
ALL$ratio = (ALL$cost)/(ALL$adap)
#Remove outlier
ALL1 <- filter(ALL, ratio < 5)

my_comparisons2 <- list( c("1Scer.1", "2Spar.1"),c("1Scer.1", "3Hybrid.1"),c("2Spar.1", "3Hybrid.1"))

ALL1 <- ALL1 %>% mutate(strain_day = as.character(interaction (ALL1$strain, ALL1$day.x)))

stat.test <- aov(ratio ~ strain_day, data = ALL1) %>%
  tukey_hsd()

  
#Fig 3D in cm
fig <- #ALL1 %>% mutate(strain_day = as.character(interaction (ALL1$strain, ALL1$day.x))) %>%
  ggplot(ALL1, aes(x=as.factor(strain_day), y=ratio, col=strain)) + 
  geom_boxplot(outlier.shape = NA) +
  theme_minimal()+ 
  xlab("Generation") + ylab("Ratio trade-off/adaptation") +
  theme(axis.title.x = element_text(size=12), axis.text.x = element_text(size = 10),
        axis.title.y = element_text(size=10), axis.text.y = element_text(size = 10)) +
  theme(legend.position="none",panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.ticks = element_line(colour="black"))+
  stat_compare_means(method = "anova", label.y = 2.8, aes(label = ..p.format..), size=2.5) +
  stat_pvalue_manual(
      stat.test, label = "p.adj", 
     y.position = c(2, 2.3, 2.6),
     size=2.5) +
  scale_color_manual(values=c("green4", "dodgerblue1", "#FF9999"), 
                     labels = toexpr(c("S. cerevisiae", "S. paradoxus", "Hybrid"))) +
  theme(legend.text.align = 0) + guides(color = guide_legend(title = "Genotype"))+
  geom_jitter(position=position_jitter(0.2), size=0.6) +
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

ggsave (plot = fig3, filename = "Fig3.jpg", units = "cm", device = "jpg",width = 20, height = 15, dpi = 300)
ggsave (plot = fig3, filename = "Fig3.pdf", units = "cm", device = "pdf",width = 20, height = 15, dpi = 300)
ggsave (plot = fig3, filename = "Fig3.tiff", units = "cm", device = "tiff",width = 20, height = 15, dpi = 300)





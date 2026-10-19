setwd("~/Activité Professionnelle/LIMNO Konstanz (2019-2023)/Collaborations/Chapter 1 Richard - genotypes Trait Variation/Fitness")

rm(list=ls())

library(car)
library(cowplot)
library(data.table)
library(DescTools)
library(deSolve)
library(dplyr)
library(foreach)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(grid)
library(gridExtra)
library(lattice)
library(lme4)
library(lmtest)
library(magrittr)
library(merTools)
library(nlme)
library(plyr)
library(plotly)
library(propagate)
library(qpcR)
library(reshape2)
library(rstatix)
library(scales)

######################################################
######################################################
##### FITNESS OF ALGAE IN DIFFERENT ENVIRONMENTS #####
######################################################
######################################################

# Import the dataset
Data=read.table("Data_F.txt", h=T, dec=",")
Data$Species=gsub("\\_","\\ ", Data$Species)
summary(Data)
names(Data)

# Specify the variables as numeric or factor
Data[,c(7:12,15)] %<>% mutate_if(is.character,as.numeric)
Data[,c(7:12,15)] %<>% mutate_if(is.integer,as.numeric)
Data[,c(7:12,15)] %<>% mutate_if(is.factor,as.numeric)

# Sort the dataset by genotype, nutrient and predator
Data$Genotype=factor(Data$Genotype, levels=unique(Data$Genotype))
Data$Species=factor(Data$Species, levels=unique(Data$Species))
Data$Nutri=factor(Data$Nutri, levels=c("Low","High"))
Data$Pred=factor(Data$Pred, levels=c("Without","With"))
Data=Data[order(Data$Genotype,Data$Nutri,Data$Pred,Data$Trial),]

# Calculate the densities
Data$DensA=Data$IXM*Data$Volu*Data$Site*Data$Dilu*Data$Cove
Data$DensR=Data$Rotifer/5

# Calculate the mean densities
Data2=setDT(Data)[, .(MeanDensA=round(mean(DensA),0), MeanDensR=round(mean(DensR),0)), by=list(Genotype,Species,Period,Nutri,Pred,Day)]
Data2=as.data.frame(Data2)
Data=as.data.frame(Data)


###############################
### Plot the fitness curves ###
###############################

tiff('Fitness Curves Alga Low.tiff', units="in", width=15, height=15, res=1000)
ggplot(subset(Data2, Nutri=="Low"), aes(Day, log(MeanDensA+1), group=Pred)) +
  geom_smooth(aes(color=Species, linetype=Pred), size=1, se=F) +
  geom_point(data=subset(Data, Nutri=="Low"), aes(Day, log(DensA+1), color=Species), alpha=0.4, size=1.5, pch=16, position=position_jitter(0.2)) +
  ylab(expression('Algal density'~'('*ln~cells~mL^-1*')')) + xlab(expression('Time (days)')) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=16)) +  
  theme(axis.text.x=element_text(face="plain", colour="black", size=16)) +  
  theme(axis.title.y=element_text(face="plain", colour="black", size=16)) +
  theme(axis.title.x=element_text(face="plain", colour="black", size=16)) +
  scale_y_continuous(labels=function(x) sprintf("%.0f", x), breaks=seq(5,14,by=3), limits=c(5,14)) +
  scale_x_continuous(labels=function(x) sprintf("%.0f", x), breaks=seq(0,10,by=2), limits=c(0,10)) +
  scale_color_manual(values=c("Chlamydomonas noctigama"="dodgerblue3", "Chlamydomonas bilatus"="firebrick3"))+
  scale_linetype_manual(values=c("Without"="11", "With"="solid")) +
  theme(axis.line=element_line(colour="black")) + theme(panel.background=element_blank()) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  theme(strip.text.x=element_text(face="bold", color="black", size=16)) +
  theme(strip.background=element_blank()) +
  facet_wrap(~Genotype, scales="free", ncol=7, nrow=7) +
  theme(plot.margin=unit(c(5.5,5.5,5.5,5.5),"pt")) +
  theme(legend.position="none")
dev.off()

tiff('Fitness Curves Alga High.tiff', units="in", width=15, height=15, res=1000)
ggplot(subset(Data2, Nutri=="High"), aes(Day, log(MeanDensA+1), group=Pred)) +
  geom_smooth(aes(color=Species, linetype=Pred), size=1, se=F) +
  geom_point(data=subset(Data, Nutri=="High"), aes(Day, log(DensA+1), color=Species), alpha=0.4, size=1.5, pch=16, position=position_jitter(0.2)) +
  ylab(expression('Algal density'~'('*ln~cells~mL^-1*')')) + xlab(expression('Time (days)')) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=16)) +  
  theme(axis.text.x=element_text(face="plain", colour="black", size=16)) +  
  theme(axis.title.y=element_text(face="plain", colour="black", size=16)) +
  theme(axis.title.x=element_text(face="plain", colour="black", size=16)) +
  scale_y_continuous(labels=function(x) sprintf("%.0f", x), breaks=seq(5,14,by=3), limits=c(5,14)) +
  scale_x_continuous(labels=function(x) sprintf("%.0f", x), breaks=seq(0,10,by=2), limits=c(0,10)) +
  scale_color_manual(values=c("Chlamydomonas noctigama"="dodgerblue3", "Chlamydomonas bilatus"="firebrick3"))+
  scale_linetype_manual(values=c("Without"="11", "With"="solid")) +
  theme(axis.line=element_line(colour="black")) + theme(panel.background=element_blank()) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  theme(strip.text.x=element_text(face="bold", color="black", size=16)) +
  theme(strip.background=element_blank()) +
  facet_wrap(~Genotype, scales="free", ncol=7, nrow=7) +
  theme(plot.margin=unit(c(5.5,5.5,5.5,5.5),"pt")) +
  theme(legend.position="none")
dev.off()

tiff('Fitness Curves Rotifer Low.tiff', units="in", width=15, height=15, res=1000)
ggplot(subset(Data2, Nutri=="Low"), aes(Day, log(MeanDensR+1), group=Pred)) +
  geom_smooth(aes(color=Species, linetype=Pred), size=1, se=F) +
  geom_point(data=subset(Data, Nutri=="Low"), aes(Day, log(DensR+1), color=Species), alpha=0.4, size=1.5, pch=16, position=position_jitter(0.2)) +
  ylab(expression('Rotifer density'~'('*ln~individuals~mL^-1*')')) + xlab(expression('Time (days)')) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=16)) +  
  theme(axis.text.x=element_text(face="plain", colour="black", size=16)) +  
  theme(axis.title.y=element_text(face="plain", colour="black", size=16)) +
  theme(axis.title.x=element_text(face="plain", colour="black", size=16)) +
  scale_y_continuous(labels=function(x) sprintf("%.0f", x), breaks=seq(0,6,by=2), limits=c(0,6)) +
  scale_x_continuous(labels=function(x) sprintf("%.0f", x), breaks=seq(0,10,by=2), limits=c(0,10)) +
  scale_color_manual(values=c("Chlamydomonas noctigama"="dodgerblue3", "Chlamydomonas bilatus"="firebrick3"))+
  scale_linetype_manual(values=c("Without"="11", "With"="solid")) +
  theme(axis.line=element_line(colour="black")) + theme(panel.background=element_blank()) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  theme(strip.text.x=element_text(face="bold", color="black", size=16)) +
  theme(strip.background=element_blank()) +
  facet_wrap(~Genotype, scales="free", ncol=7, nrow=7) +
  theme(plot.margin=unit(c(5.5,5.5,5.5,5.5),"pt")) +
  theme(legend.position="none")
dev.off()

tiff('Fitness Curves Rotifer High.tiff', units="in", width=15, height=15, res=1000)
ggplot(subset(Data2, Nutri=="High"), aes(Day, log(MeanDensR+1), group=Pred)) +
  geom_smooth(aes(color=Species, linetype=Pred), size=1, se=F) +
  geom_point(data=subset(Data, Nutri=="High"), aes(Day, log(DensR+1), color=Species), alpha=0.4, size=1.5, pch=16, position=position_jitter(0.2)) +
  ylab(expression('Rotifer density'~'('*ln~individuals~mL^-1*')')) + xlab(expression('Time (days)')) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=16)) +  
  theme(axis.text.x=element_text(face="plain", colour="black", size=16)) +  
  theme(axis.title.y=element_text(face="plain", colour="black", size=16)) +
  theme(axis.title.x=element_text(face="plain", colour="black", size=16)) +
  scale_y_continuous(labels=function(x) sprintf("%.0f", x), breaks=seq(0,6,by=2), limits=c(0,6)) +
  scale_x_continuous(labels=function(x) sprintf("%.0f", x), breaks=seq(0,10,by=2), limits=c(0,10)) +
  scale_color_manual(values=c("Chlamydomonas noctigama"="dodgerblue3", "Chlamydomonas bilatus"="firebrick3"))+
  scale_linetype_manual(values=c("Without"="11", "With"="solid")) +
  theme(axis.line=element_line(colour="black")) + theme(panel.background=element_blank()) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  theme(strip.text.x=element_text(face="bold", color="black", size=16)) +
  theme(strip.background=element_blank()) +
  facet_wrap(~Genotype, scales="free", ncol=7, nrow=7) +
  theme(plot.margin=unit(c(5.5,5.5,5.5,5.5),"pt")) +
  theme(legend.position="none")
dev.off()


###############################
### Calculating the fitness ###
###############################

### Fitness ###

# Subset the dataset
Data3=subset(Data, Day < 15)

# Calculate density time lags
Data3=Data3[order(Data3$Genotype,Data3$Nutri,Data3$Pred,Data3$Trial,Data3$Day),]
Data3=Data3 %>% group_by(Code,Genotype,Species,Period,Nutri,Pred,Trial) %>% mutate(LagA=c(rep(NA,1),diff(log(DensA),lag=1)/diff(Day,lag=1)), LagR=c(rep(NA,1),diff(log(DensR),lag=1)/diff(Day,lag=1)))
Data3=Data3 %>% mutate(DensAs=scale(DensA)[,1], DensRs=scale(DensR)[,1])

# Remove missing values
Data3$DensAs[is.na(Data3$DensAs)]=0
Data3$DensRs[is.na(Data3$DensRs)]=0
Data3$LagA[is.infinite(Data3$LagA)]=NA
Data3$LagR[is.infinite(Data3$LagR)]=NA

# Linear mixed-effect regression
ModA=lmer(LagA ~ DensAs + DensRs * Nutri * Pred + (1|Genotype) + (1|Trial), data=Data3)
summary(ModA)
Anova(ModA)

ModR=lmer(LagR ~ DensAs + DensRs * Nutri + (1|Genotype) + (1|Trial), data=Data3)
summary(ModR)
Anova(ModR)

# Model predictions
PredA=predictInterval(ModA, newdata=Data3, level=0.95, n.sims=1000, which="full", include.resid.var=T)
PredR=predictInterval(ModR, newdata=Data3, level=0.95, n.sims=1000, which="full", include.resid.var=T)
Data3$GrowA=PredA[,1]; Data3$GrowALSD=PredA[,2]; Data3$GrowAUSD=PredA[,3]
Data3$GrowR=PredR[,1]; Data3$GrowRLSD=PredR[,2]; Data3$GrowRUSD=PredR[,3]

# Replace null fitness
Data3$GrowR=ifelse(Data3$Pred=="Without", 0, Data3$GrowR)

# Calculate the fitness
Data4=setDT(Data3)[, .(GrowA=mean(GrowA), GrowR=mean(GrowR)),by=list(Genotype,Species,Period,Nutri,Pred,Trial)]
Data4=as.data.frame(Data4)

# Calculate the areas under curve
Data5=setDT(subset(Data, Day < 15))[,.(AucA=as.numeric(round(AUC(Day, log(DensA+1), from=min(Day,na.rm=T), to=max(Day,na.rm=T), method=c("spline"), absolutearea=F, subdivisions=1000, na.rm=T),0)), AucR=as.numeric(round(AUC(Day, log(DensR+1), from=min(Day,na.rm=T), to=max(Day,na.rm=T), method=c("spline"), absolutearea=F, subdivisions=1000, na.rm=T),0))), by=list(Genotype,Species,Period,Nutri,Pred,Trial)]
Data5=as.data.frame(Data5)

# Combine the datasets
Data6=cbind(Data4[,c(1:8)],Data5[,c(7:8)])

### Error bars ###

# Calculate the error bars
Data7=setDT(Data6)[, .(GrowA=round(mean(GrowA),6), GrowALSD=round(mean(GrowA)-sd(GrowA),6), GrowAUSD=round(mean(GrowA)+sd(GrowA),6), GrowR=round(mean(GrowR),6), GrowRLSD=round(mean(GrowR)-sd(GrowR),6), GrowRUSD=round(mean(GrowR)+sd(GrowR),6), AucA=round(mean(AucA),0), AucALSD=round(mean(AucA)-sd(AucA),0), AucAUSD=round(mean(AucA)+sd(AucA),0), AucR=round(mean(AucR),0), AucRLSD=round(mean(AucR)-sd(AucR),0), AucRUSD=round(mean(AucR)+sd(AucR),0)), by=list(Genotype,Species,Nutri,Pred)]
Data7=as.data.frame(Data7)

# Calculate the error bars for growth rates
Data8=as.data.frame(setDT(na.omit(Data7))[, .(GrowASpecies=round(mean(GrowA),6), GrowASpeciesL=round(mean(GrowA)-sd(GrowA),6), GrowASpeciesU=round(mean(GrowA)+sd(GrowA),6), GrowRSpecies=round(mean(GrowR),6), GrowRSpeciesL=round(mean(GrowR)-sd(GrowR),6), GrowRSpeciesU=round(mean(GrowR)+sd(GrowR),6)), by=list(Species,Nutri,Pred)])
Data7$GrowASpecies=c(rep(Data8$GrowASpecies[c(1:4)],24), rep(Data8$GrowASpecies[c(5:8)],22))
Data7$GrowASpeciesL=c(rep(Data8$GrowASpeciesL[c(1:4)],24), rep(Data8$GrowASpeciesL[c(5:8)],22))
Data7$GrowASpeciesU=c(rep(Data8$GrowASpeciesU[c(1:4)],24), rep(Data8$GrowASpeciesU[c(5:8)],22))
Data7$GrowRSpecies=c(rep(Data8$GrowRSpecies[c(1:4)],24), rep(Data8$GrowRSpecies[c(5:8)],22))
Data7$GrowRSpeciesL=c(rep(Data8$GrowRSpeciesL[c(1:4)],24), rep(Data8$GrowRSpeciesL[c(5:8)],22))
Data7$GrowRSpeciesU=c(rep(Data8$GrowRSpeciesU[c(1:4)],24), rep(Data8$GrowRSpeciesU[c(5:8)],22))

# Calculate the error bars for area under curve
Data9=as.data.frame(setDT(na.omit(Data7))[, .(AucASpecies=round(mean(AucA),0), AucASpeciesL=round(mean(AucA)-sd(AucA),0), AucASpeciesU=round(mean(AucA)+sd(AucA),0), AucRSpecies=round(mean(AucR),0), AucRSpeciesL=round(mean(AucR)-sd(AucR),0), AucRSpeciesU=round(mean(AucR)+sd(AucR),0)), by=list(Species,Nutri,Pred)])
Data7$AucASpecies=c(rep(Data9$AucASpecies[c(1:4)],24), rep(Data9$AucASpecies[c(5:8)],22))
Data7$AucASpeciesL=c(rep(Data9$AucASpeciesL[c(1:4)],24), rep(Data9$AucASpeciesL[c(5:8)],22))
Data7$AucASpeciesU=c(rep(Data9$AucASpeciesU[c(1:4)],24), rep(Data9$AucASpeciesU[c(5:8)],22))
Data7$AucRSpecies=c(rep(Data9$AucRSpecies[c(1:4)],24), rep(Data9$AucRSpecies[c(5:8)],22))
Data7$AucRSpeciesL=c(rep(Data9$AucRSpeciesL[c(1:4)],24), rep(Data9$AucRSpeciesL[c(5:8)],22))
Data7$AucRSpeciesU=c(rep(Data9$AucRSpeciesU[c(1:4)],24), rep(Data9$AucRSpeciesU[c(5:8)],22))

# Export the dataset
Data7[,c(6:10,17:22)]=round(Data7[,c(6:10,17:22)],4)
write.table(Data7[,c(1:28)], file="Data_FP.txt", sep="\t", row.names=F)


########################
### Plot the fitness ###
########################

# Combine the treatments
Data6$Treat=paste(Data6$Nutri, Data6$Pred, sep="-")
Data7$Treat=paste(Data7$Nutri, Data7$Pred, sep="-")
Data7$Treat=paste(Data7$Treat, Data7$Species, sep="-")
Data7$Treat=gsub("Chlamydomonas bilatus", "CB", Data7$Treat)
Data7$Treat=gsub("Chlamydomonas noctigama", "CN", Data7$Treat)
Data6$Treat=factor(Data6$Treat, levels=c("Low-Without","Low-With","High-Without","High-With"))
Data7$Treat=factor(Data7$Treat, levels=c("Low-Without-CN","Low-Without-CB","Low-With-CN","Low-With-CB","High-Without-CN","High-Without-CB","High-With-CN","High-With-CB"))


### Growth rate ###

# Calculate the significance between species
Stats1=subset(Data6, Pred %in% c("Without","With")) %>% group_by(Treat) %>% wilcox_test(GrowA~Species) %>% adjust_pvalue(method="bonferroni") %>% add_significance("p.adj")
Stats1=cbind(Stats1[,c(1,9:10)], x=c(1.00,3.00,5.00,7.00), xend=c(2.00,4.00,6.00,8.00), y=rep(0.3+(0.020*0.3),4), yend=rep(0.3+(0.020*0.3),4))
colnames(Stats1)[c(1:3)]=c("Trait","P","S")

# Calculate the significance between predator
Stats2=subset(Data6, Pred %in% c("Without","With")) %>% group_by(Nutri) %>% wilcox_test(GrowA~Pred) %>% adjust_pvalue(method="bonferroni") %>% add_significance("p.adj")
Stats2=cbind(Stats2[,c(1,9:10)], x=c(1.50,5.50), xend=c(3.50,7.50), y=rep(0.3+(0.145*0.3),2), yend=rep(0.3+(0.145*0.3),2))
colnames(Stats2)[c(1:3)]=c("Trait","P","S")

# Calculate the significance between nutrient
Stats3=subset(Data6, Pred %in% c("Without","With")) %>% group_by(Pred) %>% wilcox_test(GrowA~Nutri) %>% adjust_pvalue(method="bonferroni") %>% add_significance("p.adj")
Stats3=cbind(Stats3[,c(1,9:10)], x=c(2.50,2.50), xend=c(6.50,6.50), y=rep(0.3+(0.270*0.3),1), yend=rep(0.3+(0.270*0.3),1))
colnames(Stats3)[c(1:3)]=c("Trait","P","S")

tiff('Fitness Alga.tiff', units="in", width=8, height=8, res=1000)
ggplot(subset(Data7, Pred %in% c("Without","With")), aes(Treat, GrowA)) + coord_cartesian(clip="off") +
  geom_pointrange(aes(Treat, ymin=GrowALSD, ymax=GrowAUSD, color=Species), size=0.8, alpha=0.5, pch=16, linewidth=1, linetype="solid", position=position_jitter(w=0.3)) +
  geom_segment(data=Stats1, aes(x=x, xend=xend, y=y, yend=yend), color="black", size=0.8) +
  geom_segment(data=Stats2, aes(x=x, xend=xend, y=y, yend=yend), color="black", size=0.8) +
  geom_segment(data=Stats3, aes(x=x, xend=xend, y=y, yend=yend), color="black", size=0.8) +
  geom_text(data=Stats1, aes(x=x+0.50, y=0.3+(0.070*0.3), label=S), color="black", size=5) +
  geom_text(data=Stats2, aes(x=x+1.00, y=0.3+(0.195*0.3), label=S), color="black", size=5) +
  geom_text(data=Stats3, aes(x=x+2.00, y=0.3+(0.320*0.3), label=S), color="black", size=5) +
  ylab(expression('Fitness'~'('*day^-1*')')) +
  annotation_custom(grob=textGrob(label="No predation", gp=gpar(cex=1.5, color="black")), xmin=1.5, xmax=1.5, ymin=-0.1367, ymax=-0.1367) +
  annotation_custom(grob=textGrob(label="Predation", gp=gpar(cex=1.5, color="black")), xmin=3.5, xmax=3.5, ymin=-0.1367, ymax=-0.1367) +
  annotation_custom(grob=textGrob(label="No predation", gp=gpar(cex=1.5, color="black")), xmin=5.5, xmax=5.5, ymin=-0.1367, ymax=-0.1367) +
  annotation_custom(grob=textGrob(label="Predation", gp=gpar(cex=1.5, color="black")), xmin=7.5, xmax=7.5, ymin=-0.1367, ymax=-0.1367) +
  annotation_custom(grob=textGrob(label="Low nutrient", gp=gpar(cex=1.5, color="black")), xmin=2.5, xmax=2.5, ymin=-0.1570, ymax=-0.1570) +
  annotation_custom(grob=textGrob(label="High nutrient", gp=gpar(cex=1.5, color="black")), xmin=6.5, xmax=6.5, ymin=-0.1570, ymax=-0.1570) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.text.x=element_text(face="plain", colour=NA, size=18)) +  
  theme(axis.title.y=element_text(face="plain", colour="black", size=18)) +
  theme(axis.title.x=element_text(face="plain", colour=NA, size=18)) +
  scale_y_continuous(labels=function(x) sprintf("%.1f", x), breaks=seq(-0.1,0.3,by=0.1), limits=c(-0.1,0.3+(0.32*0.3))) +
  scale_x_discrete(labels=c("Low-Without-CN"="No predation", "Low-Without-CB"="No predation", "Low-With-CN"="Predation", "Low-With-CB"="Predation", "High-Without-CN"="No predation", "High-Without-CB"="No predation", "High-With-CN"="Predation", "High-With-CB"="Predation")) +
  scale_fill_manual(values=alpha(c("Chlamydomonas noctigama"="dodgerblue3", "Chlamydomonas bilatus"="firebrick3"),0.3)) +
  scale_color_manual(values=c("Chlamydomonas noctigama"="dodgerblue3", "Chlamydomonas bilatus"="firebrick3")) +
  theme(axis.line=element_line(colour="black")) + theme(panel.background=element_blank()) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  theme(plot.margin=unit(c(5.5,5.5,5.5,5.5),"pt")) +
  theme(legend.position="none")
dev.off()

# Calculate the significance between species
Stats1=subset(Data6, Pred %in% c("With")) %>% group_by(Treat) %>% wilcox_test(GrowR~Species) %>% adjust_pvalue(method="bonferroni") %>% add_significance("p.adj")
Stats1=cbind(Stats1[,c(1,9:10)], x=c(1.00,3.00), xend=c(2.00,4.00), y=rep(0.2+(0.020*0.2),2), yend=rep(0.2+(0.020*0.2),2))
colnames(Stats1)[c(1:3)]=c("Trait","P","S")

# Calculate the significance between nutrient
Stats2=subset(Data6, Pred %in% c("With")) %>% group_by(Pred) %>% wilcox_test(GrowA~Nutri) %>% adjust_pvalue(method="bonferroni") %>% add_significance("p.adj")
Stats2=cbind(Stats2[,c(1,9:10)], x=c(1.50), xend=c(3.50), y=rep(0.2+(0.145*0.2),1), yend=rep(0.2+(0.145*0.2),1))
colnames(Stats2)[c(1:3)]=c("Trait","P","S")

tiff('Fitness Rotifer.tiff', units="in", width=8, height=8, res=1000)
ggplot(subset(Data7, Pred %in% c("With")), aes(Treat, GrowR)) + coord_cartesian(clip="off") +
  geom_pointrange(aes(Treat, ymin=GrowRLSD, ymax=GrowRUSD, color=Species), size=0.8, alpha=0.5, pch=16, linewidth=1, linetype="solid", position=position_jitter(w=0.3)) +
  geom_segment(data=Stats1, aes(x=x, xend=xend, y=y, yend=yend), color="black", size=0.8) +
  geom_segment(data=Stats2, aes(x=x, xend=xend, y=y, yend=yend), color="black", size=0.8) +
  geom_text(data=Stats1, aes(x=x+0.50, y=0.2+(0.070*0.2), label=S), color="black", size=5) +
  geom_text(data=Stats2, aes(x=x+1.00, y=0.2+(0.195*0.2), label=S), color="black", size=5) +
  ylab(expression('Fitness'~'('*day^-1*')')) +
  annotation_custom(grob=textGrob(label="Low nutrient", gp=gpar(cex=1.5, color="black")), xmin=1.5, xmax=1.5, ymin=-0.2325, ymax=-0.2325) +
  annotation_custom(grob=textGrob(label="High nutrient", gp=gpar(cex=1.5, color="black")), xmin=3.5, xmax=3.5, ymin=-0.2325, ymax=-0.2325) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.text.x=element_text(face="plain", colour=NA, size=18)) +  
  theme(axis.title.y=element_text(face="plain", colour="black", size=18)) +
  theme(axis.title.x=element_text(face="plain", colour=NA, size=18)) +
  scale_y_continuous(labels=function(x) sprintf("%.1f", x), breaks=seq(-0.2,0.2,by=0.1), limits=c(-0.2,0.2+(0.195*0.2))) +
  scale_x_discrete(labels=c("Low-Without"="No predation", "Low-With"="Predation", "High-Without"="No predation", "High-With"="Predation")) +
  scale_fill_manual(values=alpha(c("Chlamydomonas noctigama"="dodgerblue3", "Chlamydomonas bilatus"="firebrick3"),0.3)) +
  scale_color_manual(values=c("Chlamydomonas noctigama"="dodgerblue3", "Chlamydomonas bilatus"="firebrick3")) +
  theme(axis.line=element_line(colour="black")) + theme(panel.background=element_blank()) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  theme(plot.margin=unit(c(5.5,5.5,5.5,5.5),"pt")) +
  theme(legend.position="none")
dev.off()

### Areas under curve ###

# Calculate the significance between species
Stats1=Data6 %>% group_by(Treat) %>% wilcox_test(AucA~Species) %>% adjust_pvalue(method="bonferroni") %>% add_significance("p.adj")
Stats1=cbind(Stats1[,c(1,9:10)], x=c(1.00,3.00,5.00,7.00), xend=c(2.00,4.00,6.00,8.00), y=rep(125+(0.020*37.5),4), yend=rep(125+(0.020*37.5),4))
colnames(Stats1)[c(1:3)]=c("Trait","P","S")

# Calculate the significance between predator
Stats2=Data6 %>% group_by(Nutri) %>% wilcox_test(AucA~Pred) %>% adjust_pvalue(method="bonferroni") %>% add_significance("p.adj")
Stats2=cbind(Stats2[,c(1,9:10)], x=c(1.50,5.50), xend=c(3.50,7.50), y=rep(125+(0.145*37.5),2), yend=rep(125+(0.145*37.5),2))
colnames(Stats2)[c(1:3)]=c("Trait","P","S")

# Calculate the significance between nutrient
Stats3=Data6 %>% group_by(Pred) %>% wilcox_test(AucA~Nutri) %>% adjust_pvalue(method="bonferroni") %>% add_significance("p.adj")
Stats3=cbind(Stats3[,c(1,9:10)], x=c(2.50), xend=c(6.50), y=rep(125+(0.270*37.5),1), yend=rep(125+(0.270*37.5),1))
colnames(Stats3)[c(1:3)]=c("Trait","P","S")

tiff('Fitness Alga.tiff', units="in", width=8, height=8, res=1000)
ggplot(subset(Data7, Pred %in% c("Without","With")), aes(Treat, AucA)) + coord_cartesian(clip="off") +
  geom_pointrange(aes(Treat, ymin=AucALSD, ymax=AucAUSD, color=Species), size=0.8, alpha=0.5, pch=16, linewidth=1, linetype="solid", position=position_jitter(w=0.3)) +
  geom_segment(data=Stats1, aes(x=x, xend=xend, y=y, yend=yend), color="black", size=0.8) +
  geom_segment(data=Stats2, aes(x=x, xend=xend, y=y, yend=yend), color="black", size=0.8) +
  geom_segment(data=Stats3, aes(x=x, xend=xend, y=y, yend=yend), color="black", size=0.8) +
  geom_text(data=Stats1, aes(x=x+0.50, y=125+(0.050*37.5), label=S), color="black", size=5) +
  geom_text(data=Stats2, aes(x=x+1.00, y=125+(0.175*37.5), label=S), color="black", size=5) +
  geom_text(data=Stats3, aes(x=x+2.00, y=125+(0.300*37.5), label=S), color="black", size=5) +
  ylab(expression('Fitness'~'('*ln~cells~day^-1*')')) +
  annotation_custom(grob=textGrob(label="No predation", gp=gpar(cex=1.5, color="black")), xmin=1.5, xmax=1.5, ymin=43.5, ymax=43.5) +
  annotation_custom(grob=textGrob(label="Predation", gp=gpar(cex=1.5, color="black")), xmin=3.5, xmax=3.5, ymin=43.5, ymax=43.5) +
  annotation_custom(grob=textGrob(label="No predation", gp=gpar(cex=1.5, color="black")), xmin=5.5, xmax=5.5, ymin=43.5, ymax=43.5) +
  annotation_custom(grob=textGrob(label="Predation", gp=gpar(cex=1.5, color="black")), xmin=7.5, xmax=7.5, ymin=43.5, ymax=43.5) +
  annotation_custom(grob=textGrob(label="Low nutrient", gp=gpar(cex=1.5, color="black")), xmin=2.5, xmax=2.5, ymin=39.5, ymax=39.5) +
  annotation_custom(grob=textGrob(label="High nutrient", gp=gpar(cex=1.5, color="black")), xmin=6.5, xmax=6.5, ymin=39.5, ymax=39.5) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.text.x=element_text(face="plain", colour=NA, size=18)) +  
  theme(axis.title.y=element_text(face="plain", colour="black", size=18)) +
  theme(axis.title.x=element_text(face="plain", colour=NA, size=18)) +
  scale_y_continuous(labels=function(x) sprintf("%.0f", x), breaks=seq(50,125,by=25), limits=c(50,125+(0.3*37.5))) +
  scale_x_discrete(labels=c("Low-Without-CN"="No predation", "Low-Without-CB"="No predation", "Low-With-CN"="Predation", "Low-With-CB"="Predation", "High-Without-CN"="No predation", "High-Without-CB"="No predation", "High-With-CN"="Predation", "High-With-CB"="Predation")) +
  scale_fill_manual(values=alpha(c("Chlamydomonas noctigama"="dodgerblue3", "Chlamydomonas bilatus"="firebrick3"),0.3)) +
  scale_color_manual(values=c("Chlamydomonas noctigama"="dodgerblue3", "Chlamydomonas bilatus"="firebrick3")) +
  theme(axis.line=element_line(colour="black")) + theme(panel.background=element_blank()) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  theme(plot.margin=unit(c(5.5,5.5,5.5,5.5),"pt")) +
  theme(legend.position="none")
dev.off()

# Calculate the significance between species
Stats1=subset(Data6, Pred %in% c("With")) %>% group_by(Treat) %>% wilcox_test(AucR~Species) %>% adjust_pvalue(method="bonferroni") %>% add_significance("p.adj")
Stats1=cbind(Stats1[,c(1,9:10)], x=c(1.00,3.00), xend=c(2.00,4.00), y=rep(125+(0.020*37.5),2), yend=rep(125+(0.020*37.5),2))
colnames(Stats1)[c(1:3)]=c("Trait","P","S")

# Calculate the significance between nutrient
Stats2=subset(Data6, Pred %in% c("With")) %>% group_by(Pred) %>% wilcox_test(AucR~Nutri) %>% adjust_pvalue(method="bonferroni") %>% add_significance("p.adj")
Stats2=cbind(Stats2[,c(1,9:10)], x=c(1.50), xend=c(3.50), y=rep(125+(0.145*37.5),1), yend=rep(125+(0.145*37.5),1))
colnames(Stats2)[c(1:3)]=c("Trait","P","S")

tiff('Fitness Rotifer.tiff', units="in", width=8, height=8, res=1000)
ggplot(subset(Data7, Pred %in% c("With")), aes(Treat, AucA)) + coord_cartesian(clip="off") +
  geom_pointrange(aes(Treat, ymin=AucALSD, ymax=AucAUSD, color=Species), size=0.8, alpha=0.5, pch=16, linewidth=1, linetype="solid", position=position_jitter(w=0.3)) +
  geom_segment(data=Stats1, aes(x=x, xend=xend, y=y, yend=yend), color="black", size=0.8) +
  geom_segment(data=Stats2, aes(x=x, xend=xend, y=y, yend=yend), color="black", size=0.8) +
  geom_text(data=Stats1, aes(x=x+0.50, y=125+(0.050*37.5), label=S), color="black", size=5) +
  geom_text(data=Stats2, aes(x=x+1.00, y=125+(0.175*37.5), label=S), color="black", size=5) +
  ylab(expression('Fitness'~'('*ln~cells~day^-1*')')) +
  annotation_custom(grob=textGrob(label="Low nutrient", gp=gpar(cex=1.5, color="black")), xmin=1.5, xmax=1.5, ymin=39.5, ymax=39.5) +
  annotation_custom(grob=textGrob(label="High nutrient", gp=gpar(cex=1.5, color="black")), xmin=3.5, xmax=3.5, ymin=39.5, ymax=39.5) +
  theme(axis.text.y=element_text(face="plain", colour="black", size=18)) +  
  theme(axis.text.x=element_text(face="plain", colour=NA, size=18)) +  
  theme(axis.title.y=element_text(face="plain", colour="black", size=18)) +
  theme(axis.title.x=element_text(face="plain", colour=NA, size=18)) +
  scale_y_continuous(labels=function(x) sprintf("%.0f", x), breaks=seq(50,125,by=25), limits=c(50,125+(0.3*37.5))) +
  scale_x_discrete(labels=c("Low-Without"="No predation", "Low-With"="Predation", "High-Without"="No predation", "High-With"="Predation")) +
  scale_fill_manual(values=alpha(c("Chlamydomonas noctigama"="dodgerblue3", "Chlamydomonas bilatus"="firebrick3"),0.3)) +
  scale_color_manual(values=c("Chlamydomonas noctigama"="dodgerblue3", "Chlamydomonas bilatus"="firebrick3")) +
  theme(axis.line=element_line(colour="black")) + theme(panel.background=element_blank()) +
  theme(panel.grid.major=element_blank(), panel.grid.minor=element_blank()) +
  theme(plot.margin=unit(c(5.5,5.5,5.5,5.5),"pt")) +
  theme(legend.position="none")
dev.off()


##################
### Statistics ###
##################

### Growth rate ###

# Rescale the data
Data4$GrowA[Data4$GrowA<=0]=0.001
Data4$GrowR[Data4$GrowR<=0]=0.001

# Subset the dataset
DataB=subset(Data4, Species=="Chlamydomonas bilatus")
DataN=subset(Data4, Species=="Chlamydomonas noctigama")

# Test for normality
shapiro.test(Data4$GrowA)
qqnorm(Data4$GrowA)
qqline(Data4$GrowA)
hist(Data4$GrowA)

# Test for homoscedasticity
bartlett.test(Data4$GrowA, Data4$Species)
bartlett.test(Data4$GrowA, Data4$Nutri)
bartlett.test(Data4$GrowA, Data4$Pred)

# Generalized linear regressions
summary(aov(glm(GrowA~Nutri + Pred + Nutri:Pred, family=Gamma(link="log"), data=Data4)))
summary(aov(glm(GrowA~Nutri + Pred + Nutri:Pred, family=Gamma(link="log"), data=DataB)))
summary(aov(glm(GrowA~Nutri + Pred + Nutri:Pred, family=Gamma(link="log"), data=DataN)))

# Test for normality
shapiro.test(Data4$GrowR)
qqnorm(Data4$GrowR)
qqline(Data4$GrowR)
hist(Data4$GrowR)

# Test for homoscedasticity
bartlett.test(Data4$GrowR, Data4$Species)
bartlett.test(Data4$GrowR, Data4$Nutri)

# Generalized linear regressions
summary(aov(glm(GrowR~Nutri, family=Gamma(link="log"), data=Data4)))
summary(aov(glm(GrowR~Nutri, family=Gamma(link="log"), data=DataB)))
summary(aov(glm(GrowR~Nutri, family=Gamma(link="log"), data=DataN)))

### Areas under curve ###

# Rescale the data
Data3$AucA=Data3$AucA+0.001
Data3$AucR=Data3$AucR+0.001

# Subset the dataset
DataB=subset(Data3, Species=="Chlamydomonas bilatus")
DataN=subset(Data3, Species=="Chlamydomonas noctigama")

# Test for normality
shapiro.test(log(Data3$AucA+1))
qqnorm(log(Data3$AucA+1))
qqline(log(Data3$AucA+1))
hist(log(Data3$AucA+1))

# Test for homoscedasticity
bartlett.test(log(Data3$AucA+1), Data3$Species)
bartlett.test(log(Data3$AucA+1), Data3$Nutri)
bartlett.test(log(Data3$AucA+1), Data3$Pred)

# Generalized linear regressions
summary(aov(glm(log(AucA+1)~Nutri + Pred + Nutri:Pred, family=Gamma(link="log"), data=Data3)))
summary(aov(glm(log(AucA+1)~Nutri + Pred + Nutri:Pred, family=Gamma(link="log"), data=DataB)))
summary(aov(glm(log(AucA+1)~Nutri + Pred + Nutri:Pred, family=Gamma(link="log"), data=DataN)))

# Test for normality
shapiro.test(log(Data3$AucR+1))
qqnorm(log(Data3$AucR+1))
qqline(log(Data3$AucR+1))
hist(log(Data3$AucR+1))

# Test for homoscedasticity
bartlett.test(log(Data3$AucR+1), Data3$Species)
bartlett.test(log(Data3$AucR+1), Data3$Nutri)

# Generalized linear regressions
summary(aov(glm(log(AucR+1)~Nutri, family=Gamma(link="log"), data=Data3)))
summary(aov(glm(log(AucR+1)~Nutri, family=Gamma(link="log"), data=DataB)))
summary(aov(glm(log(AucR+1)~Nutri, family=Gamma(link="log"), data=DataN)))

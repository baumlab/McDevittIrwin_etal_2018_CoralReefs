# Combine permadisp plots for MS: High and low heat stress together

setwd("/Users/jamiemcdevitt-irwin/Documents/Git_Repos/McDevittIrwinetal_Kiritimati16S")
getwd()
library("ggplot2")
library("vegan")
library("phyloseq")

# clear my environment
rm(list=ls()) 
theme_set(theme_bw())

# load may
load("data/secondmito/lowreads_pruned/may_mont_boxplot.Rdata")
load("data/secondmito/lowreads_pruned/may_por_boxplot.Rdata")
ls()

# load july
load("data/secondmito/lowreads_pruned/july_por_boxplot.Rdata")
load("data/secondmito/lowreads_pruned/july_mont_boxplot.Rdata")
ls()

load("data/secondmito/lowreads_pruned/otutable_top10removed_coral_866_may_rel_por.Rdata") #otu_tablef_no10_coralf_sm_866f_mayf_rel_porf
load("data/secondmito/lowreads_pruned/otutable_top10removed_coral_866_may_rel_mont.Rdata") #otu_tablef_no10_coralf_sm_866f_mayf_rel_montf
load("data/secondmito/lowreads_pruned/otutable_top10removed_coral_866_rel_july_porites.Rdata") # otu_tablef_no10_coralf_sm_866f_mayf_rel_porf
load("data/secondmito/lowreads_pruned/otutable_top10removed_coral_866_rel_july_montipora.Rdata") # otu_tablef_no10_coralf_sm_866f_mayf_rel_montf
ls()
################################################################################


################################################################################
# May Porites 

sampledf <- data.frame(sample_data(otu_tablef_no10_coralf_sm_866f_mayf_rel_porf))

# Extract mean and SD (same as Anderson paper) so we don't have to use boxplots 
maypor_mean <- tapply(maypor_beta$distances, sampledf$human_disturbance, mean) # get the means
maypor_sd <- tapply(maypor_beta$distances, sampledf$human_disturbance, sd) # get the SD
maypor_sd <- as.data.frame(maypor_sd)
maypor_sd
maypor_mean <- as.data.frame(maypor_mean)
maypor_mean

colnames(maypor_mean)[1] <- "mean"
colnames(maypor_sd)[1] <- "sd"

maypor_plot <- merge(maypor_sd, maypor_mean, by="row.names") 
maypor_plot$Row.names
maypor_plot$mean <-as.numeric(maypor_plot$mean)
maypor_plot
maypor_plot$Row.names <- as.factor(maypor_plot$Row.names)
maypor_plot$sd <-as.numeric(maypor_plot$sd)
maypor_plot

# change the level so low is first
maypor_plot$Row.names <- factor(maypor_plot$Row.names, levels = c("Very low", "Very high"))
maypor_plot

# Porites Dispersion- May
p= ggplot(maypor_plot, aes(x=Row.names, y=mean, col=Row.names, size=1)) + geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.1, size=1) + geom_point(shape=17)
p=p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none", axis.text.x=element_blank(),axis.text=element_text(size=12),axis.title=element_text(size=14)) 
p= p+ labs(x="") +labs(y="Distance to Centroid")
p= p + scale_colour_manual(values=c("darkgoldenrod1","darkorange2")) + scale_y_continuous(limits=c(0,0.9))
p
################################################################################



################################################################################
# May Montipora
sampledf <- data.frame(sample_data(otu_tablef_no10_coralf_sm_866f_mayf_rel_montf))

# Extract mean and SD (same as Anderson paper) so we don't have to use boxplots 
maymont_mean <- tapply(maymont_beta$distances, sampledf$human_disturbance, mean) # get the means
maymont_sd <- tapply(maymont_beta$distances, sampledf$human_disturbance, sd) # get the SD
maymont_sd <- as.data.frame(maymont_sd)
maymont_sd
maymont_mean <- as.data.frame(maymont_mean)
maymont_mean

colnames(maymont_mean)[1] <- "mean"
colnames(maymont_sd)[1] <- "sd"

maymont_plot <- merge(maymont_sd, maymont_mean, by="row.names") 
maymont_plot$Row.names
maymont_plot$mean <-as.numeric(maymont_plot$mean)
maymont_plot
maymont_plot$Row.names <- as.factor(maymont_plot$Row.names)
maymont_plot$sd <-as.numeric(maymont_plot$sd)

# change the level so low is first
maymont_plot$Row.names <- factor(maymont_plot$Row.names, levels = c("Very low", "Very high"))

# montipora Dispersion- May
p1= ggplot(maymont_plot, aes(x=Row.names, y=mean, col=Row.names, size=1)) + geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.1, size=1) + geom_point()
p1=p1 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none",axis.text=element_text(size=12),axis.title=element_text(size=14)) 
p1=p1+ labs(x="Disturbance") +labs(y="Distance to Centroid") 
p1= p1 + scale_colour_manual(values=c("darkgoldenrod1","darkorange2")) + scale_y_continuous(limits=c(0,0.9))+ scale_x_discrete(breaks=c("Very high", "Very low"), labels=c("High", "Low"))
p1
################################################################################


################################################################################
# July Porites
sampledf <- data.frame(sample_data(otu_tablef_no10_coralf_sm_866f_julyf_rel_porf))

# Extract mean and SD (same as Anderson paper) so we don't have to use boxplots 
julypor_mean <- tapply(julypor_beta$distances, sampledf$human_disturbance, mean) # get the means
julypor_sd <- tapply(julypor_beta$distances, sampledf$human_disturbance, sd) # get the SD
julypor_sd <- as.data.frame(julypor_sd)
julypor_sd
julypor_mean <- as.data.frame(julypor_mean)
julypor_mean

colnames(julypor_mean)[1] <- "mean"
colnames(julypor_sd)[1] <- "sd"

julypor_plot <- merge(julypor_sd, julypor_mean, by="row.names") 
julypor_plot$Row.names
julypor_plot$mean <-as.numeric(julypor_plot$mean)
julypor_plot
julypor_plot$Row.names <- as.factor(julypor_plot$Row.names)
julypor_plot$sd <-as.numeric(julypor_plot$sd)

# change the level so high is first
julypor_plot$Row.names <- factor(julypor_plot$Row.names, levels = c("Very low", "Very high"))

# July porites
p2= ggplot(julypor_plot, aes(x=Row.names, y=mean, col=Row.names, size=1)) + geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.1, size=1) + geom_point(shape=17)
p2=p2 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none",axis.text.x=element_blank(),axis.text=element_text(size=12),axis.title=element_text(size=14)) 
p2= p2+labs(x="") +labs(y="")
p2= p2 + scale_colour_manual(values=c("darkgoldenrod1","darkorange2")) + scale_y_continuous(limits=c(0,0.9))
p2
################################################################################


################################################################################
# July Montipora
sampledf <- data.frame(sample_data(otu_tablef_no10_coralf_sm_866f_julyf_rel_montf))

# Extract mean and SD (same as Anderson paper) so we don't have to use boxplots 
julymont_mean <- tapply(julymont_beta$distances, sampledf$human_disturbance, mean) # get the means
julymont_sd <- tapply(julymont_beta$distances, sampledf$human_disturbance, sd) # get the SD
julymont_sd <- as.data.frame(julymont_sd)
julymont_sd
julymont_mean <- as.data.frame(julymont_mean)
julymont_mean

colnames(julymont_mean)[1] <- "mean"
colnames(julymont_sd)[1] <- "sd"

julymont_plot <- merge(julymont_sd, julymont_mean, by="row.names") 
julymont_plot$Row.names
julymont_plot$mean <-as.numeric(julymont_plot$mean)
julymont_plot
julymont_plot$Row.names <- as.factor(julymont_plot$Row.names)
julymont_plot$sd <-as.numeric(julymont_plot$sd)

# change the level so high is first
julymont_plot$Row.names <- factor(julymont_plot$Row.names, levels = c("Very low", "Very high"))

# montimonta Dispersion- July
p3= ggplot(julymont_plot, aes(x=Row.names, y=mean, col=Row.names,size=1)) + geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.1,size=1) + geom_point()
p3=p3 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none",axis.text=element_text(size=12),axis.title=element_text(size=14)) 
p3= p3+ labs(x="Disturbance") +labs(y="")
p3= p3 + scale_colour_manual(values=c("darkgoldenrod1", "darkorange2")) + scale_y_continuous(limits=c(0,0.9)) +  scale_x_discrete(breaks=c("Very high", "Very low"), labels=c("High", "Low"))
p3
################################################################################




################################################################################
# combine the plots
p=p+ annotate("text",x=.91, y=0.88, label="(a) Porites lobata", fontface=2, size=4.5)+ annotate("text", x=2.5, y=0.88, label="*", size=6) +ggtitle("Low heat stress")
p1=p1 + annotate("text", x=1.32, y=0.88, label="(c) Montipora aequituberculata", fontface=2, size=4.5) 

p2=p2+ annotate("text", x=.55, y=0.88, label="(b)", fontface=2, size=4.5) + ggtitle("High heat stress")
p3=p3 + annotate("text",x=.55, y=0.88, label="(d)", fontface=2, size=4.5) 

# shared legend code
library(ggplot2)
library(gridExtra)
library(grid)

grid_arrange_shared_legend <- function(...) {
    plots <- list(...)
    g <- ggplotGrob(plots[[1]] + theme(legend.position="bottom"))$grobs
    legend <- g[[which(sapply(g, function(x) x$name) == "guide-box")]]
    lheight <- sum(legend$height)
    grid.arrange(
        do.call(arrangeGrob, lapply(plots, function(x)
            x + theme(legend.position="none"))),
        legend,
        ncol = 1,
        heights = unit.c(unit(1, "npc") - lheight, lheight))
}


# onefile=FALSE gets rid of the extra blank page, not sure why it was there in the first place. bug in the code?
# found onefile fix on github

pdf(file="figures/secondmito/top10removed/distance/mayandjuly_points_forms.pdf", height=8, width=8, onefile=FALSE)

grid.arrange(p, p2, p1, p3)

dev.off()





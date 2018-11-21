
# Beta Diversity: Low to high heat stress changes
# ***************************************************

setwd("/Users/jamiemcdevitt-irwin/Documents/Git_Repos/McDevittIrwinetal_Kiritimati16S")


# clear my environment
rm(list=ls()) 

# Installing & Importing  ----------------------------------------------------------
#First we need to install phyloseq
#source("http://bioconductor.org/biocLite.R")
#biocLite("phyloseq")


#now load phyloseq
library("phyloseq")
packageVersion("phyloseq")
library("ggplot2")
packageVersion("ggplot2")
library("scales")
packageVersion("scales")
library("vegan")

library("grid")
packageVersion("grid")
theme_set(theme_bw())


load("data/secondmito/lowreads_pruned/otutable_top10removed_coral_866_rel.Rdata") # otu_tablef_no10_coralf_sm_866f_rel
ls()
####################################################################################################


####################################################################################################
# SUBSET

# Montipora
otu_tablef_no10_coralf_sm_866f_rel_mont <- subset_samples(otu_tablef_no10_coralf_sm_866f_rel, field_host_genus_id=="Montipora")
print(otu_tablef_no10_coralf_sm_866f_rel_mont)
# 48 samples 
any(taxa_sums(otu_tablef_no10_coralf_sm_866f_rel_mont) == 0)
#TRUE because i filtered
otu_tablef_no10_coralf_sm_866f_rel_montf = prune_taxa(taxa_sums(otu_tablef_no10_coralf_sm_866f_rel_mont) > 0, otu_tablef_no10_coralf_sm_866f_rel_mont)
any(taxa_sums(otu_tablef_no10_coralf_sm_866f_rel_montf) == 0)
#FALSE
print(otu_tablef_no10_coralf_sm_866f_rel_montf) #48  

#save(otu_tablef_no10_coralf_sm_866f_rel_montf, file="data/secondmito/lowreads_pruned/otutable_top10removed_coral_866_rel_montipora.Rdata")



# Porites
otu_tablef_no10_coralf_sm_866f_rel_por <- subset_samples(otu_tablef_no10_coralf_sm_866f_rel, field_host_genus_id=="Porites")
print(otu_tablef_no10_coralf_sm_866f_rel_por)
# 55 samples 
any(taxa_sums(otu_tablef_no10_coralf_sm_866f_rel_por) == 0)
#TRUE because i filtered
otu_tablef_no10_coralf_sm_866f_rel_porf = prune_taxa(taxa_sums(otu_tablef_no10_coralf_sm_866f_rel_por) > 0, otu_tablef_no10_coralf_sm_866f_rel_por)
any(taxa_sums(otu_tablef_no10_coralf_sm_866f_rel_porf) == 0)
#FALSE
print(otu_tablef_no10_coralf_sm_866f_rel_porf) # 55

#save(otu_tablef_no10_coralf_sm_866f_rel_porf, file="data/secondmito/lowreads_pruned/otutable_top10removed_coral_866_rel_porites.Rdata")
##################################################



##################################################
# MONTIPORA HIGH DIST.

otu_tablef_no10_coralf_sm_866f_rel_mont_vhigh <- subset_samples(otu_tablef_no10_coralf_sm_866f_rel_mont, human_disturbance=="Very high")
print(otu_tablef_no10_coralf_sm_866f_rel_mont_vhigh)
# 24 samples
any(taxa_sums(otu_tablef_no10_coralf_sm_866f_rel_mont_vhigh) == 0)
#TRUE because i filtered
otu_tablef_no10_coralf_sm_866f_rel_mont_vhighf = prune_taxa(taxa_sums(otu_tablef_no10_coralf_sm_866f_rel_mont_vhigh) > 0, otu_tablef_no10_coralf_sm_866f_rel_mont_vhigh)
any(taxa_sums(otu_tablef_no10_coralf_sm_866f_rel_mont_vhighf) == 0)
#FALSE
print(otu_tablef_no10_coralf_sm_866f_rel_mont_vhighf) #24 

# save(otu_tablef_no10_coralf_sm_866f_rel_mont_vhighf, file="data/secondmito/lowreads_pruned/otutable_top10removed_coral_866_rel_montipora_vhigh.Rdata")

# rename the levels
as(sample_data(otu_tablef_no10_coralf_sm_866f_rel_mont_vhighf), "data.frame")[as(sample_data(otu_tablef_no10_coralf_sm_866f_rel_mont_vhighf), "data.frame")$human_disturbance %in% "Very low", "human_disturbance2"] <- "Low"
as(sample_data(otu_tablef_no10_coralf_sm_866f_rel_mont_vhighf), "data.frame")$human_disturbance2
as(sample_data(otu_tablef_no10_coralf_sm_866f_rel_mont_vhighf), "data.frame")[as(sample_data(otu_tablef_no10_coralf_sm_866f_rel_mont_vhighf), "data.frame")$human_disturbance %in% "Very high", "human_disturbance2"] <- "High"
as(sample_data(otu_tablef_no10_coralf_sm_866f_rel_mont_vhighf), "data.frame")$human_disturbance2


# Beta diversity
sample_data(otu_tablef_no10_coralf_sm_866f_rel_mont_vhighf)$human_disturbance2 # High
print(otu_tablef_no10_coralf_sm_866f_rel_mont_vhighf) # 24 samples 

monthigh_dist = phyloseq::distance(otu_tablef_no10_coralf_sm_866f_rel_mont_vhighf, "bray")
# Homogeneity of dispersion test 
sampledf <- data.frame(sample_data(otu_tablef_no10_coralf_sm_866f_rel_mont_vhighf))
monthigh_beta <- betadisper(monthigh_dist, sampledf$expedition_number)
plot(monthigh_beta, hull=FALSE, ellipse=TRUE) #showing 1 standard deviation ellipses 
t <-boxplot(monthigh_beta) # this is how you get out the lower and upper extremes of hte 'noth', this gives you almost equivalent to the 95% CI
permutest(monthigh_beta, perm=9999)

# Number of permutations: 9999

# Response: Distances
#           Df   Sum Sq    Mean Sq      F N.Perm Pr(>F)
# Groups     1 0.000043 0.00004263 0.0217   9999 0.8821
# Residuals 22 0.043159 0.00196180   



# Extract mean and SD (same as Anderson paper) so we don't have to use boxplots 
monthigh_mean <- tapply(monthigh_beta$distances, sampledf$expedition_number, mean) # get the means
monthigh_sd <- tapply(monthigh_beta$distances, sampledf$expedition_number, sd) # get the SD
monthigh_sd <- as.data.frame(monthigh_sd)
monthigh_sd
monthigh_mean <- as.data.frame(monthigh_mean)
monthigh_mean

colnames(monthigh_mean)[1] <- "mean"
colnames(monthigh_sd)[1] <- "sd"

monthigh_plot <- merge(monthigh_sd, monthigh_mean, by="row.names") 
monthigh_plot$Row.names
monthigh_plot$mean <-as.numeric(monthigh_plot$mean)
monthigh_plot
monthigh_plot$Row.names <- as.factor(monthigh_plot$Row.names)
monthigh_plot$sd <-as.numeric(monthigh_plot$sd)

# Montipora High Dispersion
p= ggplot(monthigh_plot, aes(x=Row.names, y=mean, col=Row.names, size=1)) + geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.1, size=1) + geom_point()
p=p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none",axis.text=element_text(size=12),axis.title=element_text(size=14)) + labs(x="Heat Stress") +labs(y="")
p= p + scale_colour_manual(values=c("darkblue", "darkred")) + scale_y_continuous(limits=c(0,0.9))+ scale_x_discrete(breaks=c("KI15b", "KI15c"), labels=c("Low", "High"))
p
##################################################


##################################################
# MONTIPORA LOW DIST.

otu_tablef_no10_coralf_sm_866f_rel_mont_vlow <- subset_samples(otu_tablef_no10_coralf_sm_866f_rel_mont, human_disturbance=="Very low")
print(otu_tablef_no10_coralf_sm_866f_rel_mont_vlow)
# 24 samples
any(taxa_sums(otu_tablef_no10_coralf_sm_866f_rel_mont_vlow) == 0)
#TRUE because i filtered
otu_tablef_no10_coralf_sm_866f_rel_mont_vlowf = prune_taxa(taxa_sums(otu_tablef_no10_coralf_sm_866f_rel_mont_vlow) > 0, otu_tablef_no10_coralf_sm_866f_rel_mont_vlow)
any(taxa_sums(otu_tablef_no10_coralf_sm_866f_rel_mont_vlowf) == 0)
#FALSE
print(otu_tablef_no10_coralf_sm_866f_rel_mont_vlowf) #24 

#save(otu_tablef_no10_coralf_sm_866f_rel_mont_vlowf, file="data/secondmito/lowreads_pruned/otutable_top10removed_coral_866_rel_montipora_vlow.Rdata")

# rename the levels
as(sample_data(otu_tablef_no10_coralf_sm_866f_rel_mont_vlowf), "data.frame")[as(sample_data(otu_tablef_no10_coralf_sm_866f_rel_mont_vlowf), "data.frame")$human_disturbance %in% "Very low", "human_disturbance2"] <- "Low"
as(sample_data(otu_tablef_no10_coralf_sm_866f_rel_mont_vlowf), "data.frame")$human_disturbance2
as(sample_data(otu_tablef_no10_coralf_sm_866f_rel_mont_vlowf), "data.frame")[as(sample_data(otu_tablef_no10_coralf_sm_866f_rel_mont_vlowf), "data.frame")$human_disturbance %in% "Very high", "human_disturbance2"] <- "High"
as(sample_data(otu_tablef_no10_coralf_sm_866f_rel_mont_vlowf), "data.frame")$human_disturbance2


# Montipora low disturbance
sample_data(otu_tablef_no10_coralf_sm_866f_rel_mont_vlowf)$human_disturbance2 # Low
print(otu_tablef_no10_coralf_sm_866f_rel_mont_vlowf) # 24 samples 

montlow_dist = phyloseq::distance(otu_tablef_no10_coralf_sm_866f_rel_mont_vlowf, "bray")


# Homogeneity of dispersion test 
sampledf <- data.frame(sample_data(otu_tablef_no10_coralf_sm_866f_rel_mont_vlowf))
montlow_beta <- betadisper(montlow_dist, sampledf$expedition_number)
plot(montlow_beta, hull=FALSE, ellipse=TRUE) #showing 1 standard deviation ellipses 
boxplot(montlow_beta)
permutest(montlow_beta, perm=9999)

# Number of permutations: 9999

# Response: Distances
#           Df   Sum Sq   Mean Sq      F N.Perm Pr(>F)  
# Groups     1 0.023931 0.0239306 4.4508   9999 0.0478 *
# Residuals 22 0.118288 0.0053767  


# Extract mean and SD (same as Anderson paper) so we don't have to use boxplots 
montlow_mean <- tapply(montlow_beta$distances, sampledf$expedition_number, mean) # get the means
montlow_sd <- tapply(montlow_beta$distances, sampledf$expedition_number, sd) # get the SD
montlow_sd <- as.data.frame(montlow_sd)
montlow_sd
montlow_mean <- as.data.frame(montlow_mean)
montlow_mean

colnames(montlow_mean)[1] <- "mean"
colnames(montlow_sd)[1] <- "sd"

montlow_plot <- merge(montlow_sd, montlow_mean, by="row.names") 
montlow_plot$Row.names
montlow_plot$mean <-as.numeric(montlow_plot$mean)
montlow_plot
montlow_plot$Row.names <- as.factor(montlow_plot$Row.names)
montlow_plot$sd <-as.numeric(montlow_plot$sd)

# Montipora Low Dispersion
p1= ggplot(montlow_plot, aes(x=Row.names, y=mean, col=Row.names, size=1)) + geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.1, size=1) + geom_point()
p1=p1 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none",axis.text=element_text(size=12),axis.title=element_text(size=14))
p1=p1 + labs(x="Heat Stress") +labs(y="Distance to centroid")
p1= p1 + scale_colour_manual(values=c("darkblue", "darkred"))+ scale_y_continuous(limits=c(0,0.9))+ scale_x_discrete(breaks=c("KI15b", "KI15c"), labels=c("Low", "High"))
p1
##################################################


##################################################
# PORITES HIGH DIST.

otu_tablef_no10_coralf_sm_866f_rel_por_vhigh <- subset_samples(otu_tablef_no10_coralf_sm_866f_rel_por, human_disturbance=="Very high")
print(otu_tablef_no10_coralf_sm_866f_rel_por_vhigh)
# 24 samples
any(taxa_sums(otu_tablef_no10_coralf_sm_866f_rel_por_vhigh) == 0)
#TRUE because i filtered
otu_tablef_no10_coralf_sm_866f_rel_por_vhighf = prune_taxa(taxa_sums(otu_tablef_no10_coralf_sm_866f_rel_por_vhigh) > 0, otu_tablef_no10_coralf_sm_866f_rel_por_vhigh)
any(taxa_sums(otu_tablef_no10_coralf_sm_866f_rel_por_vhighf) == 0)
#FALSE
print(otu_tablef_no10_coralf_sm_866f_rel_por_vhighf) #24 

#save(otu_tablef_no10_coralf_sm_866f_rel_por_vhighf, file="data/secondmito/lowreads_pruned/otutable_top10removed_coral_866_rel_porites_vhigh.Rdata")

# rename the levels
as(sample_data(otu_tablef_no10_coralf_sm_866f_rel_por_vhighf), "data.frame")[as(sample_data(otu_tablef_no10_coralf_sm_866f_rel_por_vhighf), "data.frame")$human_disturbance %in% "Very low", "human_disturbance2"] <- "Low"
as(sample_data(otu_tablef_no10_coralf_sm_866f_rel_por_vhighf), "data.frame")$human_disturbance2
as(sample_data(otu_tablef_no10_coralf_sm_866f_rel_por_vhighf), "data.frame")[as(sample_data(otu_tablef_no10_coralf_sm_866f_rel_por_vhighf), "data.frame")$human_disturbance %in% "Very high", "human_disturbance2"] <- "High"
as(sample_data(otu_tablef_no10_coralf_sm_866f_rel_por_vhighf), "data.frame")$human_disturbance2


sample_data(otu_tablef_no10_coralf_sm_866f_rel_por_vhighf)$human_disturbance2 # high
sample_data(otu_tablef_no10_coralf_sm_866f_rel_por_vhighf)$field_host_genus_id # Porites
print(otu_tablef_no10_coralf_sm_866f_rel_por_vhighf) # 24 samples 

# Beta-diversity test of porites very high
porhigh_dist = phyloseq::distance(otu_tablef_no10_coralf_sm_866f_rel_por_vhighf, "bray")

# Homogeneity of dispersion test 
sampledf <- data.frame(sample_data(otu_tablef_no10_coralf_sm_866f_rel_por_vhighf))
porhigh_dist_beta <- betadisper(porhigh_dist, sampledf$expedition_number)
plot(porhigh_dist_beta, hull=FALSE, ellipse=TRUE) #showing 1 standard deviation ellipses 
boxplot(porhigh_dist_beta)
permutest(porhigh_dist_beta, perm=9999)

# Number of permutations: 9999

# Response: Distances
#           Df  Sum Sq  Mean Sq      F N.Perm Pr(>F)
# Groups     1 0.03933 0.039332 0.7998   9999 0.3674
# Residuals 22 1.08185 0.049175   


# Extract mean and SD (same as Anderson paper) so we don't have to use boxplots 
porhigh_mean <- tapply(porhigh_dist_beta$distances, sampledf$expedition_number, mean) # get the means
porhigh_sd <- tapply(porhigh_dist_beta$distances, sampledf$expedition_number, sd) # get the SD
porhigh_sd <- as.data.frame(porhigh_sd)
porhigh_sd
porhigh_mean <- as.data.frame(porhigh_mean)
porhigh_mean

colnames(porhigh_mean)[1] <- "mean"
colnames(porhigh_sd)[1] <- "sd"

porhigh_plot <- merge(porhigh_sd, porhigh_mean, by="row.names") 
porhigh_plot$Row.names
porhigh_plot$mean <-as.numeric(porhigh_plot$mean)
porhigh_plot
porhigh_plot$Row.names <- as.factor(porhigh_plot$Row.names)
porhigh_plot$sd <-as.numeric(porhigh_plot$sd)

# Porites High Dispersion
p2= ggplot(porhigh_plot, aes(x=Row.names, y=mean, col=Row.names, size=1)) + geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.1, size=1) + geom_point(shape=17)
p2=p2 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none",axis.text.x=element_blank(),axis.text=element_text(size=12),axis.title=element_text(size=14)) + labs(x="") +labs(y="")
p2= p2 + scale_colour_manual(values=c("darkblue", "darkred"))+ scale_y_continuous(limits=c(0,0.9)) 
p2

##################################################


##################################################
# PORITES LOW DIST.

otu_tablef_no10_coralf_sm_866f_rel_por_vlow <- subset_samples(otu_tablef_no10_coralf_sm_866f_rel_por, human_disturbance=="Very low")
print(otu_tablef_no10_coralf_sm_866f_rel_por_vlow)
# 31 samples
any(taxa_sums(otu_tablef_no10_coralf_sm_866f_rel_por_vlow) == 0)
#TRUE because i filtered
otu_tablef_no10_coralf_sm_866f_rel_por_vlowf = prune_taxa(taxa_sums(otu_tablef_no10_coralf_sm_866f_rel_por_vlow) > 0, otu_tablef_no10_coralf_sm_866f_rel_por_vlow)
any(taxa_sums(otu_tablef_no10_coralf_sm_866f_rel_por_vlowf) == 0)
#FALSE
print(otu_tablef_no10_coralf_sm_866f_rel_por_vlowf) # 31

#save(otu_tablef_no10_coralf_sm_866f_rel_por_vlowf, file="data/secondmito/lowreads_pruned/otutable_top10removed_coral_866_rel_porites_vlow.Rdata")


# rename the levels
as(sample_data(otu_tablef_no10_coralf_sm_866f_rel_por_vlowf), "data.frame")[as(sample_data(otu_tablef_no10_coralf_sm_866f_rel_por_vlowf), "data.frame")$human_disturbance %in% "Very low", "human_disturbance2"] <- "Low"
as(sample_data(otu_tablef_no10_coralf_sm_866f_rel_por_vlowf), "data.frame")$human_disturbance2
as(sample_data(otu_tablef_no10_coralf_sm_866f_rel_por_vlowf), "data.frame")[as(sample_data(otu_tablef_no10_coralf_sm_866f_rel_por_vlowf), "data.frame")$human_disturbance %in% "Very high", "human_disturbance2"] <- "High"
as(sample_data(otu_tablef_no10_coralf_sm_866f_rel_por_vlowf), "data.frame")$human_disturbance2


sample_data(otu_tablef_no10_coralf_sm_866f_rel_por_vlowf)$human_disturbance2 # Low
print(otu_tablef_no10_coralf_sm_866f_rel_por_vlowf) # 31 samples 

porlow_dist = phyloseq::distance(otu_tablef_no10_coralf_sm_866f_rel_por_vlowf, "bray")

# Homogeneity of dispersion test 
sampledf <- data.frame(sample_data(otu_tablef_no10_coralf_sm_866f_rel_por_vlowf))
porlow_dist_beta <- betadisper(porlow_dist, sampledf$expedition_number)
plot(porlow_dist_beta, hull=FALSE, ellipse=TRUE) #showing 1 standard deviation ellipses 
boxplot(porlow_dist_beta)
permutest(porlow_dist_beta, perm=9999)
# Number of permutations: 9999

# Response: Distances
#           Df  Sum Sq Mean Sq      F N.Perm Pr(>F)    
# Groups     1 0.70563 0.70563 35.191   9999  1e-04 ***
# Residuals 29 0.58149 0.02005  


# Extract mean and SD (same as Anderson paper) so we don't have to use boxplots 
porlow_mean <- tapply(porlow_dist_beta$distances, sampledf$expedition_number, mean) # get the means
porlow_sd <- tapply(porlow_dist_beta$distances, sampledf$expedition_number, sd) # get the SD
porlow_sd <- as.data.frame(porlow_sd)
porlow_sd
porlow_mean <- as.data.frame(porlow_mean)
porlow_mean

colnames(porlow_mean)[1] <- "mean"
colnames(porlow_sd)[1] <- "sd"

porlow_plot <- merge(porlow_sd, porlow_mean, by="row.names") 
porlow_plot$Row.names
porlow_plot$mean <-as.numeric(porlow_plot$mean)
porlow_plot
porlow_plot$Row.names <- as.factor(porlow_plot$Row.names)
porlow_plot$sd <-as.numeric(porlow_plot$sd)

# Porites low Dispersion
p3= ggplot(porlow_plot, aes(x=Row.names, y=mean, col=Row.names,size=1)) + geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.1,size=1) + geom_point(shape=17)
p3=p3 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none",axis.text.x=element_blank(),axis.text=element_text(size=12),axis.title=element_text(size=14)) 
p3= p3 + labs(x="") +labs(y="Distance to centroid")
p3= p3 + scale_colour_manual(values=c("darkblue", "darkred"))+ scale_y_continuous(limits=c(0,0.9)) 
p3
##################################################


##################################################
# COMBINE PLOTS 
p=p+ annotate("text",x=.5, y=0.9, label="(d)",fontface=2, size=4.5)
p1=p1 + annotate("text", x=1.34, y=0.9, label="(c) Montipora aequituberculata",fontface=2, size=4.5) + annotate("text", x=2.5, y=0.84, label="*", size=6) 

p2=p2+ annotate("text", x=.5, y=0.9, label="(b)",fontface=2, size=4.5) + ggtitle("High disturbance")
p3=p3 + annotate("text", x=.95, y=0.9, label="(a) Porites lobata",fontface=2, size=4.5) + annotate("text", x=2.5, y=0.84, label="*", size=6)+ ggtitle("Low disturbance")

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

pdf(file="figures/secondmito/top10removed/distance/maytojuly_points_forms.pdf", height=8, width=8, onefile=FALSE)

grid.arrange(p3, p2,p1, p)

dev.off()





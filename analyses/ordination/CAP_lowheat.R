## Low Heat Stress (May) CAP

# relative abundance
# ***************************************************
# SECOND MITO
setwd("/Users/jamiemcdevitt-irwin/Documents/Git_Repos/McDevittIrwinetal_Kiritimati16S")
getwd()

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

library("grid")
packageVersion("grid")
theme_set(theme_bw())


# load normalized data 
load("data/secondmito/lowreads_pruned/otutable_top10removed_coral_866_may_rel.Rdata") # otu_tablef_no10_coralf_sm_866f_mayf_rel
ls()

sample_sums(otu_tablef_no10_coralf_sm_866f_mayf_rel) # 1

print(otu_tablef_no10_coralf_sm_866f_mayf_rel) # 39, 6037 taxa
##################################################


##################################################
# CAP
# backwards stepwise selection of variables 
otu_tablef_no10_coralf_sm_866f_mayf_rel_cap = ordinate(otu_tablef_no10_coralf_sm_866f_mayf_rel, formula=otu_tablef_no10_coralf_sm_866f_mayf_rel ~ field_host_genus_id + human_disturbance/reef_name, "CAP", "bray")
otu_tablef_no10_coralf_sm_866f_mayf_rel_cap # 27.55
otu_tablef_no10_coralf_sm_866f_mayf_rel_cap_nospecies = ordinate(otu_tablef_no10_coralf_sm_866f_mayf_rel, formula=otu_tablef_no10_coralf_sm_866f_mayf_rel ~ human_disturbance/reef_name, "CAP", "bray")
otu_tablef_no10_coralf_sm_866f_mayf_rel_cap_nospecies
otu_tablef_no10_coralf_sm_866f_mayf_rel_cap_nodist = ordinate(otu_tablef_no10_coralf_sm_866f_mayf_rel, formula=otu_tablef_no10_coralf_sm_866f_mayf_rel ~  field_host_genus_id, "CAP", "bray")
otu_tablef_no10_coralf_sm_866f_mayf_rel_cap_nodist


# Full model vs no coral species
anova(otu_tablef_no10_coralf_sm_866f_mayf_rel_cap_nospecies, otu_tablef_no10_coralf_sm_866f_mayf_rel_cap,permutations=9999)
# Model 1: OTU ~ human_disturbance/reef_name
# Model 2: OTU ~ field_host_genus_id + human_disturbance/reef_name
#   ResDf ResChiSquare Df ChiSquare      F Pr(>F)    
# 1    35       13.375                               
# 2    34       11.194  1    2.1809 6.6241  1e-04 ***
# # > 

# therefore coral species is important 

# Full model vs no human disturbance
anova(otu_tablef_no10_coralf_sm_866f_mayf_rel_cap_nodist, otu_tablef_no10_coralf_sm_866f_mayf_rel_cap, permutations=9999)
# Number of permutations: 9999

# Model 1: OTU ~ field_host_genus_id
# Model 2: OTU ~ field_host_genus_id + human_disturbance/reef_name
#   ResDf ResChiSquare Df ChiSquare      F Pr(>F)   
# 1    37       13.059                              
# 2    34       11.194  3    1.8648 1.8881 0.0046 **

# therefore disturbance is important

# therefore the best model includes: coral species and human disturbance/reef_name 
##################################################
#

##################################################
# PLOT
# rename the levels
as(sample_data(otu_tablef_no10_coralf_sm_866f_mayf_rel), "data.frame")[as(sample_data(otu_tablef_no10_coralf_sm_866f_mayf_rel), "data.frame")$human_disturbance %in% "Very low", "human_disturbance2"] <- "Low"
as(sample_data(otu_tablef_no10_coralf_sm_866f_mayf_rel), "data.frame")$human_disturbance2
as(sample_data(otu_tablef_no10_coralf_sm_866f_mayf_rel), "data.frame")[as(sample_data(otu_tablef_no10_coralf_sm_866f_mayf_rel), "data.frame")$human_disturbance %in% "Very high", "human_disturbance2"] <- "High"
as(sample_data(otu_tablef_no10_coralf_sm_866f_mayf_rel), "data.frame")$human_disturbance2


pdf(file="figures/secondmito/top10removed/ordination/CAPgenusdisturbance_reefnested_corals866_may_rel.pdf", height=6, width=6)


otu_tablef_no10_coralf_sm_866f_mayf_rel_cap = ordinate(otu_tablef_no10_coralf_sm_866f_mayf_rel, formula=otu_tablef_no10_coralf_sm_866f_mayf_rel ~ field_host_genus_id + human_disturbance2/reef_name, "CAP", "bray")
p5 = plot_ordination(otu_tablef_no10_coralf_sm_866f_mayf_rel, otu_tablef_no10_coralf_sm_866f_mayf_rel_cap, color="human_disturbance2", shape="field_host_genus_id") + scale_color_manual(values=c("darkorange2", "darkgoldenrod1"))
p5 = p5 + geom_point(size=3) + stat_ellipse()+ theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + labs(color='Local Disturbance') + labs(shape='Coral Species')
p5

otu_tablef_no10_coralf_sm_866f_mayf_rel_cap #27.55

dev.off()

save(p5, file="data/secondmito/lowreads_pruned/maycoral_CAP.RData")




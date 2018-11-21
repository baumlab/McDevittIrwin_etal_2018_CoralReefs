## PERMANOVA and Beta diversity tests/figures
# July: high heat stress
# ***************************************************

setwd("/Users/jamiemcdevitt-irwin/Documents/Git_Repos/McDevittIrwinetal_Kiritimati16S")
getwd()

# clear my environment
rm(list=ls()) 

# Installing & Importing  ----------------------------------------------------------
#source("http://bioconductor.org/biocLite.R")
#biocLite("phyloseq")

library("phyloseq")
packageVersion("phyloseq")
library("ggplot2")
packageVersion("ggplot2")
library("scales")
packageVersion("scales")
library("grid")
packageVersion("grid")
theme_set(theme_bw())
library("vegan")


load("data/secondmito/lowreads_pruned/otutable_top10removed_coral_866_july_rel.Rdata") # otu_tablef_no10_coralf_sm_866f_julyf_rel
ls()

print(otu_tablef_no10_coralf_sm_866f_julyf_rel) #64 samples, 7146 taxa 
############################################################



############################################################
# PERMANOVA

# calculate distance
julycoral_dist_track = phyloseq::distance(otu_tablef_no10_coralf_sm_866f_julyf_rel, "bray")

# permanova
adonis(julycoral_dist_track ~ field_host_genus_id + human_disturbance/reef_name, as(sample_data(otu_tablef_no10_coralf_sm_866f_julyf_rel), "data.frame"), perm=9999)

#                             Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# field_host_genus_id          1    2.5826 2.58264  6.9415 0.09591 0.0001 ***
# human_disturbance            1    0.8898 0.88976  2.3915 0.03304 0.0019 ** 
# human_disturbance:reef_name  2    1.5035 0.75173  2.0205 0.05583 0.0006 ***
# Residuals                   59   21.9513 0.37206         0.81521           
# Total                       63   26.9272                 1.00000           
############################################################



############################################################
# BETA DIVERSITY     
############################################################



############################################################
# MONTIPORA

otu_tablef_no10_coralf_sm_866f_julyf_rel_mont <- subset_samples(otu_tablef_no10_coralf_sm_866f_julyf_rel, field_host_genus_id=="Montipora")
print(otu_tablef_no10_coralf_sm_866f_julyf_rel_mont) #31 samples 

any(taxa_sums(otu_tablef_no10_coralf_sm_866f_julyf_rel_mont) == 0)
#TRUE because i filtered
otu_tablef_no10_coralf_sm_866f_julyf_rel_montf = prune_taxa(taxa_sums(otu_tablef_no10_coralf_sm_866f_julyf_rel_mont) > 0, otu_tablef_no10_coralf_sm_866f_julyf_rel_mont)
any(taxa_sums(otu_tablef_no10_coralf_sm_866f_julyf_rel_montf) == 0)
#FALSE
print(otu_tablef_no10_coralf_sm_866f_julyf_rel_montf) # 31 samples 

# save(otu_tablef_no10_coralf_sm_866f_julyf_rel_montf, file="data/secondmito/lowreads_pruned/otutable_top10removed_coral_866_rel_july_montipora.Rdata")


# rename the levels
as(sample_data(otu_tablef_no10_coralf_sm_866f_julyf_rel_montf), "data.frame")[as(sample_data(otu_tablef_no10_coralf_sm_866f_julyf_rel_montf), "data.frame")$human_disturbance %in% "Very low", "human_disturbance2"] <- "Low"
as(sample_data(otu_tablef_no10_coralf_sm_866f_julyf_rel_montf), "data.frame")$human_disturbance2
as(sample_data(otu_tablef_no10_coralf_sm_866f_julyf_rel_montf), "data.frame")[as(sample_data(otu_tablef_no10_coralf_sm_866f_julyf_rel_montf), "data.frame")$human_disturbance %in% "Very high", "human_disturbance2"] <- "High"
as(sample_data(otu_tablef_no10_coralf_sm_866f_julyf_rel_montf), "data.frame")$human_disturbance2


# calculate distance
july_mont_dist = phyloseq::distance(otu_tablef_no10_coralf_sm_866f_julyf_rel_montf, "bray")

# Homogeneity of dispersion test 
# Human Disturbance
sampledf <- data.frame(sample_data(otu_tablef_no10_coralf_sm_866f_julyf_rel_montf))
julymont_beta <- betadisper(july_mont_dist, sampledf$human_disturbance2)
plot(julymont_beta , hull=FALSE, ellipse=TRUE) #showing 1 standard deviation ellipses 
boxplot(julymont_beta )
permutest(julymont_beta , perm=9999)
# Permutation test for homogeneity of multivariate dispersions
# Permutation: free
# Number of permutations: 9999

# Response: Distances
#           Df   Sum Sq   Mean Sq      F N.Perm Pr(>F)
# Groups     1 0.003170 0.0031697 1.2228   9999 0.2774
# Residuals 29 0.075174 0.0025922  


save(julymont_beta, file="data/secondmito/lowreads_pruned/july_mont_boxplot.Rdata")
############################################################


############################################################
# PORITES

otu_tablef_no10_coralf_sm_866f_julyf_rel_por <- subset_samples(otu_tablef_no10_coralf_sm_866f_julyf_rel, field_host_genus_id=="Porites")
print(otu_tablef_no10_coralf_sm_866f_julyf_rel_por) #33 samples 

any(taxa_sums(otu_tablef_no10_coralf_sm_866f_julyf_rel_por) == 0)
#TRUE because i filtered
otu_tablef_no10_coralf_sm_866f_julyf_rel_porf = prune_taxa(taxa_sums(otu_tablef_no10_coralf_sm_866f_julyf_rel_por) > 0, otu_tablef_no10_coralf_sm_866f_julyf_rel_por)
any(taxa_sums(otu_tablef_no10_coralf_sm_866f_julyf_rel_porf) == 0)
#FALSE
print(otu_tablef_no10_coralf_sm_866f_julyf_rel_porf) # 33 samples 

# save(otu_tablef_no10_coralf_sm_866f_julyf_rel_porf, file="data/secondmito/lowreads_pruned/otutable_top10removed_coral_866_rel_july_porites.Rdata")

# rename the levels
as(sample_data(otu_tablef_no10_coralf_sm_866f_julyf_rel_porf), "data.frame")[as(sample_data(otu_tablef_no10_coralf_sm_866f_julyf_rel_porf), "data.frame")$human_disturbance %in% "Very low", "human_disturbance2"] <- "Low"
as(sample_data(otu_tablef_no10_coralf_sm_866f_julyf_rel_porf), "data.frame")$human_disturbance2
as(sample_data(otu_tablef_no10_coralf_sm_866f_julyf_rel_porf), "data.frame")[as(sample_data(otu_tablef_no10_coralf_sm_866f_julyf_rel_porf), "data.frame")$human_disturbance %in% "Very high", "human_disturbance2"] <- "High"
as(sample_data(otu_tablef_no10_coralf_sm_866f_julyf_rel_porf), "data.frame")$human_disturbance2



july_por_dist = phyloseq::distance(otu_tablef_no10_coralf_sm_866f_julyf_rel_porf, "bray")

# Homogeneity of dispersion test 
sampledf <- data.frame(sample_data(otu_tablef_no10_coralf_sm_866f_julyf_rel_porf))
julypor_beta <- betadisper(july_por_dist, sampledf$human_disturbance2)
plot(julypor_beta, hull=FALSE, ellipse=TRUE) #showing 1 standard deviation ellipses 
boxplot(julypor_beta)
permutest(julypor_beta, perm=9999)
# Number of permutations: 9999

# Response: Distances
#           Df  Sum Sq  Mean Sq     F N.Perm Pr(>F)
# Groups     1 0.03027 0.030275 0.818   9999 0.3717
# Residuals 31 1.14735 0.037011       


save(julypor_beta, file="data/secondmito/lowreads_pruned/july_por_boxplot.Rdata")
############################################################






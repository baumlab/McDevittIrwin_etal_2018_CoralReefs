
# PERMANOVA and Beta Diversity tests/figures
# May: low heat stress
# ***************************************************

setwd("/Users/jamiemcdevitt-irwin/Documents/Git_Repos/McDevittIrwinetal_Kiritimati16S")


# clear my environment
rm(list=ls()) 

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

load("data/secondmito/lowreads_pruned/otutable_top10removed_coral_866_may_rel.Rdata") # otu_tablef_no10_coralf_sm_866f_mayf_rel
ls()
sample_variables(otu_tablef_no10_coralf_sm_866f_mayf_rel)
############################################################



############################################################
# PERMANOVA
# calculate distance
maycoral_dist = phyloseq::distance(otu_tablef_no10_coralf_sm_866f_mayf_rel, "bray")


# permanova
adonis(maycoral_dist ~ field_host_genus_id + human_disturbance/reef_name, as(sample_data(otu_tablef_no10_coralf_sm_866f_mayf_rel), "data.frame"), perm=9999)
# Permutation: free
# Number of permutations: 9999

# Terms added sequentially (first to last)

#                             Df SumsOfSqs MeanSqs F.Model      R2 Pr(>F)    
# field_host_genus_id          1    2.3920 2.39205  7.2655 0.15482 0.0001 ***
# human_disturbance            1    1.2148 1.21485  3.6899 0.07863 0.0007 ***
# human_disturbance:reef_name  2    0.6500 0.32500  0.9871 0.04207 0.4241    
# Residuals                   34   11.1940 0.32924         0.72449           
# Total                       38   15.4509                 1.00000           
# ---
# Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
# > 
############################################################



############################################################
# BETA DIVERSITY     
############################################################


############################################################
# MONTIPORA

# subset out May Montipora 
otu_tablef_no10_coralf_sm_866f_mayf_rel_mont <- subset_samples(otu_tablef_no10_coralf_sm_866f_mayf_rel, field_host_genus_id=="Montipora")
print(otu_tablef_no10_coralf_sm_866f_mayf_rel_mont)
# 17 
any(taxa_sums(otu_tablef_no10_coralf_sm_866f_mayf_rel_mont) == 0)
#TRUE because i filtered
otu_tablef_no10_coralf_sm_866f_mayf_rel_montf = prune_taxa(taxa_sums(otu_tablef_no10_coralf_sm_866f_mayf_rel_mont) > 0, otu_tablef_no10_coralf_sm_866f_mayf_rel_mont)
any(taxa_sums(otu_tablef_no10_coralf_sm_866f_mayf_rel_montf) == 0)
#FALSE
print(otu_tablef_no10_coralf_sm_866f_mayf_rel_montf) #17 

# save(otu_tablef_no10_coralf_sm_866f_mayf_rel_montf, file="data/secondmito/lowreads_pruned/otutable_top10removed_coral_866_may_rel_mont.Rdata")



# rename the levels
as(sample_data(otu_tablef_no10_coralf_sm_866f_mayf_rel_montf), "data.frame")[as(sample_data(otu_tablef_no10_coralf_sm_866f_mayf_rel_montf), "data.frame")$human_disturbance %in% "Very low", "human_disturbance2"] <- "Low"
as(sample_data(otu_tablef_no10_coralf_sm_866f_mayf_rel_montf), "data.frame")$human_disturbance2
as(sample_data(otu_tablef_no10_coralf_sm_866f_mayf_rel_montf), "data.frame")[as(sample_data(otu_tablef_no10_coralf_sm_866f_mayf_rel_montf), "data.frame")$human_disturbance %in% "Very high", "human_disturbance2"] <- "High"
as(sample_data(otu_tablef_no10_coralf_sm_866f_mayf_rel_montf), "data.frame")$human_disturbance2

maymont_dist = phyloseq::distance(otu_tablef_no10_coralf_sm_866f_mayf_rel_montf, "bray")

# Homogeneity of dispersion test 
sampledf <- data.frame(sample_data(otu_tablef_no10_coralf_sm_866f_mayf_rel_montf))
maymont_beta <- betadisper(maymont_dist, sampledf$human_disturbance2)
plot(maymont_beta , hull=FALSE, ellipse=TRUE) #showing 1 standard deviation ellipses 
boxplot(maymont_beta )
permutest(maymont_beta, perm=9999)

# Permutation test for homogeneity of multivariate dispersions
# Permutation: free
# Number of permutations: 9999

# Response: Distances
#           Df   Sum Sq   Mean Sq      F N.Perm Pr(>F)
# Groups     1 0.008884 0.0088841 1.5446   9999 0.2461
# Residuals 15 0.086274 0.0057516     


save(maymont_beta, file="data/secondmito/lowreads_pruned/may_mont_boxplot.Rdata")
############################################################


############################################################
# PORITES

# subset 
otu_tablef_no10_coralf_sm_866f_mayf_rel_por <- subset_samples(otu_tablef_no10_coralf_sm_866f_mayf_rel, field_host_genus_id== "Porites")
print(otu_tablef_no10_coralf_sm_866f_mayf_rel_por)
# 22 
any(taxa_sums(otu_tablef_no10_coralf_sm_866f_mayf_rel_por) == 0)
#TRUE because i filtered
otu_tablef_no10_coralf_sm_866f_mayf_rel_porf = prune_taxa(taxa_sums(otu_tablef_no10_coralf_sm_866f_mayf_rel_por) > 0, otu_tablef_no10_coralf_sm_866f_mayf_rel_por)
any(taxa_sums(otu_tablef_no10_coralf_sm_866f_mayf_rel_porf) == 0)
#FALSE
print(otu_tablef_no10_coralf_sm_866f_mayf_rel_porf) # 22
as(sample_data(otu_tablef_no10_coralf_sm_866f_mayf_rel_porf), "data.frame")$reef_name


# save(otu_tablef_no10_coralf_sm_866f_mayf_rel_porf, file="data/secondmito/lowreads_pruned/otutable_top10removed_coral_866_may_rel_por.Rdata")


# rename the levels
as(sample_data(otu_tablef_no10_coralf_sm_866f_mayf_rel_porf), "data.frame")[as(sample_data(otu_tablef_no10_coralf_sm_866f_mayf_rel_porf), "data.frame")$human_disturbance %in% "Very low", "human_disturbance2"] <- "Low"
as(sample_data(otu_tablef_no10_coralf_sm_866f_mayf_rel_porf), "data.frame")$human_disturbance2
as(sample_data(otu_tablef_no10_coralf_sm_866f_mayf_rel_porf), "data.frame")[as(sample_data(otu_tablef_no10_coralf_sm_866f_mayf_rel_porf), "data.frame")$human_disturbance %in% "Very high", "human_disturbance2"] <- "High"
as(sample_data(otu_tablef_no10_coralf_sm_866f_mayf_rel_porf), "data.frame")$human_disturbance2


maypor_dist = phyloseq::distance(otu_tablef_no10_coralf_sm_866f_mayf_rel_porf, "bray")

# Homogeneity of dispersion test 
sampledf <- data.frame(sample_data(otu_tablef_no10_coralf_sm_866f_mayf_rel_porf))
maypor_beta <- betadisper(maypor_dist, sampledf$human_disturbance2)
boxplot(maypor_beta)
plot(maypor_beta, hull=FALSE, ellipse=TRUE) #showing 1 standard deviation ellipses 
permutest(maypor_beta, perm=9999)
# Number of permutations: 9999

# Response: Distances
#           Df  Sum Sq Mean Sq      F N.Perm Pr(>F)    
# Groups     1 0.61479 0.61479 24.719   9999  1e-04 ***
# Residuals 20 0.49743 0.02487    

save(maypor_beta, file="data/secondmito/lowreads_pruned/may_por_boxplot.Rdata")
############################################################


# May Alpha Diversity: Low heat stress
# ***************************************************

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
library("lme4")
library("MASS")
library("lsmeans")
library("MuMIn")

# load the data
load("data/secondmito/lowreads_pruned/otutable_top10removed_coral_866_may_rare.Rdata") #otu_tablef_no10_coralf_sm_866f_mayf_rare
ls()
########################################################



########################################################
# MODEL SELECTION

print(otu_tablef_no10_coralf_sm_866f_mayf_rare) # 39 samples, 2399 taxa 
sample_sums(otu_tablef_no10_coralf_sm_866f_mayf_rare) #1002


# calculate shannon diversity
may_rare_shannon <- data.frame(sample_data(otu_tablef_no10_coralf_sm_866f_mayf_rare), estimate_richness(otu_tablef_no10_coralf_sm_866f_mayf_rare), measures="Shannon" )
head(may_rare_shannon)

boxplot(Shannon ~ human_disturbance, may_rare_shannon)

qqnorm(may_rare_shannon$Shannon) # seems ok
class(may_rare_shannon$reef_name) # factor 


# Model Selection while keeping Site in as a nested effect 
# plot the contrasts of site but then do the tukey test on the human disturbance levels

# Pick the best model
model1<- lm(Shannon ~ 1, data=may_rare_shannon)
AICc(model1)

model5<- lm(Shannon ~ field_host_genus_id * human_disturbance/reef_name, data=may_rare_shannon)
summary(model5)# Adjusted R-squared:  0.479 
AICc(model5) # 113.0874

model6<- lm(Shannon ~ field_host_genus_id + human_disturbance/reef_name, data=may_rare_shannon) # same with human disturbance
AICc(model6) # 109.5598
summary(model6) # Adjusted R-squared:  0.4452 
plot(residuals(model6))
hist(residuals(model6))
plot(model6)

model7<- lm(Shannon ~ human_disturbance/reef_name, data=may_rare_shannon)
AICc(model7) # 119.0848

model8<- lm(Shannon ~ field_host_genus_id, data=may_rare_shannon)
AICc(model8) #118.4241

# so the best model is coral species + disturbance/site

model9 <- lm(Shannon ~ field_host_genus_id + human_disturbance/reef_name, data=may_rare_shannon)
AICc(model9) # 109.5598 (same when nested in HD)
summary(model9) #Adjusted R-squared:  0.4452 
hist(residuals(model9))
########################################################


########################################################
# PLOT

lsmeans.table <- function(x) {

 slm <- summary(x)
 slm_dat <- as.data.frame(slm[])
 return(slm_dat)

}


# combined human disturbance and coral species
model9.comb <- lsmeans(model9, c("human_disturbance", "field_host_genus_id"))
model9.comb

model9.comb_df <- lsmeans.table(model9.comb)
model9.comb_df



# add in human disturbance
model9.comb_df[model9.comb_df$human_disturbance %in% "Very low", "human_disturbance2"] <- "Low"
model9.comb_df
model9.comb_df[model9.comb_df$human_disturbance %in% "Very high", "human_disturbance2"] <- "High"
model9.comb_df

model9.comb_df$human_disturbance2 <- ordered(model9.comb_df$human_disturbance2, levels = c("Low","High"))
model9.comb_df

model9.comb_df$field_host_genus_id <- ordered(model9.comb_df$field_host_genus_id, levels = c("Porites","Montipora"))
model9.comb_df

# combined plot
pdf(file="figures/secondmito/top10removed/diversity/Shannon_corals_866rare_fittedmeans_may_combined.pdf", height=7, width=7)

p2= ggplot(model9.comb_df, aes(field_host_genus_id, lsmean, col=human_disturbance2, shape=field_host_genus_id)) 
p2= p2 + geom_pointrange(aes(ymin=lower.CL, ymax=upper.CL), size=1, position=position_dodge(width=0.5))
p2=p2 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text=element_text(size=12),axis.title=element_text(size=14)) 
p2= p2 + labs(x="") +labs(y="Shannon Fitted Means") + labs(colour="Local Disturbance", shape= "Coral Species")
p2= p2 + scale_colour_manual(values=c("darkgoldenrod1","darkorange2")) + scale_y_continuous(limits = c(1.3, 5))
p2= p2 + scale_shape_manual(values=c(17,19))

p2

dev.off()


# save for ms
may_p2=p2+ annotate("text", x=.5, y= 5.0, label="(a)", fontface=2, size=4.5)
may_p2
save(may_p2, file="data/secondmito/lowreads_pruned/may_coral_alphamodel_combined.Rdata")
########################################################



########################################################
# POST HOC CONTRASTS

model10 <- lm(Shannon ~ field_host_genus_id + human_disturbance, data=may_rare_shannon)
AICc(model10) # 104.6342
summary(model10) # Adjusted R-squared:  0.4689 
hist(residuals(model10))

# overall effect of human disturbance
model10.hd <- lsmeans(model10, "human_disturbance")
model10.hd
pairs(model10.hd)
# pairs(model10.hd)
#  contrast             estimate        SE df t.ratio p.value
#  Very high - Very low 1.192504 0.2761255 36   4.319  0.0001

# Results are averaged over the levels of: field_host_genus_id 
# > 

# overall effect of species
model10.gen<- lsmeans(model10, "field_host_genus_id")
model10.gen
pairs(model10.gen)
# > pairs(model10.gen)
#  contrast            estimate        SE df t.ratio p.value
#  Montipora - Porites 1.031526 0.2775982 36   3.716  0.0007

# Results are averaged over the levels of: human_disturbance 
# > 
########################################################



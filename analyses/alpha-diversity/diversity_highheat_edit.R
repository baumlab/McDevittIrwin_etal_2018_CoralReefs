# July Alpha Diversity: high heat stress
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
library("MuMIn")
library("grid")
packageVersion("grid")
theme_set(theme_bw())
#install.packages("lme4")
library("lme4")
library("MASS")
library("lsmeans")

# load the data 
load("data/secondmito/lowreads_pruned/otutable_top10removed_coral_866_july_rare.Rdata") #otu_tablef_no10_coralf_sm_866f_julyf_rare
ls()
########################################################


########################################################
# MODEL SELECTION

print(otu_tablef_no10_coralf_sm_866f_julyf_rare) # 64 samples, 3071 taxa
sample_sums(otu_tablef_no10_coralf_sm_866f_julyf_rare) # 867

# Calculate shannon diversity
july_rare_shannon <- data.frame(sample_data(otu_tablef_no10_coralf_sm_866f_julyf_rare), estimate_richness(otu_tablef_no10_coralf_sm_866f_julyf_rare), measures="Shannon" )
head(july_rare_shannon)


boxplot(Shannon ~ human_disturbance, july_rare_shannon)

qqnorm(july_rare_shannon$Shannon) # seems ok
class(july_rare_shannon$reef_name) # factor 


# Model Selection while keeping Site in as a nested effect 

# Pick the best model

model1<- lm(Shannon ~ 1, data=july_rare_shannon)
AICc(model1)

model5<- lm(Shannon ~ field_host_genus_id * human_disturbance/reef_name, data=july_rare_shannon)
summary(model5)# Adjusted R-squared:  0.1955 
AICc(model5) # 206.4857

model6<- lm(Shannon ~ field_host_genus_id + human_disturbance/reef_name, data=july_rare_shannon) # same with human disturbance
AICc(model6) #201.7589
summary(model6) # Adjusted R-squared:  0.1981 
plot(residuals(model6))
hist(residuals(model6))
plot(model6)

model7<- lm(Shannon ~ human_disturbance/reef_name, data=july_rare_shannon)
AICc(model7) # 215.2487

model8<- lm(Shannon ~ field_host_genus_id, data=july_rare_shannon)
AICc(model8) # 197.4097
summary(model8) #Adjusted R-squared:  0.2037 
plot(residuals(model8))
hist(residuals(model8))
plot(model8)

# so the best model is just coral species
########################################################


########################################################
# POST-HOC CONTRAST

model8_plot <- lsmeans(model8, c("field_host_genus_id"))
model8_plot
plot(model8_plot)
pairs(model8_plot)
# > pairs(model8_plot)
#  contrast            estimate        SE df t.ratio p.value
#  Montipora - Porites 1.131378 0.2734581 62   4.137  0.0001

########################################################

########################################################
# PLOT

# function to be able to plot the lsmeans results
lsmeans.table <- function(x) {
 slm <- summary(x)
 slm_dat <- as.data.frame(slm[])
 return(slm_dat)

}

# for plotting, include human disturbance so we can see the change
# you can't include /reefname in lsmeans
model9<- lm(Shannon ~ field_host_genus_id + human_disturbance, data=july_rare_shannon)
model9_plot <- lsmeans(model9, c("human_disturbance"))
model9_plot

# Note: this is not using the best model, this is for plotting so you can see the difference to may
model9.comb <- lsmeans(model9, c("human_disturbance", "field_host_genus_id"))
model9.comb

model9.comb_df <- lsmeans.table(model9.comb)
model9.comb_df


# Add in human disturbance
model9.comb_df[model9.comb_df$human_disturbance %in% "Very low", "human_disturbance2"] <- "Low"
model9.comb_df
model9.comb_df[model9.comb_df$human_disturbance %in% "Very high", "human_disturbance2"] <- "High"
model9.comb_df


# Reorder the levels
model9.comb_df$human_disturbance2 <- ordered(model9.comb_df$human_disturbance2, levels = c("Low","High"))
model9.comb_df

model9.comb_df$field_host_genus_id <- ordered(model9.comb_df$field_host_genus_id, levels = c("Porites","Montipora"))
model9.comb_df


# combined plot
pdf(file="figures/secondmito/top10removed/diversity/Shan_coralsf_866_rare_fittedmeans_july_combined.pdf", height=7, width=7)

p1= ggplot(model9.comb_df, aes(field_host_genus_id, lsmean, col=human_disturbance2, shape=field_host_genus_id)) 
p1= p1 + geom_pointrange(aes(ymin=lower.CL, ymax=upper.CL), size=1, position=position_dodge(width=0.5))
p1=p1 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text=element_text(size=12),axis.title.y=element_blank()) 
p1= p1 + labs(x="") +labs(y="Shannon Fitted Means") + labs(colour="Local Disturbance", shape= "Coral Species")
p1= p1 + scale_colour_manual(values=c("darkgoldenrod1","darkorange2")) + scale_y_continuous(limits = c(1.3, 5))
p1= p1+ scale_shape_manual(values=c(17,19))
p1

dev.off()

# save for ms
july_p1=p1+ annotate("text", x=.5, y=5.0, label="(b)",fontface=2, size=4.5)
july_p1
save(july_p1, file="data/secondmito/lowreads_pruned/july_coral_alphamodel_combined.Rdata")
########################################################


# Microbial Counts
############################

setwd("/Users/jamiemcdevitt-irwin/Documents/Git_Repos/McDevittIrwinetal_Kiritimati16S")
getwd()

# clear my environment
rm(list=ls()) 

library("grid")
packageVersion("grid")
theme_set(theme_bw())
library("lme4")
library("MASS")
library("ggplot2")
#install.packages("lsmeans")
library("lsmeans")
library("scales")
library("MuMIn")


cells_ml <- read.csv("data/envdata/KI_MicrobialCounts_avg_full.csv")
cells_ml
########################################################


########################################################
# MODEL SELECTION

hist(cells_ml$cells_per_ml)
min(cells_ml$cells_per_ml) # 73180.31
max(cells_ml$cells_per_ml) # 584485.8
qqnorm(cells_ml$cells_per_ml) 


cells_ml$Site <- as.factor(cells_ml$Site)
# Model Selection while keeping Site in as a nested effect 
# plot the contrasts of site but then do the tukey test on the human disturbance levels

# Pick the best model
model5<- lm(cells_per_ml ~ expedition_number * human_disturbance/Site, data=cells_ml)
summary(model5)# 0.9074 
AICc(model5) # 774.565

model6<- lm(cells_per_ml ~ expedition_number + human_disturbance/Site, data=cells_ml) # same with human disturbance
AICc(model6) # 767.2908
summary(model6) # Adjusted R-squared:  0.9074 
plot(residuals(model6))
hist(residuals(model6))
plot(model6)

model7<- lm(cells_per_ml ~ human_disturbance/Site, data=cells_ml)
AICc(model7) #768.63

model8<- lm(cells_per_ml ~ expedition_number, data=cells_ml)
AICc(model8) #836.3077

# so the best model is exp + dist/site

model9 <- lm(cells_per_ml ~ expedition_number + human_disturbance/Site, data=cells_ml)
AICc(model9) # 767.2908 
summary(model9) # Adjusted R-squared:  0.9074 
hist(residuals(model9))
########################################################

########################################################
# POSTHOC CONTRASTS

# can't use nested fixed effect for lsmeans
model11 <- lm(cells_per_ml ~ expedition_number + human_disturbance, data=cells_ml)
AICc(model11) # 780.3822 

# overall effect of human disturbance
model11.hd <- lsmeans(model11, "human_disturbance")
model11.hd
pairs(model11.hd)
# > pairs(model11.hd)
#  contrast             estimate       SE df t.ratio p.value
#  Very High - Very Low 288103.5 22974.63 28   12.54  <.0001

# Results are averaged over the levels of: expedition_number 
# > 

# overall effect of expedition
model11.ex<- lsmeans(model11, "expedition_number")
model11.ex
pairs(model11.ex)
# > pairs(model11.ex)
#  contrast      estimate       SE df t.ratio p.value
#  KI15b - KI15c 33039.63 22974.63 28   1.438  0.1615

# Results are averaged over the levels of: human_disturbance 
# > 
########################################################


########################################################
# PLOT

model10_plot <- lm(cells_per_ml ~ expedition_number + human_disturbance, data=cells_ml)

model10_dist <- lsmeans(model10_plot, "human_disturbance")
model10_dist

model10_exp <- lsmeans(model10_plot, "expedition_number")
model10_exp

# function to turn lsmeans into a dataframe
lsmeans.table <- function(x) {

 slm <- summary(x)
 slm_dat <- as.data.frame(slm[])
 return(slm_dat)

}

model10_exp_df <- lsmeans.table(model10_exp)
model10_exp_df

model10_dist_df <- lsmeans.table(model10_dist)
model10_dist_df


# change very to high and low
model10_dist_df[model10_dist_df$human_disturbance %in% "Very Low", "human_disturbance2"] <- "Low"
model10_dist_df
model10_dist_df[model10_dist_df$human_disturbance %in% "Very High", "human_disturbance2"] <- "High"
model10_dist_df

library("scales")


# Plot for human disturbance
pdf(file="figures/envdata/counts_disturbance.pdf", height=7, width=7)

p2= ggplot(model10_dist_df, aes(human_disturbance2, lsmean, col=human_disturbance2)) + geom_pointrange(aes(ymin=lower.CL, ymax=upper.CL),size=1, shape=8)
p2= p2 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.title.y=element_blank(),axis.text=element_text(size=12),axis.title=element_text(size=14))
p2= p2 + labs(x="Local Disturbance") +labs(y="Cells per ml") + labs(colour="Local Disturbance")
p2=p2 + scale_colour_manual(values=c("darkorange2","darkgoldenrod1"))+ scale_y_continuous(labels=scientific,limits = c(50000, 550000))
p2= p2 + scale_x_discrete(limits=c("Low", "High"))
p2
dev.off()

# rename the levels
model10_exp_df[model10_exp_df$expedition_number %in% "KI15b", "Hotspot"] <- "Low"
head(model10_exp_df)
model10_exp_df[model10_exp_df$expedition_number %in% "KI15c", "Hotspot"] <- "High"
head(model10_exp_df)

# change the order
model10_exp_df$Hotspot <- factor(model10_exp_df$Hotspot, levels= c("Low", "High"))


pdf(file="figures/envdata/counts_exp.pdf", height=7, width=7)

p1= ggplot(model10_exp_df, aes(Hotspot, lsmean, col=Hotspot)) + geom_pointrange(aes(ymin=lower.CL, ymax=upper.CL), size=1, shape=8)
p1=p1 + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + labs(x="Heat Stress") +labs(y="Cells per ml") + labs(colour="Human Disturbance") + 
scale_y_continuous(labels=scientific, limits = c(50000, 550000)) + 
  theme_bw() + theme(strip.background = element_blank(),legend.position="none", axis.text=element_text(size=12),axis.title=element_text(size=14),
    panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p1=p1 + scale_colour_manual(values=c("darkblue", "darkred")) 
p1

dev.off()



# combine plots: just human disturbance and expedition
p1=p1+ annotate("text", x=.5, y=550000, label="(a)", size=4.5)
p2=p2 + annotate("text", x=.5, y=550000, label="(b)", size=4.5) + annotate("text", x=2.4, y=525000, label="*", size=8)
p1
p2



# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)

  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)

  numPlots = length(plots)

  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }

 if (numPlots==1) {
    print(plots[[1]])

  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}



pdf(file="figures/envdata/cellsperml_exp_disturbancelevel.pdf", height=6, width=10)

layout <- matrix(c(1,1,1,2,2,2,2), nrow=1, byrow=TRUE)
multiplot(p1,p2, layout=layout)

dev.off()
########################################################



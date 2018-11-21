# Combine CAP plots for MS

setwd("/Users/jamiemcdevitt-irwin/Documents/Git_Repos/McDevittIrwinetal_Kiritimati16S")
getwd()
library("ggplot2")

# clear my environment
rm(list=ls()) 
theme_set(theme_bw())

##############################################################



##############################################################
# May
load("data/secondmito/lowreads_pruned/maycoral_CAP.RData") 
ls()

may_p5=p5 
may_p5=may_p5 + annotate("text", x=-1.9, y=2.8, label="(a)",fontface=2, size=4.5) +ggtitle("Low heat stress")


# July
load("data/secondmito/lowreads_pruned/july_coral_CAP.Rdata")
ls()

july_p1=p1
july_p1=july_p1 + annotate("text", x=-3.6, y=2.8, label="(b)",fontface=2, size=4.5) + ggtitle("High heat stress")


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

pdf(file="figures/secondmito/top10removed/ordination/CAP_combined_mayjuly.pdf", height=8, width=6, onefile=FALSE)

grid_arrange_shared_legend(may_p5, july_p1)

dev.off()


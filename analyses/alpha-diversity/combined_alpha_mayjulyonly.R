# Combine alpha plots MS

setwd("/Users/jamiemcdevitt-irwin/Documents/Git_Repos/McDevittIrwinetal_Kiritimati16S")
getwd()
library("ggplot2")

# clear my environment
rm(list=ls()) 
theme_set(theme_bw())

# load the data 
# May Alpha 
load("data/secondmito/lowreads_pruned/may_coral_alphamodel_combined.Rdata")

# July Alpha
load("data/secondmito/lowreads_pruned/july_coral_alphamodel_combined.Rdata")
ls()


may_p2=may_p2 + ggtitle("Low heat stress")

july_p1=july_p1 + ggtitle("High heat stress")


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


##put them all together
require(gtable)
legend = gtable_filter(ggplotGrob(may_p2), "guide-box") 
grid.draw(legend)

g_legend<-function(a.gplot){
tmp <- ggplot_gtable(ggplot_build(a.gplot))
leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
legend <- tmp$grobs[[leg]]
return(legend)}

legend <- g_legend(may_p2)
lwidth <- sum(legend$width)




pdf(file="figures/secondmito/top10removed/diversity/combined_alphadiv_for.ms_mayjulyonly_v2.pdf", height=8, width=11, onefile=FALSE)

grid.arrange(arrangeGrob(may_p2+ theme(legend.position="none") 
    + annotate("text", x=c(2.05,2.1), y=c(1.5, 1.3), label= c("Disturbance: p=0.0001", "Species: p=0.0007"), size=5), 
    july_p1 + labs(x="") + annotate("text", x=2.1, y= 1.3, label= "Species: p=0.0001", size=5)
    + theme(legend.position="none"), layout_matrix=rbind(c(1,2),c(1,2))), legend, widths=unit.c(unit(1, "npc") - lwidth, lwidth), nrow=1)

dev.off()



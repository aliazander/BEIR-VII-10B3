###############################################################################
#
# Run this script and get survival curves like [1]
#
# [1]: https://www.evernote.com/l/AAYupuVLtblNyZbvt3pwHv-U8OXiu0ufGd0B/image.png

library(plyr)
library(dplyr)
library(ggplot2)
library(scales)

setwd("/Users/g-woloschak/Downloads")
data <- read.csv('survival_graphs/data_for_survival.csv')
head(data)
# Retrieve the one unique value of x,
# raise an error if x has more than one value
only <- function(x) { 
      u <- unique(x)
      n <- length(u)
      if(n > 1) stop(paste("More than one value", paste(head(u), collapse=', ')))
      u
}

theme <- theme(
      text = element_text(size = 24, color="black"),
      axis.text = element_text(color="black"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(colour = "black", fill=NA),
      panel.background = element_blank(),
      legend.key = element_blank(),
      strip.background = element_blank()) 

# Calculate survival
data <- data %>%
      group_by(study, dose, group) %>%
      arrange(lifespan) %>%
      mutate(survival=rank(-lifespan) / length(lifespan))

# Add an observation at 100% survival for each group
h <- ddply(data, .(stratum, group), function(df) {
      top <- df[1,]
      top$survival = 1
      top$lifespan = only(df$censor)
      rbind(top, df)
})


# Reformat stratum names
h <- h %>%
      mutate(stratum = paste0(stratum_id, ". ", stratum),
             stratum = sub("???", "\n???", stratum),
             stratum = sub("???", "\n???", stratum),
             stratum = sub(" @.*", "", stratum)
      )


# Organize strata by number of animals
order <- h %>% 
      group_by(stratum) %>%
      summarize(n=length(stratum)) %>%
      arrange(-n) %>%
      select(stratum)
h$stratum <- factor(h$stratum, levels=order$stratum)

ggplot(h) +
      geom_path( data = subset(h, dose > 0),
            aes(lifespan,
                survival,
                colour= fractions >1,
                alpha = dose,
                #linetype=fractions > 1,
                group=group)) +
      geom_path( data = subset(h, dose == 0),
                 aes(lifespan,
                     survival,
                     colour= "grey30",
                     #linetype=fractions > 1,
                     group=group)) +
      facet_wrap(~ stratum_id) +
      scale_linetype(
            labels = c("acute", "protracted"),
            guide = guide_legend(title = "")) +
      scale_alpha_continuous("dose",
                           labels = c("1 Gy", "2 Gy", "3 Gy", "4 Gy"),
                           breaks = c(1, 2, 3, 4),
                           range = c(0.3, 1)) +
      scale_colour_manual(values = c("darkblue", "darkorange3", "grey30"),
                          labels = c("acute", "protracted", "control"),
                          guide = guide_legend(title = ""))+
      geom_vline(
            aes(xintercept=censor),
            color="grey",
            size=1) +
      scale_y_continuous(labels=percent) +
      expand_limits(x = -4) +
      theme(
            text = element_text(size = 10, color="black"),
            axis.text = element_text(color="black"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_rect(colour = "black", fill=NA),
            panel.background = element_blank(),
            legend.key = element_blank(),
            strip.background = element_blank()) +
      ylab("survival\n") +
      xlab("\nage (days)") 
setwd("/Users/g-woloschak/Documents")
ggsave("Bens-figs/all-graphs-survival.tiff",
       dpi = 600
)



############### individual survival graphs #########################################
#graph 1
#acute .3gy/min, total dose ~.3Gy 
#protracted 8 over 1 day, .3Gy/min
#10 over 10 days .3Gy/min

pointdata1 <- data.frame(lifespan = c(84, 84, 84, 84, 84), survival = c(.55, .51, .47, .43, .39))
pointdata2 <- data.frame(lifespan = seq(from = 84, to = 85, by = 1/7), survival = rep(.35, times = 8))
pointdata3 <- data.frame(lifespan = seq(from = 84, to = 85, by = 1/7), survival = rep(.31, times = 8))
pointdata4 <- data.frame(lifespan = seq(from = 84, to = 85, by = 1/7), survival = rep(.27, times = 8))
pointdata5 <- data.frame(lifespan = seq(from = 84, to = 93, by = 1), survival = rep(.23, times = 10))
pointdata6 <- data.frame(lifespan = seq(from = 84, to = 93, by = 1), survival = rep(.19, times = 10))
pointdata7 <- data.frame(lifespan = seq(from = 84, to = 93, by = 1), survival = rep(.15, times = 10))
pointdata8 <- data.frame(lifespan = seq(from = 84, to = 93, by = 1), survival = rep(.11, times = 10))
pointdata9 <- data.frame(lifespan = seq(from = 84, to = 93, by = 1), survival = rep(.07, times = 10))


h1 <- filter(h, stratum_id == 1)
plot1 <- ggplot(h1) +
      geom_path(data = subset(h1, dose > 0),
            aes(lifespan,
                survival,
                colour= fractions >1,
                alpha = dose,
                group=group),
            size = 1) +
      geom_path( data = subset(h1, dose == 0),
                 aes(lifespan,
                     survival,
                     colour= "grey30",
                     group=group),
                 size = 1) +
      facet_wrap(~ stratum) +
      scale_linetype(
            labels = c("acute", "protracted"),
            guide = guide_legend(title = "")) +
      scale_alpha_continuous("dose",
                             labels = c("0.25 Gy", "0.5 Gy","1.0 Gy","2.0 Gy","4.0 Gy"),
                             breaks = c(.25, .5, 1, 2, 4),
                             range = c(0.3, 1)) +
      scale_colour_manual(values = c("darkblue", "darkorange3", "grey30"),
                          labels = c("acute", "protracted", "control"),
                          guide = guide_legend(title = ""))+
      geom_vline(
            aes(xintercept=censor),
            color="grey",
            size=1) +
      scale_y_continuous(labels=percent) +
      expand_limits(x = -4) +
      theme(
            text = element_text(size = 12, color="black"),
            axis.text = element_text(color="black"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_rect(colour = "black", fill=NA),
            panel.background = element_blank(),
            legend.key = element_blank(),
            strip.background = element_blank()) +
      ylab("survival\n") +
      xlab("\nage (days)")
plot1
point1 <-      geom_point(data = pointdata1, 
                      mapping = aes(x = lifespan, y = survival),
                      pch = 24, size = 3, fill = "darkblue",
                      alpha = c(.3, .5, .6, .8, 1)) 
point2 <-      geom_point(data = pointdata2,
                 mapping = aes(x = lifespan, y = survival),
                 pch = 24, size = 3, fill = "darkorange3",
                 alpha = .6)
point3 <-      geom_point(data = pointdata3,
                  mapping = aes(x = lifespan, y = survival),
                  pch = 24, size = 3, fill = "darkorange3",
                  alpha = 0.8)
point4 <-      geom_point(data = pointdata4,
                  mapping = aes(x = lifespan, y = survival),
                  pch = 24, size = 3, fill = "darkorange3",
                  alpha = 1)
point5 <-      geom_point(data = pointdata5,
                          mapping = aes(x = lifespan, y = survival),
                          pch = 24, size = 3, fill = "darkorange3",
                          alpha = .3)
point6 <-      geom_point(data = pointdata6,
                          mapping = aes(x = lifespan, y = survival),
                          pch = 24, size = 3, fill = "darkorange3",
                          alpha = .5)
point7 <-      geom_point(data = pointdata7,
                          mapping = aes(x = lifespan, y = survival),
                          pch = 24, size = 3, fill = "darkorange3",
                          alpha = .6)
point8 <-      geom_point(data = pointdata8,
                          mapping = aes(x = lifespan, y = survival),
                          pch = 24, size = 3, fill = "darkorange3",
                          alpha = .8)
point9 <-      geom_point(data = pointdata9,
                          mapping = aes(x = lifespan, y = survival),
                          pch = 24, size = 3, fill = "darkorange3",
                          alpha = 1)
#overall suvival plot      
plot1 + point1 + point2 + point3 + point4 + point5 + 
      point6 + point7 + point8 + point9
ggsave("Bens-figs/graph1-survival.tiff",
       dpi = 600, height = 9, width = 11.5, units = "in")

#inset for survival plot
ggplot(h1)+ point1 + point2 + point3 + point4 + point5 + 
      point6 + point7 + point8 + point9+theme
ggsave("Bens-figs/graph1-survival-inset.tiff",
       dpi = 600, height = 9, width = 14.5, units = "in")


############### individual survival graphs #########################################
#graph 2
#24 over 168 days, .0037Gy/min - 4Gy
#24 over 168 days, .0018Gy/min - 2Gy
#120 continuous days - 2Gy

pointdata1 <- data.frame(lifespan = seq(from = 120, to = 287, by = 168/24), survival = rep(.35, times = 24))
pointdata2 <- data.frame(lifespan = seq(from = 120, to = 287, by = 168/24), survival = rep(.30, times = 24))
pointdata3 <- data.frame(lifespan = seq(from = 120, to = 239, by = 1), survival = rep(.25, times = 120))

h2 <- filter(h, stratum_id == 2)
plot2 <- ggplot(h2) +
      geom_path(data = subset(h2, dose > 0),
                aes(lifespan,
                    survival,
                    colour= fractions >1,
                    alpha = dose,
                    group=group),
                size = 1) +
      geom_path( data = subset(h2, dose == 0),
                 aes(lifespan,
                     survival,
                     colour= "grey30",
                     group=group),
                 size = 1) +
      facet_wrap(~ stratum) +
      scale_linetype(
            labels = c("acute", "protracted"),
            guide = guide_legend(title = "")) +
      scale_alpha_continuous("dose",
                             labels = c("2.0 Gy","4.0 Gy"),
                             breaks = c(2, 3.9),
                             range = c(0.6, 1)) +
      scale_colour_manual(values = c("darkorange3", "grey30"),
                          labels = c("protracted", "control"),
                          guide = guide_legend(title = ""))+
      geom_vline(
            aes(xintercept=censor),
            color="grey",
            size=1) +
      scale_y_continuous(labels=percent) +
      expand_limits(x = -4) +
      theme(
            text = element_text(size = 12, color="black"),
            axis.text = element_text(color="black"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_rect(colour = "black", fill=NA),
            panel.background = element_blank(),
            legend.key = element_blank(),
            strip.background = element_blank()) +
      ylab("survival\n") +
      xlab("\nage (days)")
plot2
point1 <-      geom_point(data = pointdata1, 
                          mapping = aes(x = lifespan, y = survival),
                          pch = 24, size = 3, fill = "darkorange3",
                          alpha = .6) 
point2 <-      geom_point(data = pointdata2,
                          mapping = aes(x = lifespan, y = survival),
                          pch = 24, size = 3, fill = "darkorange3",
                          alpha = 1)
point3 <-      geom_point(data = pointdata3,
                          mapping = aes(x = lifespan, y = survival),
                          pch = 24, size = 3, fill = "darkorange3",
                          alpha = .6)

#overall suvival plot      
plot2 + point1 + point2 + point3 
ggsave("Bens-figs/graph2-survival.tiff",
       dpi = 600, height = 9, width = 11.5, units = "in")

#inset for survival plot
ggplot(h2)+ point1 + point2 + point3 + theme
ggsave("Bens-figs/graph2-survival-inset.tiff",
       dpi = 600, height = 9, width = 14.5, units = "in")

############### individual survival graphs #########################################
#graph 3
#at 84 days
#acute - 4Gy/min - .25, .5, 1, 2, 4Gy
#10 over 10 days, 4Gy/min - .25, .5, 1, 2, 4Gy

pointdata1 <- data.frame(lifespan = c(84, 84, 84, 84, 84), survival = c(.39, .43, .47, .51, .55))
pointdata2 <- data.frame(lifespan = seq(from = 84, to = 93, by = 1), survival = rep(.35, times = 10))
pointdata3 <- data.frame(lifespan = seq(from = 84, to = 93, by = 1), survival = rep(.31, times = 10))
pointdata4 <- data.frame(lifespan = seq(from = 84, to = 93, by = 1), survival = rep(.27, times = 10))
pointdata5 <- data.frame(lifespan = seq(from = 84, to = 93, by = 1), survival = rep(.23, times = 10))
pointdata6 <- data.frame(lifespan = seq(from = 84, to = 93, by = 1), survival = rep(.19, times = 10))


h3 <- filter(h, stratum_id == 3)
plot3 <- ggplot(h3) +
      geom_path(data = subset(h3, dose > 0),
                aes(lifespan,
                    survival,
                    colour= fractions >1,
                    alpha = dose,
                    group=group),
                size = 1) +
      geom_path( data = subset(h3, dose == 0),
                 aes(lifespan,
                     survival,
                     colour= "grey30",
                     group=group),
                 size = 1) +
      facet_wrap(~ stratum) +
      scale_linetype(
            labels = c("acute", "protracted"),
            guide = guide_legend(title = "")) +
      scale_alpha_continuous("dose",
                             labels = c("0.25 Gy", "0.5 Gy","1.0 Gy","2.0 Gy","4.0 Gy"),
                             breaks = c(.25, .5, 1, 2, 4),
                             range = c(0.3, 1)) +
      scale_colour_manual(values = c("darkblue", "darkorange3", "grey30"),
                          labels = c("acute", "protracted", "control"),
                          guide = guide_legend(title = ""))+
      geom_vline(
            aes(xintercept=censor),
            color="grey",
            size=1) +
      scale_y_continuous(labels=percent) +
      expand_limits(x = -4) +
      theme(
            text = element_text(size = 12, color="black"),
            axis.text = element_text(color="black"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_rect(colour = "black", fill=NA),
            panel.background = element_blank(),
            legend.key = element_blank(),
            strip.background = element_blank()) +
      ylab("survival\n") +
      xlab("\nage (days)")
plot3
point1 <-      geom_point(data = pointdata1, 
                          mapping = aes(x = lifespan, y = survival),
                          pch = 24, size = 3, fill = "darkblue",
                          alpha = c(1, .8, .6, .5, .3)) 
point2 <-      geom_point(data = pointdata2,
                          mapping = aes(x = lifespan, y = survival),
                          pch = 24, size = 3, fill = "darkorange3",
                          alpha = .3)
point3 <-      geom_point(data = pointdata3,
                          mapping = aes(x = lifespan, y = survival),
                          pch = 24, size = 3, fill = "darkorange3",
                          alpha = .5)
point4 <-      geom_point(data = pointdata4,
                          mapping = aes(x = lifespan, y = survival),
                          pch = 24, size = 3, fill = "darkorange3",
                          alpha = .6)
point5 <-      geom_point(data = pointdata5,
                          mapping = aes(x = lifespan, y = survival),
                          pch = 24, size = 3, fill = "darkorange3",
                          alpha = .8)
point6 <-      geom_point(data = pointdata6,
                          mapping = aes(x = lifespan, y = survival),
                          pch = 24, size = 3, fill = "darkorange3",
                          alpha = 1)
#overall suvival plot      
plot3 + point1 + point2 + point3 + point4 + point5 + 
      point6
ggsave("Bens-figs/graph3-survival.tiff",
       dpi = 600, height = 9, width = 11.5, units = "in")

#inset for survival plot
ggplot(h3)+ point1 + point2 + point3 + point4 + point5 + 
      point6 + theme
ggsave("Bens-figs/graph3-survival-inset.tiff",
       dpi = 600, height = 9, width = 14.5, units = "in")

############### individual survival graphs #########################################
#graph 4
#acute - .133Gy/min - .3, .9, 1.5, 2.1Gy @-4 days
#acute - .133Gy/min - .5, 1, 2, 3, 4Gy @92 days
#acute - .133Gy/min - .5, 1, 2, 3, 4Gy @580 days

pointdata1 <- data.frame(lifespan = c(-4, -4, -4, -4), survival = c(.61, .57, .53, .49))
pointdata2 <- data.frame(lifespan = c(92, 92, 92, 92, 92), survival = c(.45, .41, .37, .33, .29))
pointdata3 <- data.frame(lifespan = c(580, 580, 580, 580, 580), survival = c(.25, .21, .17, .13, .09))

h4 <- filter(h, stratum_id == 4)
plot4 <- ggplot(h4) +
      geom_path(data = subset(h4, dose > 0),
                aes(lifespan,
                    survival,
                    colour= fractions >1,
                    alpha = dose,
                    group=group),
                size = 1) +
      geom_path( data = subset(h4, dose == 0),
                 aes(lifespan,
                     survival,
                     colour= "grey30",
                     group=group),
                 size = 1) +
      facet_wrap(~ stratum) +
      scale_linetype(
            labels = c("acute", "protracted"),
            guide = guide_legend(title = "")) +
      scale_alpha_continuous("dose",
                             labels = c("0.5 Gy","1.0 Gy","2.0 Gy","3.0 Gy","4.0 Gy"),
                             breaks = c(.5, 1, 2, 3, 4),
                             range = c(0.3, 1)) +
      scale_colour_manual(values = c("darkblue", "grey30"),
                          labels = c("acute", "control"),
                          guide = guide_legend(title = ""))+
      geom_vline(
            aes(xintercept=censor),
            color="grey",
            size=1) +
      scale_y_continuous(labels=percent) +
      expand_limits(x = -4) +
      theme(
            text = element_text(size = 12, color="black"),
            axis.text = element_text(color="black"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_rect(colour = "black", fill=NA),
            panel.background = element_blank(),
            legend.key = element_blank(),
            strip.background = element_blank()) +
      ylab("survival\n") +
      xlab("\nage (days)")
plot4
point1 <-      geom_point(data = pointdata1, 
                          mapping = aes(x = lifespan, y = survival),
                          pch = 24, size = 3, fill = "darkblue",
                          alpha = c(.3, .4, .5, .6)) 
point2 <-      geom_point(data = pointdata2,
                          mapping = aes(x = lifespan, y = survival),
                          pch = 24, size = 3, fill = "darkblue",
                          alpha = c(.3, .4, .6, .8, 1))
point3 <-      geom_point(data = pointdata3,
                          mapping = aes(x = lifespan, y = survival),
                          pch = 24, size = 3, fill = "darkblue",
                          alpha = c(.3, .4, .6, .8, 1))

#overall suvival plot      
plot4 + point1 + point2 + point3
ggsave("Bens-figs/graph4-survival.tiff",
       dpi = 600, height = 9, width = 11.5, units = "in")

#inset for survival plot
ggplot(h4)+ point1 + point2 + point3 + theme
ggsave("Bens-figs/graph4-survival-inset.tiff",
       dpi = 600, height = 9, width = 14.5, units = "in")

############### individual survival graphs #########################################
#graph 5
#acute - .099Gy/min - 1.9755Gy @520 days
#protracted - 60 over 420 days, .0015Gy/min 3.99Gy @100 days

pointdata1 <- data.frame(lifespan = c(520), survival = c(.3))
pointdata2 <- data.frame(lifespan = seq(from = 100, to = 519, by = 420/60), survival = rep(.25, times = 60))


h5 <- filter(h, stratum_id == 5)
plot5 <- ggplot(h5) +
      geom_path(data = subset(h5, dose > 0),
                aes(lifespan,
                    survival,
                    colour= fractions >1,
                    alpha = dose,
                    group=group),
                size = 1) +
      geom_path( data = subset(h5, dose == 0),
                 aes(lifespan,
                     survival,
                     colour= "grey30",
                     group=group),
                 size = 1) +
      facet_wrap(~ stratum) +
      scale_linetype(
            labels = c("acute", "protracted"),
            guide = guide_legend(title = "")) +
      scale_alpha_continuous("dose",
                             labels = c("2.0 Gy","4.0 Gy"),
                             breaks = c(2.1, 3.9),
                             range = c(0.6, 1)) +
      scale_colour_manual(values = c("darkblue", "darkorange3","grey30"),
                          labels = c("acute", "protracted", "control"),
                          guide = guide_legend(title = ""))+
      geom_vline(
            aes(xintercept=censor),
            color="grey",
            size=1) +
      scale_y_continuous(labels=percent) +
      expand_limits(x = -4) +
      theme(
            text = element_text(size = 12, color="black"),
            axis.text = element_text(color="black"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_rect(colour = "black", fill=NA),
            panel.background = element_blank(),
            legend.key = element_blank(),
            strip.background = element_blank()) +
      ylab("survival\n") +
      xlab("\nage (days)")
plot5
point1 <-      geom_point(data = pointdata1, 
                          mapping = aes(x = lifespan, y = survival),
                          pch = 24, size = 3, fill = "darkblue",
                          alpha = .6) 
point2 <-      geom_point(data = pointdata2,
                          mapping = aes(x = lifespan, y = survival),
                          pch = 24, size = 3, fill = "darkorange3",
                          alpha = 1)

#overall suvival plot      
plot5 + point1 + point2 
ggsave("Bens-figs/graph5-survival.tiff",
       dpi = 600, height = 9, width = 11.5, units = "in")

#inset for survival plot
ggplot(h5)+ point1 + point2 + theme
ggsave("Bens-figs/graph5-survival-inset.tiff",
       dpi = 600, height = 9, width = 14.5, units = "in")

############### individual survival graphs #########################################
#graph 6
#acute - 3.28 Gy @28 days
#protracted - 4 over 28 days, .94Gy/min -1.88, 2.82, 3.76Gy @28 days

pointdata1 <- data.frame(lifespan = c(28), survival = c(.35))
pointdata2 <- data.frame(lifespan = seq(from = 28, to = 57, by = 28/3), survival = rep(.30, times = 4))
pointdata3 <- data.frame(lifespan = seq(from = 28, to = 57, by = 28/3), survival = rep(.25, times = 4))
pointdata4 <- data.frame(lifespan = seq(from = 28, to = 57, by = 28/3), survival = rep(.20, times = 4))


h6 <- filter(h, stratum_id == 6)
plot6 <- ggplot(h6) +
      geom_path(data = subset(h6, dose > 0),
                aes(lifespan,
                    survival,
                    colour= fractions >1,
                    alpha = dose,
                    group=group),
                size = 1) +
      geom_path( data = subset(h6, dose == 0),
                 aes(lifespan,
                     survival,
                     colour= "grey30",
                     group=group),
                 size = 1) +
      facet_wrap(~ stratum) +
      scale_linetype(
            labels = c("acute", "protracted"),
            guide = guide_legend(title = "")) +
      scale_alpha_continuous("dose",
                             labels = c("2Gy", "3Gy", "4Gy"),
                             breaks = c(1.9, 2.9, 3.76),
                             range = c(0.7, 1)) +
      scale_colour_manual(values = c("darkblue", "darkorange3","grey30"),
                          labels = c("acute", "protracted", "control"),
                          guide = guide_legend(title = ""))+
      geom_vline(
            aes(xintercept=censor),
            color="grey",
            size=1) +
      scale_y_continuous(labels=percent) +
      expand_limits(x = -4) +
      theme(
            text = element_text(size = 12, color="black"),
            axis.text = element_text(color="black"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_rect(colour = "black", fill=NA),
            panel.background = element_blank(),
            legend.key = element_blank(),
            strip.background = element_blank()) +
      ylab("survival\n") +
      xlab("\nage (days)")
plot6
point1 <-      geom_point(data = pointdata1, 
                          mapping = aes(x = lifespan, y = survival),
                          pch = 24, size = 3, fill = "darkblue",
                          alpha = .85) 
point2 <-      geom_point(data = pointdata2,
                          mapping = aes(x = lifespan, y = survival),
                          pch = 24, size = 3, fill = "darkorange3",
                          alpha = .7)
point3 <-      geom_point(data = pointdata3,
                          mapping = aes(x = lifespan, y = survival),
                          pch = 24, size = 3, fill = "darkorange3",
                          alpha = .85)
point4 <-      geom_point(data = pointdata4,
                          mapping = aes(x = lifespan, y = survival),
                          pch = 24, size = 3, fill = "darkorange3",
                          alpha = 1)

#overall suvival plot      
plot6 + point1 + point2 + point3 + point4
ggsave("Bens-figs/graph6-survival.tiff",
       dpi = 600, height = 9, width = 11.5, units = "in")

#inset for survival plot
ggplot(h6)+ point1 + point2 + point3 + point4 + theme
ggsave("Bens-figs/graph6-survival-inset.tiff",
       dpi = 600, height = 9, width = 14.5, units = "in")

############### individual survival graphs #########################################
#graph 7
#acute - 1Gy/min - .5, 1, 3 Gy @ 7 days
#acute - 1Gy/min - .5, 1, 3 Gy @ 21 days


pointdata1 <- data.frame(lifespan = c(7, 7, 7), survival = c(.45, .41, .37))
pointdata2 <- data.frame(lifespan = c(21, 21, 21), survival = c(.33, .29, .25))

h7 <- filter(h, stratum_id == 7)
plot7 <- ggplot(h7) +
      geom_path(data = subset(h7, dose > 0),
                aes(lifespan,
                    survival,
                    colour= fractions >1,
                    alpha = dose,
                    group=group),
                size = 1) +
      geom_path(data = subset(h7, dose == 0),
                 aes(lifespan,
                     survival,
                     colour= "grey30",
                     group=group),
                 size = 1) +
      facet_wrap(~ stratum) +
      scale_linetype(
            labels = c("acute", "protracted"),
            guide = guide_legend(title = "")) +
      scale_alpha_continuous("dose",
                             labels = c("0.5 Gy","1.0 Gy","3.0 Gy"),
                             breaks = c(.5, 1, 3),
                             range = c(0.3, .9)) +
      scale_colour_manual(values = c("darkblue", "grey30"),
                          labels = c("acute", "control"),
                          guide = guide_legend(title = ""))+
      geom_vline(
            aes(xintercept=censor),
            color="grey",
            size=1) +
      scale_y_continuous(labels=percent) +
      expand_limits(x = -4) +
      theme(
            text = element_text(size = 12, color="black"),
            axis.text = element_text(color="black"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_rect(colour = "black", fill=NA),
            panel.background = element_blank(),
            legend.key = element_blank(),
            strip.background = element_blank()) +
      ylab("survival\n") +
      xlab("\nage (days)")
plot7
point1 <-      geom_point(data = pointdata1, 
                          mapping = aes(x = lifespan, y = survival),
                          pch = 24, size = 3, fill = "darkblue",
                          alpha = c(.3, .5, .9)) 
point2 <-      geom_point(data = pointdata2,
                          mapping = aes(x = lifespan, y = survival),
                          pch = 24, size = 3, fill = "darkblue",
                          alpha = c(.3, .5, .9))

#overall suvival plot      
plot7 + point1 + point2
ggsave("Bens-figs/graph7-survival.tiff",
       dpi = 600, height = 9, width = 11.5, units = "in")

#inset for survival plot
ggplot(h7)+ point1 + point2 + theme
ggsave("Bens-figs/graph7-survival-inset.tiff",
       dpi = 600, height = 9, width = 14.5, units = "in")

############### individual survival graphs #########################################
#graph 8
#at 56 days
#acute - .06 Gy/min - .025, .1, .4, 1, 4Gy
#5 over 35 days, .06Gy/min - .4 Gy

pointdata1 <- data.frame(lifespan = c(56, 56, 56, 56, 56), survival = c(.39, .43, .47, .51, .55))
pointdata2 <- data.frame(lifespan = seq(from = 56, to = 90, by = 35/5), survival = rep(.35, times = 10))

h8 <- filter(h, stratum_id == 8)
plot8 <- ggplot(h8) +
      geom_path(data = subset(h8, dose > 0),
                aes(lifespan,
                    survival,
                    colour= fractions >1,
                    alpha = dose,
                    group=group),
                size = 1) +
      geom_path( data = subset(h8, dose == 0),
                 aes(lifespan,
                     survival,
                     colour= "grey30",
                     group=group),
                 size = 1) +
      facet_wrap(~ stratum) +
      scale_linetype(
            labels = c("acute", "protracted"),
            guide = guide_legend(title = "")) +
      scale_alpha_continuous("dose",
                             labels = c("0.1 Gy", "0.4 Gy","1.0 Gy","4.0 Gy"),
                             breaks = c(.1, .4, 1, 4),
                             range = c(0.3, 1)) +
      scale_colour_manual(values = c("darkblue", "darkorange3", "grey30"),
                          labels = c("acute", "protracted", "control"),
                          guide = guide_legend(title = ""))+
      geom_vline(
            aes(xintercept=censor),
            color="grey",
            size=1) +
      scale_y_continuous(labels=percent) +
      expand_limits(x = -4) +
      theme(
            text = element_text(size = 12, color="black"),
            axis.text = element_text(color="black"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_rect(colour = "black", fill=NA),
            panel.background = element_blank(),
            legend.key = element_blank(),
            strip.background = element_blank()) +
      ylab("survival\n") +
      xlab("\nage (days)")
plot8
point1 <-      geom_point(data = pointdata1, 
                          mapping = aes(x = lifespan, y = survival),
                          pch = 24, size = 3, fill = "darkblue",
                          alpha = c(1, .8, .5, .4, .3)) 
point2 <-      geom_point(data = pointdata2,
                          mapping = aes(x = lifespan, y = survival),
                          pch = 24, size = 3, fill = "darkorange3",
                          alpha = .5)

plot8 + point1 + point2
ggsave("Bens-figs/graph8-survival.tiff",
       dpi = 600, height = 9, width = 11.5, units = "in")

#inset for survival plot
ggplot(h8)+ point1 + point2 + theme
ggsave("Bens-figs/graph8-survival-inset.tiff",
       dpi = 600, height = 9, width = 14.5, units = "in")

############### individual survival graphs #########################################
#graph 9
#at 56 days
#acute - .06 Gy/min - .08, .2, .4, 1.6 Gy
#20 over 28 days, .06 Gy/min - .4, 1.6 Gy

pointdata1 <- data.frame(lifespan = c(56, 56, 56, 56), survival = c(.39, .43, .47, .51))
pointdata2 <- data.frame(lifespan = seq(from = 56, to = 83, by = 28/20), survival = rep(.35, times = 20))
pointdata3 <- data.frame(lifespan = seq(from = 56, to = 83, by = 28/20), survival = rep(.31, times = 20))


h9 <- filter(h, stratum_id == 9)
plot9 <- ggplot(h9) +
      geom_path(data = subset(h9, dose > 0),
                aes(lifespan,
                    survival,
                    colour= fractions >1,
                    alpha = dose,
                    group=group),
                size = 1) +
      geom_path( data = subset(h9, dose == 0),
                 aes(lifespan,
                     survival,
                     colour= "grey30",
                     group=group),
                 size = 1) +
      facet_wrap(~ stratum) +
      scale_linetype(
            labels = c("acute", "protracted"),
            guide = guide_legend(title = "")) +
      scale_alpha_continuous("dose",
                             labels = c("0.08 Gy", "0.2 Gy","0.4 Gy","1.6 Gy"),
                             breaks = c(.08, .2, .4, 1.6),
                             range = c(0.3, .6)) +
      scale_colour_manual(values = c("darkblue", "darkorange3", "grey30"),
                          labels = c("acute", "protracted", "control"),
                          guide = guide_legend(title = ""))+
      geom_vline(
            aes(xintercept=censor),
            color="grey",
            size=1) +
      scale_y_continuous(labels=percent) +
      expand_limits(x = -4) +
      theme(
            text = element_text(size = 12, color="black"),
            axis.text = element_text(color="black"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_rect(colour = "black", fill=NA),
            panel.background = element_blank(),
            legend.key = element_blank(),
            strip.background = element_blank()) +
      ylab("survival\n") +
      xlab("\nage (days)")
plot9
point1 <-      geom_point(data = pointdata1, 
                          mapping = aes(x = lifespan, y = survival),
                          pch = 24, size = 3, fill = "darkblue",
                          alpha = c(.6, .5, .4, .3)) 
point2 <-      geom_point(data = pointdata2,
                          mapping = aes(x = lifespan, y = survival),
                          pch = 24, size = 3, fill = "darkorange3",
                          alpha = .5)
point3 <-      geom_point(data = pointdata3,
                          mapping = aes(x = lifespan, y = survival),
                          pch = 24, size = 3, fill = "darkorange3",
                          alpha = .6)

plot9 + point1 + point2 + point3
ggsave("Bens-figs/graph9-survival.tiff",
       dpi = 600, height = 9, width = 11.5, units = "in")

#inset for survival plot
ggplot(h9)+ point1 + point2 + point3 + theme
ggsave("Bens-figs/graph9-survival-inset.tiff",
       dpi = 600, height = 9, width = 14.5, units = "in")

############### individual survival graphs #########################################
#graph 10
#at 56 days
#acute - .06 Gy/min - .08, .4, 1.6 Gy
#20 over 28 days, .06 Gy/min - .4, 1.6 Gy

pointdata1 <- data.frame(lifespan = c(56, 56, 56), survival = c(.39, .43, .47))
pointdata2 <- data.frame(lifespan = seq(from = 56, to = 83, by = 28/20), survival = rep(.35, times = 20))
pointdata3 <- data.frame(lifespan = seq(from = 56, to = 83, by = 28/20), survival = rep(.31, times = 20))


h10 <- filter(h, stratum_id == 10)
plot10 <- ggplot(h10) +
      geom_path(data = subset(h10, dose > 0),
                aes(lifespan,
                    survival,
                    colour= fractions >1,
                    alpha = dose,
                    group=group),
                size = 1) +
      geom_path( data = subset(h10, dose == 0),
                 aes(lifespan,
                     survival,
                     colour= "grey30",
                     group=group),
                 size = 1) +
      facet_wrap(~ stratum) +
      scale_linetype(
            labels = c("acute", "protracted"),
            guide = guide_legend(title = "")) +
      scale_alpha_continuous("dose",
                             labels = c("0.08 Gy","0.4 Gy","1.6 Gy"),
                             breaks = c(.08, .4, 1.6),
                             range = c(0.3, .6)) +
      scale_colour_manual(values = c("darkblue", "darkorange3", "grey30"),
                          labels = c("acute", "protracted", "control"),
                          guide = guide_legend(title = ""))+
      geom_vline(
            aes(xintercept=censor),
            color="grey",
            size=1) +
      scale_y_continuous(labels=percent) +
      expand_limits(x = -4) +
      theme(
            text = element_text(size = 12, color="black"),
            axis.text = element_text(color="black"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_rect(colour = "black", fill=NA),
            panel.background = element_blank(),
            legend.key = element_blank(),
            strip.background = element_blank()) +
      ylab("survival\n") +
      xlab("\nage (days)")
plot10
point1 <-      geom_point(data = pointdata1, 
                          mapping = aes(x = lifespan, y = survival),
                          pch = 24, size = 3, fill = "darkblue",
                          alpha = c(.6, .5, .3)) 
point2 <-      geom_point(data = pointdata2,
                          mapping = aes(x = lifespan, y = survival),
                          pch = 24, size = 3, fill = "darkorange3",
                          alpha = .5)
point3 <-      geom_point(data = pointdata3,
                          mapping = aes(x = lifespan, y = survival),
                          pch = 24, size = 3, fill = "darkorange3",
                          alpha = .6)

plot10 + point1 + point2 + point3
ggsave("Bens-figs/graph10-survival.tiff",
       dpi = 600, height = 9, width = 11.5, units = "in")

#inset for survival plot
ggplot(h10)+ point1 + point2 + point3 + theme
ggsave("Bens-figs/graph10-survival-inset.tiff",
       dpi = 600, height = 9, width = 14.5, units = "in")

############### individual survival graphs #########################################
#graph 11
#at 56 days
#acute - ..9 Gy/min - 2 Gy
#10 over 305 days, .001 Gy/min - 2 Gy
#20 over 70 days, .001 Gy/min - 2 Gy
#10 over 14 days, .001 Gy/min - 2 Gy

pointdata1 <- data.frame(lifespan = c(56), survival = c(.39))
pointdata2 <- data.frame(lifespan = seq(from = 56, to = 360, by = 305/10), survival = rep(.35, times = 10))
pointdata3 <- data.frame(lifespan = seq(from = 56, to = 125, by = 70/20), survival = rep(.31, times = 20))
pointdata4 <- data.frame(lifespan = seq(from = 56, to = 69, by = 14/10), survival = rep(.27, times = 10))

h11 <- filter(h, stratum_id == 11)
plot11 <- ggplot(h11) +
      geom_path(data = subset(h11, dose > 0),
                aes(lifespan,
                    survival,
                    colour= fractions >1,
                    alpha = dose,
                    group=group),
                size = 1) +
      geom_path( data = subset(h11, dose == 0),
                 aes(lifespan,
                     survival,
                     colour= "grey30",
                     group=group),
                 size = 1) +
      facet_wrap(~ stratum) +
      scale_linetype(
            labels = c("acute", "protracted"),
            guide = guide_legend(title = "")) +
      scale_alpha_continuous("dose",
                             labels = c("2 Gy"),
                             breaks = c(2),
                             range = c(.8)) +
      scale_colour_manual(values = c("darkblue", "darkorange3", "grey30"),
                          labels = c("acute", "protracted", "control"),
                          guide = guide_legend(title = ""))+
      geom_vline(
            aes(xintercept=censor),
            color="grey",
            size=1) +
      scale_y_continuous(labels=percent) +
      expand_limits(x = -4) +
      theme(
            text = element_text(size = 12, color="black"),
            axis.text = element_text(color="black"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_rect(colour = "black", fill=NA),
            panel.background = element_blank(),
            legend.key = element_blank(),
            strip.background = element_blank()) +
      ylab("survival\n") +
      xlab("\nage (days)")
plot11
point1 <-      geom_point(data = pointdata1, 
                          mapping = aes(x = lifespan, y = survival),
                          pch = 24, size = 3, fill = "darkblue",
                          alpha = c(.8)) 
point2 <-      geom_point(data = pointdata2,
                          mapping = aes(x = lifespan, y = survival),
                          pch = 24, size = 3, fill = "darkorange3",
                          alpha = .8)
point3 <-      geom_point(data = pointdata3,
                          mapping = aes(x = lifespan, y = survival),
                          pch = 24, size = 3, fill = "darkorange3",
                          alpha = .8)
point4 <-      geom_point(data = pointdata4,
                          mapping = aes(x = lifespan, y = survival),
                          pch = 24, size = 3, fill = "darkorange3",
                          alpha = .8)

plot11 + point1 + point2 + point3 + point4
ggsave("Bens-figs/graph11-survival.tiff",
       dpi = 600, height = 9, width = 11.5, units = "in")

#inset for survival plot
ggplot(h11)+ point1 + point2 + point3 + point4 + theme
ggsave("Bens-figs/graph11-survival-inset.tiff",
       dpi = 600, height = 9, width = 14.5, units = "in")

############### individual survival graphs #########################################
#graph 12
#at 56 days
#acute - .06 Gy/min - 2 Gy
#10 over 305 days, .06 Gy/min - .2, 2 Gy
#20 over 70 days, .06 Gy/min - 2 Gy
#10 over 14 days, .06 Gy/min - 2 Gy

pointdata1 <- data.frame(lifespan = c(56), survival = c(.39))
pointdata2 <- data.frame(lifespan = seq(from = 56, to = 360, by = 305/10), survival = rep(.35, times = 10))
pointdata3 <- data.frame(lifespan = seq(from = 56, to = 360, by = 305/10), survival = rep(.31, times = 10))
pointdata4 <- data.frame(lifespan = seq(from = 56, to = 125, by = 70/20), survival = rep(.27, times = 20))
pointdata5 <- data.frame(lifespan = seq(from = 56, to = 69, by = 14/10), survival = rep(.23, times = 10))

h12 <- filter(h, stratum_id == 12)
plot12 <- ggplot(h12) +
      geom_path(data = subset(h12, dose > 0),
                aes(lifespan,
                    survival,
                    colour= fractions >1,
                    alpha = dose,
                    group=group),
                size = 1) +
      geom_path( data = subset(h12, dose == 0),
                 aes(lifespan,
                     survival,
                     colour= "grey30",
                     group=group),
                 size = 1) +
      facet_wrap(~ stratum) +
      scale_linetype(
            labels = c("acute", "protracted"),
            guide = guide_legend(title = "")) +
      scale_alpha_continuous("dose",
                             labels = c(".2 Gy","2 Gy"),
                             breaks = c(.2, 2),
                             range = c(.4, .8)) +
      scale_colour_manual(values = c("darkblue", "darkorange3", "grey30"),
                          labels = c("acute", "protracted", "control"),
                          guide = guide_legend(title = ""))+
      geom_vline(
            aes(xintercept=censor),
            color="grey",
            size=1) +
      scale_y_continuous(labels=percent) +
      expand_limits(x = -4) +
      theme(
            text = element_text(size = 12, color="black"),
            axis.text = element_text(color="black"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_rect(colour = "black", fill=NA),
            panel.background = element_blank(),
            legend.key = element_blank(),
            strip.background = element_blank()) +
      ylab("survival\n") +
      xlab("\nage (days)")
plot12
point1 <-      geom_point(data = pointdata1, 
                          mapping = aes(x = lifespan, y = survival),
                          pch = 24, size = 3, fill = "darkblue",
                          alpha = c(.8)) 
point2 <-      geom_point(data = pointdata2,
                          mapping = aes(x = lifespan, y = survival),
                          pch = 24, size = 3, fill = "darkorange3",
                          alpha = .4)
point3 <-      geom_point(data = pointdata3,
                          mapping = aes(x = lifespan, y = survival),
                          pch = 24, size = 3, fill = "darkorange3",
                          alpha = .8)
point4 <-      geom_point(data = pointdata4,
                          mapping = aes(x = lifespan, y = survival),
                          pch = 24, size = 3, fill = "darkorange3",
                          alpha = .8)
point5 <-      geom_point(data = pointdata5,
                          mapping = aes(x = lifespan, y = survival),
                          pch = 24, size = 3, fill = "darkorange3",
                          alpha = .8)

plot12 + point1 + point2 + point3 + point4 + point5
ggsave("Bens-figs/graph12-survival.tiff",
       dpi = 600, height = 9, width = 11.5, units = "in")

#inset for survival plot
ggplot(h12)+ point1 + point2 + point3 + point4 + point5 + theme
ggsave("Bens-figs/graph12-survival-inset.tiff",
       dpi = 600, height = 9, width = 14.5, units = "in")

############### individual survival graphs #########################################
#graph 13
#acute - .099Gy/min - 1.9755Gy @520 days
#protracted - 60 over 420 days, .0015Gy/min 3.99Gy @100 days

pointdata1 <- data.frame(lifespan = c(520), survival = c(.3))
pointdata2 <- data.frame(lifespan = seq(from = 100, to = 519, by = 420/60), survival = rep(.25, times = 60))


h13 <- filter(h, stratum_id == 13)
plot13 <- ggplot(h13) +
      geom_path(data = subset(h13, dose > 0),
                aes(lifespan,
                    survival,
                    colour= fractions >1,
                    alpha = dose,
                    group=group),
                size = 1) +
      geom_path( data = subset(h13, dose == 0),
                 aes(lifespan,
                     survival,
                     colour= "grey30",
                     group=group),
                 size = 1) +
      facet_wrap(~ stratum) +
      scale_linetype(
            labels = c("acute", "protracted"),
            guide = guide_legend(title = "")) +
      scale_alpha_continuous("dose",
                             labels = c("2.0 Gy","4.0 Gy"),
                             breaks = c(2.1, 3.9),
                             range = c(0.6, 1)) +
      scale_colour_manual(values = c("darkblue", "darkorange3","grey30"),
                          labels = c("acute", "protracted", "control"),
                          guide = guide_legend(title = ""))+
      geom_vline(
            aes(xintercept=censor),
            color="grey",
            size=1) +
      scale_y_continuous(labels=percent) +
      expand_limits(x = -4) +
      theme(
            text = element_text(size = 12, color="black"),
            axis.text = element_text(color="black"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_rect(colour = "black", fill=NA),
            panel.background = element_blank(),
            legend.key = element_blank(),
            strip.background = element_blank()) +
      ylab("survival\n") +
      xlab("\nage (days)")
plot13
point1 <-      geom_point(data = pointdata1, 
                          mapping = aes(x = lifespan, y = survival),
                          pch = 24, size = 3, fill = "darkblue",
                          alpha = .6) 
point2 <-      geom_point(data = pointdata2,
                          mapping = aes(x = lifespan, y = survival),
                          pch = 24, size = 3, fill = "darkorange3",
                          alpha = 1)

#overall suvival plot      
plot13 + point1 + point2 
ggsave("Bens-figs/graph13-survival.tiff",
       dpi = 600, height = 9, width = 11.5, units = "in")

#inset for survival plot
ggplot(h13)+ point1 + point2 + theme
ggsave("Bens-figs/graph13-survival-inset.tiff",
       dpi = 600, height = 9, width = 14.5, units = "in")

############### individual survival graphs #########################################
#graph 14
#acute - .06Gy/min - 2Gy @56 days
#protracted - 10 over 70 days, .06Gy/min - 2,4 Gy @56 days

pointdata1 <- data.frame(lifespan = c(56), survival = c(.3))
pointdata2 <- data.frame(lifespan = seq(from = 56, to = 125, by = 70/10), survival = rep(.25, times = 10))
pointdata3 <- data.frame(lifespan = seq(from = 56, to = 125, by = 70/10), survival = rep(.20, times = 10))

h14 <- filter(h, stratum_id == 14)
plot14 <- ggplot(h14) +
      geom_path(data = subset(h14, dose > 0),
                aes(lifespan,
                    survival,
                    colour= fractions >1,
                    alpha = dose,
                    group=group),
                size = 1) +
      geom_path( data = subset(h14, dose == 0),
                 aes(lifespan,
                     survival,
                     colour= "grey30",
                     group=group),
                 size = 1) +
      facet_wrap(~ stratum) +
      scale_linetype(
            labels = c("acute", "protracted"),
            guide = guide_legend(title = "")) +
      scale_alpha_continuous("dose",
                             labels = c("2.0 Gy","4.0 Gy"),
                             breaks = c(2, 4),
                             range = c(0.6, 1)) +
      scale_colour_manual(values = c("darkblue", "darkorange3","grey30"),
                          labels = c("acute", "protracted", "control"),
                          guide = guide_legend(title = ""))+
      geom_vline(
            aes(xintercept=censor),
            color="grey",
            size=1) +
      scale_y_continuous(labels=percent) +
      expand_limits(x = -4) +
      theme(
            text = element_text(size = 12, color="black"),
            axis.text = element_text(color="black"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_rect(colour = "black", fill=NA),
            panel.background = element_blank(),
            legend.key = element_blank(),
            strip.background = element_blank()) +
      ylab("survival\n") +
      xlab("\nage (days)")
plot14
point1 <-      geom_point(data = pointdata1, 
                          mapping = aes(x = lifespan, y = survival),
                          pch = 24, size = 3, fill = "darkblue",
                          alpha = .6) 
point2 <-      geom_point(data = pointdata2,
                          mapping = aes(x = lifespan, y = survival),
                          pch = 24, size = 3, fill = "darkorange3",
                          alpha = .6)
point3 <-      geom_point(data = pointdata3,
                          mapping = aes(x = lifespan, y = survival),
                          pch = 24, size = 3, fill = "darkorange3",
                          alpha = 1)

#overall suvival plot      
plot14 + point1 + point2 + point3
ggsave("Bens-figs/graph14-survival.tiff",
       dpi = 600, height = 9, width = 11.5, units = "in")

#inset for survival plot
ggplot(h14)+ point1 + point2 + point3 + theme
ggsave("Bens-figs/graph14-survival-inset.tiff",
       dpi = 600, height = 9, width = 14.5, units = "in")

############### individual survival graphs #########################################
#graph 15
#acute - .06Gy/min - 2Gy @56 days
#protracted - 10 over 70 days, .06Gy/min - 2,4 Gy @56 days

pointdata1 <- data.frame(lifespan = c(56), survival = c(.3))
pointdata2 <- data.frame(lifespan = seq(from = 56, to = 125, by = 70/10), survival = rep(.25, times = 10))
pointdata3 <- data.frame(lifespan = seq(from = 56, to = 125, by = 70/10), survival = rep(.20, times = 10))

h15 <- filter(h, stratum_id == 15)
plot15 <- ggplot(h15) +
      geom_path(data = subset(h15, dose > 0),
                aes(lifespan,
                    survival,
                    colour= fractions >1,
                    alpha = dose,
                    group=group),
                size = 1) +
      geom_path( data = subset(h15, dose == 0),
                 aes(lifespan,
                     survival,
                     colour= "grey30",
                     group=group),
                 size = 1) +
      facet_wrap(~ stratum) +
      scale_linetype(
            labels = c("acute", "protracted"),
            guide = guide_legend(title = "")) +
      scale_alpha_continuous("dose",
                             labels = c("2.0 Gy","4.0 Gy"),
                             breaks = c(2, 4),
                             range = c(0.6, 1)) +
      scale_colour_manual(values = c("darkblue", "darkorange3","grey30"),
                          labels = c("acute", "protracted", "control"),
                          guide = guide_legend(title = ""))+
      geom_vline(
            aes(xintercept=censor),
            color="grey",
            size=1) +
      scale_y_continuous(labels=percent) +
      expand_limits(x = -4) +
      theme(
            text = element_text(size = 12, color="black"),
            axis.text = element_text(color="black"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_rect(colour = "black", fill=NA),
            panel.background = element_blank(),
            legend.key = element_blank(),
            strip.background = element_blank()) +
      ylab("survival\n") +
      xlab("\nage (days)")
plot15
point1 <-      geom_point(data = pointdata1, 
                          mapping = aes(x = lifespan, y = survival),
                          pch = 24, size = 3, fill = "darkblue",
                          alpha = .6) 
point2 <-      geom_point(data = pointdata2,
                          mapping = aes(x = lifespan, y = survival),
                          pch = 24, size = 3, fill = "darkorange3",
                          alpha = .6)
point3 <-      geom_point(data = pointdata3,
                          mapping = aes(x = lifespan, y = survival),
                          pch = 24, size = 3, fill = "darkorange3",
                          alpha = 1)

#overall suvival plot      
plot15 + point1 + point2 + point3
ggsave("Bens-figs/graph15-survival.tiff",
       dpi = 600, height = 9, width = 11.5, units = "in")

#inset for survival plot
ggplot(h15)+ point1 + point2 + point3 + theme
ggsave("Bens-figs/graph15-survival-inset.tiff",
       dpi = 600, height = 9, width = 14.5, units = "in")

############### individual survival graphs #########################################
#graph 16
#acute .9gy/min, total dose ~.3Gy and 1.2Gy
#120 over 12 days .001Gy/min, total dose ~1.2Gy

pointdata1 <- data.frame(lifespan = c(117, 117), survival = c(.20, .15))
pointdata2 <- data.frame(lifespan = seq(from= 117, to = 128.9, by = .1), 
                         survival = rep(.10, times = 120))

h16 <- filter(h, stratum_id == 16)
plot16 <- ggplot(h16) +
      geom_path(data = subset(h16, dose > 0),
                aes(lifespan,
                    survival,
                    colour= fractions >1,
                    alpha = dose,
                    group=group),
                size = 1) +
      geom_path( data = subset(h16, dose == 0),
                 aes(lifespan,
                     survival,
                     colour= "grey30",
                     group=group),
                 size = 1) +
      facet_wrap(~ stratum) +
      scale_linetype(
            labels = c("acute", "protracted"),
            guide = guide_legend(title = "")) +
      scale_alpha_continuous("dose",
                             labels = c("0.3 Gy", "1.2 Gy"),
                             breaks = c(.3, 1.2),
                             range = c(0.3, .5)) +
      scale_colour_manual(values = c("darkblue", "darkorange3", "grey30"),
                          labels = c("acute", "protracted", "control"),
                          guide = guide_legend(title = ""))+
      geom_vline(
            aes(xintercept=censor),
            color="grey",
            size=1) +
      scale_y_continuous(labels=percent) +
      expand_limits(x = -4) +
      theme(
            text = element_text(size = 12, color="black"),
            axis.text = element_text(color="black"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_rect(colour = "black", fill=NA),
            panel.background = element_blank(),
            legend.key = element_blank(),
            strip.background = element_blank()) +
      ylab("survival\n") +
      xlab("\nage (days)") 

point1 <-      geom_point(data = pointdata1, 
                          mapping = aes(x = lifespan, y = survival),
                          pch = 24,
                          size = 3,
                          fill = "darkblue",
                          alpha = c(.3, .5))
point2 <-     geom_point(data = pointdata2, 
                         mapping = aes(x = lifespan, y = survival),
                         pch = 24,
                         size = 3,
                         fill = "darkorange3",
                         alpha = .5)
plot16 + point1 + point2
ggsave("Bens-figs/graph16-survival.tiff",
       dpi = 600, height = 9, width = 11.5, units = "in")

ggplot(h16)+ point1 + point2 + theme
ggsave("Bens-figs/graph16-survival-inset.tiff",
       dpi = 600, height = 9, width = 14.5, units = "in")

############### individual survival graphs #########################################
#graph 17
#acute .06gy/min, total dose 2Gy
#5 over 153 days .06Gy/min, total dose 4Gy

pointdata1 <- data.frame(lifespan = c(56), survival = c(.20))
pointdata2 <- data.frame(lifespan = seq(from= 56, to = 208, by = 153/5), 
                         survival = rep(.10, times = 5))

h17 <- filter(h, stratum_id == 17)
plot17 <- ggplot(h17) +
      geom_path(data = subset(h17, dose > 0),
                aes(lifespan,
                    survival,
                    colour= fractions >1,
                    alpha = dose,
                    group=group),
                size = 1) +
      geom_path( data = subset(h17, dose == 0),
                 aes(lifespan,
                     survival,
                     colour= "grey30",
                     group=group),
                 size = 1) +
      facet_wrap(~ stratum) +
      scale_linetype(
            labels = c("acute", "protracted"),
            guide = guide_legend(title = "")) +
      scale_alpha_continuous("dose",
                             labels = c("2 Gy", "4 Gy"),
                             breaks = c(2, 4),
                             range = c(0.8, 1)) +
      scale_colour_manual(values = c("darkblue", "darkorange3", "grey30"),
                          labels = c("acute", "protracted", "control"),
                          guide = guide_legend(title = ""))+
      geom_vline(
            aes(xintercept=censor),
            color="grey",
            size=1) +
      scale_y_continuous(labels=percent) +
      expand_limits(x = -4) +
      theme(
            text = element_text(size = 12, color="black"),
            axis.text = element_text(color="black"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_rect(colour = "black", fill=NA),
            panel.background = element_blank(),
            legend.key = element_blank(),
            strip.background = element_blank()) +
      ylab("survival\n") +
      xlab("\nage (days)") 
plot17
point1 <-      geom_point(data = pointdata1, 
                          mapping = aes(x = lifespan, y = survival),
                          pch = 24,
                          size = 3,
                          fill = "darkblue",
                          alpha = c(.8))
point2 <-     geom_point(data = pointdata2, 
                         mapping = aes(x = lifespan, y = survival),
                         pch = 24,
                         size = 3,
                         fill = "darkorange3",
                         alpha = 1)
plot17 + point1 + point2
ggsave("Bens-figs/graph17-survival.tiff",
       dpi = 600, height = 9, width = 11.5, units = "in")

ggplot(h17)+ point1 + point2 + theme
ggsave("Bens-figs/graph17-survival-inset.tiff",
       dpi = 600, height = 9, width = 14.5, units = "in")

############### individual survival graphs #########################################
#graph 18
#acute .9gy/min, total dose ~.3Gy and 1.2Gy
#120 over 12 days .001Gy/min, total dose ~.3Gy

pointdata1 <- data.frame(lifespan = c(56, 56), survival = c(.20, .15))
pointdata2 <- data.frame(lifespan = seq(from= 56, to = 175.9, by = 12/120), 
                         survival = rep(.10, times = 120))

h18 <- filter(h, stratum_id == 18)
plot18 <- ggplot(h18) +
      geom_path(data = subset(h18, dose > 0),
                aes(lifespan,
                    survival,
                    colour= fractions >1,
                    alpha = dose,
                    group=group),
                size = 1) +
      geom_path( data = subset(h18, dose == 0),
                 aes(lifespan,
                     survival,
                     colour= "grey30",
                     group=group),
                 size = 1) +
      facet_wrap(~ stratum) +
      scale_linetype(
            labels = c("acute", "protracted"),
            guide = guide_legend(title = "")) +
      scale_alpha_continuous("dose",
                             labels = c("0.3 Gy", "1.2 Gy"),
                             breaks = c(.3, 1.2),
                             range = c(0.3, .5)) +
      scale_colour_manual(values = c("darkblue", "darkorange3", "grey30"),
                          labels = c("acute", "protracted", "control"),
                          guide = guide_legend(title = ""))+
      geom_vline(
            aes(xintercept=censor),
            color="grey",
            size=1) +
      scale_y_continuous(labels=percent) +
      expand_limits(x = -4) +
      theme(
            text = element_text(size = 12, color="black"),
            axis.text = element_text(color="black"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_rect(colour = "black", fill=NA),
            panel.background = element_blank(),
            legend.key = element_blank(),
            strip.background = element_blank()) +
      ylab("survival\n") +
      xlab("\nage (days)") 
plot18
point1 <-      geom_point(data = pointdata1, 
                          mapping = aes(x = lifespan, y = survival),
                          pch = 24,
                          size = 3,
                          fill = "darkblue",
                          alpha = c(.3, .5))
point2 <-     geom_point(data = pointdata2, 
                         mapping = aes(x = lifespan, y = survival),
                         pch = 24,
                         size = 3,
                         fill = "darkorange3",
                         alpha = .3)
plot18 + point1 + point2
ggsave("Bens-figs/graph18-survival.tiff",
       dpi = 600, height = 9, width = 11.5, units = "in")

ggplot(h18)+ point1 + point2 + theme
ggsave("Bens-figs/graph18-survival-inset.tiff",
       dpi = 600, height = 9, width = 14.5, units = "in")


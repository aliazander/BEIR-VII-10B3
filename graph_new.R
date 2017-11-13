#############################################################################
#
# Making the main paper graph
#
# Run this file and get this graph [1]
#
# [1]: https://www.evernote.com/l/AAZC_FQggbZO3orHcapglKqdLbpzNp-1o78B/image.png

library(dplyr)
library(ggplot2)
library(scales)


setwd("/Users/g-woloschak/Downloads")
aggregated <- read.csv('figure data for alia/aggregated.csv')

# Clean up
g <- aggregated %>% 
  mutate(cluster_number = ifelse(dose == 0, "C", as.character(cluster_number)),
         type = ifelse(dose == 0, 'control', ifelse(protracted, 'protracted', 'acute')),
         stratum_full_name = paste0(stratum_id, ". ", stratum_full_name))

# Organize strata by number of animals
order <- g %>% 
  group_by(stratum_full_name) %>%
  summarize(n=sum(n)) %>%
  arrange(-n) %>%
  select(stratum_full_name)
g$stratum_full_name <- factor(g$stratum_full_name, levels=order$stratum_full_name)



ggplot(g,
       aes(dose,
           mortality,
           #label=cluster_number,
           color=type)) +
  geom_line(data = g %>% 
            filter(dose == 0 | protracted),
            aes(y=prediction), 
            color="darkorange3",
            show.legend=FALSE,
            stat="smooth",
            method="lm", 
            formula="y ~ x", 
            se=FALSE,
            size=.5) +
  geom_line(data = g %>% 
              filter(dose == 0 | !protracted),
            aes(y=prediction), 
            color="darkblue",
            show.legend=FALSE,
            stat="smooth",
            method="lm", 
            formula="y ~ x", 
            se=FALSE,
            size=.5) +
  geom_errorbar(aes(ymin=mortality - sem, ymax=mortality + sem), 
                colour="black", width=.1) +
  geom_point(size=1,
             #show.legend=FALSE,
             alpha=0.9,
             fill="white",
             pch=21) +
  #facet_wrap(~ stratum_full_name) +
  facet_wrap(~ stratum_id) +
  scale_y_continuous(labels=percent) +
  scale_size_continuous(range=c(3, 7)) +
  coord_cartesian(ylim=c(min(g$mortality) * 0.9 , max(g$mortality) * 1.1)) +
  xlab("\ndose (Gy)") +
  ylab("relative mortality\n") +
  scale_color_manual(values=c("darkblue", "grey30", "darkorange3"), 
                     name="") +
  theme(
    strip.text.x = element_text(hjust=0, size = 9),
    axis.text = element_text(color="black", size = 9),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(fill=NA),
    panel.background = element_blank(),
    legend.key = element_blank(),
    strip.background = element_blank()
  )
setwd("/Users/g-woloschak/Documents")
ggsave("Bens-figs/all-graphs-regression.tiff",
      dpi = 600
)




################## individual graphs ############################


graph.regression <- function(graph.num){
      graph <- subset(g, stratum_id == graph.num)
      ggplot(graph,
             aes(dose,
                 mortality,
                 label=cluster_number,
                 color=type)) +
            geom_line(data = graph %>% 
                            filter(dose == 0 | protracted),
                      aes(y=prediction), 
                      color="darkorange3",
                      show.legend=FALSE,
                      stat="smooth",
                      method="lm", 
                      formula="y ~ x", 
                      se=FALSE,
                      size=.5) +
            geom_line(data = graph %>% 
                            filter(dose == 0 | !protracted),
                      aes(y=prediction), 
                      color="darkblue",
                      show.legend=FALSE,
                      stat="smooth",
                      method="lm", 
                      formula="y ~ x", 
                      se=FALSE,
                      size=.5) +
            geom_errorbar(aes(ymin=mortality - sem, ymax=mortality + sem), 
                          colour="black", width=.1) +
            geom_point(size=5,
                       #show.legend=FALSE,
                       fill="white",
                       pch=21,
                       stroke = 1) +
            geom_text(size=4,
                      show.legend=FALSE) +
            scale_y_continuous(labels=percent) +
            scale_size_continuous(range=c(3, 7)) +
            coord_cartesian(ylim=c(min(graph$mortality) * 0.9 , max(graph$mortality) * 1.1)) +
            xlab("\ndose (Gy)") +
            ylab("relative mortality\n") +
            ggtitle(graph$stratum_full_name)+
            scale_color_manual(values=c("darkblue", "grey30", "darkorange3"), 
                               name="") +
            theme(
                  strip.text.x = element_text(hjust=0, size = 12),
                  axis.text = element_text(color="black", size = 12),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.border = element_rect(fill=NA),
                  panel.background = element_blank(),
                  legend.key = element_blank(),
                  strip.background = element_blank()
            )
      title <- paste("Bens-figs/graph", graph.num, "-regression.tiff", sep = "")
      ggsave(title,
             dpi = 600, height = 9, width = 11, units = "in")

}

for(i in 1:18){
      graph.regression(i)
}

#over write graph 2, colors messed up when only protracted
g2 <- subset(g, stratum_id == 2)
ggplot(g2,
       aes(dose,
           mortality,
           label=cluster_number,
           color=type)) +
      geom_line(data = g2 %>% 
                      filter(dose == 0 | protracted),
                aes(y=prediction), 
                color="darkorange3",
                show.legend=FALSE,
                stat="smooth",
                method="lm", 
                formula="y ~ x", 
                se=FALSE,
                size=.5) +
      geom_line(data = g2 %>% 
                      filter(dose == 0 | !protracted),
                aes(y=prediction), 
                color="darkblue",
                show.legend=FALSE,
                stat="smooth",
                method="lm", 
                formula="y ~ x", 
                se=FALSE,
                size=.5) +
      geom_errorbar(aes(ymin=mortality - sem, ymax=mortality + sem), 
                    colour="black", width=.1) +
      geom_point(size=5,
                 #show.legend=FALSE,
                 fill="white",
                 pch=21,
                 stroke = 1) +
      geom_text(size=4,
                show.legend=FALSE) +
      scale_y_continuous(labels=percent) +
      scale_size_continuous(range=c(3, 7)) +
      coord_cartesian(ylim=c(min(g2$mortality) * 0.9 , max(g2$mortality) * 1.1)) +
      xlab("\ndose (Gy)") +
      ylab("relative mortality\n") +
      ggtitle(g2$stratum_full_name)+
      scale_color_manual(values=c("grey30", "darkorange3"), 
                         name="") +
      theme(
            strip.text.x = element_text(hjust=0, size = 12),
            axis.text = element_text(color="black", size = 12),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_rect(fill=NA),
            panel.background = element_blank(),
            legend.key = element_blank(),
            strip.background = element_blank()
      )
ggsave("Bens-figs/graph2-regression.tiff",
       dpi = 600, height = 9, width = 11, units = "in")

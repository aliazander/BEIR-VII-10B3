library(plyr)
library(dplyr)
library(ggplot2)
library(stats4)
library(survival)
library(survminer)
library(ggpubr)
library(lattice)
library(car)
library(gridExtra)
library(ggfortify)
library(tidyr)
library(grid)
library(ggthemes)

setwd("/Users/g-woloschak/Documents")

#download data saved from website
#http://janus.northwestern.edu/janus2/data.php
macros <- read.csv('Lab_work/macro_pathologies.csv', header=TRUE, skip = 1)
macros[is.na(macros)] <- 0

demographics <- read.csv('Lab_work/demographics.csv', header=TRUE, skip = 1)

all.info <- cbind(demographics, macros)

#take lethal columns only
all.info <- select(all.info, animal_id:dose_rate, contains("_L"))

#delete JM11 (not real data set) and JM10 (peromyscus)
#delete neutron irradiated mice
#delete mice treated with radioprotectors
#delete breeder mice - not handled the same as others
#delete mice from JM 2 (different housing conditions)
#delete mice removed to another experiment - don't want to double count them
#delete mice that are from expt 3 - difference controls conditions
#delete mice that received 300 fractions, unknown stresses on mice interfering
all.info <- filter(all.info,
                   expt != 11,
                   expt != 10,
                   radn == "C" | radn == "G",
                   !(expt == 14 & !(tmt == "0S" | tmt == "C0")),
                   cause_of_death != "Grahn mice, breeder",
                   expt != 2,
                   expt != 12,
                   expt !=3,
                   fractions != 300,
                   cause_of_death != "Removal to another experiment")

all.info$expt <- as.factor(all.info$expt)
all.info$fractions <- as.factor(all.info$fractions)

#make sure each mouse only has one lethal macro --> good
lethal.columns <- grep("_L", names(all.info))
all.info$lethal_sum <- rowSums(all.info[,lethal.columns])

#add column for censoring: 0 = alive, 1 = dead
#censoring if wrong cause of death or no lethal disease listed
all.info$status <- 0
all.info$status[(all.info$cause_of_death == "Died" |
                       all.info$cause_of_death == "Sacrifice, moribund") &
                      all.info$lethal_sum == 1] <- 1

###################   Cox model for experiments #########################
#releveling to change reference level of factors to make sure
#significance doesn't change
#all.info$expt <- relevel(all.info$expt, ref = "14")

res.cox <- coxph(Surv(age, status) ~ sex + expt, 
                 data = subset(all.info, total_dose == 0))
summary(res.cox)
newd <- expand.grid(sex = c("M", "F"), expt = c("4", "7", "8", "9", "13", "14"))
newd$id <- 1:12

fortify(survfit(res.cox, newdata = newd)) %>%
      gather(strata, surv, surv.1:surv.12) %>%
      mutate(id = gsub("surv.","", strata)) %>%
      merge(newd, by = "id") %>%
      ggplot(aes(x = time, y = surv, col = expt)) +
      geom_step(size = .75) + 
      facet_grid(. ~ sex, labeller = as_labeller(c('M' = "Male", 'F' = "Female"))) +
      labs(x = "Time (days)", y = "Survival probability", colour = "Experiment") + 
      theme_bw()+
      theme(strip.background =element_rect(fill="white"))+
      theme(strip.text = element_text(colour = 'black', face = "bold", size = 12))+
      theme(text=element_text(size=12, face = "bold")) 

ggsave(file="Survival-expt-controls.png", height = 4, width = 8, units = "in")
# Test PH assumption - check residuals do not change with time
# values should not be significant if assuption is OK
test.ph <- cox.zph(res.cox)
print(test.ph)
plot2 <- ggcoxzph(test.ph,font.main = 11,
                  ggtheme = theme_bw()+
                        theme(text=element_text(size=12, face = "bold")))
plot2[[1]]
ggsave(file="Cox-proportional-expt-controls-1.png", height = 2.5, width = 4, units = "in")

plot3 <- ggcoxdiagnostics(res.cox, type = "deviance",
                          linear.predictions = FALSE, 
                          point.size = .5,
                          ggtheme = theme_bw() + 
                                theme(strip.background =element_rect(fill="white"))+
                                theme(strip.text = element_text(colour = 'black', face = "bold", size = 12))+
                                theme(text=element_text(size=12, face = "bold")))
ggsave(file="Residuals-expt-controls.png", height = 3, width = 4, units = "in")

######################### now for fractions ##########################################

all.info.no8 <- subset(all.info, expt !=8)
all.info.no8$fractions <- droplevels(all.info.no8$fractions)

#releveling to change reference level of factors to make sure
#significance doesn't change
#all.info.no8$fractions <- relevel(all.info.no8$fractions, ref = "120")

res.cox <- coxph(Surv(age, status) ~ sex + fractions, 
                 data = subset(all.info.no8, total_dose == 0))
summary(res.cox)
newd <- expand.grid(sex = c("M", "F"), fractions = levels(all.info.no8$fractions))
newd$id <- 1:10

fortify(survfit(res.cox, newdata = newd)) %>%
      gather(strata, surv, surv.1:surv.10) %>%
      mutate(id = gsub("surv.","", strata)) %>%
      merge(newd, by = "id") %>%
      ggplot(aes(x = time, y = surv, col = fractions)) +
      geom_step(size = .75) + 
      facet_grid(. ~ sex, labeller = as_labeller(c('M' = "Male", 'F' = "Female"))) +
      labs(x = "Time (days)", y = "Survival probability", colour = "Fractions") + 
      theme_bw()+
      theme(strip.background =element_rect(fill="white"))+
      theme(strip.text = element_text(colour = 'black', face = "bold", size = 12))+
      theme(text=element_text(size=12, face = "bold"))+
      labs(fill = "Experiment")
ggsave(file="Survival-fractions-controls.png", height = 4, width = 8, units = "in")
# Test PH assumption - check residuals do not change with time
# values should not be significant if assuption is OK
test.ph <- cox.zph(res.cox)
print(test.ph)
plot2 <- ggcoxzph(test.ph,font.main = 11,
                  ggtheme = theme_bw()+
                        theme(text=element_text(size=12, face = "bold")))
plot2[[1]]
ggsave(file="Cox-proportional-expt-controls-1.png", height = 2.5, width = 4, units = "in")

plot3 <- ggcoxdiagnostics(res.cox, type = "deviance",
                          linear.predictions = FALSE, 
                          point.size = .5,
                          ggtheme = theme_bw() + 
                                theme(strip.background =element_rect(fill="white"))+
                                theme(strip.text = element_text(colour = 'black', face = "bold", size = 12))+
                                theme(text=element_text(size=12, face = "bold")))
ggsave(file="Residuals-fractions-controls.png", height = 3, width = 4, units = "in")


write.csv(all.info, file = "Filtered_data.csv")
head(all.info)

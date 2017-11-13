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

##################       functions       ##############################
ordered <- function(unordered, ordered_by, logical){
      unordered[order(unordered[ordered_by], decreasing = logical),]
}

##################       end functions       ##############################

macros <- read.csv('Lab_work/macro_pathologies.csv', header=TRUE, skip = 1)
macros[is.na(macros)] <- 0

demographics <- read.csv('Lab_work/demographics.csv', header=TRUE, skip = 1)

all_info <- cbind(demographics, macros)


ncol(demographics)
ncol(macros)
#19 deomographics columns, 292 macros columns, first row of macros is animal ID

all_info$sum <- rowSums(all_info[,21:(20+290)])

mice_to_use <- filter(all_info,
                      #total_dose <= 200,
                      cause_of_death == "Died" | cause_of_death == "Sacrifice, moribund",
                      expt != 2,
                      #expt != 8,
                      expt != 10,
                      expt != 11,
                      expt != 12,
                      !(expt ==4 & (fractions == 120 | fractions == 300)),
                      !(expt == 4 & sex == "M"),
                      !(expt == 9 & sex == "M"),
                      #this should be censored instead
                      sum >= 1,
                      radn == "C" | radn == "G",
                      !(expt == 14 & !(tmt == "0S" | tmt == "C0")))
#check that the tricky filtering worked like I wanted it to
JM14 <- filter(mice_to_use,
               expt == 14)
JM14 %>%
      group_by(expt, sex, tmt) %>%
      summarise(n=n())

demographics %>%
      group_by(cause_of_death) %>%
      summarise(n=n(),
                avg.age <- mean(age),
                avg.dose <- mean(total_dose))
summary(subset(demographics, cause_of_death == "Accidental death"))
summary(demographics$age)
summary(is.na(demographics$age))
JM4 <- filter(mice_to_use,
               expt == 4)
JM4 %>%
      group_by(expt, sex, fractions) %>%
      summarise(n=n())
#to send to collaborators that want more data
JM8 <- filter(mice_to_use,
              expt ==8)

fit <- survfit(Surv(age) ~ dose_rate, data = JM8)
ggsurvplot(fit)
survdiff(Surv(age) ~ dose_rate, data = JM8)

write.csv(JM8, file = "JM8.csv",
          row.names=TRUE)

mice_to_use %>%
      group_by(dose_rate, expt) %>%
      summarise(n=n())

#10788 mice now
nrow(mice_to_use)

#make sure each mouse only has one lethal macro --> good
lethal_columns <- grep("_L", names(mice_to_use))
mice_to_use$lethal_sum <- rowSums(mice_to_use[,lethal_columns])
mice_to_use <- mice_to_use[mice_to_use$lethal_sum == 1,] 

#get lethal columns only
data_to_use <- select(mice_to_use, animal_id:dose_rate, contains("_L"))

#classify as tumor or non-tumor, this will work to use later
substr(colnames(data_to_use), 1, 1) =="T"

#create new data frame that has all types of lethal causes as rows
#first column has total number of mice with that cancer
#second calumn is if it's a tumor or not
#third column is whether it's lymphona or not.
#fourth column is total counts of tumors
#fifth column is total counts of non-tumors
#sixth column is total counts of lymphoma
age_quartiles <- summary(data_to_use$age)
age_quartiles<- age_quartiles[-4] #remove mean

quartile_1 <- filter(data_to_use, 
                     age <= age_quartiles[2])
quartile_2 <- filter(data_to_use, 
                     age > age_quartiles[2],
                     age <= age_quartiles[3])
quartile_3 <- filter(data_to_use, 
                     age > age_quartiles[3],
                     age <= age_quartiles[4])
quartile_4 <- filter(data_to_use, 
                     age > age_quartiles[4],
                     age <= age_quartiles[5])

death_counts_df <- function(quartile){
      death_counts <- data.frame(occurrences = colSums(quartile[,21:ncol(quartile)-1]))
      death_counts$tumor <- substr(rownames(death_counts), 1, 1) =="T"
      lymphoma <- c("NTYG_L","NTYL_L", "TTYG_L", "TTYL_L")
      death_counts$lymphoma <- rownames(death_counts) %in% lymphoma
      death_counts$total_tumors <- death_counts$occurrences * death_counts$tumor * !(death_counts$lymphoma)
      death_counts$total_non_tumors <- (death_counts$occurrences * !(death_counts$tumor)) * !(death_counts$lymphoma)
      death_counts$total_lymphoma <- death_counts$occurrences * death_counts$lymphoma
      death_counts <- ordered(death_counts, "occurrences", TRUE)
}
dc1 <- death_counts_df(quartile_1)
dc1

#35.48% of these deaths are from tumors
percent_tumor_deaths <- sum(dc1$total_tumors)/sum(dc1$occurrences)
percent_tumor_deaths
#21.03% of these deaths are from non-tumors
percent_non_tumor_deaths <- sum(dc1$total_non_tumors)/sum(dc1$occurrences)
percent_non_tumor_deaths
#43.48% of these deaths are from lymphomas
percent_lymphoma_deaths <- sum(death_counts$total_lymphoma)/sum(death_counts$occurrences)
percent_lymphoma_deaths

######################### male vs. female ####################################
data_to_use_f <- subset(data_to_use, sex == "F")
data_to_use_m <- subset(data_to_use, sex == "M")

death_counts_f <- data.frame(occurrences = colSums(data_to_use_f[,25:ncol(data_to_use_f)-5]))
death_counts_f$tumor <- substr(rownames(death_counts_f), 1, 1) =="T"
lymphoma <- c("NTYG_L","NTYL_L", "TTYG_L", "TTYL_L")
death_counts_f$lymphoma <- rownames(death_counts_f) %in% lymphoma
death_counts_f$total_tumors <- death_counts_f$occurrences * death_counts_f$tumor * !(death_counts_f$lymphoma)
death_counts_f$total_non_tumors <- (death_counts_f$occurrences * !(death_counts_f$tumor)) * !(death_counts_f$lymphoma)
death_counts_f$total_lymphoma <- death_counts_f$occurrences * death_counts_f$lymphoma
death_counts_f <- ordered(death_counts_f, "occurrences", TRUE)

death_counts_m <- data.frame(occurrences = colSums(data_to_use_m[,25:ncol(data_to_use_m)-5]))
death_counts_m$tumor <- substr(rownames(death_counts_m), 1, 1) =="T"
lymphoma <- c("NTYG_L","NTYL_L", "TTYG_L", "TTYL_L")
death_counts_m$lymphoma <- rownames(death_counts_m) %in% lymphoma
death_counts_m$total_tumors <- death_counts_m$occurrences * death_counts_m$tumor * !(death_counts_m$lymphoma)
death_counts_m$total_non_tumors <- (death_counts_m$occurrences * !(death_counts_m$tumor)) * !(death_counts_m$lymphoma)
death_counts_m$total_lymphoma <- death_counts_m$occurrences * death_counts_m$lymphoma
death_counts_m <- ordered(death_counts_m, "occurrences", TRUE)

######################### control vs. gamma irradiated ####################################
#check that all controls with 0Gy are listed as C in radn group
summary(data_to_use$radn)
data_to_use %>%
      group_by(radn, dose_group2) %>%
      summarize(n)

data_to_use_c <- subset(data_to_use, radn == "C")
data_to_use_g <- subset(data_to_use, radn == "G")

death_counts_c <- data.frame(occurrences = colSums(data_to_use_c[,25:ncol(data_to_use_c)-5]))
death_counts_c$tumor <- substr(rownames(death_counts_c), 1, 1) =="T"
lymphoma <- c("NTYG_L","NTYL_L", "TTYG_L", "TTYL_L")
death_counts_c$lymphoma <- rownames(death_counts_c) %in% lymphoma
death_counts_c$total_tumors <- death_counts_c$occurrences * death_counts_c$tumor * !(death_counts_c$lymphoma)
death_counts_c$total_non_tumors <- (death_counts_c$occurrences * !(death_counts_c$tumor)) * !(death_counts_c$lymphoma)
death_counts_c$total_lymphoma <- death_counts_c$occurrences * death_counts_c$lymphoma
death_counts_c <- ordered(death_counts_c, "occurrences", TRUE)

death_counts_g <- data.frame(occurrences = colSums(data_to_use_g[,25:ncol(data_to_use_g)-5]))
death_counts_g$tumor <- substr(rownames(death_counts_g), 1, 1) =="T"
lymphoma <- c("NTYG_L","NTYL_L", "TTYG_L", "TTYL_L")
death_counts_g$lymphoma <- rownames(death_counts_g) %in% lymphoma
death_counts_g$total_tumors <- death_counts_g$occurrences * death_counts_g$tumor * !(death_counts_g$lymphoma)
death_counts_g$total_non_tumors <- (death_counts_g$occurrences * !(death_counts_g$tumor)) * !(death_counts_g$lymphoma)
death_counts_g$total_lymphoma <- death_counts_g$occurrences * death_counts_g$lymphoma
death_counts_g <- ordered(death_counts_g, "occurrences", TRUE)



death_counts <- data.frame(occurrences = colSums(data_to_use[,21:ncol(data_to_use)-1]))
death_counts$tumor <- substr(rownames(death_counts), 1, 1) =="T"
lymphoma <- c("NTYG_L","NTYL_L", "TTYG_L", "TTYL_L")
death_counts$lymphoma <- rownames(death_counts) %in% lymphoma
death_counts$total_tumors <- death_counts$occurrences * death_counts$tumor * !(death_counts$lymphoma)
death_counts$total_non_tumors <- (death_counts$occurrences * !(death_counts$tumor)) * !(death_counts$lymphoma)
death_counts$total_lymphoma <- death_counts$occurrences * death_counts$lymphoma
death_counts <- ordered(death_counts, "occurrences", TRUE)

row.names(death_counts)[2]
death_counts$female <- 0

#35.48% of these deaths are from tumors
#30.26% for females
#43.78% for males
#34.94% for Controls
#35.78% for Gamma
percent_tumor_deaths <- sum(death_counts_g$total_tumors)/sum(death_counts_g$occurrences)
percent_tumor_deaths
#21.03% of these deaths are from non-tumors
#23.55% for females
#17.03% for males
#18.60% for controls
#22.38% for Gamma
percent_non_tumor_deaths <- sum(death_counts_g$total_non_tumors)/sum(death_counts_g$occurrences)
percent_non_tumor_deaths
#43.48% of these deaths are from lymphomas
#46.19% for females
#39.18% for males
#46.45% for controls
#41.83% for Gamma
percent_lymphoma_deaths <- sum(death_counts_g$total_lymphoma)/sum(death_counts_g$occurrences)
percent_lymphoma_deaths


#graph that I want
barplot(death_counts_g[1:30,1], 
        names.arg = rownames(death_counts_g[1:30,]), 
        las=2)
nrow(death_counts) #136
cancer_incidences <- death_counts_g[death_counts_g$total_tumors > 0,]
nrow(cancer_incidences) #48
noncancer_incidences <- death_counts_g[death_counts_g$total_non_tumors > 0,]
nrow(noncancer_incidences) #58
lymphoma_incidences <- death_counts_g[death_counts_g$total_lymphoma > 0,]
nrow(lymphoma_incidences) #4
#136-48-58-4=26 diseases that had no occurrences

barplot(cancer_incidences[1:25,1], 
        names.arg = rownames(cancer_incidences[1:25,]), 
        las=2)

barplot(cancer_incidences[26:nrow(cancer_incidences),1], 
        names.arg = rownames(cancer_incidences[26:nrow(cancer_incidences),]), 
        las=2)

barplot(noncancer_incidences[2:31,1], 
        names.arg = rownames(noncancer_incidences[2:31,]), 
        las=2)

barplot(noncancer_incidences[32:nrow(noncancer_incidences),1], 
        names.arg = rownames(noncancer_incidences[32:nrow(noncancer_incidences),]), 
        las=2)

barplot(lymphoma_incidences$occurrences, names.arg = rownames(lymphoma_incidences), 
        las=2)
text(x= 1:4, 
     y= lymphoma_incidences$occurrences + 100, 
     labels=as.character(lymphoma_incidences$occurrences),
     xpd=TRUE)


#6/9/17
####################################################################################
lymphoma <- c("NTYG_L","NTYL_L", "TTYG_L", "TTYL_L")
cdu <- "CDU_L"

#get list of tumor names
tumor_names <- substr(rownames(death_counts), 1, 1) == "T"
tumor_names <- rownames(death_counts)[tumor_names]
tumor_names <- tumor_names[!(tumor_names %in% lymphoma)]

#get list of nontumor names
nontumor_names <- rownames(death_counts)
nontumor_names <- nontumor_names[!(nontumor_names %in% lymphoma)]
nontumor_names <- nontumor_names[!(nontumor_names %in% tumor_names)]
nontumor_names <- nontumor_names[!(nontumor_names %in% cdu)]

#create class column to describe cause of death
data_to_use$class <- 0
location <- apply(data_to_use[,20:156], 1, function(r) which(r == "1"))
for(i in 1:length(location)){     
      column <- location[i] + 19
      if(colnames(data_to_use)[column] %in% tumor_names){
            data_to_use[i, 157] <- "tumor"
      }
      else if(colnames(data_to_use)[column] %in% nontumor_names){
            data_to_use[i, 157] <- "nontumor"
      }
      else if(colnames(data_to_use)[column] %in% lymphoma){
            data_to_use[i, 157] <- "lymphoma"
      }
      else {data_to_use[i, 157] <- "CDU"}
}

#create quartile column to describe time of death
data_to_use$quartile <- 0
age_quartiles <- summary(data_to_use$age)
age_quartiles<- age_quartiles[-4] #remove mean
for(i in 1:nrow(data_to_use)){     
      if(data_to_use$age[i] <= age_quartiles[2]){
            data_to_use$quartile[i] <- 1
      }
      else if(data_to_use$age[i] > age_quartiles[2] & 
         data_to_use$age[i] <= age_quartiles[3]){
            data_to_use$quartile[i] <- 2
      }
      else if(data_to_use$age[i] > age_quartiles[3] & 
              data_to_use$age[i] <= age_quartiles[4]){
            data_to_use$quartile[i] <- 3
      }
      else data_to_use$quartile[i] <- 4
      
}
#create lots of plots to represent data
data_to_use <- droplevels(data_to_use)
#age at death vs. cause of death & statistical analysis
boxplot(age~class, data = data_to_use, main="Q4", 
        ylab="Age at death", las = 2, subset = sex == "F")
boxplot(age~class, data = data_to_use, main="Age at death vs. cause of death", 
        ylab="Age at death", las = 2)
fit_anova <- (lm(age~factor(class), data = data_to_use))
summary(fit_anova)
anova(fit_anova)
summarized <- data_to_use %>%
      group_by(class) %>%
      summarise(avg_age = mean(age), sd_age = sd(age), n = n())
ggplot(summarized, aes(x=class, y=avg_age)) + 
      geom_bar(stat="identity") +
      geom_errorbar(aes(ymin=avg_age-sd_age, ymax=avg_age+sd_age),
                    width=.2)+
      xlab("Class") +
      ylab("Age at death") +
      ggtitle("Age at death for different causes of death") +
      theme_bw()
#age at death vs. type of treatment & statistical analysis
boxplot(age~radn, data = data_to_use, main="Age at death vs. type of treatment", 
        ylab="Age at death",las = 1)
t.test(age~radn, data = data_to_use)
boxplot(age~radn*class, data = data_to_use, main="Age at death vs. type of treatment and cause of death", 
        ylab="Age at death",las = 2)
fit_anova2 <- lm(age~radn * class, data = data_to_use)
summary(fit_anova2)
anova(fit_anova2)
boxplot(age~radn*class*quartile, data = data_to_use, main="Age at death vs. type of treatment, quartile, and class", 
        ylab="Age at death",las = 2)
class_freq_tab <- table(subset(data_to_use, class == "tumor")$quartile)
class_freq_tab <- class_freq_tab[order(class_freq_tab)]
percent_occurrence <- round(prop.table(class_freq_tab)*100,2)
barplot(class_freq_tab, ylab = "Occurrences", main = "Tumor", xlab = "Quartile")
text(x= 1:4, 
     y = class_freq_tab+20,
     labels=as.character(percent_occurrence),
     xpd=TRUE)
cdu_animals <- filter(data_to_use,
                      class == "CDU")
tumor_animals <- filter(data_to_use,
                        class == "tumor")
nontumor_animals <- filter(data_to_use,
                           class == "nontumor")
lymphoma_animals <- filter(data_to_use,
                           class == "lymphoma")


#6/19/17
####################################################################################
#add dose-group info so that I can treat them as factors and group them together

data_to_use$dose_group <- 0
for(i in 1:nrow(data_to_use)){     
      if(data_to_use$total_dose[i] <= 0){
            data_to_use$dose_group[i] <- 0
      }
      else if(data_to_use$total_dose[i] >= 190 &
              data_to_use$total_dose[i] <= 210){
            data_to_use$dose_group[i] <- 2
      }
      else if(data_to_use$total_dose[i] >= 390 &
              data_to_use$total_dose[i] <= 410){
            data_to_use$dose_group[i] <- 4
      }
      else data_to_use$dose_group[i] <- NA
      
}

#filter data for the dose groups I'm planning to look into

#get total n for all animals with these treatment conditions
df1 <- filter(data_to_use,
              #dose_group == 0 | dose_group == 2 | dose_group == 4,
              sex == "M",
              expt != 8)

df1_summary <- df1 %>%
      group_by(fractions, dose_group2) %>%
      summarise(n_group = n())

#get n for mice with these conditions that died of lung cancer or lymphoma
df2 <- filter(df1,
              #NTYG_L == 1)
              #TADN_L == 1)
              NTYG_L == 1 | NTYL_L == 1)

#filter out small sample size and exact conditions for ANOVA
df2 <- filter(df2,
              !(fractions == 60 & dose_group2 == 49.01))
              #!(fractions == 60 & dose_group == 4),
            


df2 <- filter(df2,
              !(fractions == 24 & dose_rate == 1.703780 & dose_group2 == 49.01),
              !(fractions == 24 & dose_group2 == 16.40),
              !(fractions == 24 & dose_rate == 0.182889 & dose_group2 == 4.00),
              #dose_group == 2 | dose_group == 4)
              dose_group2 != 0)

#summary to get a feel for data and use to calculate percentage later
df2_summary <- df2 %>%
      group_by(fractions, dose_group2) %>%
      summarise(n_group = n(),
                avg_age = mean(age),
                avg_dose = mean(total_dose),
                avg_rate = mean(dose_rate))

#make model and do statistical analysis on it
uw.model <- lm(age ~ factor(fractions) + factor(dose_group), subset(df2, dose_group2 != 0))
summary(uw.model)
Anova(uw.model)

ggplot(subset(df2, dose_group2 !=0),
      aes(factor(dose_group2), age, fill = factor(fractions)))+
      geom_boxplot() + 
      ggtitle("Age vs. Dose")

#create new data frame with same info as df2 that I can add on to for percentage
df3 <- data.frame(age = df2$age, dose_group = df2$dose_group, fractions = df2$fractions)
df2_summary$dose_group <- factor(df2_summary$dose_group)
#calculate percentage
df2_summary$percentage <- df2_summary$n_group/df1_summary$n_group

#plot percentages - change type of disease cuasing death each time
ggplot(df2_summary, aes(factor(dose_group2), percentage, fill = factor(fractions))) + 
      geom_bar(stat="identity", position = "dodge") +
      ggtitle("Percentage of deaths due to lung cancer")
barchart(percentage~dose_group,data=df2_summary,groups=fractions, 
         scales=list(x=list(rot=90,cex=0.8)))
boxplot(1/age~dose_group + fractions, data = df2, main="age vs. dose", 
        ylab="Age at death", las = 2, subset = quartile ==4)
ggplot(df2, aes(factor(dose_group), age, fill = factor(fractions)))+
      geom_boxplot() + 
      ggtitle("Age vs. Dose")
ggplot(subset(df2, (dose_group2 == 0 | dose_group2 == 4 | dose_group2 == 0)),
       aes(factor(rate_group), age, fill = factor(fractions)))+
      geom_boxplot() + 
      ggtitle("Age vs. Dose")

ggplot(df2,
       aes(factor(rate_group), age, fill = factor(fractions)))+
      geom_boxplot() + 
      ggtitle("Age vs. Dose")

uw.model <- lm(age ~ factor(fractions), df2)

check <-df2 %>%
      group_by(fractions, dose_group2) %>%
      summarise(n = n(), mean_age = mean(age))

check2 <-df2 %>%
      group_by(rate_group, fractions, dose_group2, dose_rate, total_dose) %>%
      summarise(n = n(), mean_age = mean(age))
cor(df2$dose_group2, df2$rate_group)

#try to create graphs for lymphoma death
nty_female <- filter(data_to_use,
              NTYG_L == 1 | NTYL_L == 1,
              expt != 8,
              sex == "F")
nty_female <- select(nty_female,
              animal_id:dose_rate,
              NTYG_L,
              NTYL_L,
              class:dose_group)

summary_nty_female <- nty_female %>%
      group_by(expt, fractions, dose_rate, total_dose, dose_group) %>%
      summarise(n=n())

data_to_use_female <- filter(data_to_use,
                             expt != 8,
                             sex == "F")

summary_female <- data_to_use_female %>%
                  group_by(expt, fractions, dose_rate, total_dose) %>%
                  summarise(n=n())

joint_nty_f <- semi_join(summary_female, summary_nty_female,
          by = c("expt", "fractions", "dose_rate", "total_dose"))

summary_nty_female$percent <- summary_nty_female$n/joint_nty_f$n
summary_nty_female_filtered <- filter(summary_nty_female,
                                      total_dose == 0 | total_dose == 197.55 | total_dose == 399.90 | total_dose == 200 | total_dose == 206)
ggboxplot(summary_nty_female_filtered, x = "total_dose", y = "percent", color = "fractions")
barplot(summary_nty_female_filtered$percent)

#7/5/17
####################################################################################
#check to see if first irradiated time is different for control vs. gamma treated mice
check <- filter(data_to_use, expt != 8) %>%
      group_by(dose_group2) %>%
      summarise(n = n(), mean = mean(first_irrad))
cor(df2$age, df2$fractions)
#add new dose groups for bigger analysis
data_to_use$dose_group2 <- 0
for(i in 1:nrow(data_to_use)){     
      if(data_to_use$total_dose[i] <= 0){
            data_to_use$dose_group2[i] <- 0
      }
      else if(data_to_use$total_dose[i] > 0 &
              data_to_use$total_dose[i] <= 150){
            data_to_use$dose_group2[i] <- 1.5
      }
      else if(data_to_use$total_dose[i] > 150 &
              data_to_use$total_dose[i] <= 400){
            data_to_use$dose_group2[i] <- 4
      }
      else if(data_to_use$total_dose[i] > 400 &
              data_to_use$total_dose[i] <= 760){
            data_to_use$dose_group2[i] <- 7.6
      }
      else if(data_to_use$total_dose[i] > 760 &
              data_to_use$total_dose[i] <= 1640){
            data_to_use$dose_group2[i] <- 16.4
      }
      else data_to_use$dose_group2[i] <- 49.01
      
}

data_to_use$rate_group <- 0
for(i in 1:nrow(data_to_use)){     
      if(data_to_use$dose_rate[i] <= 0){
            data_to_use$rate_group[i] <- 0
      }
      else if(data_to_use$dose_rate[i] > 0 &
              data_to_use$dose_rate[i] <= .5){
            data_to_use$rate_group[i] <- .5
      }
      else if(data_to_use$dose_rate[i] > .5 &
              data_to_use$dose_rate[i] <= 2){
            data_to_use$rate_group[i] <- 2
      }
      else if(data_to_use$dose_rate[i] > 2 &
              data_to_use$dose_rate[i] <= 5){
            data_to_use$rate_group[i] <- 5
      }
      else data_to_use$rate_group[i] <- 27.3
      
}

jm8 <- filter(data_to_use,
              expt == 8)

check <- jm8 %>%
      group_by(dose_rate) %>%
      summarise(n = n(), mean_age = mean(age))


summary8 <- jm8 %>% group_by(dose_group2, dose_rate) %>%
                     summarise(n = n(), mean_age = mean(age))

#5/5/17
####################################################################################



common_cancers <- names(head(sort(death_counts[1,], decreasing = TRUE), 20))

cancer_top <- select(data_to_use,
                     animal_id:dose_rate,
                     match(common_cancers,names(data_to_use)))

common_cancers_symbol <- lapply(common_cancers, as.symbol)            
cancer_top %>%
      group_by(NTYG_L, TADN_L, CDU_L, TVAS_L, ANE_L, TOVE_L) %>%
      summarise(tot = n(),
                avg_dose = mean(total_dose),
                avg_dose_rate = mean(dose_rate),
                avg_age = mean(age))
#which cancers are associated with higher doses/dose rates/age at death?
#histogram of average age if grouped by types of cancer
#histogram of average dose rate for each cancer
#histogram of average total dose for each cancer

#if == 1 for type of cancer, average age, average dose, average dose rate all get a column...
#add these rows to my death counts dataframe?

death_counts$age <- colSums(data_to_use[,21:ncol(data_to_use)-1]))
      
path_key <- read.csv('Lab_work/pathologies_key.csv', header=TRUE, skip = 1)

groups <- path_key[path_key$observation_category == "micro",]
head(groups$e_group)
head(unique(groups$e_group))













##### Reproduce Ben's figure 13 panels 5 and 13 and data from table 8

#start with janus_7, then move on to control, acute, protracted
janus_7 <- all_info[all_info$expt == 7,]
#animals cannot have died before all exposures were completed
janus_7 <- janus_7[janus_7$age >= 520,]
#animals cannot have been moved to another experiment or killed early
janus_7 <- janus_7[janus_7$cause_of_death == "Died" |
                         janus_7$cause_of_death == "Sacrifice, moribund",]

######## Ctrl - MALE
janus_7_ctrl_m <- janus_7[janus_7$total_dose == 0 & janus_7$sex == "M",]
#I have 3 more animals than Ben... not sure why. Mean is off a little, sem matches
summary(janus_7_ctrl_m)
sem <- sd(janus_7_ctrl_m$age)/sqrt(length(janus_7_ctrl_m$age))

####### Acute and protracted - MALE
janus_7_g <- janus_7[janus_7$radn == "G",]
janus_7_g_m <- janus_7_g[janus_7_g$sex == "M",]
janus_7_g_m_a <- janus_7_g_m[janus_7_g_m$fractions == 1,]
janus_7_g_m_p <- janus_7_g_m[janus_7_g_m$fractions > 1,]
#not sure why it's only this dose rate? numbers match though
janus_7_g_m_a <- janus_7_g_m_a[janus_7_g_m_a$dose_rate > 9.87 &
                                     janus_7_g_m_a$dose_rate <10,]
janus_7_g_m_p <- janus_7_g_m_p[janus_7_g_m_p$dose_rate == .1480,]

summary(janus_7_g_m_a)
summary(janus_7_g_m_p)

df_7_m <- rbind(janus_7_ctrl_m, janus_7_g_m_a, janus_7_g_m_p)
fit <- survfit(Surv(age) ~ total_dose, data = df_7_m)
#ggsurvplot(fit, linetype = c(1,2,3)) +
ggsurvplot(fit) +      
      scale_colour_brewer(palette = "Oranges")
      scale_color_manual(values=c("black", "grey30", "red"))

############# janus 7 --- FEMALE
#control
janus_7_ctrl_f <- janus_7[janus_7$total_dose == 0 & 
                                janus_7$sex == "F",]
#acute
summary(janus_7_ctrl_f)
janus_7_g_f_a <- janus_7_g[janus_7_g$sex == "F" & 
                                 janus_7_g$fractions ==1 &
                                 janus_7_g$dose_rate > 9.8 &
                                 janus_7_g$dose_rate < 10,]
summary(janus_7_g_f_a)
#protracted
janus_7_g_f_p <- janus_7_g[janus_7_g$sex == "F" & 
                                 janus_7_g$fractions > 1 &
                                 janus_7_g$dose_rate == .1480,]
summary(janus_7_g_f_p)

#combine all (ctrl, acute, protracted)
df_7_f <- rbind(janus_7_ctrl_f, janus_7_g_f_a, janus_7_g_f_p)
#survival curves
fit <- survfit(Surv(age) ~ total_dose, data = df_7_f)
ggsurvplot(fit)

################### Compare survival curves for all controls in Janus
demographics <- demographics[demographics$cause_of_death == "Died" |
                               demographics$cause_of_death == "Sacrifice, moribund",]
controls <- demographics[demographics$total_dose == 0,]

JM9 <- filter(controls,
              expt == 9)
JM9 %>%
      group_by(fractions) %>%
      summarise(n=n())

#remove JM10/11/12
controls <- controls[controls$expt != 10 &
                           controls$expt != 11 &
                           controls$expt != 12,]

#remove JM8
controls <- filter(controls, expt != 8)

#remove JM2
controls <- filter(controls, expt != 2)

#remove JM9 male mice
controls <- filter(controls, 
                   !(expt == 9 & sex == "M"))

#remove JM4 fractions == 120 and fractions == 300 and males
controls <- filter (controls,
                    !(expt ==4 & fractions == 120),
                    !(expt ==4 & fractions == 300))
                    #!(expt ==4 & sex == "M"))

controls <- filter(controls,
                   radn =="G" |radn=="C")
#get summaries of controls by expt, fractions, sex, etc.
conditions_we_will_use <- controls %>%
                              group_by(expt, fractions, total_dose, dose_rate, sex) %>%
                              summarise(tot_n = n())
mice_to_use <- filter(mice_to_use,
                      expt != 8)
conditions_we_will_use <- mice_to_use %>%
      group_by(expt, fractions, total_dose, dose_rate, sex) %>%
      summarise(tot_n = n())

write.csv(conditions_we_will_use, file = "conditions.csv",
          row.names=TRUE)
fit <- survfit(Surv(age) ~ fractions+expt, data = controls)
ggsurvplot(fit)
survdiff(Surv(age) ~ fractions+expt, data = controls)

#make lots of comparisons - start with control that excludes everything
#that I already established needs to be excluded

#try to compare JM13 to all others... had difficulty figuring out how
#http://stackoverflow.com/questions/6518456/loop-over-rows-of-dataframe-applying-function-with-if-statement

controls <- within(controls, comparison <- 1) 
controls[controls$expt==9, "comparison"] <- 0
     
comparison <- filter(controls, 
                     expt == 9 | expt == 13,
                     sex == "F")

comparison %>%
      group_by(expt, fractions, sex) %>%
      summarise(tot_n = n())


fit <- survfit(Surv(age) ~ expt, data = comparison)
ggsurvplot(fit)
survdiff(Surv(age) ~ expt, data = comparison)



#change font - ggsurvplot(fit, font.legend = c(5, "plain", "black"))

summary(controls)

controls %>%
      group_by(fractions, expt, sex) %>%
      summarise(tot_n = n(),
                mean(first_irrad),
                max(first_irrad),
                min(first_irrad))

#trick Ben showed me
#table(paste(controls$expt, controls$was_control_mock_treated))

########################################################################
#work with gamma treatments now

demographics <- demographics[demographics$cause_of_death == "Died" |
                            demographics$cause_of_death == "Sacrifice, moribund",]
gamma <- filter(demographics, 
                radn =="G", 
                expt == 3 | expt == 7 | expt == 9 | expt == 13 | expt == 14,
                total_dose <= 200)

gamma %>%
      group_by(expt, fractions, total_dose, dose_rate)%>%
      summarise(total = n(),
                avg_age = mean(age))

data_to_use <- filter(demographics,
                      radn =="G" |radn=="C",
                      expt == 3 | expt == 7 | expt == 9 | expt == 13 | expt == 14,
                      total_dose <= 200)
data_to_use %>%
      group_by(expt, fractions, dose_rate, total_dose)%>%
      summarise(n = n(),
                ave_age = mean(age))

JM13 <- filter(demographics,
              radn =="G"|radn=="C",
              expt == 13)
JM14 <- filter(demographics, 
              expt ==14,
              radn == "G"| radn=="C")
JM3 <- filter(demographics, 
              expt ==3,
              radn == "G"|radn =="C")
JM8 <- filter(demographics, 
              expt ==8,
              radn == "G"|radn =="C")
JM4 <- filter(demographics, 
              expt ==4,
              radn == "G"|radn =="C")
JM8_controls <- filter(JM8,
                       total_dose == 0)

JM7 %>%
      group_by(fractions, dose_rate, total_dose) %>%
      summarise(total = n(),
                ave_age = mean(age))

controls7 <- filter(controls,
             expt == 7)
controls3 <- filter(controls,
                    expt == 3)
controls3 %>%
      group_by(fractions) %>%
      summarise(total = n(),
                ave_age = mean(age))
JM14 %>%
      group_by(dose_rate, sex, fractions, tmt) %>%
      summarise(total = n())

table(paste(gamma$expt, gamma$fractions, gamma$total_dose))

gamma4 <- filter(demographics,
                  expt==4,
                  total_dose<=200)
table(paste(gamma4$expt, gamma4$fractions, gamma4$total_dose))

gamma_2gy <- filter(gamma,
                    total_dose <205,
                    total_dose >195)
      
fit <- survfit(Surv(age) ~ dose_rate, data = JM8)
ggsurvplot(fit)
survdiff(Surv(age) ~ fractions + expt + sex, data = gamma_2gy)









to_plot_13 <- ddply(JM13, .(total_dose, fractions, dose_rate), summarise,
                    mean_age = mean(age),
                    sem = sd(age)/sqrt(length(age)))
to_plot_3 <- ddply(JM3, .(total_dose, fractions, dose_rate), summarise,
                   mean_age = mean(age),
                   sem = sd(age)/sqrt(length(age)))

to_plot_13 <- rbind(to_plot_13, to_plot_13[rep(1, each=4),])
to_plot_3 <- rbind(to_plot_3, to_plot_3[rep(1, each=4),])

to_plot_13 <- to_plot_13 %>% 
      mutate(group = c(1, 1, 2, 3, 4, 5, 2, 3, 4, 5))
to_plot_3 <- to_plot_3 %>% 
      mutate(group = c("a1", "a1", "a2", "a3", "a4", "a5", "a2", "a3", "a4", "a5"))

acute_vs_fraction <- rbind(to_plot_13, to_plot_3)
ggplot(acute_vs_fraction,aes(total_dose,1/mean_age,color=factor(group)))+geom_point()+
      geom_smooth(data=acute_vs_fraction,
                  aes(total_dose,1/mean_age,color=factor(group)),method=lm,se=FALSE)

acute_vs_fraction %>% 
      group_by(group) %>% 
      do({
            mod = lm(1/mean_age ~ total_dose, data = .)
            data.frame(Intercept = coef(mod)[1],
                       Slope = coef(mod)[2])
      })


ggplot(aes(x=total_dose, y=1/mean_age), 
       data = to_plot_13)+
      geom_point()+
      geom_smooth(data = subset(to_plot_13, total_dose == 0 | total_dose == 100),
                  method = "lm", se = FALSE, aes(colour="Gy1"), size = .5)+
      geom_smooth(data = subset(to_plot_13, total_dose == 0 | total_dose == 200),
                  method = "lm", se = FALSE, aes(colour="Gy2"), size = .5)+
      geom_smooth(data = subset(to_plot_13, total_dose == 0 | total_dose == 300),
                  method = "lm", se = FALSE, aes(colour="Gy3"), size = .5)+
      geom_smooth(data = subset(to_plot_13, total_dose == 0 | total_dose == 450),
                  method = "lm", se = FALSE, aes(colour="Gy4.5"), size = .5)+
      geom_smooth(data = subset(to_plot_13, total_dose == 0 | total_dose == 600),
                  method = "lm", se = FALSE, aes(colour="Gy6"), size = .5)+
      scale_colour_manual(name="Line Color",
                          values=c(Gy1="red", Gy2="orange", Gy3="green", 
                                   Gy4.5="blue", Gy6="purple"))+


########################################################################

#trying to compare janus and ERA

data <- read.csv('janus/data/era/big.csv', sep='|', as.is=TRUE, na.strings = "n/a")
Janus2 <- data[grep("1003-20", data$Individual.ID),]
tail(Janus2)
nrow(Janus2)
summary(Janus2)

head(Janus2[Janus2$Treatment.Age == 106,])
head(janus_2[janus_2$first_irrad == 106,])
summary(janus_2$first_irrad)
summary(Janus2$Treatment.Age)

Janus3 <- data[grep("1003-21", data$Group.ID),]
tail(Janus3)
nrow(Janus3)
summary(Janus3)

Janus4 <- data[grep("1003-22", data$Group.ID),]
tail(Janus4)
nrow(Janus4)
summary(Janus4)

Janus7 <- data[grep("1003-24", data$Group.ID),]
tail(Janus7)
nrow(Janus7)
summary(Janus7)

Janus8 <- data[grep("1003-25", data$Group.ID),]
Janus9 <- data[grep("1003-26", data$Group.ID),]
Janus10 <- data[grep("1003-27", data$Group.ID),]
#Janus11 <- data[grep("1003-", data$Group.ID),], doesn't exist
Janus12 <- data[grep("1003-28", data$Group.ID),]

Janus13 <- data[grep("1003-29", data$Group.ID),]
Janus14 <- data[grep("1003-30", data$Group.ID),]

Janus12 %>%
      group_by(Strain) %>%
      summarise(tot_n = n())
``

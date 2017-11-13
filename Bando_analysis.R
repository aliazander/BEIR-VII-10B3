#get total n for all animals with these treatment conditions
df1 <- filter(data_to_use,
              #!(dose_group == 2 & fractions == 24),
              #sex == "M",
              #dose_group != "NA",
              expt != 8)

df1_summary <- df1 %>%
      group_by(fractions, dose_group2) %>%
      summarise(n_group = n(),
                avg_age = mean(age),
                avg_dose = mean(total_dose),
                avg_rate = mean(dose_rate))

#get n for mice with these conditions that died of lung cancer or lymphoma
df2 <- filter(df1,
              #NTYG_L == 1)
              #TADN_L == 1)
              #NTYG_L == 1 | NTYL_L == 1
              class == "tumor")
              #dose_group != "NA")

#filter out small sample size and exact conditions for ANOVA

df2 <- filter(df2,
              !(fractions == 24 & dose_group == 2))
 

#summary to get a feel for data and use to calculate percentage later
df2_summary <- df2 %>%
      group_by(fractions, dose_group2) %>%
      summarise(n_group = n(),
                avg_age = mean(age),
                avg_dose = mean(total_dose),
                avg_rate = mean(dose_rate))

#make model and do statistical analysis on it
uw.model <- lm(age ~ factor(fractions) * factor(dose_group2), subset(df2, dose_group2 != 0))
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
ggplot(df2_summary, aes(factor(dose_group), percentage, fill = factor(fractions))) + 
      geom_bar(stat="identity", position = "dodge") +
      ggtitle("Percentage of deaths due to cancer")

uw.model <- lm(age ~ factor(fractions) * factor(dose_group2), data = df2)
uw.model <- lm(age ~ factor(fractions) * factor(dose_group2), data = subset(df2, dose_group2 !=0))
uw.model <- lm(age ~ factor(fractions), subset(df2, dose_group2 == 0))

summary(uw.model)
Anova(uw.model)



ggplot(df2,
       aes(factor(dose_group2), age, fill = factor(fractions)))+
      geom_boxplot() + 
      ggtitle("Age vs. Dose")

ggplot(subset(df2, dose_group2 ==0),
       aes(factor(dose_group), age, fill = factor(fractions)))+
      geom_boxplot() + 
      ggtitle("Age vs. Dose")

boxplot(age~radn, data = subset(data_to_use, class == "tumor"), 
        main="Age at death vs. type of treatment for lethal cancers", 
        ylab="Age at death",las = 1)
t.test(age~radn, data = subset(data_to_use, class == "tumor")
       
       
survdiff(Surv(age) ~ factor(dose_group2) * factor(fractions), data = df2)
survdiff(Surv(age) ~ dose_group, data = df2)
fit <- survfit(Surv(age) ~ dose_group, data = df2)
ggsurvplot(fit)
survdiff(Surv(age) ~ dose_group, data = df2)

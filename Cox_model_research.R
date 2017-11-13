install.packages("rms")
library(OIsurv)
library(plyr)
library(dplyr)
library(ggplot2)
library(stats4)
library(survival)
library(survminer)
library(ggpubr)
library(lattice)
library(car)
library(broom)
require(broom)
library(rms)
#resources
#https://www.r-bloggers.com/cox-proportional-hazards-model/
#http://www.sthda.com/english/wiki/cox-model-assumptions


fit <- survfit(Surv(age) ~ dose_rate, data = JM8)
ggsurvplot(fit)
survdiff(Surv(age) ~ dose_rate, data = JM8)


########## tutorial and playing around with it #######
data("lung")
head(lung,30)
dim(lung)
summary(lung$sex == 1)
lung.changed1 <- filter(lung,
                       sex == 1)
lung.changed2 <- lung.changed1
#see how changing these ages changes the HR to understand meaning
lung.changed2$sex <- 2
lung.changed2$time <- lung.changed1$time*2
lung.changed <- rbind(lung.changed1, lung.changed2)

res.cox <- coxph(Surv(time, status) ~ factor(sex), data = lung.changed)
summary(res.cox)
ggsurvplot(survfit(res.cox), color = "#2E9FDF",
           ggtheme = theme_minimal())

res.cox <- coxph(Surv(time, status) ~ sex, data = lung)
summary(res.cox)

#back to their example
res.cox <- coxph(Surv(time, status) ~ age + sex + ph.ecog, data =  lung)
summary(res.cox)
ggsurvplot(survfit(res.cox), color = "#2E9FDF",
           ggtheme = theme_minimal())
# Create the new data  
sex_df <- with(lung,
            data.frame(sex = c(1, 2), 
                          age = rep(mean(age, na.rm = TRUE), 2),
                          ph.ecog = c(1, 1)
               )
)
sex_df
fit <- survfit(res.cox, newdata = sex_df)
ggsurvplot(fit, conf.int = TRUE, legend.labs=c("Sex=1", "Sex=2"),
           ggtheme = theme_minimal())
#check residuals
res.cox <- coxph(Surv(time, status) ~ age + sex + wt.loss, data =  lung)
res.cox

test.ph <- cox.zph(res.cox)
test.ph
ggcoxzph(test.ph)

#check outliers
ggcoxdiagnostics(res.cox, type = "dfbeta",
                 linear.predictions = FALSE, ggtheme = theme_bw())
ggcoxdiagnostics(res.cox, type = "deviance",
                 linear.predictions = FALSE, ggtheme = theme_bw())
#check nonlinearity
ggcoxfunctional(Surv(time, status) ~ age + log(age) + sqrt(age), data = lung)
fit <- coxph(Surv(time, status) ~ age + log(age) + sqrt(age), data = lung)
ggcoxfunctional(fit, data = lung)

#################### now with my data #################################
setwd("/Users/g-woloschak/Documents")
macros <- read.csv('Lab_work/macro_pathologies.csv', header=TRUE, skip = 1)
macros[is.na(macros)] <- 0

demographics <- read.csv('Lab_work/demographics.csv', header=TRUE, skip = 1)

all_info <- cbind(demographics, macros)
head(demographics)
# 0 = alive, 1 = dead
demographics$status <- 0
for (i in 1:nrow(demographics)){
      if(demographics$cause_of_death[i] == "Died" | demographics$cause_of_death[i] == "Sacrifice, moribund")
            demographics$status[i] <- 1
}

res.cox <- coxph(Surv(age, status) ~ sex, data = demographics)
summary(res.cox)
ggsurvplot(survfit(res.cox), color = "#2E9FDF",
           ggtheme = theme_minimal())

sex_df <- with(demographics,
               data.frame(sex = c("F", "M")
               )
)
sex_df
#new method with censoring
res.cox <- coxph(Surv(age, status) ~ radn, data = demographics)
summary(res.cox)
fit <- survfit(res.cox, newdata = sex_df)
ggsurvplot(fit, conf.int = TRUE, legend.labs=c("Sex=F", "Sex=M"),
           ggtheme = theme_minimal())
#new method without censoring
res.cox <- coxph(Surv(age) ~ radn, data = demographics)
summary(res.cox)
fit <- survfit(res.cox, newdata = sex_df)
ggsurvplot(fit, conf.int = TRUE, legend.labs=c("Sex=F", "Sex=M"),
           ggtheme = theme_minimal())

#old method
fit.old <- survfit(Surv(age) ~ sex, data = demographics)
ggsurvplot(fit.old)
survdiff(Surv(age) ~ sex, data = demographics)

#check residuals with new method with censoring
test.ph <- cox.zph(res.cox)
test.ph
ggcoxzph(test.ph)

tidy.res.cox <- tidy(res.cox)
tidy.res.cox
write.table(tidy.res.cox, file = "table1.txt", sep = ",")

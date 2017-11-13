library(plyr)
library(dplyr)
library(ggplot2)
library(stats4)
library(survival)
library(survminer)
library(ggpubr)
library(lattice)
library(car)
library(ggmap)

# https://www.r-bloggers.com/cox-proportional-hazards-model/
######################## Part 1 - Cox Models #######################

data("lung")
head(lung)
res.cox <- coxph(Surv(time, status) ~ sex, data = lung)
res.cox
summary(res.cox)


covariates <- c("age", "sex",  "ph.karno", "ph.ecog", "wt.loss")
univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(time, status)~', x)))
univ_models <- lapply( univ_formulas, function(x){coxph(x, data = lung)})
univ_results <- lapply(univ_models,
                       function(x){ 
                             x <- summary(x)
                             p.value<-signif(x$wald["pvalue"], digits=2)
                             wald.test<-signif(x$wald["test"], digits=2)
                             beta<-signif(x$coef[1], digits=2);#coeficient beta
                             HR <-signif(x$coef[2], digits=2);#exp(beta)
                             HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                             HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                             HR <- paste0(HR, " (", 
                                          HR.confint.lower, "-", HR.confint.upper, ")")
                             res<-c(beta, HR, wald.test, p.value)
                             names(res)<-c("beta", "HR (95% CI for HR)", "wald.test", 
                                           "p.value")
                             return(res)
                             #return(exp(cbind(coef(x),confint(x))))
                       })
res <- t(as.data.frame(univ_results, check.names = FALSE))
as.data.frame(res)                             
                             

#now do multivariate cox model
res.cox <- coxph(Surv(time, status) ~ age + sex + ph.ecog, data =  lung)
summary(res.cox)
ggsurvplot(survfit(res.cox), color = "#2E9FDF",
           ggtheme = theme_minimal())


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

########## Part 2 - check assumptions of Cox Model ############
#http://www.sthda.com/english/wiki/cox-model-assumptions


res.cox <- coxph(Surv(time, status) ~ age + sex + wt.loss, data =  lung)
res.cox
#tests proportional hazard assumption - not significant means they're ok for Cox
test.ph <- cox.zph(res.cox)
test.ph
#graph, non-proportional (bad for cox) if systematic departures from line
ggcoxzph(test.ph)
#check for outlier/influential data points
ggcoxdiagnostics(res.cox, type = "dfbeta",
                 linear.predictions = FALSE, ggtheme = theme_bw())
                             
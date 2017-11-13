# http://rpubs.com/alecri/258589

setwd("/Users/g-woloschak/Documents")

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
all.info <- filter(all.info,
                   expt != 11,
                   expt != 10,
                   radn == "C" | radn == "G",
                   !(expt == 14 & !(tmt == "0S" | tmt == "C0")),
                   cause_of_death != "Grahn mice, breeder",
                   expt != 2,
                   expt != 12,
                   expt !=3,
                   #!(expt == 9 & sex == "M"),
                   fractions != 300,
                   cause_of_death != "Removal to another experiment")
all.info$expt <- as.factor(all.info$expt)
all.info$fractions <- as.factor(all.info$fractions)
#add column for censoring: 0 = alive, 1 = dead
#censoring if wrong cause of death or no lethal disease listed
all.info$status <- 0
all.info$status[(all.info$cause_of_death == "Died" |
                           all.info$cause_of_death == "Sacrifice, moribund") &
                           all.info$lethal_sum == 1] <- 1
#all.info$status <- as.logical(all.info$status)
#check work with summary(all.info$status == 0) and see 1323 censored

#check numbers for censored animals
# dim(subset(all.info, cause_of_death == "Accidental death"))
# dim(subset(all.info, cause_of_death == "Escaped during irradiation"))
# dim(subset(all.info, cause_of_death == "Discard"))
# dim(subset(all.info, cause_of_death == "Improper irradiation"))
# dim(subset(all.info, cause_of_death == "Missing"))
# dim(subset(all.info, cause_of_death == "Removal to another experiment"))
# dim(subset(all.info, cause_of_death == "Sacrifice, programmed"))

#make sure each mouse only has one lethal macro --> good
lethal.columns <- grep("_L", names(all.info))
all.info$lethal_sum <- rowSums(all.info[,lethal.columns])

####################### survival analysis work for controls #################
km.survival <- function(stratified.by, df){
      f <- as.formula(paste("Surv(age) ~ ",
                          paste(stratified.by, collapse=" + ")))
  
      #had lots of weird errors... https://github.com/kassambara/survminer/issues/125
      fit <- do.call(survfit,
                     list(as.formula(paste("Surv(age) ~ ", 
                                           paste(stratified.by, collapse=" + "))),
                          data = df))
      print(survdiff(f, data = df))
      ggsurvplot(fit)
}

km.censored <- function(stratified.by, df){
      f <- as.formula(paste("Surv(age, status) ~",
                            paste(stratified.by, collapse=" + ")))
      
      #had lots of weird errors... https://github.com/kassambara/survminer/issues/125
      fit <- do.call(survfit,
                     list(as.formula(paste("Surv(age, status) ~ ", 
                                           paste(stratified.by, collapse=" + "))),
                          data = df))
      print(survdiff(f, data = df))
      ggsurvplot(fit)
}
#km.survival(stratified.by = "expt", df = subset(all.info, total_dose == 0))
#km.censored(stratified.by = "expt", df = subset(all.info, total_dose == 0))

cox.survival <- function(stratified.by, df){
      out <- NULL
      f <- as.formula(paste("Surv(age) ~ ",
                            paste(stratified.by, collapse=" + ")))
      res.cox <- coxph(f, data = df)
      print(summary(res.cox))
      expt_df <- data.frame(expt = as.factor(c(3, 4, 7, 8, 9, 13, 14)))
      fit <- survfit(res.cox, newdata = expt_df)
      plot1 <- ggsurvplot(fit, conf.int = FALSE, legend.labs=c("expt=3", "expt=4", "expt=7", "expt=8",
                                                      "expt=9", "expt=13", "expt=14"))
      # Test PH assumption - check residuals do not change with time
      # values should not be significant if assuption is OK
      test.ph <- cox.zph(res.cox)
      print(test.ph)
      plot2 <- ggcoxzph(test.ph)
      
      #Check for influential outliers
      #graphs should look symmetric around 0
      plot3 <- ggcoxdiagnostics(res.cox, type = "deviance",
                       linear.predictions = FALSE, ggtheme = theme_bw())
      
      out[[1]] <- plot1
      out[[2]] <- plot2
      out[[3]] <- plot3
      out
}
#graphs <- cox.survival(stratified.by = "expt", df = subset(all.info, total_dose == 0))
#graphs[[1]] for curve, graphs[[2]] for PH diagnostic, graphs[[3]] for outliers diagnostic

cox.censored <- function(stratified.by, df){
      out <- NULL
      f <- as.formula(paste("Surv(age, status) ~ ",
                            paste(stratified.by, collapse=" + ")))
      print(f)
      res.cox <- coxph(f, 
                       data = df)
    
      print(summary(res.cox))
      
      newd <- expand.grid(sex = c("M", "F"), expt = levels(all.info$expt))
      newd$id <- 1:16
      fit <- survfit(res.cox)
      fortify(survfit(res.cox, newdata = newd)) %>%
            gather(strata, surv, surv.1:surv.16) %>%
            mutate(id = gsub("surv.","", strata)) %>%
            merge(newd, by = "id") %>%
            ggplot(aes(x = time, y = surv, col = expt, linetype = sex)) +
            geom_step() + 
            #facet_wrap( ~ sex, nrow = 2) +
            labs(x = "Time (days)", y = "Survival probability") + theme_classic()
      
      # Test PH assumption - check residuals do not change with time
      # values should not be significant if assuption is OK
      test.ph <- cox.zph(res.cox)
      print(test.ph)
      plot2 <- ggcoxzph(test.ph)
      
      #Check for influential outliers
      #graphs should look symmetric around 0
      plot3 <- ggcoxdiagnostics(res.cox, type = "deviance",
                                linear.predictions = FALSE, ggtheme = theme_bw())
      
      out[[1]] <- plot1
      out[[2]] <- plot2
      out[[3]] <- plot3
      out
}

all.info$expt <- droplevels(all.info$expt)
all.info.no8 <- subset(all.info, expt !=8)
all.info.no8$fractions <- droplevels(all.info.no8$fractions)
#graphs <- cox.censored(stratified.by = "expt", df = subset(all.info, total_dose == 0))
#all.info$expt <- relevel(all.info$expt, ref = "14")
#all.info.no8$fractions <- relevel(all.info.no8$fractions, ref = "120")
res.cox <- coxph(Surv(age, status) ~ sex + expt, 
                 data = subset(all.info, total_dose == 0))
summary(res.cox)
newd <- expand.grid(sex = c("M", "F"), expt = c("4", "7", "8", "9", "13", "14"))
newd <- expand.grid(sex = c("M", "F"), fractions = levels(all.info.no8$fractions))
newd$id <- 1:10
# fit <- survfit(res.cox)
# res <- fortify(survfit(res.cox, newdata = newd))
# what <- fortify(survfit(res.cox, newdata = newd)) %>%
#       gather(strata, surv, surv.1:surv.16)
# huh <- fortify(survfit(res.cox, newdata = newd)) %>%
#       gather(strata, surv, surv.1:surv.16) %>%
#       mutate(id = gsub("surv.","", strata)) %>%
#       merge(newd, by = "id")
sex_names <- c(
      `M` = "Male",
      `F` = "Female")

fortify(survfit(res.cox, newdata = newd)) %>%
      gather(strata, surv, surv.1:surv.10) %>%
      mutate(id = gsub("surv.","", strata)) %>%
      merge(newd, by = "id") %>%
      ggplot(aes(x = time, y = surv, col = fractions)) +
      geom_step(size = .75) + 
      facet_grid(. ~ sex, labeller = as_labeller(sex_names)) +
      labs(x = "Time (days)", y = "Survival probability") + 
      theme_bw()+
      theme(strip.background =element_rect(fill="white"))+
      theme(strip.text = element_text(colour = 'black', face = "bold", size = 12))+
      theme(text=element_text(size=12, face = "bold"))

# Test PH assumption - check residuals do not change with time
# values should not be significant if assuption is OK
test.ph <- cox.zph(res.cox)
print(test.ph)
plot2 <- ggcoxzph(test.ph,font.main = 11,
                  ggtheme = theme_bw()+
                        theme(text=element_text(size=12, face = "bold")))

# 
# plot2 <- ggcoxzph(test.ph, font.main = 11, 
#                   text=element_text(size=11, face = "bold"),
#                   font.tickslab = 7, point.col = "red", 
#                   point.size = .5)
#Check for influential outliers
#graphs should look symmetric around 0
plot3 <- ggcoxdiagnostics(res.cox, type = "deviance",
                          linear.predictions = FALSE, 
                          ggtheme = theme_bw() + 
                                theme(strip.background =element_rect(fill="white"))+
                                theme(strip.text = element_text(colour = 'black', face = "bold", size = 12))+
                                theme(text=element_text(size=12, face = "bold")))
                                
                                
                                theme(text=element_text(size=12, face = "bold"),
                                          plot.title = element_text(size = 16)))

all.info <- subset(all.info,
                   expt != 12 & fractions!= 300)

investigate1 <- subset(all.info, expt == 3) %>%
      group_by(sex, radn, fractions, status, total_dose) %>%
      summarise(n=n(),
                mean_age = mean(age))
investigate2 <- subset(all.info, expt == 4 & total_dose == 0) %>%
      group_by(expt, sex, radn, fractions, status, total_dose) %>%
      summarise(n=n(),
                mean_age = mean(age))

df <- subset(all.info, expt ==9 & sex == "F")

what <- df %>%
      group_by(sex, radn, status, fractions) %>%
      summarise(n=n())

all.info %>%
      filter(expt == 9, sex == "F") %>%
      group_by(sex, radn, status, fractions) %>%
      summarise(n=n(),
                mean_age = mean(age))

expt_df <- with(demographics,
                data.frame(expt = c(3, 4, 7, 8, 9, 12, 13, 14)))
sex_df <- with(demographics,
               data.frame(sex = c("M", "F")))
sex_df
expt_df
#change to have or not have censoring
res.cox <- coxph(Surv(age, status) ~ factor(expt), data = subset(all.info, total_dose == 0))
summary(res.cox)
ggsurvplot(survfit(res.cox), color = "#2E9FDF",
           ggtheme = theme_minimal())
fit <- survfit(res.cox, newdata = expt_df)
ggsurvplot(fit, conf.int = FALSE, legend.labs=c("expt=3", "expt=4", "expt=7", "expt=8",
                                                "expt=9", "expt=12", "expt=13", "expt=14"))

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

cox.censored <- function(stratified.by, df){
      out <- NULL
      f <- as.formula(paste("Surv(age, status) ~ ",
                            paste(stratified.by, collapse=" + ")))
      print(f)
      res.cox <- coxph(f, 
                       data = df)
      
      print(summary(res.cox))
      # sex_df <- data.frame(sex = c("M", "F"))
      # expt_df <- data.frame(expt = as.factor(c(3, 4, 7, 8, 9, 13, 14)))
      # fractions_df <- data.frame(fractions = as.factor(c(0, 1, 24, 60, 120)))
      # fit <- survfit(res.cox, newdata = sex_df)
      # plot1 <- ggsurvplot(fit, conf.int = FALSE, legend.labs=c("Male", "Female"))
      
      # fit <- survfit(res.cox, newdata = expt_df)
      # plot1 <- ggsurvplot(fit, conf.int = FALSE, legend.labs=c("expt=3", "expt=4", "expt=7", "expt=8",
      #                                                         "expt=9", "expt=13", "expt=14"))
      # fit <- survfit(res.cox, newdata = fractions_df)
      # plot1 <- ggsurvplot(fit, conf.int = FALSE, legend.labs=c("0", "1", "24", "60",
      #                                                          "120"))
      # 
      
      
      
      
      
      # Test PH assumption - check residuals do not change with time
      # values should not be significant if assuption is OK
      test.ph <- cox.zph(res.cox)
      print(test.ph)
      plot2 <- ggcoxzph(test.ph)
      
      #Check for influential outliers
      #graphs should look symmetric around 0
      plot3 <- ggcoxdiagnostics(res.cox, type = "deviance",
                                linear.predictions = FALSE, ggtheme = theme_bw())
      
      #out[[1]] <- plot1
      out[[2]] <- plot2
      out[[3]] <- plot3
      out
}

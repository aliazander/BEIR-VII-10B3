library(plyr)
library(dplyr)
library(ggplot2)
library(RColorBrewer)

setwd("/Users/g-woloschak/Documents")
all.data <- read.csv("Filtered_data.csv")
head(all.data)
all.data$fractions <- as.factor(all.data$fractions)

#create new data frame that has all types of lethal causes of death as rows
#first column has total number of mice with that cancer
#second calumn is if it's a tumor or not
#third column is whether it's lymphona or not.
#fourth is whether it's CDU or not
#fifth column is total counts of tumors
#sixth column is total counts of non-tumors
#seventh column is total counts of lymphoma
#eight column is total counts of CDU
death.counts <- data.frame(occurrences = colSums(select(all.data, NTYG_L:TEP_L)))
death.counts$tumor <- substr(rownames(death.counts), 1, 1) =="T"
lymphoma <- c("NTYG_L","NTYL_L", "TTYG_L", "TTYL_L")
death.counts$lymphoma <- rownames(death.counts) %in% lymphoma
death.counts$cdu <- rownames(death.counts) %in% cdu
cdu <- "CDU_L"
death.counts$total.tumors <- death.counts$occurrences * death.counts$tumor * !(death.counts$lymphoma)
death.counts$total.non.tumors <- ((death.counts$occurrences * !(death.counts$tumor)) * !(death.counts$lymphoma)) * !(death.counts$cdu)
death.counts$total.lymphoma <- death.counts$occurrences * death.counts$lymphoma
death.counts$total.cdu <- death.counts$occurrences*death.counts$cdu
#death_counts <- ordered(death_counts, "occurrences", TRUE)

#calculate percentage of deaths from each cause of death
percent.deaths <- data.frame(total = colSums(select(death.counts, total.tumors:total.cdu)))
percent.deaths$percentage <- percent.deaths$total/sum(percent.deaths[,1])
percent.deaths$Type <- c("Tumors", "Non-tumors", "Lymphoma", "CDU")
percent.deaths

bp<- ggplot(percent.deaths, aes(x="", y=percentage, fill=Type))+
      geom_bar(width = 1, stat = "identity")
bp
pie <- bp + coord_polar("y", start=0)
pie + scale_fill_brewer(palette = "Purples") + 
      theme(axis.text.x=element_blank())+
      theme_void() + 
      theme(legend.title=element_text(size=12, face = "bold") , 
            legend.text=element_text(size=12, face = "bold"))+
      geom_text(aes(y = percentage/2 + c(0, cumsum(percentage)[-length(percentage)]), 
                    label = percent(percentage), x=1.25), size=4, fontface = "bold")
ggsave(file="Percent-COD.png", height = 4, width = 5, units = "in")

####################### END PIE CHART OF COD  ##########################################


#get list of tumor names
tumor.names <- substr(rownames(death.counts), 1, 1) == "T"
tumor.names <- rownames(death.counts)[tumor.names]
tumor.names <- tumor.names[!(tumor.names %in% lymphoma)]

#get list of nontumor names
nontumor.names <- rownames(death.counts)
nontumor.names <- nontumor.names[!(nontumor.names %in% lymphoma)]
nontumor.names <- nontumor.names[!(nontumor.names %in% tumor.names)]
nontumor.names <- nontumor.names[!(nontumor.names %in% cdu)]

#create class column to describe cause of death
all.data$class <- 0
#for each row, find which column has a 1 to determine COD
location <- apply(all.data[,21:157], 1, function(r) which(r == "1"))
#go through list and for each row, classify it as tumor, nontumor, etc. in new column
for(i in 1:length(location)){     
      column <- location[[i]] + 20
      #https://stackoverflow.com/questions/6451152/how-to-catch-integer0
      if(is.integer(location[[i]]) && length(location[[i]]) == 0L){
            all.data[i, 160] <- "NA"
      }
      else if(colnames(all.data)[column] %in% tumor.names){
            all.data[i, 160] <- "Tumor"
      }
      else if(colnames(all.data)[column] %in% nontumor.names){
            all.data[i, 160] <- "Non-tumor"
      }
      else if(colnames(all.data)[column] %in% lymphoma){
            all.data[i, 160] <- "Lymphoma"
      }
      else {all.data[i, 160] <- "CDU"}
}

my_palette = c(brewer.pal(4, "Purples"))
my_title <- "Age at death for different causes of death and irradiation treatments"
p <- ggplot(subset(all.data, class != "NA"), aes(x = class, y = age, fill = radn)) +
      geom_boxplot() +
      theme_bw() +
      ggtitle("Age at death for different causes of death and irradiation treatments") + 
      #ggtitle(expression(atop("A long string of text for the purpose", 
                              #paste("of illustrating my point")))) +
      #ggtitle(paste(strwrap(my_title, width = 20), collapse = "/n")) +
      scale_y_continuous(name = "Age (days)") +
      theme(text = element_text(size = 11, face = "bold"),
            axis.title.x=element_blank(),
            plot.title = element_text(hjust = 0.5)) + 
      scale_fill_manual(labels = c("Control", "Gamma"), 
                         values = c(my_palette[2], my_palette[4])) +
      labs(fill = "Treatment")
p

ggsave(file="Age-vs-cod-treatment.png", height = 4, width = 7, units = "in")
################## END BAR GRAPH OF AGE VS. COD AND TMT  ###############################


jm8 <- subset(all.data, expt ==8)
jm8$dose_rate <- as.factor(jm8$dose_rate)
res.cox <- coxph(Surv(age, status) ~ sex + dose_rate, 
                 data = jm8)
summary(res.cox)
newd <- expand.grid(sex = c("M", "F"), 
                    #dose.group = levels(all.data.no8$dose.group),
                    dose_rate = levels(jm8$dose_rate))
newd$id <- 1:8

fortify(survfit(res.cox, newdata = newd)) %>%
      gather(strata, surv, surv.1:surv.8) %>%
      mutate(id = gsub("surv.","", strata)) %>%
      merge(newd, by = "id") %>%
      ggplot(aes(x = time, y = surv, col = dose_rate)) +
      geom_step(size = .75) + 
      facet_grid(. ~ sex, labeller = as_labeller(sex_names)) +
      labs(x = "Time (days)", y = "Survival probability", colour = "Dose Rate (mGy/h)") + 
      theme_bw()+
      theme(strip.background =element_rect(fill="white"))+
      theme(strip.text = element_text(colour = 'black', face = "bold", size = 12))+
      theme(text=element_text(size=12, face = "bold"))+
      scale_color_hue(labels = c("0", "88.8", "222", "408"))

      
ggsave(file="Survival-JM8-Dose-Rate.png", height = 4, width = 8, units = "in")

test.ph <- cox.zph(res.cox)
print(test.ph)

plot3 <- ggcoxdiagnostics(res.cox, type = "deviance",
                          linear.predictions = FALSE, 
                          point.size = .5,
                          ggtheme = theme_bw() + 
                                theme(strip.background =element_rect(fill="white"))+
                                theme(strip.text = element_text(colour = 'black', face = "bold", size = 12))+
                                theme(text=element_text(size=12, face = "bold")))
ggsave(file="Residuals-dose-rate-jm8.png", height = 3, width = 4, units = "in")
################## END COX FOR JM8 - DOSE RATE  ###############################

all.data$dose.group <- 0
for(i in 1:nrow(all.data)){     
      if(all.data$total_dose[i] <= 0){
            all.data$dose.group[i] <- 0
      }
      else if(all.data$total_dose[i] > 0 &
              all.data$total_dose[i] <= 150){
            all.data$dose.group[i] <- 1.5
      }
      else if(all.data$total_dose[i] > 150 &
              all.data$total_dose[i] <= 400){
            all.data$dose.group[i] <- 4
      }
      else if(all.data$total_dose[i] > 400 &
              all.data$total_dose[i] <= 760){
            all.data$dose.group[i] <- 7.6
      }
      else if(all.data$total_dose[i] > 760 &
              all.data$total_dose[i] <= 1640){
            all.data$dose.group[i] <- 16.4
      }
      else all.data$dose.group[i] <- 49.01
}

all.data[9000:9007, 15:161]
all.data$dose.group <- as.factor(all.data$dose.group)

all.data.no8$fractions <- relevel(all.data.no8$fractions, ref = "1")
res.cox <- coxph(Surv(age, status) ~ sex + fractions * dose.group, 
                 data = all.data.no8)
summary(res.cox)
test.ph <- cox.zph(res.cox)
print(test.ph)
ggplot(subset(all.data.no8, status == 1 & class == "Tumor"),
      aes(factor(dose.group), age, fill = factor(fractions)))+
      geom_boxplot() + 
      theme_bw() +
      ggtitle("Age vs. Dose")+
      labs(x="Dose Group (cGy)",y="Age", fill = "Fractions") +
      theme(text = element_text(size = 11, face = "bold"),
            plot.title = element_text(hjust = 0.5)) +
      scale_fill_brewer(palette = "Purples") +

ggsave(file="Age-vs-dose-fractions-boxplot-no8.png", height = 4, width = 8, units = "in")

  

############### END COX AND BAR GRAPH FOR AGE VS. DOSE GROUP  ########################








all.data.no8 <- subset(all.data, expt !=8)
all.data.no8$fractions <- droplevels(all.data.no8$fractions)

res.cox <- coxph(Surv(age, status) ~ sex + fractions, 
                 data = all.data.no8)
summary(res.cox)
newd <- expand.grid(sex = c("M", "F"), 
                    #dose.group = levels(all.data.no8$dose.group),
                    fractions = levels(all.data.no8$fractions))
newd$id <- 1:10

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
ggsave(file="Survival-sex-fractions-dosegroup.png", height = 4, width = 8, units = "in")

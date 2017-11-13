library(plyr)
library(dplyr)
library(ggplot2)
library(stats4)

setwd("/Users/g-woloschak/Documents/")

##downlaod data for plot 10B-3
BEIR_VII_10B_3 <- read.csv('BEIR_VII_10B-3_data.csv')
BEIR_VII_10B_3$sem <- BEIR_VII_10B_3$sd / sqrt(BEIR_VII_10B_3$n)

#names(BEIR_VII_10B_3)[names(BEIR_VII_10B_3)=="sd"] <- "sem"
#BEIR_VII_10B_3$sd <- BEIR_VII_10B_3$sem*sqrt(BEIR_VII_10B_3$n)
##Seperate acute vs. chornic data with restrictions Ben discovered - 
## RFM mice and female only
acute_plot_10B3 <- BEIR_VII_10B_3[BEIR_VII_10B_3[,6] == "F" & 
                                        BEIR_VII_10B_3$type == "A" &  
                                        BEIR_VII_10B_3[,5] == "RFM",]
chronic_plot_10B3 <- BEIR_VII_10B_3[BEIR_VII_10B_3[,6] == "F" &
                                          BEIR_VII_10B_3$type == "C" &
                                          BEIR_VII_10B_3[,5] == "RFM",]

##replaces age with the average age when they have the same dose
##acute_plot_10B3$age <- ave(acute_plot_10B3$age, acute_plot_10B3$dose)

##combines sd when they have the same dose (sqrt(sd1^2 + sd2^2 + ...))

#returns list with the duplicate dose values (0, .5, 2 for this data set)
not_unique <- acute_plot_10B3$dose[duplicated(acute_plot_10B3$dose)]

#go through each non-unique dose
for (i in 1:length(not_unique)){
      #find row number for dose values that are duplicates
     row_numbers <- which(acute_plot_10B3$dose == not_unique[i])
     #keep the first row for averaging, all others will be deleted
     row_to_keep <- row_numbers[1]
     #keep track of total age and squared sd for all values at that dose
     total_age <- 0
     sd_sum_squared <- 0
     n <- 0
     #go through each row with the same dose
     for(k in 1:length(row_numbers)){
           row_to_use <- row_numbers[k]
           total_age <- total_age + acute_plot_10B3[row_to_use,7]
           sd_sum_squared <- sd_sum_squared + (acute_plot_10B3[row_to_use,8])^2
           n <- n + acute_plot_10B3[row_to_use, 2]
     }
     #replace current age and stdev of the top row of duplicates
     acute_plot_10B3[row_to_keep, 7] <- total_age/length(row_numbers)
     acute_plot_10B3[row_to_keep, 8] <- sqrt(sd_sum_squared)
     acute_plot_10B3[row_to_keep, 2] <- n
}
acute_plot_10B3$sem <- acute_plot_10B3$sd / sqrt(acute_plot_10B3$n)

#only use the new, non-duplicate values that were combined
acute_plot_10B3 <- acute_plot_10B3[!duplicated(acute_plot_10B3$dose),]

#calculate SEM from SD and N


##combine acute and chronic back together for plot
data_plot_10B3 <- rbind(acute_plot_10B3, chronic_plot_10B3)

#plot all of the data points 
p_10B3 <- ggplot(data_plot_10B3, aes(x=dose, y = 1/age)) + 
      geom_text(family = "Helvetica", size = 3.8, aes(label = type, group = type)) +
      scale_y_continuous(breaks = c(.0016, .0018, .0020, .0022, .0024)) +
      xlab("Dose (Gy)") +
      ylab("Reciprocal of Mean Life Length in days")

p_10B3

#filter data for linear regression - only doses below 1.5Gy

model_data <- data_plot_10B3[data_plot_10B3$dose <= 1.5,]

#create quadratic model - need the same alpha and constant for acute and linear
#models, but do not need dose squared term for chronic (divide by infinity, so equal
#to zero). Multiply dose squared term by 0 for chronic (!true = 0)
#Ben figured out - weighted by 1/(standard deviation)^2

model <- lm(I(1/age) ~ dose + I(dose^2 * (!model_data$chronic)),
            data = model_data,
            weights = 1/sem^2)
#extract constant, alpha, and theta for equations for lines
coefficients <- model$coefficients
predicted_xvalues <- as.data.frame(seq(from = 0, to = 1.5, by = .01))

acute_regression_data <- coefficients[1] +
      coefficients[2] * predicted_xvalues +
      coefficients[3] * predicted_xvalues^2

data_acute_regression_plot <- as.data.frame(cbind(predicted_xvalues, acute_regression_data))
names(data_acute_regression_plot) <- c("dose", "fit")

#chronic --> constant + alpha*D
chronic_regression_data <- coefficients[1] +
      coefficients[2] * predicted_xvalues
data_chronic_regression_plot <- as.data.frame(cbind(predicted_xvalues, chronic_regression_data))
names(data_chronic_regression_plot) <- c("dose", "fit")


#add fit lines to plot
p_10B3 <- p_10B3 + geom_path(data=data_chronic_regression_plot, aes(dose, fit)) +
      geom_path(data = data_acute_regression_plot, aes(dose,fit))

#ta-da!
p_10B3
summary(model)
DDREF <- coefficients[3]/coefficients[2] + 1
DDREF

############## 10B4, part 1  ##################################

#test to see if this way of setting up lm works by using theta that I know work and what the
# other coefficients should come out to be
o<-1.2638
model2 <- lm(I(1/age) ~ I(dose+o*I(dose^2*!model_data$chronic)), data = model_data, weights = 1/sd^2)

summary(model2)
logLik(model2)
##it worked!

#create for loop to get likelihood values for range of theta values
model_data <- BEIR_VII_10B_3[BEIR_VII_10B_3[,6] == "F" &
                                   BEIR_VII_10B_3$dose <= 1.5 &
                                   BEIR_VII_10B_3[,5] == "RFM",]

likelihood_values <- list()
o_range <- seq(-2, 6, .02)
for(i in 1:length(o_range)){
      o<- o_range[i]
      model3 <- lm(I(1/age) ~ I(dose+o*dose^2*(!model_data$chronic)), 
                   #weights = 1/sem^2,
                   data = model_data)
      likelihood_values[i] <- logLik(model3)
}
#take list of likelihood values and put in a dataframe
df <- data.frame(matrix(unlist(likelihood_values), nrow=length(o_range), byrow=T))
#add theta values to dataframe
df$theta_values <- o_range
#label columns
colnames(df) <- c("likelihoods", "thetas")
df$likelihoods <- exp(df$likelihoods)

#graph theta vs likelihood
p_10B4 <- ggplot(df, aes(x=thetas, y = likelihoods)) + 
      geom_point() +
      xlab("Possible Values of theta, in Gy-1") +
      ylab("Likelihood Value")
p_10B4

#get row number of max likelihood value
row_number <- which(df$likelihoods == max(df$likelihoods, na.rm=TRUE))
#use row number to get best theta value
best_theta <- df[row_number, 2]
best_theta


###################play with profile likelihoods ######################

y <- rnorm(10000, mean = 4, sd = 1) 
x <- seq(-3, 3, length = 10000)
y <- y + rnorm(10000, mean = 1, sd =.5)*x + rnorm(10000, mean=3, sd=1)*x^2
df2 <- data.frame(x, y)
model4 <- lm(y ~ x + I(x^2))
new_data <- data.frame(y=predict(model4), x=df2$x)
ggplot(df2, aes(x=x, y=y)) + geom_point() + geom_path(data = new_data, aes(x=x, y=y))
ggplot(df2, aes(y)) + geom_histogram()

likelihood_values <- list()
mu_range <- seq(1, 80, .05)
for(i in 1:length(mu_range)){
      mu<- mu_range[i]
      model5 <- lm(y ~ I(x+mu*x^2), data = df2)
      likelihood_values[i] <- logLik(model5)
}
#take list of likelihood values and put in a dataframe
df3 <- data.frame(matrix(unlist(likelihood_values), nrow=length(mu_range), byrow=T))
#add theta values to dataframe
df3$mu_values <- mu_range
#label columns
colnames(df3) <- c("likelihoods", "mus")


#http://stats.stackexchange.com/questions/56724/computation-of-likelihood-when-n-is-very-large-so-likelihood-gets-very-small
#divide to make numbers reasonable to calculations - R gives 0 otherwise
df3$likelihoods <- df3$likelihoods/10000
df3$likelihoods <- exp(df3$likelihoods)

row_number <- which(df3$likelihoods == max(df3$likelihoods, na.rm=TRUE))
best_theta <- df3[row_number, 2]

#graph theta vs likelihood
p_practice <- ggplot(df3, aes(x=mus, y = likelihoods)) + 
      geom_point() +
      xlab("Possible Values of theta, in Gy-1") +
      ylab("Likelihood Value")
p_practice

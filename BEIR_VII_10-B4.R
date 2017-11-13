library(plyr)
library(dplyr)
library(ggplot2)
library(stats4)
library(gridExtra)


setwd("/Users/g-woloschak/Documents/")

BEIR_VII_animal_data <- read.csv('BEIR_VII_animal_data.csv')

data_to_plot <- BEIR_VII_animal_data[BEIR_VII_animal_data$protraction == "Acute" & 
                                           !is.na(BEIR_VII_animal_data$risk) &
                                           BEIR_VII_animal_data$dose <=2,]

#create for loop to get likelihood values for range of theta values
likelihood_values <- list()
o_range <- seq(-2, 6, .02)
for(i in 1:length(o_range)){
      o<- o_range[i]
      model_one_theta <- lm(risk ~ factor(facet) * I(dose+o*dose^2), 
                            data = data_to_plot, 
                            weights = (1/error^2))
      likelihood_values[i] <- logLik(model_one_theta)
}
#take list of likelihood values and put in a dataframe
df_10b2 <- data.frame(matrix(unlist(likelihood_values), nrow=length(o_range), byrow=T))
#add theta values to dataframe
df_10b2$theta_values <- o_range
#label columns
colnames(df_10b2) <- c("likelihoods", "thetas")

#find 95% CI
bestloglik <- max(df_10b2$likelihoods, na.rm = TRUE)
minloglik <- bestloglik - qchisq(.95, 1)/2
min(df_10b2[df_10b2$likelihoods >= minloglik,][,2])
max(df_10b2[df_10b2$likelihoods >= minloglik,][,2])


#convert from log likelihood to likelihood
df_10b2$likelihoods <- exp(df_10b2$likelihoods)

#get row number of max likelihood value
row_number <- which(df_10b2$likelihoods == max(df_10b2$likelihoods, na.rm=TRUE))
#use row number to get best theta value
best_theta <- df_10b2[row_number, 2]
#plot profile likelihood
p_10B2 <- ggplot(df_10b2, aes(x=thetas, y = likelihoods)) + geom_point()

##downlaod data for plot 10B-3######################################################
BEIR_VII_10B_3 <- read.csv('BEIR_VII_10B-3_data.csv')
BEIR_VII_10B_3$sem <- BEIR_VII_10B_3$sd / sqrt(BEIR_VII_10B_3$n)

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
df_10B3 <- data.frame(matrix(unlist(likelihood_values), nrow=length(o_range), byrow=T))
#add theta values to dataframe
df_10B3$theta_values <- o_range
#label columns
colnames(df_10B3) <- c("likelihoods", "thetas")
df_10B3$likelihoods <- exp(df_10B3$likelihoods)

#graph theta vs likelihood
p_10B3 <- ggplot(df_10B3, aes(x=thetas, y = likelihoods)) + 
      geom_point() +
      xlab("Possible Values of theta, in Gy-1") +
      ylab("Likelihood Value")
p_10B3

#get row number of max likelihood value
row_number <- which(df$likelihoods == max(df_10B3$likelihoods, na.rm=TRUE))
#use row number to get best theta value
best_theta <- df_10B3[row_number, 2]
best_theta

######################## 10-B4 combo ##############################

df_10B4 <- data.frame(theta=df_10b2$thetas, 
                      like_10b2 = df_10b2$likelihoods, 
                      like_10b3 = df_10B3$likelihoods)

#https://en.wikipedia.org/wiki/Normalizing_constant
#http://stats.stackexchange.com/questions/31238/what-is-the-reason-that-a-likelihood-function-is-not-a-pdf#comment60635_31248
#http://cmd.inp.nsk.su/old/cmd2/manuals/cernlib/minuit/node46.html

df_10B4$like_10b2 <- df_10B4$like_10b2/sum(df_10B4$like_10b2)
df_10B4$like_10b3 <- df_10B4$like_10b3/sum(df_10B4$like_10b3)
df_10B4$like_10b4 <- (df_10B4$like_10b2+df_10B4$like_10b3)/2

df_10B4[,2:4] <- df_10B4[,2:4]/max(df_10B4$like_10b2)

p_10B4 <- ggplot(data = df_10B4, aes(x=theta)) +
      geom_line(aes(y=like_10b2), colour="red") +  # first layer
      geom_line(aes(y=like_10b3), colour="green") +  # second layer
      geom_line(aes(y=like_10b4), colour="blue") + # third layer
      xlab("Possible Values of Theta (1/Gy)") +
      ylab("Likelihood Value")

p_10B4

row_number <- which(df_10B4$like_10b4 == max(df_10B4$like_10b4, na.rm=TRUE))
best_theta <- df_10B4[row_number, 1]
best_theta

#find 95% CI
loglike_CI <- data.frame(theta = df_10B4$theta, loglike = log(df_10B4$like_10b4))
bestloglik <- max(loglike_CI$loglike, na.rm = TRUE)
minloglik <- bestloglik - qchisq(.95, 1)/2
min(loglike_CI[loglike_CI$loglike >= minloglik,][,1])
max(loglike_CI[loglike_CI$loglike >= minloglik,][,1])



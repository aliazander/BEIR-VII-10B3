
library(plyr)
library(dplyr)
library(ggplot2)
library(stats4)
library(gridExtra)

#after updating R I had an issue getting all 11 graphs, and 
#a warning that said "reached elapsed time limit" so I ran
#setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE)
#https://stat.ethz.ch/R-manual/R-devel/library/base/html/setTimeLimit.html
#problem fixed!

#setwd("/Users/g-woloschak/Documents/janus/data")
#BEIR_VII_animal_data <- read.csv('Edwards_1992.csv')
setwd("/Users/g-woloschak/Documents")
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
df <- data.frame(matrix(unlist(likelihood_values), nrow=length(o_range), byrow=T))
#add theta values to dataframe
df$theta_values <- o_range
#label columns
colnames(df) <- c("likelihoods", "thetas")

#find 95% CI
bestloglik <- max(df$likelihoods, na.rm = TRUE)
minloglik <- bestloglik - qchisq(.95, 1)/2
min(df[df$likelihoods >= minloglik,][,2])
max(df[df$likelihoods >= minloglik,][,2])


#convert from log likelihood to likelihood
df$likelihoods <- exp(df$likelihoods)

#get row number of max likelihood value
row_number <- which(df$likelihoods == max(df$likelihoods, na.rm=TRUE))
#use row number to get best theta value
best_theta <- df[row_number, 2]
#plot profile likelihood
ggplot(df, aes(x=thetas, y = likelihoods)) + geom_point()

#use best theta for model that keeps one theta and figure out coeffients
#for each factor
model_best_theta <- lm(risk ~ factor(facet) * I(dose + best_theta*dose^2),
                       data = data_to_plot,
                       weights = 1/error^2)
#list of 11 intercepts then list of 11 coefficients of dose part
intercept_best_theta <- (model_best_theta$coefficients[1:11])
intercept_best_theta <- c(intercept_best_theta[1], intercept_best_theta[2:11] + intercept_best_theta[1])
coefficients_best_theta <- model_best_theta$coefficients[12:22]
coefficients_best_theta <- c(coefficients_best_theta[1], coefficients_best_theta[2:11] + coefficients_best_theta[1])

facets <- unique(data_to_plot$facet)
out <- NULL
for (i in 1:length(facets)){
      facet_i <- data_to_plot[data_to_plot$facet == i,]
 
      model <- lm(risk ~ dose + I(dose^2), 
                  data = facet_i,
                  weights = 1/error^2)
      coefficients <- model$coefficients
      predicted_xvalues <- as.data.frame(seq(from = 0, to = 2.0, by = .01))
      #acute = constant + alpha*(D + theta*Dose^2)
      acute_regression_data <- coefficients[1] + 
            coefficients[2] * predicted_xvalues + 
            coefficients[3] * predicted_xvalues^2
      data_acute_regression_plot <- as.data.frame(cbind(predicted_xvalues, acute_regression_data))
      names(data_acute_regression_plot) <- c("dose", "fit")
      ylim_i <- c(min(facet_i$risk)-.1, max(facet_i$risk)+.1)
      pi <- ggplot(facet_i, 
                   aes(x=dose, y = risk))
      pi <- pi + geom_point() + coord_cartesian(xlim = c(-.02, 2.02), ylim = ylim_i) 
      limits <- aes(ymax = risk + 2*error, ymin=risk - 2*error)
      pi <- pi + geom_errorbar(limits) + geom_point()
      pi <- pi + geom_path(data=data_acute_regression_plot, aes(dose, fit)) 
      
      
      #acute = constant + alpha*(Dose + best_theta*Dose^2)
      best_theta_data <- intercept_best_theta[i] + 
            coefficients_best_theta[i] * 
            (predicted_xvalues + best_theta * predicted_xvalues^2)
      best_theta_regression_plot <- as.data.frame(cbind(predicted_xvalues, best_theta_data))
      names(best_theta_regression_plot) <- c("dose", "fit")
      pi <- pi + geom_path(data=best_theta_regression_plot, aes(dose, fit), linetype = 3) 
      
      out[[i]] <- pi
}

do.call("grid.arrange", c(out, ncol=3))

out[11]

###### things I tried and learned from, but they are no longer needed ######
#also, test area for adding new ideas

facet9 <- BEIR_VII_animal_data[BEIR_VII_animal_data$facet == 4,]
facet9_acute <- facet9[facet9$protraction == "Acute" & facet9$dose <= 2,]

p1 <- ggplot(facet9_acute, 
             aes(x=dose, y = risk))
p1 <- p1 + geom_point() + coord_cartesian(xlim = c(-.02, 2.02), ylim = c(19, 31)) 
limits <- aes(ymax = risk + 2*error, ymin=risk - 2*error)
p1 <- p1 + geom_errorbar(limits) + geom_point()

model <- lm(risk ~ dose + I(dose^2), 
            data = facet9_acute,
            weights = 1/error^2)
coefficients <- model$coefficients
predicted_xvalues <- as.data.frame(seq(from = 0, to = 2.0, by = .01))
#acute = constant + alpha*(D + theta*Dose^2)
acute_regression_data <- coefficients[1] + 
      coefficients[2] * predicted_xvalues + 
      coefficients[3] * predicted_xvalues^2
data_acute_regression_plot <- as.data.frame(cbind(predicted_xvalues, acute_regression_data))
names(data_acute_regression_plot) <- c("dose", "fit")

best_theta <- .42
model_best_theta <- lm(risk ~ I(dose + best_theta*dose^2),
                       data = facet9_acute,
                       weights = 1/error^2)
coefficients_best_theta <- model_best_theta$coefficients

#acute = constant + alpha*(Dose + best_theta*Dose^2)
best_theta_data <- coefficients_best_theta[1] + 
      coefficients_best_theta[2] * (predicted_xvalues + best_theta * predicted_xvalues^2)
best_theta_regression_plot <- as.data.frame(cbind(predicted_xvalues, best_theta_data))
names(best_theta_regression_plot) <- c("dose", "fit")

p1 <- p1 + geom_path(data=data_acute_regression_plot, aes(dose, fit)) 
p1 <- p1 + geom_path(data=best_theta_regression_plot, aes(dose, fit), linetype = 3) 
p1

facet4 <- BEIR_VII_animal_data[BEIR_VII_animal_data$facet == 4 & 
                                     BEIR_VII_animal_data$protraction == "Acute" & 
                                     BEIR_VII_animal_data$dose <= 2 ,]
p4 <- ggplot(facet4, aes(x=dose, y = risk))
limits <- aes(ymax = risk + 2*error, ymin=risk - 2*error)
p4 <- p4 + geom_errorbar(limits) + geom_point()
p4 <- p4 + geom_errorbar(limits) + geom_point() + stat_smooth(method = glm, aes(weight = (1/facet4$error^2)), formula = y ~ x + I(x^2)) 
p4 <- p4 + coord_cartesian(xlim = c(0, 2), ylim = c(20, 31))
p4

p10_b2 <- ggplot(data_to_plot, aes(x=dose, y=risk)) + geom_point()
p10_b2 <- p10_b2 + facet_wrap( ~ facet, ncol=3, scales="free_y")
p10_b2 <- p10_b2 + geom_pointrange(aes(ymin=risk - error*2, ymax=risk + error*2), size = .05)
p10_b2
limits <- aes(ymax = risk + 2*error, ymin=risk - 2*error)
p10_b2 <- p10_b2 + geom_errorbar(limits) + geom_point()
p10_b2 <- p10_b2 + geom_smooth(data = data_acute_regression_plot, aes(dose,fit))
p10_b2

#factors work by just seperating them into different categories
data_test <- data_to_plot[data_to_plot$facet ==1,]
test_model <- lm(risk ~ I(dose+best_theta*dose^2), data = data_test, weights = 1/error^2)

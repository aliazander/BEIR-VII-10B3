library(plyr)
library(dplyr)
library(ggplot2)
library(stats4)

setwd("/Users/g-woloschak/Documents/")

##downlaod data for plot 10B-3
BEIR_VII_10B_3 <- read.csv('BEIR_VII_10B-3_data.csv')

##Seperate acute vs. chornic data with restrictions Ben discovered - 
## RFM mice and female only
acute_plot_10B3 <- BEIR_VII_10B_3[BEIR_VII_10B_3[,6] == "F" & 
                                        BEIR_VII_10B_3$type == "A" &  
                                        BEIR_VII_10B_3[,5] == "RFM",]
chronic_plot_10B3 <- BEIR_VII_10B_3[BEIR_VII_10B_3[,6] == "F" &
                                          BEIR_VII_10B_3$type == "C" &
                                          BEIR_VII_10B_3[,5] == "RFM",]

##replaces age with the average age when they have the same dose
acute_plot_10B3$age <- ave(acute_plot_10B3$age, acute_plot_10B3$dose, FUN = mean)
##delete duplicate dose points
acute_plot_10B3 <- acute_plot_10B3[!duplicated(acute_plot_10B3$dose),]

##combine acute and chronic back together for plot
data_plot_10B3 <- rbind(acute_plot_10B3, chronic_plot_10B3)

#plot all of the data points 
p_10B3 <- ggplot(data_plot_10B3, aes(x=dose, y = 1/age)) + 
      geom_text(family = "Helvetica", size = 3.8, aes(label = type)) +
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
            weights = 1/sd^2)
#extract constant, alpha, and theta for equations for lines
coefficients <- model$coefficients
predicted_xvalues <- as.data.frame(seq(from = 0, to = 1.5, by = .01))
#acute - constant + alpha*(D + theta*Dose^2)
acute_regression_data <- coefficients[1] + 
      coefficients[2] * predicted_xvalues + 
      coefficients[3] * predicted_xvalues^2

data_acute_regression_plot <- as.data.frame(cbind(predicted_xvalues, acute_regression_data))
names(data_acute_regression_plot) <- c("dose", "fit")

##chronic --> constant + alpha*D
chronic_regression_data <- coefficients[1] + 
      coefficients[2] * predicted_xvalues
data_chronic_regression_plot <- as.data.frame(cbind(predicted_xvalues, chronic_regression_data))
names(data_chronic_regression_plot) <- c("dose", "fit")

#add fit lines to plot
p_10B3 <- p_10B3 + geom_path(data=data_chronic_regression_plot, aes(dose, fit)) +
      geom_path(data = data_acute_regression_plot, aes(dose,fit))

#ta-da!
p_10B3  

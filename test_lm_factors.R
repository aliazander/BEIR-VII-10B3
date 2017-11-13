x<- c(0, 1, 2, 3, 0, 1, 2, 3)
y<- c(0, 1, 2, 3, 0, 0, 0, 0)
facet <- c(1, 1, 1, 1, 2, 2, 2 ,2)
test_data <- data.frame(x,y,facet)

plot<-ggplot(aes(x=x, y=y), data = test_data) + geom_point()
plot

test_data_1 <- test_data[test_data$facet == 1,]
model_test1 <- lm(y ~ x, data = test_data_1)
fit_line1 <- data.frame(predict(model_test1), test_data_1$x)
names(fit_line1)<-c("predicted","x")
plot_1 <- ggplot(aes(x=x, y=y), data = test_data_1) + geom_point() + 
      geom_path(aes(x=x, y=predicted), data=fit_line1)

test_data_2 <- test_data[test_data$facet == 2,]
model_test2 <- lm(y ~ x, data = test_data_2)
predict(model_test2)
fit_line2 <- data.frame(predict(model_test2), test_data_2$x)
names(fit_line2)<-c("predicted","x")
plot_2 <- ggplot(aes(x=x, y=y), data = test_data_2) + geom_point() + 
      geom_path(aes(x=x, y=predicted), data=fit_line2)
plot_2

combined_model <- lm(y ~ factor(facet)*x, data = test_data)
summary(combined_model)
coeff <- combined_model$coefficients
intercepts <- coeff[1:2]
coeffs <- coeff[3:4]
predicted_yvalues_1 <- intercepts[1] + coeffs[1] * test_data_1$x
predicted_yvalues_2 <- intercepts[1] + intercepts[2] + (coeffs[1] + coeffs[2]) * test_data_2$x
combined_fit_1 <- data.frame(predict_1 = predicted_yvalues_1, x_values_1 =test_data_1$x)
combined_fit_2 <- data.frame(predict_2 = predicted_yvalues_2, x_values_2 = test_data_2$x)
plot_1 <- plot_1 + geom_path(aes(x=x_values_1, y=predict_1), data=combined_fit_1)
plot_1
plot_2 <- plot_2 + geom_path(aes(x=x_values_2, y=predict_2), data=combined_fit_2)
plot_2


plots <- list()
unique(test_data$facet)
for(i in 1:length(unique(test_data$facet))){
      test_data_i <- test_data[test_data$facet == i,]
      model <- lm(y~x, data=test_data_i)
      fitline <- data.frame(predicted = predict(model), x = test_data_i$x)
      plot_i <- ggplot(aes(x=x, y=y), data = test_data_i) + geom_point() +
            geom_path(aes(x=x, y=predicted), data = fitline)
      plots[[i]] <- plot_i
      
}
do.call("grid.arrange", c(plots, ncol=2))



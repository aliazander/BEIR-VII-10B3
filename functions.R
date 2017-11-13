#Functions!

#get data frame of theta and loglikelihood values with linear quadratic
loglike_df <- function(o_range, data_for_function) {
      likelihood_values <- list()
      #create for loop to get likelihood values for range of theta values
      for(i in 1:length(o_range)){
            o<- o_range[i]
            model_one_theta <- lm(risk ~ factor(facet) * I(dose+o*dose^2), 
                                  data = data_for_function, 
                                  weights = (1/error^2))
            likelihood_values[i] <- logLik(model_one_theta)
      }
      
      #take list of likelihood values and put in a dataframe
      df <- data.frame(matrix(unlist(likelihood_values), nrow=length(o_range), byrow=T))
      #add theta values to dataframe
      df$theta_values <- o_range
      #label columns
      colnames(df) <- c("likelihoods", "thetas")
      df
}

#find the best log likelihood value (max)
find_bestlog <- function(likelihood_df){
      bestloglik <- max(df$likelihoods, na.rm = TRUE)
      bestloglik
}

#get 95% confidene interval for likelihood dataframe
confidence_interval <- function(df, level=.95){
      #find 95% CI
      
      bestloglik <- max(df$likelihoods, na.rm = TRUE)
      minloglik <- bestloglik - qchisq(.95, 1)/2
      confinterval <- data.frame(min = min(df[df$likelihoods >= minloglik,][,2]))
      confinterval$max <- max(df[df$likelihoods >= minloglik,][,2])
      confinterval
}

get_best_theta <- function(df){
      #get row number of max likelihood value
      row_number <- which(df$likelihoods == max(df$likelihoods, na.rm=TRUE))
      #use row number to get best theta value
      best_theta <- df[row_number, 2]
}

lq_best_theta <- function(best_theta, data_to_plot){
      #use best theta for model that keeps one theta
      model_best_theta <- lm(risk ~ factor(facet) * I(dose + best_theta*dose^2),
                             data = data_to_plot,
                             weights = 1/error^2)
}

ordered_df <- function(unordered, ordered_by, logical){
      unordered[order(unordered[ordered_by], decreasing = logical),]
}



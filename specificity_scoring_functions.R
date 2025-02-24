#contains the functions used to calculate the specificity scores 

z_score <- function(expr){
  sapply(seq_along(expr), function(i) {
    sd_val <- sd(expr[-i])
    if (sd_val == 0) { #nog naar kijken, 
      return((expr[i] - mean(expr[-i])) / 1)  
    } else {
      return((expr[i] - mean(expr[-i])) / sd_val)
    }
  })
}

proportion <- function(expr) {
  sapply(seq_along(expr), function(i){
    total <- sum(expr)
    return(expr[i] / total)
  })
}

calc_tau <- function(x) {
  max_x <- max(x)
  tau <- sum(1-(x/max_x)) / (length(x)-1)
  return(tau)
}



posterior_probs <- function(celltype_expression, mean, var) {
  likelihood <- dnorm(celltype_expression, mean = mean, sd = sqrt(var))
  prior <- mean / sum(mean)
  posterior <- likelihood * prior
  return(posterior)
}


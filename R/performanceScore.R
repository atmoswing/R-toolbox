CRPS <- function(x,x0,a,b){
  # r: rank
  # n: size of the sample
  # x : vector of length n
  # x0 : observed value
  # a, b : constants, that have a law-dependent optimum from which the samples are derived
  x.s <- sort(x)
  n <- length(x)
  r <- 1:n
  Fx <- (r-a)/(n+b)
  res <- 0
  # add rectangle on right side if observed value is on the right
  if (x0 > x.s[n]) {
    res <- res + x0 - x.s[n]
  }
  # add rectangle on the left side if observed value is on the left
  if (x0 < x.s[1]) {
    res <- res + x.s[1] - x0
  }
  if(n>1){
    for (i in 1:(n-1)) {
      # if value on left side of observed value
      if (x.s[i] < x0) {
        # case 1: next value also on left side of observed value
        if (x.s[i+1] <= x0) {
          res <- res + (Fx[i]^2+Fx[i+1]^2)*(x.s[i+1]-x.s[i])/2
        }
        # case 2: x0 before next value
        else { 
          F0 <- (Fx[i+1]-Fx[i]) * (x0-x.s[i])/(x.s[i+1]-x.s[i]) + Fx[i]
          res <- res + (F0^2+Fx[i]^2)*(x0-x.s[i])/2
          res <- res + ((F0-1)^2+(Fx[i+1]-1)^2)*(x.s[i+1]-x0)/2
        }
      }
      # if value on the right side of observed value
      else {
        res <- res + ((Fx[i]-1)^2+(Fx[i+1]-1)^2)*(x.s[i+1]-x.s[i])/2
      }
    }
  }
  return(res)
} 
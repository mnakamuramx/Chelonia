rm(list = ls())
library(ggplot2)
library(nlstools)

##  ............................................................................
##  Box-Cox transformation                                                  ####

bc_trans = function(y, lambda) {
  if (lambda == 0) {
    return(log(y))
  } else {
    return((y ^ lambda - 1) / lambda)
  }
}

##  ............................................................................
##  Inverse of Box-Cox transformation                                       ####

bc_trans_inv = function(x, lambda) {
  if (lambda == 0) {
    return(exp(x))
  } else {
    return((lambda * x + 1) ^ (1 / lambda))
  }
}

##  ............................................................................
##  Box-Cox transformation least squares diagnostics                        ####

diagnostics_nls = function(fit) {
  resids = residuals(fit)
  preds = fitted.values(fit)
  obs = resids + preds
  n = length(resids)
  fitdata = data.frame(Predicted = preds, Residuals = resids)
  nlsresids = nlsResiduals(fit)
  plot(nlsresids)
  print(test.nlsResiduals(nlsresids))
  
  # Sequential plot of residuals
  
  p1 = ggplot(fitdata, aes(y = resids, x = 1:n)) +
    geom_point() +
    ylab("Residuals") +
    xlab("Index") +
    ggtitle("Residuals") +
    geom_hline(yintercept = 0, col = "red")
  print(p1)
  
  # qq normal plot residuals
  
  p2 = ggplot(fitdata, aes(sample = resids)) +
    geom_qq() +
    stat_qq_line(col = "red") +
    ggtitle("Q-Q plot")
  print(p2)
  
  # Histogram of residuals
  
  p3 = ggplot(fitdata, aes(x = resids)) +
    geom_histogram(aes(y = after_stat(density)),
                   bins = log(n) / log(2),
                   fill = "darkgray") +
    stat_function(fun = dnorm, args = list(mean = mean(resids),
                                           sd = sd(resids))) +
    ggtitle("Residuals")
  print(p3)
  
  # Observed vs. predicted
  
  p4 = ggplot(fitdata, aes(x = preds, y = obs)) +
    geom_point() +
    geom_abline(slope = 1,
                intercept = 0,
                col = "red") +
    xlab("Predicted") +
    ylab("Observed") +
    ggtitle("Observed vs. predicted")
  print(p4)
  
  # Predicted vs. residuals
  
  p5 = ggplot(fitdata, aes(x = preds, y = resids)) +
    geom_point() +
    geom_hline(yintercept = 0, col = "red") +
    xlab("Predicted") +
    ylab("Residuals") +
    ggtitle("Residuals vs. predicted")
  print(p5)
  
  # Kolmogorov-Smirnov test of normality
  
  print(ks.test(scale(resids), "pnorm"))
  
  # Autocorrelation test
  
  acf(resids, main = "Residual autocorrelation")
  
}

##  ............................................................................

##  Box-Cox transformation fit and diagnostics with AMO.26                  ####

fit_bc_trans = function(datos,
                        plotflag,
                        diagnosticsflag,
                        printflag,
                        lambda) {
  nls.control(maxiter = 5000, warnOnly = TRUE)
  
  fit = NULL
  try(fit <- nls(
    bc_trans(Hatchlings, lambda) ~
      bc_trans((beta0) * Females ^ (gamma0 + gamma1 * AMO.26),
               lambda),
    data = datos,
    
    # Provide starting values for least squares here:
    
    start = c(
      beta0 = 264.29,
      gamma0 = 0.998,
      gamma1 = -0.61
    ),
    algorithm = "port",
    lower = c(0, -Inf, -Inf),
    upper = c(Inf, Inf, Inf)
  ))
  
  if (is.null(fit)) {
    return(fit)
  }
  if (printflag) {
    cat("Summary fit:", "\n")
    print(summary(fit))
  }
  
  if (diagnosticsflag == TRUE) {
    diagnostics_nls(fit)
  }
  
  coefs4 = coefficients(fit)
  if (printflag) {
    print(summary(fit)$parameters[, 1:2])
  }
  sigma = summary(fit)$sigma
  covariance = vcov(fit)
  
  a = (coefs4[1])
  b = coefs4[2] + coefs4[3] * datos$AMO.26
  hatchlings_med = a * datos$Females ^ b
  
  hatchlings_025 = bc_trans_inv(bc_trans(hatchlings_med, lambda) +
                                  qnorm(.025, mean = 0, sd = sigma),
                                lambda)
  hatchlings_975 = bc_trans_inv(bc_trans(hatchlings_med, lambda) +
                                  qnorm(0.975, mean = 0, sd = sigma),
                                lambda)
  densodep = a * b * datos$Females ^ (b - 1)
  if (printflag) {
    cat("Correlation (sst, densodep): ",
        cor(densodep, datos$AMO.26),
        "\n")
    print(cor.test(densodep, datos$AMO.26))
    
    cat("Correlation (sst, b): ",
        cor(b, datos$AMO.26),
        "\n")
    print(cor.test(b, datos$AMO.26))
    
    cat(
      "Correlation (observed, predicted): ",
      cor(datos$Hatchlings, hatchlings_med),
      "\n"
    )
    print(cor.test(datos$Hatchlings, hatchlings_med))
    
    cat("sigma parameter: ", sigma, "\n")
    
  }
  databio <<- data.frame(
    datos,
    Hatchlings_pred = hatchlings_med,
    DensoDep = densodep,
    A = a,
    B = b,
    Q025 = hatchlings_025,
    Q975 = hatchlings_975,
    row.names = NULL
  )
  if (printflag) {
    write.csv(databio, "databio.csv")
  }
  
  if (plotflag == TRUE) {
    p1 = ggplot(databio, aes(y = Hatchlings_pred, x = Year)) +
      geom_line(col = "red") +
      geom_point(aes(x = Year, y = Hatchlings), col = "blue") +
      xlab("Year") +
      ylab("Hatchlings") +
      ggtitle("Predicted and observed") +
      geom_line(aes(x = Year, y = Q025)) +
      geom_line(aes(x = Year, y = Q975))
    print(p1)
    
    p2 = ggplot(databio, aes(x = Year, ymin = Q025, ymax = Q975)) +
      geom_errorbar() +
      geom_point(aes(x = Year, y = hatchlings_med)) +
      geom_point(aes(x = Year, y = Hatchlings), col = "blue")
    print(p2)
    
    p3 = ggplot(databio, aes(x = Year, y = AMO.26)) +
      geom_line() +
      ggtitle("AMO.26")
    print(p3)
    
    p4 = ggplot(databio, aes(x = Year, y = a)) +
      geom_line() +
      ggtitle("A parameter")
    print(p4)
    
    p5 = ggplot(databio, aes(x = Year, y = b)) +
      geom_line() +
      ggtitle("B parameter")
    print(p5)
    
    p6 = ggplot(databio, aes(x = Year, y = DensoDep)) +
      geom_line() +
      ggtitle("Hatchlings/Female")
    print(p6)
    
    p7 = ggplot(databio, aes(x = AMO.26, a)) +
      geom_line() +
      ggtitle("A parameter")
    print(p7)
    
    p8 = ggplot(databio, aes(x = AMO.26, b)) +
      geom_line() +
      ggtitle("B parameter")
    print(p8)
  }
  return(fit)
}

#   ____________________________________________________________________________
#   Read raw data and prepare for regression                                ####

datos = read.csv("Chelonia_data.csv")
datos = datos[, c("Year", "Females", "Hatchlings", "SEA.LEVEL", "AMO.26")]

datosreg = data.frame(
  datos,
  v1 = datos$SEA.LEVEL * log(datos$Females),
  v2 = datos$AMO.26 * log(datos$Females)
)

#   ____________________________________________________________________________
#   Run Box-Cox transformation fit                                          ####

### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
### Main call                                                               ####

fit4 = fit_bc_trans(
  datosreg,
  plotflag = TRUE,
  diagnosticsflag = TRUE,
  printflag = TRUE,
  lambda = 0.5
)

### . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ..
### Extract fitted values                                                   ####

coefs4 = coefficients(fit4)
summary(fit4)
confint(fit4)

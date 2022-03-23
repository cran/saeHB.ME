#' @title Small Area Estimation with Measurement Error using Hierarchical Bayesian Method under Normal Distribution
#' @description This function is implemented to variable of interest \eqn{(y)} that assumed to be a Normal Distribution when auxiliary variable is measured with error.
#' @param formula  an object of class \code{\link[stats]{formula}} (or one that can be coerced to that class): a symbolic description of the model to be fitted. The variables included \code{formula} must have a length equal to the number of domains \code{m}. This formula can provide auxiliary variable either measured with error or combination between measured with \code{error} and \code{without error}. If the auxiliary variable are combination between \code{error} and \code{without error}, input the \code{error} variable first followed by \code{without error} variable.
#' @param vardir vector containing the \code{m} sampling variances of direct estimators for each domain. The values must be sorted as the \code{Y}.
#' @param var.x vector containing mean squared error of \code{X}. The values must be sorted as the \code{X}.
#' @param coef a vector contains prior initial value of Coefficient of Regression Model for fixed effect with default vector of \code{0} with the length of the number of regression coefficients.
#' @param var.coef a vector contains prior initial value of variance of Coefficient of Regression Model with default vector of \code{1} with the length of the number of regression coefficients.
#' @param tau.u  prior initial value of inverse of Variance of area random effect with default \code{1}.
#' @param iter.update number of updates with default \code{3}.
#' @param iter.mcmc number of total iterations per chain with default \code{10000}.
#' @param thin thinning rate, must be a positive integer with default \code{2}.
#' @param burn.in number of iterations to discard at the beginning with default \code{2000}.
#' @param data the data frame.
#'
#' @return This function returns a list with the following objects:
#'    \item{Est}{A vector with the values of Small Area mean Estimates using Hierarchical bayesian method }
#'    \item{refVar}{Estimated random effect variances}
#'    \item{coefficient}{A data frame with the estimated model coefficient}
#'    \item{plot}{Trace, Dencity, Autocorrelation Function Plot of MCMC samples}
#'
#' @export meHBNormal
#'
#' @examples
#' ## Load dataset
#' data(dataHBME)
#'
#' ## Auxiliary variables only contains variable with error
#' example <- meHBNormal(Y~x1+x2, vardir = "vardir",
#'                    var.x = c("v.x1","v.x2"), iter.update = 3, iter.mcmc = 10000,
#'                    thin = 5, burn.in = 1000, data = dataHBME)
#'
#' ## Auxiliary variables contains variable with error and without error
#' example_mix <- meHBNormal(Y~x1+x2+x3, vardir = "vardir",
#'                    var.x = c("v.x1","v.x2"), iter.update = 3, iter.mcmc = 10000,
#'                    thin = 5, burn.in = 1000, data = dataHBME)
#'
#' ## Create dataset with nonsampled area
#' dataHBMEns <- dataHBME
#' dataHBMEns[c(1,10,20,30),"Y"] <- NA
#'
#' ## For data with nonsampled area use dataHBMEns
#'
meHBNormal <- function(formula, vardir, var.x, coef, var.coef,
                       iter.update=3, iter.mcmc=10000, thin = 2, tau.u = 1, burn.in =2000, data){
  result <- list(Est = NA, refVar = NA, coefficient = NA,
                 plot=NA)
  formuladata <- model.frame(formula, na.action = NULL, data)
  n <- nrow(formuladata)
  if (any(is.na(formuladata[,-1]))){
    stop("Auxiliary Variables contains NA values.")
  }
  X <- model.matrix(formula, data)
  p <- dim(X)[2]

  if(!missing(var.coef)){
    if (length(var.coef) != p) {
      stop("length of vector var.coef does not match the number of regression coefficients, the length must be ",
           p)
    }
    tau.b <- 1/var.coef
  } else {
    tau.b <- 1/rep(1, p)
  }
  if(!missing(coef)){
    if (length(coef) != p) {
      stop("length of vector coef does not match the number of regression coefficients, the length must be ",
           p)
    }
    mu.b <- coef
  } else {
    mu.b <- rep(0, p)
  }
  if (iter.update < 3){
    stop("the number of iteration updates at least 3 times")
  }

  tau.ua = tau.ub = 1
  c <- data[, var.x]
  c <- as.data.frame(c)
  n_c <- dim(c)[2]
  psi <- data[, vardir]

  if (!any(is.na(formuladata[,1]))){
    formuladata <- model.frame(formula, na.action = na.omit, data)
    X <- model.matrix(formula, data)
    y <- formuladata[,1]
    m <- length(y)

    if(n_c < p-1){
      for (iter in 1:iter.update) {
        dat <- list(m = m, y = y, Xme=as.matrix(formuladata[,2:(n_c+1)]), x=as.matrix(formuladata[,(n_c+2):p]),
                    mu.b=mu.b, p = p , var.x=c, n_c = n_c,
                    tau.b=tau.b,
                    tau.ua=tau.ua, tau.ub=tau.ub, vardir=psi)  # names list of numbers
        inits <- list(u = rep(0,m), b = mu.b, tau.u = tau.u)

        cat("model {
					for (i in 1:m) {
							y[i] ~ dnorm(mu[i],tau[i])

              for (j in 1:(n_c)){
                em[i,j] ~ dnorm(Xme[i,j],te[i,j])
							}

							for(I in 1:(n_c)){
                kjk[i,I] <- b[I+1]*em[i,I]
							}

							mu[i] <- b[1] + sum(kjk[i,]) + sum(b[(n_c+2):p]*x[i,])+u[i]

							for (k in 1:(n_c)){
							  te[i,k] <- 1/var.x[i,k]
							}
							tau[i] <- 1/vardir[i]

							u[i] ~ dnorm(0,tau.u)
					}

					for (l in 1:p){
					  b[l]~dnorm(mu.b[l], tau.b[l])
					}

					tau.u ~ dgamma(tau.ua,tau.ub)
					a.var <- 1 / tau.u

			}", file="ME.txt")

        jags.m <- jags.model( file = "ME.txt", data=dat, inits=inits, n.chains=1, n.adapt=500)
        file.remove("ME.txt")
        params <- c("mu","a.var","b", "tau.u")
        samps1 <- coda.samples( jags.m, params, n.iter=iter.mcmc, thin=thin)
        samps11 <- window(samps1, start=burn.in + 1, end=iter.mcmc)

        result_samps=summary(samps11)
        a.var=result_samps$statistics[1]
        beta=result_samps$statistics[2:(p+1),1:2]
        for (i in 1:p){
          mu.b[i]  = beta[i,1]
          tau.b[i] = 1/(beta[i,2]^2)
        }

        tau.ua = result_samps$statistics[2+p+m,1]^2/result_samps$statistics[2+p+m,2]^2
        tau.ub = result_samps$statistics[2+p+m,1]/result_samps$statistics[2+p+m,2]^2
      }
    } else {
      for (iter in 1:iter.update) {
        dat <- list(m = m, y = y, x=as.matrix(formuladata[,2:p]), mu.b=mu.b, p = p , var.x=c,
                    tau.b=tau.b,
                    tau.ua=tau.ua, tau.ub=tau.ub, vardir=psi)  # names list of numbers
        inits <- list(u = rep(0,m), b = mu.b, tau.u = tau.u)

        cat("model {
					for (i in 1:m) {
							y[i] ~ dnorm(mu[i],tau[i])

              for (j in 1:(p-1)){
                em[i,j] ~ dnorm(x[i,j],te[i,j])
							}

							for(I in 1:(p-1)){
                kjk[i,I] <- b[I+1]*em[i,I]
              }
							mu[i] <- b[1] + sum(kjk[i,])+u[i]

							for (k in 1:(p-1)){
							  te[i,k] <- 1/var.x[i,k]
							}
							tau[i] <- 1/vardir[i]

							u[i] ~ dnorm(0,tau.u)
					}

					for (l in 1:p){
					  b[l]~dnorm(mu.b[l], tau.b[l])
					}

					tau.u ~ dgamma(tau.ua,tau.ub)
					a.var <- 1 / tau.u

			}", file="ME.txt")

        jags.m <- jags.model( file = "ME.txt", data=dat, inits=inits, n.chains=1, n.adapt=500)
        file.remove("ME.txt")
        params <- c("mu","a.var","b", "tau.u")
        samps1 <- coda.samples( jags.m, params, n.iter=iter.mcmc, thin=thin)
        samps11 <- window(samps1, start=burn.in + 1, end=iter.mcmc)

        result_samps=summary(samps11)
        a.var=result_samps$statistics[1]
        beta=result_samps$statistics[2:(p+1),1:2]
        for (i in 1:p){
          mu.b[i]  = beta[i,1]
          tau.b[i] = 1/(beta[i,2]^2)
        }

        tau.ua = result_samps$statistics[2+p+m,1]^2/result_samps$statistics[2+p+m,2]^2
        tau.ub = result_samps$statistics[2+p+m,1]/result_samps$statistics[2+p+m,2]^2
      }
    }
    result_samps=summary(samps1)
    b.varnames <- list()
    for (i in 1:(p)) {
      idx.b.varnames <- as.character(i-1)
      b.varnames[i] <-str_replace_all(paste("b[",idx.b.varnames,"]"),pattern=" ", replacement="")
    }

    result_mcmc <- samps1[,c(2:(p+1))]
    colnames(result_mcmc[[1]]) <- b.varnames

    a.var=result_samps$statistics[1]

    beta=result_samps$statistics[2:(p+1),1:2]
    rownames(beta) <- b.varnames

    mu=result_samps$statistics[(p+2):(1+p+m),1:2]

    Estimation=data.frame(mu)
    colnames(Estimation)=c("mean","sd")

    Quantiles <- as.data.frame(result_samps$quantiles[1:(2+p+m),])
    q_beta <- (Quantiles[2:(p+1),])
    q_Estimation <- (Quantiles[(p+2):(p+1+m),])
    rownames(q_beta) <- b.varnames
    beta <- cbind(beta,q_beta)
  } else{

    formuladata <- as.data.frame(formuladata)
    X <- model.matrix(formula, data)
    y <- formuladata[,1]
    m <- length(y)

    formuladata$idx <- rep(1:m)
    formuladata[,2+p] <- psi
    formuladata <- cbind(formuladata, c)
    data_sampled <- na.omit(formuladata)
    psi_sampled <- data_sampled[, 2+p]
    c_sampled <- data_sampled[, var.x]
    data_nonsampled <- formuladata[-data_sampled$idx,]
    c_nonsampled <- data_nonsampled[, var.x]
    r <- data_nonsampled$idx
    m1 <- nrow(data_sampled)
    m2 <- nrow(data_nonsampled)

    if(n_c < p-1){
      for (iter in 1:iter.update) {
        dat <- list(m1 = m1, m2 = m2, y_sampled = data_sampled[,1],
                    Xme_sampled = as.matrix(data_sampled[,2:(n_c+1)]), Xme_nonsampled = as.matrix(data_nonsampled[,2:(n_c+1)]),
                    x_sampled = as.matrix(data_sampled[,(n_c+2):p]), x_nonsampled = as.matrix(data_nonsampled[,(n_c+2):p]),
                    mu.b = mu.b, p = p, n_c = n_c,
                    tau.b = tau.b,
                    tau.ua = tau.ua, tau.ub = tau.ub, vardir = psi_sampled, var.x = c_sampled,
                    var.x.T = c_nonsampled)  # names list of numbers
        inits <- list(u = rep(0,m1), b = mu.b, tau.u = tau.u)

        cat("model {
					for (i in 1:m1) {
							y_sampled[i] ~ dnorm(mu[i],tau[i])

              for (j in 1:(n_c)){
                em[i,j] ~ dnorm(Xme_sampled[i,j],te[i,j])
							}

              for(I in 1:(n_c)){
                kjk[i,I] <- b[I+1]*em[i,I]
              }

							mu[i] <- b[1] + sum(kjk[i,]) + sum(b[(n_c+2):p]*x_sampled[i,]) + u[i]

							for (k in 1:(n_c)){
							  te[i,k] <- 1/var.x[i,k]
							}
							tau[i] <- 1/vardir[i]

							u[i] ~ dnorm(0,tau.u)
					}

					for (l in 1:p){
					  b[l]~dnorm(mu.b[l], tau.b[l])
					}

					tau.u ~ dgamma(tau.ua,tau.ub)
					a.var <- 1 / tau.u

					for(j in 1:m2){
					  v[j]~dnorm(0,1/a.var)

					  for(J in 1:(n_c)){
					    em.T[j,J] ~ dnorm(Xme_nonsampled[j,J], te.T[j,J])
					  }
					  for(I in 1:(n_c)){
                kjk.T[j,I] <- mu.b[I+1]*em.T[j,I]
					  }

					  mu.T[j]<- mu.b[1] + sum(kjk.T[j,]) + sum(mu.b[(n_c+2):p]*x_nonsampled[j,])+ v[j]
					   for(G in 1:(n_c)){
					    te.T[j,G] <-  1/var.x.T[j,G]
					  }
					}

			}", file="ME.txt")

        jags.m <- jags.model( file = "ME.txt", data=dat, inits=inits, n.chains=1, n.adapt=500)
        file.remove("ME.txt")
        params <- c("mu","mu.T","a.var","b", "tau.u")
        samps1 <- coda.samples( jags.m, params, n.iter=iter.mcmc, thin=thin)
        samps11 <- window(samps1, start=burn.in + 1, end=iter.mcmc)

        result_samps=summary(samps11)
        a.var=result_samps$statistics[1]
        beta=result_samps$statistics[2:(p+1),1:2]
        for (i in 1:p){
          mu.b[i]  = beta[i,1]
          tau.b[i] = 1/(beta[i,2]^2)
        }

        tau.ua = result_samps$statistics[2+p+m,1]^2/result_samps$statistics[2+p+m,2]^2
        tau.ub = result_samps$statistics[2+p+m,1]/result_samps$statistics[2+p+m,2]^2
      }
    } else {
      for (iter in 1:iter.update) {
        dat <- list(m1 = m1, m2 = m2, y_sampled = data_sampled[,1],
                    x_sampled = as.matrix(data_sampled[,2:p]), x_nonsampled = as.matrix(data_nonsampled[,2:p]),
                    mu.b = mu.b, p = p,
                    tau.b = tau.b,
                    tau.ua = tau.ua, tau.ub = tau.ub, vardir = psi_sampled, var.x = c_sampled,
                    var.x.T = c_nonsampled)  # names list of numbers
        inits <- list(u = rep(0,m1), b = mu.b, tau.u = tau.u)

        cat("model {
					for (i in 1:m1) {
							y_sampled[i] ~ dnorm(mu[i],tau[i])

              for (j in 1:(p-1)){
                em[i,j] ~ dnorm(x_sampled[i,j],te[i,j])
							}

              for(I in 1:(p-1)){
                kjk[i,I] <- b[I+1]*em[i,I]
              }
							mu[i] <- b[1] + sum(kjk[i,])+u[i]

							for (k in 1:(p-1)){
							  te[i,k] <- 1/var.x[i,k]
							}
							tau[i] <- 1/vardir[i]

							u[i] ~ dnorm(0,tau.u)
					}

					for (l in 1:p){
					  b[l]~dnorm(mu.b[l], tau.b[l])
					}

					tau.u ~ dgamma(tau.ua,tau.ub)
					a.var <- 1 / tau.u

					for(j in 1:m2){
					  mu.T[j]<- mu.b[1]+sum(kjk.T[j,])+v[j]
					  v[j]~dnorm(0,1/a.var)
					  for(J in 1:(p-1)){
					    em.T[j,J] ~ dnorm(x_nonsampled[j,J], te.T[j,J])
					  }
					  for(I in 1:(p-1)){
                kjk.T[j,I] <- mu.b[I+1]*em.T[j,I]
              }
					   for(G in 1:(p-1)){
					    te.T[j,G] <-  1/var.x.T[j,G]
					  }
					}

			}", file="ME.txt")

        jags.m <- jags.model( file = "ME.txt", data=dat, inits=inits, n.chains=1, n.adapt=500)
        file.remove("ME.txt")
        params <- c("mu","mu.T","a.var","b", "tau.u")
        samps1 <- coda.samples( jags.m, params, n.iter=iter.mcmc, thin=thin)
        samps11 <- window(samps1, start=burn.in+1, end=iter.mcmc)

        result_samps=summary(samps11)
        a.var=result_samps$statistics[1]
        beta=result_samps$statistics[2:(p+1),1:2]
        for (i in 1:p){
          mu.b[i]  = beta[i,1]
          tau.b[i] = 1/(beta[i,2]^2)
        }

        tau.ua = result_samps$statistics[2+p+m,1]^2/result_samps$statistics[2+p+m,2]^2
        tau.ub = result_samps$statistics[2+p+m,1]/result_samps$statistics[2+p+m,2]^2
      }
    }

    result_samps=summary(samps1)
    b.varnames <- list()
    for (i in 1:(p)) {
      idx.b.varnames <- as.character(i-1)
      b.varnames[i] <-str_replace_all(paste("b[",idx.b.varnames,"]"),pattern=" ", replacement="")
    }

    result_mcmc <- samps1[,c(2:(p+1))]
    colnames(result_mcmc[[1]]) <- b.varnames

    a.var <- result_samps$statistics[1]

    beta <- result_samps$statistics[2:(p+1),1:2]
    rownames(beta) <- b.varnames

    mu <- result_samps$statistics[(p+2):(1+p+m1),1:2]
    mu.nonsampled <- result_samps$statistics[(2+p+m1):(1+p+m),1:2]

    Estimation <- matrix(rep(0,m),nrow = m, ncol = 2)
    Estimation[r,] <- mu.nonsampled
    Estimation[-r,] <- mu
    Estimation <- as.data.frame(Estimation)
    colnames(Estimation) <- c("mean","sd")

    Quantiles <- as.data.frame(result_samps$quantiles[1:(2+p+m),])
    q_beta <- (Quantiles[2:(p+1),])
    q_mu <- (Quantiles[(p+2):(p+1+m1),])
    q_mu.nonsampled <- (Quantiles[(2+p+m1):(1+p+m),])
    q_Estimation <- matrix(0,m,5)
    for (i in 1:5) {
      q_Estimation[r,i] <- q_mu.nonsampled[,i]
      q_Estimation[-r,i] <- q_mu[,i]
    }
    rownames(q_beta) <- b.varnames
    beta <- cbind(beta,q_beta)
  }

  result$Est = cbind(Estimation, q_Estimation)
  result$refVar      = a.var
  result$coefficient = beta
  result$plot       = list(graphics.off(),
                           par(mar=c(2,2,2,2)),
                           autocorr.plot(result_mcmc,col="brown2",lwd=2),
                           plot(result_mcmc,col="brown2",lwd=2))
  on.exit(result)
  return(result)
}

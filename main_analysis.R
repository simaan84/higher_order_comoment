library(quantmod)
library(lubridate)
library(stargazer)
library(xtable)
library(mixtools)
library(parallel)
library(data.table)
library(plyr)
library(moments)
library(boot)
library(nortest)
library(PerformanceAnalytics)
library(forecast)
library(rvest)
library(MASS)


rm(list = ls())
gc()

output.dir <- "latex_folder"
output.dir.fig <- paste(output.dir,"/Figures/",sep = "")
output.dir.tables <- paste(output.dir,"/Tables/",sep = "")


###############################################
############### LOAD DATA FROM GITHUB #########

crsp_file <- "https://raw.githubusercontent.com/simaan84/higher_order_comoment/main/snp_crsp.csv"
ff_file <- "https://raw.githubusercontent.com/simaan84/higher_order_comoment/main/FF_BTM.csv"

xts_all <- read.csv(crsp_file,stringsAsFactors = F)
rownames(xts_all) <- xts_all$date
xts_all$date <- NULL
xts_all <- as.xts(xts_all)

FF_data <- read.csv(ff_file,stringsAsFactors = F)
rownames(FF_data) <- FF_data$date
FF_data$date <- NULL
FF_data <- as.xts(FF_data)

all_data_list <- list(xts_all,FF_data)

###############################################################################################################

##################################################
############# ESTIMATION FUNCTION ################
##################################################

# the following is the main estimation function
# for a given r time series T times d, it returns the model estimates
# N stands for the number of times we perform the estimation
# elapse_second denotes the number of seconds to halt after no convergence

theta_f <- function(r,N,elapse_second) {
  r <- na.omit(r)
  
  MMNM_est <- function(r) {
    setTimeLimit(elapse=elapse_second, trans=T)
    mix1 <- mvnormalmixEM(r,k = 2,arbvar =F,verb = TRUE)
    states <- 1:2
    state1 <- which.max(mix1$lambda)
    state2 <- states[!states %in% state1] 
    
    SIGMA <- mix1$sigma
    PI <- min(mix1$lambda)
    ETA <- mix1$mu[[state1]]
    ETA_ALPHA <- mix1$mu[[state2]]
    GAMMA_tilde <- ETA_ALPHA - ETA
    post_prob <- mix1$posterior[,c(state1,state2)]
    
    list(PI = PI, SIGMA = SIGMA,  ETA = ETA, ETA_ALPHA = ETA_ALPHA,GAMMA_tilde = GAMMA_tilde,post_prob = post_prob, loglik = mix1$loglik)
  }
  
  results1 <- mclapply(1:N,function(i) {set.seed(i); return(try(MMNM_est(r),silent = T))},mc.cores = detectCores() )
  results1 <- results1[!sapply(results1,function(x) inherits(x,"try-error")  )]
  
  keep_results1 <- sapply(results1, function(x) x$PI) > 0.05 &  sapply(results1, function(x) x$PI) < 0.5
  results1_keep <- results1[keep_results1]
  
  results1_log <- sapply(results1_keep, function(x) x$loglik )
  
  # keep results with maximum likelihood
  keep_results2 <- results1_log   == max(results1_log)
  
  results1_keep <- results1_keep[keep_results2]
  
  PI <- mean(sapply(results1_keep, function(x) x$PI))
  
  SIGMA_LIST <- lapply(results1_keep, function(x) x$SIGMA)
  SIGMA <- Reduce("+", SIGMA_LIST) / length(SIGMA_LIST)
  
  ETA_LIST <-  lapply(results1_keep, function(x) x$ETA)
  ETA <- Reduce("+", ETA_LIST) / length(ETA_LIST)
  
  ETA_ALPHA_LIST <-  lapply(results1_keep, function(x) x$ETA_ALPHA)
  ETA_ALPHA <- Reduce("+", ETA_ALPHA_LIST) / length(ETA_ALPHA_LIST)
  
  GAMMA_tilde_LIST <-  lapply(results1_keep, function(x) x$GAMMA_tilde)
  GAMMA_tilde <- Reduce("+", GAMMA_tilde_LIST) / length(GAMMA_tilde_LIST)
  
  post_prob_list <- lapply(results1_keep, function(x) x$post_prob[,1] )
  post_prob <- Reduce("+", post_prob_list) / length(post_prob_list)
  post_prob <- cbind(post_prob,1-post_prob)
  
  list(PI = PI, SIGMA = SIGMA,  ETA = ETA, ETA_ALPHA = ETA_ALPHA,GAMMA_tilde = GAMMA_tilde,post_prob = post_prob)
  
}                                       


###############################################
########### NUMERICAL EXAMPLE #################
###############################################

# The following is an illustration of the portfolio solution for simple data

{
  d <- ncol(xts_all)
  set.seed(13)
  SELECT_STOCKS <- sort(sample(1:d,10))
  R_choose <- xts_all["2005-01-01/2010-01-01",SELECT_STOCKS]
  THETA <- theta_f(R_choose,100,20)
  PI <- THETA$PI
  SIGMA <- THETA$SIGMA
  GAMMA_tilde <- THETA$GAMMA_tilde
  ETA <- THETA$ETA
  ETA_ALPHA <- ETA + GAMMA_tilde
  E_R <- ETA + PI*GAMMA_tilde
  V_R <- SIGMA + PI*(1-PI)*GAMMA_tilde%*%t(GAMMA_tilde)
  
  # consider the optimal portfolio with risk-free asset
  A <- 1
  k_function <- function(s)  PI/(((1-PI))*exp(A*s) + PI )
  
  s_seq <- seq(-10,10,length = 10^3)
  # function k
  k_seq<- sapply(s_seq,k_function)
  plot(k_seq~s_seq,type = "l")
  
  # let's write the optimal portfolio as a function of s
  xi_f <- function(s) {
    GAMMA_tilde <- GAMMA_tilde
    port <- (1/A)*solve(SIGMA)%*%(ETA + k_function(s)*GAMMA_tilde )
    gamma_p <- t(port)%*%GAMMA_tilde
    list(port,(gamma_p - s)*10^6,gamma_p )
  }
  
  s <- uniroot(function(s) xi_f(s)[[2]], c(-100,100),tol = 10^-16 )$root
  cat("the gamma_p is ",s*10^6, "\n" )
  
  s_approx <- (1/A)*  (c(t(GAMMA_tilde)%*%solve(SIGMA)%*%ETA) +  PI*c(t(GAMMA_tilde)%*%solve(SIGMA)%*%GAMMA_tilde) ) /  (  1 + PI*(1-PI)* c(t(GAMMA_tilde)%*%solve(SIGMA)%*%GAMMA_tilde)   )
  s_approx*10^6
  k_approx <- PI - A*PI*(1-PI)*s_approx
  k_function(s)
  
  # as soon we find k solution, we can determine the two funds portfolio
  opt_port_approx <- (1/A)*solve(SIGMA)%*%(ETA + k_approx*GAMMA_tilde  )
  
  opt_port <- xi_f(s)[[1]]
  plot(opt_port~opt_port_approx)
  abline(a = 0,b = 1)
  
  # create a function for opt portfolio as a function of A and shock 
  my_port_f <- function(A,stock = 1,shock = 0) { 
    k_function <- function(s)  PI/(((1-PI))*exp(A*s) + PI )
    
    # find s
    xi_f <- function(s) {
      GAMMA_tilde[stock] <- GAMMA_tilde[stock]+shock
      port <- (1/A)*solve(SIGMA)%*%(ETA + k_function(s)*GAMMA_tilde )
      gamma_p <- t(port)%*%GAMMA_tilde
      list(port,(gamma_p - s)*10^6,gamma_p )
    }
    
    s <- uniroot(function(s) xi_f(s)[[2]], c(-100,100),tol = 10^-16 )$root
    opt_port <- xi_f(s)[[1]]
    return(opt_port)
    
  }
  
  # solve for different portfolios
  shock_seq <- seq(-0.01,0.01,length = 100)
  stock <- 1
  cor_mat <- cor(R_choose)
  cor_mat <- cor_mat[stock,]
  cor_mat <- cor_mat[cor_mat!=1]
  max_cor_stock <- which(names(R_choose) ==  names(which.max(cor_mat)))
  min_cor_stock <- which(names(R_choose) ==  names(which.min(cor_mat)))
  cor(R_choose[,c(stock,max_cor_stock,min_cor_stock)])[,1]
  
  opt_port_list_1 <- sapply(shock_seq, function(shock)  my_port_f(1,stock,shock)  )
  opt_port_list_2 <- sapply(shock_seq, function(shock)  my_port_f(2,stock,shock)  )
  opt_port_list_100 <- sapply(shock_seq, function(shock)  my_port_f(100,stock,shock)  )
  slopes <- apply(opt_port_list_1, 1,  function(x) lm(x  ~ shock_seq)[[1]][[2]] )
  
  # this plot portfolio weight sensitivity to shocks in gamma tilde estimation
  plot(opt_port_list_1[stock,] ~ shock_seq,type = "l", ylab = expression(xi[i]),xlab = expression(gamma[i]) )
  for (i in which(slopes > 0)) {
    lines(opt_port_list_1[i,] ~ shock_seq,type = "l", col = i) 
  }
  grid(10)
  
  
  plot(opt_port_list_1[stock,] ~ shock_seq,type = "l", ylab = expression(xi[i]),xlab = expression(gamma[i]) )
  lines(opt_port_list_1[max_cor_stock,] ~ shock_seq,type = "l", col = 3) 
  lines(opt_port_list_1[min_cor_stock,] ~ shock_seq,type = "l", col = 2) 
  grid(10)
  
}


###########################################
#### SIMULATED DATA BASED ON ESTIMATES ####
###########################################


## The following code estimates the estimation risk of each input using simulation study

scale_par <- 1

sim_f <- function(n,sims = sims) {
  cat("This is trial ",n,"\n")
  
  set.seed(n)
  B_sim <- rbinom(sims,1,PI)
  set.seed(n)
  U <- mvtnorm::rmvnorm(sims,mean = ETA,sigma = SIGMA)
  GAMMA_tilde_mat <- matrix(GAMMA_tilde,sims,length(GAMMA_tilde),byrow = T) 
  
  R_sim_jump <- GAMMA_tilde_mat*B_sim
  R_sim_gauss <-  U
  
  # the simualted return is the sum of the two
  R_sim <-  R_sim_gauss +  R_sim_jump
  
  THETA_hat <- theta_f(R_sim,100,20)
  PI_hat <- THETA_hat$PI
  SIGMA_hat <- THETA_hat$SIGMA
  GAMMA_tilde_hat <- THETA_hat$GAMMA_tilde
  ETA_hat <- THETA_hat$ETA
  ETA_ALPHA_hat <- ETA_hat + GAMMA_tilde_hat
  
  list(ETA_hat = ETA_hat,SIGMA_hat = SIGMA_hat, GAMMA_tilde_hat = GAMMA_tilde_hat,PI_hat = PI_hat)
  
}

n <- 17
sim_f_100 <- sim_f(n,10^2) 
sim_f_1000 <- sim_f(n,10^3) 
sim_f_10000 <- sim_f(n,10^4) 

cex_size <- 2

{
  file.i <- paste(output.dir.fig,"EM_PI.pdf",sep = "")
  pdf(file.i)
  
  plot(I(sim_f_100$PI_hat*100) ~ I(PI*100), 
       ylab = expression( hat(pi)), xlab = expression(pi), pch = 20, cex = cex_size, 
       ylim = c(20,40), xlim = c(20,40)   )
  points(I(sim_f_1000$PI_hat*100) ~  I(PI*100), pch = 15, col = 2,cex = cex_size)
  points(I(sim_f_10000$PI_hat*100) ~  I(PI*100), pch = 18, col = 3,cex = cex_size)
  abline(a=0,b = 1,lty = 2)
  legend("topleft",c("T = 100","T = 1,000","T = 10,000"), col = 1:3, pch = c(20,15,18), cex = cex_size)
  grid(10)  
  dev.off()
}

{
  file.i <- paste(output.dir.fig,"EM_ETA.pdf",sep = "")
  pdf(file.i)
  plot(I(sim_f_100$ETA_hat*100^2) ~ I(ETA*100^2), 
       ylab = expression( hat(eta)), xlab = expression(eta), pch = 20,cex = cex_size )
  points(I(sim_f_1000$ETA_hat*100^2) ~  I(ETA*100^2), pch = 15, col = 2, cex = cex_size)
  points(I(sim_f_10000$ETA_hat*100^2) ~  I(ETA*100^2), pch = 18, col = 3, cex = cex_size)
  abline(a=0,b = 1,lty = 2)
  legend("topleft",c("T = 100","T = 1,000","T = 10,000"), col = 1:3, pch = c(20,15,18), cex = cex_size)
  grid(10)
  dev.off()
}

{
  file.i <- paste(output.dir.fig,"EM_SIGMA.pdf",sep = "")
  pdf(file.i)
  plot(I(sim_f_100$SIGMA_hat*100^2) ~ I(SIGMA*100^2), 
       ylab = expression( hat(Sigma)), xlab = expression(Sigma), pch = 20, cex = cex_size )
  points(I(sim_f_1000$SIGMA_hat*100^2) ~  I(SIGMA*100^2), pch = 15, col = 2,cex = cex_size)
  points(I(sim_f_10000$SIGMA_hat*100^2) ~  I(SIGMA*100^2), pch = 18, col = 3,cex = cex_size)
  abline(a=0,b = 1,lty = 2)
  legend("topleft",c("T = 100","T = 1,000","T = 10,000"), col = 1:3, pch = c(20,15,18), cex = cex_size)
  grid(10)
  dev.off()
}


{
  file.i <- paste(output.dir.fig,"EM_GAMMA.pdf",sep = "")
  pdf(file.i)
  plot(I(sim_f_100$GAMMA_tilde_hat*100^2) ~ I(GAMMA_tilde*100^2), 
       ylab = expression( hat(tilde(gamma))), xlab = expression(tilde(gamma)), pch = 20, cex = cex_size, 
       ylim = (100^2)*range(c(sim_f_100$GAMMA_tilde_hat,sim_f_1000$GAMMA_tilde_hat,sim_f_10000$GAMMA_tilde_hat))   )
  points(I(sim_f_1000$GAMMA_tilde_hat*100^2) ~  I(GAMMA_tilde*100^2), pch = 15, col = 2,cex = cex_size)
  points(I(sim_f_10000$GAMMA_tilde_hat*100^2) ~  I(GAMMA_tilde*100^2), pch = 18, col = 3,cex = cex_size)
  abline(a=0,b = 1,lty = 2)
  legend("topleft",c("T = 100","T = 1,000","T = 10,000"), col = 1:3, pch = c(20,15,18), cex = cex_size)
  grid(10)  
  dev.off()
  
}



sims <- 10^3
N <- 10^3

gc()
norm_diff <- function(x) {
  which.norm <- "2"
  a1 <-  norm(matrix(x$ETA_hat  - ETA)/norm(ETA,which.norm),which.norm)^2
  a2 <- norm(matrix(x$SIGMA_hat  - SIGMA)/norm(SIGMA,which.norm),which.norm)^2
  a3 <-  norm(matrix(x$GAMMA_tilde_hat  - GAMMA_tilde)/norm(GAMMA_tilde,which.norm),which.norm)^2
  a4 <-  norm(matrix(x$PI_hat  - PI)/norm(matrix(PI),which.norm),which.norm)
  return(c(a1,a2,a3,a4))
}


NORM_ds <- t(sapply(sim_f_list, norm_diff))
NORM_ds <- data.frame(NORM_ds)
names(NORM_ds) <- c("eta","Sigma","gamma","pi")
stargazer(NORM_ds,digits = 3)

NORM_ds <- lapply(1:ncol(NORM_ds), function(i) data.frame(MSE = NORM_ds[,i],   Parameter = names(NORM_ds)[i]   )  )
NORM_ds <- Reduce(rbind,NORM_ds)
NORM_ds$Parameter <- as.factor(NORM_ds$Parameter)

my.labs <- list(bquote(eta),bquote(Sigma),bquote(gamma),bquote(pi))

p <- ggplot(NORM_ds, aes(x= Parameter, y= MSE),labels = my.labs) + 
  geom_boxplot()
p <- p + ylab(expression(Delta) ) + xlab("")
p <- p + coord_flip(expand = T)
p <- p + theme_gray()


file.i <- paste(output.dir.fig,"EM_MSE.pdf",sep = "")
pdf(file.i)
print(p)
dev.off()




############################################################################################################

#############################################
############ PORTFOLIO FUNCTIONS ############
#############################################

# the following are optimization functions for decision rules


OPT_PORT_function <- function(ETA,SIGMA,PI,GAMMA_tilde,kappa,which.BC) { 
  
  d <- length(ETA)
  U <- function(X) {
    u1 <- t(X)%*%(ETA)
    u2 <- t(X)%*%SIGMA%*%X
    u3 <- log(1 - PI + PI*exp(-kappa*sum(X*GAMMA_tilde)) )
    total <- u1 - (kappa/2)*u2 - (1/kappa)*u3
    return(c(total))
  }
  
  
  G <- function(X) {
    g1 <- ETA
    g2 <- kappa*SIGMA%*%X
    gamma_p <- t(GAMMA_tilde)%*%X
    k <- c((PI*exp( -kappa*gamma_p))/(1 - PI + PI*exp( -kappa*gamma_p  )   ))
    g3 <- as.matrix(k*GAMMA_tilde)
    total <- g1 - g2 + g3
    return(total)
  }
  
  # add constraints
  BC_f <- function(d) {
    # sum to one constraint
    A <- matrix(1,1,d)
    A <- rbind(A,-A)
    B <- c(0.999,-1.001)
    
    # short-sales constraints
    A2 <- diag(rep(1,d))
    B2 <- rep(0,d)
    A2 <- rbind(A,A2)
    B2 <- c(B,B2)
    
    # stack altogether in a list
    BC1 <- list(A,B)
    BC2 <- list(A2,B2)
    list(BC1,BC2)
  }
  
  BC <- BC_f(d)[[which.BC]] # 2 for no short-sales and 1 for yes
  A <- BC[[1]]
  B <- BC[[2]]
  
  # solving this should give the minimum variance portfolio (GMV)
  X0 <- rep(1/d,d)
  X_opt <- constrOptim(X0,function(X) -U(X),grad = function(X) -G(X) ,ui = A,ci = B)
  X1 <- X_opt$par
  X1 <- X1/sum(X1)
  
  
  list(X_opt = X_opt, U = U(X1))
}

# portfolio 2 ignores higher-order moments
OPT_PORT2_function <- function(ETA,SIGMA,PI,GAMMA_tilde,kappa,which.BC) { 
  d <- length(ETA)
  U <- function(X) {
    u1 <- t(X)%*%(ETA)
    u2 <- t(X)%*%SIGMA%*%X
    total <- u1/u2
    
    return(c(total))
  }
  
  G <- function(X) {
    M <- ETA
    S <- SIGMA
    total <-  ((1/c(t(X)%*%S%*%X)^2)*( c(t(X)%*%S%*%X)*M - 2*c(t(X)%*%M)*S%*%X  ))
    return(total)
  }
  
  # add constraints
  BC_f <- function(d) {
    # sum to one constraint
    A <- matrix(1,1,d)
    A <- rbind(A,-A)
    B <- c(0.999,-1.001)
    
    # short-sales constraints
    A2 <- diag(rep(1,d))
    B2 <- rep(0,d)
    A2 <- rbind(A,A2)
    B2 <- c(B,B2)
    
    # stack altogether in a list
    BC1 <- list(A,B)
    BC2 <- list(A2,B2)
    list(BC1,BC2)
  }
  
  BC <- BC_f(d)[[which.BC]] # 2 for no short-sales and 1 for yes
  A <- BC[[1]]
  B <- BC[[2]]
  
  # solving this should give the minimum variance portfolio (GMV)
  X0 <- rep(1/d,d)
  X_opt <- constrOptim(X0,function(X) -U(X),grad = function(X) -G(X) ,ui = A,ci = B)
  X1 <- X_opt$par
  X1 <- X1/sum(X1)
  
  list(X_opt = X_opt, U = U(X1))
}

### add MV as benchmark portfolio
OPT_PORT_MV_function <- function(R_sub,kappa,which.BC) { 
  
  MU <- apply(R_sub, 2, mean)
  SIGMA <- var(R_sub)
  
  d <- length(MU)
  
  U <- function(X) {
    u1 <- t(X)%*%(MU)
    u2 <- t(X)%*%SIGMA%*%X
    total <- u1 - (kappa/2)*u2
    return(c(total))
  }
  
  G <- function(X) {
    g1 <- MU
    g2 <- kappa*SIGMA%*%X
    total <- g1 - g2 
    return(total)
  }
  
  # add constraints
  BC_f <- function(d) {
    # sum to one constraint
    A <- matrix(1,1,d)
    A <- rbind(A,-A)
    B <- c(0.999,-1.001)
    
    # short-sales constraints
    A2 <- diag(rep(1,d))
    B2 <- rep(0,d)
    A2 <- rbind(A,A2)
    B2 <- c(B,B2)
    
    # stack altogether in a list
    BC1 <- list(A,B)
    BC2 <- list(A2,B2)
    list(BC1,BC2)
  }
  
  BC <- BC_f(d)[[which.BC]] # 2 for no short-sales and 1 for yes
  A <- BC[[1]]
  B <- BC[[2]]
  
  # solving this should give the minimum variance portfolio (GMV)
  X0 <- rep(1/d,d)
  X_opt <- constrOptim(X0,function(X) -U(X),grad = function(X) -G(X) ,ui = A,ci = B)
  X1 <- X_opt$par
  X1 <- X1/sum(X1)
  
  list(X_opt = X_opt, U = U(X1) )
  
}

# this corresponds to decision rule 3

OPT_PORT2_MV_function <- function(R_sub,kappa,which.BC) { 
  
  MU <- apply(R_sub, 2, mean)
  SIGMA <- var(R_sub)
  
  d <- length(MU)
  
  U <- function(X) {
    u1 <- t(X)%*%(MU)
    u2 <- t(X)%*%SIGMA%*%X
    total <- u1/u2
    return(c(total))
  }
  
  G <- function(X) {
    M <- MU
    S <- SIGMA
    total <-  ((1/c(t(X)%*%S%*%X)^2)*( c(t(X)%*%S%*%X)*M - 2*c(t(X)%*%M)*S%*%X  ))
    return(total)
  }
  
  
  # add constraints
  BC_f <- function(d) {
    # sum to one constraint
    A <- matrix(1,1,d)
    A <- rbind(A,-A)
    B <- c(0.999,-1.001)
    
    # short-sales constraints
    A2 <- diag(rep(1,d))
    B2 <- rep(0,d)
    A2 <- rbind(A,A2)
    B2 <- c(B,B2)
    
    # stack altogether in a list
    BC1 <- list(A,B)
    BC2 <- list(A2,B2)
    list(BC1,BC2)
  }
  
  BC <- BC_f(d)[[which.BC]] # 2 for no short-sales and 1 for yes
  A <- BC[[1]]
  B <- BC[[2]]
  
  # solving this should give the minimum variance portfolio (GMV)
  X0 <- rep(1/d,d)
  X_opt <- constrOptim(X0,function(X) -U(X),grad = function(X) -G(X) ,ui = A,ci = B)
  X1 <- X_opt$par
  X1 <- X1/sum(X1)
  
  list(X_opt = X_opt,  U = U(X1) )
  
}

# this corresponds to decision rule 5

OPT_GMV_function <- function(SIGMA,which.BC) { 
  
  d <- nrow(SIGMA)
  U <- function(X) {
    u2 <- t(X)%*%SIGMA%*%X
    return(c(u2))
  }
  
  
  G <- function(X) {
    g2 <- SIGMA%*%X
    total <- g2
    return(total)
  }
  
  # add constraints
  BC_f <- function(d) {
    # sum to one constraint
    A <- matrix(1,1,d)
    A <- rbind(A,-A)
    B <- c(0.999,-1.001)
    
    # short-sales constraints
    A2 <- diag(rep(1,d))
    B2 <- rep(0,d)
    A2 <- rbind(A,A2)
    B2 <- c(B,B2)
    
    # stack altogether in a list
    BC1 <- list(A,B)
    BC2 <- list(A2,B2)
    list(BC1,BC2)
  }
  
  BC <- BC_f(d)[[which.BC]] # 2 for no short-sales and 1 for yes
  A <- BC[[1]]
  B <- BC[[2]]
  
  # solving this should give the minimum variance portfolio (GMV)
  X0 <- rep(1/d,d)
  X_opt <- constrOptim(X0,function(X) U(X),grad = function(X) G(X) ,ui = A,ci = B)
  X1 <- X_opt$par
  X1 <- X1/sum(X1)
  
  
  list(X_opt = X_opt, U = U(X1))
}




########################################################################################################################################


###################################################
############# PERFORMANCE METRICS #################
###################################################

my_sum <- function(x) {
  m1 <- mean(x)*12
  s1 <- sd(x)*sqrt(12)
  s2 <- sd(x[x<0])*sqrt(12)
  sr1 <- m1/s1
  sr2 <- m1/s2
  skew <- skewness(x)
  kurt <- kurtosis(x)
  VaR_1 <- m1 - quantile(x,0.05)
  CVAR <- mean(x[x < quantile(x,0.05)])
  result <- c(Mean = m1*100, Std = s1*100, Std_neg = s2*100, Sharpe = sr1, Sortino = sr2, Skewness = skew, Kurtosis = kurt)
  return(result)
}

my_sum_d <- function(x) {
  m1 <- mean(x)*252
  s1 <- sd(x)*sqrt(252)
  s2 <- sd(x[x<0])*sqrt(252)
  sr1 <- m1/s1
  sr2 <- m1/s2
  skew <- skewness(x)
  kurt <- kurtosis(x)
  VaR_1 <- m1 - quantile(x,0.05)
  CVAR <- mean(x[x < quantile(x,0.05)])
  result <- c(Mean = m1, Std = s1, Std_neg = s2, Sharpe = sr1, Sortino = sr2, Skewness = skew, Kurtosis = kurt)
  return(result)
}



########################################################################################################################################

###################################################
############# ESITMATION ON ROLLING WINDOW  #######
###################################################

select_period1 <- "2005-01-01/2010-12-31"
select_period2 <- "2011-01-01/2015-12-31"
select_period_full <- "1990-01-01/2020-12-31"

# xts_all determines the data of interest
# 1 stands for CRSP
# 2 stands for Fama-French BTM

xts_all <- all_data_list[[2]] # <---------------------------- IMPORTANT TO CHOOSE WHICH DATA

ks_f <- function(x) {
  ks_x <- ks.test(x, "pnorm", mean=mean(x), sd=sd(x))
  return(ks_x$statistic)
}

run_variation_i <- function(choose_group) {
  
  cat(" ##################### THIS VARIATION ",choose_group, " ###########","\n")
  
  R <- xts_all
  total_months <- unique(ceiling_date(date(R),"m") - 1)
  d <- ncol(R)
  start_month <- "1994-12-31"
  month_seq <- total_months[total_months >= start_month]
  i_seq <- 1:(length(month_seq)-1)
  
  
  GAMMA_tilde_mat <- PI_mat <- post_prob_mat <- PORT_RET_mat <- R_next_mat <- data.frame()
  X_MVS_mat <- X_MV_mat <- X_MVS2_mat <- X_MV2_mat <- X_MVS3_mat <- X_MV3_mat <-   X_naive_mat <- data.frame()
  
  
  for (i in i_seq) {
    
    m_i <- month_seq[i]
    m_i_next <- month_seq[i+1]
    cat("This is month ",as.character(m_i), "\n")
    
    # choose sub data
    R_sub <- tail(R[date(R) <= m_i,],window_size)
    ks_seq <- apply(R_sub,2, ks_f)
    order_stocks  <- order(ks_seq)
    R_order <- R_sub[,order_stocks]
    # split into 10 groups
    split_group <- ncol(R)/10
    group_list <- lapply(0:9,function(k)  1:split_group + split_group*k )
    
    # split stocks based on 
    R_list <- lapply(group_list,function(x) R_order[,x] )
    R_sub <- R_list[[choose_group]]
    d <- ncol(R_sub)
    R_choose <- R_sub 
    
    R_next <- R[date(R) <= m_i_next & date(R) > m_i,names(R_sub)]
    
    { # estimation
      try_theta <- try(THETA <- theta_f(R_sub,10,10),silent = T) # try theta. if error keep using the previous one
      if(inherits(try_theta,"try-error")) {
        cat("THETA ERROR AT ", i, " MONTH: ", as.character(m_i),"\n")
      }
      
      PI_i <- THETA$PI
      SIGMA_i <- THETA$SIGMA
      GAMMA_tilde_i <- THETA$GAMMA_tilde
      
      # stack data in a data.frame
      GAMMA_tilde_mat <- rbind(GAMMA_tilde_mat,data.frame(date = m_i,t(GAMMA_tilde_i)) )
      PI_mat <- rbind(PI_mat,data.frame(date = m_i,PI = PI_i) )
      
      post_prob <- THETA$post_prob
      rownames(post_prob) <- as.character(date(R_sub))
      post_prob <- as.xts(post_prob)
      keep_last_month <- floor_date(date(post_prob),"m") + months(1) - 1
      keep_last_month <- keep_last_month == m_i
      post_prob <- post_prob[keep_last_month,]
      post_prob <- data.frame(date = date(post_prob),post_prob)
      rownames(post_prob) <- NULL
      post_prob_mat <- rbind(post_prob_mat,post_prob)
      
      ETA_i <- THETA$ETA
      ETA_ALPHA <- ETA_i + GAMMA_tilde_i
    }
    
    kappa_seq <- exp(seq(1,10,length = 10^2))
    port_list <- mclapply(kappa_seq, function(k) try(OPT_PORT_function(ETA_i,SIGMA_i,PI_i,GAMMA_tilde_i,k,which.BC),silent = T),mc.cores = 10 )
    port_list2 <- mclapply(kappa_seq, function(k) try(OPT_PORT_MV_function(R_sub,k,which.BC),silent = T),mc.cores = 10 )
    
    keep_list1 <- !sapply(port_list, function(x) inherits(x,what = "try-error") )
    keep_list2 <- !sapply(port_list2, function(x) inherits(x,what = "try-error") )
    
    keep_list <- keep_list1&keep_list2
    
    port_list <- port_list[keep_list]
    port_list2 <- port_list2[keep_list]
    kappa_seq <- kappa_seq[keep_list]
    
    X_max1 <- sapply(port_list, function(x) max(x$X_opt$par)  )
    X_max2 <- sapply(port_list2, function(x) max(x$X_opt$par)  )
    
    R_seq1 <- lapply(port_list, function(x) R_sub %*% x$X_opt$par  )
    R_seq2 <- lapply(port_list2, function(x) R_sub %*% x$X_opt$par  )
    
    Util_approx <- function(choose_R_seq, k_i) {
      X <- choose_R_seq[[k_i]]
      kappa_i <- kappa_seq[k_i]
      mu_seq <- mean(X)
      sig_seq <- var(X)
      sr_seq <- mu_seq/sig_seq
      skew_seq <- skewness(X)
      kurt_seq <- kurtosis(X)
      result <- mu_seq- 0.5*kappa_i*sig_seq  + ((kappa_i^2)/6)*skew_seq 
      result <- mu_seq- 0.5*kappa_i*sig_seq  + ((kappa_i^2)/6)*skew_seq - ((kappa_i^3)/24)*kurt_seq
      return(result)
    }
    
    Util_approx1 <- sapply(1:length(kappa_seq), function(k_i) Util_approx(R_seq1,k_i)  )
    Util_approx2 <- sapply(1:length(kappa_seq), function(k_i) Util_approx(R_seq2,k_i)  )
    U_diff <- (Util_approx1 - Util_approx2)
    
    U_diff_roll <- na.omit(rollapply(U_diff,3,mean,align = "right",fill = NA))
    kappa_seq_roll <- na.omit(rollapply(kappa_seq,3,mean,align = "right",fill = NA))
    choose_kappa <- kappa_seq_roll[which.max(U_diff_roll)]
    
    port_list <- mclapply(choose_kappa, function(k) OPT_PORT_function(ETA_i,SIGMA_i,PI_i,GAMMA_tilde_i,k,which.BC),mc.cores = 1 )
    port_list2 <- mclapply(choose_kappa, function(k) OPT_PORT_MV_function(R_sub,k,which.BC),mc.cores = 1 )
    
    # consider the GMV portfolio A
    port_list5 <- mclapply(max(kappa_seq), function(k) OPT_GMV_function(SIGMA_i,which.BC),mc.cores = 1 )
    port_list6 <- mclapply(max(kappa_seq), function(k) OPT_GMV_function(var(R_sub),which.BC),mc.cores = 1 )
    
    # consider the portfolio without the coskewness
    port_list3 <- mclapply(1, function(k) OPT_PORT2_function(ETA_i,SIGMA_i,PI_i,GAMMA_tilde_i,k,which.BC),mc.cores = 10 )
    port_list4 <- mclapply(1, function(k) OPT_PORT2_MV_function(R_sub,k,which.BC),mc.cores = 10 )
    
    SR_approx <- function(R_seq) {
      mu_seq <- apply(R_seq,2,mean)
      sig_seq <-  apply(R_seq,2,var)
      sr_seq <- mu_seq/sig_seq
      return(sr_seq)
    }
    
    # # consider 4 portfolio rules
    X_MVS <- port_list[[1]]$X_opt$par
    X_MVS <- X_MVS/sum(X_MVS)
    X_MV <- port_list2[[1]]$X_opt$par
    X_MV <- X_MV/sum(X_MV)
    # consdier those with max A
    X_MVS3 <- port_list5[[1]]$X_opt$par
    X_MVS3 <- X_MVS3/sum(X_MVS3)
    X_MV3 <- port_list6[[1]]$X_opt$par
    X_MV3 <- X_MV3/sum(X_MV3)
    
    # cosnider thoe without skewness that return max SR
    X_MVS2 <- port_list3[[1]]$X_opt$par
    X_MVS2 <- X_MVS2/sum(X_MVS2)
    X_MV2 <- port_list4[[1]]$X_opt$par
    X_MV2 <- X_MV2/sum(X_MV2)
    
    # finally consider the naive portfolio
    X_naive <- rep(1/d,d)
    
    names(X_MVS) <- names(X_MV) <- names(X_MVS2) <- names(X_MV2) <- names(X_MVS3) <- names(X_MV3) <- names(X_naive) <- names(R_sub)
    
    X_MVS_mat <- rbind.fill(X_MVS_mat,data.frame(t(X_MVS)) )
    X_MV_mat <- rbind.fill(X_MV_mat,data.frame(t(X_MV)))
    
    X_MVS3_mat <- rbind.fill(X_MVS3_mat,data.frame(t(X_MVS3)) )
    X_MV3_mat <- rbind.fill(X_MV3_mat,data.frame(t(X_MV3)))
    
    X_MVS2_mat <- rbind.fill(X_MVS2_mat,data.frame(t(X_MVS2)) )
    X_MV2_mat <- rbind.fill(X_MV2_mat,data.frame(t(X_MV2)) )
    
    X_naive_mat <- rbind.fill(X_naive_mat,data.frame(t(X_naive)) )
    
    X_list <- list(MVS = X_MVS_mat,MV = X_MV_mat,  X_MVS3_mat = X_MVS3_mat,X_MV3_mat = X_MV3_mat,X_MVS2_mat = X_MVS2_mat,X_MV2_mat = X_MV2_mat, Naive = X_naive_mat)
    
    # return_next
    R_next <- apply(R_next,2,function(x) prod(x + 1) - 1)
    R_next_mat <- rbind(R_next_mat,R_next)
    
    port_ret_MVS <- R_next%*%X_MVS # our proposed utility
    port_ret_MV <- R_next%*%X_MV # conventional MV
    
    port_ret_MVS3 <- R_next%*%X_MVS3 # max A
    port_ret_MV3 <- R_next%*%X_MV3 # max A which corresponds to GMV
    
    port_ret_MVS2 <- R_next%*%X_MVS2 # without skewness
    port_ret_MV2 <- R_next%*%X_MV2 # counterpart to the previous
    
    port_ret_naive <- R_next%*%X_naive
    
    PORT_RET <- data.frame(date = m_i_next , MVS = port_ret_MVS,  MV = port_ret_MV,
                           MVS3 = port_ret_MVS3, MV3 = port_ret_MV3, 
                           MVS2 = port_ret_MVS2, MV2 = port_ret_MV2, 
                           Naive =  port_ret_naive)
    
    PORT_RET_mat <- rbind(PORT_RET_mat,PORT_RET)
    
    
  }
  
  save_list <- list(GAMMA_tilde_mat = GAMMA_tilde_mat, PI_mat = PI_mat, post_prob_mat = post_prob_mat,PORT_RET_mat = PORT_RET_mat, X_list = X_list, R_next_mat = R_next_mat)
  return(save_list)
  
}


which.BC <- 2
window_size <- 3*252

out_of_sample_results <- mclapply(1:10,run_variation_i,mc.cores = 10)


#######################################################################################################


#### summary of results and bootstrap

get.asterisk <- function(pv) {
  star <- rep("",length(pv))
  star[pv < 0.1] <- "*"
  star[pv < 0.05] <- "**"
  star[pv < 0.01] <- "***"
  return(star)
}

out_of_sample_results_i <- out_of_sample_results[[1]]

main_sum_perf_function <- function(out_of_sample_results_i,select_period,select_port) {
  
  PORT_RET_mat <- out_of_sample_results_i$PORT_RET_mat
  rownames(PORT_RET_mat) <- PORT_RET_mat$date
  PORT_RET_mat$date <- NULL
  PORT_RET_mat <- as.xts(PORT_RET_mat)
  
  all_dates <- date(PORT_RET_mat) 
  
  PORT_RET_mat <- PORT_RET_mat[select_period,c(select_port,7)]
  
  date_index <- all_dates %in% date(PORT_RET_mat)
  
  port_r_boot <- function(port_r_i) { 
    
    boot.f <- function(port_r_i) {
      A2 <- apply(port_r_i,2,my_sum)
      A2 <- cbind(A2,A2[,1] - A2[,2],A2[,1] - A2[,3])
      return( A2[,4:5]  )
    }
    
    boot.metric <- tsboot(port_r_i,statistic =  boot.f,R = 10^3, l = 10, sim = "fixed")
    Z <- t(apply( boot.metric[[2]],1,function(v) v - apply(boot.metric[[2]],2,mean)  ) )
    p.v.A2 <- sapply(1:length(boot.metric[[1]]), function(i)  2*(1-ecdf(Z[,i])(abs(boot.metric[[1]][i] )))    )
    p.v.A2 <- get.asterisk(p.v.A2)
    p.v.A2 <- matrix(p.v.A2,7)
    
    A3 <- apply(port_r_i,2,my_sum)
    A3 <- cbind(A3,A3[,1] - A3[,2],A3[,1] - A3[,3])
    
    A3 <- apply(A3,2, function(x) sprintf("%.2f", round(x,3)) )
    A3[,4:5] <- paste(A3[,4:5],p.v.A2,sep = "")
    return(A3)
  }
  
  M1 <- port_r_boot(PORT_RET_mat)
  result <- M1
  
  ## repeat the same for TO
  X_list <- out_of_sample_results_i$X_list
  X_list <- X_list[c(select_port,7)]
  
  TO_boot <- function(X_list) { 
    
    TO_function <- function(X_i) {
      X_i[is.na(X_i)] <- 0
      TO <- abs(X_i[-1,] - X_i[-nrow(X_i),])
      TO <- (apply(TO, 1,sum))
      return(TO)
    }
    
    TO <- sapply(X_list, TO_function)
    TO <- TO[date_index[-1],]
    
    boot.f <- function(TO) {
      TO_mean <- apply(TO, 2, mean)
      TO_mean <- c(TO_mean,TO_mean[1] - TO_mean[2],TO_mean[1] - TO_mean[3])
      return( TO_mean[4:5]  )
    }
    
    boot.metric <- tsboot(TO,statistic =  boot.f,R = 10^3, l = 10, sim = "fixed")
    Z <- t(apply( boot.metric[[2]],1,function(v) v - apply(boot.metric[[2]],2,mean)  ) )
    p.v.A2 <- sapply(1:length(boot.metric[[1]]), function(i)  2*(1-ecdf(Z[,i])(abs(boot.metric[[1]][i] )))    )
    p.v.A2 <- get.asterisk(p.v.A2)
    
    A3 <- apply(TO,2,mean)
    A3 <- c(A3,A3[1] - A3[2],A3[1] - A3[3])
    
    A3 <- sprintf("%.2f", round(A3,3))
    
    A3[4:5] <- paste(A3[4:5],p.v.A2,sep = "")
    
    return(A3)
  }
  
  M2 <- TO_boot(X_list)
  result <- rbind(result,M2)
  
  return(result)
}

main_sum_perf_diff_function <- function(out_of_sample_results_i,out_of_sample_results_j,select_period,select_port) {
  
  PORT_RET_mat <- out_of_sample_results_i$PORT_RET_mat
  rownames(PORT_RET_mat) <- PORT_RET_mat$date
  PORT_RET_mat$date <- NULL
  PORT_RET_mat <- as.xts(PORT_RET_mat)
  
  all_dates <- date(PORT_RET_mat) 
  
  PORT_RET_mat <- PORT_RET_mat[select_period,c(select_port,7)]
  
  date_index <- all_dates %in% date(PORT_RET_mat)
  
  PORT_RET_mat2 <- out_of_sample_results_j$PORT_RET_mat
  rownames(PORT_RET_mat2) <- PORT_RET_mat2$date
  PORT_RET_mat2$date <- NULL
  PORT_RET_mat2 <- as.xts(PORT_RET_mat2)
  PORT_RET_mat2 <- PORT_RET_mat2[select_period,c(select_port,7)]
  
  PORT_RET_mat12 <- merge(PORT_RET_mat,PORT_RET_mat2)
  
  port_r_boot <- function(port_r_i) { 
    
    boot.f <- function(port_r_i) {
      A1 <- port_r_i[,1:3]
      A1 <- apply(A1,2,my_sum)
      A1 <- cbind(A1,A1[,1] - A1[,2],A1[,1] - A1[,3])
      
      A2 <- port_r_i[,4:6]
      A2 <- apply(A2,2,my_sum)
      A2 <- cbind(A2,A2[,1] - A2[,2],A2[,1] - A2[,3])
      
      A12 <- A2 - A1
      
      return( A12[,4:5]  )
    }
    
    boot.metric <- tsboot(port_r_i,statistic =  boot.f,R = 10^3, l = 10, sim = "fixed")
    Z <- t(apply( boot.metric[[2]],1,function(v) v - apply(boot.metric[[2]],2,mean)  ) )
    p.v.A2 <- sapply(1:length(boot.metric[[1]]), function(i)  2*(1-ecdf(Z[,i])(abs(boot.metric[[1]][i] )))    )
    p.v.A2 <- get.asterisk(p.v.A2)
    p.v.A2 <- matrix(p.v.A2,7)
    
    A3_1 <- apply(port_r_i[,1:3],2,my_sum)
    A3_1 <- cbind(A3_1,A3_1[,1] - A3_1[,2],A3_1[,1] - A3_1[,3])
    
    A3_2 <- apply(port_r_i[,4:6],2,my_sum)
    A3_2 <- cbind(A3_2,A3_2[,1] - A3_2[,2],A3_2[,1] - A3_2[,3])
    
    A3 <- A3_2 - A3_1
    A3 <- apply(A3,2, function(x) sprintf("%.2f", round(x,3)) )
    A3[,4:5] <- paste(A3[,4:5],p.v.A2,sep = "")
    
    return(A3)
  }
  
  M1 <- port_r_boot(PORT_RET_mat12)
  result <- M1
  
  ## repeat the same for TO
  X_list <- out_of_sample_results_i$X_list
  X_list <- X_list[c(select_port,7)]
  
  X_list2 <- out_of_sample_results_j$X_list
  X_list2 <- X_list2[c(select_port,7)]
  
  X_list12 <- c(X_list,X_list2)
  
  TO_function <- function(X_i) {
    X_i[is.na(X_i)] <- 0
    TO <- abs(X_i[-1,] - X_i[-nrow(X_i),])
    TO <- (apply(TO, 1,sum))
    return(TO)
  }
  
  TO_matrix <- sapply(X_list12, TO_function)
  TO_matrix <- TO_matrix[date_index[-1],]
  
  TO_boot <- function(TO_matrix) { 
    
    boot.f <- function(TO_matrix) {
      TO_mean <- apply(TO_matrix[,1:3], 2, mean)
      TO_mean <- c(TO_mean,TO_mean[1] - TO_mean[2],TO_mean[1] - TO_mean[3])
      
      TO_mean2 <- apply(TO_matrix[,4:6], 2, mean)
      TO_mean2 <- c(TO_mean2,TO_mean2[1] - TO_mean2[2],TO_mean2[1] - TO_mean2[3])
      
      return( TO_mean2[4:5] -TO_mean[4:5] )
    }
    
    # there seems to be too many outliers for Theta_m_n
    
    boot.metric <- tsboot(TO_matrix,statistic =  boot.f,R = 10^3, l = 10, sim = "fixed")
    Z <- t(apply( boot.metric[[2]],1,function(v) v - apply(boot.metric[[2]],2,mean)  ) )
    p.v.A2 <- sapply(1:length(boot.metric[[1]]), function(i)  2*(1-ecdf(Z[,i])(abs(boot.metric[[1]][i] )))    )
    p.v.A2 <- get.asterisk(p.v.A2)
    p.v.A2
    
    
    TO_matrix <- sapply(X_list12, TO_function)
    TO_matrix <- TO_matrix[date_index[-1],]
    
    A3 <- apply(TO_matrix,2,mean)
    
    A3_1 <- A3[1:3]
    A3_1 <- c(A3_1,A3_1[1] - A3_1[2],A3_1[1] - A3_1[3])
    
    A3_2 <- A3[4:6]
    A3_2 <- c(A3_2,A3_2[1] - A3_2[2],A3_2[1] - A3_2[3])
    
    
    A3_12 <- A3_2 - A3_1
    A3_12 <- sprintf("%.2f", round(A3_12,3))
    
    A3_12[4:5] <- paste(A3_12[4:5],p.v.A2,sep = "")
    
    # need to add transaction cost to the above table
    return(A3_12)
  }
  
  M2 <- TO_boot(TO_matrix)
  result <- rbind(result,M2)
  return(result)
}

choose_port <- c(1,2)

for(i in 1:10) {
  cat("------------> this is ",i,"\n")
  P_i <- out_of_sample_results[[i]]
  M1 <- main_sum_perf_function(P_i,select_period_full,choose_port)
  M2 <- main_sum_perf_function(P_i,select_period1,choose_port)
  M3 <- main_sum_perf_function(P_i,select_period2,choose_port)
  
  group_i_sum <- cbind(M1,M2,M3)
  rownames(group_i_sum) <-   c("Mean", "Std", "Std_neg", "Sharpe", "Sortino", "Skewness" , "Kurtosis","TO")
  group_i_sum <- strsplit(print(xtable(group_i_sum)),"\n")[[1]][9:16]
  file.i <- paste( output.dir.tables,"table_out_sample_g_",which.BC,"_",i,".txt",sep ="")
  cat(group_i_sum,file = file.i, sep = "\n")
}

i <- 0
{
  
  P1 <- out_of_sample_results[[1]]
  P10 <- out_of_sample_results[[10]]
  
  M1 <- main_sum_perf_diff_function(P1,P10,select_period_full,choose_port)
  M2 <- main_sum_perf_diff_function(P1,P10,select_period1,choose_port)
  M3 <- main_sum_perf_diff_function(P1,P10,select_period2,choose_port)
  
  group_i_sum <- cbind(M1,M2,M3)
  rownames(group_i_sum) <-   c("Mean", "Std", "Std_neg", "Sharpe", "Sortino", "Skewness" , "Kurtosis","TO")
  group_i_sum <- strsplit(print(xtable(group_i_sum)),"\n")[[1]][9:16]
  file.i <- paste( output.dir.tables,"table_out_sample_g_",which.BC,"_",i,".txt",sep ="")
  cat(group_i_sum,file = file.i, sep = "\n")
}


##### GMV portfolio 
choose_port <- c(3,4)

for(i in 1:10) {
  cat("------------> this is ",i,"\n")
  
  P_i <- out_of_sample_results[[i]]
  M1 <- main_sum_perf_function(P_i,select_period_full,choose_port)
  M2 <- main_sum_perf_function(P_i,select_period1,choose_port)
  M3 <- main_sum_perf_function(P_i,select_period2,choose_port)
  
  group_i_sum <- cbind(M1,M2,M3)
  rownames(group_i_sum) <-   c("Mean", "Std", "Std_neg", "Sharpe", "Sortino", "Skewness" , "Kurtosis","TO")
  group_i_sum <- strsplit(print(xtable(group_i_sum)),"\n")[[1]][9:16]
  file.i <- paste( output.dir.tables,"table_out_sample_g_gmv_",which.BC,"_",i,".txt",sep ="")
  cat(group_i_sum,file = file.i, sep = "\n")
}

i <- 0
{
  
  P1 <- out_of_sample_results[[1]]
  P10 <- out_of_sample_results[[10]]
  
  M1 <- main_sum_perf_diff_function(P1,P10,select_period_full,choose_port)
  M2 <- main_sum_perf_diff_function(P1,P10,select_period1,choose_port)
  M3 <- main_sum_perf_diff_function(P1,P10,select_period2,choose_port)
  
  group_i_sum <- cbind(M1,M2,M3)
  rownames(group_i_sum) <-   c("Mean", "Std", "Std_neg", "Sharpe", "Sortino", "Skewness" , "Kurtosis","TO")
  group_i_sum <- strsplit(print(xtable(group_i_sum)),"\n")[[1]][9:16]
  file.i <- paste( output.dir.tables,"table_out_sample_g_gmv_",which.BC,"_",i,".txt",sep ="")
  cat(group_i_sum,file = file.i, sep = "\n")
}

### report the shrinkage portfolio
choose_port <- c(5,6)

for(i in 1:10) {
  cat("------------> this is ",i,"\n")
  
  P_i <- out_of_sample_results[[i]]
  M1 <- main_sum_perf_function(P_i,select_period_full,choose_port)
  M2 <- main_sum_perf_function(P_i,select_period1,choose_port)
  M3 <- main_sum_perf_function(P_i,select_period2,choose_port)
  
  group_i_sum <- cbind(M1,M2,M3)
  rownames(group_i_sum) <-   c("Mean", "Std", "Std_neg", "Sharpe", "Sortino", "Skewness" , "Kurtosis","TO")
  group_i_sum <- strsplit(print(xtable(group_i_sum)),"\n")[[1]][9:16]
  file.i <- paste( output.dir.tables,"table_out_sample_g_shrink_",which.BC,"_",i,".txt",sep ="")
  cat(group_i_sum,file = file.i, sep = "\n")
}

i <- 0
{
  
  P1 <- out_of_sample_results[[1]]
  P10 <- out_of_sample_results[[10]]
  
  M1 <- main_sum_perf_diff_function(P1,P10,select_period_full,choose_port)
  M2 <- main_sum_perf_diff_function(P1,P10,select_period1,choose_port)
  M3 <- main_sum_perf_diff_function(P1,P10,select_period2,choose_port)
  
  group_i_sum <- cbind(M1,M2,M3)
  rownames(group_i_sum) <-   c("Mean", "Std", "Std_neg", "Sharpe", "Sortino", "Skewness" , "Kurtosis","TO")
  group_i_sum <- strsplit(print(xtable(group_i_sum)),"\n")[[1]][9:16]
  file.i <- paste( output.dir.tables,"table_out_sample_g_shrink_",which.BC,"_",i,".txt",sep ="")
  cat(group_i_sum,file = file.i, sep = "\n")
  
  
}



###########################################################3
### summary statistics
Sk <- data.frame(t(apply(xts_all, 2, my_sum_d)))
Sk <- Sk[,c("Mean","Std","Skewness",  "Kurtosis")]
stargazer(Sk)
Mean <- apply(xts_all, 2, mean)*252
Std <- apply(xts_all, 2, sd)*sqrt(252)



###############################################
#### plot latent factor exposure over time ####
bar.col <- "lightblue"
ggplot_recession1 <- geom_rect(fill = bar.col,col = bar.col,aes(xmin=date("1990-07-31"), xmax=date("1991-03-31"), ymin=-Inf, ymax=Inf))
ggplot_recession2 <- geom_rect(fill = bar.col,col = bar.col,aes(xmin=date("2001-03-31"), xmax=date("2001-11-30"), ymin=-Inf, ymax=Inf))
ggplot_recession3 <- geom_rect(fill = bar.col,col = bar.col,aes(xmin=date("2007-12-31"), xmax=date("2009-06-30"), ymin=-Inf, ymax=Inf))

for(g in c(1,2,5,9,10)) {
  ds_plot_g <- ds_plot_all[ds_plot_all$Group  == g,]
  
  p <- ggplot(ds_plot_g, aes(Date, Median),colour = Group)
  p <- p + ggplot_recession2 +ggplot_recession3
  p <- p + geom_line(data=ds_plot_g)+
    geom_ribbon(data=ds_plot_g,aes(ymin=P10,ymax=P90),alpha=0.3)
  p <- p + geom_hline(yintercept = 0, linetype="dashed")
  p <- p + ylim(-2.5,2.5) + ylab("")
  
  file.i <- paste(output.dir.fig,"gamma_tilde_median_",g,".pdf",sep = "")
  pdf(file.i)
  print(p)
  dev.off()
}






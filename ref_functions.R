# Functions
# Moran's I estimated ( ~ Pearson correlation)
#########################################################
# function to convert images to unit square
#########################################################
image_convert = function(X_coordinate, Y_coordinate){
  h_min = min(Y_coordinate)
  h_max = max(Y_coordinate)
  l_min = min(X_coordinate)
  l_max = max(X_coordinate)
  L = max(l_max - l_min, h_max - h_min)
  
  x_new = (X_coordinate - l_min)/L
  y_new = (Y_coordinate - h_min)/L
  return(list("X_new" = x_new, "Y_new" = y_new))
}

#########################################################
# fix r_max at 1/4 for now
# to make sure all patients have Moran's I correlation evaluated at the same distance
moran_cross = function(X, Y, m, d, r){
  # X, Y: marked subpatterns
  # index of marks, (if there's multiple)
  #calculate cross distance between X and Y
  D12 = crossdist(X, Y)
  # d: divider of max distance
  #r_max = #1/4#max(D12)/d
  #breaks = handle.r.b.args(r = NULL, breaks = NULL, window, pixeps = NULL, r_max)
  #r_seq = breaks$r
  r_seq = r
  mi = marks(X)[m]
  mj = marks(Y)[m]
  moran_res = c()
  for (k in 1: length(r_seq)){
    idx_pair = which(D12 < r_seq[k], arr.ind = TRUE)
    if (!is.empty(idx_pair)){
      mi_k = mi[idx_pair[,1],]
      mj_k = mj[idx_pair[,2],]
      mi_bar = mean(mi_k)#mean(mi[unique(idx_pair[,1]),])
      mj_bar = mean(mj_k)#mean(mj[unique(idx_pair[,2]),])
      num = crossprod(mi_k - mi_bar, mj_k - mj_bar)
      m_bar = mean(mi_bar, mj_bar)
      deno = crossprod(c(mi_k,mj_k) - m_bar, c(mi_k, mj_k) - m_bar)
      moran_res[k] = num/(deno+ 0.0001) #num/deno
    }
    else {moran_res[k] = 0}
  }
  return(list("r" = r_seq, "moran_est" = moran_res))
}



moran_cross1 = function(X, Y, m, d, r){
  # X, Y: marked subpatterns
  # index of marks, (if there's multiple)
  #calculate cross distance between X and Y
  D12 = crossdist(X, Y)
  # d: divider of max distance
  #r_max = #1/4#max(D12)/d
  #breaks = handle.r.b.args(r = NULL, breaks = NULL, window, pixeps = NULL, r_max)
  #r_seq = breaks$r
  r_seq = r
  mi = marks(X)[m]
  mj = marks(Y)[m]
  moran_res = c()
  for (k in 1: length(r_seq)){
    idx_pair = which(D12 <= r_seq[k], arr.ind = TRUE)
    if (!is.empty(idx_pair)){
      mi_k = mi[idx_pair[,1],]
      mj_k = mj[idx_pair[,2],]
      mi_bar = mean(mi[unique(idx_pair[,1]),])
      mj_bar = mean(mj[unique(idx_pair[,2]),])
      num = crossprod(mi_k - mi_bar, mj_k - mj_bar)
      m_bar = mean(mi_bar, mj_bar)
      deno = crossprod(c(mi_k,mj_k) - m_bar, c(mi_k, mj_k) - m_bar)
      moran_res[k] = num/deno
    }
    else {moran_res[k] = 0}
  }
  return(list("r" = r_seq, "moran_est" = moran_res))
}

# cal_grid for LFCM and AFCM models

## construct the grid (100*100 matrix by default) with specified ranges
### function "cal_grid()" returns grid (grid_coef) with pre-specified range of x and s
#### grid_length: length of grid in each direction, 100 by default
#### range_x: pre-specified upper bound of x, 1 by default
#### range_s: pre-specified upper bound of s, 1 by default
cal_grid_old <- function(grid_length = 100, range_x = 1, range_s = 1){
  tind_pred <- seq(0, range_s, length.out = grid_length)
  xind_pred <- seq(-1, range_x, length.out = grid_length)
  grid_coef <- data.frame(t_int = rep(tind_pred, grid_length),
                          func_cov = rep(xind_pred, each = grid_length),
                          l_int=1)
  return(grid_coef)
}

#grid generating
#cal_grid()

## construct the grid (100*100 matrix by default) with specified ranges
### function "cal_grid()" returns grid (grid_coef) with pre-specified range of x and s
#### grid_length: length of grid in each direction, 100 by default
#### range_x: pre-specified upper bound of x, 1 by default
#### range_s: pre-specified upper bound of s, 1 by default
cal_grid <- function(grid_length = 100, range_x = c(0, 0.25), range_s = c(0,295)){
  tind_pred <- seq(range_s[1], range_s[2], length.out = grid_length)
  xind_pred <- seq(range_x[1], range_x[2], length.out = grid_length)
  grid_coef <- data.frame(t_int = rep(tind_pred, grid_length),
                          func_cov = rep(xind_pred, each = grid_length),
                          l_int=1)
  return(grid_coef)
}

# following Figure 5 in (Cui et al. 2021)
### function "vis_surf()" returns the visualization of functional surface using ggplot2 package
#### coef: the matrix to visualize
#### simulation: indicate whether the plot is for application (FALSE) or simulation (TRUE)
#### transformation: indicate whether the functional covariates are transformed (TRUE) or not (FALSE)
#### range_x: the upper bound of functional value domain
#### range_s: the upper bound of functional domain
#### trans_name: name of the transformation
vis_surf_old <- function(coef, simulation = FALSE, transformation = TRUE, range_x = 1, range_s = 1, 
                     trans_name = "smoothed LAC"){
  coef_plot <- melt(coef) # re-organize to plot form
  grid_length <- nrow(coef) ## length of grid in each direction
  ## when the surface is from real data, x axis is time of a day, y axis is transformed AC
  if(simulation == FALSE){
    if(transformation == TRUE){
      name_y <- expression(h[is](log(1 + AC)))
    }else{
      name_y <- expression(log(1 + AC))
    }
    surf <- ggplot(coef_plot, aes(Var1,Var2, fill = value)) + 
      geom_raster() +
      scale_fill_gradientn(colours=c("blue3", "#FFFFFF", "#F67F7F", "#FF0000"), name = expression(hat(F)(.,.))) +
      ggtitle(trans_name) +
      scale_x_continuous(name = "Time of day",
                         breaks = c(1, floor(0.25 * grid_length), floor(0.5 * grid_length), floor(0.75 * grid_length), grid_length),
                         labels = c("00:00","06:00", "12:00", "18:00", "24:00")) +
      scale_y_continuous(name = name_y, 
                         breaks = floor(seq(0, floor(range_x), length.out = 5) / range_x * grid_length), 
                         labels = seq(0, floor(range_x), length.out = 5)) +
      theme_minimal() + 
      theme(plot.title = element_text(face = "bold", size = 15, hjust = 0.5))
  }else{ ## when the surface is from simulated data, x axis is the functional domain, y axis is the functional value domain,
    surf <- ggplot(coef_plot, aes(Var1,Var2, fill = value)) + 
      geom_raster() +
      scale_fill_gradientn(colours = c("blue3", "#FFFFFF", "red2"), name = expression(hat(F)(.,.))) +
      ggtitle(trans_name) +
      scale_x_continuous(name = "Distance (r)",
                         breaks = floor(seq(0, range_s, length.out = 6) / range_s * grid_length), 
                         labels = seq(0, range_s, length.out = 6)) +
      scale_y_continuous(name = "Moran's I cross correlation",
                         #breaks = floor(seq(0, floor(range_x), length.out = 6) / range_x * grid_length), 
                         breaks = floor(seq(0, floor(range_x), length.out = 6) / range_x * grid_length), 
                         labels = seq(-1, floor(range_x), length.out = 6)) +
      theme_minimal() + 
      theme(plot.title = element_text(face = "bold", size = 15, hjust = 0.5))
  }
  return(surf)
}


vis_surf <- function(coef, simulation = FALSE, transformation = TRUE, range_x = c(0,1), range_s = c(0,1), 
                     trans_name = "smoothed LAC", x_name = "r (microns)", y_name = "mcf" ){
  coef_plot <- melt(coef) # re-organize to plot form
  grid_length <- nrow(coef) ## length of grid in each direction
  ## when the surface is from real data, x axis is time of a day, y axis is transformed AC
  if(simulation == FALSE){
    if(transformation == TRUE){
      name_y <- expression(h[is](log(1 + AC)))
    }else{
      name_y <- expression(log(1 + AC))
    }
    surf <- ggplot(coef_plot, aes(Var1,Var2, fill = value)) + 
      geom_raster() +
      scale_fill_gradientn(colours=c("blue3", "#FFFFFF", "#F67F7F", "#FF0000"), name = expression(hat(F)(.,.))) +
      ggtitle(trans_name) +
      scale_x_continuous(name = "Time of day",
                         breaks = c(1, floor(0.25 * grid_length), floor(0.5 * grid_length), floor(0.75 * grid_length), grid_length),
                         labels = c("00:00","06:00", "12:00", "18:00", "24:00")) +
      scale_y_continuous(name = name_y, 
                         breaks = floor(seq(0, floor(range_x), length.out = 5) / range_x * grid_length), 
                         labels = seq(0, floor(range_x), length.out = 5)) +
      theme_minimal() + 
      theme(plot.title = element_text(face = "bold", size = 15, hjust = 0.5))
  }else{ ## when the surface is from simulated data, x axis is the functional domain, y axis is the functional value domain,
    surf <- ggplot(coef_plot, aes(Var1,Var2, fill = value)) + 
      geom_raster() +
      scale_fill_gradientn(colours = c("blue3", "#FFFFFF", "red2"), name = expression(hat(F)(.,.))) +
      ggtitle(trans_name) +
      scale_x_continuous(name = x_name,
                         breaks = floor(seq(0, grid_length, length.out = 6)) , 
                         labels = floor(seq(range_s[1], range_s[2], length.out = 6))) +
      scale_y_continuous(name = y_name,
                         breaks = floor(seq(0, range_x[2], length.out = 6) / range_x[2] * grid_length), 
                         labels = round(seq(range_x[1], range_x[2], length.out = 6),2))  +
      #geom_hline(yintercept=1, linetype="dashed", color = "black")+
      theme_minimal() + 
      theme(plot.title = element_text(face = "bold", size = 15, hjust = 0.5))
  }
  return(surf)
}

#vis_surf_new

vis_surf_new <- function(coef, x_lo = -1, x_up = 1, s_lo = 0, s_up = 2, name, fontsize = 12){
  coef_plot <- melt(coef) # re-organize to plot form
  grid_length <- nrow(coef)
  
  color <- c("blue3", "#FFFFFF", "red2")
  surf <- ggplot(coef_plot, aes(Var1,Var2, fill = value)) + 
    geom_raster() +
    scale_fill_gradientn(colours=color, name = expression(hat(F)(.,.))) +
    ggtitle(name) +
    scale_x_continuous(name = "Functional domain",
                       breaks = floor(seq(0, grid_length, length.out = 6)), 
                       labels = floor(seq(s_lo, s_up, length.out = 6))) +
    scale_y_continuous(name = "Value of functional covariate",
                       breaks = floor(seq(0, grid_length, length.out = 6)), 
                       labels = round(seq(x_lo, x_up, length.out = 6), 2)) +
    #geom_hline(yintercept=1, linetype="dashed", color = "black")+
    theme_minimal() +
    theme(plot.title = element_text(size = fontsize, face = "bold", hjust = 0.5))
  
  return(surf)
}

#===============================================================================================#
# SIMULATION
#===============================================================================================#

## simulate functional covariate from real data
### function "sim_func_cov" simulate functional covariates using FPCA
#### real_func_cov: the matrix to simulate functional covariates from 
#### N: number of subjects to simulate
sim_func_cov <- function(real_func_cov, N){
  ## FPCA
  fit_fpca <- fpca.face(unclass(real_func_cov))
  
  ## extract estimated eigenfunctions/eigenvalues
  phi <- fit_fpca$efunctions
  lambda <- fit_fpca$evalues
  K <- length(fit_fpca$evalues)
  
  ## simulate random coefficients 
  theta <- matrix(rnorm(N*K), nrow=N, ncol=K)  # generate independent standard normals 
  theta <- theta %*% diag(sqrt(lambda))        # scale to have appropriate variance
  X_new <- theta %*% t(phi)                    # simulate new functions
  X_new <- X_new + t(fit_fpca$mu %o% rep(1,N))      # add back in the mean function
  return(X_new)
}



### integrate F
cal_ISE_F <- function(sim_coef, real_coef){
  res = sum((sim_coef - real_coef)^2) / (ncol(real_coef) * nrow(real_coef))
  return(res)
}



## simulate survival data with pre-specified F(.,.)
### function "sim_surv_F()" simulate survival data using pre-specified F(.,.)
#### fit: the model used to estimate cumulative hazard
#### N: number of subjects to simulate
#### tmax: upper bound of time points to simulate from
#### time_mort: time to event from real data
#### event: event status from real data
#### nt_pred: number of potential survival times
#### model: name of the specified F(.,.)
#### range_x: range of functional covariate value
#### range_s: range of functional domain
sim_surv_F1 <- function(fit, N, func_data, scalar_var, tmax, time_mort, event, nt_pred = 1000, model = "linear", range_x = c(0,1), range_s = c(0,1)){
  ## calculate estimated baseline hazard
  t0 <- rev(fit$family$data$tr)
  H0_hat <- rev(fit$family$data$h) ## estimate for the baseline cumulative hazard function \int_0^t \lambda_0(s)ds
  H0_fit <- scam(H0_hat ~ s(t0, bs = "mpi")-1)
  ## evaluate smoothed H0 on fine grid
  tgrid_sim <- seq(0, tmax, len=nt_pred) ## potential grid of time points to simulate from
  H0_prd <- pmax(0, predict(H0_fit, newdata = data.frame(t0 = tgrid_sim))) ## truncate cumulative hazard at 0
  H0_prd[1] <- 0
  
  ## define pre-specified functional coefficient surface F(.,.), labeled as beta
  if(model == "linear"){
    beta <- function(X, s) X*s
  }else if(model == "const"){
    beta <- function(X, s) X*2
  }else if(model == "nonlinear1"){
    beta <- function(X, s) (X^3)*s
  }else if(model == "nonlinear2"){
    beta <- function(X, s) sin(X*s)
  }
  
  ## simulate scalar covariate (eg, age)
  X_scalar = rnorm(N, mean(scalar_var), sd(scalar_var)/1) #/2
  #X_scalar <- scale(X_scalar, center = T, scale = F)
  
  ## simulate functional covariate from FPCA
  X_new <- sim_func_cov(func_data, N)
  #X_new <- I(apply(X_new, 2,function(y) ecdf(y)(y))) * range_x[2] * 2 #this is for quantile transformation, we don't need that?
  #X_new <- sweep(X_new, 2, colMeans(X_new)) ## mean center simulated functional covariates
  
  ## combine functional covariates with necessary elements in model fitting
  nt <- ncol(X_new)
  tind <- seq(range_s[1], range_s[2], length.out = nt)
  data_sim <- data.frame(func_cov = I(X_new), 
                         l_int = I(matrix(1/nt, ncol = nt, nrow = N)),
                         t_int = I(matrix(tind, ncol = nt, nrow = N, byrow = TRUE)))
  
  # calculate eta_i using numerical intergration given pre-specified F
  # eta_i should also consider scalar covariate (i.e age) 
  # scalar linear predictor using the coefficient from coxph model with only age
  scalar_lp = X_scalar * fit$coefficients["Age"] #(X_scalar - fit_scalar$means) * coef(fit_scalar)
  func_lp =  as.matrix(apply(X_new, 1, function(x) sum(beta(X=x, s=tind))/nt), ncol = 1)
  eta_i <-   scalar_lp + func_lp
  
  ## calcualte survival times
  S_i <- exp(-(exp(eta_i) %*% H0_prd))
  U   <- runif(N)
  Ti  <- rep(NA, N)
  for(i in 1:N){
    if(all(S_i[i,] > U[i])){
      Ti[i] <- max(tgrid_sim) + 1
    } else {
      Ti[i] <- tgrid_sim[min(which(S_i[i,] < U[i]))]
    }
  }
  
  ## simulate censoring times from empirical distribution
  C_i <- sample(time_mort[event == 0], size=N, replace=TRUE)
  
  ## Create simulated censoring indicator based on "true" event time S_i and censoring time C_i
  Ti_obs <- pmin(C_i, Ti)
  di <- as.numeric(Ti <= C_i)
  
  ## combine the simulated survival data in a data frame
  data_sim <- data.frame(data_sim, Age = X_scalar, time_mort = Ti_obs, event=di, lp = eta_i)
  #data_sim$act_log_mat_ave <- rowSums(data_sim$func_cov) / ncol(data_sim$func_cov)
  
  #return(list("data_sim" = data_sim, "lambda0" = H0_prd))
  return(data_sim)
}




sim_surv_F2 <- function(fit, N, func_data, scalar_var, tmax, time_mort, event, nt_pred = 1000, model = "linear", range_x = c(0,1), range_s = c(0,1)){
  ## calculate estimated baseline hazard
  t0 <- rev(fit$family$data$tr)
  H0_hat <- rev(fit$family$data$h) ## estimate for the baseline cumulative hazard function \int_0^t \lambda_0(s)ds
  H0_fit <- scam(H0_hat ~ s(t0, bs = "mpi")-1)
  ## evaluate smoothed H0 on fine grid
  tgrid_sim <- seq(0, tmax, len=nt_pred) ## potential grid of time points to simulate from
  H0_prd <- pmax(0, predict(H0_fit, newdata = data.frame(t0 = tgrid_sim))) ## truncate cumulative hazard at 0
  H0_prd[1] <- 0
  
  ## define pre-specified functional coefficient surface F(.,.), labeled as beta
  if(model == "linear"){
    beta <- function(X, s) X*s
  }else if(model == "const"){
    beta <- function(X, s) X*2
  }else if(model == "nonlinear1"){
    beta <- function(X, s) (X^2)*s
  }else if(model == "nonlinear2"){
    beta <- function(X, s) sin(X*s)
  }
  
  ## simulate scalar covariate (eg, age)
  X_scalar = rnorm(N, mean(scalar_var), sd(scalar_var)/1) #/2
  #X_scalar <- scale(X_scalar, center = T, scale = F)
  
  ## simulate functional covariate from FPCA
  X_new <- sim_func_cov(func_data, N)
  #X_new <- I(apply(X_new, 2,function(y) ecdf(y)(y))) * range_x[2] * 2 #this is for quantile transformation, we don't need that?
  #X_new <- sweep(X_new, 2, colMeans(X_new)) ## mean center simulated functional covariates
  
  ## combine functional covariates with necessary elements in model fitting
  nt <- ncol(X_new)
  tind <- seq(range_s[1], range_s[2], length.out = nt)
  data_sim <- data.frame(func_cov = I(X_new), 
                         l_int = I(matrix(1/nt, ncol = nt, nrow = N)),
                         t_int = I(matrix(tind, ncol = nt, nrow = N, byrow = TRUE)))
  
  # calculate eta_i using numerical intergration given pre-specified F
  # eta_i should also consider scalar covariate (i.e age) 
  # scalar linear predictor using the coefficient from coxph model with only age
  scalar_lp = X_scalar  #(X_scalar - fit_scalar$means) * coef(fit_scalar)
  func_lp =  as.matrix(apply(X_new, 1, function(x) sum(beta(X=x, s=tind))/nt), ncol = 1)
  eta_i <-   scalar_lp + func_lp
  
  ## calcualte survival times
  S_i <- exp(-(exp(eta_i) %*% H0_prd))
  U   <- runif(N)
  Ti  <- rep(NA, N)
  for(i in 1:N){
    if(all(S_i[i,] > U[i])){
      Ti[i] <- max(tgrid_sim) + 1
    } else {
      Ti[i] <- tgrid_sim[min(which(S_i[i,] < U[i]))]
    }
  }
  
  ## simulate censoring times from empirical distribution
  C_i <- sample(time_mort[event == 0], size=N, replace=TRUE)
  
  ## Create simulated censoring indicator based on "true" event time S_i and censoring time C_i
  Ti_obs <- pmin(C_i, Ti)
  di <- as.numeric(Ti <= C_i)
  
  ## combine the simulated survival data in a data frame
  data_sim <- data.frame(data_sim, scal_avg = X_scalar, time_mort = Ti_obs, event=di, lp = eta_i)
  #data_sim$act_log_mat_ave <- rowSums(data_sim$func_cov) / ncol(data_sim$func_cov)
  
  #return(list("data_sim" = data_sim, "lambda0" = H0_prd))
  return(data_sim)
}


sim_surv_F3 <- function(fit, N, func_data, scalar_var, tmax, time_mort, event, nt_pred = 1000, model = "linear", range_x = c(0,1), range_s = c(0,1)){
  ## calculate estimated baseline hazard
  t0 <- rev(fit$family$data$tr)
  H0_hat <- rev(fit$family$data$h) ## estimate for the baseline cumulative hazard function \int_0^t \lambda_0(s)ds
  H0_fit <- scam(H0_hat ~ s(t0, bs = "mpi")-1)
  ## evaluate smoothed H0 on fine grid
  tgrid_sim <- seq(0, tmax, len=nt_pred) ## potential grid of time points to simulate from
  H0_prd <- pmax(0, predict(H0_fit, newdata = data.frame(t0 = tgrid_sim))) ## truncate cumulative hazard at 0
  H0_prd[1] <- 0
  
  ## define pre-specified functional coefficient surface F(.,.), labeled as beta
  if(model == "linear"){
    beta <- function(X, s) X*s
  }else if(model == "const"){
    beta <- function(X, s) X*2
  }else if(model == "nonlinear1"){
    beta <- function(X, s) (X^3)*s
  }else if(model == "nonlinear2"){
    beta <- function(X, s) sin(X*s)
  }
  
  ## simulate scalar covariate (eg, age)
  X_scalar = rnorm(N, mean(scalar_var), sd(scalar_var)/1) #/2
  X_scalar <- scale(X_scalar, center = T, scale = T)
  
  ## simulate functional covariate from FPCA
  X_new <- sim_func_cov(func_data, N)
  #X_new <- I(apply(X_new, 2,function(y) ecdf(y)(y))) * range_x[2] * 2 #this is for quantile transformation, we don't need that?
  X_new = sweep(X_new, 2, colMeans(X_new))/apply(X_new, 2, sd)
  
  ## combine functional covariates with necessary elements in model fitting
  nt <- ncol(X_new)
  tind <- seq(range_s[1], range_s[2], length.out = nt)
  data_sim <- data.frame(func_cov = I(X_new), 
                         l_int = I(matrix(1/nt, ncol = nt, nrow = N)),
                         t_int = I(matrix(tind, ncol = nt, nrow = N, byrow = TRUE)))
  
  # calculate eta_i using numerical intergration given pre-specified F
  # eta_i should also consider scalar covariate (i.e age) 
  # scalar linear predictor using the coefficient from coxph model with only age
  scalar_lp = X_scalar # * fit$coefficients["Age"]
  func_lp =  as.matrix(apply(X_new, 1, function(x) sum(beta(X=x, s=tind))/nt), ncol = 1)
  eta_i <-   scalar_lp + func_lp
  
  ## calcualte survival times
  S_i <- exp(-(exp(eta_i) %*% H0_prd))
  U   <- runif(N)
  Ti  <- rep(NA, N)
  for(i in 1:N){
    if(all(S_i[i,] > U[i])){
      Ti[i] <- max(tgrid_sim) + 1
    } else {
      Ti[i] <- tgrid_sim[min(which(S_i[i,] < U[i]))]
    }
  }
  
  ## simulate censoring times from empirical distribution
  C_i <- sample(time_mort[event == 0], size=N, replace=TRUE)
  
  ## Create simulated censoring indicator based on "true" event time S_i and censoring time C_i
  Ti_obs <- pmin(C_i, Ti)
  di <- as.numeric(Ti <= C_i)
  
  ## combine the simulated survival data in a data frame
  data_sim <- data.frame(data_sim, X_scalar = X_scalar, time_mort = Ti_obs, event=di, lp = eta_i)
  #data_sim$act_log_mat_ave <- rowSums(data_sim$func_cov) / ncol(data_sim$func_cov)
  
  #return(list("data_sim" = data_sim, "lambda0" = H0_prd))
  return(data_sim)
}



#' Calculate knot sequence for a B-spline
#' See Eilers & Marx (2010) for details, the code was
#' taken from their paper.
#' @param xl minimum of X
#' @param xr maximum of X
#' @param ndx number of intervals
#' @param bdeg degree of basis function (2 for cubic)
get_knots <- function(xl, xr, ndx, bdeg) {
  dx <- (xr - xl) / ndx
  knots <- seq(xl - bdeg * dx, xr + bdeg * dx, by = dx)
  return(knots)
}

#' @title 
#' Simulate pupil data to test papss
#'
#' @description 
#' Uses the pupil response function described by Hoeks & Levelt (1993)
#' with the refinements proposed by Wierda e t al. (2012) to simulate pupil
#' dilation time-course data with three sources of deviation:
#' \enumerate{
#' \item Simulates systematic per-subject deviation from a 'shared' population trend in demand
#' \item Simulates per-trial deviation from each subject's individual 'true demand'
#' \item Assumes random noise (N(0, sigma), with constant sigma) for each trial
#' }
#' 
#' @details
#' Demand is modelled according to a simple B-spline (see Eilers & Marx, 2010)
#' with equally spaced knots with associated coefficients of which only a small
#' percentage will be different from zero (to ensure that the simulated demand
#' trajectory is sparse).
#' 
#' @param nk Number of knots for the B-spline basis.
#' @param n_sub How many subjects to simulate.
#' @param n_trials How many trials to simulate.
#' @param pulse_loc_diff Assume a pulse every 'pulse_loc_diff' samples (one sample = 20 ms)
#' @param n_diff Maximum number of spline basis coefficients with systematic per-subject variation
#' @param spars_deg The degree of sparsity enforced in spline basis coefficient vector
#' @param sub_dev Standard deviation of normal distribution used to sample by-subject variation
#' @param slope_dev Standard deviation of normal distribution used to sample by-subject slope variation
#' @param trial_dev Standard deviation of normal distribution used to sample by-trial variation (in spike weights and coefficients)
#' @param residual_dev Standard deviation of normal distribution used to sample residuals per trial
#' @param n Parameter defined by Hoeks & Levelt (number of laters)
#' @param t_max Parameter defined by Hoeks & Levelt (response maximum in ms)
#' @param f Parameter defined by Wierda et al. (scaling factor)
#' @param should_plot Whether the generated data should be visualized.
#' @param seed For replicability
#' 
#' @examples 
#' pupil_sim <- additive_pupil_sim(n_sub=5)
#' sim_dat <- pupil_sim$data
#' @export
additive_pupil_sim <- function(nk=20,
                               n_sub=10,
                               n_trials=250,
                               pulse_loc_diff=1,
                               n_diff=4,
                               spars_deg=0.5,
                               sub_dev=0.15,
                               slope_dev = 1.5,
                               trial_dev=0.05,
                               residual_dev=15.0,
                               n=10.1,
                               t_max=930,
                               f=1/(10^24),
                               should_plot=T,
                               seed=124){
  
  # Response function from Hoeks and Levelt
  # + scale parameter introduced by Wierda et al.(2012).
  # Code was taken from their example code.
  h_pupil <- function(t,n,t_max,f)
  {
    h<-f*(t^n)*exp(-n*t/t_max)
    h[0] <- 0
    h
  }
  
  # Optionally set seed
  if(!is.null(seed)){
    set.seed(seed)
  }
  
  # Now create time variable
  X <- seq(0,3000,by=20)
  
  # To simulate demand as a B-spline, we first need to
  # calculate the knot locations over the time variable.
  # Taken from Eilers & Marx (2010)
  
  k <- get_knots(min(X),max(X),(nk-2),2)
  
  # Now calculate the basis functions (here deg + 1 is passed for
  # cubic B-spline basis, see Wood (2017) for details).
  B <- splines::spline.des(k,X,3,0*X,outer.ok = T)$design
  
  # Start with demand on population level
  coef <- runif(ncol(B)) # Get coef for each B-spline basis
  coef[1:2] <- 0 # First two should always be zero so that demand increases 'after stimulus onset'.
  coef[sample(1:ncol(B),round(nk*spars_deg))] <- 0 # Induce sparsity.
  coef[(length(coef) - 4):length(coef)] <- 0 # Let demand decay towards the end.
  
  # Now this is the actual demand change over time for the population.
  # i.e., these are our demand spike weights/coefficients, that we try to recover.
  pop_truth_demand <- B %*% matrix(coef,ncol = 1) 
  
  # Optionally restrict demand spikes to every nth. sample step
  pulse_locations <- seq(1,length(pop_truth_demand),pulse_loc_diff) # Take only every nth. pulse
  real_locations <- X[pulse_locations] # Get corresponding time points
  pop_truth_demand <- pop_truth_demand[pulse_locations] # Reflect this in demand
  
  # Embed demand again in vector matching X/time length.
  pop_spikes <- rep(0,length(X))
  pop_spikes[pulse_locations] <- pop_truth_demand
  
  if(should_plot){
    plot(X,pop_spikes,type = "l",lwd=3,
         xlab="Time",
         ylab="Demand",
         main="Demand over time on the population level")
  }
  
  # Get predicted pupil dilation time-course on population level.
  # Based on code from Wierda et al. (2012)
  hc <- h_pupil(X,n,t_max,f)
  o <- convolve(pop_spikes,rev(hc),type="open") # 'o' notation follows Wierda et al.
  o <- o[1:length(X)]
  
  if(should_plot){
    plot(X,o,type="l",lwd=3,
         xlab="Time",
         ylab="Pupil dilation",
         main="Pupil dilation time-course on the population level")
  }
  
  ### By-Subject Variation ###
  
  # Select non-negative B-spline weights
  nn_ind_all <- (1:ncol(B))[coef > 0]
  
  # Select weights that should be different for subjects.
  nn_ind <- sample(nn_ind_all,min(n_diff,length(nn_ind_all)))
  
  # Now form sub. deviation matrix
  sub_dev_mat <- matrix(coef,nrow=length(coef),ncol=n_sub)
  for (i in 1:length(nn_ind)) {
    dev <- rnorm(n_sub,sd=sub_dev)
    sub_dev_mat[nn_ind[i],] <- sub_dev_mat[nn_ind[i],] + dev
  }
  
  # Now calculate the true demand weights for each subject.
  sub_truth_demand <- B %*% sub_dev_mat
  
  # Re-enforce non-negativity constrain that might have been violated because of the dev.
  sub_truth_demand[sub_truth_demand < 0] <- 0
  
  # Again optionally restrict demand spikes to every nth. sample step.
  sub_truth_demand <- sub_truth_demand[pulse_locations,]
  
  # Again create matrix that matches actual time dimension
  sub_spikes <- matrix(0,nrow = length(X),ncol = n_sub)
  sub_spikes[pulse_locations,] <- sub_truth_demand
  
  if(should_plot){
    plot(X,pop_spikes,type = "l",lwd=3,
         ylim = c(0,
                  (max(pop_spikes) + 4*sub_dev*max(pop_spikes))),
         xlab="Time",
         ylab="Demand",
         main="Demand over time (with subject variation)")
    
    for(si in 1:n_sub){
      lines(X,sub_spikes[,si],lty=2,col="red")
    }
  }
  
  # Get trend/slope deviation for each subject
  slopes <- rnorm(n_sub,sd=slope_dev)
  
  # Get predicted pupil dilation time-course for each subject
  o_subs <- matrix(0,nrow = length(X),ncol=n_sub)
  
  if(should_plot){
    plot(X,o,type="l",lwd=3,ylim = c((min(o) - 4*sub_dev*max(o)),
                                     (max(o) + 4*sub_dev*max(o))),
         xlab="Time",
         ylab="Pupil dilation",
         main="Pupil dilation time-course (with subject variation)")
  }
  
  for(si in 1:n_sub){
    o_sub <- convolve(sub_spikes[,si],rev(hc),type="open")
    o_sub <- o_sub[1:length(X)]
    o_sub <- o_sub + (slopes[si] * 1:length(X))
    o_subs[,si] <- o_sub
    
    if(should_plot){
      lines(X,o_sub,lty=2,col="red")
    }
  }
  
  # Now simulate trials for each subject
  trial <- rep(rep(1:n_trials,each=length(X)),times=n_sub)
  time <- rep(rep(X,times=n_trials),times=n_sub)
  pupil <- rnorm(length(trial),sd=residual_dev) # Add some white noise on each trial.
  subject <- rep(0,lenght.out=length(trial))
  
  ci <- 1
  for (si in 1:n_sub) {
    
    # Calculate trial deviation for subject (for spline and slope coef)
    trial_dev_mat <- matrix(0,nrow = n_trials,ncol = length(nn_ind_all))
    trial_dev_mat_slope <- matrix((rnorm(n_trials,sd=trial_dev)+slopes[si]),nrow=n_trials,ncol=1)
    
    # Take a new sample for each difference
    for(nni in 1:length(nn_ind_all)){
      dev_t <- rnorm(n_trials,sd=trial_dev)
      trial_dev_mat[,nni] <- dev_t
    }
    
    # Now calculate all trials for subject
    for(ti in 1:n_trials){
      # Now add by-trial deviation to each subject ground-truth.
      sub_trial_coef <- sub_dev_mat[,si]
      sub_trial_coef[nn_ind_all] <- sub_trial_coef[nn_ind_all] + trial_dev_mat[ti,]
      sub_trial_truth <- B %*% sub_trial_coef
      
      # Re-enforce constraints again:
      sub_trial_truth[sub_trial_truth < 0] <- 0
      
      # Optionally restrict demand spikes to every nth. sample step.
      sub_trial_truth <- sub_trial_truth[pulse_locations]
      
      # And finally get spike vector again on time domain.
      sub_trial_spikes <- matrix(0,nrow = length(X),ncol = 1)
      sub_trial_spikes[pulse_locations,] <- sub_trial_truth
      
      # Now get trial-level pupil time-course prediction.
      o_sub_trial <- convolve(sub_trial_spikes,rev(hc),type="open")
      o_sub_trial <- o_sub_trial[1:length(X)]
      o_sub_trial <- o_sub_trial + (trial_dev_mat_slope[ti] * 1:length(X))
      
      # Embed in variables
      pupil[ci:(ci+length(o_sub_trial)-1)] <- pupil[ci:(ci+length(o_sub_trial)-1)] + o_sub_trial
      subject[ci:(ci+length(o_sub_trial)-1)] <- si
      ci <- ci + length(o_sub_trial)
      
    }
  }
  
  # Collect data
  sim_dat <- data.frame(pupil,trial,subject,time)
  sim_dat$subject <- as.factor(sim_dat$subject)
  sim_dat$trial <- as.factor(sim_dat$trial)
  
  if(should_plot){
    plot(0,0,type="n",
         xlim=c(min(X),max(X)),
         ylim=c(min(sim_dat$pupil),max(sim_dat$pupil)),
         xlab="Time",
         ylab="Pupil dilation",
         main="Pupil dilation time-course (with subject & trial variation)")
    
    for(si in 1:n_sub){
      sub_dat <- sim_dat[sim_dat$subject ==si,]
      for(t in unique(sub_dat$trial)){
        lines(sub_dat$time[sub_dat$trial == t],sub_dat$pupil[sub_dat$trial == t],col=si)
      }
    }
  }
  
  # Return simulated data frame and the ground truth objects to allow for comparisons
  return(list("pop_truth_demand"=pop_spikes,
              "pop_truth_pupil"=o,
              "sub_truth_demand"=sub_spikes,
              "sub_truth_pupil"=o_subs,
              "time"=X,
              "data"=sim_dat))
}

#' @title 
#' Plot simulated pupil data against recovered weights.
#' 
#' @description
#' Only works for model setup like the one by Wierda et al. (2012).
#' 
#' @param n_sub How many subjects were simulate.
#' @param aggr_dat Aggregated simulation data
#' @param sim_obj The list returned by additive_pupil_sim
#' @param recovered_coef The coefficients returned by papss:pupil_solve()
#' @param pulse_locations The index vector pointing at pulse locations passed to papss::pupil_solve()
#' @param real_locations The vector with the time-points at which pulses are assumed passed to papss::pupil_solve()
#' @param pulses_in_time Boolean vector indicating which pulse was within the time-window. Returned by papss::pupil_solve()
#' @param expanded_time The expanded time series returned by papss:pupil_solve()
#' @param expanded_by Expansion time in ms passed to papss::pupil_solve(expand_by=) divided by sample length in ms
#' @param f_est f parameter value
#' @param scaling_factor A scaling factor to scale up or down the demand signal
#' @param se NULL or list with standard errors recovered for each subject
#' @param plot_avg Whether to plot average or not
#' @export
plot_sim_vs_recovered <- function(n_sub,
                                  aggr_dat,
                                  sim_obj,
                                  recovered_coef,
                                  pulse_locations,
                                  real_locations,
                                  pulses_in_time,
                                  expanded_time,
                                  expanded_by,
                                  f_est=1/(10^24),
                                  scaling_factor=1,
                                  se=NULL,
                                  plot_avg=T){

  # First re-create model matrix, assuming Wierda et al. (2012) like setup.
  slopePredX <- papss::create_slope_term(unique(aggr_dat$time),1)
  
  semiPredX <- papss::create_spike_matrix_term(unique(expanded_time),
                                               expanded_by,
                                               unique(aggr_dat$time),
                                               pulse_locations,
                                               n=10.1,
                                               t_max=930,
                                               f=f_est)
  # Combine slope and spline terms
  predMat <- cbind(slopePredX,semiPredX)
  
  # First n coefficients are the slopes for each subject
  slopes <- recovered_coef[1:n_sub]
  
  # Placeholder to later collect the population estimate, i.e., a simple average.
  pop_spike_est <- rep(0,length.out=length(unique(aggr_dat$time)))
  
  # Get estimates for each subject
  for(i in 1:n_sub){
    sub_i <- i
    
    # Get corresponding spline spike weights
    splineCoefSub <- recovered_coef[((n_sub + 1) + ((i - 1) * ncol(semiPredX))):(n_sub+(i * ncol(semiPredX)))]
    
    if(!is.null(se)){
      standardErrorSub <- se[((n_sub + 1) + ((i - 1) * ncol(semiPredX))):(n_sub+(i * ncol(semiPredX)))]
    }
    
    # combine all subj. coefficients
    allCoefSub <- c(slopes[i],splineCoefSub)
    subCoefMat <- matrix(0,nrow = length(allCoefSub),ncol=1)
    subCoefMat[,1] <- allCoefSub
    
    # Get prediction
    predSub <- predMat %*% subCoefMat
    
    splineCoefSub <- splineCoefSub * scaling_factor
    
    # Plot
    plot(aggr_dat$time[aggr_dat$subject == sub_i],
         aggr_dat$pupil[aggr_dat$subject == sub_i],
         type="l",
         main=sub_i,
         xlab = "time",
         ylab="Pupil dilation (base-lined)",
         lwd=3)
    lines(unique(aggr_dat$time),predSub,lty=2,col="red",lwd=3)
    
    # Spike plot
    true_peaks_plot <- sim_obj$sub_truth_demand[,sub_i]
    
    
    est_peaks_plot <- rep(0,length.out = length(unique(aggr_dat$time)))
    est_peaks_plot[unique(aggr_dat$time) %in%
                     real_locations[pulses_in_time]] <- splineCoefSub[pulses_in_time]
    
    # Calculate sum for average
    pop_spike_est <- pop_spike_est + est_peaks_plot
    
    plot(unique(aggr_dat$time),true_peaks_plot,type="l",
         ylim = c(0,1.5),
         main = sub_i,
         xlab = "time",
         ylab="Spike strength",lwd=3)
    lines(unique(aggr_dat$time),est_peaks_plot,col="red",lwd=3)
    
    if(!is.null(se)){
      lower_se <- est_peaks_plot
      upper_se <- est_peaks_plot
      
      lower_se[unique(aggr_dat$time) %in%
          real_locations[pulses_in_time]] <- lower_se[unique(aggr_dat$time) %in%
                                                        real_locations[pulses_in_time]] - standardErrorSub[pulses_in_time]
      
      upper_se[unique(aggr_dat$time) %in%
                 real_locations[pulses_in_time]] <- upper_se[unique(aggr_dat$time) %in%
                                                               real_locations[pulses_in_time]] + standardErrorSub[pulses_in_time]
      
      lines(unique(aggr_dat$time),lower_se,col="red",lwd=2,lty=2)
      lines(unique(aggr_dat$time),upper_se,col="red",lwd=2,lty=2)
    }
    
    
  }
  par(mfrow=c(1,1))
  # Calculate average to approximate population estimate
  pop_spike_est <- pop_spike_est/n_sub
  if(plot_avg){
    # Now plot population estimate
    plot(unique(aggr_dat$time),sim_obj$pop_truth_demand,type="l",
         ylim = c(0,1.5),
         main = "Population level",
         xlab = "time",
         ylab="Spike strength",lwd=3)
    lines(unique(aggr_dat$time),pop_spike_est,col="red",lwd=3)
  }
}

#'@title
#' Extracts demand curves
#'
#' @description
#' Extracts demand curves for each level of the factor passed to the pupil_solve
#' call.
#' @param aggr_dat Aggregated data passed to solver
#' @param recovered_coef The coefficients returned by papss::pupil_solve()
#' @param pulse_locations The index vector pointing at pulse locations passed to papss::pupil_solve()
#' @param real_locations The vector with the time-points at which pulses are assumed passed to papss::pupil_solve()
#' @param pulses_in_time Boolean vector indicating which pulse was within the time-window. Returned by papss::pupil_solve()
#' @param expanded_time The expanded time series returned by papss::pupil_solve()
#' @param expanded_by Expansion time in ms passed to papss::pupil_solve(expand_by=) divided by sample length in ms
#' @param factor_id Name of factor column in aggr_dat
#' @param n_fact Choice for n parameter
#' @param t_max_fact Choice for t_max parameter
#' @param f_fact Choice for f parameter
#' @export
extract_demand_for_fact <- function(aggr_dat,
                                    recovered_coef,
                                    pulse_locations,
                                    real_locations,
                                    pulses_in_time,
                                    expanded_time,
                                    expanded_by,
                                    factor_id="subject",
                                    n_fact=10.1,
                                    t_max_fact = 930,
                                    f_fact=1/(10^24)){
  
  factor <- aggr_dat[,colnames(aggr_dat) == factor_id]
  unq_factor <- unique(factor)
  n_fact <- length(unq_factor)
  
  # First re-create demand-related parts of model matrix.
  semiPredX <- papss::create_spike_matrix_term(unique(expanded_time),
                                               expanded_by,
                                               unique(aggr_dat$time),
                                               pulse_locations,
                                               n=n_fact,
                                               t_max=t_max_fact,
                                               f=f_fact)
  
  demand_dat <- NULL
  
  # Get estimates for each factor level
  for(fi in 1:n_fact){
    
    # Get corresponding spline spike weights
    splineCoef <- recovered_coef[((n_fact + 1) + ((fi - 1) * ncol(semiPredX))):(n_fact+(fi * ncol(semiPredX)))]
    
    # Align with original time variable
    demand_trajectory <- rep(0,length.out = length(unique(aggr_dat$time)))
    demand_trajectory[unique(aggr_dat$time) %in%
                     real_locations[pulses_in_time]] <- splineCoef[pulses_in_time]
    
    # Collect demand data
    demand_fact_dat <- data.frame("demand"=demand_trajectory,
                                  "time"=unique(aggr_dat$time),
                                  "factor"=rep(unq_factor[fi],
                                               length.out=length(demand_trajectory)))
    
    demand_dat <- rbind(demand_dat,demand_fact_dat)
  }
  return(demand_dat)
}
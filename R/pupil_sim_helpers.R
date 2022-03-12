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


#' Simulate pupil data to test the functionality offered by this package
#'
#' Uses the pupil response function described by Hoeks & Levelt (1993)
#' with the refinements proposed by Wierda e t al. (2012) to simulate pupil
#' dilation time-course data with three sources of deviation:
#' 
#' 1) Simulates systematic per-subject deviation from a 'shared' population trend in demand
#' 2) Simulates per-trial deviation from each subject's individual 'true demand'
#' 3) Assumes random noise (N(0, sigma), with constant sigma) for each trial
#' 
#' Demand is modelled according to a simple B-spline (see Eilers & Marx, 2010)
#' with equally spaced knots with associated coefficients of which only a small
#' percentage will be different from zero (to ensure that the simulated demand
#' trajectory is sparse).
#' 
#' @param nk Number of knots for the B-spline basis.
#' @param n_subs How many subjects to simulate.
#' @param n_trials How many trials to simulate.
#' @param pulse_loc_diff Assume a pulse every 'pulse_loc_diff' samples (one sample = 20 ms)
#' @param n_diff Maximum number of spline basis coefficients with systematic per-subject variation
#' @param spars_deg The degree of sparsity enforced in spline basis coefficient vector
#' @param sub_dev Standard deviation of normal distribution used to sample by-subject variation
#' @param trial_dev Standard deviation of normal distribution used to sample by-trial variation (in spike weights and coefficients)
#' @param residual_dev Standard deviation of normal distribution used to sample residuals per trial
#' @param n Parameter defined by Hoeks & Levelt (number of laters)
#' @param t_max Parameter defined by Hoeks & Levelt (response maximum in ms)
#' @param f Parameter defined by Wierda et al. (scaling factor)
#' @param should_plot Whether the generated data should be visualized.
#' @param seed For replicability
#' @export
additive_pupil_sim <- function(nk=20,
                               n_subs=10,
                               n_trials=250,
                               pulse_loc_dev=1,
                               n_diff=4,
                               spars_deg=0.5,
                               sub_dev=0.25,
                               trial_dev=0.25,
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
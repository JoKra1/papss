#' Create "basis-functions" that make up the pupil spline.
#' 
#' The basis is calculated using the pupil response function originally
#' described by Hoeks & levelt (1993). Wierda et al. (2012) refined it
#' by adding an additional scaling paramter (f). The code here makes use
#' of the convolution operation between a spike unit vector and said pupil
#' response function and is based on the code provided by Wierda et al (2012)
#' in their supplementary materials.
#' 
#' See: Wierda, S. M., van Rijn, H., Taatgen, N. A., & Martens, S. (2012).
#' Pupil dilation deconvolution reveals the dynamics of attention at high
#' temporal resolution. Proceedings of the National Academy of Sciences of
#' the United States of America, 109(22), 8456–8460.
#' https://doi.org/10.1073/pnas.1201858109
#' 
#' See: Hoeks, B., & Levelt, W. (1993). Pupillary dilation as a measure of
#' attention: A quantitative system analysis.
#' Behav. Res. Meth. Ins. C., 25, 16–26.
#' 
#' @param i An integer representing index of a basis function.
#' @param time A numeric vector containing positive time values in ms
#' @param pulse_locations A numeric vector containing index values of pulse loc.
#' @param n Parameter defined by Hoeks & Levelt (number of laters)
#' @param t_max Parameter defined by Hoeks & Levelt (response maximum in ms)
#' @param f Parameter defined by Wierda et al. (scaling factor)
#' @export
h_basis <- function(i,expanded_time,expand_by,time,pulse_locations,n,t_max,f) {
  
  # We need to make sure that bases towards the end do not 'contaminate'
  # the data for the next subject. Since we assume that each subject
  # contributes the same time-steps we form the basis over the unique
  # expanded time-values and then just repeat the basis multiple times (for each subject).
  unq_time <- unique(expanded_time)
  
  # Response function from Hoeks and Levelt
  # + scale parameter introduced by Wierda et al.
  # n+1 = number of laters
  # t_max = response maximum
  # f = scaling factor
  h<-f*(unq_time^n)*exp(-n*unq_time/t_max)
  
  # Set weight of current basis to 1
  strengths <- rep(0,length(pulse_locations))
  strengths[i] <- 1
  
  # Create spike vector
  peaks <- rep(0,max(pulse_locations)+1)
  peaks[pulse_locations] <- strengths
  
  # Convolve "spike" defined by peaks with h
  o <- convolve(peaks,rev(h),type="open")
  
  # Keep only the realization of the response function
  # within the un-expanded time window
  o_restr <- o[(expand_by + 1):length(unq_time)]
  
  # Now repeat the basis function for each subject (i.e., until the dimension
  # matches the dimension of time).
  o_restr <- rep(o_restr,length.out=length(time))
  
  return(o_restr)
  
}

#' Creates a single intercept term (at the population level).
#' Can be manipulated via "by" argument as done in mgcv.
#' 
#' @param time A numeric vector containing positive time values in ms
#' @export
create_constant_term <- function(time) {
  intercept <- rep(1,length(time))
  return(intercept)
}

#' Creates a single slope term (at the population level)
#' Can be manipulated via "by" argument as done in mgcv.
#' Note: Use of unique values that are just repeated requires all levels
#' of the category factor (usually subjects) to have the same values on the
#' covariate time.
#' 
#' @param time_unq A numeric vector containing unique positive time values in ms
#' @param n_cat Integer, number of categories for which slope should be repeated.
#' @export
create_slope_term <- function(time_unq,n_cat) {
  slope <- 1:length(time_unq)
  slope <- rep(slope,n_cat)
  return(slope)
}

#' Create basis part of the model matrix as described by Wood (2017).
#' Essentially, this means creating a matrix with a column for each
#' pulse location and to assign a pupil basis (convolved with that location)
#' to each column.
#' 
#' See: Wood, S. N. (2017). Generalized Additive Models: An Introduction with R,
#' Second Edition (2nd ed.). Chapman and Hall/CRC.
#' 
#' @param time A numeric vector containing positive time values in ms
#' @param pulse_locations A numeric vector containing index values of pulse loc.
#' @param n Parameter defined by Hoeks & Levelt (number of laters)
#' @param t_max Parameter defined by Hoeks & Levelt (response maximum in ms)
#' @param f Parameter defined by Wierda et al. (scaling factor)
#' @export
create_spike_matrix_term <- function(expanded_time,expand_by,time,pulse_locations,n,t_max,f) {
  
  spike_matrix <- matrix(nrow = length(time),ncol = (length(pulse_locations)))
  
  for (i in 1:length(pulse_locations)) {
    spike_matrix[,i] <- h_basis(i,expanded_time,expand_by,time,pulse_locations,n,t_max,f)
  }
  
  return(spike_matrix)
}

#' Achieves "by" keyword functionality used in mgcv to enable by-factor
#' smooths (see. Wood, 2017).
#' 
#' See: Wood, S. N. (2017). Generalized Additive Models: An Introduction with R,
#' Second Edition (2nd ed.). Chapman and Hall/CRC.
#' 
#' See 'by' at: https://www.rdocumentation.org/packages/mgcv/versions/1.8-38/topics/s
#' 
#' Code is based on: https://stats.stackexchange.com/questions/110472/how-exactly-is-the-sum-or-mean-centering-constraint-for-splines-also-w-r-t-g
#' 
#' @param term A slope, splike_matrix, or intercept.
#' @param fact A factor variable.
#' @export
term_by_factor <- function(term,fact) {
  
  term_by <- NULL
  for(l in unique(fact)) {
    term_by <- cbind(term_by, term * (fact == l))
  }
  return(term_by)
}

### PAPSS model templates ###

#' Creates trainings matrix and penalties for a fully penalized
#' version of the original Wierda et al. (2012) model.
#' Like by Wierda et al. separate sets of coefficients are estimated per
#' subject with the h_basis terms. All those coefficients are again constrained
#' to be non-negative. Model also includes the slope terms for each subject that
#' Wierda et al. (2012) introduced. These are not constrained - since their
#' primary purpose was to account for drift in the pupil averages of individual
#' subjects.
#' 
#' This model differs from the one by Wierda et al. (2012) in that it penalizes
#' the coefficients corresponding to the h_basis terms. Specifically, it
#' enforces a single penalty term shared by all subjects (e.g., this is
#' similar to the 'fs' basis in mgcv or can be achieved by using the 'id' keyword
#' in a 'by' factor smooth.). The form of the penalty expressed on all of the
#' basis functions is a simple identity matrix. We here also penalize all slope
#' terms, again with a single penalty (this time applied to a single matrix though).
#' 
#' @param time A numeric vector containing positive time values in ms
#' @param subs A factor vector containing subject identifiers
#' @param pulse_locations A numeric vector containing index values of pulse loc.
#' @param n Parameter defined by Hoeks & Levelt (number of laters)
#' @param t_max Parameter defined by Hoeks & Levelt (response maximum in ms)
#' @param f Parameter defined by Wierda et al. (scaling factor)
WIER_SHARED_NNLS_model_setup <- function(expanded_time,expand_by,time,subs,pulse_locations,n,t_max,f) {
  
  # Extract number of subjects
  n_subs <- length(unique(subs))
  # Setup model matrix
  slope <- create_slope_term(unique(time),n_subs)
  spike_matrix <- create_spike_matrix_term(expanded_time,expand_by,time,pulse_locations,n,t_max,f)
  slope_matrix <- term_by_factor(slope,subs)
  spike_matrix_by <- term_by_factor(spike_matrix,subs)
  
  trainingsMatrix <- cbind(slope_matrix,spike_matrix_by)
  
  # Setup Penalty definition to be implemented by c++
  # one individual slope penaltiy plus one
  # shared penalty for all subjects (expressed on bases)
  freq <- c(1, n_subs)
  # individual penalty is of size n_subs*n_subs. Also, Each of the shared
  # penalties is of size length(pulse_locations)*length(pulse_locations)
  size <- c(n_subs,length(pulse_locations))
  
  
  # Define positive constraints
  constraints <- rep("c",ncol(trainingsMatrix))
  constraints[1:n_subs] <- "u" # All parametric terms are unconstrained
  
  return(list("X"=trainingsMatrix,
              "Penalties"=list("freq"=freq,
                               "size"=size,
                               "startIndex"=0),
              "constraints"=constraints))
}

#' Creates trainings matrix and penalties for a fully penalized
#' version of the original Wierda et al. (2012) model.
#' Like by Wierda et al. separate sets of coefficients are estimated per
#' subject with the h_basis terms. All those coefficients are again constrained
#' to be non-negative. Model also includes the slope terms for each subject that
#' Wierda et al. (2012) introduced. These are not constrained - since their
#' primary purpose was to account for drift in the pupil averages of individual
#' subjects.
#' 
#' This model differs from the one by Wierda et al. (2012) in that it penalizes
#' the coefficients corresponding to the h_basis terms. Specifically, it
#' enforces a single penalty term shared by all subjects (e.g., this is
#' similar to the 'fs' basis in mgcv or can be achieved by using the 'id' keyword
#' in a 'by' factor smooth.). The form of the penalty expressed on all of the
#' basis functions is a simple identity matrix. We here also penalize all slope
#' terms, again with a single penalty (this time applied to a single matrix though).
#' 
#' @param time A numeric vector containing positive time values in ms
#' @param subs A factor vector containing subject identifiers
#' @param pulse_locations A numeric vector containing index values of pulse loc.
#' @param n Parameter defined by Hoeks & Levelt (number of laters)
#' @param t_max Parameter defined by Hoeks & Levelt (response maximum in ms)
#' @param f Parameter defined by Wierda et al. (scaling factor)
WIER_IND_NNLS_model_setup <- function(expanded_time,expand_by,time,subs,pulse_locations,n,t_max,f) {
  
  # Extract number of subjects
  n_subs <- length(unique(subs))
  # Setup model matrix
  slope <- create_slope_term(unique(time),n_subs)
  spike_matrix <- create_spike_matrix_term(expanded_time,expand_by,time,pulse_locations,n,t_max,f)
  slope_matrix <- term_by_factor(slope,subs)
  spike_matrix_by <- term_by_factor(spike_matrix,subs)
  
  trainingsMatrix <- cbind(slope_matrix,spike_matrix_by)
  
  # Setup Penalty definition to be implemented by c++
  # one individual slope penaltiy plus one
  # shared penalty for all subjects (expressed on bases)
  freq <- c(1, rep(1,length.out=n_subs))
  # individual penalty is of size n_subs*n_subs. Also, Each of the shared
  # penalties is of size length(pulse_locations)*length(pulse_locations)
  size <- c(n_subs,rep(length(pulse_locations),length.out=n_subs))
  
  
  # Define positive constraints
  constraints <- rep("c",ncol(trainingsMatrix))
  constraints[1:n_subs] <- "u" # All parametric terms are unconstrained
  
  return(list("X"=trainingsMatrix,
              "Penalties"=list("freq"=freq,
                               "size"=size,
                               "startIndex"=0),
              "constraints"=constraints))
}

#' Creates trainings matrix and penalties for a fully penalized
#' model inspired by the one used by Denison et al., (2020).
#' Like by Wierda et al. (2012) separate sets of coefficients are estimated per
#' subject with the h_basis terms. All those coefficients are again constrained
#' to be non-negative. Model also includes the intercept terms for each subject that
#' Denison et al. (2012) introduced. These are not constrained - since their
#' primary purpose was to account for pupil trajectories that were overall just
#' containing very negative samples (relative to baseline).
#' 
#' This model differs from the one by Wierda et al. (2012) and Denison et al., (2020)
#' in that it penalizes the coefficients corresponding to the h_basis terms.
#' Specifically, it enforces a single penalty term shared by all subjects (e.g., this is
#' similar to the 'fs' basis in mgcv or can be achieved by using the 'id' keyword
#' in a 'by' factor smooth.). The form of the penalty expressed on all of the
#' basis functions is a simple identity matrix. We here also penalize all intercept
#' terms, again with a single penalty (this time applied to a single matrix though).
#' 
#' @param time A numeric vector containing positive time values in ms
#' @param subs A factor vector containing subject identifiers
#' @param pulse_locations A numeric vector containing index values of pulse loc.
#' @param n Parameter defined by Hoeks & Levelt (number of laters)
#' @param t_max Parameter defined by Hoeks & Levelt (response maximum in ms)
#' @param f Parameter defined by Wierda et al. (scaling factor)
DEN_SHARED_NNLS_model_setup <- function(expanded_time,expand_by,time,subs,pulse_locations,n,t_max,f) {
  
  # Extract number of subjects
  n_subs <- length(unique(subs))
  # Setup model matrix
  intercept <- create_constant_term(time)
  spike_matrix <- create_spike_matrix_term(expanded_time,expand_by,time,pulse_locations,n,t_max,f)
  intercept_matrix <- term_by_factor(intercept,subs)
  spike_matrix_by <- term_by_factor(spike_matrix,subs)
  
  trainingsMatrix <- cbind(intercept_matrix,spike_matrix_by)
  
  # Setup Penalty definition to be implemented by c++
  # one individual intercept penalty plus one
  # shared penalty for all subjects (expressed on bases)
  freq <- c(1, n_subs)
  # individual penalty is of size of size n_subs*n_subs. Also, Each of the shared
  # penalties is of size length(pulse_locations)*length(pulse_locations)
  size <- c(n_subs,length(pulse_locations))
  
  
  # Define positive constraints
  constraints <- rep("c",ncol(trainingsMatrix))
  constraints[1:n_subs] <- "u" # All parametric terms are unconstrained
  
  return(list("X"=trainingsMatrix,
              "Penalties"=list("freq"=freq,
                               "size"=size,
                               "startIndex"=0),
              "constraints"=constraints))
}

#' Creates trainings matrix and penalties for a fully penalized
#' model that combines elements of both the models by Denison et al., (2020)
#' and Wierda et al. (2012).
#'
#' Like by Wierda et al. (2012) separate sets of coefficients are estimated per
#' subject with the h_basis terms. All those coefficients are again constrained
#' to be non-negative. Model also includes the intercept terms for each subject that
#' Denison et al. (2012) introduced as well as the slope term introduced by
#' Wierda et al., (2012). These are not constrained - since their
#' primary purpose was to account for pupil trajectories that were overall just
#' containing very negative samples (relative to baseline) and drifts in the
#' trajectories respectively.
#' 
#' This model differs from the one by Wierda et al. (2012) and Denison et al., (2020)
#' in that it penalizes the coefficients corresponding to the h_basis terms.
#' Specifically, it enforces a single penalty term shared by all subjects (e.g., this is
#' similar to the 'fs' basis in mgcv or can be achieved by using the 'id' keyword
#' in a 'by' factor smooth.). The form of the penalty expressed on all of the
#' basis functions is a simple identity matrix. We here also penalize all intercept
#' and slope terms, again with a single penalty
#' (this time applied to a single matrix though).
#' 
#' @param time A numeric vector containing positive time values in ms
#' @param subs A factor vector containing subject identifiers
#' @param pulse_locations A numeric vector containing index values of pulse loc.
#' @param n Parameter defined by Hoeks & Levelt (number of laters)
#' @param t_max Parameter defined by Hoeks & Levelt (response maximum in ms)
#' @param f Parameter defined by Wierda et al. (scaling factor)
WIER_DEN_SHARED_NNLS_model_setup <- function(expanded_time,expand_by,time,subs,pulse_locations,n,t_max,f) {
  
  # Extract number of subjects
  n_subs <- length(unique(subs))
  # Setup model matrix
  intercept <- create_constant_term(time)
  slope <- create_slope_term(unique(time),n_subs)
  spike_matrix <- create_spike_matrix_term(expanded_time,expand_by,time,pulse_locations,n,t_max,f)
  intercept_matrix <- term_by_factor(intercept,subs)
  slope_matrix <- term_by_factor(slope,subs)
  spike_matrix_by <- term_by_factor(spike_matrix,subs)
  
  trainingsMatrix <- cbind(intercept_matrix,slope_matrix,spike_matrix_by)
  
  # Setup Penalty definition to be implemented by c++
  # one individual intercept  and slope penalty plus one
  # shared penalty for all subjects (expressed on bases)
  freq <- c(1, 1, n_subs)
  # individual penalties are of size n_subs*n_subs. Also, Each of the shared
  # penalties is of size length(pulse_locations)*length(pulse_locations)
  size <- c(n_subs, n_subs, length(pulse_locations))
  
  
  # Define positive constraints
  constraints <- rep("c",ncol(trainingsMatrix))
  constraints[1:(2 * n_subs)] <- "u" # All parametric terms are unconstrained
  
  return(list("X"=trainingsMatrix,
              "Penalties"=list("freq"=freq,
                               "size"=size,
                               "startIndex"=0),
              "constraints"=constraints))
}


#' The main wrapper function that is also exposed. Fits the desired model.
#' @export
pupil_solve <- function(pulse_spacing,
                        data,
                        model="WIER_SHARED",
                        n=10.1,
                        t_max=930,f=1/(10^24),
                        pulse_dropping_factor=5,
                        maxiter_inner=10000,
                        maxiter_outer=25,
                        convergence_tol=1e-08,
                        should_collect_progress=F,
                        start_lambda=0.1,
                        should_accum_H=T,
                        init_cf = NULL,
                        expand_by = 800,
                        sample_length = 20) {
  
  
  # Create y vector
  y <- matrix(nrow = length(data$pupil),ncol=1)
  y[,1] <- data$pupil
  
  # Expand time for pulses that happened before the the time-window
  # that is considered.
  expanded_time <- data$time
  
  if(expand_by > 0){
    time_expansion <- seq(0,(expand_by - sample_length), by = sample_length)
    expanded_time <- rep(c(time_expansion,
                           (unique(data$time) + expand_by)),
                         times=length(unique(data$subject)))
  }
  
  # Create pulse location vector
  last_pulse <- length(unique(expanded_time)) - (pulse_dropping_factor * round(t_max/100)) - 1
  pulse_locations <- seq(1,last_pulse,by=pulse_spacing)
  real_locations <- unique(expanded_time)[pulse_locations] - expand_by
  
  # Create model setup
  if (model == "WIER_SHARED") {
    # Wierda et al. (2012) model, but with shared penalty!
    setup <- WIER_SHARED_NNLS_model_setup(expanded_time,
                                          (expand_by/sample_length),
                                          data$time,data$subject,
                                          pulse_locations,
                                          n,t_max,f)
  } else if (model == "WIER_IND") {
    # Wierda et al. (2012) model, but with individual penalties!
    setup <- WIER_IND_NNLS_model_setup(expanded_time,
                                      (expand_by/sample_length),
                                      data$time,data$subject,
                                      pulse_locations,
                                      n,t_max,f)
  } else if (model == "DEN_SHARED") {
    # Denison et al. (2012) model, but with shared penalty!
    setup <- DEN_SHARED_NNLS_model_setup(expanded_time,
                                         (expand_by/sample_length),
                                         data$time,data$subject,
                                         pulse_locations,
                                         n,t_max,f)
  } else if (model == "WIER_DEN_SHARED") {
    # Combined model with shared penalty!
    setup <- WIER_DEN_SHARED_NNLS_model_setup(expanded_time,
                                              (expand_by/sample_length),
                                              data$time,data$subject,
                                              pulse_locations,
                                              n,t_max,f)
  } else {
    stop("Model not specified.")
  }
  
  # Basic identifiable check
  if(ncol(setup$X) > nrow(setup$X)) {
    stop("Model is not identifiable, reduce `expand_by` or increase `pulse_dropping_factor`.")
  }
  
  # Initialize coefficients if none are provided
  if(is.null(init_cf)){
    init_cf <- matrix(nrow=ncol(setup$X),ncol=1)
    init_cf[,1] <- runif(n=ncol(setup$X))
  }

  solved_pup <- wrapAmSolve(setup$X,y,init_cf,setup$constraints,
                            setup$Penalties$freq,setup$Penalties$size,
                            setup$Penalties$startIndex,maxiter_outer,
                            maxiter_inner,convergence_tol,
                            should_collect_progress,start_lambda,
                            should_accum_H)
  
  return(list("coef" = solved_pup$coefficients,
              "modelmat" = setup$X,
              "lambdas" = solved_pup$finalLambdas,
              "sigma" = solved_pup$sigma,
              "convergence" = solved_pup$convergence,
              "coefHistory" = solved_pup$coefChanges,
              "expandedTime" = expanded_time,
              "pulseLocations" = pulse_locations,
              "realLocations" = real_locations,
              "pulsesInTime"= real_locations >= 0,
              "setup"=setup,
              "resid"= y - setup$X %*% solved_pup$coefficients,
              "fitted"= setup$X %*% solved_pup$coefficients))
}

#' Bootstrap estimate of the standard error in the pulse weights. Calculation
#' is based on the "bootstrapping residuals" approach presented in chapter 9
#' of "An introduction to the Bootstrap" by Efron & Tribshirani. Algorithm 6.1
#' for calculating standard error estimates based on bootstrap samples
#' from chapter 6 is used for the actual calculations.
#' @export
bootstrap_papss_standard_error <- function(cf,
                                           pupil_var,
                                           setup,
                                           warm_start=T,
                                           N=100,n=10.1,
                                           t_max=930,f=1/(10^24),
                                           maxiter_inner=10000,
                                           maxiter_outer=25,
                                           convergence_tol=1e-08,
                                           start_lambda=0.1,
                                           should_accum_H=T) {
  
  # Get pupil into required format
  y <- matrix(nrow = length(pupil_var),ncol=1)
  y[,1] <- pupil_var
  
  # Get prediction
  pred <- setup$X %*% cf
  
  # Calculate residual vector
  resid <- y - pred
  
  # Get N * length(time_var) bootstrap samples from residual vector
  resid_b <- sample(resid,size=(N * length(resid)),replace = T)
  
  # Prepare bootstrap sets
  R_b <- matrix(resid_b,nrow = N,ncol = length(resid))
  
  # Prepare storage for bootstrap parameter estimates
  B_b <- matrix(0,nrow=N,ncol=length(cf)) # individual samples
  B_b_m <- rep(0,length.out=length(cf)) # bootstrap mean
  
  for (b_iter in 1:N) {
    
    # Prepare y*
    y_star <- pred + R_b[b_iter,]
    
    # Re-use coefficients as starting values or estimate from scratch
    init_cf <- cf
    
    if(!warm_start){
      init_cf <- matrix(nrow=ncol(setup$X),ncol=1)
      init_cf[,1] <- runif(n=ncol(setup$X))
    }
    
    # Solve sample
    solved_pup <- wrapAmSolve(setup$X,y_star,init_cf,setup$constraints,
                              setup$Penalties$freq,setup$Penalties$size,
                              setup$Penalties$startIndex,maxiter_outer,
                              maxiter_inner,convergence_tol,
                              F,start_lambda,should_accum_H)
    
    # Update estimated parameters
    B_b[b_iter,] <- solved_pup$coefficients
    B_b_m <- B_b_m + solved_pup$coefficients
  }
  
  # Finalize bootstrap mean
  B_b_m <- B_b_m/N
  
  # Subtract mean row-wise from bootstrap parameters
  # See: https://stackoverflow.com/questions/24520720/
  # subtract-a-constant-vector-from-each-row-in-a-matrix-in-r
  B_b_diff <- sweep(B_b,2,B_b_m)
  
  # Square all differences
  B_b_diff_pow <- B_b_diff**2
  
  # Calculate col sum
  B_b_sum <- colSums(B_b_diff_pow)
  
  # Divide by N-1 and calculate root to get standard error
  B_b_se <- sqrt(B_b_sum/(N-1))
  
  return(list("standardError"=B_b_se,
              "mean"=B_b_m,
              "individualParams"=B_b))
}
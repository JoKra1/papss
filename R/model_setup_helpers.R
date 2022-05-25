#' @title
#' Create "basis-functions" that make up the pupil spline.
#' 
#' @description
#' The basis is calculated using the pupil response function originally
#' described by Hoeks & levelt (1993). Wierda et al. (2012) refined it
#' by adding an additional scaling paramter (f). The code here makes use
#' of the convolution operation between a spike unit vector and said pupil
#' response function and is based on the code provided by Wierda et al (2012)
#' in their supplementary materials.
#' 
#' @details
#' See: Wierda, S. M., van Rijn, H., Taatgen, N. A., & Martens, S. (2012).
#' Pupil dilation deconvolution reveals the dynamics of attention at high
#' temporal resolution. Proceedings of the National Academy of Sciences of
#' the United States of America, 109(22), 8456–8460.
#' https://doi.org/10.1073/pnas.1201858109
#' See: Hoeks, B., & Levelt, W. (1993). Pupillary dilation as a measure of
#' attention: A quantitative system analysis.
#' Behav. Res. Meth. Ins. C., 25, 16–26.
#' 
#' @param i An integer representing index of a basis function.
#' @param time A numeric vector containing positive time values in ms
#' @param expanded_time A numeric vector containing positive time values in ms, expanded by a certain amount of ms
#' @param expand_by Expansion time in ms passed to papss::pupil_solve(expand_by=) divided by sample length in ms
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

#' @title
#' Intercept term creation
#' 
#' @description
#' Creates a single intercept term (at the population level).
#' Can be manipulated via "by" argument as done in mgcv (Wood, 2017).
#' 
#' @param time A numeric vector containing positive time values in ms
#' @export
create_constant_term <- function(time) {
  intercept <- rep(1,length(time))
  return(intercept)
}

#' @title
#' Slope term creation
#' 
#' @description
#' Creates a single slope term (at the population level)
#' Can be manipulated via "by" argument as done in mgcv (Wood, 2017).
#' Note: Use of unique values that are just repeated requires all levels
#' of the category factor (usually subjects) to have the same values on the
#' variable time.
#' 
#' @param time_unq A numeric vector containing unique positive time values in ms
#' @param n_cat Integer, number of categories for which slope should be repeated.
#' @export
create_slope_term <- function(time_unq,n_cat) {
  slope <- 1:length(time_unq)
  slope <- rep(slope,n_cat)
  return(slope)
}

#' @title
#' Spike matrix creation
#' 
#' @description
#' Create basis part of the model matrix as described by Wood (2017).
#' Essentially, this means creating a matrix with a column for each
#' pulse location and to assign a pupil basis (convolved with that location)
#' to each column.
#' 
#' @details
#' See: Wood, S. N. (2017). Generalized Additive Models: An Introduction with R,
#' Second Edition (2nd ed.). Chapman and Hall/CRC.
#' 
#' @param expanded_time A numeric vector containing positive time values in ms, expanded by a certain amount of ms
#' @param expand_by Expansion time in ms passed to papss::pupil_solve(expand_by=) divided by sample length in ms
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

#' @title
#' Expand term by factor
#'
#'@description
#' Achieves "by" keyword functionality available in mgcv 
#' to enable by-factor smooths (Wood, 2017)!
#'
#' @details
#' See: Wood, S. N. (2017). Generalized Additive Models: An Introduction with R,
#' Second Edition (2nd ed.). Chapman and Hall/CRC.
#' See 'by' at: https://www.rdocumentation.org/packages/mgcv/versions/1.8-38/topics/s
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

#' @title
#' Model template based on Wierda et al.(2012)'s model
#'
#' @description
#' Creates trainings matrix and penalties for a fully penalized
#' version of the original Wierda et al. (2012) model.
#' Like by Wierda et al. separate sets of coefficients are estimated per
#' subject with the h_basis terms. All those coefficients are again constrained
#' to be non-negative. Model also includes the slope terms for each subject that
#' Wierda et al. (2012) introduced. These are not constrained - since their
#' primary purpose was to account for drift in the pupil averages of individual
#' subjects.
#' This model differs from the one by Wierda et al. (2012) in that it penalizes
#' the coefficients corresponding to the h_basis terms. Specifically, it
#' enforces a single penalty term shared by all subjects (e.g., this is
#' similar to the 'fs' basis in mgcv). The form of the penalty expressed on all of the
#' basis functions is a simple identity matrix. We here also penalize all slope
#' terms, again with a single penalty.
#' 
#' @param expanded_time A numeric vector containing positive time values in ms, expanded by a certain amount of ms
#' @param expand_by Expansion time in ms passed to papss::pupil_solve(expand_by=) divided by sample length in ms
#' @param time A numeric vector containing positive time values in ms
#' @param fact A factor vector containing factor level identifiers
#' @param pulse_locations A numeric vector containing index values of pulse loc.
#' @param n Parameter defined by Hoeks & Levelt (number of laters)
#' @param t_max Parameter defined by Hoeks & Levelt (response maximum in ms)
#' @param f Parameter defined by Wierda et al. (scaling factor)
#' @export
WIER_SHARED_NNLS_model_setup <- function(expanded_time,expand_by,time,fact,pulse_locations,n,t_max,f) {
  
  # Extract number of factor levels
  n_fact <- length(unique(fact))
  # Setup model matrix
  slope <- create_slope_term(unique(time),n_fact)
  spike_matrix <- create_spike_matrix_term(expanded_time,expand_by,time,pulse_locations,n,t_max,f)
  slope_matrix <- term_by_factor(slope,fact)
  spike_matrix_by <- term_by_factor(spike_matrix,fact)
  
  trainingsMatrix <- cbind(slope_matrix,spike_matrix_by)
  
  # Setup Penalty definition to be implemented by c++
  # one individual slope penalty plus one
  # shared penalty for all factor levels (expressed on bases)
  freq <- c(1, n_fact)
  # individual penalty is of size n_fact*n_fact Also, Each of the shared
  # penalties is of size length(pulse_locations)*length(pulse_locations)
  size <- c(n_fact,length(pulse_locations))
  
  
  # Define positive constraints
  constraints <- rep("c",ncol(trainingsMatrix))
  constraints[1:n_fact] <- "u" # All parametric terms are unconstrained
  
  return(list("X"=trainingsMatrix,
              "Penalties"=list("freq"=freq,
                               "size"=size,
                               "startIndex"=0),
              "constraints"=constraints))
}

#' @title
#' Alternative model template based on Wierda et al.(2012)'s model
#'
#' @description
#' Creates trainings matrix and penalties for a fully penalized
#' version of the original Wierda et al. (2012) model.
#' Like by Wierda et al. separate sets of coefficients are estimated per
#' subject with the h_basis terms. All those coefficients are again constrained
#' to be non-negative. Model also includes the slope terms for each subject that
#' Wierda et al. (2012) introduced. These are not constrained - since their
#' primary purpose was to account for drift in the pupil averages of individual
#' subjects.
#' This model differs from the one by Wierda et al. (2012) in that it penalizes
#' the coefficients corresponding to the h_basis terms. Specifically, it
#' enforces a **penalty term for each individual subject** (e.g., this is
#' similar to using the 'by' keyword in mgcv). The form of the penalty expressed on all of the
#' basis functions is a simple identity matrix. We here also penalize all slope
#' terms, again with a single penalty.
#' 
#' @param expanded_time A numeric vector containing positive time values in ms, expanded by a certain amount of ms
#' @param expand_by Expansion time in ms passed to papss::pupil_solve(expand_by=) divided by sample length in ms
#' @param time A numeric vector containing positive time values in ms
#' @param fact A factor vector containing factor level identifiers
#' @param pulse_locations A numeric vector containing index values of pulse loc.
#' @param n Parameter defined by Hoeks & Levelt (number of laters)
#' @param t_max Parameter defined by Hoeks & Levelt (response maximum in ms)
#' @param f Parameter defined by Wierda et al. (scaling factor)
#' @export
WIER_IND_NNLS_model_setup <- function(expanded_time,expand_by,time,fact,pulse_locations,n,t_max,f) {
  
  # Extract number of factor levels
  n_fact <- length(unique(fact))
  # Setup model matrix
  slope <- create_slope_term(unique(time),n_fact)
  spike_matrix <- create_spike_matrix_term(expanded_time,expand_by,time,pulse_locations,n,t_max,f)
  slope_matrix <- term_by_factor(slope,fact)
  spike_matrix_by <- term_by_factor(spike_matrix,fact)
  
  trainingsMatrix <- cbind(slope_matrix,spike_matrix_by)
  
  # Setup Penalty definition to be implemented by c++
  # one individual slope penaltiy plus one
  # shared penalty for all factor levels (expressed on bases)
  freq <- c(1, rep(1,length.out=n_fact))
  # individual penalty is of size n_fact*n_fact Also, Each of the shared
  # penalties is of size length(pulse_locations)*length(pulse_locations)
  size <- c(n_fact,rep(length(pulse_locations),length.out=n_fact))
  
  
  # Define positive constraints
  constraints <- rep("c",ncol(trainingsMatrix))
  constraints[1:n_fact] <- "u" # All parametric terms are unconstrained
  
  return(list("X"=trainingsMatrix,
              "Penalties"=list("freq"=freq,
                               "size"=size,
                               "startIndex"=0),
              "constraints"=constraints))
}

#' @title
#' Model template based on Denison et al.(2020)'s model investigation
#'
#' @description
#' Creates trainings matrix and penalties for a fully penalized
#' version inspired by the models investigated by Denison et al. (2020).
#' Similarly, separate sets of coefficients are estimated per
#' subject with the h_basis terms. All those coefficients are again constrained
#' to be non-negative. The model also includes an intercept/offset term for each
#' subject as introduced by Denison et al., (2020. These are not constrained.
#' This model again penalizes
#' the coefficients corresponding to the h_basis terms. Specifically, it
#' enforces a single penalty term shared by all subjects (e.g., this is
#' similar to the 'fs' basis in mgcv). The form of the penalty expressed on all of the
#' basis functions is a simple identity matrix. We here also penalize all intercept
#' terms, again with a single penalty.
#' 
#' @param expanded_time A numeric vector containing positive time values in ms, expanded by a certain amount of ms
#' @param expand_by Expansion time in ms passed to papss::pupil_solve(expand_by=) divided by sample length in ms
#' @param time A numeric vector containing positive time values in ms
#' @param fact A factor vector containing factor level identifiers
#' @param pulse_locations A numeric vector containing index values of pulse loc.
#' @param n Parameter defined by Hoeks & Levelt (number of laters)
#' @param t_max Parameter defined by Hoeks & Levelt (response maximum in ms)
#' @param f Parameter defined by Wierda et al. (scaling factor)
#' @export
DEN_SHARED_NNLS_model_setup <- function(expanded_time,expand_by,time,fact,pulse_locations,n,t_max,f) {
  
  # Extract number of factor levels
  n_fact <- length(unique(fact))
  # Setup model matrix
  intercept <- create_constant_term(time)
  spike_matrix <- create_spike_matrix_term(expanded_time,expand_by,time,pulse_locations,n,t_max,f)
  intercept_matrix <- term_by_factor(intercept,fact)
  spike_matrix_by <- term_by_factor(spike_matrix,fact)
  
  trainingsMatrix <- cbind(intercept_matrix,spike_matrix_by)
  
  # Setup Penalty definition to be implemented by c++
  # one individual intercept penalty plus one
  # shared penalty for all factor levels (expressed on bases)
  freq <- c(1, n_fact)
  # individual penalty is of size of size n_fact*n_fact Also, Each of the shared
  # penalties is of size length(pulse_locations)*length(pulse_locations)
  size <- c(n_fact,length(pulse_locations))
  
  
  # Define positive constraints
  constraints <- rep("c",ncol(trainingsMatrix))
  constraints[1:n_fact] <- "u" # All parametric terms are unconstrained
  
  return(list("X"=trainingsMatrix,
              "Penalties"=list("freq"=freq,
                               "size"=size,
                               "startIndex"=0),
              "constraints"=constraints))
}

#' @title
#' Model template based on Denison et al.(2020)'s model investigation and Wierda et al. (2012)'s model
#'
#' @description
#' Like by Wierda et al. (2012) separate sets of coefficients are estimated per
#' subject with the h_basis terms. All those coefficients are again constrained
#' to be non-negative. The moodel also includes the intercept terms for each subject that
#' Denison et al. (2012) introduced as well as the slope term introduced by
#' Wierda et al., (2012). These are not constrained.
#' 
#' This model again penalizes the coefficients corresponding to the h_basis terms.
#' Specifically, it enforces a single penalty term shared by all subjects (e.g., this is
#' similar to the 'fs' basis in mgcv). The form of the penalty expressed on all of the
#' basis functions is a simple identity matrix. We here also penalize all intercept
#' and slope terms, again with a single penalty.
#' 
#' @param expanded_time A numeric vector containing positive time values in ms, expanded by a certain amount of ms
#' @param expand_by Expansion time in ms passed to papss::pupil_solve(expand_by=) divided by sample length in ms
#' @param time A numeric vector containing positive time values in ms
#' @param fact A factor vector containing factor level identifiers
#' @param pulse_locations A numeric vector containing index values of pulse loc.
#' @param n Parameter defined by Hoeks & Levelt (number of laters)
#' @param t_max Parameter defined by Hoeks & Levelt (response maximum in ms)
#' @param f Parameter defined by Wierda et al. (scaling factor)
#' @export
WIER_DEN_SHARED_NNLS_model_setup <- function(expanded_time,expand_by,time,fact,pulse_locations,n,t_max,f) {
  
  # Extract number of factor levels
  n_fact <- length(unique(fact))
  
  # Setup model matrix
  intercept <- create_constant_term(time)
  slope <- create_slope_term(unique(time),n_fact)
  spike_matrix <- create_spike_matrix_term(expanded_time,expand_by,time,pulse_locations,n,t_max,f)
  intercept_matrix <- term_by_factor(intercept,fact)
  slope_matrix <- term_by_factor(slope,fact)
  spike_matrix_by <- term_by_factor(spike_matrix,fact)
  
  trainingsMatrix <- cbind(intercept_matrix,slope_matrix,spike_matrix_by)
  
  # Setup Penalty definition to be implemented by c++
  # one individual intercept  and slope penalty plus one
  # shared penalty for all factor levels (expressed on bases)
  freq <- c(1, 1, n_fact)
  # individual penalties are of size n_fact*n_fact Also, Each of the shared
  # penalties is of size length(pulse_locations)*length(pulse_locations)
  size <- c(n_fact, n_fact, length(pulse_locations))
  
  
  # Define positive constraints
  constraints <- rep("c",ncol(trainingsMatrix))
  constraints[1:(2 * n_fact)] <- "u" # All parametric terms are unconstrained
  
  return(list("X"=trainingsMatrix,
              "Penalties"=list("freq"=freq,
                               "size"=size,
                               "startIndex"=0),
              "constraints"=constraints))
}


#' @title
#' Main function that fits the penalized additive pupil model.
#'
#' @description
#' This function fits one of the penalized pupil models that are available with this package.
#' @details
#' See the artificial_data_analysis vignette for details and usage examples.
#'
#' @param pulse_spacing Model pulses every 'pulse_spacing' samples. Setting this to 1
#' ensures 1 pulse every sample
#' @param data Aggregated data with a time and pupil column. Also needs a factor column
#' @param factor_id Name of the factor column. Model will estimate demand trajectory for each level of this factor
#' @param model Model template.
#' @param n Choice for parameter defined by Hoeks & Levelt (number of laters)
#' @param t_max Choice for parameter defined by Hoeks & Levelt (response maximum in ms)
#' @param f Choice for parameter defined by Wierda et al. (scaling factor)
#' @param pulse_dropping_factor Last pulse is modelled at index corresponding to: length(unique(data$time)) - (pulse_dropping_factor * round(t_max/100)) - 1
#' @param maxiter_inner Maximum steps taken by inner optimizer
#' @param maxiter_outer Maximum steps taken by outer optimizer
#' @param convergence_tol Convergence check to terminate early
#' @param should_collect_progress If T, then the entire coefficient update history is collected and returned. VERY COSTLY.
#' @param start_lambda Initial lambda value. Must be > 0 if a penalty should be used! Setting this to 0 and maxiter_outer=1, leads to estimation of an un-penalized additive model, i.e., recovers the traditional NNLS estimate used by Wierda et al. (2012) and Denison et al. (2012).
#' @param should_accum_H Whether Hessian should be approximated using BFGS rule or not (see Fletcher, R. (2000). Practical Methods of Optimization). If not, then least squares Hessian matrix is used. With the BFGS rule models ended up being much smoother in our simulations. So this should be set to true if under-smoothing is observed. However, the BFGS update is much more costly and takes much more time!
#' @param init_cf NULL or vector with initial coefficient estimate
#' @param expand_by Time in ms by which to expand the time-series in the past. Then pulses that happened before the recorded time-window can still be approximated! See artificial_data_analysis vignette for details.
#' @param sample_length Duration in ms of a single sample. If pupil dilation time-course was down-sampled to 50HZ, set this to 20
#' @export
pupil_solve <- function(pulse_spacing,
                        data,
                        factor_id="subject",
                        model="WIER_SHARED",
                        n=10.1,
                        t_max=930,f=1/(10^24),
                        pulse_dropping_factor=5,
                        maxiter_inner=10000,
                        maxiter_outer=25,
                        convergence_tol=1e-08,
                        should_collect_progress=F,
                        start_lambda=0.1,
                        should_accum_H=F,
                        init_cf = NULL,
                        expand_by = 800,
                        sample_length = 20) {
  
  
  # Create y vector
  y <- matrix(nrow = length(data$pupil),ncol=1)
  y[,1] <- data$pupil
  
  # Extract factor variable
  fact <- data[,colnames(data) == factor_id]
  
  # Expand time for pulses that happened before the the time-window
  # that is considered.
  expanded_time <- data$time
  
  if(expand_by > 0){
    time_expansion <- seq(0,(expand_by - sample_length), by = sample_length)
    expanded_time <- rep(c(time_expansion,
                           (unique(data$time) + expand_by)),
                         times=length(unique(fact)))
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
                                          data$time,fact,
                                          pulse_locations,
                                          n,t_max,f)
  } else if (model == "WIER_IND") {
    # Wierda et al. (2012) model, but with individual penalties!
    setup <- WIER_IND_NNLS_model_setup(expanded_time,
                                      (expand_by/sample_length),
                                      data$time,fact,
                                      pulse_locations,
                                      n,t_max,f)
  } else if (model == "DEN_SHARED") {
    # Denison et al. (2012) model, but with shared penalty!
    setup <- DEN_SHARED_NNLS_model_setup(expanded_time,
                                         (expand_by/sample_length),
                                         data$time,fact,
                                         pulse_locations,
                                         n,t_max,f)
  } else if (model == "WIER_DEN_SHARED") {
    # Combined model with shared penalty!
    setup <- WIER_DEN_SHARED_NNLS_model_setup(expanded_time,
                                              (expand_by/sample_length),
                                              data$time,fact,
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

#' @title
#' Boostrap estimate of standard error
#'
#' @description
#' Bootstrap estimate of the standard error in the pulse weights. Calculation
#' is based on the "bootstrapping residuals" approach presented in chapter 9
#' of "An introduction to the Bootstrap" by Efron & Tribshirani. Algorithm 6.1
#' for calculating standard error estimates based on bootstrap samples
#' from chapter 6 is used for the actual calculations.
#' 
#' @param cf Coefficients recovered by papss::pupil_solve()
#' @param pupil_var Pupil column from aggregated data frame
#' @param setup Model setup object returned by papss::pupil_solve()
#' @param warm_start Whether or not to re-use cf - if F then coefficients are estimated from scratch for every bootstrap sample
#' @param N Number of repetitions for boot-strapping
#' @param n Choice for parameter defined by Hoeks & Levelt (number of laters)
#' @param t_max Choice for parameter defined by Hoeks & Levelt (response maximum in ms)
#' @param f Choice for parameter defined by Wierda et al. (scaling factor)
#' @param maxiter_inner Maximum steps taken by inner optimizer
#' @param maxiter_outer Maximum steps taken by outer optimizer
#' @param convergence_tol Convergence check to terminate early
#' @param start_lambda Initial lambda value. Must be > 0 if a penalty should be used! Setting this to 0 and maxiter_outer=1, leads to estimation of an un-penalized additive model, i.e., recovers the traditional NNLS estimate used by Wierda et al. (2012) and Denison et al. (2012).
#' @param should_accum_H Whether Hessian should be approximated using BFGS rule or not. If not, then least squares Hessian matrix is used. With the BFGS rule models ended up being much smoother in our simulations. So this should be set to true if under-smoothing is observed. However, the BFGS update is much more costly and takes much more time!
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
                                           should_accum_H=F) {
  
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

#' @title
#' Cross-validation for t_max estimation
#'
#' @description
#' Denison et al. (2020) report large variance in the optimal t_max parameters
#' between subjects. This function can thus be used to recover the optimal parameter
#' for subjects using cross-validation. The same procedure utilized by Denison et
#' al. (2020) is adopted here: 1/n trials are held-out and the model is fitted
#' on the remaining trials. The error between the average on the held-out set and
#' the model prediction is then taken as the cross-validation error. The held-out set
#' cycles through the entire data-set resulting in n repetitions and n cross-validation errors.
#' The average cross-validation error for a specific t_max is then reported.
#' 
#' @details
#' Note that different forms of cross-validation are possible depending on the
#' experimental design and one's assumptions. It is possible to optimize t_max for
#' each subject for each condition individually (then only data from one subject and
#' one condition should be passed to the function) or across conditions (then data from
#' all conditions should be passed to the function). Based on the findings by Denison
#' et al. (2020), the latter is likely sufficient and more appropriate.
#'
#' @param cand_tmax vector with all t_max values to be considered
#' @param folds list of vectors, each vector corresponds to fold and contains trial values to be held-out in that fold!
#' @param pulse_spacing Model pulses every 'pulse_spacing' samples. Setting this to 1
#' ensures 1 pulse every sample
#' @param trial_data trial-level data with a time and pupil column. Also needs a factor column
#' @param factor_id Name of the factor column. Model will estimate demand trajectory for each level of this factor
#' @param model Model template.
#' @param n Choice for parameter defined by Hoeks & Levelt (number of laters)
#' @param t_max Choice for parameter defined by Hoeks & Levelt (response maximum in ms)
#' @param f Choice for parameter defined by Wierda et al. (scaling factor), can also be a vector with values for each t_max candidate
#' @param pulse_dropping_factor Last pulse is modelled at index corresponding to: length(unique(data$time)) - (pulse_dropping_factor * round(t_max/100)) - 1
#' @param maxiter_inner Maximum steps taken by inner optimizer
#' @param maxiter_outer Maximum steps taken by outer optimizer
#' @param convergence_tol Convergence check to terminate early
#' @param start_lambda Initial lambda value. Must be > 0 if a penalty should be used! Setting this to 0 and maxiter_outer=1, leads to estimation of an un-penalized additive model, i.e., recovers the traditional NNLS estimate used by Wierda et al. (2012) and Denison et al. (2012).
#' @param should_accum_H Whether Hessian should be approximated using BFGS rule or not. If not, then least squares Hessian matrix is used. With the BFGS rule models ended up being much smoother in our simulations. So this should be set to true if under-smoothing is observed. However, the BFGS update is much more costly and takes much more time!
#' @param init_cf NULL or vector with initial coefficient estimate
#' @param expand_by Time in ms by which to expand the time-series in the past. Then pulses that happened before the recorded time-window can still be approximated! See artificial_data_analysis vignette for details.
#' @param sample_length Duration in ms of a single sample. If pupil dilation time-course was down-sampled to 50HZ, set this to 20
#' @param should_plot Whether or not fit plots should be generated as well.
#' @export
cross_val_tmax <- function(cand_tmax,
                        folds,
                        pulse_spacing,
                        trial_data,
                        factor_id="subject",
                        model="WIER_SHARED",
                        n=10.1,
                        f=1/(10^24),
                        pulse_dropping_factor=5,
                        maxiter_inner=10000,
                        maxiter_outer=25,
                        convergence_tol=1e-08,
                        start_lambda=0.1,
                        should_accum_H=F,
                        init_cf = NULL,
                        expand_by = 800,
                        sample_length = 20,
                        should_plot=T){
  
  # Collect cross-validation errors
  errs <- c()
  
  # Loop over all t_max candidates
  for (tmc_i in 1:length(cand_tmax)){
    
    sqrt_err <- 0
    
    tmc <- cand_tmax[tmc_i]
    fc <- f
    
    if(length(f) > 1) {
      fc <- f[tmc_i]
    }
    
    # Now handle each fold
    for (fold in folds) {

      # These are the trials excluded from the model training set.
      held_out_dat <- trial_data[trial_data$num_trial %in% fold,]
      
      # Calculate average for held-out set.
      aggr_held_out <- aggregate(list("pupil"=held_out_dat$pupil),
                                 by=list("factor"=held_out_dat[,colnames(held_out_dat) == factor_id],
                                         "time"=held_out_dat$time),FUN=mean)
      
      # This is the remaining data on which the model is fitted.
      remaining_dat <- trial_data[!(trial_data$num_trial %in% fold),]
      
      # This are the actual averages passed to the model
      aggr_remaining <- aggregate(list("pupil"=remaining_dat$pupil),
                                  by=list("factor"=remaining_dat[,colnames(remaining_dat) == factor_id],
                                          "time"=remaining_dat$time),FUN=mean)
      
      # Check that all data-points exist in both sets.
      aggr_held_out$timeCond <- interaction(aggr_held_out$time,aggr_held_out$factor)
      aggr_remaining$timeCond <- interaction(aggr_remaining$time,aggr_remaining$factor)
      
      # If one is missing completely in held-out set, just set that one to the 
      # from remaining set.
      if(nrow(aggr_held_out) < nrow(aggr_remaining)){
        cat("Warning: held-out set has ", nrow(aggr_remaining)- nrow(aggr_held_out), "time-points less. Assigning missing ones from remaining set.\n")
        aggr_held_out <- rbind(aggr_held_out,aggr_remaining[!(aggr_remaining$timeCond %in% aggr_held_out$timeCond),])
        cat("Now, nrow(held-out) == ", nrow(aggr_held_out), " nrow(remaining) == ", nrow(aggr_remaining), ".\n")
      }
      
      # Sort data correctly
      aggr_held_out <- aggr_held_out[order(aggr_held_out$factor),]
      aggr_remaining <- aggr_remaining[order(aggr_remaining$factor),]
      
      # Solve pupil model
      solvedPupil <- pupil_solve(pulse_spacing,
                                 aggr_remaining,
                                 "factor",
                                 model,
                                 n,
                                 t_max=tmc,
                                 fc,
                                 pulse_dropping_factor,
                                 maxiter_inner,
                                 maxiter_outer,
                                 convergence_tol,
                                 F,
                                 start_lambda,
                                 should_accum_H,
                                 init_cf,
                                 expand_by,
                                 sample_length)
      
      recovered_coef <- solvedPupil$coef
      model_mat <- solvedPupil$modelmat
      
      # Calculate held-out residuals
      pred <- model_mat %*% recovered_coef
      res <- aggr_held_out$pupil - pred
      
      # And squared error
      sum_sqrt_res <- sum(res**2)
      
      # Update average error
      sqrt_err <- sqrt_err + ((nrow(held_out_dat)/nrow(trial_data)) * sum_sqrt_res)
      
      if(should_plot){
        plot(1:nrow(model_mat),
             aggr_remaining$pupil,
             type="l",lwd=3,
             xlab="Index",
             ylab="Pupil dilation",
             main=paste0("t_max: ", tmc))
        lines(1:nrow(model_mat),
              model_mat %*% recovered_coef,
              lwd=3,
              col= "red",
              lty=2)
        lines(1:nrow(model_mat),
              aggr_held_out$pupil,
              lty=3,
              col="blue")
        legend("topleft",c("Remaining pupil","Predicted pupil","Held-out pupil"),
               lty = c(1,2,3),col=c("black","red","blue"),lwd = 3)
      }
      

    }
    
    # Collect error for this particular t_max
    errs <- c(errs,sqrt_err)
  }
  
  # Plot cross-validation curve
  plot(cand_tmax,errs,ylim=c(min(errs) - 0.01*min(errs),
                             max(errs) + 0.01*max(errs)),
       ylab="Cross-validation error",
       xlab="t_max candidates")
  lines(cand_tmax,errs)
  
  return(errs)
}

#'@title
#' Extracts demand curves
#'
#' @description
#' Extracts demand curves for each level of the factor passed to the pupil_solve
#' call.
#' @param aggr_dat Aggregated data passed to solver
#' @param recovered_coef The coefficients returned by papss::pupil_solve()
#' @param pulse_locations The index vector pointing at pulse locations returned by papss::pupil_solve()
#' @param real_locations The vector with the time-points at which pulses are modelled, returned by papss::pupil_solve()
#' @param pulses_in_time Boolean vector indicating which pulse was within the time-window. Returned by papss::pupil_solve()
#' @param expanded_time The expanded time series returned by papss::pupil_solve()
#' @param expanded_by Expansion time in ms passed to papss::pupil_solve(expand_by=) divided by sample length in ms
#' @param factor_id Name of factor column in aggr_dat
#' @param model Model_template name
#' @param n Choice for n parameter
#' @param t Choice for t_max parameter
#' @param f Choice for f parameter
#' @export
extract_demand_for_fact <- function(aggr_dat,
                                    recovered_coef,
                                    pulse_locations,
                                    real_locations,
                                    pulses_in_time,
                                    expanded_time,
                                    expanded_by,
                                    factor_id="subject",
                                    model="WIER_SHARED",
                                    n = 10.1,
                                    t = 930,
                                    f = 1/(10^24)){
  
  factor <- aggr_dat[,colnames(aggr_dat) == factor_id]
  unq_factor <- unique(factor)
  n_fact <- length(unq_factor)
  
  # First re-create demand-related parts of model matrix.
  semiPredX <- papss::create_spike_matrix_term(unique(expanded_time),
                                               expanded_by,
                                               unique(aggr_dat$time),
                                               pulse_locations,
                                               n,
                                               t,
                                               f)
  
  demand_dat <- NULL
  
  # Start indexing after parametric coefficients
  # i.e., intercept/slope or both in case of WIER_DEN model
  n_param_coef <- n_fact

  if(model == "WIER_DEN_SHARED"){
    n_param_coef <- 2 * n_fact
  }
  
  # Get estimates for each factor level
  for(fi in 1:n_fact){
    
    # Get corresponding spline spike weights
    splineCoef <- recovered_coef[((n_param_coef + 1) + ((fi - 1) * ncol(semiPredX))):(n_param_coef+(fi * ncol(semiPredX)))]
    
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
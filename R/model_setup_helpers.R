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
h_basis <- function(i,time,pulse_locations,n,t_max,f) {
  
  # We need to make sure that bases towards the end do not 'contaminate'
  # the data for the next subject. Since we assume that each subject
  # contributes the same time-steps we form the basis over the unique
  # time-values and then just repeat the basis multiple times (for each subject).
  unq_time <- unique(time)
  
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
  o_restr <- o[1:length(unq_time)]
  
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
create_spike_matrix_term <- function(time,pulse_locations,n,t_max,f) {
  
  spike_matrix <- matrix(nrow = length(time),ncol = (length(pulse_locations)))
  
  for (i in 1:length(pulse_locations)) {
    spike_matrix[,i] <- h_basis(i,time,pulse_locations,n,t_max,f)
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

expand_var_by_rep <- function(var, fact, by) {
  expanded_var <- c()
  for (l in unique(fact)) {
    fact_var <- var[fact == l]
    expanded_var <- c(expanded_var,rep(fact_var[1],times=by),fact_var)
  }
  return(expanded_var)
}

expand_var_by_seq <- function(var, fact, from, step=20) {
  expanded_var <- c()
  for (l in unique(fact)) {
    fact_var <- var[fact == l]
    expanded_var <- c(expanded_var,
                      (seq(-from,(fact_var[1] - step),by=step) + from),
                      (fact_var + from))
  }
  return(expanded_var)
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
WIER_SHARED_NNLS_model_setup <- function(time,subs,pulse_locations,n,t_max,f) {
  
  # Extract number of subjects
  n_subs <- length(unique(subs))
  # Setup model matrix
  slope <- create_slope_term(unique(time),n_subs)
  spike_matrix <- create_spike_matrix_term(time,pulse_locations,n,t_max,f)
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
DEN_SHARED_NNLS_model_setup <- function(time,subs,pulse_locations,n,t_max,f) {
  
  # Extract number of subjects
  n_subs <- length(unique(subs))
  # Setup model matrix
  intercept <- create_constant_term(time)
  spike_matrix <- create_spike_matrix_term(time,pulse_locations,n,t_max,f)
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
WIER_DEN_SHARED_NNLS_model_setup <- function(time,subs,pulse_locations,n,t_max,f) {
  
  # Extract number of subjects
  n_subs <- length(unique(subs))
  # Setup model matrix
  intercept <- create_constant_term(time)
  slope <- create_slope_term(unique(time),n_subs)
  spike_matrix <- create_spike_matrix_term(time,pulse_locations,n,t_max,f)
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
WIER_DEN_SHARED_NNLS_EXP_model_setup <- function(expanded_time,time,subs,pulse_locations,n,t_max,f,expand_by,step) {
  
  # Extract number of subjects
  n_subs <- length(unique(subs))
  
  # Create expansion index vector
  original_index <- expanded_time >= expand_by

  # Setup model matrix
  
  # Two intercepts, one for original period, one for expansion
  intercept <- create_constant_term(expanded_time)
  intercept_expansion <- intercept
  intercept_original <- intercept
  
  intercept_expansion[original_index] = 0
  intercept_original[!original_index] = 0
  
  intercept_expansion_by = term_by_factor(intercept_expansion,subs)
  intercept_original_by = term_by_factor(intercept_original,subs)
  
  
  # Slope, only over original period
  slope <- create_slope_term(unique(expanded_time),n_subs)
  slope[slope <= (expand_by/step)] <- 0
  slope[slope > (expand_by/step)] <- slope[slope > (expand_by/step)] - (expand_by/step)
  slope_matrix <- term_by_factor(slope,subs)
  
  # Spike matrix, over entire period
  spike_matrix <- create_spike_matrix_term(expanded_time,pulse_locations,n,t_max,f)
  spike_matrix_by <- term_by_factor(spike_matrix,subs)
  
  trainingsMatrix <- cbind(intercept_expansion_by,
                           intercept_original_by,
                           slope_matrix,
                           spike_matrix_by)
  
  # Setup Penalty definition to be implemented by c++
  # two individual intercept  and slope penalties plus one
  # shared penalty for all subjects (expressed on bases)
  freq <- c(1, 1, 1, n_subs)
  # individual penalties are of size n_subs*n_subs. Also, Each of the shared
  # penalties is of size length(pulse_locations)*length(pulse_locations)
  size <- c(n_subs, n_subs, n_subs, length(pulse_locations))
  
  
  # Define positive constraints
  constraints <- rep("c",ncol(trainingsMatrix))
  constraints[1:(3 * n_subs)] <- "u" # All parametric terms are unconstrained
  
  return(list("X"=trainingsMatrix,
              "Penalties"=list("freq"=freq,
                               "size"=size,
                               "startIndex"=0),
              "constraints"=constraints))
}

#' The main wrapper function that is also exposed. Fits the desired model.
#' @export
pupil_solve <- function(pulse_locations,real_locations,
                        data,model="WIER_SHARED",n=10.1,
                        t_max=930,f=1/(10^24),
                        maxiter_inner=10000,
                        maxiter_outer=25,
                        convergence_tol=1e-08,
                        should_collect_progress=F,
                        start_lambda=0.1,
                        init_cf = NULL,
                        expand_by = 1000,
                        sample_len=20) {
  
  
  # create optional variables for expansion
  expanded_time <- NULL
  expanded_subs <- NULL
  expanded_pupil <- NULL
  
  # Create y vector
  y <- matrix(nrow = length(data$pupil),ncol=1)
  y[,1] <- data$pupil
  
  # Create model setup
  if (model == "WIER_SHARED") {
    # Wierda et al. (2012) model, but with shared penalty!
    setup <- WIER_SHARED_NNLS_model_setup(data$time,data$subject,
                                          pulse_locations,
                                          n,t_max,f)
  } else if (model == "DEN_SHARED") {
    # Denison et al. (2012) model, but with shared penalty!
    setup <- DEN_SHARED_NNLS_model_setup(data$time,data$subject,
                                         pulse_locations,
                                         n,t_max,f)
  } else if (model == "WIER_DEN_SHARED") {
    # Combined model with shared penalty!
    setup <- WIER_DEN_SHARED_NNLS_model_setup(data$time,data$subject,
                                              pulse_locations,
                                              n,t_max,f)
  } else if (model == "WIER_DEN_SHARED_EXP") {
    # Combined model with shared penalty and baseline expansion!
    
    # Expand time and subject
    expanded_time <- expand_var_by_seq(data$time,data$subject,expand_by,sample_len)
    
    expanded_subs <- as.factor(expand_var_by_rep(as.numeric(data$subject),data$subject,(expand_by/sample_len)))
    
    # Re-calculate y vector
    expanded_pupil <- expand_var_by_rep(data$pupil,data$subject,(expand_by/sample_len))
    y <- matrix(nrow = length(expanded_pupil),ncol=1)
    y[,1] <- expanded_pupil
    
    # Re-calculate pulse locations
    pulse_locations_diff <- pulse_locations[2] - pulse_locations[1]
    last_pulse <- length(unique(expanded_time)) - (3 * round(930/100)) - 1
    pulse_locations <- seq(1,last_pulse,pulse_locations_diff)
    real_locations <- unique(expanded_time)[pulse_locations]
    
    # Now proceed with setup
    setup <- WIER_DEN_SHARED_NNLS_EXP_model_setup(expanded_time,
                                                  data$time,expanded_subs,
                                                  pulse_locations,
                                                  n,t_max,f,expand_by,
                                                  sample_len)
  } else {
    stop("Model not specified.")
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
                            should_collect_progress,start_lambda)
  
  return(list("coef" = solved_pup$coefficients,
              "modelmat" = setup$X,
              "lambdas" = solved_pup$finalLambdas,
              "sigma" = solved_pup$sigma,
              "convergence" = solved_pup$convergence,
              "coefHistory" = solved_pup$coefChanges,
              "usedPulseLocations"=pulse_locations,
              "usedRealLocations"=real_locations,
              "expandedPupil"=expanded_pupil,
              "expandedTime"=expanded_time,
              "expandedSubjects"=expanded_subs))
}
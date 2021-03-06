% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pupil_sim_helpers.R
\name{additive_pupil_sim}
\alias{additive_pupil_sim}
\title{Simulate pupil data to test papss}
\usage{
additive_pupil_sim(
  nk = 20,
  n_sub = 10,
  n_trials = 250,
  pulse_loc_diff = 1,
  n_diff = 4,
  spars_deg = 0.5,
  sub_dev = 0.15,
  slope_dev = 1.5,
  trial_dev = 0.05,
  residual_dev = 15,
  n = 10.1,
  t_max = 930,
  f = 1/(10^24),
  should_plot = T,
  seed = 124
)
}
\arguments{
\item{nk}{Number of knots for the B-spline basis.}

\item{n_sub}{How many subjects to simulate.}

\item{n_trials}{How many trials to simulate.}

\item{pulse_loc_diff}{Assume a pulse every 'pulse_loc_diff' samples (one sample = 20 ms)}

\item{n_diff}{Maximum number of spline basis coefficients with systematic per-subject variation}

\item{spars_deg}{The degree of sparsity enforced in spline basis coefficient vector}

\item{sub_dev}{Standard deviation of normal distribution used to sample by-subject variation}

\item{slope_dev}{Standard deviation of normal distribution used to sample by-subject slope variation}

\item{trial_dev}{Standard deviation of normal distribution used to sample by-trial variation (in spike weights and coefficients)}

\item{residual_dev}{Standard deviation of normal distribution used to sample residuals per trial}

\item{n}{Parameter defined by Hoeks & Levelt (number of laters)}

\item{t_max}{Parameter defined by Hoeks & Levelt (response maximum in ms)}

\item{f}{Parameter defined by Wierda et al. (scaling factor)}

\item{should_plot}{Whether the generated data should be visualized.}

\item{seed}{For replicability}
}
\description{
Uses the pupil response function described by Hoeks & Levelt (1993)
with the refinements proposed by Wierda e t al. (2012) to simulate pupil
dilation time-course data with three sources of deviation:
\enumerate{
\item Simulates systematic per-subject deviation from a 'shared' population trend in demand
\item Simulates per-trial deviation from each subject's individual 'true demand'
\item Assumes random noise (N(0, sigma), with constant sigma) for each trial
}
}
\details{
Demand is modelled according to a simple B-spline (see Eilers & Marx, 2010)
with equally spaced knots with associated coefficients of which only a small
percentage will be different from zero (to ensure that the simulated demand
trajectory is sparse).
}
\examples{
pupil_sim <- additive_pupil_sim(n_sub=5)
sim_dat <- pupil_sim$data
}

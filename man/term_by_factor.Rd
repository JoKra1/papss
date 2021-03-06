% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/model_setup_helpers.R
\name{term_by_factor}
\alias{term_by_factor}
\title{Expand term by factor}
\usage{
term_by_factor(term, fact)
}
\arguments{
\item{term}{A slope, splike_matrix, or intercept.}

\item{fact}{A factor variable.}
}
\description{
Achieves "by" keyword functionality available in mgcv 
to enable by-factor smooths (Wood, 2017)!
}
\details{
See: Wood, S. N. (2017). Generalized Additive Models: An Introduction with R,
Second Edition (2nd ed.). Chapman and Hall/CRC.
See 'by' at: https://www.rdocumentation.org/packages/mgcv/versions/1.8-38/topics/s
Code is based on: https://stats.stackexchange.com/questions/110472/how-exactly-is-the-sum-or-mean-centering-constraint-for-splines-also-w-r-t-g
}

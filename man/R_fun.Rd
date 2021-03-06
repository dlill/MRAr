% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/response_matrices.R
\name{R_fun}
\alias{R_fun}
\title{Compute R from fractional changes of steady states}
\usage{
R_fun(pars_opt, perturbation_prediction, obs_fun, p_fun, pars)
}
\arguments{
\item{pars_opt}{Free parameters a_ij to control the influence of the complexes}

\item{perturbation_prediction}{Simulated steady state perturbation data}

\item{obs_fun}{Observation function which returns the communicating species}

\item{p_fun}{Parameter transformation which returns different conditions}

\item{pars}{Parameters of the ODE model}
}
\description{
Compute R from fractional changes of steady states
}

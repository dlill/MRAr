% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dose_response.R
\name{r_opt_dose_response}
\alias{r_opt_dose_response}
\title{Compute dose responses for one parameter}
\usage{
r_opt_dose_response(which_par, dosages, alpha_scan_range = seq(-9, 9, by = 2),
  pars_opt = pars_opt_0, obs_fun = g, p_fun = (p_log * p_pert),
  pars = pars_0)
}
\arguments{
\item{which_par}{Which dynamical parameter shall be scanned?}

\item{dosages}{Scan range of this parameter}

\item{alpha_scan_range}{Which values of free parameters a shall be scanned?}

\item{pars_opt}{Free parameters a_ij to control the influence of the complexes}

\item{obs_fun}{Observation function which returns the communicating species}

\item{p_fun}{Parameter transformation which returns different conditions}

\item{pars}{Parameters of the ODE model}
}
\description{
Scan free parameters separately
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/condition_generation.R
\name{p_ic_list_fun}
\alias{p_ic_list_fun}
\title{Generate a list of parameter trafos.}
\usage{
p_ic_list_fun(ic_list, pars, mixed = FALSE)
}
\arguments{
\item{ic_list}{Named list of parameter values to be set.}

\item{pars}{The original set of parameters to be replaced.}

\item{mixed}{Should all combinations of the pars be gone through (TRUE)
or should the other pars be set to their first value supplied, when scanning through a certain parameter?}
}
\description{
This function is similar to p_ic_fun, but it returns a list of parfns.
`p_ic_list` needs to be used like `lapply(p_ic_list, function(p_ic) {...(g*x*p_ic)})`.
This allows to do optimization over more than one condition (e.g. replace the dots by `obj`).
}

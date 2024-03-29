% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kmean_util.R
\name{kmean_infer}
\alias{kmean_infer}
\title{Inference method for the partially functional linear model.}
\usage{
kmean_infer(data_list, group, num_pca = 5, tau_1 = c(0.005, 0.05, 0.1),
  tau_2 = 0.1, boot_size = 5000)
}
\arguments{
\item{data_list}{A list of data. Several elements must be present in the list. The reponse \code{y},
the functional covariate \code{x_recv}, the scalar covariates \code{z}, and an index matrix \code{index}. 
The functional covariate \code{x_recv} must be generated from the \code{fda} package by, e.g., spline smoothing. 
The scalar covariates \code{z} is a matrix. The index matrix \code{index} is a data.frame recording the structure of the
data. The first column of \code{index} is the family number, the second column is the within family index. The column names of \code{index}
must be \code{ind_b} and \code{ind_w}.}

\item{group}{The estimated group membership from the function \code{bic_kmean_est}.}

\item{num_pca}{The number of principal components.}

\item{tau_1}{The confidence level is \code{1-tau_1}, can be a vector.}

\item{tau_2}{The proportion of the index of the functional coefficient that might be not covered by the band.}

\item{boot_size}{The size of bootstrap sampling.}
}
\description{
This function outputs the standard error of the scarlar estimates and the confidence band of the functional coefficient.
}

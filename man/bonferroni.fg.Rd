% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/core.R
\name{bonferroni.fg}
\alias{bonferroni.fg}
\title{Bonferroni's correction with fixed \eqn{\Gamma}}
\usage{
bonferroni.fg(d, gamma = 1, mm = c(2, 2, 2), two.sided = TRUE)
}
\arguments{
\item{d}{a matrix of treatment-minus-control differences.}

\item{gamma}{sensitivity parameter (maximum odds different from a randomized experiment).}

\item{mm}{test statistic, either a vector of length 3 or a matrix of three rows where each column corresponds to a U-statistic. Default is the (approximate) Wilcoxon's signed rank test.}

\item{two.sided}{whether a two-sided test should be used. If FALSE, test the one-sided alternative that the center of \code{d} is positive.}
}
\value{
a vector of sensitivity values for each column of \code{d}
}
\description{
Bonferroni's correction with fixed \eqn{\Gamma}
}
\details{
If \code{mm} is a matrix, this function computes a one-sided or two-sided p-value with each statistic (i.e. there is a p-value for every column of \code{d} and every column of $mm$), then does a Bonferroni correction over all the p-values.
}
\author{
Qingyuan Zhao
}


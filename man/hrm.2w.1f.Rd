% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/f3_sub1.R
\name{hrm.2w.1f}
\alias{hrm.2w.1f}
\title{Test for influence of factor A}
\usage{
hrm.2w.1f(X, alpha, group, subgroup, factor, subject, data, H, text = "",
  nonparametric, ranked, varQGlobal, np.correction, tmpQ1g, tmpQ2g)
}
\arguments{
\item{X}{dataframe containing the data in the long table format}

\item{alpha}{alpha level used for the test}

\item{group}{column name of the data frame X specifying the groups}

\item{subgroup}{column name of the subgroups (crossed with groups)}

\item{factor}{column name of the data frame X of within-subject factor}

\item{subject}{column name of the data frame X identifying the subjects}

\item{data}{column name of the data frame X containing the measurement data}

\item{H}{string specifying the hypothesis}

\item{text}{a string, which will be printed in the output}
}
\value{
Returns a data frame consisting of the degrees of freedom, the test value, the critical value and the p-value
}
\description{
Test for influence of factor A
}
\keyword{internal}

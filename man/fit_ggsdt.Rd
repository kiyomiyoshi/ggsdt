\name{fit_ggsdt}
\alias{fit_ggsdt}
\title{Fitting GGSDT model}
\usage{
fit_ggsdt(nR_S1, nR_S2, add_constant = T)
}
\arguments{
  \item{nR_S1}{response frequency for S1 stimuli ordered from highest confidence S1 to highest confidence S2, such as (170, 120, 50, 25, 45, 40)}
  \item{nR_S2}{response frequency for S2 stimuli ordered from highest confidence S1 to highest confidence S2, such as (40, 50, 30, 20, 70, 240)}
  \item{add_constant}{if true, adds small value to response frequency to avoid zero-cell-related issues (default value is TRUE)}
}
\description{
Estimates GGSDT parameters
}
\examples{
fit_ggsdt(nR_S1 = c(170, 120, 50, 25, 45, 40),
          nR_S2 = c(40, 50, 30, 20, 70, 240),
          add_constant = F)
}
\value{
  data.frame with following columns
  \item{mu2}{estimated mean value for S2 distribution}
  \item{alpha2}{estimated alpha value for S2 distribution}
  \item{beta}{estimated beta value shared with S1 and S2 distributions}
  \item{LogLike}{log likelihood}
  \item{sigma1}{SD for S1 distribution calculated with estimated beta value (alpha1 is fixed at 1)}
  \item{sigma2}{SD for S2 distribution calculated with estimated alpha2 and beta values}
  \item{kurtosis}{kurtosis for S1 and S2 distributions caluculated with estimated beta value}
  \item{Xi}{ith decision criterion (from lowest to highest)}
}

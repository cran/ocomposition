\name{negbinom}
\alias{dtnegbin}
\alias{rtnegbin}


\title{Truncated negative binomial distribution function and sampling.}

\description{Distribution function of truncated negative binomial distribution and random draws from the distribution (uses rejection sampling).}

\usage{
dtnegbin(x, mu, dispersion, l.bound)
rtnegbin(N, mu, dispersion, l.bound)
}

\arguments{
  \item{x}{
value at which density is evaluated
}
  \item{mu}{
mean of the distribution
}
  \item{dispersion}{
dispersion parameter
}
  \item{l.bound}{
the lower bound of truncation
}
  \item{N}{
Number of draws
}
}
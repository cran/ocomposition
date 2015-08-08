\name{predict.composition}
\alias{predict.composition}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{
Predicted compositional vector
}

\description{
Computes the predicted composition given a vector of covariates \code{newdata}.
}
\usage{
\method{predict}{composition}(object, newdata, n.method = "median", l.bound = NULL, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{Object of class \code{composition} from the previous call of\code{fitcomp}.}
  \item{newdata}{\code{data.frame} of predictors.}
  \item{n.method}{If \code{"median"} (default), then the number of components is set to the median of the truncated negative binomial distribution, conditional on the estimated parameters; alternatively, it can be set to "mode."}
  \item{l.bound}{Set the lower bound for predicted number of components. If \code{NULL} (default), then \code{l.bound} is taken from the \code{composition} object.}
  \item{\dots}{additional arguments}
}

\value{
Matrix of predicted probabilities of for each component. Rows are MCMC iterations and columns are the components.
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
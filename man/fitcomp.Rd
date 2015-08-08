\name{fitcomp}
\alias{fitcomp}
%- Also NEED an '\alias' for EACH other topic documented here.

\title{Gibbs sampler for parameter estimation}

\description{The main regression function for compositional rank-index data. For units \eqn{i = 1, ..., n}, the response variable is vector \eqn{(y_{i1}, ..., y_{in})}, where \eqn{\sum_j y_{ij} = 1} and \eqn{y_{i1} \geq y_{i2} \geq ... \geq y_{in}} for all \eqn{i} and \eqn{y_{ij} \in [0, 1]} for all \eqn{i} and \eqn{j}. The regression model has two parts: a truncated negative binomial model for the count of non-zero components and a set of seemingly unrelated t regressions for the compositions. See \code{References} for further details.}

\usage{
fitcomp(data.v, data.x, n.formula, v.formula, l.bound = 1, 
n.sample = 100, burn = 0, thin = 1, init = NULL)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data.v}{
Matrix of compositional data: rows for units and columns for components. Rows must add up to 1; if not, they are automatically rescaled. NA values turned into 0 automatically. Ordering done automatically.
}
\item{data.x}{
Data frame with covariates, missing values not allowed.
}
  \item{n.formula}{
formula for the number of components: e.g., \code{ ~ x1 + x2 + factor(z)}.
}
  \item{v.formula}{
 formula for the size of components: e.g., \code{ ~ x1 + x2}.
}

  \item{l.bound}{
lower bound for the negative binomial regression; must be greater or equal to 1; \code{default = 1}.
}
  \item{n.sample}{
number of samples you want to have after burn-in and thinning; default 100
}
  \item{burn}{
number of burn-in samples; default 0
}
  \item{thin}{
thinning of the MCMC chain; default 1
}
  \item{init}{
initial parameters; not required
}
}

\value{
\item{g}{samples of \code{gamma} coefficients for the multivariate regression model}
\item{b}{posterior samples of the coefficients for the negative binomial regression}
\item{mu}{hyperparameters for \code{gamma} coefficients}
\item{rho}{shrinkage hyperparameters for \code{gamma} coefficients}
\item{Sigma}{posterior samples of the covariance matrix}
\item{nu}{degrees of freedom for the Student's \code{t} distribution}
}

\examples{
	
data(data)	
out <- fitcomp(data$v, data$m, ~ log(m), ~ log(m) + log(n), n.sample = 50)
summary(out)

# predict distribution of votes in a country with 5-member median district

v.hat <- predict(out, data.frame(m=5)) 
plot(v.hat)		
}

\references{Rozenas, Arturas (2012) 'A Statistical Model for Party Systems Analysis', \emph{Political Analysis}, 2(20), p.235-247.}

\keyword{ compositional }
\keyword{ multivariate response }
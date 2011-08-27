\name{orcomp}
\alias{orcomp}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Gibbs sampler}
\description{Implements Gibbs sampler for the ordered compositional data model. Has two parts: negative binomial model for the count of components, and a latent seemingly unrelated t regressions for compositions.}
\usage{
orcomp(v, n.formula, v.formula, l.bound = 1, 
data, n.sample = 100, burn = 1, thin = 1, init = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{v}{
Matrix of vote-shares: rows for districts/countries; columns for parties. Rows must add up to 1; if not, they are automatically rescaled. NA values turned into 0 automatically. Ordering done automatically.
}
  \item{n.formula}{
\code{glm} formula for number of parties: \code{n ~ x1 + x2}
}
  \item{v.formula}{
\code{glm} formula for vote-shares: \code{v ~ x1+x2}
}
  \item{l.bound}{
lower bound for the negative binomial regression. default = 0.
}
  \item{data}{
data list: a list of vote-shares and a list of covariates
}
  \item{n.sample}{
number of samples you want to have after burn-in and thinning
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
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
\item{g}{samples of \code{gamma} coefficients for the vote-share model}
\item{b}{posterior samples of the negative binomial regression}
\item{mu}{vote-share model hyperparameters}
\item{rho}{shrinkage hyperparameters}
\item{Sigma}{posterior samples of the covariance matrix}
\item{tau}{latent \code{tau} parameters for t regressions}
\item{nu}{degrees of freedom}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
data(data)
data$n=apply(data$v,1,function(x) sum(x>0))
fit=orcomp(data$v, n.formula=n ~ log(m), v.formula=v ~ log(m),l.bound=1,data=data,n.sample=500)
}

\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

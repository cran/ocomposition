\name{dtlnegbin}
\alias{dtlnegbin}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{neg bin density}
\description{truncated negative binomial log-density}
\usage{
dtlnegbin(x, mu, dispersion, l.bound)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
%%     ~~Describe \code{x} here~~
}
  \item{mu}{
%%     ~~Describe \code{mu} here~~
}
  \item{dispersion}{
%%     ~~Describe \code{dispersion} here~~
}
  \item{l.bound}{
%%     ~~Describe \code{l.bound} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function(x,mu,dispersion,l.bound){
N=length(x)
mu=rep(mu,length=N)
l.bound=rep(l.bound,length=N)
ldensity=ifelse(x>=l.bound,dnbinom(x,mu=mu, size=dispersion,log=TRUE)-
pnbinom(l.bound-.1,mu=mu,size=dispersion,lower.tail=FALSE,log=TRUE),log(1e-300))
return(ldensity)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
\keyword{internal}
\name{rtnegbin}
\alias{rtnegbin}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{sample from tr neg bin}
\description{Rejection sampling from truncated negative binomial distribution.}
\usage{
rtnegbin(N, mu, dispersion, l.bound)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{N}{
%%     ~~Describe \code{N} here~~
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
function(N,mu,dispersion,l.bound){
sample=rep(NA,N)
mu=rep(mu,length=N)
l.bound=rep(l.bound,length=N)
n=1:N
ind=n
while(length(ind)!=0){ ##rejection sampling
sample[ind]=rnbinom(length(ind),mu=mu[ind],size=dispersion)
ind=which(sample<l.bound)
  }
return(sample)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
\keyword{internal}
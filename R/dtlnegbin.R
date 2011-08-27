dtlnegbin <-
function(x,mu,dispersion,l.bound){
N=length(x)
mu=rep(mu,length=N)
l.bound=rep(l.bound,length=N)
ldensity=ifelse(x>=l.bound,dnbinom(x,mu=mu, size=dispersion,log=TRUE)-
pnbinom(l.bound-.1,mu=mu,size=dispersion,lower.tail=FALSE,log=TRUE),log(1e-300))
return(ldensity)
}


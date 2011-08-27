orcomp <-
function(v, n.formula, v.formula, l.bound=1, data=data, n.sample=100, burn=1, thin=1, init=NULL){

X=model.matrix(n.formula, data=data)
G=model.matrix(v.formula, data=data)

v=as.matrix(v)
v[is.na(v)]=0
if (sum(round(apply(v,1,sum),5))/nrow(v)!=1) stop("Vote-shares don't add up to 1. Check your data or rescale")
v=t(apply(v,1,sort,decreasing=TRUE))

y=t(apply(v,1,v.y)) # transform vote-share
n=apply(v,1,function(x)(sum(x>0)))

N=nrow(X); 
L=D=max(n)-1; # DIMENSION OF THE VOTE_SHARE MODEL
K=ncol(X);
T=ncol(G);

#=====================================================
#  COUNT MODEL PARAMETERS 
#====================================================

# L-LIKELIHOOD
n.lik=function(b){
sum=sum(dtlnegbin(n,mu=exp(X%*%b[1:K]),dispersion=exp(b[K+1]),l.bound=l.bound))
return(sum)
}

cat("----  MLE estimates for neg binomial ---- ", '\n')
# MLE ESTIMATES FOR TRUNCATED N BIN MODEL
n.mle=optim(c(rep(0,K),1), n.lik, method="L-BFGS-B",hessian=TRUE, control=list(fnscale=-1),lower=c(rep(-Inf,K),-10),upper=c(rep(Inf,K),10))
H=-.4*solve(n.mle$hessian)



# LOG POSTERIOR
pois.post=function(b){
sum=sum(dtlnegbin(n,mu=exp(X%*%b[1:K]),dispersion=exp(b[K+1]),l.bound=l.bound)) + dnorm(b[K+1],0,10,log=TRUE)
return(sum)
}



# SAMPLING FROM LOG-POSTERIOR
update.b=function(b){
cand=c(rmnorm(1,mean=b,varcov=H))
prob=min(0,pois.post(cand)-pois.post(b))
if (runif(1) < exp(prob)){b = cand}
return(b)
}


#========================================================# 		VOTE-SHARE MODEL PARAMETERS
#========================================================

# Z, gamma, mu, rho

#========================================================
#---     Sweep operator for conditional MVN sampling
#========================================================

cov.f=function(x,Sigma){
if (x[L+1]==L){return(x[-(L+1)])}
if (x[L+1]<L){
l=1:x[L+1]
t=(length(l)+1):L
S=Sigma[t,t]-Sigma[t,l]%*%solve(Sigma[l,l])%*%Sigma[l,t]
m=c(x[l],x[t]%*%chol(S))
return(m)
}
}

mvm=function(x,Sigma,g){
if (x[T+L+1]==L){return(rep(0,L))}
if (x[T+L+1]<L){
l=1:x[T+L+1]
t=(length(l)+1):L
m=x[1:T]%*%g
M=m[t]+Sigma[t,l]%*%solve(Sigma[l,l])%*%(x[(T+1):(T+length(l))]-m[l])
return(c(rep(0,length(l)),M))
}
}


# SAMPLE LATENT PARAMETER Z
update.Z=function(y,n,g,Sigma,tau){
G.n=cbind(G,y,n-1)
Z=rmnorm(N,0,varcov=diag(L))/sqrt(tau)
X.n=cbind(Z,n-1)
M=t(apply(G.n,1,mvm,Sigma=Sigma,g=g))
Z=t(apply(X.n,1,cov.f,Sigma=Sigma))+M
Z[!is.na(y)]=y[!is.na(y)]
return(Z)
}

# SAMPLE REGRESSION COEFFICIENTS FOR Y MODEL
update.g=function(Z,Sigma,mu,rho,tau){
sum.xx=kronecker(solve(Sigma),t(G*tau)%*%G)
sum.xy=kronecker(solve(Sigma),t(G*tau))%*%c(Z)
G.prior=rep(rho,each=L)*diag(T*L)
g.prior=rep(mu,each=L)
V=solve(sum.xx+solve(G.prior))
M=V%*%(sum.xy+solve(G.prior)%*%g.prior)
g=rmnorm(1,mean=M,varcov=V)
return(matrix(g,T,L))
}


# HYPERPARAMETERS FOR COEFFICIENTS
update.mu=function(g,rho){
mean=apply(g,1,mean)
sd=sqrt(rho/D)
mu=rnorm(T,mean=mean,sd=sd)
return(mu)
}

update.rho=function(g,mu,al,be){
W=apply((g-mu)^2,1,sum)
return(1/rgamma(T,shape=D/2+al,rate=W/2+be))
}

# SAMPLE COVARIANCE MATRIX
update.Sigma=function(Z,g,psi,tau){
e=sqrt(tau)*(Z-G%*%g)	
v_0=L+1+exp(psi[1]);  B=diag(L)*exp(sum(psi))+t(e)%*%(e)
return(riwish(v_0+N, B))
}

# SAMPLE HYPERPARAMETERS OF THE COV MATRIX
wpost=function(psi,Sigma,sd){
v_0=L+1+exp(psi[1])
IW=diag(D)*(exp(sum(psi)))
lpost=lndIWishart(v_0, IW,Sigma)+sum(dnorm(c(psi),0,sd,log=TRUE))
return(lpost)
}

update.psi=function(step,psi,Sigma,sd){
cand=c(psi+rnorm(2,0,step))
prob=min(0,wpost(cand,Sigma,sd)-wpost(psi,Sigma,sd))
if(runif(1)<exp(prob)){psi=cand}
return(psi)
}

# SAMPLE LATENT MIXING PARAMETER TAU
m.dist=function(x,Sigma){t(x)%*%solve(Sigma)%*%x}
update.tau=function(Z,nu,g,Sigma){
rgamma(N,shape=(nu+L)/2,rate=apply(Z-G%*%g,1,m.dist,Sigma=Sigma)/2+nu/2)
}


# SAMPLE DEGREES OF FREEDOM PARAMETER NU

lpost.nu=function(nu,tau){sum(dgamma(tau, shape=nu/2, rate=nu/2,log=TRUE))-2*log(nu)}
update.nu=function(nu,tau,step){
cand=1+rgamma(1,shape=nu*step,rate=step)
prob=min(0,lpost.nu(cand,tau)-lpost.nu(nu,tau)-dgamma(cand-1,shape=nu*step,rate=step,log=TRUE)+dgamma(nu,shape=nu*step,rate=step,log=TRUE))
if (runif(1) < exp(prob)){nu = cand}
return(nu)
}

#### STARTING VALUES
if(is.null(init)){
cat("---- initializing starting values ----- ", '\n')
out=list(
Z=NULL,
b=n.mle$par, 
g=rmnorm(T,0,.5*diag(L)),
mu=rep(0,T),
rho=rep(1,T),
Sigma=diag(L),
psi=c(0,0),
tau=rep(1,N),
nu=15
)
}

if(!is.null(init)) out=init

samples=list(
g=array(NA, c(T,L,n.sample)),
b=matrix(NA, nrow=n.sample, ncol=K+1),
mu=matrix(NA, nrow=n.sample, ncol=T),
rho=matrix(NA, nrow=n.sample, ncol=T),
psi=matrix(NA,nrow=n.sample, ncol=2),
Sigma=array(NA,c(L,L,n.sample)),
tau=matrix(NA,n.sample,N),
nu=rep(NA,n.sample)
)

colnames(samples$b)=c(colnames(X),"ln(w)")
colnames(samples$mu)=colnames(samples$rho)=colnames(G)
colnames(samples$g)=c(paste("y",1:L,sep="-"))
rownames(samples$g)=colnames(G)

gibbs=function(out){
out$Z=update.Z(y,n,out$g,out$Sigma,out$tau)
out$tau=update.tau(out$Z,out$nu,out$g,out$Sigma)
out$nu=update.nu(out$nu,out$tau,step=1.5)
out$g=update.g(out$Z,out$Sigma,out$mu,out$rho,out$tau)
out$b=update.b(out$b)
out$mu=update.mu(out$g,out$rho)
out$rho=update.rho(out$g,out$mu,al=1e-3,be=1e-3)
out$Sigma=update.Sigma(out$Z,out$g,out$psi,out$tau)
out$psi=update.psi(.4,out$psi,out$Sigma,sd=1)
return(out)
}

cat("---- burn-in samples ------", '\n')

t=1
while(t<=burn) { ### burn-in chain
out=gibbs(out) 
if (t%%(burn/10)==0)  cat(100*t/burn, "%  ", sep="")
t=t+1
}

cat("                          ", '\n')
cat("---- saving samples ------", '\n')

t=1
while(t<=n.sample){
for (i in 1:thin){
out=gibbs(out)
samples$Sigma[,,t]=out$Sigma
samples$g[,,t]=out$g
samples$b[t,]=unlist(out$b)
samples$mu[t,]=out$mu
samples$rho[t,]=out$rho
samples$psi[t,]=out$psi
samples$tau[t,]=out$tau
samples$nu[t]=out$nu
}
if (t%%(n.sample/10)==0){
cat(100*t/n.sample, "%  ", sep="")
}
t=t+1
}

cat("",'\n')

print.res=function(x){
x=as.mcmc(x)
R=cbind(apply(x,2,mean),apply(x,2,sd),HPDinterval(x),batchSE(x,n.sample/10),heidel.diag(x)[,4])
R=round(R,3)
rownames(R)=colnames(x)
colnames(R)=c("Estimate", "SE", "0.025-HPD", "0.975-HPD","MCMC error", "Heidel-conv")
return(R)
}

cat("",'\n')
cat("---------  PARTY COUNT MODEL PARAMETERS --------- ", '\n')
print(print.res(samples$b))
cat(" ",'\n')

cat("---------  V0TE-SHARE MODEL HYPER-PARAMETERS ---- ", '\n')
print(print.res(samples$m))
cat(" ",'\n')

cat("---------  V0TE-SHARE MODEL PARAMETERS ---------- ", '\n')
cat(" ",'\n')

for (i in 1:L){
cat("gamma",i,sep="-",'\n')
print(print.res(t(samples$g[,i,])))
cat(" ",'\n')
}

return(invisible(samples))

}


y.v <-
function(y){
y=y[!is.na(y)]
A=-a(y)*diag(length(y))
A=cbind(0,A)
A=rbind(A,1)
diag(A)=1
y.star=solve(A)[,length(y)+1]
return(y.star)
}


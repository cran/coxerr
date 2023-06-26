"coxerr" <- function(t,dlt,wuz,method,initbt=rep(0,dim(as.matrix(wuz))[2]-1),
                     derr=1e-6)
{  size <- length(t)
   npred <- dim(as.matrix(wuz))[2]-1

   fit <- .Fortran("coxerr", as.double(t), as.integer(dlt), as.double(t(wuz)),
          as.integer(size), as.integer(npred), as.integer(method),
          as.double(derr), 
          bt=as.double(initbt), va=double(npred^2), succ=integer(1),
          integer(size), double(size), double(npred), double(npred),
	  double(npred), double(npred^2), double(npred+1), double(npred^2),
	  PACKAGE="coxerr")

   list(bt=fit$bt, va=matrix(fit$va,ncol=npred),
        succ=ifelse(fit$succ==1,T,F))
}

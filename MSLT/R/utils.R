##' logLik.mslt 
##'
##' @param object fitted model retuned by \code{\link{fitMSLT}} 
##' @param ... Extra argument
##' @details Returns the log-likelihood and the number of parameters 
##' @export
logLik.mslt<-function(object, ...){
  ret<- -object$opt$objective

  mmppProfiled = 0
  if(object$conf$mmpp==1)mmppProfiled = 1
  
  attr(ret,"df")<-length(object$opt$par) + mmppProfiled
  class(ret)<-"logLik"
  ret
}


##' Print mslt object
##' 
##' @param  x Fitted model retuned by \code{\link{fitMSLT}}
##' @param ... Extra argument
##' @details Print log-likelihood and the main convergence criteria
##' @export
print.mslt<-function(x, ...){
  cat("mslt model: log likelihood is", logLik.mslt(x,...),"Convergence", ifelse(0==x$opt$convergence, "OK\n", "failed\n"))
}


##' print.jit
##' 
##' @param  jit List with the jitter runs
##' @param ... Extra argument
##' @details Print results of the jitter analysis; Difference in parameters, abundance and log likelihood
##' @export
print.jit<-function(jit,...){
  j2 = list()
  counter =1
  for(i in 1:length(jit)){
    if(!is.na(jit[[i]])[1]){
      if(jit[[i]]$opt$convergence==0){
        j2[[counter]] = jit[[i]]
        counter = counter + 1
      }
    }
  }

  run = attributes(jit)$run

  scale = 1/((4*3.14159265)*exp(run$pl$log_kappa*2))

  maxabsdiff <- apply(abs(do.call(cbind, lapply(j2, function(f)unlist(f$pl)-unlist(run$pl)))),1,max)
  maxlist <- utils::relist(maxabsdiff, run$pl)
  ret <- as.data.frame(unlist(lapply(maxlist,function(x)if(length(x)>0)max(x) else NULL)))
  logLik <- max(abs(unlist(lapply(j2, logLik))-logLik(run)))
  abundance <- max(unlist(lapply(j2, function(f)abs(f$rl$logAbundance-run$rl$logAbundance))))
  ret <- rbind(ret,  logLik=logLik, abundance = abundance)
  names(ret) <- "max(|delta|)"
  ret[which(rownames(ret)=="x_intensity"),] = ret[which(rownames(ret)=="x_intensity"),]/sqrt(scale[1]) *exp(run$pl$log_sigma[1])
  ret[which(rownames(ret)=="x_size"),] = ret[which(rownames(ret)=="x_size"),]/sqrt(scale[2]) *exp(run$pl$log_sigma[2])

  print(ret)
  print(paste0("Proportion convergence: ", length(j2)/length(jit)))
}



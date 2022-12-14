## used by stats::AIC()
##'
##' @param object
##' @details
##' @export
logLik.lineTransectLGPP<-function(object, ...){
  ret<- -object$opt$objective

  attr(ret,"df")<-length(object$opt$par)
  class(ret)<-"logLik"
  ret
}


##' Print lineTransectLGPP object
##' @method print lineTransectLGPP
##' @param  x
##' @details Print log-likelihood and the main convergence criteria
##' @export
print.lineTransectLGPP<-function(x, ...){
  cat("lineTransectLGPP model: log likelihood is", logLik.lineTransectLGPP(x,...),"Convergence", ifelse(0==x$opt$convergence, "OK\n", "failed\n"))
}


##' Print jit
##' @method print jit
##' @param  jj
##' @details Print jit
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
  maxlist <- relist(maxabsdiff, run$pl)
  ret <- as.data.frame(unlist(lapply(maxlist,function(x)if(length(x)>0)max(x) else NULL)))
  logLik <- max(abs(unlist(lapply(j2, logLik))-logLik(run)))
  abundance <- max(unlist(lapply(j2, function(f)abs(logAbundance(f)-logAbundance(run)))))
  ret <- rbind(ret,  logLik=logLik, abundance = abundance) #TODO: Add abundance
  names(ret) <- "max(|delta|)"
  ret[which(rownames(ret)=="x_intensity"),] = ret[which(rownames(ret)=="x_intensity"),]/sqrt(scale[1]) *exp(run$pl$log_sigma[1])
  ret[which(rownames(ret)=="x_size"),] = ret[which(rownames(ret)=="x_size"),]/sqrt(scale[2]) *exp(run$pl$log_sigma[2])

  print(ret)
  print(paste0("Proportion convergence: ", length(j2)/length(jit)))
}



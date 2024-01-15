#' Fit the LGPP model
#' @description Fit the LGPP model
#' @param data data needed
#' @param par data needed
#' @param conf data needed
#' @param rel.tol data needed
#' @param map data needed
#' @param ... sent to sdreport
#' @return fitted LGPP model
#' @export
#' @examples
fitLGPP = function(data,par,conf,rel.tol=1e-10,map = setMap(conf, par),...){

  #Estimating the model and extract results-------------
  startTime <- Sys.time()
  if(conf$mmpp==1){
    obj <- RTMB::MakeADFun(mslt, par, random=c("x_intensity","x_size"), profile = c("log_c_mmpp"), map = map)	
  }else{
    obj <- RTMB::MakeADFun(mslt, par, random=c("x_intensity","x_size"), map = map)	
  }
  lower = list()
  lower$log_c_mmpp = -10#NB, set lower boundary on jump in MMPP

  fit = tryCatch({
    opt <- nlminb(obj$par, obj$fn, obj$gr,control = list(rel.tol = rel.tol,trace = 1,iter.max = 1000,eval.max = 2000),lower = lower)
    rep = sdreport(obj, ...)
    pl = as.list(rep,"Est")
    plSd = as.list(rep,"Std")
    list(obj = obj, opt = opt,rep = rep, map = map, conf = conf, data = data,pl = pl, plSd = plSd)},
  error=function(cond) {
    print("NB!!! Optization failed!")
    NA
  })

  endTime = Sys.time()
  timeUsed = endTime - startTime
  print(timeUsed)

  class(fit) = "lineTransectLGPP"
  return(fit)
}




#' @param
#' @return fitted model
#' @export
#' @examples
jit = function(run, sd = 0.2, nojit = 50,ncores = 2){
  par = setPar(run$data,run$conf)
  parTmp = setPar(run$data,run$conf)
  parv <- unlist(par)
  parsTmp <- lapply(1:nojit, function(i)relist(parv+rnorm(length(parv),sd=sd), par))
  pars =  lapply(parsTmp, function(p){
    p
  })


  if(ncores>1){
    cl <- makeCluster(ncores) #set up nodes
    on.exit(stopCluster(cl)) #shut it down
    lib.ver <- dirname(path.package("mslt"))
    clusterExport(cl, varlist=c("lib.ver","run"), envir=environment())
    clusterEvalQ(cl, {library(mslt, lib.loc=lib.ver)})
    runs <- parLapply(cl, pars, function(p){
      rr = fitLGPP(run$data,  p,run$conf)
      FreeADFun(rr$obj)#Free memory from C-side
      rr
    })
  } else {
    runs <- lapply(pars, function(p){
      rr = fitLGPP(run$data,  p, run$conf)
      FreeADFun(rr$obj)#Free memory from C-side
      rr
      })
  }
  attr(runs,"run") <- run
  class(runs) = "jit"

  return(runs)
}


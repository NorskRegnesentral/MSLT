#' @param
#' @return fitted model
#' @useDynLib mslt
#' @export
#' @examples
fitLGPP = function(data,par,conf,rel.tol=1e-10,...){

  map = setMap(conf, par);

  #Estimating the model and extract results-------------
  startTime <- Sys.time()
  obj <- MakeADFun(data, par, random=c("x_intensity","x_size"), DLL="mslt", map = map)

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
    scale =1/((4*3.14159265)*exp(p$log_kappa)*exp(p$log_kappa))
    p$x_intensity = p$x_intensity*sqrt(scale[1])
    p$x_size = p$x_size*sqrt(scale[2])
    p
  })


  if(ncores>1){
    cl <- makeCluster(ncores) #set up nodes
    on.exit(stopCluster(cl)) #shut it down
    lib.ver <- dirname(path.package("niseEst"))
    clusterExport(cl, varlist=c("lib.ver","run"), envir=environment())
    clusterEvalQ(cl, {library(niseEst, lib.loc=lib.ver)})
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







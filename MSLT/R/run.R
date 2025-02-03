#' fitMSLT
#' 
#' @param data data list returned by \code{\link{setData}}
#' @param par parameter list  returned by \code{\link{setPar}}
#' @param conf model configurations list  returned by \code{\link{defConf}}
#' @param rel.tol rel.tol in nlminb
#' @param map map object used by TMB to couple and remove parameters
#' @param ... additional parameters passed to sdreport 
#' @return fitted model
#' @useDynLib mslt
#' @export
fitMSLT = function(data,par,conf,rel.tol=1e-10,map = setMap(conf, par),...){

  #Estimating the model and extract results-------------
  startTime <- Sys.time()
  if(conf$mmpp==1){
    obj <- TMB::MakeADFun(data, par, random=c("x_intensity","x_size"), profile = c("log_c_mmpp"),DLL="mslt", map = map)	
  }else{
    obj <- TMB::MakeADFun(data, par, random=c("x_intensity","x_size"),DLL="mslt", map = map)	
  }
  lower = list()
  lower$log_c_mmpp = -10#NB, set lower boundary on jump in MMPP

  fit = tryCatch({
    opt <- stats::nlminb(obj$par, obj$fn, obj$gr,control = list(rel.tol = rel.tol,trace = 1,iter.max = 1000,eval.max = 2000),lower = lower)
    rep = TMB::sdreport(obj, ...)
    pl = as.list(rep,"Est")
    plSd = as.list(rep,"Std")
    rl = as.list(rep, what = "Est", report = TRUE)
    rlSd = as.list(rep, what = "Std", report = TRUE)
    list(obj = obj, opt = opt,rep = rep, map = map, conf = conf, data = data,pl = pl, plSd = plSd, rl = rl, rlSd = rlSd)},
  error=function(cond) {
    print("NB!!! Optization failed!")
    NA
  })

  endTime = Sys.time()
  timeUsed = endTime - startTime
  print(timeUsed)

  class(fit) = "mslt"
  return(fit)
}




#' jit
#' 
#' @param run Fitted model returned by \code{\link{fitMSLT}}
#' @param sd Standard deviation in jitter analysis
#' @param nojit Number of sets of starting values
#' @param ncores Number of cores used in jitter analysis
#' @return Fitted model
#' @export
jit = function(run, sd = 0.2, nojit = 50,ncores = 2){
  par = setPar(run$data,run$conf)
  parTmp = setPar(run$data,run$conf)
  parv <- unlist(par)
  parsTmp <- lapply(1:nojit, function(i)utils::relist(parv+stats::rnorm(length(parv),sd=sd), par))
  pars =  lapply(parsTmp, function(p){
    scale =1/((4*3.14159265)*exp(p$log_kappa)*exp(p$log_kappa))
    p$x_intensity = p$x_intensity*sqrt(scale[1])
    p$x_size = p$x_size*sqrt(scale[2])
    p
  })


  if(ncores>1){
    cl <- parallel::makeCluster(ncores) #set up nodes
    on.exit(parallel::stopCluster(cl)) #shut it down
    lib.ver <- dirname(path.package("mslt"))
    parallel::clusterExport(cl, varlist=c("lib.ver","run"), envir=environment())
    parallel::clusterEvalQ(cl, {library(mslt)})
    runs <- parallel::parLapply(cl, pars, function(p){
      rr = fitMSLT(run$data,  p,run$conf)
      TMB::FreeADFun(rr$obj)#Free memory from C-side
      rr
    })
  } else {
    runs <- lapply(pars, function(p){
      rr = fitMSLT(run$data,  p, run$conf)
      TMB::FreeADFun(rr$obj)#Free memory from C-side
      rr
      })
  }
  attr(runs,"run") <- run
  class(runs) = "jit"

  return(runs)
}







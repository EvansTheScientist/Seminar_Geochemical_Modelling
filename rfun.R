## This script contains the evaporated rainwater as initial solution and the Northern zone as the target solution

######################################################################################
##' This function performs a combinatorial analysis, creating all
##' possible PHREEQC models with a number "len" of equilibrium phases
##' across a pool of primary and secondary minerals.
##'
##' Requires packages doParallel and foreach (for parallelization),
##' but they should already be installed as dependency of RedModRphree
##' @title Perform combinatorial analysis of a given pool of
##'     equilibrium minerals
##' @param initsol the base PHREEQC script
##' @param primary character vector containing minerals which are
##'     inserted in the models as primary minerals
##' @param secondary character vector containing minerals which are
##'     inserted in the models as secondary minerals (i.e., initial 0
##'     amount)
##' @param len number of phases in each model
##' @param procs integer, how many CPUs can we use? Defaults to 4
##' @return a list containing the parsed results (list of many blocks)
##' @author MDL
DoCombSim <- function(initsol, db, primary, secondary, len, procs=4L) {
  require(parallel)
  require(doParallel)
  require(foreach)
  
  ## create all combinations of primary and secondary phases of length len
  phases <- c(primary, secondary)
  combs <- combn(phases, len, FUN = NULL, simplify=FALSE)
  
  cat(":: Going to do ", length(combs), " simulations\n")
  
  ## create the phreeqc scripts
  .addPhaseComb <- function(x) {
    inp <- initsol
    for (phase in x) {
      inp <- AddProp(inp, name=phase, values=ifelse(phase %in% primary, "0.0 2", "0.0 0"), cat="pphases")
    }
    return(inp)
  }
  
  ## apply this function to all entries
  biginp <- lapply(combs, .addPhaseComb)
  
  ## workhorse function to run simulations
  .runPQC <- function(input) {
    phreeqc::phrSetOutputStringsOn(TRUE)
    phreeqc::phrRunString(input)
    tmpout <- phreeqc::phrGetOutputStrings()
    res <- RedModRphree::ReadOut(tmpout)[[1]]
    return(res)
  }
  
  if (procs > 1) {
    if (Sys.info()[["sysname"]]=="Windows") {
      ThisRunCluster <- parallel::makePSOCKcluster(procs)
    } else {
      ThisRunCluster <- parallel::makeForkCluster(procs)
    }
    
    doParallel::registerDoParallel(ThisRunCluster)
    cat(":: Registered default doParallel cluster with ", procs, "nodes")
    parallel::clusterCall(cl=ThisRunCluster, phreeqc::phrLoadDatabase, db)
    msg(":: Database loaded on each worker")
    
    res <- foreach::foreach(i = seq_along(biginp)) %dopar% 
      .runPQC(biginp[[i]])
    cat("[ DONE ]\n")
    parallel::stopCluster(ThisRunCluster)
    
  } else {
    ## revert to sequential computation
    cat(":: Firing up PHREEQC onsingle CPU...")
    res <- lapply(biginp, .runPQC)
    cat("[ DONE ]\n")
    
  }
  
  return(res)
}

##' Computes a specific metric allowing for selection of the
##' components to be included
##'
##' @title Compute metric selecting components
##' @param data the matrix or data.frame containing all the results
##'     from the PHREEQC simulations. Its columns need to be named!
##' @param target the named vector with the target concentrations
##' @param FUN the name of the metric function
##' @param comp optional, a char vector with the names of the
##'     components. If unspecified, all components are selected
##' @param ... further parameter passed to FUN, such as "na.rm"
##' @return numeric vector with the computed metric
##' @author Marco
ComputeMetric <- function(data, target, FUN="rmse", comp=colnames(data), ...){
  ## find the metric function
  .Fun <- match.fun(FUN)
  
  ## retain only the columns given as argument
  tmp <- subset(data, select = comp)
  
  cvec <- target[comp]
  ## compute stuff using apply
  res <- apply(tmp, 1, function(x) .Fun(cvec, x, ...))
  return(res)
}

##### Filtering

Filter <- function(lin, delta=0.5) {
  excluded <- sapply(lin, function(x) any(abs(x$pphases$delta)>delta))
  out <- lin[which(!excluded)]
  return(out)
}

Filter2 <- function(lin, delta=0.5) {
  excluded <- sapply(lin, function(x) any(abs(x$pphases$delta)>delta))
  return(which(excluded))
}

FilterAll <- function(lin, delta=0.5) {
  retain <- sapply(lin, function(x) all(abs(x$pphases$delta)<delta))
  out <- which(!retain)
  return(out)
}

## Some metrics
rrmse <- function(y_true, y_pred, na.rm=TRUE)
  sqrt(mean(((y_true - y_pred)/y_true)^2, na.rm = na.rm))

## mean absolute percent error
rmape <- function(y_true, y_pred, na.rm=TRUE)
  mean(abs((y_true- y_pred)/y_true), na.rm = na.rm) * 100

##Additional errors relative Mean absolute error (mae)
rmae <- function(y_true, y_pred, na.rm=TRUE)
  mean(abs((y_true- y_pred)/y_true), na.rm = na.rm)


##
PlotComb <- function(res, samples, comp, ...) {
  if (missing(comp)) {
    comp <- intersect(colnames(res), colnames(samples))
  }
  tmp <- subset(res, select=comp)
  sam <- subset(samples, select=comp)
  
  mins <- apply(sam, 2, min, na.rm=TRUE)
  maxs <- apply(sam, 2, max, na.rm=TRUE)
  meds <- apply(sam, 2, median, na.rm=TRUE)
  meas <- apply(sam, 2, mean, na.rm=TRUE)
  
  colors <- heat.colors(nrow(tmp))
  
  out <- barplot(tmp, beside=TRUE, ylab="", log="y",
                 col=colors, las=1, ...)
  for (i in seq_along(comp)) {
    rect(out[1,i]-0.6, mins[i],out[nrow(out),i]+0.6, maxs[i],col=rgb(0,0,1.0,alpha=0.5))
    segments(out[1,i]-0.6, meds[i], out[nrow(out),i]+0.6, meds[i],col="red", lwd=2, lty="dashed")
    segments(out[1,i]-0.6, meas[i], out[nrow(out),i]+0.6, meas[i],col="grey", lwd=2, lty="dotted")
  }
}




# the operating model was modified from R package "fishdynr"
# Two major changes were done compared with the original package: 
# 1) the age of new individual is 0 instead of 1 (original setting);
# 2) observation process for length-at-age data is added 
# R version 3.5.3
library(fishdynr)

#tincr = 1/12;K.mu = 0.5; K.cv = 0.1;Linf.mu = 80; Linf.cv = 0.1;ts = 0; C = 0.85;LWa = 0.01; LWb = 3;Lmat = 0.5*Linf.mu; wmat = Lmat*0.2;rmax = 10000; beta = 1;repro_wt = c(0,0,0,1,0,0,0,0,0,0,0,0);M = 0.7; harvest_rate = M;L50 = 0.25*Linf.mu; wqs = L50*0.2;bin.size = 1;timemin = 0; timemax = 20; timemin.date = as.Date("1980-01-01");N0 = 10000;fished_t = seq(17+1/12,19,tincr);lfqFrac = 1;progressBar = TRUE


lfqSim <- function(
  tincr = 1/12,
  K.mu = 0.5, K.cv = 0.1,
  Linf.mu = 80, Linf.cv = 0.1,
  ts = 0, C = 0.85,
  LWa = 0.01, LWb = 3,
  Lmat = 0.5*Linf.mu, wmat = Lmat*0.2,
  rmax = 10000, beta = 1,
  repro_wt = c(0,0,0,1,0,0,0,0,0,0,0,0),
  M = 0.7, harvest_rate = M, 
  L50 = 0.25*Linf.mu, wqs = L50*0.2,
  bin.size = 1,
  timemin = 0, timemax = 20, timemin.date = as.Date("1980-01-01"),
  N0 = 10000,
  fished_t = seq(17+1/12,19,tincr),
  lfqFrac = 1,
  progressBar = TRUE
){
  
  # times
  timeseq = seq(from=timemin+tincr, to=timemax, by=tincr)
  if(!zapsmall(1/tincr) == length(repro_wt)) stop("length of repro_wt must equal the number of tincr in one year")
  repro_wt <- repro_wt/sum(repro_wt)
  repro_t <- rep(repro_wt, length=length(timeseq)) 
  # repro_t <- seq(timemin+repro_toy, timemax+repro_toy, by=1)
  
  # make empty lfq object and length-at-age boject
  lfq <- vector(mode="list", length(timeseq))
  names(lfq) <- timeseq
  LA<-vector(mode="list", length(timeseq))
  names(LA) <- timeseq
  # Estimate tmaxrecr
  tmaxrecr <- (which.max(repro_wt)-1)*tincr
  
  # mean phiprime
  phiprime.mu = log10(K.mu) + 2*log10(Linf.mu)
  
  
  
  # required functions ------------------------------------------------------
  date2yeardec <- function(date){as.POSIXlt(date)$year+1900 + (as.POSIXlt(date)$yday)/365}
  yeardec2date <- function(yeardec){as.Date(strptime(paste(yeardec%/%1, ceiling(yeardec%%1*365+1), sep="-"), format="%Y-%j"))}
  
  make.inds <- function(
    id=NaN, A = 0, L = 0, W=NaN, mat=0,
    K = K.mu, Winf=NaN, Linf=NaN, phiprime=NaN,
    F=NaN, Z=NaN, 
    Fd=0, alive=1
  ){
    inds <- data.frame(
      id = id,
      A = A,
      L = L,
      W = W,
      Lmat=NaN,
      mat = mat,
      K = K,
      Linf = Linf,
      Winf = Winf,
      phiprime = phiprime,
      F = F,
      Z = Z,
      Fd = Fd,
      alive = alive
    )
    lastID <<- max(inds$id)
    return(inds)
  }
  
  
  express.inds <- function(inds){
    inds$Linf <- Linf.mu * rlnorm(nrow(inds), 0, Linf.cv)
    inds$Winf <- LWa*inds$Linf^LWb
    # inds$K <- 10^(phiprime.mu - 2*log10(inds$Linf)) * rlnorm(nrow(inds), 0, K.cv)
    inds$K <- K.mu * rlnorm(nrow(inds), 0, K.cv)
    inds$W <- LWa*inds$L^LWb
    inds$phiprime <- log10(inds$K) + 2*log10(inds$Linf)
    inds$Lmat <- rnorm(nrow(inds), mean=Lmat, sd=wmat/diff(qnorm(c(0.25, 0.75))))
    return(inds)
  }
  
  
  
  grow.inds <- function(inds){
    # grow
    L2 <- dt_growth_soVB(Linf = inds$Linf, K = inds$K, ts = ts, C = C, L1 = inds$L, t1 = t-tincr, t2 = t)
    # update length and weight
    inds$L <- L2
    inds$W <- LWa*inds$L^LWb
    # age inds
    inds$A <- inds$A + tincr
    return(inds)
  }
  
  mature.inds <- function(inds){
    # p <- pmat_w(inds$L, Lmat, wmat) # probability of being mature at length L
    # p1t <- 1-((1-p)^tincr)
    # inds$mat <- ifelse(runif(nrow(inds)) < p1t | inds$mat == 1, 1, 0)
    inds$mat <- ifelse((inds$L > inds$Lmat | inds$mat == 1), 1, 0)
    return(inds)
  }
  
  death.inds <- function(inds){
    pSel <- logisticSelect(inds$L, L50, wqs)
    inds$F <- pSel * Fmax
    inds$Z <- M + inds$F
    pDeath <- 1 - exp(-inds$Z*tincr)
    dead <- which(runif(nrow(inds)) < pDeath)
    # determine if natural or fished
    if(length(dead) > 0){
      inds$alive[dead] <- 0
      tmp <- cbind(inds$F[dead], inds$Z[dead])
      # Fd=1 for fished individuals; Fd=0, for those that died naturally
      Fd <- apply(tmp, 1, FUN=function(x){sample(c(0,1), size=1, prob=c(M/x[2], x[1]/x[2]) )})
      inds$Fd[dead] <- Fd
      rm(tmp)
    }
    return(inds)
  }
  
  remove.inds <- function(inds){
    dead <- which(inds$alive == 0)
    if(length(dead)>0) {inds <- inds[-dead,]}
    return(inds)
  }
  
  reproduce.inds <- function(inds){
    # reproduction can only occur of population contains >1 mature individual
    if(repro > 0 & sum(inds$mat) > 0){
      #calc. SSB
      SSB <- sum(inds$W*inds$mat)
      n.recruits <- ceiling(srrBH(rmax, beta, SSB) * repro)
      # make recruits 
      offspring <- make.inds(
        id = seq(lastID+1, length.out=n.recruits)
      )
      # express genes in recruits
      offspring <- express.inds(offspring)	
      #combine all individuals
      inds <- rbind(inds, offspring)
    }	
    return(inds)	
  }
  
  record.inds <- function(inds, ids=1:10, rec=NULL){
    if(is.null(rec)) {
      rec <- vector(mode="list", length(ids))
      names(rec) <- ids
      inds <- inds
    } else {
      ids <- as.numeric(names(rec))
    }
    if(length(rec) > 0) {
      inds.rows.rec <- which(!is.na(match(inds$id, ids)))
      if(length(inds.rows.rec) > 0){
        for(ii in inds.rows.rec){
          match.id <- match(inds$id[ii], ids)
          if(is.null(rec[[match.id]])) {
            rec[[match.id]] <- inds[ii,]
          } else {
            rec[[match.id]] <- rbind(rec[[match.id]], inds[ii,])
          }
        }
      }
    }
    rec
  }
  
  
  
  
  # run model ---------------------------------------------------------------
  
  # Initial population
  lastID <- 0
  inds <- make.inds(
    id=seq(N0)
  )
  inds <- express.inds(inds)
  
  # results object
  res <- list()
  res$pop <- list(
    dates = yeardec2date( date2yeardec(timemin.date) + (timeseq - timemin) ),
    N = NaN*timeseq,
    B = NaN*timeseq,
    SSB = NaN*timeseq,
    catch = NaN*timeseq
  )
  
  # simulation

  if(progressBar) pb <- txtProgressBar(min=1, max=length(timeseq), style=3)
  for(j in seq(timeseq)){
    t <- timeseq[j]
    
    # harvest rate applied? lfq sampled?
    if(length(fished_t) == 0){
      Fmax <- 0
      lfqSamp <- 0
    } else {
      if(min(sqrt((t-fished_t)^2)) < 1e-8){
        Fmax <- harvest_rate
        lfqSamp <- 1
      } else {
        Fmax <- 0
        lfqSamp <- 0
      }
    }
    
    repro <- repro_t[j]
    
    # population processes
    inds <- grow.inds(inds)
    inds <- mature.inds(inds)
    inds <- reproduce.inds(inds)
    inds <- death.inds(inds)
    if(lfqSamp){
      samp <- try( sample(nrow(inds), ceiling(sum(inds$Fd)*lfqFrac), prob = inds$Fd), silent = TRUE)
      if(class(samp) != "try-error"){lfq[[j]] <-inds$L[samp]}
      if(class(samp) != "try-error"){LA[[j]] <- cbind(inds$L[samp],inds$A[samp])}
      if(class(samp) != "try-error"){res$pop$catch[j] <- sum(inds$W[samp])}
      rm(samp)
    }
    
    inds <- remove.inds(inds)
    
    # update results
    res$pop$N[j] <- nrow(inds)
    res$pop$B[j] <- sum(inds$W)
    res$pop$SSB[j] <- sum(inds$W*inds$mat)
    if(progressBar) setTxtProgressBar(pb, j)
    
  }
  if(progressBar) close(pb)
  
  
  
  # Export data -------------------------------------------------------------
  
  # Trim and Export 'lfq'
  lfq2 <- lfq[which(sapply(lfq, length) > 0)]
  
  
  # binned version of lfq
  dates <- yeardec2date( date2yeardec(timemin.date) + (as.numeric(names(lfq2)) - timemin) )
  Lran <- range(unlist(lfq2))
  Lran[1] <- floor(Lran[1])
  Lran[2] <- (ceiling(Lran[2])%/%bin.size + ceiling(Lran[2])%%bin.size + 1) * bin.size
  bin.breaks <- seq(Lran[1], Lran[2], by=bin.size)
  bin.mids <- bin.breaks[-length(bin.breaks)] + bin.size/2
  res$lfqbin <- list(
    sample.no = seq(bin.mids),
    midLengths = bin.mids,
    dates = dates,
    catch = sapply(lfq2, FUN = function(x){
      hist(x, breaks=bin.breaks, plot = FALSE, include.lowest = TRUE)$counts
    })
  )
  res$LA<-LA
  # record mean parameters
  res$growthpars <- list(
    K = K.mu,
    Linf = Linf.mu,
    C = C,
    ts = ts,
    phiprime = phiprime.mu,
    tmaxrecr = tmaxrecr
  )
  
  return(res)
  
} # end of function




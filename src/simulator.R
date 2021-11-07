#assume fit.ntv is 

getFitLTV <- function( prok.df ){
  #fit.ltv <- lm( lm ~ lr + g, data = prok.df )
  fit.ltv <- loess( lm ~ lr , data=prok.df )
  fit.ltv
}

assignCountsForVal <- function(lrval, cl.id, b.id, start.counter=1, dalph=1, nk=NULL, fit.ltv=NULL,... ){
  require(gtools)
  if(lrval==-Inf){
    ids <- paste0( "cluster>", cl.id, ">", NA)
    data.frame( ids=ids, abun=NA )
  } else {
    
    if(is.null(nk)){
      stopifnot(!is.null( fit.ltv ))
      #mod <- do.call( rbind,  lapply( seq_along(lrval), function(x) fit.ltv$model[1,,drop=FALSE] ) )
      #mod[,2] <- lrval
      #ltv <- predict( fit.ltv, mod )
      ltv <- predict( fit.ltv, newdata= lrval )
      
      nk <- ceiling( exp( ltv ) )
    }
    rd <- rdirichlet( n =1, alpha =  dalph/seq(nk)  )
    rmul <- rmultinom( n=1, size=exp(lrval), prob=rd)
    print("---")
    print(start.counter)
    print(length(rmul))
    ids <- paste0( "cluster>", cl.id, ">", b.id, ">", seq(start.counter, start.counter + length(rmul)-1))
    data.frame( ids=ids, abun=rmul )
    
  }
  
}

assignCountsForVec <- function(lrvec, fit.ltv, b.id, cl.id=1, ... ){
  #mod <- do.call( rbind,  lapply( seq_along(lrvec), function(x) fit.ltv$model[1,,drop=FALSE] ) )
  #mod[,2] <- lrvec
  #ltv <- predict( fit.ltv, mod )
  ltv <- predict( fit.ltv, newdata = lrvec )
  ltv[is.na(ltv)] <- 1
  ntv <- ceiling(exp(ltv))
  cltv <- cumsum(ntv)
  lst <- do.call( rbind, lapply( seq_along(lrvec), function(i){
    if( i==1 ){
      sc <- 1
    } else {
      sc <- cltv[i-1]+1
    }
    assignCountsForVal( lrval=lrvec[i], nk=ntv[i], cl.id=cl.id, b.id=i, start.counter = sc, ...  )
  } ) )
  
  lst[ which(lst$abun>0), ]
  
}

assignTaxFor <- function(cl.id, lr, fit.ltv, nb=1, ... ){
  
  #bio
  u <- runif( n=nb, min = (1e-5), max=1-(1e-5) ); u <- u/sum(u)
  lrumat <- do.call( cbind, lapply( lr, function(lr.inst){
    log( rmultinom(n=1, size=exp( lr.inst ), prob=u) )
  } ) ) 
  
  #tech
  swiseTaxForFeature <- apply( lrumat, 2, function(lruvec){
    assignCountsForVec( lruvec, fit.ltv = fit.ltv, cl.id = cl.id, ...  ) 
  } )
  
  Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "ids", all.x = TRUE, all.y=TRUE), swiseTaxForFeature )
  
}

genGroupData <- function( lrmat, nbvec=rep(1, nrow(lrmat)), fit.ltv, cl.ids=seq(nrow(lrmat)), ncore=4, fnvec = c( "assignTaxFor", "assignCountsForVal", "assignCountsForVec", "getFitLTV" ), ...  ){
  require(doParallel)
  cl <- makeCluster( ncore )
  registerDoParallel(cl)
  
  ivec <- seq( nrow( lrmat ) )
  dat <- foreach(zid=ivec, .combine = rbind, .export = fnvec) %dopar% {
    print( zid )
    assignTaxFor( cl.id = zid, lr=lrmat[zid,], fit.ltv=fit.ltv, nb = nbvec[zid]  )
  }
  dat[is.na(dat)] <- 0
  
  stopCluster(cl)
  
  dat
}


gen.lrobj <- function( df, K=50, ng=20 , 
                       med.thresh=median( df$lr ), 
                       lr.cpi=.05, 
                       ud.thresh=.5, 
                       upwards=TRUE, symm=FALSE ){
  
  z <- df$lr
  z<- z[!is.na(z)]
  z <- z[ z!=0 ]
  med.thresh <- median(z)
  df.lo <- z[ which( z < med.thresh ) ]
  df.hi <- z[ which( z>med.thresh ) ]
  contr.lr.mu <- sample( df.lo, size = K )
  case.lr.mu <- sample( df.hi, size = K )
  nc <- floor( lr.cpi*K )
  ru <- NA
  
  if(nc >0){
    toch.i <- sample( seq(K), nc, replace=FALSE )
    
    if(upwards){
      case.lr.mu <- contr.lr.mu
      case.lr.mu[toch.i] <- sample( df.hi, nc )  
    } else {
      contr.lr.mu[toch.i] <- case.lr.mu[toch.i]
      case.lr.mu <- contr.lr.mu
      case.lr.mu[toch.i] <- sample( df.lo, nc )  
    }
    print( case.lr.mu)
    print( contr.lr.mu)
    print( case.lr.mu - contr.lr.mu )
    if(symm){
      #choose downers, swap means.
      ru <- runif( nc )
      down <- which( ru < ud.thresh )
      tmp <- contr.lr.mu
      contr.lr.mu[toch.i[down]] <- case.lr.mu[ toch.i[down] ]
      case.lr.mu[toch.i[down]] <- tmp[ toch.i[down] ]  
    }
  }
  
  mu <- aggregate( df$lr, list( df$samples ), mean, na.rm=TRUE )$x
  sd <- aggregate( df$lr, list( df$samples ), sd, na.rm=TRUE )$x
  fit <- lm( sd ~ mu )
  
  mu <- .5*(contr.lr.mu  + case.lr.mu)
  
  
  contr.dat <- do.call( cbind, lapply( seq(ng),function(x) contr.lr.mu )  )
  sds.contr <- predict( fit, newdata=data.frame( "mu"= contr.lr.mu) )
  sds.contr[!is.finite(sds.contr)] <- 1
  contr.noise <- do.call( rbind, lapply( seq(K),function(g) rnorm(  n=(ng), 
                                                                    mean=0, sd=sds.contr[g]  ) )  ) 
  contr.dat <- contr.dat + contr.noise
  
  case.dat <- do.call( cbind, lapply( seq(ng),function(x) case.lr.mu )  ) 
  sds.case <- predict( fit, newdata=data.frame( "mu"=case.lr.mu) )
  sds.case[!is.finite(sds.case)] <- 1
  case.noise <- do.call( rbind, lapply( seq(K),function(g) rnorm(  n=(ng), 
                                                                    mean=0, sd=sds.case[g]  ) )  ) 
  case.dat <- case.dat + case.noise
  all.mat <- cbind( contr.dat, case.dat )
  
  mat <- all.mat
  mdlr <- max( df$lr )
  indx <- which( mat > mdlr ) #those that exceed what we have observed ; very few
  print( paste( "Total indx out of range: ", length( indx ) ) )
  if( length( indx ) >0 ){
    sds <- predict( fit, newdata=data.frame("mu"=all.mat[indx]) ); sds[!is.finite(sds)] <- 1
    mat[indx ] <- truncnorm::rtruncnorm( length(indx), mean = all.mat[indx], sd = sds, a = -Inf, b=mdlr   )   
  }
  colnames( mat ) <- paste0( "sample", seq(ncol(mat))  )
  group <- factor( rep( c(1,2), each=ng ) )
  pd <- data.frame( group ); rownames(pd) <- colnames(mat)
  
  obj <- ExpressionSet( assayData = mat, phenoData = AnnotatedDataFrame(pd)  )
  
  list( "obj"=obj, 
        "sig.dat"=list( "contr.dat"=contr.dat, "case.dat"=case.dat ), 
        "noise.dat"=list("contr.noise"=contr.noise, "case.noise"=case.noise ),
        "fi"=indx, #those resampled
        "sds.contr"=sds.contr,
        "sds.case"=sds.case,
        "mdlr"=mdlr,
        "lrfc"=case.lr.mu - contr.lr.mu,
        "contr.lr.mu"=contr.lr.mu,
        "case.lr.mu"=case.lr.mu,
        "ru"=ru,
        "ud.thresh"=ud.thresh
  )
  
}

simDatProk <- function( prok.df, lr.cpi=.05, nb.cpi=0, 
                        nb.sample.vals=seq(2,10), 
                        force.same=FALSE, force.diff=TRUE, ... ){
  
  if(force.same){
    stopifnot( lr.cpi == nb.cpi ) #force same means lr.cpi must equal nb.cpi cannot be true! 
  }
  fit.ltv <- getFitLTV( prok.df )
  lrobj <- gen.lrobj( df=prok.df, lr.cpi = lr.cpi, ...  )
  obj <- lrobj$obj
  group <- pData(obj)$group
  contr.gr <- which( group==1 )
  nbvec.contr <- rep( 1, nrow( obj) )
  inf.contr <- genGroupData( lrmat = exprs( obj[,contr.gr,drop=FALSE] ), 
                             nbvec = nbvec.contr, 
                             fit.ltv = fit.ltv  )
  gd <- genGroupData( lrmat = exprs( obj ), nbvec = rep( 1, nrow( obj ) ), fit.ltv = fit.ltv  )
  
  
  #nb change
  nbvec.case <- rep( 1,nrow(obj) )
  nc <- ceiling( nb.cpi*nrow( obj ) )
  if( nc > 0 ){
    if( force.same ){
      toch <- which( lrobj$lrfc != 0 )
      nc <- sum( lrobj$lrfc!=0 )
    } else {
      if(force.diff){
        toch <- sample( seq(nrow(obj))[-which( lrobj$lrfc != 0 )], nc )  
      } else {
        toch <- sample( nrow(obj), nc )    
      }
      
    }
    
    nbvec.case[toch] <- sample( nb.sample.vals, nc ,replace=TRUE )
  } 
  
  inf.case <- genGroupData( lrmat = exprs( obj[,-contr.gr,drop=FALSE] ), 
                            nbvec = nbvec.case, 
                            fit.ltv = fit.ltv  )
  inf.mat <- merge( inf.contr, inf.case, by="ids", all.x = TRUE, all.y=TRUE ) 
  inf.mat[is.na(inf.mat) ] <- 0
  rownames( inf.mat ) <- inf.mat$ids
  inf.mat <- inf.mat[,-1] #remove ids column
  colnames( inf.mat ) <- paste0( "sample", seq(ncol(inf.mat)) )
  
  fd <- as.data.frame( do.call( rbind, strsplit( rownames(inf.mat), ">" ) )[,-1])
  colnames( fd ) <- c( "cluster", "subtype", "error.variant")
  rownames(fd ) <- rownames( inf.mat )
  
  inf.obj <- ExpressionSet( assayData = as.matrix( inf.mat), 
                            phenoData = phenoData( obj ) , 
                            featureData = AnnotatedDataFrame(fd) )
  
  list(
    "lrobj"=lrobj, 
    "inf.obj"=inf.obj,
    "nbvec.case"=nbvec.case, 
    "nbvec.contr"=nbvec.contr,
    "nbvec.diff"=nbvec.case - nbvec.contr,
    "fit.ltv"=fit.ltv
  )
  
}


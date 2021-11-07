chiuchao16 <- function(y){
  
  #Eqn. 6 in file:///Users/senthil/Zotero/storage/Y6SWTR8K/1634.html#methods
  # we found the estimator to be undefined in microbiomes whenever f[3]=0 or f[4]=0.
  # we found it to correlate with Chao1 extremely highly. We expect our conclusions about Chao1 to carry over here.
  
  Sobs <- sum( y>0 )
  n <- sum(y)
  tab <- table( y )
  f <- tab[c("1","2","3","4")]
  f[is.na(f)] <- 0
  f1h <- (2*(f[2]^2)/(3*f[3])) + 2*f[2]*( (f[2]/(3*f[3])) - (f[3]/(4*f[4])) )
  
  if( f[2] >0){
    res <- Sobs - f[1] +  f1h + ((n-1)/n)*((f[1]^2)/(2*f[2]))
  } else {
    res <- Sobs - f[1] +  f1h + (f1h*(f1h-1))/(2*(f[2]+1))
  }
  
  res
}


ccobj <- function(eobj){
  mat <- exprs(eobj)
  sn <- colnames(mat); names(sn) <- sn
  sapply( sn, function(x) chiuchao16( mat[,x] ) )
}

getlfc <- function(mat, sindx.case, sindx.cont){
  mat[!is.finite(mat)] <- NA
  log(rowMeans(mat[,sindx.case], na.rm=TRUE) ) - log( rowMeans(mat[,sindx.cont], na.rm=TRUE) ) 
}


getBreakaway <- function(obj, ... ){
  library(breakaway)
  gbr.sub <- function(obj){
    tryCatch({
      frequencytablelist <- build_frequency_count_tables(exprs(obj))
      lapply(frequencytablelist, function(x){
       tryCatch( {
         breakaway( x )
       }, error=function(e){
         print(e)
         list( "est"=NA,
               "seest"=NA
               )
       } ) 
      })  
    },error=function(e){
      as.list( rep(NA, ncol(obj)) )
    })
  }
  ss <- colSums( exprs(obj) )
  gbbr <- as.list( rep(NA, length(ss)) )
  indx <- which( ss==0 )
  if(length(indx)>0){
    gbbr[-indx] <- gbr.sub( obj[,-indx,drop=FALSE] ) 
  } else {
    gbbr <- gbr.sub( obj ) 
  }
  
  gb.mu <- do.call( cbind, lapply( gbbr, function(x){
    if( is.na(x)){
      0
    } else {
      x$est  
    }
  }  ) )
  gb.se <- do.call( cbind, lapply( gbbr, function(x) {
    if( is.na(x)){
      0
    } else {
      x$seest
    }
  }  ) )
  
  list( 
    "gbbr"=gbbr, 
    "gb.mu"=gb.mu, 
    "gb.se"=gb.se
  )
  
}

getBetta <- function(chats.mat, ses.mat, X, adj=TRUE, ... ){
  library(breakaway)
  res <- lapply( seq(nrow(chats.mat)), function(i){
    tryCatch({
      y <- chats.mat[i,]
      y.se <- ses.mat[i,]
      if( any(y.se==0, na.rm = TRUE) ){
        y.se[y.se==0] <- min(y.se[y.se>0])  
      }
      betta(chats = y,ses = y.se, X = X)  
      },error=function(e){
        print(e)
        list("NA")
    })
  } )
  names(res) <- rownames(chats.mat)
  res
  
}

getWLSAsympt <- function(chats.mat, ses.mat, X, adj=TRUE, ... ){
    library(breakaway)
    res <- lapply( seq(nrow(chats.mat)), function(i){
        tryCatch({
            y <- chats.mat[i,]
            y.se <- ses.mat[i,]
            se0 <- which( y.se==0 )
            if( length(se0)>0 ){
                y.se[se0] <- min(y.se[y.se>0])
            }
            sena <- which( !is.finite(y.se) )
            if(length(sena)>0){
                y.se[ sena ] <- max(y.se, na.rm=TRUE)
            }
            w <- 1/(y.se^2); w <- w/sum(w)
            summary( lm( y ~ -1 + X, weights=w ) )
        },error=function(e){
            print(e)
            list("NA")
        })
    } )
    names(res) <- rownames(chats.mat)
    res
    
}

getSwideStatsMat <- function(obj, ... ){
    sws <- getSwideStats( obj, ...  )
    chats.mat <- rbind( sws$chao1, sws$ace, sws$gb.mu )
    rownames( chats.mat ) <- c("chao1", "ace","breakaway")
    chats.se <- rbind( sws$chao1.se, sws$ace.se, sws$gb.se )
    rownames( chats.se ) <- c("chao1", "ace","breakaway")
    list( "chats.mat"=chats.mat, "chats.se"=chats.se )
}

getSwideInferences <- function(obj, X, ... ){
  sws <- getSwideStats( obj, ...  )
  chats.mat <- rbind( sws$chao1, sws$ace, sws$gb.mu )
  rownames( chats.mat ) <- c("chao1", "ace","breakaway")
  chats.se <- rbind( sws$chao1.se, sws$ace.se, sws$gb.se )
  rownames( chats.se ) <- c("chao1", "ace","breakaway")
  bet <- getBetta( chats.mat=chats.mat, ses.mat = chats.se, X=X, ... )
  tab <- list( "chao1"=bet$chao1$table, 
               "ace"=bet$ace$table, 
               "breakaway"=bet$breakaway$table 
               )
  list(
    "chats.mat"=chats.mat,
    "chats.se"=chats.se,
    "tab"=tab,
    "bet"=bet
  )
}


getGwideInferences <- function(obj, g, X, ncore=4, fns2export=c("getSwideStats", "getGwideStats", 
                                                                "getSwideInferences", "ccobj", "getBetta", "getBreakaway","chiuchao16"), ... ){
  cl <- makeCluster( ncore )
  registerDoParallel(cl)
  
  g <- as.character(g)
  ug <- unique( g )
  
  res <- foreach( gi=ug, .export = fns2export ) %dopar% {
    fi <- which( g==gi )
    require(metagenomeSeq)
    sws <- getSwideStats( obj[fi,,drop=FALSE], ...  )
    chats.mat <- rbind( sws$chao1, sws$ace, sws$gb.mu )
    rownames( chats.mat ) <- c("chao1", "ace","breakaway")
    chats.se <- rbind( sws$chao1.se, sws$ace.se, sws$gb.se )
    rownames( chats.se ) <- c("chao1", "ace","breakaway")
    bet <- getBetta( chats.mat=chats.mat, ses.mat = chats.se, X=X, ... )
    tab <- list( "chao1"=bet$chao1$table, 
                 "ace"=bet$ace$table, 
                 "breakaway"=bet$breakaway$table 
    )
    list(
      "chats.mat"=chats.mat,
      "chats.se"=chats.se,
      "tab"=tab,
      "bet"=bet
    ) 
  }
  
  
  stopCluster(cl)
  
  res
  
}



getSwideStats <- function(obj, ... ){
  obj <- obj[ which(rowSums( exprs(obj) ) > 0 ), ]
  getGwideStats( obj, g=rep( 1, nrow(obj) ), ...  )
}

getGwideStats <- function(obj, levl=NULL, g=as.character(fData(obj)[[levl]]), ... ){
  
  obj <- obj[which(rowSums( exprs(obj) )>0),]
  
  
  gabun <- do.call(rbind, lapply(as.character(unique(g)), function(x) colSums( exprs(obj[which( g ==x ), , drop=FALSE]) ) ))
  
  res <- lapply( unique(g), function(x){
    oo <- obj[ which( g==x ), ,drop=FALSE ]
    vegan:::estimateR( t(exprs(oo)) )
  })  
  chao1 <- do.call( rbind, lapply( res, function(x) x[2,] ) )
  chao1.se <- do.call( rbind, lapply( res, function(x) x[3,] ) )
  ace <-  do.call( rbind, lapply( res, function(x) x[4,] ) )
  ace.se <-  do.call( rbind, lapply( res, function(x) x[5,] ) )
  
  gbs <- lapply( unique(g), function(x){
    oo <- obj[ which( g==x ), ,drop=FALSE ]
    getBreakaway( oo ) 
  })  
  gb.mu <- do.call( rbind, lapply(gbs, function(x) x$gb.mu) ) 
  gb.se <- do.call( rbind, lapply(gbs, function(x) x$gb.se) ) 
  #getBetta( c(gb.mu), c(gb.se), X =  )
  
  sh <- do.call(rbind, lapply( unique(g), function(x){
    oo <- obj[ which( g==x ), ,drop=FALSE ]
    vegan:::diversity( t(exprs(oo)), index="shannon" )
  })  )
  
  res <- lapply( unique(g), function(x){
    oo <- obj[ which( g==x ), ,drop=FALSE ]
    ccobj( oo ) 
  })
  cc16 <- do.call(rbind, res)
  
  
  list( gabun=gabun, 
        chao1=chao1,
        chao1.se=chao1.se,
        ace = ace,
        ace.se = ace.se, 
        gbs=gbs,
        gb.mu=gb.mu,
        gb.se=gb.se,
        sh = sh,
        cc16 = cc16
  )
  
}


genFig <- function(gws, sindx.case, sindx.cont, plt.sh=FALSE, plt.cc16=FALSE){ #shannon, not "exactly" richness, but still used so let's show it in the supplementary.
  afc <- getlfc( gws$gabun, sindx.case = sindx.case, sindx.cont = sindx.cont  )
  chao1.fc <- getlfc( gws$chao1, sindx.case = sindx.case, sindx.cont = sindx.cont  )
  ace.fc <- getlfc( gws$ace, sindx.case = sindx.case, sindx.cont = sindx.cont  )
  sh.fc <- getlfc( gws$sh, sindx.case = sindx.case, sindx.cont = sindx.cont  )
  cc16.fc <- getlfc( gws$cc16, sindx.case = sindx.case, sindx.cont = sindx.cont  )
  gb.fc <- getlfc( gws$gb.mu, sindx.case = sindx.case, sindx.cont = sindx.cont  )
  
  op <- par(mfrow=c(1,ifelse(plt.sh, 4, 3)), pty='s')
  plot( chao1.fc ~ afc, xlab = "RelAb.LFC", ylab = "Chao1 LFC" ); abline( h=0 ); abline( v=0 )
  plot( ace.fc ~ afc, xlab = "RelAb.LFC", ylab = "ACE LFC" ); abline( h=0 ); abline( v=0 )
  plot( gb.fc ~ afc, xlab = "RelAb.LFC", ylab = "Breakaway" ); abline( h=0 ); abline( v=0 )
  if( plt.sh ){
    plot( sh.fc ~ afc, xlab = "RelAb.LFC", ylab = "ACE LFC" ); abline( h=0 ); abline( v=0 )
  }
  if( plt.cc16 ){
    tryCatch({
      plot( cc16.fc ~ afc, xlab = "RelAb.LFC", ylab = "CC16 LFC" ); abline( h=0 ); abline( v=0 )
    },error=function(e){
      print(e)
      print(cc16.fc)
    })
    
  }
  #par(op)
  
}

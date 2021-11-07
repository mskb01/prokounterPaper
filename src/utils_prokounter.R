#' Converst a metagenomeSeq object to an ExpressionSet object.
#' We need such functions because all all the functions in this project
#' utilize ExpressionSet objects.
#'
#'@param mrexpobj a MRexperiment object
convertMRtoES <- function(mrexpobj){
  require(Biobase)
  ExpressionSet( assayData = MRcounts(mrexpobj, norm=F),
                 phenoData = phenoData(mrexpobj), 
                 featureData = featureData(mrexpobj) )
}

#' With some matrix and group information, creates an ESet. 
#' @param mat matrix
#' @param grp vector of group information
createESet <- function(mat, grp){
  grp <- as.matrix(grp)
  colnames(grp) <- "group"
  if(ncol(mat) != nrow(grp)){
    stop("Not enough group information.")
  }
  nmes <- paste("sample", seq(ncol(mat)))
  colnames(mat) <- nmes
  row.names(grp) <- nmes
  ExpressionSet( assayData = mat,
                 phenoData = AnnotatedDataFrame(as.data.frame(grp))
  )
}

#' Creates an edgeR DGEList object
#' @param obj an ExpressionSet object
#' @param factr a factor label, whose corresponding sample entries are 
#'                are listed in \code{pData(obj)[[factr]]} 
createEdgeRObj <-function(obj, factr){
  require(edgeR)
  pd <- as.matrix(pData(obj)[[factr]])
  na.ind <- is.na(pd)
  if(any(na.ind)){
    rownames(pd)[na.ind] <- colnames(exprs(obj))[na.ind]
  }
  test.obj <- DGEList(counts = exprs(obj),
                      group = pd
  )
}



#' Plots boxplots in increasing median order
#' @param values a vector 
#' @param factr categories within which the values are to be plotted
boxplotByOrder <- function( values, factr, cols=NULL, ... ){
  values[is.na(values)] <- 0
  df <- data.frame( "values"=values )
  factr <- as.character(factr)
  df <- cbind(df, "factr"=factr)
  bymed <- with( df, reorder( factr, values, median  )  )
  if(is.null(cols)){
    cols <- getColsForFactors( df$factr )
  }
  colors <- sapply( levels(bymed), function(x) cols$map[which(cols$map[,2]==x),1]  )
  boxplot( values~bymed, data=df, las=2, 
           col=colors, ...
  )
}

#' Simply wraps RColorBrewer and brewer.pal
getCols <- function(set="Set2", m=10, n=1, ... ){
  require(RColorBrewer)
  colorRampPalette(brewer.pal(m,set))(n)
}

#' Given a pData column, returns colors for each group.
#' @param gr the pData column name
#' @param ... parameters that can be passed to getCols function, which in turn
#'   feeds to RColorBrewer function.
getColsForFactors <- function(gr,...){
  #
  # This function is useful for building colors for data points
  #  when each datapoint is associated with some group.
  # Any extra arguments passed are passed to the getCols function
  #  above.
  # For a given group information, builds the colors
  # palette.
  #
  stopifnot(is.factor(gr))
  
  if(FALSE){
    cols <- c(matrix(0, length(gr), 1))
    levs <- unique(sort(gr))
    cs <- getCols(n=length(levs),...)
    for(l in seq(1, length(levs))){
      cols[gr==levs[l]] <- cs[l]
    }
  }
  
  levs <- levels(gr)
  cols <- c(matrix(0, length(levs), 1))
  cs <- getCols(n=length(levs),...)
  for(l in seq(1, length(levs))){
    cols[ gr==levs[l] ] <- cs[l]
  }
  list("cols"=cols, "map"=cbind( cs, levs ) )
}

getShapesForFactors <- function(gr, ...){
  shapes <- c( 15, 16, 17, 18, 25, seq(14, 0, by=-1) )
  cols <- c(matrix(0, length(gr), 1))
  levs <- unique(sort(gr))
  nlevs <- length(levs)
  if(nlevs > length(shapes)){
    stop("Not enough shapes available.")
  }
  cs <- shapes[seq(1,length(levs))]
  for(l in seq(1, length(levs))){
    cols[gr==levs[l]] <- cs[l]
  }
  mapdf <- data.frame(cbind( cs, levs ))
  mapdf[,2] <- levs
  list("shapes"=cols, "map"=mapdf )
  
}





extrGLMSummary.series <- function( glmfits.list, type=NULL, row.indx=c(2,5), col.indx=c(1,2,3,4), ...  ) {
  require(abind)
  res <- lapply( glmfits.list, function(ft){
    if(is.null(type)){
      smry <- extrGLMSummary(glmfit = ft, row.indx=row.indx, col.indx=col.indx, ... )
    } else {
      smry <- extrGLMSummary.poisson(glmfit = ft, row.indx=row.indx, col.indx=col.indx, ... )
      }
    
    
  } )
  bound <- abind( res, along = 3 )
  apply(bound, 1, function(x) data.frame( t(x) ) ) 
}

extrGLMAIC <- function(glmfit){
  res <- NA
  if( !is.na(glmfit) ){
    if(glmfit$converged ){
      res <- AIC( glmfit )
    }
  }
  res  
}
extrGLMAIC.series <- function(glmfit.s){
  sapply( glmfit.s, extrGLMAIC)
}

extrGLMDeviance <- function(glmfit){
  res <- NA
  if( !is.na(glmfit) ){
    if(glmfit$converged ){
      res <- deviance( glmfit )
    }
  }
  res  
}
extrGLMDeviance.series <- function(glmfit.s){
  sapply( glmfit.s, extrGLMDeviance)
}



extr.disp <- function(glmfit){
  res <- NA
  if( !is.na(glmfit) ){
    if(glmfit$converged ){
      res <- summary( glmfit )$dispersion
    }
  }
  res  
}
extr.disp.series <- function(glmfit.s){
  sapply( glmfit.s, extr.disp)
}

pltImage <-function( omeg, nofcol=1000, axis.plt=FALSE, ... ){
  
  require(RColorBrewer)
  icols <- rev(colorRampPalette( brewer.pal( n = 10, name="RdYlBu" ) )(nofcol))
  
  x <- ncol(omeg)
  y <- nrow(omeg)
  image( seq(x),seq(y), t(omeg), col=icols, xaxt="n", yaxt="n", xlab="", ylab="", ... )
  if( axis.plt ){
    axis( 1, at=seq(x), colnames(omeg), las=2, pos=-0.25)
    axis( 2, at=seq(y), rownames(omeg), las=2, cex.axis=.5  )  
    #axis(2, at=seq(5), labels=colnames(omeg), 
    #     lwd=0, pos=-0.25)
  }
  
  
}

otuGrowthToGenusIntensity <- function( obj1, levv="Genus", legend.on=FALSE, plt.sdepth.dep=FALSE, plt.trend=FALSE, plt.ExceptDom=FALSE, ...   ){
  
  obj1 <- obj1[ which(rowMeans( exprs(obj1) )>0), ]
  
  if(plt.sdepth.dep){
    otu.per.sample <- colSums( exprs( obj1 )>0 )
    sdepth <- colSums( exprs(obj1) )
    plot( as.numeric(log2(otu.per.sample)) ~ log2(sdepth), 
          xlab = "Log2 Sample Depth", ylab = "Log2 Num. Taxa"  )
  }
  
  
  
  otu.cnt <- table( fData(obj1)[,levv] )
  gi <- which(  !( names(otu.cnt) %in% c( "", "g__" )) )
  otu.cnt <- otu.cnt[ gi ] 
  
  mr <- newMRexperiment( exprs(obj1), 
                         phenoData = phenoData(obj1), 
                         featureData = featureData(obj1)  )
  obj1.ge <- convertMRtoES(metagenomeSeq:::aggregateByTaxonomy( mr, lvl=levv ) ) #genus only
  
  fi <- intersect( rownames(obj1.ge), names(otu.cnt)  )  
  otu.cnt <- otu.cnt[fi]
  obj1.ge <- obj1.ge[fi,]
  
  ge.names <- names( otu.cnt)
  ge.sum <- rowSums( exprs(obj1.ge)[ge.names,,drop=FALSE] )
  ge.mu <- rowMeans( exprs(obj1.ge)[ge.names,,drop=FALSE] )
  
  
  #op <- par(mfrow=c(1,1), pty='s', mar=c(4,4,2,1))
  y <-  as.numeric(log2(otu.cnt))
  x <- as.numeric(log2(ge.mu))
  x <- as.numeric(log2(ge.sum))
  
  #my.fit <- lm( y ~ x )
  plot( y ~ x, ylab = "Log2 Num. Taxa", xlab = "Log2 Dataset-Wide Genus Recov. Abundance", ...  )
  #abline( my.fit$coeff[1], my.fit$coeff[2], lty=2, col = 'red' )
  if(legend.on)
      legend( "topright", legend= paste( 2^my.fit$coefficients[2] ) )
  
  df <- data.frame("x"=x, "y"=y)
  if(plt.trend){
    
    myf <- ssanova( y ~ x, data=df )
    df <- df[ order( df$x), ]
    #points( predict( myf, newdata=df ) ~ df$x , col = 'red', cex=.75 )
    lines( predict( myf, newdata=df ) ~ df$x , col = 'red', lwd=2 )
  }
  
  
  glvs <- as.character(fData(obj1)[[levv]])
  ug <- unique(glvs)
  if(plt.ExceptDom){
    z <- do.call(rbind, lapply( ug, function(g){
      fi <- which( glvs == g )
      sm <- exprs(obj1)[fi,,drop=FALSE]
      first2doms <- sort( rowSums( sm ) , decreasing=TRUE)[1:2] # for both strand clusters
      data.frame( "recov.ab"=sum(rowSums( sm )), 
                  "first.dom.cnt"=first2doms[1], 
                  "sec.dom.cnt"=first2doms[2],
                  "num.ftrs"=length(fi) )
    } ) )
    x <- z$recov.ab
    y <- z$recov.ab - z$first.dom.cnt - z$sec.dom.cnt
    
    plot( log2(y) ~ log2(x) )
    #summary(lm(log2(y) ~ log2(x)))
    z
  } else {
    df
  }
  
  
  
  #par(op)
  
}
pltMemsVsCountForWithinGenusClustersOnly <- function(mat, mems, categ.genus, bplot=TRUE, ... ){
  
  otugr <-do.call( rbind, lapply( mems, function(fi){
    ug <- unique( categ.genus[fi] )
    stopifnot( length(ug)==1 )
    genus <- categ.genus[ fi ][1]
    data.frame( "genus"=genus, "cluster.cnt"=sum( colSums( mat[fi,,drop=FALSE]) ) )
  } ) )
  
  genus.nummems <- table( otugr$genus )
  genus.tot <- sapply( names(genus.nummems), function(gt) 
    sum( otugr$cluster.cnt[which(otugr$genus==gt)] )
  )
  genus.nummems <- c(genus.nummems)
  
  y <- log2(genus.nummems)
  x <- log2( genus.tot )
  plot( y ~ x, xlab = "Log2 Genus Total", ylab = "Log2 Num. Taxa", ... )
  #print( summary( lm(y~x) ) )
  
  if(bplot){
    boxplot( y ~ cut( x, breaks = 15  ) )
  }
  
  list( "y"=y, "x"=x )
}

pltTaxAndRelFreqForGenus <- function(obj, levv, genus.name, rel.freq=TRUE, cc=NULL, legend.on=FALSE){
  
  fi <- grep( genus.name ,fData(obj)[[levv]] )
  
  tmp.mat <- exprs(obj[fi,])
  xl <- bquote( log[10]~.(genus.name)~"Recovered Abundance" )
  yl <- bquote( log[10]~"Num."~.(genus.name)~"Taxa")
  
  l10nf <- log10(colSums(tmp.mat>0))
  l10ra <- log10(colSums(tmp.mat))
  df <- data.frame(l10nf, l10ra)
  plot( l10nf ~ l10ra, xlab = xl, ylab=yl, 
        col = cc$cols, data=df )
  if(legend.on){
    legend( "topleft", cc$map[,2], fill = cc$map[,1] , bty="n" )
  }
  
  if(rel.freq){
    rf <- rowSums(tmp.mat)/sum(rowSums(tmp.mat))
    xl <- bquote(log[10]~"Net Rel. Freq. "~.(genus.name)~"Taxa" )
    hist(log10(rf), xlab = xl, main ="")
    legend( "topright", c( paste("Total num. taxa: ", nrow(tmp.mat) ), 
                           paste("Median num. taxa per sample: ", median(colSums(tmp.mat>0)))
                           ), bty="n" )
  }
  
  
  df
}



pltMat <- function(x=1:nrow(mat),mat, yl, xl, facs=NULL, legend.on=TRUE, legend.pos="topleft", ccset="Paired", legend.tit=NULL, ... ){
  if(is.null(facs)){
    facs <- factor(colnames(mat))  
  }
  cc <- getColsForFactors( facs, set=ccset  )
  matplot( x,mat, col = cc$map[,1], type='l', lty=1, ylab = yl, xlab = xl)
  if(legend.on){
    legend(legend.pos, legend=cc$map[,2], fill = cc$map[,1], bty="n", title=legend.tit )
  }
  
}

boxplotByOrder <- function( values, factr, cols=NULL, pnts=FALSE, ... ){
  values[is.na(values)] <- 0
  df <- data.frame( "values"=values )
  factr <- as.character(factr)
  df <- cbind(df, "factr"=factr)
  bymed <- with( df, reorder( factr, values, median  )  )
  if(is.null(cols)){
    cols <- getColsForFactors( df$factr )
  }
  colors <- sapply( levels(bymed), function(x) cols$map[which(cols$map[,2]==x),1]  )
  b <- boxplot( values~bymed, data=df, las=2, 
           col=colors, ...
  )
  if(pnts)
    points( values~bymed, data=df )
  b
}


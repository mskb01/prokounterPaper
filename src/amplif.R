
get_cond_s_for_n <- function(sub, spar, L, pm, n){
  if( (spar > L) || (sub > L) ){
    stop("Bad spar or sub argument: spar > L, or sub > L.")
  }
  tip <- min(spar,sub)
  fnm1 <- (1-exp(n*log(1-4*pm)))
  sum( sapply(seq( 0, tip  ), function(k){
    exp( dbinom( x = k,  size = spar, prob = .75*(1-pm)*fnm1, log = TRUE ) + 
           dbinom( x = sub-k,  size = (L-spar), prob = 3*pm, log = TRUE )
    )
  }) )
}

getCondMat_for_n <- function(L, pm, n){ #one step prob
  #matrices to store intermediate results
  hmat <-matrix( 0, L+1, L+1 )
  for( spar in seq(0,L) ){
    for( sub in seq(0, L)  ){
      hmat[ spar+1, sub+1  ] <- get_cond_s_for_n( sub, spar, L, pm, n )
    }
  }
  hmat
}

getCondMat.array <- function(L, pm, N, nc=4){ #one step prob
  #matrices to store intermediate results
  library(doParallel)
  library(abind)
  registerDoParallel(cores = nc)
  hmat <-array( 0, c(N, L+1, L+1) )
  arr.comb <- function(...) abind(..., along=3)
  hmat <- foreach( n=seq(1,N), .combine = 'arr.comb', .multicombine = TRUE ) %dopar% {
    getCondMat_for_n( L, pm, n )
  }
  hmat
  
}

getMarg <- function(condMatArr){ #yields the marginals for all n along the third axis (i.e., conditional of Sn | S0=0).
  Pn.giv.0 <- array(0, dim = dim( condMatArr  )) 
  Pn.giv.0[,,1]<- condMatArr[,,1]
  for( n in seq(2,dim(condMatArr)[3]) ){
    Pn.giv.0[,,n] <- Pn.giv.0[,,n-1]%*%condMatArr[,,n]
  } #transiting to states in n^th step when the initial state is 0.  
  Pn.giv.0
}

#matrices to store intermediate results
#hmat <- getCondMat( L=LENGTH, pm=pM )
#margMat <- getMargMat( N=40, LENGTH, condMat = hmat, s0=0 )

#LENGTH <- 250
#PM <- .00012 #3pm is substiution generation rate; pm is sub to truth rate; reference for rate: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5218489/

#pM.kappa <- 1.6 * (10^(-5)) # reference: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5218489/ , Table 3
#pM.Q5 <- 5.3*(10^(-7))
#CYCLES <- 50
#N <- 14 # = CYCLES-1 from simulation for comparison purposes.  number of cycles, bridge PCR included; N=50 computationally infeasible to sample from binomials with size 2^50.

#hmat <- getCondMat.array( L=LENGTH, pm=pM, N=50 )
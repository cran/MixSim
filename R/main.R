
overlap <- function(Pi, Mu, S, eps = 1e-06, acc = 1e-06, lim = 1e06){

	K <- dim(Mu)[1]
	p <- dim(Mu)[2]

#	if (!is.loaded("runExactOverlap")) {
#		dyn.load("libMixSim.so")
#	}

        Mu1 <- as.vector(t(Mu))
        S1 <- as.vector(S)
        OmegaMap1 <- rep(0, K*K)

        pars <- c(eps, acc)
        rcMax <- c(0, 0)
        

	Q <- .C("runExactOverlap", p1 = as.integer(p), K1 = as.integer(K), Pi = as.double(Pi), Mu1 = as.double(Mu1), S1 = as.double(S1), pars = as.double(pars), lim1 = as.integer(lim), OmegaMap1 = as.double(OmegaMap1), BarOmega = as.double(1), MaxOmega = as.double(1), rcMax = as.integer(rcMax), PACKAGE = MixSim)

        return(list(OmegaMap = matrix(Q$OmegaMap1, byrow = TRUE, ncol = K), BarOmega = Q$BarOmega, MaxOmega = Q$MaxOmega, rcMax = Q$rcMax + 1))

}



MixSim <- function(BarOmega = NULL, MaxOmega = NULL, K, p, sph = 0, ecc = 0.90, PiLow = 1.0, Ubound = 1.0, resN = 100, eps = 1e-06, acc = 1e-06, lim = 1e06){

        Pi <- rep(0, K)
        Mu1 <- rep(0, K*p)
        S1 <- rep(0, K*p*p)
        OmegaMap1 <- rep(0, K*K)
        rcMax <- c(0, 0)
        pars <- c(eps, acc)

 
        if ((is.null(MaxOmega)) & (!is.null(BarOmega))){
               method <- 0
               Omega <- BarOmega
        }
        if ((is.null(BarOmega)) & (!is.null(MaxOmega))){
               method <- 1
               Omega <- MaxOmega
             }
        if ((!is.null(BarOmega)) & (!is.null(MaxOmega)))  method <- 2   
		if ((is.null(BarOmega)) & (is.null(MaxOmega)))  method <- -1   

        if ((method == 0) | (method == 1)){
        
#              if (!is.loaded("runOmegaClust")) {
#			dyn.load("libMixSim.so")
#              }
     
              Q <- .C("runOmegaClust", Omega1 = as.double(Omega), method1 = as.integer(method), p1 = as.integer(p), K1 = as.integer(K), PiLow1 = as.double(PiLow), Ubound1 = as.double(Ubound), emax1 = as.double(ecc), pars = as.double(pars), lim1 = as.integer(lim), resN1 = as.integer(resN), sph1 = as.integer(sph), Pi = as.double(Pi), Mu1 = as.double(Mu1), S1 = as.double(S1), OmegaMap1 = as.double(OmegaMap1), BarOmega = as.double(1), MaxOmega = as.double(1), rcMax = as.integer(rcMax), fail = as.integer(1), PACKAGE = MixSim)


        }

        if (method == 2){
              
#              if (!is.loaded("runOmegaBarOmegaMax")) {
#			  		dyn.load("libMixSim.so")
#              }

              Q <- .C("runOmegaBarOmegaMax", p1 = as.integer(p), K1 = as.integer(K), PiLow1 = as.double(PiLow), Ubound1 = as.double(Ubound), emax1 = as.double(ecc), pars = as.double(pars), lim1 = as.integer(lim), resN1 = as.integer(resN), sph1 = as.integer(sph), Pi = as.double(Pi), Mu1 = as.double(Mu1), S1 = as.double(S1), OmegaMap1 = as.double(OmegaMap1), BarOmega = as.double(BarOmega), MaxOmega = as.double(MaxOmega), rcMax = as.integer(rcMax), fail = as.integer(1), PACKAGE = MixSim)

        }

        if (method != -1){

		if (Q$fail == 0){
		       	return(list(Pi = Q$Pi, Mu = matrix(Q$Mu1, byrow = TRUE, ncol = p), S = array(Q$S1, c(p, p, K)), OmegaMap = matrix(Q$OmegaMap1, byrow = TRUE, ncol = K), BarOmega = Q$BarOmega, MaxOmega = Q$MaxOmega, rcMax = Q$rcMax + 1, fail = Q$fail))
		}

	} else {
		cat("Error: at least one overlap characteristic should be specified...\n")		
	}
		
}



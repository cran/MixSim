
simdataset <- function(n, Pi, Mu, S){

	K <- dim(Mu)[1]
	p <- dim(Mu)[2]

	Sample <- NULL
	Nk <- drop(rmultinom(1, n, Pi))
	id <- NULL

	for (k in 1:K){
		id <- c(id, rep(k, Nk[k]))
		Sample <- rbind(Sample, mvrnorm(n = Nk[k], mu =
		Mu[k,], Sigma = S[,,k]))
	}

	return(list(x = Sample, id = id))
}


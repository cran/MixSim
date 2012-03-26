### Example of Section 4

library("mclust")
library("cluster")

p <- 4
K <- 6
n <- 200
n.files <- 100
omega <- c(0.4, 0.3, 0.25, 0.2, 0.15, 0.1, 0.05, 0.01, 0.005, 0.001)

set.seed(1234)

AR.Mclust <- rep(NA, n.files)
AR.Kmeans <- rep(NA, n.files)
AR.PAM <- rep(NA, n.files)
AR.Ward <- rep(NA, n.files)
	
P.Mclust <- rep(NA, n.files)
P.Kmeans <- rep(NA, n.files)
P.PAM <- rep(NA, n.files)
P.Ward <- rep(NA, n.files)
	
VI.Mclust <- rep(NA, n.files)
VI.Kmeans <- rep(NA, n.files)
VI.PAM <- rep(NA, n.files)
VI.Ward <- rep(NA, n.files)

for (bo in omega){
	
	for (i in 1:n.files){

		Q <- MixSim(BarOmega = bo, K = K, p = p)
		A <- simdataset(n = n, Pi = Q$Pi, Mu = Q$Mu, S = Q$S)

		# Use try(...) to catch Mclust failure instances
		id.Mclust <- try(Mclust(A$X, G = K, model = "VVV")$class)
		# Provide error ids if a failure occurs
		if (inherits(id.Mclust, what = "try-error") || (is.null(id.Mclust)))
			id.Mclust <- rep(1, n)

		id.Kmeans <- kmeans(A$X, K, iter.max = 1000, nstart = 10)$cluster
		id.PAM <- pam(A$X, K)$clustering
		d <- dist(A$X, method = "euclidean")
		id.Ward <- cutree(hclust(d = d, method = "ward"), k = K)

		AR.Mclust[i] <- RandIndex(A$id, id.Mclust)$AR
		AR.Kmeans[i] <- RandIndex(A$id, id.Kmeans)$AR
		AR.PAM[i] <- RandIndex(A$id, id.PAM)$AR
		AR.Ward[i] <- RandIndex(A$id, id.Ward)$AR

		P.Mclust[i] <- ClassProp(A$id, id.Mclust)
		P.Kmeans[i] <- ClassProp(A$id, id.Kmeans)
		P.PAM[i] <- ClassProp(A$id, id.PAM)
		P.Ward[i] <- ClassProp(A$id, id.Ward)

		VI.Mclust[i] <- VarInf(A$id, id.Mclust)
		VI.Kmeans[i] <- VarInf(A$id, id.Kmeans)
		VI.PAM[i] <- VarInf(A$id, id.PAM)
		VI.Ward[i] <- VarInf(A$id, id.Ward)

	}

}


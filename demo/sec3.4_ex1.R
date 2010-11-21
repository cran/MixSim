### Example of Section 3.4 plot (a)
data(iris)
K <- 3
p <- dim(iris)[2] - 1
n <- dim(iris)[1]
id <- as.numeric(iris[, 5])

Pi <- NULL
Mu <- NULL
S <- array(rep(0, p * p * K), c(p, p, K))
for (k in 1:K){
	Pi <- c(Pi, sum(id == k) / n)
	Mu <- rbind(Mu, apply(iris[id == k, -5], 2, mean))
	S[, , k] <- var(iris[id == k, -5])
}
pdplot(Pi = Pi, Mu = Mu, S = S)

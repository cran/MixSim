### Example for Section 3.2
p <- dim(iris)[2] - 1
n <- dim(iris)[1]
K <- 3
id <- as.numeric(iris[, 5])

Pi <- NULL
Mu <- NULL
S <- array(rep(0, p * p * K), c(p, p, K))

# estimate mixture parameters
for (k in 1:K){
	Pi <- c(Pi, sum(id == k) / n)
	Mu <- rbind(Mu, apply(iris[id == k, -5], 2, mean))
	S[, , k] <- var(iris[id == k, -5])
}

overlap(Pi = Pi, Mu = Mu, S = S)

### Example of Section 3.4 plot (a)
data(iris)
p <- dim(iris)[2] - 1
K <- 3
id <- as.numeric(iris[, 5])

Pi <- sapply(1:K, function(k){ sum(id == k) }) / dim(iris)[1]
Mu <- t(sapply(1:K, function(k){ colMeans(iris[id == k, -5]) }))
S <- sapply(1:K, function(k){ var(iris[id == k, -5]) })
dim(S) <- c(p, p, K)

pdplot(Pi = Pi, Mu = Mu, S = S)

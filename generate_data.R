### ENSURE THAT THE SAME DATA SETS ARE CREATED EACH TIME THIS IS RUN
set.seed(12345)

### MAKE THE DIRECTORY TO HOLD THE DATA, SPLITS, AND RNG SEEDS
dir.create("problems", showWarnings=FALSE)
dir.create("problems/data", showWarnings=FALSE)
dir.create("problems/splits", showWarnings=FALSE)
dir.create("problems/seeds", showWarnings=FALSE)


### MAKE THE "MODIFIED BISHOP" DATA
mbishop <- function(x, a=0) 2.3 * (x - a) + sin(2 * pi * (x - a)^2)

n.samples <- 100
n.per.sample <- 25
n.test <- 1000
n <- n.samples * n.per.sample + n.test
p <- 1
X <- matrix(runif(n * p), ncol=p)

for (a in c(1, 5, 10)) {
    h  <- mbishop(X, a)
    t <- h + rnorm(n, sd=0.3)

    write.table(cbind(X, h + rnorm(n, sd=0.3)), sprintf("problems/data/mbishop-%02d", a), row.names=FALSE, col.names=FALSE)

    writeBin(as.raw(as.integer(sample(0:255, 256 * n.samples, replace=TRUE))), sprintf("problems/seeds/mbishop-%02d", a))

    idx <- sapply(seq(n.samples), function(i, nsamples, ntrain, ntest) {
        idx <- rep(c("#", "1"), times=c(nsamples*ntrain, ntest))
        idx[seq(ntrain) + (i - 1) * ntrain] <- "0"
        paste(idx, collapse="")
    }, nsamples=n.samples, ntrain=n.per.sample, ntest=n.test)

    write(idx, sprintf("problems/splits/mbishop-%02d", a), sep="\n")
}


### MAKE THE FRIEDMAN1 DATA
friedman1 <- function(X) 10 * sin(pi * X[, 1] * X[, 2]) + 20 * (X[, 3] - 0.5)^2 + 10 * X[, 4] + 5 * X[, 5]

n.samples <- 100
n.per.sample <- 250
n.test <- 1000
n <- n.samples * n.per.sample + n.test
p <- 10
X <- matrix(runif(n * p), ncol=p)
h <- friedman1(X)
t <- h + rnorm(n)

write.table(cbind(X, t), "problems/data/friedman1", row.names=FALSE, col.names=FALSE)

writeBin(as.raw(as.integer(sample(0:255, 256 * n.samples, replace=TRUE))), "problems/seeds/friedman1")

idx <- sapply(seq(n.samples), function(i, nsamples, ntrain, ntest) {
    idx <- rep(c("#", "1"), times=c(nsamples*ntrain, ntest))
    idx[seq(ntrain) + (i - 1) * ntrain] <- "0"
    paste(idx, collapse="")
}, nsamples=n.samples, ntrain=n.per.sample, ntest=n.test)

write(idx, "problems/splits/friedman1", sep="\n")

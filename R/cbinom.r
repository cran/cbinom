dcbinom <- function(x, size, prob, log = FALSE){
  smalli <- which(x < 0)
  nsm <- length(smalli)
  for (i in smalli) warning(paste0("x = ", x[i], " < 0"), call. = T)
  bigi <- which(x > size + 1)
  for (i in bigi) warning(paste0("x = ", x[i], " > size + 1"), call. =T)
  fmat <- dcblp(x, size, prob)

  # columns of fmat are:
  #  log(pcbinom(x + h)), log(pcbinom(x - h)), 2*h if h < x < size + 1 - h
  #  log(pcbinom(x + h)), log(pcbinom(x)), h if x <= h
  #  log(pcbinom(x), log(pcbinom(x - h), h if x > size + 1 - h
  #   where h is small
  suppressWarnings({
    if (!log) {
      ans <- (exp(fmat[, 1]) - exp(fmat[, 2]))/fmat[,3]
      ans[ans < 0] <- 0
    } else {
      ans <- log((fmat[, 1] - fmat[, 2])/fmat[, 3]) +
        pcbinom(x, size, prob, log.p = T)
      ans[is.na(ans)] <- -Inf
    }
  })
  return(ans)
}
pcbinom <- function(q, size, prob, lower.tail = TRUE, log.p = FALSE){
  q[q < 0] <- 0
  q[q > size + 1] <- size + 1
  pbeta(prob, q, size - q + 1, lower.tail = !lower.tail, log.p = log.p)
}
qcbinom <- function(p, size, prob, lower.tail = TRUE, log.p = FALSE){
  qcbinomC(p, size, prob, lower.tail, log.p)
}
rcbinom <- function(n, size, prob){
  u <- runif(n, 0, 1)
  qcbinomC(log(u), size, prob, lowertail = T, logp = T, rcb = T)
}



calcLik <- function (x) {
  
  fe <-  apply(visit_tab, 1, negBinLik, x, rep(0, nrow(visit_tab)), c(0, 0, 0), c(0, 0, 0))
  post <- log(prod(fe))
              
  return (post)
              
}

liks <- lapply(seq(-10, -1, by = 0.01), calcLik)

plot(seq(-10, -1, by = 0.01), liks)

m <- seq(-10, -1, by = 0.01)[which.max(liks)]

abline(v = m], col = "red")

mean(rnbinom(1000, 1, exp(m)))



calcLik <- function (x) {
  
  fe <- apply(cd, 1, calcFe, sig.y, sig.s, sig.cc, sig.wc, sha2, x, sig.ss)
  u <- integrate(f.gamma.function, 0, 1000, x, sha2)$value/1000
  post <- log(prod(fe/u)) 
  
  return (post)
  
}

liks <- lapply(200:400, calcLik)

plot(200:400, liks)
abline(v = 318, col = "red")

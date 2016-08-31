
calcLik <- function (x) {
  
  fe <- apply(cd, 1, calcFe, sig.y, sig.s, sig.cc, sig.wc, sha2, x, sig.ss)
  u <- integrate(f.gamma.function, 0, 1000, x, sha2)$value/1000
  post <- log(prod((fe/u)/2)) 
  
  return (post)
  
}

liks <- lapply(100:400, calcLik)

plot(100:400, liks)
abline(v = 318, col = "red")

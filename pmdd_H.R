mirrorstat = function(y1, Z1, y2, Z2) {
  n1 = nrow(Z1)
  p = ncol(Z1)
  n2 = nrow(Z2)
  n = n1 + n2
  
  a1= diag( c(1/sqrt(1+ Z1^2)) )
  a2= diag( c(1/sqrt(1+ Z2^2)) )
  a=a1%*%(1+Z1%*%t(Z2))%*%a2
  
  a[a > 1] = 1
  a[a < -1] = -1
  A = asin(a)
  
  H1 = t(c(y1 - mean(y1)) * A) * c(y2 - mean(y2))
  H2 = t(c(y1 - mean(y2)) * A) * c(y2 - mean(y2))
  return(list(H1,H2))
}

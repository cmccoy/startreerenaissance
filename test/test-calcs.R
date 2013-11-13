m <- matrix(1, nrow=4, ncol=4)
m
diag(m) <- -3
m
s <- eigen(m)

eig <- eigen(m)
eig
inv(eig$vector)
solve(eig$vector)
solve(eig$vector) * diag(exp(eig$values, 0.02)) * eig$vector
solve(eig$vector) * diag(exp(eig$values * 0.02)) * eig$vector
eig$vector * diag(exp(eig$values * 0.02)) * eig$vector
eig$vector * diag(exp(eig$values * 0.02))
eig$vector * diag(exp(eig$values * 0.02)) * solve(eig$vector)
eig$vector %*% diag(exp(eig$values * 0.02)) %*% solve(eig$vector)
solve(eig$vector) %*% diag(exp(eig$values * 0.02)) %*% eig$vector
eig$vector %*% diag(exp(eig$values * 0.02)) %*% solve(eig$vector)
eig <- eigen(m)
m
m <- m / 3
m
eig <- eigen(m)
eig$vector %*% diag(exp(eig$values * 0.02)) %*% solve(eig$vector)
log(eig) * diag(rep(1, 4))
rep(1, 4)
diag(rep(1, 4))
p <- eig$vector %*% diag(exp(eig$values * 0.02)) %*% solve(eig$vector)
log(o) * diag(rep(1, 4))
log(p) * diag(rep(1, 4))
sum(log(p) * diag(rep(1, 4)))
m <- matrix(c(94, 3, 2, 1, 2, 95, 2, 1, 2, 4, 89,5, 1,3,2, 94), nrow=4, byrow=TRUE)
m
sum(log(p) * m)
writehistory('test-calcs.R'
)
savehistory('test-calcs.R'
)

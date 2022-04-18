debug <- 1
set.seed(1)



## load dataset
bcb:::load_example(eg = "dnet", method = "bcb-bma", network = "asia")
jpt <- bcb:::bn.fit2jpt(bn.fit = bn.fit)
nodes <- c("V8", "V6", "V2", "V4")
levels <- c(2, 2)
p <- bcb:::jpt2p(jpt, nodes, levels)



## generate data
set.seed(1)
n <- 1e3
nlarge <- 1e4
X <- rmultinom(n = nlarge, size = n, prob = p)
N <- t(X)
dim(N) <- c(ncol(X), nrow(p), 3)
P <- N / n
Q <- M <- W <- matrix(0, nrow = ncol(X), ncol = nrow(p))
M <- N[,,1] * (N[,,1] + N[,,2] + N[,,3]) / n^2
W <- (N[,,1] + N[,,2]) / n
Q <- M / W
Pr <- rowSums(Q)



## misc quantities
i <- 1
j <- ip1 <- (i + 1) %% 3
prime <- 1  # prime+1
phat <- P[1,,]



##################################################
## Test Pr: Done!
##################################################

## Var(Pr)
var(Pr)
# summary(temp1 <- apply(P, 1, bcb:::boot_Var_Pr, n = n))
summary(temp2 <- apply(P, 1, bcb:::Var_Pr, n = n, M_plus = 0, W_plus = 0))
# p_true <- apply(P, c(2, 3), mean)  # large sample estimate
# dim(p_true) <- c(1, dim(p_true))
# summary(temp3 <- apply(p_true[rep(1, nlarge),,], 1, bcb:::boot_Var_Pr, n = n))



##################################################
## Test Denominator W: Done!
##################################################

## E(W_i)
mean(W[,i])
summary(temp <- (N[,i,1] + N[,i,2]) / n)
summary(temp <- apply(P, 1, bcb:::E_W, n = n, i = i))

## Var(W_i)
var(W[,i])
summary(temp <- (N[,i,1] + N[,i,2]) / n * (1 - (N[,i,1] + N[,i,2]) / n)) / n
summary(temp <- apply(P, 1, bcb:::Var_W, n = n, i = i))

## Cov(W_i, W_j)
cov(W[,i], W[,j])
summary(temp <- -((N[,i,1] + N[,i,2]) * (N[,j,1] + N[,j,2]))) / n^3
summary(temp <- apply(P, 1, bcb:::Cov_W, n = n, i = i, j = j))



##################################################
## Test Numerator M: Done!
##################################################

## E(M_i)
mean(M[,i])
summary(temp <- apply(P, 1, bcb:::E_M, n = n, i = i))

## Var(N_i^2)
var(N[,i,1]^2)
summary(temp <- apply(P, 1, bcb:::Var_aa, n = n, i = i))

## Var(N_i pN_i)
var(N[,i,1] * N[,i,2])
summary(temp <- apply(P, 1, bcb:::Var_ab, n = n, i = i, prime = 1))
var(N[,i,1] * N[,i,3])
summary(temp <- apply(P, 1, bcb:::Var_ab, n = n, i = i, prime = 2))



## Cov(N_i^2, N_i pN_i)
cov(N[,i,1]^2, N[,i,1] * N[,i,2])
summary(temp <- apply(P, 1, bcb:::Cov_aa_ab, n = n, i = i, prime = 1))
cov(N[,i,1]^2, N[,i,1] * N[,i,3])
summary(temp <- apply(P, 1, bcb:::Cov_aa_ab, n = n, i = i, prime = 2))

## Cov(N_i pN_i, N_i ppN_i)
cov(N[,i,1] * N[,i,2], N[,i,1] * N[,i,3])
summary(temp <- apply(P, 1, bcb:::Cov_ab_ac, n = n, i = i))
cov(N[,1,1] * N[,1,2], N[,1,1] * N[,1,3])
summary(temp <- apply(P, 1, bcb:::Cov_ab_ac, n = n, i = 1))
cov(N[,2,1] * N[,2,2], N[,2,1] * N[,2,3])
summary(temp <- apply(P, 1, bcb:::Cov_ab_ac, n = n, i = 2))
cov(N[,3,1] * N[,3,2], N[,3,1] * N[,3,3])
summary(temp <- apply(P, 1, bcb:::Cov_ab_ac, n = n, i = 3))
cov(N[,4,1] * N[,4,2], N[,4,1] * N[,4,3])
summary(temp <- apply(P, 1, bcb:::Cov_ab_ac, n = n, i = 4))

## Var(M_i)
var(M[,i])
summary(temp <- apply(P, 1, bcb:::Var_M, n = n, i = i))



## Cov(N_i^2, N_j^2)
cov(N[,i,1]^2, N[,j,1]^2)
summary(temp <- apply(P, 1, bcb:::Cov_ii_jj, n = n, i = i, j = j))

## Cov(N_i^2, N_j pN_j)
cov(N[,i,1]^2, N[,j,1] * N[,j,2])
summary(temp <- apply(P, 1, bcb:::Cov_ii_jajb, n = n, i = i, j = j, prime = 1))
cov(N[,i,1]^2, N[,j,1] * N[,j,3])
summary(temp <- apply(P, 1, bcb:::Cov_ii_jajb, n = n, i = i, j = j, prime = 2))
cov(N[,j,1]^2, N[,i,1] * N[,i,2])
summary(temp <- apply(P, 1, bcb:::Cov_ii_jajb, n = n, i = j, j = i, prime = 1))
cov(N[,j,1]^2, N[,i,1] * N[,i,3])
summary(temp <- apply(P, 1, bcb:::Cov_ii_jajb, n = n, i = j, j = i, prime = 2))

## Cov(N_i pN_i, N_j pN_j)
cov(N[,i,1] * N[,i,2], N[,j,1] * N[,j,2])
summary(temp <- apply(P, 1, bcb:::Cov_iaib_jajc, n = n, i = i, j = j, prime_i = 1, prime_j = 1))
cov(N[,i,1] * N[,i,2], N[,j,1] * N[,j,3])
summary(temp <- apply(P, 1, bcb:::Cov_iaib_jajc, n = n, i = i, j = j, prime_i = 1, prime_j = 2))
cov(N[,i,1] * N[,i,3], N[,j,1] * N[,j,2])
summary(temp <- apply(P, 1, bcb:::Cov_iaib_jajc, n = n, i = i, j = j, prime_i = 2, prime_j = 1))
cov(N[,i,1] * N[,i,3], N[,j,1] * N[,j,3])
summary(temp <- apply(P, 1, bcb:::Cov_iaib_jajc, n = n, i = i, j = j, prime_i = 2, prime_j = 2))

## Cov(M_i, M_j)
cov(M[,i], M[,j])
summary(temp <- apply(P, 1, bcb:::Cov_M, n = n, i = i, j = j))



##################################################
## Test M and W:
##################################################

## Cov(N_i^2, N_i)
cov(N[,i,1]^2, N[,i,1])
summary(temp <- apply(P, 1, bcb:::Cov_ii_i, n = n, i = i))

## Cov(N_i^2, pN_i)
cov(N[,i,1]^2, N[,i,2])
summary(temp <- apply(P, 1, bcb:::Cov_ii_ja, n = n, i = i, j = i, prime = 1))

## Cov(N_i pN_i, N_i)
cov(N[,i,1] * N[,i,2], N[,i,1])
summary(temp <- apply(P, 1, bcb:::Cov_ab_a, n = n, i = i, prime_a = 0, prime_b = 1))
cov(N[,i,1] * N[,i,2], N[,i,2])
summary(temp <- apply(P, 1, bcb:::Cov_ab_a, n = n, i = i, prime_a = 1, prime_b = 0))
cov(N[,i,1] * N[,i,3], N[,i,1])
summary(temp <- apply(P, 1, bcb:::Cov_ab_a, n = n, i = i, prime_a = 0, prime_b = 2))

## Cov(N_i ppN_i, pN_i)
cov(N[,i,1] * N[,i,3], N[,j,2])
summary(temp <- apply(P, 1, bcb:::Cov_iaib_jc, n = n, i = i, j = j, prime_i = 2, prime_j = 1))

## Cov(M_i, W_i)
cov(M[,i], W[,i])
summary(temp <- apply(P, 1, bcb:::Cov_M_W, n = n, i = i))
cov(M[,1], W[,1])
summary(temp <- apply(P, 1, bcb:::Cov_M_W, n = n, i = 1))
cov(M[,2], W[,2])
summary(temp <- apply(P, 1, bcb:::Cov_M_W, n = n, i = 2))
cov(M[,3], W[,3])
summary(temp <- apply(P, 1, bcb:::Cov_M_W, n = n, i = 3))

## Cov(N_i^2, (p)N_j)
cov(N[,i,1]^2, N[,j,1])
summary(temp <- apply(P, 1, bcb:::Cov_ii_ja, n = n, i = i, j = j, prime = 0))
cov(N[,i,1]^2, N[,j,2])
summary(temp <- apply(P, 1, bcb:::Cov_ii_ja, n = n, i = i, j = j, prime = 1))

## Cov(N_i pN_i, (p)N_j)
cov(N[,i,1] * N[,i,2], N[,j,1])
summary(temp <- apply(P, 1, bcb:::Cov_iaib_jc, n = n, i = i, j = j, prime_i = 1, prime_j = 0))
cov(N[,i,1] * N[,i,2], N[,j,2])
summary(temp <- apply(P, 1, bcb:::Cov_iaib_jc, n = n, i = i, j = j, prime_i = 1, prime_j = 1))
cov(N[,i,1] * N[,i,3], N[,j,1])
summary(temp <- apply(P, 1, bcb:::Cov_iaib_jc, n = n, i = i, j = j, prime_i = 2, prime_j = 0))
cov(N[,i,1] * N[,i,3], N[,j,2])
summary(temp <- apply(P, 1, bcb:::Cov_iaib_jc, n = n, i = i, j = j, prime_i = 2, prime_j = 1))

## Cov(M_i, W_j)
cov(M[,i], W[,j])
summary(temp <- apply(P, 1, bcb:::Cov_Mi_Wj, n = n, i = i, j = j))
cov(M[,j], W[,i])
summary(temp <- apply(P, 1, bcb:::Cov_Mi_Wj, n = n, i = j, j = i))
cov(M[,1], W[,2])
summary(temp <- apply(P, 1, bcb:::Cov_Mi_Wj, n = n, i = 1, j = 2))
cov(M[,1], W[,3])
summary(temp <- apply(P, 1, bcb:::Cov_Mi_Wj, n = n, i = 1, j = 3))
cov(M[,1], W[,4])
summary(temp <- apply(P, 1, bcb:::Cov_Mi_Wj, n = n, i = 1, j = 4))
cov(M[,2], W[,3])
summary(temp <- apply(P, 1, bcb:::Cov_Mi_Wj, n = n, i = 2, j = 3))
cov(M[,2], W[,4])
summary(temp <- apply(P, 1, bcb:::Cov_Mi_Wj, n = n, i = 2, j = 4))
cov(M[,3], W[,4])
summary(temp <- apply(P, 1, bcb:::Cov_Mi_Wj, n = n, i = 3, j = 4))
## and more...



##################################################
## Test Q: Done!
##################################################

## E(Q_i)
mean(Q[,i])
summary(temp <- apply(P, 1, bcb:::E_Q, n = n, i = i))
mean(Q[,1])
summary(temp <- apply(P, 1, bcb:::E_Q, n = n, i = 1))
mean(Q[,2])
summary(temp <- apply(P, 1, bcb:::E_Q, n = n, i = 2))
mean(Q[,3])
summary(temp <- apply(P, 1, bcb:::E_Q, n = n, i = 3))
mean(Q[,4])
summary(temp <- apply(P, 1, bcb:::E_Q, n = n, i = 4))

## Var(Q_i)
var(Q[,i])
summary(temp <- apply(P, 1, bcb:::Var_Q, n = n, i = i))
var(Q[,1])
summary(temp <- apply(P, 1, bcb:::Var_Q, n = n, i = 1))
var(Q[,2])
summary(temp <- apply(P, 1, bcb:::Var_Q, n = n, i = 2))
var(Q[,3])
summary(temp <- apply(P, 1, bcb:::Var_Q, n = n, i = 3))
var(Q[,4])
summary(temp <- apply(P, 1, bcb:::Var_Q, n = n, i = 4))

## Cov(Q_i, Q_j)
cov(Q[,i], Q[,j])
summary(temp <- apply(P, 1, bcb:::Cov_Q, n = n, i = i, j = j))
## and more...




##################################################
## Test Multinomial Derivations: Done!
##################################################

## E(N_u^3)
mean(N[,i,prime]^3)
summary(temp <- sapply(N[,i,prime], function(x) E_k(n = n, p = x / n, k = 3)))
# all.equal(temp, sapply(N[,i,1], function(x) mnmoment(n = n, p = x / n, moment = 3)))

## E(N_u^4)
mean(N[,i,prime]^4)
summary(temp <- sapply(N[,i,prime], function(x) E_k(n = n, p = x / n, k = 4)))
# all.equal(temp, sapply(N[,i,prime], function(x) mnmoment(n = n, p = x / n, moment = 4)))

## E(Nu^2 Nv^2)
mean(N[,i,1]^2 * N[,i,2]^2)
summary(temp <- sapply(seq_len(nrow(N)),
                       function(x) E_2_2(n = n, p1 = N[x,i,1] / n, p2 = N[x,i,2] / n)))

## E(N_u^3 N_v)
mean(N[,i,1]^3 * N[,i,2])
summary(temp <- sapply(seq_len(nrow(N)),
                       function(x) E_3_1(n = n, p1 = N[x,i,1] / n, p2 = N[x,i,2] / n)))

## E(N_u^2 N_v N_w)
mean(N[,i,1]^2 * N[,i,2] * N[,i,3])
summary(temp <- sapply(seq_len(nrow(N)), function(x){

  E_2_1_1(n = n, p1 = N[x,i,1] / n,
          p2 = N[x,i,2] / n, p3 = N[x,i,3] / n)
}))

## E(N_u N_v N_w N_x)
mean(N[,i,1] * N[,i,2] * N[,i,3] * N[,ip1,1])
summary(temp <- sapply(seq_len(nrow(N)), function(x){

  E_1_1_1_1(n = n, p1 = N[x,i,1] / n, p2 = N[x,i,2] / n,
            p3 = N[x,i,3] / n, p4 = N[x,ip1,1] / n)
}))

## E(N_u^2 N_v)
mean(N[,i,1]^2 * N[,i,2])
summary(temp <- sapply(seq_len(nrow(N)), function(x){

  E_2_1(n = n, p1 = N[x,i,1] / n, p2 = N[x,i,2] / n)
}))

## E(N_u N_v N_w)
mean(N[,i,1] * N[,i,2] * N[,i,3])
summary(temp <- sapply(seq_len(nrow(N)), function(x){

  E_1_1_1(n = n, p1 = N[x,i,1] / n,
          p2 = N[x,i,2] / n, p3 = N[x,i,3] / n)
}))

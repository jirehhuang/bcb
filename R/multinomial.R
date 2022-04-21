######################################################################
## Multinomial expectations
######################################################################

# E(N_u^k)

E_k <- function(n, p, k = 1){

  if (k == 1){

    n * p

  } else if (k == 2){

    # n * (n-1) * p^2 + E_k(n = n, p = p, k = 1)

    n * p * (1 + (n-1) * p)

  } else if (k == 3){

    # n * (n-1) * (
    #   (n-2) * p^3 + 2 * p^2
    # ) + E_k(n = n, p = p, k = 2)

    n * p * (1 + (n-1) * p * (3 + (n-2) * p))

  } else if (k == 4){

    # n * (n-1) * (n-2) * (
    #   (n-3) * p^4 + 3 * p^3
    # ) + 2 * n * (n-1) * (
    #   (n-2) * p^3 + 2 * p^2
    # ) + E_k(n = n, p = p, k = 3)

    n * p * (1 + (n-1) * p * (7 + (n-2) * p * (6 + (n-3) * p)))

  } else{

    debug_cli(TRUE, cli::cli_abort,
              "k = {k} not supported")
  }
}



# E(N_u^2, N_v^2)

E_2_2 <- function(n, p1, p2){

  # n*(n-1)*p1*p2 + n*(n-1)*(n-2)*(n-3)*p1^2*p2^2 + n*(n-1)*(n-2)*p1*p2^2 + n*(n-1)*(n-2)*p1^2*p2  # verified

  n * (n-1) * p1 * p2 * (1 + (n-2) * (p1 + p2 + (n-3) * p1 * p2))
}



# E(N_u^3, N_v)

E_3_1 <- function(n, p1, p2){

  # n*(n-1)*p1*p2+n*(n-1)*(n-2)*(n-3)*p1^3*p2+3*n*(n-1)*(n-2)*p1^2*p2  # verified

  n * (n-1) * p1 * p2 * (1 + (n-2) * p1 * (3 + (n-3) * p1))
}



# E(N_u^2, N_v, N_w)

E_2_1_1 <- function(n, p1, p2, p3){

  # n*(n-1)*(n-2)*p1*p2*p3+n*(n-1)*(n-2)*(n-3)*p1^2*p2*p3  # verified

  n * (n-1) * (n-2) * p1 * p2 * p3 * (1 + (n-3) * p1)
}



# E(N_u, N_v, N_w, N_x)

E_1_1_1_1 <- function(n, p1, p2, p3, p4){

  n * (n-1) * (n-2) * (n-3) * p1 * p2 * p3 * p4
}



# E(N_u^2, N_v)

E_2_1 <- function(n, p1, p2){

  # n*(n-1)*p1*p2+n*(n-1)*(n-2)*p1^2*p2  # verified

  n * (n-1) * p1 * p2 * (1 + (n-2) * p1)
}



# E(N_u, N_v, N_w)

E_1_1_1 <- function(n, p1, p2, p3){

  n * (n-1) * (n-2) * p1 * p2 * p3
}



# E(N_u, N_v)

E_1_1 <- function(n, p1, p2){

  n * (n-1) * p1 * p2
}



######################################################################
## Denominator W
######################################################################



# E(W_i)

E_W <- function(n, p, i){

  (p[i,1] + p[i,2])
}



# Var(W_i)

Var_W <- function(n, p, i){

  (p[i,1] + p[i,2]) * (1 - (p[i,1] + p[i,2])) / n
}



# Cov(W_i, W_j)

Cov_W <- function(n, p, i, j){

  -(p[i,1] + p[i,2]) * (p[j,1] + p[j,2]) / n
}



######################################################################
## Numerator M
######################################################################



# E(M_i)

E_M <- function(n, p, i){

  # (E_k(n, p[i,1], 2) + E_1_1(n, p[i,1], p[i,2]) + E_1_1(n, p[i,1], p[i,3])) / n^2  # verified

  p[i,1] * (1 + (n-1) * sum(p[i,])) / n
}



# Var(N_i^2)

Var_aa <- function(n, p, i){

  # E_k(n, p[i,1], 4) - E_k(n, p[i,1], 2)^2  # verified

  n * p[i,1] * (1 + (n-1) * p[i,1] * (7 + (n-2) * p[i,1] * (6 + (n-3) * p[i,1])) -
                  n * p[i,1] * (1 + (n-1) * p[i,1])^2)
}



# Var(N_i pN_i)

Var_ab <- function(n, p, i, prime = 1){

  # E_2_2(n, p[i,1], p[i,1+prime]) - E_1_1(n, p[i,1], p[i,1+prime])^2  # verified

  n * (n-1) * p[i,1] * p[i,1+prime] * (

    1 + (n-2) * (

      p[i,1] + p[i,1+prime] + (n-3) * p[i,1] * p[i,1+prime]

    ) - n * (n-1) * p[i,1] * p[i,1+prime]
  )
}



# Cov(N_i^2, N_i pN_i)

Cov_aa_ab <- function(n, p, i, prime = 1){

  # E_3_1(n, p[i,1], p[i,1+prime]) -
  #   E_k(n, p[i,1], 2) * E_1_1(n, p[i,1], p[i,1+prime])  # verified

  n * (n-1) * p[i,1] * p[i,1+prime] * (

    1 + (n-2) * (

      3 * p[i,1] + (n-3) * p[i,1]^2

    ) - n * p[i,1] * (1 + (n-1) * p[i,1])
  )
}



# Cov(N_i pN_i, N_i ppN_i)

Cov_ab_ac <- function(n, p, i){

  # E_2_1_1(n, p[i,1], p[i,2], p[i,3]) -
  #   E_1_1(n, p[i,1], p[i,2]) * E_1_1(n, p[i,1], p[i,3])  # verified

  n * (n-1) * p[i,1] * p[i,2] * p[i,3] * (

    (n-2) * (1 + (n-3) * p[i,1]) - n * (n-1) * p[i,1]
  )
}



# Var(M_i)

Var_M <- function(n, p, i){

  (Var_aa(n, p, i) + Var_ab(n, p, i, 1) + Var_ab(n, p, i, 2) +
    2 * (Cov_aa_ab(n, p, i, 1) + Cov_aa_ab(n, p, i, 2) + Cov_ab_ac(n, p, i))) / n^4
}



# Cov(N_i^2, N_j^2)

Cov_ii_jj <- function(n, p, i, j){

  # E_2_2(n, p[i,1], p[j,1]) - E_k(n, p[i,1], 2) * E_k(n, p[j,1], 2)  # verified

  n * p[i,1] * p[j,1] * (

    (n-1) * (

      1 + (n-2) * (p[i,1] + p[j,1] + (n-3) * p[i,1] * p[j,1])

    ) - n * (1 + (n-1) * p[i,1]) * (1 + (n-1) * p[j,1])
  )
}



# Cov(N_i^2, N_j pN_j)

Cov_ii_jajb <- function(n, p, i, j, prime = 1){

  # E_2_1_1(n, p[i,1], p[j,1], p[j,1+prime]) -
  #   E_k(n, p[i,1], 2) * E_1_1(n, p[j,1], p[j,1+prime])  # verified

  n * (n-1) * p[i,1] * p[j,1] * p[j,1+prime] * (

    (n-2) * (

      1 + (n-3) * p[i,1]

    ) - n * (1 + (n-1) * p[i,1])
  )
}



# Cov(N_i pN_i, N_j pN_j)

Cov_iaib_jajc <- function(n, p, i, j, prime_i = 1, prime_j = 1){

  # E_1_1_1_1(n, p[i,1], p[i,1+prime_i], p[j,1], p[j,1+prime_j]) -
  #   E_1_1(n, p[i,1], p[i,1+prime_i]) * E_1_1(n, p[j,1], p[j,1+prime_j])  # verified

  n * (n-1) * p[i,1] * p[i,1+prime_i] * p[j,1] * p[j,1+prime_j] * (

    (n-2) * (n-3) - n * (n-1)
  )
}



# Cov(M_i, M_j)

Cov_M <- function(n, p, i, j){

  (Cov_ii_jj(n, p, i, j) +
    Cov_ii_jajb(n, p, i, j, 1) + Cov_ii_jajb(n, p, i, j, 2) +
    Cov_ii_jajb(n, p, j, i, 1) + Cov_ii_jajb(n, p, j, i, 2) +
    Cov_iaib_jajc(n, p, i, j, 1, 1) + Cov_iaib_jajc(n, p, i, j, 1, 2) +
    Cov_iaib_jajc(n, p, i, j, 2, 1) + Cov_iaib_jajc(n, p, i, j, 2, 2)) / n^4
}



######################################################################
## M and W
######################################################################



# Cov(N_i^2, N_i)

Cov_ii_i <- function(n, p, i){

  # E_k(n, p[i,1], 3) - E_k(n, p[i,1], 2) * E_k(n, p[i,1], 1)  # verified

  n * p[i,1] * (

    1 + p[i,1] * ((n-1) * (3 - 2 * p[i,1]) - n)
  )
}



# Cov(N_i^2, pN_i)

Cov_ii_ja <- function(n, p, i, j, prime = 0){

  # E_2_1(n, p[i,1], p[j,1+prime]) -
  #   E_k(n, p[i,1], 2) * E_k(n, p[j,1+prime], 1)  # verified

  n * p[i,1] * p[j,1+prime] * (

    (n-1) * (

      1 + (n-2) * p[i,1]

    ) - n * (1 + (n-1) * p[i,1])
  )
}



# Cov(N_i pN_i, N_i)

Cov_ab_a <- function(n, p, i, prime_a = 0, prime_b = 1){

  # E_2_1(n, p[i,1+prime_a], p[i,1+prime_b]) -
  #   E_1_1(n, p[i,1+prime_a], p[i,1+prime_b]) * E_k(n, p[i,1+prime_a])  # verified

  n * (n-1) * p[i,1+prime_a] * p[i,1+prime_b] * (1 - 2 * p[i,1+prime_a])
}



# Cov(N_i ppN_i, pN_i)

Cov_iaib_jc <- function(n, p, i, j, prime_i = 1, prime_j = 0){

  # E_1_1_1(n, p[i,1], p[i,1+prime_i], p[j,1+prime_j]) -
  #   E_1_1(n, p[i,1], p[i,1+prime_i]) * E_k(n, p[j,1+prime_j], 1)  # verified

  -2 * n * (n-1) * p[i,1] * p[i,1+prime_i] * p[j,1+prime_j]
}



# Cov(M_i, W_i)

Cov_M_W <- function(n, p, i){

  (Cov_ii_i(n, p, i) + Cov_ii_ja(n, p, i, i, 1) +
    Cov_ab_a(n, p, i, 0, 1) + Cov_ab_a(n, p, i, 1, 0) + Cov_ab_a(n, p, i, 0, 2) +
    Cov_iaib_jc(n, p, i, i, 2, 1)) / n^3
}



# Cov(M_i, W_j)

Cov_Mi_Wj <- function(n, p, i, j){

  # (Cov_ii_ja(n, p, i, j, 0) + Cov_ii_ja(n, p, i, j, 1) +
  #   Cov_iaib_jc(n, p, i, j, 1, 0) + Cov_iaib_jc(n, p, i, j, 1, 1) +
  #   Cov_iaib_jc(n, p, i, j, 2, 0) + Cov_iaib_jc(n, p, i, j, 2, 1)) / n^3  # verified

  (p[i,1] * (p[j,1] + p[j,2]) * (

    (n-1) * (1 + (n-2) * p[i,1]) - n * (1 + (n-1) * p[i,1])

  ) - 2 * (n-1) * p[i,1] * ((p[i,2] + p[i,3]) * (p[j,1] + p[j,2]))) / n^2
}



######################################################################
## Q
######################################################################



# Q = q(M, W)

q <- function(x, M_plus = 0, W_plus = 0){

  (x[1] + M_plus) / (x[2] + W_plus)
}
dqdM <- function(x, M_plus = 0, W_plus = 0){

  1 / (x[2] + W_plus)
}
dqdW <- function(x, M_plus = 0, W_plus = 0){

  -(x[1] + M_plus) / (x[2] + W_plus)^2
}
d2qdM2 <- function(x, M_plus = 0, W_plus = 0){

  0
}
d2qdW2 <- function(x, M_plus = 0, W_plus = 0){

  2 * (x[1] + M_plus) / (x[2] + W_plus)^3
}
d2qdMW <- function(x, M_plus = 0, W_plus = 0){

  -1 / (x[2] + W_plus)^2
}



# E(Q_i)

E_Q <- function(n, p, i, M_plus = 0, W_plus = 0){

  mu <- c(E_M(n, p, i), E_W(n, p, i))

  q(mu, M_plus, W_plus) +
    # d2fdM2(mu) * Var_M(n, p, i) / 2 +  # d2fdM2(mu) = 0
    d2qdW2(mu, M_plus, W_plus) * Var_W(n, p, i) / 2 +
    d2qdMW(mu, M_plus, W_plus) * Cov_M_W(n, p, i)
}



# Var(Q_i)

Var_Q <- function(n, p, i, M_plus = 0, W_plus = 0){

  mu <- c(E_M(n, p, i), E_W(n, p, i))

  dqdM(mu, M_plus, W_plus)^2 * Var_M(n, p, i) +
    dqdW(mu, M_plus, W_plus)^2 * Var_W(n, p, i) +
    2 * dqdM(mu, M_plus, W_plus) * dqdW(mu, M_plus, W_plus) * Cov_M_W(n, p, i)
}



# E(Q_i, Q_j)

E_Qi_Qj <- function(n, p, i, j, M_plus = 0, W_plus = 0){

  mu_i <- c(E_M(n, p, i), E_W(n, p, i))
  mu_j <- c(E_M(n, p, j), E_W(n, p, j))

  q(mu_i, M_plus, W_plus) * q(mu_j, M_plus, W_plus) +
    dqdM(mu_i, M_plus, W_plus) * dqdM(mu_j, M_plus, W_plus) * Cov_M(n, p, i, j) +
    dqdM(mu_i, M_plus, W_plus) * dqdW(mu_j, M_plus, W_plus) * Cov_Mi_Wj(n, p, i, j) +
    dqdW(mu_i, M_plus, W_plus) * dqdM(mu_j, M_plus, W_plus) * Cov_Mi_Wj(n, p, j, i) +
    dqdW(mu_i, M_plus, W_plus) * dqdW(mu_j, M_plus, W_plus) * Cov_W(n, p, i, j)
}



# Cov(Q_i, Q_j)

Cov_Q <- function(n, p, i, j, M_plus = 0, W_plus = 0){

  mu_i <- c(E_M(n, p, i), E_W(n, p, i))
  mu_j <- c(E_M(n, p, j), E_W(n, p, j))

  dqdM(mu_i, M_plus, W_plus) * dqdM(mu_j, M_plus, W_plus) * Cov_M(n, p, i, j) +
    dqdM(mu_i, M_plus, W_plus) * dqdW(mu_j, M_plus, W_plus) * Cov_Mi_Wj(n, p, i, j) +
    dqdW(mu_i, M_plus, W_plus) * dqdM(mu_j, M_plus, W_plus) * Cov_Mi_Wj(n, p, j, i) +
    dqdW(mu_i, M_plus, W_plus) * dqdW(mu_j, M_plus, W_plus) * Cov_W(n, p, i, j)
}



# Compute Var(Pr) from a matrix of generated Q

Q2Var_Pr <- function(Q, n){

  seq_r <- seq_len(ncol(Q))
  sum(sapply(seq_r, function(i){

    var(Q[,i], na.rm = TRUE) +
      2 * sum(unlist(sapply(seq_r[seq_r > i],
                            function(j) cov(Q[,i], Q[,j],
                                            use = "complete.obs"))))
  }))
}



# Compute taylor approximation of Q using M and W

MW2taylorQ <- function(M, W, i, M_plus = 0, W_plus = 0, order = 1){

  mu <- c(mean(M[,i]), mean(W[,i]))

  temp <- q(mu, M_plus, W_plus) +
    (M[,i] - mu[1]) * dqdM(mu, M_plus, W_plus) +
    (W[,i] - mu[2]) * dqdW(mu, M_plus, W_plus)

  if (order > 1){

    temp <- temp +
      # (M[,i] - mu[1])^2 * d2qdM2(mu, M_plus, W_plus) / 2 +  # zero
      (W[,i] - mu[2])^2 * d2qdW2(mu, M_plus, W_plus) / 2 +
      (M[,i] - mu[1]) * (W[,i] - mu[2]) * d2qdMW(mu, M_plus, W_plus)
  }
  return(temp)
}



# Compute Var(Pr) from generated M and W with different order approximation

MW2Var_Pr <- function(M, W, n, M_plus = 0, W_plus = 0){

  seq_r <- seq_len(ncol(M))
  covQ <- sapply(seq_r, function(i){

    sapply(seq_r, function(j){

      if (i == j){

        mu <- c(mean(M[,i]), mean(W[,i]))

        dqdM(mu, M_plus, W_plus)^2 * var(M[,i]) +
          dqdW(mu, M_plus, W_plus)^2 * var(W[,i]) +
          2 * dqdM(mu, M_plus, W_plus) * dqdW(mu, M_plus, W_plus) * cov(M[,i], W[,i])

      } else{

        mu_i <- c(mean(M[,i]), mean(W[,i]))
        mu_j <- c(mean(M[,j]), mean(W[,j]))

        dqdM(mu_i, M_plus, W_plus) * dqdM(mu_j, M_plus, W_plus) * cov(M[,i], M[,j]) +
          ## using more complex approximation for E[Q_i]
          # - (d2qdW2(mu_i, M_plus, W_plus) * var(W[,i]) / 2 +
          #    d2qdMW(mu_i, M_plus, W_plus) * cov(M[,i], W[,i])) *
          # (d2qdW2(mu_j, M_plus, W_plus) * var(W[,j]) / 2 +
          #    d2qdMW(mu_j, M_plus, W_plus) * cov(M[,j], W[,j]))+
          dqdM(mu_i, M_plus, W_plus) * dqdW(mu_j, M_plus, W_plus) * cov(M[,i], W[,j]) +
          dqdW(mu_i, M_plus, W_plus) * dqdM(mu_j, M_plus, W_plus) * cov(M[,j], W[,i]) +
          dqdW(mu_i, M_plus, W_plus) * dqdW(mu_j, M_plus, W_plus) * cov(W[,i], W[,j])
      }
    })
  })
  sum(sapply(seq_r, function(i){

    covQ[i, i] +
      2 * sum(unlist(sapply(seq_r[seq_r > i],
                            function(j) covQ[i, j])))
  })) / n^2
}



# Var(Pr)

Var_Pr <- function(n,
                   p,
                   M_plus = 0,
                   W_plus = 0){

  seq_r <- seq_len(nrow(p))
  sum(sapply(seq_r, function(i){

    Var_Q(n, p, i, M_plus, W_plus) +
      2 * sum(unlist(sapply(seq_r[seq_r > i],
                            function(j) Cov_Q(n, p, i, j, M_plus, W_plus))))
  }))
}



# Empirical boostrap version

boot_Var_Pr <- function(n,
                        p,
                        nrep = 1e3,
                        nprior = 0){

  p <- (p + nprior / (n * prod(dim(p))))
  p <- p / sum(p)

  samples <- rmultinom(n = nrep,
                       size = n, prob = p)

  estimates <- apply(samples, 2, function(x){

    dim(x) <- dim(p)
    rS <- rowSums(x)
    sum(x[,1] * rS / (rS - x[,3]) / n)
  })
  return(var(estimates, na.rm = TRUE))
}



######################################################################
## Miscellaneous
######################################################################



# Convert string to wolfram

string2wolfram <- function(string){

  string <- gsub("pp\\{p\\_i\\}", "z", string)
  string <- gsub("p\\{p\\_i\\}", "y", string)
  string <- gsub("p\\_i", "x", string)

  string <- gsub("pp\\{p\\_j\\}", "w", string)
  string <- gsub("p\\{p\\_j\\}", "v", string)
  string <- gsub("p\\_j", "u", string)

  string <- gsub("left|right|quad|\\&|\\\n| ", "", string)
  string <- gsub("\\[", "(", string)
  string <- gsub("\\]", ")", string)

  return(string)
}



# Convert jpt to p

jpt2p <- function(jpt,
                  nodes,
                  levels){

  jpt <- query_jpt(jpt, target = nodes)
  jpt <- aperm(a = jpt, perm = match(nodes, names(dimnames(jpt))))

  seq_nnodes <- seq_len(length(nodes))
  zmat <- expand.grid(lapply(dim(jpt)[-seq_len(length(levels))], seq_len))
  if (!length(zmat)){

    zmat <- data.frame(1)
    dim(jpt) <- c(dim(jpt), 1)
    seq_nnodes <- seq_len(length(nodes) + 1)
  }
  p <- t(apply(zmat, 1, function(z){

    ## y, x, z
    idx1 <- idx2 <- idx3 <- as.list(c(levels, unname(z)))

    ## =/=y, x, z
    idx2[[1]] <- setdiff(seq_len(dim(jpt)[1]), idx2[[1]])

    ## =/=x, z
    idx3[[1]] <- seq_len(dim(jpt)[1])
    idx3[[2]] <- setdiff(seq_len(dim(jpt)[2]), idx2[[2]])

    sapply(list(idx1, idx2, idx3), function(idx){

      sum(abind::asub(x = jpt, idx = idx, dims = seq_nnodes))
    })
  }))
  return(p)
}



# Test Var_Pr()

test_Var_Pr <- function(eg,  # grid of scenarios with seed, r, and n
                        path,
                        nq = 1e4,
                        nboot = 1e3,
                        nrep = 1e3,
                        clear = FALSE,
                        debug = 1){


  nq <- min(nq, 1e6)
  nboot <- min(nboot, 1e6)
  debug_cli(debug, cli::cli_alert,
            c("simulating {nrow(eg)} scenarios with ",
              "nq = {nq}, nboot = {nboot}, and nrep = {nrep}"))
  folder <- file.path(path, sprintf("%s_%s_%s", nq, nboot, nrep))
  dir_check(folder)


  ## combine all saved results
  compile_fn <- function(){

    tryCatch({

      files <- list.files(folder)
      files <- files[grepl(".rds", files)]
      files <- files[!grepl("test_Var_Pr", files)]
      files <- file.path(folder, files)
      df <- do.call(rbind, lapply(files, readRDS))
      saveRDS(object = df,
              file.path(path, sprintf("test_Var_Pr_%s_%s_%s.rds",
                                      nq, nboot, nrep)))
    },
    error = function(err){

      debug_cli(debug, cli::cli_alert_danger,
                "error: {as.character(err)}",
                .envir = environment())
    })
  }
  on.exit(compile_fn(), add = TRUE)


  ## function to use on.exit()
  fn <- function(i){


    ## file setup
    rds <- file.path(folder, sprintf("%s.rds", paste(eg[i,], collapse = "_")))
    if (file.exists(rds)){

      if (clear && is.null(readRDS(file = rds))){

        debug_cli(debug, cli::cli_alert,
                  "deleting file {i}")
        file.remove(rds)
      }
      return(NULL)
    }
    saveRDS(object = NULL, file = rds)
    on.exit(expr = {  # delete if failed before completion

      if (is.null(readRDS(file = rds))){

        file.remove(rds)
      }
    }, add = TRUE)
    start_time <- Sys.time()


    ## simulate
    list2env(eg[i,], envir = environment())
    debug_cli(debug, cli::cli_alert,
              "{i}. seed = {seed}, r = {r}, n = {n}",
              .envir = environment())
    set.seed(seed)


    ## generate network
    # nnodes <- 10
    # nodes <- sprintf("V%g", seq_len(nnodes))
    # target <- nodes[nnodes]
    # repeat{
    #
    #   ## generate network
    #   bn <- bcb:::dkpar_bn(p1 = 3, p2 = 6, d = max(3, r))
    #
    #   ## should not have any disconnected components
    #   ig <- igraph::graph.adjacency(bnlearn::amat(bn))
    #   if (length(igraph::decompose.graph(ig)) > 1){
    #
    #     next
    #   }
    #   ## detect nodes with r parents and select one
    #   if (length(int <- which(colSums(bnlearn::amat(bn))[-nnodes] == r))){
    #
    #     if (length(int) > 1)
    #       int <- sample(int, size = 1)
    #     int <- sprintf("V%g", int)
    #
    #   } else{
    #
    #     next
    #   }
    #   ## convert to dnet
    #   dnet <- bn2dnet(bn = bn, # seed = seed,
    #                   min_levels = 2,
    #                   max_levels = 2,
    #                   marginal_lb = 1e-2,
    #                   ce_lb = 0.05,
    #                   n_attempts = 100,
    #                   time_limit = Inf,
    #                   debug = debug)
    #   break
    # }
    # nodes <- union(c(target, int), nodes)
    # dnet <- bcb:::reorder_bn.fit(dnet, ordering = match(nodes, names(dnet)))
    # amat <- bnlearn::amat(dnet)
    # nodes <- c(target, int, dnet[[int]]$parents)
    ######################################################################
    bn <- bnlearn::empty.graph(nodes = sprintf("V%s", seq_len(r + 2)))
    amat <- bnlearn::amat(bn)
    amat[-seq_len(2), 2] <- amat[2, 1] <- 1L
    amat[-seq_len(2), 1] <- sample(c(0, 1), r,  # confounding
                                   replace = TRUE)
    bnlearn::amat(bn) <- amat
    dnet <- bn2dnet(bn = bn, # seed = seed,
                    min_levels = 2,
                    max_levels = 2,
                    marginal_lb = 1e-2,
                    ce_lb = 0.05,
                    n_attempts = 100,
                    time_limit = Inf,
                    debug = debug)
    nodes <- names(dnet)
    ######################################################################
    jpt <- get_jpt(bn.fit = dnet,
                   nodes = nodes)
    p <- jpt2p(jpt = jpt,
               levels = c(1, 1),
               nodes = nodes)
    # true_Pr <- query_jpt(jpt,
    #                      target = nodes[1],
    #                      given = nodes[2],
    #                      adjust = nodes[-c(1, 2)])[1, 1]
    true_Pr <- sum(p[,1] * rowSums(p) /
                     rowSums(p[,seq_len(2),drop=F]))


    ## generate data
    X <- rmultinom(n = 1e6, size = n, prob = p)
    N <- t(X)
    dim(N) <- c(ncol(X), nrow(p), 3)
    Q <- M <- W <- matrix(0, nrow = ncol(X), ncol = nrow(p))
    M <- N[,,1,drop=F] * (N[,,1,drop=F] + N[,,2,drop=F] + N[,,3,drop=F])
    W <- N[,,1,drop=F] + N[,,2,drop=F]
    dim(M) <- dim(W) <- dim(M)[seq_len(2)]
    Q <- M / W
    Pr <- rowSums(Q) / n
    Phat <- N / n

    n_prior <- 1
    N_prior <- array(1, dim = dim(N))
    N_prior[,,3] <- 2
    N_prior <- N_prior * n_prior / (4 * dim(N)[2])
    Phat1 <- (N + N_prior) / (n + n_prior)


    ## for conciseness
    coverage_fn <- function(se_vec){

      mean(sapply(seq_len(length(se_vec)), function(i){

        (Pr[i] - 2 * se_vec[i]) < true_Pr &&
          true_Pr < (Pr[i] + 2 * se_vec[i])

      }), na.rm = TRUE)
    }


    ## compute estimates
    estimates <- list()

    ## simulated true sampling distribution
    estimates$sampling <- apply(aperm(replicate(n = nboot, p), c(3, 1, 2)), 1,
                                boot_Var_Pr, n = n, nrep = nrep, nprior = 0)

    ## proposed
    estimates$proposed0 <- apply(Phat[seq_len(nq),,,drop=F], 1, Var_Pr,
                                 n = n, M_plus = 0, W_plus = 0)
    estimates$proposed1 <- apply(Phat1[seq_len(nq),,,drop=F], 1, Var_Pr,
                                 n = n, M_plus = 0, W_plus = 0)
    estimates$proposed_ <- ifelse(is.na(estimates$proposed0),
                                  estimates$proposed1, estimates$proposed0)

    ## bootstrap from p-hat estimates
    estimates$bootstrap0 <- apply(Phat[seq_len(nboot),,,drop=F], 1, boot_Var_Pr,
                                  n = n, nrep = nrep, nprior = 0)
    estimates$bootstrap1 <- apply(Phat1[seq_len(nboot),,,drop=F], 1, boot_Var_Pr,
                                  n = n, nrep = nrep)
    estimates$bootstrap_ <- ifelse(is.na(estimates$bootstrap0),
                                   estimates$bootstrap1, estimates$bootstrap0)

    ## treat as a conditional proportion, which is only true if r = 0
    estimates$naive <- Pr * (1 - Pr) / (n * query_jpt(jpt,
                                                      target = nodes[2])[1])


    ## compile and save results
    results <- do.call(rbind, lapply(estimates, function(x){

      len <- length(x)
      x <- x[!is.na(x)]
      temp <- data.frame(matrix(quantile(x, probs = seq(0, 1, 0.25)), nrow = 1))
      names(temp) <- sprintf("quantile%s", c("000", "025", "050", "075", "100"))
      temp <- cbind(data.frame(mean = mean(x), sd = sd(x),
                               coverage = coverage_fn(sqrt(x))),
                    temp, data.frame(na_method = (len - length(x)) / len))
      return(temp)
    }))
    results <- cbind(
      eg[rep(i, nrow(results)),],
      true_Pr = true_Pr,
      min_p = min(p), max_p = max(p),
      mean_p1 = mean(p[,1]), mean_p2 = mean(p[,2]), mean_p3 = mean(p[,3]),
      na_Pr = mean(is.na(Pr)),
      mean_Pr = mean(Pr, na.rm = TRUE),
      var_Pr = var(Pr, na.rm = TRUE),
      mean_sampling = mean(estimates$sampling, na.rm = TRUE),
      sd_sampling = sd(estimates$sampling, na.rm = TRUE),
      method = names(estimates),
      results
    )
    end_time <- Sys.time()
    run_time <- as.numeric(end_time - start_time, units = 'secs')
    results$time <- run_time
    rownames(results) <- NULL


    saveRDS(object = results, file = rds)
    debug_cli(debug, cli::cli_alert_success,
              "completed {i} in {prettyunits::pretty_sec(run_time)}",
              .envir = environment())


    if (i %% 1e2 == 0 || i == nrow(eg))
      compile_fn()
  }
  null <- sapply(seq_len(nrow(eg)),
                 fn)
  return(NULL)
}

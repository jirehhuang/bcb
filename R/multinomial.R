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

  n * (p[i,1] + p[i,2])
}



# Var(W_i)

Var_W <- function(n, p, i){

  n * (p[i,1] + p[i,2]) * (1 - (p[i,1] + p[i,2]))
}



# Cov(W_i, W_j)

Cov_W <- function(n, p, i, j){

  -n * (p[i,1] + p[i,2]) * (p[j,1] + p[j,2])
}



######################################################################
## Numerator M
######################################################################



# E(M_i)

E_M <- function(n, p, i){

  # E_k(n, p[i,1], 2) + E_1_1(n, p[i,1], p[i,2]) + E_1_1(n, p[i,1], p[i,3])  # verified

  n * p[i,1] * (1 + (n-1) * sum(p[i,]))
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

  Var_aa(n, p, i) + Var_ab(n, p, i, 1) + Var_ab(n, p, i, 2) +
    2 * (Cov_aa_ab(n, p, i, 1) + Cov_aa_ab(n, p, i, 2) + Cov_ab_ac(n, p, i))  # verified

  # Var_aa(n, p, i) + Var_ab(n, p, i, 1) + Var_ab(n, p, i, 2) +
  #   2 * (Cov_aa_ab(n, p, i, 1) / p[i,2] * (p[i,2] + p[i,3]) + Cov_ab_ac(n, p, i))  # can result in NaN
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

  Cov_ii_jj(n, p, i, j) +
    Cov_ii_jajb(n, p, i, j, 1) + Cov_ii_jajb(n, p, i, j, 2) +
    Cov_ii_jajb(n, p, j, i, 1) + Cov_ii_jajb(n, p, j, i, 2) +
    Cov_iaib_jajc(n, p, i, j, 1, 1) + Cov_iaib_jajc(n, p, i, j, 1, 2) +
    Cov_iaib_jajc(n, p, i, j, 2, 1) + Cov_iaib_jajc(n, p, i, j, 2, 2)  # verified

  # Cov_ii_jj(n, p, i, j) +
  #   Cov_ii_jajb(n, p, i, j, 1) / p[j,2] * (p[j,2] + p[j,3]) +
  #   Cov_ii_jajb(n, p, j, i, 1) / p[i,2] * (p[i,2] + p[i,3]) +
  #   Cov_iaib_jajc(n, p, i, j, 1, 1) / (p[i,2] * p[j,2]) *
  #   (p[i,2] + p[i,3]) * (p[j,2] + p[j,3])  # can result in NaN
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

  Cov_ii_i(n, p, i) + Cov_ii_ja(n, p, i, i, 1) +
    Cov_ab_a(n, p, i, 0, 1) + Cov_ab_a(n, p, i, 1, 0) + Cov_ab_a(n, p, i, 0, 2) +
    Cov_iaib_jc(n, p, i, i, 2, 1)  # verified

  # Cov_ii_i(n, p, i) + Cov_ii_ja(n, p, i, i, 1) +
  #   Cov_ab_a(n, p, i, 0, 1) / p[i,2] * (p[i,2] + p[i,3]) +
  #   Cov_ab_a(n, p, i, 1, 0) +
  #   Cov_iaib_jc(n, p, i, i, 2, 1)  # can result in NaN
}



# Cov(M_i, W_j)

Cov_Mi_Wj <- function(n, p, i, j){

  # Cov_ii_ja(n, p, i, j, 0) + Cov_ii_ja(n, p, i, j, 1) +
  #   Cov_iaib_jc(n, p, i, j, 1, 0) + Cov_iaib_jc(n, p, i, j, 1, 1) +
  #   Cov_iaib_jc(n, p, i, j, 2, 0) + Cov_iaib_jc(n, p, i, j, 2, 1)  # verified

  n * p[i,1] * (p[j,1] + p[j,2]) * (

    (n-1) * (1 + (n-2) * p[i,1]) - n * (1 + (n-1) * p[i,1])

  ) - 2 * n * (n-1) * p[i,1] * ((p[i,2] + p[i,3]) * (p[j,1] + p[j,2]))
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



# Cov(Q_i, Q_j)

Cov_Q <- function(n, p, i, j, M_plus = 0, W_plus = 0){

  mu_i <- c(E_M(n, p, i), E_W(n, p, i))
  mu_j <- c(E_M(n, p, j), E_W(n, p, j))

  dqdM(mu_i, M_plus, W_plus) * dqdM(mu_j, M_plus, W_plus) * Cov_M(n, p, i, j) +
    dqdM(mu_i, M_plus, W_plus) * dqdW(mu_j, M_plus, W_plus) * Cov_Mi_Wj(n, p, i, j) +
    dqdW(mu_i, M_plus, W_plus) * dqdM(mu_j, M_plus, W_plus) * Cov_Mi_Wj(n, p, j, i) +
    dqdW(mu_i, M_plus, W_plus) * dqdW(mu_j, M_plus, W_plus) * Cov_W(n, p, i, j)
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
  })) / n^2
}



# Empirical boostrap version

boot_Var_Pr <- function(n,
                        p,
                        nreps = 1e3,
                        nprior = 0){

  p <- (p + nprior / n)
  p <- p / sum(p)

  samples <- rmultinom(n = nreps,
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

  jpt <- bcb:::query_jpt(jpt, target = nodes)
  jpt <- aperm(a = jpt, perm = match(nodes, names(dimnames(jpt))))

  seq_nnodes <- seq_len(length(nodes))
  zmat <- expand.grid(lapply(dim(jpt)[-seq_len(length(levels))], seq_len))

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
                        nlarge = 1e5,
                        nreps = 1e3,
                        clear = FALSE,
                        debug = 1){


  debug_cli(debug, cli::cli_alert,
            c("simulating {nrow(eg)} scenarios with ",
              "nlarge = {nlarge} and nreps = {nreps}"))


  ## function to use on.exit()
  fn <- function(i){


    ## file setup
    rds <- file.path(path, sprintf("%s_%s_%s.rds", paste(eg[i,], collapse = "_"),
                                   nlarge, nreps))
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
    }, add = FALSE)
    start_time <- Sys.time()


    ## simulate
    list2env(eg[i,], envir = environment())
    debug_cli(debug, cli::cli_alert,
              "{i}. seed = {seed}, r = {r}, n = {n}",
              .envir = environment())
    set.seed(seed)

    p <- runif(r * 3)
    p <- p / sum(p)
    dim(p) <- c(r, 3)

    X <- rmultinom(n = nlarge, size = n, prob = p)
    N <- t(X)
    dim(N) <- c(ncol(X), nrow(p), 3)
    Phat <- N / n
    Q <- M <- W <- matrix(0, nrow = ncol(X), ncol = nrow(p))
    M <- N[,,1] * (N[,,1] + N[,,2] + N[,,3])
    W <- N[,,1] + N[,,2]
    Q <- M / W
    Pr <- rowSums(Q) / n


    ## compute estimates
    estimates <- list()

    ## simulated true sampling distribution
    estimates$sampling <- apply(aperm(replicate(n = nlarge, p), c(3, 1, 2)), 1,
                                bcb:::boot_Var_Pr,
                                n = n, nreps = nreps, nprior = 1)

    ## proposed
    estimates$proposed00 <- apply(Phat, 1, bcb:::Var_Pr,
                                  n = n, M_plus = 0, W_plus = 0)
    estimates$proposed01 <- apply(Phat, 1, bcb:::Var_Pr,
                                  n = n, M_plus = 0, W_plus = 1)
    estimates$proposed11 <- apply(Phat, 1, bcb:::Var_Pr,
                                  n = n, M_plus = 1, W_plus = 1)

    ## bootstrap from p-hat estimates
    estimates$bootstrap0 <- apply(Phat, 1, bcb:::boot_Var_Pr,
                                  n = n, nreps = nreps, nprior = 0)
    estimates$bootstrap1 <- apply(Phat, 1, bcb:::boot_Var_Pr,
                                  n = n, nreps = nreps, nprior = 1)


    ## compile and save results
    results <- do.call(rbind, lapply(estimates, function(x){

      x <- x[!is.na(x)]
      temp <- data.frame(matrix(quantile(x, probs = seq(0, 1, 0.25)), nrow = 1))
      names(temp) <- sprintf("quantile%s", c("000", "025", "050", "075", "100"))
      temp <- cbind(data.frame(mean = mean(x), sd = sd(x)),
                    temp, data.frame(na_method = (nlarge - length(x)) / nlarge))
      return(temp)
    }))
    results <- cbind(eg[rep(i,6),],
                     sampling = mean(estimates$sampling, na.rm = TRUE),
                     large = var(Pr, na.rm = TRUE), na_large = mean(is.na(Pr)),
                     method = names(estimates), results)
    rownames(results) <- NULL
    saveRDS(object = results, file = rds)


    end_time <- Sys.time()
    run_time <- as.numeric(end_time - start_time, units = 'secs')
    debug_cli(debug, cli::cli_alert_success,
              "completed {i} in {prettyunits::pretty_sec(run_time)}",
              .envir = environment())


    ## combine all saved results
    tryCatch({

      files <- list.files(path)
      files <- files[grepl(".rds", files)]
      files <- files[!grepl("test_Var_Pr", files)]
      files <- file.path(path, files)
      df <- do.call(rbind, lapply(files, readRDS))
      saveRDS(object = df, file.path(path, "test_Var_Pr.rds"))
    },
    error = function(err){

      debug_cli(debug, cli::cli_alert_danger,
                "error: {as.character(err)}",
                .envir = environment())
    })
  }
  sapply(seq_len(nrow(eg)), fn)

  return(NULL)
}

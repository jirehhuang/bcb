library(bcb)
get_projects_dir()
path0 <- file.path(projects_dir, "current", "simulations")
path <- file.path(path0, sprintf("test-var-pr"))
dir_check(path)


## simulation settings
eg <- expand.grid(seed = seq_len(1e4),
                  r = seq(2, 10, 2),
                  n = c(25, 50, 100, 250, 500, 1000))
nlarge <- 1e5
nreps <- 1e3

for (i in seq_len(nrow(eg))){


  ## file setup
  rds <- file.path(path, sprintf("%s.rds", paste(eg[i,], collapse = "_")))
  if (file.exists(rds)){

    next
  }
  saveRDS(object = NULL, file = rds)
  on.exit(expr = {  # delete if failed before completion

    if (is.null(is.null(readRDS(file = rds)))){

      file.remove(rds)
    }
  }, add = FALSE)


  ## simulate
  list2env(eg[i,], envir = environment())
  debug_cli(debug, cli::cli_alert,
            "{i}. seed = {seed}, r = {r}, n = {n}")
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
}

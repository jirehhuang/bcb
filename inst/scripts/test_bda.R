######################################################################
## Test backdoor adjustment
######################################################################


library(bcb)
path0 <- ifelse(get_projects_dir(envir = environment()) == getwd(),
                projects_dir, file.path(projects_dir, "current",
                                        "simulations"))
debug <- 1


######################################################################
## Discrete
######################################################################

path <- file.path(path0, "test_bda-d")
eg <- expand.grid(seed = seq_len(1e3),
                  r = c(0, 1, 2, 3),
                  n = 100 * 2^(0:5))

## Execute simulation
bcb:::test_Var_Pr(eg = eg,
                  path = path,
                  nest = 1e3,
                  nrep = 1e3,
                  clear = FALSE,
                  debug = debug)


######################################################################
## Gaussian
######################################################################

# Test Var(beta)

test_Var_beta <- function(eg,  # grid of scenarios with seed, r, and n
                          path,
                          nlarge = 1e5,
                          coefs = c(0.5, 1),
                          vars = c(0.5, 1),
                          clear = FALSE,
                          debug = 1){


  nlarge <- min(nlarge, 1e6)
  debug_cli(debug, cli::cli_alert,
            "simulating {nrow(eg)} scenarios with nlarge = {nlarge}",
            .envir = environment())
  folder <- file.path(path, nlarge)
  dir_check(folder)


  ## combine all saved results
  compile_fn <- function(){

    tryCatch({

      files <- list.files(folder)
      files <- files[grepl(".rds", files)]
      files <- files[!grepl("test_Var_beta", files)]
      files <- file.path(folder, files)
      df <- do.call(rbind, lapply(files, readRDS))
      saveRDS(object = df,
              file.path(path, sprintf("test_Var_beta_%s.rds", nlarge)))
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
                  "deleting file {i}",
                  .envir = environment())
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
    bn <- bnlearn::empty.graph(nodes = sprintf("V%s", seq_len(r + 2)))
    amat <- bnlearn::amat(bn)
    amat[-seq_len(2), 2] <- amat[2, 1] <- 1L
    amat[-seq_len(2), 1] <- sample(c(0, 1), r,  # confounding
                                   replace = TRUE)
    bnlearn::amat(bn) <- amat
    gnet <- bcb:::bn2gnet(bn = bn, # seed = seed,
                          coefs = coefs, vars = vars,
                          normalize = TRUE, intercept = FALSE)
    nodes <- names(gnet)
    others <- nodes[-c(1, 2)]
    true_beta <- bcb:::bn.fit2effects(bn.fit = gnet)[2,1,1]


    ## generate data
    yX_array_obs <- replicate(nlarge, as.matrix(bnlearn::rbn(x = gnet,
                                                             n = n)))
    n_types <- 1 + length(others) * 2  # obs and {-1, 1} for each node
    prob <- rep(1 / n_types, n_types)
    n_split <- rmultinom(n = 1, size = n * nlarge,
                         prob = prob)
    intervene <- c(
      list(list(n = 0)),
      do.call(c, lapply(others, function(node){

        lapply(c(-1, 1), function(x){

          temp <- list(n = 0, x)
          names(temp) <- c("n",
                           node)
          return(temp)
        })
      }))
    )
    for (j in seq_len(length(n_split))){

      intervene[[j]]$n <- n_split[j]
    }
    yX_array_mix <- as.matrix(bcb::ribn(x = gnet, intervene = intervene))
    yX_array_mix <- yX_array_mix[sample.int(nrow(yX_array_mix)),,
                                 drop = FALSE]
    seq_n <- seq_len(n)
    yX_array_mix <- sapply(seq_len(nlarge), function(j){

      yX_array_mix[seq_n + (j-1) * n,]
    })
    dim(yX_array_mix) <- dim(yX_array_obs)


    ## compute estimates and compile results
    t_star <- qt(p = 0.975, df = n - r - 1)
    data2results <- function(yX_array, method){

      values <- t(apply(yX_array, 3, function(yX){

        values <- numeric(5)
        bcb:::lm_cpp(X = yX[, c(2, match(gnet[[2]]$parents,
                                         nodes)), drop=F],
                     y = yX[, 1], values = values)

        return(values[seq_len(2)])
      }))
      beta <- values[,1]
      se_beta <- values[,2]
      Var_beta <- se_beta^2
      coverage <- mean(sapply(seq_len(length(beta)), function(j){

        (beta[j] - t_star * se_beta[j]) < true_beta &&
          true_beta < (beta[j] + t_star * se_beta[j])

      }), na.rm = TRUE)
      results <- data.frame(matrix(quantile(Var_beta, probs = seq(0, 1, 0.25)),
                                   nrow = 1))
      names(results) <- sprintf("quantile%s", c("000", "025", "050", "075", "100"))

      results <- cbind(eg[i,], results,
                       data.frame(method = method,
                                  true_beta = true_beta,
                                  extra = sum(amat[-seq_len(2), 1]),
                                  mean_beta = mean(beta),
                                  sd_beta = sd(beta),
                                  var_beta = var(beta),
                                  mean_Var_beta = mean(Var_beta),
                                  sd_Var_beta = sd(Var_beta),
                                  coverage = coverage))
      return(results)
    }
    results <- rbind(data2results(yX_array_obs, "obs"),
                     data2results(yX_array_mix, "mix"))
    end_time <- Sys.time()
    run_time <- as.numeric(end_time - start_time, units = 'secs')
    results$time <- run_time
    rownames(results) <- NULL


    ## save results
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

## simulation settings
path <- file.path(path0, sprintf("test_bda-g"))
eg <- expand.grid(seed = seq_len(1e3),
                  r = c(1, 2, 3, 4),
                  n = 10 * 2^(0:5))

## Execute simulation
test_Var_beta(eg = eg,
              path = path,
              nlarge = 1e5,
              coefs = c(0.5, 1),
              vars = c(0.5, 1),
              clear = FALSE,
              debug = debug)

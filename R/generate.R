#' Build data grid.
#'
#' Builds and formats a data grid of simulations.
#'
#' @param network character value indicating network name.
#' @param debug logical value that activates printing debugging output.
#' @return data.frame with grid of network settings.
#' @author Jireh Huang (\email{jirehhuang@@ucla.edu})
#' @examples
#' ## Default data_grid
#' build_data_grid()
#' @export

build_data_grid <- function(network = "survey",
                            data_type = "gaussian",
                            n_dat = 0,
                            n_obs = 0,
                            var_lb = 0.5,
                            var_ub = 1,
                            coef_lb = 0.5,
                            coef_ub = 1,
                            normalize = TRUE,
                            copies = 1,
                            debug = FALSE){

  ## TODO:
  # check functions
  # manual network structures
  # random network structures
  # data_type: discrete
  # max_levels: merge discrete levels
  # manipulate cpts
  # add and remove discrete edges
  # max_in_deg and max_out_deg
  # tiling

  ## dormant parameters
  k <- 1
  avg_deg <- 4
  max_in_deg <- Inf
  max_out_deg <- Inf
  max_levels <- Inf

  data_grid <- expand.grid(normalize = normalize,
                           coef_ub = coef_ub,
                           coef_lb = coef_lb,
                           var_ub = var_ub,
                           var_lb = var_lb,
                           max_levels = max_levels,
                           max_out_deg = max_out_deg,
                           max_in_deg = max_out_deg,
                           avg_deg = avg_deg,
                           n_obs = n_obs,
                           n_dat = n_dat,
                           k = k,
                           data_type = data_type,
                           network = network,
                           id = seq_len(copies),
                           stringsAsFactors = FALSE)

  data_grid <- check_data_grid(data_grid)

  return(data_grid)
}



check_data_grid <- function(data_grid){

  ## TODO:
  # check values

  ## remove duplicates
  data_grid <- data_grid[! duplicated(data_grid), , drop = FALSE]

  ## add missing columns
  data_grid$index <- seq_len(nrow(data_grid))
  if (is.null(data_grid$id))
    data_grid$id <- 1
  nms <- c("index", "id", "seed", "network", "data_type", "k", "n_dat",
           "n_obs", "avg_deg", "max_in_deg", "max_out_deg", "max_levels",
           "var_lb", "var_ub", "coef_lb", "coef_ub", "normalize",
           "n_node", "n_edge", "n_within", "n_between",
           "n_compelled", "n_reversible", "n_params")
  data_grid[setdiff(nms, names(data_grid))] <- 0

  ## rearrange columns
  data_grid <- data_grid[, nms, drop = FALSE]

  return(data_grid)
}



#' @export

generate_data_grid <- function(data_grid = build_data_grid(),
                               path = NULL,
                               n_dat = NULL,
                               n_cores = 1,
                               seed0 = 0,
                               debug = TRUE){

  ## initialize output directory
  if (is.null(path)){

    ## default directory
    # path <- file.path(getwd(), "simulations", Sys.time())
    path <- file.path(path.expand("~"),  # TODO: change this back
                      "Documents/ucla/research/projects/current",
                      "simulations", Sys.time())

  } else if (! grepl(path.expand("~"), path)){

    ## append dir to getwd() if home directory not included
    path <- file.path(getwd(), path)
  }
  dir_check(path)

  ## data_grid
  dg_path <- file.path(path, "data_grid.txt")
  if (!is.null(n_dat) && is.numeric(n_dat)){

    if (file.exists(dg_path)){

      data_grid <- read.table(dg_path, stringsAsFactors = FALSE)
    }
    data_grid$n_dat <- n_dat
  }
  data_grid <- check_data_grid(data_grid)
  data_grid$seed <- seed0 + data_grid$id
  write.table(data_grid, dg_path)

  ## set up parallel execution
  if (n_cores < 1)
    n_cores <- min(parallel::detectCores(), nrow(data_grid))
  n_cores <- round(n_cores)

  debug_cli_sprintf(debug, "", "Generating %g datasets using %g cores",
                    sum(data_grid$n_dat), n_cores)

  mclapply <- if (FALSE && ncores > 1 &&
                  Sys.info()[["sysname"]] %in% c("Windows")){

    ## TODO: eventually support windows
    ## windows workaround with https://github.com/nathanvan/parallelsugar
    # parallelsugar::mclapply

  } else{

    ## reduces to lapply() when ncores = 1
    parallel::mclapply
  }

  ## function for preparing networks
  net_fn <- function(i){

    error <- tryCatch(
      {
        ## prepare data row and directory
        data_row <- data_grid[i, , drop = FALSE]
        data_dir <- file.path(path,
                              sprintf("%g_%g_%s_%g",
                                      data_row$index, data_row$id,
                                      data_row$network, data_row$n_obs))
        dir_check(data_dir)

        ## files already completed
        if (all(c("bn.fit.rds", "effects_list.rds",
                  "true_dag.txt", "true_cpdag.txt",
                  "effects_mat.txt", "order_mat.txt") %in% list.files(data_dir))){

          debug_cli_sprintf(debug, "success",
                            "%g Previously completed preparing network %s",
                            i, data_row$network)

          return(data_row)
        }

        debug_cli_sprintf(debug, "", "%g Preparing network %s",
                          i, data_row$network)

        network <- data_row$network

        if (data_row$data_type == "gaussian"){

          bn.fit <- generate_gnet(x = network,
                                  coefs = c(data_row$coef_lb, data_row$coef_ub),
                                  vars = c(data_row$var_lb, data_row$var_ub),
                                  seed = data_row$seed,
                                  normalize = data_row$normalize,
                                  reorder = TRUE, rename = TRUE)

        } else if (data_row$data_type == "discrete"){

          bn.fit <- load_bn.fit(x = network,
                                reorder = TRUE, rename = TRUE)
        }

        ## true graphs
        true_dag <- bnlearn::amat(bn.fit)
        true_cpdag <- bnlearn::amat(bnlearn::cpdag(bn.fit))

        ## get true effect sizes
        effects_list <- bn.fit2effects(bn.fit, debug = debug)
        effects_mat <- effects_list2mat(effects_list, level = 1)

        ## random orderings
        set.seed(data_row$seed)
        order_mat <- do.call(cbind,
          lapply(seq_len(max(1, data_row$n_dat)), function(x){
            sample(seq_len(data_row$n_node))
          })
        )

        ## update data_row
        data_row <- bn.fit2data_row(bn.fit, data_row)

        ## write files
        write.table(data_row, file.path(data_dir, "data_row.txt"))
        saveRDS(bn.fit, file.path(data_dir, "bn.fit.rds"))
        saveRDS(effects_list, file.path(data_dir, "effects_list.rds"))
        write.table(true_dag, file.path(data_dir, "true_dag.txt"))
        write.table(true_cpdag, file.path(data_dir, "true_cpdag.txt"))
        write.table(effects_mat, file.path(data_dir, "effects_mat.txt"))
        write.table(order_mat, file.path(data_dir, "order_mat.txt"))

        debug_cli_sprintf(debug, "success",
                          "%g Completed preparing network %s",
                          i, data_row$network)

        return(data_row)
      }
      , error = function(err){

        debug_cli_sprintf(TRUE, "danger", "Error in %g: %s", i, err)
        browser()
      }
    )
  }
  data_rows <- mclapply(seq_len(nrow(data_grid)), mc.cores = n_cores,
                        mc.preschedule = FALSE, net_fn)
  data_grid <- as.data.frame(data.table::rbindlist(data_rows))
  write.table(data_grid, dg_path)

  ## function for generating datasets
  gen_fn <- function(i){

    error <- tryCatch(
      {
        ## prepare data row and directory
        data_row <- data_grid[i, , drop = FALSE]
        data_dir <- file.path(path,
                              sprintf("%g_%g_%s_%g",
                                      data_row$index, data_row$id,
                                      data_row$network, data_row$n_obs))
        dir_check(data_dir)

        debug_cli_sprintf(debug, "", "%g Generating %g datasets for network %s",
                          i, data_row$n_dat, data_row$network)

        if (data_row$n_obs <= 0)
          return(NULL)

        true_scores <- data.frame(loglik = numeric(data_row$n_dat),
                                  aic = 0, bic = 0)

        ## generate and write datasets
        for (j in seq_len(data_row$n_dat)){

          if (file.exists(data_file <- file.path(data_dir,
                                                 sprintf("data%g.txt", j)))){

            data <- read.table(data_file)

          } else{

            ## read bn.fit object
            bn.fit <- readRDS(file.path(data_dir, "bn.fit.rds"))

            set.seed(data_row$seed + j)

            if ("bn.fit.gnet" %in% class(bn.fit)){

              data <- bnlearn::rbn(x = bn.fit, n = data_row$n_obs)

            } else if ("bn.fit.dnet" %in% class(bn.fit)){

              attempt <- 1
              repeat{

                data <- bnlearn::rbn(x = bn.fit, n = data_row$n_obs)

                ## check if any with only one discrete level
                invalid <- sum(sapply(data, function(x) var(as.integer(x))) == 0)

                debug_cli_sprintf(debug, ifelse(invalid, "danger", "success"),
                                  "Dataset %g has %g/%g invalid variables on attempt %g",
                                  j, invalid, data_row$n_node, attempt)

                if (invalid == 0)
                  break
                attempt <- attempt + 1
              }
            }
          }
          ## write data
          write.table(data, data_file)

          ## compute log-likelihood
          true_scores[j, "loglik"] <- bnlearn:::logLik.bn.fit(bn.fit, data)
        }
        ## AIC
        true_scores$aic <- true_scores$loglik -
          1 * bnlearn::nparams(bn.fit)

        ## BIC
        true_scores$bic <- true_scores$loglik -
          log(nrow(data)) / 2 * bnlearn::nparams(bn.fit)

        ## write scores
        write.table(true_scores, file.path(data_dir,
                                           sprintf("true_scores%g.txt", j)))

        return(NULL)
      }
      , error = function(err){

        debug_cli_sprintf(TRUE, "danger", "Error: %s", err)
      }
    )
  }
  null <- mclapply(seq_len(nrow(data_grid)), mc.cores = n_cores,
                   mc.preschedule = FALSE, gen_fn)
}



# Generate Gaussian network parameters
#' @export

generate_gnet <- function(x,
                          coefs = c(0.5, 1),
                          vars = c(0.5, 1),
                          seed = 1,
                          normalize = TRUE,
                          reorder = TRUE,
                          rename = TRUE){

  ## TODO: check arguments

  set.seed(seed)

  net <- load_bn.fit(x = x, reorder = reorder, rename = rename)
  gnet <- bnlearn::empty.graph(nodes = bnlearn::nodes(net))
  bnlearn::amat(gnet) <- bnlearn::amat(net)

  ## generate parameters for gnet
  dist <- lapply(gnet$nodes, function(node){

    ## sample coefficients and standard deviations
    params <- list(coef = c(0,  # zero mean
                            sample(c(-1, 1),  # negative or positive
                                   length(node$parents), replace = TRUE) *
                              runif(length(node$parents),  # magnitudes
                                    coefs[1], coefs[2])),
                   sd = runif(1, sqrt(vars[1]), sqrt(vars[2])))

    names(params$coef) <- c("(Intercept)", node$parents)

    return(params)
  })

  ## normalize variances
  if (normalize){

    beta <- wamat(bnlearn::custom.fit(gnet, dist = dist))  # coefficient matrix
    I <- diag(length(dist))  # identity matrix
    Omega <- diag(sapply(dist, `[[`, "sd")^2)  # error variances

    rownames(I) <- colnames(I) <-
      rownames(Omega) <- colnames(Omega) <- names(dist)

    ## visit nodes topologically
    for (node in bnlearn::node.ordering(net)){

      if (length(net[[node]]$parents) == 0){

        ## if no parents, variance 1
        Omega[node, node] <- 1
        dist[[node]]$sd <- 1

      } else{

        ## if there are parents, scale coefficients and error
        ## variances such that the variable variance is 1

        ## TODO: other normalizing strategies

        ## estimate covariance matrix
        Sigma <- solve(t(I - beta)) %*% Omega %*% solve(I - beta)

        ## scale coefficients
        beta[, node] <- beta[, node] / sqrt(Sigma[node, node])
        dist[[node]]$coef[-1] <-
          dist[[node]]$coef[-1] / sqrt(Sigma[node, node])

        ## scale error variance
        Omega[node, node] <- Omega[node, node] / Sigma[node, node]
        dist[[node]]$sd <- dist[[node]]$sd / sqrt(Sigma[node, node])
      }
    }
  }
  gnet <- bnlearn::custom.fit(gnet, dist = dist)

  return(gnet)
}



# Simulate random interventional Bayesian network data
#' @export

ribn <- function(x,
                 n = 0,
                 fix = TRUE,
                 intervene = list(),
                 seed = NULL,
                 debug = FALSE){

  if (!is.null(seed))
    set.seed(seed)

  bnlearn:::check.bn.or.fit(x)

  ## generate observational data if no intervention specified
  if (missing(intervene) || is.null(intervene) || length(intervene) == 0){

    return(bnlearn:::rbn.backend(x = x, n = n, fix = fix, debug = debug))
  }

  debug_cli_sprintf(! class(x)[2] %in% c("bn.fit.gnet", "bn.fit.dnet"),
                    "abort", "Currently only bn.fit.gnet and bn.fit.dnet supported for generating interventional data")

  ## TODO: check intervene

  nodes <- bnlearn::nodes(x)
  n_int <- sum(sapply(intervene, `[[`, "n"))
  n_obs <- ifelse(n > n_int, n - n_int, 0)
  intervene[[length(intervene) + 1]] <-
    list(n = n_obs)  # observational intervention

  ## for each intervention
  data <- lapply(intervene, function(int){

    debug_cli_sprintf(debug, "",
                      "Generating %g samples with %g interventions",
                      int$n, sum(names(int) %in% nodes))

    ## convert to bn_list
    xi <- lapply(x, function(y) lapply(y, function(z) z))

    ## for each intervened node
    for (node in intersect(names(int), nodes)){

      if (class(x)[2] == "bn.fit.gnet"){

        ## mutilate xi, cutting off parents to node
        xi[[node]]$parents <- character(0)

        ## set the intercept to the intervention value
        xi[[node]]$coefficients <- c(`(Intercept)` = int[[node]])

        ## remove error variance
        xi[[node]]$sd <- 0

      } else if (class(x)[2] == "bn.fit.dnet"){

        ## mutilate xi, cutting off parents to node
        xi[[node]]$parents <- character(0)

        ## initialize discrete probabilities as 0
        dim_nms <- dimnames(xi[[node]]$prob)[1]
        xi[[node]]$prob <- rep(0, r <- dim(xi[[node]]$prob)[1])

        if (is.character(int[[node]]) &&
            int[[node]] %in% dim_nms[[1]]){

          ## intervention is a character value indicating category
          xi[[node]]$prob[match(int[[node]], dim_nms[[1]])] <- 1

        } else if (is.numeric(int[[node]])  &&
                   int[[node]] %% 1 == 0 &&
                   int[[node]] >= 1 && int[[node]] <= r){

          ## intervention is a numeric value indicating index
          xi[[node]]$prob[int[[node]]] <- 1

        } else if (int[[node]] == "dirichlet"){

          ## dirichlet intervention
          xi[[node]]$prob <- c(DirichletReg::rdirichlet(1, alpha = rep(1/r, r)))
        }
        ## update dimensions and names
        dim(xi[[node]]$prob) <- length(xi[[node]]$prob)
        dimnames(xi[[node]]$prob) <- dim_nms
      }
    }
    xi <- bn_list2bn.fit(xi)

    ## return interventional data
    bnlearn:::rbn.backend(x = xi, n = int$n, fix = fix, debug = debug)
  })

  ## generally faster than do.call(rbind, )
  data <- as.data.frame(data.table::rbindlist(data))

  return(data)
}

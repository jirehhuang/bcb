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
  # n_levels: merge discrete levels
  # manipulate cpts
  # add and remove discrete edges
  # max_in_degree and max_out_degree
  # tiling

  ## dormant parameters
  data_type <- "continuous"
  k <- 1
  n_dat <- 0
  n_obs <- 0
  max_in_degree <- Inf
  max_out_degree <- Inf
  n_levels <- Inf

  data_grid <- expand.grid(normalize = normalize,
                           coef_ub = coef_ub,
                           coef_lb = coef_lb,
                           var_ub = var_ub,
                           var_lb = var_lb,
                           n_levels = n_levels,
                           max_out_degree = max_out_degree,
                           max_in_degree = max_out_degree,
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
           "n_obs", "max_in_degree", "max_out_degree", "n_levels",
           "var_lb", "var_ub", "coef_lb", "coef_ub", "normalize",
           "n_node", "n_edge", "sparsity", "n_within", "n_between",
           "n_reversible", "n_compelled", "n_params")
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
  if (!is.null(n_dat)){

    if (file.exists(dg_path)){

      data_grid <- read.table(dg_path, stringsAsFactors = FALSE)
      data_grid$n_dat <- n_dat
    }
  }
  data_grid <- check_data_grid(data_grid)
  data_grid$seed <- seed0 + data_grid$id
  debug_cli_sprintf(debug, "", "Writing data_grid to output directory")
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
    ## windows workaround with parLapply()
    # parallelsugar::mclapply
  } else{

    ## reduces to lapply() when ncores = 1
    parallel::mclapply
  }

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

        debug_cli_sprintf(debug, "", "(%g / %g) Preparing network %s",
                          i, nrow(data_grid), data_row$network)

        ## TODO:
        network <- data_row$network
        net <- generate_gnet(x = network,
                             coefs = c(data_row$coef_lb, data_row$coef_ub),
                             vars = c(data_row$var_lb, data_row$var_ub),
                             seed = data_row$seed, normalize = data_row$normalize,
                             reorder = TRUE, rename = TRUE)

        browser()
      }
      , error = function(err){

        debug_cli_sprintf(TRUE, "danger", "Error in %g: %s", i, err)
      }
    )
  }

  gen_fn <- function(i){

    error <- tryCatch(
      {
        ## prepare data row and directory
        data_row <- data_grid[i, , drop = FALSE]
        data_dir <- file.path(path,
                              sprintf("%g_%g_%s_n%g",
                                      data_row$index, data_row$id,
                                      data_row$network, data_row$n_obs))
        dir_check(data_dir)

        debug_cli_sprintf(debug, "", "(%g / %g) Generating $g datasets for network %s",
                          i, nrow(data_grid), data_row$n_dat, data_row$network)

        browser()

        ## TODO: generate datasets
      }
      , error = function(err){

        debug_cli_sprintf(TRUE, "danger", "Error: %s", err)
      }
    )
  }

  null <- mclapply(seq_len(nrow(data_grid)), mc.cores = n_cores,
                   mc.preschedule = FALSE, net_fn)

  ## TODO: rewrite data_grid and data_row

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

# Generate data grid
#' @export

gen_data_grid <- function(data_grid = build_data_grid(),
                          path = NULL,
                          n_dat = NULL,
                          n_cores = 1,
                          seed0 = 0,
                          regenerate = FALSE,
                          cache = 0,
                          debug = 1){

  ## initialize output directory
  path <- check_path(path)
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
  write.table(x = data_grid,  # write initial data_grid settings
              file = file.path(path, "data_grid0.txt"))

  ## set up parallel execution
  if (n_cores < 1)
    n_cores <- min(parallel::detectCores(), nrow(data_grid))
  n_cores <- round(n_cores)

  debug_cli(debug, cli::cli_alert_info,
            "generating {sum(data_grid$n_dat)} datasets using {n_cores} core(s)")

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
                              sprintf("%s%g_%s_%g",
                                      data_row$index, data_row$id,
                                      data_row$network, data_row$n_obs))
        dir_check(data_dir)
        write.table(data_row, file.path(data_dir, "data_row0.txt"))

        ## files already completed
        if (all(c("bn.fit.rds", "effects_array.rds",
                  "true_dag.txt", "true_cpdag.txt",
                  "effects_mat.txt", "order_mat.txt") %in% list.files(data_dir))){

          debug_cli(debug, cli::cli_alert_success,
                    "{i} previously prepared network {data_row$network}",
                    .envir = environment())

          data_row <- read.table(file = file.path(data_dir, "data_row.txt"))

          return(data_row)
        }

        debug_cli(debug, "",
                  "{i} preparing network {data_row$network}",
                  .envir = environment())

        set.seed(data_row$seed)
        rename <- !grepl("bnrpar",
                         data_row$network)

        ## load bn structure
        if (grepl("rand|dkpar", data_row$network)){

          repeat{

            ## generate random graph
            bn.fit <- load_bn.fit(x = data_row$network,
                                  reorder = TRUE, rename = rename)
            nodes <- names(bn.fit)
            if (grepl("randpar", data_row$network)){

              nodes <- nodes[-length(nodes)]
            }
            ## too many parents
            if (any(sapply(lapply(bn.fit[nodes], `[[`, "parents"),
                           length) > data_row$max_in_deg)){

              debug_cli(debug, "",
                        "random structure has too many parents",
                        .envir = environment())
              next
            }
            ## too many children
            if (any(sapply(lapply(bn.fit[nodes], `[[`, "children"),
                           length) > data_row$max_out_deg)){

              debug_cli(debug, "",
                        "random structure has too many parents",
                        .envir = environment())
              next
            }
            ## incorrect number of disconnected components
            ig <- igraph::graph.adjacency(bnlearn::amat(bn.fit))
            if (length(igraph::decompose.graph(ig)) != data_row$k){

              debug_cli(debug, "",
                        "random structure has too many disconnected components",
                        .envir = environment())
              next
            }
            ## success
            break
          }
        } else{

          bn.fit <- load_bn.fit(x = data_row$network,
                                reorder = TRUE, rename = rename)
        }
        if (data_row$data_type == "gaussian"){

          attempt <- 1
          repeat{

            seed <- data_row$seed + (attempt - 1) * nrow(data_grid)
            gnet <- bn2gnet(bn = bn.fit,
                            seed = seed,
                            coefs = c(data_row$coef_lb, data_row$coef_ub),
                            vars = c(data_row$var_lb, data_row$var_ub),
                            normalize = data_row$normalize,
                            intercept = FALSE)

            ## check if invalid
            temp_row <- bn.fit2data_row(gnet, data_row)
            invalid <- temp_row$reg_lb < data_row$reg_lb ||
              temp_row$ri_lb < data_row$ri_lb

            debug_cli(debug, ifelse(invalid, cli::cli_alert_danger, cli::cli_alert_success),
                      c("gnet {ifelse(invalid, 'violates', 'satisfies')} ",
                        "constraints on attempt {attempt} with ",
                        "ri_lb = {format(temp_row$ri_lb, digits = 3, nsmall = 3)}, ",
                        "reg_lb = {format(temp_row$reg_lb, digits = 3, nsmall = 3)}"),
                      .envir = environment())

            if (! invalid){

              bn.fit <- gnet
              data_row$seed <- seed
              break
            }
            attempt <- attempt + 1
          }
        } else if (data_row$data_type == "discrete"){

          if (data_row$normalize ||
              grepl("parallel|chain|sink|rand|par", data_row$network)){

            bn.fit <- process_dnet(bn.fit,
                                   min_levels = 1,
                                   max_levels = Inf,
                                   max_in_deg = data_row$max_in_deg,
                                   max_out_deg = data_row$max_out_deg,
                                   remove_order = "decreasing",
                                   min_cp = 0,
                                   ce_lb = data_row$ce_lb,
                                   rename = rename,
                                   debug = debug)
            ## generate cpts
            attempt <- 1
            repeat{

              seed <- data_row$seed + (attempt - 1) * nrow(data_grid)
              dnet <- bn2dnet(bn = bn.fit,
                              seed = seed,
                              min_levels = data_row$var_lb,
                              max_levels = data_row$var_ub,
                              marginal_lb = data_row$coef_lb,
                              ce_lb = data_row$ce_lb,
                              n_attempts = 10000,
                              time_limit = 120,
                              debug = debug)

              ## check if invalid
              temp_row <- bn.fit2data_row(dnet, data_row)
              invalid <- temp_row$reg_lb < data_row$reg_lb ||
                temp_row$ri_lb < data_row$ri_lb

              debug_cli(debug, ifelse(invalid, cli::cli_alert_danger, cli::cli_alert_success),
                        c("dnet {ifelse(invalid, 'violates', 'satisfies')} ",
                          "constraints on attempt {attempt} with ",
                          "ri_lb = {format(temp_row$ri_lb, digits = 3, nsmall = 3)}, ",
                          "reg_lb = {format(temp_row$reg_lb, digits = 3, nsmall = 3)}"),
                        .envir = environment())

              if (! invalid){

                bn.fit <- dnet
                data_row$seed <- seed
                break
              }
              attempt <- attempt + 1
            }
          } else{

            bn.fit <- process_dnet(bn.fit,
                                   min_levels = data_row$var_lb,
                                   max_levels = data_row$var_ub,
                                   merge_order = "increasing",
                                   max_in_deg = data_row$max_in_deg,
                                   max_out_deg = data_row$max_out_deg,
                                   remove_order = "decreasing",
                                   # min_cp = .Machine$double.eps,
                                   ce_lb = data_row$ce_lb,
                                   rename = rename,
                                   debug = debug)
          }  # end if else normalize
        }  # end if gaussian else if discrete

        ## sample from data
        if (data_row$network == "cytometry"){

          ## load from sparsebn
          require(sparsebn, quietly = TRUE)

          target <- unlist(lapply(cytometry$ivn,
                                  function(x) ifelse(is.null(x), NA, x)))
          target <- unname(c(pka = "PKA", akt = "Akt", pkc = "PKC", pip2 = "PIP2", mek = "Mek")[target])

          if (data_row$data_type == "gaussian"){

            data(cytometryContinuous)
            cytometry <- cytometryContinuous

            data <- cytometry$data
            names(data) <- cytometry_nodes

            # data <- as.data.frame(lapply(data, function(x) x - mean(x)))

            not_na <- !is.na(target)
            obs_means <- colMeans(data[!not_na,])

            ## TODO: Gaussian implementation
            browser()

          } else if (data_row$data_type == "discrete"){

            data(cytometryDiscrete)
            cytometry <- cytometryDiscrete

            data <- cytometry$data
            names(data) <- cytometry_nodes
            # data <- as.data.frame(lapply(data, pmax, 1)) - 1
            data <- as.data.frame(lapply(data, pmin, 1))
            data <- as.data.frame(lapply(data, as.factor), stringsAsFactors = TRUE)
          }

          if (!identical(names(bn.fit), names(data))){

            if (setequal(names(bn.fit), names(data))){

              data <- data[names(bn.fit)]

            } else{

              new_names <- names(bn.fit)
              names(new_names) <- names(data)

              names(data) <- unname(new_names[names(data)])
              target <- unname(new_names[target])
            }
          }
          attr(bn.fit, "data") <- data
          attr(bn.fit, "target") <- target
        }

        ## true graphs
        true_dag <- bnlearn::amat(bn.fit)
        true_cpdag <- bnlearn::amat(bnlearn::cpdag(bn.fit))

        ## get true effects
        effects_array <- bn.fit2effects(bn.fit = bn.fit)
        effects_mat <- effects_array[,,1]

        ## random orderings
        set.seed(data_row$seed)
        order_mat <- do.call(cbind,
          lapply(seq_len(max(1, data_row$n_dat)), function(x){
            sample(seq_len(data_row$n_node))
          })
        )
        ## update data_row
        data_row <- bn.fit2data_row(bn.fit = bn.fit,
                                    data_row = data_row,
                                    effects_array = effects_array)

        ## write files
        write.table(data_row, file.path(data_dir, "data_row.txt"))
        saveRDS(bn.fit, file.path(data_dir, "bn.fit.rds"))
        saveRDS(effects_array, file.path(data_dir, "effects_array.rds"))
        write.table(true_dag, file.path(data_dir, "true_dag.txt"))
        write.table(true_cpdag, file.path(data_dir, "true_cpdag.txt"))
        write.table(effects_mat, file.path(data_dir, "effects_mat.txt"))
        write.table(order_mat, file.path(data_dir, "order_mat.txt"))

        debug_cli(debug, cli::cli_alert_success,
                  "{i} successfully prepared network {data_row$network}",
                  .envir = environment())

        return(data_row)
      }
      , error = function(err){

        debug_cli(TRUE, cli::cli_alert_danger, "error in {i}: {err}",
                  .envir = environment())
        browser()
      }
    )
  }
  data_rows <- mclapply(seq_len(nrow(data_grid)), mc.cores = n_cores,
                        mc.preschedule = FALSE, net_fn)
  data_grid <- as.data.frame(data.table::rbindlist(data_rows))
  write.table(x = data_grid, file = dg_path)

  ## expand data.grid so one dataset per thread
  dataset_grid <- do.call(
    rbind,
    lapply(seq_len(nrow(data_grid)), function(i){

      cbind(dataset = seq_len(data_grid[i, "n_dat"]),
            data_grid[rep(i, data_grid[i, "n_dat"]), ])
    })
  )
  rownames(dataset_grid) <- NULL

  ## TODO: dataset_grid for gen_fn() (tricky with true_scores)

  ## function for generating datasets
  gen_fn <- function(i){

    error <- tryCatch(
      {
        ## prepare data row and directory
        data_row <- data_grid[i, , drop = FALSE]
        data_dir <- file.path(path,
                              sprintf("%s%g_%s_%g",
                                      data_row$index, data_row$id,
                                      data_row$network, data_row$n_obs))
        dir_check(data_dir)

        if (data_row$n_obs <= 0 ||
            data_row$n_dat <= 0)
          return(NULL)

        if (!regenerate &&
            file.exists(file.path(data_dir, sprintf("true_scores.txt"))) &&
            all(sapply(seq_len(data_row$n_dat), function(j)
              file.exists(file.path(data_dir, sprintf("data%g.txt", j)))))){

          debug_cli(debug, cli::cli_alert_success,
                    "{i} already generated {data_row$n_dat} datasets for network {data_row$network}",
                    .envir = environment())

          return(NULL)
        }
        debug_cli(debug, "",
                  "{i} generating {data_row$n_dat} datasets for network {data_row$network}",
                  .envir = environment())

        true_scores <- data.frame(loglik = numeric(data_row$n_dat),
                                  aic = 0, bic = 0)

        ## generate and write datasets
        for (j in seq_len(data_row$n_dat)){

          ## read bn.fit object
          bn.fit <- readRDS(file.path(data_dir, "bn.fit.rds"))

          if (!regenerate &  # to create data_file
              file.exists(data_file <- file.path(data_dir,
                                                 sprintf("data%g.txt", j)))){

            data <- read.table(data_file)

            if ("bn.fit.dnet" %in% class(bn.fit))
              data <- as.data.frame(lapply(data, as.factor))

          } else{

            set.seed(data_row$seed + j)

            if ("bn.fit.gnet" %in% class(bn.fit)){

              data <- ribn(x = bn.fit, n = data_row$n_obs)

            } else if ("bn.fit.dnet" %in% class(bn.fit)){

              attempt <- 1
              repeat{

                data <- ribn(x = bn.fit, n = data_row$n_obs)

                ## check if any with only one discrete level
                invalid <- sum(sapply(data, function(x) var(as.integer(x))) == 0)

                debug_cli(debug, ifelse(invalid, cli::cli_alert_danger, cli::cli_alert_success),
                          "{i} dataset {j} has {invalid} of {data_row$n_node} invalid variables on attempt {attempt}",
                          .envir = environment())

                if (invalid == 0)
                  break
                attempt <- attempt + 1
              }
            }
          }
          ## write data
          write.table(x = data, file = data_file)

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
        write.table(x = true_scores, file = file.path(data_dir,
                                                      sprintf("true_scores.txt")))

        return(NULL)
      }
      , error = function(err){

        debug_cli(TRUE, cli::cli_alert_danger, "error in {i}: {err}",
                  .envir = environment())
        browser()
      }
    )
  }
  null <- mclapply(seq_len(nrow(data_grid)), mc.cores = n_cores,
                   mc.preschedule = FALSE, gen_fn)

  if (cache){

    settings <- list(method = "cache",
                     n_obs = max(data_grid$n_obs), n_int = 0)
    simulate_method(method_num = "",
                    settings = settings,
                    path = path,
                    n_cores = n_cores,
                    resimulate = cache > 1,
                    debug = debug)
  }
}



# Simulate random interventional Bayesian network data
#' @export

ribn <- function(x,
                 n = 0,
                 fix = TRUE,
                 intervene = list(),
                 seed = NULL,
                 debug = 0){

  if (!is.null(seed))
    set.seed(seed)

  ## sample from data
  if (!is.null(data <- attr(x, "data")) &&
      !is.null(target <- attr(x, "target"))){

    if (length(intervene) == 0){

      debug_cli(!any(bool_obs <- is.na(target)), cli::cli_abort,
                "no observational data available to resample")

      data <- data[sample(x = which(bool_obs),
                          size = n, replace = TRUE),, drop = FALSE]
      return(data)

    } else{

      ## TODO: generalize to multiple interventions
      if (length(intervene) > 1 ||
          length(intervene[[1]]) > 2){

        browser()
      }
      int <- intervene[[1]]
      node <- setdiff(names(int), "n")
      value <- ifelse(is.numeric(int[[node]]),
                      levels(data[[node]])[int[[node]]], int[[node]])
      data <- data[sample(x = which(data[[node]] == value & target == node),
                          size = int$n, replace = TRUE),, drop = FALSE]
      return(data)
    }
  }

  bnlearn:::check.bn.or.fit(x)

  ## generate observational data if no intervention specified
  if (missing(intervene) || is.null(intervene) || length(intervene) == 0){

    return(bnlearn:::rbn.backend(x = x, n = n, fix = fix,
                                 debug = debug >= 4))
  }

  debug_cli(! class(x)[2] %in% c("bn.fit.gnet", "bn.fit.dnet"), cli::cli_abort,
            c("currently only bn.fit.gnet and bn.fit.dnet supported for ",
              "generating interventional data"))

  ## TODO: check intervene

  nodes <- bnlearn::nodes(x)
  n_int <- sum(sapply(intervene, `[[`, "n"))
  n_obs <- ifelse(n > n_int, n - n_int, 0)
  intervene[[length(intervene) + 1]] <-
    list(n = n_obs)  # observational intervention

  ## for each intervention
  data <- lapply(intervene, function(int){

    debug_cli(debug, "",
              c("generating {int$n} samples with ",
                "{sum(names(int) %in% nodes)} interventions"),
              .envir = environment())

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
    bnlearn:::rbn.backend(x = xi, n = int$n, fix = fix,
                          debug = debug >= 4)
  })

  ## generally faster than do.call(rbind, )
  data <- as.data.frame(data.table::rbindlist(data))

  return(data)
}



######################################################################
## Initialize and check
######################################################################



# Build data grid.
#
# Builds and formats a data grid of simulations.
#
# @param network character value indicating network name.
# @param debug logical value that activates printing debugging output.
# @return data.frame with grid of network settings.
# @examples
# ## Default data_grid
# build_data_grid()
#' @export

build_data_grid <- function(network = "survey",
                            data_type = "gaussian",
                            n_dat = 0,
                            n_obs = 0,
                            max_in_deg = Inf,
                            max_out_deg = Inf,
                            target = "",
                            ce_lb = 1e-2,  # causal effect lower bound
                            ri_lb = 1e-2,  # reward identifiability lower bound
                            reg_lb = 0,  # regret difference lower bound
                            var_lb = ifelse(data_type == "gaussian", 0.1, 2),
                            var_ub = ifelse(data_type == "gaussian", 0.2, 2),
                            coef_lb = ifelse(data_type == "gaussian", 0.5, 1e-2),
                            coef_ub = 1,
                            normalize = TRUE,
                            copies = 1,
                            debug = 0){

  ## TODO:
  # check functions
  # manual network structures
  # random network structures
  # data_type: discrete
  # var_ub: merge discrete levels
  # manipulate cpts
  # add and remove discrete edges
  # max_in_deg and max_out_deg
  # tiling

  ## dormant parameters
  k <- 1
  avg_deg <- 4

  data_grid <- expand.grid(normalize = normalize,
                           coef_ub = coef_ub,
                           coef_lb = coef_lb,
                           var_ub = var_ub,
                           var_lb = var_lb,
                           reg_lb = reg_lb,
                           ri_lb = ri_lb,
                           ce_lb = ce_lb,
                           target = target,
                           max_out_deg = max_out_deg,
                           max_in_deg = max_in_deg,
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



# Check data grid

check_data_grid <- function(data_grid){

  ## TODO: check values

  ## column names
  nms <- c("index", "id", "seed", "network", "data_type", "k",
           "n_dat", "n_obs", "avg_deg", "max_in_deg", "max_out_deg",
           "target", "ce_lb", "ri_lb", "reg_lb", "var_lb", "var_ub",
           "coef_lb", "coef_ub", "normalize", "n_node", "n_edge", "n_within",
           "n_between", "n_compelled", "n_reversible", "n_params")

  ## remove extra columns
  data_grid <- data_grid[, intersect(names(data_grid), nms)]

  ## remove duplicates
  data_grid <- data_grid[! duplicated(data_grid), , drop = FALSE]

  ## add missing columns
  data_grid$index <- stringr::str_pad(string = seq_len(nrow(data_grid)),
                                      width = nchar(nrow(data_grid)),
                                      side = "left", pad = "0")
  data_grid$index <- sprintf("%s_", data_grid$index)
  if (is.null(data_grid$id))
    data_grid$id <- 1
  data_grid[setdiff(nms, names(data_grid))] <- 0

  ## gaussian coefficients must be at least ce_lb
  data_grid$coef_lb <- ifelse(data_grid$data_type == "gaussian",
                              pmax(data_grid$ce_lb, data_grid$coef_lb),
                              data_grid$coef_lb)

  ## rearrange columns
  data_grid <- data_grid[, nms, drop = FALSE]

  return(data_grid)
}

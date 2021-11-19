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
                            target = "",
                            reg_lb = 0,
                            reg_ub = 1,
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
  # var_ub: merge discrete levels
  # manipulate cpts
  # add and remove discrete edges
  # max_in_deg and max_out_deg
  # tiling

  ## dormant parameters
  k <- 1
  avg_deg <- 4
  max_in_deg <- Inf
  max_out_deg <- Inf

  data_grid <- expand.grid(normalize = normalize,
                           coef_ub = coef_ub,
                           coef_lb = coef_lb,
                           var_ub = var_ub,
                           var_lb = var_lb,
                           reg_ub = reg_ub,
                           reg_lb = reg_lb,
                           target = target,
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

  ## TODO: check values

  ## column names
  nms <- c("index", "id", "seed", "network", "data_type", "k",
           "n_dat", "n_obs", "avg_deg", "max_in_deg", "max_out_deg",
           "target", "reg_lb", "reg_ub", "var_lb", "var_ub",
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

  ## rearrange columns
  data_grid <- data_grid[, nms, drop = FALSE]

  return(data_grid)
}



#' @export

gen_data_grid <- function(data_grid = build_data_grid(),
                          path = NULL,
                          n_dat = NULL,
                          n_cores = 1,
                          seed0 = 0,
                          regenerate = FALSE,
                          recache = FALSE,
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
  write.table(data_grid,  # write initial data_grid settings
              file.path(path, "data_grid0.txt"))

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

        ## files already completed
        if (all(c("bn.fit.rds", "effects_list.rds",
                  "true_dag.txt", "true_cpdag.txt",
                  "effects_mat.txt", "order_mat.txt") %in% list.files(data_dir))){

          debug_cli(debug, cli::cli_alert_success,
                    "{i} previously prepared network {data_row$network}",
                    .envir = environment())

          write.table(data_row, file.path(data_dir, "data_row.txt"))

          return(data_row)
        }

        debug_cli(debug, "",
                  "{i} preparing network {data_row$network}",
                  .envir = environment())

        ## load bn structure
        bn.fit <- load_bn.fit(x = data_row$network,
                              reorder = TRUE, rename = TRUE)

        if (data_row$data_type == "gaussian"){

          attempt <- 1
          repeat{

            seed <- data_row$seed + (attempt - 1) * nrow(data_grid)
            gnet <- bn2gnet(bn = bn.fit,
                            seed = seed,
                            coefs = c(data_row$coef_lb, data_row$coef_ub),
                            vars = c(data_row$var_lb, data_row$var_ub),
                            normalize = data_row$normalize)

            ## check if invalid
            temp_mat <- effects_list2mat(bn.fit2effects(gnet, debug = debug))
            temp_row <- bn.fit2data_row(gnet, data_row, temp_mat)
            invalid <- temp_row$reg_lb < data_row$reg_lb ||
              temp_row$reg_ub > data_row$reg_ub

            debug_cli(debug, ifelse(invalid, cli::cli_alert_danger, cli::cli_alert_success),
                      "gnet {ifelse(invalid, 'violates', 'satisfies')} regret constraints on attempt {attempt}",
                      .envir = environment())

            if (! invalid){

              bn.fit <- gnet
              data_row$seed <- seed
              break
            }
            attempt <- attempt + 1
          }
        } else if (data_row$data_type == "discrete"){

          ## TODO: generate cpts
        }

        ## true graphs
        true_dag <- bnlearn::amat(bn.fit)
        true_cpdag <- bnlearn::amat(bnlearn::cpdag(bn.fit))

        ## get true effect sizes
        effects_list <- bn.fit2effects(bn.fit, debug = debug)
        effects_mat <- effects_list2mat(effects_list,
                                        level = 2)  # only applicable for discrete

        ## random orderings
        set.seed(data_row$seed)
        order_mat <- do.call(cbind,
          lapply(seq_len(max(1, data_row$n_dat)), function(x){
            sample(seq_len(data_row$n_node))
          })
        )

        ## update data_row
        data_row <- bn.fit2data_row(bn.fit, data_row, effects_mat)

        ## write files
        write.table(data_row, file.path(data_dir, "data_row.txt"))
        saveRDS(bn.fit, file.path(data_dir, "bn.fit.rds"))
        saveRDS(effects_list, file.path(data_dir, "effects_list.rds"))
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
  write.table(data_grid, dg_path)

  ## expand data.grid so one dataset per thread
  dataset_grid <- do.call(
    rbind,
    lapply(seq_len(nrow(data_grid)), function(i){

      cbind(dataset = seq_len(data_grid[i, "n_dat"]),
            data_grid[rep(i, data_grid[i, "n_dat"]), ])
    })
  )
  rownames(dataset_grid) <- NULL

  ## TODO: dataset_grid for gen_fn() (difficult with true_scores)

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

          if (!regenerate &&
              file.exists(data_file <- file.path(data_dir,
                                                 sprintf("data%g.txt", j)))){

            data <- read.table(data_file)

            if ("bn.fit.dnet" %in% class(bn.fit))
              data <- as.data.frame(lapply(data, as.factor))

          } else{

            set.seed(data_row$seed + j)

            if ("bn.fit.gnet" %in% class(bn.fit)){

              data <- bnlearn::rbn(x = bn.fit, n = data_row$n_obs)

            } else if ("bn.fit.dnet" %in% class(bn.fit)){

              attempt <- 1
              repeat{

                data <- bnlearn::rbn(x = bn.fit, n = data_row$n_obs)

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

  ## function for caching observational rounds
  obs_fn <- function(i){

    error <- tryCatch(
      {
        ## prepare data row and directory
        data_row <- dataset_grid[i, , drop = FALSE]
        data_dir <- file.path(path,
                              sprintf("%s%g_%s_%g",
                                      data_row$index, data_row$id,
                                      data_row$network, data_row$n_obs))
        dir_check(data_dir)

        if (data_row$n_obs <= 0 ||
            data_row$n_dat <= 0)
          return(NULL)

        j <- data_row$dataset
        if (!recache &&
            file.exists(file.path(data_dir,
                                  sprintf("rounds%g.rds", j)))){

          debug_cli(debug, cli::cli_alert_success,
                    "{i} already cached rounds {j} of {data_row$n_dat} for network {data_row$network}",
                    .envir = environment())

          next
        }
        debug_cli(debug, "",
                  "{i} caching rounds {j} of {data_row$n_dat} for network {data_row$network}",
                  .envir = environment())

        ## read bn.fit object
        bn.fit <- readRDS(file.path(data_dir, "bn.fit.rds"))

        settings <- list(method = "cache", target = data_row$target,
                         run = j, n_run = data_row$n_dat, n_obs = data_row$n_obs, n_int = 0)
        settings <- check_settings(bn.fit = bn.fit, settings = settings, debug = debug > 1)

        ## execute bandit
        roundsj <- bandit(bn.fit = bn.fit, settings = settings, debug = debug)

        ## write results in folder roundsj and as roundsj.rds
        rounds_dir <- file.path(data_dir, sprintf("rounds%g", settings$run))
        write_rounds(rounds = roundsj, where = rounds_dir)
        write_rounds(rounds = roundsj, where = sprintf("%s.rds", rounds_dir))

        return(NULL)
      }
      , error = function(err){

        debug_cli(TRUE, cli::cli_alert_danger, "error in {i}: {err}",
                  .envir = environment())
        browser()
      }
    )
  }
  null <- mclapply(seq_len(nrow(dataset_grid)), mc.cores = n_cores,
                   mc.preschedule = FALSE, obs_fn)
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
                "{sum(names(int) %in% nodes} interventions"))

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

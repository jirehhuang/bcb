# Compile path
#' @export

compile_path <- function(path,
                         concise = TRUE,
                         n_cores = 1,
                         debug = 1){

  path <- check_path(path)

  compile_dir <- file.path(path, ifelse(concise, "concise", "complete"))
  dir_check(compile_dir)

  data_grid <- read.table(file.path(path, "data_grid.txt"))
  write.table(data_grid, file.path(compile_dir, "data_grid.txt"))

  if (file.exists(file.path(path, "method_grid.txt"))){

    method_grid <- read.table(file.path(path, "method_grid.txt"))
    write.table(method_grid, file.path(compile_dir, "method_grid.txt"))
  }
  methods <- list.files(path)
  methods <- methods[grepl(paste(avail_methods,
                                 collapse = "|"), methods)]

  ## set up parallel execution
  if (n_cores < 1)
    n_cores <- min(parallel::detectCores(), nrow(data_grid))
  n_cores <- round(n_cores)

  debug_cli(debug && length(methods), cli::cli_alert_info,
            c("compiling {ifelse(concise, 'concise', 'complete')} .rds files ",
              "for {length(methods)} method(s) using {n_cores} core(s)"))

  mclapply <- if (FALSE && ncores > 1 &&
                  Sys.info()[["sysname"]] %in% c("Windows")){

    ## TODO: eventually support windows
    ## windows workaround with https://github.com/nathanvan/parallelsugar
    # parallelsugar::mclapply

  } else{

    ## reduces to lapply() when ncores = 1
    parallel::mclapply
  }

  ## function for writing rds file
  write_fn <- function(method){

    debug_cli(debug, cli::cli_alert_info,
              "compiling {method}", .envir = environment())

    method_dir <- file.path(path, method)

    compiled <- do.call(c, lapply(seq_len(nrow(data_grid)), function(i){

      data_row <- data_grid[i, , drop = FALSE]
      method_data_dir <- file.path(path, method,
                                   nm <- sprintf("%s%g_%s_%g",
                                                 data_row$index, data_row$id,
                                                 data_row$network, data_row$n_obs))

      method_data_rds <- list.files(file.path(method_data_dir, "rds"))
      # method_data_rds <- sprintf("rounds%g.rds", seq_len(data_row$n_dat))
      net_rounds <- sapply(method_data_rds, function(x){

        if (!file.exists(file.path(method_data_dir, "rds", x)))
          return(NULL)

        roundsj <- readRDS(file.path(method_data_dir, "rds", x))
        if (concise){

          cp_dag <- unlist(lapply(avail_bda[-seq_len(2)],
                                  function(x) sprintf("%s_%s", c("dag", "cpdag"), x)))
          roundsj$ji <- as.data.frame(sapply(cp_dag,
                                             function(x) roundsj[[x]]$JI, simplify = FALSE))

          roundsj$selected <-
            roundsj$selected[intersect(names(roundsj$selected),
                                       c("arm", "interventions", "reward", "estimate", "criteria",
                                         "simple_reward", "expected", "correct", "simple_regret",
                                         "cumulative", "mu_est"))]
          roundsj <- roundsj[c("arms", "selected", "beta_true", "mu_true",
                               "mu_est", "criteria", "ji")]
        }
        return(roundsj)

      }, simplify = FALSE)
      names(net_rounds) <- gsub(".rds", "", names(net_rounds))
      net_rounds <- net_rounds[sprintf("rounds%g", seq_len(length(net_rounds)))]

      net_rounds <- list(net_rounds)
      names(net_rounds) <- nm

      return(net_rounds)
    }))
    # attr(compiled, "concise") <- concise
    # attr(compiled, "average") <- 0
    # attr(compiled, "format") <- 0
    if (all(sapply(compiled, length) > 0))
      saveRDS(compiled, file.path(compile_dir, sprintf("%s.rds", method)))
  }
  null <- mclapply(methods, mc.cores = n_cores,
                   mc.preschedule = FALSE, write_fn)
}



# Average compiled rds written by compile_rds()
#' @export

average_compiled <- function(compiled,
                             across_networks = FALSE){

  # concise <- attr(compiled, "concise")
  # average <- attr(compiled, "average")
  # format <- attr(compiled, "format")
  # debug_cli(average > 0 || format > 0, cli::cli_abort,
  #           "compiled must not already be averaged or formatted")

  nms <- names(compiled[[1]][[1]])
  averaged <- sapply(compiled, function(net_rounds){

    net_rounds <- net_rounds[!sapply(net_rounds, is.null)]
    net_averaged <- sapply(nms, function(nm){

      is_numeric <- sapply(compiled[[1]][[1]][[nm]],
                           is.numeric)
      reduced <- Reduce(`+`, lapply(net_rounds, function(roundsj){

        roundsj[[nm]][is_numeric]

      })) / length(net_rounds)

      if (is.matrix(net_rounds[[1]][[nm]])){

        dim(reduced) <- dim(net_rounds[[1]][[nm]])
        dimnames(reduced) <- dimnames(net_rounds[[1]][[nm]])
      }
      return(reduced)

    }, simplify = FALSE)
  }, simplify = FALSE)
  # attr(averaged, "average") <- 1
  # attr(averaged, "format") <- 1

  if (across_networks){

    averaged <- sapply(nms, function(nm){

      is_numeric <- sapply(averaged[[1]][[nm]],
                           is.numeric)
      reduced <- Reduce(`+`, lapply(averaged, function(roundsj){

        roundsj[[nm]][is_numeric]

      })) / length(averaged)

      if (is.matrix(averaged[[1]][[nm]])){

        dim(reduced) <- dim(averaged[[1]][[nm]])
        dimnames(reduced) <- dimnames(averaged[[1]][[nm]])
      }
      return(reduced)

    }, simplify = FALSE)
    # attr(averaged, "average") <- 2
    # attr(averaged, "format") <- 2
  }
  return(averaged)
}



# Convert rounds to df
# Also works for averaged
#' @export

rounds2df <- function(rounds){

  ## TODO: extend beyond concise format

  colnames(rounds$mu_est) <- sprintf("estimate%s", seq_len(ncol(rounds$mu_est)))
  colnames(rounds$criteria) <- sprintf("criteria%s", seq_len(ncol(rounds$criteria)))

  df <- cbind(
    rounds$selected[, setdiff(names(rounds$selected), "interventions")],
    rounds$mu_est, rounds$criteria, rounds$ji
  )
  return(df)
}

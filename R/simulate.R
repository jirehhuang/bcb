# Simulate method
#' @export

simulate_method <- function(method_num,
                            settings,
                            path,
                            n_cores = 1,
                            resimulate = FALSE,
                            debug = 0){

  ## method
  method <- settings$method
  method_dir <- file.path(path, sprintf("%s%s", method, method_num))
  settings$temp_dir <- file.path(method_dir, "progress")
  dir_check(file.path(method_dir, "progress"))  # create method_dir and add progress

  ## data_grid
  data_grid <- read.table(file.path(path, "data_grid.txt"))
  write.table(data_grid, file.path(method_dir, "data_grid.txt"))

  ## expand data.grid so one dataset per thread
  dataset_grid <- do.call(
    rbind,
    lapply(seq_len(nrow(data_grid)), function(i){

      cbind(dataset = seq_len(data_grid[i, "n_dat"]),
            data_grid[rep(i, data_grid[i, "n_dat"]), ])
    })
  )
  rownames(dataset_grid) <- NULL

  ## set up parallel execution
  if (n_cores < 1)
    n_cores <- min(parallel::detectCores(), nrow(data_grid))
  n_cores <- round(n_cores)

  debug_cli(debug, cli::cli_alert_info,
            "executing {method}{method_num} on {sum(data_grid$n_dat)} datasets using {n_cores} core(s)")

  mclapply <- if (FALSE && ncores > 1 &&
                  Sys.info()[["sysname"]] %in% c("Windows")){

    ## TODO: eventually support windows
    ## windows workaround with https://github.com/nathanvan/parallelsugar
    # parallelsugar::mclapply

  } else{

    ## reduces to lapply() when ncores = 1
    parallel::mclapply
  }

  ## function for simulating
  sim_fn <- function(i){

    error <- tryCatch(
      {
        ## prepare data row and directory
        data_row <- dataset_grid[i, , drop = FALSE]
        data_dir <- file.path(path,
                              nm <- sprintf("%s%g_%s_%g",
                                            data_row$index, data_row$id,
                                            data_row$network, data_row$n_obs))
        method_data_dir <- file.path(method_dir, nm)

        j <- data_row$dataset
        if (!resimulate &  # to create rounds_dir
            file.exists(rounds_rds <- file.path(method_data_dir, "rds",
                                                sprintf("rounds%g.rds", j)))){

          debug_cli(debug, cli::cli_alert_success,
                    "{i} already executed {method}{method_num} on dataset {j} of {data_row$n_dat} for network {data_row$network}",
                    .envir = environment())

          return(NULL)
        }

        ## keep track of whether in progress
        progressi <- file.path(method_dir, "progress",
                               sprintf("progress%g", i))
        if (file.exists(progressi))
          return(NULL)  # skip if in progress or complete
        write.table(x = 0, file = progressi,  # mark as in progress
                    row.names = FALSE, col.names = FALSE)
        ## TODO: delete doesn't always work when manually stop multi-core
        on.exit(expr = {  # delete if failed before completion

          if (!file.exists(rounds_rds)){

            file.remove(progressi)

          } else{

            write.table(x = 1, file = progressi,  # mark as complete
                        row.names = FALSE, col.names = FALSE)
          }
        })
        debug_cli(debug, "",
                  "{i} executing {method}{method_num} on dataset {j} of {data_row$n_dat} for network {data_row$network}",
                  .envir = environment())

        ## settings
        bn.fit <- readRDS(file.path(data_dir, "bn.fit.rds"))
        if (method == "cache"){

          settings$n_obs <- min(settings$n_obs, data_row$n_obs)

        } else{

          cache_data_dir <- gsub(sprintf("%s%s", method, method_num),
                                 "cache", method_data_dir)
          cache_file <- file.path(cache_data_dir, "rds",
                                  sprintf("rounds%g.rds", j))
          if (file.exists(cache_file))
            settings$rounds0 <- readRDS(cache_file)
        }
        settings$target <- data_row$target
        settings$run <- j
        settings$data_obs <- file.path(data_dir,
                                       sprintf("data%g.txt", j))
        # set <- check_settings(settings = settings, bn.fit = bn.fit)

        ## execute bandit
        roundsj <- bandit(bn.fit = bn.fit, settings = settings,
                          seed0 = data_row$seed, debug = debug)

        ## write results in folder roundsj and as roundsj.rds
        write_rounds(rounds = roundsj, where = file.path(method_data_dir, "txt",
                                                         sprintf("rounds%g", j)))
        dir_check(file.path(method_data_dir, "rds"))
        write_rounds(rounds = roundsj, where = rounds_rds)

        # browser()

        ## TODO: progress

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
                   mc.preschedule = FALSE, sim_fn)
}



# Simulate method grid
#' @export

sim_method_grid <- function(method,
                            method_nums,
                            method_grid,
                            path,
                            n_cores = 1,
                            resimulate = FALSE,
                            debug = 0){

  ## write method_grid
  write.table(method_grid, file.path(path, "method_grid.txt"))

  ## execute
  null <- lapply(method_nums, function(method_num){

    settings <- as.list(method_grid[method_num,])
    settings$method <- method

    simulate_method(method_num = method_num,
                    settings = settings,
                    path = path,
                    n_cores = n_cores,
                    resimulate = resimulate,
                    debug = debug)
  })
}




######################################################################
## Initialize and check
######################################################################



# Check method grid
#' @export

check_method_grid <- function(method_grid){

  ## TODO: check values

  ## column names
  nms <- c("target", "n_run", "n_obs", "n_int", "n_ess", "n_t", "int_parents",
           "epsilon", "c", "score", "max_parents", "eta", "borrow")

  ## remove extra columns
  method_grid <- method_grid[, intersect(names(method_grid), nms)]

  ## remove duplicates
  method_grid <- method_grid[! duplicated(method_grid), , drop = FALSE]

  ## add missing columns
  method_grid$index <- stringr::str_pad(string = seq_len(nrow(method_grid)),
                                        width = nchar(nrow(method_grid)),
                                        side = "left", pad = "0")
  data_grid$index <- sprintf("%s_", data_grid$index)
  if (is.null(method_grid$target))
    method_grid$target <- ""
  if (is.null(method_grid$c))
    method_grid$c <- 1
  if (is.null(method_grid$int_parents))
    method_grid$int_parents <- 1
  if (is.null(method_grid$score))
    method_grid$score <- ""
  if (is.null(method_grid$max_parents))
    method_grid$max_parents <- 5
  method_grid[setdiff(nms, names(method_grid))] <- 0

  ## rearrange columns
  method_grid <- method_grid[, nms, drop = FALSE]

  return(method_grid)
}

# Simulate method
#' @export

simulate_method <- function(method_num,
                            settings,
                            path = NULL,
                            n_cores = 1,
                            resimulate = FALSE,
                            debug = 0){

  ## initialize output directory
  path <- check_path(path)
  dir_check(path)

  ## method
  method <- settings$method
  method_dir <- file.path(path, sprintf("%s%s", method, method_num))
  dir_check(file.path(method_dir, "progress"))  # create method_dir and add progress
  settings$temp_dir <- file.path(method_dir, "temp")
  dir_check(settings$temp_dir)

  ## data_grid
  data_grid <- read.table(file.path(path, "data_grid.txt"))
  write.table(x = data_grid, file = file.path(method_dir, "data_grid.txt"))

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
        rounds_rds <- file.path(method_data_dir, "rds",
                                sprintf("rounds%g.rds", j))

        ## keep track of whether in progress
        progressi <- file.path(method_dir, "progress",
                               sprintf("progress%g", i))
        if (!resimulate &&
            file.exists(progressi)){

          debug_cli(debug, cli::cli_alert_success,
                    "{i} already executed {method}{method_num} on dataset {j} of {data_row$n_dat} for network {data_row$network}",
                    .envir = environment())

          return(NULL)  # skip if in progress or complete
        }
        write.table(x = 0, file = progressi,  # mark as in progress
                    row.names = FALSE, col.names = FALSE)
        ## TODO: delete doesn't always work when manually stop multi-core
        on.exit(expr = {  # delete if failed before completion

          if (!file.exists(rounds_rds)){

            if (file.exists(progressi))
              file.remove(progressi)

          } else{

            write.table(x = 1, file = progressi,  # mark as complete
                        row.names = FALSE, col.names = FALSE)
          }
        }, add = TRUE)
        debug_cli(debug, "",
                  "{i} executing {method}{method_num} on dataset {j} of {data_row$n_dat} for network {data_row$network}",
                  .envir = environment())

        ## settings
        bn.fit <- readRDS(file.path(data_dir, "bn.fit.rds"))
        seed0 <- data_row$seed * data_row$n_dat
        if (method == "cache"){

          settings$n_obs <- min(settings$n_obs, data_row$n_obs)
          seed0 <- seed0 - settings$n_obs

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
        ## ensure unique id; reset in bandit()
        set.seed(round(as.numeric(Sys.time()) * 1e6) %%
                   .Machine$integer.max)
        settings$id <- random_id(n = 12)

        ## execute bandit
        roundsj <- bandit(bn.fit = bn.fit, settings = settings,
                          seed0 = seed0, debug = debug)

        ## write results in folder roundsj and as roundsj.rds
        # write_rounds(rounds = roundsj, where = file.path(method_data_dir, "txt",
        #                                                  sprintf("rounds%g", j)))
        dir_check(file.path(method_data_dir, "rds"))
        write_rounds(rounds = roundsj, where = rounds_rds)

        ## write progress
        progress <- get_progress(path = path, data_grid = data_grid)
        write.table(x = progress, file = file.path(path, "progress.txt"),
                    quote = FALSE, row.names = TRUE, col.names = TRUE)

        ## print progress
        k <- trimws(rownames(progress)) == sprintf("%s%s", method,
                                                   method_num)
        last2 <- ncol(progress) - c(1, 0)
        tot <- paste(progress[k, last2], collapse = " + ")
        prog <- unlist(progress[k, -last2])
        names(prog) <- seq_len(length(prog))
        prog <- prog[prog != "1.000"]
        prog <- paste(sapply(names(prog), function(x){

          sprintf("%s. %s", x, prog[x])

        }), collapse = " | ")
        prog <- ifelse(nchar(prog) == 0, tot,
                       sprintf("%s || %s", prog, tot))
        debug_cli(TRUE, cli::cli_alert,
                  "{sprintf('%s%s || %s', method, method_num, prog)}",
                  .envir = environment())

        return(NULL)

      }, error = function(err){

        debug_cli(TRUE, cli::cli_alert_danger, "error in {i}: {err}",
                  .envir = environment())
        browser()
      }
    )
  }
  null <- mclapply(seq_len(nrow(dataset_grid)), mc.cores = n_cores,
                   mc.preschedule = FALSE, sim_fn)

  ## write progress
  progress <- get_progress(path = path, data_grid = data_grid)
  write.table(x = progress, file = file.path(path, "progress.txt"),
              quote = FALSE, row.names = TRUE, col.names = TRUE)
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
  method_grid <- check_method_grid(method_grid = method_grid)
  write.table(x = method_grid, file = file.path(path, "method_grid.txt"))

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



# Get progress

get_progress <- function(path,
                         method = "all",
                         data_grid){

  ## TODO: remove; temporary for development
  # if (!dir.exists(path)){
  #
  #   path <- file.path(get_projects_dir(debug = 0),
  #                     "current","simulations", path)
  # }
  path <- check_path(path)
  methods <- list.files(path)
  methods <- methods[grepl(paste(avail_methods,
                                 collapse = "|"), methods)]
  if (method != "all"){

    methods <- match.arg(method, methods)
  }
  progress <- sapply(methods, function(method){

    temp <- lapply(seq_len(nrow(data_grid)), function(i){

      data_row <- data_grid[i,]
      method_data_dir <- file.path(path, method,
                                   sprintf("%s%g_%s_%g",
                                           data_row$index, data_row$id,
                                           data_row$network, data_row$n_obs))
      ## proportion of completed rds files
      length(list.files(file.path(method_data_dir,
                                  "rds"))) / data_row$n_dat
    })
    temp <- do.call(rbind,
                    temp)
    return(temp)
  })
  ## add average progress for each method (Total)
  progress <- rbind(progress, colMeans(progress))

  ## add proportion of incomplete (in progress) executions
  progress <- rbind(progress, sapply(methods, function(method){

    progress_files <- list.files(file.path(path, method, "progress"))
    progress_files <- progress_files[grepl("progress", progress_files)]
    sum(sapply(progress_files, function(x)
      read.table(file.path(path, method, "progress", x))) == 0)

  }) / sum(data_grid$n_dat))

  ## add average progress for all methods
  progress <- cbind(progress, rowMeans(progress))

  ## format
  width <- 5
  progress <- apply(format(progress, digits = 3, nsmall = 3), 1,
                    function(x) stringr::str_pad(x, width = width,
                                                 side = "right"))
  progress <- as.data.frame(progress)

  nc <- max(nchar(c(methods, "total")))
  rownames(progress) <- stringr::str_pad(c(methods, "total"),
                                         width = nc, side = "right")
  colnames(progress) <- c(
    stringr::str_pad(stringr::str_pad("1", width = nc+2,
                                      side = "left"),
                     width = nc + width + 1, side = "right"),
    stringr::str_pad(seq_len(ncol(progress)-2)[-1],
                     width = width, side = "right"),
    "total", "incomplete"
  )
  return(progress)
}



# Clear path
#' @export

clear_path <- function(path,
                       method = "all",
                       clear_type = c("incomplete", "all"),
                       match_type = c("multiple", "single"),
                       debug = 1){

  clear_type <- match.arg(clear_type)
  match_type <- match.arg(match_type)

  ## TODO: remove; temporary for development
  # if (!dir.exists(path)){
  #
  #   path <- file.path(get_projects_dir(debug = 0),
  #                     "current","simulations", path)
  # }
  path <- check_path(path)
  methods <- list.files(path)
  methods <- methods[grepl(paste(avail_methods,
                                 collapse = "|"), methods)]
  if (method != "all"){

    methods <- switch(match_type,
                      multiple = methods[grepl(method, methods)],
                      single = match.arg(method, methods))
  }
  debug_cli(debug && length(methods), cli::cli_alert_info,
            c("clearing {clear_type} files ",
              "for method(s): {paste(methods, collapse = ', ')}"))

  for (method in methods){

    progress_dir <- file.path(path, method, "progress")
    files <- list.files(progress_dir)
    files <- files[grepl("progress", files)]

    for (file in files){

      progressi <- file.path(progress_dir, file)

      if ((method != "cache" && clear_type == "all") ||
          any(read.table(progressi) == 0)){

        debug_cli(debug, cli::cli_alert,
                  "deleting `{file}` for `{method}`")
        file.remove(progressi)
      }
    }
  }
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
           "epsilon", "c", "mu_0", "nu_0", "b_0", "a_0",
           "score", "restrict", "alpha", "max.sx",
           "max_parents", "eta", "unique_make")

  ## remove extra columns
  method_grid <- method_grid[, intersect(names(method_grid), nms)]

  ## remove duplicates
  method_grid <- method_grid[! duplicated(method_grid), , drop = FALSE]

  ## add missing columns
  method_grid$index <- stringr::str_pad(string = seq_len(nrow(method_grid)),
                                        width = nchar(nrow(method_grid)),
                                        side = "left", pad = "0")
  method_grid$index <- sprintf("%s_", method_grid$index)
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

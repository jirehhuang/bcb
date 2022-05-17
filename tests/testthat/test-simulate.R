debug <- 3



## simulate_method()

## cache
settings <- list(method = "cache",
                 restrict = "pc",
                 alpha = 0.1,
                 max.sx = 1,
                 n_obs = 20,
                 n_int = 10)  # should get changed to 0

networks <- c("asia")

for (network in networks){

  path <- file.path(gsub("/tests.*", "", getwd()),
                    "tests", "temp", sprintf("test-%s", network))

  start_time <- Sys.time()
  simulate_method(method_num = "",
                  settings = settings,
                  path = path,
                  n_cores = 1,
                  resimulate = TRUE,
                  debug = debug)
  end_time <- Sys.time()
  run_time <- as.numeric(end_time - start_time, unit = "secs")
  bcb:::debug_cli(TRUE, cli::cli_alert_success,
                  "Completed in {prettyunits::pretty_sec(run_time)}")
}

browser()

## other methods
settings <- list(method = "random",
                 n_obs = 10, n_int = 10)

networks <- c("asia")

for (network in networks){

  path <- file.path(gsub("/tests.*", "", getwd()),
                    "tests", "temp", sprintf("test-%s", network))

  simulate_method(method_num = 1,
                  settings = settings,
                  path = path,
                  n_cores = 1,
                  resimulate = TRUE,
                  debug = debug)
}

# ## check_method_grid()
# bcb:::avail_methods[-1]
# method_grid <- data.frame()
#
# nms <- c("method", "target", "run", "n_obs", "n_int",
#          "n_ess", "n_t", "int_parents", "optimistic", "epsilon",
#          "c", "score", "max_parents", "threshold", "eta")
# settings <- settings[union(nms, c("nodes", "nnodes", "type", "temp_dir",
#                                   "aps_dir", "mds_dir", "id", "data_obs",
#                                   "bn.fit"))]

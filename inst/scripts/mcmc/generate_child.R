######################################################################
## Generate dkpar datasets
######################################################################


library(bcb)
path0 <- ifelse(get_projects_dir(envir = environment()) == getwd(),
                projects_dir, file.path(projects_dir, "current",
                                        "simulations", "mcmc"))
network <- "child"
n_cores <- -1
n_dat <- 50
debug <- 1
start_time <- Sys.time()


######################################################################
## Discrete
######################################################################

data_grid <- build_data_grid(network = network,
                             data_type = "discrete",
                             n_dat = n_dat,
                             n_obs = 1e5,
                             var_lb = 2,  # min_levels
                             var_ub = 2,  # max_levels
                             normalize = FALSE,
                             copies = 1)

path <- file.path(path0, sprintf("%s-d_0", network))

gen_data_grid(data_grid = data_grid,
              path = path,
              n_dat = n_dat,
              seed0 = 0,
              regenerate = FALSE,
              cache = 0,
              n_cores = n_cores,
              debug = debug)


## Execution time
time <- prettyunits::pretty_sec(as.numeric(Sys.time() - start_time,
                                           unit = "secs"))
debug_cli_sprintf(TRUE, "",
                  "executed in {time}")

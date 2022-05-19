######################################################################
## Generate dkpar datasets
######################################################################


library(bcb)
path0 <- ifelse(get_projects_dir(envir = environment()) == getwd(),
                projects_dir, file.path(projects_dir, "current",
                                        "simulations"))
network <- "dkpar_3_6_3"
n_cores <- 1
debug <- 1
start_time <- Sys.time()


######################################################################
## Discrete
######################################################################

data_grid <- build_data_grid(network = network,
                             data_type = "discrete",
                             n_dat = 10,
                             n_obs = 1e4,
                             max_in_deg = 3,
                             ce_lb = 0.05,
                             ri_lb = 0.05,
                             var_lb = 2,  # min_levels
                             var_ub = 2,  # max_levels
                             coef_lb = 1e-2,  # min_marginal
                             normalize = TRUE,
                             copies = 100)

path <- file.path(path0, sprintf("%s-d_0", network))

gen_data_grid(data_grid = data_grid,
              path = path,
              seed0 = 0,
              regenerate = FALSE,
              cache = 0,
              n_cores = n_cores,
              debug = debug)


######################################################################
## Gaussian
######################################################################

data_grid <- build_data_grid(network = network,
                             data_type = "gaussian",
                             n_dat = 10,
                             n_obs = 1e3,
                             max_in_deg = 3,
                             ce_lb = 0.05,
                             ri_lb = 0.05,
                             var_lb = 0.5,
                             var_ub = 1,
                             coef_lb = 0.5,
                             coef_ub = 1,
                             normalize = TRUE,
                             copies = 100)

path <- file.path(path0, sprintf("%s-g_0", network))

gen_data_grid(data_grid = data_grid,
              path = path,
              seed0 = 0,
              regenerate = FALSE,
              cache = 0,
              n_cores = n_cores,
              debug = debug)


## Execution time
time <- prettyunits::pretty_sec(as.numeric(Sys.time() - start_time,
                                           unit = "secs"))
debug_cli_sprintf(TRUE, "",
                  "executed in {time}")  # 40m 9.2s on i5-9600k

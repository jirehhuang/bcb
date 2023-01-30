######################################################################
## Compile dkpar results
######################################################################


library(bcb)
devtools::load_all(path = file.path(get_projects_dir(debug = 0), "current/bcb"))
path0 <- ifelse(get_projects_dir(envir = environment(), debug = 0) == getwd(),
                projects_dir, file.path(projects_dir, "current",
                                        "simulations", "mcmc"))
network <- "dkpar_3_6_3"
n_cores <- 1
debug <- 1
start_time <- Sys.time()


## Discrete

path <- file.path(path0, sprintf("%s-d_done", network))

compile_path(path = path,
             concise = 2,
             n_cores = n_cores,
             debug = debug)

compiled2results(path = path,
                 concise = 2,
                 format = 2,
                 as_df = TRUE,
                 n_cores = n_cores,
                 debug = debug)

compiled2results(path = path,
                 concise = 2,
                 format = 0,
                 as_df = TRUE,
                 n_cores = n_cores,
                 debug = debug)


## Execution time
time <- prettyunits::pretty_sec(as.numeric(Sys.time() - start_time,
                                           unit = "secs"))
debug_cli_sprintf(TRUE, "",
                  "executed in {time}")

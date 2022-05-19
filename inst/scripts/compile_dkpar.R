######################################################################
## Compile dkpar results
######################################################################


library(bcb)
path0 <- ifelse(get_projects_dir(envir = environment()) == getwd(),
                projects_dir, file.path(projects_dir, "current",
                                        "simulations"))
network <- "dkpar_3_6_3"
n_cores <- 4
debug <- 1
start_time <- Sys.time()


## Discrete and Gaussian

for (x in c("d", "g")){

  path <- file.path(path0, sprintf("%s-%s", network, x))

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
}


## Execution time
time <- prettyunits::pretty_sec(as.numeric(Sys.time() - start_time,
                                           unit = "secs"))
debug_cli_sprintf(TRUE, "",
                  "executed in {time}")  # 2h 55m 56.3s with 4 cores on cluster

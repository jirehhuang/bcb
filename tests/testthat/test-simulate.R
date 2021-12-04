debug <- 3



## simulate_method()
settings <- list(method = "random",
                 n_obs = 10, n_int = 10)

networks <- c("asia")

for (network in networks){

  path <- file.path(gsub("/tests.*", "", getwd()),
                    "tests", "temp", sprintf("test-%s", network))

  simulate_method(method_num = 1,
                  settings = settings,
                  path = path,
                  n_cores = 4,
                  resimulate = TRUE,
                  debug = 3)
}

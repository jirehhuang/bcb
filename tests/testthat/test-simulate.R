debug <- 3



settings <- list(method = "random",
                 n_obs = 10, n_int = 10)
path <- file.path(gsub("/tests.*", "", getwd()),
                  "tests", "temp", "asia")

result <- simulate_method(method_num = 1,
                          settings = settings,
                          path = path,
                          n_cores = 1,
                          seed0 = 0,
                          resimulate = FALSE,
                          debug = 3)

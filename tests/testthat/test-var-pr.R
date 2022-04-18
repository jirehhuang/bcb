library(bcb)
get_projects_dir()
path0 <- file.path(projects_dir, "current", "simulations")
path <- file.path(path0, sprintf("test-var-pr"))
dir_check(path)


## simulation settings
eg <- expand.grid(seed = seq_len(1e3),
                  r = c(1, 2, 4, 8, 16),
                  n = c(25, 50, 100, 250, 500, 1000))


## execute simulation
bcb:::test_Var_Pr(eg = eg,
                  path = path,
                  nq = 1e4,
                  nboot = 1e3,
                  nrep = 1e3,
                  clear = FALSE,
                  debug = 1)

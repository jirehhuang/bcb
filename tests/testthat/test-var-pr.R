library(bcb)
get_projects_dir()
path0 <- file.path(projects_dir, "current", "simulations")
path <- file.path(path0, sprintf("test-var-pr"))
dir_check(path)


## simulation settings
eg <- expand.grid(seed = seq_len(1e3),
                  r = c(0, 1, 2, 3),
                  n = 100 * 2^(0:5))


## execute simulation
bcb:::test_Var_Pr(eg = eg,
                  path = path,
                  nest = 1e3,
                  nrep = 1e3,
                  clear = FALSE,
                  debug = 1)

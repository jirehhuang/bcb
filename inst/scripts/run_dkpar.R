######################################################################
## Execute algorithms on dkpar datasets
######################################################################


## Remove _0 from the generated directories from generate_dkpar.R,
## resulting in folder names dkpar_3_6_3-d and dkpar_3_6_3-g


library(bcb)
network <- "dkpar_3_6_3"
path0 <- ifelse(get_projects_dir(envir = environment()) == getwd(),
                projects_dir, file.path(projects_dir, "current",
                                        "simulations", network))
n_cores <- 1
debug <- 1


######################################################################
## Discrete
######################################################################

method_grid <- data.frame(
  n_int = 5000,
  n_obs = rep(10^2 * 2^seq(0, 5), each = 10),
  initial_n_ess = NA,
  n_t = 1,
  max_cache = 1,
  int_parents = 1,
  epsilon = seq(0.01, 0.1, length.out = 10),
  c = c(sqrt(2)^seq(-6, 0), sqrt(3/2), sqrt(2)^seq(1, 2)),
  mu_0 = 0,
  delta = log(5000)^seq(0, 9),
  ucb_criteria = "c",
  bcb_combine = "conjugate",
  bcb_criteria = "bucb",
  score = "bde",
  restrict = "none",
  alpha = 0.1,
  max.sx = 1,
  max_parents = 3,
  threshold = 1 - 1e-6,
  eta = 0.1,
  minimal = TRUE,
  unique_make = TRUE,
  stringsAsFactors = FALSE
)

path <- file.path(path0, sprintf("%s-d", network))
ARGS <- list()
ARGS <- c(ARGS, list(
  list("bucb", c(1, 6), path),  # 1
  list("ts", c(1), path),
  list("ucb", c(1:10), path),  # 4
  list("cn-bucb", c(1, 6), path),  # 1
  list("cn-ts", c(1), path),
  list("cn-ucb", c(1:10), path),  # 4
  list("bcb-bucb", c(seq(1, 60, 10)), path),
  list("bcb-ts", c(seq(1, 60, 10)), path),
  list("bcb-ucb", c(seq(4, 60, 10)), path)
))


## Execute
lapply(ARGS, function(args){

  if (!length(args)) return(NULL)

  n_dat <- ifelse(grepl("bcb", args[[1]]), 5, 10)

  sim_method_grid(method = args[[1]],
                  method_nums = args[[2]],
                  method_grid = method_grid,
                  path = args[[3]],
                  n_dat = n_dat,
                  n_cores = n_cores,
                  resimulate = FALSE,
                  debug = debug)
})


######################################################################
## Gaussian
######################################################################

method_grid$n_obs <- rep(10 * 2^seq(0, 5), each = 10)
method_grid$score <- "bge"

path <- file.path(path0, sprintf("%s-g", network))
ARGS <- list()
ARGS <- c(ARGS, list(
  list("bucb", c(1, 6), path),  # 1
  list("ts", c(1), path),
  list("ucb", c(1:10), path),  # 5
  list("cn-bucb", c(1, 6), path),  # 1
  list("cn-ts", c(1), path),
  list("cn-ucb", c(1:10), path),  # 6
  list("bcb-bucb", c(seq(1, 60, 10)), path),
  list("bcb-ts", c(seq(1, 60, 10)), path),
  list("bcb-ucb", c(seq(6, 60, 10)), path)
))


## Execute
method_grid$unique_make <- method_grid$unique_make || n_cores > 1
lapply(ARGS, function(args){

  if (!length(args)) return(NULL)

  n_dat <- ifelse(grepl("bcb", args[[1]]), 5, 10)

  sim_method_grid(method = args[[1]],
                  method_nums = args[[2]],
                  method_grid = method_grid,
                  path = args[[3]],
                  n_dat = n_dat,
                  n_cores = n_cores,
                  resimulate = FALSE,
                  debug = debug)
})

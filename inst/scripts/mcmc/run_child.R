######################################################################
## Execute algorithms on large datasets
######################################################################


## Remove _0 from child-d_0 generated by generate_child.R,
## resulting in folder name child-d


library(bcb)
library(bidag2)
path0 <- ifelse(get_projects_dir(envir = environment()) == getwd(),
                projects_dir, file.path(projects_dir, "current",
                                        "simulations", "mcmc"))
network <- "child"
n_dat <- 20
n_cores <- 1
debug <- 1


######################################################################
## Discrete
######################################################################

method_grid <- data.frame(
  n_int = 5000,
  n_obs = rep(10^2 * 2^seq(0, 9), each = 10),
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
  bcb_engine = "exact",
  plus1every = Inf,
  plus1post = 0.05,
  plus1it = 2,
  bidag_type = "order",
  min_iterations = 1e4,
  max_iterations = Inf,
  stepsave = 0,
  burnin = 0.2,
  score = "bde",
  restrict = rep(c("none", "pc"), each = 100),
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
  list("bucb", c(1, 6), path),
  list("ts", c(1), path),
  list("ucb", c(1:10), path),
  list("cn-bucb", c(1, 6), path),
  list("cn-ts", c(1), path),
  list("cn-ucb", c(1:10), path),
  list("bcb-mcmc-bucb", c(151), path),
  list("bcb-mcmc-ts", c(151), path),
  list("bcb-mcmc-ucb", c(154), path)
))


## Execute
lapply(ARGS, function(args){

  if (!length(args)) return(NULL)

  sim_method_grid(method = args[[1]],
                  method_nums = args[[2]],
                  method_grid = method_grid,
                  path = args[[3]],
                  n_dat = n_dat,
                  n_cores = n_cores,
                  resimulate = FALSE,
                  debug = debug)
})

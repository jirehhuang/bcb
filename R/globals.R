# Available Bayesian networks from bnrepository
# See \code{\link[phsl]{bnrepository}}

avail_bnrepository <- c(
  # Discrete Bayesian networks
  "asia", "cancer", "earthquake", "sachs", "survey",
  "alarm", "barley", "child", "insurance", "mildew", "water",
  "hailfinder", "hepar2", "win95pts",
  "munin", "munin1", "munin2", "munin3", "munin4",
  "andes", "diabetes", "link", "pathfinder", "pigs",
  # Gaussian Bayesian networks
  "ecoli70", "magic-niab", "magic-irri", "arth150",
  # Conditional linear Gaussian Bayesian networks
  "healthcare", "sangiovese", "mehra"
)



# Available bandit methods (i.e. algorithms; policies)

avail_methods <- c("cache",     # cache observational rounds
                   "random",    # random interventions
                   "greedy",    # (epsilon)-greedy interventions
                   "ucb",       # upper-confidence bounds
                   "ts",        # thompson sampling
                   "bcb",       # bayesian causal bandit; default bma
                   "bcb-star",  # true underlying dag
                   "bcb-bma",   # bayesian model averaging
                   "bcb-mpg",   # median probability graph
                   "bcb-mds",   # thompson sampling via modular dag sampling
                   "bcb-gies")  # greedy interventional equivalence search



# Available back-door adjustment methods

avail_bda <- c("star",  # true underlying dag
               "bma",   # bayesian model averaging
               "mpg",   # median probability graph
               "mds",   # thompson sampling via modular dag sampling
               "gies")  # greedy interventional equivalence search



# Width of function portion of debugging output

DEBUG_WIDTH <- 20



# Symbols from cli package

green_tick <- "{cli::col_green({cli::symbol$tick})}"
red_cross <- "{cli::col_red({cli::symbol$cross})}"

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

avail_methods <- c(
  "cache",          # cache observational rounds
  "random",         # random interventions
  "bucb",           # bayes-ucb
  "ts",             # thompson sampling
  "greedy",         # epsilon-greedy interventions
  "ucb",            # upper-confidence bounds
  "cn-bucb",        # bucb + int_parents = 2
  "cn-ts",          # ts + int_parents = 2
  "cn-greedy",      # epsilon-greedy + int_parents = 2
  "cn-ucb",         # ucb + int_parents = 2
  "bcb-bucb",       # bma + bucb with delta
  "bcb-ts",         # bma + ts
  "bcb-mds-ts",     # mds + ts
  "bcb-greedy",     # bma + greedy with epsilon
  "bcb-ucb",        # bma + tuned with c
  "bcb-star",       # true underlying dag
  "bcb-bma",        # bayesian model averaging
  "bcb-mpg",        # median probability graph
  "bcb-mds",        # thompson sampling via modular dag sampling
  "bcb-gies",       # greedy interventional equivalence search
  "bcb-eg",         # empty graph
  "bcb-mcmc-bucb",  # bma + bucb with delta
  "bcb-mcmc-ts",    # bma + ts
  "bcb-mcmc-ucb"    # bma + ts
)



# Available back-door adjustment methods

avail_bda <- c("star",  # true underlying dag
               "bma",   # bayesian model averaging
               "mpg",   # median probability graph
               "mds",   # thompson sampling via modular dag sampling; TODO: change to sample
               "gies",  # greedy interventional equivalence search
               "eg")    # empty graph



# Available restrict methods

avail_restrict <- c("none",  # no restriction
                    "star",  # skeleton of true underlying skeleton
                    "pc",    # phsl::ppc(max_groups = 1)
                    "ppc",   # phsl::ppc()
                    "gies",  # pcalg::gies()
                    bnlearn:::constraint.based.algorithms)



# Available ucb methods for determining criteria

avail_ucb_criteria <- c("c",      # default unknown distribution with c
                        "tuned")  # ucb-tuned



# Available bcb methods for combining obs and int estimates

avail_bcb_combine <- c("conjugate",  # default conjugate
                       "average")    # TODO: deprecate weighted average



# Available bcb methods for determining criteria

avail_bcb_criteria <- c("bucb",    # default bayes-ucb with delta
                        "ts",      # thompson sampling
                        "greedy",  # greedy with epsilon
                        "c",       # ucb unknown distribution with c
                        "tuned",   # ucb using variance with c
                        "csd")     # +c*sd



# Available bcb methods for

avail_bcb_engine <- c("exact",  # using bida or mds
                      "mcmc")   # using MCMC with bidag



# Flow cytometry node names

cytometry_nodes <- c("Raf",  # praf
                     "Mek",  # pmek
                     "Plcg",  # plcg
                     "PIP2", "PIP3",
                     "Erk",  # p44/42
                     "Akt",  # pakts472
                     "PKA", "PKC",
                     "P38",  # P38
                     "Jnk")  # pjnk



# Width of function portion of debugging output

DEBUG_WIDTH <- 20



# Symbols from cli package

green_tick <- "{cli::col_green({cli::symbol$tick})}"
red_cross <- "{cli::col_red({cli::symbol$cross})}"

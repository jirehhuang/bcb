debug <- 1
set.seed(1)



bcb:::load_example(eg = "gnet", method = "bcb-bma", network = "asia")
settings$n_obs <- 1e1
settings$n_int <- 0
interventions <- rep("", settings$n_obs)
skel <- 1 * (bnlearn::amat(bn.fit) | t(bnlearn::amat(bn.fit)))



mse <- numeric(100)
for (i in seq_len(length(mse))){

  debug_cli(TRUE, cli::cli_alert, "{i}")

  data <- ribn(x = bn.fit, n = settings$n_obs, debug = debug)
  rounds <- list(arms = bcb:::build_arms(bn.fit, settings, debug),
                 data = data,
                 selected = data.frame(arm = rep(0, nrow(data)),
                                       interventions = interventions),
                 ps = list(),
                 bda = list(),
                 beta_true = bcb:::bn.fit2effects(bn.fit = bn.fit)[, , 1],
                 node_values = bcb:::bn.fit2values(bn.fit))
  rounds$mu_true <- mu_true <- sapply(rounds$arms, function(arm){

    arm$value * rounds$beta_true[arm$node, settings$target]
  })

  scores <- bcb:::compute_scores(data = data, settings = settings,
                                 blmat = 1 - skel,
                                 # blmat = diag(settings$nnodes),
                                 interventions = interventions,
                                 output = TRUE, debug = debug)

  rounds$ps <- bcb:::compute_ps(settings = settings,
                                interventions = interventions, debug = debug)

  rounds$bda <- bcb:::compute_bda(data = data,
                                  settings = settings, rounds = rounds,
                                  target = NULL, debug = debug)

  dag <- bnlearn::amat(bn.fit)
  # dag <- NULL
  type <- "bda"
  mu_bda <- sapply(rounds$arms, function(arm){

    int_index <- match(arm$value, rounds$node_values[[arm$node]])
    bcb:::expect_post(rounds = rounds, dag = dag,
                      from = arm$node, to = settings$target,
                      metric = sprintf("mu%g_%s", int_index, type))
  })
  mse[i] <- mean((mu_true - mu_bda)^2)
}
# mse_recenter <- mse
t.test(x = mse, y = mse_recenter, paired = TRUE)

## conclusion: indistinguishable for larger sample sizes (e.g. n >= 1e2),
##             but non-centered significantly better for smaller sample sizes

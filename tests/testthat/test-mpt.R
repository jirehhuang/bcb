debug <- 3



avail_disc_net <- bcb:::avail_bnrepository[seq_len(24)]

for (network in avail_disc_net){

  network <- "alarm"  # TODO: temporary for debugging

  debug_cli(TRUE, cli::cli_alert,
            "testing `{network}`")

  ## unordered and not renamed
  bn.fit <- load_bn.fit(x = network, reorder = FALSE, rename = TRUE)
  mpt <- bcb:::compute_mpt(bn.fit = bn.fit, reorder = FALSE, ordered = FALSE, debug = debug)
  data <- ribn(x = bn.fit, n = 1e6, seed = 1)

  success <- sapply(mpt, function(x) names(x$marginal)[1])

  probs <- sapply(names(bn.fit), function(x) unname(mpt[[x]]$marginal[success[x]]))
  props <- sapply(names(bn.fit), function(x) mean(data[[x]] == success[x]))

  if (mean((probs - props)^2) > 1e-5){

    worst <- order(abs(probs - props), decreasing = TRUE)
    print(abs(probs - props)[worst][seq_len(5)])

    bnlearn::graphviz.plot(bn.fit)
    browser()
  }
}

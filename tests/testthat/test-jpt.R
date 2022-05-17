debug <- 3



# bn.fit <- load_bn.fit(x = "alarm", reorder = FALSE, rename = TRUE)
# bcb:::get_jpt(bn.fit, nodes = "V16")
#
# jpt <- bcb:::get_jpt(bn.fit, nodes = c("V16", "V32", "V33", "V31"))
# bcb:::query_jpt(jpt, target = "V16")
#
# data <- ribn(bn.fit, n = 1e6, seed = 1)
# table(data$V16) / 1e6
#
# ig <- igraph::graph.adjacency(bnlearn::amat(bn.fit))
#
# igraph::all_simple_paths(graph = ig, from = "V31", to = "V33")
#
# bcb:::get_jpt(bn.fit,
#               nodes = "V9")
# props$V9
#
# table(data$V9, data$V8, data$V35) / 1e-6
# jpt
#
# jpt <- bcb:::get_jpt(bn.fit, nodes = c("V9", "V8", "V35"))
# bcb:::query_jpt(jpt, target = "V9")
#
# bn.fit <- load_bn.fit(x = "barley", reorder = FALSE, rename = TRUE)
# bnlearn::graphviz.plot(bn.fit)
# bcb:::get_jpt(bn.fit, nodes = "V2")
# bcb:::get_jpt(bn.fit, nodes = "V40")
# bn_list <- bcb:::add_j_m_pt(bn.fit)



avail_disc_net <- bcb:::avail_bnrepository[seq_len(24)]

for (network in avail_disc_net[-c(seq_len(5), 7)]){

  network <- "insurance"

  debug_cli(TRUE, cli::cli_alert,
            "testing `{network}`")

  bn.fit <- load_bn.fit(x = network, reorder = FALSE, rename = TRUE)
  bnlearn::graphviz.plot(bn.fit)

  data <- ribn(x = bn.fit, n = 1e6, seed = 1)
  props <- lapply(data, function(x) as.vector(table(x)) / length(x))

  bn_list <- lapply(bn.fit[13], function(node){

    print(node$node)

    node$jpt <- bcb:::get_jpt(bn.fit,
                              nodes = union(node$node,
                                            bnlearn::parents(bn.fit, node = node$node)))

    ## mpt
    node$mpt <- bcb:::get_jpt(bn.fit,
                              nodes = node$node)

    diffs <- (props[[node$node]] - node$mpt)

    if (sum(diffs^2) > 1e-4){

      print(node$mpt)
      print(props[[node$node]])
      browser()
    }
    ## cpt
    prob <- bcb:::query_jpt(node$jpt, target = node$node,
                            given = node$parents)
    if (length(names(dimnames(prob)))){

      prob <- aperm(prob, match(names(dimnames(prob)),
                                names(dimnames(node$prob))))
    }
    ae <- all.equal(prob, node$prob,
                    check.attributes = FALSE)
    if (ae != TRUE){

      print(ae)
      browser()
    }
    return(node)
  })
}

## full jpt
for (network in avail_disc_net[seq_len(5)]){

  debug_cli(TRUE, cli::cli_alert,
            "testing `{network}`")

  bn.fit <- load_bn.fit(x = network, reorder = FALSE, rename = TRUE)
  bnlearn::graphviz.plot(bn.fit)

  debug_cli(TRUE, cli::cli_alert,
            "computing jpt for {length(bn.fit)} nodes")
  jpt <- bcb:::bn.fit2jpt(bn.fit = bn.fit)

  ## mpt
  debug_cli(TRUE, cli::cli_alert,
            "computing mpts for {length(bn.fit)} nodes")
  data <- ribn(x = bn.fit, n = 1e6, seed = 1)
  props <- lapply(data, function(x) as.vector(table(x)) / length(x))
  probs <- sapply(names(bn.fit), function(x){

    bcb:::query_jpt(jpt = jpt, target = x, given = character(0))

  }, simplify = FALSE)
  diffs <- unlist(lapply(names(dimnames(jpt)), function(node){

    props[[node]] - probs[[node]]
  }))
  if (sum(diffs^2) > 1e-5){

    worst <- sort(abs(diffs), decreasing = TRUE)
    print(worst[seq_len(5)])

    browser()
  }
  ## cpt
  for (node in names(dimnames(jpt))){

    debug_cli(TRUE, cli::cli_alert,
              "computing cpt for `{node}`")

    prob <- bcb:::query_jpt(jpt = jpt,
                            target = node,
                            given = bn.fit[[node]]$parents)

    if (length(names(dimnames(prob)))){

      prob <- aperm(prob, match(names(dimnames(prob)),
                                names(dimnames(bn.fit[[node]]$prob))))
    }
    ae <- all.equal(prob, bn.fit[[node]]$prob,
                    check.attributes = FALSE)
    if (ae != TRUE){

      print(ae)
      browser()
    }
  }
}


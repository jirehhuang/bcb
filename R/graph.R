######################################################################
## Functions for manipulating graphs
######################################################################



# Function for summarizing (C)(P)DAGs

eval_graph <- function(est, true, cp = FALSE){

  if (any(class(est) %in% c("bn", "bn.fit")))
    est <- bnlearn::amat(est)

  if (any(class(true) %in% c("bn", "bn.fit")))
    true <- bnlearn::amat(true)

  ## convert to cpdags
  if (cp){

    true <- amat2cpdag(amat = true,
                       nodes = colnames(true), debug = 0)
    est <- amat2cpdag(amat = est,
                      nodes = colnames(est), debug = 0)
  }

  est_undir <- est * t(est)
  est_dir <- est - est_undir
  true_undir <- true * t(true)
  true_dir <- true - true_undir

  TE <- sum(true | t(true)) / 2  # true number of edges
  P <- sum(est) - sum(est_undir) / 2  # predicted

  debug_cli_sprintf(P != (sum(est) - sum(est_undir)/2),
                    "abort", "Inconsistent est an est_cp")

  ## summarize (p)dag
  D <- sum(est_dir)  # directed
  U <- sum(est_undir) / 2  # undirected
  TP <- sum(est_dir * true_dir) +
    sum(est_undir * true_undir) / 2  # true positive
  R_D <- sum(true_dir == 1 & (t(est_dir) == true_dir))  # reversed directed
  R_DU <- sum(true_dir == 1 & true == est_undir)  # reversed dir -> undir
  R_UD <- sum(true_undir == 1 & true_undir == est_dir)  # reversed undir -> dir
  R <- R_D + R_DU + R_UD  # reversed
  FP <- P - TP - R  # false positive

  M <- TE - TP - R  # missing
  SHD <- R + FP + M  # structural hamming distance
  JI <- TP / (TE + P - TP)  # jaccard index

  df <- data.frame(`T`=TE, P=P, D=D, U=U, TP=TP, RD=R_D, RDU=R_DU,
                    RUD=R_UD, FP=FP, SHD = SHD, JI=JI)
  df[which(is.na(df))] <- 0

  if (df[["T"]] == 0 && df$P == 0)
    df$JI <- 1

  return(df)
}



# Convert adjacency matrix to CPDAG
# Circumvent converting amat to a bn, which fails
# if amat contains cycles

amat2cpdag <- function(amat,
                       nodes = colnames(amat),
                       debug = 0){

  if (sum(amat) == 0){

    return(amat)
  }
  mode(amat) <- "integer"
  is_acyclic <-
    bnlearn:::is.acyclic(arcs = bnlearn:::amat2arcs(a = amat, nodes = nodes),
                         nodes = nodes, directed = TRUE, debug = debug >= 3)
  ## phsl version, for remove_invalid
  if (!is_acyclic){

    ## get v-structures
    vsmat <- amat * 0
    vs <- phsl:::vstructs_cpp(amat = amat, vsmat = vsmat,
                              debug = debug >= 3)

    ## apply v-structures
    pdag <- (amat | t(amat)) -  # undirected graph
      (vsmat | t(vsmat)) +  # remove v-structures
      vsmat  # orient v-structures

    ## apply cpdag rules
    ## remove_invalid: remove arc if conflicting R2
    mode(pdag) <- "integer"
    cpdag <- phsl:::apply_cpdag_rules(pdag = pdag, nodes = nodes,
                                      remove_invalid = TRUE,
                                      debug = debug >= 3)
    rownames(cpdag) <-
      colnames(cpdag) <- nodes

    ## TODO: remove; temporary for debugging
    is_acyclic <-
      bnlearn:::is.acyclic(arcs = bnlearn:::amat2arcs(a = cpdag, nodes = nodes),
                           nodes = nodes, directed = TRUE, debug = debug >= 3)
    if (!is_acyclic){

      debug_cli(TRUE, cli::cli_alert_danger,
                "CPDAG still contains cycles")
      browser()
    }
  } else{

    ## bnlearn version
    mode(amat) <- "integer"
    arcs <- bnlearn:::amat2arcs(a = amat, nodes = nodes)
    cpdag = .Call(bnlearn:::call_cpdag, arcs = arcs, nodes = nodes,
                  moral = FALSE, fix = FALSE, wlbl = FALSE,
                  whitelist = NULL, blacklist = NULL,
                  illegal = NULL, debug = debug >= 3)
  }
  return(cpdag)
}

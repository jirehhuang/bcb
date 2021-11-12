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

    bn <- bnlearn::empty.graph(nodes = colnames(true))

    bnlearn::amat(bn) <- true
    true <- bnlearn::amat(bnlearn::cpdag(bn))

    bnlearn::amat(bn) <- est
    est <- bnlearn::amat(bnlearn::cpdag(bn))
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

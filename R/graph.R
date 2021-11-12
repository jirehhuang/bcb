######################################################################
## Functions for manipulating graphs
######################################################################



# Function for summarizing (CP)DAGs

summarize_cp_dag <- function(est_dag, true_dag){

  ## get cpdags
  bn <- bnlearn::empty.graph(nodes = colnames(true_dag))

  bnlearn::amat(bn) <- true_dag
  true_cp <- bnlearn::amat(bnlearn::cpdag(bn))

  bnlearn::amat(bn) <- est_dag
  est_cp <- bnlearn::amat(bnlearn::cpdag(bn))

  est_dag_undir <- est_dag * t(est_dag)
  est_dag_dir <- est_dag - est_dag_undir
  est_cp_undir <- est_cp * t(est_cp)
  est_cp_dir <- est_cp - est_cp_undir
  true_cp_undir <- true_cp * t(true_cp)
  true_cp_dir <- true_cp - true_cp_undir

  T <- sum(true_dag | t(true_dag)) / 2  # true
  P <- sum(est_cp) - sum(est_cp_undir)/2  # predicted

  debug_cli_sprintf(P != (sum(est_dag) - sum(est_dag_undir)/2),
                    "abort", "Inconsistent est_dag an est_cp")

  ## summarize (p)dag
  D <- sum(est_dag_dir)  # directed
  U <- sum(est_dag_undir) / 2  # undirected
  TP <- sum(est_dag_dir * true_dag)  # true positives
  R_D <- sum(true_dag == 1 & (t(est_dag_dir) == true_dag))  # reversed directed
  R_DU <- sum(true_dag == 1 & true_dag == est_dag_undir)  # reversed dir -> undir
  R_UD <- sum(true_dag_undir == 1 & true_dag_undir == est_dag_dir)  # reversed undir -> dir
  R <- R_D + R_DU + R_UD  # reversed
  FP <- P - TP - R  # false positives

  M <- T - TP - R  # missing
  SHD <- R + FP + M  # structural hamming distance
  JI <- TP/(T+P-TP)  # jaccard index

  dag <- data.frame(T=T, P=P, D=D, U=U, TP=TP, RD=R_D, RDU=R_DU,
                    RUD=R_UD, FP=FP, SHD = SHD, JI=JI)
  dag[which(is.na(dag))] <- 0

  ## summarize cpdag
  D <- sum(est_cp_dir)
  U <- sum(est_cp_undir) / 2
  TP <- sum(est_cp_dir * true_cp_dir) + sum(est_cp_undir * true_cp_undir) / 2
  R_D <- sum(true_cp_dir == 1 & (t(est_cp_dir) == true_cp_dir))
  R_DU <- sum(true_cp_dir == 1 & true_cp_dir == est_cp_undir)
  R_UD <- sum(true_cp_undir == 1 & true_cp_undir == est_cp_dir)
  R <- R_D + R_DU + R_UD
  FP <- P-TP-R

  M <- T-TP-R
  SHD <- R+FP+M
  JI <- TP/(T+P-TP)

  cpdag <- data.frame(T=T, P=P, D=D, U=U, TP=TP, RD=R_D, RDU=R_DU,
                      RUD=R_UD, FP=FP, SHD = SHD, JI=JI)
  cpdag[which(is.na(cpdag))] <- 0
  if (cpdag[["T"]] == 0 && cpdag$P == 0) cpdag$JI <- 1

  return(list(dag = dag,
              cpdag = cpdag))
}

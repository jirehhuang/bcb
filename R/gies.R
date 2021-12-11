######################################################################
## Functions for efficiently executing gies
######################################################################



# Estimate adjacency matrix using gies

estimate_gies <- function(rounds,
                          settings,
                          interventions,
                          dag = TRUE,
                          debug = 0){

  debug_cli(debug >= 2, cli::cli_alert_info,
            "estimating graph with gies")

  nodes <- settings$nodes
  if (any(rounds$blmat == 1)){

    data <- rounds$data[seq_len(settings$n_obs),]
    score <- new("bnlearn_score", data = data, interventions = interventions,
                 nodes = nodes, score = settings$score, extra.args =
                   bnlearn:::check.score.args(score = settings$score,
                                              network = settings$bn.fit,
                                              data = data, extra.args = list()))
  } else{

    score <- new("lookup_score", ps = rounds$ps,
                 interventions = interventions,
                 nodes = nodes)
  }
  ## TODO: score that tries to lookup and uses bnlearn if failed

  gies <- pcalg::gies(score = score, maxDegree = settings$max_parents,
                      fixedGaps = rounds$blmat, iterate = TRUE,
                      adaptive = ifelse(any(rounds$blmat == 1),
                                        "vstructures", "none"),
                      phase = c("forward", "backward", "turning"),
                      verbose = max(0, debug - 2))

  amat <- 1 * as(gies$essgraph, "matrix")
  rownames(amat) <- colnames(amat) <- nodes

  if (dag)
    amat <- phsl:::pdag2dag_cpp(amat, nodes = nodes)$graph

  return(amat)
}



# Convert interventions to targets and target.index as in pcalg

int2targets <- function(interventions,
                        nodes){

  uint <- unique(interventions)

  targets <- lapply(uint, function(x){
    int <- match(x, nodes)
    if (is.na(int))
      int <- integer(0)
    return(int)
  })
  target.index <- match(interventions, uint)

  return(list(targets = targets,
              target.index = target.index))
}

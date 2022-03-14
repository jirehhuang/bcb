######################################################################
## Functions for efficiently executing gies
######################################################################



# Estimate adjacency matrix using gies

estimate_gies <- function(rounds,
                          # blmat,
                          settings,
                          interventions,
                          dag = FALSE,
                          debug = 0){

  ## if minimal, return empty graph
  nodes <- settings$nodes
  if (settings$minimal &&
      settings$method != "bcb-gies" &&
      settings$restrict != "gies"){

    return(bnlearn::amat(bnlearn::empty.graph(nodes = nodes)))
  }
  debug_cli(debug >= 2, cli::cli_alert_info,
            "estimating graph with gies")

  # if (is.null(blmat)){
  #
  #   blmat <- 1 - diag(settings$nnodes)
  #
  # } else if (is.null(dim(blmat))){
  #
  #   blmat <- row2mat(row = blmat, nodes = nodes)
  # }
  data <- rounds$data[seq_len(length(interventions)),]
  score <- new("bnlearn_score", data = data, interventions = interventions,
               nodes = nodes, score = settings$score, extra.args =
                 bnlearn:::check.score.args(score = settings$score,
                                            network = settings$bn.fit,
                                            data = data, extra.args = list()))

  gies <- pcalg::gies(score = score, maxDegree = settings$max_parents,
                      # fixedGaps = blmat, adaptive = ifelse(any(blmat == 1),
                      #                                      "vstructures", "none"),
                      iterate = TRUE, phase = c("forward", "backward", "turning"),
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

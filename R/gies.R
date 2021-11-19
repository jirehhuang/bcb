######################################################################
## Functions for efficiently executing gies
######################################################################



# Estimate adjacency matrix using gies

estimate_gies <- function(ps,
                          settings,
                          interventions,
                          dag = TRUE,
                          debug = 0){

  debug_cli(debug >= 2, cli::cli_alert_info,
            "estimating graph with gies")

  nodes <- settings$nodes
  score <- new("lookup_score", ps = ps,
               interventions = interventions,
               nodes = nodes)

  ## TODO: fixedGaps and adaptive with blmat

  gies <- pcalg::gies(score = score, maxDegree = settings$max_parents,
                      phase = c("forward", "backward", "turning"),
                      iterate = TRUE, verbose = max(0, debug - 2))

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

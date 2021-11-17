# Estimate adjacency matrix using gies

estimate_gies <- function(cache,
                          settings,
                          dag = TRUE,
                          debug = FALSE){

  nodes <- names(cache)

  score <- new("lookup_score", cache = cache,
               nodes = nodes)

  ## TODO: fixedGaps and adaptive with blmat

  gies <- pcalg::gies(score = score, maxDegree = settings$max_parents,
                      phase = c("forward", "backward", "turning"),
                      iterate = TRUE, verbose = debug)

  amat <- 1 * as(gies$essgraph, "matrix")
  rownames(amat) <- colnames(amat) <- nodes

  if (dag)
    amat <- phsl:::pdag2dag_cpp(amat, nodes = nodes)$graph

  return(amat)
}

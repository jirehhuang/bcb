debug <- 3



######################################################################
## Gaussian
######################################################################

bcb:::load_example(eg = "gnet", network = "asia")



## estimate_gies(), which uses lookup_score class and lookup()
dag <- bcb:::estimate_gies(ps = ps, settings = settings,
                           interventions = interventions,
                           dag = FALSE, debug = debug)
testthat::expect_true(all(dag == bnlearn::amat(bn.fit)))

## gies using bnlearn_score class
targets <- lapply(unique(interventions), function(x){
  int <- match(x, settings$nodes)
  if (is.na(int))
    int <- integer(0)
  return(int)
})
target.index <- match(interventions,
                      unique(interventions))
score <- new("bnlearn_score", data = data,
             interventions = interventions,
             # targets = targets,
             # target.index = target.index,
             nodes = settings$nodes, score = "bge")
gies <- pcalg::gies(score = score, maxDegree = settings$max_parents,
                    phase = c("forward", "backward", "turning"),
                    iterate = TRUE, verbose = debug)
dag <- as(gies$essgraph, "matrix")
testthat::expect_true(all(dag == bnlearn::amat(bn.fit)))

debug <- 3



######################################################################
## Gaussian
######################################################################

bcb:::load_example(eg = "gnet", network = "asia")



## compute_scores()
scores <- bcb:::compute_scores(data = data, settings = settings,
                               interventions = interventions,
                               output = TRUE, debug = debug)

## compute_ps()
ps <- bcb:::compute_ps(settings = settings,
                       interventions = interventions, debug = debug)

## compute_arp()
arp <- bcb:::compute_arp(settings = settings,
                         interventions = interventions, debug = debug)

## ps2es()
bma <- bcb:::ps2es(ps = ps, settings = settings)

## es2mpg()
mpg <- bcb:::es2mpg(es = bma)

## convert_ps()
testthat::expect_identical(
  ps,
  ## list -> data.frame -> list
  bcb:::convert_ps(bcb:::convert_ps(ps, "data.frame"), "list")
)

## lookup(), lookup_scores(), lookup_scores_cpp()
lu <- all(unlist(lapply(settings$nodes, function(node){

  apply(ps[[node]], 1, function(row){

    parents <- row[seq_len(settings$max_parents)]

    row[settings$max_parents + 1] ==
      c(
        ps[[node]][bcb:::lookup(parents = parents[!is.na(parents)],
                             ps_i = ps[[node]]), settings$max_parents + 1],
        bcb:::lookup_score(target = node,
                           parents = parents,
                           ps = ps),
        bcb:::lookup_score_cpp(parents = parents[!is.na(parents)],
                               ps_i = ps[[node]]),
        bcb:::lookup_score_cpp(parents = parents[!is.na(parents)],
                               ps_i = ps[[node]],
                               score_col = match("score",
                                                 names(ps[[node]]))-1)
      )
  })
})))
testthat::expect_true(lu)

## execute_mds()
graph <- bcb:::execute_mds(ps = ps, settings = settings,
                           seed = 1, debug = debug)
testthat::expect_identical(dim(graph), dim(bnlearn::amat(bn.fit)))



######################################################################
## General
######################################################################



## test bida github examples
bcb:::test_bida(debug = debug)

## test mds github examples
bcb:::test_mds(debug = debug)


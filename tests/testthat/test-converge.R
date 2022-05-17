debug <- 3
set.seed(1)



######################################################################
## Gaussian
######################################################################

bcb:::load_example(eg = "gnet", network = "asia")
true <- bnlearn::amat(bn.fit)

bcb:::eval_graph(est = true, true = true, cp = FALSE)

## large sample observational data
data

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

bcb:::eval_graph(est = mpg, true = true, cp = FALSE)

bn.fit <- bcb::load_bn.fit(x = "survey")
bcb:::load_example(network = "asia")


test <- bcb:::remove_arcs(bn.fit, arcs = bnlearn::arcs(bn.fit))




## test
bn.fit <- bcb::load_bn.fit(x = "sachs", reorder = TRUE, rename = TRUE)
bnlearn::graphviz.plot(bn.fit)

test <- bcb:::add_j_m_pt(bn.fit = bn.fit, nodes = c("V1", "V11"), ignore_if_present = TRUE)
lapply(test, `[[`, "mpt")
lapply(test, `[[`, "jpt")

bn.fit2 <- bcb:::solve_cpt_dnet(bn.fit = bn.fit, target = "V1", parents = c("V7", "V3"), ordered = TRUE)


bn <- bnlearn::empty.graph(names(bn.fit))
bnlearn::amat(bn) <- bnlearn::amat(bn.fit)

dnet <- bcb:::bn2dnet(bn = bn, seed = 1, min_marginal = 0.1)


source('~/Documents/ucla/research/projects/current/code/generate_data.grid.R')
bn.fit2 <- solve_new.cpt(child = "V1", parents = c("V7", "V3"),
                         bn.fit = bn.fit, ordered = TRUE, check.acyclic = FALSE)



## test process_dnet()
bn.fit <- bcb::load_bn.fit(x = "munin1")
bnlearn::graphviz.plot(bn.fit)

dnet <- bcb:::process_dnet(bn.fit = bn.fit,
                           min_levels = 2, max_levels = 2,
                           max_in_deg = 2,
                           max_out_deg = 10,
                           debug = 2)

bcb:::remove_arcs(bn.fit = bn.fit, arcs = cbind("V175", "V186"), debug = 2)


length(dnet)
sapply(dnet, function(x) dim(x$prob)[1])
sapply(dnet, function(x) dimnames(x$prob)[[1]])

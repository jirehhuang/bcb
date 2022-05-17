debug <- 3



## test bn_list2bn.fit
bn.fit0 <- phsl::bnrepository(x = "asia")
bn_list <- lapply(bn.fit0, function(x) x)
bn.fit <- bcb:::bn_list2bn.fit(bn_list = bn_list)

testthat::expect_true(all.equal(bn.fit, bn.fit0))

## test parallel graph in load_bn.fit
bn.fit <- load_bn.fit(x = "parallel_10")

testthat::expect_true(length(bnlearn::nodes(bn.fit)) == 10)

## test chain graph in load_bn.fit
bn.fit <- load_bn.fit(x = "chain_10")

testthat::expect_true(length(bnlearn::nodes(bn.fit)) == 10)

## test random graph in load_bn.fit
bn.fit <- load_bn.fit(x = "random_10_2_1")

testthat::expect_true(length(bnlearn::nodes(bn.fit)) == 10)

##
# bn.fit <- load_bn.fit(x = "barley", reorder = FALSE, rename = TRUE)
mpt <- bcb:::compute_mpt(bn.fit = bn.fit, ordered = FALSE, debug = debug)


bn.fit <- bcb:::restrict_bn.fit(bn.fit = bn.fit,
                                max_in_deg = 2, max_out_deg = 2,
                                max_levels = 3, min_levels = 2,
                                reorder = FALSE, debug = debug)

## TODO: test jpt

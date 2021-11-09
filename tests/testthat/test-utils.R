## test bn_list2bn.fit

bn.fit0 <- phsl::bnrepository(x = "asia")
bn_list <- lapply(bn.fit0, function(x) x)
bn.fit <- bn_list2bn.fit(bn_list = bn_list)

testthat::expect_true(all.equal(bn.fit, bn.fit0))

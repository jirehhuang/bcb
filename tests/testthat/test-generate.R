## test generate_gnet()

bn.fit <- load_bn.fit(x = "asia")
gnet <- bn2gnet(bn = bn.fit, seed = 1,
                coefs = c(0.5, 1),
                vars = c(0.5, 1),
                normalize = TRUE)



## test ribn() for observational data, and normalize

data <- ribn(x = gnet, n = 1e6)
testthat::expect_true(all(abs(sapply(data, sd) - 1) < 1e-2))



## test ribn() for gnet

intervene <- list(list(V1 = 1, n = 10),
                  list(V2 = 1, V4 = -1, n = 20),
                  list(n = 50))
data <- ribn(x = gnet, n = 100, intervene = intervene, seed = 1, debug = FALSE)

testthat::expect_true(nrow(data) == 100)
testthat::expect_true(sum(data$V1 == 1) == 10)
testthat::expect_true(sum(data$V2 == 1 & data$V4 == -1) == 20)



## test ribn() for dnet

dnet <- load_bn.fit(x = "asia", reorder = TRUE, rename = TRUE)

intervene <- list(list(V1 = "0", n = 10),
                  list(V2 = "0", V4 = 1, n = 20),
                  list(V3 = "dirichlet", n = 50))
data <- ribn(x = dnet, n = 100, intervene = intervene, seed = 1, debug = FALSE)

testthat::expect_true(nrow(data) == 100)
testthat::expect_true(sum(data$V1 == "0") >= 10)
testthat::expect_true(sum(data$V2 == "0" & data$V4 == "0") >= 20)

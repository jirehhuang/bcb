debug <- 3



## bn2gnet()
bn.fit <- load_bn.fit(x = "asia")
gnet <- bcb:::bn2gnet(bn = bn.fit, seed = 1,
                      coefs = c(0.5, 1),
                      vars = c(0.5, 1),
                      normalize = TRUE)

## ribn() for observational data, and normalize
data <- ribn(x = gnet, n = 1e6)
testthat::expect_true(all(abs(sapply(data, sd) - 1) < 1e-2))

## ribn() for gnet
intervene <- list(list(V1 = 1, n = 10),
                  list(V2 = 1, V4 = -1, n = 20),
                  list(n = 50))
data <- ribn(x = gnet, n = 100, intervene = intervene, seed = 1, debug = debug)

testthat::expect_true(nrow(data) == 100)
testthat::expect_true(sum(data$V1 == 1) == 10)
testthat::expect_true(sum(data$V2 == 1 & data$V4 == -1) == 20)

##  ribn() for dnet
dnet <- load_bn.fit(x = "asia", reorder = TRUE, rename = TRUE)

intervene <- list(list(V1 = "0", n = 10),
                  list(V2 = "0", V4 = 1, n = 20),
                  list(V3 = "dirichlet", n = 50))
data <- ribn(x = dnet, n = 100, intervene = intervene, seed = 1, debug = debug)

testthat::expect_true(nrow(data) == 100)
testthat::expect_true(sum(data$V1 == "0") >= 10)
testthat::expect_true(sum(data$V2 == "0" & data$V4 == "0") >= 20)

## gen_data_grid()
networks <- c("asia")

for (network in networks){

  data_grid <- build_data_grid(network = network,
                               data_type = "gaussian",
                               n_obs = 20,
                               reg_lb = 0.2,
                               copies = 2)

  path <- file.path(gsub("/tests.*", "", getwd()),
                    "tests", "temp", sprintf("test-%s", network))

  gen_data_grid(data_grid = data_grid,
                path = path,
                n_dat = 2,
                seed0 = 0,
                regenerate = TRUE,
                recache = TRUE,
                n_cores = 1,
                debug = max(1, debug - 2))
}

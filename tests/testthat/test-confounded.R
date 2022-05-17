## example with confounded graph

debug <- 0

nds <- c("Z", "X", "Y")
bn <- bnlearn::empty.graph(nodes = nds)
bnlearn::arcs(bn) <- matrix(c("Z", "X",
                              "Z", "Y",
                              "X", "Y"),
                            ncol = 2, byrow = TRUE)
bnlearn::graphviz.plot(bn)
# bn <- bcb:::rename_bn.fit(bn)

gnet <- bcb:::bn2gnet(bn = bn,
                      seed = 5,
                      coefs = c(0.5, 1),
                      vars = c(0.5, 1),
                      normalize = TRUE,
                      intercept = FALSE)
bn.fit <- gnet




settings0 <- bcb:::check_settings(settings = list(method = "bcb-ts",
                                                  n_obs = 0,
                                                  n_int = 0,
                                                  eta = 0.1,
                                                  max_cache = 1),
                                  bn.fit = bn.fit,
                                  debug = debug)
settings <- settings0

settings2rounds <- function(settings){

  rounds <- bcb::bandit(bn.fit = bn.fit,
                        settings = settings,
                        seed0 = 1,
                        debug = debug)
  return(rounds)
}
rounds2table <- function(rounds){

  tbl <- cbind(rounds$ps$X[,c(1, 2, 4)],
               psi_bda = rounds$bda$X$Y$beta_bda,
               psi_est = rounds$bda$X$Y$beta_est,
               psi_est_avg = rounds$arms$estimate[4],
               psi_true = rounds$mu_true[4])
  return(round(tbl, digits = 2))
}



settings$n_obs <- 10
settings$n_int <- 0
print(rounds2table(rounds <- settings2rounds(settings)))
bcb:::ps2es(ps = rounds$ps, settings = settings)


settings$n_obs <- 50
settings$n_int <- 0
print(rounds2table(rounds <- settings2rounds(settings)))
bcb:::ps2es(ps = rounds$ps, settings = settings)


settings$n_obs <- 10
settings$n_int <- 40
settings$eta <- 0
print(rounds2table(rounds <- settings2rounds(settings)))
bcb:::ps2es(ps = rounds$ps, settings = settings)
rounds$arms


settings$n_obs <- 10
settings$n_int <- 40
settings$eta <- 0.1
print(rounds2table(rounds <- settings2rounds(settings)))
bcb:::ps2es(ps = rounds$ps, settings = settings)
rounds$arms


rounds$bda$X$Y

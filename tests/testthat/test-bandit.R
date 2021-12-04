debug <- 3
set.seed(1)



######################################################################
## Gaussian
######################################################################

bcb:::load_example(eg = "gnet", network = "asia")

if (FALSE){

  bnlearn::graphviz.plot(bn.fit)

  set.seed(1)
  data <- bnlearn::rbn(x = bn.fit, n = 1e6)
  data4 <- ribn(x = bn.fit, intervene = list(list(V4 = mean(data$V4) + 1, n = 1e6)))
  mean(data4$V8)
  sd(data4$V8)

  # data$V3 <- data$V3 * 10

  bn.fit$V4
  bn.fit$V8

  effects_true <- bcb:::bn.fit2effects(bn.fit = bn.fit)
  effects_true["V4", "V8", 1]
  mean(data4$V8) - mean(data$V8)
  effects_true[,,1] + effects_true[,,2]

  m <- lm(V8 ~ V4 + V3, data = data)
  mean(data$V8 - m$coefficients[2] * data$V4)
  effects_true["V4", "V8", 2]

  sd(data$V8 - m$coefficients[1] - m$coefficients[2] * data$V4)
  sd(data$V8 - m$coefficients[2] * data$V4)
  sd(data4$V8)

  m <- lm(V8 ~ V4 + V3, data = data)
  phi <- matrix(c(m$coefficients[1], 1, mean(data$V3)), ncol = 1)
  y <- as.matrix(data$V8, ncol = 1)
  X <- as.matrix(cbind(V0 = rep(1, nrow(data)), data[,c("V4", "V3")]))

  var(m$residuals) * t(phi) %*% solve(t(X) %*% X) %*% phi

  1 / t(phi) %*% solve(t(X) %*% X) %*% phi
}

## bandit()
settings$n_obs <- 10
settings$n_int <- 10
bndt <- bandit(bn.fit = bn.fit, settings = settings, debug = debug)

## write_rounds()
rounds_dir <- file.path(settings$temp_dir,
                        sprintf("rounds%g", settings$run))
bcb:::write_rounds(rounds = bndt,
                   where = rounds_dir)

## read_rounds()
bndt_dir <- bcb:::read_rounds(where = rounds_dir)

testthat::expect_identical(names(bndt), names(bndt_dir))

nms <- names(bndt)
nms <- nms[sapply(nms, function(nm) ! identical(bndt[[nm]], bndt_dir[[nm]]))]
nms <- setdiff(nms, c("settings"))

for (nm in nms){

  ae <- all.equal(bndt[[nm]], bndt_dir[[nm]],
                  check.attributes = FALSE)

  testthat::expect_true(ae)
}

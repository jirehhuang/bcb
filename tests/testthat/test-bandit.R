debug <- 3
set.seed(1)



######################################################################
## Gaussian
######################################################################

bcb:::load_example(eg = "gnet", network = "asia")

## bandit()
settings$method <- "bcb-star"
settings$n_obs <- 10
settings$n_int <- 10
settings$restrict <- "pc"
bndt <- bandit(bn.fit = bn.fit, settings = settings, debug = debug)

## write_rounds()
settings <- bndt$settings
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

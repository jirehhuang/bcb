######################################################################
## Analyze compiled dkpar results
######################################################################


## Append _done from the generated directories from test_bda.R,
## resulting in folder names test_bda-d_done and test_bda-g_done


library(bcb)
library(ggplot2)
path0 <- ifelse(get_projects_dir(envir = environment()) == getwd(),
                projects_dir, file.path(projects_dir, "current",
                                        "simulations", "test_bda"))
wh <- c(1.3, 0.6) * 139.68


######################################################################
## Discrete
######################################################################

path <- file.path(path0, sprintf("test_bda-d"))
df <- readRDS(file.path(path, "test_Var_Pr_1000_1000.rds"))

df$Method <- ifelse(df$method == "naive", "Naive",
                    ifelse(grepl("sampling", df$method), "Sampling",
                           ifelse(grepl("bootstrap", df$method),
                                  "Bootstrap", "Proposed")))
df$Method <- factor(df$Method, levels = c("Naive", "Bootstrap",
                                          "Proposed", "Sampling"))

## Coverage grid
coverage_d <- ggplot(data = df, aes(x = Method, y = coverage)) +
  geom_hline(yintercept = 0.95, linetype = "dashed") +
  geom_boxplot(outlier.shape = NA, width = 0.5, aes(fill = Method)) +
  coord_cartesian(ylim = c(0.85, 1)) +
  theme_fixed(bool_axis_text = FALSE,
              bool_legend_title = FALSE,
              legend.position = "bottom") +
  facet_wrap(~ n + r, ncol = 4) +
  labs(x = "Estimator", y = "Coverage per Scenario")

sapply(c("eps", "png"), function(x){

  ggsave(filename = sprintf("%s/bda_coverage-d.%s", path, x),
         plot = coverage_d, device = x, dpi = 1600,
         width = wh[1], height = wh[2] * 2.5, units = "mm")
})

## Discrete execution time
time_d <- sum(df$time[df$method == "sampling"])


######################################################################
## Gaussian
######################################################################

path <- file.path(path0, sprintf("test_bda-g"))
df <- readRDS(file.path(path, "test_Var_beta_1e+05.rds"))

df$n <- factor(df$n, levels = sort(unique(df$n)))
df$Data <- factor(ifelse(df$method == "mix", "Ensemble", "Observational"),
                  levels = c("Observational", "Ensemble"))

## Coverage
coverage_g <- ggplot(data = df, aes(x = n, y = coverage)) +
  geom_hline(yintercept = 0.95, linetype = "dashed") +
  geom_boxplot(aes(fill = Data)) +
  coord_cartesian(ylim = c(0.953, 0.947)) +
  theme_fixed(bool_axis_text = TRUE,
              bool_legend_title = FALSE,
              legend.position = "bottom") +
  labs(x = "Sample Size", y = "Coverage per Scenario")

sapply(c("eps", "png"), function(x){

  ggsave(filename = sprintf("%s/bda_coverage-g.%s", path, x),
         plot = coverage_g, device = x, dpi = 1600,
         width = wh[1], height = wh[2], units = "mm")
})

## Gaussian execution time
time_g <- sum(df$time) / 2


## Total execution time (conservative)
time_dg <- prettyunits::pretty_sec(time_d + time_g)

######################################################################
## Analyze compiled dkpar results
######################################################################


## Append _done from the compiled directories from compile_dkpar.R,
## resulting in folder names dkpar_3_6_3-d_done and dkpar_3_6_3-g_done


library(bcb)
library(ggplot2)
library(ggpubr)
library(dplyr)
path0 <- ifelse(get_projects_dir(envir = environment(), debug = 0) == getwd(),
                projects_dir, file.path(projects_dir, "current",
                                        "simulations", "mcmc"))
network <- "dkpar_3_6_3"
wh <- c(1.3, 0.6) * 139.68


######################################################################
## Discrete
######################################################################

path <- file.path(path0, sprintf("%s-d_done", network))
df0 <- readRDS(file.path(path, "concise", "df_0.rds")) %>%
  filter(t %in% seq(500, 5000, 1500))
df2 <- readRDS(file.path(path, "concise", "df_2.rds"))
method_grid <- read.table(file.path(path, "concise", "method_grid.txt"))

unique(df0$method)
methods <- c(sprintf("bcb-mcmc-bucb%s", c(51, 151)),
             sprintf("bcb-mcmc-ts%s", c(51, 151)),
             sprintf("bcb-mcmc-ucb%s", c(54, 154)),
             sprintf("bcb-bucb%s", c(51)),
             sprintf("bcb-ts%s", c(51)),
             sprintf("bcb-ucb%s", c(54)))

df <- df0 %>%
  filter(method %in% methods) %>%
  mutate(t = as.factor(t),
         Alg = case_when(grepl("bucb", method) ~ "Bayes-UCB",
                         grepl("ts", method) ~ "TS",
                         TRUE ~ "UCB")) %>%
  mutate(Alg = sprintf("BBB-%s", Alg),
         type = case_when(grepl("15", method) ~ "Hybrid MCMC",
                          grepl("mcmc", method) ~ "MCMC",
                          TRUE ~ "Exact")) %>%
  mutate(type = factor(type, levels = c("Exact", "MCMC", "Hybrid MCMC")))

colors <- scales::hue_pal()(8)[c(8, 6, 4)]
ltys <- scales::linetype_pal()(8)
ltys <- c(ltys[1], ltys)[c(8, 6, 4)]

df %>%
  group_by(Alg, type, t) %>%
  summarize(Q1 = quantile(expected_cumulative, probs = c(0.25)),
            Q3 = quantile(expected_cumulative, probs = c(0.75)),
            median = median(expected_cumulative),
            gt60 = sum(expected_cumulative > 60)) %>%
  mutate(IQR = Q3 - Q1) %>%
  mutate(upper = Q3 + 1.5 * IQR)

## Discrete cumulative regret boxplot
cumulative_bp <- df %>%
  filter(expected_cumulative <= 60) %>%
  ggplot(aes(x = t, y = expected_cumulative,
             fill = type)) +
  geom_boxplot(color = "black",
               outlier.shape = 21,
               outlier.size = 2,
               position = position_dodge(width = 0.9)) +
  geom_text(data = df %>%
              group_by(Alg, type, t) %>%
              summarize(gt = sum(expected_cumulative > 60)),
            aes(x = t, y = 65, group = type, label = gt),
            position = position_dodge(width = 0.9),
            size = 2.5) +
  scale_y_continuous(breaks = c(seq(0, 60, 20), 65),
                     labels = c(seq(0, 60, 20), ">60")) +
  scale_fill_manual(values = colors,
                    breaks = levels(df$type)) +
  coord_cartesian(ylim = c(0, 65)) +
  facet_wrap(~Alg) +
  theme_fixed(bool_axis_text = TRUE) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  labs(x = "Time Step", y = "Cumulative Regret")
cumulative_bp

sapply(c("eps", "png"), function(x){

  ggsave(filename = sprintf("%s/cumulative_mcmc_bp.%s", path, x),
         plot = cumulative_bp, device = x, dpi = 1600,
         width = wh[1], height = wh[2] * 0.9, units = "mm")
})


######################################################################
## Computation time
######################################################################

## Total (conservative)
sum(df0$time * 200) %>%
  prettyunits::pretty_sec()

## Per iteration for different methods
df %>%
  mutate(alg_group = case_when(grepl("mcmc", method) &
                                 grepl("15", method) ~ "Hybrid MCMC",
                               grepl("mcmc", method) ~ "MCMC",
                               TRUE ~ "Exact")) %>%
  group_by(alg_group) %>%
  summarize(mean = mean(time))

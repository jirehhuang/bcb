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
                                        "simulations", "mcmc", "analyze"))
network <- "child"
wh <- c(1.3, 0.6) * 139.68


######################################################################
## Discrete
######################################################################

path <- file.path(path0, sprintf("%s-d_done", network))
df0 <- readRDS(file.path(path, "concise", "df_2.rds"))
method_grid <- read.table(file.path(path, "concise", "method_grid.txt"))


## parameter tuning
# methods <- unique(df0$method)
# methods <- methods[grepl("ucb", methods) & !grepl("cn", methods) & !grepl("bucb", methods)]
# methods <- methods[grepl("cn-ucb", methods)]
# methods <- sprintf("ucb%g", seq_len(10))
methods <- sprintf("cn-ucb%g", seq_len(5))
df <- df0[df0$method %in% methods,]
ggplot(data = df,
       aes(x = t, y = expected_cumulative, group = method,
           color = method, lty = method)) +
  geom_line(size = 1) +
  theme_fixed(bool_axis_text = TRUE) +
  labs(x = "Time Step", y = "Cumulative Regret") +
  theme(legend.key.width = unit(2, "cm"))


unique(df0$method)
shift <- switch(network, dkpar_3_6_3 = 3, 0)
methods <- c("bcb-mcmc-bucb151",
             "bcb-mcmc-ts151",
             "bcb-mcmc-ucb154",
             "bucb1", "ts1", "ucb1",
             "cn-bucb1", "cn-ts1", "cn-ucb4")

df <- df0[df0$method %in% methods,]
df$Alg <- ifelse(grepl("bucb", df$method), "Bayes-UCB",
                 ifelse(grepl("ts", df$method), "TS", "UCB"))
df$Alg <- sprintf("Alg: %s", df$Alg)
df$Method <- ifelse(grepl("bcb", df$method),
                    sprintf("BBB-Alg(%g)",
                            method_grid$n_obs[as.numeric(gsub("[^\\d]+", "",
                                                              df$method, perl=TRUE))]),
                    ifelse(grepl("cn", df$method), "Alg*", "Alg"))

labels <- if (is.factor(df$Method)) levels(df$Method) else unique(df$Method)
# labels <- c(labels[7:8], labels[-(7:8)])
labels <- c("Alg", "Alg*", sprintf("BBB-Alg(%g)", 100 * 2^seq(0, 5)))
df$Method <- factor(df$Method, levels = labels)
colors <- scales::hue_pal()(length(labels))
ltys <- scales::linetype_pal()(length(labels) - 1)
ltys <- c(ltys[1], ltys)

## For combined plots
df$k <- as.character(log2(method_grid$n_obs[as.numeric(gsub("[^\\d]+", "",
                                                            df$method,
                                                            perl=TRUE))] / 1e2))
df$k <- ifelse(grepl("bcb-", df$method), df$k, as.character(df$Method))
df$Method_k <- sapply(seq_len(nrow(df)), function(x)
  gsub("\\(.*\\)", sprintf("(k=%s)", df$k[x]), df$Method[x]))
# labels2 <- unique(df$Method_k)[c(7:8, 1:6)]
labels2 <- c("Alg", "Alg*", sprintf("BBB-Alg(k=%g)", seq(0, 5)))
# labels2 <- c("Alg", "Alg*", "BBB-Alg")
# df$Method_k <- ifelse(grepl("BBB", df$Method_k), "BBB-Alg", df$Method_k)

## Discrete cumulative regret
cumulative_d <- ggplot(data = df,
                       aes(x = t, y = expected_cumulative, group = Method_k,
                           color = Method_k, lty = Method_k)) +
  geom_line(size = 1) +
  theme_fixed(bool_axis_text = TRUE,
              legend.position = "bottom") +
  theme(legend.margin = margin(c(-5, 0, 5, 0))) +
  labs(x = element_blank(), y = "Average Cumulative") +
  scale_fill_manual(values = colors, breaks = labels2) +
  scale_color_manual(values = colors, breaks = labels2) +
  scale_linetype_manual(values = ltys, breaks = labels2) +
  facet_wrap(~Alg) +
  theme(legend.key.width = unit(2, "cm"))
cumulative_d

## Discrete cumulative regret with error bars
Method_k <- c("Alg", "Alg*", sprintf("BBB-Alg(k=%g)", c(5)))
width0 <- 0.7
dodge_width0 <- 0.9
shapes <- c(25, 24, 23)
cumulative_err_d <-
  ggplot(data = df[df$t %in% c(seq(500, 5000, 1500)) &
                     df$Method_k %in% Method_k,],
         aes(x = as.factor(t), y = expected_cumulative_500,
             ymin = expected_cumulative_025, ymax = expected_cumulative_975,
             group = Method_k, color = Method_k,
             fill = Method_k, shape = Method_k)) +
  geom_line(aes(lty = Method_k),
            size = 0.5,
            position = position_dodge(width = dodge_width0)) +
  geom_errorbar(size = 1, width = width0,
                position = position_dodge(width = dodge_width0)) +
  geom_point(size = 1.5, stroke = 1, fill = "white",
             position = position_dodge(width = dodge_width0)) +
  theme_fixed(bool_axis_text = TRUE) +
  # theme(legend.margin = margin(c(0, 0, 0, 0))) +
  labs(x = "Time Step", y = "Median Cumulative (95%)") +
  scale_fill_manual(values = colors, breaks = labels2) +
  scale_color_manual(values = colors, breaks = labels2) +
  scale_linetype_manual(values = ltys, breaks = labels2) +
  scale_shape_manual(values = shapes, breaks = Method_k) +
  facet_wrap(~Alg) +
  theme(legend.key.width = unit(1, "cm"))
cumulative_err_d

## Discrete simple regret
simple_d <- ggplot(data = df,
                   aes(x = t, y = greedy_regret, group = Method_k,
                       color = Method_k, lty = Method_k)) +
  geom_line(size = 1) +
  theme_fixed(bool_axis_text = TRUE) +
  labs(x = "Time Step", y = "Discrete Simple") +
  scale_fill_manual(values = colors, breaks = labels2) +
  scale_color_manual(values = colors, breaks = labels2) +
  scale_linetype_manual(values = ltys, breaks = labels2) +
  facet_wrap(~Alg) +
  # coord_cartesian(ylim = c(0, 0.05)) +
  theme(legend.key.width = unit(2, "cm"))
simple_d

## Discrete ESSAE
essae_d <- ggplot(data = df[grepl("bcb", df$method),],
                  aes(x = t, y = sae_bma, group = Method_k,
                      color = Method_k, lty = Method_k)) +
  geom_line(size = 1) +
  theme_fixed(bool_axis_text = TRUE) +
  labs(x = "Time Step", y = "Discrete ESSAE") +
  scale_fill_manual(values = colors, breaks = labels2) +
  scale_color_manual(values = colors, breaks = labels2) +
  scale_linetype_manual(values = ltys, breaks = labels2) +
  facet_wrap(~Alg) +
  # scale_y_continuous(breaks = seq(0, 18, 3)) +
  # coord_cartesian(xlim = c(120, 5000-120), ylim = c(0.6, 15-0.6)) +
  theme(legend.key.width = unit(2, "cm"))
essae_d

## Discrete ESSAE for optimal arm
essae_arm0_d <- ggplot(data = df[grepl("bcb", df$method),],
                       aes(x = t, y = sae_arm0, group = Method_k,
                           color = Method_k, lty = Method_k)) +
  geom_line(size = 1) +
  theme_fixed(bool_axis_text = TRUE) +
  labs(x = "Time Step", y = "Discrete ESSAE") +
  scale_fill_manual(values = colors, breaks = labels2) +
  scale_color_manual(values = colors, breaks = labels2) +
  scale_linetype_manual(values = ltys, breaks = labels2) +
  facet_wrap(~Alg) +
  # scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
  # coord_cartesian(xlim = c(120, 5000-120), ylim = c(0, 4)) +
  theme(legend.key.width = unit(2, "cm"))
essae_arm0_d

## Combined cumulative regret
cumulative_dd <- ggpubr::ggarrange(plotlist = list(cumulative_d, cumulative_err_d), ncol = 1,
                                   common.legend = FALSE, legend = "bottom")

sapply(c("eps", "png"), function(x){

  ggsave(filename = sprintf("%s/cumulative_%s.%s", path, network, x),
         plot = cumulative_dd, device = x, dpi = 1600,
         width = wh[1], height = wh[2] * 1.5, units = "mm")
})


######################################################################
## Computation time
######################################################################

## Total (conservative)
sum(df0$time * 20) %>%
  prettyunits::pretty_sec()

## Per iteration for different methods
df %>%
  mutate(alg_group = case_when(grepl("cn-", method) ~ "Alg*",
                               grepl("bcb-", method) ~ "BBB-Alg",
                               TRUE ~ "Alg")) %>%
  group_by(alg_group) %>%
  summarize(mean = mean(time))

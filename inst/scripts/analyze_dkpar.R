######################################################################
## Analyze compiled dkpar results
######################################################################


## Append _done from the compiled directories from compile_dkpar.R,
## resulting in folder names dkpar_3_6_3-d_done and dkpar_3_6_3-g_done


library(bcb)
library(ggplot2)
library(ggpubr)
network <- "dkpar_3_6_3"
path0 <- ifelse(get_projects_dir(envir = environment()) == getwd(),
                projects_dir, file.path(projects_dir, "current",
                                        "simulations", network))
wh <- c(1.3, 0.6) * 139.68


######################################################################
## Discrete
######################################################################

path <- file.path(path0, sprintf("%s-d_done", network))
df0 <- readRDS(file.path(path, "concise", "df_2.rds"))
method_grid <- read.table(file.path(path, "concise", "method_grid.txt"))

methods <- c(sprintf("bcb-bucb%s", c(seq(1, 60, 10))),
             sprintf("bcb-ts%s", c(seq(1, 60, 10))),
             sprintf("bcb-ucb%s", c(seq(4, 60, 10))),
             "bucb1", "ts1", "ucb4",
             "cn-bucb1", "cn-ts1", "cn-ucb4")

df <- df0[df0$method %in% methods,]
df$Alg <- ifelse(grepl("bucb", df$method), "Bayes-UCB",
                 ifelse(grepl("ts", df$method), "TS", "UCB"))
df$Method <- ifelse(grepl("bcb", df$method),
                    sprintf("BBB-Alg(%g)",
                            method_grid$n_obs[as.numeric(gsub("[^\\d]+", "",
                                                              df$method, perl=TRUE))]),
                    ifelse(grepl("cn", df$method), "Alg*", "Alg"))

labels <- if (is.factor(df$Method)) levels(df$Method) else unique(df$Method)
labels <- c(labels[7:8], labels[-(7:8)])
df$Method <- factor(df$Method, levels = labels)
colors <- scales::hue_pal()(length(labels))
ltys <- scales::linetype_pal()(length(labels) - 1)
ltys <- c(ltys[1], ltys)

## Discrete cumulative regret with head start (hs)
colors_3 <- colors[c(1, 2, 5)]
Alg_hs_list <- lapply(c("Bayes-UCB", "TS", "UCB"), function(Alg){

  start <- ifelse(Alg == "UCB", 4, 1)
  df_Alg <- do.call(rbind, lapply(seq(start, 60, 10), function(x){

    methods <- switch(
      Alg,
      `Bayes-UCB` = c("bucb1", "cn-bucb1", sprintf("bcb-bucb%s", x)),
      TS = c("ts1", "cn-ts1", sprintf("bcb-ts%s", x)),
      UCB = c("ucb4", "cn-ucb4", sprintf("bcb-ucb%s", x))
    )
    temp_mg <- method_grid
    temp_mg$n_obs <- temp_mg$n_obs[x]
    temp_df <- df[df$method %in% methods,]
    temp_df <- headstart_df(temp_df,
                            headstart = temp_mg)
    temp_df$n_obs <- temp_mg$n_obs[x]
    temp_df$Method <- gsub("\\s*\\([^\\)]+\\)", "", temp_df$Method)
    temp_df$Alg <- temp_df$Method

    return(temp_df)
  }))
  labels_Alg <- sapply(c("%s", "%s*", "BBB-%s"), sprintf, "Alg")

  p <- ggplot(data = df_Alg,
              aes(x = t, y = expected_cumulative, group = Alg,
                  color = Alg, lty = Alg)) +
    geom_line(size = 1) +
    theme_fixed(bool_axis_text = FALSE) +
    ylab(sprintf("%s", Alg)) +
    scale_fill_manual(values = colors_3, breaks = labels_Alg) +
    scale_color_manual(values = colors_3, breaks = labels_Alg) +
    scale_linetype_manual(values = ltys, breaks = labels_Alg) +
    facet_wrap(~n_obs, scales = "free", nrow = 1) +
    theme(legend.key.width = unit(2, "cm"),
          axis.text.y = element_blank(),
          axis.ticks = element_blank())

  p <- p + if (Alg == "UCB"){
    xlab("Time Step")
  } else{
    xlab(element_blank())
  }
  return(p)
})
hs_grid <- ggpubr::ggarrange(plotlist = Alg_hs_list,
                             nrow = 3, ncol = 1,
                             common.legend = TRUE,
                             legend = "bottom")

sapply(c("eps", "png"), function(x){

  ggsave(filename = sprintf("%s/hs-d.%s", path, x),
         plot = hs_grid, device = x, dpi = 1600,
         width = wh[1], height = wh[2] * 1.2, units = "mm")
})

## For combined plots
df$k <- as.character(log2(method_grid$n_obs[as.numeric(gsub("[^\\d]+", "",
                                                            df$method,
                                                            perl=TRUE))] / 1e2))
df$k <- ifelse(grepl("bcb-", df$method), df$k, as.character(df$Method))
df$Method_k <- sapply(seq_len(nrow(df)), function(x)
  gsub("\\(.*\\)", sprintf("(k=%s)", df$k[x]), df$Method[x]))
labels2 <- unique(df$Method_k)[c(7:8, 1:6)]

## Discrete cumulative regret
cumulative_d <- ggplot(data = df,
                       aes(x = t, y = expected_cumulative, group = Method_k,
                           color = Method_k, lty = Method_k)) +
  geom_line(size = 1) +
  theme_fixed(bool_axis_text = FALSE) +
  labs(x = element_blank(), y = "Discrete Cumulative") +
  scale_fill_manual(values = colors, breaks = labels2) +
  scale_color_manual(values = colors, breaks = labels2) +
  scale_linetype_manual(values = ltys, breaks = labels2) +
  facet_wrap(~Alg) +
  theme(legend.key.width = unit(2, "cm"))

## Discrete cumulative regret with error bars
Method_k <- c("Alg", "Alg*", sprintf("BBB-Alg(k=%g)", c(0, 3, 5)))
width0 <- 0.7
dodge_width0 <- 0.9
shapes <- c(25, 24, 21, 22, 23)
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
  theme_fixed(bool_axis_text = FALSE) +
  labs(x = element_blank(), y = "Discrete Cumulative") +
  scale_fill_manual(values = colors, breaks = labels2) +
  scale_color_manual(values = colors, breaks = labels2) +
  scale_linetype_manual(values = ltys, breaks = labels2) +
  scale_shape_manual(values = shapes, breaks = Method_k) +
  facet_wrap(~Alg) +
  theme(legend.key.width = unit(1, "cm"))

## Discrete simple regret
simple_d <- ggplot(data = df,
                   aes(x = t, y = greedy_regret, group = Method_k,
                       color = Method_k, lty = Method_k)) +
  geom_line(size = 1) +
  theme_fixed(bool_axis_text = FALSE) +
  labs(x = element_blank(), y = "Discrete Simple") +
  scale_fill_manual(values = colors, breaks = labels2) +
  scale_color_manual(values = colors, breaks = labels2) +
  scale_linetype_manual(values = ltys, breaks = labels2) +
  facet_wrap(~Alg) +
  theme(legend.key.width = unit(2, "cm")) +
  coord_cartesian(ylim = c(0, 0.05))

## Discrete ESSAE
essae_d <- ggplot(data = df[grepl("bcb", df$method),],
                  aes(x = t, y = sae_bma, group = Method_k,
                      color = Method_k, lty = Method_k)) +
  geom_line(size = 1) +
  theme_fixed(bool_axis_text = FALSE) +
  labs(x = element_blank(), y = "Discrete ESSAE") +
  scale_fill_manual(values = colors, breaks = labels2) +
  scale_color_manual(values = colors, breaks = labels2) +
  scale_linetype_manual(values = ltys, breaks = labels2) +
  facet_wrap(~Alg) +
  theme(legend.key.width = unit(2, "cm")) +
  scale_y_continuous(breaks = seq(0, 18, 3)) +
  coord_cartesian(xlim = c(120, 5000-120), ylim = c(0.6, 15-0.6))

## Discrete ESSAE for optimal arm
essae_arm0_d <- ggplot(data = df[grepl("bcb", df$method),],
                       aes(x = t, y = sae_arm0, group = Method_k,
                           color = Method_k, lty = Method_k)) +
  geom_line(size = 1) +
  theme_fixed(bool_axis_text = FALSE) +
  labs(x = element_blank(), y = "Discrete ESSAE") +
  scale_fill_manual(values = colors, breaks = labels2) +
  scale_color_manual(values = colors, breaks = labels2) +
  scale_linetype_manual(values = ltys, breaks = labels2) +
  facet_wrap(~Alg) +
  theme(legend.key.width = unit(2, "cm")) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
  coord_cartesian(xlim = c(120, 5000-120), ylim = c(0, 4))

## Execution time (conservative)
time_d <- sum(df0$time * 100 *
                ifelse(grepl("bcb", df0$method), 5, 10))


######################################################################
## Gaussian
######################################################################

path <- file.path(path0, sprintf("%s-g_done", network))
df0 <- readRDS(file.path(path, "concise", "df_2.rds"))
method_grid <- read.table(file.path(path, "concise", "method_grid.txt"))

methods <- c(sprintf("bcb-bucb%s", c(seq(1, 60, 10))),
             sprintf("bcb-ts%s", c(seq(1, 60, 10))),
             sprintf("bcb-ucb%s", c(seq(6, 60, 10))),
             "bucb1", "ts1", "ucb5",
             "cn-bucb1", "cn-ts1", "cn-ucb6")

df <- df0[df0$method %in% methods,]
df$Alg <- ifelse(grepl("bucb", df$method), "Bayes-UCB",
                 ifelse(grepl("ts", df$method), "TS", "UCB"))
df$Method <- ifelse(grepl("bcb", df$method),
                    sprintf("BBB-Alg(%g)",
                            method_grid$n_obs[as.numeric(gsub("[^\\d]+", "",
                                                              df$method, perl=TRUE))]),
                    ifelse(grepl("cn", df$method), "Alg*", "Alg"))

labels <- if (is.factor(df$Method)) levels(df$Method) else unique(df$Method)
labels <- c(labels[7:8], labels[-(7:8)])
df$Method <- factor(df$Method, levels = labels)
colors <- scales::hue_pal()(length(labels))
ltys <- scales::linetype_pal()(length(labels) - 1)
ltys <- c(ltys[1], ltys)

## For combined plots
df$k <- as.character(log2(method_grid$n_obs[as.numeric(gsub("[^\\d]+", "",
                                                            df$method,
                                                            perl=TRUE))] / 1e1))
df$Method_k <- sapply(seq_len(nrow(df)), function(x)
  gsub("\\(.*\\)", sprintf("(k=%s)", df$k[x]), df$Method[x]))
labels2 <- unique(df$Method_k)[c(7:8, 1:6)]

## Gaussian cumulative regret
cumulative_g <- ggplot(data = df,
                       aes(x = t, y = expected_cumulative, group = Method_k,
                           color = Method_k, lty = Method_k)) +
  geom_line(size = 1) +
  theme_fixed() +
  labs(x = "Time Step", y = "Gaussian Cumulative") +
  scale_fill_manual(values = colors, breaks = labels2) +
  scale_color_manual(values = colors, breaks = labels2) +
  scale_linetype_manual(values = ltys, breaks = labels2) +
  facet_wrap(~Alg) +
  theme(legend.key.width = unit(2, "cm"))

## Combined cumulative regret
cumulative_dg <- ggpubr::ggarrange(plotlist = list(cumulative_d, cumulative_g), ncol = 1,
                                   common.legend = TRUE, legend = "bottom")

sapply(c("eps", "png"), function(x){

  ggsave(filename = sprintf("%s/cumulative.%s", path, x),
         plot = cumulative_dg, device = x, dpi = 1600,
         width = wh[1], height = wh[2] * 1.5, units = "mm")
})

## Gaussian cumulative regret with error bars
cumulative_err_g <-
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
  theme_fixed() +
  labs(x = "Time Step", y = "Gaussian Cumulative") +
  scale_fill_manual(values = colors, breaks = labels2) +
  scale_color_manual(values = colors, breaks = labels2) +
  scale_linetype_manual(values = ltys, breaks = labels2) +
  scale_shape_manual(values = shapes, breaks = Method_k) +
  facet_wrap(~Alg) +
  theme(legend.key.width = unit(1, "cm"))

## Combined cumulative regret
cumulative_err_dg <- ggpubr::ggarrange(plotlist = list(cumulative_err_d, cumulative_err_g), ncol = 1,
                                   common.legend = TRUE, legend = "bottom")

sapply(c("eps", "png"), function(x){

  ggsave(filename = sprintf("%s/cumulative_err.%s", path0, x),
         plot = cumulative_err_dg, device = x, dpi = 1600,
         width = wh[1], height = wh[2] * 2, units = "mm")
})

## Gaussian simple regret
simple_g <- ggplot(data = df,
                   aes(x = t, y = greedy_regret, group = Method_k,
                       color = Method_k, lty = Method_k)) +
  geom_line(size = 1) +
  theme_fixed() +
  labs(x = "Time Step", y = "Gaussian Simple") +
  scale_fill_manual(values = colors, breaks = labels2) +
  scale_color_manual(values = colors, breaks = labels2) +
  scale_linetype_manual(values = ltys, breaks = labels2) +
  facet_wrap(~Alg) +
  theme(legend.key.width = unit(2, "cm")) +
  coord_cartesian(ylim = c(0, 0.1))

## Combined simple regret
simple_dg <- ggpubr::ggarrange(plotlist = list(simple_d, simple_g), ncol = 1,
                               common.legend = TRUE, legend = "bottom")

## Gaussian ESSAE
essae_g <- ggplot(data = df[grepl("bcb", df$method),],
                  aes(x = t, y = sae_bma, group = Method_k,
                      color = Method_k, lty = Method_k)) +
  geom_line(size = 1) +
  theme_fixed() +
  labs(x = "Time Step", y = "Gaussian ESSAE") +
  scale_fill_manual(values = colors, breaks = labels2) +
  scale_color_manual(values = colors, breaks = labels2) +
  scale_linetype_manual(values = ltys, breaks = labels2) +
  facet_wrap(~Alg) +
  theme(legend.key.width = unit(2, "cm"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_y_continuous(breaks = seq(0, 12, 2)) +
  coord_cartesian(xlim = c(120, 5000-120), ylim = c(0.4, 10-0.4))

## Combined ESSAE
essae_dg <- ggpubr::ggarrange(plotlist = list(essae_d, essae_g), ncol = 1,
                              common.legend = TRUE, legend = "bottom")

sapply(c("eps", "png"), function(x){

  ggsave(filename = sprintf("%s/essae.%s", path, x),
         plot = essae_dg, device = x, dpi = 1600,
         width = wh[1], height = wh[2] * 1.5, units = "mm")
})

## Gaussian ESSAE for optimal arm
essae_arm0_g <- ggplot(data = df[grepl("bcb", df$method),],
                       aes(x = t, y = sae_arm0, group = Method_k,
                           color = Method_k, lty = Method_k)) +
  geom_line(size = 1) +
  theme_fixed() +
  labs(x = "Time Step", y = "Gaussian ESSAE") +
  scale_fill_manual(values = colors, breaks = labels2) +
  scale_color_manual(values = colors, breaks = labels2) +
  scale_linetype_manual(values = ltys, breaks = labels2) +
  facet_wrap(~Alg) +
  theme(legend.key.width = unit(2, "cm"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.01)) +
  coord_cartesian(xlim = c(120, 5000-120), ylim = c(0, 1))

## Combined ESSAE for optimal arm
essae_arm0_g <- ggpubr::ggarrange(plotlist = list(essae_arm0_d, essae_arm0_g), ncol = 1,
                                  common.legend = TRUE, legend = "bottom")

## Gaussian execution time (conservative)
time_g <- sum(df0$time * 100 *
                ifelse(grepl("bcb", df0$method), 5, 10))


######################################################################
## Computation time
######################################################################

## Total (conservative)
prettyunits::pretty_sec(time_d + time_g)

## Per iteration for different methods
df %>%
  mutate(dist_group = "Gaussian") %>%
  full_join(readRDS(file.path(file.path(path0, sprintf("%s-d_done", network)),
                              "concise", "df_2.rds")) %>%
              mutate(dist_group = "Discrete")) %>%
  mutate(alg_group = case_when(grepl("cn-", method) ~ "Alg*",
                               grepl("bcb-", method) ~ "BBB-Alg",
                               TRUE ~ "Alg")) %>%
  group_by(dist_group, alg_group) %>%

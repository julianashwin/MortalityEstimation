setwd("/Users/julianashwin/Documents/GitHub/MortalityEstimation/")
rm(list=ls())

require(ggplot2)
require(ggpubr)
require(tidyverse)


decomp_df <- read.csv("figures/benchmark/siler_i2drift_decomp_ber.csv")

params_df <- read.csv("figures/benchmark/siler_i2drift_preds.csv")

"
Siler params
"
plot_df <- params_df %>% 
  filter(year >1900 & year < 2020 & parameter %in% c("b", "B", "c", "C", "d", "σ")) %>%
  mutate(parameter = factor(parameter, levels =c("b", "c", "d", "B", "C", "σ"), ordered = T))

ggplot(plot_df) + 
  theme_bw() + facet_wrap(~parameter, scales = "free_y") + 
  geom_ribbon(aes(x = year, ymin=pc025, ymax=pc975), alpha = 0.2) +
  geom_ribbon(aes(x = year, ymin=pc15, ymax=pc85), alpha = 0.3) + 
  geom_line(aes(x = year, y = median)) + 
  labs(x = "Year", y = "")
ggsave("figures/benchmark/siler_i2drift_params.pdf", width = 8, height = 4, device = cairo_pdf)


"
Drift terms params
"
params_df$parameter <- str_replace_all(params_df$parameter, "α_", "Drift ")
plot_df <- params_df %>% 
  filter(year >1900 & year < 2020 & parameter %in% c("Drift b", "Drift B", "Drift c", "Drift C", 
                                                     "Drift d", "Drift σ")) %>%
  mutate(parameter = factor(parameter, levels =c("Drift b", "Drift c", "Drift d", "Drift B", 
                                                 "Drift C", "Drift σ"), ordered = T))
ggplot(plot_df) + 
  theme_bw() + facet_wrap(~parameter, scales = "free_y") + 
  geom_ribbon(aes(x = year, ymin=pc025, ymax=pc975), alpha = 0.2) +
  geom_ribbon(aes(x = year, ymin=pc15, ymax=pc85), alpha = 0.3) + 
  geom_line(aes(x = year, y = median)) + 
  labs(x = "Year", y = "")
ggsave("figures/benchmark/siler_i2drift_drifts.pdf", width = 8, height = 4, device = cairo_pdf)


"
Variance params
"
params_df$parameter <- str_replace_all(params_df$parameter, "σ", "sigma")
params_df$parameter <- str_replace_all(params_df$parameter, "α", "alpha")
plot_df <- params_df %>% 
  mutate(parameter = str_replace_all(parameter, "σ", "sigma"),
         parameter = str_replace_all(parameter, "alpha", "alpha")) %>%
  filter(str_detect(parameter, "sigma_")) %>%
  mutate(across(where(is.numeric), round, digits=3)) %>%
  select(parameter, median, pc025, pc15, pc85, pc975)
  
stargazer(as.matrix(plot_df), table.placement = "H", label = "tab:variance_terms",
          title = "Posterior distributions of variance terms")



ggplot(decomp_df) + theme_bw() + 
  
  geom_bar

# Get long versions of the changes in parameters and derivatives
change_df <- decomp_df %>%
  select(year, Δb, ΔB, Δc, ΔC, Δd) %>%
  pivot_longer( cols = -year, names_to = "parameter", values_to = "change") %>%
  mutate(parameter = str_remove_all(parameter,"Δ"))
LEgrad_df <- decomp_df %>%
  select(year, LE_b, LE_B, LE_c, LE_C, LE_d) %>%
  pivot_longer( cols = -year, names_to = "parameter", values_to = "LE_grad") %>%
  mutate(parameter = str_remove_all(parameter,"LE_"))
hgrad_df <- decomp_df %>%
  select(year, h_b, h_B, h_c, h_C, h_d) %>%
  pivot_longer( cols = -year, names_to = "parameter", values_to = "h_grad") %>%
  mutate(parameter = str_remove_all(parameter,"h_"))
Lstargrad_df <- decomp_df %>%
  select(year, Lstar_b, Lstar_B, Lstar_c, Lstar_C, Lstar_d) %>%
  pivot_longer( cols = -year, names_to = "parameter", values_to = "Lstar_grad") %>%
  mutate(parameter = str_remove_all(parameter,"Lstar_"))
# Combine into long dataframe
decomp_long <- decomp_df %>%
  select(year, LE_mod, h_mod, Lstar_mod, ΔLE_mod, Δh_mod, ΔLstar_mod) %>%
  full_join(change_df, by = "year") %>%
  left_join(LEgrad_df, by = c("year", "parameter")) %>%
  left_join(hgrad_df, by = c("year", "parameter")) %>%
  left_join(Lstargrad_df, by = c("year", "parameter"))
# Plot the decomposition of changes over time
LE_plt <- ggplot(decomp_long) + theme_bw() + 
  geom_bar(aes(x = year, y = change*LE_grad, fill = parameter), 
           stat = "identity", position = "stack") + 
  geom_line(aes(x = year, y = ΔLE_mod)) +
  labs(title = "Life Expectancy", x = "Year", y = "ΔLE", fill = "Parameter")
Lstar_plt <- ggplot(decomp_long) + theme_bw() + 
  geom_bar(aes(x = year, y = change*Lstar_grad, fill = parameter), 
           stat = "identity", position = "stack") + 
  geom_line(aes(x = year, y = ΔLstar_mod)) +
  labs(title = "Lifespan", x = "Year", y = "ΔL*", fill = "Parameter")
h_plt <- ggplot(decomp_long) + theme_bw() + 
  geom_bar(aes(x = year, y = change*h_grad, fill = parameter), 
           stat = "identity", position = "stack") + 
  geom_line(aes(x = year, y = Δh_mod)) +
  labs(title = "Lifespan Equality", x = "Year", y = "Δh", fill = "Parameter")
ggarrange(LE_plt, Lstar_plt, h_plt, nrow = 1, common.legend = TRUE, legend = "right")
ggsave("figures/benchmark/siler_i2drift_decomp.pdf", width = 12, height = 4, device = cairo_pdf)


"
LE
"
plot_df <- filter(params_df, parameter %in% c("LE"))
ggplot(plot_df) + theme_bw() +
  geom_ribbon(aes(x = year, ymin=pc025, ymax=pc975), alpha = 0.2) +
  geom_ribbon(aes(x = year, ymin=pc15, ymax=pc85), alpha = 0.3) + 
  geom_line(aes(x = year, y = median)) + 
  labs(x = "Year", y = "BPLE") + 
  geom_vline(aes(xintercept = 2018), linetype = "dashed")
ggsave("figures/benchmark/siler_i2drift_bple_fcast.pdf", width = 5, height = 3, device = cairo_pdf)

"
alpha_C
"
plot_df <- filter(params_df, parameter %in% c("α_C"))
ggplot(plot_df) + theme_bw() +
  geom_ribbon(aes(x = year, ymin=pc025, ymax=pc975), alpha = 0.2) +
  geom_ribbon(aes(x = year, ymin=pc15, ymax=pc85), alpha = 0.3) + 
  geom_line(aes(x = year, y = median)) + 
  labs(x = "Year", y = "")






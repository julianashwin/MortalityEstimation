setwd("/Users/julianashwin/Documents/GitHub/MortalityEstimation/")
rm(list=ls())

library(ggplot2)
library(ggpubr)
library(tidyverse)
library(stargazer)
library(lfe)


"
Define color scheme
"
col_scheme <- c("Australia" = "darkolivegreen4", 
                "Canada" = "pink",
                "France" = "blue3", "United Kingdom" = "darkgoldenrod4", 
                "Hong Kong" = "lightgoldenrod", 
                "Italy" = "forestgreen", 
                "Japan" = "red","New Zealand" = "black", "Russia" = "firebrick",
                "Sweden" = "yellow", "United States of America" = "cornflowerblue",
                "Other" = "gray")

"
Import data and results
"
## Import mortality data
mort_df <- read.csv("data/clean/all_lifetab_5y.csv")
bp_df <- read.csv("data/clean/bp_5y.csv")
pop_df <- as_tibble(mort_df) %>%
  group_by(code, year) %>%
  summarise(pop = sum(Total, na.rm = T))

## Import BPLE results
bp_parests_df <- read.csv("figures/benchmark/siler_i2_params_ber.csv")
bp_params_df <- read.csv("figures/benchmark/siler_i2drift_preds.csv")
bp_decomp_df <- read.csv("figures/benchmark/siler_i2drift_decomp_ber.csv")
bp_LEgrad_df <- read.csv("figures/benchmark/siler_i2drift_LEgrads.csv")

## Import parameter estimates for each country
import_files <- dir("figures/countries/")
import_files <- import_files[which(str_detect(import_files, "_preds.csv"))]
countries_df <- data.frame(matrix(NA,nrow=0,ncol = 17))
names(countries_df) <- c("parameter", "code", "forecast", "year", "mean", "min", "median", "max", "nmissing", "eltype",
                         "std", "pc975", "pc025", "pc85", "pc15", "pc75", "pc25")
for (ii in 1:length(import_files)){
  filename <- import_files[ii]
  country_df <- read.csv(paste0("figures/countries/", filename), stringsAsFactors = F)
  country_df$code <- str_remove(filename, "_i2_preds.csv")
  country_df <- country_df[,names(countries_df)]
  countries_df <- rbind(countries_df, country_df)
}
# Combine into dataframe for all countries
all_pars_df <- merge(countries_df, 
                unique(mort_df[which(mort_df$age == 0),c("code", "name")]), by = "code")
all_pars_df <- merge(all_pars_df, pop_df, by = c("code", "year"), all.x = T)
all_pars_df <- all_pars_df[which(all_pars_df$year > 0),]
all_pars_df$code <- str_replace(all_pars_df$code, "NZL_NM", "NZL")
#Forecast label for pretty legends
all_pars_df$Forecast <- "Estimate"
all_pars_df$Forecast[which(all_pars_df$forecast == 1)] <- "Forecast"
all_pars_df <- all_pars_df[which(all_pars_df$year > 1900),]
rm(countries_df, country_df)

## Import decomposition for each country
import_files <- dir("figures/countries/")
import_files <- import_files[which(str_detect(import_files, "decomp_pred.csv"))]
all_decomp_df <- data.frame(matrix(NA,nrow=0,ncol = 50))
names(all_decomp_df) <- c("year", "code", "forecast", "B", "b", "C", "c", "d", "σ", 
                         "Lmed", "LE_mod", "H_mod", "h_mod", "Lstar_mod", "Lmed_mod",
                         "LE_b", "LE_B", "LE_c", "LE_C", "LE_d", 
                         "H_b", "H_B", "H_c", "H_C", "H_d",
                         "h_b", "h_B", "h_c", "h_C", "h_d", 
                         "Lstar_b", "Lstar_B", "Lstar_c", "Lstar_C", "Lstar_d", 
                         "Lmed_b", "Lmed_B", "Lmed_c", "Lmed_C", "Lmed_d", 
                         "Δb", "ΔB", "Δc", "ΔC", "Δd", 
                         "ΔLE_mod", "ΔH_mod", "Δh_mod", "ΔLstar_mod", "ΔLmed_mod")
for (ii in 1:length(import_files)){
  filename <- import_files[ii]
  country_df <- read.csv(paste0("figures/countries/", filename), stringsAsFactors = F)
  country_df$code <- str_remove(filename, "_i2_decomp_pred.csv")
  country_df <- country_df[,names(all_decomp_df)]
  all_decomp_df <- rbind(all_decomp_df, country_df)
}
all_decomp_df <- all_decomp_df %>%
  inner_join(unique(mort_df[,c("code", "name")]), by = "code") %>%
  filter(year > 1918) %>%
  mutate(DeltaLE_c = Δc*LE_c, DeltaLE_C = ΔC*LE_C,
         year_group = case_when(
           year <= 1943 ~ "1918-1943",
           year > 1943 & year <= 1968 ~ "1943-1968",
           year > 1968 & year <= 1993 ~ "1968-1993",
           year > 1993 & year <= 2018 ~ "1993-2018",
           TRUE ~ "pred"))
table(all_decomp_df$year_group)
rm(country_df)

## Import panel data with econ variables
all_panel_df <- read.csv("data/results/siler_econ_panel.csv")
all_panel_df[is.na(all_panel_df)] <- NA
all_panel_df <- as_tibble(all_panel_df) %>%
  select(-pop) %>%
  left_join(pop_df, by = c("code", "year")) %>%
  relocate(pop, .after = "code")

## Forecasts for BPLE 
forecasts_df <- read.csv("data/results/bp_forecasts.csv")

## Forecasts for countries 
int_forecasts_df <- read.csv("data/results/int_forecasts.csv")


"
Figure 1: Best practice mortality and life expectancy (1900-2019)
"
# log mortality
p1 <- bp_df %>%
  filter(year > 1900) %>%
  filter(mx_f > 0) %>%
  ggplot() + theme_bw() +
  geom_line(aes(x = age, y = log(mx_f), group = year, color = year)) + 
  scale_color_gradientn(colours = rainbow(5), name = "Year",
                        breaks=c(1900,1925,1950,1975,2000,2019),
                        labels=c(1900,1925,1950,1975,2000,2019),
                        limits=c(1900,2019)) +
  xlab("Age") + ylab("log Mortality rate") + ggtitle("Mortality")
# Survival
p2 <- bp_df %>%
  filter(year > 1900) %>%
  ggplot() + theme_bw() +
  geom_line(aes(x = age, y = lx_f, group = year, color = year)) + 
  scale_color_gradientn(colours = rainbow(5), name = "Year",
                        breaks=c(1900,1925,1950,1975,2000,2019),
                        labels=c(1900,1925,1950,1975,2000,2019),
                        limits=c(1900,2019)) +
  ylim(c(0,1)) + labs(x = "Age", y = "Survival rate", title = "Survival")
# Life expectancy
p3 <- bp_df %>%
  filter(year > 1900 & age == 0) %>%
  ggplot() + theme_bw() +
  geom_smooth(aes(x = year, y = ex_f), method = "lm") +
  geom_point(aes(x = year, y = ex_f, color = name)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  scale_color_discrete(name = "Country") + ggtitle("Life Expectancy") +
  xlab("Year") + ylab("Life Expectancy at birth")
# Lifespan equality
p4 <- bp_df %>%
  filter(year > 1900 & age == 0) %>%
  ggplot() + theme_bw() +
  geom_smooth(aes(x = year, y = -log(Hx_f)), method = "loess") +
  geom_point(aes(x = year, y = -log(Hx_f), color = name)) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  scale_color_discrete(name = "Country") + ggtitle("Lifespan Equality") +
  xlab("Year") + ylab("Lifespan Equality at birth")
# Lifespan
p5 <-  bp_df %>%
  filter(year > 1900 & age == 0) %>%
  ggplot() + theme_bw() +
  geom_smooth(aes(x = year, y = l_99p9_f), method = "loess") +
  geom_point(aes(x = year, y = l_99p9_f, color = name)) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  scale_color_discrete(name = "Country") + ggtitle("UBAD (99.9%)") +
  xlab("Year") + ylab("UBAD (99.9%)")

# Combine panels
p_upper <- ggarrange(p1, p2, nrow = 1, ncol=2, common.legend = TRUE,
                  legend = "right")
p_lower <- ggarrange(p3, p4, p5, nrow = 1, ncol=3, common.legend = TRUE,
                  legend = "right")

ggarrange(p_upper, p_lower, nrow = 2, ncol=1, common.legend = FALSE)
ggsave("figures_paper/best_practice_5y_summary.pdf", width = 9, height = 6)
rm(p1,p2,p3,p4,p5,p_upper, p_lower)

"
Figure 2: Siler parameter estimates for best practice mortality (1900-2019)
"
bp_params_df %>% 
  filter(year >1900 & year < 2020 & parameter %in% c("b", "B", "c", "C", "d", "σ")) %>%
  mutate(parameter = factor(parameter, levels =c("b", "c", "d", "B", "C", "σ"), ordered = T)) %>%
  ggplot() + theme_bw() + facet_wrap(~parameter, scales = "free_y") + 
  geom_ribbon(aes(x = year, ymin=pc025, ymax=pc975), alpha = 0.2) +
  geom_ribbon(aes(x = year, ymin=pc15, ymax=pc85), alpha = 0.3) + 
  geom_line(aes(x = year, y = median)) + 
  labs(x = "Year", y = "")
ggsave("figures_paper/siler_i2drift_params.pdf", width = 8, height = 4, device = cairo_pdf)


"
Figure 3: Historical decomposition of changes (1900-2019)
"
# Get changes in parameters and gradients
change_df <- bp_decomp_df %>%
  select(year, Δb, ΔB, Δc, ΔC, Δd) %>%
  pivot_longer( cols = -year, names_to = "parameter", values_to = "change") %>%
  mutate(parameter = str_remove_all(parameter,"Δ"))
LEgrad_df <- bp_decomp_df %>%
  select(year, LE_b, LE_B, LE_c, LE_C, LE_d) %>%
  pivot_longer( cols = -year, names_to = "parameter", values_to = "LE_grad") %>%
  mutate(parameter = str_remove_all(parameter,"LE_"))
hgrad_df <- bp_decomp_df %>%
  select(year, h_b, h_B, h_c, h_C, h_d) %>%
  pivot_longer( cols = -year, names_to = "parameter", values_to = "h_grad") %>%
  mutate(parameter = str_remove_all(parameter,"h_"))
Lstargrad_df <- bp_decomp_df %>%
  select(year, Lstar_b, Lstar_B, Lstar_c, Lstar_C, Lstar_d) %>%
  pivot_longer( cols = -year, names_to = "parameter", values_to = "Lstar_grad") %>%
  mutate(parameter = str_remove_all(parameter,"Lstar_"))
# Combine into long df
decomp_long <- bp_decomp_df %>%
  select(year, LE_mod, h_mod, Lstar_mod, ΔLE_mod, Δh_mod, ΔLstar_mod) %>%
  full_join(change_df, by = "year") %>%
  left_join(LEgrad_df, by = c("year", "parameter")) %>%
  left_join(hgrad_df, by = c("year", "parameter")) %>%
  left_join(Lstargrad_df, by = c("year", "parameter"))
# Plot decomps
p1 <- decomp_long %>%
  ggplot() + theme_bw() + 
  geom_bar(aes(x = year, y = change*LE_grad, fill = parameter), 
           stat = "identity", position = "stack") + 
  geom_line(aes(x = year, y = ΔLE_mod)) +
  labs(title = "Life Expectancy", x = "Year", y = "ΔLE", fill = "Parameter")
p2 <- decomp_long %>%
  ggplot() + theme_bw() + 
  geom_bar(aes(x = year, y = change*Lstar_grad, fill = parameter), 
           stat = "identity", position = "stack") + 
  geom_line(aes(x = year, y = ΔLstar_mod)) +
  labs(title = "Lifespan", x = "Year", y = "ΔL*", fill = "Parameter")
p3 <- decomp_long %>%
  ggplot() + theme_bw() + 
  geom_bar(aes(x = year, y = change*h_grad, fill = parameter), 
           stat = "identity", position = "stack") + 
  geom_line(aes(x = year, y = Δh_mod)) +
  labs(title = "Lifespan Equality", x = "Year", y = "Δh", fill = "Parameter")
ggarrange(p1, p2, p3, nrow = 1, common.legend = TRUE, legend = "right")
ggsave("figures_paper/siler_i2drift_decomp.pdf", width = 12, height = 4, device = cairo_pdf)

rm(p1,p2,p3,decomp_long,change_df,LEgrad_df,Lstargrad_df,hgrad_df)


"
Figure 4: 
"
decomp_table <- all_decomp_df[which(!is.na(all_decomp_df$DeltaLE_c)),]
decomp_table <- aggregate(decomp_table[,c("DeltaLE_c", "DeltaLE_C", "ΔLE_mod")],
                          by = list(name = decomp_table$name, year = decomp_table$year_group), 
                          FUN = sum, drop = FALSE)
decomp_table$propLE_c <- round(decomp_table$DeltaLE_c/decomp_table$ΔLE_mod,3)
decomp_table$propLE_C <- round(decomp_table$DeltaLE_C/decomp_table$ΔLE_mod,3)
decomp_table <- pivot_wider(decomp_table, id_cols = c(name), names_from = year, 
                            names_glue = "{.value}_{year}",
                            values_from = c(propLE_c, propLE_C))
decomp_table <- decomp_table[,c("name", "propLE_c_1918-1943",  "propLE_C_1918-1943",
                                "propLE_c_1943-1968", "propLE_C_1943-1968",
                                "propLE_c_1968-1993", "propLE_C_1968-1993",
                                "propLE_c_1993-2018", "propLE_C_1993-2018")]
decomp_table <- merge(decomp_table, all_decomp_df[which(all_decomp_df$year == 2018),c("name", "C")],
                      by = "name")
decomp_table$C <- round(decomp_table$C, 1)

# Export table as tex
stargazer(as.matrix(decomp_table), table.placement = "H", column.sep.width = "2")
# Or plot as heatmap
as_tibble(decomp_table) %>%
  rename(country = name) %>%
  pivot_longer(cols = -country) %>%
  mutate(country = factor(country)) %>%
  mutate(parameter = str_c("Proportion~of~LE~gains~due~to~", str_sub(name, 8, 8), "[t]")) %>%
  mutate(interval = str_replace(str_sub(name, 10, 18), "\\.", "-")) %>%
  mutate(interval = case_when(interval == "" ~ "Latest", TRUE ~ interval),
         parameter = case_when(parameter == "Proportion~of~LE~gains~due~to~[t]" ~ "Latest~C[t]", TRUE ~ parameter)) %>%
  mutate(parameter = factor(parameter, levels = c("Proportion~of~LE~gains~due~to~c[t]", "Proportion~of~LE~gains~due~to~C[t]",
                                                  "Latest~C[t]"))) %>%
  mutate(value_latest = case_when(name == "C" ~ value, TRUE ~ NA_real_)) %>%
  mutate(value = case_when(name == "C" ~ NA_real_, TRUE ~ value)) %>%
  ggplot(aes(y =fct_rev(country), x = interval)) + theme_bw() + 
  facet_grid(~parameter, labeller = label_parsed, scales = "free_x", space="free_x") +
  geom_tile(aes(fill = value), color = "black") +
  geom_text(aes(label = value), color = "black", size = 2.5) +
  geom_text(aes(label = value_latest), color = "black", size = 2.5) +
  scale_fill_gradient2(low = "white", high = "red", na.value = 'white') +  #, midpoint = 0)
  theme(axis.text.x=element_text(angle=45,hjust=1)) +
  labs(x = "Interval", y = "Country", fill = "Proportion of LE gains \ndue to parameter") 

ggsave("figures_paper/decomp_table.pdf", width = 8, height = 8)
rm(decomp_table)





"
Figure 5: Evolution of $c$ and $C$ over time
"
plot_df <- data.frame(pivot_wider(all_pars_df, id_cols = c(name, code, year, Forecast), 
                                  names_from = parameter, values_from = median))
extra_obs <- plot_df[which((plot_df$year == 2018 & !(plot_df$code %in% c("NZL","RUS")))|
                             (plot_df$year == 2013 & plot_df$code %in% c("NZL","RUS"))),]
extra_obs$Forecast <- "Estimate"
plot_df <- rbind(plot_df, extra_obs)
plot_df$h <- -log(plot_df$H)
plot_df$name[which(!(plot_df$name %in% names(col_scheme)))] <- "Other"
plot_df$name <- factor(plot_df$name, levels = names(col_scheme), ordered = TRUE)
# Include some years as labels
plot_df$year_label <- NA
plot_df$year_label[which(plot_df$year==1903)] <- 1903
plot_df$year_label[which(plot_df$year==1933)] <- 1933
plot_df$year_label[which(plot_df$year==1963)] <- 1963
plot_df$year_label[which(plot_df$year==1993)] <- 1993
plot_df$year_label[which(plot_df$year==2023)] <- 2023



# Cc path
plot_df %>%
  ggplot() + theme_bw() + 
  scale_color_manual("Country", values = col_scheme) + 
  geom_path(data = plot_df[which(plot_df$name == "Other"),], alpha = 0.5,
            aes(x = C, y = c, color = name, linetype = Forecast)) + 
  geom_path(data = plot_df[which(plot_df$name != "Other"),],
            aes(x = C, y = c, color = name, linetype = Forecast)) + 
  geom_text(data = plot_df[which(plot_df$name != "Other"),],
            aes(x = C, y = c, color = name, label = year_label), show.legend=FALSE, size = 2) +
  xlab("Timing of aging (C)") + ylab("Speed of aging (c)") + 
  guides(color=guide_legend(ncol=2))
ggsave("figures_paper/Cc_paths.pdf", width = 10, height = 4)

rm(plot_df, extra_obs)


"
Figure 6: Evolution of L*, C and LE over time
"
plot_df <- data.frame(pivot_wider(all_pars_df, id_cols = c(name, code, year, Forecast), 
                                  names_from = parameter, values_from = median))
extra_obs <- plot_df[which((plot_df$year == 2018 & !(plot_df$code %in% c("NZL","RUS")))|
                             (plot_df$year == 2013 & plot_df$code %in% c("NZL","RUS"))),]
extra_obs$Forecast <- "Estimate"
plot_df <- rbind(plot_df, extra_obs)
plot_df$h <- -log(plot_df$H)
plot_df$name[which(!(plot_df$name %in% names(col_scheme)))] <- "Other"
plot_df$name <- factor(plot_df$name, levels = names(col_scheme), ordered = TRUE)
# Include some years as labels
plot_df$year_label <- NA
plot_df$year_label[which(plot_df$year==1903)] <- 1903
plot_df$year_label[which(plot_df$year==1933)] <- 1933
plot_df$year_label[which(plot_df$year==1963)] <- 1963
plot_df$year_label[which(plot_df$year==1993)] <- 1993
plot_df$year_label[which(plot_df$year==2023)] <- 2023

p1 <- ggplot(plot_df) + theme_bw() + 
  scale_color_manual("Country", values = col_scheme) + 
  geom_path(data = plot_df[which(plot_df$name == "Other"),], alpha = 0.5,
            aes(x = LE, y = Lstar_99p9, color = name, linetype = Forecast)) + 
  geom_path(data = plot_df[which(plot_df$name != "Other"),],
            aes(x = LE, y = Lstar_99p9, color = name, linetype = Forecast)) + 
  geom_text(data = plot_df[which(plot_df$name != "Other"),],
            aes(x = LE, y = Lstar_99p9, color = name, label = year_label), show.legend=FALSE, size = 2) +
  xlab("Life expectancy at birth") + ylab("UBAD") + 
  guides(color=guide_legend(ncol=2)) + 
  scale_x_continuous(limits = c(46,93),breaks=seq(0,130,10)) + 
  scale_y_continuous(limits = c(96,118),breaks=seq(0,130,10)) 
p2 <- ggplot(plot_df) + theme_bw() + 
  scale_color_manual("Country", values = col_scheme) + 
  geom_path(data = plot_df[which(plot_df$name == "Other"),], alpha = 0.5,
            aes(x = LE, y = Lstar_99p9-C, color = name, linetype = Forecast)) + 
  geom_path(data = plot_df[which(plot_df$name != "Other"),],
            aes(x = LE, y = Lstar_99p9-C, color = name, linetype = Forecast)) + 
  geom_text(data = plot_df[which(plot_df$name != "Other"),],
            aes(x = LE, y = Lstar_99p9-C, color = name, label = year_label), show.legend=FALSE, size = 2) +
  xlab("Life expectancy at birth") + ylab("UBAD - C") + 
  guides(color=guide_legend(ncol=2))
ggarrange(p1, p2, nrow = 1, common.legend = TRUE, legend = "right")
ggsave("figures_paper/Lstar_C_LE_paths.pdf", width = 10, height = 4)

rm(p1, p2, plot_df, extra_obs)


"
Figure 7: Forecast LE for Best Practice case
"

bp_params_df %>%
  filter(parameter %in% c("LE")) %>%
  ggplot() + theme_bw() +
  geom_ribbon(aes(x = year, ymin=pc025, ymax=pc975), alpha = 0.2) +
  geom_ribbon(aes(x = year, ymin=pc15, ymax=pc85), alpha = 0.3) + 
  geom_line(aes(x = year, y = median)) + 
  labs(x = "Year", y = "BPLE") + 
  geom_vline(aes(xintercept = 2018), linetype = "dashed")
ggsave("figures_paper/siler_i2drift_bple_fcast.pdf", width = 5, height = 3, device = cairo_pdf)



"
Figure 8: Gradients
"
p1 <- bp_LEgrad_df %>%
  mutate(Forecast = case_when(year < 2020 ~ "Estimate",
                              year >= 2020 ~ "Forecast")) %>%
  ggplot() + theme_bw() +
  geom_line(aes(x = age, y = LE_Cs, group = year, color = year, linetype = Forecast)) +
  geom_hline(yintercept=0, linetype="dashed", color = "black") +
  scale_color_gradientn(colours = rainbow(5), name = "Year",
                        breaks=c(1900, 1950, 2000,2049),
                        labels=c(1900,1950,2000,2049),
                        limits=c(1900,2049)) +
  xlab("Age") + ylab("Gradient") + 
  ggtitle(expression(Gradient~of~LE~wrt~C[t]))
p2 <- bp_LEgrad_df %>%
  mutate(Forecast = case_when(year < 2020 ~ "Estimate",
                              year >= 2020 ~ "Forecast")) %>%
  ggplot() + theme_bw() +
  geom_line(aes(x = age, y = LE_cs, group = year, color = year, linetype = Forecast)) +
  geom_hline(yintercept=0, linetype="dashed", color = "black") +
  scale_color_gradientn(colours = rainbow(5), name = "Year",
                        breaks=c(1900, 1950, 2000,2049),
                        labels=c(1900,1950,2000,2049),
                        limits=c(1900,2049)) +
  xlab("Age") + ylab("Gradient") + 
  ggtitle(expression(Gradient~of~LE~wrt~c[t]))
ggarrange(p1, p2, nrow = 1, ncol = 2, common.legend = TRUE, legend = "right")
ggsave(paste0("figures_paper/LEgrads_cC.pdf"), width = 10, height = 4)

rm(p1,p2)


"
Figure 9: Evolution of $h(0)$ over time
"
plot_df <- all_pars_df[which(all_pars_df$parameter == "H"),]
extra_obs <- plot_df[which((plot_df$year == 2018 & !(plot_df$code %in% c("NZL","RUS")))|
                             (plot_df$year == 2013 & plot_df$code %in% c("NZL","RUS"))),]
extra_obs$Forecast <- "Estimate"
plot_df <- rbind(plot_df, extra_obs)
plot_df$name[which(!(plot_df$name %in% names(col_scheme)))] <- "Other"
plot_df$name <- factor(plot_df$name, levels = names(col_scheme), ordered = TRUE)


ggplot(plot_df) + theme_bw() + 
  scale_color_manual("Country", values = col_scheme) + guides(color=guide_legend(ncol=2)) + 
  geom_line(data = plot_df[which(plot_df$name == "Other"),], alpha = 0.5,
            aes(x = year, y = -log(median), color = name, 
                group = interaction(code, Forecast), linetype = Forecast)) +
  geom_line(data = plot_df[which(plot_df$name != "Other"),],
            aes(x = year, y = -log(median), color = name, 
                group = interaction(code, Forecast), linetype = Forecast)) +
  xlab("Year") + ylab("Lifespan equality at birth")
ggsave("figures_paper/h_international.pdf", width = 10, height = 4)

rm(plot_df, extra_obs)




"
Appendix A: Data availability by country, see clean_econ_data.R
"


"
Appendix B: 1) Convergence of HMC sampler for each Best Practice estimated parameter
"
ggplot(bp_parests_df) + theme_bw() +
  geom_density(aes(x = rhat)) + 
  ylab("Density") + xlab(expression(hat(R)))
ggsave("figures_paper/rhat_convergence.pdf", width = 4, height = 4)



"
Appendix D: 1) Siler parameter drift estimates for best practice mortality
"
bp_params_df %>% 
  mutate(parameter = str_replace_all(parameter, "α_", "Drift ")) %>%
  filter(year >1900 & year < 2020 & parameter %in% c("Drift b", "Drift B", "Drift c", "Drift C", 
                                                     "Drift d", "Drift σ")) %>%
  mutate(parameter = factor(parameter, levels =c("Drift b", "Drift c", "Drift d", "Drift B", 
                                                 "Drift C", "Drift σ"), ordered = T)) %>%
  ggplot() + theme_bw() + facet_wrap(~parameter, scales = "free_y") + 
  geom_ribbon(aes(x = year, ymin=pc025, ymax=pc975), alpha = 0.2) +
  geom_ribbon(aes(x = year, ymin=pc15, ymax=pc85), alpha = 0.3) + 
  geom_line(aes(x = year, y = median)) + 
  labs(x = "Year", y = "")
ggsave("figures_paper/siler_i2drift_drifts.pdf", width = 8, height = 4, device = cairo_pdf)

"
Appendix D: Table 2) Posterior distributions of variance terms
"
plot_df <- bp_params_df %>% 
  mutate(parameter = str_replace_all(parameter, "σ", "sigma"),
         parameter = str_replace_all(parameter, "α", "alpha")) %>%
  filter(str_detect(parameter, "sigma_")) %>%
  mutate(across(where(is.numeric), round, digits=3)) %>%
  select(parameter, median, pc025, pc15, pc85, pc975)

stargazer(as.matrix(plot_df), table.placement = "H", label = "tab:variance_terms",
          title = "Posterior distributions of variance terms")

rm(plot_df)

"
Appendix D: Table 3) Economic Determinants of Senescent Mortality
"

reg_df <- all_panel_df %>%
  select(code, name, year, c, C, income_gini, gdp_pc,pop, health_share) %>%
  drop_na %>% group_by(year) %>% 
  mutate(pop_prop = pop/sum(pop), US = as.numeric(code == "USA"))
reg_lag_df <- all_panel_df %>%
  select(code, name, year, c, C, income_gini_lag, gdp_pc_lag, pop, health_share_lag) %>%
  drop_na %>% group_by(year) %>% 
  mutate(pop_prop = pop/sum(pop), US = as.numeric(code == "USA"))


model1 <- (felm(c ~ income_gini + log(gdp_pc)-US |
                  name + year, data = reg_df, weights = reg_df$pop_prop))
model2 <- (felm(C ~ income_gini + log(gdp_pc)-US |
                  name + year, data = reg_df, weights = reg_df$pop_prop))
model3 <- (felm(c ~ income_gini + log(gdp_pc) + health_share|
                  name + year, data = reg_df, weights = reg_df$pop_prop))
model4 <- (felm(C ~ income_gini + log(gdp_pc) + health_share|
                  name + year, data = reg_df, weights = reg_df$pop_prop))
model5 <- (felm(c ~ income_gini_lag + log(gdp_pc_lag) + health_share_lag|
                  name + year, data = reg_lag_df, weights = reg_lag_df$pop_prop))
model6 <- (felm(C ~ income_gini_lag + log(gdp_pc_lag) + health_share_lag|
                  name + year, data = reg_lag_df, weights = reg_lag_df$pop_prop))
models <- list(model1, model2, model3, model4, model5, model6)
stargazer(models, table.placement = "H",
          df = FALSE, title = "C and c on economic variables", label = "tab:cC_econ",
          font.size = "scriptsize")

rm(reg_df, reg_lag_df, models, model1, model2, model3, model4, model5, model6, model7, model8)



"
Appendix E: see figs_forecasts.R
"
model_cols <- c("Siler" = "purple", "Lee-Carter (demo)" = "red", "Lee-Carter (SMM)" = "deeppink", 
                "FDM" = "forestgreen", "CBD" = "darkorange1",
                "Lee-Carter (dt)" = "gold", "Lee-Carter (dxt)" = "green", 
                "Lee-Carter (e0)" = "blue")


# Plot the errors
fe_full_plt <-
  tibble(forecasts_df) %>%
  filter(age == 0 & n_ahead > 0) %>%
  ggplot() + theme_bw() + 
  scale_fill_manual("Models", values = model_cols) +
  geom_bar(aes(x = n_ahead-1.5, y = sqrt(siler_fe2), fill = "Siler"), 
           stat = "summary", fun = mean, width = 0.5) +
  geom_bar(aes(x = n_ahead-1.0, y = sqrt(LC_fe2), fill = "Lee-Carter (demo)"), 
           stat = "summary", fun = mean, width = 0.5, alpha = 0.5) +
  geom_bar(aes(x = n_ahead-0.5, y = sqrt(LC_SMM_fe2), fill = "Lee-Carter (SMM)"), 
           stat = "summary", fun = mean, width = 0.5, alpha = 0.5) +
  geom_bar(aes(x = n_ahead+0, y = sqrt(LC_dt_fe2), fill = "Lee-Carter (dt)"), 
           stat = "summary", fun = mean, width = 0.5, alpha = 0.5) +
  geom_bar(aes(x = n_ahead+0.5, y = sqrt(LC_dxt_fe2), fill = "Lee-Carter (dxt)"), 
           stat = "summary", fun = mean, width = 0.5, alpha = 0.5) +
  geom_bar(aes(x = n_ahead+1.0, y = sqrt(LC_e0_fe2), fill = "Lee-Carter (e0)"), 
           stat = "summary", fun = mean, width = 0.5, alpha = 0.5) +
  geom_bar(aes(x = n_ahead+1.5, y = sqrt(FDM_fe2), fill = "FDM"), 
           stat = "summary", fun = mean, width = 0.5, alpha = 0.5) +
  geom_bar(aes(x = n_ahead+2.0, y = sqrt(CBD_SMM_fe2), fill = "CBD"), 
           stat = "summary", fun = mean, width = 0.5, alpha = 0.5) +
  xlab("Years ahead") + ylab("RMSE")+ ggtitle("Full sample")


fe_rolling_plt <-
  tibble(forecasts_df) %>%
  filter(age == 0 & n_ahead > 0) %>%
  ggplot() + theme_bw() + 
  scale_fill_manual("Models", values = model_cols) +
  geom_bar(aes(x = n_ahead-1.5, y = sqrt(siler_fe2), fill = "Siler"), 
           stat = "summary", fun = mean, width = 0.5) +
  geom_bar(aes(x = n_ahead-1.0, y = sqrt(LC_r_fe2), fill = "Lee-Carter (demo)"), 
           stat = "summary", fun = mean, width = 0.5, alpha = 0.5) +
  geom_bar(aes(x = n_ahead-0.5, y = sqrt(LC_SMM_r_fe2), fill = "Lee-Carter (SMM)"), 
           stat = "summary", fun = mean, width = 0.5, alpha = 0.5) +
  geom_bar(aes(x = n_ahead+0, y = sqrt(LC_dt_r_fe2), fill = "Lee-Carter (dt)"), 
           stat = "summary", fun = mean, width = 0.5, alpha = 0.5) +
  geom_bar(aes(x = n_ahead+0.5, y = sqrt(LC_dxt_r_fe2), fill = "Lee-Carter (dxt)"), 
           stat = "summary", fun = mean, width = 0.5, alpha = 0.5) +
  geom_bar(aes(x = n_ahead+1.0, y = sqrt(LC_e0_r_fe2), fill = "Lee-Carter (e0)"), 
           stat = "summary", fun = mean, width = 0.5, alpha = 0.5) +
  geom_bar(aes(x = n_ahead+1.5, y = sqrt(FDM_r_fe2), fill = "FDM"), 
           stat = "summary", fun = mean, width = 0.5, alpha = 0.5) +
  geom_bar(aes(x = n_ahead+2.0, y = sqrt(CBD_SMM_r_fe2), fill = "CBD"), 
           stat = "summary", fun = mean, width = 0.5, alpha = 0.5) +
  xlab("Years ahead") + ylab("RMSE")+ ggtitle("Rolling window")

ggarrange(fe_full_plt, fe_rolling_plt, nrow = 1, common.legend = T, legend = "right")
ggsave("figures_paper/BP_forecast_errors.pdf", width = 8, height = 4, device = cairo_pdf)




"
Appendix E: Table 4) Increases in LE (years per decade) 
"
increase_df <- all_pars_df %>%
  filter(year %in% c(1988, 1998, 2008, 2018, 2028, 2038, 2048) & 
           parameter == "LE" & 
           name != "Russia" & name != "Ukraine") %>%
  mutate(median = round(median, 1)) %>%
  select(name, year, median) %>%
  pivot_wider(id_cols = name, names_from = year, values_from = median) %>%
  mutate(`1988-1998` = round(`1998`-`1988`, 3), 
         `1998-2008` = round(`2008`-`1998`, 3), 
         `2008-2018` = round(`2018`-`2008`, 3), 
         `2018-2028` = round(`2028`-`2018`, 3), 
         `2028-2038` = round(`2038`-`2028`, 3), 
         `2038-2048` = round(`2048`-`2038`, 3),
         `2018-2048` = (`2048`-`2018`)/3,
         `1988-2018` = (`2018`-`1988`)/3,
         trend_change = (`2048`-`2018`)/3 - (`2018`-`1988`)/3) %>%
  select(-`1988`, -`1998`, -`2008`, -`2018`, -`2028`, -`2038`, -`2048`) %>%
  na.omit() %>%
  pivot_longer(cols = -name, names_to = "Period", values_to = "Change")  %>%
  mutate(Change_bucket = cut(Change, breaks = c(-Inf, 0, 0.5, 1, 1.5, 2, 2.5, Inf)))

increase_table <- table(select(increase_df, Period, Change_bucket))
increase_tab_df <- data.frame(matrix(increase_table, ncol = 7))
colnames(increase_tab_df) <- colnames(increase_table)
rownames(increase_tab_df) <- rownames(increase_table)

stargazer(as.matrix(increase_tab_df), table.placement = "H", 
          label = "intnlforecast", title = "Increases in LE (years per decade)")

rm(increase_df, increase_table, increase_tab_df)




"
Appendix E: see figs_forecasts.R
"
model_cols <- c("Siler" = "purple", "Lee-Carter (demo)" = "red", "Lee-Carter (SMM)" = "deeppink", 
                "FDM" = "forestgreen", "CBD" = "darkorange1",
                "Lee-Carter (dt)" = "gold", "Lee-Carter (dxt)" = "green", 
                "Lee-Carter (e0)" = "blue")


# Plot the errors
fe_mean_plt <-
  tibble(int_forecasts_df) %>%
  filter(age == 0 & n_ahead > 0) %>%
  ggplot() + theme_bw() + 
  scale_fill_manual("Models", values = model_cols) +
  geom_bar(aes(x = n_ahead-1.5, y = siler_fe, fill = "Siler"), 
           stat = "summary", fun = mean, width = 0.5) +
  geom_bar(aes(x = n_ahead-1.0, y = LC_r_fe, fill = "Lee-Carter (demo)"), 
           stat = "summary", fun = mean, width = 0.5, alpha = 0.5) +
  geom_bar(aes(x = n_ahead-0.5, y = (LC_SMM_r_fe), fill = "Lee-Carter (SMM)"), 
           stat = "summary", fun = mean, width = 0.5, alpha = 0.5) +
  geom_bar(aes(x = n_ahead+0, y = (LC_dt_r_fe), fill = "Lee-Carter (dt)"), 
           stat = "summary", fun = mean, width = 0.5, alpha = 0.5) +
  geom_bar(aes(x = n_ahead+0.5, y = (LC_dxt_r_fe), fill = "Lee-Carter (dxt)"), 
           stat = "summary", fun = mean, width = 0.5, alpha = 0.5) +
  geom_bar(aes(x = n_ahead+1.0, y = (LC_e0_r_fe), fill = "Lee-Carter (e0)"), 
           stat = "summary", fun = mean, width = 0.5, alpha = 0.5) +
  geom_bar(aes(x = n_ahead+1.5, y = (FDM_r_fe), fill = "FDM"), 
           stat = "summary", fun = mean, width = 0.5, alpha = 0.5) +
  geom_bar(aes(x = n_ahead+2.0, y = (CBD_SMM_r_fe), fill = "CBD"), 
           stat = "summary", fun = mean, width = 0.5, alpha = 0.5) +
  xlab("Years ahead") + ylab("Mean Forecast error")+ ggtitle("Bias")


fe_rmse_plt <-
  tibble(int_forecasts_df) %>%
  filter(age == 0 & n_ahead > 0) %>%
  ggplot() + theme_bw() + 
  scale_fill_manual("Models", values = model_cols) +
  geom_bar(aes(x = n_ahead-1.5, y = sqrt(siler_fe2), fill = "Siler"), 
           stat = "summary", fun = mean, width = 0.5) +
  geom_bar(aes(x = n_ahead-1.0, y = sqrt(LC_r_fe2), fill = "Lee-Carter (demo)"), 
           stat = "summary", fun = mean, width = 0.5, alpha = 0.5) +
  geom_bar(aes(x = n_ahead-0.5, y = sqrt(LC_SMM_r_fe2), fill = "Lee-Carter (SMM)"), 
           stat = "summary", fun = mean, width = 0.5, alpha = 0.5) +
  geom_bar(aes(x = n_ahead+0, y = sqrt(LC_dt_r_fe2), fill = "Lee-Carter (dt)"), 
           stat = "summary", fun = mean, width = 0.5, alpha = 0.5) +
  geom_bar(aes(x = n_ahead+0.5, y = sqrt(LC_dxt_r_fe2), fill = "Lee-Carter (dxt)"), 
           stat = "summary", fun = mean, width = 0.5, alpha = 0.5) +
  geom_bar(aes(x = n_ahead+1.0, y = sqrt(LC_e0_r_fe2), fill = "Lee-Carter (e0)"), 
           stat = "summary", fun = mean, width = 0.5, alpha = 0.5) +
  geom_bar(aes(x = n_ahead+1.5, y = sqrt(FDM_r_fe2), fill = "FDM"), 
           stat = "summary", fun = mean, width = 0.5, alpha = 0.5) +
  geom_bar(aes(x = n_ahead+2.0, y = sqrt(CBD_SMM_r_fe2), fill = "CBD"), 
           stat = "summary", fun = mean, width = 0.5, alpha = 0.5) +
  xlab("Years ahead") + ylab("RMSE")+ ggtitle("RMSE")

ggarrange(fe_mean_plt, fe_rmse_plt, nrow = 1, common.legend = T, legend = "right")
ggsave("figures_paper/country_forecast_errors.pdf", width = 8, height = 4, device = cairo_pdf)





"
Appendix E: Table 5)
"
obs <- which(!is.na(int_forecasts_df$siler_fe) & int_forecasts_df$age == 0 & 
               int_forecasts_df$est_year > 1960)
mse_tab <- aggregate(int_forecasts_df[obs,c("siler_fe2","LC_r_fe2", "LC_SMM_r_fe2", 
                                            "LC_dt_r_fe2", "LC_dxt_r_fe2", 
                                            "LC_e0_r_fe2", "FDM_r_fe2", "CBD_SMM_r_fe2")], 
                     by = list(name = int_forecasts_df$name[obs]), FUN = mean,)
rownames(mse_tab) <- NULL
mse_tab[,2:ncol(mse_tab)] <- sqrt(mse_tab[,2:ncol(mse_tab)])
stargazer(mse_tab, summary = FALSE, title = "Mean Squared Forecast error by country",
          label = "tab:country_mse", table.placement = "H")









"
End of script
"


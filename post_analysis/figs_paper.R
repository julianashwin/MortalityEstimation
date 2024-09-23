setwd("/Users/julianashwin/Documents/GitHub/MortalityEstimation/")
rm(list=ls())

library(ggplot2)
library(ggpubr)
library(tidyverse)
library(stargazer)
library(lfe)
library(janitor)


"
Define color scheme
"
col_scheme <- c("Australia" = "darkolivegreen4", 
                "Canada" = "pink",
                "France" = "blue3", "United Kingdom" = "darkgoldenrod4", 
                "Hong Kong" = "lightgoldenrod", 
                "Italy" = "forestgreen", 
                "Japan" = "red","New Zealand" = "black", "Russia" = "firebrick",
                "Sweden" = "yellow", "USA" = "cornflowerblue",
                "Other" = "gray")


col_scheme <- c("Average" = "black", 
                "Australia" = "darkolivegreen4", 
                "China" = "lightgoldenrod", 
                "Ethiopia" = "yellow", 
                "France" = "blue3", 
                "India" = "pink",
                "Japan" = "red", 
                "Mexico" = "darkgoldenrod4",
                "Nigeria" = "forestgreen", 
                "Russia" = "firebrick",
                "USA" = "cornflowerblue",
                "Other" = "gray")

model_cols <- c("Siler" = "purple", "Lee-Carter (demo)" = "red", "Lee-Carter (SMM)" = "deeppink", 
                "FDM" = "forestgreen", "CBD" = "darkorange1",
                "Lee-Carter (dt)" = "gold", "Lee-Carter (dxt)" = "green", 
                "Lee-Carter (e0)" = "blue")


keep_codes <- c("AUS", "BEL", "CAN", "DNK", "FRA", "ITA", "NLD", "NZL_NM", "NOR",
                "PRT", "RUS", "ESP", "SWE", "CHE", "GBR", "USA", "JPN", "DEU",
                "IND", "CHN", "IDN", "PAK", "NGA", "BRA", "BGD", "MEX", "ETH",
                "PHL", "EGY", "COG", "VNM", "IRN", "TUR", "THA")

high_income <- c("Australia", "Belgium", "Canada", "Denmark", "France", "Germany", "Italy", "Japan", "Netherlands", "New Zealand", 
                 "Norway", "Portugal", "Russia", "Spain", "Sweden", "Switzerland", "UK", "USA")

other_income <- c("Bangladesh", "Belarus", "Brazil", "China", "Congo", "Egypt", "Ethiopia", "India", "Indonesia", 
                  "Iran", "Mexico", "Nigeria", "Pakistan", "Philippines", "Thailand", "Türkiye", "Viet Nam")

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
  mutate(DeltaLE_b = Δb*LE_b, DeltaLE_B = ΔB*LE_B, DeltaLE_d = Δd*LE_d,
         DeltaLE_c = Δc*LE_c, DeltaLE_C = ΔC*LE_C,
         year_group = case_when(
           year <= 1943 ~ "1918-1943",
           year > 1943 & year <= 1968 ~ "1943-1968",
           year > 1968 & year <= 1993 ~ "1968-1993",
           year > 1993 & year <= 2018 ~ "1993-2018",
           TRUE ~ "2018-2048"))
table(all_decomp_df$year_group)
rm(country_df)


## Import LE gradients for each country
import_files <- dir("figures/countries/")
import_files <- import_files[which(str_detect(import_files, "_LEgrads.csv"))]
all_LEgrads_df <- data.frame(matrix(NA,nrow=0,ncol = 40))
names(all_LEgrads_df) <- c("code", "age", "year", "LE", "LE_Bs", "LE_bs", "LE_Cs", "LE_cs", "LE_ds", "LE_cC", "Lstar", 
                           "Lstar_Bs", "Lstar_bs", "Lstar_Cs", "Lstar_cs", "Lstar_ds", "Lstar_cC", "Lmed", "Lmed_Bs", "Lmed_bs", "Lmed_Cs", 
                           "Lmed_cs", "Lmed_ds", "Lmed_cC", "h", "h_Bs", "h_bs", "h_Cs", "h_cs", "h_ds", "h_cC", 
                           "H", "H_Bs", "H_bs", "H_Cs", "H_cs", "H_ds", "H_cC", "mortality", "survival")
for (ii in 1:length(import_files)){
  filename <- import_files[ii]
  country_df <- read.csv(paste0("figures/countries/", filename), stringsAsFactors = F)
  country_df$code <- str_remove(filename, "_i2_decomp_pred.csv")
  country_df <- country_df[,names(all_LEgrads_df)]
  all_LEgrads_df <- rbind(all_LEgrads_df, country_df)
}

rm(country_df)


## Import panel data with econ variables
all_panel_df <- read_csv("data/clean/siler_econ_panel.csv")
all_panel_df <- as_tibble(all_panel_df) %>%
  select(-pop.x) %>%
  left_join(pop_df, by = c("code", "year")) %>%
  relocate(pop, .after = "code")

## Forecasts for BPLE 
forecasts_df <- read.csv("data/results/bp_forecasts.csv")

## Forecasts for countries 
int_forecasts_df <- read.csv("data/results/int_forecasts.csv")

####### Export data for Andrew to play with 
tib_list <- all_pars_df %>%
  filter(parameter %in% c("LE", "B", "b", "C", "c", "d", "h", "Lstar_99p9", "Lstar_99", "Lstar_95", "Lstar_90")) %>%
  mutate(parameter = str_replace(parameter, "B", "big_B"),
         parameter = str_replace(parameter, "b", "small_b"),
         parameter = str_replace(parameter, "C", "big_C"),
         parameter = str_replace(parameter, "c", "small_c")) %>%
  pivot_wider(id_cols = c(parameter, name), names_from = year, values_from = median) %>%
  split(f = as.factor(.$parameter))

library(openxlsx)
blank_excel <- createWorkbook()
Map(function(df, tab_name){     
  addWorksheet(blank_excel, tab_name)
  writeData(blank_excel, tab_name, df)
}, 
tib_list, names(tib_list)
)
saveWorkbook(blank_excel, file = "data/Andrew_siler_data.xlsx", overwrite = TRUE)


"
Convergence correlation table
"
conv_corr_tab <- all_pars_df %>%
  filter(code %in% keep_codes) %>%
  rbind(mutate(filter(all_pars_df, code %in% keep_codes), name = str_c(name, " (for all)"))) %>%
  mutate(name = case_when(name == "United States of America" ~ "USA", 
                             name == "Iran (Islamic Republic of)" ~ "Iran", 
                             name == "United Kingdom" ~ "UK", 
                             TRUE ~ name)) %>%
  mutate(income = case_when(name %in% high_income ~ "High~Income", 
                            name %in% other_income ~ "Large~Emerging~Economies",
                            TRUE ~ "All")) %>%
  mutate(income = factor(income, levels = c("High~Income", "Large~Emerging~Economies", "All"))) %>%
  filter(parameter %in% c("b", "B", "c", "C", "d", "LE", "h"),
         year %in% c(1948, 2018)) %>%
  mutate(parameter = factor(parameter, levels = c("b", "B", "c", "C", "d", "LE", "h"), ordered = T)) %>%
  group_by(income, name, parameter) %>%
  mutate(value = median, 
         initial_value = lag(median, n = 1, order_by = year),
         change = value - initial_value) %>%
  select(income, name, parameter, year, value, initial_value, change) %>%
  filter(year == 2018) %>%
  group_by(income, parameter) %>%
  summarise(cor = cor.test(initial_value, change)$estimate, 
            pval = cor.test(initial_value, change)$p.value)


conv_corr_tab %>%
  mutate(coef_text = sprintf(cor, fmt = '%#.2f')) %>%
  mutate(coef_text = case_when(pval <= 0.01 ~ str_c(coef_text, "***"),
                               pval <= 0.05 ~ str_c(coef_text, "**"),
                               pval <= 0.1 ~ str_c(coef_text, "*"),
                               TRUE ~ str_c(coef_text, ""))) %>%
  ggplot(aes(y =fct_rev(income), x = 1)) + theme_bw() + 
  facet_wrap(~ parameter, nrow = 1) +
  geom_tile(aes(fill = as.numeric(cor)), color = "black") +
  geom_text(aes(label = coef_text), color = "black", size = 3) +
  scale_fill_gradient2(low = "firebrick", high = "forestgreen", na.value = "white", 
                       mid = "white", midpoint = 0) +  #, midpoint = 0)
  theme(axis.text.x = element_blank(), axis.ticks = element_blank()) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  theme(legend.position="none") +
  labs(x = "", y = " ", fill = "Correlation of initial value and change \nfrom 1950 to present") 
ggsave("figures_paper/con_corr_table.pdf", width = 7, height = 1.5)

"
Mortality derivatives
"

mort_df %>%
  tibble() %>%
  arrange(code, year, age) %>%
  filter(code %in% c("USA", "FRA", "JPN", "GBR", "NGA", "CHN")) %>% 
  #filter(year %in% c(1948, 1968, 1988,2008, 2028, 2048)) %>%
  group_by(code, year) %>%
  filter(age != 100) %>%
  mutate(age = case_when(age == 99 ~ 100, TRUE ~ age)) %>%
  filter(age %in% c(60, 70, 80, 90, 100)) %>%
  mutate(Forecast = case_when(year < 2020 ~ "Estimate",
                              year >= 2020 ~ "Forecast")) %>%
  filter(age <= 110, year <= 2020, year >= 1950) %>%
  group_by(code, age) %>%
  mutate(mx_first = first(mx_f), 
         mort_growth = (mx_f) - (lag(mx_f, n = 1)),
         mx_norm = mx_f - mx_first) %>%
  ungroup() %>%
  ggplot() + theme_bw() + 
  facet_wrap(~name) +
  geom_smooth(aes(x = year, y = mx_norm, color = age, group = age), se = F) + 
  geom_point(aes(x = year, y = mx_norm, color = age, group = age), alpha = 0.2, size = 1) + 
  geom_hline(yintercept=0, linetype="dashed", color = "black") +
  scale_color_gradientn(colours = rainbow(5), name = "Age") +
  labs(x = "Year", y = "Change in mortality rate relative to 1950")



mort_df %>%
  tibble() %>%
  arrange(code, year, age) %>%
  filter(code %in% c("USA", "FRA", "JPN", "GBR", "NGA", "CHN")) %>% 
  #filter(year %in% c(1948, 1968, 1988,2008, 2028, 2048)) %>%
  group_by(code, year) %>%
  filter(age != 100) %>%
  mutate(age = case_when(age == 99 ~ 100, TRUE ~ age)) %>%
  filter(age %in% c(60, 70, 80, 90, 100)) %>%
  mutate(Forecast = case_when(year < 2020 ~ "Estimate",
                              year >= 2020 ~ "Forecast")) %>%
  filter(age <= 110, year <= 2020, year >= 1950) %>%
  group_by(code, age) %>%
  mutate(mx_first = first(mx_f), 
         mort_growth = (mx_f) - (lag(mx_f, n = 1)),
         mx_norm = mx_f - mx_first) %>%
  ungroup() %>%
  ggplot() + theme_bw() + 
  facet_wrap(~name) +
  geom_smooth(aes(x = year, y = mx_f, color = age, group = age), se = F) + 
  geom_point(aes(x = year, y = mx_f, color = age, group = age), alpha = 0.2, size = 1) + 
  geom_hline(yintercept=0, linetype="dashed", color = "black") +
  scale_color_gradientn(colours = rainbow(5), name = "Age") +
  labs(x = "Year", y = "Mortality rate")

bp_LEgrad_df %>%
  tibble() %>%
  arrange(age, year) %>%
  filter(year %in% c(1948, 1968, 1988,2008, 2028, 2048)) %>%
  left_join(select(bp_decomp_df, year, B, b, C, c, d)) %>%
  mutate(mu_C = -c^2 * exp(c*(age-C))) %>%
  select(age, year,  B, b, C, c, d, mu_C) %>%
  mutate(Forecast = case_when(year < 2020 ~ "Estimate",
                              year >= 2020 ~ "Forecast")) %>%
  filter(age <= 110) %>%
  group_by(age) %>%
  mutate(mu_C_growth = log(mu_C/lag(mu_C, n = 1))) %>%
  ggplot() + theme_bw() + 
  geom_line(aes(x = age, y = mu_C, color = year, group = year, linetype = Forecast)) + 
  geom_hline(yintercept=0, linetype="dashed", color = "black") +
  scale_color_gradientn(colours = rainbow(5), name = "Year",
                        breaks=c(1900, 1950, 2000,2049),
                        labels=c(1900,1950,2000,2049),
                        limits=c(1900,2049)) +
  labs(x = "Age", y = "Gradient")



"
Figure 1: Survival rate examples
"
survc_example <- read_csv("figures/interpret/survival_c_example.csv")
survc_example %>%
  ggplot(aes(x = age)) + theme_bw() +
  geom_line(aes(y = S_base, color = "Baseline"), size = 1) + 
  geom_line(aes(y = S_cchange, color = "Change to c"), linetype = "dashed", size = 1) +
  geom_line(aes(y = S_Cchange, color = "Change to C"), linetype = "dashed", size = 1) + 
  scale_color_manual(values = c("blue3", "firebrick", "forestgreen")) +
  labs(x = "Age", y = "Survival Rate", color = "Parameters")
ggsave("figures_paper/survival_example.pdf", width = 5, height = 3)
rm(survc_example)

"
Figure 2: Best practice mortality and life expectancy (1900-2019)
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


bp_df %>%
  filter(year > 1900 & age == 0) %>%
  ggplot() + theme_bw() +
  geom_smooth(aes(x = year, y = l_90_f), method = "loess") +
  geom_point(aes(x = year, y = l_90_f, color = name)) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  scale_color_discrete(name = "Country") + ggtitle("UBAD (90%)") +
  xlab("Year") + ylab("UBAD (90%)")

"
Figure 3: Siler parameter estimates for best practice mortality (1900-2019)
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
Figure 4: Historical decomposition of changes (1900-2019)
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
Figure 4: Prop LE gains across countries
"
decomp_table <- all_decomp_df %>%
  filter(!is.na(DeltaLE_c)) %>%
  group_by(code, name, year_group) %>%
  summarise(DeltaLE_bB = sum(DeltaLE_b) + sum(DeltaLE_B), DeltaLE_d = sum(DeltaLE_d),
            DeltaLE_c = sum(DeltaLE_c), DeltaLE_C = sum(DeltaLE_C), ΔLE_mod = sum(ΔLE_mod)) %>%
  ungroup() %>%
  mutate(propLE_bB = DeltaLE_bB/ΔLE_mod,
         propLE_d = DeltaLE_d/ΔLE_mod,
         propLE_c = DeltaLE_c/ΔLE_mod,
         propLE_C = DeltaLE_C/ΔLE_mod) %>%
  pivot_wider(id_cols = c(code, name), names_from = year_group, 
              names_glue = "{.value}_{year_group}",
              values_from = c(propLE_bB, propLE_d, propLE_c, propLE_C)) %>%
  select(code, name, 
         `propLE_bB_1918-1943`, `propLE_bB_1943-1968`, `propLE_bB_1968-1993`, `propLE_bB_1993-2018`,
         #`propLE_B_1918-1943`, `propLE_B_1943-1968`, `propLE_B_1968-1993`, `propLE_B_1993-2018`,
         `propLE_d_1918-1943`, `propLE_d_1943-1968`, `propLE_d_1968-1993`, `propLE_d_1993-2018`,
         `propLE_c_1918-1943`, `propLE_c_1943-1968`, `propLE_c_1968-1993`, `propLE_c_1993-2018`,
         `propLE_C_1918-1943`, `propLE_C_1943-1968`, `propLE_C_1968-1993`, `propLE_C_1993-2018`) %>%
  left_join(all_decomp_df[which(all_decomp_df$year == 2018),c("name", "C")]) %>%
  filter(code %in% keep_codes) %>%
  rename(country = name) %>%
  mutate(country = case_when(country == "United States of America" ~ "USA", 
                             country == "Iran (Islamic Republic of)" ~ "Iran", 
                             country == "United Kingdom" ~ "UK", 
                             TRUE ~ country)) %>%
  mutate(income = case_when(country %in% high_income ~ "High~Income", 
                            country %in% other_income ~ "Large~Emerging~Economies"))

# Export table as tex
stargazer(as.matrix(decomp_table), table.placement = "H", column.sep.width = "2")
# Or plot as heatmap
decomp_table %>%
  rbind(mutate(decomp_table, country = "Average")) %>%
  dplyr::select(-code) %>%
  group_by(income, country) %>%
  summarise(across(where(is.numeric), ~ mean(.x, na.rm = T))) %>%
  pivot_longer(cols = -c(income, country)) %>%
  mutate(country = factor(country, levels = c(high_income, other_income, "Average"), ordered = T)) %>%
  arrange(country) %>%
  mutate(parameter = case_when(nchar(name) >= 8 ~ str_c("Prop~LE~gains~", str_sub(name, 8, 8), "[t]"),
                               TRUE ~ "Latest~C[t]")) %>%
  mutate(interval = str_replace(str_sub(str_replace(name, "B", ""), 10, 18), "\\.", "-")) %>%
  mutate(interval = case_when(interval == "" ~ "Latest", TRUE ~ interval),
         parameter = case_when(parameter == "Prop~LE~gains~[t]" ~ "Latest~C[t]", 
                               parameter == "Prop~LE~gains~b[t]" ~ "Prop~LE~gains~b[t]~and~B[t]", 
                               TRUE ~ parameter)) %>%
  mutate(parameter = factor(parameter, levels = c("Prop~LE~gains~b[t]~and~B[t]",
                                                  "Prop~LE~gains~d[t]", 
                                                  "Prop~LE~gains~c[t]", "Prop~LE~gains~C[t]",
                                                  "Latest~C[t]"))) %>%
  mutate(value_latest = sprintf(value, fmt = '%#.1f'),
         value_latest = case_when(name == "C" ~ value_latest, TRUE ~ "")) %>%
  mutate(value = sprintf(value, fmt = '%#.2f'),
         value = case_when(name == "C" ~ "", 
                           value == "NaN" ~ "",
                           TRUE ~ value)) %>%
  ggplot(aes(y =fct_rev(country), x = interval)) + theme_minimal() + 
  facet_grid(vars(income), vars(parameter), labeller = label_parsed,
             scales = "free", space="free") +
  theme(axis.text.x=element_text(angle=45,hjust=1)) +
  geom_tile(aes(fill = as.numeric(value)), color = "black") +
  geom_text(aes(label = value), color = "black", size = 3) +
  geom_text(aes(label = value_latest), color = "black", size = 3) +
  scale_fill_gradient2(low = "firebrick", high = "forestgreen", na.value = "white", 
                       mid = "white", midpoint = 0) +  #, midpoint = 0)
  labs(x = "Interval", y = "", fill = "Proportion \nLE gains \ndue to \nparameter") 

ggsave("figures_paper/decomp_table.pdf", width = 12.2, height = 8)
rm(decomp_table)




"
Figure xxx: Prop LE gains across countries including the forecast
"
decomp_table_preds <- all_decomp_df %>%
  filter(!is.na(DeltaLE_c)) %>%
  group_by(code, name, year_group) %>%
  summarise(DeltaLE_bB = sum(DeltaLE_b) + sum(DeltaLE_B), DeltaLE_d = sum(DeltaLE_d),
            DeltaLE_c = sum(DeltaLE_c), DeltaLE_C = sum(DeltaLE_C), ΔLE_mod = sum(ΔLE_mod)) %>%
  ungroup() %>%
  mutate(propLE_bB = DeltaLE_bB/ΔLE_mod,
         propLE_d = DeltaLE_d/ΔLE_mod,
         propLE_c = DeltaLE_c/ΔLE_mod,
         propLE_C = DeltaLE_C/ΔLE_mod) %>%
  pivot_wider(id_cols = c(code, name), names_from = year_group, 
              names_glue = "{.value}_{year_group}",
              values_from = c(propLE_bB, propLE_d, propLE_c, propLE_C)) %>%
  select(code, name, 
         `propLE_bB_1918-1943`, `propLE_bB_1943-1968`, `propLE_bB_1968-1993`, `propLE_bB_1993-2018`, `propLE_bB_2018-2048`,
         #`propLE_B_1918-1943`, `propLE_B_1943-1968`, `propLE_B_1968-1993`, `propLE_B_1993-2018`,
         `propLE_d_1918-1943`, `propLE_d_1943-1968`, `propLE_d_1968-1993`, `propLE_d_1993-2018`, `propLE_d_2018-2048`,
         `propLE_c_1918-1943`, `propLE_c_1943-1968`, `propLE_c_1968-1993`, `propLE_c_1993-2018`, `propLE_c_2018-2048`,
         `propLE_C_1918-1943`, `propLE_C_1943-1968`, `propLE_C_1968-1993`, `propLE_C_1993-2018`, `propLE_C_2018-2048` ) %>%
  left_join(all_decomp_df[which(all_decomp_df$year == 2018),c("name", "C")]) %>%
  filter(code %in% keep_codes) %>%
  rename(country = name) %>%
  mutate(country = case_when(country == "United States of America" ~ "USA", 
                             country == "Iran (Islamic Republic of)" ~ "Iran", 
                             country == "United Kingdom" ~ "UK", 
                             TRUE ~ country)) %>%
  mutate(income = case_when(country %in% high_income ~ "High~Income", 
                            country %in% other_income ~ "Large~Emerging~Economies"))

# Export table as tex
stargazer(as.matrix(decomp_table), table.placement = "H", column.sep.width = "2")
# Or plot as heatmap
decomp_table_preds %>%
  rbind(mutate(decomp_table_preds, country = "Average")) %>%
  dplyr::select(-code) %>%
  group_by(income, country) %>%
  summarise(across(where(is.numeric), ~ mean(.x, na.rm = T))) %>%
  pivot_longer(cols = -c(income, country)) %>%
  mutate(country = factor(country, levels = c(high_income, other_income, "Average"), ordered = T)) %>%
  arrange(country) %>%
  mutate(parameter = case_when(nchar(name) >= 8 ~ str_c("Prop~LE~gains~", str_sub(name, 8, 8), "[t]"),
                               TRUE ~ "Latest~C[t]")) %>%
  mutate(interval = str_replace(str_sub(str_replace(name, "B", ""), 10, 18), "\\.", "-")) %>%
  mutate(interval = case_when(interval == "" ~ "Latest", TRUE ~ interval),
         parameter = case_when(parameter == "Prop~LE~gains~[t]" ~ "Latest~C[t]", 
                               parameter == "Prop~LE~gains~b[t]" ~ "Prop~LE~gains~b[t]~and~B[t]", 
                               TRUE ~ parameter)) %>%
  mutate(parameter = factor(parameter, levels = c("Prop~LE~gains~b[t]~and~B[t]",
                                                  "Prop~LE~gains~d[t]", 
                                                  "Prop~LE~gains~c[t]", "Prop~LE~gains~C[t]",
                                                  "Latest~C[t]"))) %>%
  mutate(value_latest = sprintf(value, fmt = '%#.1f'),
         value_latest = case_when(name == "C" ~ value_latest, TRUE ~ "")) %>%
  mutate(value = sprintf(value, fmt = '%#.2f'),
         value = case_when(name == "C" ~ "", 
                           value == "NaN" ~ "",
                           TRUE ~ value)) %>%
  ggplot(aes(y =fct_rev(country), x = interval)) + theme_minimal() + 
  facet_grid(vars(income), vars(parameter), labeller = label_parsed,
             scales = "free", space="free") +
  theme(axis.text.x=element_text(angle=45,hjust=1)) +
  geom_tile(aes(fill = as.numeric(value)), color = "black") +
  geom_text(aes(label = value), color = "black", size = 3) +
  geom_text(aes(label = value_latest), color = "black", size = 3) +
  scale_fill_gradient2(low = "firebrick", high = "forestgreen", na.value = "white", 
                       mid = "white", midpoint = 0) +  #, midpoint = 0)
  labs(x = "Interval", y = "", fill = "Proportion \nLE gains \ndue to \nparameter") 

ggsave("figures_paper/decomp_table_preds.pdf", width = 12.2, height = 8)
rm(decomp_table_preds)





"
Figure xxx: Absolute LE gains across countries including the forecast
"
LEchange_table_preds <- all_decomp_df %>%
  filter(!is.na(DeltaLE_c)) %>%
  group_by(code, name, year_group) %>%
  summarise(DeltaLE_bB = sum(DeltaLE_b) + sum(DeltaLE_B), DeltaLE_d = sum(DeltaLE_d),
            DeltaLE_c = sum(DeltaLE_c), DeltaLE_C = sum(DeltaLE_C), ΔLE_mod = sum(ΔLE_mod)) %>%
  ungroup() %>%
  mutate(propLE_bB = DeltaLE_bB,
         propLE_d = DeltaLE_d,
         propLE_c = DeltaLE_c,
         propLE_C = DeltaLE_C) %>%
  pivot_wider(id_cols = c(code, name), names_from = year_group, 
              names_glue = "{.value}_{year_group}",
              values_from = c(propLE_bB, propLE_d, propLE_c, propLE_C)) %>%
  select(code, name, 
         `propLE_bB_1918-1943`, `propLE_bB_1943-1968`, `propLE_bB_1968-1993`, `propLE_bB_1993-2018`, `propLE_bB_2018-2048`,
         #`propLE_B_1918-1943`, `propLE_B_1943-1968`, `propLE_B_1968-1993`, `propLE_B_1993-2018`,
         `propLE_d_1918-1943`, `propLE_d_1943-1968`, `propLE_d_1968-1993`, `propLE_d_1993-2018`, `propLE_d_2018-2048`,
         `propLE_c_1918-1943`, `propLE_c_1943-1968`, `propLE_c_1968-1993`, `propLE_c_1993-2018`, `propLE_c_2018-2048`,
         `propLE_C_1918-1943`, `propLE_C_1943-1968`, `propLE_C_1968-1993`, `propLE_C_1993-2018`, `propLE_C_2018-2048` ) %>%
  filter(code %in% keep_codes) %>%
  rename(country = name) %>%
  mutate(country = case_when(country == "United States of America" ~ "USA", 
                             country == "Iran (Islamic Republic of)" ~ "Iran", 
                             country == "United Kingdom" ~ "UK", 
                             TRUE ~ country)) %>%
  mutate(income = case_when(country %in% high_income ~ "High~Income", 
                            country %in% other_income ~ "Large~Emerging~Economies"))

# Export table as tex
stargazer(as.matrix(decomp_table), table.placement = "H", column.sep.width = "2")
# Or plot as heatmap
LEchange_table_preds %>%
  rbind(mutate(LEchange_table_preds, country = "Average")) %>%
  dplyr::select(-code) %>%
  group_by(income, country) %>%
  summarise(across(where(is.numeric), ~ mean(.x, na.rm = T))) %>%
  pivot_longer(cols = -c(income, country)) %>%
  mutate(country = factor(country, levels = c(high_income, other_income, "Average"), ordered = T)) %>%
  arrange(country) %>%
  mutate(parameter = case_when(nchar(name) >= 8 ~ str_c("LE~gains~", str_sub(name, 8, 8), "[t]"))) %>%
  mutate(interval = str_replace(str_sub(str_replace(name, "B", ""), 10, 18), "\\.", "-")) %>%
  mutate(interval = case_when(interval == "" ~ "Latest", TRUE ~ interval),
         parameter = case_when(parameter == "LE~gains~[t]" ~ "Latest~C[t]", 
                               parameter == "LE~gains~b[t]" ~ "LE~gains~b[t]~and~B[t]", 
                               TRUE ~ parameter)) %>%
  mutate(parameter = factor(parameter, levels = c("LE~gains~b[t]~and~B[t]",
                                                  "LE~gains~d[t]", 
                                                  "LE~gains~c[t]", "LE~gains~C[t]"))) %>%
  mutate(value = sprintf(value, fmt = '%#.2f'),
         value = case_when(name == "C" ~ "", 
                           value == "NaN" ~ "",
                           TRUE ~ value)) %>%
  ggplot(aes(y =fct_rev(country), x = interval)) + theme_minimal() + 
  facet_grid(vars(income), vars(parameter), labeller = label_parsed,
             scales = "free", space="free") +
  theme(axis.text.x=element_text(angle=45,hjust=1)) +
  geom_tile(aes(fill = as.numeric(value)), color = "black") +
  geom_text(aes(label = value), color = "black", size = 3) +
  scale_fill_gradient2(low = "firebrick", high = "forestgreen", na.value = "white", 
                       mid = "white", midpoint = 0) +  #, midpoint = 0)
  labs(x = "Interval", y = "", fill = "Absolute \nLE gains \ndue to \nparameter") 

ggsave("figures_paper/abs_LEgains_table_preds.pdf", width = 12.2, height = 8)
rm(decomp_table_preds)





"
Figure 6: Evolution of $c$ and $C$ over time
"
plot_df <- data.frame(pivot_wider(all_pars_df, id_cols = c(name, code, year, Forecast), 
                                  names_from = parameter, values_from = median)) %>%
  filter(code %in% keep_codes) %>%
  mutate(name = case_when(name == "United States of America" ~ "USA", 
                          name == "Iran (Islamic Republic of)" ~ "Iran", 
                          name == "United Kingdom" ~ "UK", 
                          TRUE ~ name))
extra_obs <- plot_df[which((plot_df$year == 2018 & !(plot_df$code %in% c("NZL","RUS")))|
                             (plot_df$year == 2013 & plot_df$code %in% c("NZL","RUS"))),]
extra_obs$Forecast <- "Estimate"
plot_df <- rbind(plot_df, extra_obs)
plot_df$h <- -log(plot_df$H)
plot_df$plot_name <- plot_df$name
plot_df$plot_name[which(!(plot_df$name %in% names(col_scheme)))] <- "Other"
plot_df$plot_name <- factor(plot_df$plot_name, levels = names(col_scheme), ordered = TRUE)
# Include some years as labels
plot_df$year_label <- NA
plot_df$year_label[which(plot_df$year==1903)] <- 1903
plot_df$year_label[which(plot_df$year==1933)] <- 1933
plot_df$year_label[which(plot_df$year==1963)] <- 1963
plot_df$year_label[which(plot_df$year==1993)] <- 1993
plot_df$year_label[which(plot_df$year==2023)] <- 2023

plot_df <- plot_df %>%
  rbind(mutate(filter(plot_df, Forecast != "Forecast" | year != 2013), 
               name = "Average", plot_name = "Average", code = "Average")) %>%
  group_by(name, plot_name, code, Forecast, year, year_label) %>%
  summarise(across(where(is.numeric), ~ mean(.x, na.rm = T))) %>%
  mutate(plot_name = factor(plot_name, levels = names(col_scheme), ordered = T))

plot_df %>%
  mutate(code = factor(code, levels = c(keep_codes, "Average"), ordered = T)) %>%
  mutate(other_alph = case_when(plot_name == "Other" ~ "0", TRUE ~ "1"),
         line_size =  case_when(plot_name == "Average" ~ "2", TRUE ~ "1"),
         year_label = case_when(plot_name == "Average" ~ year_label, TRUE ~ NA_real_)) %>%
  ggplot() + theme_bw() + 
  scale_color_manual("Country", values = col_scheme) + 
  scale_alpha_manual(guide = "none", values = c("0" = 0.4, "1" = 1)) +
  scale_size_manual(guide = "none", values = c("1" = 1, "2" = 2)) +
  geom_path(aes(x = C, y = c, color = plot_name, linetype = Forecast, size = line_size,
                group = interaction(code, Forecast), alpha = other_alph)) + 
  #geom_text(aes(x = C, y = c, color = plot_name, label = year_label), show.legend=FALSE, size = 5) +
  xlab("Timing of aging (C)") + ylab("Speed of aging (c)") + 
  guides(color=guide_legend(ncol=2))
ggsave("figures_paper/Cc_paths.pdf", width = 10, height = 6)


### For the Appendix: each country on a separate panel
plot_df %>%
  ggplot() + theme_bw() + 
  facet_wrap(~name) +
  scale_color_manual("Country", values = col_scheme) + 
  geom_path(aes(x = C, y = c, linetype = Forecast, group = interaction(code, Forecast))) + 
  geom_text(aes(x = C, y = c, label = year_label), show.legend=FALSE, size = 2) +
  xlab("Timing of aging (C)") + ylab("Speed of aging (c)")
ggsave("figures_paper/Cc_paths_facets.pdf", width = 12, height = 10)

rm(plot_df, extra_obs)


"
Figure 7: Evolution of L*, C and LE over time
"
Lstar_df <- all_pars_df %>%
  tibble() %>%
  filter(code %in% keep_codes) %>%
  filter(str_detect(parameter, "Lstar"))
Lstar_df %>%
  pivot_wider(id_cols = c("code", "year"), names_from = "parameter", values_from = "median") %>%
  filter(year == 2018) %>%
  arrange(-Lstar_90) %>%
  mutate(Lstar_90 = sprintf(Lstar_90, fmt = '%#.2f'),
         Lstar_95 = sprintf(Lstar_95, fmt = '%#.2f'),
         Lstar_99 = sprintf(Lstar_99, fmt = '%#.2f'),
         Lstar_99p9 = sprintf(Lstar_99p9, fmt = '%#.2f'))

plot_df <- data.frame(pivot_wider(all_pars_df, id_cols = c(name, code, year, Forecast), 
                                  names_from = parameter, values_from = median)) %>%
  filter(code %in% keep_codes) %>%
  mutate(name = case_when(name == "United States of America" ~ "USA", 
                          name == "Iran (Islamic Republic of)" ~ "Iran", 
                          name == "United Kingdom" ~ "UK", 
                          TRUE ~ name))
extra_obs <- plot_df[which((plot_df$year == 2018 & !(plot_df$code %in% c("NZL","RUS")))|
                             (plot_df$year == 2013 & plot_df$code %in% c("NZL","RUS"))),]
extra_obs$Forecast <- "Estimate"
plot_df <- rbind(plot_df, extra_obs)
plot_df$h <- -log(plot_df$H)
plot_df$plot_name <- plot_df$name
plot_df$plot_name[which(!(plot_df$name %in% names(col_scheme)))] <- "Other"
plot_df$plot_name <- factor(plot_df$plot_name, levels = names(col_scheme), ordered = TRUE)

# Include some years as labels
plot_df$year_label <- NA
plot_df$year_label[which(plot_df$year==1903)] <- 1903
plot_df$year_label[which(plot_df$year==1933)] <- 1933
plot_df$year_label[which(plot_df$year==1963)] <- 1963
plot_df$year_label[which(plot_df$year==1993)] <- 1993
plot_df$year_label[which(plot_df$year==2023)] <- 2023

plot_df <- plot_df %>%
  rbind(mutate(filter(plot_df, Forecast != "Forecast" | year != 2013), 
               name = "Average", plot_name = "Average", code = "Average")) %>%
  group_by(name, plot_name, code, Forecast, year, year_label) %>%
  summarise(across(where(is.numeric), ~ mean(.x, na.rm = T))) %>%
  mutate(plot_name = factor(plot_name, levels = names(col_scheme), ordered = T))

p1 <- plot_df %>%
  mutate(code = factor(code, levels = c(keep_codes, "Average"), ordered = T)) %>%
  mutate(other_alph = case_when(plot_name == "Other" ~ "0", TRUE ~ "1"),
         line_size =  case_when(plot_name == "Average" ~ "2", TRUE ~ "1"),
         year_label = case_when(plot_name == "Average" ~ year_label, TRUE ~ NA_real_)) %>%
  ggplot() + theme_bw() + 
  scale_color_manual("Country", values = col_scheme) + 
  scale_alpha_manual(guide = "none", values = c("0" = 0.4, "1" = 1)) +
  scale_size_manual(guide = "none", values = c("1" = 1, "2" = 2)) +
  geom_path(aes(x = LE, y = Lstar_99p9, color = plot_name, linetype = Forecast, size = line_size,
                group = interaction(code, Forecast), alpha = other_alph)) + 
  #geom_text(aes(x = LE, y = Lstar_99p9, color = plot_name, label = year_label), show.legend=FALSE, size = 5) +
  xlab("Life expectancy at birth") + ylab("UBAD") + 
  guides(color=guide_legend(ncol=2)) + 
  scale_x_continuous(limits = c(46,93),breaks=seq(0,130,10)) + 
  scale_y_continuous(limits = c(96,118),breaks=seq(0,130,10)) 
p1

### For the Appendix: each country on a separate panel
plot_df %>%
  ggplot() + theme_bw() + 
  facet_wrap(~name) +
  scale_color_manual("Country", values = col_scheme) + 
  geom_path(aes(x = LE, y = Lstar_99p9, linetype = Forecast, group = interaction(code, Forecast))) + 
  geom_text(aes(x = LE, y = Lstar_99p9, label = year_label), show.legend=FALSE, size = 2) +
  xlab("Life expectancy at birth") + ylab("UBAD")
ggsave("figures_paper/LE_Lstar_paths_facets.pdf", width = 12, height = 10)


p2 <- plot_df %>%
  mutate(code = factor(code, levels = c(keep_codes, "Average"), ordered = T)) %>%
  mutate(other_alph = case_when(plot_name == "Other" ~ "0", TRUE ~ "1"),
         line_size =  case_when(plot_name == "Average" ~ "2", TRUE ~ "1"),
         year_label = case_when(plot_name == "Average" ~ year_label, TRUE ~ NA_real_)) %>%
  ggplot() + theme_bw() + 
  scale_color_manual("Country", values = col_scheme) + 
  scale_alpha_manual(guide = "none", values = c("0" = 0.4, "1" = 1)) +
  scale_size_manual(guide = "none", values = c("1" = 1, "2" = 2)) +
  geom_path(aes(x = LE, y = Lstar_99p9-C, color = plot_name, linetype = Forecast, size = line_size,
                group = interaction(code, Forecast), alpha = other_alph)) + 
  #geom_text(aes(x = LE, y = Lstar_99p9-C, color = plot_name, label = year_label), show.legend=FALSE, size = 2) +
  xlab("Life expectancy at birth") + ylab("UBAD - C") + 
  guides(color=guide_legend(ncol=2))
p2

### For the Appendix: each country on a separate panel
plot_df %>%
  ggplot() + theme_bw() + 
  facet_wrap(~name) +
  scale_color_manual("Country", values = col_scheme) + 
  geom_path(aes(x = LE, y = Lstar_99p9-C, linetype = Forecast, group = interaction(code, Forecast))) + 
  geom_text(aes(x = LE, y = Lstar_99p9-C, label = year_label), show.legend=FALSE, size = 2) +
  xlab("Life expectancy at birth") + ylab("UBAD - C")
ggsave("figures_paper/LE_Lstar_C_paths_facets.pdf", width = 12, height = 10)


ggarrange(p1, p2, nrow = 1, common.legend = TRUE, legend = "right")
ggsave("figures_paper/Lstar_C_LE_paths.pdf", width = 10, height = 4)

  

rm(p1, p2, plot_df, extra_obs)




"
Figure 8: Forecast LE for Best Practice case
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
Figure 9: Gradients
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
Figure 10: Evolution of $h(0)$ over time
"
plot_df <- data.frame(pivot_wider(all_pars_df, id_cols = c(name, code, year, Forecast), 
                                  names_from = parameter, values_from = median)) %>%
  filter(code %in% keep_codes) %>%
  mutate(name = case_when(name == "United States of America" ~ "USA", 
                          name == "Iran (Islamic Republic of)" ~ "Iran", 
                          name == "United Kingdom" ~ "UK", 
                          TRUE ~ name))
extra_obs <- plot_df[which((plot_df$year == 2018 & !(plot_df$code %in% c("NZL","RUS")))|
                             (plot_df$year == 2013 & plot_df$code %in% c("NZL","RUS"))),]
extra_obs$Forecast <- "Estimate"
plot_df <- rbind(plot_df, extra_obs)
plot_df$h <- -log(plot_df$H)
plot_df$plot_name <- plot_df$name
plot_df$plot_name[which(!(plot_df$name %in% names(col_scheme)))] <- "Other"
plot_df$plot_name <- factor(plot_df$plot_name, levels = names(col_scheme), ordered = TRUE)

# Include some years as labels
plot_df$year_label <- NA
plot_df$year_label[which(plot_df$year==1903)] <- 1903
plot_df$year_label[which(plot_df$year==1933)] <- 1933
plot_df$year_label[which(plot_df$year==1963)] <- 1963
plot_df$year_label[which(plot_df$year==1993)] <- 1993
plot_df$year_label[which(plot_df$year==2023)] <- 2023

plot_df <- plot_df %>%
  rbind(mutate(filter(plot_df, Forecast != "Forecast" | year != 2013), 
               name = "Average", plot_name = "Average", code = "Average")) %>%
  group_by(name, plot_name, code, Forecast, year, year_label) %>%
  summarise(across(where(is.numeric), ~ mean(.x, na.rm = T))) %>%
  mutate(plot_name = factor(plot_name, levels = names(col_scheme), ordered = T))

plot_df %>%
  mutate(code = factor(code, levels = c(keep_codes, "Average"), ordered = T)) %>%
  mutate(other_alph = case_when(plot_name == "Other" ~ "0", TRUE ~ "1"),
         line_size =  case_when(plot_name == "Average" ~ "2", TRUE ~ "1"),
         year_label = case_when(plot_name == "Average" ~ year_label, TRUE ~ NA_real_)) %>%
  ggplot() + theme_bw() + 
  scale_color_manual("Country", values = col_scheme) + 
  scale_alpha_manual(guide = "none", values = c("0" = 0.4, "1" = 1)) +
  scale_size_manual(guide = "none", values = c("1" = 1, "2" = 2)) +
  guides(color=guide_legend(ncol=2)) + 
  geom_line(aes(x = year, y = h, color = plot_name, linetype = Forecast, size = line_size,
                group = interaction(code, Forecast), alpha = other_alph)) +
  xlab("Year") + ylab("Lifespan equality at birth")
ggsave("figures_paper/h_international.pdf", width = 10, height = 4)

rm(plot_df, extra_obs)



plot_df %>%
  ggplot() + theme_bw() + 
  facet_wrap(~name) +
  scale_color_manual("Country", values = col_scheme) + 
  geom_line(aes(x = year, y = h, linetype = Forecast, group = interaction(code, Forecast))) +
  xlab("Year") + ylab("Lifespan equality at birth")
ggsave("figures_paper/h_international_facets.pdf", width = 12, height = 10)


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
Appendix C: Compare different parameterised mortality functions
"
param_gomp = list(C = 85, c = 0.09, d = 0.0004)
gompertz_example_plot <- tibble(ages = 0:109) %>%
  mutate(`Age independent (λ)` = param_gomp$d,
         `Later life (α & β)` = param_gomp$c*exp(param_gomp$c*(ages - param_gomp$C)),
         total = (`Age independent (λ)` + `Later life (α & β)`)) %>%
  pivot_longer(cols = c(`Age independent (λ)`, `Later life (α & β)`)) %>%
  mutate(name = factor(name, levels = c("Age independent (λ)", "Later life (α & β)"), ordered = T)) %>%
  ggplot(aes(x = ages)) + theme_bw() + 
  geom_line(aes(y = log(total)), color = "black") + 
  geom_line(aes(y = log(value), color = name), linetype = "dashed") +
  coord_cartesian(ylim = c(-10.2, 0)) +
  scale_color_manual(values = c("forestgreen", "firebrick")) + 
  labs(x = "Age", y = expression(log(mu(a))), color = "Parts", title = "Gompertz-Makeham") + 
  theme(legend.position = c(0.33, 0.82), axis.text.y=element_blank(), axis.ticks.y=element_blank())
#gompertz_example_plot

param_siler = list(B = 3.5, b = 1.5, C = 85, c = 0.09, d = 0.0004)
siler_example_plot <- tibble(ages = 0:109) %>%
  mutate(`Infant (b & B)` = param_siler$b*exp(-param_siler$b*(ages + param_siler$B)),
         `Age independent (d)` = param_siler$d,
         `Later life (c & C)` = param_siler$c*exp(param_siler$c*(ages - param_siler$C)),
         total = (`Infant (b & B)` + `Age independent (d)` + `Later life (c & C)`)) %>%
  pivot_longer(cols = c(`Infant (b & B)`, `Age independent (d)`, `Later life (c & C)`)) %>%
  mutate(name = factor(name, levels = c("Infant (b & B)", "Age independent (d)", 
                                        "Later life (c & C)"), ordered = T)) %>%
  ggplot(aes(x = ages)) + theme_bw() + 
  geom_line(aes(y = log(total)), color = "black") + 
  geom_line(aes(y = log(value), color = name), linetype = "dashed") +
  coord_cartesian(ylim = c(-10.2, 0)) +
  scale_color_manual(values = c("blue3", "forestgreen", "firebrick")) + 
  labs(x = "Age", y = expression(log(mu(a))), color = "Parts", title = "Siler") + 
  theme(legend.position = c(0.33, 0.82), axis.text.y=element_blank(), axis.ticks.y=element_blank())
#siler_example_plot

param_hp = list(A = 0.0009, B = 0.04, C = 0.1, D = 0.0005, E = 6.0,
                F = 25.0, G = 0.00009, H = 1.08)
hp_example_plot <- tibble(ages = 0:109) %>%
  mutate(`Infant (A, B & C)` = param_hp$A^((ages + param_hp$B)^param_hp$C) ,
         `Accident hump (D, E & F)` = param_hp$D*exp(-param_hp$E*(log(ages) - log(param_hp$F))^2) ,
         `Later life (G & H)` = param_hp$G*param_hp$H^(ages),
         total = (`Infant (A, B & C)` + `Accident hump (D, E & F)` + `Later life (G & H)`)) %>%
  pivot_longer(cols = c(`Infant (A, B & C)`, `Accident hump (D, E & F)`, `Later life (G & H)`)) %>%
  mutate(name = factor(name, levels = c("Infant (A, B & C)", "Accident hump (D, E & F)", 
                                        "Later life (G & H)"), ordered = T)) %>%
  ggplot(aes(x = ages)) + theme_bw() + 
  geom_line(aes(y = log(total/(1-total))), color = "black") + 
  geom_line(aes(y = log(value/(1-value)), color = name), linetype = "dashed") +
  coord_cartesian(ylim = c(-10.2, 0)) +
  scale_color_manual(values = c("blue3", "forestgreen", "firebrick")) + 
  labs(x = "Age", y = expression(log(mu(a))), color = "Parts", title = "Heligman-Pollard") + 
  theme(legend.position = c(0.33, 0.82), axis.text.y=element_blank(), axis.ticks.y=element_blank())
# Alternative y label : expression(paste("log", bgroup("(", frac(mu(a),1-mu(a)), ")")))
#hp_example_plot

ggarrange(gompertz_example_plot, siler_example_plot, hp_example_plot, nrow = 1)
ggsave("figures_paper/mort_function_comparison.pdf", width = 10, height = 4, device = cairo_pdf)



    



"
Appendix D: Figure D.1 Siler parameter drift estimates for best practice mortality
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
  geom_hline(aes(yintercept = 0), linetype = "dashed") + 
  geom_line(aes(x = year, y = median)) + 
  labs(x = "Year", y = "")
ggsave("figures_paper/siler_i2drift_drifts.pdf", width = 8, height = 4, device = cairo_pdf)

"
Appendix D: Table D.2 Posterior distributions of variance terms
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
Appendix D: Figure D.3 independent periods
"

parests_compare <- read_csv("figures/benchmark/siler_indep_params_ber.csv") %>%
  mutate(Model = "Independent", forecast = 0) %>%
  select(any_of(c(names(bp_params_df), "Model"))) %>%
  rbind(mutate(bp_params_df, Model = "Dynamic"))


parests_compare %>% 
  filter(year >1900 & year < 2020 & parameter %in% c("b", "B", "c", "C", "d", "σ")) %>%
  mutate(parameter = factor(parameter, levels =c("b", "B", "c", "C", "d", "σ"), ordered = T)) %>%
  group_by(parameter) %>%
  mutate(min_par = min(pc025), max_par = max(pc975)) %>%
  ungroup() %>%
  ggplot() + theme_bw() + facet_wrap(~Model+parameter, scales = "free_y", nrow = 2) + 
  geom_line(aes(x = year, y = median, color = Model)) + 
  geom_blank(aes(x = year, y = min_par)) + 
  geom_blank(aes(x = year, y = max_par)) + 
  geom_ribbon(aes(x = year, ymin=pc025, ymax=pc975, fill = Model), alpha = 0.2) +
  geom_ribbon(aes(x = year, ymin=pc15, ymax=pc85, fill = Model), alpha = 0.3) + 
  labs(x = "Year", y = "")
#ggsave("figures_paper/siler_comp_params_sep.pdf", width = 10, height = 5, device = cairo_pdf)



parests_compare %>% 
  filter(year >1900 & year < 2020 & parameter %in% c("b", "B", "c", "C", "d", "σ")) %>%
  mutate(parameter = factor(parameter, levels =c("b", "c", "d", "B", "C", "σ"), ordered = T)) %>%
  group_by(parameter) %>%
  mutate(min_par = min(pc025), max_par = max(pc975)) %>%
  ungroup() %>%
  ggplot() + theme_bw() + facet_wrap(~parameter, scales = "free_y", nrow = 2) + 
  geom_line(aes(x = year, y = median, color = Model)) + 
  geom_ribbon(aes(x = year, ymin=pc025, ymax=pc975, fill = Model), alpha = 0.2) +
  geom_ribbon(aes(x = year, ymin=pc15, ymax=pc85, fill = Model), alpha = 0.3) + 
  labs(x = "Year", y = "")
ggsave("figures_paper/siler_comp_params_overlay.pdf", width = 8, height = 4, device = cairo_pdf)



"
Appendix D: Normality assumption
"
library(rstatix)
library(ggpubr)
library(tseries)


errors_df <- tibble(forecasts_df) %>%
  filter(est_year == 2018, year <= 2018, !is.na(mx_f)) %>%
  select(age, year, mortality, mx_f) %>%
  left_join(select(bp_decomp_df, year, σ)) %>%
  mutate(lmort_error = log(mx_f) - log(mortality)) 

errors_df %>%
  group_by(year, σ) %>%
  shapiro_test(lmort_error) %>%
  ggplot() + 
  geom_histogram(aes(x = p))
  arrange(-p)
  
errors_df %>%
  left_join(select(bp_decomp_df, year, σ)) %>%
  group_by(year, σ) %>%
  summarise(p = jarque.bera.test(lmort_error)$p.value)  %>%
  ggplot() + 
  geom_point(aes(x = p, fill = year))

errors_df %>%
  left_join(select(bp_decomp_df, year, σ)) %>%
  group_by(year, σ) %>%
  summarise(p = ks.test(lmort_error, rnorm(100,mean(lmort_error), σ))$p.value) %>%
  arrange(-p)
    
errors_df %>%
  ggqqplot(x = "lmort_error", facet.by = "year")
ggsave("figures_paper/residual_qqplot.pdf", width = 12, height = 10, device = cairo_pdf)

errors_df %>%
  ggplot(aes(x = lmort_error)) + theme_bw() + 
  geom_density()


int_fore_full <- read_csv("data/results/int_forecasts_full.csv")
int_errors_df <- int_fore_full %>%
  filter(est_year == 2018, year <= 2018, !is.na(mx), !is.na(mx_siler)) %>%
  select(code, name, year, age, mx_siler, mx) %>%
  left_join(select(all_decomp_df, code, year, σ)) %>%
  mutate(lmort_error = log(mx) - log(mx_siler)) 

int_errors_df %>%
  group_by(year, name, σ) %>%
  shapiro_test(lmort_error) %>%
  ggplot() + 
  geom_histogram(aes(x = p))

errors_df %>%
  left_join(select(bp_decomp_df, year, σ)) %>%
  group_by(year, σ) %>%
  summarise(p = ks.test(lmort_error, rnorm(100,mean(lmort_error), σ))$p.value) %>%
  arrange(-p)

errors_df %>%
  ggqqplot(x = "lmort_error", facet.by = "year")
ggsave("figures_paper/residual_qqplot.pdf", width = 12, height = 10, device = cairo_pdf)

errors_df %>%
  ggplot(aes(x = lmort_error)) + theme_bw() + 
  geom_density()

mean(int_errors_df$lmort_error)/(sd(int_errors_df$lmort_error)/sqrt(nrow(int_errors_df)))
ks.test(int_errors_df$lmort_error, "pnorm")


"
Appendix D: 1 year periods
"


"
Appendix E: All country level parameters
"
plot_df <- data.frame(pivot_wider(all_pars_df, id_cols = c(name, code, year, Forecast), 
                                  names_from = parameter, values_from = median)) %>%
  filter(code %in% keep_codes) %>%
  mutate(name = case_when(name == "United States of America" ~ "USA", 
                          name == "Iran (Islamic Republic of)" ~ "Iran", 
                          name == "United Kingdom" ~ "UK", 
                          TRUE ~ name))
extra_obs <- plot_df[which((plot_df$year == 2018 & !(plot_df$code %in% c("NZL","RUS")))|
                             (plot_df$year == 2013 & plot_df$code %in% c("NZL","RUS"))),]
extra_obs$Forecast <- "Estimate"
plot_df <- rbind(plot_df, extra_obs)
plot_df$h <- -log(plot_df$H)
plot_df$plot_name <- plot_df$name
plot_df$plot_name[which(!(plot_df$name %in% names(col_scheme)))] <- "Other"
plot_df$plot_name <- factor(plot_df$plot_name, levels = names(col_scheme), ordered = TRUE)

plot_df <- plot_df %>%
  rbind(mutate(filter(plot_df, Forecast != "Forecast" | year != 2013), 
               name = "Average", plot_name = "Average", code = "Average")) %>%
  group_by(name, plot_name, code, Forecast, year) %>%
  summarise(across(where(is.numeric), ~ mean(.x, na.rm = T))) %>%
  mutate(plot_name = factor(plot_name, levels = names(col_scheme), ordered = T))


plot_df %>%
  rename(Country = name) %>%
  mutate(code = factor(code, levels = c(keep_codes, "Average"), ordered = T)) %>%
  mutate(other_alph = case_when(plot_name == "Other" ~ "0", TRUE ~ "1"),
         line_size =  case_when(plot_name == "Average" ~ "2", TRUE ~ "1")) %>%
  pivot_longer(cols = c(b, B, c, C, d, σ)) %>%
  ggplot() + theme_bw() + 
  facet_wrap(~ name, scales = 'free_y', nrow = 3) +
  scale_color_manual("Country", values = col_scheme) + 
  scale_alpha_manual(guide = "none", values = c("0" = 0.4, "1" = 1)) +
  scale_size_manual(guide = "none", values = c("1" = 1, "2" = 2)) +
  geom_line(aes(x = year, y = value, color = plot_name, linetype = Forecast, size = line_size,
                group = interaction(code, Forecast), alpha = other_alph)) + 
  xlab("Year") + ylab("Values")
ggsave("figures_paper/country_params.pdf", width = 8, height = 10, device = cairo_pdf)


plot_df %>%
  rename(Country = name) %>%
  mutate(code = factor(code, levels = c(keep_codes, "Average"), ordered = T)) %>%
  mutate(other_alph = case_when(plot_name == "Other" ~ "0", TRUE ~ "1"),
         line_size =  case_when(plot_name == "Average" ~ "2", TRUE ~ "1")) %>%
  pivot_longer(cols = c(α_b, α_B, α_c, α_C, α_d, α_σ)) %>%
  mutate(name = str_c(str_replace_all(name, "_", "["), "]")) %>%
  ggplot() + theme_bw() + 
  facet_wrap(~ name, scales = 'free_y', nrow = 3, labeller = label_parsed) +
  scale_color_manual("Country", values = col_scheme) + 
  scale_alpha_manual(guide = "none", values = c("0" = 0.4, "1" = 1)) +
  scale_size_manual(guide = "none", values = c("1" = 1, "2" = 2)) +
  geom_line(aes(x = year, y = value, color = plot_name, linetype = Forecast, size = line_size,
                group = interaction(code, Forecast), alpha = other_alph)) + 
  xlab("Year") + ylab("Values")
ggsave("figures_paper/country_drifts.pdf", width = 8, height = 10, device = cairo_pdf)




tibble(plot_df) %>%
  group_by(code) %>%
  mutate(incl_1950 = max(year <= 1953)) %>%
  filter(incl_1950 == 1) %>%
  filter(year >= 1953) %>%
  rename(Country = name) %>%
  pivot_longer(cols = c(b, B, c, C, d, LE)) %>%
  group_by(name, year, Forecast) %>%
  arrange(name, year) %>%
  summarise(mean = mean(value), lmean = mean(log(value)), sd = sd(value), lsd = sd(log(value))) %>%
  ggplot() + theme_bw() + 
  facet_wrap(~ name, scales = 'free_y') +
  scale_color_manual("Country", values = col_scheme) + 
  geom_line(aes(x = year, y = sd, linetype = Forecast))



"
Appendix E: Table 3) Economic Determinants of Senescent Mortality
"


reg_df <- all_panel_df %>%
  select(code, name, year, c, C, gini_income, gdp_pc_usd, pop, health_exp_share) %>%
  drop_na(pop, gini_income) %>% 
  group_by(year) %>% 
  mutate(pop_prop = pop/sum(pop), US = as.numeric(code == "USA"))
reg_lag_df <- all_panel_df %>%
  select(code, name, year, c, C, gini_income_lag, gdp_pc_usd_lag, pop, health_exp_share_lag) %>%
  drop_na(pop, gini_income_lag) %>% group_by(year) %>% 
  mutate(pop_prop = pop/sum(pop), US = as.numeric(code == "USA"))


model1 <- (felm(c ~ gini_income + log(gdp_pc_usd) |
                  name + year, data = reg_df, weights = reg_df$pop_prop))
model2 <- (felm(C ~ gini_income + log(gdp_pc_usd) |
                  name + year, data = reg_df, weights = reg_df$pop_prop))
model3 <- (felm(c ~ gini_income + log(gdp_pc_usd) + health_exp_share|
                  name + year, data = reg_df, weights = reg_df$pop_prop))
model4 <- (felm(C ~ gini_income + log(gdp_pc_usd) + health_exp_share|
                  name + year, data = reg_df, weights = reg_df$pop_prop))
model5 <- (felm(c ~ gini_income_lag + log(gdp_pc_usd_lag) |
                  name + year, data = reg_lag_df, weights = reg_lag_df$pop_prop))
model6 <- (felm(C ~ gini_income_lag + log(gdp_pc_usd_lag) |
                  name + year, data = reg_lag_df, weights = reg_lag_df$pop_prop))
model7 <- (felm(c ~ gini_income_lag + log(gdp_pc_usd_lag) + health_exp_share_lag|
                  name + year, data = reg_lag_df, weights = reg_lag_df$pop_prop))
model8 <- (felm(C ~ gini_income_lag + log(gdp_pc_usd_lag) + health_exp_share_lag|
                  name + year, data = reg_lag_df, weights = reg_lag_df$pop_prop))

models <- list(model1, model2, model3, model4, model5, model6, model7, model8)
stargazer(models, table.placement = "H",
          df = FALSE, title = "C and c on economic variables", label = "tab:cC_econ",
          font.size = "scriptsize")

model1 <- felm(c ~ gini_income + log(gdp_pc_usd) |
                  name + year, data = reg_df)
model2 <- felm(C ~ gini_income + log(gdp_pc_usd) |
                  name + year, data = reg_df)
model3 <- felm(c ~ gini_income + log(gdp_pc_usd) + health_exp_share|
                  name + year, data = reg_df)
model4 <- felm(C ~ gini_income + log(gdp_pc_usd) + health_exp_share|
                  name + year, data = reg_df)
model5 <- felm(c ~ gini_income_lag + log(gdp_pc_usd_lag) + health_exp_share_lag|
                  name + year, data = reg_lag_df)
model6 <- felm(C ~ gini_income_lag + log(gdp_pc_usd_lag) + health_exp_share_lag|
                  name + year, data = reg_lag_df)
models <- list(model1, model2, model3, model4, model5, model6)
stargazer(models, table.placement = "H",
          df = FALSE, title = "C and c on economic variables", label = "tab:cC_econ",
          font.size = "scriptsize")


model1 <- (felm(c ~ gini_income + log(gdp_pc_usd) |
                  year, data = reg_df, weights = reg_df$pop_prop))
model2 <- (felm(C ~ gini_income + log(gdp_pc_usd) |
                  year, data = reg_df, weights = reg_df$pop_prop))
model3 <- (felm(c ~ gini_income + log(gdp_pc_usd) + health_exp_share|
                  year, data = reg_df, weights = reg_df$pop_prop))
model4 <- (felm(C ~ gini_income + log(gdp_pc_usd) + health_exp_share|
                  year, data = reg_df, weights = reg_df$pop_prop))
model5 <- (felm(c ~ gini_income_lag + log(gdp_pc_usd_lag) + health_exp_share_lag|
                  year, data = reg_lag_df, weights = reg_lag_df$pop_prop))
model6 <- (felm(C ~ gini_income_lag + log(gdp_pc_usd_lag) + health_exp_share_lag|
                  year, data = reg_lag_df, weights = reg_lag_df$pop_prop))
models <- list(model1, model2, model3, model4, model5, model6)
stargazer(models, table.placement = "H",
          df = FALSE, title = "C and c on economic variables", label = "tab:cC_econ_nocountryfe",
          font.size = "scriptsize")

rm(reg_df, reg_lag_df, models, model1, model2, model3, model4, model5, model6, model7, model8)





"
Appendix F: Comparing static and dynamic model forecasts
"



forecasts_dist_df <- read.csv("figures/benchmark/held-out/siler_i2drift_preds_all.csv", stringsAsFactors = F)


indep_pars <- parests_compare %>%
  filter(Model == "Independent") %>%
  select(Model, year, forecast, parameter, median, pc85, pc15) %>%
  full_join(crossing(Model = "Independent", est_year = unique(forecasts_df$est_year), 
                     year = unique(forecasts_df$year), parameter = c("B","b","C","c","d","σ"))) %>%
  filter(year - est_year <= 30) %>%
  mutate(median = case_when(year > est_year ~ NA_real_, TRUE ~ median)) %>%
  filter(est_year - year <= 50) %>%
  group_by(parameter, est_year) %>%
  mutate(intercept = lm(log(median) ~ year)$coefficients[1],
         slope = lm(log(median) ~ year)$coefficients[2],
         Forecast = case_when(year > est_year ~ "Forecast", TRUE ~ "Estimate"))  %>%
  mutate(median = case_when(Forecast == "Forecast" ~ exp(intercept + slope*year), TRUE ~ median)) %>%
  select(-intercept, -slope) %>%
  ungroup() %>% arrange(est_year, parameter, year)
indep_pars <- indep_pars %>% 
  rbind(mutate(filter(indep_pars, year == est_year), Forecast = "Forecast"))
  
  
indep_pars %>%
  ggplot() + theme_bw() + facet_wrap(~parameter, scales = "free_y", nrow = 2) + 
  geom_ribbon(aes(x = year, ymin=pc15, ymax=pc85), alpha = 0.3) + 
  geom_line(aes(x = year, y = median, color = est_year, 
                group = interaction(est_year, Forecast), linetype = Forecast)) + 
  labs(x = "Year", y = "")
#ggsave("figures_paper/indep_param_forecasts.pdf", width = 8, height = 4, device = cairo_pdf)





forecasts_dist_df %>%
  mutate(Forecast = case_when(year > est_year ~ "Forecast", TRUE ~ "Estimate")) %>%
  rbind(mutate(filter(forecasts_dist_df, year == est_year), Forecast = "Forecast")) %>%
  mutate(Model = "Dynamic") %>%
  filter(parameter %in% c("B","b","C","c","d","σ")) %>%
  select(Model, year, forecast, parameter, median, pc85, pc15, est_year, Forecast) %>%
  rbind(indep_pars) %>%
  ggplot() + theme_bw() + facet_wrap(~Model+parameter, scales = "free_y", nrow = 2) + 
  geom_ribbon(aes(x = year, ymin=pc15, ymax=pc85, group = est_year), alpha = 0.3) + 
  geom_line(aes(x = year, y = median, color = Model, 
                group = interaction(est_year, Forecast), linetype = Forecast)) + 
  labs(x = "Year", y = "") + 
  theme(legend.position = "bottom")
ggsave("figures_paper/indep_param_forecasts.pdf", width = 10, height = 6, device = cairo_pdf)



indep_forecast_df <- indep_pars %>%
  pivot_wider(id_cols = c(est_year, year, Forecast), names_from = parameter, values_from = median) %>%
  crossing(age = 0:100) %>%
  mutate(Sa = exp(-d*age + (exp(-b*(age + B)) - exp(-b*B)) - (exp(c*(age - C)) - exp(-c*C)))) %>%
  group_by(est_year, year, Forecast) %>%
  summarise(age = min(age), b = mean(b), B = mean(B), c = mean(c), C = mean(C), d = mean(d), σ = mean(σ), LE = sum(Sa))
  
indep_forecast_df %>%
  mutate(n_ahead = case_when(year > est_year ~ year - est_year, TRUE ~ NA_real_)) %>%
  left_join(select(forecasts_df, est_year, year, age, ex_f, ex_siler, siler_fe, siler_fe2)) %>%
  filter(n_ahead > 0 ) %>%
  mutate(siler_indep_fe = ex_f - LE, 
         siler_indep_fe2 = siler_indep_fe^2) %>%
  ggplot() + theme_bw() + 
  geom_bar(aes(x = n_ahead-0.25, y = (siler_indep_fe2), fill = "Independent"), 
           stat = "summary", fun = mean, width = 0.5) + 
  geom_bar(aes(x = n_ahead+0.25, y = (siler_fe2), fill = "Dynamic"), 
           stat = "summary", fun = mean, width = 0.5) +
  labs(x = "Years ahead", y = "RMSE", fill = "Model")
ggsave("figures_paper/indep_param_rmse.pdf", width =6, height = 4, device = cairo_pdf)




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
Appendix F: Table 4) Increases in LE (years per decade) 
"
increase_df <- all_pars_df %>%
  filter(code %in% keep_codes) %>%
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
         `2038-2048` = round(`2048`-`2038`, 3)#,
         #`2018-2048` = (`2048`-`2018`)/3,
         #`1988-2018` = (`2018`-`1988`)/3,
         #trend_change = (`2048`-`2018`)/3 - (`2018`-`1988`)/3
         ) %>%
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


increase_df %>%
  mutate(income = case_when(name %in% high_income ~ "High~Income", 
                            name %in% other_income ~ "Large~Emerging~Economies")) %>%
  group_by(income, Period, Change_bucket) %>%
  summarise(value = n()) %>%
  ungroup() %>%
  complete(income, Period, Change_bucket) %>%
  replace_na(list(value = 0)) %>%
  ggplot(aes(y =fct_rev(Period), x = Change_bucket)) + theme_bw() + 
  facet_wrap(vars(income), labeller = label_parsed, nrow = 2) +
  theme(axis.text.x=element_text(angle=45,hjust=1)) +
  geom_tile(aes(fill = as.numeric(value)), color = "black") +
  geom_text(aes(label = value), color = "black", size = 3) +
  scale_fill_gradient2(low = "firebrick", high = "forestgreen", na.value = "white", 
                       mid = "white", midpoint = 0, guide = "none") +  #, midpoint = 0)
  labs(x = "Increase in LE (years per decade)", y = "Period") 
ggsave("figures_paper/LE_increase_table.pdf", width = 6, height = 4, device = cairo_pdf)

rm(increase_df, increase_table, increase_tab_df)







"
Appendix F: forecast errors
"

# Plot the errors
fe_mean_plt <-
  tibble(int_forecasts_df) %>%
  filter(code %in% keep_codes) %>%
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
  filter(code %in% keep_codes) %>%
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




### Now do the same but split by income 
fe_mean_plt <-
  tibble(int_forecasts_df) %>%
  filter(code %in% keep_codes) %>%
  left_join(select(all_decomp_df, code, name)) %>%
  mutate(income = case_when(name %in% high_income ~ "High~Income", 
                            name %in% other_income ~ "Large~Emerging~Economies")) %>%
  filter(age == 0 & n_ahead > 0) %>%
  ggplot() + theme_bw() + 
  facet_wrap(~income, ncol = 1) +
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
  filter(code %in% keep_codes) %>%
  left_join(select(all_decomp_df, code, name)) %>%
  mutate(income = case_when(name %in% high_income ~ "High~Income", 
                            name %in% other_income ~ "Large~Emerging~Economies")) %>%
  filter(age == 0 & n_ahead > 0) %>%
  ggplot() + theme_bw() + 
  facet_wrap(~income, ncol = 1) +
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
ggsave("figures_paper/country_incgrps_forecast_errors.pdf", width = 8, height = 6.5, device = cairo_pdf)




"
Appendix F: Table 5)
"

mse_tab <- int_forecasts_df %>%
  filter(code %in% keep_codes) %>%
  mutate(name = case_when(name == "United States of America" ~ "USA", 
                          name == "Iran (Islamic Republic of)" ~ "Iran", 
                          name == "United Kingdom" ~ "UK", 
                          TRUE ~ name)) %>%
  filter(!is.na(siler_fe) & age == 0 & est_year > 1960) %>%
  mutate(income = case_when(name %in% high_income ~ "High~Income", 
                            name %in% other_income ~ "Large~Emerging~Economies")) %>%
  group_by(income, name) %>%
  summarise(Siler = sqrt(mean(siler_fe2)), `LC (demo)` = sqrt(mean(LC_r_fe2)), 
            `LC (SMM)` = sqrt(mean(LC_SMM_r_fe2)), `LC (dt)` = sqrt(mean(LC_dt_r_fe2)),
            `LC (dxt)` = sqrt(mean(LC_dxt_r_fe2)), `LC (e0)` = sqrt(mean(LC_e0_r_fe2)), 
            `FDM` = sqrt(mean(FDM_r_fe2)), `CBD` = sqrt(mean(CBD_SMM_r_fe2))) %>%
  arrange(income, name)

mse_tab %>%
  rename(`Country` = name) %>%
  mutate(Country = case_when(Country == "United States of America" ~ "USA", 
                             Country == "Iran (Islamic Republic of)" ~ "Iran", 
                             Country == "United Kingdom" ~ "UK", 
                             TRUE ~ Country)) %>%
  pivot_longer(cols = c(-income, -Country)) %>%
  mutate(`Model` = factor(name, levels = c("Siler", "LC (demo)", "LC (SMM)", "LC (dt)", "LC (dxt)", "LC (e0)", "FDM", "CBD"), ordered = T)) %>%
  group_by(Country) %>%
  mutate(lowest = as.character(value == min(value))) %>%
  mutate(value_text =  sprintf(value, fmt = '%#.2f')) %>%
  ggplot(aes(y =fct_rev(Country), x = Model)) + theme_minimal() + 
  facet_grid(vars(income), labeller = label_parsed,
             scales = "free", space="free") +
  theme(axis.text.x=element_text(angle=45,hjust=1)) +
  geom_tile(aes(fill = lowest), color = "black") +
  geom_text(aes(label = value_text), color = "black", size = 3) +
  scale_fill_manual(guide = "none", values = c("FALSE" = "white", "TRUE" = "firebrick")) +  #, midpoint = 0)
  labs(x = "Model", y = "Country") 
ggsave("figures_paper/rmse_by_country.pdf", width = 5, height = 7, device = cairo_pdf)


ggplot(aes(y =fct_rev(country), x = interval)) + theme_minimal() + 
  facet_grid(vars(income), vars(parameter), labeller = label_parsed,
             scales = "free", space="free") +
  theme(axis.text.x=element_text(angle=45,hjust=1)) +
  geom_tile(aes(fill = as.numeric(value)), color = "black") +
  geom_text(aes(label = value), color = "black", size = 3) +
  geom_text(aes(label = value_latest), color = "black", size = 3) +
  scale_fill_gradient2(low = "firebrick", high = "forestgreen", na.value = "white", 
                       mid = "white", midpoint = 0) +  #, midpoint = 0)
  labs(x = "Interval", y = "", fill = "Proportion \nLE gains \ndue to \nparameter") 


"
End of script
"


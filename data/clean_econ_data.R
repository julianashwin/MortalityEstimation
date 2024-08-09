setwd("/Users/julianashwin/Documents/GitHub/MortalityEstimation/")
rm(list=ls())


require(tidyverse)
require(ggpubr)
require(lfe)
require(stargazer)
require(readxl)
require(janitor)

"
Set colour scheme
"
col_scheme <- c("Other" = "gray",
                "Australia" = "darkolivegreen4", 
                "Canada" = "pink",
                "France" = "blue3", "United Kingdom" = "darkgoldenrod4", 
                "Hong Kong" = "lightgoldenrod", 
                "Italy" = "forestgreen", 
                "Japan" = "red","New Zealand" = "black", "Russia" = "firebrick",
                "Sweden" = "yellow", "United States of America" = "cornflowerblue",
                "Best Practice" = "darkmagenta")

col_scheme <- c("Australia" = "darkolivegreen4", 
                "India" = "pink",
                "France" = "blue3", "Mexico" = "darkgoldenrod4", 
                "Hong Kong" = "lightgoldenrod", 
                "Nigeria" = "forestgreen", 
                "Japan" = "red","China" = "black", "Russia" = "firebrick",
                "Ethiopia" = "yellow", "United States of America" = "cornflowerblue",
                "Other" = "gray")


keep_codes <- c("AUS", "BEL", "CAN", "DNK", "FRA", "ITA", "NLD", "NZL_NM", "NOR",
                "PRT", "RUS", "ESP", "SWE", "CHE", "GBR", "USA", "JPN", "DEU",
                "IND", "CHN", "IDN", "PAK", "NGA", "BRA", "BGD", "MEX", "ETH",
                "PHL", "EGY", "COG", "VNM", "IRN", "TUR", "THA")

un_data_codes <- c("IND", "CHN", "IDN", "PAK", "NGA", "BRA", "BGD", "MEX", "ETH",
                   "PHL", "EGY", "COG", "VNM", "IRN", "TUR", "THA")

#mort_df <- read.csv("data/clean/all_lifetab_1y.csv", stringsAsFactors = FALSE)


siler_df <- read.csv("data/results/siler_panel.csv", stringsAsFactors = FALSE)
siler_df$plot_name <- siler_df$name
siler_df$plot_name[which(!(siler_df$code %in% keep_codes))] <- "Other"

siler_df %>%
  filter(plot_name != "Other") %>%
  ggplot() + theme_bw() + 
  #scale_color_manual("Country", values = col_scheme) + guides(color=guide_legend(ncol=2)) + 
  geom_line(aes(x = C, y = Lstar_90, color = plot_name, group = code))



"
data/national_indicators.csv
"
wb_df <- read_csv("data/national_indicators.csv") %>% 
  mutate(code = str_replace(code, "NZL", "NZL_NM")) %>%
  filter(code %in% keep_codes,
         year < 2020 & year > 1900) %>%
  mutate(years = as.character(cut(year, seq(1900, 2020, 5)) )) %>%
  relocate(years, .after = year)


"
Table for Data Appendix
"
mort_sum <- siler_df %>% 
  filter(code %in% keep_codes) %>%
  arrange(code, name, year) %>%
  group_by(name, code) %>%
  mutate(forecast = case_when(year < 2020 ~ 0,
                              TRUE ~ 1)) %>%
  filter(forecast == 0) %>%
  summarise(mort_start = min(year), mort_end = max(year)) %>%
  mutate(mort_start = paste0("{",mort_start - 2, ",", mort_start +2,"}")) %>%
  mutate(mort_end = paste0("{",mort_end - 2, ",", mort_end +2,"}")) %>%
  mutate(mort_sample = paste0(mort_start,":",mort_end),
         mort_sample = str_replace(mort_sample, "2020", "2019")) %>%
  mutate(source = case_when(code %in% un_data_codes ~ "UN WPP", 
         TRUE ~ "HMD")) %>%
  select(name, code, mort_sample, source)


gdp_sum <- wb_df %>%
  group_by(code, region) %>%
  filter(!is.na(gdp_pc_usd) & year < 2019) %>%
  summarise(gdp_start = min(year), gdp_end = max(year)) %>%
  mutate(gdp_sample = paste0(gdp_start,":",gdp_end)) %>%
  select(code, gdp_sample)

gini_sum <- wb_df %>%
  group_by(code, region) %>%
  filter(!is.na(gini_income) & year < 2019) %>%
  summarise(gini_start = min(year), gini_end = max(year)) %>%
  mutate(gini_sample = paste0(gini_start,":",gini_end)) %>%
  select(code, gini_sample)

health_sum <- wb_df %>%
  group_by(code, region) %>%
  filter(!is.na(health_exp_share) & year < 2019) %>%
  summarise(health_share_start = min(year), health_share_end = max(year)) %>%
  mutate(health_share_sample = paste0(health_share_start,":",health_share_end)) %>%
  select(code, health_share_sample)


summary_df <- mort_sum %>%
  full_join(gdp_sum) %>%
  full_join(gini_sum) %>%
  full_join(health_sum) %>%
  select(-code) %>%
  rename(Country = name, Mortality = mort_sample, `Mortality Source` = source, 
         GDP = gdp_sample, `Income Inequality` = gini_sample,
         `Health Share` = health_share_sample) 
stargazer(as.matrix(summary_df), table.placement = "H", label = "tab:sample_starts",
          title = "Data availability by country")



"
Merge together and aggregate
"
econ_5y_df <- aggregate(subset(wb_df, select=-c(region,code,years)), 
                        FUN = mean, na.rm = TRUE, na.action = na.omit,
                        by = list(region = wb_df$region, code = wb_df$code, 
                                  years = wb_df$years))
econ_5y_df <- tibble(econ_5y_df) %>%
  arrange(code, year) %>% 
  mutate(year = as.integer(case_when(year == 2017.5 ~ 2018, TRUE ~ year)))


# Test plot
ggplot(econ_5y_df) + theme_bw() + 
  scale_color_manual("Country", values = col_scheme) + guides(color=guide_legend(ncol=2)) + 
  geom_point(aes(x = health_exp_share, y = health_exp_pc_ppp, color = region, group = region))


all_panel_df <- tibble(siler_df) %>%
  left_join(econ_5y_df, by = c("code", "year")) %>%
  arrange(year) %>%
  mutate(year_num = as.numeric(as.factor(year))) %>%
  relocate(year_num, .after = year) %>%
  arrange(code, year) %>%
  group_by(code) %>%
  mutate(C_lag = lag(C),
         C_diff = C - C_lag,
         c_lag = lag(c),
         c_diff = c - c_lag,
         gdp_pc_usd_lag = lag(gdp_pc_usd),
         gdp_pc_usd_diff = gdp_pc_usd - lag(gdp_pc_usd),
         log_gdp_pc_usd_diff = log(gdp_pc_usd) - log(lag(gdp_pc_usd)),
         gini_income_lag = lag(gini_income),
         gini_income_diff = gini_income - lag(gini_income),
         health_exp_share_lag = lag(health_exp_share),
         health_exp_share_diff = health_exp_share - lag(health_exp_share)) 

write_csv(all_panel_df, "data/clean/siler_econ_panel.csv")



ggplot(all_panel_df) + theme_bw() + 
  scale_color_manual("Country", values = col_scheme) + guides(color=guide_legend(ncol=2)) + 
  geom_line(aes(x = year, y = gini_income, color = plot_name, group = name))

### Income Gini ###
xvar <- "income_gini"
xlabel <- "Income Gini Index"
plot_cC_econ <- function(all_panel_df, xvar, xlabel, reg_method = "lm"){
  command <- paste0(
    "C_plt <- ggplot(all_panel_df) + theme_bw() + 
      scale_color_manual(\"Country\", values = col_scheme) + guides(color=guide_legend(ncol=2)) + 
      geom_point(data = all_panel_df[which(all_panel_df$plot_name == \"Other\"),], alpha = 0.5,
                 aes(x = ",xvar,", y = C, color = plot_name)) +
      geom_point(data = all_panel_df[which(all_panel_df$plot_name != \"Other\"),],
                 aes(x = ",xvar,", y = C, color = plot_name)) +
      geom_smooth(data = data.frame(all_panel_df),
                  aes(x = ",xvar,", y = C), method = \"",reg_method,"\", color = \"black\") +
      xlab(\"",xlabel,"\") + ylab(\"C\") + ggtitle(\"C\")
    c_plt <- ggplot(all_panel_df) + theme_bw() + 
      scale_color_manual(\"Country\", values = col_scheme) + guides(color=guide_legend(ncol=2)) + 
      geom_point(data = all_panel_df[which(all_panel_df$plot_name == \"Other\"),], alpha = 0.5,
                 aes(x = ",xvar,", y = c, color = plot_name)) +
      geom_point(data = all_panel_df[which(all_panel_df$plot_name != \"Other\"),],
                 aes(x = ",xvar,", y = c, color = plot_name)) +
      geom_smooth(data = data.frame(all_panel_df),
                  aes(x = ",xvar,", y = c), method = \"",reg_method,"\", color = \"black\") +
      xlab(\"",xlabel,"\") + ylab(\"c\") + ggtitle(\"c\")
    ggarrange(c_plt, C_plt, nrow = 1, ncol=2, common.legend = TRUE,
              legend = \"right\")"
  )
  plt <- eval(parse(text=command))
  return(plt)
  
}

plot_cC_econ(all_panel_df, "gini_income", "Income Gini Coefficient")
#ggsave("figures/countries/summary/Cc_income_gini.pdf", width = 10, height = 4)
plot_cC_econ(all_panel_df, "inc_share_bottom20", "inc_share_bottom20")
plot_cC_econ(all_panel_df, "pv_rate215", "pv_rate215")
plot_cC_econ(all_panel_df, "gdp_pc_growth", "")
plot_cC_econ(all_panel_df, "log_gdp_pc_usd_diff", "")
plot_cC_econ(all_panel_df, "health_exp_share_diff", "health_exp_share_diff")







### For paper
all_panel_df$US <- as.numeric(all_panel_df$name == "United States of America")

reg_df <- all_panel_df[complete.cases(all_panel_df[,c("c","C","income_gini","gdp_pc","pop", "health_share")]),] %>%
  group_by(year) %>%
  mutate(pop = pop/sum(pop))

reg_lag_df <- all_panel_df[complete.cases(all_panel_df[,c("c","C","pop","income_gini_lag", "gdp_pc_lag", "health_share_lag")]),] %>%
  group_by(year) %>%
  mutate(pop = pop/sum(pop))




model1 <- (felm(c ~ income_gini + log(gdp_pc) |
                            name + year, data = reg_df, weights = reg_df$pop))
model2 <- (felm(C ~ income_gini + log(gdp_pc) |
                  name + year, data = reg_df, weights = reg_df$pop))
model3 <- (felm(c ~ income_gini + log(gdp_pc) + health_share|
                  name + year, data = reg_df, weights = reg_df$pop))
model4 <- (felm(C ~ income_gini + log(gdp_pc) + health_share|
                  name + year, data = reg_df, weights = reg_df$pop))
model5 <- (felm(c ~ income_gini_lag + log(gdp_pc_lag) + health_share_lag|
                  name + year, data = reg_lag_df, weights = reg_lag_df$pop))
model6 <- (felm(C ~ income_gini_lag + log(gdp_pc_lag) + health_share_lag|
                  name + year, data = reg_lag_df, weights = reg_lag_df$pop))
models <- list(model1, model2, model3, model4, model5, model6)
stargazer(models, table.placement = "H",
          df = FALSE, title = "C and c on economic variables", label = "tab:cC_econ",
          font.size = "scriptsize")

# OLS
model1 <- (lm(c ~ income_gini + log(gdp_pc) , data = all_panel_df))
model2 <- (lm(C ~ income_gini + log(gdp_pc) , data = all_panel_df))
# Period FE
model3 <- (felm(c ~ income_gini + log(gdp_pc) |
                  year, data = all_panel_df))
model4 <- (felm(C ~ income_gini + log(gdp_pc) |
                  year, data = all_panel_df))
# Country and period FE
model5 <- (felm(c ~ income_gini + log(gdp_pc) |
                  name + year, data = all_panel_df))
model6 <- (felm(C ~ income_gini + log(gdp_pc) |
                  name + year, data = all_panel_df))
stargazer(model1, model2, model3, model4, model5, model6, table.placement = "H",
          df = FALSE, title = "C and c on economic variables", label = "tab:cC_econ",
          column.labels = c("OLS", "OLS", "t FE", "t FE", "i,t FE", "i,t FE"),
          font.size = "scriptsize")



# OLS
model1 <- (lm(c ~ income_gini + log(gdp_pc) + health_share + school_avg, data = all_panel_df))
model2 <- (lm(C ~ income_gini + log(gdp_pc) + health_share + school_avg, data = all_panel_df))
# Period FE
model3 <- (felm(c ~ income_gini + log(gdp_pc) + health_share + school_avg|
               year, data = all_panel_df))
model4 <- (felm(C ~ income_gini + log(gdp_pc) + health_share + school_avg|
               year, data = all_panel_df))
# Country and period FE
model5 <- (felm(c ~ income_gini + log(gdp_pc) + health_share + school_avg|
               name + year, data = all_panel_df))
model6 <- (felm(C ~ income_gini + log(gdp_pc) + health_share + school_avg|
             name + year, data = all_panel_df))
stargazer(model1, model2, model3, model4, model5, model6, table.placement = "H",
          df = FALSE, title = "C and c on economic variables, health and schooling", 
          label = "tab:cC_econ_health_edu",
          column.labels = c("OLS", "OLS", "t FE", "t FE", "i,t FE", "i,t FE"),
          font.size = "scriptsize")



model1 <- (felm(c_diff ~ c_lag + income_gini_lag + log(gdp_pc_lag) +
               income_gini_diff + log_gdp_pc_diff , data = all_panel_df))
model2 <- (felm(C_diff ~ C_lag + income_gini_lag + log(gdp_pc_lag) + 
               income_gini_diff + log_gdp_pc_diff , data = all_panel_df))

model3 <- (felm(c_diff ~ c_lag + income_gini_lag + log(gdp_pc_lag) +
               income_gini_diff + log_gdp_pc_diff | year, data = all_panel_df))
model4 <- (felm(C_diff ~ C_lag + income_gini_lag + log(gdp_pc_lag) + 
               income_gini_diff + log_gdp_pc_diff | year, data = all_panel_df))


model5 <- (felm(c_diff ~ c_lag + income_gini_lag + log(gdp_pc_lag) + health_share_lag + school_avg_lag + 
               income_gini_diff + log_gdp_pc_diff + health_share_diff + school_avg_diff | year, data = all_panel_df))
model6 <- (felm(C_diff ~ C_lag + income_gini_lag + log(gdp_pc_lag) + health_share_lag + school_avg + 
               income_gini_diff + log_gdp_pc_diff + health_share_diff + school_avg_diff | year, data = all_panel_df))

stargazer(model1, model2, model3, model4, model5, model6, table.placement = "H",
          df = FALSE, title = "C and c ECM type models", 
          label = "tab:cC_econ_ecm",
          column.labels = c("OLS", "OLS", "t FE", "t FE", "t FE", "t FE"),
          font.size = "scriptsize")


summary(felm(Lstar_99 ~ income_gini + log(gdp_pc)| year, data =all_panel_df))























############## Everything below is redundant, but kept in case we want to go back. 


"
Per capita income, wealth and GDP
"
wid_pc <- data.frame(read_xls("data/economic-inequality/WID_pc.xls", 
                              sheet = "Data", col_names = F))
names(wid_pc) <- c("country", "variable", "perc", "year", "value")
# Rename countries to match HMD
wid_pc$name <- wid_pc$country
wid_pc$name <- str_replace(wid_pc$name, "USA", "United States of America")
wid_pc$name <- str_replace(wid_pc$name, "Russian Federation", "Russia")
wid_pc$name <- str_replace(wid_pc$name, "Czech Republic", "Czechia")
wid_pc$name <- str_replace(wid_pc$name, "Korea", "South Korea")
wid_pc <- wid_pc[which(wid_pc$name != "Czechoslovakia"),]

# Simplify variable names
wid_pc$variable_full <- wid_pc$variable
wid_pc$variable[which(str_detect(wid_pc$variable, "average") & 
                        str_detect(wid_pc$variable, "national wealth"))] <- "wealth_pc"
wid_pc$variable[which(str_detect(wid_pc$variable, "average") & 
                        str_detect(wid_pc$variable, "National income"))] <- "income_pc"
wid_pc$variable[which(str_detect(wid_pc$variable, "average") & 
                        str_detect(wid_pc$variable, "domestic product"))] <- "gdp_pc"
# Keep only necessary columns
wid_pc <- wid_pc[,c("name", "variable", "year", "value")]
wid_pc <- data.frame(pivot_wider(wid_pc, id_cols = c(name, year), names_from = variable, 
                                 values_from = value))
# Create name for plotting
wid_pc$plot_name <- wid_pc$name
wid_pc$plot_name[which(!(wid_pc$plot_name %in% names(col_scheme)))] <- "Other"
# Test plot
wid_pc %>%
  ggplot() + theme_bw() + 
  #scale_color_manual("Country", values = col_scheme) + guides(color=guide_legend(ncol=2)) + 
  geom_line(aes(x = year, y = income_pc, color = plot_name, group = name))



"
Income inequality
"
wid_inc <- data.frame(read_xls("data/economic-inequality/WID_income_ineq.xls", 
                               sheet = "Data", col_names = F))
names(wid_inc) <- c("country", "variable", "perc", "year", "value")
# Rename countries to match HMD
wid_inc$name <- wid_inc$country
wid_inc$name <- str_replace(wid_inc$name, "USA", "United States of America")
wid_inc$name <- str_replace(wid_inc$name, "Russian Federation", "Russia")
wid_inc$name <- str_replace(wid_inc$name, "Czech Republic", "Czechia")
wid_inc$name <- str_replace(wid_inc$name, "Korea", "South Korea")

# Simplify variable names
wid_inc$variable_full <- wid_inc$variable
wid_inc$variable[which(str_detect(wid_inc$variable, "Top 1%"))] <- "income_top_1pc"
wid_inc$variable[which(str_detect(wid_inc$variable, "Top 10%"))] <- "income_top_10pc"
wid_inc$variable[which(str_detect(wid_inc$variable, "Middle 40%"))] <- "income_mid_40pc"
wid_inc$variable[which(str_detect(wid_inc$variable, "Bottom 50%"))] <- "income_bot_50pc"
wid_inc$variable[which(str_detect(wid_inc$variable, "Gini"))] <- "income_gini"
# Keep only necessary columns
wid_inc <- wid_inc[,c("name", "variable", "year", "value")]
wid_inc <- data.frame(pivot_wider(wid_inc, id_cols = c(name, year), names_from = variable, 
                                  values_from = value))
# Create name for plotting
wid_inc$plot_name <- wid_inc$name
wid_inc$plot_name[which(!(wid_inc$plot_name %in% names(col_scheme)))] <- "Other"
# Test plot
ggplot(wid_inc) + theme_bw() + 
  scale_color_manual("Country", values = col_scheme) + guides(color=guide_legend(ncol=2)) + 
  geom_line(aes(x = year, y = income_top_1pc, color = plot_name, group = name))




"
Wealth inequality
"
wid_wealth <- data.frame(read_xls("data/economic-inequality/WID_wealth_ineq.xls", 
                                  sheet = "Data", col_names = F))
names(wid_wealth) <- c("country", "variable", "perc", "year", "value")
# Rename countries to match HMD
wid_wealth$name <- wid_wealth$country
wid_wealth$name <- str_replace(wid_wealth$name, "USA", "United States of America")
wid_wealth$name <- str_replace(wid_wealth$name, "Russian Federation", "Russia")
wid_wealth$name <- str_replace(wid_wealth$name, "Czech Republic", "Czechia")
wid_wealth$name <- str_replace(wid_wealth$name, "Korea", "South Korea")

# Simplify variable names
wid_wealth$variable_full <- wid_wealth$variable
wid_wealth$variable[which(str_detect(wid_wealth$variable, "Top 1%"))] <- "wealth_top_1pc"
wid_wealth$variable[which(str_detect(wid_wealth$variable, "Top 10%"))] <- "wealth_top_10pc"
wid_wealth$variable[which(str_detect(wid_wealth$variable, "Middle 40%"))] <- "wealth_mid_40pc"
wid_wealth$variable[which(str_detect(wid_wealth$variable, "Bottom 50%"))] <- "wealth_bot_50pc"
wid_wealth$variable[which(str_detect(wid_wealth$variable, "Gini"))] <- "wealth_gini"
# Keep only necessary columns
wid_wealth <- wid_wealth[,c("name", "variable", "year", "value")]
wid_wealth <- data.frame(pivot_wider(wid_wealth, id_cols = c(name, year), names_from = variable, 
                                     values_from = value))
# Create name for plotting
wid_wealth$plot_name <- wid_wealth$name
wid_wealth$plot_name[which(!(wid_wealth$plot_name %in% names(col_scheme)))] <- "Other"
# Test plot
ggplot(wid_wealth) + theme_bw() + 
  scale_color_manual("Country", values = col_scheme) + guides(color=guide_legend(ncol=2)) + 
  geom_line(aes(x = year, y = wealth_top_1pc, color = plot_name, group = name))



"
Average years education 
"
owid_educ <- read.csv("data/economic-inequality/mean-years-of-schooling-long-run.csv", 
                      stringsAsFactors = FALSE)
names(owid_educ) <- c("name", "code", "year", "school_avg")
unique(siler_df$code)[!(unique(siler_df$code) %in% unique(owid_educ$code))]
# Merge in Siler names 
owid_educ <- merge(owid_educ[,c("code", "year", "school_avg")], 
                   unique(siler_df[,c("code", "name")]), by = "code")
owid_educ <- owid_educ[,c("name", "year", "school_avg")]


"
Govt healthcare spend
"
owid_health <- read.csv("data/economic-inequality/public-health-expenditure-share-GDP-OWID.csv", 
                        stringsAsFactors = FALSE)
names(owid_health) <- c("name", "code", "year", "health_share")
unique(siler_df$code)[!(unique(siler_df$code) %in% unique(owid_health$code))]
# Merge in Siler names 
owid_health <- merge(owid_health[,c("code", "year", "health_share")], 
                     unique(siler_df[,c("code", "name")]), by = "code")
owid_health <- owid_health[,c("name", "year", "health_share")]


"
Combine old availability table
"
gdp_sum <- wid_pc %>%
  group_by(name) %>%
  filter(!is.na(gdp_pc) & year <= 2020) %>%
  summarise(gdp_start = min(year), gdp_end = max(year)) %>%
  mutate(gdp_sample = paste0(gdp_start,":",gdp_end)) %>%
  select(name, gdp_sample)
gini_sum <- wid_inc %>%
  group_by(name) %>%
  filter(!is.na(income_gini) & year <= 2020) %>%
  summarise(gini_start = min(year), gini_end = max(year)) %>%
  mutate(gini_sample = paste0(gini_start,":",gini_end)) %>%
  select(name, gini_sample)
health_sum <- owid_health %>%
  group_by(name) %>%
  filter(!is.na(health_share) & year >= 1990 & year <= 2020) %>%
  summarise(health_start = min(year), health_end = max(year)) %>%
  mutate(health_sample = paste0(health_start,":",health_end)) %>%
  select(name, health_sample)
mort_sum <- mort_df %>%
  group_by(name) %>%
  filter(year <= 2020) %>%
  summarise(mort_start = min(year), mort_end = max(year)) %>%
  mutate(mort_sample = paste0(mort_start,":",mort_end)) %>%
  select(name, mort_sample)





summary_df <- mort_sum %>%
  full_join(gdp_sum) %>%
  full_join(gini_sum) %>%
  full_join(health_sum) %>%
  rename(Country = name, Mortality = mort_sample, GDP = gdp_sample, `Income Inequality` = gini_sample,
         `Health Share` = health_sample)
stargazer(as.matrix(summary_df), table.placement = "H", label = "tab:sample_starts",
          title = "Data availability by country")


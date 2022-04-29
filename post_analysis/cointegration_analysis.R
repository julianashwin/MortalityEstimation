setwd("/Users/julianashwin/Documents/GitHub/MortalityEstimation/")
rm(list=ls())

require(ggplot2)
require(ggpubr)
require(stringr)
require(plm)
require(lfe)
require(tidyr)


"
Set colour scheme
"
col_scheme <- c("Australia" = "darkolivegreen4", "Belgium" = "gold3", 
                "Canada" = "pink", "Switzerland" = "darkorchid3", 
                "Spain" = "darkorange", "Finland" = "darkslategray1", 
                "France" = "blue3", "United Kingdom" = "gray", 
                "Greece" = "lemonchiffon2", "Hong Kong" = "lightgoldenrod", 
                "Iceland" = "cornsilk3", "Italy" = "forestgreen", 
                "Japan" = "red", "South Korea" = "blue",
                "Norway" = "deeppink", "New Zealand" = "black", 
                "Portugal" = "green", "Sweden" = "yellow", 
                "United States of America" = "cornflowerblue",
                "Best Practice" = "darkmagenta")

"
Import data and results
"
# Import parameter estimates
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

# Merge in the actual data to get country names
mort_df <- read.csv("data/clean/all_lifetab_5y.csv", stringsAsFactors = FALSE)
mort_df <- mort_df[which(mort_df$age == 0),]
all_df <- merge(countries_df, unique(mort_df[,c("code", "name")]), by = "code")
all_df <- all_df[which(all_df$year > 0),]

# Merge in the best practice results
bp_df <- read.csv("figures/benchmark/siler_i2_preds.csv", stringsAsFactors = F)
bp_df$code <- "BP"
bp_df$name <- "Best Practice"
bp_df$best_practice <- 2
bp_df <- bp_df[,names(all_df)]
all_df <- rbind(all_df, bp_df)

"
Plot estimates and forecasts
"
## Forecast label for pretty legends
all_df$Forecast <- "Estimate"
all_df$Forecast[which(all_df$forecast == 1)] <- "Forecast"
## Life expectancy
plot_df <- all_df[which(all_df$parameter == "LE"),]
extra_obs <- plot_df[which(plot_df$year == 2018),]
extra_obs$Forecast <- "Estimate"
plot_df <- rbind(plot_df, extra_obs)
ggplot(plot_df) + theme_bw() + 
  scale_color_manual("Country", values = col_scheme) + 
  geom_line(aes(x = year, y = median, color = name, linetype = Forecast)) +
  xlab("Year") + ylab("Life expectancy at birth")
ggsave("figures/countries/summary/LE_international.pdf", width = 8, height = 4)
## Equality
plot_df <- all_df[which(all_df$parameter == "H"),]
extra_obs <- plot_df[which(plot_df$year == 2018),]
extra_obs$Forecast <- "Estimate"
plot_df <- rbind(plot_df, extra_obs)
ggplot(plot_df) + theme_bw() + 
  scale_color_manual("Country", values = col_scheme) + 
  geom_line(aes(x = year, y = -log(median), color = name, linetype = Forecast)) +
  xlab("Year") + ylab("Lifespan equality at birth")
ggsave("figures/countries/summary/h_international.pdf", width = 8, height = 4)
## c
plot_df <- all_df[which(all_df$parameter == "c"),]
extra_obs <- plot_df[which(plot_df$year == 2018),]
extra_obs$Forecast <- "Estimate"
plot_df <- rbind(plot_df, extra_obs)
ggplot(plot_df) + theme_bw() + 
  scale_color_manual("Country", values = col_scheme) + 
  geom_line(aes(x = year, y = median, color = name, linetype = Forecast)) +
  xlab("Year") + ylab("Senescent rectangularity (c)")
ggsave("figures/countries/summary/c_rect_international.pdf", width = 8, height = 4)
## C
plot_df <- all_df[which(all_df$parameter == "C"),]
extra_obs <- plot_df[which(plot_df$year == 2018),]
extra_obs$Forecast <- "Estimate"
plot_df <- rbind(plot_df, extra_obs)
ggplot(plot_df) + theme_bw() + 
  scale_color_manual("Country", values = col_scheme) + 
  geom_line(aes(x = year, y = median, color = name, linetype = Forecast)) +
  xlab("Year") + ylab("Senescent elongation (C)")
ggsave("figures/countries/summary/C_elg_international.pdf", width = 8, height = 4)
## b
plot_df <- all_df[which(all_df$parameter == "b"),]
extra_obs <- plot_df[which(plot_df$year == 2018),]
extra_obs$Forecast <- "Estimate"
plot_df <- rbind(plot_df, extra_obs)
ggplot(plot_df) + theme_bw() + 
  scale_color_manual("Country", values = col_scheme) + 
  geom_line(aes(x = year, y = median, color = name, linetype = Forecast)) +
  xlab("Year") + ylab("Infant rectangularity (b)")
ggsave("figures/countries/summary/b_rect_international.pdf", width = 8, height = 4)
## C
plot_df <- all_df[which(all_df$parameter == "B"),]
extra_obs <- plot_df[which(plot_df$year == 2018),]
extra_obs$Forecast <- "Estimate"
plot_df <- rbind(plot_df, extra_obs)
ggplot(plot_df) + theme_bw() + 
  scale_color_manual("Country", values = col_scheme) + 
  geom_line(aes(x = year, y = median, color = name, linetype = Forecast)) +
  xlab("Year") + ylab("Infant elongation (B)")
ggsave("figures/countries/summary/B_elg_international.pdf", width = 8, height = 4)
## d
plot_df <- all_df[which(all_df$parameter == "d"),]
extra_obs <- plot_df[which(plot_df$year == 2018),]
extra_obs$Forecast <- "Estimate"
plot_df <- rbind(plot_df, extra_obs)
ggplot(plot_df) + theme_bw() + 
  scale_color_manual("Country", values = col_scheme) + 
  geom_line(aes(x = year, y = median, color = name, linetype = Forecast)) +
  xlab("Year") + ylab("Age-independent mortality (d)")
ggsave("figures/countries/summary/d_international.pdf", width = 8, height = 4)


"
Cointegration tests
"
library(urca)
# Wide df for c
c_df <- all_df[which(all_df$parameter == "c" & all_df$year < 2022 
                       & all_df$year > 1900),]
c_df <- c_df[,c("name", "year", "median")]
rownames(c_df) <- NULL
c_wide_df <- spread(c_df, key = name, value = median)
# Wide df for C 
C_df <- all_df[which(all_df$parameter == "C" & all_df$year < 2022 
                     & all_df$year > 1900),]
C_df <- C_df[,c("name", "year", "median")]
rownames(C_df) <- NULL
C_wide_df <- spread(C_df, key = name, value = median)


# Unit root test on country - BP
coint_table <- data.frame(country = unique(all_df$name), c = NA, C = NA, N_obs = NA)
for (country in coint_table$country){
  if (country != "Best Practice"){
    # Test for c
    diff_c_bp = c_wide_df[,country] - c_wide_df[,"Best Practice"]
    diff_c_bp <- diff_c_bp[which(!is.na(diff_c_bp))]
    urtest_c = ur.df(diff_c_bp, type = "none", lags = 0)
    coint_table[which(coint_table$country == country), "c"] <- round(urtest_c@teststat[1],3)
    
    # Test for C
    diff_C_bp = C_wide_df[,country] - C_wide_df[,"Best Practice"]
    diff_C_bp <- diff_C_bp[which(!is.na(diff_C_bp))]
    urtest_C = ur.df(diff_C_bp, type = "none", lags = 0)
    coint_table[which(coint_table$country == country), "C"] <- round(urtest_C@teststat[1],3)
    
    coint_table[which(coint_table$country == country), "N_obs"] <- length(diff_C_bp)
  }
  
}

stargazer(as.matrix(coint_table), title = "Pairwise cointegration tests", 
          table.placement = "H", label = "tab:pair_coint")


"
Some panel regressions
"
# Convert to date and factor to facilitate panel analysis
full_parests_df <- all_parests_df[which(all_parests_df$code %in% full_codes),]
full_parests_df <- full_parests_df[which(!is.na(full_parests_df$year)),]
full_parests_df$year <- as.Date(paste0(full_parests_df$year, "-01-01"))
full_parests_df$time <- as.numeric(as.factor(full_parests_df$year))
full_parests_df$code <- as.factor(full_parests_df$code)

# Get separate df for each parameter
full_B_df <- pdata.frame(full_parests_df[which(full_parests_df$parameter == "B"),], index = c("code", "year"))
full_B_df$lB <- log(full_B_df$median)
full_B_df$lB_1diff <- full_B_df$lB - plm::lag(full_B_df$lB)

full_b_df <- pdata.frame(full_parests_df[which(full_parests_df$parameter == "b"),], index = c("code", "year"))
full_b_df$lb <- log(full_b_df$median)
full_b_df$lb_1diff <- full_b_df$lb - plm::lag(full_b_df$lb)

full_C_df <- pdata.frame(full_parests_df[which(full_parests_df$parameter == "C"),], index = c("code", "year"))
full_C_df$lC <- log(full_C_df$median)
full_C_df$lC_1diff <- full_C_df$lC - plm::lag(full_C_df$lC)

full_c_df <- pdata.frame(full_parests_df[which(full_parests_df$parameter == "c"),], index = c("code", "year"))
full_c_df$lc <- log(full_c_df$median)
full_c_df$lc_1diff <- full_c_df$lc - plm::lag(full_c_df$lc)

full_d_df <- pdata.frame(full_parests_df[which(full_parests_df$parameter == "d"),], index = c("code", "year"))
full_d_df$ld <- log(full_d_df$median)
full_d_df$ld_1diff <- full_d_df$ld - plm::lag(full_d_df$ld)

# Regs for B
model1 <- felm(lB ~ plm::lag(lB, 1) | code, data = full_B_df)
summary(model1)
model2 <- felm(lB ~ plm::lag(lB, 1) + time | code, data = full_B_df)
summary(model2)
model3 <- felm(lB_1diff ~ plm::lag(lB_1diff, 1) | code, data = full_B_df)
summary(model3)
model4 <- felm(lB_1diff ~ time + plm::lag(lB_1diff, 1) | code, data = full_B_df)
summary(model4)
model5 <- felm(lB_1diff ~  plm::lag(lB_1diff, 1) + plm::lag(lB_1diff, 2) | code, data = full_B_df)
summary(model5)
stargazer(model1,model2,model3, model4, model5, table.placement = "H", df = FALSE,
          title = "Panel model for B", font.size = "small")

# Regs for b
model1 <- felm(lb ~ plm::lag(lb, 1) | code, data = full_b_df)
summary(model1)
model2 <- felm(lb ~ plm::lag(lb, 1) + time| code, data = full_b_df)
summary(model2)
model3 <- felm(lb_1diff ~ plm::lag(lb_1diff, 1) | code, data = full_b_df)
summary(model3)
model4 <- felm(lb_1diff ~ time + plm::lag(lb_1diff, 1) | code, data = full_b_df)
summary(model4)
model5 <- felm(lb_1diff ~  plm::lag(lb_1diff, 1) + plm::lag(lb_1diff, 2) | code, data = full_b_df)
summary(model5)
stargazer(model1,model2,model3, model4, model5, table.placement = "H", df = FALSE,
          title = "Panel model for b", font.size = "small")


# Regs for C
model1 <- felm(lC ~ plm::lag(lC, 1) | code, data = full_C_df)
summary(model1)
model2 <- felm(lC ~ plm::lag(lC, 1) + time| code, data = full_C_df)
summary(model2)
model3 <- felm(lC_1diff ~ plm::lag(lC_1diff, 1) | code, data = full_C_df)
summary(model3)
model4 <- felm(lC_1diff ~ time + plm::lag(lC_1diff, 1) | code, data = full_C_df)
summary(model4)
model5 <- felm(lC_1diff ~  plm::lag(lC_1diff, 1) + plm::lag(lC_1diff, 2) | code, data = full_C_df)
summary(model5)
stargazer(model1,model2,model3, model4, model5, table.placement = "H", df = FALSE,
          title = "Panel model for C", font.size = "small")



# Regs for c
model1 <- felm(lc ~ plm::lag(lc, 1) | code, data = full_c_df)
summary(model1)
model2 <- felm(lc ~ plm::lag(lc, 1) + time| code, data = full_c_df)
summary(model2)
model3 <- felm(lc_1diff ~ plm::lag(lc_1diff, 1) | code, data = full_c_df)
summary(model3)
model4 <- felm(lc_1diff ~ time + plm::lag(lc_1diff, 1) | code, data = full_c_df)
summary(model4)
model5 <- felm(lc_1diff ~  plm::lag(lc_1diff, 1) + plm::lag(lc_1diff, 2) | code, data = full_c_df)
summary(model5)
stargazer(model1,model2,model3, model4, model5, table.placement = "H", df = FALSE,
          title = "Panel model for c", font.size = "small")


model1 <- felm(ld ~ plm::lag(ld, 1) | code, data = full_d_df)
summary(model1)
model2 <- felm(ld ~ plm::lag(ld, 1) + time| code, data = full_d_df)
summary(model2)
model3 <- felm(ld_1diff ~ plm::lag(ld_1diff, 1) | code, data = full_d_df)
summary(model3)
model4 <- felm(ld_1diff ~ time + plm::lag(ld_1diff, 1) | code, data = full_d_df)
summary(model4)
model5 <- felm(ld_1diff ~  plm::lag(ld_1diff, 1) + plm::lag(ld_1diff, 2) | code, data = full_d_df)
summary(model5)
stargazer(model1,model2,model3, model4, model5, table.placement = "H", df = FALSE,
          title = "Panel model for d", font.size = "small")



"
Plot the parameter estimates across all countries
"
B_plt <- ggplot(all_parests_df[which(all_parests_df$parameter == "B"),]) + theme_bw() + 
  theme(axis.title.y = element_text(angle = 0, vjust = 0.5)) +
  geom_line(aes(x = year, y = median, color = name)) + 
  xlab("Year") + ylab("B") + ggtitle("B parameter from dynamic Siler model")
b_plt <- ggplot(all_parests_df[which(all_parests_df$parameter == "b"),]) + theme_bw() + 
  theme(axis.title.y = element_text(angle = 0, vjust = 0.5)) +
  geom_line(aes(x = year, y = median, color = name)) + 
  xlab("Year") + ylab("b") + ggtitle("b parameter from dynamic Siler model")
C_plt <- ggplot(all_parests_df[which(all_parests_df$parameter == "C"),]) + theme_bw() + 
  theme(axis.title.y = element_text(angle = 0, vjust = 0.5)) +
  geom_line(aes(x = year, y = median, color = name)) + 
  xlab("Year") + ylab("C") + ggtitle("C parameter from dynamic Siler model")
C_plt
c_plt <- ggplot(all_parests_df[which(all_parests_df$parameter == "c"),]) + theme_bw() + 
  theme(axis.title.y = element_text(angle = 0, vjust = 0.5)) +
  geom_line(aes(x = year, y = median, color = name)) + 
  xlab("Year") + ylab("c") + ggtitle("c parameter from dynamic Siler model")
d_plt <- ggplot(all_parests_df[which(all_parests_df$parameter == "d"),]) + theme_bw() + 
  theme(axis.title.y = element_text(angle = 0, vjust = 0.5)) +
  geom_line(aes(x = year, y = median, color = name)) + 
  xlab("Year") + ylab("d") + ggtitle("d parameter from dynamic Siler model")
sigma_plt <- ggplot(all_parests_df[which(all_parests_df$parameter == "σ"),]) + theme_bw() + 
  theme(axis.title.y = element_text(angle = 0, vjust = 0.5)) +
  geom_line(aes(x = year, y = median, color = name)) + 
  xlab("Year") + ylab(expression(sigma)) + 
  ggtitle(expression(sigma~"parameter from dynamic Siler model"))


# Export plots
ggarrange(B_plt, b_plt, nrow = 1, ncol = 2, common.legend = TRUE, legend = "right")
ggsave("figures/all_HMD/infant_compare.pdf", width = 12, height = 4)
ggarrange(C_plt, c_plt, nrow = 1, ncol = 2, common.legend = TRUE, legend = "right")
ggsave("figures/all_HMD/elderly_compare.pdf", width = 12, height = 4)
ggarrange(d_plt, sigma_plt, nrow = 1, ncol = 2, common.legend = TRUE, legend = "right")
ggsave("figures/all_HMD/base_compare.pdf", width = 12, height = 4)


# Plot the estimates time series parameters






"
Dispersion across countries
"
years <- unique(full_parests_df$year[which(!is.na(full_parests_df$year))])
disp_df <- data.frame(year = years)
for (year in years){
  temp_df <- full_parests_df[which(full_parests_df$year == year),]  
  for (par in c("B", "b", "C", "c", "d", "σ")){
    disp_df[which(disp_df$year == year), paste0("log_",par, "_var")] <- 
      var(log(temp_df$median[which(temp_df$parameter == par)]), na.rm = TRUE)
    disp_df[which(disp_df$year == year), paste0(par, "_var")] <- 
      var((temp_df$median[which(temp_df$parameter == par)]), na.rm = TRUE)
    disp_df[which(disp_df$year == year), paste0(par, "_mean")] <- 
      mean((temp_df$median[which(temp_df$parameter == par)]), na.rm = TRUE)
  }
}


B_plt <- ggplot(disp_df) + theme_bw() + ylab("Var(log(B))") + xlab("") +
  geom_line(aes(x = year, y = log_B_var)) 
b_plt <- ggplot(disp_df) + theme_bw() + ylab("Var(log(b))") + xlab("") +
  geom_line(aes(x = year, y = log_b_var)) 
C_plt <- ggplot(disp_df) + theme_bw() + ylab("Var(log(C))") + xlab("") +
  geom_line(aes(x = year, y = log_C_var)) 
c_plt <- ggplot(disp_df) + theme_bw() + ylab("Var(log(c))") + xlab("") +
  geom_line(aes(x = year, y = log_c_var)) 
d_plt <- ggplot(disp_df) + theme_bw() + ylab("Var(log(d))") + xlab("") +
  geom_line(aes(x = year, y = log_d_var)) 
sigma_plt <- ggplot(disp_df) + theme_bw() + 
  ylab(expression("Var(log("~sigma~"))")) + xlab("Year") +
  geom_line(aes(x = year, y = log_σ_var))
ggarrange(B_plt, b_plt, C_plt, c_plt, d_plt, sigma_plt,
          nrow = 3, ncol = 2, common.legend = TRUE, legend = "right")
ggsave("figures/cross_country/dispersion.pdf", width = 5, height = 3)
  

















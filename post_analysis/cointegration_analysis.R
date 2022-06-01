setwd("/Users/julianashwin/Documents/GitHub/MortalityEstimation/")
rm(list=ls())

require(ggplot2)
require(ggpubr)
require(stringr)
require(plm)
require(lfe)
require(tidyr)
require(stargazer)

"
Set colour scheme
"
col_scheme <- c("Australia" = "darkolivegreen4", "Belgium" = "gold3", 
                "Canada" = "pink", "Switzerland" = "darkorchid3", 
                "Spain" = "firebrick", "Finland" = "darkslategray1", 
                "France" = "blue3", "United Kingdom" = "gray", 
                "Greece" = "lemonchiffon2", "Hong Kong" = "lightgoldenrod", 
                "Iceland" = "cornsilk3", "Italy" = "forestgreen", 
                "Japan" = "red", "South Korea" = "blue",
                "Norway" = "deeppink", "Netherlands" = "orange",
                "New Zealand" = "black", "Portugal" = "green", 
                "Sweden" = "yellow", "United States of America" = "cornflowerblue",
                "Best Practice" = "darkmagenta")

col_scheme <- c("Other" = "gray",
                "Australia" = "darkolivegreen4", 
                "Canada" = "pink",
                "France" = "blue3", "United Kingdom" = "darkgoldenrod4", 
                "Hong Kong" = "lightgoldenrod", 
                "Italy" = "forestgreen", 
                "Japan" = "red","New Zealand" = "black", "Russia" = "firebrick",
                "Sweden" = "yellow", "United States of America" = "cornflowerblue",
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

all_df$code <- str_replace(all_df$code, "NZL_NM", "NZL")

## Forecast label for pretty legends
all_df$Forecast <- "Estimate"
all_df$Forecast[which(all_df$forecast == 1)] <- "Forecast"


export_df <- all_df[which(all_df$year > 1900),]
export_df <- data.frame(pivot_wider(export_df, id_cols = c(name, code, year), names_from = parameter, 
                                  values_from = median))
write.csv(export_df, "data/results/siler_panel.csv", row.names = FALSE)
rm(export_df)


# Merge in the best practice results
bp_df <- read.csv("figures/benchmark/siler_i2drift_preds.csv", stringsAsFactors = F)
bp_df$code <- "BP"
bp_df$name <- "Best Practice"
bp_df$best_practice <- 2
bp_df <- bp_df[,names(all_df)]
#all_df <- rbind(all_df, bp_df)


cC_df <- all_df[which(all_df$parameter %in% c("C","c")),c("code", "name", "year","parameter", "median")]
cC_df <- data.frame(pivot_wider(cC_df, id_cols = c(code, year), names_from = parameter, 
                                values_from = c(median)))
cC_df <- data.frame(pivot_wider(cC_df, id_cols = c(year), names_from = code, 
                                       names_glue = "{.value}_{code}", values_from = c(c, C)))
cC_df <- cC_df <- cC_df[order(cC_df$year),]
write.csv(cC_df, "data/results/cC_panel.csv", row.names = FALSE)
rm(cC_df)



# Economic inequality data
ineq_df1 <- read.csv("data/economic-inequality/economic-inequality-gini-index.csv", 
                    stringsAsFactors = FALSE)
ineq_df2 <- read.csv("data/economic-inequality/gini-coefficient-equivalized-income-chartbook.csv",
                    stringsAsFactors = FALSE)
names(ineq_df2)[4] <- "Gini.coeff.att"
econ_df <- merge(ineq_df1, ineq_df2, by = c("Entity", "Code", "Year" ), all.x=T, all.y=T)
rm(ineq_df1,ineq_df2)
econ_df <- econ_df[which(econ_df$Code != ""),]
unique(all_df$code)[which(!(unique(all_df$code) %in% unique(econ_df$Code)))]

# GDP data
gdp_df <- read.csv("data/economic-inequality/real-gdp-per-capita-PennWT.csv", 
                   stringsAsFactors = FALSE)
names(gdp_df)[4] <- "real.gdp.pc"
econ_df <- merge(econ_df, gdp_df, by = c("Entity", "Code", "Year" ), all.x=T, all.y=T)
rm(gdp_df)

# Health expenditure data
health_df <- read.csv("data/economic-inequality/health-expenditure-and-financing-per-capita.csv",
                      stringsAsFactors = FALSE)
names(health_df)[4] <- "health.exp.pc"
health_df <- health_df[which(health_df$health.exp.pc != ".."),]
health_df$health.exp.pc <- as.numeric(health_df$health.exp.pc)*1.044 # deflator to match gdp
health_df <- health_df[which(health_df$Code != ""),]
unique(all_df$code)[which(!(unique(all_df$code) %in% unique(health_df$Code)))]
econ_df <- merge(econ_df, health_df, by = c("Entity", "Code", "Year" ), all.x=T, all.y=T)
rm(health_df)
econ_df <- econ_df[which(econ_df$Code != ""),]
unique(all_df$code)[which(!(unique(all_df$code) %in% unique(econ_df$Code)))]




# Fill in the blanks with previous values where possible
names(econ_df) <- tolower(names(econ_df))
econ_df <- econ_df[which(econ_df$code %in% unique(all_df$code)),]
econ_df <- econ_df[order(econ_df$code, econ_df$year),]
for (ii in 2:nrow(econ_df)){
  if (is.na(econ_df$gini.index[ii]) & (econ_df$code[ii] == econ_df$code[ii-1])){
    econ_df$gini.index[ii] <- econ_df$gini.index[ii-1]
  }
}

"
Plot estimates and forecasts
"
## Life expectancy
plot_df <- all_df[which(all_df$parameter == "LE"),]
extra_obs <- plot_df[which((plot_df$year == 2018 & !(plot_df$code %in% c("NZL","RUS")))|
                             (plot_df$year == 2013 & plot_df$code %in% c("NZL","RUS"))),]
extra_obs$Forecast <- "Estimate"
plot_df <- rbind(plot_df, extra_obs)
plot_df$name[which(!(plot_df$name %in% names(col_scheme)))] <- "Other"
ggplot(plot_df) + theme_bw() + 
  scale_color_manual("Country", values = col_scheme) + guides(color=guide_legend(ncol=2)) + 
  geom_line(data = plot_df[which(plot_df$name == "Other"),], alpha = 0.5,
            aes(x = year, y = median, color = name, 
                group = interaction(code, Forecast), linetype = Forecast)) +
  geom_line(data = plot_df[which(plot_df$name != "Other"),],
            aes(x = year, y = median, color = name, 
                group = interaction(code, Forecast), linetype = Forecast)) +
  xlab("Year") + ylab("Life expectancy at birth")
ggsave("figures/countries/summary/LE_international.pdf", width = 8, height = 4)
## Equality
plot_df <- all_df[which(all_df$parameter == "H"),]
extra_obs <- plot_df[which((plot_df$year == 2018 & !(plot_df$code %in% c("NZL","RUS")))|
                             (plot_df$year == 2013 & plot_df$code %in% c("NZL","RUS"))),]
extra_obs$Forecast <- "Estimate"
plot_df <- rbind(plot_df, extra_obs)
plot_df$name[which(!(plot_df$name %in% names(col_scheme)))] <- "Other"
ggplot(plot_df) + theme_bw() + 
  scale_color_manual("Country", values = col_scheme) + guides(color=guide_legend(ncol=2)) + 
  geom_line(data = plot_df[which(plot_df$name == "Other"),], alpha = 0.5,
            aes(x = year, y = -log(median), color = name, 
                group = interaction(code, Forecast), linetype = Forecast)) +
  geom_line(data = plot_df[which(plot_df$name != "Other"),],
            aes(x = year, y = -log(median), color = name, 
                group = interaction(code, Forecast), linetype = Forecast)) +
  xlab("Year") + ylab("Lifespan equality at birth")
ggsave("figures/countries/summary/h_international.pdf", width = 8, height = 4)
## Lifespan
plot_df <- all_df[which(all_df$parameter == "Lstar"),]
extra_obs <- plot_df[which((plot_df$year == 2018 & !(plot_df$code %in% c("NZL","RUS")))|
                             (plot_df$year == 2013 & plot_df$code %in% c("NZL","RUS"))),]
extra_obs$Forecast <- "Estimate"
plot_df <- rbind(plot_df, extra_obs)
plot_df$name[which(!(plot_df$name %in% names(col_scheme)))] <- "Other"
ggplot(plot_df) + theme_bw() + 
  scale_color_manual("Country", values = col_scheme) + guides(color=guide_legend(ncol=2)) + 
  geom_line(data = plot_df[which(plot_df$name == "Other"),], alpha = 0.5,
            aes(x = year, y = median, color = name, 
                group = interaction(code, Forecast), linetype = Forecast)) +
  geom_line(data = plot_df[which(plot_df$name != "Other"),],
            aes(x = year, y = median, color = name, 
                group = interaction(code, Forecast), linetype = Forecast)) +
  xlab("Year") + ylab("Lifespan")
ggsave("figures/countries/summary/Lstar_international.pdf", width = 8, height = 4)
## c
plot_df <- all_df[which(all_df$parameter == "c"),]
extra_obs <- plot_df[which((plot_df$year == 2018 & !(plot_df$code %in% c("NZL","RUS")))|
                             (plot_df$year == 2013 & plot_df$code %in% c("NZL","RUS"))),]
extra_obs$Forecast <- "Estimate"
plot_df <- rbind(plot_df, extra_obs)
plot_df$name[which(!(plot_df$name %in% names(col_scheme)))] <- "Other"
ggplot(plot_df) + theme_bw() + 
  scale_color_manual("Country", values = col_scheme) + guides(color=guide_legend(ncol=2)) + 
  geom_line(data = plot_df[which(plot_df$name == "Other"),], alpha = 0.5,
            aes(x = year, y = median, color = name, 
                group = interaction(code, Forecast), linetype = Forecast)) +
  geom_line(data = plot_df[which(plot_df$name != "Other"),],
            aes(x = year, y = median, color = name, 
                group = interaction(code, Forecast), linetype = Forecast)) +
  xlab("Year") + ylab("Senescent rectangularity (c)")
ggsave("figures/countries/summary/c_rect_international.pdf", width = 8, height = 4)
## C
plot_df <- all_df[which(all_df$parameter == "C"),]
extra_obs <- plot_df[which((plot_df$year == 2018 & !(plot_df$code %in% c("NZL","RUS")))|
                             (plot_df$year == 2013 & plot_df$code %in% c("NZL","RUS"))),]
extra_obs$Forecast <- "Estimate"
plot_df <- rbind(plot_df, extra_obs)
plot_df$name[which(!(plot_df$name %in% names(col_scheme)))] <- "Other"
ggplot(plot_df) + theme_bw() + 
  scale_color_manual("Country", values = col_scheme) + guides(color=guide_legend(ncol=2)) + 
  geom_line(data = plot_df[which(plot_df$name == "Other"),], alpha = 0.5,
            aes(x = year, y = median, color = name, 
                group = interaction(code, Forecast), linetype = Forecast)) +
  geom_line(data = plot_df[which(plot_df$name != "Other"),],
            aes(x = year, y = median, color = name, 
                group = interaction(code, Forecast), linetype = Forecast)) +
  xlab("Year") + ylab("Senescent elongation (C)")
ggsave("figures/countries/summary/C_elg_international.pdf", width = 8, height = 4)
## b
plot_df <- all_df[which(all_df$parameter == "b"),]
extra_obs <- plot_df[which((plot_df$year == 2018 & !(plot_df$code %in% c("NZL","RUS")))|
                             (plot_df$year == 2013 & plot_df$code %in% c("NZL","RUS"))),]
extra_obs$Forecast <- "Estimate"
plot_df <- rbind(plot_df, extra_obs)
plot_df$name[which(!(plot_df$name %in% names(col_scheme)))] <- "Other"
ggplot(plot_df) + theme_bw() + 
  scale_color_manual("Country", values = col_scheme) + guides(color=guide_legend(ncol=2)) + 
  geom_line(data = plot_df[which(plot_df$name == "Other"),], alpha = 0.5,
            aes(x = year, y = median, color = name, 
                group = interaction(code, Forecast), linetype = Forecast)) +
  geom_line(data = plot_df[which(plot_df$name != "Other"),],
            aes(x = year, y = median, color = name, 
                group = interaction(code, Forecast), linetype = Forecast)) +
  xlab("Year") + ylab("Infant rectangularity (b)")
ggsave("figures/countries/summary/b_rect_international.pdf", width = 8, height = 4)
## B
plot_df <- all_df[which(all_df$parameter == "B"),]
extra_obs <- plot_df[which((plot_df$year == 2018 & !(plot_df$code %in% c("NZL","RUS")))|
                             (plot_df$year == 2013 & plot_df$code %in% c("NZL","RUS"))),]
extra_obs$Forecast <- "Estimate"
plot_df <- rbind(plot_df, extra_obs)
plot_df$name[which(!(plot_df$name %in% names(col_scheme)))] <- "Other"
ggplot(plot_df) + theme_bw() + 
  scale_color_manual("Country", values = col_scheme) + guides(color=guide_legend(ncol=2)) + 
  geom_line(data = plot_df[which(plot_df$name == "Other"),], alpha = 0.5,
            aes(x = year, y = median, color = name, 
                group = interaction(code, Forecast), linetype = Forecast)) +
  geom_line(data = plot_df[which(plot_df$name != "Other"),],
            aes(x = year, y = median, color = name, 
                group = interaction(code, Forecast), linetype = Forecast)) +
  xlab("Year") + ylab("Infant elongation (B)")
ggsave("figures/countries/summary/B_elg_international.pdf", width = 8, height = 4)
## d
plot_df <- all_df[which(all_df$parameter == "d"),]
extra_obs <- plot_df[which((plot_df$year == 2018 & !(plot_df$code %in% c("NZL","RUS")))|
                             (plot_df$year == 2013 & plot_df$code %in% c("NZL","RUS"))),]
extra_obs$Forecast <- "Estimate"
plot_df <- rbind(plot_df, extra_obs)
plot_df$name[which(!(plot_df$name %in% names(col_scheme)))] <- "Other"
ggplot(plot_df) + theme_bw() + 
  scale_color_manual("Country", values = col_scheme) + guides(color=guide_legend(ncol=2)) + 
  geom_line(data = plot_df[which(plot_df$name == "Other"),], alpha = 0.5,
            aes(x = year, y = median, color = name, 
                group = interaction(code, Forecast), linetype = Forecast)) +
  geom_line(data = plot_df[which(plot_df$name != "Other"),],
            aes(x = year, y = median, color = name, 
                group = interaction(code, Forecast), linetype = Forecast)) +
  xlab("Year") + ylab("Age-independent mortality (d)")
ggsave("figures/countries/summary/d_international.pdf", width = 8, height = 4)





"
Income vs lifespan inequality
"
plot_df <- all_df[which(all_df$year > 1900),]
plot_df <- data.frame(pivot_wider(plot_df, id_cols = c(name, code, year), names_from = parameter, 
                                  values_from = median))
plot_df <- merge(plot_df, econ_df, by = c("code", "year"), all.x = T)
plot_df$real.gdp.pc <- plot_df$real.gdp.pc/10000
plot_df$health.exp.pc <- plot_df$health.exp.pc/10000
plot_df$gini.index <- plot_df$gini.index/100
plot_df$name[which(!(plot_df$name %in% names(col_scheme)))] <- "Other"
plot_df$Period <- "pred"
plot_df$Period[which(plot_df$year <= 1943)] <- "1918-1943"
plot_df$Period[which(plot_df$year > 1943 & plot_df$year <= 1968)] <- "1943-1968"
plot_df$Period[which(plot_df$year > 1968 & plot_df$year <= 1993)] <- "1968-1993"
plot_df$Period[which(plot_df$year > 1993 & plot_df$year <= 2018)] <- "1993-2018"
table(plot_df$Period)
#plot_df <- plot_df[which(!is.na(plot_df$gini.index) | !is.na(plot_df$gini.coeff.att)),]


## C and c vs GDP per capita
plot_df1 <- plot_df[!is.na(plot_df$real.gdp.pc),]
C_plt <- ggplot(plot_df1) + theme_bw() + 
  scale_color_manual("Country", values = col_scheme) + guides(color=guide_legend(ncol=2)) + 
  geom_point(data = plot_df1[which(plot_df1$name == "Other"),], alpha = 0.5,
             aes(x = real.gdp.pc, y = C, color = name), size = 1) +
  geom_path(data = plot_df1[which(plot_df1$name == "Other"),], alpha = 0.2,
            aes(x = real.gdp.pc, y = C, color = name)) +
  geom_point(data = plot_df1[which(plot_df1$name != "Other"),],
             aes(x = real.gdp.pc, y = C, color = name), size = 1) +
  geom_path(data = plot_df1[which(plot_df1$name != "Other"),],
            aes(x = real.gdp.pc, y = C, color = name), alpha = 0.3) +
  geom_smooth(aes(x = real.gdp.pc, y = C), method = "lm", color = "black") +
  xlab("Real GDP p.c.") + ylab("C") + ggtitle("C")
c_plt <- ggplot(plot_df1) + theme_bw() + 
  scale_color_manual("Country", values = col_scheme) + guides(color=guide_legend(ncol=2)) + 
  geom_point(data = plot_df1[which(plot_df1$name == "Other"),], alpha = 0.5,
             aes(x = real.gdp.pc, y = c, color = name)) +
  geom_path(data = plot_df1[which(plot_df1$name == "Other"),], alpha = 0.2,
            aes(x = real.gdp.pc, y = c, color = name)) + 
  geom_point(data = plot_df1[which(plot_df1$name != "Other"),],
             aes(x = real.gdp.pc, y = c, color = name)) +
  geom_path(data = plot_df1[which(plot_df1$name != "Other"),],
            aes(x = real.gdp.pc, y = c, color = name), alpha = 0.3) + 
  geom_smooth(aes(x = real.gdp.pc, y = c), method = "lm", color = "black") +
  xlab("Real GDP p.c.") + ylab("c") + ggtitle("c")
ggarrange(c_plt, C_plt, nrow = 1, ncol=2, common.legend = TRUE,
          legend = "right")
ggsave("figures/countries/summary/Cc_gdp.pdf", width = 8, height = 4)


C_plt <- ggplot(plot_df1) + theme_bw() + 
  geom_point(aes(x = real.gdp.pc, y = C, color = Period), size = 1, alpha = 0.5) +
  geom_smooth(aes(x = real.gdp.pc, y = C, color = Period), method = "lm", se = FALSE) +
  xlab("Real GDP p.c.") + ylab("C") + ggtitle("C")
c_plt <- ggplot(plot_df1) + theme_bw() + 
  geom_point(aes(x = real.gdp.pc, y = c, color = Period), size = 1, alpha = 0.5) +
  geom_smooth(aes(x = real.gdp.pc, y = c, color = Period), method = "lm", se = FALSE) +
  xlab("Real GDP p.c.") + ylab("c") + ggtitle("c")
ggarrange(c_plt, C_plt, nrow = 1, ncol=2, common.legend = TRUE,
          legend = "right")
ggsave("figures/countries/summary/Cc_gdp_byperiod.pdf", width = 8, height = 4)




## C and c vs Health expenditure
plot_df1 <- plot_df[!is.na(plot_df$health.exp.pc),]
C_plt <- ggplot(plot_df1) + theme_bw() + 
  scale_color_manual("Country", values = col_scheme) + guides(color=guide_legend(ncol=2)) + 
  geom_point(data = plot_df1[which(plot_df1$name == "Other"),], alpha = 0.5,
             aes(x = health.exp.pc, y = C, color = name), size = 1) +
  geom_path(data = plot_df1[which(plot_df1$name == "Other"),], alpha = 0.2,
            aes(x = health.exp.pc, y = C, color = name)) +
  geom_point(data = plot_df1[which(plot_df1$name != "Other"),],
             aes(x = health.exp.pc, y = C, color = name), size = 1) +
  geom_path(data = plot_df1[which(plot_df1$name != "Other"),],
            aes(x = health.exp.pc, y = C, color = name), alpha = 0.3) +
  geom_smooth(aes(x = health.exp.pc, y = C), method = "lm", color = "black") +
  xlab("Health expenditure p.c.") + ylab("C") + ggtitle("C")
c_plt <- ggplot(plot_df1) + theme_bw() + 
  scale_color_manual("Country", values = col_scheme) + guides(color=guide_legend(ncol=2)) + 
  geom_point(data = plot_df1[which(plot_df1$name == "Other"),], alpha = 0.5,
             aes(x = health.exp.pc, y = c, color = name)) +
  geom_path(data = plot_df1[which(plot_df1$name == "Other"),], alpha = 0.2,
            aes(x = health.exp.pc, y = c, color = name)) + 
  geom_point(data = plot_df1[which(plot_df1$name != "Other"),],
             aes(x = health.exp.pc, y = c, color = name)) +
  geom_path(data = plot_df1[which(plot_df1$name != "Other"),],
            aes(x = health.exp.pc, y = c, color = name), alpha = 0.3) + 
  geom_smooth(aes(x = health.exp.pc, y = c), method = "lm", color = "black") +
  xlab("Health expenditure p.c.") + ylab("c") + ggtitle("c")
ggarrange(c_plt, C_plt, nrow = 1, ncol=2, common.legend = TRUE,
          legend = "right")
ggsave("figures/countries/summary/Cc_hexp.pdf", width = 8, height = 4)

C_plt <- ggplot(plot_df1) + theme_bw() + 
  geom_point(aes(x = health.exp.pc, y = C, color = Period), size = 1, alpha = 0.5) +
  geom_smooth(aes(x = health.exp.pc, y = C, color = Period), method = "lm", se = FALSE) +
  xlab("Health expenditure p.c.") + ylab("C") + ggtitle("C")
c_plt <- ggplot(plot_df1) + theme_bw() + 
  geom_point(aes(x = health.exp.pc, y = c, color = Period), size = 1, alpha = 0.5) +
  geom_smooth(aes(x = health.exp.pc, y = c, color = Period), method = "lm", se = FALSE) +
  xlab("Health expenditure p.c.") + ylab("c") + ggtitle("c")
ggarrange(c_plt, C_plt, nrow = 1, ncol=2, common.legend = TRUE,
          legend = "right")
ggsave("figures/countries/summary/Cc_hexp_byperiod.pdf", width = 8, height = 4)


## C and c vs Health expenditure/GDP
plot_df1 <- plot_df[!is.na(plot_df$health.exp.pc/plot_df$real.gdp.pc),]
C_plt <- ggplot(plot_df1) + theme_bw() + 
  scale_color_manual("Country", values = col_scheme) + guides(color=guide_legend(ncol=2)) + 
  geom_point(data = plot_df1[which(plot_df1$name == "Other"),], alpha = 0.5,
             aes(x = health.exp.pc/real.gdp.pc, y = C, color = name), size = 1) +
  geom_path(data = plot_df1[which(plot_df1$name == "Other"),], alpha = 0.2,
            aes(x = health.exp.pc/real.gdp.pc, y = C, color = name)) +
  geom_point(data = plot_df1[which(plot_df1$name != "Other"),],
             aes(x = health.exp.pc/real.gdp.pc, y = C, color = name), size = 1) +
  geom_path(data = plot_df1[which(plot_df1$name != "Other"),],
            aes(x = health.exp.pc/real.gdp.pc, y = C, color = name), alpha = 0.3) +
  geom_smooth(aes(x = health.exp.pc/real.gdp.pc, y = C), method = "lm", color = "black") +
  xlab("Health expenditure share of GDP") + ylab("C") + ggtitle("C")
c_plt <- ggplot(plot_df1) + theme_bw() + 
  scale_color_manual("Country", values = col_scheme) + guides(color=guide_legend(ncol=2)) + 
  geom_point(data = plot_df1[which(plot_df1$name == "Other"),], alpha = 0.5,
             aes(x = health.exp.pc/real.gdp.pc, y = c, color = name)) +
  geom_path(data = plot_df1[which(plot_df1$name == "Other"),], alpha = 0.2,
            aes(x = health.exp.pc/real.gdp.pc, y = c, color = name)) + 
  geom_point(data = plot_df1[which(plot_df1$name != "Other"),],
             aes(x = health.exp.pc/real.gdp.pc, y = c, color = name)) +
  geom_path(data = plot_df1[which(plot_df1$name != "Other"),],
            aes(x = health.exp.pc/real.gdp.pc, y = c, color = name), alpha = 0.3) + 
  geom_smooth(aes(x = health.exp.pc/real.gdp.pc, y = c), method = "lm", color = "black") +
  xlab("Health expenditure share of GDP") + ylab("c") + ggtitle("c")
ggarrange(c_plt, C_plt, nrow = 1, ncol=2, common.legend = TRUE,
          legend = "right")
ggsave("figures/countries/summary/Cc_hexp_share.pdf", width = 8, height = 4)

C_plt <- ggplot(plot_df1) + theme_bw() + 
  geom_point(aes(x = health.exp.pc/real.gdp.pc, y = C, color = Period), size = 1, alpha = 0.5) +
  geom_smooth(aes(x = health.exp.pc/real.gdp.pc, y = C, color = Period), method = "lm", se = FALSE) +
  xlab("Health expenditure share of GDP") + ylab("C") + ggtitle("C")
c_plt <- ggplot(plot_df1) + theme_bw() + 
  geom_point(aes(x = health.exp.pc/real.gdp.pc, y = c, color = Period), size = 1, alpha = 0.5) +
  geom_smooth(aes(x = health.exp.pc/real.gdp.pc, y = c, color = Period), method = "lm", se = FALSE) +
  xlab("Health expenditure share of GDP") + ylab("c") + ggtitle("c")
ggarrange(c_plt, C_plt, nrow = 1, ncol=2, common.legend = TRUE,
          legend = "right")
ggsave("figures/countries/summary/Cc_hexp_share_byperiod.pdf", width = 8, height = 4)


C_plt <- ggplot(plot_df1) + theme_bw() + 
  geom_point(aes(x = health.exp.pc/real.gdp.pc, y = LE, color = Period), size = 1, alpha = 0.5) +
  geom_smooth(aes(x = health.exp.pc/real.gdp.pc, y = LE, color = Period), method = "lm", se = FALSE) +
  xlab("Health expenditure share of GDP") + ylab("LE") + ggtitle("LE")
c_plt <- ggplot(plot_df1) + theme_bw() + 
  geom_point(aes(x = health.exp.pc/real.gdp.pc, y = h, color = Period), size = 1, alpha = 0.5) +
  geom_smooth(aes(x = health.exp.pc/real.gdp.pc, y = h, color = Period), method = "lm", se = FALSE) +
  xlab("Health expenditure share of GDP") + ylab("h") + ggtitle("h")
ggarrange(c_plt, C_plt, nrow = 1, ncol=2, common.legend = TRUE,
          legend = "right")
ggsave("figures/countries/summary/LEh_hexp_share_byperiod.pdf", width = 8, height = 4)



## C and c vs inequality
plot_df1 <- plot_df[!is.na(plot_df$gini.index),]
C_plt <- ggplot(plot_df1) + theme_bw() + 
  scale_color_manual("Country", values = col_scheme) + guides(color=guide_legend(ncol=2)) + 
  geom_point(data = plot_df1[which(plot_df1$name == "Other"),], alpha = 0.5,
             aes(x = gini.index, y = C, color = name)) +
  geom_path(data = plot_df1[which(plot_df1$name == "Other"),], alpha = 0.2,
            aes(x = gini.index, y = C, color = name)) + 
  geom_point(data = plot_df1[which(plot_df1$name != "Other"),],
             aes(x = gini.index, y = C, color = name)) +
  geom_path(data = plot_df1[which(plot_df1$name != "Other"),], 
            aes(x = gini.index, y = C, color = name)) + 
  geom_smooth(aes(x = gini.index, y = C), method = "lm", color = "black") +
  xlab("Income Gini Index") + ylab("C") + ggtitle("C")
c_plt <- ggplot(plot_df1) + theme_bw() + 
  scale_color_manual("Country", values = col_scheme) + guides(color=guide_legend(ncol=2)) + 
  geom_point(data = plot_df1[which(plot_df1$name == "Other"),], alpha = 0.5,
             aes(x = gini.index, y = c, color = name)) +
  geom_path(data = plot_df1[which(plot_df1$name == "Other"),], alpha = 0.2,
            aes(x = gini.index, y = c, color = name)) + 
  geom_point(data = plot_df1[which(plot_df1$name != "Other"),],
             aes(x = gini.index, y = c, color = name)) +
  geom_path(data = plot_df1[which(plot_df1$name != "Other"),],
            aes(x = gini.index, y = c, color = name)) + 
  geom_smooth(aes(x = gini.index, y = c), method = "lm", color = "black") +
  xlab("Income Gini Index") + ylab("c") + ggtitle("c")
ggarrange(c_plt, C_plt, nrow = 1, ncol=2, common.legend = TRUE,
          legend = "right")
ggsave("figures/countries/summary/Cc_gini.pdf", width = 8, height = 4)


C_plt <- ggplot(plot_df1) + theme_bw() + 
  geom_point(aes(x = gini.index, y = C, color = Period), size = 1, alpha = 0.5) +
  geom_smooth(aes(x = gini.index, y = C, color = Period), method = "lm", se = FALSE) +
  xlab("Income Gini Index") + ylab("C") + ggtitle("C")
c_plt <- ggplot(plot_df1) + theme_bw() + 
  geom_point(aes(x = gini.index, y = c, color = Period), size = 1, alpha = 0.5) +
  geom_smooth(aes(x = gini.index, y = c, color = Period), method = "lm", se = FALSE) +
  xlab("Income Gini Index") + ylab("c") + ggtitle("c")
ggarrange(c_plt, C_plt, nrow = 1, ncol=2, common.legend = TRUE,
          legend = "right")
ggsave("figures/countries/summary/Cc_gini_byperiod.pdf", width = 8, height = 4)



## Regressions
plot_df$health.share <- plot_df$health.exp.pc/plot_df$real.gdp.pc

plot_df <- plot_df[order(plot_df$code, plot_df$year),]
plot_df$year_num <- as.numeric(as.factor(plot_df$year))

write.csv(plot_df, "data/results/econ_panel.csv", row.names = FALSE)


panel_df <- pdata.frame(data.frame(plot_df), index = c("code", "year_num"))
panel_df <- panel_df[which(panel_df$year < 2020),]

# Lags
panel_df$gini.index_1lag <- plm::lag(panel_df$gini.index, 1)
panel_df$real.gdp.pc_1lag <- plm::lag(panel_df$real.gdp.pc, 1)
panel_df$health.share_1lag <- plm::lag(panel_df$health.share, 1)
panel_df$c_1lag <- plm::lag(panel_df$c, 1)
panel_df$C_1lag <- plm::lag(panel_df$C, 1)
# First diffs
panel_df$gini.index_diff <- panel_df$gini.index - panel_df$gini.index_1lag
panel_df$real.gdp.pc_diff <- panel_df$real.gdp.pc - panel_df$real.gdp.pc_1lag
panel_df$health.share_diff <- panel_df$health.share - panel_df$health.share_1lag
panel_df$c_diff <- panel_df$c - panel_df$c_1lag
panel_df$C_diff <- panel_df$C - panel_df$C_1lag

names(panel_df) <- str_replace(names(panel_df), "α", "alpha")
names(panel_df) <- str_replace(names(panel_df), "σ", "sigma")

write.csv(panel_df, "data/results/econ_panel.csv", row.names = FALSE)



require(plm)

# c on econ
model1 <- felm(c ~ gini.index + real.gdp.pc + health.share, data = panel_df)
summary(model1)
model1a <- felm(c ~ gini.index_1lag + real.gdp.pc_1lag + health.share_1lag, data = panel_df)
summary(model1a)
model1b <- felm(c ~ gini.index + real.gdp.pc + health.share + gini.index_1lag + 
                  real.gdp.pc_1lag + health.share_1lag, data = panel_df)
summary(model1b)
# c with period f.e.
model2 <- felm(c ~ gini.index + real.gdp.pc + health.share| year, data = panel_df)
summary(model2)
model2a <- felm(c ~ gini.index_1lag + real.gdp.pc_1lag + health.share_1lag| year, data = panel_df)
summary(model2a)
model2b <- felm(c ~ gini.index + real.gdp.pc + health.share + gini.index_1lag + 
                  real.gdp.pc_1lag + health.share_1lag| year, data = panel_df)
summary(model2b)
# c with period and country f.e.
model3 <- felm(c ~ gini.index + real.gdp.pc + health.share| code + year, data = panel_df)
summary(model3)
model3a <- felm(c ~ gini.index_1lag + real.gdp.pc_1lag + health.share_1lag| code + year, data = panel_df)
summary(model3a)
model3b <- felm(c ~ gini.index + real.gdp.pc + health.share + gini.index_1lag + 
                  real.gdp.pc_1lag + health.share_1lag| code + year, data = panel_df)
summary(model3b)

stargazer(model1, model1a, model1b, model2, model2a, model2b, model3, model3a, model3b,
          table.placement = "H", df = FALSE, title = "c (rectangularisation) on econ variables")

# C on econ
model1 <- felm(C ~ gini.index + real.gdp.pc + health.share, data = panel_df)
summary(model1)
model1a <- felm(C ~ gini.index_1lag + real.gdp.pc_1lag + health.share_1lag, data = panel_df)
summary(model1a)
model1b <- felm(C ~ gini.index + real.gdp.pc + health.share + gini.index_1lag + 
                  real.gdp.pc_1lag + health.share_1lag, data = panel_df)
summary(model1b)
# C with period f.e.
model2 <- felm(C ~ gini.index + real.gdp.pc + health.share| year, data = panel_df)
summary(model2)
model2a <- felm(C ~ gini.index_1lag + real.gdp.pc_1lag + health.share_1lag| year, data = panel_df)
summary(model2a)
model2b <- felm(C ~ gini.index + real.gdp.pc + health.share + gini.index_1lag + 
                  real.gdp.pc_1lag + health.share_1lag| year, data = panel_df)
summary(model2b)
# C with period and country f.e.
model3 <- felm(C ~ gini.index + real.gdp.pc + health.share| code + year, data = panel_df)
summary(model3)
model3a <- felm(C ~ gini.index_1lag + real.gdp.pc_1lag + health.share_1lag| code + year, data = panel_df)
summary(model3a)
model3b <- felm(C ~ gini.index + real.gdp.pc + health.share + gini.index_1lag + 
                  real.gdp.pc_1lag + health.share_1lag| code + year, data = panel_df)
summary(model3b)

stargazer(model1, model1a, model1b, model2, model2a, model2b, model3, model3a, model3b,
          table.placement = "H", df = FALSE, title = "C (elongation) on econ variables")


xs_df2013 <- plot_df[which(plot_df$year == 2013),]
model1 <- lm(c ~ gini.index + real.gdp.pc + health.share, data = xs_df2013)
summary(model1)
model2 <- lm(C ~ gini.index + real.gdp.pc + health.share, data = xs_df2013)
summary(model2)
stargazer(model1, model2, table.placement = "H", df = FALSE, 
          title = "2013 Cross-section")


### ECM ###
model1 <- felm(c_diff ~ c_1lag + gini.index_diff + real.gdp.pc_diff + health.share_diff + 
                 gini.index_1lag + real.gdp.pc_1lag + health.share_1lag, data = panel_df)
summary(model1)
model1a <- felm(c_diff ~ c_1lag + gini.index_diff + real.gdp.pc_diff + health.share_diff + 
                 gini.index_1lag + real.gdp.pc_1lag + health.share_1lag | year, data = panel_df)
summary(model1a)
model1b <- felm(c_diff ~ c_1lag + gini.index_diff + real.gdp.pc_diff + health.share_diff + 
                  gini.index_1lag + real.gdp.pc_1lag + health.share_1lag | code + year, data = panel_df)
summary(model1b)


model2 <- felm(C_diff ~ C_1lag + gini.index_diff + real.gdp.pc_diff + health.share_diff + 
                 gini.index_1lag + real.gdp.pc_1lag + health.share_1lag, data = panel_df)
summary(model2)
model2a <- felm(C_diff ~ C_1lag + gini.index_diff + real.gdp.pc_diff + health.share_diff + 
                  gini.index_1lag + real.gdp.pc_1lag + health.share_1lag | year, data = panel_df)
summary(model2a)
model2b <- felm(C_diff ~ C_1lag + gini.index_diff + real.gdp.pc_diff + health.share_diff + 
                  gini.index_1lag + real.gdp.pc_1lag + health.share_1lag | code + year, data = panel_df)
summary(model2b)

stargazer(model1, model1a, model1b, model2, model2a, model2b,
          table.placement = "H", df = FALSE, title = "ECM type regressions for c and C econ variables")










model1 <- felm(LE ~ gini.index + real.gdp.pc + health.share | year, data = plot_df)
summary(model1)
model2 <- felm(Lmed ~ gini.index + real.gdp.pc + health.share | year, data = plot_df)
summary(model2)
model3 <- felm(Lstar ~ gini.index + real.gdp.pc + health.share | year, data = plot_df)
summary(model3)
model4 <- felm(H ~ gini.index + real.gdp.pc + health.share | year, data = plot_df)
summary(model4)
model5 <- felm(h ~ gini.index + real.gdp.pc + health.share | year, data = plot_df)
summary(model5)

stargazer(model1, model2, model3, model4, model5,
          table.placement = "H", df = FALSE)





"
C-c path over time
"
plot_df <- all_df[which(all_df$year > 1900),]
plot_df <- data.frame(pivot_wider(plot_df, id_cols = c(name, code, year, Forecast), 
                                  names_from = parameter, values_from = median))
extra_obs <- plot_df[which((plot_df$year == 2018 & !(plot_df$code %in% c("NZL","RUS")))|
                             (plot_df$year == 2013 & plot_df$code %in% c("NZL","RUS"))),]
extra_obs$Forecast <- "Estimate"
plot_df <- rbind(plot_df, extra_obs)
plot_df$h <- -log(plot_df$H)
plot_df$name[which(!(plot_df$name %in% names(col_scheme)))] <- "Other"
# Include some years as labels
plot_df$year_label <- NA
plot_df$year_label[which(plot_df$year==1903)] <- 1903
plot_df$year_label[which(plot_df$year==1933)] <- 1933
plot_df$year_label[which(plot_df$year==1963)] <- 1963
plot_df$year_label[which(plot_df$year==1993)] <- 1993
plot_df$year_label[which(plot_df$year==2023)] <- 2023


# Cc path
ggplot(plot_df) + theme_bw() + 
  scale_color_manual("Country", values = col_scheme) + 
  geom_path(data = plot_df[which(plot_df$name == "Other"),], alpha = 0.5,
            aes(x = C, y = c, color = name, linetype = Forecast)) + 
  geom_path(data = plot_df[which(plot_df$name != "Other"),],
            aes(x = C, y = c, color = name, linetype = Forecast)) + 
  geom_text(data = plot_df[which(plot_df$name != "Other"),],
            aes(x = C, y = c, color = name, label = year_label), show.legend=FALSE, size = 2) +
  xlab("Senescent elongation (C)") + ylab("Senescent rectangularity (c)") + 
  guides(color=guide_legend(ncol=2))
ggsave("figures/countries/summary/Cc_paths.pdf", width = 10, height = 5)

# Bb path
ggplot(plot_df) + theme_bw() + 
  scale_color_manual("Country", values = col_scheme) + 
  geom_path(data = plot_df[which(plot_df$name == "Other"),], alpha = 0.5,
            aes(x = B, y = b, color = name, linetype = Forecast)) + 
  geom_path(data = plot_df[which(plot_df$name != "Other"),],
            aes(x = B, y = b, color = name, linetype = Forecast)) + 
  geom_text(data = plot_df[which(plot_df$name != "Other"),],
            aes(x = B, y = b, color = name, label = year_label), show.legend=FALSE, size = 2) +
  xlab("Infant elongation (B)") + ylab("Infant rectangularity (b)") + 
  guides(color=guide_legend(ncol=2))
ggsave("figures/countries/summary/Bb_paths.pdf", width = 10, height = 5)

# LEh path
ggplot(plot_df) + theme_bw() + 
  scale_color_manual("Country", values = col_scheme) + 
  geom_path(data = plot_df[which(plot_df$name == "Other"),], alpha = 0.5,
            aes(x = LE, y = h, color = name, linetype = Forecast)) + 
  geom_path(data = plot_df[which(plot_df$name != "Other"),],
            aes(x = LE, y = h, color = name, linetype = Forecast)) + 
  geom_text(data = plot_df[which(plot_df$name != "Other"),],
            aes(x = LE, y = h, color = name, label = year_label), show.legend=FALSE, size = 2) +
  xlab("Life expectancy at birth") + ylab("Lifespan equality at birth") + 
  guides(color=guide_legend(ncol=2))
ggsave("figures/countries/summary/Leh_paths.pdf", width = 10, height = 5)

# LstarLE path
ggplot(plot_df) + theme_bw() + 
  scale_color_manual("Country", values = col_scheme) + 
  geom_path(data = plot_df[which(plot_df$name == "Other"),], alpha = 0.5,
            aes(x = LE, y = Lstar, color = name, linetype = Forecast)) + 
  geom_path(data = plot_df[which(plot_df$name != "Other"),],
            aes(x = LE, y = Lstar, color = name, linetype = Forecast)) + 
  geom_text(data = plot_df[which(plot_df$name != "Other"),],
            aes(x = LE, y = Lstar, color = name, label = year_label), show.legend=FALSE, size = 2) +
  xlab("Life expectancy at birth") + ylab("Lifespan") + 
  guides(color=guide_legend(ncol=2)) + 
  scale_x_continuous(limits = c(46,93),breaks=seq(0,130,5)) + 
  scale_y_continuous(limits = c(96,118),breaks=seq(0,130,5)) 
ggsave("figures/countries/summary/LstarLE_paths.pdf", width = 10, height = 4)


# LstarC path
ggplot(plot_df) + theme_bw() + 
  scale_color_manual("Country", values = col_scheme) + 
  geom_path(data = plot_df[which(plot_df$name == "Other"),], alpha = 0.5,
            aes(x = C, y = Lstar, color = name, linetype = Forecast)) + 
  geom_path(data = plot_df[which(plot_df$name != "Other"),],
            aes(x = C, y = Lstar, color = name, linetype = Forecast)) + 
  geom_text(data = plot_df[which(plot_df$name != "Other"),],
            aes(x = C, y = Lstar, color = name, label = year_label), show.legend=FALSE, size = 2) +
  xlab("C") + ylab("Lifespan") + 
  guides(color=guide_legend(ncol=2))
ggsave("figures/countries/summary/LstarC_paths.pdf", width = 10, height = 5)

# Lstar-C vs LE path
ggplot(plot_df) + theme_bw() + 
  scale_color_manual("Country", values = col_scheme) + 
  geom_path(data = plot_df[which(plot_df$name == "Other"),], alpha = 0.5,
            aes(x = LE, y = Lstar-C, color = name, linetype = Forecast)) + 
  geom_path(data = plot_df[which(plot_df$name != "Other"),],
            aes(x = LE, y = Lstar-C, color = name, linetype = Forecast)) + 
  geom_text(data = plot_df[which(plot_df$name != "Other"),],
            aes(x = LE, y = Lstar-C, color = name, label = year_label), show.legend=FALSE, size = 2) +
  xlab("Life expectancy at birth") + ylab("Lifespan - C") + 
  guides(color=guide_legend(ncol=2))
ggsave("figures/countries/summary/LE_vs_Lstar-C_paths.pdf", width = 10, height = 5)




"
Decomposition analysis
"

import_files <- dir("figures/countries/")
import_files <- import_files[which(str_detect(import_files, "decomp_pred.csv"))]

countries_df <- data.frame(matrix(NA,nrow=0,ncol = 50))
names(countries_df) <- c("year", "code", "forecast", "B", "b", "C", "c", "d", "σ", 
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
  country_df <- country_df[,names(countries_df)]
  countries_df <- rbind(countries_df, country_df)
}

countries_df <- merge(countries_df, unique(mort_df[,c("code", "name")]), by = "code")
countries_df <- countries_df[which(countries_df$year > 1918),]

countries_df$DeltaLE_c <- countries_df$Δc*countries_df$LE_c
countries_df$DeltaLE_C <- countries_df$ΔC*countries_df$LE_C

countries_df$year_group <- "pred"
countries_df$year_group[which(countries_df$year <= 1943)] <- "1918-1943"
countries_df$year_group[which(countries_df$year > 1943 & countries_df$year <= 1968)] <- "1943-1968"
countries_df$year_group[which(countries_df$year > 1968 & countries_df$year <= 1993)] <- "1968-1993"
countries_df$year_group[which(countries_df$year > 1993 & countries_df$year <= 2018)] <- "1993-2018"
table(countries_df$year_group)

export_table <- countries_df[which(!is.na(countries_df$DeltaLE_c)),]
export_table <- aggregate(export_table[,c("DeltaLE_c", "DeltaLE_C", "ΔLE_mod")],
                          by = list(name = export_table$name, year = export_table$year_group), 
                          FUN = sum, drop = FALSE)
export_table$propLE_c <- round(export_table$DeltaLE_c/export_table$ΔLE_mod,3)
export_table$propLE_C <- round(export_table$DeltaLE_C/export_table$ΔLE_mod,3)

export_table <- data.frame(pivot_wider(export_table, 
                                       id_cols = c(name), names_from = year, 
                                       names_glue = "{.value}_{year}",
                                       values_from = c(propLE_c, propLE_C)))

export_table <- export_table[,c("name", "propLE_c_1918.1943",  "propLE_C_1918.1943",
                                "propLE_c_1943.1968", "propLE_C_1943.1968",
                                "propLE_c_1968.1993", "propLE_C_1968.1993",
                                "propLE_c_1993.2018", "propLE_C_1993.2018")]

export_table <- merge(export_table, countries_df[which(countries_df$year == 2018),c("name", "C")],
                      by = "name")
export_table$C <- round(export_table$C, 1)

stargazer(as.matrix(export_table), table.placement = "H", column.sep.width = "2")

















# Lstar vs LE
ggplot(plot_df[which(plot_df$code == "USA"),]) + theme_bw() + 
  scale_color_manual("Legend", values = c("red","blue")) + 
  geom_line(aes(x = year, y = Lstar, color = "Lifespan", linetype = Forecast)) + 
  geom_line(aes(x = year, y = LE, color = "Life expectancy", linetype = Forecast)) +
  geom_ribbon(aes(x = year, ymin = LE, ymax = Lstar), fill = "purple", alpha = 0.2) + 
  xlab("Years from birth") + ylab("Year") + 
  guides(color=guide_legend(ncol=1)) 
#ggsave("figures/benchmark/USA_Lstar-LE.pdf", width = 6, height = 4)


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
  

















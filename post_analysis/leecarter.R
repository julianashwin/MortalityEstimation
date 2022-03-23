setwd("/Users/julianashwin/Documents/GitHub/MortalityEstimation/")
rm(list=ls())

require(ggplot2)
require(ggpubr)
require(stringr)
require(tidyr)
require(reshape)
require(reshape2)
require(plm)
require(lfe)
require(stargazer)
require(demography)

"
Define some useful color schemes
"
# Color scheme
col_scheme <- c("Canada" = "pink", "France" = "blue3", "Italy" =  "forestgreen", 
                "United States" = "cornflowerblue", "West Germany" = "darkgoldenrod2", 
                "United Kingdom" = "gray", "Japan"= "red","Best Practice"= "black",
                "Belgium" = "gold3", "Denmark" = "firebrick4", "Finland" = "darkslategray1",
                "Netherlands" = "darkorange1", "Norway" = "deeppink", "Sweden" = "yellow",
                "Switzerland" = "darkorchid3", "Iceland" = "cornsilk3")
line_colors <- scale_color_manual("Country", values = col_scheme)
fill_colors <- scale_fill_manual("Country", values = col_scheme)


par_cols <- c("LE Change" = "black", "H Change" = "black", "b" = "blue", "B" =  "purple", 
              "c" = "red", "C" = "orange", "d" = "green")
par_lines <- scale_color_manual("Parameter", values = par_cols)
par_fills <- scale_fill_manual("Parameter", values = par_cols)



"
Import data and results
"
# Import parameter estimates
parests_df <- read.csv("results_firstdiff/select_country_siler_est_results.csv", stringsAsFactors = FALSE)
#other_parests_df <- read.csv("results_justrw/other_country_siler_est_results.csv", stringsAsFactors = FALSE)
# Shorten USA to United States
parests_df$name[which(parests_df$code == "USA")] <- "United States"
# year = 0 means year should be NA
parests_df$year[which(parests_df$year == 0)] <- NA
#other_parests_df$year[which(other_parests_df$year == 0)] <- NA
# Combine
all_parests_df <- rbind(parests_df)#, other_parests_df)

# Include only countries that have a full sample (plus UK, US and Japan)
full_codes <- names(table(all_parests_df$code)[which(table(all_parests_df$code) == 162)])
full_codes <- c(full_codes, "GBR", "USA", "JPN")
full_parests_df <- all_parests_df[which(all_parests_df$code %in% full_codes),]


# Import mortality data
mort_df <- read.csv("data/clean/all_lifetab.csv", stringsAsFactors = FALSE)

# Create data matrix
code <- "SWE"
country_df <- mort_df[which(mort_df$code == code),]

mx_matrix <- cast(country_df[,c("year", "age", "mx")], age ~ year, value = "mx")
lx_matrix <- cast(country_df[,c("year", "age", "lx")], age ~ year, value = "lx")

country_demdata <- demogdata(data=mx_matrix[,2:ncol(mx_matrix)], 
                             pop = 1000000*lx_matrix[,2:ncol(lx_matrix)],
                             ages=mx_matrix$age, as.numeric(names(mx_matrix[,2:ncol(mx_matrix)])), 
                             type= "mortality",label=code, name="both")


country_LC <- lca(country_demdata)
plot(country_LC)
par(mfrow=c(1,2))
plot(country_demdata,ylim=c(-11,1))
plot(forecast(country_LC,jumpchoice="actual"),ylim=c(-17,1))

country_BMS <- bms(country_demdata, breakmethod="bai")
fcast_BMS <- forecast(country_BMS)
par(mfrow=c(1,1))
plot(fcast_BMS$kt)



names(mx_matrix)






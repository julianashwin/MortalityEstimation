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


"
Import data and results
"
# Import parameter estimates
bp_df <- read.csv("data/clean/bp_5y.csv", stringsAsFactors = FALSE)
bp_df <- bp_df[which(bp_df$year > 1900),]
bp_df$mx_f[which(bp_df$mx_f==0)] <- min(bp_df$mx_f[which(bp_df$mx_f>0)])
# The Siler model oos forecasts
forecasts_dist_df <- read.csv("figures/benchmark/held-out/siler_i2drift_preds_all.csv", stringsAsFactors = F)
forecasts_dist_df <- forecasts_dist_df[which(forecasts_dist_df$parameter == "LE"),]
forecasts_dist_df$Forecast <- "Estimate"
forecasts_dist_df$Forecast[which(forecasts_dist_df$forecasts_dist_df == 1)] <- "Forecast"
forecasts_dist_df$`Estimation Year` <- as.character(forecasts_dist_df$est_year)


forecasts_dist_df <- merge(forecasts_dist_df, bp_df[which(bp_df$age == 0),c("year", "ex_f")],
                      by = "year", all.x = T)

ggplot(forecasts_dist_df) + theme_bw() + 
  geom_ribbon(aes(x = year, ymin = pc15, ymax = pc85, fill = `Estimation Year`,
                  group = interaction(`Estimation Year`, Forecast)), alpha = 0.3) + 
  geom_point(aes(x = year, y = ex_f), shape = 3) + 
  xlab("Year") + ylab("Life expectancy at birth")




forecasts_df <- read.csv("figures/benchmark/held-out/siler_i2drift_LEgrad_oos.csv", stringsAsFactors = F)
forecasts_df$Forecast <- "Estimate"
forecasts_df$Forecast[which(forecasts_df$forecast == 1)] <- "Forecast"
forecasts_df$`Estimation Year` <- as.character(forecasts_df$est_year)

forecasts_df <- merge(forecasts_df, bp_df[,c("year", "age", "mx_f", "lx_f", "ex_f")],
                      by = c("year", "age"), all.x = T)

ggplot(forecasts_df[which(forecasts_df$year %in% c(1903,1963,2018) &
                            forecasts_df$est_year == 2018),]) + theme_bw() + 
  scale_color_gradientn(colours = rainbow(5), name = "Year") +
  geom_line(aes(x = age, y = log(mortality), group = year, color = year)) +
  geom_point(aes(x = age, y = log(mx_f), color = year), shape = 3) 
  
  
  xlab("Year") + ylab("Life expectancy at birth")



"
Cycle through to get out of sample results
"
LC_df <- data.frame(matrix(NA, nrow = 0, ncol = 6))
names(LC_df) <-c("year", "est_year", "age", "mx_LC", "lx_LC", "ex_LC")

est_years <- sort(unique(forecasts_df$est_year))

yy <- est_years[1]
for (yy in est_years){
  mx_matrix <- cast(bp_df[which(bp_df$year <= yy),c("year", "age", "mx_f")], 
                    age ~ year, value = "mx_f")
  mx_matrix[(mx_matrix == 0)]
  pop_matrix <- cast(bp_df[which(bp_df$year <= yy),c("year", "age", "Female")], 
                     age ~ year, value = "Female")
  
  bp_demdata <- demogdata(data=mx_matrix[,2:ncol(mx_matrix)], 
                          pop = pop_matrix[2:ncol(mx_matrix)],
                          ages=mx_matrix$age, 
                          years = as.numeric(names(mx_matrix[,2:ncol(mx_matrix)])), 
                          type= "mortality",label="BP", name="female")
  #plot(bp_demdata)
  bp_lt <- lifetable(bp_demdata, ages = bp_demdata$age)
  # Historical ex
  bp_ex <- data.frame(bp_lt$ex)
  bp_ex$age <- rownames(bp_ex)
  bp_ex <- reshape2::melt(bp_ex, id = "age", variable.name = "year",value.name = "ex_LC")
  bp_ex$year <- as.numeric(str_remove(bp_ex$year, "X"))
  # Historical lx
  bp_lx <- data.frame(bp_lt$lx)
  bp_lx$age <- rownames(bp_lx)
  bp_lx <- reshape2::melt(bp_lx, id = "age", variable.name = "year",value.name = "mx_LC")
  bp_lx$year <- as.numeric(str_remove(bp_lx$year, "X"))
  # Historical mx
  bp_mx <- data.frame(bp_lt$mx)
  bp_mx$age <- rownames(bp_mx)
  bp_mx <- reshape2::melt(bp_mx, id = "age", variable.name = "year",value.name = "mx_LC")
  bp_mx$year <- as.numeric(str_remove(bp_mx$year, "X"))
  
  # Estimate model and forecast
  bp_LC <- lca(bp_demdata)
  bp_LC_fc <- forecast(bp_LC,jumpchoice="fit", h = 6)
  bp_lt_fc <- lifetable(bp_LC_fc)
  # Historical ex
  bp_ex_fc <- data.frame(bp_lt_fc$ex)
  bp_ex_fc$age <- rownames(bp_ex_fc)
  bp_ex_fc <- reshape2::melt(bp_ex_fc, id = "age", variable.name = "year",value.name = "ex_LC")
  bp_ex_fc$year <- as.numeric(str_remove(bp_ex_fc$year, "X"))
  # Historical lx
  bp_lx_fc <- data.frame(bp_lt_fc$lx)
  bp_lx_fc$age <- rownames(bp_lx_fc)
  bp_lx_fc <- reshape2::melt(bp_lx_fc, id = "age", variable.name = "year",value.name = "mx_LC")
  bp_lx_fc$year <- as.numeric(str_remove(bp_lx_fc$year, "X"))
  # Historical mx
  bp_mx_fc <- data.frame(bp_lt_fc$mx)
  bp_mx_fc$age <- rownames(bp_mx_fc)
  bp_mx_fc <- reshape2::melt(bp_mx_fc, id = "age", variable.name = "year",value.name = "mx_LC")
  bp_mx_fc$year <- as.numeric(str_remove(bp_mx_fc$year, "X"))
  
  rbind(bp_mx,bp_mx_fc)
  
  temp_df <- data.frame(year = as.numeric(c(colnames(bp_lt$ex), colnames(bp_lt_fc$ex))),
                        est_year = yy, ex_LC= c(bp_lt$ex[1,], bp_lt_fc$ex[1,]))

  LC_df <- rbind(LC_df, temp_df)

}

forecasts_df <- merge(forecasts_df, LC_df,by = c("year", "est_year"), all.x = T)


"
Compare Lee-Carter and siler model forecasts
"

ggplot(forecasts_df) + theme_bw() + 
  #geom_ribbon(aes(x = year, ymin = pc15, ymax = pc85, fill = `Estimation Year`,
  #                group = interaction(`Estimation Year`, Forecast)), alpha = 0.3) + 
  geom_point(aes(x = year, y = ex_f), shape = 3) +
  geom_line(aes(x = year, y = median, color = `Estimation Year`, linetype = "Siler")) + 
  geom_line(aes(x = year, y = ex_LC, color = `Estimation Year`, linetype = "LC")) + 
  xlab("Year") + ylab("Life expectancy at birth")


ggplot(forecasts_df[which(forecasts_df$year > 1958 & forecasts_df$forecast == 1),]) + theme_bw() + 
  #geom_ribbon(aes(x = year, ymin = pc15, ymax = pc85, fill = `Estimation Year`,
  #                group = interaction(`Estimation Year`, Forecast)), alpha = 0.3) + 
  geom_point(aes(x = year, y = ex_f), shape = 3) +
  geom_line(aes(x = year, y = median + 0.5, color = `Estimation Year`, linetype = "Siler")) + 
  geom_line(aes(x = year, y = ex_LC, color = `Estimation Year`, linetype = "LC")) + 
  xlab("Year") + ylab("Life expectancy at birth")




"
End
"
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
require(qlcMatrix)

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

forecasts_df <- merge(forecasts_df, bp_df[,c("year", "age", "mx_f", "lx_f", "ex_f", "Female")],
                      by = c("year", "age"), all.x = T)

ggplot(forecasts_df[which(forecasts_df$year %in% c(1903,1963,2018) &
                            forecasts_df$est_year == 2018),]) + theme_bw() + 
  scale_color_gradientn(colours = rainbow(5), name = "Year") +
  geom_line(aes(x = age, y = log(mortality), group = year, color = year)) +
  geom_point(aes(x = age, y = log(mx_f), color = year), shape = 3) +
  xlab("Year") + ylab("Life expectancy at birth")



### Create a lifetable from the siler mortality estimates. 
est_years <- sort(unique(forecasts_df$est_year))
siler_df <- data.frame(matrix(NA, nrow = 0, ncol = 4))
names(siler_df) <-c("year", "est_year", "age", "ex_siler")

yy <- est_years[1]
for (yy in est_years){
  mx_matrix <- cast(forecasts_df[which(forecasts_df$est_year == yy),c("year", "age", "mortality")], 
                  age ~ year, value = "mortality")
  pop_matrix <- cast(forecasts_df[which(forecasts_df$est_year == yy),c("year", "age", "Female")], 
                   age ~ year, value = "Female")
  siler_demdata <- demogdata(data=mx_matrix[,2:ncol(mx_matrix)], 
                        pop = pop_matrix[2:ncol(mx_matrix)],
                        ages=mx_matrix$age, 
                        years = as.numeric(names(mx_matrix[,2:ncol(mx_matrix)])), 
                        type= "mortality",label="Siler", name="female")
  siler_lt <- lifetable(siler_demdata, ages = siler_demdata$age)
  # Extract the life expectancy
  siler_ex <- data.frame(siler_lt$ex)
  siler_ex$age <- rownames(siler_ex)
  siler_ex <- reshape2::melt(siler_ex, id = "age", variable.name = "year",value.name = "ex_siler")
  siler_ex$year <- as.numeric(str_remove(siler_ex$year, "X"))
  
  # Add the recalculated life expectancy to the bottom of the df
  siler_ex$est_year <- yy
  siler_df <- rbind(siler_df, siler_ex[names(siler_df)])
}
forecasts_df <- merge(forecasts_df, siler_df,by = c("year", "est_year", "age"), all.x = T)



"
Cycle through to get out of sample results
"
LC_df <- data.frame(matrix(NA, nrow = 0, ncol = 6))
names(LC_df) <-c("year", "est_year", "age", "mx_LC", "lx_LC", "ex_LC")

yy <- est_years[1]
for (yy in est_years){
  mx_matrix <- cast(bp_df[which(bp_df$year <= yy),c("year", "age", "mx_f")], 
                    age ~ year, value = "mx_f")
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
  bp_lx <- reshape2::melt(bp_lx, id = "age", variable.name = "year",value.name = "lx_LC")
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
  bp_lx_fc <- reshape2::melt(bp_lx_fc, id = "age", variable.name = "year",value.name = "lx_LC")
  bp_lx_fc$year <- as.numeric(str_remove(bp_lx_fc$year, "X"))
  # Historical mx
  bp_mx_fc <- data.frame(bp_lt_fc$mx)
  bp_mx_fc$age <- rownames(bp_mx_fc)
  bp_mx_fc <- reshape2::melt(bp_mx_fc, id = "age", variable.name = "year",value.name = "mx_LC")
  bp_mx_fc$year <- as.numeric(str_remove(bp_mx_fc$year, "X"))
  
  temp_df <- merge(rbind(bp_mx,bp_mx_fc), rbind(bp_lx,bp_lx_fc), 
                   by = c("year", "age"))
  temp_df <- merge(temp_df, rbind(bp_ex,bp_ex_fc), 
                   by = c("year", "age"))
  temp_df$est_year <- yy
  
  
  LC_df <- rbind(LC_df, temp_df[names(LC_df)])

}

forecasts_df <- merge(forecasts_df, LC_df,by = c("year", "est_year", "age"), all.x = T)


## Add a corresponding adjustment to the forecast distributions
for (ii in 1:nrow(forecasts_dist_df)){
  yy <- forecasts_dist_df$year[ii]
  yy_est <- forecasts_dist_df$est_year[ii]
  obs <- which(forecasts_df$year == yy & forecasts_df$est_year == yy_est & 
                 forecasts_df$age == 0)
  adjustment <- forecasts_df$LE[obs] - forecasts_df$ex_siler[obs]
  
  forecasts_dist_df$median_adj[ii] <- forecasts_dist_df$median[ii] - adjustment
  forecasts_dist_df$pc975_adj[ii] <- forecasts_dist_df$pc975[ii] - adjustment
  forecasts_dist_df$pc025_adj[ii] <- forecasts_dist_df$pc025[ii] - adjustment
  forecasts_dist_df$pc85_adj[ii] <- forecasts_dist_df$pc85[ii] - adjustment
  forecasts_dist_df$pc15_adj[ii] <- forecasts_dist_df$pc15[ii] - adjustment
  forecasts_dist_df$pc75_adj[ii] <- forecasts_dist_df$pc75[ii] - adjustment
  forecasts_dist_df$pc25_adj[ii] <- forecasts_dist_df$pc25[ii] - adjustment
  
}



"
Compute the forecast errors
"

forecasts_df$n_ahead <- forecasts_df$year - forecasts_df$est_year
forecasts_df$n_ahead[which(forecasts_df$n_ahead <= 0)] <- NA
# Errors for genuinely out-of-sample predictions
oos_obs <- which(forecasts_df$year > forecasts_df$est_year)
forecasts_df$LC_fe <- NA
forecasts_df$LC_fe[oos_obs] <- forecasts_df$ex_f[oos_obs] - forecasts_df$ex_LC[oos_obs]
forecasts_df$siler_fe[oos_obs] <- forecasts_df$ex_f[oos_obs] - forecasts_df$ex_siler[oos_obs]
# Squared errors
forecasts_df$siler_fe2 <- forecasts_df$siler_fe^2
forecasts_df$LC_fe2 <- forecasts_df$LC_fe^2

ggplot(forecasts_df[which(forecasts_df$age == 0),]) + theme_bw() + 
  geom_bar(aes(x = n_ahead-1, y = siler_fe2, fill = "Siler"), 
           stat = "summary", fun = mean, width = 2) +
  geom_bar(aes(x = n_ahead+1, y = LC_fe2, fill = "Lee-Carter"), 
           stat = "summary", fun = mean, width = 2) +
  xlab("Years ahead") + ylab("MSE")



ggplot(forecasts_df[which(forecasts_df$age == 0),]) + theme_bw() + 
  geom_line(aes(x = year, y = siler_fe, color = `Estimation Year`, linetype = "Siler")) + 
  geom_line(aes(x = year, y = LC_fe, color = `Estimation Year`, linetype = "LC")) + 
  geom_hline(yintercept=0, linetype="dashed", color = "black") +
  xlab("Year") + ylab("Forecast errors")


t.test(forecasts_df$siler_fe[which(forecasts_df$age == 0)], na.rm = T)
t.test(forecasts_df$LC_fe[which(forecasts_df$age == 0)], na.rm = T)
mean(forecasts_df$LC_fe[which(forecasts_df$age == 0)]^2, na.rm = T)
mean(forecasts_df$siler_fe[which(forecasts_df$age == 0)]^2, na.rm = T)

lm(siler_fe^2 ~ n_ahead, forecasts_df)
lm(LC_fe^2 ~ n_ahead, forecasts_df)

"
Compare Lee-Carter and siler model forecasts
"

sum((forecasts_df$mx_f[which(forecasts_df$est_year == 2018)] - 
  forecasts_df$mx_LC[which(forecasts_df$est_year == 2018)])^2, na.rm = T)
sum((forecasts_df$mx_f[which(forecasts_df$est_year == 2018)] - 
       forecasts_df$mortality[which(forecasts_df$est_year == 2018)])^2, na.rm = T)

## Compare some actual mortality curves
ggplot(forecasts_df[which(forecasts_df$year %in% c(1903,1963,2018, 2048) &
                            forecasts_df$est_year == 2018),]) + theme_bw() + 
  scale_color_gradientn(colours = rainbow(5), name = "Year") +
  geom_line(aes(x = age, y = log(mortality), group = year, color = year, linetype = "Siler")) +
  geom_line(aes(x = age, y = log(mx_LC), group = year, color = year, linetype = "Lee-Carter")) +
  geom_point(aes(x = age, y = log(mx_f), color = year), shape = 3) +
  xlab("Year") + ylab("Life expectancy at birth")


## Out of sample LE forecasts
# Unadjusted
ggplot(forecasts_df[which(forecasts_df$age == 0),]) + theme_bw() + 
  geom_ribbon(data = forecasts_dist_df[which(forecasts_dist_df$forecast ==1),], 
              aes(x = year, ymin = pc15, ymax = pc85, fill = `Estimation Year`,
                  group = interaction(`Estimation Year`, Forecast)), alpha = 0.3) +
  geom_line(aes(x = year, y = LE, color = `Estimation Year`)) + 
  geom_point(aes(x = year, y = ex_f), shape = 3) +
  xlab("Year") + ylab("Life expectancy at birth")


ggplot(forecasts_df[which(forecasts_df$age == 0),]) + theme_bw() + 
  geom_ribbon(data = forecasts_dist_df[which(forecasts_dist_df$forecast ==1),], 
              aes(x = year, ymin = pc15_adj, ymax = pc85_adj, fill = `Estimation Year`,
                  group = interaction(`Estimation Year`, Forecast)), alpha = 0.3) +
  geom_line(aes(x = year, y = ex_siler, color = `Estimation Year`)) + 
  geom_point(aes(x = year, y = ex_f), shape = 3) +
  xlab("Year") + ylab("Life expectancy at birth")


ggplot(forecasts_df[which(forecasts_df$year > forecasts_df$est_year  & 
                            forecasts_df$age == 0),]) + theme_bw() + 
  #geom_ribbon(aes(x = year, ymin = pc15, ymax = pc85, fill = `Estimation Year`,
  #                group = interaction(`Estimation Year`, Forecast)), alpha = 0.3) + 
  geom_line(aes(x = year, y = ex_siler, color = `Estimation Year`, linetype = "Siler")) + 
  geom_line(aes(x = year, y = ex_LC, color = `Estimation Year`, linetype = "LC")) + 
  geom_point(aes(x = year, y = ex_f), shape = 3) +
  xlab("Year") + ylab("Life expectancy at birth")

ggplot() 




"
End
"
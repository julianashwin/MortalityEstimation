setwd("/Users/julianashwin/Documents/GitHub/MortalityEstimation/")
rm(list=ls())

require(ggplot2)
require(ggpubr)
require(HMDHFDplus)
require(stringr)
require(plyr)
require(zoo)
require(reshape2)
require(tidyr)



# Colour scheme
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
Import and prepare data
"
# Import Best practice gradients
folder <- "figures/benchmark/"
LEgrad_df <- read.csv(paste0(folder,"siler_i2drift_LEgrads.csv"), stringsAsFactors = FALSE)
LEgrad_df$Forecast <- "Estimate"
LEgrad_df$Forecast[which(LEgrad_df$year > 2020)] <- "Forecast"
LEgrad_df$code <- "Best Practice"

"
Quickly add some parameter estimates to best practice
"
bp_df <- LEgrad_df
bp_pars_df <- read.csv("figures/benchmark/siler_i2drift_preds.csv", stringsAsFactors = F)
bp_pars_df <- bp_pars_df[which(bp_pars_df$parameter %in% c("b","B","c","C","d","Lstar_99p9", 
                                                           "Lstar_99","Lstar_95","Lstar_90")),]
bp_pars_df <- data.frame(pivot_wider(bp_pars_df, id_cols = c(year), names_from = parameter, 
                                     values_from = c(median)))
bp_df <- merge(bp_df, bp_pars_df, by = "year")

bp_mort_df <- read.csv("data/clean/bp_5y.csv", stringsAsFactors = FALSE)
bp_df <- merge(bp_df, bp_mort_df[,c("year", "age", "mx_f", "lx_f", "ex_f", "Hx_f",
                                   "l_90_f","l_95_f","l_99_f","l_99p9_f",
                                   "Female")], 
                by = c("year", "age"), all.x = T)

rm(bp_pars_df,bp_mort_df)

# Import gradients for each country
import_files <- dir("figures/countries/")
import_files <- import_files[which(str_detect(import_files, "_LEgrads.csv"))]
for (ii in 1:length(import_files)){
  filename <- import_files[ii]
  country_df <- read.csv(paste0("figures/countries/", filename), stringsAsFactors = F)
  country_df$Forecast <- "Estimate"
  country_df$Forecast[which(country_df$year > 2020)] <- "Forecast"
  country_df$code <- str_remove(filename, "_i2_LEgrads.csv")
  
  if (is.null(country_df$mortality)){
    country_df$mortality <- NA
    country_df$survival <- NA
  }
  
  country_df <- country_df[,names(LEgrad_df)]
  LEgrad_df <- rbind(LEgrad_df, country_df)
}
rm(country_df)



# Merge in full names
mort_df <- read.csv("data/clean/all_lifetab_5y.csv", stringsAsFactors = FALSE)
all_df <- merge(LEgrad_df, unique(mort_df[,c("code","name")]), by = c("code"), all.x = T)
all_df <- merge(all_df, mort_df[,c("code", "year", "age", "mx_f", "lx_f", "ex_f", "Hx_f",
                                   "l_90_f","l_95_f","l_99_f","l_99p9_f",
                                   "Female", "best_practice")], 
                   by = c("code", "year", "age"), all.x = T)

# Tidy up a little bit
all_df$code <- str_replace(all_df$code, "NZL_NM", "NZL")
all_df <- all_df[which(all_df$year > 1900),]
all_df$plot_name <- all_df$name
all_df$plot_name[which(!(all_df$name %in% names(col_scheme)))] <- "Other"

all_df$year_label <- NA
all_df$year_label[which(all_df$year==1903)] <- 1903
all_df$year_label[which(all_df$year==1933)] <- 1933
all_df$year_label[which(all_df$year==1963)] <- 1963
all_df$year_label[which(all_df$year==1993)] <- 1993
all_df$year_label[which(all_df$year==2023)] <- 2023

## Add in parameter estimates
econ_panel_df <- read.csv("data/results/siler_econ_panel.csv", stringsAsFactors = FALSE)

all_df <- merge(all_df, econ_panel_df[,c("code", "year", "B", "b", "C", "c", "d",
                                         "Lstar_99p9","Lstar_99","Lstar_95","Lstar_90")],
                by = c("code", "year"), all.x = T)
all_df <- all_df[order(all_df$code, all_df$year, all_df$age),]
all_df <- all_df[which(all_df$code != "Best Practice"),]


rm(econ_panel_df, mort_df, LEgrad_df)


"
Gradients for C and c
"

bp_leC_plt <- ggplot(bp_df) + theme_bw() +
  geom_line(aes(x = age, y = LE_Cs, group = year, color = year, linetype = Forecast)) +
  geom_hline(yintercept=0, linetype="dashed", color = "black") +
  scale_color_gradientn(colours = rainbow(5), name = "Year") + 
  xlab("Age") + ylab("Gradient") + 
  ggtitle(expression(Gradient~of~LE~wrt~C[t]))
bp_lec_plt <- ggplot(bp_df) + theme_bw() +
  geom_line(aes(x = age, y = LE_cs, group = year, color = year, linetype = Forecast)) +
  geom_hline(yintercept=0, linetype="dashed", color = "black") +
  scale_color_gradientn(colours = rainbow(5), name = "Year") + 
  xlab("Age") + ylab("Gradient") + 
  ggtitle(expression(Gradient~of~LE~wrt~c[t]))
bp_lecC_plt <- ggplot(bp_df) + theme_bw() +
  geom_line(aes(x = age, y = LE_cC, group = year, color = year, linetype = Forecast)) +
  geom_hline(yintercept=0, linetype="dashed", color = "black") +
  scale_color_gradientn(colours = rainbow(5), name = "Year") + 
  xlab("Age") + ylab("Cross-derivative") + 
  ggtitle(expression(Cross~derivative~of~LE~wrt~c[t]~and~C[t]))
ggarrange(bp_lec_plt, bp_leC_plt, nrow = 1, ncol = 2, common.legend = TRUE, 
          legend = "right")
ggsave(paste0("figures/benchmark/LEgrads_cC.pdf"), width = 10, height = 4)
rm(bp_leC_plt,bp_lec_plt,bp_lecC_plt)


bp_hC_plt <- ggplot(bp_df[bp_df$age <111,]) + theme_bw() +
  geom_line(aes(x = age, y = h_Cs, group = year, color = year, linetype = Forecast)) +
  geom_hline(yintercept=0, linetype="dashed", color = "black") +
  scale_color_gradientn(colours = rainbow(5), name = "Year") + 
  xlab("Age") + ylab("Gradient") + 
  ggtitle(expression(Gradient~of~h~wrt~C[t]))
bp_hc_plt <- ggplot(bp_df[bp_df$age <111,]) + theme_bw() +
  geom_line(aes(x = age, y = h_cs, group = year, color = year, linetype = Forecast)) +
  geom_hline(yintercept=0, linetype="dashed", color = "black") +
  scale_color_gradientn(colours = rainbow(5), name = "Year") + 
  xlab("Age") + ylab("Gradient") + 
  ggtitle(expression(Gradient~of~h~wrt~c[t]))
bp_hcC_plt <- ggplot(bp_df[bp_df$age <111,]) + theme_bw() +
  geom_line(aes(x = age, y = h_cC, group = year, color = year, linetype = Forecast)) +
  geom_hline(yintercept=0, linetype="dashed", color = "black") +
  scale_color_gradientn(colours = rainbow(5), name = "Year") + 
  xlab("Age") + ylab("Cross-derivative") + ylim(c(-3,0.36)) +
  ggtitle(expression(Cross~derivative~of~h~wrt~c[t]~and~C[t]))
ggarrange(bp_hc_plt, bp_hC_plt, nrow = 1, ncol = 2, common.legend = TRUE, 
          legend = "right")
ggsave(paste0("figures/benchmark/hgrads_cC.pdf"), width = 10, height = 4)
rm(bp_hC_plt,bp_hc_plt,bp_hcC_plt)



"
Best practice gradients of LE(0) and h(0) over time
"

plot_df <- bp_df[which(bp_df$age == 0 & bp_df$Forecast == "Estimate"),]
## Gradients
LEgrad_plt <- ggplot(plot_df) + theme_bw() + xlab("Year") +
  geom_hline(yintercept=0, color = "gray") +
  geom_line(aes(x = year, y = LE_Cs), linetype = "solid") +
  geom_line(aes(x = year, y = LE_cs/max(LE_cs)), linetype = "dashed") +
  scale_y_continuous(expression(LE[C]~at~birth~(solid)), sec.axis = sec_axis(
    ~ . *max(plot_df$LE_cs), name = expression(LE[c]~at~birth~(dashed))))
hgrad_plt <- ggplot(plot_df) + theme_bw() + xlab("Year") +
  geom_hline(yintercept=0, color = "gray") +
  geom_line(aes(x = year, y = h_Cs), linetype = "solid") +
  geom_line(aes(x = year, y = h_cs*(max(h_Cs)/max(h_cs))), linetype = "dashed") +
  scale_y_continuous(expression(h[C]~at~birth~(solid)), sec.axis = sec_axis(
    ~ . *(max(plot_df$h_Cs)/max(plot_df$h_cs)), name = expression(h[c]~at~birth~(dashed))))
ggarrange(LEgrad_plt, hgrad_plt, nrow = 1, ncol=2, legend = FALSE)
ggsave("figures/benchmark/LEgrad_hgrad_birth.pdf", width = 8, height = 4)
rm(plot_df, LEgrad_plt,hgrad_plt)
## Elasticities for C and c: e_C = (dLE/dc)(c/LE)
plot_df <- bp_df[which(bp_df$age == 0 & bp_df$Forecast == "Estimate"),]
LEelas_plt <- ggplot(plot_df) + theme_bw() + xlab("Year") +
  scale_linetype_manual("Elasticity wrt ", values = c("C" = "solid",
                                             "c" = "dashed")) + 
  geom_hline(yintercept=0, color = "gray") +
  geom_line(aes(x = year, y = LE_Cs*(C/LE), linetype = "C")) +
  geom_line(aes(x = year, y = LE_cs*(c/LE), linetype = "c")) + 
  ylab("Elasticity of Life Expectancy at birth")
# h elasticity
helas_plt <- ggplot(plot_df) + theme_bw() + xlab("Year") +
  scale_linetype_manual("Elasticity wrt ", values = c("C" = "solid",
                                                      "c" = "dashed")) + 
  geom_hline(yintercept=0, color = "gray") +
  geom_line(aes(x = year, y = h_Cs*(C/h), linetype = "C")) +
  geom_line(aes(x = year, y = h_cs*(c/h), linetype = "c")) + 
  ylab("Elasticity of Lifespan Equality at birth") 
ggarrange(LEelas_plt, helas_plt, nrow = 1, ncol=2, common.legend = TRUE, legend = "right")
ggsave("figures/benchmark/LEelas_helas_birth.pdf", width = 8, height = 4)
rm(plot_df,LEelas_plt,helas_plt)




"
LE(0) and LE(80) versus L*-C
"
extra_obs <- all_df[which((all_df$year == 2018 & !(all_df$code %in% c("NZL","RUS")))|
                            (all_df$year == 2013 & all_df$code %in% c("NZL","RUS"))),]
extra_obs$Forecast <- "Forecast"
plot_df <- rbind(all_df[which(all_df$code != "Best Practice"),], extra_obs)
plot_df <- plot_df[order(plot_df$code, plot_df$year),]
plot_df <- plot_df[which(plot_df$age %in% c(0, 80)),
                     c("code", "name", "plot_name", "year", "year_label", "age", "Forecast", "LE", "Lstar", "c", "C")]
plot_df <- data.frame(pivot_wider(plot_df, id_cols = c(code, name, plot_name, year, year_label, Forecast, Lstar, c, C) , 
                                  names_from = age, names_glue = "{.value}_{age}",values_from = c(LE)))
### L*-C vs LE
LE0_plt <- ggplot(plot_df) + theme_bw() + 
  scale_color_manual("Country", values = col_scheme) + 
  geom_path(data = plot_df[which(plot_df$plot_name == "Other"),], alpha = 0.5,
            aes(x = LE_0, y = Lstar-C, group = interaction(code, Forecast), 
                color = plot_name, linetype = Forecast)) + 
  geom_path(data = plot_df[which(plot_df$plot_name != "Other"),],
            aes(x = LE_0, y = Lstar-C, group = interaction(code, Forecast), 
                color = plot_name, linetype = Forecast)) + 
  geom_text(data = plot_df[which(plot_df$plot_name != "Other"),],
            aes(x = LE_0, y = Lstar-C, color = plot_name, label = year_label), 
            show.legend=FALSE, size = 2) +
  xlab("LE(0)") + ylab("L*-C") + guides(color=guide_legend(ncol=1)) + 
  scale_x_continuous(limits = c(45,95),breaks=seq(0,100,5))+ 
  scale_y_continuous(limits = c(12,30),breaks=seq(0,100,5))
LE80_plt <- ggplot(plot_df) + theme_bw() + 
  scale_color_manual("Country", values = col_scheme) + 
  geom_path(data = plot_df[which(plot_df$plot_name == "Other"),], alpha = 0.5,
            aes(x = LE_80, y = Lstar-C, group = interaction(code, Forecast), 
                color = plot_name, linetype = Forecast)) + 
  geom_path(data = plot_df[which(plot_df$plot_name != "Other"),],
            aes(x = LE_80, y = Lstar-C, group = interaction(code, Forecast), 
                color = plot_name, linetype = Forecast)) + 
  geom_text(data = plot_df[which(plot_df$plot_name != "Other"),],
            aes(x = LE_80, y = Lstar-C, color = plot_name, label = year_label), 
            show.legend=FALSE, size = 2) +
  xlab("LE(80)") + ylab("L*-C") + guides(color=guide_legend(ncol=2)) + 
  scale_x_continuous(limits = c(4,17.5),breaks=seq(0,100,5)) + 
  scale_y_continuous(limits = c(12,30),breaks=seq(0,100,5))
Lstar_plt <- ggplot(plot_df) + theme_bw() + 
  scale_color_manual("Country", values = col_scheme) + 
  geom_path(data = plot_df[which(plot_df$plot_name == "Other"),], alpha = 0.5,
            aes(x = LE_80, y = Lstar-LE_80-80, group = interaction(code, Forecast), 
                color = plot_name, linetype = Forecast)) + 
  geom_path(data = plot_df[which(plot_df$plot_name != "Other"),],
            aes(x = LE_80, y = Lstar-LE_80-80, group = interaction(code, Forecast), 
                color = plot_name, linetype = Forecast)) + 
  geom_text(data = plot_df[which(plot_df$plot_name != "Other"),],
            aes(x = LE_80, y = Lstar-LE_80-80, color = plot_name, label = year_label), 
            show.legend=FALSE, size = 2) +
  xlab("LE(80)") + ylab("L*-LE(80)-80") + guides(color=guide_legend(ncol=2)) + 
  scale_x_continuous(limits = c(4,17.5),breaks=seq(0,100,5)) + 
  scale_y_continuous(limits = c(12,30),breaks=seq(0,100,5)) 
ggarrange(LE0_plt, LE80_plt, Lstar_plt, nrow = 1, ncol=3, common.legend = TRUE, 
          legend = "right")
ggsave("figures/benchmark/LE0LE80_Lstar_C.pdf", width = 10, height = 4)
rm(plot_df,extra_obs,LE0_plt,LE80_plt,Lstar_plt)





"
Model lifespan vs data lifespan
"

plot_df <- bp_df[which(bp_df$age == 0),]

ggplot(plot_df) + theme_bw() + 
  geom_line(aes(x = year, y = Lstar_99p9, color = "1) 99.9%")) + 
  geom_point(aes(x = year, y = l_99p9_f, color = "1) 99.9%")) + 
  geom_line(aes(x = year, y = Lstar_99, color = "2) 99%")) + 
  geom_point(aes(x = year, y = l_99_f, color = "2) 99%")) + 
  geom_line(aes(x = year, y = Lstar_95, color = "3) 95%")) + 
  geom_point(aes(x = year, y = l_95_f, color = "3) 95%")) + 
  geom_line(aes(x = year, y = Lstar_90, color = "4) 90%")) + 
  geom_point(aes(x = year, y = l_90_f, color = "4) 90%")) + 
  ylab("Lifespan") + xlab("Year")
ggsave("figures/benchmark/Lstar_defs_modelfit.pdf", width = 6, height = 4)
rm(plot_df)
  






"
Out of sample forecasts
"

forecasts_df <- read.csv("figures/benchmark/held-out/siler_i2drift_preds_all.csv", stringsAsFactors = F)
forecasts_df <- forecasts_df[which(forecasts_df$parameter == "LE"),]
forecasts_df$Forecast <- "Estimate"
forecasts_df$Forecast[which(forecasts_df$forecast == 1)] <- "Forecast"
forecasts_df$`Estimation Year` <- as.character(forecasts_df$est_year)


bp_data_df <- mort_df[which(mort_df$best_practice == 1 & mort_df$age == 0 & 
                              mort_df$year > 1900),]
plot_df <- bp_df[which(bp_df$age == 0),]

ggplot() + theme_bw() + 
  geom_ribbon(data=forecasts_df, aes(x = year, ymin = pc15, ymax = pc85, fill = `Estimation Year`,
              group = interaction(`Estimation Year`, Forecast)), alpha = 0.3) + 
  geom_point(data=plot_df, aes(x = year, y = ex_f), shape = 3) + 
  xlab("Year") + ylab("Life expectancy at birth")
ggsave("figures/benchmark/LE_oos_forecasts.pdf", width = 6, height = 3)
rm(forecasts_df, plot_df)


"
Remaining LE by age group
"

plot_df <- bp_df[which(bp_df$age%% 10 ==0),]
plot_df$Age <- factor(plot_df$age, 
                     levels = unique(plot_df$age)[order(as.numeric(plot_df$age))])
extra_obs <- plot_df[which(plot_df$year == 2018),]
extra_obs$Forecast <- "Forecast"
plot_df <- rbind(plot_df, extra_obs)

# Plot
ggplot(plot_df[which(plot_df$age < 100),]) + theme_bw() + 
  geom_point(aes(x = year, y = ex_f, color = Age)) + 
  geom_line(aes(x = year, y = LE, color = Age, linetype = Forecast)) + 
  xlab("Year") + ylab("Remaining Life Expectancy") +guides(color=guide_legend(ncol=2))
ggsave("figures/benchmark/rLE_byage.pdf", width = 6, height = 3)
rm(plot_df, extra_obs)




"
Jean Calment exercise
"
JC_age <- 122+(164/365)



for (code in unique(all_df$code)){
  country_df <- all_df[which(all_df$code == code & all_df$year >= 2018),]
  current_df <- country_df[which(country_df$year == 2018),c("age", "Female")]
  current_df$Female[which(is.na(current_df$Female))] <- 0
  current_df$pop_2022 <- current_df$Female
  years <- 2023:max(country_df$year)
  
  N_age <- nrow(current_df)
  
  temp_df <- data.frame(year = years)
  temp_df$closest <- c(rep(2023,3), rep(2028,5), rep(2033,5), rep(2038,5),
                       rep(2043,5), rep(2048,3))
  
  for(yy in 1:nrow(temp_df)){
    yy_closest <- temp_df$closest[yy]
    
    current_mort <- country_df[which(country_df$year == yy_closest),c("age", "mortality")]
    names(current_mort)[2] <- paste0("m_", temp_df$year[yy])
    current_df <- merge(current_df, current_mort, by = "age")
    current_df[,paste0("pop_", temp_df$year[yy])] <- c(0, 
      current_df[1:(N_age-1),paste0("pop_", temp_df$year[yy]-1)]*
      (1 - current_df[1:(N_age-1),paste0("m_", temp_df$year[yy])]))
    
  }
  
  
  
  
}


























"
End of script
"
























"
Spare
"

"
Save data on Lstar, C and LEa
"
plot_df <- LEgrad_df[which(LEgrad_df$year > 1900),]
extra_obs <- plot_df[which((plot_df$year == 2018 & !(plot_df$code %in% c("NZL","RUS")))|
                             (plot_df$year == 2013 & plot_df$code %in% c("NZL","RUS"))),]
extra_obs$Forecast <- "Forecast"
plot_df <- rbind(plot_df, extra_obs)
plot_df$name[which(!(plot_df$name %in% names(col_scheme)))] <- "Other"
# Include some years as labels
plot_df$year_label <- NA
plot_df$year_label[which(plot_df$year==1903)] <- 1903
plot_df$year_label[which(plot_df$year==1933)] <- 1933
plot_df$year_label[which(plot_df$year==1963)] <- 1963
plot_df$year_label[which(plot_df$year==1993)] <- 1993
plot_df$year_label[which(plot_df$year==2023)] <- 2023
plot_df <- plot_df[order(plot_df$code, plot_df$year),]

plot_df$LE0 <- plot_df$LE
plot_df$LE0[which(plot_df$age != 0)] <- NA
plot_df <- plot_df[order(plot_df$code, plot_df$year, plot_df$age),]
plot_df$LE0 <- na.locf(plot_df$LE0, fromLast = FALSE)

temp_df <- plot_df[plot_df$age %in% c(0, 80),]
temp_df <- merge(temp_df, econ_panel_df[,c("code", "year", "c", "C")], by = c("code", "year"),
                 all.x = T)
temp_df <- temp_df[,c("code", "year", "age", "LE", "Lstar", "c", "C")]
temp_df <- temp_df[!duplicated(temp_df),]
temp_df <- data.frame(pivot_wider(temp_df, id_cols = c(code, year,Lstar, c, C) , 
                                  names_from = age, names_glue = "{.value}_{age}",
                                  values_from = c(LE)))
temp_df <- data.frame(pivot_wider(temp_df, id_cols = c( year) , 
                                  names_from = code, names_glue = "{code}_{.value}",
                                  values_from = c( Lstar, c, C, LE_0, LE_80)))
write.csv(temp_df, "data/results/Lstar_C_LEa.csv", row.names = FALSE)
rm(temp_df)










"
Path of L vs h at different ages
"
plot_df <- rbind(plot_df, extra_obs)
# Include some years as labels
plot_df$LE0 <- plot_df$LE
plot_df$LE0[which(plot_df$age != 0)] <- NA
plot_df <- plot_df[order(plot_df$code, plot_df$year, plot_df$age),]
plot_df$LE0 <- na.locf(plot_df$LE0, fromLast = FALSE)


# LE(0) vs h(a) for different ages 
plot_df1 <- plot_df[plot_df$age %in% c(0, 50,70,80),]
ggplot(plot_df1) + theme_bw() + 
  scale_color_manual("Country", values = col_scheme) + 
  geom_path(data = plot_df1[which(plot_df1$name == "Other"),], alpha = 0.5,
            aes(x = LE0, y = h, group = interaction(code, Forecast), 
                color = name, linetype = Forecast)) + 
  geom_path(data = plot_df1[which(plot_df1$name != "Other"),],
            aes(x = LE0, y = h, group = interaction(code, Forecast), 
                color = name, linetype = Forecast)) + 
  geom_text(data = plot_df1[which(plot_df1$name != "Other"),],
            aes(x = LE0, y = h, color = name, label = year_label), show.legend=FALSE, size = 1) +
  xlab("LE(0)") + ylab("h(age)") + guides(color=guide_legend(ncol=2)) +
  facet_grid(.~age)
ggsave("figures/countries/summary/LE0_vs_ha_international.pdf", width = 12, height = 4)

# LE(a) vs h(a) for different ages 
plot_df1 <- plot_df[plot_df$age %in% c(0, 50,70,80),]
ggplot(plot_df1) + theme_bw() + 
  scale_color_manual("Country", values = col_scheme) + 
  geom_path(data = plot_df1[which(plot_df1$name == "Other"),], alpha = 0.5,
            aes(x = LE, y = h, group = interaction(code, Forecast), 
                color = name, linetype = Forecast)) + 
  geom_path(data = plot_df1[which(plot_df1$name != "Other"),],
            aes(x = LE, y = h, group = interaction(code, Forecast), 
                color = name, linetype = Forecast)) + 
  geom_text(data = plot_df1[which(plot_df1$name != "Other"),],
            aes(x = LE, y = h, color = name, label = year_label), show.legend=FALSE, size = 1) +
  xlab("LE(age)") + ylab("h(age)") + guides(color=guide_legend(ncol=2)) +
  facet_grid(.~age, scales="free")
ggsave("figures/countries/summary/LEa_vs_ha_international.pdf", width = 12, height = 4)


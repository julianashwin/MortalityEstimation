setwd("/Users/julianashwin/Documents/GitHub/MortalityEstimation/")
rm(list=ls())

require(ggplot2)
require(HMDHFDplus)
require(stringr)
require(plyr)


"
Import data for all available countries
"
## Mortality, survival and life expectancy
# Initialise dataframe to fill with each country
var_names <- c("code", "year", "age", "mx", "lx", "ex")
lifetab_df <- as.data.frame(matrix(NA, nrow = 0, ncol = length(var_names)), var_names)
# Identify which directory to import from
import_dir <- "data/raw/bltper_1x1/"
files <- dir(import_dir)
for (file in files){
  print(file)
  country_code <- str_split(file, "\\.")[[1]][1]
  country_df <- as.data.frame(read.table(paste0(import_dir,file), header = TRUE, skip = 2))
  country_df <- cbind(code = country_code, country_df)  
  country_df$year <- country_df$Year
  country_df$age <- country_df$Age
  country_df$mx[which(country_df$mx == ".")] <- NA
  country_df$mx <- as.numeric(country_df$mx)
  country_df$lx[which(country_df$lx == ".")] <- NA
  country_df$lx <- as.numeric(country_df$lx)/100000
  country_df$ex[which(country_df$ex == ".")] <- NA
  country_df$ex <- as.numeric(country_df$ex)
  
  lifetab_df <- rbind(lifetab_df, country_df[,var_names])
}
rm(country_df)

lifetab_df <- lifetab_df[which(lifetab_df$age != "110+"),]
lifetab_df$age <- as.numeric(lifetab_df$age)
lifetab_df <- lifetab_df[which(!is.na(lifetab_df$mx)),]

## Population
# Initialise dataframe to fill with each country
pop_names <- c("code", "year", "age", "Male", "Female", "Total")
pop_df <- as.data.frame(matrix(NA, nrow = 0, ncol = length(pop_names)), pop_names)
# Identify which directory to import from
import_dir <- "data/raw/Population/"
files <- dir(import_dir)
for (file in files){
  print(file)
  country_code <- str_split(file, "\\.")[[1]][1]
  country_df <- as.data.frame(read.table(paste0(import_dir,file), header = TRUE, skip = 2))
  country_df <- cbind(code = country_code, country_df)  
  country_df$year <- country_df$Year
  country_df$age <- country_df$Age
  country_df$Female[which(country_df$Female == ".")] <- NA
  country_df$Female <- as.numeric(country_df$Female)
  country_df$Male[which(country_df$Male == ".")] <- NA
  country_df$Male <- as.numeric(country_df$Male)
  country_df$Total[which(country_df$Total == ".")] <- NA
  country_df$Total <- as.numeric(country_df$Total)
  
  pop_df <- rbind(pop_df, country_df[,pop_names])
}
rm(country_df, file, files, country_code)

pop_df <- pop_df[which(pop_df$age != "110+"),]
pop_df$age <- as.numeric(pop_df$age)
pop_df <- pop_df[which(!is.na(pop_df$Total)),]



lifetab_df <- merge(lifetab_df, pop_df, by = c("code", "year", "age"))
lifetab_df <- lifetab_df[with(lifetab_df, order(code, year, age)), ]











# Plot for most recent period
plt_codes <- c("BEL", "DNK", "ENW", "FRA", "ITA", "NOR", "NZL_NM", "CHE", "FIN", "ISL", "NLD", "SWE")
all2010_plt <- ggplot(lifetab_df[which(lifetab_df$year == "2010" & 
                                         lifetab_df$code %in% plt_codes),], aes(x = age, y = mx)) + 
  geom_line(aes(group = code, color = code)) + scale_color_discrete(name = "Country") + 
  xlab("Age") + ylab("") + theme_bw() + ylim(c(0,0.85)) +
  ggtitle("2010")
all1901_plt <- ggplot(lifetab_df[which(lifetab_df$year == "1901" & 
                                         lifetab_df$code %in% plt_codes),], aes(x = age, y = mx)) + 
  geom_line(aes(group = code, color = code)) + scale_color_discrete(name = "Country") + 
  xlab("Age") + ylab("Mortality Rate") + theme_bw() + ylim(c(0,0.85)) +
  ggtitle("1901")
ggarrange(all1901_plt,all2010_plt, nrow = 1, ncol=2, common.legend = TRUE)
ggsave("figures/data/all_1990_2010_mortality_rates.pdf", width = 6, height = 4)

all2010_plt <- ggplot(lifetab_df[which(lifetab_df$year == "2010" & 
                                         lifetab_df$code %in% plt_codes),], aes(x = age, y = Total)) + 
  geom_line(aes(group = code, color = code)) + scale_color_discrete(name = "Country") + 
  xlab("Age") + ylab("") + theme_bw() + 
  ggtitle("2010")
all1901_plt <- ggplot(lifetab_df[which(lifetab_df$year == "1901" & 
                                         lifetab_df$code %in% plt_codes),], aes(x = age, y = Total)) + 
  geom_line(aes(group = code, color = code)) + scale_color_discrete(name = "Country") + 
  xlab("Age") + ylab("Population") + theme_bw() + 
  ggtitle("1901")
ggarrange(all1901_plt,all2010_plt, nrow = 1, ncol=2, common.legend = TRUE)
ggsave("figures/data/all_1990_2010_population.pdf", width = 6, height = 4)

rm(all1901_plt, all2010_plt)

"
Get five year average
"
# Collapse down to averages over 5 year intervals
lifetab_5y_df <- lifetab_df
lifetab_5y_df$years <- as.character(cut(lifetab_5y_df$year, seq(1750, 2020, 5)) )
lifetab_5y_df <- ddply(lifetab_5y_df, .(code, years, age), numcolwise(mean))

# Make sure year represents the mid-point of the window
lifetab_5y_df[,c("year_lower", "year_upper")] <- 
  do.call(rbind, str_split(str_remove(str_remove(lifetab_5y_df$years, "\\("), "\\]"), ","))
lifetab_5y_df$year <- as.numeric(lifetab_5y_df$year_lower) + 3

# Plot an example
enw_plt <- ggplot(lifetab_5y_df[which(lifetab_5y_df$code == "ENW" & lifetab_5y_df$year > 1900),], 
       aes(x = age, y = mx)) + scale_color_continuous(name = "Year") + theme_bw() +
  geom_line(aes(group = year, color = year)) + ylim(c(0,0.85)) +
  xlab("Age") + ylab("Mortality Rate") + ggtitle("England and Wales")
jpn_plt <- ggplot(lifetab_5y_df[which(lifetab_5y_df$code == "JPN" & lifetab_5y_df$year > 1900),], 
       aes(x = age, y = mx)) + scale_color_continuous(name = "Year") + theme_bw() +
  geom_line(aes(group = year, color = year)) + ylim(c(0,0.85)) +
  xlab("Age") + ylab("") + ggtitle("Japan")
usa_plt <- ggplot(lifetab_5y_df[which(lifetab_5y_df$code == "USA" & lifetab_5y_df$year > 1900),], 
                  aes(x = age, y = mx)) + scale_color_continuous(name = "Year") + theme_bw() +
  geom_line(aes(group = year, color = year)) + ylim(c(0,0.85)) +
  xlab("Age") + ylab("") + ggtitle("USA")
nzl_plt <- ggplot(lifetab_5y_df[which(lifetab_5y_df$code == "NZL_NM" & lifetab_5y_df$year > 1900),], 
                  aes(x = age, y = mx)) + scale_color_continuous(name = "Year") + theme_bw() +
  geom_line(aes(group = year, color = year)) + ylim(c(0,0.85)) +
  xlab("Age") + ylab("") + ggtitle("New Zealand")
ggarrange(enw_plt,jpn_plt,usa_plt, nrow = 1, ncol=3, common.legend = FALSE)
ggsave("figures/data/mortality_rates_time.pdf", width = 12, height = 3)
rm(enw_plt, jpn_plt, usa_plt, nzl_plt)
     
"
Calculate best practice
"
lifetab_5y_df$best_practice <- 0
lifetab_5y_df$best_practice_alt <- 0
all_years <- sort(unique(lifetab_5y_df$years))
all_ages <- sort(unique(lifetab_5y_df$age))
for (yy in all_years ){
  print(yy)
  # Find which country had best LE at birth 
  year_df <- lifetab_5y_df[which(lifetab_5y_df$years == yy & 
                                   lifetab_5y_df$age == 0),]
  max_le <- max(year_df$ex, na.rm = T)
  bp_country <- year_df$code[which(year_df$ex == max_le)]
  # Assign best practice variable
  lifetab_5y_df$best_practice[which(lifetab_5y_df$years == yy & 
                        lifetab_5y_df$code == bp_country)] <- 1
  # Alternative, find the lowest mortality rate at each age for that year
  for (aa in all_ages){
    age_df <- lifetab_5y_df[which(lifetab_5y_df$years == yy & 
                                     lifetab_5y_df$age == aa),]
    min_mx <- min(age_df$mx, na.rm = T)
    bp_alt_country <- age_df$code[which(age_df$mx == min_mx)]
    if (length(bp_alt_country) > 1){
      print(paste(c("for",yy,"age",aa,"multiple countries:",bp_alt_country), collapse = " "))
      lifetab_5y_df$best_practice_alt[
        which(lifetab_5y_df$years == yy & lifetab_5y_df$age == aa &
                lifetab_5y_df$code == bp_alt_country[1])] <- 1
      lifetab_5y_df$best_practice_alt[
        which(lifetab_5y_df$years == yy & lifetab_5y_df$age == aa &
                lifetab_5y_df$code %in% bp_alt_country[2:length(bp_alt_country)])] <- 0.5
      
    } else{
      lifetab_5y_df$best_practice_alt[which(lifetab_5y_df$years == yy & 
                                          lifetab_5y_df$age == aa &
                                          lifetab_5y_df$code == bp_alt_country)] <- 1
    }
    
  }
  
}
rm(year_df, age_df, bp_country, bp_alt_country, min_mx, max_le, aa, yy, all_ages, all_years)

# Plot BP life expectancy over time, survival and mortality curves
bp_le_plt <- ggplot(lifetab_5y_df[which(lifetab_5y_df$age == 0 & lifetab_5y_df$year > 1900 &
                                          lifetab_5y_df$best_practice == 1),]) + theme_bw() +
  geom_smooth(aes(x = year, y = ex), method = "lm") +
  geom_point(aes(x = year, y = ex, color = code)) +
  scale_color_discrete(name = "Country") + ggtitle("Life Expectancy") +
  xlab("Year") + ylab("Life Expectancy at birth")
bp_s_plt <- ggplot(lifetab_5y_df[which(lifetab_5y_df$best_practice == 1 & 
                                   lifetab_5y_df$year > 1900),]) + theme_bw() +
  geom_line(aes(x = age, y = lx, group = year, color = year)) + ylim(c(0,1)) +
  scale_color_continuous(name = "Year") + xlab("Age") + ylab("Survival Rate") +
  ggtitle("Survival Rate")
bp_m_plt <- ggplot(lifetab_5y_df[which(lifetab_5y_df$best_practice == 1 & 
                                   lifetab_5y_df$year > 1900),]) + theme_bw() +
  geom_line(aes(x = age, y = mx, group = year, color = year)) + ylim(c(0,1)) +
  scale_color_continuous(name = "Year") + xlab("Age") + ylab("Mortality  Rate") +
  ggtitle("Mortality Rate")
ggarrange(bp_m_plt,bp_s_plt,bp_le_plt, nrow = 1, ncol=3, common.legend = FALSE)
ggsave("figures/data/best_practice_data.pdf", width = 15, height = 4)
rm(bp_le_plt, bp_s_plt, bp_m_plt)




# Plot BP life expectancy over time, survival and mortality curves
m2S <- function(mx){
  S <- rep(1, length(mx))
  for (aa in 2:length(mx)){
    S[aa] <- S[aa-1]*(1-mx[aa-1])
  }
  return(S)
}
S2LE <- function(S){
  LE <- rep(NA, length(S))
  for (aa in 1:length(mx)){
    S_aa <- S[aa:length(mx)]/S[aa]
    LE[aa] <- sum(S_aa)
  }
  return(LE) 
}

bp_df_alt <- lifetab_5y_df[which(lifetab_5y_df$best_practice_alt ==1 & 
                                   lifetab_5y_df$year > 1900),]
bp_df_alt <- bp_df_alt[with(bp_df_alt, order(year, age)),]
bp_df_alt[,c("lx_alt", "ex_alt")] <- NA
for (yy in unique(bp_df_alt$year)){
  mx_alt <- bp_df_alt$mx[which(bp_df_alt$year == yy)]
  lx_alt <- m2S(mx_alt)
  ex_alt <- S2LE(lx_alt)
  bp_df_alt$lx_alt[which(bp_df_alt$year == yy)] <- lx_alt
  bp_df_alt$ex_alt[which(bp_df_alt$year == yy)] <- ex_alt
}
rm(mx_alt, lx_alt, ex_alt)


bp_le_plt <- ggplot(bp_df_alt[which(lifetab_5y_df$age == 0),]) + theme_bw() +
  geom_smooth(aes(x = year, y = ex_alt), method = "lm") +
  geom_point(aes(x = year, y = ex_alt)) +
  scale_color_discrete(name = "Country") + ggtitle("Life Expectancy") +
  xlab("Year") + ylab("Life Expectancy at birth")
bp_s_plt <- ggplot(bp_df_alt) + theme_bw() +
  geom_line(aes(x = age, y = lx_alt, group = year, color = year)) + ylim(c(0,1)) +
  scale_color_continuous(name = "Year") + xlab("Age") + ylab("Survival Rate") +
  ggtitle("Survival Rate")
bp_m_plt <- ggplot(bp_df_alt) + theme_bw() +
  geom_line(aes(x = age, y = mx, group = year, color = year)) + ylim(c(0,1)) +
  scale_color_continuous(name = "Year") + xlab("Age") + ylab("Mortality  Rate") +
  ggtitle("Mortality Rate")
ggarrange(bp_m_plt,bp_s_plt,bp_le_plt, nrow = 1, ncol=3, common.legend = FALSE)
ggsave("figures/data/best_practice_data_alt.pdf", width = 15, height = 4)
rm(bp_le_plt, bp_s_plt, bp_m_plt)

## Zoom in on a couple years
# 2018 
bp_2018_plt <- ggplot(bp_df_alt[which(bp_df_alt$year == 2018),]) + theme_bw() +
  geom_point(aes(x = age, y = mx, color = code)) +
  xlab("Age") + ylab("Mortality Rate") +
  ggtitle("BP Mortality Rate 2018")
# 1983
bp_1983_plt <- ggplot(bp_df_alt[which(bp_df_alt$year == 1983),]) + theme_bw() +
  geom_point(aes(x = age, y = mx, color = code)) +
  xlab("Age") + ylab("Mortality Rate") +
  ggtitle("BP Mortality Rate 1983")
# 1943
bp_1943_plt <- ggplot(bp_df_alt[which(bp_df_alt$year == 1943),]) + theme_bw() +
  geom_point(aes(x = age, y = mx, color = code)) +
  xlab("Age") + ylab("Mortality Rate") +
  ggtitle("BP Mortality Rate 1943")
# 1903
bp_1903_plt <- ggplot(bp_df_alt[which(bp_df_alt$year == 1903),]) + theme_bw() +
  geom_point(aes(x = age, y = mx, color = code)) +
  xlab("Age") + ylab("Mortality Rate") +
  ggtitle("BP Mortality Rate 1903")
ggarrange(bp_2018_plt,bp_1983_plt,bp_1943_plt,bp_1903_plt, nrow = 2, ncol=2, common.legend = FALSE)
ggsave("figures/data/best_practice_examples_alt.pdf", width = 15, height = 8)
rm(bp_df_alt, bp_2018_plt, bp_1983_plt, bp_1943_plt, bp_1903_plt)




"
Add names as well as codes
"

all_codes <- unique(lifetab_5y_df$code)
convert_df <- data.frame(code = all_codes, name = NA)

convert_df$name[which(convert_df$code== "AUS")] <- "Australia"
convert_df$name[which(convert_df$code== "AUT")] <- "Austria"
convert_df$name[which(convert_df$code== "BEL")] <- "Belgium"
convert_df$name[which(convert_df$code== "BGR")] <- "Bulgaria"
convert_df$name[which(convert_df$code== "BLR")] <- "Belarus"
convert_df$name[which(convert_df$code== "CAN")] <- "Canada"
convert_df$name[which(convert_df$code== "CHE")] <- "Switzerland"
convert_df$name[which(convert_df$code== "CHL")] <- "Chile"
convert_df$name[which(convert_df$code== "CZE")] <- "Czechia"
convert_df$name[which(convert_df$code== "DEU")] <- "Germany"
convert_df$name[which(convert_df$code== "DEUTE")] <- "East Germany"
convert_df$name[which(convert_df$code== "DEUTW")] <- "West Germany"
convert_df$name[which(convert_df$code== "DNK")] <- "Denmark"
convert_df$name[which(convert_df$code== "ENW")] <- "England and Wales"
convert_df$name[which(convert_df$code== "ESP")] <- "Spain"
convert_df$name[which(convert_df$code== "EST")] <- "Estonia"
convert_df$name[which(convert_df$code== "FIN")] <- "Finland"
convert_df$name[which(convert_df$code== "FRA")] <- "France"
convert_df$name[which(convert_df$code== "GBR")] <- "United Kingdom"
convert_df$name[which(convert_df$code== "GRC")] <- "Greece"
convert_df$name[which(convert_df$code== "HKG")] <- "Hong Kong"
convert_df$name[which(convert_df$code== "HRV")] <- "Croatia"
convert_df$name[which(convert_df$code== "HUN")] <- "Hungary"
convert_df$name[which(convert_df$code== "IRL")] <- "Ireland"
convert_df$name[which(convert_df$code== "ISL")] <- "Iceland"
convert_df$name[which(convert_df$code== "ISR")] <- "Israel"
convert_df$name[which(convert_df$code== "ITA")] <- "Italy"
convert_df$name[which(convert_df$code== "JPN")] <- "Japan"
convert_df$name[which(convert_df$code== "KOR")] <- "South Korea"
convert_df$name[which(convert_df$code== "LTU")] <- "Lithuania"
convert_df$name[which(convert_df$code== "LUX")] <- "Luxembourg"
convert_df$name[which(convert_df$code== "LVA")] <- "Latvia"
convert_df$name[which(convert_df$code== "NLD")] <- "Netherlands"
convert_df$name[which(convert_df$code== "NOR")] <- "Norway"
#convert_df$name[which(convert_df$code== "NZL")] <- "New Zealand"
convert_df$name[which(convert_df$code== "NZL_NM")] <- "New Zealand (non-Maori)"
convert_df$name[which(convert_df$code== "POL")] <- "Poland"
convert_df$name[which(convert_df$code== "PRT")] <- "Portugal"
convert_df$name[which(convert_df$code== "RUS")] <- "Russia"
convert_df$name[which(convert_df$code== "SCO")] <- "Scotland"
convert_df$name[which(convert_df$code== "SVK")] <- "Slovakia"
convert_df$name[which(convert_df$code== "SVN")] <- "Slovenia"
convert_df$name[which(convert_df$code== "SWE")] <- "Sweden"
convert_df$name[which(convert_df$code== "TWN")] <- "Taiwan"
convert_df$name[which(convert_df$code== "UKR")] <- "Ukraine"
convert_df$name[which(convert_df$code== "USA")] <- "United States of America"


lifetab_5y_export <- merge(lifetab_5y_df, convert_df, by = "code", all.x = TRUE)

lifetab_5y_export <- lifetab_5y_export[,c("code", "name", "years", "year", "age",
                                          "mx", "lx", "ex", "Total", 
                                          "best_practice")]

write.csv(lifetab_5y_export, "data/clean/all_lifetab.csv", row.names = FALSE)





"
Back out a mortality rate from raw deaths and population data
"
deaths_df <- read.csv("data/raw/USA_deaths.csv", stringsAsFactors = FALSE)

pop_df <- as.data.frame(read.table("data/raw/USA_population.txt", header = TRUE))

# Pre-cleaningh
pop_df$Age[which(pop_df$Age == "110+")] <- "110"
pop_df$Age <- as.numeric(pop_df$Age)
pop_df$Year[which(pop_df$Year == "1959+")] <- "1959"
pop_df <- pop_df[which(pop_df$Year != "1959-"),] 
pop_df$Year <- as.numeric(pop_df$Year)
names(pop_df)[1:2] <- c("year", "age")



years <- c(1933, 1940, 1950, 1960, 1970, 1980, 1990, 2000, 2010, 2020)
ggplot(pop_df[which(pop_df$year %in% years),], aes(x = age, y = Total)) + 
  geom_line(aes(group = year, color = as.factor(year))) +
  xlab("Age") + ylab("Population") + theme_bw() +
  ggtitle("US Population over time")
ggsave("figures/USA/population.pdf", width = 6, height = 4)


write.csv(pop_df, "data/clean/USA_population.csv", row.names = FALSE)





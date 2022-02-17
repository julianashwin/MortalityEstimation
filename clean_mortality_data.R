setwd("/Users/julianashwin/Documents/GitHub/MortalityEstimation/")
rm(list=ls())

require(ggplot2)
require(HMDHFDplus)
require(stringr)
require(plyr)


"
Import data for all available countries
"
# Initialise a dataframe to fill with each country
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

table(lifetab_df$year)

# Plot for most recent period
plt_codes <- c("BEL", "DNK", "FRA", "ITA", "NOR", "CHE", "FIN", "ISL", "NLD", "SWE")
all2010_plt <- ggplot(lifetab_df[which(lifetab_df$year == "2010" & 
                                         lifetab_df$code %in% plt_codes),], aes(x = age, y = mx)) + 
  geom_line(aes(group = code, color = code)) + scale_color_discrete(name = "Country") + 
  xlab("Age") + ylab("") + theme_bw() + ylim(c(0,0.85)) +
  ggtitle("2010")
all1900_plt <- ggplot(lifetab_df[which(lifetab_df$year == "1900" & 
                                         lifetab_df$code %in% plt_codes),], aes(x = age, y = mx)) + 
  geom_line(aes(group = code, color = code)) + scale_color_discrete(name = "Country") + 
  xlab("Age") + ylab("Mortality Rate") + theme_bw() + ylim(c(0,0.85)) +
  ggtitle("1900")
ggarrange(all1900_plt,all2010_plt, nrow = 1, ncol=2, common.legend = TRUE)
ggsave("figures/all/all_1990_2010_mortality_rates.pdf", width = 6, height = 4)
rm(all1900_plt, all2010_plt)

"
Get five year average
"
# Collapse down to averages over 5 year intervals
lifetab_5y_df <- lifetab_df
lifetab_5y_df$years <- as.character(cut(lifetab_5y_df$year, seq(1750, 2020, 5)) )
lifetab_5y_df <- ddply(lifetab_5y_df, .(code, years, age), numcolwise(mean))

# Plot an example
ita_plt <- ggplot(lifetab_5y_df[which(lifetab_5y_df$code == "ITA" & lifetab_5y_df$year > 1900),], 
       aes(x = age, y = mx)) + scale_color_continuous(name = "Year") + theme_bw() +
  geom_line(aes(group = year, color = year)) + ylim(c(0,0.85)) +
  xlab("Age") + ylab("Mortality Rate") + ggtitle("Italy")
jpn_plt <- ggplot(lifetab_5y_df[which(lifetab_5y_df$code == "JPN" & lifetab_5y_df$year > 1900),], 
       aes(x = age, y = mx)) + scale_color_continuous(name = "Year") + theme_bw() +
  geom_line(aes(group = year, color = year)) + ylim(c(0,0.85)) +
  xlab("Age") + ylab("") + ggtitle("Japan")
usa_plt <- ggplot(lifetab_5y_df[which(lifetab_5y_df$code == "USA" & lifetab_5y_df$year > 1900),], 
                  aes(x = age, y = mx)) + scale_color_continuous(name = "Year") + theme_bw() +
  geom_line(aes(group = year, color = year)) + ylim(c(0,0.85)) +
  xlab("Age") + ylab("") + ggtitle("USA")
ggarrange(ita_plt,jpn_plt,usa_plt, nrow = 1, ncol=3, common.legend = FALSE)
ggsave("figures/all/mortality_rates_time.pdf", width = 12, height = 3)
rm(ita_plt, jpn_plt, usa_plt)
     
"
Calculate best practice
"
lifetab_5y_df$best_practice <- 0
all_years <- sort(unique(lifetab_5y_df$years))
for (yy in all_years ){
  # Find which country had best LE at birth 
  year_df <- lifetab_5y_df[which(lifetab_5y_df$years == yy & 
                                   lifetab_5y_df$age == 0),]
  max_le <- max(year_df$ex, na.rm = T)
  bp_country <- year_df$code[which(year_df$ex == max_le)]
  # Assign best practice variable
  lifetab_5y_df$best_practice[which(lifetab_5y_df$years == yy & 
                        lifetab_5y_df$code == bp_country)] <- 1
}
rm(year_df)

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
ggsave("figures/all/best_practice_data.pdf", width = 15, height = 4)
rm(bp_le_plt, bp_s_plt, bp_m_plt)





write.csv(lifetab_5y_df, "data/clean/all_lifetab.csv", row.names = FALSE)

 
"
Sweden
"
# Import life table from txct file (copied from HMD website)
lifetab_df <- as.data.frame(read.table("data/raw/SWE_lifetab.txt", header = TRUE))

# Some quick pre-cleaning
lifetab_df <- lifetab_df[which(lifetab_df$Age != "110+"),]
lifetab_df$Age <- as.numeric(lifetab_df$Age)
lifetab_df <- lifetab_df[which(lifetab_df$Age <= 100),]

ggplot(lifetab_df, aes(x = Age, y = mx)) + 
  geom_line(aes(group = Year, color = Year)) +
  xlab("Age") + ylab("Mortality Rate") + theme_bw() +
  ggtitle("Swedish Mortality Rates over time")
ggsave("figures/SWE/mortality_rates.pdf", width = 6, height = 4)

# Save the cleaned up data
write.csv(lifetab_df, "data/clean/SWE_life.csv", row.names = FALSE)



"
USA
"
# Import life table from txt file (copied from HMD website)
lifetab_df <- as.data.frame(read.table("data/raw/USA_lifetab.txt", header = TRUE))

# Some quick pre-cleaning
lifetab_df <- lifetab_df[which(lifetab_df$Age != "110+"),]
lifetab_df$Age <- as.numeric(lifetab_df$Age)
lifetab_df <- lifetab_df[which(lifetab_df$Age <= 100),]

ggplot(lifetab_df[lifetab_df$Year > 1900,], aes(x = Age, y = mx)) + 
  geom_line(aes(group = Year, color = Year)) +
  xlab("Age") + ylab("Mortality Rate") + theme_bw() +
  ggtitle("USA Mortality Rates over time")
ggsave("figures/USA/mortality_rates.pdf", width = 6, height = 4)


# Save the cleaned up data

names(lifetab_df)[1:2] <- c("year", "age")
write.csv(lifetab_df, "data/clean/USA_life.csv", row.names = FALSE)



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





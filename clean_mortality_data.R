setwd("/Users/julianashwin/Documents/GitHub/MortalityEstimation/")
rm(list=ls())

require(ggplot2)
require(HMDHFDplus)

"
Import data directly from HMB with readHMD
"
#readHMDweb("USA", "julianashwin@gmail.com", "pp@G..Jr7Rfa..i")

#US_lt <- read.hmd("julianashwin@gmail.com", "pp@G..Jr7Rfa..i",
#                  country = ("U.S.A"), sex = c("m", "f", "b" ), 
#                  HMDurl='http://www.mortality.org/hmd', dataType = "lt", 
#                  ltCol = c('m', 'q', 'a', 'l', 'd', 'L', 'T', 'e'))

      
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
# Import life table from txct file (copied from HMD website)
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
pop_df <- read.csv("data/raw/USA_population.csv", stringsAsFactors = FALSE)


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





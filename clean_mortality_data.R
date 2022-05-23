setwd("/Users/julianashwin/Documents/GitHub/MortalityEstimation/")
rm(list=ls())

require(ggplot2)
require(ggpubr)
require(HMDHFDplus)
require(stringr)
require(plyr)


"
Import data for all available countries
"
## Mortality, survival and life expectancy

# Function to import hmd lifetable data from a given directory
import_hmd_data <- function(import_dir){
  
  var_names <- c("code", "year", "age", "mx", "lx", "ex")
  # Initialise dataframe to fill with each country
  lifetab_df <- as.data.frame(matrix(NA, nrow = 0, ncol = length(var_names)), var_names)
  
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
  
  lifetab_df <- lifetab_df[which(lifetab_df$age != "110+"),]
  lifetab_df$age <- as.numeric(lifetab_df$age)
  lifetab_df <- lifetab_df[which(!is.na(lifetab_df$mx)),]
  
  return(lifetab_df)
}
lifetab_df_female <- import_hmd_data("data/raw/fltper_1x1/")
lifetab_df_male <- import_hmd_data("data/raw/mltper_1x1/")
lifetab_df_both <- import_hmd_data("data/raw/bltper_1x1/")


# Combine these mortality curves for the separate genders
names(lifetab_df_female)[4:6] <- paste0(c("mx","lx", "ex"),"_f")
names(lifetab_df_male)[4:6] <- paste0(c("mx","lx", "ex"),"_m")

lifetab_df <- merge(lifetab_df_both, lifetab_df_female, by = c("code", "year", "age"))
lifetab_df <- merge(lifetab_df, lifetab_df_male, by = c("code", "year", "age"))
rm(lifetab_df_both, lifetab_df_female, lifetab_df_male)

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



lifetab_df <- merge(lifetab_df, pop_df, by = c("code", "year", "age"), all.x = T)
lifetab_df <- lifetab_df[with(lifetab_df, order(code, year, age)), ]
rm(pop_df)

# Fix some codes up 
unique(lifetab_df$code)
lifetab_df$code[which(lifetab_df$code == "DEUTNP")] <- "DEU"
lifetab_df$code[which(lifetab_df$code == "FRATNP")] <- "FRA"
lifetab_df$code[which(lifetab_df$code == "GBR_NIR")] <- "NIR"
lifetab_df$code[which(lifetab_df$code == "GBR_NP")] <- "GBR"
lifetab_df$code[which(lifetab_df$code == "GBR_SCO")] <- "SCO"
lifetab_df$code[which(lifetab_df$code == "GBRTENW")] <- "ENW"
lifetab_df$code[which(lifetab_df$code == "NZL_NP")] <- "NZL"

extra_nzl_years <- unique(lifetab_df$year[which(lifetab_df$code == "NZL")])[
  which(!unique(lifetab_df$year[which(lifetab_df$code == "NZL")]) %in%
          unique(lifetab_df$year[which(lifetab_df$code == "NZL_NM")]))]
extra_obs <- lifetab_df[which(lifetab_df$code == "NZL" & lifetab_df$year %in% extra_nzl_years),]
extra_obs$code <- "NZL_NM"

lifetab_df <- rbind(lifetab_df, extra_obs)
rm(extra_obs, extra_nzl_years)

unique(lifetab_df$code)
lifetab_df <- lifetab_df[which(lifetab_df$code != "FRACNP"),]
lifetab_df <- lifetab_df[which(lifetab_df$code != "GBRCENW"),]


"
Some comparison plots
"
# Mortality rates
plt_codes <- c("BEL", "DNK", "ENW", "FRA", "ITA", "NOR", "NZL_NM", "CHE", "FIN", "ISL", "NLD", "SWE")
all2010_plt <- ggplot(lifetab_df[which(lifetab_df$year == "2010" & 
                                         lifetab_df$code %in% plt_codes),], aes(x = age, y = mx_f)) + 
  geom_line(aes(group = code, color = code)) + scale_color_discrete(name = "Country") + 
  xlab("Age") + ylab("") + theme_bw() + ylim(c(0,0.85)) +
  ggtitle("2010")
all1901_plt <- ggplot(lifetab_df[which(lifetab_df$year == "1901" & 
                                         lifetab_df$code %in% plt_codes),], aes(x = age, y = mx_f)) + 
  geom_line(aes(group = code, color = code)) + scale_color_discrete(name = "Country") + 
  xlab("Age") + ylab("Mortality Rate") + theme_bw() + ylim(c(0,0.85)) +
  ggtitle("1901")
ggarrange(all1901_plt,all2010_plt, nrow = 1, ncol=2, common.legend = TRUE)
ggsave("figures/data/all_1990_2010_mortality_rates.pdf", width = 6, height = 4)

# Life expectancy
plt_codes <- c("BEL", "DNK", "ENW", "FRA", "ITA", "NOR", "NZL_NM", "CHE", "FIN", "ISL", "NLD", "SWE")
all2010_plt <- ggplot(lifetab_df[which(lifetab_df$year == "2010" & 
                                         lifetab_df$code %in% plt_codes),], aes(x = age, y = ex_f)) + 
  geom_line(aes(group = code, color = code)) + scale_color_discrete(name = "Country") + 
  xlab("Age") + ylab("") + theme_bw() + ylim(c(0,85)) +
  ggtitle("2010")
all1901_plt <- ggplot(lifetab_df[which(lifetab_df$year == "1901" & 
                                         lifetab_df$code %in% plt_codes),], aes(x = age, y = ex_f)) + 
  geom_line(aes(group = code, color = code)) + scale_color_discrete(name = "Country") + 
  xlab("Age") + ylab("Remaining life expectancy") + theme_bw() + ylim(c(0,85)) +
  ggtitle("1901")
ggarrange(all1901_plt,all2010_plt, nrow = 1, ncol=2, common.legend = TRUE)
ggsave("figures/data/all_1990_2010_life_expectancies.pdf", width = 6, height = 4)

# Population
all2010_plt <- ggplot(lifetab_df[which(lifetab_df$year == "2010" & 
                                         lifetab_df$code %in% plt_codes),], aes(x = age, y = Female)) + 
  geom_line(aes(group = code, color = code)) + scale_color_discrete(name = "Country") + 
  xlab("Age") + ylab("") + theme_bw() + 
  ggtitle("2010")
all1901_plt <- ggplot(lifetab_df[which(lifetab_df$year == "1901" & 
                                         lifetab_df$code %in% plt_codes),], aes(x = age, y = Female)) + 
  geom_line(aes(group = code, color = code)) + scale_color_discrete(name = "Country") + 
  xlab("Age") + ylab("Population") + theme_bw() + 
  ggtitle("1901")
ggarrange(all1901_plt,all2010_plt, nrow = 1, ncol=2, common.legend = TRUE)
ggsave("figures/data/all_1990_2010_population.pdf", width = 6, height = 4)

rm(all1901_plt, all2010_plt)





"
Who's looking good in the 21st century?
"
ex_df <- lifetab_df[which(lifetab_df$year >= 1900 & lifetab_df$age == 0 & 
                                  lifetab_df$ex_f >= 60),]
ggplot(ex_df, aes(x = year, y = ex_f)) + 
  geom_line(aes(group = code, color = code)) + scale_color_discrete(name = "Country") + 
  xlab("Year") + ylab("LE(0)") + theme_bw() +
  ggtitle("2010")

current_ex_df <- lifetab_df[which(lifetab_df$year == 2018 & lifetab_df$age == 0 & 
                                    lifetab_df$ex_f >= 84),]

rm(ex_df,current_ex_df)

"
Get five year average
"
# Collapse down to averages over 5 year intervals
lifetab_5y_df <- lifetab_df[which(lifetab_df$year <= 2020),]
lifetab_5y_df$years <- as.character(cut(lifetab_5y_df$year, seq(1750, 2020, 5)) )
lifetab_5y_df <- ddply(lifetab_5y_df, .(code, years, age), numcolwise(mean))

# Make sure year represents the mid-point of the window
lifetab_5y_df[,c("year_lower", "year_upper")] <- 
  do.call(rbind, str_split(str_remove(str_remove(lifetab_5y_df$years, "\\("), "\\]"), ","))
lifetab_5y_df$year <- as.numeric(lifetab_5y_df$year_lower) + 3

## Plot some examples
# ENW
enw_plt_m <- ggplot(lifetab_5y_df[which(lifetab_5y_df$code == "ENW" & lifetab_5y_df$year > 1900),], 
                    aes(x = age, y = mx_m)) + scale_color_continuous(name = "Year") + theme_bw() +
  geom_line(aes(group = year, color = year)) + ylim(c(0,0.85)) +
  xlab("Age") + ylab("Mortality Rate") + ggtitle("England and Wales (Male)")
enw_plt_f <- ggplot(lifetab_5y_df[which(lifetab_5y_df$code == "ENW" & lifetab_5y_df$year > 1900),], 
                    aes(x = age, y = mx_f)) + scale_color_continuous(name = "Year") + theme_bw() +
  geom_line(aes(group = year, color = year)) + ylim(c(0,0.85)) +
  xlab("Age") + ylab("Mortality Rate") + ggtitle("England and Wales (Female)")
# JAP
jpn_plt_m <- ggplot(lifetab_5y_df[which(lifetab_5y_df$code == "JPN" & lifetab_5y_df$year > 1900),], 
                    aes(x = age, y = mx_m)) + scale_color_continuous(name = "Year") + theme_bw() +
  geom_line(aes(group = year, color = year)) + ylim(c(0,0.85)) +
  xlab("Age") + ylab("") + ggtitle("Japan (Male)")
jpn_plt_f <- ggplot(lifetab_5y_df[which(lifetab_5y_df$code == "JPN" & lifetab_5y_df$year > 1900),], 
                    aes(x = age, y = mx_f)) + scale_color_continuous(name = "Year") + theme_bw() +
  geom_line(aes(group = year, color = year)) + ylim(c(0,0.85)) +
  xlab("Age") + ylab("") + ggtitle("Japan (Female)")
# USA
usa_plt_m <- ggplot(lifetab_5y_df[which(lifetab_5y_df$code == "USA" & lifetab_5y_df$year > 1900),], 
                    aes(x = age, y = mx_m)) + scale_color_continuous(name = "Year") + theme_bw() +
  geom_line(aes(group = year, color = year)) + ylim(c(0,0.85)) +
  xlab("Age") + ylab("") + ggtitle("USA (Male)")
usa_plt_f <- ggplot(lifetab_5y_df[which(lifetab_5y_df$code == "USA" & lifetab_5y_df$year > 1900),], 
                    aes(x = age, y = mx_f)) + scale_color_continuous(name = "Year") + theme_bw() +
  geom_line(aes(group = year, color = year)) + ylim(c(0,0.85)) +
  xlab("Age") + ylab("") + ggtitle("USA (Female)")
# NZL_NM
nzl_plt_m <- ggplot(lifetab_5y_df[which(lifetab_5y_df$code == "NZL_NM" & lifetab_5y_df$year > 1900),], 
                    aes(x = age, y = mx_m)) + scale_color_continuous(name = "Year") + theme_bw() +
  geom_line(aes(group = year, color = year)) + ylim(c(0,0.85)) +
  xlab("Age") + ylab("") + ggtitle("New Zealand (Male)")
nzl_plt_f <- ggplot(lifetab_5y_df[which(lifetab_5y_df$code == "NZL_NM" & lifetab_5y_df$year > 1900),], 
                    aes(x = age, y = mx_f)) + scale_color_continuous(name = "Year") + theme_bw() +
  geom_line(aes(group = year, color = year)) + ylim(c(0,0.85)) +
  xlab("Age") + ylab("") + ggtitle("New Zealand (Female)")

ggarrange(enw_plt_f,jpn_plt_f,usa_plt_f, nrow = 1, ncol=3, common.legend = FALSE)
ggsave("figures/data/mortality_rates_time.pdf", width = 12, height = 3)
rm(enw_plt_m, enw_plt_f, jpn_plt_m, jpn_plt_f, usa_plt_m, usa_plt_f, nzl_plt_m, nzl_plt_f)





"
Compute empirical lifespan inequality from the mortality data
" 
compute_ineq <- function(lifetab_5y_df){
  lifetab_5y_df$Hx <- NA
  lifetab_5y_df$Hx_f <- NA
  lifetab_5y_df$Hx_m <- NA
  for (code in unique(lifetab_5y_df$code)){
    print(code)
    years <- unique(lifetab_5y_df$year[which(lifetab_5y_df$code == code)])
    for (year in years){
      data_obs <- which(lifetab_5y_df$code == code & lifetab_5y_df$year == year)
      ages <- lifetab_5y_df$age[data_obs]
      LEs <- lifetab_5y_df$ex[data_obs]
      Ss <- lifetab_5y_df$lx[data_obs]
      LE_fs <- lifetab_5y_df$ex_f[data_obs]
      S_fs <- lifetab_5y_df$lx_f[data_obs]
      LE_ms <- lifetab_5y_df$ex_m[data_obs]
      S_ms <- lifetab_5y_df$lx_m[data_obs]
      
      Hs <- rep(NA, length(ages))
      H_fs <- rep(NA, length(ages))
      H_ms <- rep(NA, length(ages))
      for (ii in 1:length(ages)){
        # Extract LEs
        LE <- LEs[ii]
        LE_f <- LE_fs[ii]
        LE_m <- LE_ms[ii]
        # Extract Ss
        Sa_s <- Ss[ii:length(ages)]/Ss[ii]
        Sa_s <- Sa_s[Sa_s>0]
        Sa_fs <- S_fs[ii:length(ages)]/S_fs[ii]
        Sa_fs <- Sa_fs[Sa_fs>0]
        Sa_ms <- S_ms[ii:length(ages)]/S_ms[ii]
        Sa_ms <- Sa_ms[Sa_ms>0]
        # Compute Hs
        if (length(Sa_s) > 0){
          Hs[ii] <- -sum(Sa_s*log(Sa_s))/LE
        } else{
          Hs[ii] <- 0
        }
        if (length(Sa_fs) > 0){
          H_fs[ii] <- -sum(Sa_fs*log(Sa_fs))/LE_f
        } else{
          H_fs[ii] <- 0
        }
        if (length(Sa_ms) > 0){
          H_ms[ii] <- -sum(Sa_ms*log(Sa_ms))/LE_m
        } else{
          H_ms[ii] <- 0
        }
      }
      lifetab_5y_df$Hx[data_obs] <- Hs
      lifetab_5y_df$Hx_f[data_obs] <- H_fs
      lifetab_5y_df$Hx_m[data_obs] <- H_ms
    }
  }
  return(lifetab_5y_df)
}
lifetab_5y_df <- compute_ineq(lifetab_5y_df)
lifetab_df <- compute_ineq(lifetab_df)


"
Calculate best practice
"
identify_best_practice <- function(lifetab_df){
  # Initialise indicators
  lifetab_df$best_practice <- 0
  lifetab_df$best_practice_gender <- NA
  lifetab_df$best_practice_alt <- 0
  # Identify the set of years and ages 
  all_years <- sort(unique(lifetab_df$year))
  all_ages <- sort(unique(lifetab_df$age))
  for (yy in all_years ){
    print(yy)
    # Find which country had best LE at birth 
    year_df <- lifetab_df[which(lifetab_df$year == yy & 
                                     lifetab_df$age == 0),]
    max_le_f <- max(year_df$ex_f, na.rm = T)
    max_le_m <- max(year_df$ex_m, na.rm = T)
    if (max_le_f >= max_le_m){
      bp_country <- year_df$code[which(year_df$ex_f == max_le_f)]
      bp_gender <- "F"
    } else {
      bp_country <- year_df$code[which(year_df$ex_m == max_le_m)]
      bp_gender <- "M"
    }
    
    # Assign best practice variable
    lifetab_df$best_practice[which(lifetab_df$year == yy & 
                                        lifetab_df$code == bp_country)] <- 1
    lifetab_df$best_practice_gender[which(lifetab_df$year == yy & 
                                     lifetab_df$code == bp_country)] <- bp_gender
    # Alternative, find the lowest mortality rate at each age for that year
    for (aa in all_ages){
      age_df <- lifetab_df[which(lifetab_df$year == yy & 
                                      lifetab_df$age == aa),]
      min_mx <- min(age_df$mx_f, na.rm = T)
      bp_alt_country <- age_df$code[which(age_df$mx_f == min_mx)]
      if (length(bp_alt_country) > 1){
        print(paste(c("for",yy,"age",aa,"multiple countries:",bp_alt_country), collapse = " "))
        lifetab_df$best_practice_alt[
          which(lifetab_df$year == yy & lifetab_df$age == aa &
                  lifetab_df$code == bp_alt_country[1])] <- 1
        lifetab_df$best_practice_alt[
          which(lifetab_df$year == yy & lifetab_df$age == aa &
                  lifetab_df$code %in% bp_alt_country[2:length(bp_alt_country)])] <- 0.5
        
      } else{
        lifetab_df$best_practice_alt[which(lifetab_df$year == yy & 
                                             lifetab_df$age == aa &
                                             lifetab_df$code == bp_alt_country)] <- 1
      }
    }
  }
  return(lifetab_df)
}

lifetab_df <- identify_best_practice(lifetab_df)
lifetab_5y_df <- identify_best_practice(lifetab_5y_df)




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
convert_df$name[which(convert_df$code== "NZL_NM")] <- "New Zealand"# (non-Maori)"
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
lifetab_export <- merge(lifetab_df, convert_df, by = "code", all.x = TRUE)


" 
Plot BP life expectancy over time, survival and mortality curves
"
# Five year intervals
bp_le_plt <- ggplot(lifetab_5y_export[which(lifetab_5y_export$age == 0 & lifetab_5y_export$year > 1900 &
                                              lifetab_5y_export$best_practice == 1),]) + theme_bw() +
  geom_smooth(aes(x = year, y = ex_f), method = "lm") +
  geom_point(aes(x = year, y = ex_f, color = name)) +
  scale_color_discrete(name = "Country") + ggtitle("Life Expectancy") +
  xlab("Year") + ylab("Life Expectancy at birth")
bp_h_plt <- ggplot(lifetab_5y_export[which(lifetab_5y_export$age == 0 & lifetab_5y_export$year > 1900 &
                                 lifetab_5y_export$best_practice == 1),]) + theme_bw() +
  geom_smooth(aes(x = year, y = -log(Hx_f)), method = "loess") +
  geom_point(aes(x = year, y = -log(Hx_f), color = name)) +
  scale_color_discrete(name = "Country") + ggtitle("Lifespan Equality") +
  xlab("Year") + ylab("Lifespan Equality at birth")
bp_H_plt <- ggplot(lifetab_5y_export[which(lifetab_5y_export$age == 0 & lifetab_5y_export$year > 1900 &
                                             lifetab_5y_export$best_practice == 1),]) + theme_bw() +
  geom_smooth(aes(x = year, y = Hx_f), method = "loess") +
  geom_point(aes(x = year, y = Hx_f, color = name)) +
  scale_color_discrete(name = "Country") + ggtitle("Lifespan Inequality") +
  xlab("Year") + ylab("Lifespan Inequality at birth")
bp_s_plt <- ggplot(lifetab_5y_export[which(lifetab_5y_export$best_practice == 1 & 
                                             lifetab_5y_export$year > 1900),]) + theme_bw() +
  geom_line(aes(x = age, y = lx_f, group = year, color = year)) + ylim(c(0,1)) +
  scale_color_gradientn(colours = rainbow(5), name = "Year") + 
  xlab("Age") + ylab("Survival Rate") + ggtitle("Survival Rate")
bp_m_plt <- ggplot(lifetab_5y_export[which(lifetab_5y_export$best_practice == 1 & 
                                             lifetab_5y_export$year > 1900),]) + theme_bw() +
  geom_line(aes(x = age, y = mx_f, group = year, color = year)) + ylim(c(0,1)) +
  scale_color_gradientn(colours = rainbow(5), name = "Year") + 
  xlab("Age") + ylab("Mortality  Rate") + ggtitle("Mortality Rate")
ggarrange(bp_s_plt,bp_le_plt, bp_h_plt, nrow = 1, ncol=3, common.legend = FALSE)
ggsave("figures/data/best_practice_5y_data.pdf", width = 15, height = 4)

ggarrange(bp_H_plt, bp_h_plt, nrow = 1, ncol=2, common.legend = TRUE, legend = "right")
ggsave("figures/data/best_practice_5y_Hvh.pdf", width = 10, height = 4)


# One year intervals
bp_le_plt <- ggplot(lifetab_export[which(lifetab_export$age == 0 & lifetab_export$year > 1900 &
                                       lifetab_export$best_practice == 1),]) + theme_bw() +
  geom_smooth(aes(x = year, y = ex_f), method = "lm") +
  geom_point(aes(x = year, y = ex_f, color = name)) +
  scale_color_discrete(name = "Country") + ggtitle("Life Expectancy") +
  xlab("Year") + ylab("Life Expectancy at birth")
bp_h_plt <- ggplot(lifetab_export[which(lifetab_export$age == 0 & lifetab_export$year > 1900 &
                                      lifetab_export$best_practice == 1),]) + theme_bw() +
  geom_smooth(aes(x = year, y = -log(Hx_f)), method = "loess") +
  geom_point(aes(x = year, y = -log(Hx_f), color = name)) +
  scale_color_discrete(name = "Country") + ggtitle("Lifespan Equality") +
  xlab("Year") + ylab("Lifespan Equality at birth")
bp_s_plt <- ggplot(lifetab_export[which(lifetab_export$best_practice == 1 & 
                                      lifetab_export$year > 1900),]) + theme_bw() +
  geom_line(aes(x = age, y = lx_f, group = year, color = year)) + ylim(c(0,1)) +
  scale_color_gradientn(colours = rainbow(5), name = "Year") + 
  xlab("Age") + ylab("Survival Rate") + ggtitle("Survival Rate")
bp_m_plt <- ggplot(lifetab_export[which(lifetab_df$lifetab_export == 1 & 
                                      lifetab_export$year > 1900),]) + theme_bw() +
  geom_line(aes(x = age, y = mx_f, group = year, color = year)) + ylim(c(0,1)) +
  scale_color_gradientn(colours = rainbow(5), name = "Year") + 
  xlab("Age") + ylab("Mortality  Rate") + ggtitle("Mortality Rate")
ggarrange(bp_s_plt,bp_le_plt, bp_h_plt, nrow = 1, ncol=3, common.legend = FALSE)
ggsave("figures/data/best_practice_data.pdf", width = 15, height = 4)

rm(bp_le_plt, bp_s_plt, bp_m_plt)



if (FALSE){
  
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
    for (aa in 1:length(S)){
      S_aa <- S[aa:length(S)]/S[aa]
      LE[aa] <- sum(S_aa)
    }
    return(LE) 
  }
  
  create_bp_df <- function(lifetab_df){
    bp_df_alt <- lifetab_df[which(lifetab_df$best_practice_alt ==1 & 
                                       lifetab_df$year > 1900),]
    bp_df_alt <- bp_df_alt[with(bp_df_alt, order(year, age)),]
    bp_df_alt[,c("lx_alt", "ex_alt")] <- NA
    for (yy in unique(bp_df_alt$year)){
      mx_alt <- bp_df_alt$mx_f[which(bp_df_alt$year == yy)]
      lx_alt <- m2S(mx_alt)
      ex_alt <- S2LE(lx_alt)
      bp_df_alt$lx_alt[which(bp_df_alt$year == yy)] <- lx_alt
      bp_df_alt$ex_alt[which(bp_df_alt$year == yy)] <- ex_alt
    }
    return(bp_df_alt)
  }
  bp_alt_df <- create_bp_df(lifetab_df)
  bp_alt_5y_df <- create_bp_df(lifetab_5y_df)
  
  
  # 5 year intervals
  bp_le_plt <- ggplot(bp_alt_5y_df[which(bp_alt_5y_df$age == 0),]) + theme_bw() +
    geom_smooth(aes(x = year, y = ex_alt), method = "lm") +
    geom_point(aes(x = year, y = ex_alt)) +
    scale_color_discrete(name = "Country") + ggtitle("Life Expectancy") +
    xlab("Year") + ylab("Life Expectancy at birth")
  bp_s_plt <- ggplot(bp_alt_5y_df) + theme_bw() +
    geom_line(aes(x = age, y = lx_alt, group = year, color = year)) + ylim(c(0,1)) +
    scale_color_continuous(name = "Year") + xlab("Age") + ylab("Survival Rate") +
    ggtitle("Survival Rate")
  bp_m_plt <- ggplot(bp_alt_5y_df) + theme_bw() +
    geom_line(aes(x = age, y = mx, group = year, color = year)) + ylim(c(0,1)) +
    scale_color_continuous(name = "Year") + xlab("Age") + ylab("Mortality  Rate") +
    ggtitle("Mortality Rate")
  ggarrange(bp_m_plt,bp_s_plt,bp_le_plt, nrow = 1, ncol=3, common.legend = FALSE)
  ggsave("figures/data/best_practice_5y_data_alt.pdf", width = 15, height = 4)
  rm(bp_le_plt, bp_s_plt, bp_m_plt)
  # One year intervals
  bp_le_plt <- ggplot(bp_alt_df[which(bp_alt_df$age == 0),]) + theme_bw() +
    geom_smooth(aes(x = year, y = ex_alt), method = "lm") +
    geom_point(aes(x = year, y = ex_alt)) +
    scale_color_discrete(name = "Country") + ggtitle("Life Expectancy") +
    xlab("Year") + ylab("Life Expectancy at birth")
  bp_s_plt <- ggplot(bp_alt_df) + theme_bw() +
    geom_line(aes(x = age, y = lx_alt, group = year, color = year)) + ylim(c(0,1)) +
    scale_color_continuous(name = "Year") + xlab("Age") + ylab("Survival Rate") +
    ggtitle("Survival Rate")
  bp_m_plt <- ggplot(bp_alt_df) + theme_bw() +
    geom_line(aes(x = age, y = mx, group = year, color = year)) + ylim(c(0,1)) +
    scale_color_continuous(name = "Year") + xlab("Age") + ylab("Mortality  Rate") +
    ggtitle("Mortality Rate")
  ggarrange(bp_m_plt,bp_s_plt,bp_le_plt, nrow = 1, ncol=3, common.legend = FALSE)
  ggsave("figures/data/best_practice_data_alt.pdf", width = 15, height = 4)
  rm(bp_le_plt, bp_s_plt, bp_m_plt)
  
  ## Zoom in on a couple years
  # 2018 
  bp_2018_plt <- ggplot(bp_alt_df[which(bp_alt_df$year == 2018),]) + theme_bw() +
    geom_point(aes(x = age, y = mx, color = code)) +
    xlab("Age") + ylab("Mortality Rate") +
    ggtitle("BP Mortality Rate 2018")
  # 1983
  bp_1983_plt <- ggplot(bp_alt_df[which(bp_alt_df$year == 1983),]) + theme_bw() +
    geom_point(aes(x = age, y = mx, color = code)) +
    xlab("Age") + ylab("Mortality Rate") +
    ggtitle("BP Mortality Rate 1983")
  # 1943
  bp_1943_plt <- ggplot(bp_alt_df[which(bp_alt_df$year == 1943),]) + theme_bw() +
    geom_point(aes(x = age, y = mx, color = code)) +
    xlab("Age") + ylab("Mortality Rate") +
    ggtitle("BP Mortality Rate 1943")
  # 1903
  bp_1903_plt <- ggplot(bp_alt_df[which(bp_alt_df$year == 1903),]) + theme_bw() +
    geom_point(aes(x = age, y = mx, color = code)) +
    xlab("Age") + ylab("Mortality Rate") +
    ggtitle("BP Mortality Rate 1903")
  ggarrange(bp_2018_plt,bp_1983_plt,bp_1943_plt,bp_1903_plt, nrow = 2, ncol=2, common.legend = FALSE)
  ggsave("figures/data/best_practice_examples_alt.pdf", width = 15, height = 8)
  rm(bp_df_alt, bp_2018_plt, bp_1983_plt, bp_1943_plt, bp_1903_plt)
  
}  
  



lifetab_export <- lifetab_export[,c("code", "name", "year", "age",
                                    "mx", "mx_f", "mx_m", "lx", "lx_f", "lx_m", 
                                    "ex", "ex_f", "ex_m", "Hx", "Hx_f", "Hx_m",
                                    "Total", "Male", "Female", 
                                    "best_practice", "best_practice_alt")]
lifetab_5y_export <- lifetab_5y_export[,c("code", "name", "years", "year", "age",
                                          "mx", "mx_f", "mx_m", "lx", "lx_f", "lx_m", 
                                          "ex", "ex_f", "ex_m", "Hx", "Hx_f", "Hx_m",
                                          "Total", "Male", "Female",
                                          "best_practice", "best_practice_alt")]

write.csv(lifetab_5y_export, "data/clean/all_lifetab_5y.csv", row.names = FALSE)
write.csv(lifetab_export, "data/clean/all_lifetab_1y.csv", row.names = FALSE)
#lifetab_5y_export <- read.csv("data/clean/all_lifetab_5y.csv", stringsAsFactors = FALSE)

bp_df <-  lifetab_export[which(lifetab_export$best_practice == 1),]
write.csv(bp_df, "data/clean/bp.csv", row.names = FALSE)

bp_alt_df <-  lifetab_export[which(lifetab_export$best_practice_alt == 1),]
write.csv(bp_alt_df, "data/clean/bp_alt.csv", row.names = FALSE)

bp_5y_df <-  lifetab_5y_export[which(lifetab_5y_export$best_practice == 1),]
write.csv(bp_5y_df, "data/clean/bp_5y.csv", row.names = FALSE)
#bp_5y_df <- read.csv("data/clean/bp_5y.csv", stringsAsFactors = FALSE)

bp_alt_5y_df <-  lifetab_5y_export[which(lifetab_5y_export$best_practice_alt == 1),]
write.csv(bp_alt_5y_df, "data/clean/bp_alt_5y.csv", row.names = FALSE)

"
Look at life expectancy at older ages
"
plt_df <- bp_df[which(bp_df$age %in% c(0,10,20,30,40,50,60,70,80,90)),]
plt_df <- plt_df[order(plt_df$year, plt_df$age),]
bp_rle_plt <- ggplot(plt_df, aes(x = year)) + theme_bw() + 
  geom_line(aes(y = ex_f, color = as.factor(age))) + 
  scale_color_discrete(name = "Age") + 
  xlab("Year") + ylab("Remaining LE") + ggtitle("Best practice")
  

plt_df <- bp_alt_df[which(bp_alt_df$age %in% c(0,10,20,30,40,50,60,70,80,90)),]
plt_df <- plt_df[order(plt_df$year, plt_df$age),]
bp_alt_rle_plt <- ggplot(plt_df, aes(x = year)) + theme_bw() + 
  geom_line(aes(y = ex_f, color = as.factor(age))) + 
  scale_color_discrete(name = "Age") + 
  xlab("Year") + ylab("Remaining LE") + ggtitle("Best practice (alternative)")




plt_df <- bp_5y_df[which(bp_5y_df$age %in% c(0,10,20,30,40,50,60,70,80,90)),]
plt_df <- plt_df[order(plt_df$year, plt_df$age),]
bp_rle_plt <- ggplot(plt_df, aes(x = year)) + theme_bw() + 
  geom_line(aes(y = ex_f, color = as.factor(age))) + 
  scale_color_discrete(name = "Age") + 
  xlab("Year") + ylab("Remaining LE") + ggtitle("Best practice")


plt_df <- bp_alt_5y_df[which(bp_alt_5y_df$age %in% c(0,10,20,30,40,50,60,70,80,90)),]
plt_df <- plt_df[order(plt_df$year, plt_df$age),]
bp_alt_rle_plt <- ggplot(plt_df, aes(x = year)) + theme_bw() + 
  geom_line(aes(y = ex_f, color = as.factor(age))) + 
  scale_color_discrete(name = "Age") + 
  xlab("Year") + ylab("Remaining LE") + ggtitle("Best practice (alternative)")


ggarrange(bp_rle_plt, bp_alt_rle_plt, nrow = 1, ncol=2, common.legend = FALSE)
ggsave("figures/data/best_practice_rle.pdf", width = 10, height = 5)
rm(bp_rle_plt, bp_alt_rle_plt)



# Mortality 
temp_df <- bp_5y_df[which(bp_5y_df$year > 1900),]
temp_df <- temp_df[order(temp_df$year, temp_df$age),which(!str_detect(names(temp_df), "_m"))]
temp_df <- temp_df[,which(!str_detect(names(temp_df), "Ma|Fe|To|best"))]
temp_df <- temp_df[,which(!(names(temp_df) %in% c("ex", "lx", "mx", "Hx")))]

ggplot(temp_df) + theme_bw() +
  geom_line(aes(x = age, y = -log(Hx_f), group = year, color = year)) + 
  scale_color_gradientn(colours = rainbow(5), name = "Year")

bp_m_plt <- ggplot(bp_5y_df[which(bp_5y_df$year > 1900),]) + theme_bw() +
  geom_line(aes(x = age, y = mx_f, group = year, color = year)) + 
  scale_color_gradientn(colours = rainbow(5), name = "Year") + ylim(c(0,1)) + 
  xlab("Age") + ylab("Mortality rate") + ggtitle("Mortality")
# Survival
bp_s_plt <- ggplot(bp_5y_df[which(bp_5y_df$year > 1900),]) + theme_bw() +
  geom_line(aes(x = age, y = lx_f, group = year, color = year)) + 
  scale_color_gradientn(colours = rainbow(5), name = "Year") + ylim(c(0,1)) + 
  xlab("Age") + ylab("Survival rate") + ggtitle("Survival")
# Life expectancy
bp_le_plt <- ggplot(bp_5y_df[which(bp_5y_df$year > 1900),]) + theme_bw() +
  geom_line(aes(x = age, y = ex_f, group = year, color = year)) + 
  scale_color_gradientn(colours = rainbow(5), name = "Year") + 
  xlab("Age") + ylab("Remaining life expectancy") + ggtitle("Life expectancy")
# Lifespan inequality
bp_h_plt <- ggplot(bp_5y_df[which(bp_5y_df$year > 1900),]) + theme_bw() +
  geom_line(aes(x = age, y = -log(Hx_f), group = year, color = year)) + 
  scale_color_gradientn(colours = rainbow(5), name = "Year") + 
  xlab("Age") + ylab("Remaining lifespan equality") + ggtitle("Lifespan equality")

ggarrange(bp_m_plt, bp_s_plt,bp_le_plt, bp_h_plt, nrow = 1, ncol=4, common.legend = TRUE,
          legend = "right")
ggsave("figures/data/best_practice_5y_data.pdf", width = 16, height = 4)

plt1 <- ggarrange(bp_m_plt, bp_s_plt, nrow = 1, ncol=2, common.legend = TRUE,
          legend = "left")



bp_le_plt <- ggplot(bp_5y_df[which(bp_5y_df$age == 0 & bp_5y_df$year > 1900),]) + theme_bw() +
  geom_smooth(aes(x = year, y = ex_f), method = "lm") +
  geom_point(aes(x = year, y = ex_f, color = name)) +
  scale_color_discrete(name = "Country") + ggtitle("Life Expectancy") +
  xlab("Year") + ylab("Life Expectancy at birth")
bp_h_plt <- ggplot(bp_5y_df[which(bp_5y_df$age == 0 & bp_5y_df$year > 1900),]) + theme_bw() +
  geom_smooth(aes(x = year, y = -log(Hx_f)), method = "loess") +
  geom_point(aes(x = year, y = -log(Hx_f), color = name)) +
  scale_color_discrete(name = "Country") + ggtitle("Lifespan Equality") +
  xlab("Year") + ylab("Lifespan Equality at birth")
plt2 <- ggarrange(bp_le_plt, bp_h_plt, nrow = 1, ncol=2, common.legend = TRUE,
          legend = "right")

ggarrange(plt1, plt2, nrow = 1, ncol=2, common.legend = FALSE)
ggsave("figures/data/best_practice_5y_summary.pdf", width = 16, height = 4)





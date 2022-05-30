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
ggplot(lifetab_5y_df[which(lifetab_5y_df$code == "ENW" & lifetab_5y_df$year > 1900),], 
                    aes(x = age, y = mx_m)) + scale_color_continuous(name = "Year") + theme_bw() +
  geom_line(aes(group = year, color = year)) + ylim(c(0,0.85)) +
  xlab("Age") + ylab("Mortality Rate") + ggtitle("England and Wales (Male)")
ggplot(lifetab_5y_df[which(lifetab_5y_df$code == "ENW" & lifetab_5y_df$year > 1900),], 
                    aes(x = age, y = mx_f)) + scale_color_continuous(name = "Year") + theme_bw() +
  geom_line(aes(group = year, color = year)) + ylim(c(0,0.85)) +
  xlab("Age") + ylab("Mortality Rate") + ggtitle("England and Wales (Female)")





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
Compute 99.9% lifespan
" 
compute_lifespan <- function(lifetab_5y_df){
  lifetab_5y_df$l_99p9 <- NA
  lifetab_5y_df$l_99p9_f <- NA
  lifetab_5y_df$l_99p9_m <- NA
  for (code in unique(lifetab_5y_df$code)){
    print(code)
    years <- unique(lifetab_5y_df$year[which(lifetab_5y_df$code == code)])
    for (year in years){
      data_obs <- which(lifetab_5y_df$code == code & lifetab_5y_df$year == year)
      ages <- lifetab_5y_df$age[data_obs]
      Ss <- lifetab_5y_df$lx[data_obs]
      S_fs <- lifetab_5y_df$lx_f[data_obs]
      S_ms <- lifetab_5y_df$lx_m[data_obs]
      
      # Both genders
      if (all(Ss > 0.001)){
        Lstar <- 110
      } else {
        above <- min(ages[which(Ss < 0.001)])
        below <- max(ages[which(Ss >= 0.001)])
        Lstar <- below + (Ss[below+1] - 0.001)/(Ss[below+1] - Ss[above+1])
      }
      # Only Male
      if (all(S_fs > 0.001)){
        Lstar_f <- 110
      } else {
        above <- min(ages[which(S_fs < 0.001)])
        below <- max(ages[which(S_fs >= 0.001)])
        Lstar_f <- below + (S_fs[below+1] - 0.001)/(S_fs[below+1] - S_fs[above+1])
      }
      # Only Female
      if (all(S_ms > 0.001)){
        Lstar_m <- 110
      } else {
        above <- min(ages[which(S_ms < 0.001)])
        below <- max(ages[which(S_ms >= 0.001)])
        Lstar_m <- below + (S_ms[below+1] - 0.001)/(S_ms[below+1] - S_ms[above+1])
      }
      
      
      lifetab_5y_df$l_99p9[data_obs] <- Lstar
      lifetab_5y_df$l_99p9_f[data_obs] <- Lstar_f
      lifetab_5y_df$l_99p9_m[data_obs] <- Lstar_m
    }
  }
  return(lifetab_5y_df)
}
lifetab_5y_df <- compute_lifespan(lifetab_5y_df)
lifetab_df <- compute_lifespan(lifetab_df)



"
Calculate best practice
"
identify_best_practice <- function(lifetab_df){
  # Initialise indicators
  lifetab_df$best_practice <- 0
  lifetab_df$best_practice_gender <- NA
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

rm(convert_df)

" 
Plot BP life expectancy over time, survival and mortality curves
"
# Five year intervals
ggplot(lifetab_5y_export[which(lifetab_5y_export$age == 0 & lifetab_5y_export$year > 1900 &
                                              lifetab_5y_export$best_practice == 1),]) + theme_bw() +
  geom_smooth(aes(x = year, y = ex_f), method = "lm") +
  geom_point(aes(x = year, y = ex_f, color = name)) +
  scale_color_discrete(name = "Country") + ggtitle("Life Expectancy") +
  xlab("Year") + ylab("Life Expectancy at birth")
ggplot(lifetab_5y_export[which(lifetab_5y_export$age == 0 & lifetab_5y_export$year > 1900 &
                                 lifetab_5y_export$best_practice == 1),]) + theme_bw() +
  geom_smooth(aes(x = year, y = -log(Hx_f)), method = "loess") +
  geom_point(aes(x = year, y = -log(Hx_f), color = name)) +
  scale_color_discrete(name = "Country") + ggtitle("Lifespan Equality") +
  xlab("Year") + ylab("Lifespan Equality at birth")
ggplot(lifetab_5y_export[which(lifetab_5y_export$age == 0 & lifetab_5y_export$year > 1900 &
                                 lifetab_5y_export$best_practice == 1),]) + theme_bw() +
  geom_smooth(aes(x = year, y = l_99p9_f), method = "loess") +
  geom_point(aes(x = year, y = l_99p9_f, color = name)) +
  scale_color_discrete(name = "Country") + ggtitle("Lifespan (99.9%)") +
  xlab("Year") + ylab("Lifespan")
ggplot(lifetab_5y_export[which(lifetab_5y_export$best_practice == 1 & 
                                             lifetab_5y_export$year > 1900),]) + theme_bw() +
  geom_line(aes(x = age, y = lx_f, group = year, color = year)) + ylim(c(0,1)) +
  scale_color_gradientn(colours = rainbow(5), name = "Year") + 
  xlab("Age") + ylab("Survival Rate") + ggtitle("Survival Rate")
ggplot(lifetab_5y_export[which(lifetab_5y_export$best_practice == 1 & 
                                             lifetab_5y_export$year > 1900),]) + theme_bw() +
  geom_line(aes(x = age, y = mx_f, group = year, color = year)) + ylim(c(0,1)) +
  scale_color_gradientn(colours = rainbow(5), name = "Year") + 
  xlab("Age") + ylab("Mortality  Rate") + ggtitle("Mortality Rate")




"
Export
"
lifetab_export <- lifetab_export[,c("code", "name", "year", "age",
                                    "mx", "mx_f", "mx_m", "lx", "lx_f", "lx_m", 
                                    "ex", "ex_f", "ex_m", "Hx", "Hx_f", "Hx_m",
                                    "l_99p9", "l_99p9_f", "l_99p9_m",
                                    "Total", "Male", "Female", "best_practice")]
lifetab_5y_export <- lifetab_5y_export[,c("code", "name", "years", "year", "age",
                                          "mx", "mx_f", "mx_m", "lx", "lx_f", "lx_m", 
                                          "ex", "ex_f", "ex_m", "Hx", "Hx_f", "Hx_m",
                                          "l_99p9", "l_99p9_f", "l_99p9_m",
                                          "Total", "Male", "Female","best_practice")]

write.csv(lifetab_5y_export, "data/clean/all_lifetab_5y.csv", row.names = FALSE)
write.csv(lifetab_export, "data/clean/all_lifetab_1y.csv", row.names = FALSE)
#lifetab_5y_export <- read.csv("data/clean/all_lifetab_5y.csv", stringsAsFactors = FALSE)

bp_df <-  lifetab_export[which(lifetab_export$best_practice == 1),]
write.csv(bp_df, "data/clean/bp.csv", row.names = FALSE)

bp_5y_df <-  lifetab_5y_export[which(lifetab_5y_export$best_practice == 1),]
write.csv(bp_5y_df, "data/clean/bp_5y.csv", row.names = FALSE)
#bp_5y_df <- read.csv("data/clean/bp_5y.csv", stringsAsFactors = FALSE)


"
Look at life expectancy at older ages
"

plt_df <- bp_5y_df[which(bp_5y_df$age %in% c(0,10,20,30,40,50,60,70,80,90)),]
plt_df <- plt_df[order(plt_df$year, plt_df$age),]
ggplot(plt_df, aes(x = year)) + theme_bw() + 
  geom_line(aes(y = ex_f, color = as.factor(age))) + 
  scale_color_discrete(name = "Age") + 
  xlab("Year") + ylab("Remaining LE") + ggtitle("Best practice")

rm(plt_df)



# Mortality 
temp_df <- bp_5y_df[which(bp_5y_df$year > 1900),]
temp_df <- temp_df[order(temp_df$year, temp_df$age),which(!str_detect(names(temp_df), "_m"))]
temp_df <- temp_df[,which(!str_detect(names(temp_df), "Ma|Fe|To|best"))]
temp_df <- temp_df[,which(!(names(temp_df) %in% c("ex", "lx", "mx", "Hx")))]

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
# Lifespan
bp_lstar_plt <- ggplot(bp_5y_df[which(bp_5y_df$year > 1900),]) + theme_bw() +
  geom_line(aes(x = age, y = l_99p9_f, group = year, color = year)) + 
  scale_color_gradientn(colours = rainbow(5), name = "Year") + 
  xlab("Age") + ylab("Lifespan") + ggtitle("Lifespan (99.9%)")

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
bp_lstar_plt <- ggplot(bp_5y_df[which(bp_5y_df$age == 0 & bp_5y_df$year > 1900),]) + theme_bw() +
  geom_smooth(aes(x = year, y = l_99p9_f), method = "loess") +
  geom_point(aes(x = year, y = l_99p9_f, color = name)) +
  scale_color_discrete(name = "Country") + ggtitle("Lifespan (99.9%)") +
  xlab("Year") + ylab("Lifespan")
plt2 <- ggarrange(bp_le_plt, bp_h_plt, bp_lstar_plt, nrow = 1, ncol=3, common.legend = TRUE,
          legend = "right")

ggarrange(plt1, plt2, nrow = 1, ncol=2, common.legend = FALSE)
ggsave("figures/data/best_practice_5y_summary.pdf", width = 16, height = 4)





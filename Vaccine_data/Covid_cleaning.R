library(tidyverse)
covid <- read.delim("vaccinations.txt", sep = ",")
#Subsetting the data from the 15 countries used for this analysis:
subset <- covid[covid$location %in% c("Australia", "Belarus", "China", "Denmark", 
                                      "Germany", "Greece", "South Korea", "Russia", 
                                      "Austria", "Turkey", "England", "Ukraine", 
                                      "United States", "Saudi Arabia", "Oman"), ]


#The three dates used for this analysis - it was the most recent date where
#all the countries had sufficient vaccine data
subsetdate <- subset[subset$date %in% c("2022-02-18" ,"2022-02-14", "2022-02-13"), ]

library(dplyr)
finaldf <- subsetdate %>% select(location, date, people_vaccinated_per_hundred)
finaldf <- na.omit(finaldf)

newsub <- subset_redDat[subset_redDat$city %in% c("Muscat", "Riyadh"),]
length(unique(newsub$groupid))
length(unique(subset_redDat$groupid))

# Count the number of unique group IDs for each nation
table <- count(subset_redDat, nation, groupid, city)

# Add a new column to the table that contains the count of unique group IDs
table <- mutate(table, nrgroup = n_distinct(groupid))

# Count the number of unique group IDs for each nation
table_new <- subset_redDat %>%
  group_by(nation) %>%
  summarize(nrgroup = n_distinct(groupid))

#Plotting the 15 countries vaccine rate
subset$date <- as.Date(subset$date)
ggplot(subset, aes(x = date, y = people_vaccinated_per_hundred, color = location)) +
  geom_line() +
  geom_vline(xintercept = as.numeric(as.Date("2022-02-16")), color = "red", linetype = "dashed") +
  theme_classic() +
  ggtitle("Number of people vaccinated per 100")

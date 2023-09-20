library(dplyr)
library(tidyr)
library(ggplot2)



######################## Load plate map and set controls #######################

# Empty vectors for control locations
positive_wells <- vector()
negative_wells <- vector()

# File name for the plate map
filepath <- "DegTopHits_RXSLC1C4_Platemap.csv"


plateMap <- data.frame()
#colnames(plateMap) <- c("Well", "Motif", "Substrate")

#Reading the file and creating the plate map
con = file(filepath, "r")
while ( TRUE ) {
  # Read a line and check its not empty
  line = readLines(con, n = 1)
  if ( length(line) == 0 ) {
    break
  }
  print(line)
  test <- unlist(strsplit(line, ','))
  row_letter <- test[1]
  # Get information for the plate map
  for (i in 2:length(test)){ 
    # Find the positive controls
    well_location <- paste0(row_letter,
                            sprintf("%02d", i - 1))
    if (test[i] == "WT"){
      positive_wells <- append(positive_wells, well_location)
    } else if (test[i] == "Mut"){ 
      # Find the negative controls
      negative_wells <- append(negative_wells, well_location)
    } else if (test[i] != ""){ 
      # Find the well locations and add to plateMap
      plateMap <- rbind(plateMap, append(well_location, test[i]))
    }
  }
}
close(con)

colnames(plateMap) <- c("Well", "Unique_ID")

plateMap <- separate_wider_delim(plateMap, Unique_ID, "_", names = c("Library", "Substitution", "Substrate"))


################################## Graph data ##################################

data <- read.csv("Figure_7.csv")
data2 <- merge(plateMap, data2) %>%
  mutate(Substrate = factor(Substrate, levels = c("Ala", "Arg", "Asn", "Gly", "Gln", "Glu", "His", "fhOrn", "Phe", "Ser", "Thr", "Tyr", "Val", "Xle"))) %>%
  mutate(Substitution = factor(Substitution, levels = c('C1-A10', 'C4-A10')))


# Summarise by substitution type to get average and standard deviation
data_summary <- data2 %>%
  group_by(Substitution, Library, Substrate, Well) %>%
  summarise(average = mean(Yield_per_litre), stDev = sd(Yield_per_litre))
  
# Draw graph
ggplot(data_summary, aes(Substrate, average, fill = Substitution)) +
  geom_col(position = position_dodge(0.6), width=0.5) +  
  geom_errorbar(aes(ymin = average - stDev, ymax = average + stDev), position = position_dodge(0.6), width = 0.2) +
  geom_point(data = data3, aes(Substrate, Yield_per_litre, fill = Substitution), position = position_dodge(0.6), alpha = 0.3, show.legend = FALSE) + 
  theme_classic() +
  facet_wrap(vars(Library)) +
  labs(y = "Yield (mg/l)", x = "") +
  theme(legend.title = element_blank(),
        legend.position = "none",
        text = element_text(size = 10),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_manual(values=c("#FF7171", "#EDCDB7"))




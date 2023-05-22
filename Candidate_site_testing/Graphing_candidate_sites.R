library(dplyr)
library(ggpubr)

library(ggplot2)



######################## Load plate map and set controls #######################

# Empty vectors for control locations
positive_wells <- vector()
negative_wells <- vector()

# File name for the plate map
filepath <- "Motif_recombination_plate_map.csv"


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
      plateMap <- rbind(plateMap, append(well_location,
             unlist(strsplit(test[i], '_'))))
    }
  }
}
close(con)

colnames(plateMap) <- c("Well", "Substrate", "Substitution")


############################# Load data and rearrange ##########################

# Function for reading each plate
read_plates <- function(replicate, positive_wells, negative_wells){
  # Read the file
  filename = paste0("Replicate_", replicate, ".csv")
  file_input <- read.csv(filename, header = FALSE, skip = 12, nrows = 384)[c(TRUE, FALSE, FALSE, FALSE), ]
  
  # Get the raw data and label columns
  output <- data.frame(file_input$V1,
                       replicate,
                       as.numeric(file_input$V2))
  colnames(output) <- c("Well", "Replicate", "Raw")

  # Get the average for the negative controls
  data_neg <- output %>%
    filter(Well %in% negative_wells) %>%
    select(Raw) %>%
    colMeans()
  
  # Get the average for the positive controls
  data_pos <- output %>%
    filter(Well %in% c(positive_wells)) %>%
    select(Raw) %>%
    colMeans()
  
  # Calculate the percentage relative to positives
  output$Percentage <- (output$Raw-data_neg)/(data_pos-data_neg)*100

  output
}

# Empty data frame to add each plate to
unlabelledData <- data.frame()

# Go through each file and add to the dataframe
for(replicate in 1:3){
  unlabelledData <- rbind(unlabelledData, read_plates(replicate, positive_wells, negative_wells))
}


################## Combine plateMap with the unlabelled data ###################

# Merge platemap with data and add factoring to set the correct ordering
data <- merge(plateMap, unlabelledData) %>%
  mutate(Substrate = factor(Substrate, levels = c('Gly', 'fhOrn', 'Ser', 'Ala', 'Arg1', 'Arg2', 'Glu', 'Phe'))) %>%
  mutate(Substitution = factor(Substitution, levels = c('C1-A10', 'C3-A10', 'C4-A10', 'C7-A10', '1-A10', '2-A10', 'A2-A10', 'A3-A10',
                                                        'X-A6', 'X-A7', 'X-A8a', 'X-A8b', 'X-A8c', 'X-A10')))


################################## Graph data ##################################

# A function to make the graphs consistent with each other
standard_graph <- function(data, legend_x){

  # Summarise by substitution type to get average and standard deviation
  data_summary <- data %>%
    group_by(Substitution, Substrate) %>%
    summarise(average = mean(Percentage), stDev = sd(Percentage))
  
  # Draw graph
  graph <- ggplot(data_summary, aes(Substitution, average, fill = Substrate)) +
    geom_col(position = position_dodge(0.6), width=0.5, alpha = 0.8) +  
    geom_errorbar(aes(ymin = average - stDev, ymax = average + stDev), position = position_dodge(0.6), width = 0.2) +
    geom_point(data = data, aes(Substitution, Percentage, fill = Substrate), position = position_dodge(0.6), alpha = 0.3, show.legend = FALSE) + 
    theme_classic() +
    scale_y_continuous(breaks=seq(0, 100, 25)) +
    labs(y = "Yield relative to wildtype (%)", x = "Substituted region") +
    theme(legend.title = element_blank(),
          legend.position = c(legend_x, 0.8),
          text = element_text(size = 20)) +
    coord_cartesian(ylim = c(-5, 110)) +
    scale_fill_manual(breaks = c('Gly', 'fhOrn', 'Ser', 'Ala', 'Arg1', 'Arg2', 'Glu', 'Phe'), values=c("#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD", "#8C564B", "#E377C2", "#7F7F7F"))
  
  graph
  
}

# Make a graph of the first substitutions
graph_1 <- data %>%
  filter(Substitution %in% c('C1-A10', 'C3-A10', 'C4-A10', 'C7-A10', '1-A10', '2-A10', 'A2-A10', 'A3-A10')) %>%
  filter(Substrate %in% c('Gly', 'fhOrn', 'Ser')) %>%
  standard_graph(legend_x = 0.1)

png(filename = "graph_1.png", width = 760, height = 560)
graph_1 # Print to file
dev.off()

# Make a graph of the second substitutions
graph_2 <- data %>%
  filter(Substitution %in% c('X-A6', 'X-A7', 'X-A8a', 'X-A8b', 'X-A8c', 'X-A10')) %>%
  standard_graph(legend_x = 0.128)

png(filename = "graph_2.png", width = 590, height = 560)
graph_2 # Print to file
dev.off()

# Make a graph of the third substitutions
graph_3 <- data %>%
  filter(Substitution %in% c('C1-A10', 'C4-A10')) %>%
  filter(!Substrate %in% c('Gly', 'fhOrn', 'Ser')) %>%
  standard_graph(legend_x = 0.18)

png(filename = "graph_3.png", width = 420, height = 560)
graph_3 # Print to file
dev.off()

################################# Run t tests ##################################

# Comparisons between 'C1-A10' and 'C4-A10' for each substrate

sapply(c('Gly', 'fhOrn', 'Ser', 'Ala', 'Arg1', 'Arg2', 'Glu', 'Phe'), function(x, control_site = 'C1-A10', test_site = 'C4-A10'){
  control <- data$Percentage[data$Substrate == x & data$Substitution == control_site]
  tester <- data$Percentage[data$Substrate == x & data$Substitution == test_site]
  t.test(control, tester)$p.value < 0.05
})





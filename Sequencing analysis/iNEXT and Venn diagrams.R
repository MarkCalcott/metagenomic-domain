
library("iNEXT")
library(ggplot2)
library(dplyr)
library(gplots)
library(webr)

################################ Load OTU data #################################

# Import OTU table
otu_tab <- read.csv(file = 'Clustering/processed_reverse_sequences_otutab.txt', header = TRUE, 
                    sep = '\t', row.names = 1)

# Import blast data
otu_blast <-  read.csv(file = 'OTU blast/OTU_Blast.csv', header = TRUE, row.names = 1)

# merge the two data sets
otu_complete_data <- merge(otu_tab, otu_blast, by = 'row.names')

######################## Adding columns for comparisons ########################

# Add columns for comparisons
otu_complete_data <- cbind(otu_complete_data, RX = rowSums(select(otu_complete_data, matches("U.RX")))) # Unselected RX sequences
otu_complete_data <- cbind(otu_complete_data, "RX hits" = rowSums(select(otu_complete_data, matches("H.RX")))) # Selected RX sequences
otu_complete_data <- cbind(otu_complete_data, SL = rowSums(select(otu_complete_data, matches("U.SL")))) # Unselected SL sequences
otu_complete_data <- cbind(otu_complete_data, "SL hits" = rowSums(select(otu_complete_data, matches("H.SL")))) # Selected SL sequences
otu_complete_data <- cbind(otu_complete_data, "C1-A10" = rowSums(select(otu_complete_data, matches("UA+")))) # Unselected C1-A10 sequences
otu_complete_data <- cbind(otu_complete_data, "C4-A10" = rowSums(select(otu_complete_data, matches("U4+")))) # Unselected C4-A10 sequences

# Unselected RX1 - RX4 sequences
otu_complete_data <- cbind(otu_complete_data, RX1 = rowSums(select(otu_complete_data, matches("U.RX1"))))
otu_complete_data <- cbind(otu_complete_data, RX2 = rowSums(select(otu_complete_data, matches("U.RX2"))))
otu_complete_data <- cbind(otu_complete_data, RX3 = rowSums(select(otu_complete_data, matches("U.RX3"))))
otu_complete_data <- cbind(otu_complete_data, RX4 = rowSums(select(otu_complete_data, matches("U.RX4"))))

# Unselected SL1 - SL4 sequences
otu_complete_data <- cbind(otu_complete_data, SL1 = rowSums(select(otu_complete_data, matches("U.SL1"))))
otu_complete_data <- cbind(otu_complete_data, SL2 = rowSums(select(otu_complete_data, matches("U.SL2"))))
otu_complete_data <- cbind(otu_complete_data, SL3 = rowSums(select(otu_complete_data, matches("U.SL3"))))
otu_complete_data <- cbind(otu_complete_data, SL4 = rowSums(select(otu_complete_data, matches("U.SL4"))))


######################## Small stats on OTUs for paper ########################

sum(colSums(select(otu_tab, matches("U.")))) # Sum of the total sequences unselected
min(colSums(select(otu_tab, matches("U.")))) # Min OTUs in a sublibrary
max(colSums(select(otu_tab, matches("U.")))) # Max OTUs in a sublibrary
sum(rowSums(select(otu_tab, matches("U."))) > 0) # Count OTUs for unselected
sum(rowSums(select(otu_tab, matches("H."))) > 0) # Count OTUs for selected

sum(colSums(select(otu_tab, matches("H."))))

##################### Estimate OTU diversity and draw graph ####################

# Calculate species diversity for the SL and RX data
out <- otu_complete_data %>%
  select(SL, RX) %>%
  iNEXT(., q = 0, datatype="abundance", endpoint=3000)
out

# Settings for drawing graph
df <- fortify(out, type=1)
df.point <- df[which(df$Method=="Observed"),]
df.line <- df[which(df$Method!="Observed"),]
df.line$Method <- factor(df.line$Method, levels = c("Rarefaction", "Extrapolation"), labels = c("Interpolation", "Extrapolation"))

# Draw graph
graph <- ggplot(df, aes(x = x, y = y, colour = Assemblage)) +
  geom_point(aes(shape = Assemblage), size = 5, data = df.point) +
  geom_ribbon(aes(ymin = y.lwr, ymax = y.upr, fill = Assemblage, colour = NULL), alpha=0.2) +
  geom_line(aes(linetype = Method), lwd = 1.5, data = df.line) +
  scale_colour_manual(values = c("#009900", "#0031CC")) + # Colour line
  scale_fill_manual(values = c("#00FE00", "blue")) + # Colour ribbon
  labs(x = "Sample size", y = "OTU diversity") +
  theme_classic() +
  theme(legend.position = c(0.8, 0.27),
        text = element_text(size = 30))

# Make a graph of OTU diversity
png(filename = "OTU_diversity.png", width = 760, height = 560)
graph # Print to file
dev.off()

############################## Make Venn diagrams ##############################

### Compare unselected sequences from RX and SL
# Draw graph
comparison = c("RX", "SL")
venn(sapply(comparison, function(x) {x=rownames(otu_complete_data)[otu_complete_data[, x]>0]}))



### Compare unselected sequences from CA and C4
# Draw graph
comparison = c("C1-A10", "C4-A10")
venn(sapply(comparison, function(x) {x=rownames(otu_complete_data)[otu_complete_data[, x]>0]}))


### Compare unselected sequences from RX sublibraries
# Draw graph
comparison = c("RX1", "RX2", "RX3", "RX4")
venn(sapply(comparison, function(x) {x=rownames(otu_complete_data)[otu_complete_data[, x]>0]}))


### Compare unselected sequences from SL sublibraries
# Draw graph
comparison = c("SL1", "SL2", "SL3", "SL4")
venn(sapply(comparison, function(x) {x=rownames(otu_complete_data)[otu_complete_data[, x]>0]}))


#### Compare unselected and selected otus for RX
# Draw graph
comparison = c("RX hits", "RX")
venn(sapply(comparison, function(x) {x=rownames(otu_complete_data)[otu_complete_data[, x]>0]}))


#### Compare unselected and selected otus for SL
# Draw graph
comparison = c("SL hits", "SL")
venn(sapply(comparison, function(x) {x=rownames(otu_complete_data)[otu_complete_data[, x]>0]}))


######################## Donut graphs of taxonomic data ########################

# Rearrange to highlight Phyla and genus for Streptomyces and Pseudomonas
otu_complete_data$Phyla_PS <- factor(otu_complete_data$Phyla, levels = c("Actinomycetota", "Pseudomonadota", "Other"))
otu_complete_data$Genus_PS <- factor(otu_complete_data$Genus, levels = c("Streptomyces", "Pseudomonas", ""))

otu_complete_data$Phyla_PS[is.na(otu_complete_data$Phyla_PS)] <- "Other"
otu_complete_data$Genus_PS[is.na(otu_complete_data$Genus_PS)] <- ""

# Make graph for SL
otu_complete_data %>%
  filter(SL > 0) %>%
  PieDonut(aes(Phyla_PS, Genus_PS), ratioByGroup = FALSE)

# Make graph for RX
otu_complete_data %>%
  filter(RX > 0) %>%
  PieDonut(aes(Phyla_PS, Genus_PS), ratioByGroup = FALSE)


##### Some small descriptive states on the unselected data
SL_RX_unselected <- otu_complete_data %>%
  filter(SL > 0 | RX > 0)

length(unique(SL_RX_unselected$Genus)) - 1 # How many genera in the unselected, minus one to remove ""
length(unique(SL_RX_unselected$Phyla)) # How many Phyla in the unselected
sum(SL_RX_unselected$Phyla == "Pseudomonadota") # How many Pseudomonadota
sum(SL_RX_unselected$Phyla == "Actinomycetota") # How many Actinomycetota



       
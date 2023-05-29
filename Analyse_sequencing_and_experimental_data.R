library(dplyr)
library(ggplot2)


# Load experimental data
experimental_data <- read.csv("Screening/Absorbance_and_MALDI/Screen_data.csv")

# Remove samples with no substrate prediction
experimental_data <- experimental_data %>%
  filter(!is.na(Substrate), !is.na(Library))

# Load sequence data
sequence_data <- read.csv("Sequencing analysis/Sequence_data.csv")

# Combine sequence and experimental data
combined_data <- merge(experimental_data, sequence_data)



################### Frequency of each hit during screening #####################

# Create a graph for showing the average for each hit
graph_points <- ggplot(experimental_data, aes(Substrate, Abs, color = Substitution)) +
  geom_point(position = position_dodge(0.6), size = 3, alpha = 0.7) +
  facet_wrap(vars(Library)) +
  theme_bw() +
  labs(y = "Average yield for each hit (%)", x = element_blank()) + 
  theme(legend.position = "none",
        text = element_text(size = 20),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_color_manual(values=c("#FF7171", "#EDCDB7"))

png(filename = "graph_points.png", width = 820, height = 560)
graph_points # Print to file
dev.off()

# Create a graph showing boxplot of each
graph_boxplot <- combined_data %>%
  filter(!Substrate == "Val") %>%
  ggplot(aes(Substitution, Abs, color = Substitution)) +
  geom_boxplot(position = position_dodge(0.6)) +
  facet_wrap(vars(Library)) +
  labs(y = "Average yield for each hit (%)", x = element_blank()) + 
  theme(legend.position = "none", text = element_text(size = 20),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme_bw() +
  scale_color_manual(values=c("#FF7171", "#EDCDB7"))

png(filename = "graph_boxplot.png", width = 350, height = 560)
graph_boxplot # Print to file
dev.off()

# Create data to graph relationships
graph_rev_PID <- combined_data %>%
  ggplot(aes(rev_PID, Abs)) +
  geom_smooth(method='lm') +
  geom_point(aes(color = Substrate)) +
  labs(y = "Average yield for each hit (%)", 
       x = "Reverse sequence identity (DNA)") +
  theme_bw() + 
  theme(legend.position = "none", text = element_text(size = 15),
        axis.text.x = element_text(angle = 45, hjust = 1)) 
png(filename = "graph_rev_PID.png", width = 400, height = 400)
graph_rev_PID # Print to file
dev.off()

# Create of rev_PID
graph_fwd_PID <- combined_data %>%
  ggplot(aes(fwd_PID, Abs)) +
  geom_smooth(method='lm') +
  geom_point(aes(color = Substrate)) +
  labs(y = "Average yield for each hit (%)", 
       x = "Forward sequence identity (DNA)") +
  theme_bw() + 
  theme(legend.position = "none", text = element_text(size = 15),
        axis.text.x = element_text(angle = 45, hjust = 1)) 
png(filename = "graph_fwd_PID.png", width = 400, height = 400)
graph_fwd_PID # Print to file
dev.off()

graph_rev_aaPID <- combined_data %>%
  ggplot(aes(rev_aaPID, Abs)) +
  geom_smooth(method='lm') +
  geom_point(aes(color = Substrate)) +
  labs(y = "Average yield for each hit (%)", 
       x = "Reverse sequence identity (amino acids)") +
  theme_bw() + 
  theme(legend.position = "none", text = element_text(size = 15),
        axis.text.x = element_text(angle = 45, hjust = 1)) 
png(filename = "graph_rev_aaPID.png", width = 400, height = 400)
graph_rev_aaPID # Print to file
dev.off()

graph_fwd_aaPID <- combined_data %>%
  ggplot(aes(fwd_aaPID, Abs)) +
  geom_smooth(method='lm') +
  geom_point(aes(color = Substrate)) +
  labs(y = "Average yield for each hit (%)", 
       x = "Forward sequence identity (amino acids)") +
  theme_bw() + 
  theme(legend.position = "none", text = element_text(size = 15),
        axis.text.x = element_text(angle = 45, hjust = 1)) 
png(filename = "graph_fwd_aaPID.png", width = 400, height = 400)
graph_fwd_aaPID # Print to file
dev.off()

graph_rev_HEG <- combined_data %>%
  ggplot(aes(rev_HEG, Abs)) +
  geom_smooth(method='lm') +
  geom_point(aes(color = Substrate)) +
  labs(y = "Average yield for each hit (%)", 
       x = "Reverse codon adaptiveness (HEG)") +
  theme_bw() + 
  theme(legend.position = "none", text = element_text(size = 15),
        axis.text.x = element_text(angle = 45, hjust = 1)) 
png(filename = "graph_rev_HEG.png", width = 400, height = 400)
graph_rev_HEG # Print to file
dev.off()

# Create of rev_HEG
graph_fwd_HEG <- combined_data %>%
  ggplot(aes(fwd_HEG, Abs)) +
  geom_smooth(method='lm') +
  geom_point(aes(color = Substrate)) +
  labs(y = "Average yield for each hit (%)", 
       x = "Forward codon adaptiveness (HEG)") +
  theme_bw() + 
  theme(legend.position = "none", text = element_text(size = 15),
        axis.text.x = element_text(angle = 45, hjust = 1)) 
png(filename = "graph_fwd_HEG.png", width = 400, height = 400)
graph_fwd_HEG # Print to file
dev.off()

graph_rev_RPG <- combined_data %>%
  ggplot(aes(rev_RPG, Abs)) +
  geom_smooth(method='lm') +
  geom_point(aes(color = Substrate)) +
  labs(y = "Average yield for each hit (%)", 
       x = "Reverse codon adaptiveness (RPG)") +
  theme_bw() + 
  theme(legend.position = "none", text = element_text(size = 15),
        axis.text.x = element_text(angle = 45, hjust = 1)) 
png(filename = "graph_rev_RPG.png", width = 400, height = 400)
graph_rev_RPG # Print to file
dev.off()

graph_fwd_RPG <- combined_data %>%
  ggplot(aes(fwd_RPG, Abs)) +
  geom_smooth(method='lm') +
  geom_point(aes(color = Substrate)) +
  labs(y = "Average yield for each hit (%)", 
       x = "Forward codon adaptiveness (RPG)") +
  theme_bw() + 
  theme(legend.position = "none", text = element_text(size = 15),
        axis.text.x = element_text(angle = 45, hjust = 1)) 
png(filename = "graph_fwd_RPG.png", width = 400, height = 400)
graph_fwd_RPG # Print to file
dev.off()

graph_rev_GC <- combined_data %>%
  ggplot(aes(rev_GC, Abs)) +
  geom_smooth(method='lm') +
  geom_point(aes(color = Substrate)) +
  labs(y = "Average yield for each hit (%)", 
       x = "Reverse GC content") +
  theme_bw() + 
  theme(legend.position = "none", text = element_text(size = 15),
        axis.text.x = element_text(angle = 45, hjust = 1)) 
png(filename = "graph_rev_GC.png", width = 400, height = 400)
graph_rev_GC # Print to file
dev.off()

graph_fwd_GC <- combined_data %>%
  ggplot(aes(fwd_GC, Abs)) +
  geom_smooth(method='lm') +
  geom_point(aes(color = Substrate)) +
  labs(y = "Average yield for each hit (%)", 
       x = "Forward GC content") +
  theme_bw() + 
  theme(legend.position = "none", text = element_text(size = 15),
        axis.text.x = element_text(angle = 45, hjust = 1)) 
png(filename = "graph_fwd_GC.png", width = 400, height = 400)
graph_fwd_GC # Print to file
dev.off()






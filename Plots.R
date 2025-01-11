# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(ggpubr)
# Load the RMSD data
# setwd()
data <- read.csv("project_data/RMSD_results.csv", sep=',')  # Replace with your actual file path

# Add classification for A and B
data <- data %>%
  mutate(
    Closer_To_A = ifelse(HOLO_vs_APO_ALPHA0 < APO_vs_APO_ALPHA0, "Holo", "Apo"),
    Closer_To_B = ifelse(HOLO_vs_HOLO_ALPHA0 < APO_vs_HOLO_ALPHA0, "Holo", "Apo")
  )

# FIG 1
figure1 <-  
  ggplot() +
  geom_density(data = data,
               aes(x = APO_vs_HOLO), 
               alpha = 0.3, color = 'black',
               fill = 'black',
               adjust = 1) + 
  xlab(expression(RMSD[Apo~vs~Holo]~(Å))) + ylab('Density') +
  
  scale_x_continuous(limits = c(0,15.5),breaks=seq(0, 15, 1), expand = c(0,0)) +
  scale_y_continuous(limits = c(0,0.25), breaks = seq(0.000, 0.25, 0.05), expand = c(0,0)) +
  
  theme(legend.position = "none") +
  
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) +
  
  theme(panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_blank()) +
  
  theme(axis.title.x = element_text(face ='plain',size = 13, margin = margin(t = 0, r = 8, b = 0, l = 0)),
        axis.title.y = element_text(face ='plain',size = 13, margin = margin(t = 0, r = 8, b = 0, l = 0))) +
  
  theme(axis.text=element_text(size = 10))+
  theme(plot.margin = unit(c(.2,.2,.2,.2), "cm")) 

summary(data$APO_vs_HOLO)

print(figure1)
ggsave("Stats/fig0_density.png", plot = figure1, width = 8, height = 6)

plddt_data <- read.csv("project_data/plddt_results.tsv", sep='\t')  # Replace with your actual file path

# plDDT scores
plddt_mean = mean(plddt_data$PLDDT)
print(plddt_mean)
print(count(plddt_data))
print((count(plddt_data[plddt_data$PLDDT >= 85,]) / count(plddt_data)) * 100)

# FIG 2A - APO_ALPHA TO APO AND HOLO
figure2a <- data %>% 
  ggplot(aes(x = APO_vs_APO_ALPHA0, y = HOLO_vs_APO_ALPHA0, color = Closer_To_A)) +
  geom_point(size = 4, alpha = 0.7) +
  scale_color_manual(values = c("Holo" = "cyan", "Apo" = "red")) +
  labs(x = expression(RMSD[~vs~Apo]~(Å)),
       y = expression(RMSD[~vs~Holo]~(Å))) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = 'red') +
  scale_x_continuous(breaks = seq(0, 15, 1), limits = c(0, 15), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 15, 1), limits = c(0, 15), expand = c(0, 0)) +
  theme_minimal() +
  ggtitle('A')

print(figure2a)

ggsave("Stats/fig1_Alpha_vs_Holo_vs_Apo.png", plot = figure2a, width = 8, height = 6)

# FIG 2B - HOLO_ALPHA TO APO AND HOLO
figure2b <- data %>% 
  ggplot(aes(x = APO_vs_HOLO_ALPHA0, y = HOLO_vs_HOLO_ALPHA0, color = Closer_To_B)) +
  geom_point(size = 4, alpha = 0.7) +
  scale_color_manual(values = c("Holo" = "cyan", "Apo" = "red")) +
  labs(x = expression(RMSD[~vs~Apo]~(Å)),
       y = expression(RMSD[~vs~Holo]~(Å))) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = 'red') +
  scale_x_continuous(breaks = seq(0, 15, 1), limits = c(0, 15), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 15, 1), limits = c(0, 15), expand = c(0, 0)) +
  theme_minimal() +
  ggtitle('B')

print(figure2b)
ggsave("Stats/fig2_HOLO_ALPHA_VS_HOLO_VS_APO.png", plot = figure2b, width = 8, height = 6)

# Wilcoxon Signed-Rank Test (Right-tailed)
# Ar HOLO su apo_alpha yra reiksmingai mažesnis RMSD
wilcoxon_result <- wilcox.test(
  c(data$HOLO_vs_APO_ALPHA0), 
  c(data$APO_vs_APO_ALPHA0),
  alternative = "less"
)
# W = 109, p-value = 0.006597
# Print the results
cat("Wilcoxon Signed-Rank Test:\n")
print(wilcoxon_result)

# Ar HOLO su holo_alpha yra reiksmingai mazesnis RMSD
# Pasirodo nera.
wilcoxon_result <- wilcox.test(
  c(data$HOLO_vs_HOLO_ALPHA0), 
  c(data$APO_vs_HOLO_ALPHA0),
  alternative = "less"
)
# W = 141, p-value = 0.05713
# Print the results
cat("Wilcoxon Signed-Rank Test:\n")
print(wilcoxon_result)

# FIG 3A
data$Average_RMSD <- rowMeans(data[, c("HOLO_vs_APO_ALPHA0", "APO_vs_APO_ALPHA0")], na.rm = TRUE)
# Preprocess for Holo group (Figure A)
holo_data <- data %>%
  filter(Closer_To_A == "Holo") %>%  # Filter for Holo
  select(HOLO_vs_HOLO_ALPHA0, APO_vs_HOLO_ALPHA0, Closer_To_A) %>%  # Use exact column names
  pivot_longer(
    cols = c(HOLO_vs_HOLO_ALPHA0, APO_vs_HOLO_ALPHA0),
    names_to = "RMSD_Type",
    values_to = "RMSD_Value"
  ) %>%
  mutate(RMSD_Type = ifelse(RMSD_Type == "HOLO_vs_HOLO_ALPHA0", "Holo vs HOLO Alpha", "APO vs HOLO Alpha"))

# Plot A
figure_holo <- ggplot(holo_data, aes(x = RMSD_Type, y = RMSD_Value, fill = RMSD_Type)) +
  geom_boxplot() +
  labs(
    title = "Proteins Closer to Holo",
    x = "",
    y = expression(RMSD~(Å))
  ) +
  theme_minimal() +
  scale_y_continuous(breaks = seq(0, 15, 2), limits = c(0, 15), expand = c(0, 0)) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", size = 20, hjust = 0.5),
    axis.title.y = element_text(face = "plain", size = 20, margin = margin(t = 0, r = 10, b = 0, l = 0)),
    axis.text = element_text(size = 13),
    axis.ticks = element_line(size = 0.5),
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
    strip.text = element_text(size = 15, face = "plain")
  ) +
  scale_fill_manual(values = c("orange", "green"))+
  theme_minimal() +
  ggtitle('A')
print(figure_holo)
ggsave("Stats/closer_holo_holo_alpha.png", plot = figure_holo, width = 8, height = 6)

# FIG 3B
holo_data <- data %>%
  filter(Closer_To_A == "Holo") %>%  # Filter for Holo
  select(HOLO_vs_APO_ALPHA0, APO_vs_APO_ALPHA0, Closer_To_A) %>%  # Use exact column names
  pivot_longer(
    cols = c(HOLO_vs_APO_ALPHA0, APO_vs_APO_ALPHA0),
    names_to = "RMSD_Type",
    values_to = "RMSD_Value"
  ) %>%
  mutate(RMSD_Type = ifelse(RMSD_Type == "HOLO_vs_APO_ALPHA0", "Holo vs APO Alpha", "APO vs APO Alpha"))

# B plot
figure_holo <- ggplot(holo_data, aes(x = RMSD_Type, y = RMSD_Value, fill = RMSD_Type)) +
  geom_boxplot() +
  labs(
    title = "Proteins Closer to Holo",
    x = "",
    y = expression(RMSD~(Å))
  ) +
  theme_minimal() +
  scale_y_continuous(breaks = seq(0, 15, 2), limits = c(0, 15), expand = c(0, 0)) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", size = 20, hjust = 0.5),
    axis.title.y = element_text(face = "plain", size = 20, margin = margin(t = 0, r = 10, b = 0, l = 0)),
    axis.text = element_text(size = 13),
    axis.ticks = element_line(size = 0.5),
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
    strip.text = element_text(size = 15, face = "plain")
  ) +
  scale_fill_manual(values = c("orange", "green"))+
  theme_minimal() +
  ggtitle('B')
print(figure_holo)
ggsave("Stats/closer_holo_apo_alpha.png", plot = figure_holo, width = 8, height = 6)

# FIG 3C
apo_data <- data %>%
  filter(Closer_To_A == "Apo") %>%  # Filter for Apo
  select(HOLO_vs_APO_ALPHA0, APO_vs_APO_ALPHA0, Closer_To_A) %>%  # Use exact column names
  pivot_longer(
    cols = c(HOLO_vs_APO_ALPHA0, APO_vs_APO_ALPHA0),
    names_to = "RMSD_Type",
    values_to = "RMSD_Value"
  ) %>%
  mutate(RMSD_Type = ifelse(RMSD_Type == "HOLO_vs_APO_ALPHA0", "Holo vs APO Alpha", "APO vs APO Alpha"))

# Plot C
figure_apo <- ggplot(apo_data, aes(x = RMSD_Type, y = RMSD_Value, fill = RMSD_Type)) +
  geom_boxplot() +
  labs(
    title = "Proteins Closer to Apo",
    x = "",
    y = expression(RMSD~(Å))
  ) +
  theme_minimal() +
  scale_y_continuous(breaks = seq(0, 15, 2), limits = c(0, 15), expand = c(0, 0)) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", size = 20, hjust = 0.5),
    axis.title.y = element_text(face = "plain", size = 20, margin = margin(t = 0, r = 10, b = 0, l = 0)),
    axis.text = element_text(size = 13),
    axis.ticks = element_line(size = 0.5),
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
    strip.text = element_text(size = 15, face = "plain")
  ) +
  scale_fill_manual(values = c("orange", "green"))+
  theme_minimal() +
  ggtitle('C')
# Print both figures
print(figure_apo)   # Figure for Apo group
ggsave("Stats/closer_apo_apo_alpha.png", plot = figure_holo, width = 8, height = 6)

# FIG 3D
apo_data <- data %>%
  filter(Closer_To_A == "Apo") %>%  # Filter for Apo
  select(HOLO_vs_HOLO_ALPHA0, APO_vs_HOLO_ALPHA0, Closer_To_A) %>%  # Use exact column names
  pivot_longer(
    cols = c(HOLO_vs_HOLO_ALPHA0, APO_vs_HOLO_ALPHA0),
    names_to = "RMSD_Type",
    values_to = "RMSD_Value"
  ) %>%
  mutate(RMSD_Type = ifelse(RMSD_Type == "HOLO_vs_HOLO_ALPHA0", "Holo vs HOLO Alpha", "APO vs HOLO Alpha"))

# Plot D
figure_apo <- ggplot(apo_data, aes(x = RMSD_Type, y = RMSD_Value, fill = RMSD_Type)) +
  geom_boxplot() +
  labs(
    title = "Proteins Closer to Apo",
    x = "",
    y = expression(RMSD~(Å))
  ) +
  theme_minimal() +
  scale_y_continuous(breaks = seq(0, 15, 2), limits = c(0, 15), expand = c(0, 0)) +
  theme(
    legend.position = "none",
    plot.title = element_text(face = "bold", size = 20, hjust = 0.5),
    axis.title.y = element_text(face = "plain", size = 20, margin = margin(t = 0, r = 10, b = 0, l = 0)),
    axis.text = element_text(size = 13),
    axis.ticks = element_line(size = 0.5),
    plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"),
    strip.text = element_text(size = 15, face = "plain")
  ) +
  scale_fill_manual(values = c("orange", "green"))+
  theme_minimal() +
  ggtitle('D')

print(figure_apo)
ggsave("Stats/closer_apo_holo_alpha.png", plot = figure_holo, width = 8, height = 6)

# FIG 4
# Reshape data to have PDB_ID and corresponding APO_vs_HOLO
reshaped_rmsd <- data %>%
  select(Apo_PDB_ID, Holo_PDB_ID, APO_vs_HOLO) %>% 
  mutate(APO_vs_HOLO = as.numeric(gsub(",", ".", APO_vs_HOLO))) %>%  # Ensure numeric format
  pivot_longer(
    cols = c(Apo_PDB_ID, Holo_PDB_ID),
    names_to = "PDB_Type",
    values_to = "PDB_ID"
  )

# Merge with plDDT data
merged_data <- reshaped_rmsd %>%
  inner_join(plddt_data, by = c("PDB_ID" = "PDB_NAME"))

# Create the scatter plot
figure <- ggplot(merged_data, aes(x = APO_vs_HOLO, y = PLDDT)) +
  geom_point(fill = 'black', color = 'black', alpha = 0.3, size = 4) +
  geom_smooth(method = "lm", color = "blue") +
  
  # Axis labels and titles
  labs(
    x = expression(RMSD[Apo~vs~Holo]~(Å)),
    y = expression(plDDT)
  ) +
  
  # Axis scales
  scale_x_continuous(breaks = seq(0, 15.5, 1), limits = c(0, 15.5), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 100, 10), limits = c(0, 100), expand = c(0, 0)) +
  
  # Theme customization
  theme(
    axis.line = element_line(colour = "black"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    legend.position = "none",
    axis.title.x = element_text(face = 'plain', size = 25, margin = margin(t = 10, r = 8, b = 0, l = 0)),
    axis.title.y = element_text(face = 'plain', size = 25, margin = margin(t = 0, r = 8, b = 0, l = 0)),
    axis.text = element_text(size = 20),
    plot.margin = unit(c(.2, .2, .2, .2), "cm"),
    plot.title = element_text(face = 'bold', vjust = 0, size = 20)
  ) +
  
  ggtitle("Scatter Plot of plDDT vs. RMSD (Apo vs. Holo)")

print(figure)
ggsave("Stats/fig5_plDDT_RMSD.png", plot = figure, width = 8, height = 6)

# p-value = 0.6583, -0.1157285
cor.test(merged_data$APO_vs_HOLO, merged_data$PLDDT,
         method = 'pearson')

# FIG 5A
data_families <- read.csv("project_data/family_RMSD.csv", sep=';')  # Replace with your actual file path

# Fix columns with commas as decimal separators
data_families <- data_families %>%
  mutate(
    # Apo_PDB_ID = gsub("-", "_", Apo_PDB_ID),
    # Holo_PDB_ID = gsub("-", "_", Holo_PDB_ID),
    low_RMSD = as.numeric(gsub(",", ".", low_RMSD)),
    avg_RMSD = as.numeric(gsub(",", ".", avg_RMSD)),
    high_RMSD = as.numeric(gsub(",", ".", high_RMSD)))

data_families <- data_families %>% mutate(family = ifelse((high_RMSD-low_RMSD) < 4, "Homogeneous", "Heterogeneous"))

merged_family_data <- merge(data, data_families, by = c("Apo_PDB_ID", "Holo_PDB_ID"), all=TRUE)

# Calculate the minimum of HOLO_vs_HOLO_ALPHA0 and APO_vs_HOLO_ALPHA0
merged_family_data <- merged_family_data %>%
  rowwise() %>%
  mutate(Min_RMSD_Value = min(HOLO_vs_HOLO_ALPHA0, APO_vs_HOLO_ALPHA0, na.rm = TRUE))
y_range <- range(merged_family_data$Min_RMSD_Value, na.rm = TRUE)
# Create the box plot
figure5A <- ggplot(merged_family_data, aes(x = family, y = Min_RMSD_Value)) +
  geom_boxplot() +
  labs(
    title = 'A',
    # title = "Box Plot of lowest RMSD Values by Family Group",
    # x = "Family Group",
    y = "lowest RMSD"
  ) +
  scale_y_continuous(limits = y_range) +
  theme_minimal()

ggsave("Stats/RMSD_vs_families_holo.png", plot = figure5A, width = 8, height = 6)


wilcox_test <- wilcox.test(
  Min_RMSD_Value ~ family, 
  data = merged_family_data,
  na.action = na.omit,  # Exclude NA values
  alternative = "greater"
)

# Print the Wilcoxon test results
print("Wilcoxon Rank-Sum Test Results:")
print(wilcox_test)


# FIG 5B
merged_family_data <- merged_family_data %>%
  rowwise() %>%
  mutate(Min_RMSD_Value = min(HOLO_vs_APO_ALPHA0, APO_vs_APO_ALPHA0, na.rm = TRUE))

# Create the box plot
figure5b <- ggplot(merged_family_data, aes(x = family, y = Min_RMSD_Value)) +
  geom_boxplot() +
  labs(
    title = 'B',
    # title = "Box Plot of lowest RMSD Values by Family Group(Apo Alpha)",
    # x = "Family Group",
    y = "lowest RMSD"
  ) +
  scale_y_continuous(limits = y_range) +
  theme_minimal()

ggsave("Stats/RMSD_vs_families_apo.png", plot = figure5B, width = 8, height = 6)

wilcox_test <- wilcox.test(
  Min_RMSD_Value ~ family, 
  data = merged_family_data,
  na.action = na.omit,  # Exclude NA values
  alternative = "greater"
)

# Print the Wilcoxon test results
print("Wilcoxon Rank-Sum Test Results:")
print(wilcox_test)

## CODE FOR CREATING UNIFIED AND FACETED TRIANGLE PLOTS USING ADMIXTURE (q) AND INTERSOURCE ANCESTRY (Q)

## Based on code created by Amanda Meuser
## Last edited by Jenna LeBlanc, February 23, 2026
## Note: portions of this code have been supplemented by ChatGPT and/or Seqera AI. Jenna has vetted the whole 
## script and confirmed that it functions as intended with the data provided. 

# Clear environment and start fresh
rm(list = ls())

# Set working directory
setwd("C:/Users/jenna/OneDrive - York University/PhD Program/Data/Hybridization data/Entropy/Fixed_MPGL/big_Q")

# Load required packages
library(rhdf5)
library(abind)
library(dplyr)

# Load updated metadata with genetic species assignments
# This file contains sample names, genotypic species assignments (based on q-score confidence interval overlaps),
# and metadata like sample site and zone etc. 
# 
admixture <- read.csv("fescue_reassignment_CIs.csv")

admixture


# Load sample names 
# This is just a list of sample names in the order that they show up in the ENTROPY output files 
# (so that the Q scores line up properly)

names <- read.table("sample_names_list.txt", header = F)

names

# Read in HDF5 files
# These are ENTROPY output files for big and little q (admixture & intersource ancestry)
# There's 3 files because you do three replicates for accuracy (each hdf5 file contains the big and little q scores)

q1 <- h5read("fescue_bigQ_k2_rep1.hdf5", "q")
q2 <- h5read("fescue_bigQ_k2_rep2.hdf5", "q")
q3 <- h5read("fescue_bigQ_k2_rep3.hdf5", "q")

Q1 <- h5read("fescue_bigQ_k2_rep1.hdf5", "Q")
Q2 <- h5read("fescue_bigQ_k2_rep2.hdf5", "Q")
Q3 <- h5read("fescue_bigQ_k2_rep3.hdf5", "Q")

# Calculate means across replicates
allq <- abind(q1, q2, q3, along=1) 
q <- t(apply(allq, 2:3, mean))

allQ <- abind(Q1, Q2, Q3, along=1)
Q <- t(apply(allQ, 2:3, mean))

# Calculate confidence intervals
qci <- apply(allq, 2:3, quantile, probs=c(0.025, 0.975))
Qci <- apply(allQ, 2:3, quantile, probs=c(0.025, 0.975))

# Merge sample names with updated metadata
triangle_data <- merge(names, admixture, by.x = "V1", by.y = "sample_name")

### Basic triangle plot

# Define colours for the plot
# Use the same colors as the faceted plots
species_colors_alpha <- c("FC" = "#779be7A0", "FH" = "#f56476A0", "HY" = "#5bc8afA0")
species_colors_solid <- c("FC" = "#779be7", "FH" = "#f56476", "HY" = "#5bc8af")


## Single triangle plot

plot(q[,2], Q[,2], xlim=c(0,1), ylim=c(0,1), axes = F, 
     xlab="Admixture Proportion (q)", ylab="Inter-source Ancestry (Q)", 
     type="n",
     cex.lab = 1.3) 

axis(1, at = seq(0, 1, by = 0.25))
axis(2, at = seq(0, 1, by = 0.25))
box()

# Triangle marks 
arrows(0,0,0.5,1, length=0, col="gray45", lwd=3)  # Changed from "gray" to "gray45" and added lwd=3
arrows(0.5,1,1,0, length=0, col="gray45", lwd=3)  # Changed from "gray" to "gray45" and added lwd=3

# Confidence intervals 
arrows(qci[2,2,], Q[,2], qci[1,2,], Q[,2], length=0, col="gray45", lwd=1) # Added lwd=1
arrows(q[,2], Qci[1,2,], q[,2], Qci[2,2,], length=0, col="gray45", lwd=1) # Added lwd=1

# Plot points with SEMI-TRANSPARENT colors and larger size
points(q[,2], Q[,2], 
       col = species_colors_alpha[triangle_data$ci_geno_species], 
       pch = 19, cex = 2.5)  # Changed from pch=16 to pch=19, cex=1.5 to cex=2.5

# Add legend 
legend("topright", 
       legend = c("F. campestris", "F. hallii", "Hybrid"), 
       fill = species_colors_solid,  # Use the solid colors that match
       cex = 1.0,
       title = "Species",
       title.cex = 1.2,  # Makes the title bigger than the legend text
       title.font = 2)

print("Basic triangle plot completed!")



## Faceted Plots ---------------------------------------------------------------

# This facets the plots according to the zone the samples were collected from. 
# A bunch of this code is just annoying workarounds to make the plots look the way I want them to. The code 
# works well for the fescue samples but it has been supplemented by Seqera, just as a heads up!  


library(gridExtra)
library(grid)

# Replace zone names FIRST, before creating the zones variable
triangle_data$Zone[triangle_data$Zone == "Sympatric 1"] <- "Sympatric (Glenbow Ranch)"
triangle_data$Zone[triangle_data$Zone == "Sympatric 2"] <- "Sympatric (Cypress Hills)"

# NOW create the zones variable with the updated names
zones <- unique(triangle_data$Zone)
zones <- zones[!is.na(zones)]

# Create semi-transparent colors
species_colors_alpha <- c("FC" = "#779be7A0", "FH" = "#f56476A0", "HY" = "#5bc8afA0")
species_colors_solid <- c("FC" = "#779be7", "FH" = "#f56476", "HY" = "#5bc8af")

# Function to create individual triangle plot
create_triangle_plot <- function(zone) {
  zone_indices <- which(triangle_data$Zone == zone)
  
  # Set up individual plot with normal margins
  par(mar = c(5, 5, 3, 2))
  
  plot(q[zone_indices, 2], Q[zone_indices, 2], 
       xlim = c(0, 1), ylim = c(0, 1), axes = F,
       xlab = "Admixture Proportion (q)", 
       ylab = "Inter-source Ancestry (Q)",
       main = zone, cex.main = 1.2,
       cex.lab = 1.5)
  
  axis(1, at = seq(0, 1, by = 0.25))
  axis(2, at = seq(0, 1, by = 0.25))
  box()
  
  # Triangle marks
  arrows(0, 0, 0.5, 1, length = 0, col = "gray", lwd = 3)
  arrows(0.5, 1, 1, 0, length = 0, col = "gray", lwd = 3)
  
  # Confidence intervals
  arrows(qci[2, 2, zone_indices], Q[zone_indices, 2], 
         qci[1, 2, zone_indices], Q[zone_indices, 2], 
         length = 0, col = "gray45", lwd = 1)
  arrows(q[zone_indices, 2], Qci[1, 2, zone_indices], 
         q[zone_indices, 2], Qci[2, 2, zone_indices], 
         length = 0, col = "gray45", lwd = 1)
  
  # Plot points
  points(q[zone_indices, 2], Q[zone_indices, 2],
         col = species_colors_alpha[triangle_data$ci_geno_species[zone_indices]],
         pch = 19, cex = 2.5)
  
  # Add sample size
  mtext(paste0("n = ", length(zone_indices)), side = 3, line = 0, cex = 0.8)
}

# Set figure margins BEFORE creating layout - INCREASED BOTTOM MARGIN
par(mai = c(1.4, 1, 0.6, 0.6))  # Increased bottom margin from 1 to 1.4

# Create a layout with plots and legend
# Set up the layout: 2x2 grid for plots + space for legend
layout_matrix <- rbind(c(1, 2, 5),
                       c(3, 4, 5))

layout(layout_matrix, widths = c(1, 1, 0.4))

# Create each plot
for(i in 1:length(zones)) {
  create_triangle_plot(zones[i])
}

# Create the shared legend in the 5th panel
par(mar = c(0, 0, 0, 0))
plot.new()
legend("center", 
       legend = c("F. campestris", "F. hallii", "Hybrid"), 
       fill = species_colors_solid,
       cex = 1.2, 
       title = "Species",
       title.cex = 1.3)

# Reset layout AND margins
layout(1)
par(mar = c(5, 5, 4, 2) + 0.1, mai = c(1.02, 1, 0.82, 0.42))

print("Combined triangle plots with shared legend completed!")
getwd()



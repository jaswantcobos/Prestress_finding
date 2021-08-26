# ------------------------------------------------------------------------------
# Integral Feasible Prestress of Tensegrity Structures 
#
# Date: 27-07-2021
# Author:Jaswant Cobos, Marlon E. Cobos
# Contact: Jaswant Cobos (jaswant.cobos@gmail.com)
# ------------------------------------------------------------------------------

# Code to determine the Integral Feasible Prestress of Tensegrity Structures
# using the Doble Singular Value Decomposition (DSVD) method developed by 
# X.F.Yuan and S.L.Dong (https://doi.org/10.1016/S0045-7949(03)00254-2) 

# ------------------------------------------------------------------------------

# load functions and package
source("R_code/prestress_finding_functions.R")

#set working directory
#setwd("C:/Users") #only an example

# ------------------------------------------------------------------------------

# Data

# This code read the information from the csv documents (connectivity.csv,
# coordinates.csv, and free_nodes.csv)

# The user should insert data in the documents this way:

# Connectivity matrix:
# |Element|Initial_Node|Final_Node|Symmetry_Group| 
# Symmetry groups must be ordered from least to greatest

# Coordinate matrix:
# |Node|x_Coordinates|y_Coordinates|y_Coordinates|

# Free nodes vector:
# Insert data in the first column of this document
con <- read.csv("R_code/connectivity.csv") # connectivity matrix
coor <- read.csv("R_code/coordinates.csv") # matrix with coordinates
fn <- read.csv("R_code/free_nodes.csv") # free nodes of the structure

# ------------------------------------------------------------------------------

# calculations
p_stress <- find_prestress(coordinates = coor, connectivity = con, 
                            free_nodes = fn)

# ------------------------------------------------------------------------------

## save analysis results
save(p_stress, file = "R_code/Prestress_complete_results.RData")

write.csv(p_stress$presstress_matrix, "R_code/prestress_results.csv",
          row.names = FALSE)

# ------------------------------------------------------------------------------

# plotting
plot_prestress(p_stress)


#writeWebGL(filename = "plot_prestress.html")

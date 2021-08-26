# ------------------------------------------------------------------------------
# R package and functions to calculate and plot Integral Feasible Prestress for 
# Tensegrity Structures
#
# Date: 27-07-2021
# Author:Jaswant Cobos, Marlon E. Cobos
# ------------------------------------------------------------------------------

# r package for plotting
if (!require(rgl)) {
  install.packages("rgl")
  library(rgl)
}

# ------------------------------------------------------------------------------

# useful for 3d plots in UNIX like OS
if (.Platform$OS.type == "unix") {
  options(rgl.printRglwidget = TRUE)
}

# ------------------------------------------------------------------------------

# function for prestress calculations
find_prestress <- function(coordinates, connectivity, free_nodes) {
  
  # dimensions
  b <- nrow(connectivity)
  
  # connectivity matrix
  cs <- matrix(0, nrow = b, ncol = nrow(free_nodes))
  
  for (i in 1:b) {
    cs[i, min(connectivity[i, 2:3])] <- 1
    cs[i, max(connectivity[i, 2:3])] <- -1
  }
  
  # equilibrium matrix
  ## element length
  l <- list()
  
  for (i in 1:b) {
    ### finding coordinates
    ci <- coordinates[coordinates[, 1] == connectivity[i, 2], 2:4]
    cf <- coordinates[coordinates[, 1] == connectivity[i, 3], 2:4]
    
    ### finding lengths
    l[[i]] <- sqrt((cf[1] - ci[1])^2 + (cf[2] - ci[2])^2 + (cf[3] - ci[3])^2)
  }
  
  l <- unlist(l)
  
  l <- diag(l, nrow = b, ncol = b)
  
  ## equilibrium matrix calculation
  il <- solve(l)
  tcs <- t(cs)
  
  A <- t(cbind(t(tcs %*% diag(c(cs %*% coordinates[, 2]), b, b) %*% il),
               t(tcs %*% diag(c(cs %*% coordinates[, 3]), b, b) %*% il),
               t(tcs %*% diag(c(cs %*% coordinates[, 4]), b, b) %*% il)))
  
  # pre-stress independent mode finding
  r1 <- qr(A)$rank
  s1 <- b - r1
  v <- svd(A, nu = nrow(A), nv = ncol(A))$v
  mpi <- as.matrix(v[, (r1 + 1):b]) # check if it is V or what?
  
  # pre-stress integral mode finding
  a <- nrow(mpi)
  n <- max(connectivity[, 4])
  
  ## e vectors from the method
  E <- matrix(0, nrow = a, ncol = n)
  
  ### ensemble of e vectors
  p <- 1
  q <- 1
  
  for (h in 1:n) {
    while ((connectivity[p, 4] - h) == 0) {
      E[p, q] <- -1
      p <- p + 1
      
      if (p == (a + 1)) {
        p <- 1
      }
    }
    q <- q + 1
  }
  
  
  # matrix of method
  X <- cbind(mpi, E)
  v2 <- svd(X, nu = nrow(X), nv = ncol(X))$v
  r2 <- qr(X)$rank
  b2 <- ncol(X)
  s2 <- s1 + n - r2
  
  mpi2 <- as.matrix(v2[, (r2 + 1):b2])
  w1 <- E %*% mpi2[(ncol(mpi) + 1):b2, ]
  
  ## ordering coefficients according to symmetry group
  W <- as.matrix(w1[!duplicated(w1[, 1]), ])
  W <- W / max(abs(c(W)))
  
  # results
  return(list(coordinates = coordinates, connectivity = connectivity, 
              free_nodes = free_nodes, presstress_matrix = W))
}

# ------------------------------------------------------------------------------

# function for plotting 
plot_prestress <- function(preestress_list, 
                           cols = c("#D7D7D7", "#8E8E8E", "#000000")) {
  
  # plotting
  for (i in 1:nrow(preestress_list[[2]])) {
    ## coordinates
    x <- c(preestress_list[[1]][preestress_list[[2]][i, 2], 2], 
           preestress_list[[1]][preestress_list[[2]][i, 3], 2])
    y <- c(preestress_list[[1]][preestress_list[[2]][i, 2], 3], 
           preestress_list[[1]][preestress_list[[2]][i, 3], 3])
    z <- c(preestress_list[[1]][preestress_list[[2]][i, 2], 4], 
           preestress_list[[1]][preestress_list[[2]][i, 3], 4])
    
    ## color and width coding
    el <- preestress_list[[4]][preestress_list[[2]][i, 4]]
    
    if (el == 0) {
      coli <- cols[1]
      lwdi <- 0.5
    } else {
      coli <- ifelse(el > 0, cols[3], cols[2]) 
      
      lwdi <- round(abs(el) * 10, 2)
      lwdi <- ifelse(lwdi == 0, 0.8, lwdi)
    } 
    
    ## plotting
    if (i == 1) {
      plot3d(x, y, z, xlab = "", ylab = "", zlab = "", type = "l", 
             col = coli, lwd = lwdi, aspect = "iso", axes = FALSE) 
    } else {
      plot3d(x, y, z, type = "l", aspect = "iso",
             col = coli, lwd = lwdi, add = TRUE)
    }
  }
  
  if (.Platform$OS.type == "unix") {
    bg3d("white")
  }
}

# ------------------------------------------------------------------------------
# This R script takes in a csv file containing xyz coordinates of spots and return a 3D contour field by kernel density of dense regions
# code is adoped from https://stackoverflow.com/questions/60001481/create-contour-in-a-3d-kernel-density-and-find-which-points-are-within-that-co

library(MASS)
library(misc3d)
library(rgl)
library(sp)
library(oce)

mv = read.csv("filepath")

# Create kernel density
dens3d <- kde3d(mv[,1], mv[,2], mv[,3], n = 40)
dens3d$d
dens3d$x
# Find the estimated density at each observed point
datadensity <- approx3d(dens3d$x, dens3d$y, dens3d$z, dens3d$d, 
                        mv[,1], mv[,2], mv[,3])
# Find the contours
prob <-  .5
levels <- quantile(datadensity, probs = prob, na.rm = TRUE)
levels
# Plot it
colours <- c("gray", "orange")
cuts <- cut(datadensity, c(0, levels, Inf))
cuts
for (i in seq_along(levels(cuts))) {
  gp <- as.numeric(cuts) == i
  spheres3d(mv[gp,1], mv[gp,2], mv[gp,3], col = colours[i], radius = 5)
}
box3d(col = "gray")
contour3d(dens3d$d, level = levels, x = dens3d$x, y = dens3d$y, z = dens3d$z, #exp(-12)
          alpha = 0.1, color = "gray", color2 = "gray", add = TRUE)
rgl.snapshot( "filepath", fmt = "png", top = TRUE )

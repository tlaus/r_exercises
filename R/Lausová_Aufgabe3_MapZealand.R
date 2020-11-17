## Aufgabe 3 - New Zealand map
## Tereza Lausová
## 14.9.2020
## Code from https://www.stat.auckland.ac.nz/~paul/RGraphics/examples-arden.R

# Load libraries
# install.packages("mapproj")
library(maps)

# Plot
pdf(file = here("plots/Lausová_Aufgabe3_MapZealand.pdf"))
par(mar=rep(0, 4))
map("nz", fill=TRUE, col="cornflowerblue")
points(174.75, -36.87, pch=16, cex=2, add = TRUE)
arrows(172, -36.87, 174, -36.87, lwd=3)
text(172, -36.87, "Heidelberg", adj=1, cex=2)
# mini world map as guide
maplocs <- map(projection="sp_mercator", wrap=TRUE, lwd=0.1, 
               col="grey", ylim=c(-60, 75),
               interior=FALSE, orientation=c(90, 180, 0), add=TRUE,
               plot=FALSE)
xrange <- range(maplocs$x, na.rm=TRUE)
yrange <- range(maplocs$y, na.rm=TRUE)
aspect <- abs(diff(yrange))/abs(diff(xrange))
# customised to 6.5 by 4.5 figure size
par(fig=c(0.99 - 0.5, 0.99, 0.01, 0.01 + 0.5*aspect*4.5/6.5), 
    mar=rep(0, 4), new=TRUE)
plot.new()
plot.window(xlim=xrange,
            ylim=yrange)
map(projection="sp_mercator", wrap=TRUE, lwd=0.1, ylim=c(-60, 75),
    interior=FALSE, orientation=c(90, 180, 0), add=TRUE)
symbols(-.13, -0.8, circles=1, inches=0.1, add=TRUE)

dev.off()
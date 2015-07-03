library(compositions)
data(ArcticLake)
ilrDefinition <- function(sbp, side="+-") {
  
  if (nrow(sbp) != (ncol(sbp)-1)) stop("SBP not valid")
  
  ilrDef <- vector()
  for (n in 1:nrow(sbp)) {
    pos <- names(sbp[n,][which(sbp[n,] == 1)])
    neg <- names(sbp[n,][which(sbp[n,] == -1)])
    if (side=="-+") {
      pos <- rev(pos)
      neg <- rev(neg)
    }
    pos.group <- character()
    neg.group <- character()
    for (i in 1:length(pos)) {
      if (i == 1) {
        pos.group <- paste(pos.group,pos[i], sep="")
      } else {
        pos.group <- paste(pos.group,pos[i], sep=",")
      }
    }
    for (i in 1:length(neg)) {
      if (i == 1) {
        neg.group <- paste(neg.group,neg[i], sep="")
      } else {
        neg.group <- paste(neg.group,neg[i], sep=",")
      }
    }
    if (side=="+-") {
      ilrDef[n] <- paste("[",pos.group," | ",neg.group,"]", sep="")
    } else if (side=="-+") {
      ilrDef[n] <- paste("[",neg.group," | ",pos.group,"]", sep="")
    }
    
  }
  
  ilrDef
}
rotate <- function(xy, deg, center = c(0, 0)) {
  xy <- unclass(xy)
  alpha <- -deg*pi/180 # rotation angle
  rotm <- matrix(c(cos(alpha),sin(alpha),-sin(alpha),cos(alpha)),ncol=2) #rotation matrix
  rotated <- t(rotm %*% (t(xy) - center) + center)
  return(rotated)
}
circle <- function(rad = 1, center = c(0, 0), npoints = 100) {
  seqc <- seq(0, 2*pi, length=npoints)
  circ <- t(rbind(center[1] + sin(seqc)*rad, center[2] + cos(seqc)*rad)) 
  return(circ)
}
arrowBal <- function(labels, deg = 0, rad = 1, center = c(0,0), pos = c(1, 1), color = "black") {
  x0 <- center[1]
  x1 <- cos(deg*pi/180)*rad + center[1]
  x2 <- cos((deg+90)*pi/180)*rad + center[1]
  y0 <- center[2]
  y1 <- sin(deg*pi/180)*rad + center[2]
  y2 <- sin((deg+90)*pi/180)*rad + center[2]
  arrows(x0, y0, x1, y1, col = color)
  text(x1, y1, labels = labels[1], pos = pos[1], col = color)
  arrows(x0, y0, x2, y2, col = color)
  text(x2, y2, labels = labels[2], pos = pos[2], col = color)
}


parts <- ArcticLake[, 1:3]
comp <- acomp(parts)

sbp1 <- matrix(c( 1,-1,-1,
                  0, 1,-1),
               byrow = TRUE,
               ncol = 3)
sbp2 <- matrix(c(-1,-1, 1,
                  1,-1, 0),
               byrow = TRUE,
               ncol = 3)
sbp3 <- matrix(c( 1,-1, 1,
                  1, 0,-1),
               byrow = TRUE,
               ncol = 3)
colnames(sbp1) <- colnames(parts)
colnames(sbp2) <- colnames(parts)
colnames(sbp3) <- colnames(parts)

bal1 <- ilr(comp, V = gsi.buildilrBase(t(sbp1)))
bal2 <- ilr(comp, V = gsi.buildilrBase(t(sbp2)))
bal3 <- ilr(comp, V = gsi.buildilrBase(t(sbp3)))
colnames(bal1) <- ilrDefinition(sbp1, side = "-+")
colnames(bal2) <- ilrDefinition(sbp2, side = "-+")
colnames(bal3) <- ilrDefinition(sbp3, side = "-+")

par(mar = c(4,4,1,1), mfrow = c(1, 2))
plot(unclass(bal1), xlim = c(-4, 4), ylim = c(-4, 4), col = "black", cex=0.5,
     xlab = "ilr 1", ylab = "ilr 2")
abline(v=0, col = "grey70")
abline(h=0, col = "grey70")
points(unclass(bal2), col = "red", cex = 0.75)
points(unclass(bal3), col = "blue", cex = 1)

plot(unclass(bal1), xlim = c(-4, 4), ylim = c(-4, 4), col = "black", cex = 0.5) 
points(rotate(bal2, deg = 120), col = "red", cex = 0.75) # rotate 240 deg
points(rotate(bal3, deg = 60), col = "blue", cex = 1) # rotate 300 deg
rad = 2
lines(circle(rad = rad))
arrowBal(labels = colnames(bal2), deg = 0, rad = rad, pos = c(4,3), color = "black")
arrowBal(labels = colnames(bal2), deg = 120, rad = rad, pos = c(2,2), color = "red")
arrowBal(labels = colnames(bal3), deg = 60, rad = rad, pos = c(4,2), color = "blue")

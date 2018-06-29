#' class for ddgv
#'
#' @description  Instan
#' @export
#' @useDynLib moddicomV2
#' @import misc3d



morphologicalFeatures <- function(imgObj,px,py,pz){

#dyn.load("/home/kbocalc/Desktop/J/DresdenFeatures/lateralSurface.so");


  nVoxel <- dim(imgObj)[1]*dim(imgObj)[2]*dim(imgObj)[3] - sum(is.na(imgObj))
  n <- numeric()
  n <- table(imgObj)
  p <- n/nVoxel

  #Initialise data table for storing Morphological features
  featNames <- c("F_morph.surface", "F_morph.volume", "F_morph.av",
                 "F_morph.comp.1", "F_morph.comp.2", "F_morph.sph.dispr",
                 "F_morph.sphericity", "F_morph.asphericity","F_morph.com", "F_morph.diam",
                 "L_major","L_minor","L_least","F_morph.pca.elongation","F_morph.pca.flatness")
  F_morph <- data.frame(matrix(NA, ncol=length(featNames)))
  colnames(F_morph) <- featNames


  #Generate roi mask consisting of 0s and 1s
  roiObj    <- is.finite(imgObj) * 1

  #Pad roi mask with an additional boundary of 0s - otherwise the marching cubes algorithm in the misc3d library will not create an appropriate mesh.
  roiObj.dim <- dim(roiObj)
  roiObj.pad <- array(data=0, dim=roiObj.dim+2)

  #Place the original roi mask into the centre of the padded roi
  roiObj.pad[2:(roiObj.dim[1] + 1), 2:(roiObj.dim[2] + 1), 2:(roiObj.dim[3] + 1) ] <- roiObj

###
objS <- services()
ppx <-seq(0,px*(dim(roiObj.pad)[1]-1),by=px)
ppy <-seq(0,py*(dim(roiObj.pad)[2]-1),by=py)
ppz <-seq(0,pz*(dim(roiObj.pad)[3]-1),by=pz)

mesh.triangle <- contour3d(roiObj.pad,level=0.5, x = ppx,y = ppy, z = ppz, engine = "none")
mesh<-objS$triangle2mesh(x = mesh.triangle)
F_morph$F_morph.surface <- StructureSurface(mesh = mesh,measure.unit = "mm2")

#F_morph.surface <- vcgArea(mesh)

###
voxelCoords <-which(!is.na(imgObj),arr.ind = TRUE)
voxelCoordsLenghts <- matrix(nrow=dim(voxelCoords)[1],ncol=3)
voxelCoordsLenghts[,1] <- px * voxelCoords[,1]
voxelCoordsLenghts[,2] <- py * voxelCoords[,2]
voxelCoordsLenghts[,3] <- pz * voxelCoords[,3]
###

F_morph$F_morph.volume <- StructureVolume(mesh,measure.unit = "mm3")


F_morph$F_morph.av <- F_morph$F_morph.surface / F_morph$F_morph.volume
F_morph$F_morph.comp.1 <-F_morph$F_morph.volume / (sqrt(pi) * (sqrt(F_morph$F_morph.surface)^3))
F_morph$F_morph.comp.2 <-  36 * pi * F_morph$F_morph.volume^2 / F_morph$F_morph.surface^3
F_morph$F_morph.sph.dispr <- F_morph$F_morph.surface / (36*pi*F_morph$F_morph.volume^2)^(1/3)
F_morph$F_morph.sphericity <- (36*pi*F_morph$F_morph.volume^2)^(1/3) / F_morph$F_morph.surface
F_morph$F_morph.asphericity <- (F_morph$F_morph.surface^3 / (36*pi*F_morph$F_morph.volume^2))^(1/3) - 1


## stuff to calculate CENTER OF MASS SHIFT
CoMgeom <- numeric()
CoMgeom[1] <- sum(voxelCoords[,1])/nVoxel
CoMgeom[2] <- sum(voxelCoords[,2])/nVoxel
CoMgeom[3] <- sum(voxelCoords[,3])/nVoxel




grayLevels <- numeric()
for (i in seq(1,dim(voxelCoords)[1])){
  grayLevels[i] <- imgObj[voxelCoords[i,1],voxelCoords[i,2],voxelCoords[i,3]]
}

voxelCoordsGl <- cbind(voxelCoords,grayLevels)
voxelWeightedCoords <- voxelCoordsGl[,1:3] * voxelCoordsGl[,4]


CoMgl <- numeric()
CoMgl[1] <- sum(voxelWeightedCoords[,1])/sum(voxelCoordsGl[,4])
CoMgl[2] <- sum(voxelWeightedCoords[,2])/sum(voxelCoordsGl[,4])
CoMgl[3] <- sum(voxelWeightedCoords[,3])/sum(voxelCoordsGl[,4])
###



F_morph$F_morph.com <- sqrt((px*(CoMgeom[1]-CoMgl[1]))^2 + (py*(CoMgeom[2]-CoMgl[2]))^2 + (pz*(CoMgeom[3]-CoMgl[3]))^2)



edgPoints <- voxelCoords[chull(voxelCoords),]
distanceMatrixHull <- dist(edgPoints)
F_morph$F_morph.diam <- max(distanceMatrixHull)


#Principal Compnent Analysis
imgObj.data.frame <- which(x = !is.na(imgObj), arr.ind = T)
imgObj.data.frame_pca <- c(px,py,pz)*imgObj.data.frame
pca <- prcomp(imgObj.data.frame_pca,
                 center = FALSE,
                 scale. = FALSE)


pca$sdev <- sort(pca$sdev,decreasing=T)

F_morph$L_major <- 4*sqrt(pca$sdev[1])
F_morph$L_minor <- 4*sqrt(pca$sdev[2])
F_morph$L_least <- 4*sqrt(pca$sdev[3])

F_morph$F_morph.pca.elongation <- sqrt(F_morph$L_minor/F_morph$L_major)
F_morph$F_morph.pca.flatness <- sqrt(F_morph$L_least/F_morph$L_major)

F_morph <- as.numeric(F_morph)
names(F_morph) <- featNames
return(F_morph)

}

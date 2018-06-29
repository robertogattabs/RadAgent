#' class for viewing
#' 
#' @description  allows to view what's happening
#' @useDynLib moddicomV2
#' @import Rvcg rgl 
#' @export
viewer<-function() {
  plotROIs<-function(obj , arrayROINames = NA , sliceNumber = 1, ps.x = NA, ps.y = NA, ps.z = NA,paletteColor = c()) {
    
    
    # controlli formali
    if( length(arrayROINames)>20 )  stop("Ahahahahahahah!!!! Ciccio, al massimo 20 ROI, forse meno.....");
    # definisce una palette di circa 20 color (bo'... son 20?)
    if(length(paletteColor)==0)
      paletteColor<-c('#38D747','#F366FD','#C9360E','#ADCDE7','#391C1F','#F3BC6E','#058557','#5E601A','#FCACA9','#F7556B','#E7EDBE','#759D16','#551E50','#A76A19','#E994F2','#2DC3CB','#B3EE5A','#FB8F5B','#40B85C','#04C1A9')
    # verifica che ce ne sia almeno una
    if(length(arrayROINames)==0) return;
    align<-list();  cubeDim<-c();
    # browser()
    # dataStorage <- obj$getAttribute(attribute = "dataStorage")
    # initial.point <- objService$get3DPosFromNxNy(Nx = 1,Ny = 1,oM = dataStorage$info[[SeriesInstanceUID]][[1]]$orientationMatrix)
    # final.point <- objService$get3DPosFromNxNy(Nx = 100,Ny = 100,oM = dataStorage$info[[SeriesInstanceUID]][[numberOfSlices]]$orientationMatrix)
    # directions <- final.point - initial.point
    
    
    # frulla sull'array di ROINames per allinearle una per una
    for( ROIName in seq(1,length(arrayROINames))) {
      tmpAlign <- obj$getAlignedStructureAndVoxelCube(ROIName = arrayROINames[ROIName] , ps.x = ps.x, ps.y = ps.y,ps.z = ps.z)
      # prendi i punti ROI allineati
      align[[ arrayROINames[ ROIName] ]]<-tmpAlign$ROI
      # prendi il voxelCube (anche se riscrive non è importante, tanto è sempre lo stesso)
      vc<-tmpAlign$voxelCube
      # prendi le dimensioni del voxelCube (sovrascrive? No problem, vedi sopra)
      #cubeDim<-tmpAlign$cubeDim
      cubeDim<-dim(tmpAlign$voxelCube)
    }
    
    # plotta l'immagine  
    image(vc[,,sliceNumber],col = grey.colors(255,start=0),axes = FALSE)
    # prendi in numero di slice
    whichROILine<-sliceNumber;
    colorIncrement<-1
    
    # Cicla per plottare ogni ROI nell'array delle ROI
    for( ROIName in seq(1,length(arrayROINames))) {
      # scorri nella lista delle ROIPointList per cercare la fetta corretta
      for( internalROIPlanes in seq(1,length( align[[   arrayROINames[ ROIName]    ]] ))    ) {
        # 
        if(align[[   arrayROINames[ ROIName]    ]][[internalROIPlanes]][[1]][1,3] == whichROILine  ) {
          
          for(subROIOnTheSamePlane in seq(1, length(align[[   arrayROINames[ ROIName]    ]][[internalROIPlanes]])    )){
            ROI<-align[[   arrayROINames[ ROIName]    ]][[internalROIPlanes]][[subROIOnTheSamePlane]]
            
            ROI[,1]<-ROI[,1] / cubeDim[1]
            ROI[,2]<-1-ROI[,2] / cubeDim[2]
            points(x = ROI[,1],y = ROI[,2],type='l', lwd = 3, col = paletteColor[colorIncrement])  
          }
          colorIncrement<-colorIncrement+1;
        }
      }
    }
  }
  return(
    list( 
      "plotROIs"=plotROIs 
    )
  )
}
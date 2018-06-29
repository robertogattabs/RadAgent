getEmptyStructure<-function( objGLT , mask ) {
  objS<-services();

  res <-list("DOM"=list(),"final.array"=list(),"masked.images"=list(),"geometricalInformationOfImages"=list(),"resamplingInformation"=list())

  res$final.array <- objGLT$getImageVoxelCube();
  res$masked.images <- mask * res$final.array;
  
  res$geometricalInformationOfImages<-list()
  res$geometricalInformationOfImages$pixelSpacing <- objGLT$getPixelSpacing()[1:2]
  res$geometricalInformationOfImages$Rows <- dim(res$final.array)[1]
  res$geometricalInformationOfImages$Columns <- dim(res$final.array)[2]
  res$geometricalInformationOfImages$supposedNumberOfSlices <- dim(res$final.array)[3]
  res$resamplingInformation<-list()
  res$resamplingInformation$resampled<-FALSE
  
  croppedRes<-list()
  croppedRes$DOM<-res$DOM
  croppedRes$geometricalInformationOfImages<-res$geometricalInformationOfImages

  croppedRes$masked.images<-objS$cropCube( bigCube = res$masked.images)

  croppedRes$masked.images$location$fe<-dim(res$masked.images)[1]
  croppedRes$masked.images$location$se<-dim(res$masked.images)[2]
  croppedRes$masked.images$location$te<-dim(res$masked.images)[3]
  croppedRes$geometricalInformationOfImages$koc<-"littleCube"
  croppedRes$resamplingInformation<-res$resamplingInformation
  class(croppedRes)<-"geoLetStructureVoxelList"
  
  return( croppedRes )
}
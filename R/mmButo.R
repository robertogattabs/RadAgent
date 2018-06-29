#' class for loading multiple sets of DICOM studies from filesystem
#' 
#' @description  Instantiate an object of the class \code{new.mmButo}. This represents just the classname,
#'               for each instantiated object many methods are available not i S3 or S4 but by closures method.
#'               The available methods are:
#' @export
mmButo<-function(
                  folderCleanUp = FALSE,
                  loadXMLInCache = FALSE,
                  loadRAWInCache = TRUE) {
  list_geoLet<-list()
  attr_folderCleanUp<-FALSE                              # force to re-dump DICOM files 
  attr_ROIVoxelMemoryCache<-FALSE                        # force to cache ROI Voxel
  attr_loadXMLInCache<-loadXMLInCache
  attr_loadRAWInCache<-loadRAWInCache  
  # ========================================================================================
  # loadCollection: load a set of subfolders
  # ========================================================================================
  loadCollection<-function( Path, collectionID = "default") {
    objService<-services();
    listaFolders<-list.dirs( Path )
    ct<-1
    if( length(list_geoLet[[collectionID]]) == 0 ) list_geoLet[[collectionID]]<<-list()
    for( folderName in listaFolders[2:length(listaFolders)] ) {
      
      if ( Sys.info()["sysname"] == "Windows") {
        folderNameFS<-chartr("\\","/",folderName)
      }
      else folderNameFS<-folderName;
      
      list_geoLet[[collectionID]][[ folderNameFS ]]<<-geoLet( ROIVoxelMemoryCache = TRUE , 
                                                              folderCleanUp = attr_folderCleanUp,
                                                              loadXMLInCache = attr_loadXMLInCache,
                                                              loadRAWInCache = attr_loadRAWInCache   )
      list_geoLet[[collectionID]][[ folderNameFS ]]$openDICOMFolder( folderNameFS )
    }
  }  
  # ========================================================================================
  # extractROIs: extract one or more ROIs voxels
  # ======================================================================================== 
  getROIVoxel<-function(  ROIName,  new.pixelSpacing=c() , collectionID = "default") {
    objS<-services();
    singleROI<-ROIName
    print("=================================================================");
    print( paste( c("getROIVoxel for ROI: ",ROIName)   , collapse='') );
    print("=================================================================");
    list_extractROIVoxel<-list();
    list_extractROIVoxel[[collectionID]]<-list();    
    
    for(folderName in names(list_geoLet[[collectionID]])) {
      list_extractROIVoxel[[collectionID]][[ folderName ]]<-list();
      
      print( paste( c("Now processing=",folderName)   , collapse='') );
      
      a <- list_geoLet[[collectionID]][[folderName]]$getROIVoxels(Structure = singleROI, new.pixelSpacing = new.pixelSpacing );
      
      invisible(gc())
      
      if((TRUE %in% is.na(a)) == FALSE  ) {
        list_extractROIVoxel[[collectionID]][[ folderName ]][[ singleROI ]]<-a;
      }
      else {
        list_extractROIVoxel[[collectionID]][[ folderName ]][[ singleROI ]]<-NA
      }
      
    }
    arr2Return<-list();
    for( folderName in names(list_extractROIVoxel[[collectionID]])) {
      arr2Return[[folderName]]<-list_extractROIVoxel[[collectionID]][[folderName]][[singleROI]]
    }
    class(arr2Return)<-"mmButoStructureVoxelList"
    invisible(gc())
    return(arr2Return)
  }  
  # ========================================================================================
  # getImageVoxel: give back the entire image voxel cubes
  # ========================================================================================   
  getImageVoxel<-function( ps.x = NA, ps.y = NA, ps.z = NA  ) {
    collectionID<-"default"
    imgVC<-list();
    
    for(i in names(list_geoLet[[collectionID]])) {
      imgVC[[i]]<-list_geoLet[[collectionID]][[i]]$getImageVoxelCube(ps.x = ps.x, ps.y = ps.y, ps.z = ps.z)
    }
    return(imgVC)
  }
  # ========================================================================================
  # getAttribute: the usual 'getAttribute' method
  # ======================================================================================== 
  getAttribute<-function( attributeName ) {
    if( attributeName == 'list_geoLet' ) return( list_geoLet );
    logObj$handle( "error" , "The attribute does not exist"  );
  }  
  # ========================================================================================
  # getAllPixelSpacing: return the pixelSpacings of all the loaded geoLet object
  # ========================================================================================   
  getAllPixelSpacing<-function( collectionID='default' ) {
    elArr<-list();
    for(patName in names(list_geoLet[[ collectionID ]]) ) {
      elArr[[patName]]<-list_geoLet[[ collectionID ]][[patName]]$getPixelSpacing();
    }
    return(elArr);    
  }  
  # ========================================================================================
  # mmButoLittleCube.expand: expand a cropped ROI in order to satisfy compatibility
  # with older releases of moddicom
  # ========================================================================================   
  expandCube<-function( ROIVoxelElement ) {
    # get the needed parameters
    pc<-ROIVoxelElement
    x<-pc$masked.images$location$min.x; y<-pc$masked.images$location$min.y; z<-pc$masked.images$location$min.z
    fe<-pc$masked.images$location$fe; se<-pc$masked.images$location$se;  te<-pc$masked.images$location$te
    # invocke the procedure from class Services. The procedure is in Services because it could also be used 
    # for different issues, in perspective, from classes different than mmButo.
    objS<-services()

    bigVoxelCube<-objS$expandCube(littleCube = pc$masked.images$voxelCube, x.start = x, y.start=y, z.start=z, fe = fe, se = se, te = te )    
    return(bigVoxelCube)
  }  
  class<-function() {
    return("mmButo");
  }  
  # ========================================================================================
  # getROIVoxelStats
  # ========================================================================================
  getROIVoxelStats<-function( ROIVoxelList ) {
    # define the empty arrays
    dataInfo<-list()
    min.arr<-c(); max.arr<-c(); mean.arr<-c(); sd.arr<-c(); median.arr<-c()
    # loop in order to calcualte min, max, mean, medians, sd
    for(i in names(ROIVoxelList)) {
      if((TRUE %in% is.na(ROIVoxelList[[i]])) == FALSE  ) {
        # consider only the voxel which are NOT ZERO
        listaGrigiDaConsiderare<-ROIVoxelList[[i]]$masked.images$voxelCube[which( !is.na(ROIVoxelList[[i]]$masked.images$voxelCube))]
        dataInfo[[i]]<-list()
        # collect the details
        dataInfo[[i]]$mean<-mean(listaGrigiDaConsiderare)
        dataInfo[[i]]$min<-min(listaGrigiDaConsiderare)
        dataInfo[[i]]$max<-max(listaGrigiDaConsiderare)
        dataInfo[[i]]$sd<-sd(listaGrigiDaConsiderare)
        dataInfo[[i]]$median<-median(listaGrigiDaConsiderare)
        # and get the summary
        min.arr<-c( min.arr, dataInfo[[i]]$min )
        max.arr<-c( max.arr, dataInfo[[i]]$max )
        mean.arr<-c( mean.arr, dataInfo[[i]]$mean )
        sd.arr<-c( sd.arr, dataInfo[[i]]$sd )
        median.arr<-c( median.arr, dataInfo[[i]]$median )
      }
      else {
        dataInfo[[i]]<-NA
      }
    }
    return(list(
      "details"=dataInfo,
      "summary"=list(
        "min"=min.arr,
        "max"=max.arr,
        "mean"=mean.arr,
        "sd"=sd.arr,
        "median"=median.arr)
    ))
  }
  
  # ========================================================================================
  # constructor
  # ========================================================================================
  constructor<-function( folderCleanUp , loadXMLInCache , loadRAWInCache ) {
    list_geoLet<<-list()
    attr_folderCleanUp<<-folderCleanUp                      # force to re-dump DICOM files
    attr_loadXMLInCache<<-loadXMLInCache
    attr_loadRAWInCache<<-loadRAWInCache
  }
  constructor( folderCleanUp = folderCleanUp, 
               loadXMLInCache = loadXMLInCache,
               loadRAWInCache = loadRAWInCache)
  
  return( list( "loadCollection"=loadCollection,
                "getAttribute"=getAttribute,
                "getROIVoxel"=getROIVoxel,
                "getImageVoxel"=getImageVoxel,
                "getAllPixelSpacing"=getAllPixelSpacing,
                "expandCube"=expandCube,
                "class"=class,
                "getROIVoxelStats"=getROIVoxelStats
  ) )
}
#' class for loading and presenting DICOM data
#' 
#' @description  Instantiate an object of the class \code{geoLet}.This represents just the classname, 
#'               methods are exposed with the technique of 'closure'.
#'               In order to see manuals for the single mathods, consider the vignette or use the 
#'               available for the following wrapping functions:
#'               \itemize{
#'               \item \code{GLT.openDICOMFolder( );} : to load a DICOM series into an geoLet object
#'               \item \code{GLT.getImageVoxelCube( );} : to get the ImageVoxelCube stored into a geoLet object
#'               \item \code{GLT.getPixelSpacing( );} : to get the pixelSpacing (x,y,z) of the main ImageVoxelCube stored into a geoLet object
#'               \item \code{GLT.getROIList( );} : to get the list of the ROI defined in a geoLet object
#'               \item \code{GLT.getTag( );} : to get a single DICOM-tag of a DICOM file loaded into a geoLet object
#'               \item \code{GLT.getROIVoxels( );} : to get the IMAGE Voxels geometrically located into a ROI, for a given geoLet object
#'               \item \code{GLT.extractDoseVoxels( );} : to get the DOSE Voxels geometrically located into a ROI, for a given geoLet object
#'               \item \code{GLT.calculateDVH( );} : to get the DVH calculated from a geoLet object
#'               }
#'               The original methods for the class geoLet can also be invocked using the same name without the previx 'GTL.', i.e.:
#' @examples \dontrun{
#' # first of all create an object obj.1 and obj.2
#' obj.1<-geoLet()
#' obj.2<-geoLet()
#' 
#' # now load a DICOM serie using the methods for obj.1 and the wrapping function for obj.2
#' # The result is the same: the only difference is some froceries available in the second way (i.e.: help online)
#' 
#' obj.1$openDICOMFolder(pathToOpen='./DICOMSeries/pat001' );
#' obj.2<-GLT.openDICOMFolder(obj = obj.2, pathToOpen='./DICOMSeries/pat001' );#' 
#' 
#' }
#' @export
#' @useDynLib moddicomV2 
#' @import stringr XML misc3d rgl Rvcg oce rmarkdown
geoLet<-function(ROIVoxelMemoryCache = TRUE , 
                 folderCleanUp = FALSE,
                 loadXMLInCache = FALSE,
                 loadRAWInCache = TRUE,
                 defaultExtension = "*.dcm"
) {
  # Attributes - set by user
  attr_folderCleanUp<-FALSE                              # force to re-dump DICOM files
  attr_ROIVoxelMemoryCache<-FALSE                        # force to cache ROI Voxel
  attr_ROIVoxelMemoryCacheArray<-list();
  attr_loadXMLInCache<-loadXMLInCache
  attr_loadRAWInCache<-loadRAWInCache
  attr_arrayXMLCache<-list()                              # array containing XML files
  attr_arrayRAWCache<-list()                              # array containing RAW files
  attr_defaultExtension<-c()
  # Internal Attributes
  attr_mainFrameOfReferenceUID<-NA                       # frameOfReference (geometry)
  logObj<-logHandler()                                   # log/error handler Object 
  dataStorage<-list()                                    # memory data structure
  attr_dataChache<-list();                               # Cache per i voxel delle ROI
  attr_attributeList<-''                                 # contiene attributi generici 
  SOPClassUIDList<-''                                    # Lista delle SOP Classes UID
  attr_ROI.non.compl<-''                                 # array delle ROI non complanari alle imgs
  #=================================================================================
  # openDICOMFolder
  # Loads a Folder containing one or more DICOM Studies
  #=================================================================================
  # Open a folder and load the content
  openDICOMFolder<-function(pathToOpen) {
    
    if(!dir.exists(pathToOpen)) logObj$handle( "error" , "The indicate Path does not exist"  );
    defaultExtension<-attr_defaultExtension;
    
    # get the dcm file type
    SOPClassUIDList<<-getFolderContent( pathToOpen , defaultExtension );   
    
    # Load RTPLan (if exists)
    logObj$sendLog("---------------------------------------")
    logObj$sendLog("Load RTPlan")
    logObj$sendLog("---------------------------------------")
    loadRTPlan( SOPClassUIDList = SOPClassUIDList  );
    
    # Load CT/RMN Scans
    logObj$sendLog("---------------------------------------")
    logObj$sendLog("Load Images")
    logObj$sendLog("---------------------------------------")
    loadCTRMNRDScans( SOPClassUIDList = SOPClassUIDList  );   
    
    # Load RTStruct Files
    logObj$sendLog("---------------------------------------")
    logObj$sendLog("Load RTStruct")
    logObj$sendLog("---------------------------------------")
    dataStorage[["structures"]]<<-loadRTStructFiles(SOPClassUIDList);    
    # Associate ROI and Images
    if(is.list(dataStorage[["structures"]])) {
      associateROIandImageSlices();
    }    
    # set the internal attribute indicating the path
    attr_attributeList[["path"]]<<-pathToOpen
    changeDVHROIIDInROINames();
  }
  #=================================================================================
  # NAME: loadRTPlan
  # load the RTPLAN
  #=================================================================================  
  loadRTPlan<-function( SOPClassUIDList , setValidCTRMNSeriesInstanceUID = 'any') {
    logObj$sendLog("... to be implemented properly")
  }
  #=================================================================================
  # getROIList
  # restituisce la lista delle ROI
  #=================================================================================  
  getROIList<-function() {
    mat2Ret<-matrix( c(seq(1,length(names(dataStorage$structures))),names(dataStorage$structures)),nrow=2 ,byrow=T )
    return(mat2Ret)
  }  
  #=================================================================================
  # loadRTStructFiles
  # Loads a DICOM RT Struct (one x folder)
  #=================================================================================
  loadRTStructFiles<-function(SOPClassUIDList) {
    
    imageSerie<-list()
    listaPuntiROI<-list()
    # loop over the list    
    # even if the assumption is that only one RTStruct is admitted for a CT scan serie
    TMP<-list()
    for(i in names(SOPClassUIDList)) {
      if(  SOPClassUIDList[[i]]$kind=="RTStructureSetStorage") {
        logObj$sendLog(i) 
        TMP[[i]]<-getStructuresFromXML( i );
      }
    }
    if(length(TMP)==0) {
      return( TMP );
    }
    
    # now let me use some more easy to handle variable names
    matrice2<-c(); matrice3<-c(); FORUID.m<-NA;
    for(i in names(TMP)) {
      matrice2<-cbind(matrice2,TMP[[i]]$IDROINameAssociation)
      matrice3<-rbind(matrice3,TMP[[i]]$tableROIPointList)
      # Aggiungi le informazioni relative al FrameOfReferenceUID delle ROI caricate
      # ed il ReferencedROINumber!!!!! (xè non è il numero della ROI di moddicom, possono essere diversi)
      if(!is.list(dataStorage$info[["structures"]])) dataStorage$info[["structures"]]<-list();
      for( nomeROI in TMP[[i]]$IDROINameAssociation[2,] ) {
        dataStorage$info[["structures"]][[nomeROI]]<<-list();
        dataStorage$info[["structures"]][[nomeROI]]$FrameOfReferenceUID<<-TMP[[i]]$FORUID.m
        dataStorage$info[["structures"]][[nomeROI]]$SeriesInstanceUID<<-TMP[[i]]$RTStructSeriesInstanceUID
      }      
    }
    listaROI<-list()
    
    # for each ROI
    for(i in matrice2[2,]) {
      
      # get the points
      subMatrix<-matrice3[which(matrice3[,2]==i,arr.ind = TRUE),]
      # if some points exist
      quantiElementiTrovati<--1
      
      if(is.list(subMatrix) & !is.array(subMatrix)) quantiElementiTrovati<-1
      if(is.matrix(subMatrix) & is.array(subMatrix)) quantiElementiTrovati<-dim(subMatrix)[1]
      if(quantiElementiTrovati==-1) {
        logObj$handle( "error" , "Unexpected error in loading slices. No slices found."  );
      }
      
      if( quantiElementiTrovati >0 ) {
        listaROI[[i]]<-list()
        
        # add properly the points to the 'listaROI' structure
        for(contatore in seq(1,quantiElementiTrovati) ) {
          
          if( quantiElementiTrovati == 1) {
            ROIPointStringList<-subMatrix[[3]]
            SOPInstance<-subMatrix[[4]]
          }
          else {
            ROIPointStringList<-subMatrix[contatore,3][[1]]
            SOPInstance<-subMatrix[contatore,4][[1]]
          }
          listaCoords<-strsplit(ROIPointStringList,"\\\\");
          listaCoords<-as.numeric(listaCoords[[1]])
          # if a ROI already exists for the slice, well append it to the list
          if( !( SOPInstance  %in% names(listaROI[[i]])  ) )  listaROI[[i]][[   SOPInstance  ]]<-list()          
          listaROI[[i]][[   SOPInstance  ]][[ length(listaROI[[i]][[   SOPInstance  ]])+1  ]]<-matrix(listaCoords,ncol=3,byrow=T)
          # Add the first one as last (close the loop)
          listaROI[[i]][[   SOPInstance  ]][[ length(listaROI[[i]][[   SOPInstance  ]]) ]]<-rbind(listaROI[[i]][[   SOPInstance  ]][[length(listaROI[[i]][[   SOPInstance  ]])]],listaROI[[i]][[   SOPInstance  ]][[length(listaROI[[i]][[   SOPInstance  ]])]][1,])
        }   
      } else {
        listaROI[[i]]<-NA
      }
    }
    check.not.coplanar.ROI(listaROI)
    return(listaROI); 
  }  
  #=================================================================================
  # check.not.coplanar.ROI
  # cerca di verificare la coplanarietà di tutte le ROI rispetto alle immagini. 
  # Quelle che non lo sono verranno iscritte in un array globale 'attr_ROI.non.compl'
  # l'eventuale calcolo dei ROIVoxel passerà dall'utilizzo delle mesh 
  # IPOTESI: le ROI sono definite su piani paralleli (sennò ecchecazzo!)
  #=================================================================================   
  check.not.coplanar.ROI<-function(listaROI) {
    #     # Per ogni ROI
    #     for ( ROIName in names(listaROI)) {
    #       # prendi la fetta con più punti censiti
    #       for(SOPInstance.ROI in names(listaROI[[ROIName]])) {
    #         for(closedSequence in names(listaROI[[ROIName]][[closedSequence]]))
    #         listaROI[[ROIName]][[SOPInstance.ROI]]
    #         browser()
    #       }
    #     }
    attr_ROI.non.compl<<-c()
  }
  #=================================================================================
  # changeDVHROIIDInROINames
  # at the end of the computation, change the ROIId in the ROINames. This cannot be 
  # done at the beginning because DicomRT object can be loaded in any order
  #=================================================================================   
  changeDVHROIIDInROINames<-function() {
    matriceNomiROI<-getROIList();
    if(!is.list(dataStorage$info$DVHs)) return;
    if(!is.list(dataStorage$info$DVHs[[1]])) return;
    for( SOPInstanceUID in names(dataStorage$info$DVHs) ) {
      listaDaSostituire<-list();
      for(indColonna in seq(1,dim(matriceNomiROI)[2] )) {
        for(numericID in names( dataStorage$info$DVHs[[SOPInstanceUID]]$DVHFromFile )) {
          nuovoNome<-matriceNomiROI[2,which(matriceNomiROI[1,]==numericID,arr.ind = T)]
          if(length(nuovoNome)>0) {
            listaDaSostituire[[nuovoNome]]<-dataStorage$info$DVHs[[SOPInstanceUID]]$DVHFromFile[[numericID]]
          }
        }
      }
      dataStorage$info$DVHs[[SOPInstanceUID]]$DVHFromFile<<-list();
      dataStorage$info$DVHs[[SOPInstanceUID]]$DVHFromFile<<-listaDaSostituire
    }
  }  
  #=================================================================================
  # getDoseVoxelCube
  # restituisce il voxel cube della dose
  #=================================================================================    
  getDoseVoxelCube<-function() {
    doseVC<-dataStorage$dose[[1]]
    doseInfo<-dataStorage$info$doses[[1]]
    return( list("voxelCube"=doseVC,"info"=doseInfo)  )
  }  
  #=================================================================================
  # getXYZFromNxNyNzOfImageVolume
  #=================================================================================    
  getXYZFromNxNyNzOfImageVolume<-function(Nx, Ny, Nz, startFromZero=FALSE) {
    objS<-services();
    
    serieInstanceUID<-giveBackImageSeriesInstanceUID();
    istanceNumbers<-as.character(sort(as.numeric(names(dataStorage$info[[serieInstanceUID]]))))
    if(startFromZero==FALSE) ct<-1;
    if(startFromZero==TRUE) ct<-0;
    for(slice in istanceNumbers) {
      if(ct==Nz) {
        imagePosition<-dataStorage$info[[serieInstanceUID]][[slice]]$ImagePositionPatient;
        ImageOrientationPatient<-dataStorage$info[[serieInstanceUID]][[slice]]$ImageOrientationPatient;
        pixelSpacing<-dataStorage$info[[serieInstanceUID]][[slice]]$pixelSpacing;
        orientationMatrix<-dataStorage$info[[serieInstanceUID]][[slice]]$orientationMatrix;
        punto<-objS$get3DPosFromNxNy(Nx,Ny,orientationMatrix);
        #        print(punto);
        return(punto);
      }
      ct<-ct+1;
    }
    stop();
    print(c(Nx,Ny,Nz))
  }  
  #=================================================================================
  # NAME: get.MESH.from.ROI
  # prendi una mesh dalla roi
  #=================================================================================    
  getMESHfromROI<-function( ROIName ,   
                            decimation=FALSE, decimation.percentage=0.8, 
                            smoothing=FALSE, smoothing.iterations = 10 ) {
    objS<-services();
    
    pixelSpacing<-getPixelSpacing();
    # ---------------------------------------------------
    # Prendi la ROI e informazioni preliminari
    # ---------------------------------------------------  
    ROIVoxels<-getROIVoxels(Structure = ROIName)
    ROILocations<-ROIVoxels$masked.images$location
    pip.arr<-ROIVoxels$masked.images$voxelCube
    pip.arr.exploded<-objS$expandCube(littleCube = ROIVoxels$masked.images$voxelCube,
                                      x.start = ROIVoxels$masked.images$location$min.x,
                                      y.start = ROIVoxels$masked.images$location$min.y,
                                      z.start = ROIVoxels$masked.images$location$min.z,
                                      fe = ROIVoxels$masked.images$location$fe,
                                      se = ROIVoxels$masked.images$location$se,te = ROIVoxels$masked.images$location$te)
    geomInfo<-getGeometricalInformationOfImage();

    # fai la mesh (senza visualizzare)
    x.coor <- seq(0, ((dim(pip.arr.exploded)[1]-1)*pixelSpacing[1]), by=pixelSpacing[1])
    y.coor <- seq(0, ((dim(pip.arr.exploded)[2]-1)*pixelSpacing[2]), by=pixelSpacing[2])
    z.coor <- seq(0, ((dim(pip.arr.exploded)[3]-1)*pixelSpacing[3]), by=pixelSpacing[3])
    
    threshold<-min(pip.arr.exploded[which(!is.na(pip.arr.exploded))])-1000  # prendi il minimo e togli ancora 1000
    pip.arr.exploded[which(is.na(pip.arr.exploded),arr.ind = T )]<-threshold
    #mesh.triangle<-contour3d(f = pip.arr.exploded, level = 1, x = x.coor, y = y.coor, rev(z.coor), engine = "none")
    mesh.triangle<-contour3d(f = pip.arr.exploded, level = threshold+10, x = x.coor, y = y.coor, rev(z.coor), engine = "none")
    
    # passa dai triangoli alle mesh
    mesh<-objS$triangle2mesh(x = mesh.triangle) 
    # pulizia
    mesh<-vcgClean(mesh = mesh, sel = c(0,0,1,1,2,2,3,3,4,4,5,5,6,6,7,7,0,0,7,7))  
    if(decimation == TRUE) mesh<-vcgQEdecim(mesh = mesh, percent = decimation.percentage);
    if(smoothing == TRUE) mesh<-vcgSmooth(mesh = mesh, iteration = smoothing.iterations)
    
    return(mesh)
    
  }
  #=================================================================================
  # NAME: extractDoseVoxels
  # estrae i voxel di dose interni ad una ROI
  #=================================================================================    
  extractDoseVoxels<-function(ROIName ,newPixelSpacing=NA, plotIT = FALSE, 
                              verbose=FALSE, forceReCalculus=FALSE, fastEngine = TRUE,
                              decimation=FALSE, decimation.percentage=0.8, 
                              smoothing=FALSE, smoothing.iterations = 10,
                              interpolate.dose = TRUE) {
    objS<-services();
    cat("\n --------------------------------------------- ")
    cat("\n Attenzione: il calcolo della dose deve essere ancora")
    cat("\n validato ")
    cat("\n --------------------------------------------- ")
    # verifica se è in cache
    if(forceReCalculus==FALSE) {
      if(!is.list(attr_dataChache)) attr_dataChache<<-list();
      if(!is.list(attr_dataChache$calculatedDVH)) attr_dataChache$calculatedDVH<<-list();
      if(!is.list(attr_dataChache$calculatedDVH[[ROIName]])) attr_dataChache$calculatedDVH[[ROIName]]<<-list() 
      else return (attr_dataChache$calculatedDVH[[ROIName]])
    }
    # non è in cache... elabora!
    if(length(newPixelSpacing)==1) newPixelSpacing<-getPixelSpacing();
    pixelSpacing<-getPixelSpacing();
    # ---------------------------------------------------
    # Prendi la ROI e informazioni preliminari
    # ---------------------------------------------------  
    
    ROIVoxels<-getROIVoxels(Structure = ROIName)
    ROILocations<-ROIVoxels$masked.images$location
    pip.arr<-ROIVoxels$masked.images$voxelCube
    pip.arr.exploded<-objS$expandCube(littleCube = ROIVoxels$masked.images$voxelCube,
                                      x.start = ROIVoxels$masked.images$location$min.x,
                                      y.start = ROIVoxels$masked.images$location$min.y,
                                      z.start = ROIVoxels$masked.images$location$min.z,
                                      fe = ROIVoxels$masked.images$location$fe,
                                      se = ROIVoxels$masked.images$location$se,te = ROIVoxels$masked.images$location$te)
    geomInfo<-getGeometricalInformationOfImage();
    
    # prendi il DoseVoxelCube
    doseList<-getDoseVoxelCube();
    doseVoxelCube<-doseList$voxelCube
    doseInfo<-doseList$info
    CTVC<-getImageVoxelCube();
    
    # ---------------------------------------------------
    # costruisci gli assi con le coordinate (si tratta delle coordinate dei voxel della CT interni al
    # voxelCube della ROI di interesse, ovviamente CALCOLATI CON IL PIXELSPACING della CT)
    # ---------------------------------------------------    
    # ROIBoundingBox 
    bBox.min.x<-getXYZFromNxNyNzOfImageVolume(Nx = ROILocations$min.x, Ny = ROILocations$min.y, Nz = ROILocations$min.z)[1]
    bBox.min.y<-getXYZFromNxNyNzOfImageVolume(Nx = ROILocations$min.x, Ny = ROILocations$min.y, Nz = ROILocations$min.z)[2]
    bBox.min.z<-getXYZFromNxNyNzOfImageVolume(Nx = ROILocations$min.x, Ny = ROILocations$min.y, Nz = ROILocations$min.z)[3]
    bBox.max.x<-getXYZFromNxNyNzOfImageVolume(Nx = ROILocations$max.x, Ny = ROILocations$max.y, Nz = ROILocations$max.z)[1]
    bBox.max.y<-getXYZFromNxNyNzOfImageVolume(Nx = ROILocations$max.x, Ny = ROILocations$max.y, Nz = ROILocations$max.z)[2]
    bBox.max.z<-getXYZFromNxNyNzOfImageVolume(Nx = ROILocations$max.x, Ny = ROILocations$max.y, Nz = ROILocations$max.z)[3]
    
    min.x<-0; max.x<-as.numeric(geomInfo$Columns)-1;
    min.y<-0; max.y<-as.numeric(geomInfo$Rows)-1;
    min.z<-1; max.z<-geomInfo$supposedNumberOfSlices;
    
    fromX<-getXYZFromNxNyNzOfImageVolume(Nx = min.x, Ny = 0, Nz = 1)[1];
    toX<-getXYZFromNxNyNzOfImageVolume(Nx = max.x, Ny = 0, Nz = 1)[1];
    fromY<-getXYZFromNxNyNzOfImageVolume(Nx = 0, Ny = min.y, Nz = 1)[2];
    toY<-getXYZFromNxNyNzOfImageVolume(Nx = 0, Ny = max.y, Nz = 1)[2];
    fromZ<-getXYZFromNxNyNzOfImageVolume(Nx = 0, Ny = 0, Nz = min.z)[3];
    toZ<-getXYZFromNxNyNzOfImageVolume(Nx = 0, Ny = 0, Nz = max.z)[3];
    
    if( (toX-fromX)>0 ) x.coor<-seq(fromX, toX, by = pixelSpacing[1]    )
    else x.coor<-seq(toX, fromX, by = pixelSpacing[1]    )
    if( (toY-fromY)>0 ) y.coor<-seq(fromY, toY, by =  pixelSpacing[2]    )
    else y.coor<-seq(toY, fromY, by =  pixelSpacing[2]    )
    if( (toZ-fromZ)>0 ) z.coor<-seq(fromZ, toZ, by =  as.numeric(geomInfo$SliceThickness)    )
    else z.coor<-seq(toZ, fromZ, by =  as.numeric(geomInfo$SliceThickness)   )
    
    # ---------------------------------------------------
    # calcola le MESH
    # dalla struttura della ROI, nell'intero spazio della CT
    # (fin qui lavoro ancora solo su CT e ROI, le dosi non compaiono)
    # ---------------------------------------------------
    # fai la mesh (senza visualizzare)
    
    threshold<-min(pip.arr.exploded[which(!is.na(pip.arr.exploded))])-1000  # prendi il minimo e togli ancora 1000
    pip.arr.exploded[which(is.na(pip.arr.exploded),arr.ind = T )]<-threshold
    #mesh.triangle<-contour3d(f = pip.arr.exploded, level = 1, x = x.coor, y = y.coor, rev(z.coor), engine = "none")
    mesh.triangle<-contour3d(f = pip.arr.exploded, level = threshold+10, x = x.coor, y = y.coor, rev(z.coor), engine = "none")
    
    # passa dai triangoli alle mesh
    mesh<-objS$triangle2mesh(x = mesh.triangle) 
    # pulizia
    mesh<-vcgClean(mesh = mesh, sel = c(0,0,1,1,2,2,3,3,4,4,5,5,6,6,7,7,0,0,7,7))  
    if(decimation == TRUE) mesh<-vcgQEdecim(mesh = mesh, percent = decimation.percentage);
    if(smoothing == TRUE) mesh<-vcgSmooth(mesh = mesh, iteration = smoothing.iterations)
    
    # ---------------------------------------------------
    # interpola la DOSE
    # ---------------------------------------------------  
    # imposta le dinamiche lungo i 3 assi con i NUOVI PIXELSPACING
    # (che se non specificati saranno gli stessi di default della CT)
    resampling.coords.x<-seq(min(x.coor),max(x.coor), by = newPixelSpacing[1])
    resampling.coords.y<-seq(min(y.coor),max(y.coor), by = newPixelSpacing[2])
    resampling.coords.z<-seq(min(z.coor),max(z.coor), by = newPixelSpacing[3])
    
    # fai l'expand.grid
    resampling.coords.grid<-expand.grid(resampling.coords.x,resampling.coords.y,resampling.coords.z)
    
    pptmp1<-doseInfo$pixelSpacing
    pptmp2<-doseInfo$GridFrameOffsetVector
    dose.px<-pptmp1[1];  dose.py<-pptmp1[2];  dose.pz<-abs(pptmp2[1]-pptmp2[2])
    
    interpolatedValues<-array(0,nrow(resampling.coords.grid))
    ps.x<-doseInfo$pixelSpacing[1]; ps.y<-doseInfo$pixelSpacing[2];
    
    f.x<-getXYZFromNxNyNzOfDoseVolume(Nx = 0,Ny = 0,Nz = 1)[1]
    f.y<-getXYZFromNxNyNzOfDoseVolume(Nx = 0,Ny = 0,Nz = 1)[2]
    t.x<-getXYZFromNxNyNzOfDoseVolume(Nx = as.numeric(doseInfo$Columns)-1 ,Ny = 0,Nz = 1)[1];
    t.y<-getXYZFromNxNyNzOfDoseVolume(Nx = 0 ,Ny = as.numeric(doseInfo$Rows)-1,Nz = 1)[2];
    
    if(t.x>f.x) {  xDosePointsCoord<-seq(f.x, t.x,by=ps.x ) }
    else { xDosePointsCoord<-seq(t.x, f.x,by=ps.x )   }
    if( t.y> f.y ){  yDosePointsCoord<-seq(f.y, t.y,by=ps.y  )  }
    else {  yDosePointsCoord<-seq(t.y, f.y,by=ps.y  )  }
    zDosePointsCoord<-doseInfo$GridFrameOffsetVector+doseInfo$imagePositionPatient[3];
    # shifta il valore della Y rispetto al valore centrale
    # (raffinatissimo problema di ribaltamento della y non percepibile dato il passaggio
    # alla gestione in punti nello spazio)
    delta1<-min(yDosePointsCoord) - min(resampling.coords.grid[,2])
    delta2<-max(resampling.coords.grid[,2]) - max(yDosePointsCoord)
    deltaT<-delta1-delta2;
    yDosePointsCoord<-yDosePointsCoord-deltaT
    xDosePointsCoord<-xDosePointsCoord-ps.x
    if(interpolate.dose == TRUE) {
      interpolata <- approx3d(x = xDosePointsCoord, y = yDosePointsCoord, z = zDosePointsCoord, 
                              f = doseVoxelCube, 
                              xout = as.array(resampling.coords.grid[,1]),
                              yout = as.array(resampling.coords.grid[,2]), 
                              zout = as.array(resampling.coords.grid[,3]))  
    } else {
      interpolata <- getCloserDoseValue(x = xDosePointsCoord, y = yDosePointsCoord, z = zDosePointsCoord, 
                                        f = doseVoxelCube, 
                                        xout = as.array(resampling.coords.grid[,1]),
                                        yout = as.array(resampling.coords.grid[,2]), 
                                        zout = as.array(resampling.coords.grid[,3])) 
      browser();
    }  
    
    # costruisci la matrice 3D
    matriciona<-array(interpolata,dim=c(length(resampling.coords.x),length(resampling.coords.y),length(resampling.coords.z)   ))
    # e ribalta l'asse Z
    matriciona[,,seq(1,dim(matriciona)[3])]<-matriciona[,,seq(dim(matriciona)[3],1,by=-1)]
    # ---------------------------------------------------
    # ora calcola PUNTI INTERNI/PUNTI ESTERNI
    # ---------------------------------------------------   
    # prendi la CT originale
    # ricampionala
    CTVC.interpolata <- approx3d(x = x.coor, y = y.coor, z = rev(z.coor), 
                                 f = CTVC, 
                                 xout = as.array(resampling.coords.grid[,1]),
                                 yout = as.array(resampling.coords.grid[,2]), 
                                 zout = as.array(resampling.coords.grid[,3]))  
    # costruisci la matrice 3D
    CTInterpolata<-array(CTVC.interpolata,dim=c(length(resampling.coords.x),length(resampling.coords.y),length(resampling.coords.z)   ))
    # e ribalta l'asse Z
    CTInterpolata[,,seq(1,dim(CTInterpolata)[3])]<-CTInterpolata[,,seq(dim(CTInterpolata)[3],1,by=-1)]  
    
    voxelCube.CTInterpolata<-array(0,dim=c(  length(resampling.coords.x),length(resampling.coords.y),length(resampling.coords.z)   ))
    #voxelCube.DoseInterpolata<-array(0,dim=c(  length(resampling.coords.x),length(resampling.coords.y),length(resampling.coords.z)   ))
    voxelCube.DoseInterpolata<-array(NA,dim=c(  length(resampling.coords.x),length(resampling.coords.y),length(resampling.coords.z)   ))
    
    casted.coords.x<-resampling.coords.x[resampling.coords.x>=bBox.min.x & resampling.coords.x<=bBox.max.x]
    casted.coords.y<-resampling.coords.y[resampling.coords.y>=bBox.min.y & resampling.coords.y<=bBox.max.y]
    casted.coords.z<-resampling.coords.z[resampling.coords.z<=bBox.min.z & resampling.coords.z>=bBox.max.z]
    
    firstX<-which(resampling.coords.x>=bBox.min.x & resampling.coords.x<=bBox.max.x)[1]
    firstY<-which(resampling.coords.y>=bBox.min.y & resampling.coords.y<=bBox.max.y)[1]
    if(verbose==TRUE) {cat(  paste(  c("\n|",rep("-",length(resampling.coords.z)),"|\n|"),collapse='')    ) }
    
    # loop per ogni slice lungo la cranio-caudale
    ct<-1;
    for(sliceRunner in seq(1,length(resampling.coords.z) )) {
      
      sliceZLocation<-rev(resampling.coords.z)[sliceRunner]
      
      # se la coordinata z è interna al campo di interesse
      # allora preparati a riflettere sui punti interni/esterni
      if(bBox.max.z<rev(resampling.coords.z)[sliceRunner] &
         bBox.min.z>rev(resampling.coords.z)[sliceRunner] 
      ) {
        if(verbose==TRUE) {cat("*")}
        casted.coords.x<-resampling.coords.x[resampling.coords.x>=bBox.min.x & resampling.coords.x<=bBox.max.x]
        casted.coords.y<-resampling.coords.y[resampling.coords.y>=bBox.min.y & resampling.coords.y<=bBox.max.y]
        casted.coords.z<-resampling.coords.z[resampling.coords.z<=bBox.min.z & resampling.coords.z>=bBox.max.z]
        
        firstX<-which(resampling.coords.x>=bBox.min.x & resampling.coords.x<=bBox.max.x)[1]
        firstY<-which(resampling.coords.y>=bBox.min.y & resampling.coords.y<=bBox.max.y)[1]
        
        resampling.coords.grid.single.slice<-expand.grid(casted.coords.x,casted.coords.y,rev(resampling.coords.z)[sliceRunner])
        
        puntoInCulonia.00<-c(3*bBox.min.x,3*bBox.min.y,rev(resampling.coords.z)[sliceRunner][1]);
        puntoInCulonia.01<-c(3*bBox.max.x,3*bBox.max.y,rev(resampling.coords.z)[sliceRunner][1]);
        puntoInCulonia.02<-c(3*bBox.min.x,3*bBox.max.y,rev(resampling.coords.z)[sliceRunner][1]);
        puntoInCulonia.03<-c(3*bBox.max.x,3*bBox.min.y,rev(resampling.coords.z)[sliceRunner][1]);
        
        matricePunti<-rbind(puntoInCulonia.00,puntoInCulonia.01,puntoInCulonia.01,puntoInCulonia.01,as.matrix(resampling.coords.grid.single.slice))
        resampling.coords.grid.single.slice<-matricePunti
        
        if ( fastEngine == TRUE ) clost <- vcgClostKD(as.matrix(resampling.coords.grid.single.slice), mesh)
        else clost <- vcgClost(as.matrix(resampling.coords.grid.single.slice), mesh)
        
        qualita<-clost$quality
        qualita<-qualita[5: (length(qualita)) ]
        
        arrayOne<-array(0,length(resampling.coords.grid.single.slice));
        arrayOne[which(qualita<0)]<-1
        #arrayOne[which(is.na(arrayOne))]<-0
        arrayOne[which(is.na(arrayOne))]<-NA
        
        bigMaskCube<-array(arrayOne,dim=c( length(casted.coords.x),length(casted.coords.y),length(resampling.coords.z[sliceRunner])  ))
        bigMaskCubeTMP<-array(NA,dim=c( length(resampling.coords.x),length(resampling.coords.y), 1  ))
        bigMaskCubeTMP[firstX:(firstX+dim(bigMaskCube)[1]-1), firstY:(firstY+dim(bigMaskCube)[2]-1),1 ] <-bigMaskCube[,,1]
        
        if(plotIT == TRUE) {
          image(CTInterpolata[,,sliceRunner], col = grey.colors(255*255))  
          image(matriciona[,,sliceRunner]*bigMaskCubeTMP[,,1], col = heat.colors(20, alpha = 0.5),add = TRUE)  
        }
        
        voxelCube.CTInterpolata[,,ct]<-CTInterpolata[,,sliceRunner]
        #voxelCube.DoseInterpolata[,,ct]<-matriciona[,,sliceRunner]*bigMaskCubeTMP[,,1]
        bigMaskCubeTMP[ which( bigMaskCubeTMP ==0,arr.ind = TRUE)   ]<-NA
        voxelCube.DoseInterpolata[,,ct]<-matriciona[,,sliceRunner]*bigMaskCubeTMP[,,1]
        ct<-ct+1;
      }
      else { if(verbose==TRUE) {cat(".")} }
    }
    
    if(verbose==TRUE) {cat("|");}
    if(!is.list(attr_dataChache)) attr_dataChache<<-list();
    if(!is.list(attr_dataChache$calculatedDVH)) attr_dataChache$calculatedDVH<<-list();
    if(!is.list(attr_dataChache$calculatedDVH[[ROIName]])) attr_dataChache$calculatedDVH[[ROIName]]<<-list();
    res<-list("voxelCube.CT"=voxelCube.CTInterpolata,"voxelCube.Dose"=voxelCube.DoseInterpolata,"mesh"=mesh);
    attr_dataChache$calculatedDVH[[ROIName]]<<-res
    return( res ) ;
  }    
  #=================================================================================
  # NAME: extractDoseVoxels
  # estrae i voxel di dose interni ad una ROI
  #=================================================================================    
  extractDoseVoxels.part<-function(ROIName ,newPixelSpacing=NA, plotIT = FALSE, 
                                   verbose=FALSE, forceReCalculus=FALSE, fastEngine = TRUE,
                                   decimation=FALSE, decimation.percentage=0.8, 
                                   smoothing=FALSE, smoothing.iterations = 10,
                                   interpolate.dose = TRUE) {
    objS<-services();
    stop("no more implemented: nemmeno più ricordo che doveva fare, sto' metodo...")
    
    # verifica se è in cache
    if(forceReCalculus==FALSE) {
      if(!is.list(attr_dataChache)) attr_dataChache<<-list();
      if(!is.list(attr_dataChache$calculatedDVH)) attr_dataChache$calculatedDVH<<-list();
      if(!is.list(attr_dataChache$calculatedDVH[[ROIName]])) attr_dataChache$calculatedDVH[[ROIName]]<<-list() 
      else return (attr_dataChache$calculatedDVH[[ROIName]])
    }
    # non è in cache... elabora!
    if(length(newPixelSpacing)==1) newPixelSpacing<-getPixelSpacing();
    pixelSpacing<-getPixelSpacing();
    # ---------------------------------------------------
    # Prendi la ROI e informazioni preliminari
    # ---------------------------------------------------  
    
    ROIVoxels<-getROIVoxels(Structure = ROIName)
    ROILocations<-ROIVoxels$masked.images$location
    pip.arr<-ROIVoxels$masked.images$voxelCube
    pip.arr.exploded<-objS$expandCube(littleCube = ROIVoxels$masked.images$voxelCube,
                                      x.start = ROIVoxels$masked.images$location$min.x,
                                      y.start = ROIVoxels$masked.images$location$min.y,
                                      z.start = ROIVoxels$masked.images$location$min.z,
                                      fe = ROIVoxels$masked.images$location$fe,
                                      se = ROIVoxels$masked.images$location$se,te = ROIVoxels$masked.images$location$te)
    geomInfo<-getGeometricalInformationOfImage();
    
    # prendi il DoseVoxelCube
    #     doseList<-getDoseVoxelCube();
    #     doseVoxelCube<-doseList$voxelCube
    #     doseInfo<-doseList$info
    CTVC<-getImageVoxelCube();
    
    # ---------------------------------------------------
    # costruisci gli assi con le coordinate (si tratta delle coordinate dei voxel della CT interni al
    # voxelCube della ROI di interesse, ovviamente CALCOLATI CON IL PIXELSPACING della CT)
    # ---------------------------------------------------    
    # ROIBoundingBox 
    bBox.min.x<-getXYZFromNxNyNzOfImageVolume(Nx = ROILocations$min.x, Ny = ROILocations$min.y, Nz = ROILocations$min.z)[1]
    bBox.min.y<-getXYZFromNxNyNzOfImageVolume(Nx = ROILocations$min.x, Ny = ROILocations$min.y, Nz = ROILocations$min.z)[2]
    bBox.min.z<-getXYZFromNxNyNzOfImageVolume(Nx = ROILocations$min.x, Ny = ROILocations$min.y, Nz = ROILocations$min.z)[3]
    bBox.max.x<-getXYZFromNxNyNzOfImageVolume(Nx = ROILocations$max.x, Ny = ROILocations$max.y, Nz = ROILocations$max.z)[1]
    bBox.max.y<-getXYZFromNxNyNzOfImageVolume(Nx = ROILocations$max.x, Ny = ROILocations$max.y, Nz = ROILocations$max.z)[2]
    bBox.max.z<-getXYZFromNxNyNzOfImageVolume(Nx = ROILocations$max.x, Ny = ROILocations$max.y, Nz = ROILocations$max.z)[3]
    
    min.x<-0; max.x<-as.numeric(geomInfo$Columns)-1;
    min.y<-0; max.y<-as.numeric(geomInfo$Rows)-1;
    min.z<-1; max.z<-geomInfo$supposedNumberOfSlices;
    
    fromX<-getXYZFromNxNyNzOfImageVolume(Nx = min.x, Ny = 0, Nz = 1)[1];
    toX<-getXYZFromNxNyNzOfImageVolume(Nx = max.x, Ny = 0, Nz = 1)[1];
    fromY<-getXYZFromNxNyNzOfImageVolume(Nx = 0, Ny = min.y, Nz = 1)[2];
    toY<-getXYZFromNxNyNzOfImageVolume(Nx = 0, Ny = max.y, Nz = 1)[2];
    fromZ<-getXYZFromNxNyNzOfImageVolume(Nx = 0, Ny = 0, Nz = min.z)[3];
    toZ<-getXYZFromNxNyNzOfImageVolume(Nx = 0, Ny = 0, Nz = max.z)[3];
    
    if( (toX-fromX)>0 ) x.coor<-seq(fromX, toX, by = pixelSpacing[1]    )
    else x.coor<-seq(toX, fromX, by = pixelSpacing[1]    )
    if( (toY-fromY)>0 ) y.coor<-seq(fromY, toY, by =  pixelSpacing[2]    )
    else y.coor<-seq(toY, fromY, by =  pixelSpacing[2]    )
    if( (toZ-fromZ)>0 ) z.coor<-seq(fromZ, toZ, by =  as.numeric(geomInfo$SliceThickness)    )
    else z.coor<-seq(toZ, fromZ, by =  as.numeric(geomInfo$SliceThickness)   )
    
    # ---------------------------------------------------
    # calcola le MESH
    # dalla struttura della ROI, nell'intero spazio della CT
    # (fin qui lavoro ancora solo su CT e ROI, le dosi non compaiono)
    # ---------------------------------------------------
    # fai la mesh (senza visualizzare)
    
    threshold<-min(pip.arr.exploded[which(!is.na(pip.arr.exploded))])-1000  # prendi il minimo e togli ancora 1000
    pip.arr.exploded[which(is.na(pip.arr.exploded),arr.ind = T )]<-threshold
    #mesh.triangle<-contour3d(f = pip.arr.exploded, level = 1, x = x.coor, y = y.coor, rev(z.coor), engine = "none")
    x.coor<- seq(fromX,   fromX+pixelSpacing[1]*(dim(pip.arr.exploded)[1]-1),by=pixelSpacing[1])
    y.coor<- seq(fromY,   fromY+pixelSpacing[2]*(dim(pip.arr.exploded)[2]-1),by=pixelSpacing[2])
    z.coor<- seq(fromZ,   fromZ+as.numeric(geomInfo$SliceThickness)*(dim(pip.arr.exploded)[3]-1),by=as.numeric(geomInfo$SliceThickness))
    mesh.triangle<-contour3d(f = pip.arr.exploded, level = threshold+10, x = x.coor, y = y.coor, rev(z.coor), engine = "none")
    
    # passa dai triangoli alle mesh
    mesh<-objS$triangle2mesh(x = mesh.triangle) 
    # pulizia
    mesh<-vcgClean(mesh = mesh, sel = c(0,0,1,1,2,2,3,3,4,4,5,5,6,6,7,7,0,0,7,7))  
    if(decimation == TRUE) mesh<-vcgQEdecim(mesh = mesh, percent = decimation.percentage);
    if(smoothing == TRUE) mesh<-vcgSmooth(mesh = mesh, iteration = smoothing.iterations)
    return(mesh)
  }   
  #=================================================================================
  # getXYZFromNxNyNzOfDoseVolume
  # Create the association between ROI and images
  #=================================================================================   
  getXYZFromNxNyNzOfDoseVolume<-function(Nx, Ny, Nz, startFromZero=FALSE) {
    objS<-services();
    if(length(dataStorage$info$doses)==1) {seriesInstanceUID<-names(dataStorage$info$doses)[[1]];}
    else {stop(" caso non ancora previsto ( più volumi di dose a questo livello #njj9)");} 
    pixelSpacing<-dataStorage$info$doses[[seriesInstanceUID]]$pixelSpacing
    imagePositionPatient<-dataStorage$info$doses[[seriesInstanceUID]]$imagePositionPatient
    ImageOrientationPatient<-dataStorage$info$doses[[seriesInstanceUID]]$ImageOrientationPatient
    mat<-array(dataStorage$info$doses[[seriesInstanceUID]]$ImageOrientationPatient,dim=c(3,2))
    mat<-rbind(mat,c(0,0));    mat<-cbind(mat,c(0,0,0,0));
    mat[,1]<-mat[,1]*pixelSpacing[1];   mat[,2]<-mat[,2]*pixelSpacing[2];
    mat<-cbind(mat,c(0,0,0,1)); mat[1:3,4]<-imagePositionPatient;
    mat[3,4]<-dataStorage$info$doses[[seriesInstanceUID]]$GridFrameOffsetVector[Nz]
    punto<-objS$get3DPosFromNxNy(Nx,Ny,mat);
    return(punto);
  }  
  #=================================================================================
  # calculateDVH
  # Create the association between ROI and images
  #=================================================================================  
  calculateDVH<-function(ROIName, newPixelSpacing=NA , justTheDVH=TRUE,
                         verbose=FALSE, forceReCalculus=FALSE, fastEngine = TRUE,
                         decimation=FALSE, decimation.percentage=0.8, 
                         smoothing=FALSE, smoothing.iterations = 10) {
    if(length(newPixelSpacing)==1) newPixelSpacing<-getPixelSpacing();
    voxelVolume<-newPixelSpacing[1]*newPixelSpacing[2]*newPixelSpacing[3];
    vv<-extractDoseVoxels(ROIName = ROIName, newPixelSpacing = newPixelSpacing,
                          verbose = verbose, forceReCalculus = forceReCalculus,
                          fastEngine = fastEngine, decimation = decimation, 
                          decimation.percentage = decimation.percentage,
                          smoothing = smoothing, 
                          smoothing.iterations = smoothing.iterations);
    
    voxelCube.Dose<-vv$voxelCube.Dose
    voxelCube.Dose<-array(voxelCube.Dose);
    voxelCube.Dose<-voxelCube.Dose[ which(!is.na(voxelCube.Dose) ) ]
    voxelCube.Dose<-voxelCube.Dose[ which( voxelCube.Dose!=0) ] 
    max.voxelCube.Dose <- 1.05 * max(voxelCube.Dose);
    voxelVolume<-voxelVolume/1000
    dose.bin<-0.25
    if( (max.voxelCube.Dose / dose.bin)<100  ) dose.bin<-max.voxelCube.Dose/100;
    a<-DVH.extract(x = voxelCube.Dose,dvh.type = 'cumulative',vol.distr = 'absolute',createObj = TRUE,voxel.volume = voxelVolume,max.dose = max.voxelCube.Dose,dose.bin = dose.bin)
    if( justTheDVH == TRUE) return(a);
    return( list("DVHobj"=a, "voxelCube.CT"=vv$voxelCube.CT, "voxelCube.Dose"=vv$voxelCube.Dose)  );
  }  
  #=================================================================================
  # associateROIandImageSlices
  # Create the association between ROI and images
  #=================================================================================
  associateROIandImageSlices<-function( relaxCoPlanarity = FALSE) {
    
    numeroROIDaAssegnare<-0;
    numeroROIAssegnate<-0;
    ROINames<-getROIList();
    objServ<-services();
    
    for(nomeROI in ROINames[2,]) {
      FrameOfReferenceUID<-dataStorage$info$structures[[nomeROI]]$FrameOfReferenceUID;
      # prendi la seriesInstanceUID con quel FrameOfReferenceUID
      # (per sapere quale serie di immagini avrà associata la ROI)
      
      seriesInstanceUID<-giveBackImageSeriesInstanceUID(FrameOfReferenceUID=FrameOfReferenceUID)
      # bene, ora per ogni ROI scorri i contorni sui vari piani assiali
      listaPuntiROI<-getROIPointList(ROINumber = nomeROI );
      
      # frulla su ogni assiale
      if(is.list(listaPuntiROI)) {
        for( indiciROI in seq(1,length(listaPuntiROI) )) {
          if(length(listaPuntiROI[[indiciROI]])>0) {
            numeroROIDaAssegnare<-numeroROIDaAssegnare+1
            # prendi un punto campione
            sampleROIPoint<-listaPuntiROI[[indiciROI]][[1]][1,]
            # ora cicla sulla serie di immagini identificata come pertinente 
            # l'indice è l'instance number
            for(  imgInstanceNumber  in names(dataStorage$img[[seriesInstanceUID]] ) ) {
              # prendi l'equazione del piano di quella fetta
              planeEquation<-dataStorage$info[[seriesInstanceUID]][[imgInstanceNumber]]$planeEquation
              # calcola la distanza fra il piano dell'immagine in esame ed il punto campione della ROI
              distanza<-objServ$getPointPlaneDistance(Punto=sampleROIPoint, Piano=planeEquation)
              # prendi la slice Thickness
              sliceThickness<-as.numeric(dataStorage$info[[seriesInstanceUID]][[imgInstanceNumber]]$SliceThickness) 
              # RSTruct SeriesInstanceUID
              STRUCTSeriesInstanceUID<-names(listaPuntiROI)[indiciROI]
              # verifica se la distanza è inferiore ad un dato delta, se sì
              # ciò per decidere se la ROI giace sulla slice
              if( abs(distanza)<0.1 |  ( abs(distanza)<=(sliceThickness/2) & relaxCoPlanarity==TRUE   ) ) {
                if( nomeROI %in% dataStorage$info[[seriesInstanceUID]][[imgInstanceNumber]][["ROIList"]][,1] ) {
                  logObj$handle( "error" , "the ROI seems to be associated two times to the same slice.. do we have a problem in recognizing points in space?"  );
                } else {
                  # se la tabella manco c'era
                  if(length(dataStorage$info[[seriesInstanceUID]][[imgInstanceNumber]][["ROIList"]])==0) {
                    
                    dataStorage$info[[seriesInstanceUID]][[imgInstanceNumber]][["ROIList"]]<<-matrix(0,ncol=4,nrow=1)
                    dataStorage$info[[seriesInstanceUID]][[imgInstanceNumber]][["ROIList"]][1,1]<<-nomeROI
                    dataStorage$info[[seriesInstanceUID]][[imgInstanceNumber]][["ROIList"]][1,2]<<-seriesInstanceUID
                    dataStorage$info[[seriesInstanceUID]][[imgInstanceNumber]][["ROIList"]][1,3]<<-FrameOfReferenceUID
                    dataStorage$info[[seriesInstanceUID]][[imgInstanceNumber]][["ROIList"]][1,4]<<-STRUCTSeriesInstanceUID
                    numeroROIAssegnate<-numeroROIAssegnate+1
                  } else {
                    # se invece c'era serve solo aggiungere una riga
                    riga<-nrow(dataStorage$info[[seriesInstanceUID]][[imgInstanceNumber]][["ROIList"]])+1
                    dataStorage$info[[seriesInstanceUID]][[imgInstanceNumber]][["ROIList"]]<<-rbind(dataStorage$info[[seriesInstanceUID]][[imgInstanceNumber]][["ROIList"]],c("","","",""))
                    dataStorage$info[[seriesInstanceUID]][[imgInstanceNumber]][["ROIList"]][riga,1]<<-nomeROI
                    dataStorage$info[[seriesInstanceUID]][[imgInstanceNumber]][["ROIList"]][riga,2]<<-seriesInstanceUID
                    dataStorage$info[[seriesInstanceUID]][[imgInstanceNumber]][["ROIList"]][riga,3]<<-FrameOfReferenceUID
                    dataStorage$info[[seriesInstanceUID]][[imgInstanceNumber]][["ROIList"]][riga,4]<<-STRUCTSeriesInstanceUID                
                    numeroROIAssegnate<-numeroROIAssegnate+1
                  }
                }
              }
            }
          }
        }
      }
    }
    if(numeroROIAssegnate!=numeroROIDaAssegnare) {
      logObj$handle(type="warning", msg=c("Some ROIs are not coplanar with the loaded image slices. Further ROIVoxel extraction will be performed by meshes (they could take a lot of time and memory)"))
    }
  }  
  #=================================================================================
  # getROIPointList
  # restituisce la lista delle coordinate dei vertici di una data ROI
  #=================================================================================    
  getROIPointList<-function(ROINumber) {    
    return(dataStorage$structures[ROINumber][[names(dataStorage$structures[ROINumber])]])
  }  
  #####################################################################################
  # getStructuresFromXML: carica il file xml del RT struct
  #   
  # INPUT:
  #   - fileName: il nome del file del RT struct
  # OUTPUT:
  #   -  
  #################################################################################  
  getStructuresFromXML<-function(fileName) {
    obj.S<-services();
    massimo<-0
    
    # Load the XML file if not in cache
    doc<-getXMLStructureFromDICOMFile(fileName = fileName, folderCleanUp = folderCleanUp)
    
    # prima di tutto controlla che la FrameOfReferenceUID sia la stessa OVUNQUE e che punti
    # ad una serie di immagini ESISTENTE!
    # E' un chiodo ma .... ragionevole, almeno per ora
    # Estrae la seriesIstanceUID, la frame of reference UID e il Referenced Frame of Reference UID
    RTStructSeriesInstanceUID<-xpathApply(doc,'//element[@tag="0020,000e" and @name="SeriesInstanceUID"]',xmlValue)[[1]]
    FORUID.m<-xpathApply(doc,'//element[@tag="0020,0052" and @name="FrameOfReferenceUID"]',xmlValue)[[1]]
    FORUID.d<-xpathApply(doc,'//element[@tag="3006,0024" and @name="ReferencedFrameOfReferenceUID"]',xmlValue)
    # Analizza tutti i valori del Referenced frame of reference UID e controlla se tutti i valori coincidono con il 
    # frame of reference UID
    for(FORUID.d_index in seq(1,length(FORUID.d))) {
      if( FORUID.d[[ FORUID.d_index ]] !=  FORUID.m ) {
        logObj$handle( "error" , "FrameOfReferenceUID not aligned in RTStruct file"  );
      }
    }
    if(is.na(giveBackImageSeriesInstanceUID(FrameOfReferenceUID = FORUID.m))) {
      logObj$handle( "error" , "the FrameOfReferenceUID of the RTStruct is not associated to an image"  );
    }
    
    # SEQUENCES: the one with the attribute  tag="3006,0020"  and name="StructureSetROISequence" 
    # is the one with association NAME<->ID
    # Estrazione della parte di xml contenente il nome,numero delle varie ROI
    n2XML<-getNodeSet(doc,'/file-format/data-set/sequence[@tag="3006,0020" and @name="StructureSetROISequence"]/item')
    # SEQUENCES: now get the true coords
    # Estrazione delle coordinate di ciascuna ROI
    n3XML<-getNodeSet(doc,'/file-format/data-set/sequence[@tag="3006,0039" and @name="ROIContourSequence"]/item')
    
    # ROI Names
    matrice2<-c()
    # Per ciascuna ROI viene estratto il numero, il nome e organizzati in una matrice
    for(i in n2XML) {
      ROINumber<-xpathApply(xmlDoc(i),'/item/element[@tag="3006,0022"]',xmlValue)[[1]]
      ROIName<-xpathApply(xmlDoc(i),'/item/element[@tag="3006,0026"]',xmlValue)[[1]]
      matrice2<-rbind(matrice2,c(ROINumber,ROIName))
    }
    
    matrice2<-t(matrice2)
    # ROI Point list
    massimo<-0
    
    matrice3<-c()
    
    # esegue una interazione per ogni ROI 
    # Non fa altro che estrarre da 'i' la parte che contiene i vertici della ROI.
    for(i in n3XML) {
      
      ROINumber<-xpathApply(xmlDoc(i),'/item/element[@tag="3006,0084"]',xmlValue)[[1]]
      ROIName<-matrice2[2,which(matrice2[1,]==ROINumber)]
      # browser()
      
      
      # la funzione getNodeSet() restituisce l'insieme delle coordinate dei vertici di una sola ROI
      # (sto intendendo per ROI l'insieme di tutti i contorni con lo stesso ROIName)
      listaPuntiDaRavanare<-getNodeSet(xmlDoc(i),'/item/sequence/item')
      
      numero.Punti.semiperimetro.massimo<-0
      
      # Estrae solo un punto (fPoint.x, fPoint.y, fPoint.z) da listaPuntiDaRavanare in quanto viene assunto
      # che la ROI rispetto alla immagine sono tra loro paralleli (non sempre necessariamente vero)
      for(i2 in listaPuntiDaRavanare)   {
        
        # ReferencedSOPInstanceUID<-xpathApply(xmlDoc(i2),'//element[@tag="0008,1155"]',xmlValue)[[1]]
        # salva le coordinate delle ROI in una lista
        ROIPointList<-xpathApply(xmlDoc(i2),'/item/element[@tag="3006,0050"]',xmlValue)        
        # organizza le coordinate in un vettore
        splittedROIPointList<-as.numeric(strsplit(ROIPointList[[1]],split = "\\\\")[[1]])
        # Separa in vettori diversi le coordinate dei punti x, y, z
        fPoint.x<-splittedROIPointList[1]
        fPoint.y<-splittedROIPointList[2]
        fPoint.z<-splittedROIPointList[3]
        
        # Calcola di quanti punti consta il "semiperimetro" della ROI
        # (questo serve per avere un'idea di dove prendere 3 punti abbstanza distanti per calcolare
        # il piano su cui giace)
        numero.Punti.semiperimetro <- as.integer((length(splittedROIPointList)/3)/2)
        
        # Se sei in presenza del numero di punti massimo, visto finora, prendi 3 punti
        if(numero.Punti.semiperimetro>numero.Punti.semiperimetro.massimo & 
           numero.Punti.semiperimetro>10) {
          
          
          # prendi i tre punti sperabilmente più "distanti"
          # E' una stima, lo sa il cielo quali siano in realtà: dovrei
          # calcolare tutte le distanze reciproche! (ma anche no...)
          # Considero il primo punto della sequenza, quello ad un quarto e quello a metà
          p1 <- 1;    p2 <- numero.Punti.semiperimetro/2;      p3 <- numero.Punti.semiperimetro;
          # Costruisci la matrice di tutti i punti (così è più facile estrarre i 3 punti)
          matrice.punti <- matrix(splittedROIPointList,ncol=3,byrow = TRUE)
          p1 <- matrice.punti[p1,];
          p2 <- matrice.punti[p2,];
          p3 <- matrice.punti[p3,];
          #           p1 <- c(splittedROIPointList[(p1)],splittedROIPointList[(p1)+1],splittedROIPointList[(p1)+2])
          #           p2 <- c(splittedROIPointList[(p2*3)+1],splittedROIPointList[(p2*3)+2],splittedROIPointList[(p2*3)+3])
          #           p3 <- c(splittedROIPointList[(p3*3)+1],splittedROIPointList[(p3*3)+2],splittedROIPointList[(p3*3)+3])
          
          # memorizza l'equazione del piano ed il numero di punti massimo della ROI
          numero.Punti.semiperimetro.massimo<-numero.Punti.semiperimetro
        }
        assegnato<-FALSE
        ReferencedSOPInstanceUID<-''
        list.index<-giveBackImageSeriesInstanceUID()
        
        if(list.index=='')  {
          logObj$handle( "error" , "giveBackImageSeriesInstanceUID(); gave back nothing"  );
        }
        # Cerca di assegnarlo ad una slice di immagine
        # Considerando un punto nella matrice della ROI calcola per ogni slice la distanza punto piano 
        # in modo che se la distanza risulta minore di 0.2 assegna la ROI a quella determinata slice
        for(slice.index in seq(1,length(dataStorage$info[[list.index]]))) {
          distanza<-obj.S$getPointPlaneDistance(c(fPoint.x,fPoint.y,fPoint.z),dataStorage$info[[list.index]][[slice.index]]$planeEquation)
          if( abs(distanza)<0.2 ) {
            ReferencedSOPInstanceUID<-dataStorage$info[[list.index]][[slice.index]]$SOPInstanceUID
            assegnato<-TRUE
          }
        }
        # Assumento le slice di immagine fra loro parallele, vedi se è su un piano parallelo ad esse
        # se 'numero.Punti.semiperimetro.massimo'>0 significa che i tre punti della ROI sono stati estratti
        if(numero.Punti.semiperimetro.massimo>0) {
          
          distanza1<-obj.S$getPointPlaneDistance(p1,dataStorage$info[[list.index]][[slice.index]]$planeEquation)
          distanza2<-obj.S$getPointPlaneDistance(p2,dataStorage$info[[list.index]][[slice.index]]$planeEquation)
          distanza3<-obj.S$getPointPlaneDistance(p3,dataStorage$info[[list.index]][[slice.index]]$planeEquation)
          # calcola il massimo gap fra i 3 punti. 
          tot <- c(distanza1,distanza2,distanza3)
          tot <- max(tot)-min(tot)
          
          # Se è maggiore di un errore indicato ( .1 mm ) dichiara non paralleli i piani
          # (in realtà punti troppo vicini potrebbero fregarmi. Speriamo di no!)
          if(tot > .1 & !(ROIName %in% attr_ROI.non.compl) ) {
            
            logObj$handle( type="warning",  msg = c("ROI ",ROIName," is not coplanar with images)")  );
            attr_ROI.non.compl<<-c(attr_ROI.non.compl,ROIName)
          }
        }
        if( assegnato == FALSE ) {
          
          logObj$handle( type="warning",  msg = c("ROI ",ROIName,": the point (",fPoint.x,",",fPoint.y,",",fPoint.z,") has no image slice! ")  );
        }
        matrice3<-rbind(matrice3,c(ROINumber,ROIName,ROIPointList,ReferencedSOPInstanceUID))
      }
      # browser()
    }
    # browser()
    return(list("IDROINameAssociation"=matrice2,"tableROIPointList"=matrice3,"FORUID.m"=FORUID.m,"RTStructSeriesInstanceUID"=RTStructSeriesInstanceUID))
    
  }   
  #=================================================================================
  # loadCTRMNRDScans.SingleSlice
  # Si occupa del caricamento di una singola slice. E' stato scorporato in quanto deve 
  # essere potenzialmente evocato da più punti del programma ( gestione cache )
  #=================================================================================   
  loadCTRMNRDScans.SingleSlice<-function( fileName, seriesInstanceUID , instanceNumber) {
    
    objServ<-services()
    
    # get the image data
    immagine<-getDICOMTag(tag = "7fe0,0010", fileName = fileName);
    # browser()
    # apply rescaleSlope and rescaleIntercept, if needed
    rescale.intercept<-as.numeric(getDICOMTag(fileName = fileName,tag ="0028,1052" )); # rescale Intercept
    rescale.slope<-as.numeric(getDICOMTag(fileName = fileName,tag ="0028,1053" )); # rescale Slope
    rescale.type<-getDICOMTag(fileName = fileName,tag ="0028,1054" ); # rescale Type   
    
    if(is.na(rescale.intercept)) rescale.intercept = 0;
    if(is.na(rescale.slope)) rescale.slope = 1;
    immagine<-immagine * rescale.slope + rescale.intercept   
    # browser()
    # Do I have to rotate the image?
    imageToBeRotated <- -1
    if( dataStorage[["info"]][[seriesInstanceUID]][[instanceNumber]][["PatientPosition"]] == "HFS" ) imageToBeRotated<-1
    if( dataStorage[["info"]][[seriesInstanceUID]][[instanceNumber]][["PatientPosition"]] == "FFS" ) imageToBeRotated<-1
    if( dataStorage[["info"]][[seriesInstanceUID]][[instanceNumber]][["PatientPosition"]] == "HFP" ) imageToBeRotated<-1
    if( dataStorage[["info"]][[seriesInstanceUID]][[instanceNumber]][["PatientPosition"]] == "FFP" ) imageToBeRotated<-1
    # browser()
    if ( imageToBeRotated != -1 ) {
      # verifica a mano la corretta rotazione: servono queste righe?????
      immagine <- objServ$rotateMatrix( immagine, imageToBeRotated )
    }
    if(SOPClassUIDList[[fileName]]$kind=="PositronEmissionTomographyImageStorage") {
      SUVCoefficient.BW<-calculate.SUVCoefficient.BW( fileName = fileName )     
      # browser()
      immagine <- SUVCoefficient.BW * immagine 
    }
    return(immagine)
  }
  #=================================================================================
  # calculate.SUVCoefficient.BW
  # Calcola il coefficiente moltiplicativo per il SUV
  #=================================================================================    
  calculate.SUVCoefficient.BW<-function(fileName) {
    AcquisitionTime<-getDICOMTag(fileName = fileName,tag ="0008,0032" );
    RadiopharmaceuticalStartTime<-getDICOMTag(fileName = fileName,tag ="0018,1072" );
    PatientWeight<-as.numeric(getDICOMTag(fileName = fileName,tag ="0010,1030" ));
    RadionuclideTotalDose<-as.numeric(getDICOMTag(fileName = fileName,tag ="0018,1074" )); # RadionuclideTotalDose
    RadionuclideHalfLife<-as.numeric(getDICOMTag(fileName = fileName,tag ="0018,1075" )); # RadionuclideHalfLife      
    rescale.type<-getDICOMTag(fileName = fileName,tag ="0028,1054" ); # rescale Type)
    rescale.type <- 'BQML'
    deltaT<-as.numeric(difftime(as.POSIXct(AcquisitionTime, format = "%H:%M:%S"),as.POSIXct(RadiopharmaceuticalStartTime, format = "%H:%M:%S"),units = 'secs'))
    deltaT<-as.numeric(difftime(as.POSIXct(AcquisitionTime, format = "%H%M%S"),as.POSIXct(str_sub(RadiopharmaceuticalStartTime,0,6), format = "%H%M%S"),units = 'secs'))
    # chiodo
    # browser()
    rescaleDueToUM<-1
    if(rescale.type=='BQML' || rescale.type=='CNTS') rescaleDueToUM<-1;
    #SUVCoefficient.BW<-PatientWeight/( RadionuclideTotalDose * exp( -deltaT *log(2)/(RadionuclideHalfLife) ) )
    SUVCoefficient.BW<-PatientWeight/( RadionuclideTotalDose * 2^( -deltaT / RadionuclideHalfLife ) ) * 1000
    SUVCoefficient.BW<-SUVCoefficient.BW*rescaleDueToUM    
    return(SUVCoefficient.BW);
  }
  #=================================================================================
  # loadCTRMRDNScans
  # Loads a DICOM CT/MR Scans
  #=================================================================================  
  loadCTRMNRDScans<-function( SOPClassUIDList , setValidCTRMNSeriesInstanceUID='any' ) {   
    imageSerie<-list()
    objServ<-services()
    if(length(SOPClassUIDList)==0) stop("SOPClassUIDList Not Built, please check DICOM SOPClasses...")
    # loop over the list    
    for(i in names(SOPClassUIDList)) {
      if(  (
        SOPClassUIDList[[i]]$kind=="CTImageStorage" ||
        SOPClassUIDList[[i]]$kind=="MRImageStorage"  ||
        SOPClassUIDList[[i]]$kind=="PositronEmissionTomographyImageStorage")  ) {
        # get the Series number
        #              seriesInstanceUID<-getDICOMTag(i,"0020,000e")
        seriesInstanceUID<-getDICOMTag(tag = "0020,000e", fileName = i)
        if(seriesInstanceUID == setValidCTRMNSeriesInstanceUID || 
           setValidCTRMNSeriesInstanceUID=='any') {
          
          # carica il FrameOfReferenceUID e cerca di capire se è uguale a quello evenutalmente
          # già caricato: se no le geometrie son troppo diverse! (e skippa la serie)
          FrameOfReferenceUID<-getDICOMTag(tag = "0020,0052", fileName = i)
          if(  is.na(attr_mainFrameOfReferenceUID) 
               || ( !is.na(attr_mainFrameOfReferenceUID) & FrameOfReferenceUID==attr_mainFrameOfReferenceUID )  ) {
            
            # se son qui significa che il FrameOfReferenceUID  è compatibile e le geometrie sono coerenti
            
            # get the Instance number (number of CT in the serie)
            instanceNumber<-as.character( getDICOMTag( tag = "0020,0013", fileName = i) )
            
            # if istance number vector is empty
            if(instanceNumber==''){
              
              instanceNumber<- SOPClassUIDList[[i]]$calculatedIstanceNumber
            }
            
            dataStorage[["info"]][[seriesInstanceUID]][[instanceNumber]]<<-list()
            # get the Patient Position
            PatientPosition<-getDICOMTag(fileName = i,tag = "0018,5100")  
            dataStorage[["info"]][[seriesInstanceUID]][[instanceNumber]][["PatientPosition"]]<<-PatientPosition
            
            # get the image data
            immagine<-loadCTRMNRDScans.SingleSlice(fileName = i, seriesInstanceUID = seriesInstanceUID, instanceNumber = instanceNumber)
            # browser()
            # now update the structure in memory
            setImageSlice( SeriesInstanceUID = seriesInstanceUID, instanceNumber = instanceNumber, imageSlice = immagine)
            
            # if the unity of measures are different, track an error!
            # ( I mean, in rescale slope/intercept and capted value in image)
            if(SOPClassUIDList[[i]]$kind=="PositronEmissionTomographyImageStorage") {
              UM<-getDICOMTag(fileName = i,tag ="0054,1001" ); # UM of voxel Cube
              CountsSource<-getDICOMTag(fileName = i,tag ="0054,1002" ); # CountsSource
              DecayCorrection<-getDICOMTag(fileName = i,tag ="0054,1102" ); # DecayCorrection
              
              Radiopharmaceutical<-getDICOMTag(fileName = i,tag ="0018,0031" ); # Radiopharmaceutical
              RadiopharmaceuticalStartTime<-getDICOMTag(fileName = i,tag ="0018,1072" ); # RadiopharmaceuticalStartTime
              RadionuclideTotalDose<-as.numeric(getDICOMTag(fileName = i,tag ="0018,1074" )); # RadionuclideTotalDose
              RadionuclideHalfLife<-as.numeric(getDICOMTag(fileName = i,tag ="0018,1075" )); # RadionuclideHalfLife
              rescale.type<-getDICOMTag(fileName = i,tag ="0028,1054" ); # rescale Type)
              # browser();
              rescale.type<-'BQML'
              # if( UM != rescale.type) {
              #   logObj$handle( "error" , "in PET image the rescale slope/intercept have different UM than the one used in the image (0054,1001) vs (0028,1054)"  );
              # }
              if(CountsSource!='EMISSION' | DecayCorrection!='START') {
                logObj$handle( "error" , c("\n ERROR: CountsSource!='EMISSION' or DecayCorrection!='START' ! This modality is not yet supported") );
              }
              if(rescale.type!='BQML' & rescale.type!='CNTS') { 
                logObj$handle( "error" , 'Rescale Type is not BQML, for the moment no other UM are supported.' );
              }  
            } else UM<-'';
            
            dataStorage[["info"]][[seriesInstanceUID]][[instanceNumber]][["SOPClassUID"]]<<-SOPClassUIDList[[i]]$kind
            dataStorage[["info"]][[seriesInstanceUID]][[instanceNumber]][["UM"]]<<-UM
            dataStorage[["info"]][[seriesInstanceUID]][[instanceNumber]][["SOPInstanceUID"]]<<-getDICOMTag( tag = "0008,0018", fileName = i)
            dataStorage[["info"]][[seriesInstanceUID]][[instanceNumber]][["FrameOfReferenceUID"]]<<-FrameOfReferenceUID
            dataStorage[["info"]][[seriesInstanceUID]][[instanceNumber]][["fileName"]]<<-i
            pixelSpacing<-splittaTAG(getDICOMTag(fileName = i,tag = "0028,0030"))
            
            ImagePositionPatient<-splittaTAG(getDICOMTag(fileName = i, tag = "0020,0032")) 
            ImageOrientationPatient<-splittaTAG(getDICOMTag(fileName = i , tag = "0020,0037"))
            iPP<-ImagePositionPatient
            iOP<-ImageOrientationPatient
            oM<-matrix(c(iOP[1],iOP[2],iOP[3],0,iOP[4],iOP[5],iOP[6],0,0,0,0,0,iPP[1],iPP[2],iPP[3],1),ncol=4); 
            
            dataStorage[["info"]][[seriesInstanceUID]][[instanceNumber]][["ImagePositionPatient"]]<<-ImagePositionPatient
            dataStorage[["info"]][[seriesInstanceUID]][[instanceNumber]][["ImageOrientationPatient"]]<<-ImageOrientationPatient
            oM[1,1]<-oM[1,1]*pixelSpacing[1]
            oM[2,1]<-oM[2,1]*pixelSpacing[1]
            oM[3,1]<-oM[3,1]*pixelSpacing[1]
            oM[1,2]<-oM[1,2]*pixelSpacing[2]
            oM[2,2]<-oM[2,2]*pixelSpacing[2]
            oM[3,2]<-oM[3,2]*pixelSpacing[2]
            dataStorage[["info"]][[seriesInstanceUID]][[instanceNumber]][["orientationMatrix"]]<<-oM
            dataStorage[["info"]][[seriesInstanceUID]][[instanceNumber]][["pixelSpacing"]]<<-pixelSpacing
            dataStorage[["info"]][[seriesInstanceUID]][[instanceNumber]][["Rows"]]<<-getDICOMTag(tag = "0028,0010", fileName = i)
            dataStorage[["info"]][[seriesInstanceUID]][[instanceNumber]][["Columns"]]<<-getDICOMTag(tag = "0028,0011", fileName = i)
            dataStorage[["info"]][[seriesInstanceUID]][[instanceNumber]][["SliceThickness"]]<<-getDICOMTag( tag = "0018,0050", fileName = i)
            if(SOPClassUIDList[[i]]$kind=="MRImageStorage") {
              dataStorage[["info"]][[seriesInstanceUID]][[instanceNumber]][["RepetitionTime"]]<<-getDICOMTag( tag = "0018,0080", fileName = i)
              dataStorage[["info"]][[seriesInstanceUID]][[instanceNumber]][["EchoTime"]]<<-getDICOMTag(tag = "0018,0081", fileName = i)
            }
            if(SOPClassUIDList[[i]]$kind=="PositronEmissionTomographyImageStorage") {
              PatientWeight<-as.numeric(getDICOMTag(fileName = i,tag ="0010,1030" ));
              AcquisitionTime<-getDICOMTag(fileName = i,tag ="0008,0032" );
              dataStorage[["info"]][[seriesInstanceUID]][[instanceNumber]][["CountsSource"]]<<-CountsSource
              dataStorage[["info"]][[seriesInstanceUID]][[instanceNumber]][["DecayCorrection"]]<<-DecayCorrection
              dataStorage[["info"]][[seriesInstanceUID]][[instanceNumber]][["Radiopharmaceutical"]]<<-Radiopharmaceutical
              dataStorage[["info"]][[seriesInstanceUID]][[instanceNumber]][["RadiopharmaceuticalStartTime"]]<<-RadiopharmaceuticalStartTime
              dataStorage[["info"]][[seriesInstanceUID]][[instanceNumber]][["AcquisitionTime"]]<<-AcquisitionTime
              dataStorage[["info"]][[seriesInstanceUID]][[instanceNumber]][["RadionuclideTotalDose"]]<<-RadionuclideTotalDose
              dataStorage[["info"]][[seriesInstanceUID]][[instanceNumber]][["RadionuclideHalfLife"]]<<-RadionuclideHalfLife
              dataStorage[["info"]][[seriesInstanceUID]][[instanceNumber]][["PatientWeight"]]<<-PatientWeight
              
              #deltaT<-as.numeric(difftime(as.POSIXct(RadiopharmaceuticalStartTime, format = "%H:%M:%S"),as.POSIXct(AcquisitionTime, format = "%H:%M:%S"),units = 'secs'))
              deltaT<-as.numeric(difftime(as.POSIXct(AcquisitionTime, format = "%H:%M:%S"),as.POSIXct(RadiopharmaceuticalStartTime, format = "%H:%M:%S"),units = 'secs'))
              if((is.na(as.POSIXct(AcquisitionTime, format = "%H:%M:%S")) & is.na(as.POSIXct(RadiopharmaceuticalStartTime, format = "%H:%M:%S",units = 'secs')))) {
                deltaT<-(as.numeric(AcquisitionTime) - as.numeric(RadiopharmaceuticalStartTime))/(60)
              }
              
              
              if(is.na(deltaT) | is.null(deltaT) | deltaT==0) {
                # browser()
                logObj$handle( "error" , "\n Error: deltaT between RadiopharmaceuticalStartTime and AcquisitionTime seems to be invalid"  );
              }
              #               rescaleDueToUM<-1
              #               #if(rescale.type=='BQML') rescaleDueToUM<-10^(-3);
              #               if(rescale.type=='BQML') rescaleDueToUM<-1;
              #               #SUVCoefficient.BW<-PatientWeight/( RadionuclideTotalDose * exp( -deltaT *log(2)/(RadionuclideHalfLife) ) )
              #               SUVCoefficient.BW<-PatientWeight/( RadionuclideTotalDose * 2^( -deltaT / RadionuclideHalfLife ) ) * 1000
              #               SUVCoefficient.BW<-SUVCoefficient.BW*rescaleDueToUM
              #               
              #               dataStorage[["info"]][[seriesInstanceUID]][[instanceNumber]][["SUVCoefficient.BW"]]<<-SUVCoefficient.BW
              #               
              dataStorage[["info"]][[seriesInstanceUID]][[instanceNumber]][["SUVCoefficient.BW"]]<<-calculate.SUVCoefficient.BW( i )
              # browser()
              
              #              dataStorage[["img"]][[seriesInstanceUID]][[instanceNumber]]<<-dataStorage[["img"]][[seriesInstanceUID]][[instanceNumber]] * SUVCoefficient.BW
              #               tmpCorrectedImg <- SUVCoefficient.BW * getImageSlice( SeriesInstanceUID = seriesInstanceUID, instanceNumber = instanceNumber) 
              #               setImageSlice(SeriesInstanceUID = seriesInstanceUID, 
              #                             instanceNumber = instanceNumber,
              #                             imageSlice = tmpCorrectedImg)              
              #               
              
            }
            # -im 
            # instanceNumber<-getDICOMTag( tag = "0020,0013", fileName = i)
            # -fm
            
            # three points to find out plane equation
            Pa<-c(oM[1,4],oM[2,4],oM[3,4])  
            Pb<-objServ$get3DPosFromNxNy(1000,0,oM)
            Pc<-objServ$get3DPosFromNxNy(0,1000,oM)
            
            abcd<-objServ$getPlaneEquationBetween3Points(Pa,Pb,Pc) 
            piano<-matrix(abcd,nrow=1)
            colnames(piano)<-c("a","b","c","d")
            dataStorage[["info"]][[seriesInstanceUID]][[instanceNumber]][["planeEquation"]]<<-piano
            # send the log
            logObj$sendLog(i)
          }
        }
      }
      
      if(  SOPClassUIDList[[i]]$kind=="RTDoseStorage") {
        
        # get the mainFrameOfReferenceUID
        FrameOfReferenceUID<-getDICOMTag( tag = "0020,0052", fileName = i)
        
        # carica il FrameOfReferenceUID e cerca di capire se è uguale a quello evenutalmente
        # già caricato: se no le geometrie son troppo diverse! (e skippa la serie)        
        if(  is.na(attr_mainFrameOfReferenceUID) 
             || ( !is.na(attr_mainFrameOfReferenceUID) & FrameOfReferenceUID==attr_mainFrameOfReferenceUID )  ) {
          
          # carica l'XML ed estraine i campi
          datiDiDose<-getDoseFromXML( i )
          
          # verifica prima di tutto se la dose caricata fa riferimento all'RTPLAN
          # presente in memoria (se c'è!)
          if( is.null(dataStorage$info$plan$SOPInstanceUID) || ( 
            !is.null(dataStorage$info$plan$SOPInstanceUID) &
            dataStorage$info$plan$SOPInstanceUID == datiDiDose$ReferencedRTPlanSequence_ReferencedSOPInstanceUID )) {
            
            if(length(dataStorage[["info"]])==0) dataStorage[["info"]]<<-list();
            if(length(dataStorage[["info"]][["doses"]])==0) dataStorage[["info"]][["doses"]]<<-list();
            SOPInstanceUID<-datiDiDose$SOPInstanceUID;
            dataStorage[["info"]][["doses"]][[SOPInstanceUID]]<<-list();
            dataStorage[["info"]][["doses"]][[SOPInstanceUID]][["imagePositionPatient"]]<<-splittaTAG(datiDiDose$ImagePositionPatient)
            dataStorage[["info"]][["doses"]][[SOPInstanceUID]][["seriesInstanceUID"]]<<-datiDiDose$SeriesInstanceUID
            dataStorage[["info"]][["doses"]][[SOPInstanceUID]][["FrameOfReferenceUID"]]<<-FrameOfReferenceUID
            dataStorage[["info"]][["doses"]][[SOPInstanceUID]][["ImageOrientationPatient"]]<<-splittaTAG(datiDiDose$ImageOrientationPatient)
            dataStorage[["info"]][["doses"]][[SOPInstanceUID]][["Rows"]]<<-datiDiDose$Rows
            dataStorage[["info"]][["doses"]][[SOPInstanceUID]][["Columns"]]<<-datiDiDose$Columns
            dataStorage[["info"]][["doses"]][[SOPInstanceUID]][["doseType"]]<<-datiDiDose$DoseType
            dataStorage[["info"]][["doses"]][[SOPInstanceUID]][["GridFrameOffsetVector"]]<<-splittaTAG(datiDiDose$GridFrameOffsetVector)
            dataStorage[["info"]][["doses"]][[SOPInstanceUID]][["DoseGridScaling"]]<<-datiDiDose$DoseGridScaling
            dataStorage[["info"]][["doses"]][[SOPInstanceUID]][["ReferencedSOPInstanceUID"]]<<-datiDiDose$ReferencedRTPlanSequence_ReferencedSOPInstanceUID 
            dataStorage[["info"]][["doses"]][[SOPInstanceUID]][["SOPClassUID"]]<<-SOPClassUIDList[[i]]$kind
            dataStorage[["info"]][["doses"]][[SOPInstanceUID]][["pixelSpacing"]]<<-splittaTAG(datiDiDose$PixelSpacing)
            dataStorage[["info"]][["DVHs"]][[SOPInstanceUID]][["DVHFromFile"]]<<-datiDiDose$DVHList
            
            # estrai l'immagine
            immagine<-getDICOMTag( tag = "7fe0,0010", fileName = i);
            if(length(dataStorage[["dose"]])==0) dataStorage[["dose"]]<<-list();
            dataStorage[["dose"]][[SOPInstanceUID]]<<-immagine * as.numeric( dataStorage[["info"]][["doses"]][[SOPInstanceUID]][["DoseGridScaling"]] )
            
            # controlli 'formali' di congruità
            # Verifica che le DOSI abbiano tutte la stessa geometria
            if(length(dataStorage[["info"]][["doses"]])>1) {
              if(!all(dataStorage[["info"]][["doses"]][[SOPInstanceUID]][["imagePositionPatient"]] %in% dataStorage[["info"]][["doses"]][[1]][["imagePositionPatient"]]))
                logObj$handle( "error" , "two RTDoses have different 'imagePositionPatient'"  );
              if(dataStorage[["info"]][["doses"]][[SOPInstanceUID]][["FrameOfReferenceUID"]]!=dataStorage[["info"]][["doses"]][[1]][["FrameOfReferenceUID"]])
                logObj$handle( "error" , "two RTDoses have different 'FrameOfReferenceUID'"  );
              if(!all(dataStorage[["info"]][["doses"]][[SOPInstanceUID]][["ImageOrientationPatient"]] %in% dataStorage[["info"]][["doses"]][[1]][["ImageOrientationPatient"]]))
                logObj$handle( "error" , "two RTDoses have different 'ImageOrientationPatient'"  );
              if(dataStorage[["info"]][["doses"]][[SOPInstanceUID]][["Rows"]]!=dataStorage[["info"]][["doses"]][[1]][["Rows"]])
                logObj$handle( "error" , "two RTDoses have different 'Rows'"  );
              if(dataStorage[["info"]][["doses"]][[SOPInstanceUID]][["Columns"]]!=dataStorage[["info"]][["doses"]][[1]][["Columns"]])
                logObj$handle( "error" , "two RTDoses have different 'Columns'"  );
              if(dataStorage[["info"]][["doses"]][[SOPInstanceUID]][["pixelSpacing"]][1]!=dataStorage[["info"]][["doses"]][[1]][["pixelSpacing"]][1])
                logObj$handle( "error" , "two RTDoses have different 'x-dim'"  );
              if(dataStorage[["info"]][["doses"]][[SOPInstanceUID]][["pixelSpacing"]][2]!=dataStorage[["info"]][["doses"]][[1]][["pixelSpacing"]][2])
                logObj$handle( "error" , "two RTDoses have different 'y-dim'"  );
              if(!all(dataStorage[["info"]][["doses"]][[SOPInstanceUID]][["GridFrameOffsetVector"]] %in% dataStorage[["info"]][["doses"]][[1]][["GridFrameOffsetVector"]]))
                logObj$handle( "error" , "two RTDoses have different 'GridFrameOffsetVector'"  );
              if(dataStorage[["info"]][["doses"]][[SOPInstanceUID]][["doseType"]]!="PHYSICAL")
                logObj$handle( "error" , "rhe dose is not 'PHYSICAL' (attribute 'DoseType')"  );
            }
          }
        }
      }      
    }
  }  
  #=================================================================================
  # NAME: getROIVoxels
  # restituisce i voxel interni ad una data ROI
  #=================================================================================    
  getROIVoxels<-function( Structure = Structure , new.pixelSpacing=c() ) {
    objS<-services();
    
    # cerca nella cache di memoria, se attivata, la presenza della ROI già estratta
    if(attr_ROIVoxelMemoryCache==TRUE  & !is.null(attr_ROIVoxelMemoryCacheArray[[Structure]]) & length(new.pixelSpacing)==0) {
      return(attr_ROIVoxelMemoryCacheArray[[Structure]]);
    }
    
    if(!(Structure %in% getROIList()[2,])) logObj$handle( "error" , c( Structure," not present."  )  );
    
    # try to find out which Series is the CT/MR serie
    SeriesInstanceUID<-giveBackImageSeriesInstanceUID()
    if(SeriesInstanceUID == '' ) {
      logObj$handle( "error" , "missing CT/MR series"  );
    }
    # browser()
    res<-getROIVoxelsFromCTRMN( Structure = Structure, SeriesInstanceUID = SeriesInstanceUID,
                                new.pixelSpacing = new.pixelSpacing)
    # browser()
    croppedRes<-list()
    croppedRes$DOM<-res$DOM
    croppedRes$geometricalInformationOfImages<-res$geometricalInformationOfImages
    # browser()
    # image(res$masked.images[,,25])
    croppedRes$masked.images<-objS$cropCube( bigCube = res$masked.images)
    croppedRes$masked.images$location$fe<-dim(res$masked.images)[1]
    croppedRes$masked.images$location$se<-dim(res$masked.images)[2]
    croppedRes$masked.images$location$te<-dim(res$masked.images)[3]
    croppedRes$geometricalInformationOfImages$koc<-"littleCube"
    croppedRes$resamplingInformation<-res$resamplingInformation
    # croppedRes$originalUnCropped<-res$immagineNonMascherata
    class(croppedRes)<-"geoLetStructureVoxelList"
    
    # se la cache in memoria è attiva, salvane una copia
    if(attr_ROIVoxelMemoryCache==TRUE  & length(new.pixelSpacing)==0) attr_ROIVoxelMemoryCacheArray[[Structure]]<<-croppedRes;
    invisible(gc())
    return( croppedRes )
  }  
  #=================================================================================
  # NAME: getEmptyStructure
  #=================================================================================    
  getEmptyStructure<-function(  ) {
    objS<-services();
    
    # # cerca nella cache di memoria, se attivata, la presenza della ROI già estratta
    # if(attr_ROIVoxelMemoryCache==TRUE  & !is.null(attr_ROIVoxelMemoryCacheArray[[Structure]]) & length(new.pixelSpacing)==0) {
    #   return(attr_ROIVoxelMemoryCacheArray[[Structure]]);
    # }
    # 
    # if(!(Structure %in% getROIList()[2,])) logObj$handle( "error" , c( Structure," not present."  )  );
    
    # try to find out which Series is the CT/MR serie
    SeriesInstanceUID<-giveBackImageSeriesInstanceUID()
    if(SeriesInstanceUID == '' ) {
      logObj$handle( "error" , "missing CT/MR series"  );
    }
    # browser()
    res<-getROIVoxelsFromCTRMN( Structure = Structure, SeriesInstanceUID = SeriesInstanceUID,
                                new.pixelSpacing = new.pixelSpacing)
    # browser()
    croppedRes<-list()
    croppedRes$DOM<-res$DOM
    croppedRes$geometricalInformationOfImages<-res$geometricalInformationOfImages
    # browser()
    # image(res$masked.images[,,25])
    croppedRes$masked.images<-objS$cropCube( bigCube = res$masked.images)
    croppedRes$masked.images$location$fe<-dim(res$masked.images)[1]
    croppedRes$masked.images$location$se<-dim(res$masked.images)[2]
    croppedRes$masked.images$location$te<-dim(res$masked.images)[3]
    croppedRes$geometricalInformationOfImages$koc<-"littleCube"
    croppedRes$resamplingInformation<-res$resamplingInformation
    # croppedRes$originalUnCropped<-res$immagineNonMascherata
    class(croppedRes)<-"geoLetStructureVoxelList"
    
    # se la cache in memoria è attiva, salvane una copia
    # if(attr_ROIVoxelMemoryCache==TRUE  & length(new.pixelSpacing)==0) attr_ROIVoxelMemoryCacheArray[[Structure]]<<-croppedRes;
    # invisible(gc())
    return( croppedRes )
  }   
  #=================================================================================
  # getGeometricalInformationOfImage
  # give back pixelspacing and other little stuff about the CT/MR
  #=================================================================================  
  getGeometricalInformationOfImage<-function() {
    
    serieInstanceUID<-giveBackImageSeriesInstanceUID();
    ddd<-dataStorage$info[[ serieInstanceUID ]][[1]]
    return(list(
      "PatientPosition"=ddd$PatientPosition,
      "SOPClassUID"=ddd$MRImageStorage,
      "pixelSpacing"=ddd$pixelSpacing,
      "ImagePositionPatient"=ddd$ImagePositionPatient,
      "Rows"=ddd$Rows,
      "Columns"=ddd$Columns,
      "SliceThickness"=ddd$SliceThickness,
      "supposedNumberOfSlices"=length(dataStorage$info[[ serieInstanceUID ]]),
      "randomSliceImageOrientationPatient"=ddd$ImageOrientationPatient,
      "randomSlicePlaneEquation"=ddd$planeEquation
    ))
  }  
  #=================================================================================
  # NAME: getImageSlice
  # E' un wrapper della dataStorage per restituire la slice dell'immagine. E' necessaria
  # questa funzione di wrap per poter gestire la cache in maniera indolore ed efficace.
  #=================================================================================    
  getImageSlice<-function( SeriesInstanceUID, instanceNumber, typeOfImage="CT/MRI/PET"  ) {
    if(typeOfImage=="CT/MRI/PET" & attr_loadRAWInCache == TRUE)   return( dataStorage$img[[SeriesInstanceUID]][[instanceNumber]]  );
    if(typeOfImage=="CT/MRI/PET" & attr_loadRAWInCache == FALSE) {
      fileName<-dataStorage$info[[SeriesInstanceUID]][[instanceNumber]]$fileName
      immagine<-loadCTRMNRDScans.SingleSlice( fileName  = fileName, seriesInstanceUID = SeriesInstanceUID , instanceNumber = instanceNumber)
      return(immagine);
    }
  }
  setImageSlice<-function( SeriesInstanceUID, instanceNumber, imageSlice, typeOfImage="CT/MRI/PET"  ) {
    if(length(dataStorage$img)==0) dataStorage$img<<-list()
    if(length(dataStorage$img[[SeriesInstanceUID]])==0) dataStorage$img[[SeriesInstanceUID]]<<-list()
    if(typeOfImage=="CT/MRI/PET" & attr_loadRAWInCache == TRUE ) dataStorage$img[[SeriesInstanceUID]][[instanceNumber]]<<-imageSlice
    if(typeOfImage=="CT/MRI/PET" & attr_loadRAWInCache == FALSE ) dataStorage$img[[SeriesInstanceUID]][[instanceNumber]]<<-NA
  }  
  #=================================================================================
  # NAME: getROIVoxelsFromCTRMN
  # Estrae i voxel da scansioni CT,MR
  #=================================================================================   
  getROIVoxelsFromCTRMN<-function( Structure = Structure, SeriesInstanceUID = SeriesInstanceUID,
                                   new.pixelSpacing=c() ) {
    objService<-services()  
    
    if ( (Structure %in% getROIList()[2,]) == FALSE )  return(NA)
    # define some variables to make more clear the code
    numberOfRows<-as.numeric(dataStorage$info[[SeriesInstanceUID]][[1]]$Rows);
    #numberOfRows<-as.numeric(dataStorage$info[[1]][[1]]$Rows);
    numberOfColumns<-as.numeric(dataStorage$info[[SeriesInstanceUID]][[1]]$Columns);
    #numberOfRows<-as.numeric(dataStorage$info[[1]][[1]]$Columns);
    numberOfSlices<-length(dataStorage$img[[SeriesInstanceUID]]);
    
    old.ps <- getPixelSpacing()
    if( length(new.pixelSpacing)>0 ) {
      if(new.pixelSpacing[1] == old.ps[1] & new.pixelSpacing[2]==old.ps[2]) new.pixelSpacing<-c()
    }
    if(  length(new.pixelSpacing)>0 ){
      new.pixelSpacing <- c(new.pixelSpacing,old.ps[3])
      ratio.ps.x <- old.ps[1] / new.pixelSpacing[1] 
      ratio.ps.y <- old.ps[2] / new.pixelSpacing[2]
      # new.numberOfRows <- numberOfRows * delta.ps.x
      # new.numberOfColumns <- numberOfColumns * delta.ps.y
    } else {
      new.numberOfRows <- numberOfRows
      new.numberOfColumns <- numberOfColumns
    }
    
    # initialize the image array with the right dimension
    # -im
    # image.arr<-array( data = -1, dim = c(numberOfRows, numberOfColumns, numberOfSlices ) )
    # image.arr<-array( data = -1, dim = c(numberOfColumns, numberOfRows, numberOfSlices ) )
    image.arr<-array( data = -1, dim = c(numberOfColumns, numberOfRows, numberOfSlices ) )
    # -fm
    
    # index values listed as characters: creates empty array of DICOM orientation matrices
    index<-as.character(sort(as.numeric(names( dataStorage$img[[SeriesInstanceUID]]) )))  
    
    # create and fill the vector DOM, of DICOM orientation matrices
    DOM<-c();  nn<-0
    for (n in index) {
      nn<-nn+1
      #      image.arr[,,nn]<-dataStorage$img[[SeriesInstanceUID]][[n]]   
      image.arr[,,nn] <- getImageSlice( SeriesInstanceUID = SeriesInstanceUID, instanceNumber = n)
    }  
    
    if(all(new.pixelSpacing==old.ps)==FALSE) {
      # ooo <- image.arr
      image.arr<-objService$new.trilinearInterpolator(voxelCube = image.arr,
                                                      pixelSpacing.new = new.pixelSpacing,
                                                      pixelSpacing.old = old.ps)
      
      new.pixelSpacing <- c( numberOfColumns/dim(image.arr)[1] * old.ps[1], 
                             numberOfRows/dim(image.arr)[2] * old.ps[2], 
                             new.pixelSpacing[3])      
      
      # new.pixelSpacing <- c( numberOfColumns/dim(image.arr)[1], numberOfRows/dim(image.arr)[2] , new.pixelSpacing[3])
      new.numberOfColumns <- dim(image.arr)[1]
      new.numberOfRows <- dim(image.arr)[2]
      resampled <- TRUE
    } else { resampled <- FALSE  }
    
    for (n in index) {
      # browser()
      Dic.Or.Mat <- dataStorage$info[[SeriesInstanceUID]][[n]]$orientationMatrix[c(1:3,5:7,13:15)]
      if(all(new.pixelSpacing==old.ps)==FALSE) {
        Dic.Or.Mat[1:3] <- Dic.Or.Mat[1:3]* (new.pixelSpacing[1]/old.ps[1])
        Dic.Or.Mat[4:6] <- Dic.Or.Mat[4:6]* (new.pixelSpacing[2]/old.ps[2])
      }
      DOM<-c(DOM, Dic.Or.Mat)
    }  
    
    # fills the vectors of X and Y coordinates 
    # and other Vectors 'associatedInstanceNumberVect' and 'arrayInstanceNumberWithROI'
    TotalX<- -10000;  TotalY<- -10000;  arrayAssociationROIandSlice<- -10000;
    OriginX<- -10000;   OriginY<- -10000
    associatedInstanceNumberVect<- -10000
    contatoreROI<-1; indiceDOM<-1;
    
    # for each instance number
    for (n in index) {
      # check if there is a ROI for such slice
      # browser()
      for (m in which(dataStorage$info[[SeriesInstanceUID]][[n]]$ROIList[,1]==Structure)) {
        # browser()
        # find the slice and gets the key for accessing at coordinates vectors
        #key<-dataStorage$info[[SeriesInstanceUID]][[n]]$ROIList[m,2]
        key<-dataStorage$info[[SeriesInstanceUID]][[n]]$ROIList[m,4]
        
        # calculate how many ROIs are co-planar
        numeroROIComplanari<-length(dataStorage$structures[[Structure]][[key]])
        # for each one of them concat the array
        for(indiceROI in seq(1,numeroROIComplanari)) {    
          #          browser();
          TotalX<-c(TotalX, dataStorage$structures[[Structure]][[key]][[indiceROI]][,1])
          TotalY<-c(TotalY, dataStorage$structures[[Structure]][[key]][[indiceROI]][,2])
          
          # calculate how many points compose the ROI
          numeroPunti<-length(dataStorage$structures[[Structure]][[key]][[indiceROI]][,1])
          
          # for each point write which is the related Slice in the cube-matrix
          arrayAssociationROIandSlice<-c(arrayAssociationROIandSlice,rep(   which( index == n ) -1  , numeroPunti ))
          
          # Usa OriginX and OriginY as terminator
          TotalX<-c(TotalX, OriginX)
          TotalY<-c(TotalY, OriginY)      
          arrayAssociationROIandSlice<-c(arrayAssociationROIandSlice,OriginX)
          
          contatoreROI<-contatoreROI+1      
        } 
      }      
      indiceDOM<-indiceDOM+1;
    }
    
    # ok, call the Wrapper!
    final.array<-NewMultiPointInPolyObl(
      # array of DICOM Orientation Matrices
      DICOMOrientationVector = DOM, 
      # X and Y vector Points
      totalX = TotalX, totalY = TotalY, 
      # association between ROIs and Slices in the 3D Matrix
      arrayAssociationROIandSlice = arrayAssociationROIandSlice,
      # matrices dimensions (rows and columns)
      nX = new.numberOfColumns, 
      nY = new.numberOfRows,
      nZ = numberOfSlices
      # nX = numberOfColumns, 
      # nY = numberOfRows,
      # nZ = numberOfSlices
    )
    # browser();
    final.array<-array(data = final.array, dim = c(   new.numberOfColumns, new.numberOfRows, numberOfSlices )   )
    # print( dim(final.array))
    # In Example: 
    #
    # > TotalX[0:70]
    # [1] -10000.00     -5.16     -3.28     -1.41      0.47      2.34      4.22      5.12      6.09      7.61      7.97      9.48      9.84     10.96     11.72
    # [16]     12.14     12.93     13.59     13.61     14.19     14.70     14.85     14.23     13.59     13.45     12.68     11.72     11.11      9.84      8.78
    # [31]      7.97      6.09      4.22      2.84      2.34      0.85      0.47      0.14     -1.41     -3.28     -3.92     -5.16     -5.73     -6.56     -7.03
    # [46]     -8.91    -10.26    -10.78    -11.25    -11.86    -12.01    -11.83    -11.29    -10.78    -10.64     -9.88     -8.91     -8.81     -7.51     -7.03
    # [61]     -5.46     -5.16 -10000.00     -5.16     -3.28     -1.41      0.47      2.34      4.22      5.38
    # 
    # > arrayAssociationROIandSlice[0:70]
    # [1] -10000      8      8      8      8      8      8      8      8      8      8      8      8      8      8      8      8      8      8      8      8
    # [22]      8      8      8      8      8      8      8      8      8      8      8      8      8      8      8      8      8      8      8      8      8
    # [43]      8      8      8      8      8      8      8      8      8      8      8      8      8      8      8      8      8      8      8      9 -10000
    # [64]      9      9      9      9      9      9      9
    
    # ROTATE THE MATRIX
    
    for ( i in seq(1,dim(image.arr)[3] )) {
      #      image.arr[,,i]<-objService$SV.rotateMatrix(image.arr[,,i])
      final.array[,,i]<-t(objService$rotateMatrix(final.array[,,i],rotations=3))
    }
    
    # -im
    # vedi se e come è necessario CAPPOTTARE il voxelCube dell'immagine per allinearlo alla ROI
    Dic.Or.Mat.1 <- dataStorage$info[[SeriesInstanceUID]][[1]]$orientationMatrix
    Dic.Or.Mat.last <- dataStorage$info[[SeriesInstanceUID]][[numberOfSlices]]$orientationMatrix
    if(all(new.pixelSpacing==old.ps)==FALSE) {
      Dic.Or.Mat.1[1:3,1] <- Dic.Or.Mat.1[1:3,1]* (new.pixelSpacing[1]/old.ps[1])
      Dic.Or.Mat.1[1:3,2] <- Dic.Or.Mat.1[1:3,2]* (new.pixelSpacing[2]/old.ps[2])
    }
    
    initial.point <- objService$get3DPosFromNxNy(Nx = 1,Ny = 1,oM = Dic.Or.Mat.1)
    final.point <- objService$get3DPosFromNxNy(Nx = 100,Ny = 100,oM = Dic.Or.Mat.last)
    directions <- final.point - initial.point
    
    # cat("\nDirections: ",directions[1:3],"\n")
    # Sarosh E' correttamente orientato per <+,+,+>'
    # browser(directions)
    # Se devo ribaltare lungo la Z
    # if(directions[3] < 0) {
    #   # browser()
    #   # final.array[,,i]<-t(objService$rotateMatrix(final.array[,,i],rotations=1))
    # }
    # if(directions[1] < 0) {
    #   cat("\n ------------------------------------------------------")
    #   cat("\n ATTENZIONE: orientamento supportato ma non ancora testato")
    #   cat("\n si prega di contattare il manutentore del pacchetto")
    #   cat("\n ------------------------------------------------------")
    #   # stop()
    # }
    # if(directions[2] < 0) {
    #   cat("\n ------------------------------------------------------")
    #   cat("\n ATTENZIONE: orientamento supportato ma non ancora testato")
    #   cat("\n si prega di contattare il manutentore del pacchetto")
    #   cat("\n ------------------------------------------------------")
    #   # stop()
    # }    
    # -fm      
    
    #    return(list(TotalX=TotalX, TotalY=TotalY, FullZ=FullZ, Offset=Offset, 
    #                DOM=array(DOM, dim = c(3,3,length(index))), final.array=final.array, masked.images=final.array*image.arr))
    
    # browser()
    # immagineNonMascherata <- final.array
    immagineMascherata<- array(NA, dim=c(  dim(image.arr)[1],dim(image.arr)[2],dim(image.arr)[3]  ))
    immagineMascherata[which(final.array==1,arr.ind = TRUE)]<-image.arr[ which(final.array==1,arr.ind = TRUE) ]
    # image(immagineMascherata[,,25],col = grey.colors(258))
    # image(immagineMascherata[,,25])
    return(list(
      "DOM"=array(DOM, dim = c(3,3,length(index))), 
      "final.array"=final.array, 
      "masked.images"=immagineMascherata,
      "geometricalInformationOfImages"=getGeometricalInformationOfImage(),
      # "immagineNonMascherata"=immagineNonMascherata,
      "resamplingInformation"=list(
        "px" = new.pixelSpacing[1], "py" = new.pixelSpacing[1], "pz" = new.pixelSpacing[3], "resampled" = resampled, 
        "numberOfRows" = new.numberOfRows, "numberOfColumns" = new.numberOfColumns)
    )
    )
  }  
  old.getROIVoxelsFromCTRMN<-function( Structure = Structure, SeriesInstanceUID = SeriesInstanceUID,
                                       new.pixelSpacing=c() ) {
    objService<-services()  
    
    if ( (Structure %in% getROIList()[2,]) == FALSE )  return(NA)
    # define some variables to make more clear the code
    numberOfRows<-as.numeric(dataStorage$info[[SeriesInstanceUID]][[1]]$Rows);
    #numberOfRows<-as.numeric(dataStorage$info[[1]][[1]]$Rows);
    numberOfColumns<-as.numeric(dataStorage$info[[SeriesInstanceUID]][[1]]$Columns);
    #numberOfRows<-as.numeric(dataStorage$info[[1]][[1]]$Columns);
    numberOfSlices<-length(dataStorage$img[[SeriesInstanceUID]]);
    
    # initialize the image array with the right dimension
    # -im
    # image.arr<-array( data = -1, dim = c(numberOfRows, numberOfColumns, numberOfSlices ) )
    image.arr<-array( data = -1, dim = c(numberOfColumns, numberOfRows, numberOfSlices ) )
    # -fm
    
    # index values listed as characters: creates empty array of DICOM orientation matrices
    index<-as.character(sort(as.numeric(names( dataStorage$img[[SeriesInstanceUID]]) )))  
    
    # create and fill the vector DOM, of DICOM orientation matrices
    DOM<-c();  nn<-0
    for (n in index) {
      DOM<-c(DOM, dataStorage$info[[SeriesInstanceUID]][[n]]$orientationMatrix[c(1:3,5:7,13:15)])
      nn<-nn+1
      #      image.arr[,,nn]<-dataStorage$img[[SeriesInstanceUID]][[n]]   
      image.arr[,,nn] <- getImageSlice( SeriesInstanceUID = SeriesInstanceUID, instanceNumber = n)
    }  
    
    # fills the vectors of X and Y coordinates 
    # and other Vectors 'associatedInstanceNumberVect' and 'arrayInstanceNumberWithROI'
    TotalX<- -10000;  TotalY<- -10000;  arrayAssociationROIandSlice<- -10000;
    OriginX<- -10000;   OriginY<- -10000
    associatedInstanceNumberVect<- -10000
    contatoreROI<-1; indiceDOM<-1;
    
    # for each instance number
    for (n in index) {
      # check if there is a ROI for such slice
      # browser()
      for (m in which(dataStorage$info[[SeriesInstanceUID]][[n]]$ROIList[,1]==Structure)) {
        # browser()
        # find the slice and gets the key for accessing at coordinates vectors
        #key<-dataStorage$info[[SeriesInstanceUID]][[n]]$ROIList[m,2]
        key<-dataStorage$info[[SeriesInstanceUID]][[n]]$ROIList[m,4]
        
        # calculate how many ROIs are co-planar
        numeroROIComplanari<-length(dataStorage$structures[[Structure]][[key]])
        # for each one of them concat the array
        for(indiceROI in seq(1,numeroROIComplanari)) {    
          #          browser();
          TotalX<-c(TotalX, dataStorage$structures[[Structure]][[key]][[indiceROI]][,1])
          TotalY<-c(TotalY, dataStorage$structures[[Structure]][[key]][[indiceROI]][,2])
          
          # calculate how many points compose the ROI
          numeroPunti<-length(dataStorage$structures[[Structure]][[key]][[indiceROI]][,1])
          
          # for each point write which is the related Slice in the cube-matrix
          arrayAssociationROIandSlice<-c(arrayAssociationROIandSlice,rep(   which( index == n ) -1  , numeroPunti ))
          
          # Usa OriginX and OriginY as terminator
          TotalX<-c(TotalX, OriginX)
          TotalY<-c(TotalY, OriginY)      
          arrayAssociationROIandSlice<-c(arrayAssociationROIandSlice,OriginX)
          
          contatoreROI<-contatoreROI+1      
        } 
      }      
      indiceDOM<-indiceDOM+1;
    }
    # browser();
    # ok, call the Wrapper!
    final.array<-NewMultiPointInPolyObl(
      # array of DICOM Orientation Matrices
      DICOMOrientationVector = DOM, 
      # X and Y vector Points
      totalX = TotalX, totalY = TotalY, 
      # association between ROIs and Slices in the 3D Matrix
      arrayAssociationROIandSlice = arrayAssociationROIandSlice,
      # matrices dimensions (rows and columns)
      nX = numberOfColumns, 
      nY = numberOfRows,
      nZ = numberOfSlices
    )
    
    final.array<-array(data = final.array, dim = c(   numberOfColumns, numberOfRows, numberOfSlices )   )
    # In Example: 
    #
    # > TotalX[0:70]
    # [1] -10000.00     -5.16     -3.28     -1.41      0.47      2.34      4.22      5.12      6.09      7.61      7.97      9.48      9.84     10.96     11.72
    # [16]     12.14     12.93     13.59     13.61     14.19     14.70     14.85     14.23     13.59     13.45     12.68     11.72     11.11      9.84      8.78
    # [31]      7.97      6.09      4.22      2.84      2.34      0.85      0.47      0.14     -1.41     -3.28     -3.92     -5.16     -5.73     -6.56     -7.03
    # [46]     -8.91    -10.26    -10.78    -11.25    -11.86    -12.01    -11.83    -11.29    -10.78    -10.64     -9.88     -8.91     -8.81     -7.51     -7.03
    # [61]     -5.46     -5.16 -10000.00     -5.16     -3.28     -1.41      0.47      2.34      4.22      5.38
    # 
    # > arrayAssociationROIandSlice[0:70]
    # [1] -10000      8      8      8      8      8      8      8      8      8      8      8      8      8      8      8      8      8      8      8      8
    # [22]      8      8      8      8      8      8      8      8      8      8      8      8      8      8      8      8      8      8      8      8      8
    # [43]      8      8      8      8      8      8      8      8      8      8      8      8      8      8      8      8      8      8      8      9 -10000
    # [64]      9      9      9      9      9      9      9
    
    # ROTATE THE MATRIX
    
    for ( i in seq(1,dim(image.arr)[3] )) {
      #      image.arr[,,i]<-objService$SV.rotateMatrix(image.arr[,,i])
      final.array[,,i]<-t(objService$rotateMatrix(final.array[,,i],rotations=3))
    }
    
    # -im
    # vedi se e come è necessario CAPPOTTARE il voxelCube dell'immagine per allinearlo alla ROI
    
    initial.point <- objService$get3DPosFromNxNy(Nx = 1,Ny = 1,oM = dataStorage$info[[SeriesInstanceUID]][[1]]$orientationMatrix)
    final.point <- objService$get3DPosFromNxNy(Nx = 100,Ny = 100,oM = dataStorage$info[[SeriesInstanceUID]][[numberOfSlices]]$orientationMatrix)
    directions <- final.point - initial.point
    
    cat("\nDirections: ",directions[1:3],"\n\n")
    # Sarosh E' correttamente orientato per <+,+,+>'
    # browser(directions)
    # Se devo ribaltare lungo la Z
    # if(directions[3] < 0) {
    #   # browser()
    #   # final.array[,,i]<-t(objService$rotateMatrix(final.array[,,i],rotations=1))
    # }
    if(directions[1] < 0) {
      cat("\n ------------------------------------------------------")
      cat("\n ATTENZIONE: orientamento supportato ma non ancora testato")
      cat("\n si prega di contattare il manutentore del pacchetto")
      cat("\n ------------------------------------------------------")
      # stop()
    }
    if(directions[2] < 0) {
      cat("\n ------------------------------------------------------")
      cat("\n ATTENZIONE: orientamento supportato ma non ancora testato")
      cat("\n si prega di contattare il manutentore del pacchetto")
      cat("\n ------------------------------------------------------")
      # stop()
    }    
    # -fm      
    
    #    return(list(TotalX=TotalX, TotalY=TotalY, FullZ=FullZ, Offset=Offset, 
    #                DOM=array(DOM, dim = c(3,3,length(index))), final.array=final.array, masked.images=final.array*image.arr))
    
    # browser()
    immagineMascherata<- array(NA, dim=c(  dim(image.arr)[1],dim(image.arr)[2],dim(image.arr)[3]  ))
    immagineMascherata[which(final.array==1,arr.ind = TRUE)]<-image.arr[ which(final.array==1,arr.ind = TRUE) ]
    
    return(list(
      "DOM"=array(DOM, dim = c(3,3,length(index))), 
      "final.array"=final.array, 
      "masked.images"=immagineMascherata,
      "geometricalInformationOfImages"=getGeometricalInformationOfImage()
    )
    )
  }   
  NewMultiPointInPolyObl<-function(DICOMOrientationVector,totalX,totalY,arrayAssociationROIandSlice,nX,nY,nZ ) {  
    
    maxX<-max(totalX)
    minX<-min(totalX[which(totalX>-10000)])
    maxY<-max(totalY)
    minY<-min(totalY[which(totalY>-10000)])
    
    # creates the PIPvector
    PIPvector<-rep.int(x = 0, times = nX * nY * nZ)  
    numberOfPoints<-length(totalX);
    result<-.C("NewMultiPIPObl", 
               as.integer(PIPvector), as.double(totalX), as.double(totalY), as.integer(numberOfPoints), 
               as.integer(nX), as.integer(nY), as.integer(nZ),             
               as.integer(arrayAssociationROIandSlice), 
               as.double(DICOMOrientationVector),as.double(minX),as.double(maxX),as.double(minY),as.double(maxY))  
    return(result[[1]])
  }  
  #=================================================================================
  # getImageVoxelCube
  # give back the greyLevel voxel cube. If no ps.x/y/z are specified it gives back 
  # the voxelCube of the original dimensions, otherwise it gives back the interpolated
  # voxelCube according to the wished pixelSpacing along x,y or z
  #=================================================================================     
  getImageVoxelCube<-function( ps.x=NA, ps.y=NA, ps.z=NA) {
    objS<-services();
    # prendi il cubone
    voxelCube<-createImageVoxelCube()
    # se non  server interpolare
    if(is.na(ps.x) && is.na(ps.y) && is.na(ps.z) ) return(voxelCube)
    
    # se invece serve interpolare: prendi i pixelSpacing lungo la X, la Y e la Z (slice thickness)
    oldPixelSpacing<-getPixelSpacing();
    
    if(is.na(ps.x))  ps.x <- oldPixelSpacing[1];
    if(is.na(ps.y))  ps.y <- oldPixelSpacing[2];
    if(is.na(ps.z))  ps.z <- oldPixelSpacing[3];
    
    voxelCube<-objS$new.trilinearInterpolator(
      voxelCube = voxelCube,
      pixelSpacing.new = c(ps.x,ps.y,ps.z),
      pixelSpacing.old = oldPixelSpacing )    
    
    invisible(gc())
    
    return( voxelCube )
  }  
  #=================================================================================
  # getPixelSpacing
  # una funzione specifica per il pixelSpacing, visto quanto è usato!
  #================================================================================= 
  getPixelSpacing<-function( seriesInstanceUID  = NA) {
    if(is.na(seriesInstanceUID))  seriesInstanceUID<-giveBackImageSeriesInstanceUID();   
    ps<-c();
    ps[1]<-as.numeric(dataStorage$info[[seriesInstanceUID]][[1]]$pixelSpacing[1])
    ps[2]<-as.numeric(dataStorage$info[[seriesInstanceUID]][[1]]$pixelSpacing[2])
    ps[3]<-as.numeric(dataStorage$info[[seriesInstanceUID]][[1]]$SliceThickness)
    return(ps);
  }   
  #=================================================================================
  # createImageVoxelCube
  # create the imageVoxelCube for the current obj and for the image stored
  #=================================================================================   
  createImageVoxelCube<-function() {
    # ge the series Instance UID of the images
    seriesInstanceUID<-giveBackImageSeriesInstanceUID();
    # order them according with the Instance Number
    listaSeqImages<-as.character(sort(as.numeric(names( dataStorage$img[[seriesInstanceUID]]) )))
    # get dimensional data
    Rows<-dataStorage$info[[seriesInstanceUID]][[1]]$Rows
    Columns<-dataStorage$info[[seriesInstanceUID]][[1]]$Columns
    Slices<-length(listaSeqImages)
    
    cubone<-array(data = 0,dim = c(Columns,Rows,Slices))
    numSlice<-1
    # add the slices and build the cube
    for(i in listaSeqImages) {
      #      cubone[,,numSlice]<-dataStorage$img[[seriesInstanceUID]][[as.character(i)]]
      cubone[,,numSlice]<-getImageSlice( SeriesInstanceUID = seriesInstanceUID, instanceNumber = as.character(i))
      numSlice<-numSlice+1
    }
    return(cubone)
  }    
  #=================================================================================
  # getPlanFromXML
  # get the interesting tag via XML instead of via normal dump (more robust)
  #=================================================================================   
  getDoseFromXML<-function(fileName) { 
    obj.S<-services();
    # browser()
    if(!file.exists(fileName)) logObj$handle( "error" , "The indicated file does not exist ( geoLet::getDoseFromXML() )"  );
    
    # build the XML file and get the XML structure
    doc<-getXMLStructureFromDICOMFile(fileName = fileName, folderCleanUp = folderCleanUp)
    
    ImagePositionPatient<-xpathApply(doc,'/file-format/data-set/element[@tag="0020,0032" and @name="ImagePositionPatient"]',xmlValue)[[1]]
    ImageOrientationPatient<-xpathApply(doc,'/file-format/data-set/element[@tag="0020,0037" and @name="ImageOrientationPatient"]',xmlValue)[[1]]
    Rows<-xpathApply(doc,'/file-format/data-set/element[@tag="0028,0010" and @name="Rows"]',xmlValue)[[1]]
    Columns<-xpathApply(doc,'/file-format/data-set/element[@tag="0028,0011" and @name="Columns"]',xmlValue)[[1]]
    PixelSpacing<-xpathApply(doc,'/file-format/data-set/element[@tag="0028,0030" and @name="PixelSpacing"]',xmlValue)[[1]]
    PixelRepresentation<-xpathApply(doc,'/file-format/data-set/element[@tag="0028,0103" and @name="PixelRepresentation"]',xmlValue)[[1]]
    SOPInstanceUID<-xpathApply(doc,'/file-format/data-set/element[@tag="0008,0018" and @name="SOPInstanceUID"]',xmlValue)[[1]]
    SOPInstanceUID<-xpathApply(doc,'/file-format/data-set/element[@tag="0008,0018" and @name="SOPInstanceUID"]',xmlValue)[[1]]
    SeriesInstanceUID<-xpathApply(doc,'/file-format/data-set/element[@tag="0020,000e" and @name="SeriesInstanceUID"]',xmlValue)[[1]]
    GridFrameOffsetVector<-xpathApply(doc,'/file-format/data-set/element[@tag="3004,000c" and @name="GridFrameOffsetVector"]',xmlValue)[[1]]
    DoseUnits<-xpathApply(doc,'/file-format/data-set/element[@tag="3004,0002" and @name="DoseUnits"]',xmlValue)[[1]]
    DoseType<-xpathApply(doc,'/file-format/data-set/element[@tag="3004,0004" and @name="DoseType"]',xmlValue)[[1]]
    DoseGridScaling<-xpathApply(doc,'/file-format/data-set/element[@tag="3004,000e" and @name="DoseGridScaling"]',xmlValue)[[1]]
    ReferencedRTPlanSequence_ReferencedSOPInstanceUID<-xpathApply(doc,'/file-format/data-set/sequence[@tag="300c,0002" and @name="ReferencedRTPlanSequence"]//element[@tag="0008,1155" and @name="ReferencedSOPInstanceUID"]',xmlValue)[[1]]    
    
    # now look for some DVHs
    #a<-xpathApply(doc,'/file-format/data-set/sequence[@tag="3004,0050" and @name="DVHSequence"]',xmlValue)
    a<-getNodeSet(doc,'/file-format/data-set/sequence[@tag="3004,0050" and @name="DVHSequence"]/item')
    
    DVHList<-list();
    ct<-1
    if(length(a)>0) {
      for(ct in seq(1,length(a))) {
        ReferencedROINumber<-xpathApply(a[[ct]],'//item/sequence[@tag="3004,0060"]//element[@tag="3006,0084"]',xmlValue)[[ct]]
        
        ROIName<-as.character(ReferencedROINumber);
        
        DVHList[[ROIName]]<-list();
        
        DVHList[[ROIName]][["DVHType"]]<-xpathApply(a[[ct]],'//item/element[@tag="3004,0001"]',xmlValue)[[ct]]
        DVHList[[ROIName]][["DoseUnits"]]<-xpathApply(a[[ct]],'//item/element[@tag="3004,0002"]',xmlValue)[[ct]]
        DVHList[[ROIName]][["DoseType"]]<-xpathApply(a[[ct]],'//item/element[@tag="3004,0004"]',xmlValue)[[ct]]
        DVHList[[ROIName]][["DVHDoseScaling"]]<-xpathApply(a[[ct]],'//item/element[@tag="3004,0052"]',xmlValue)[[ct]]
        DVHList[[ROIName]][["DVHVolumeUnits"]]<-xpathApply(a[[ct]],'//item/element[@tag="3004,0054"]',xmlValue)[[ct]]
        DVHList[[ROIName]][["DVHNumberOfBins"]]<-xpathApply(a[[ct]],'//item/element[@tag="3004,0056"]',xmlValue)[[ct]]
        DVHList[[ROIName]][["DVHMeanDose"]]<-xpathApply(a[[ct]],'//item/element[@tag="3004,0074"]',xmlValue)
        if(length(DVHList[[ROIName]][["DVHMeanDose"]])>1) DVHList[[ROIName]][["DVHMeanDose"]]<-DVHList[[ROIName]][["DVHMeanDose"]][[ct]]
        DVHList[[ROIName]][["DVHMaximumDose"]]<-xpathApply(a[[ct]],'//item/element[@tag="3004,0072"]',xmlValue)
        if(length(DVHList[[ROIName]][["DVHMaximumDose"]])>1) DVHList[[ROIName]][["DVHMaximumDose"]]<-DVHList[[ROIName]][["DVHMaximumDose"]][[ct]]
        DVHList[[ROIName]][["ReferencedROINumber"]]<-xpathApply(a[[ct]],'//item/sequence[@tag="3004,0060"]//element[@tag="3006,0084"]',xmlValue)[[ct]]
        
        dvhString<-xpathApply(a[[ct]],'//item/element[@tag="3004,0058"]',xmlValue)[[ct]];
        dvhArr<-strsplit(dvhString,"\\\\");
        DVHList[[ROIName]][["DVHData.volume"]]<-as.numeric(dvhArr[[1]][seq(2,length(dvhArr[[1]]),by=2 )])
        
        DVHList[[ROIName]][["DVHData.dose"]]<-cumsum(  as.numeric(dvhArr[[1]][seq(1,length(dvhArr[[1]]),by=2 )])  )
        DVHList[[ROIName]][["DVHData.dose"]]<- DVHList[[ROIName]][["DVHData.dose"]] - ( (DVHList[[ROIName]][["DVHData.dose"]][2]-DVHList[[ROIName]][["DVHData.dose"]][1])/2 )
        dvh.type<-''
        final.matrix<-as.matrix(  cbind( DVHList[[ROIName]][["DVHData.volume"]],DVHList[[ROIName]][["DVHData.dose"]] ) )
        if(DVHList[[ROIName]][["DVHType"]]=="CUMULATIVE") dvh.type<-"cumulative";
        if(DVHList[[ROIName]][["DVHType"]]=="DIFFERENTIAL") dvh.type<-"differential";
        if(dvh.type=='') 
          logObj$handle( "error" , "type of DVH not yet supported"  );
        if(DVHList[[ROIName]][["DoseUnits"]]!='GY') 
          logObj$handle( "error" , "only 'GY' are supported as DoseUnit"  );
        if(DVHList[[ROIName]][["DoseType"]]!='PHYSICAL') 
          logObj$handle( "error" , "only DoseType 'physical' is supported"  );
        if(as.numeric(DVHList[[ROIName]][["DVHDoseScaling"]])!=1) 
          logObj$handle( "error" , "only DVHDoseScaling equal to 1 is supported"  );
        if(DVHList[[ROIName]][["DVHVolumeUnits"]]!='CM3') 
          logObj$handle( "error" , "only DVHVolumeUnits equal to 'CM3' is supported  ");
        final.matrix<-as.matrix(cbind(DVHList[[ROIName]][["DVHData.dose"]],DVHList[[ROIName]][["DVHData.volume"]]))
        
        DVHObj<-new("dvhmatrix", dvh=final.matrix, dvh.type=dvh.type, vol.distr='absolute', volume=final.matrix[1,1])
        DVHList[[ROIName]][["DVHObj"]]<-DVHObj
      }
    }    
    
    return(list(
      "ImagePositionPatient"=ImagePositionPatient,"ImageOrientationPatient"=ImageOrientationPatient,
      "Rows"=Rows,"Columns"=Columns,"PixelSpacing"=PixelSpacing,"PixelRepresentation"=PixelRepresentation,
      "SOPInstanceUID"=SOPInstanceUID,"SeriesInstanceUID"=SeriesInstanceUID,"GridFrameOffsetVector"=GridFrameOffsetVector,
      "DoseUnits"=DoseUnits,"DoseType"=DoseType,"DoseGridScaling"=DoseGridScaling,
      "ReferencedRTPlanSequence_ReferencedSOPInstanceUID"=ReferencedRTPlanSequence_ReferencedSOPInstanceUID,
      "DVHList"=DVHList
    ))
  }  
  #=================================================================================
  # getImageFromRAW
  # build a row data from a DICOM file stored on filesystem and load it 
  # into memory (using DCMTK)
  #=================================================================================
  getImageFromRAW<-function(fileName) {
    
    objSV<-services()
    fileNameRAW<-paste(fileName,".0.raw")    
    fileNameRAW<-str_replace_all(string = fileNameRAW , pattern = " .0.raw",replacement = ".0.raw")
    
    if(!file.exists(fileName)) logObj$handle( "error" , " the fileName is missing in geoLet::getImageFromRAW()"  );
    
    pathToStore<-substr(fileName,1,tail(which(strsplit(fileName, '')[[1]]=='/'),1)-1)
    if(!file.exists( fileNameRAW )  | folderCleanUp==TRUE) {
      stringa1<-"dcmdump";
      if ( Sys.info()["sysname"] == "Windows") {
        fileNameFS<-chartr("\\","/",fileName);
        stringa2<-chartr("/","\\\\",stringa1)
      }
      else fileNameFS<-fileName;
      stringa2<-paste(" +W  ",pathToStore,fileNameFS,collapse='')
      options(warn=-1)
      stringone<-as.character(paste( c(stringa1," ",stringa2),collapse=''))
      # gestisci le system call in maniera diversa in funzione che sia WINDOWS o LINUX
      if ( Sys.info()["sysname"] == "Windows") {
        res<-.C("executeCMDLine",  as.character(stringone), as.integer(str_length(stringone))  )
      }
      else {
        system2(stringa1,stringa2,stdout=NULL)
      }
      options(warn=0)
    }
    rowsDICOM<-as.numeric(getDICOMTag(fileName = fileName,tag = '0028,0010'))
    columnsDICOM<-as.numeric(getDICOMTag(fileName = fileName,tag = '0028,0011'))
    bitsAllocated<-as.numeric(getDICOMTag(fileName = fileName,tag = '0028,0100'))
    if(SOPClassUIDList[[fileName]]$kind!="RTDoseStorage"){
      if(bitsAllocated!=16) 
        logObj$handle( "error" , "16bit pixel are allowed only for non-RTDoseStorage"  );
      if ( Sys.info()["sysname"] == "Windows") {
        fileNameRAWFS<-chartr("\\","/",fileNameRAW);
        fileNameRAWFS<-chartr("/","\\\\",fileNameRAWFS);
      }
      else fileNameRAWFS<-fileNameRAW;
      
      if(!file.exists(fileNameRAWFS)) logObj$handle( "error" , "problem in creating image binary file in geoLet::getImageFromRAW()"  );
      
      rn<-readBin(con = fileNameRAWFS, what="integer", size=2, endian="little",n=rowsDICOM*columnsDICOM)    
      rn<-matrix(rn,ncol=columnsDICOM, byrow = TRUE)
    }
    if(SOPClassUIDList[[fileName]]$kind=="RTDoseStorage"){
      if(bitsAllocated==32) {
        if(SOPClassUIDList[[fileName]]$kind!="RTDoseStorage") 
          logObj$handle( "error" , "32bit pixel are allowed only for RTDoseStorage"  );
        numberOfFrames<-as.numeric(getDICOMTag(fileName = fileName, tag = '0028,0008'))
        if ( Sys.info()["sysname"] == "Windows") {
          fileNameRAWFS<-chartr("\\","/",fileNameRAW);
          fileNameRAWFS<-chartr("/","\\\\",fileNameRAWFS);
        }
        else fileNameRAWFS<-fileNameRAW;
        
        if(!file.exists(fileNameRAWFS)) logObj$handle( "error" , "problem in creating image binary file in geoLet::getImageFromRAW()"  );
        
        rn<-readBin(con = fileNameRAWFS, what="integer", size=4, endian="little",n=rowsDICOM*columnsDICOM*numberOfFrames) 
        # per ora va via come ciclo FOR, poi ci ragioniamo.... (per le performances)        
        matRN<-array(0,c(rowsDICOM,columnsDICOM,numberOfFrames))
        ct<-1
        for( z in seq(1,numberOfFrames)) {
          for(x in seq(1,rowsDICOM)) {
            for(y in seq(1,columnsDICOM)) {
              matRN[x,columnsDICOM-y,z]<-rn[ct]
              ct<-ct+1 
            }
          }
        }        
        new_atRN<-array(0,c(columnsDICOM,rowsDICOM,numberOfFrames))
        for(ct in seq(1:dim(matRN)[3]  )) {
          new_atRN[,,ct]<-t(objSV$rotateMatrix( matRN[,,ct], rotations=2 ))
        }
        rn<-new_atRN
      } else  {
        numberOfFrames<-as.numeric(getDICOMTag(fileName = fileName, tag = '0028,0008'))
        
        if ( Sys.info()["sysname"] == "Windows") {
          fileNameRAWFS<-chartr("\\","/",fileNameRAW);
          fileNameRAWFS<-chartr("/","\\\\",fileNameRAWFS)
        }
        else fileNameRAWFS<-fileNameRAW;
        
        if(!file.exists(fileNameRAWFS)) logObj$handle( "error" , "problem in creating image binary file in geoLet::getImageFromRAW()"  );
        
        rn<-readBin(con = fileNameRAWFS, what="integer", size=2, endian="little",n=rowsDICOM*columnsDICOM*numberOfFrames)
        matRN<-array(0,c(rowsDICOM,columnsDICOM,numberOfFrames))
        ct<-1
        for( z in seq(1,numberOfFrames)) {
          for(x in seq(1,rowsDICOM)) {
            for(y in seq(1,columnsDICOM)) {
              matRN[x,columnsDICOM-y,z]<-rn[ct]
              ct<-ct+1 
            }
          }
        }        
        new_atRN<-array(0,c(columnsDICOM,rowsDICOM,numberOfFrames))
        for(ct in seq(1:dim(matRN)[3]  )) {
          new_atRN[,,ct]<-t(objSV$rotateMatrix( matRN[,,ct], rotations=2 ))
        }
        rn<-new_atRN        
      }
    }    
    return(rn)
  }  
  #####################################################################################
  # getFolderContent: check a folder and find out the content in terms of DICOM objects
  #   
  # INPUT:
  #   - pathToOpen: il path da aprire
  #   - defaultExtension: estensione di default è il .dicom
  #       occured
  # OUTPUT:
  #   - lista di SOP Class UID per ciascun DICOM presente nel path indicato 
  #################################################################################
  getFolderContent<-function(pathToOpen, defaultExtension) {
    
    # if no path is given, use the set one
    if(!dir.exists(pathToOpen)) logObj$handle( "error" , "The indicate Path does not exist"  );
    
    # salva in un array tutti i DICOM presenti nella cartella
    DCMFilenameArray<-list.files(pathToOpen,defaultExtension)    
    
    # lista con la SOP Class UID di ciascun DICOM
    SOPClassUIDList<-list()
    ImagingPositionArray <- c()
    Iteration <- 0
    ImageNumber <- 0
    # browser()
    for(i in 1:length(DCMFilenameArray) ) {
      fileNameWithPath<-paste(pathToOpen,"/",DCMFilenameArray[i] , sep="");
      if( substr(fileNameWithPath,nchar(fileNameWithPath)-3,nchar(fileNameWithPath))=='.dcm' ) {
        valore<-getDICOMTag( fileName = fileNameWithPath, tag = "0008,0016")
        # do the system call
        SOPClassUIDList[[fileNameWithPath]]<-list();
        SOPClassUIDList[[fileNameWithPath]]$tag<-valore
        SOPClassUIDList[[fileNameWithPath]]$kind<-"Unknown"
        # browser()
        if( SOPClassUIDList[[fileNameWithPath]]$tag == "1.2.840.10008.5.1.4.1.1.2" ) SOPClassUIDList[[fileNameWithPath]]$kind<-"CTImageStorage"
        if( SOPClassUIDList[[fileNameWithPath]]$tag == "1.2.840.10008.5.1.4.1.1.481.2" ) SOPClassUIDList[[fileNameWithPath]]$kind<-"RTDoseStorage"
        if( SOPClassUIDList[[fileNameWithPath]]$tag == "1.2.840.10008.5.1.4.1.1.481.3" ) SOPClassUIDList[[fileNameWithPath]]$kind<-"RTStructureSetStorage"
        if( SOPClassUIDList[[fileNameWithPath]]$tag == "1.2.840.10008.5.1.4.1.1.481.5" ) SOPClassUIDList[[fileNameWithPath]]$kind<-"RTPlanStorage"
        if( SOPClassUIDList[[fileNameWithPath]]$tag == "1.2.840.10008.5.1.4.1.1.4" ) SOPClassUIDList[[fileNameWithPath]]$kind<-"MRImageStorage"
        if( SOPClassUIDList[[fileNameWithPath]]$tag == "1.2.840.10008.5.1.4.1.1.128" ) SOPClassUIDList[[fileNameWithPath]]$kind<-"PositronEmissionTomographyImageStorage"
        if( SOPClassUIDList[[fileNameWithPath]]$tag == "1.2.840.10008.5.1.4.1.1.2.1" ) SOPClassUIDList[[fileNameWithPath]]$kind<-"CTImageStorage"
        
        # if it is an image check istance number TAG availability
        if (SOPClassUIDList[[fileNameWithPath]]$kind=="CTImageStorage" | SOPClassUIDList[[fileNameWithPath]]$kind=="MRImageStorage"|
            SOPClassUIDList[[fileNameWithPath]]$kind=="PositronEmissionTomographyImageStorage") {
          ImageNumber <- ImageNumber+1 
          
          # verify istance number TAG
          InstanceNumberI<-getDICOMTag( fileName = fileNameWithPath, tag = "0020,0013")
          # if istance number not available
          if (identical(InstanceNumberI, "")){
            
            # count number of images with no information about the istance number
            Iteration <- Iteration+1
            # load z imaging position patient
            ImagingPosition<-splittaTAG(getDICOMTag( fileName = fileNameWithPath, tag = "0020,0032"))[3]
            # save z imaging position patient in a vector and into the SOPClassUIDList
            ImagingPositionArray <- c(ImagingPositionArray, ImagingPosition)
            SOPClassUIDList[[fileNameWithPath]]$zImagingPosition <- ImagingPosition
          }
        }
      }
    }
    
    if ( Iteration!=0 ){
      # check if number of patients (excluding RT struct, RT dose e RT plan) is equal to the number of
      # images with no information about the Istance Number
      if( ImageNumber!=Iteration ){
        logObj$handle(type="error",msg="Attenzione! C'e' almeno un Immagine che presenta la TAG Istance Number")
      }
      # sort ascending z Imaging Position Array
      SortedImagingPositionArray <- sort(ImagingPositionArray)
      # for each SOPClassUIDList
      for(i in 1:length(SOPClassUIDList)){
        # only for the images, calculatedIstanceNumber as position number based on imaging position patient 
        if(SOPClassUIDList[[i]]$kind=="CTImageStorage" | SOPClassUIDList[[i]]$kind=="MRImageStorage"|
          SOPClassUIDList[[i]]$kind=="PositronEmissionTomographyImageStorage"){
          
          SOPClassUIDList[[i]]$calculatedIstanceNumber <- which(SortedImagingPositionArray==SOPClassUIDList[[i]]$zImagingPosition)
          SOPClassUIDList[[i]]$calculatedIstanceNumber <- as.character(SOPClassUIDList[[i]]$calculatedIstanceNumber)
        }
      }
    }
    
    # applica un primo filtro basato sull'esperienza di Tarducci
    SOPClassUIDList<-filter_01(SOPClassUIDList)
    
    return(SOPClassUIDList);
  }
  ####################################################################################################################
  # filter_01: controlla che il FrameOfReferenceUID dell'RTstruct coincide con il FrameOfReferenceUID
  # delle SOPclassUID e dell'RTdose  
  #
  # INPUT:
  #   - SOPClassUIDList: lista di SOP Class UID per ciascun DICOM presente nel path indicato
  #   
  # OUTPUT:
  #   - oggetti DICOM che presentano il FrameOfReferenceUID uguale all'RTstruct e il corrispettivo FrameOfReferenceUID
  #####################################################################################################################
  filter_01<-function(SOPClassUIDList){
    new.SOPClassUIDList<-list()
    
    # Per ogni file indicato in SOPClassUIDList
    for(fileName in names(SOPClassUIDList)) {
      
      # Se RTstructure è presente e il FrameOfReferenceUID coincide con il 
      # ReferencedFrameOfReferenceUID, salva il valore del FrameOfReferenceUID
      if (SOPClassUIDList[[fileName]]$kind == "RTStructureSetStorage"){
        
        # Check sul numero di RTstruct presenti nel file se piu' di uno restituisce errore
        numRTstruct <- lapply(SOPClassUIDList[[fileName]], function(x) x[[1]])
        if (length((numRTstruct)$kind)==1){
          document <- getXMLStructureFromDICOMFile( fileName = fileName )
          FrameReference <- xpathApply(document,'//element[@tag="0020,0052" and @name="FrameOfReferenceUID"]',xmlValue)[[1]]
          arr.ReferencedFrameReference <- unlist(xpathApply(document, '//element[@tag="3006,0024" and @name="ReferencedFrameOfReferenceUID"]',xmlValue))
          if(sum(arr.ReferencedFrameReference %in% FrameReference) == length(arr.ReferencedFrameReference)) FrameReferenceFin <- FrameReference
          else {
            logObj$handle(type="err",msg="Attenzione! C'e' almeno un ReferencedFrameOfReference diverso dal FrameOfReference principale dell RTStruct")
          }
        }
        else{
          logObj$handle(type="err",msg="Attenzione! C'e' piu' di un RTstruct")
        }
      }
      else new.SOPClassUIDList[[fileName]] <- SOPClassUIDList[[fileName]]
    }
    
    new.new.SOPClassUIDList<-list()
    
    # Per ogni file indicato in SOPClassUIDList
    for(fileName in names(SOPClassUIDList)) {
      
      # Se il FrameOfReferenceUID della SOP class UID trovata coincide con il FrameOfReferenceUID
      # dell'RTstruct, salva il corrispettivo oggetto dicom altrimenti lo salta
      if (SOPClassUIDList[[fileName]]$kind == "CTImageStorage" |
          SOPClassUIDList[[fileName]]$kind == "RTDoseStorage" |
          SOPClassUIDList[[fileName]]$kind == "MRImageStorage" | 
          SOPClassUIDList[[fileName]]$kind == "PositronEmissionTomographyImageStorage"  ){
        document <- getXMLStructureFromDICOMFile( fileName = fileName )
        FrameReference <- xpathApply(document,'//element[@tag="0020,0052" and @name="FrameOfReferenceUID"]',xmlValue)[[1]]
        if(FrameReference==FrameReferenceFin)  {
          new.new.SOPClassUIDList[[fileName]] <- SOPClassUIDList[[fileName]]
        } 
        else {
          logObj$handle(type="warning",msg=c("Attenzione! l'immagine ",fileName," non ha un FOR associabile all'RTStruct"))
        }
      } 
      else new.new.SOPClassUIDList[[fileName]] <- SOPClassUIDList[[fileName]]
    }
    
    return(new.new.SOPClassUIDList);
  }
  #=================================================================================
  # NAME: getTag
  # restituisce all'utente una TAG prenddola dall'XML
  #=================================================================================    
  getTag<-function(tag=tag, fileName="",whichFile='first', whichSerie='firstImageSerie') {   
    
    if(fileName!="") return(getDICOMTag(tag = tag , fileName = fileName) )
    if(whichSerie!='firstImageSerie') return(getDICOMTag(tag = tag , whichSerie = whichSerie, whichFile = whichFile)  )
    
    serieName<-giveBackImageSeriesInstanceUID();
    fileName<-dataStorage$info[[serieName]][[1]]$fileName
    return(getDICOMTag(tag = tag , fileName = fileName , whichFile = '' , whichSerie = serieName))
  }
  #####################################################################################
  # getDICOMTag: estrae una tag dall'XML
  #   
  # INPUT:
  #   - tag: tag considerato
  #   - fileName: nome del oggetto DICOM
  #   - whichFile: quale file
  #   - whichSerie: quale serie
  #     
  # OUTPUT:
  #   - Valore della tag richiesta
  ################################################################################# 
  getDICOMTag<-function(tag=tag, fileName="",whichFile='first', whichSerie='firstImageSerie') {    
    obj.S<-services();
    # cat("\n ** ",fileName)
    # exemption: you want an Image!
    if(tag == "7fe0,0010") return( getImageFromRAW(fileName)  ); 
    
    return(obj.S$getDICOMTag(tag = tag,fileName = fileName, folderCleanUp = attr_folderCleanUp ))
  } 
  #=================================================================================
  # getROIList
  # restituisce la lista delle ROI
  #=================================================================================  
  getROIList<-function() {
    mat2Ret<-matrix( c(seq(1,length(names(dataStorage$structures))),names(dataStorage$structures)),nrow=2 ,byrow=T )
    return(mat2Ret)
  } 
  #=================================================================================
  # getXMLStructureFromDICOMFile
  # this function is a wrapper for the method from Services with the same name
  # the reason of the wrap is due to the need to handle a cache
  #=================================================================================  
  getXMLStructureFromDICOMFile<-function(fileName, folderCleanUp = FALSE) {
    obj.S<-services();
    # Carica l'XML se non in cache
    if( attr_loadXMLInCache == TRUE &  length(attr_arrayXMLCache[[fileName]])!=0 ) {
      doc<-attr_arrayXMLCache[[fileName]]
    }
    else {  
      doc<-obj.S$getXMLStructureFromDICOMFile(fileName = fileName, folderCleanUp = folderCleanUp)
      # mettilo in cache se è previsto di farlo
      if(attr_loadXMLInCache == TRUE) attr_arrayXMLCache[[fileName]]<<-doc      
    }
    return(doc)
  }
  #=================================================================================
  # getAlignedStructureAndVoxelCube
  # rimappa le coordinate delle ROI rispetto alla geometria bitmap del voxelCuba
  # valori reali sono ovviamente ammessi anche se la geometria dei voxelcube
  # sarebbe per sua natura a domino N^3... (orbo se beccherai mai un vertice esattamente
  # posizionato nel centroide di un voxel, dopo le rotazioni/interpolazioni :) )
  #=================================================================================    
  getAlignedStructureAndVoxelCube<-function(  ps.x=NA, ps.y=NA, ps.z=NA, ROIName ) {
    objService <- services()
    
    if(is.na(ps.x) & is.na(ps.y) & is.na(ps.z)) {
      ps.x<-getPixelSpacing()[1]
      ps.y<-getPixelSpacing()[2]
      ps.z<-getPixelSpacing()[3]
    }
    voxelCube<-getImageVoxelCube( ps.x = ps.x, ps.y = ps.y, ps.z = ps.z ) 
    ROI<-rotateToAlign(ROIName = ROIName)
    old.ps<-getPixelSpacing();
    delta.x<-old.ps[1] / ps.x;
    delta.y<-old.ps[2] / ps.y;
    delta.z<-old.ps[3] / ps.z;
    #    browser();
    for(i in names(ROI$pointList)) {
      for( ct in seq(1,length(ROI$pointList[[i]]) ) ) {
        ROI$pointList[[i]][[ct]][,3]<-ROI$pointList[[i]][[ct]][,3] * delta.z;
        ROI$pointList[[i]][[ct]][,2]<-ROI$pointList[[i]][[ct]][,2] * delta.y;
        ROI$pointList[[i]][[ct]][,1]<-ROI$pointList[[i]][[ct]][,1] * delta.x;
      }
    }
    invisible(gc())
    
    # -im 
    SeriesInstanceUID<-giveBackImageSeriesInstanceUID()
    if(SeriesInstanceUID == '' ) {
      logObj$handle( "error" , "missing CT/MR series"  );
    }
    numberOfSlices<-length(dataStorage$img[[SeriesInstanceUID]]);
    
    initial.point <- objService$get3DPosFromNxNy(Nx = 1,Ny = 1,oM = dataStorage$info[[SeriesInstanceUID]][[1]]$orientationMatrix)
    final.point <- objService$get3DPosFromNxNy(Nx = 100,Ny = 100,oM = dataStorage$info[[SeriesInstanceUID]][[numberOfSlices]]$orientationMatrix)
    directions <- final.point - initial.point
    
    if(directions[3] < 0) {
      cat("\n -----------------------------------------")
      cat("\n Attenzione, rotazione non testata: verificare l'outcome")
      cat("\n -----------------------------------------")
      aaa <- voxelCube
      aaa[,, seq(1,numberOfSlices)] <- voxelCube[,, seq(numberOfSlices,1,by=-1)]
      voxelCube<-aaa
      # voxelCube[,,i]<-t(objService$rotateMatrix(voxelCube[,,i],rotations=1))
    }
    if(directions[1] < 0) {
      cat("\n ------------------------------------------------------")
      cat("\n ATTENZIONE: orientamento supportato ma non ancora testato")
      cat("\n si prega di contattare il manutentore del pacchetto")
      cat("\n ------------------------------------------------------")
      # stop()
    }
    if(directions[2] < 0) {
      cat("\n ------------------------------------------------------")
      cat("\n ATTENZIONE: orientamento supportato ma non ancora testato")
      cat("\n si prega di contattare il manutentore del pacchetto")
      cat("\n ------------------------------------------------------")
      # stop()
    }
    
    # -fm
    return( list( "voxelCube"=voxelCube, "ROI"=ROI$pointList )  )
  }  
  rotateToAlign<-function(ROIName) {
    SeriesInstanceUID<-giveBackImageSeriesInstanceUID();
    tabella1<-getAssociationTable( "SOPInstance_vs_SliceLocation" , ROIName )
    pointList<-getROIPointList( ROIName )
    newPointList<-pointList
    iterazione<-1
    
    for(index in names( pointList )) {
      for(internalIndex in seq(1,length(pointList[[index]]))) {
        IMGsliceInfo<-dataStorage$info[[SeriesInstanceUID]][[ tabella1[ which(tabella1[,4]==index)  ,5  ]  ]]
        sliceLocation<-as.numeric(tabella1[which(tabella1[,4]==index),5])
        DOM<-IMGsliceInfo$orientationMatrix
        
        m<-pointList[[index]][[internalIndex]];
        
        numeroRighe<-dim(m)[1]
        for(i in seq(1,numeroRighe)) {
          valori<-calcolaNX(m[i,],DOM )
          newPointList[[index]][[internalIndex]][i,1]<-valori$Nx
          newPointList[[index]][[internalIndex]][i,2]<-valori$Ny
          newPointList[[index]][[internalIndex]][i,3]<-sliceLocation
        }
      }
      iterazione<-iterazione+1
    }
    return( list("pointList"=newPointList) ) 
  }  
  getAssociationTable<-function( tipoTabella="SOPInstance_vs_SliceLocation", ROIName ) {
    SeriesInstanceUID<-giveBackImageSeriesInstanceUID();
    
    if(tipoTabella=="SOPInstance_vs_SliceLocation") {
      matrice<-c()
      involvedCT<-names(dataStorage$structures[[ROIName]]);
      
      for(index in names(dataStorage$info[[SeriesInstanceUID]]) ) {
        if(!is.null(dataStorage$info[[SeriesInstanceUID]][[index]]$ROIList)) {
          matrice<-rbind(matrice,cbind(dataStorage$info[[SeriesInstanceUID]][[index]]$ROIList,index) );
        }
      }
      matrice<-matrice[which(matrice[,1]==ROIName),]
      return(matrice)
    }
  }  
  calcolaNX<-function(riga,DOM) {
    Px<-riga[1];    Py<-riga[2];
    a11<-DOM[1,1];    a21<-DOM[2,1];    a31<-DOM[3,1];    a12<-DOM[1,2];
    a22<-DOM[2,2];    a32<-DOM[3,2];    Sx<-DOM[1,4];     Sy<-DOM[2,4];     Sz<-DOM[3,4]; 
    Nx<-(a22*Px-a12*Py-a22*Sx+a12*Sy)/(a11*a22-a21*a12);
    Ny<-(a11*Py-a21*Px+a21*Sx-a11*Sy)/(a22*a11-a21*a12);
    Nz<-Sz;
    return(list("Nx"=Nx,"Ny"=Ny,"Nz"=Nz))
  }   
  #=================================================================================
  # giveBackImageSeriesInstanceUID
  # from dataStorage it gives back the SOPInstanceUID of the series which has a 
  # SOPClassUID as 'CTImageStorage' or 'MRImageStorage' or 'PET'. 
  #=================================================================================
  giveBackImageSeriesInstanceUID<-function(FrameOfReferenceUID=NA, ROIName=NA) {
    
    # se è per FrameOfReferenceUID
    if(!is.na(FrameOfReferenceUID)) return(searchIMGSeriesForFrameOfReferenceUID(FrameOfReferenceUID=FrameOfReferenceUID));
    # se le vuoi tutte...
    list.index<-''
    for(whichIdentifier in seq(1,length(dataStorage$info) )) {
      if(names(dataStorage$info)[whichIdentifier]!="structures" &
         names(dataStorage$info)[whichIdentifier]!="plan"  &
         names(dataStorage$info)[whichIdentifier]!="DVHs"  &
         names(dataStorage$info)[whichIdentifier]!="doses" ) {
        if(dataStorage$info[[whichIdentifier]][[1]]$SOPClassUID == 'CTImageStorage' ||
           dataStorage$info[[whichIdentifier]][[1]]$SOPClassUID == 'MRImageStorage' ||
           dataStorage$info[[whichIdentifier]][[1]]$SOPClassUID == 'PositronEmissionTomographyImageStorage'
        ) {
          #list.index<-whichIdentifier
          list.index<-names(dataStorage$info)[whichIdentifier]
        }
      }
    }
    return(list.index)
  } 
  #=================================================================================
  # searchIMGSeriesForFrameOfReferenceUID
  # restituisce la seriesInstanceUID della serie di immagini che fa riferimento 
  # ad un dato FrameOfReferenceUID
  #=================================================================================
  searchIMGSeriesForFrameOfReferenceUID<-function(FrameOfReferenceUID) {
    listaSeriesInstanceUID<-names(dataStorage$img)
    for( index4Info in names(dataStorage$img)) {
      if(index4Info!="structures") {
        if(dataStorage$info[[index4Info]][[1]]$FrameOfReferenceUID == FrameOfReferenceUID) {
          return(index4Info);
        }
      }
    }
    return(NA)
  }   
  getAttribute<-function( attribute ) {
    if( attribute == "ROIVoxelMemoryCache") return (attr_ROIVoxelMemoryCache)
    if( attribute == "folderCleanUp") return (attr_folderCleanUp)
    if( attribute == "needPhysicalFile") return (attr_needPhysicalFile)
    if( attribute == "dataStorage") return( dataStorage )
  }
  #   #=================================================================================
  #   # buildOrientationMatrix
  #   # internal and stupid function to build Orientation Matrix
  #   #=================================================================================
  #   buildOrientationMatrix<-function(fileName) {
  #     iPP<-getAttribute(attribute="ImagePositionPatient",fileName=fileName)
  #     iOP<-getAttribute(attribute="ImageOrientationPatient",fileName=fileName)
  #     matrice<-matrix(c(iOP[1],iOP[2],iOP[3],0,iOP[4],iOP[5],iOP[6],0,0,0,0,0,iPP[1],iPP[2],iPP[3],1),ncol=4); 
  #     return(matrice)
  #   }  
  #=================================================================================
  # splittaTAG
  # internal and stupid function useful to kill artifacts from 
  # a text taken from DCMTK files
  #=================================================================================
  splittaTAG<-function(stringa) {
    return( as.numeric(strsplit(stringa,split = "\\\\")[[1]])   )  
  }  
  class<-function() {
    return("geoLet");
  }
  #=================================================================================
  # Constructor
  #=================================================================================
  constructor<-function( ROIVoxelMemoryCache = TRUE , 
                         folderCleanUp = FALSE,
                         loadXMLInCache = FALSE,
                         loadRAWInCache = TRUE,
                         defaultExtension = ".dcm") {
    # Attributes - set by user
    attr_folderCleanUp<<-folderCleanUp                      # force to re-dump DICOM files
    attr_ROIVoxelMemoryCache<<-ROIVoxelMemoryCache          # force to cache ROI Voxel
    attr_ROIVoxelMemoryCacheArray<<-list();
    attr_arrayXMLCache<<-list();                             # array containint XML files
    
    # Internal Attributes
    attr_mainFrameOfReferenceUID<<-NA                       # frameOfReference (geometry)
    logObj<<-logHandler()                                   # log/error handler Object 
    dataStorage<<-list()                                    # memory data structure
    attr_dataChache<<-list();
    attr_attributeList<<-''
    attr_loadXMLInCache<<-loadXMLInCache
    attr_loadRAWInCache<<-loadRAWInCache
    attr_arrayXMLCache<<-list();
    attr_arrayRAWCache<<-list();
    attr_ROI.non.compl<<-c();
    attr_defaultExtension<<-defaultExtension
    #    cat("\n folderCleanUp=",folderCleanUp, "  loadXMLInCache =",loadXMLInCache," loadRAWInCache=",loadRAWInCache)
  }  
  constructor( ROIVoxelMemoryCache = ROIVoxelMemoryCache , 
               folderCleanUp = folderCleanUp, 
               loadXMLInCache = loadXMLInCache,
               loadRAWInCache = loadRAWInCache,
               defaultExtension = defaultExtension)
  return(list(
    "openDICOMFolder" = openDICOMFolder,
    "getMESHfromROI"=getMESHfromROI,
    "getImageVoxelCube" = getImageVoxelCube,
    "getPixelSpacing"= getPixelSpacing,
    "getROIList"=getROIList,
    "getROIPointList"=getROIPointList,
    "getTag"=getTag,
    "getROIVoxels"=getROIVoxels,
    "getEmptyStructure"=getEmptyStructure,
    "extractDoseVoxels" = extractDoseVoxels,
    "calculateDVH" = calculateDVH,
    "getAttribute"= getAttribute,
    "getAlignedStructureAndVoxelCube"=getAlignedStructureAndVoxelCube,
    "rotateToAlign"=rotateToAlign,
    "calcolaNX"=calcolaNX,
    "class"=class,
    "extractDoseVoxels.part"=extractDoseVoxels.part,
    "getFolderContent"=getFolderContent
  ))
}
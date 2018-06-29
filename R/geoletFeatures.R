#' f.extractor.pluri.LoG.par
#'
#' @description  Instantiate an object of the class \code{geoLet}.This represents just the classname (f.extractor.pluri.LoG.par)
#' @export
f.extractor.pluri.LoG.par<- function( path, ROIName ,
                                 feature.family=c("stat","morph","glcm","rlm","szm","fractal"),
                                 interpolate = FALSE, px="", py="", error = .1,
                                 fileName = "tmp.f.extractor.pluri.par.RData" ,
                                 forceRecalculus = TRUE, from.sigma=1, to.sigma=2, def.by = .1, strategy="fixed",
                                 sigma.array=c()) {
  
  # sigma.array <-c()
  big.matrice<-list()
  if(strategy=="fixed") sigma.array <- seq( from.sigma, to.sigma, by = def.by)
# browser()
  
  if(forceRecalculus==FALSE) {
    if(file.exists(fileName)==TRUE) {
      load(file = fileName)  
    }
  } else {

    if(file.exists(fileName)==TRUE) {
      file.remove(fileName)
    }
  }  
  
  for( iterazione in sigma.array) {
    
    if( !(iterazione %in% names(big.matrice))) {
    
      filterPipeline<- list()
      filterPipeline[[1]]<-list("kernel.type"="LoG", "sigma"=iterazione)
  
      nomeFile <- paste(c( "f.extractor.single.par_",iterazione,".RData" ),collapse='')
      
      tmp.matrice <- f.extractor.sing.par(path = path,ROIName = ROIName,feature.family = feature.family,
                           filterPipeline = filterPipeline,interpolate = interpolate,px = px,py = py,
                           fileName = nomeFile, forceRecalculus = forceRecalculus)
      
      big.matrice[[as.character(iterazione)]]<-tmp.matrice
      save(big.matrice, file = fileName)
      file.remove(nomeFile)
      
    }
  }
  return(big.matrice);
}

#' f.extractor.sing.par
#'
#' @description  Instantiate an object of the class \code{geoLet}. This represents just the classname (f.extractor.sing.par)
#' @param path The path to be browsed for folders containing DICOM images and Structure Set of patients
#' @param ROIName A \code{character} object containing the name(s) of ROIs to be analyzed. If there is more than one ROI (multiple references) or no ROI is 
#' listed the analysis is stopped
#' @param feature.family The list of features to be extracted, by default all
#' @param filterPipeline A \code{list} containint the filter pipeline to be applied to images before feature extraction
#' @param interpolate A \code{logical} value, by default \code{FALSE}, for interpolating sigma values in filters using sigma for computation as LoG
#' @param px The pixel spacing in \emph{x} direction for interpolating sigma
#' @param py The pixel spacing in \emph{y} direction for interpolating sigma
#' @param error A \code{numeric} value setting the threshold for applying the interpolation
#' @param forceRecalculus A \code{logical} value, by defaul \code{TRUE}, for forcing recalulus
#' @details Using an object of class \code{ROImap} it is possible to list all the ROIs contained in an explored directory for the analysis. It is possible to analyze
#' different ROIs names referencing to the same structure in the different studies by putting al their names a a vector of \code{character} objects in the argument \code{ROIName}
#' @export
f.extractor.sing.par<- function(path, ROIName ,
                          feature.family=c("stat","morph","glcm","rlm","szm","fractal"),
                          filterPipeline=list(), interpolate = FALSE, px="", py="", error = .001,
                          fileName = "tmp.f.extractor.sing.par.RData" ,
                          forceRecalculus = TRUE) {
  
  if((px == "" | py=="") & !(px == "" & py == "" ))  {
    stop("ERRORE: px e py devono essere specificati entrambi!")
  }

  patList <- list.dirs(path = path,recursive = FALSE)
  matrice <- c()
  
  if(forceRecalculus==FALSE) {
    if(file.exists(fileName)==TRUE) {
      load(file = fileName)  
    }   
  } else {
    if(file.exists(fileName)==TRUE) {
      file.remove(fileName)
    }    
  }
  
  for( patID in patList) {
    
    if(!(patID %in% matrice[,1])) {
      temp.ROIName <- ROIName
      ogg.geoLet <- geoLet()
      ogg.geoLet$openDICOMFolder(pathToOpen = patID)    
      #browser();
      cat("\n ---------------------------------\n Computing Patient: ",patID,"\n")
      
      if (length(which((ROIName %in% ogg.geoLet$getROIList()[2,]) == TRUE)) > 1) {  # check for multiple ROI mapped
        str <- ogg.geoLet$getROIList()[2, ][which((ROIName %in% ogg.geoLet$getROIList()[2,]) == TRUE)] # multiple ROI names mapped
        stop(paste('ERROR: more than one ROI mapped in a single structure set:', paste(str, collapse = ', ')))
      }
      if (length(which((ROIName %in% ogg.geoLet$getROIList()[2,]) == TRUE)) == 0) 
        stop(paste('ROI', paste(ROIName, collapse = ', '), 'not mapped into structure set')) # check for not mapped ROI
      ROIName <- ROIName[which((ROIName %in% ogg.geoLet$getROIList()[2,]) == TRUE)]
      riga <- computeFeatures.geoLet( obj.geoLet = ogg.geoLet, ROIName = ROIName ,
                                      feature.family=feature.family, filterPipeline=filterPipeline,
                                      px = px, py = py) 
      # aggiungi il nome del paziente
      riga <- c( patID , riga)
      matrice <- rbind(matrice,riga)    
      # salva la matrice
      save(matrice,file = fileName)
      ROIName <- temp.ROIName
    }
    
  }
  return(matrice)
}

#' map.ROI
#' @description Function for mapping all ROI names in a series a studies
#' @param path A \code{character} value containing the path to browse for DICOM studies
#' @details DICOM RT structure set files often contain multiple \emph{ROI names} referencing to the same structure. \code{moddicom} is desinged to analyze only
#' one structure at a time. Using this function a \code{ROImap} class object is created, listing all ROI names for each case inside the borwsed directory.
#' Using the method \code{print.ROImap} a table showing all ROI names and their frequency in the explored path is shown.
#' @return An object of class \code{ROImap} that is a list containing all mapped ROIs in a series of studies
#' @export
#' @examples ## NOT RUN
#' ListOfROI <- map.ROI(path = '/my_path')
#' ## to get the table of ROIs in the explored path
#' print(ListOfROI)
map.ROI <- function(path) {
  patList <- list.dirs(path = path,recursive = FALSE)
  ROI.map <- list()
  n <- 0
  for(patID in patList) {
    n <- n + 1
    cat("\n ---------------------------------\n Computing Patient: ", patID, "\n")
    ogg.geoLet <- geoLet()
    ogg.geoLet$openDICOMFolder(pathToOpen = patID)  
    ROI.map[[n]] <- GLT.getROIList(obj.geoLet = ogg.geoLet)
    rm(ogg.geoLet)
    gc()
  }
  attr(ROI.map, "class") <- "ROImap"
  return(ROI.map)
}

#' print.ROImap
#' @export print.ROImap
print.ROImap <- function(obj) {
  cat('\nNumber of patients:', length(x = obj))
  ROI.Table <- c()
  # browser()
  for (N in obj) {
    ROI.Table <- c(ROI.Table, N[2, ])
  }
  cat('\n')
  print(table(ROI.Table))
}

#' computeFeatures.geoLet lavora su piu
#'
#' @description  Intiate an object of the class \code{geoLet}.This represents just the classname (f.extractor.sing.par)
#' @export
computeFeatures.geoLet<- function( obj.geoLet, ROIName , feature.family=c("stat","morph","glcm","rlm","szm","fractal"), filterPipeline=c(),
                                   px = "", py ="", error = .1) {
  # browser()
  objS <- services()
  if( (px == "" | py=="") & !(px == "" & py == "" ))  {
    stop("ERRORE: px e py devono essere specificati entrambi!")
  }
  if(px == "" & py == "" )
    ROI <- obj.geoLet$getROIVoxels(Structure = ROIName)
  else{
    pixelSpacing <- obj.geoLet$getPixelSpacing()
    if( (abs(pixelSpacing[1]-as.numeric(px)))>error | (abs(pixelSpacing[2]-as.numeric(py))>error) ){
      ROI <- obj.geoLet$getROIVoxels(Structure = ROIName, new.pixelSpacing=c(px,py))
    }
    else{
      ROI <- obj.geoLet$getROIVoxels(Structure = ROIName)
    }
  }

  # Ragionamento sul pixelspaxing legato alle esigenze del filtro di convoluzione
  old.px <- ROI$geometricalInformationOfImages$pixelSpacing[1]
  old.py <- ROI$geometricalInformationOfImages$pixelSpacing[2]
  old.pz <- as.numeric(ROI$geometricalInformationOfImages$SliceThickness)

  if( px == "") px <- ROI$geometricalInformationOfImages$pixelSpacing[1]
  if( py == "") py <- ROI$geometricalInformationOfImages$pixelSpacing[2]
  pz <-  as.numeric(ROI$geometricalInformationOfImages$SliceThickness)

  if(length(filterPipeline)!=0) {
    res <- FIL.2D.conv.Filter.geoLet(obj.geoLet = obj.geoLet,ROIName = ROIName,
                                     filter.pipeline = filterPipeline,
                                     scaleFactor = 'space', px = px, py = py)

  } else {
    res <- ROI$masked.images
  }

  naVALUE <- as.integer(min(voxelCube = res$voxelCube,na.rm = T)- 10)
  res.noNA <- res$voxelCube
  res.noNA[ which( is.na(res.noNA),arr.ind = T  )  ] <- naVALUE

  # stat.df<-c()
  def <- c()
  if("stat" %in% feature.family) {
    cat("computing stat features...")
    stat.df <- statisticalFeatures(res$voxelCube)
    def <- c(def,stat.df)
  }

  if("morph" %in% feature.family) {
    cat("computing morph features...")
    morph.df <- morphologicalFeatures(res$voxelCube,px=px,py=py,pz=old.pz)
    def <- c(def,morph.df)
  }

  if("glcm" %in% feature.family){
    cat("computing glcm features...")
    F_cm <-  glcmTexturalFeatures(res$voxelCube)
    cm.df <- colMeans(F_cm)
    def <- c(def,cm.df)
  }

  if("rlm" %in% feature.family) {
    cat("computing glrm features...")
    F_rlm <-glrlmTexturalFeatures(res$voxelCube)
    rlm.df <- colMeans(F_rlm)
    def <- c(def,rlm.df)
  }

  if("szm" %in% feature.family) {
    cat("computing szm features...")
    F_szm <- glszmTexturalFeatures(res$voxelCube)
    szm.df <- colMeans(F_szm)
    def <- c(def,szm.df)
  }

  if("fractal" %in% feature.family) {
    cat("computing fractal features...")
    F_fractal <- fractalFeatures(ROI)
    fractal.df <- colMeans(F_fractal)
    def <- c(def,fractal.df )
  }


  return(def)
}

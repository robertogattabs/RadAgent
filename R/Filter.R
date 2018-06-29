# ================================================================
#' Apply a filter to a 2d or a 3d array wow
#' 
#' @description  It applies filter to a 2d or 3d array. '...' is banned, so shit cannot flies ( ;) ).
#' @param arr2App is the array containing the image that should be filtered
#' @param kernel.type is the kernel you want to use: \code{gaussian}, \code{laplacian}, \code{LoG},\code{shrpen}, \code{emboss}, \code{sobel}
#' @param sigma is the \code{sigma} value for gaussian filter
#' @return the filteres array (same geometry of the one in input)
#' @examples \dontrun{
#' # Create an instante of new.mmButo and load some cases
#' a<-array( 0, dim=c(10,10,10))
#' a[5,5,5]<-177; a[5,6,5]<-198; a[6,6,5]<-45; a[6,5,5]<-12;
#' b<-FIL.applyFilter( a , kernel.type="gaussian", sigma=1.4)
#' }#' 
#' @export
#' @useDynLib moddicomV2
#' @import spatialfil
FIL.applyFilter<-function( arr2App, kernel.type, sigma = 1.4 ) {
  if( is.null(sigma) ) sigma = 1.4;
  # inp.arr <- arr2App
  
  # memorizza le dimensioni originali, in caso debba fare un reshape
  # per rendere la matrice quadrata 
  original.dim <- c(  dim(arr2App)[1] , dim(arr2App)[2]   ) 
  
  if(dim(arr2App)[1] != dim(arr2App)[2]) {
    # prendi la dimensione massima fra X e Y, dell'immagine
    max.dim <- max(dim(arr2App)[1:2])
    # costruisci il nuovo voxel cube quadrato rispetto ad X ed Y
    sq.m <- array( 0 , dim = c( max.dim , max.dim, dim(arr2App)[3] ))
    # popola il nuovo voxelcube(quadrato) con i dati del vecchio
    sq.m[ 1:(dim(arr2App)[1]) , 1:(dim(arr2App)[2]), 1:(dim(arr2App)[3])  ] <- arr2App
    
    # per ogni slice, riempi i voxel aggiunti, ora a zero, con i valori del confine più prossimo
    for( z.slice in 1:dim(arr2App)[3]) {
      if(dim(arr2App)[1] > dim(arr2App)[2]) {
        sq.m[ , (dim(arr2App)[2]):max.dim ,  z.slice  ] <- rep(  arr2App[ , dim(arr2App)[2], z.slice ] , length((dim(arr2App)[2]):max.dim))
      }
      else {
        sq.m[ (dim(arr2App)[1]):max.dim ,  ,  z.slice  ] <- rep(  arr2App[ dim(arr2App)[1] ,, z.slice ] , length((dim(arr2App)[1]):max.dim))
      }
    }
    arr2App <- sq.m
  }
  
  kern2Apply<-convKernel(  sigma = sigma, k = kernel.type);
  ret <- applyFilter(x =arr2App, kernel = kern2Apply);
  if(original.dim[1] != original.dim[2]) {
    ret <- ret[ 1:original.dim[1] , 1:original.dim[2] , ]
  }
  return(ret)
}
#' DEPRECATED!!!!:  Apply a filter to mmButo object
#' 
#' @description  It applies filter to an mButo object to simplify RADIOMICS issues
#' @param obj.mmButo is the \code{new.mmButo} object
#' @param ROINameForNormalization the name of the ROI used to normalize the signal
#' @param valueForNormalization the highest value used to normalize the 'ROINameForNormalization' voxel values
#' @param ROIName the name of the ROI to extract
#' @param filter.pipeline is a \code{list} where is indicated the pipeline of filtering to be applied (in sequence)
#' @param collection is the collection: the default is '\code{default}'
#' @param cropResult is a boolean (\code{TRUE} or \code{FALSE}) which indicates if the result should be cropped or not, in order to save memory. Default is \code{TRUE}
#' @param scaleFactor can be 'voxel' or 'space'. If 'voxel' (the default) it consider the sam sigma values (if specified) for all the geoLet object, if 'space' it normalize the sigma according to the different pixelSpacing.
#' @return a list containing the filtered images (cropped or not)
#' @useDynLib moddicomV2
#' @examples \dontrun{
#' # DEPRECATED!!!!!
#' } 
FIL.applyFilterToStudy<-function(obj.mmButo, ROINameForNormalization=NA, valueForNormalization = NA, ROIName, filter.pipeline ,collection="default",cropResult=TRUE, scaleFactor="voxel") {
  objS<-services();
  FUNMap<-list();
  cat("\n\n ------------------------------------------------- ")
  cat("\n funzione da testare a causa della nuova interpolazione del getROIVoxels()")
  cat("\n ------------------------------------------------- \n")
  stop();
  if(scaleFactor != "voxel" & scaleFactor != "space") stop("\n Error: 'scaleFactor' can be only 'voxel' or 'space'.");
  # prendi i pixelSpacing
  pixelSpacingArr<-obj.mmButo$getAllPixelSpacing();
  
  # prendi la ROI
  ROIVoxelData<-obj.mmButo$getROIVoxel(ROIName = ROIName)
  
  if(!is.na(ROINameForNormalization)) {
    # prendi le ROI per la normalizzazione e calcolane le statistiche
    ROIVoxelDataForNormalization<-obj.mmButo$getROIVoxel(ROIName = ROINameForNormalization)
    statsForNormalization<-obj.mmButo$getROIVoxelStats(ROIVoxelList = ROIVoxelDataForNormalization)

    # se non è stato passato un valore di normalizzazione prendi il massimo delle medie della ROI
    # utilizzata per normalizzare
    if(is.na(valueForNormalization)) {
      arrayNormalizationROIVoxelValue<-c();
      for ( patient in names(statsForNormalization$details) ) {
        arrayNormalizationROIVoxelValue<-c(arrayNormalizationROIVoxelValue,statsForNormalization$details[[ patient ]]$mean)
      }  
      valueForNormalization<-max(arrayNormalizationROIVoxelValue)
    }
    # calcola l'array dei fattori moltiplicativi
    multiplier<-list();
    for ( patient in names(statsForNormalization$details) ) {
      multiplier[[ patient ]] <- valueForNormalization / statsForNormalization$details[[ patient ]]$mean 
    }  
  } 
  else {
    # calcola l'array dei fattori moltiplicativi nel caso in cui non si voglia normalizzazione 
    # (piazzali tutti a uno)
    multiplier<-list();
    for ( patient in names(ROIVoxelData)) {
      multiplier[[ patient ]] <- 1
    }      
  } 

  # prendi la lista di oggetti geoLet
  list_geoLet<-obj.mmButo$getAttribute("list_geoLet");
  
  # inizia a ciclare su tutti i pazienti
  for(patient in names(ROIVoxelData)) {
    
    if((TRUE %in% is.na(ROIVoxelData[[patient]])) == FALSE  ) {
      print( paste("FUN - Now filtering:",patient),collapse='' )
      # prendi i voxelData ed estendili
      voxelData.ready<-ROIVoxelData[[patient]]
      voxelData.ready.espanso <- obj.mmButo$expandCube(   voxelData.ready   ) 
      voxelDataDaRiapplicare <- voxelData.ready.espanso;
      
      # e crea la relativa maschera
      #voxelData.ready.espanso[voxelData.ready.espanso!=0]<-1
      voxelData.ready.espanso[!is.na(voxelData.ready.espanso)]<-1
      
      # prendi l'original voxelCute non tagliato dalla ROI
      # è su questo che dovrò applicare il filtro
      originalMR <- list_geoLet[[collection]][[patient]]$getImageVoxelCube()
      # ora moltiplica l'MR originale per il valore di normalizzazione
      originalMR <- originalMR * multiplier[[ patient ]] 
      for(i in seq(1,length(filter.pipeline))) {
        print( paste("     => applying",filter.pipeline[[i]]$kernel.type),collapse='' )
        if (scaleFactor == "space") {
          normalizedSigma<-sqrt((filter.pipeline[[i]]$sigma^2)/(pixelSpacingArr[[patient]][1]*pixelSpacingArr[[patient]][2]))
        } else {
          normalizedSigma<-filter.pipeline[[i]]$sigma
        }
        originalMR<-FIL.applyFilter( originalMR, 
                                     kernel.type = filter.pipeline[[i]]$kernel.type,
                                     sigma = normalizedSigma)
                                     #sigma = filter.pipeline[[i]]$sigma)
      }
      # l'output deve essere mascherato
      u<-voxelData.ready.espanso
      
      u[which(!is.na(u),arr.ind = TRUE)]<-1
      u[which(is.na(u),arr.ind = TRUE)]<-0
      
      #FUNMap[[patient]]<-originalMR * voxelData.ready.espanso;
      # Applica la maschera

      FUNMap[[patient]]<-originalMR * u;
      FUNMap[[patient]][which(is.na(voxelData.ready.espanso),arr.ind = TRUE)]<-NA
      
      # se richiesto, CROPPA!
      if ( cropResult == TRUE )  FUNMap[[patient]]<-objS$cropCube( FUNMap[[patient]] )   
    }
    else {
      FUNMap[[patient]]<-NA
    }
  }
  return(FUNMap)
}
#' a new kind of filtering
#' 
#' @description  It applies filter to an mButo object to simplify RADIOMICS issues
#' @param obj.mmButo is the \code{new.mmButo} object
#' @param ROINameForNormalization the name of the ROI used to normalize the signal
#' @param valueForNormalization the highest value used to normalize the 'ROINameForNormalization' voxel values
#' @param ROIName the name of the ROI to extract
#' @param filter.pipeline is a \code{list} where is indicated the pipeline of filtering to be applied (in sequence)
#' @param collection is the collection: the default is '\code{default}'
#' @param cropResult is a boolean (\code{TRUE} or \code{FALSE}) which indicates if the result should be cropped or not, in order to save memory. Default is \code{TRUE}
#' @param scaleFactor can be 'voxel' or 'space'. If 'voxel' (the default) it consider the sam sigma values (if specified) for all the geoLet object, if 'space' it normalize the sigma according to the different pixelSpacing.
#' @param ROINameForNormalization Optional, default \code{NA}. If you want to normalize the voxelCubes accoring to some specific ROIs (i.e. Urin) you can indicate here the ROIName. By this way, the other parameter \code{upperValueForNormalization}. Signals into such ROIS will be aligned at the biggest.
#' @return a list containing the filtered images (cropped or not)
#' @export
#' @useDynLib moddicomV2 
#' @examples \dontrun{
#' # Create an instante of new.mmButo and load some cases
#' obj<-new.mmButo()
#' obj$loadCollection("/progetti/immagini/urinaEasy")
#' 
#' # build the filtering pipeline
#' filterPipeline<- list()
#' filterPipeline[[1]]<-list("kernel.type"="LoG", "sigma"=1.5)
#'
#' # filter the images (normalizing the GTV signal with Urin) cropping the result
#' a<-FIL.2D.conv.Filter(obj.mmButo = obj,ROIName = "GTV",ROINameForNormalization='Urina', valueForNormalization=10000,filter.pipeline = filterPipeline )
#'
#' # if you want to normalize with the signal in urins, for example:
#' a<-FIL.2D.conv.Filter(obj.mmButo = obj,ROIName = "GTV",ROINameForNormalization='Urina', valueForNormalization=10000,filter.pipeline = filterPipeline , ROINameForNormalization = "Urina")
#' }
FIL.2D.conv.Filter<-function(obj.mmButo,ROIName, filter.pipeline ,collection="default",cropResult=TRUE, scaleFactor="voxel", cropImageVoxelCube_inProcessing = TRUE, ROINameForNormalization=NA) {
  objS<-services();
  FUNMap<-list();
  upperValueForNormalization<-NA
  cat("\n\n ------------------------------------------------- ")
  cat("\n funzione da testare a causa della nuova interpolazione del getROIVoxels()")
  cat("\n ------------------------------------------------- \n")
  stop();  
  listaValoriMedi<-list();
  if(scaleFactor != "voxel" & scaleFactor != "space") stop("\n Error: 'scaleFactor' can be only 'voxel' or 'space'.");
  # prendi i pixelSpacing
  pixelSpacingArr<-obj.mmButo$getAllPixelSpacing();
  
  # prendi la ROI
  ROIVoxelData<-obj.mmButo$getROIVoxel(ROIName = ROIName)
  # se bisogna normalizzare, prendi anche le ROI di normalizzazione
  if(!is.na(ROINameForNormalization)) {
    norma.ROI<-obj.mmButo$getROIVoxel(ROIName = ROINameForNormalization)
    # cicla su tutte 
    for(i in names(norma.ROI)) {
      listaValoriMedi[[i]]<-mean(norma.ROI[[i]]$masked.images$voxelCube[ which(!is.na(  norma.ROI[[i]]$masked.images$voxelCube ),arr.ind = TRUE ) ])
    }
    upperValueForNormalization<-max(  unlist( listaValoriMedi ) )
    for(i in names(listaValoriMedi)) {
      listaValoriMedi[[i]]<- upperValueForNormalization / listaValoriMedi[[i]]
    }
  }

  # prendi la lista di oggetti geoLet
  list_geoLet<-obj.mmButo$getAttribute("list_geoLet");
  
  # inizia a ciclare su tutti i pazienti
  for(patient in names(ROIVoxelData)) {
    
    if((TRUE %in% is.na(ROIVoxelData[[patient]])) == FALSE  ) {
      print( paste("FUN - Now filtering:",patient),collapse='' )
      
      # prendi i voxelData ed estendili
      voxelData.ready<-ROIVoxelData[[patient]]
      voxelData.ready.espanso <- obj.mmButo$expandCube(   voxelData.ready   ) 
      voxelDataDaRiapplicare <- voxelData.ready.espanso;
      
      # e crea la relativa maschera
      #voxelData.ready.espanso[voxelData.ready.espanso!=0]<-1
      voxelData.ready.espanso[!is.na(voxelData.ready.espanso)]<-1
      
      # prendi l'original voxelCute non tagliato dalla ROI
      # è su questo che dovrò applicare il filtro
      originalMR <- list_geoLet[[collection]][[patient]]$getImageVoxelCube()
      
      # se è stato chiesto di normalizzare, normalizza
      if(!is.na(ROINameForNormalization)) {
        if(min(originalMR)<0) cat("\n WARNING: negative values... are you sure you want to normalize???")
        originalMR<- originalMR * listaValoriMedi[[patient]]
      }      
      
      minZ<-voxelData.ready$masked.images$location$min.z
      maxZ<-voxelData.ready$masked.images$location$max.z
      
      # se e' stato chiesto di croppare il voxelcube per migliorere le pervormance
      # (ha senso sul 2d!)
      if(cropImageVoxelCube_inProcessing == TRUE) {
        croppedSubMatrix<- originalMR[,,minZ:maxZ]
        OLDoriginalMR<-originalMR
        originalMR<-croppedSubMatrix
      }
      
      for(i in seq(1,length(filter.pipeline))) {
        print( paste("     => applying",filter.pipeline[[i]]$kernel.type),collapse='' )
        if (scaleFactor == "space") {
          normalizedSigma<-sqrt((filter.pipeline[[i]]$sigma^2)/(pixelSpacingArr[[patient]][1]*pixelSpacingArr[[patient]][2]))
        } else {
          normalizedSigma<-filter.pipeline[[i]]$sigma
        }
        originalMR<-FIL.applyFilter( originalMR, 
                                     kernel.type = filter.pipeline[[i]]$kernel.type,
                                     sigma = normalizedSigma)
        #sigma = filter.pipeline[[i]]$sigma)
      }
      # se era stato croppato, ripristina!
      if(cropImageVoxelCube_inProcessing == TRUE) {
        for(zPos in seq(minZ,maxZ)) {
          OLDoriginalMR[,,zPos]<-originalMR[,,(zPos-minZ+1)]
        }
        originalMR<-OLDoriginalMR
      }
      
      # l'output deve essere mascherato
      u<-voxelData.ready.espanso
      
      u[which(!is.na(u),arr.ind = TRUE)]<-1
      u[which(is.na(u),arr.ind = TRUE)]<-0
      
      # Applica la maschera
      FUNMap[[patient]]<-originalMR * u;
      FUNMap[[patient]][which(is.na(voxelData.ready.espanso),arr.ind = TRUE)]<-NA
      
      # se richiesto, CROPPA!
      if ( cropResult == TRUE )  FUNMap[[patient]]<-objS$cropCube( FUNMap[[patient]] )   
    }
    else {
      FUNMap[[patient]]<-NA
    }
  }
  return(FUNMap)
}
#' a new kind of filtering able to work on geoLet objects
#' 
#' @description  It applies filter to an geoLet object to simplify RADIOMICS issues
#' @param obj.geoLet is the \code{new.mmButo} object
#' @param ROINameForNormalization the name of the ROI used to normalize the signal
#' @param valueForNormalization the highest value used to normalize the 'ROINameForNormalization' voxel values
#' @param ROIName the name of the ROI to extract
#' @param filter.pipeline is a \code{list} where is indicated the pipeline of filtering to be applied (in sequence)
#' @param collection is the collection: the default is '\code{default}'
#' @param cropResult is a boolean (\code{TRUE} or \code{FALSE}) which indicates if the result should be cropped or not, in order to save memory. Default is \code{TRUE}
#' @param scaleFactor can be 'voxel' or 'space'. If 'voxel' (the default) it consider the sam sigma values (if specified) for all the geoLet object, if 'space' it normalize the sigma according to the different pixelSpacing.
#' @param ROINameForNormalization Optional, default \code{NA}. If you want to normalize the voxelCubes accoring to some specific ROIs (i.e. Urin) you can indicate here the ROIName. By this way, the other parameter \code{upperValueForNormalization}. Signals into such ROIS will be aligned at the biggest.
#' @return a list containing the filtered images (cropped or not)
#' @export
#' @useDynLib moddicomV2 
#' @examples \dontrun{
#' # Create an instante of new.mmButo and load some cases
#' obj<-geoLet()
#' obj$openDICOMFolder("/progetti/immagini/urinaEasy")
#' 
#' # build the filtering pipeline
#' filterPipeline<- list()
#' filterPipeline[[1]]<-list("kernel.type"="LoG", "sigma"=1.5)
#'
#' # filter the images (normalizing the GTV signal with Urin) cropping the result
#' a<-FIL.2D.conv.Filter(obj.geoLet = obj,ROIName = "GTV",ROINameForNormalization='Urina', valueForNormalization=10000,filter.pipeline = filterPipeline )
#'
#' # if you want to normalize with the signal in urins, for example:
#' a<-FIL.2D.conv.Filter(obj.geoLet = obj,ROIName = "GTV",ROINameForNormalization='Urina', valueForNormalization=10000,filter.pipeline = filterPipeline , ROINameForNormalization = "Urina")
#' } 
FIL.2D.conv.Filter.geoLet<-function(obj.geoLet,ROIName, filter.pipeline ,collection="default",cropResult=TRUE, 
                                    scaleFactor="voxel", cropImageVoxelCube_inProcessing = TRUE, 
                                    ROINameForNormalization=NA, px ="", py="", returnBigCube = FALSE) {
  objS<-services();
  FUNMap<-list();

  upperValueForNormalization<-NA
  listaValoriMedi<-list();
  if(scaleFactor != "voxel" & scaleFactor != "space") stop("\n Error: 'scaleFactor' can be only 'voxel' or 'space'.");
  # prendi i pixelSpacing
#  pixelSpacingArr<-obj.mmButo$getAllPixelSpacing();

  pixelSpacingArr<-obj.geoLet$getPixelSpacing();
  if(px=="") px <- pixelSpacingArr[1]
  if(py=="") py <- pixelSpacingArr[2]
  
#   # prendi la ROI
  # browser()
  new.pixelSpacing <- c( px , py)
  ROIVoxelData<-obj.geoLet$getROIVoxel(Structure = ROIName, new.pixelSpacing = new.pixelSpacing)

#   # se bisogna normalizzare, prendi anche le ROI di normalizzazione
#   if(!is.na(ROINameForNormalization)) {
#     norma.ROI<-obj.mmButo$getROIVoxel(ROIName = ROINameForNormalization)
#     # cicla su tutte 
#     for(i in names(norma.ROI)) {
#       listaValoriMedi[[i]]<-mean(norma.ROI[[i]]$masked.images$voxelCube[ which(!is.na(  norma.ROI[[i]]$masked.images$voxelCube ),arr.ind = TRUE ) ])
#     }
#     upperValueForNormalization<-max(  unlist( listaValoriMedi ) )
#     for(i in names(listaValoriMedi)) {
#       listaValoriMedi[[i]]<- upperValueForNormalization / listaValoriMedi[[i]]
#     }
#   }
  
  # prendi la lista di oggetti geoLet
#  list_geoLet<-obj.mmButo$getAttribute("list_geoLet");
  
  # inizia a ciclare su tutti i pazienti
#  for(patient in names(ROIVoxelData)) {
    
    if((TRUE %in% is.na(ROIVoxelData)) == FALSE  ) {
      print( paste("FUN - Now filtering:"),collapse='' )
      
      # prendi i voxelData ed estendili
      voxelData.ready<-ROIVoxelData
      pc<-ROIVoxelData
      x<-pc$masked.images$location$min.x; y<-pc$masked.images$location$min.y; z<-pc$masked.images$location$min.z
      fe<-pc$masked.images$location$fe; se<-pc$masked.images$location$se;  te<-pc$masked.images$location$te      
#      voxelData.ready.espanso <- obj.mmButo$expandCube(   voxelData.ready   ) 
      # browser()
      voxelData.ready.espanso<-objS$expandCube(littleCube = pc$masked.images$voxelCube, x.start = x, y.start=y, z.start=z, fe = fe, se = se, te = te )    
      voxelDataDaRiapplicare <- voxelData.ready.espanso;

      # e crea la relativa maschera
      #voxelData.ready.espanso[voxelData.ready.espanso!=0]<-1
      voxelData.ready.espanso[!is.na(voxelData.ready.espanso)]<-1
      
      # prendi l'original voxelCute non tagliato dalla ROI
      # è su questo che dovrò applicare il filtro
      # originalMR <- obj.geoLet$getImageVoxelCube()  (  new.pixelSpacing  )
      originalMR <- obj.geoLet$getImageVoxelCube( 
        ps.x = px, 
        ps.y = py,
        ps.z = pixelSpacingArr[3]
        )
      

      minZ<-voxelData.ready$masked.images$location$min.z
      maxZ<-voxelData.ready$masked.images$location$max.z
      
      # se e' stato chiesto di croppare il voxelcube per migliorere le pervormance
      # (ha senso sul 2d!)
      if(cropImageVoxelCube_inProcessing == TRUE) {
        croppedSubMatrix<- originalMR[,,minZ:maxZ]
        OLDoriginalMR<-originalMR
        originalMR<-croppedSubMatrix
      }
      
      for(i in seq(1,length(filter.pipeline))) {
        print( paste("     => applying",filter.pipeline[[i]]$kernel.type),collapse='' )
        if (scaleFactor == "space") {
          normalizedSigma<-sqrt((filter.pipeline[[i]]$sigma^2)/(pixelSpacingArr[1]*pixelSpacingArr[2]))
        } else {
          normalizedSigma<-filter.pipeline[[i]]$sigma
        }
        originalMR<-FIL.applyFilter( originalMR, 
                                     kernel.type = filter.pipeline[[i]]$kernel.type,
                                     sigma = normalizedSigma)
        #sigma = filter.pipeline[[i]]$sigma)
      }
      # se era stato croppato, ripristina!
      if(cropImageVoxelCube_inProcessing == TRUE) {
        # browser()
        for(zPos in seq(minZ,maxZ)) {
          OLDoriginalMR[,,zPos]<-originalMR[,,(zPos-minZ+1)]
        }
        originalMR<-OLDoriginalMR
      }
      
      # l'output deve essere mascherato
      u<-voxelData.ready.espanso
      
      u[which(!is.na(u),arr.ind = TRUE)]<-1
      u[which(is.na(u),arr.ind = TRUE)]<-0
      
      # Applica la maschera
      FUNMap<-originalMR * u;
      FUNMap[which(is.na(voxelData.ready.espanso),arr.ind = TRUE)]<-NA
      
      # browser()
      # se richiesto, CROPPA!
      if ( cropResult == TRUE )  { FUNMap<-objS$cropCube( FUNMap ) }
    }
    else {
      FUNMap<-NA
    }
#  }
  if(returnBigCube==TRUE) FUNMap$originalBigCube<-originalMR
  return(FUNMap)
}
old.FIL.2D.conv.Filter.geoLet<-function(obj.geoLet,ROIName, filter.pipeline ,collection="default",cropResult=TRUE, 
                                        scaleFactor="voxel", cropImageVoxelCube_inProcessing = TRUE, 
                                        ROINameForNormalization=NA, px ="", py="") {
  objS<-services();
  FUNMap<-list();
  cat("\n\n ------------------------------------------------- ")
  cat("\n funzione da testare a causa della nuova interpolazione del getROIVoxels()")
  cat("\n ------------------------------------------------- \n")
  # stop();  
  upperValueForNormalization<-NA
  listaValoriMedi<-list();
  if(scaleFactor != "voxel" & scaleFactor != "space") stop("\n Error: 'scaleFactor' can be only 'voxel' or 'space'.");
  # prendi i pixelSpacing
  #  pixelSpacingArr<-obj.mmButo$getAllPixelSpacing();
  browser()
  
  pixelSpacingArr<-obj.geoLet$getPixelSpacing();
  
  
  #   # prendi la ROI
  # browser()
  ROIVoxelData<-obj.geoLet$getROIVoxel(Structure = ROIName)
  #   # se bisogna normalizzare, prendi anche le ROI di normalizzazione
  #   if(!is.na(ROINameForNormalization)) {
  #     norma.ROI<-obj.mmButo$getROIVoxel(ROIName = ROINameForNormalization)
  #     # cicla su tutte 
  #     for(i in names(norma.ROI)) {
  #       listaValoriMedi[[i]]<-mean(norma.ROI[[i]]$masked.images$voxelCube[ which(!is.na(  norma.ROI[[i]]$masked.images$voxelCube ),arr.ind = TRUE ) ])
  #     }
  #     upperValueForNormalization<-max(  unlist( listaValoriMedi ) )
  #     for(i in names(listaValoriMedi)) {
  #       listaValoriMedi[[i]]<- upperValueForNormalization / listaValoriMedi[[i]]
  #     }
  #   }
  
  # prendi la lista di oggetti geoLet
  #  list_geoLet<-obj.mmButo$getAttribute("list_geoLet");
  
  # inizia a ciclare su tutti i pazienti
  #  for(patient in names(ROIVoxelData)) {
  
  if((TRUE %in% is.na(ROIVoxelData)) == FALSE  ) {
    print( paste("FUN - Now filtering:"),collapse='' )
    
    # prendi i voxelData ed estendili
    voxelData.ready<-ROIVoxelData
    pc<-ROIVoxelData
    x<-pc$masked.images$location$min.x; y<-pc$masked.images$location$min.y; z<-pc$masked.images$location$min.z
    fe<-pc$masked.images$location$fe; se<-pc$masked.images$location$se;  te<-pc$masked.images$location$te      
    #      voxelData.ready.espanso <- obj.mmButo$expandCube(   voxelData.ready   ) 
    voxelData.ready.espanso<-objS$expandCube(littleCube = pc$masked.images$voxelCube, x.start = x, y.start=y, z.start=z, fe = fe, se = se, te = te )    
    voxelDataDaRiapplicare <- voxelData.ready.espanso;
    
    # e crea la relativa maschera
    #voxelData.ready.espanso[voxelData.ready.espanso!=0]<-1
    voxelData.ready.espanso[!is.na(voxelData.ready.espanso)]<-1
    
    # prendi l'original voxelCute non tagliato dalla ROI
    # è su questo che dovrò applicare il filtro
    originalMR <- obj.geoLet$getImageVoxelCube()
    
    
    minZ<-voxelData.ready$masked.images$location$min.z
    maxZ<-voxelData.ready$masked.images$location$max.z
    
    # se e' stato chiesto di croppare il voxelcube per migliorere le pervormance
    # (ha senso sul 2d!)
    if(cropImageVoxelCube_inProcessing == TRUE) {
      croppedSubMatrix<- originalMR[,,minZ:maxZ]
      OLDoriginalMR<-originalMR
      originalMR<-croppedSubMatrix
    }
    
    for(i in seq(1,length(filter.pipeline))) {
      print( paste("     => applying",filter.pipeline[[i]]$kernel.type),collapse='' )
      if (scaleFactor == "space") {
        normalizedSigma<-sqrt((filter.pipeline[[i]]$sigma^2)/(pixelSpacingArr[1]*pixelSpacingArr[2]))
      } else {
        normalizedSigma<-filter.pipeline[[i]]$sigma
      }
      originalMR<-FIL.applyFilter( originalMR, 
                                   kernel.type = filter.pipeline[[i]]$kernel.type,
                                   sigma = normalizedSigma)
      #sigma = filter.pipeline[[i]]$sigma)
    }
    # se era stato croppato, ripristina!
    if(cropImageVoxelCube_inProcessing == TRUE) {
      for(zPos in seq(minZ,maxZ)) {
        OLDoriginalMR[,,zPos]<-originalMR[,,(zPos-minZ+1)]
      }
      originalMR<-OLDoriginalMR
    }
    
    # l'output deve essere mascherato
    u<-voxelData.ready.espanso
    
    u[which(!is.na(u),arr.ind = TRUE)]<-1
    u[which(is.na(u),arr.ind = TRUE)]<-0
    
    # Applica la maschera
    FUNMap<-originalMR * u;
    FUNMap[which(is.na(voxelData.ready.espanso),arr.ind = TRUE)]<-NA
    
    # browser()
    # se richiesto, CROPPA!
    if ( cropResult == TRUE )  { FUNMap<-objS$cropCube( FUNMap ) }
  }
  else {
    FUNMap<-NA
  }
  #  }
  return(FUNMap)
}
#' a new kind of filtering able to work on geoLet objects
#' 
#' @description  It applies filter to an geoLet object to simplify RADIOMICS issues
#' @param obj.geoLet is the \code{new.mmButo} object
#' @param ROINameForNormalization the name of the ROI used to normalize the signal
#' @param valueForNormalization the highest value used to normalize the 'ROINameForNormalization' voxel values
#' @param ROIName the name of the ROI to extract
#' @param filter.pipeline is a \code{list} where is indicated the pipeline of filtering to be applied (in sequence)
#' @param collection is the collection: the default is '\code{default}'
#' @param cropResult is a boolean (\code{TRUE} or \code{FALSE}) which indicates if the result should be cropped or not, in order to save memory. Default is \code{TRUE}
#' @param scaleFactor can be 'voxel' or 'space'. If 'voxel' (the default) it consider the sam sigma values (if specified) for all the geoLet object, if 'space' it normalize the sigma according to the different pixelSpacing.
#' @param ROINameForNormalization Optional, default \code{NA}. If you want to normalize the voxelCubes accoring to some specific ROIs (i.e. Urin) you can indicate here the ROIName. By this way, the other parameter \code{upperValueForNormalization}. Signals into such ROIS will be aligned at the biggest.
#' @return a list containing the filtered images (cropped or not)
#' @export
#' @useDynLib moddicomV2 
#' @examples \dontrun{
#' # Create an instante of new.mmButo and load some cases
#' obj<-geoLet()
#' obj$openDICOMFolder("/progetti/immagini/urinaEasy")
#' 
#' # build the filtering pipeline
#' filterPipeline<- list()
#' filterPipeline[[1]]<-list("kernel.type"="LoG", "sigma"=1.5)
#'
#' # filter the images (normalizing the GTV signal with Urin) cropping the result
#' a<-FIL.2D.conv.Filter(obj.geoLet = obj,ROIName = "GTV",ROINameForNormalization='Urina', valueForNormalization=10000,filter.pipeline = filterPipeline )
#'
#' # if you want to normalize with the signal in urins, for example:
#' a<-FIL.2D.conv.Filter(obj.geoLet = obj,ROIName = "GTV",ROINameForNormalization='Urina', valueForNormalization=10000,filter.pipeline = filterPipeline , ROINameForNormalization = "Urina")
#' }
FIL.2D.conv.Filter.voxelVolumes<-function(study.VV, ROI.VV, filter.pipeline , pixelSpacingArr , collection="default",cropResult=TRUE, scaleFactor="voxel", cropImageVoxelCube_inProcessing = TRUE, ROINameForNormalization=NA) {
  objS<-services();
  if(scaleFactor != "voxel" & scaleFactor != "space") stop("\n Error: 'scaleFactor' can be only 'voxel' or 'space'.");
  
  cat("\n\n ------------------------------------------------- ")
  cat("\n funzione da testare a causa della nuova interpolazione del getROIVoxels()")
  cat("\n ------------------------------------------------- \n")
  stop();
  
  print( paste("FUN - Now filtering:"),collapse='' )
  margine<-5
  
  min.x<-ROI.VV$masked.images$location$min.x-margine; min.y<-ROI.VV$masked.images$location$min.y-margine; min.z<-ROI.VV$masked.images$location$min.z-margine
  max.x<-ROI.VV$masked.images$location$max.x+margine; max.y<-ROI.VV$masked.images$location$max.y+margine; max.z<-ROI.VV$masked.images$location$max.z+margine
  delta.x<-ROI.VV$masked.images$location$max.x - ROI.VV$masked.images$location$min.x
  delta.y<-ROI.VV$masked.images$location$max.y - ROI.VV$masked.images$location$min.y
  delta.z<-ROI.VV$masked.images$location$max.z - ROI.VV$masked.images$location$min.z
  def.from.x<-margine; def.to.x<-(max.x-min.x);
  def.from.y<-margine; def.to.y<-(max.y-min.y);

  study.VV <- study.VV[ min.x: max.x, min.y:max.y, min.z:max.z ]

  for(i in seq(1,length(filter.pipeline))) {
    print( paste("     => applying",filter.pipeline[[i]]$kernel.type),collapse='' )
    if (scaleFactor == "space") {
      normalizedSigma<-sqrt((filter.pipeline[[i]]$sigma^2)/(pixelSpacingArr[1]*pixelSpacingArr[2]))
    } else {
      normalizedSigma<-filter.pipeline[[i]]$sigma
    }
    originalMR<-FIL.applyFilter( study.VV, 
                                 kernel.type = filter.pipeline[[i]]$kernel.type,
                                 sigma = normalizedSigma)
  }
  
  originalMR<-originalMR[ margine:(margine+delta.x),margine:(margine+delta.y),margine:(margine+delta.z) ]

  maschera <- ROI.VV$masked.images$voxelCube
  # l'output deve essere mascherato
  maschera[which(!is.na(maschera),arr.ind = TRUE)]<-1
  maschera[which(is.na(maschera),arr.ind = TRUE)]<-0
  
  # Applica la maschera
  originalMR<-originalMR * maschera;
  originalMR[which(is.na(ROI.VV$masked.images$voxelCube),arr.ind = TRUE)]<-NA
  
  # se richiesto, CROPPA!
  if ( cropResult == TRUE )  FUNMap<-objS$cropCube( originalMR )   

  return(originalMR)
}
#' Return a filtered VoxelCube
#' 
#' @description  returns a filtered voxelcube (WIP)
#' @return an list of items
#' @export
FIL.getFilteredVoxelCube<-function( obj.geoLet, ROIName, kind.of.filter, sigma=0.8 ) {
  if( kind.of.filter != "LoG" & kind.of.filter!="none") stop("Error: only 'LoG' and 'none' are admitted, in this release")

  cat("\n\n ------------------------------------------------- ")
  cat("\n funzione da testare a causa della nuova interpolazione del getROIVoxels()")
  cat("\n ------------------------------------------------- \n")
  stop();
  
  if(kind.of.filter == "LoG") {
    filterPipeline<- list()
    filterPipeline[[1]]<-list("kernel.type"="LoG", "sigma"=sigma)
    voxelCube <- FIL.2D.conv.Filter.geoLet(obj.geoLet = obj.geoLet,ROIName = ROIName,
                                           filter.pipeline = filterPipeline,
                                           cropResult = TRUE,scaleFactor = 'space',
                                           cropImageVoxelCube_inProcessing = TRUE)
  }
  if(kind.of.filter == "none") {
    extracted.VC <-  obj.geoLet$getROIVoxels(Structure = ROIName)
    voxelCube<-extracted.VC$masked.images
  }
  return( list( "voxelCube"=voxelCube  )   )
}
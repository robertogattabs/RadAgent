#' class for doing Radiomics
#' 
#' @description  A class to do Radiomics
#' @export
radiomics<-function(){
  DICOMFolder <- ""
  #=================================================================================
  # setDICOMFolder
  # set the folder where DICOM studies are stored
  #=================================================================================  
  setDICOMFolder<-function(dicom.path) {
    DICOMFolder<<-dicom.path
  }
  #=================================================================================
  # run
  # run a Radiomics Analysis
  #=================================================================================  
  run<-function() {
    
  }  
  class<-function() {
    return("radiomics");
  }
  #=================================================================================
  # Constructor
  #=================================================================================
  constructor<-function( ) {
    DICOMFolder <<- ""
  }  
  constructor( )
  return(list(
    "setDICOMFolder" = setDICOMFolder
  ))  
} 





#' function to calculate the first order features
#' 
#' @description  Calculates Shannon entropy, kursosis, Skewness, mean, standard deviation and energy of a given list of arrays of voxels
#' @param inputData a list where each element is an array of the voxel of the image. Each element of the list normally refers to a patient.             
#' @return six lists: the first list contains the entropies, the second the kurtosis, the third the skewness, the fourth the mean, the fifth the standard deviation and the sixth the energy
#' @export
#' @examples \dontrun{
#' # Create an instante of mmButo and load some cases
#' obj<-mmButo()
#' obj$loadCollection(Path = '/progetti/immagini/urinaEasy')
#' 
#' # get the three ROIs
#' Retto<-obj$getROIVoxel(ROIName="Retto")  
#' 
#' # get the possible biopsy
#' aa<-RAD.firstOrderFeatureImage(inputData = Retto )
#' aa$entropy
#' }#' #' 
#' @import entropy moments 
RAD.firstOrderFeatureImage <- function ( inputData ) { 
  # NON USATA, DA RANZARE?
  if(class(inputData) == "geoLetStructureVoxelList") 
    return(RAD.firstOrderFeatureImage.geoLet( inputData))
  if(class(inputData) == "mmButoStructureVoxelList")  
    return(RAD.firstOrderFeatureImage.mmButo( inputData))
}
RAD.firstOrderFeatureImage.geoLet<-function(inputData) {
  # NON USATA, DA RANZARE?
  numPatient<-1
  ImageEntropy <- array(data = c(0), dim = c(numPatient))
  ImageKurtosis <- array(data = c(0), dim = c(numPatient))
  ImageSkewness <- array(data = c(0), dim = c(numPatient))
  ImageMean <- array(data = c(0), dim = c(numPatient))
  ImageStandDeviat <- array(data = c(0), dim = c(numPatient))

  istogr <- c();    freq <- c();  i<-1;
  if(is.list(inputData) == TRUE  ) {

    histSamples<-500

    #voxelCube.values<-inputData$masked.images$voxelCube[ inputData$masked.images$voxelCube!=0  ]
    voxelCube.values<-inputData$masked.images$voxelCube[ !is.na(inputData$masked.images$voxelCube)  ] 
    
    maxVoxelValue<-max(voxelCube.values);  minVoxelValue<-min(voxelCube.values);
    histSamples.array<-seq( from = minVoxelValue, to=maxVoxelValue, by = (maxVoxelValue-minVoxelValue)/histSamples   )
    
    # Calcola l'istogramma dei grigi
    #istogr <- hist(voxelCube.values, breaks = histSamples,  plot=FALSE)
    istogr <- discretize(x = voxelCube.values[which(!is.na(voxelCube.values),arr.ind = TRUE ) ], numBins = histSamples)        
    # Calcola le frequenze da dover utilizzare nel calcolo dell'entropia
    #freq <- freqs(y = istogr$counts)
    # Calcola l'entropia di Shannon per ogni paziente
    #ImageEntropy[i] <- entropy.plugin(freqs = freq, unit = c("log2"))
    ImageEntropy[i] <- entropy.plugin(freqs = istogr, unit = c("log2"))
    # Calcola la Kurtosis per ogni paziente (Kurtosis=0 distribuzione normale, Kurtosis > 0 distribuzione leptocurtica
    # cioè stretta, Kurtosis < 0 distribuzione platicurtica cioè larga)
    ImageKurtosis[i] <- kurtosis (x = voxelCube.values[which(!is.na(voxelCube.values),arr.ind = TRUE ) ])
    # Calcola la Skewness per ogni paziente (Skewness = 0 simmetria perfetta, Skewness > 0 asimmetrica verso destra
    # Skewness < 0 asimmetrica verso sinistra)
    ImageSkewness[i] <- skewness(x = voxelCube.values[which(!is.na(voxelCube.values),arr.ind = TRUE ) ])
    #Calcola media dei grigi dei singoli pazienti
    ImageMean[i] <- mean(x = voxelCube.values[which(!is.na(voxelCube.values),arr.ind = TRUE ) ])
    #Calcola deviazione standard
    ImageStandDeviat[i] <- sd (x = voxelCube.values[which(!is.na(voxelCube.values),arr.ind = TRUE ) ])
    
  } else {
    ImageEntropy[i] <-NA;
    ImageKurtosis[i]<-NA;
    ImageSkewness[i]<-NA;
    ImageMean[i]<-NA;
    ImageStandDeviat[i]<-NA;
  }  
  #Restituisce una lista di array; ciascun valore corrisponde al singolo paziente ()
  return(list ("entropy"=ImageEntropy, "kurtosis"=ImageKurtosis, "skewness"=ImageSkewness, "mean"=ImageMean, 
               "standardDeviation"=ImageStandDeviat))   
}
RAD.firstOrderFeatureImage.mmButo <- function ( inputData ) {
  
  # set some variables;
  numPatient<-length(inputData)
  obj.mButo<-mmButo()
  ImageEntropy <- array(data = c(0), dim = c(numPatient))
  ImageKurtosis <- array(data = c(0), dim = c(numPatient))
  ImageSkewness <- array(data = c(0), dim = c(numPatient))
  ImageMean <- array(data = c(0), dim = c(numPatient))
  ImageStandDeviat <- array(data = c(0), dim = c(numPatient))

  histSamples<-500
  
  voxel.stats<-obj.mButo$getROIVoxelStats( inputData )
  # stop("NOT YET SUPPORTED #dnf9df89")
  maxVoxelValue<-max(voxel.stats$summary$max)
  minVoxelValue<-min(voxel.stats$summary$min)
  histSamples.array<-seq( from = minVoxelValue, to=maxVoxelValue, by = (maxVoxelValue-minVoxelValue)/histSamples   )
  
  # loop on each patient
  for (i in 1:numPatient)  {
    istogr <- c();    freq <- c()
    if((TRUE %in% is.na(inputData[[i]])) == FALSE  ) {
      voxelCube.values<-unlist(inputData[[i]]$masked.images$voxelCube)
      voxelCube.values<-voxelCube.values[ !is.na(voxelCube.values)  ] 
      # Calcola l'istogramma dei grigi
      #istogr <- hist(voxelCube.values, breaks = histSamples,  plot=FALSE)
      
      #istogr <- discretize(x = voxelCube.values, numBins = histSamples)
      istogr <- discretize(x = voxelCube.values[which(!is.na(voxelCube.values),arr.ind = TRUE ) ], numBins = histSamples)        
      # Calcola le frequenze da dover utilizzare nel calcolo dell'entropia
      #freq <- freqs(y = istogr$counts)
      # Calcola l'entropia di Shannon per ogni paziente
      #ImageEntropy[i] <- entropy.plugin(freqs = freq, unit = c("log2"))
      ImageEntropy[i] <- entropy.plugin(freqs = istogr, unit = c("log2"))
      # Calcola la Kurtosis per ogni paziente (Kurtosis=0 distribuzione normale, Kurtosis > 0 distribuzione leptocurtica
      # cioè stretta, Kurtosis < 0 distribuzione platicurtica cioè larga)
      ImageKurtosis[i] <- kurtosis (x = voxelCube.values[which(!is.na(voxelCube.values),arr.ind = TRUE ) ])
      # Calcola la Skewness per ogni paziente (Skewness = 0 simmetria perfetta, Skewness > 0 asimmetrica verso destra
      # Skewness < 0 asimmetrica verso sinistra)
      ImageSkewness[i] <- skewness(x = voxelCube.values[which(!is.na(voxelCube.values),arr.ind = TRUE ) ])
      #Calcola media dei grigi dei singoli pazienti
      ImageMean[i] <- mean(x = voxelCube.values[which(!is.na(voxelCube.values),arr.ind = TRUE ) ])
      #Calcola deviazione standard
      ImageStandDeviat[i] <- sd (x = voxelCube.values[which(!is.na(voxelCube.values),arr.ind = TRUE ) ])
      

    } else {
      ImageEntropy[i] <-NA;
      ImageKurtosis[i]<-NA;
      ImageSkewness[i]<-NA;
      ImageMean[i]<-NA;
      ImageStandDeviat[i]<-NA;
    }
  }
  #Restituisce una lista di array; ciascun valore corrisponde al singolo paziente ()
  return(list ("entropy"=ImageEntropy, "kurtosis"=ImageKurtosis, "skewness"=ImageSkewness, "mean"=ImageMean, 
               "standardDeviation"=ImageStandDeviat)) 
}


#' function to calculate Are/Volume and related measures
#' 
#' @description  calculates Are, Volume, Area/Volume Ratio and equivolumetric Spherical Area Ratio
#' @param listaROIVoxels an output of a \code{obj$getROIVoxel()} method
#' @return a list containing, for each patient the indicated measures
#' @export
#' @examples \dontrun{
#' # Create an instante of mmButo and load some cases
#' obj<-mmButo()
#' obj$loadCollection(Path = '/progetti/immagini/urinaEasy')
#' 
#' # get the three ROIs
#' Retto<-obj$getROIVoxel(ROIName="Retto")  
#' 
#' # get the possible biopsy
#' uu<-RAD.areaVolume(listaROIVoxels = Retto)
#' }#' #' 
RAD.areaVolume <- function ( listaROIVoxels ) {
  if(class(listaROIVoxels) == "geoLetStructureVoxelList") 
    return(RAD.areaVolume.geoLet( listaROIVoxels ))
  if(class(listaROIVoxels) == "mmButoStructureVoxelList")  
    return(RAD.areaVolume.mmButo( listaROIVoxels ))
}
RAD.areaVolume.geoLet<-function( listaROIVoxels ) {
  objS<-services();  arrayAV<-list(); obj.mmButo<-mmButo(); i<-1;
  if(is.list(listaROIVoxels) == TRUE  ) {
    arrayAV[[ i ]]<-list();
    geometry<-listaROIVoxels$geometricalInformationOfImages;
    pSX<-geometry$pixelSpacing[1]
    pSY<-geometry$pixelSpacing[2]
    pSZ<-as.numeric(geometry$SliceThickness  )
    # expand the cropped voxelCube
    voxelCube <- listaROIVoxels$masked.images$voxelCube
    arrayAV[[ i ]]$Area<-objS$rawSurface(voxelMatrix = voxelCube, pSX = pSX, pSY=pSY,pSZ=pSZ)    
    if ( arrayAV[[ i ]]$Area == -1 ) {
      arrayAV[[ i ]]$Volume<- -1
      arrayAV[[ i ]]$equivolumetricSphericAreaRatio<- -1
    }
    else {
      arrayAV[[ i ]]$Volume<-length(which(voxelCube!=0))*pSX*pSY*pSZ
      arrayAV[[ i ]]$equivolumetricSphericAreaRatio<- ( 4*pi* (   (3/(4*pi))*arrayAV[[ i ]]$Volume   )^(2/3) ) / arrayAV[[ i ]]$Area
    }
  } else {
    arrayAV[[ i ]]$Volume<-NA
    arrayAV[[ i ]]$equivolumetricSphericAreaRatio<-NA
    arrayAV[[ i ]]$Area<-NA
  }
  
  return(arrayAV[[i]])  
}
RAD.areaVolume.mmButo<-function( listaROIVoxels ) {
  objS<-services();
  obj.mmButo<-mmButo();
  arrayAV<-list()
  
  # progression bar  
  iterazione<-0
  pb = txtProgressBar(min = 0, max = length(listaROIVoxels), initial = 0)
  
  for ( i in names(listaROIVoxels) ) {
    if((TRUE %in% is.na(listaROIVoxels[[i]])) == FALSE  ) {
      geometry<-listaROIVoxels[[ i ]]$geometricalInformationOfImages;
      pSX<-geometry$pixelSpacing[1]
      pSY<-geometry$pixelSpacing[2]
      pSZ<-as.numeric(geometry$SliceThickness  )
      # expand the cropped voxelCube
      voxelCube <- obj.mmButo$mmButoLittleCube.expand(   listaROIVoxels[[i]] )
      arrayAV[[ i ]]$Area<-objS$rawSurface(voxelMatrix = voxelCube, pSX = pSX, pSY=pSY,pSZ=pSZ)    
      if ( arrayAV[[ i ]]$Area == -1 ) {
        arrayAV[[ i ]]$Volume<- -1
        arrayAV[[ i ]]$equivolumetricSphericAreaRatio<- -1
      }
      else {
        arrayAV[[ i ]]$Volume<-length(which(voxelCube!=0))*pSX*pSY*pSZ
        arrayAV[[ i ]]$equivolumetricSphericAreaRatio<- ( 4*pi* (   (3/(4*pi))*arrayAV[[ i ]]$Volume   )^(2/3) ) / arrayAV[[ i ]]$Area
      }
    } else {
      arrayAV[[ i ]]$Volume<-NA
      arrayAV[[ i ]]$equivolumetricSphericAreaRatio<-NA
      arrayAV[[ i ]]$Area<-NA
    }
    setTxtProgressBar(pb,iterazione)
    iterazione<-iterazione+1
    
  }
  close(pb)
  
  return(arrayAV)
}
#' Return the list of the available features
#' 
#' @description  returns the list of the available features
#' @return an array of strings
#' @export
RAD.getFeaturesListNames<-function(){
  return(
    c("Mean","Variance","Standard Deviation","Entropy","Skewness","Kurtosys")
    )
}
#' Extract a set of features 
#' 
#' @description  extract a set of features from a given voxel cube and geoLet obj
#' @return many stuff
#' @export
RAD.ExtractFeatures<-function( whichFeatures=c(), obj.geoLet = obj , voxelCube = vc1  ) {
  if( length(whichFeatures)==0){
    whichFeatures<-RAD.getFeaturesListNames()
  }
  arr.features <- length(whichFeatures)
  
  res<-list()
  for(i in seq(1,arr.features)) {
    res[[ whichFeatures[[ i ]] ]] <- rnorm(n = 1)
  }
  return( list("result"=res)   )
}
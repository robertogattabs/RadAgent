#' class for loading and presenting DICOM data
#' 
#' @description  Instantiate an object of the class \code{geoLet}.This represents just the classname, 
#'               methods are exposed with the technique of 'closure'.
#'               In order to see manuals for the single mathods, consider the vignette or use the 
#'               available for the following wrapping functions:
#'               \itemize{
#'               \item \code{GLT.openDICOMFolder( );} : to load a DICOM series into an geoLet object
#'               }
#'               The original methods for the class geoLet can also be invocked using the same name without the previx 'GTL.', i.e.:
#' @examples \dontrun{
#' # first of all create an object aa
#' aa <- RAD.scoutV2(path = './RETTO/pazienti182/Pre/')
#' 
#' 
#' aa$scoutSigma(CSpatient = datasetALL$Codice.sanitario, outcome = datasetALL$TRG, ROIName = 'GTV', feature.family = c("stat"), interpolate = F, from.sigma = 0.3, to.sigma = 0.5, def.by = 0.1, forceRecalculus = F)
#' aa$pixelSpacing()
#' 
#' 
#' 
#' 
#' }
#' @export
#' @useDynLib moddicomV2 
#' @import rmarkdown
RAD.scoutV2 <- function (path){
  attr.G.def.RMark.dir<-''
  attr.G.def.output.report<-''
  attr.G.def.DICOM.pattern.ext<-''
  attr.G.def.report.fileName.scout<-''
  attr.G.def.report.fileName.sigma<-''
  attr.G.path.final<-''
  attr.G.fileName<-''
  
  
  
  # --------------------------------------------------
  # setPath
  # -------------------------------------------------- 
  setPath<-function(def.RMark.dir = NA, def.output.report = NA,
                    def.DICOM.pattern.ext = NA, def.report.fileName.scout=NA, def.report.fileName.sigma=NA,
                    fileName=NA) {
    
    if(!is.na(def.RMark.dir)) attr.G.def.RMark.dir<<-def.RMark.dir
    if(!is.na(def.output.report)) attr.G.def.output.report<<-def.output.report
    if(!is.na(def.DICOM.pattern.ext)) attr.G.def.DICOM.pattern.ext<<-def.DICOM.pattern.ext
    if(!is.na(def.report.fileName.scout)) attr.G.def.report.fileName.scout<<-def.report.fileName.scout
    if(!is.na(def.report.fileName.sigma)) attr.G.def.report.fileName.sigma<<-def.report.fileName.sigma
    if(!is.na(fileName)) attr.G.fileName<<-fileName
    
  }
  # --------------------------------------------------
  # getPath
  # --------------------------------------------------
  getPath<-function() {
    return(list("def.RMark.dir"=attr.G.def.RMark.dir,"def.output.report"=attr.G.def.output.report,"def.DICOM.pattern.ext"=attr.G.def.DICOM.pattern.ext,
                "def.report.fileName.scout"=attr.G.def.report.fileName.scout, 
                "def.report.fileName.sigma"=attr.G.def.report.fileName.sigma, "fileName"=attr.G.fileName))
  }
  # --------------------------------------------------
  # pixelSpacing
  # --------------------------------------------------  
  pixelSpacing <- function(){
    
    patList <- list.dirs(path = attr.G.path.final,recursive = FALSE)
    pixelX <- c()
    pixelY <- c()
    
    for( patID in patList) {
      objServ<-services()
      pixelSpacing <- c()
      ogg.geoLet <- geoLet()
      dicomList <- list.files(path = patID, pattern = attr.G.def.DICOM.pattern.ext)
      devo.uscire <- FALSE
      ct<-1
      while( devo.uscire == FALSE ) {
       
        TAG <- objServ$getDICOMTag(tag = "0008,0016", fileName = paste(patID,dicomList[ct], sep = "/"))
        
        if(TAG=="1.2.840.10008.5.1.4.1.1.2" | TAG=="1.2.840.10008.5.1.4.1.1.4" | TAG=="1.2.840.10008.5.1.4.1.1.128"){
          pixelSpacing <- objServ$getDICOMTag(tag = "0028,0030", fileName = paste(patID,dicomList[ct], sep = "/"))
          pixelSpacing <- as.numeric(strsplit(pixelSpacing,split = "\\\\")[[1]])
          devo.uscire <- TRUE
        }
        ct <- ct +1
      }
      
      pixelX <- c(pixelX, pixelSpacing[1])
      pixelY <- c(pixelY, pixelSpacing[2])
      rm(ogg.geoLet)
      gc()
    }
    
    dir.x.RMARKDOWN <- paste( c(attr.G.def.RMark.dir,"scoutPixelSpacing.Rmd")  , collapse='')
    
    
    render(dir.x.RMARKDOWN, "pdf_document", output_file=attr.G.def.report.fileName.scout, output_dir=attr.G.def.output.report)
    
    
  }
  
  
  # --------------------------------------------------
  # pixelSpacing
  # --------------------------------------------------
  
  scoutSigma <- function(CSpatient, outcome, ROIName, feature.family, interpolate=F, px="", py="", error,
                         from.sigma, to.sigma, def.by, forceRecalculus){
    
    Dataset <- data.frame("CS"= as.character(CSpatient), "outcome"=as.logical(outcome))
    Dataset$CS <- as.character(Dataset$CS)
    
    CalcoloSigma <- f.extractor.pluri.LoG.par(path = attr.G.path.final, ROIName = ROIName, 
                                              feature.family = feature.family,
                                              interpolate = interpolate, px=px, py=py, error = error,
                                              fileName = attr.G.fileName ,
                                              forceRecalculus = forceRecalculus, from.sigma=from.sigma, 
                                              to.sigma=to.sigma, def.by = def.by)
    
    CalcoloSigmaNEW <- list()
    CalcoloSigmaID <- list()
    for (i in names(CalcoloSigma)){
      
      ID <- c()
      for (indice in 1:nrow(CalcoloSigma[[1]])){
        
        CS <- strsplit(CalcoloSigma[[i]][indice],split = "//")[[1]][2]
        CS <- as.character(as.numeric(CS))
        
        ID <- c(ID, CS)
        
      }
      
      CalcoloSigmaID[[i]] <- cbind(ID,  CalcoloSigma[[i]])
      colnames(CalcoloSigmaID[[i]])[1] <- "CS"
      CalcoloSigmaNEW[[i]] <- merge(Dataset, CalcoloSigmaID[[i]], by = "CS")
      
    }
    
    NomeColonne <- (colnames(CalcoloSigmaNEW[[1]]))[-c(1:3)]
    SignificantPValue <- list()
    NotSignificantPValue <- list()
    
    for (features in NomeColonne){
      pValue <- c()
      for (i in names(CalcoloSigmaNEW)){
        print(i)
        pValue <- c(pValue, round((wilcox.test( as.numeric(as.character(CalcoloSigmaNEW[[i]][which(CalcoloSigmaNEW[[i]]$outcome==T),features])), 
                                                as.numeric(as.character(CalcoloSigmaNEW[[i]][which(CalcoloSigmaNEW[[i]]$outcome==F),features]))))$p.value,
                                  digits = 4))
        
        
      }
      if(length(which(pValue<=0.05) )!=0){
        SignificantPValue[[features]] <- pValue
      }
      else{
        NotSignificantPValue[[features]] <- pValue
      }
    }
    
    dir.x.RMARKDOWN2 <- paste( c(attr.G.def.RMark.dir,"scoutSigma.Rmd")  , collapse='')
    
    
    render(dir.x.RMARKDOWN2, "pdf_document", output_file=attr.G.def.report.fileName.sigma, output_dir=attr.G.def.output.report)
    
  }
  
  
  constructor<-function(path) {
    attr.G.def.RMark.dir<<-'./'
    attr.G.def.output.report<<-'./'
    attr.G.def.DICOM.pattern.ext <<-'\\.dcm$'
    attr.G.def.report.fileName.scout<<-'scoutPixelSpacing.pdf'
    attr.G.def.report.fileName.sigma<<-'scoutSigma.pdf'
    attr.G.path.final<<-path
    attr.G.fileName<<-"tmp.f.extractor.pluri.par.RData"
  }
  constructor(path = path)
  
  return(list("pixelSpacing"=pixelSpacing,"getPath"=getPath,"setPath"=setPath,"scoutSigma"=scoutSigma))
  
}
  
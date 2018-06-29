#' Function that creates a \code{dvhmatrix} class object
#' 
#' @param dvh.number The number of DVHs to be generated
#' @param type The type of DVH to be generated: \code{random} generates all possible variation in dose distribution, 
#'        \code{convex} for DVH similar to PTV structures with good dose homogeneity to high levels,
#'        \code{concave} for DVH with low level of dose delivered to the volume as a spared structure,
#'        \code{mix} for DVH where dose distribution is averaged in the whole volume spectrum.
#' @param dvh.type The type of volume distribution to be created: \code{differential} or \code{cumulative}.
#' @param vol.distr Defines if the volume bins have to be divided by the total volume of the structure (\code{relative}) or not (\code{absolute}).
#' @param max.dose The upper bound of the dose distribution to be simulated.
#' @param dose.bin The dose bin in Gy to compute the value of volume parts in the final DVH.
#' @param volbin.side The value of the side of each cube (in mm) that builds the final volume of the simulated structure.
#' @param min.vol The minimum volume of the range to be simulated.
#' @param max.vol The maximum volume of the range to be simulated.
#' @description  Creates an object of class \code{dvhmatrix}.
#' @details Dose Volume Histograms (DVH) are the basic representation of how the radiation doses distribute inside structures.
#'          Usually they are shown in two forms: \code{differential} or \code{cumulative} being related each other because
#'          \code{cumultive} DVHs are the integral form of \code{"differential"} ones. This function provides a siumlated series
#'          of structures with the features defined by function parameters. The simulated series are useful for producing simulated
#'          Dose/Response models using the modeling functions defined in this package.
#'          class objects.
#' @references Van den Heuvela F. \emph{Decomposition analysis of differential dose volume histograms.} Med Phys. 2006 Feb;33(2):297-307. PubMed PMID: 16532934.
#' @return An object of \code{dvhmatrix} class.
#' @examples ## creates a dvhmatrix object with 200 differential - relative histograms
#' a<-DVH.generate(dvh.number=200, dvh.type="differential", vol.distr="relative")
#' 
#' ## creates a dvhmatrix object containing 150 cumulative - absolute histograms
#' ## maximum dose in simulated series is 60 Gy
#' b<-DVH.generate(dvh.number=150, dvh.type="cumulative", vol.distr="absolute", max.dose=60)
#' @export
DVH.generate<-function(dvh.number, type=c("random","convex","concave","mix"), 
                       dvh.type=c("differential", "cumulative"), vol.distr=c("relative", "absolute"),
                       max.dose = 75, dose.bin = 0.5, volbin.side = 2.5, min.vol=180, max.vol=220) {
  if (min.vol<2) {
    warning("Minimum volumes under 2 cc aren't allowed, min.vol set to 2.")
    min.vol <- 2
  }
  if (max.dose <= 10) stop("max.dose must by higher than 10 Gy")
  # create the vector of volumes  
  volumes <- runif(n = dvh.number, min = min.vol, max = max.vol)  
  volbin.num <- round(volumes/((volbin.side/10)^3)) # number of bins for each volume
  # function for generating convex DVHs voxels series
  convex.dvh <- function(n) {
    mean.dose <- runif(n = 1, min = max.dose - (max.dose/3), max = max.dose - max.dose/6)
    sd.dose <- (max.dose - mean.dose)/2.5
    return(rnorm(n = n, mean = mean.dose, sd = sd.dose))
  }
  # function for generating concave DVHs voxels series
  concave.dvh <- function(n) {
    mean.dose <- runif(n = 1, min = 2, max = max(10, max.dose/2.5))
    sd.dose <- mean.dose/3
    result<-rnorm(n = n, mean = mean.dose, sd = sd.dose)
    # takes only the voxels with dose>=0
    result<-result[which(result>=0)]
    # adds more voxels to reach the expected number with random uniform distribution
    return(c(result, runif(n = n - length(result), min = 0, max = 2*mean.dose)))
  }
  # function for creating mix DVHs voxels series
  mix.dvh <- function(n) {
    contrib<-c(runif(n = 1, min = .05, max = .3), runif(n = 1, min = .05, max = .3))
    # proportions in contributions to final DVH
    contrib<-c(contrib[1], 1 - sum(contrib), contrib[2])
    part1<-concave.dvh(n = round(contrib[1] * n))
    part3<-convex.dvh(n = round(contrib[3] * n))
    part2.mean<-runif(n = 1, min = max.dose/20, max = max.dose - max.dose/5)
    part2.sd  <-(max.dose - part2.mean)/(8 * 1/runif(n = 1, min = 1.5, max = 3.5))
    part2<-rnorm(n = n - (length(part1) + length(part3)), 
                 mean = part2.mean, 
                 sd = part2.sd)
    result<-c(part1, part2, part3)
    result<-result[which(result>=0)]    
    # compensate the negative values when available
    return(c(result, runif(n = n - length(result), min = 0, max = max.dose)))    
  }
  # function that creates random DVHs according the previous three given functions
  random.dvh <- function(n) {
    FUN <- sample(x = c(convex.dvh, concave.dvh, mix.dvh), size = 1, replace = T, prob = c(.25,.25,.5))
    return(FUN[[1]](n))
  }
  
  # voxels creation
  type <- match.arg(type)
  if (type=="convex") dose.voxels<-sapply(X = volbin.num, FUN = convex.dvh)
  if (type=="concave") dose.voxels<-sapply(X = volbin.num, FUN = concave.dvh)
  if (type=="mix") dose.voxels<-sapply(X = volbin.num, FUN = mix.dvh)
  if (type=="random") dose.voxels<-sapply(X = volbin.num, FUN = random.dvh)
  VolBin<-volbin.side^3/1000 # Volume Bin in cc
  # delete dose.voxels with exceeding value over 15% max
  dose.voxels<-sapply(X = dose.voxels, FUN = function(x) return(x[which(x<(max.dose+max.dose*.15))]), USE.NAMES = FALSE)
  # creates the vector of structures volumes
  if (dvh.number > 1) volume<-sapply(X = dose.voxels, FUN = function(x) VolBin * length(x), simplify = TRUE, USE.NAMES = FALSE) else volume <- VolBin * length(dose.voxels)
  result<-new("dvhmatrix")
  # creates the dvhmatrix object
  dvh.type<-match.arg(arg = dvh.type)
  vol.distr<-match.arg(arg = vol.distr)
  result@dvh.type<-"differential" # default value, corrected by DVH.diff.to.cum if dvh.type = "cumulative"
  result@vol.distr<-"absolute"    # default DVH vol.distr type
  result@volume<-volume
  # creates the list of differential histograms
  
  if (dvh.number > 1) {
    hlist<-lapply(X = dose.voxels, FUN = hist, 
                  breaks = seq(from = 0, to = max(c(unlist(lapply(X = dose.voxels, FUN = max))) + dose.bin * 4, max.dose), by = dose.bin), plot = FALSE)
    if (dvh.type=="differential") {
      result@dvh<-cbind(hlist[[1]]$mids, sapply(X = hlist, FUN = function(x) cbind(x$counts * VolBin), simplify = TRUE, USE.NAMES = FALSE))
      #if (vol.distr=="relative") for (n in 2:ncol(result@dvh)) result@dvh[,n]<-result@dvh[,n]/result@volume[n-1]
    }
    if (dvh.type=="cumulative") {
      result@dvh<-cbind(hlist[[1]]$mids, sapply(X = hlist, FUN = function(x) cbind(x$counts * VolBin), simplify = TRUE, USE.NAMES = FALSE))
      result<-DVH.diff.to.cum(dvh = result)
      #if (vol.distr=="relative") result<-DVH.diff.to.cum(dvh = result, relative = TRUE) else result<-DVH.diff.to.cum(dvh = result, relative = FALSE)
    }    
  } else {    
    hlist<-hist(x = dose.voxels, breaks = seq(from = 0, to = max(max(dose.voxels) + dose.bin * 4, max.dose), by = dose.bin), plot = FALSE)
    if (dvh.type=="differential") {
      result@dvh <- cbind(hlist$mids, hlist$counts * VolBin)
      #if (vol.distr=="relative") result@dvh[,2] <- result@dvh[,2]/result@volume
    }
    if (dvh.type=="cumulative") {
      result@dvh <- cbind(hlist$mids, hlist$counts * VolBin)
      result<-DVH.diff.to.cum(dvh = result)
      #if (vol.distr=="relative") result<-DVH.diff.to.cum(dvh = result, relative = TRUE) else result<-DVH.diff.to.cum(dvh = result, relative = FALSE)
    }
  }
  if (vol.distr == "relative") result<-DVH.relative(dvh = result)
  return(result)
}

#' Function that converts differential DVHs into cumulative ones
#' 
#' @param dvh Either an object of class \code{dvhmatrix} or \code{matrix} type.
#' @description Function that converts an object of class \code{dvhmatrix} or a simple \code{matrix} that 
#'              represents a differential DVH into a cumulative one.
#' @return Either an object of class \code{dvhmatrix} or a \code{matrix} according the argument \code{dvh}.
#' @export
DVH.diff.to.cum <- function(dvh) {
  if ((!is.matrix(dvh))&&(class(dvh)!="dvhmatrix")) stop("dvh MUST be either an object of class dvhmatrix or a matrix")
  if (class(dvh)=="dvhmatrix") dvh.matrix<-dvh@dvh else dvh.matrix<-dvh  
  if (class(dvh)=="dvhmatrix") if (dvh@dvh.type=="cumulative") return(dvh)
  
  dvh.size <- dim(dvh.matrix)
  DVHList <- matrix(nrow=dvh.size[1] + 1, ncol=dvh.size[2])   # create the matrix of cumulative DVHs     
  for (m in 2:dvh.size[2]) {                                  
    total.volume <- sum(dvh.matrix[,m])    
    for (n in 1:dvh.size[1])                                     # loop for rows
      DVHList[n+1, m] <- total.volume - sum(dvh.matrix[c(1:n),m]) # elements of the matrix as relative volume
    DVHList[1,m] <- total.volume  # first element is total volume by default    
  }
  DVHList[1,1]<-0
  
  for (n in 2:nrow(DVHList)) DVHList[n,1]<-2*dvh.matrix[n-1,1]-DVHList[n-1,1] # set doses
  if (class(dvh)=="dvhmatrix") {
    dvh@dvh<-DVHList
    dvh@dvh.type<-"cumulative"    
    return(dvh)
  } else return(DVHList)
}


#' Function that converts cumulative DVHs into differential ones
#' 
#' @param dvh Either an object of class \code{dvhmatrix} or \code{matrix} type.
#' @description Function that converts an object of class \code{dvhmatrix} or a simple \code{matrix} that 
#'              represents a cumulative DVH into a differential one.
#' @return Either an object of class \code{dvhmatrix} or a \code{matrix} according the argument \code{dvh}.
#' @export
DVH.cum.to.diff <- function(dvh) {
  if ((!is.matrix(dvh))&&(class(dvh)!="dvhmatrix")) stop("dvh MUST be either an object of class dvhmatrix or a matrix")
  if (class(dvh)=="dvhmatrix") dvh.matrix<-dvh@dvh else dvh.matrix<-dvh
  if (class(dvh)=="dvhmatrix") if (dvh@dvh.type=="differential") return(dvh)
  
  dvh.size <- dim(dvh.matrix)
  DVHList <- matrix(nrow=dvh.size[1] - 1, ncol=dvh.size[2])   # create the matrix of differential DVHs
  for (m in 2:dvh.size[2]) {                                  # loop for columns (volumes)
    total.volume<-dvh.matrix[1,m]                             # set total volume to 1st element of the cumulative DVH
    for (n in 2:dvh.size[1])                                  # loop for rows
      DVHList[n-1,m] <- (dvh.matrix[n-1,m]-dvh.matrix[n,m])   # differential volume for absolute distribution
  }       
  for (m in 2:dvh.size[1]) {                                   # loop for column of dose values
    DVHList[m-1,1]<-dvh.matrix[m-1,1]+(dvh.matrix[m,1]-dvh.matrix[m-1,1])/2 # sets the dose values
  }
  if (class(dvh)=="dvhmatrix") {
    dvh@dvh<-DVHList
    dvh@dvh.type<-"differential"
    return(dvh)
  } else return(DVHList)
}

#' Extracts a DVH from a vector of dose bins
#' @param x vector of dose bins.
#' @param max.dose Upper dose bound for limiting the DVH computation.
#' @param dose.bin The dose bin for DVH computation in Gy.
#' @param dvh.type The type of volume distribution to be created: \code{differential} or \code{cumulative}.
#' @param vol.distr Defines if the volume bins have to be divided by the total volume of the structure (\code{relative}) or not (\code{absolute}).
#' @param createObj if \code{TRUE} returns a \code{dvhmatrix} class object.
#' @param volbin.side The value of the side of each cube (in cm) that builds the final volume of the structure.
#' @param voxel.volume The value of the volume of a single voxel (in cc) to compute the ROI total volume
#' @param total.volume The value of the total volume of the ROI (in cc)
#' @description Function that given a vector of dose bins (either a vector got from sampling a 3D mesh or a vector of
#' simple values of dose) extracts the DVH that summarizes that vector. Either \code{volbin.side}, \code{voxel.volume} or \code{total.volume} have to be
#' not \code{NULL} values. The priority list for computing the volume (if more than one argument is not \code{NULL}) is total.volume \code{->} voxel.volume \code{->} volbin.side.
#' @return Either an object of class \code{dvhmatrix} or a \code{matrix} according the argument \code{dvh} type.
#' @examples ## simulate a vector of dose bins
#' doses <- c(rnorm(n = 10000, mean = 45, sd = 3), rnorm(n = 7000, mean = 65, sd = 2.5))
#' 
#' ## creates a dvhmatrix class object
#' DVH<-DVH.extract(x = doses)
#' @export
DVH.extract<-function(x, max.dose=NULL, dose.bin=.25, dvh.type=c("differential","cumulative"),vol.distr=c("relative","absolute"), createObj=TRUE,  volbin.side=NULL, voxel.volume = NULL, total.volume = NULL) {
  # default VolBin is given in cm3 
  if (!is.null(total.volume)) VolBin <- total.volume / length(x)
  else if (!is.null(voxel.volume)) VolBin <- voxel.volume
  else if (is.null(volbin.side)) stop('Either total.volume, voxel.volume or volbin.side MUST be a not NULL value')
  else VolBin <- volbin.side^3 
  TotalVol<-length(x)*VolBin
  dvh.type=match.arg(dvh.type)
  vol.distr=match.arg(vol.distr)
  if (is.null(x=max.dose)) max.dose<-max(x)+dose.bin*4
  h<-hist(x=x, breaks=seq(from=0, to=max.dose,by=dose.bin), plot=FALSE)
  diff<-cbind(h$mids, h$density/sum(h$density)*TotalVol)
  # function for converting differential DVHs into
  # relative differential DVHS
  rel.diff.dvh <- function(dvh.matrix) {
    dvh.size <- dim(dvh.matrix)
    DVHList <- matrix(nrow=dvh.size[1], ncol=dvh.size[2])       # create the matrix of  DVHs
    for (m in 2:dvh.size[2]) {                                  # loop for columns (volumes) 
      tot.volume <- sum(dvh.matrix[,m])                       # calculate the total volume
      for (n in 1:dvh.size[1]) {                                # loop for rows
        DVHList[n, m] <- dvh.matrix[n, m]/tot.volume          # elements of the matrix as relative volume
      }    
    }
    DVHList[,1]<-dvh.matrix[,1]                                 # sets the dose column  
    return(DVHList)
  }
  # function for converting a matrix of differential DVH 
  # into matrix of cumulative DVH (relative volume)
  cum.dvh <- function(dvh.matrix, relative=TRUE)
  {
    dvh.size <- dim(dvh.matrix)
    DVHList <- matrix(nrow=dvh.size[1] + 1, ncol=dvh.size[2])   # create the matrix of cumulative DVHs     
    for (m in 2:dvh.size[2]) {                                  # loop for columns (volumes) 
      tot.volume <- sum(dvh.matrix[,m])                       # calculate the total volume
      if (relative==TRUE) {
        for (n in 1:dvh.size[1]) {                                    # loop for rows
          DVHList[n+1, m] <- (tot.volume - sum(dvh.matrix[c(1:n),m]))/tot.volume # elements of the matrix as relative volume
        }
        DVHList[1,m] <- 1 # first element is 1 by default if relative==TRUE
      } else {
        for (n in 1:dvh.size[1]) {                                    # loop for rows
          DVHList[n+1, m] <- tot.volume - sum(dvh.matrix[c(1:n),m]) # elements of the matrix as relative volume
        }
        DVHList[1,m] <- tot.volume  # first element is total volume by default
      }                                       
    }
    DVHList[1,1]<-0
    #  DVHList[2:nrow(DVHList),1]<-dvh.matrix[,1]
    for (n in 2:nrow(DVHList)) DVHList[n,1]<-2*dvh.matrix[n-1,1]-DVHList[n-1,1]
    return(DVHList)
  }
  # return matrix without dvhmatrix class structure
  if ((dvh.type=="differential") && (vol.distr=="absolute"))  final.matrix<-diff
  if ((dvh.type=="differential") && (vol.distr=="relative"))  final.matrix<-rel.diff.dvh(diff)
  if ((dvh.type=="cumulative")   && (vol.distr=="absolute"))  final.matrix<-cum.dvh(dvh.matrix=diff, relative=FALSE)
  if ((dvh.type=="cumulative")   && (vol.distr=="relative"))  final.matrix<-cum.dvh(dvh.matrix=diff, relative=TRUE)
  if (createObj==FALSE) return(final.matrix)
  # return matrix within a dvhmatrix object
  if (createObj==TRUE) return(new("dvhmatrix", dvh=final.matrix, dvh.type=dvh.type, vol.distr=vol.distr, volume=TotalVol))  
}

#' mean of DVHs
#' @param dvh A \code{dvhmatrix} object
#' @description Function that gives the value of the mean dose of a \code{dvhmatrix} class object
#' @return a vector with the means of doses in the \code{dvhmatrix} object
#' @export
#' @examples ## generate a dataset of DVHs
#' a<-DVH.generate(dvh.number = 100)
#' m<-DVH.mean(a)
DVH.mean<-function(dvh)  {
  if (class(dvh)!="dvhmatrix") stop("dvh MUST be a dvhmatrix class object")
  if (dvh@dvh.type=="cumulative") dvh<-DVH.cum.to.diff(dvh = dvh)
  return(DVH.eud(dvh = dvh, a = 1)) # use a = 1 that corresponds to mean dose
}

#' Converts absolute \code{dvhmatrix} class objects to relative
#' @param dvh A \code{dvhmatrix} class object
#' @description Function that converts an object of \code{dvhmatrix} class where DVH are stored in absoulte mode into
#' another \code{dvhmatrix} class object where DVH are stored in relative mode.
#' @return A \code{dvhmatrix} class object.
#' @export
#' @examples ## generate a dataset of absolute DVHs
#' a<-DVH.generate(dvh.number = 10, vol.distr = "absolute")
#' b<-DVH.relative(dvh = a)
DVH.relative<-function(dvh) {
  if (dvh@vol.distr=="absolute")
    for (n in 1:(ncol(dvh@dvh) - 1)) dvh@dvh[,n+1]<-dvh@dvh[,n+1]/dvh@volume[n]
    dvh@vol.distr<-"relative"
    return(dvh)
}

#' Converts relative \code{dvhmatrix} class objects to absolute
#' @param dvh A \code{dvhmatrix} class object
#' @description Function that converts an object of \code{dvhmatrix} class where DVH are stored in relative mode into
#' another \code{dvhmatrix} class object where DVH are stored in absolute mode.
#' @return A \code{dvhmatrix} class object.
#' @export
#' @examples ## generate a dataset of relative DVHs
#' a<-DVH.generate(dvh.number = 10, vol.distr = "relative")
#' b<-DVH.absolute(dvh = a)
DVH.absolute<-function(dvh) {
  if (dvh@vol.distr=="relative")
    for (n in 1:(ncol(dvh@dvh) - 1)) dvh@dvh[,n+1]<-dvh@dvh[,n+1]*dvh@volume[n]
    dvh@vol.distr<-"absolute"
    return(dvh)
}

#' Calculates Equivalent Uniform Dose for a \code{dvhmatrix} object
#' @param dvh A \code{dvhmatrix} class object
#' @param a Value for parallel-serial correlation in radiobiological response
#' @description Function that calculates the value of Equivalent Uniform Dose (EUD) for a \code{dvhmatrix} object. The
#' formula to compute the equivalent uniform dose
#' is \deqn{EUD=\left (\sum_{i=1}^{ }D_{i}^{a}\cdot \frac{V_{i}}{V_{tot}}  \right )^{1/a}}{EUD=(Sum(Di^a * Vi/Vtot))^(1/a)}
#' where \emph{i} is the counter of volume bins in the structure and \emph{a} is the parameter that relates
#' the organs response to the \emph{serial} or \emph{parallel} physiological function of the organ. Usually \emph{parallel}
#' structures are organs that can suffer depletion of functional subunits up to a given threshold, as liver, lungs or kidneys.
#' They are more sentitive to mean dose rather than hot spots and maximum dose levels delivered to the structure. On the other hand
#' \emph{serial} structures are organs more sensitive to hot-spots delivered in the volume, such as spinal cord.
#' @return A vector containing the values of EUD(s) for the given DVH(s)
#' @export
#' @references Niemierko A. \emph{Reporting and analyzing dose distributions: a concept of equivalent uniform dose.} Med Phys. 1997 Jan;24(1):103-10. PubMed PMID: 9029544.
#' @useDynLib moddicomV2
DVH.eud<-function(dvh, a = 1) {
  dvh<-DVH.cum.to.diff(dvh = dvh)
  dvh<-DVH.relative(dvh = dvh)
  Ncol<-ncol(dvh@dvh) - 1
  Nrow<-nrow(dvh@dvh)
  ceud<-rep.int(x = 0, times = Ncol) 
  dosebin<-dvh@dvh[,1]
  volumebin<-dvh@dvh[,2:(Ncol + 1)]
  result<-.C("cEUD", as.double(dosebin), as.double(volumebin), as.double(a), as.integer(Nrow), 
             as.integer(Ncol), as.double(ceud))
  return(result[[6]])
}


#' Merge two different \code{dvhmatrix} class objects into one
#' @param receiver The \code{dvhmatrix} object that will receive the \code{addendum} object.
#' @param addendum The \code{dvhmatrix} object to be merged with \code{addendum}.
#' @description This function can e used to create a \code{dvhmatrix} class object from two different
#'              ones. The slots of the resulting object will be the same in the \code{receiver} one,
#'              so some convertions can be automatically realized to create the homogeneous final result.
#' @return A \code{dvhmatrix} class object.
#' @export
#' @examples ## creates two different dvhmatrx objects
#' a<-DVH.generate(dvh.number = 100, dvh.type="differential", vol.distr = "relative")
#' b<-DVH.generate(dvh.number = 100, dvh.type="cumulative", vol.distr = "absolute")
#' ab<-DVH.merge(receiver = a, addendum = b)
DVH.merge<-function(receiver=NULL, addendum=NULL) {
  # check dvhmatrix class
  if ((class(receiver)!="dvhmatrix") || (class(addendum)!="dvhmatrix")) 
    stop("BOTH receiver AND addendum MUST be dvhmatrix class objects")
  # creates list of dose bins
  dbin<-list(receiver=receiver@dvh[3,1]-receiver@dvh[2,1], addendum=addendum@dvh[3,1]-addendum@dvh[2,1])
  
  if (receiver@vol.distr=="relative") { # relative state of receiver
    rel<-TRUE
    vol.distr.fun.name <- "DVH.relative"
  } 
  else {
    rel<-FALSE
    vol.distr.fun.name <- "DVH.absolute"
  }  
  # STEP (0): check the relative/absolute distribution
  if (receiver@vol.distr!=addendum@vol.distr) {
    vol.distr.fun <- match.fun(FUN = vol.distr.fun.name)
    addendum <- vol.distr.fun(addendum)
  }
  # STEP (1): identify conversion function if DVH types are different and convert addendum according receiver type
  if (receiver@dvh.type!=addendum@dvh.type) {
    if (receiver@dvh.type=="cumulative") conv.funct.name<-"DVH.diff.to.cum"
    if (receiver@dvh.type=="differential") conv.funct.name<-"DVH.cum.to.diff"
    conv.funct<-match.fun(FUN=conv.funct.name)  
    # correct the addendum if receiver has absolute distribution
    addendum@dvh<-conv.funct(dvh=addendum@dvh)
  }
  # STEP (2): check the dose-bin for interpolating addendum or receiver according the lowest dbin
  if (dbin$receiver!=dbin$addendum) {
    warning("Different dose bins between receiver and addendum: linear interpolation performed.")
    dose.seq<-seq(from = min(receiver@dvh[,1], addendum@dvh[,1]),
                  to   = max(receiver@dvh[,1], addendum@dvh[,1]),
                  by   = min(diff(x = receiver@dvh[,1]), diff(x = addendum@dvh[,1])))
    if (ncol(receiver@dvh) == 2) { # single DVH receiver
      apf.receiver  <- approxfun(x = receiver@dvh[,1], y = receiver@dvh[,2])
      temp.receiver <- cbind(dose.seq, apf.receiver(dose.seq))
    } else {  # multiple DVH receiver
      apf.receiver  <- apply(X = receiver@dvh[,2:ncol(receiver@dvh)], MARGIN = 2, FUN = approxfun, x = receiver@dvh[,1])  # generate list of approxfun
      temp.receiver <- cbind(dose.seq, sapply(X = apf.receiver, FUN = function(x) return(x(dose.seq)), USE.NAMES = FALSE) )  # generate new dvh matrix
    }
    if (ncol(addendum@dvh) == 2) { # single DVH addendum
      apf.addendum  <- approxfun(x = addendum@dvh[,1], y = addendum@dvh[,2])
      temp.addendum <- cbind(dose.seq, apf.addendum(dose.seq))
    } else {  # multiple DVH addendum
      apf.addendum  <- apply(X = addendum@dvh[,2:ncol(addendum@dvh)], MARGIN = 2, FUN = approxfun, x = addendum@dvh[,1])  # generate list of approxfun
      temp.addendum <- cbind(dose.seq, sapply(X = apf.addendum, FUN = function(x) return(x(dose.seq)), USE.NAMES = FALSE) )  # generate new dvh matrix
    }
    temp.receiver[which(is.na(temp.receiver))]<-0
    temp.addendum[which(is.na(temp.addendum))]<-0
  } else {
    temp.receiver<-receiver@dvh
    temp.addendum<-addendum@dvh
  }  
  
  # STEP (3): check if the two dvhs have equal length
  if (nrow(temp.receiver)>nrow(temp.addendum)) {
    for (n in 1:(nrow(temp.receiver)-nrow(temp.addendum))) temp.addendum<-rbind(temp.addendum, 0)
    temp.addendum[,1]<-temp.receiver[,1]
  } 
  if (nrow(temp.receiver)<nrow(temp.addendum)) {
    for (n in 1:(nrow(temp.addendum)-nrow(temp.receiver))) temp.receiver<-rbind(temp.receiver, 0)
    temp.receiver[,1]<-temp.addendum[,1]
  }  
  
  # STEP (4): joins the two dvhs
  result.DVH<-unname(cbind(temp.receiver, temp.addendum[,2:ncol(temp.addendum)]))
  return(new("dvhmatrix", dvh=result.DVH, dvh.type=receiver@dvh.type, vol.distr=receiver@vol.distr, 
             volume=c(receiver@volume, addendum@volume)))
}

#' method for plotting the Log-Likelihood values in \code{DoseVolumeModel} class object
#' @description This method plots the results of fitting process for a \emph{Dose/Volume} model. The \emph{x} axis represents either the Dose,
#'              or the Volume, the \emph{y} axis represents the Log-Likelihood result of the modeing process.
#' @param x The \code{DoseVolumeModel} class object
#' @param ShowFittingValue A logical value for showing the fitting point in the Log-Likelihood plot
#' @param ... Other parameters to be passed to \code{plot} function
#' @details The \code{DoseVolumeModel} class object contains the results of a modeling process across different values of Dose or 
#'          Volumes (depending on which kind of fitting was chosen during fitting process). The plot doesn't show the impact of other covariates
#'          that could be part of the model itself.
#' @exportMethod plot          
setMethod('plot', signature(x = 'DoseVolumeModel'),
          function(x, ShowFittingValue = TRUE, ...) {
            M  <- slot(object = x, name = 'output.matrix')
            tp <- slot(object = x, name = 'fitted.parameter')
            value <- slot(object = x, name = 'fitted.value')
            if (tp == 'Vdose')   xlbl <- 'V-Dose [Gy]'
            if (tp == 'Dvolume') xlbl <- 'D-Volume [cc]'
            plot(x = M[,1], y = M[,2], xlab = xlbl, ylab = 'Log-Likelihood', type = 'l', ...)
            if (ShowFittingValue == TRUE) {
              abline(h = min(M[,2]), lty = 2)
              abline(v = value, lty = 2)
              points(x = value, y = min(M[,2]), ...)
              axis(side = 1, at = round(x = value, digits = 2))
            }
          }
)

# method for class dvhmatrix
# enhance: is a vector of number of DVHs to be enhanced in the plot
# mean: plots the mean DVH overlapped to the background
# elements: if equal to a string "all" (default value) plots all the DVHs, if it is a vector plots only the
#   chosen DVHs
#' method for plotting \code{dvhmatrix} class objects
#' @description This method overrides the default method \code{plot} and handles the plots of \code{dvhmatrix}
#'              class objects.
#' @param x The \code{dvhmatrix} object to be plotted.
#' @param elements A vector containing the indeces of the DVHs to be plotted.
#' @param enhance A vector containing the indeces of the DVHs curves to be enhanced in the plot by changing de default color.
#' @param mean.dvh A \code{logical} value, if \code{TRUE} the mean dvh is plotted.
#' @param median.dvh A \code{logical} value, if \code{TRUE} the median dvh is plotted.
#' @param mean.median.alone A \code{logical} value, if \code{TRUE} only the mean and/or median DVH(s) is/are plotted.
#'        Mean or median DVHs are plotted according the values in \code{mean.dvh} and \code{median.dvh}.
#' @param el.color The color for plotting DVHs. The default value is \code{"black"}.
#' @param en.color The color for plotting the enhanced DVHs. The default value is \code{"red"}.
#' @param mean.color The color for the mean dvh plot. The default value is \code{"blue"}.
#' @param median.color The color for the median dvh plot. The default value is \code{"red"}.
#' @param lwd The size of the lines in the plot. Default value is 1.
#' @param C.I.dvh Plot the confidence interval of mean and/or median DVH(s) according the values in \code{mean.dvh} and
#'        \code{median.dvh}.
#' @param C.I.dvh.width The width of confidence interval of mean, median DVH(s) and overall dvh series.
#' @param C.I.dvh.range Plot the confidence interval of the whole DVHs series according the value in \code{C.I.dvh.width}.
#' @param n.boot The number of bootstrap repetitions for calculating the confidence interval of mean and median DVHs.
#' @param ... Other parameters to be passed to \code{plot} function.
#' @details The dvh are stored in the slot \code{dvh} of a \code{dvhmatrix} class object. Of course it is possible for7
#'          users to use these matrices for plotting the dvhs considering that the structure of the matrix is referenced
#'          under \code{\link{dvhmatrix-class}}. In order to simplify the plotting operations this method has been coded
#'          allowing users to automatically obtain some useful features in the graphs, as the mean and median dvh curves,
#'          and the confidence intervals for mean, median and overall dvh matrix. The confidence intervals of mean and
#'          median are calculated by bootstrapping the whole \code{dvhmatrix} and calculating a given number of mean and
#'          median dvh values that finally are reduced by using a quantile function to create the final plot.
#' @exportMethod plot
#' @examples ## create a dvhmatrix object
#' b <- DVH.generate(dvh.number = 100, dvh.type = "cumulative", vol.distr = "absolute")
#' plot(x = b, mean.dvh = TRUE, median.dvh = TRUE, mean.median.alone = TRUE, 
#'      C.I.dvh = TRUE, C.I.dvh.range = TRUE, C.I.dvh.width = .67, lwd = 2)
setMethod("plot", signature(x="dvhmatrix"), 
          function(x, elements=NULL, enhance=NULL, mean.dvh=FALSE, median.dvh=FALSE, 
                   mean.median.alone=FALSE, el.color="black", en.color="red",  mean.color="blue", median.color="red",
                   lwd=1, C.I.dvh=FALSE, C.I.dvh.width=.95, C.I.dvh.range=FALSE, n.boot=2000, ...){
            dvh <- slot(x,"dvh")
            dvh.type <- slot(x, "dvh.type")
            vol.distr <- slot(x, "vol.distr")
            # define the lables for the y axis
            if ((dvh.type=="cumulative") && (vol.distr=="relative")) ylab<-"Volume [*100%]"
            if ((dvh.type=="cumulative") && (vol.distr=="absolute")) ylab<-"Volume [cc]"
            if ((dvh.type=="differential") && (vol.distr=="relative")) ylab<-"dVolume/dDose [*100%/Gy]"
            if ((dvh.type=="differential") && (vol.distr=="absolute")) ylab<-"dVolume/dDose [cc/Gy]"
            # max value of y axes
            ymax<-max(dvh[,2:ncol(dvh)])
            # direct plot for a single DVH
            if (ncol(dvh)==2) plot(x=dvh[,1], y=dvh[,2], type="l", ylab=ylab, xlab="Dose [Gy]", lwd=lwd)
            # plot for more than one DVH
            if ((ncol(dvh)>2) && (mean.median.alone==FALSE)) {
              # creates the vector of elements to be plotted and enhanced
              if (is.null(elements)) elements<-c(2:ncol(dvh)) else elements<-elements+1
              if (!is.null(enhance)) for (n in 1:length(enhance)) elements<-elements[elements!=(enhance[n]+1)]              
              # open the plot panel and raws the 1st plot
              if (length(elements)>0) plot(x=dvh[,1], y=dvh[,elements[1]], ylim=c(0, ymax), col=el.color, 
                                           lwd=lwd, type="l", xlab="Dose [Gy]", ylab=ylab)
              else
                if (length(elements)>0) lines(x=dvh[,1], y=dvh[,elements[1]], ylim=c(0, ymax), col=el.color, 
                                              lwd=lwd)              
              # plots the remaining DVHs into "elements"
              if (length(elements)==0) {
                start.enhance<-2
                plot(x=dvh[,1], y=dvh[,enhance[1]+1], col=en.color, lwd=lwd, type="l")
              } else start.enhance<-1
              if (length(elements)>1)
                for (n in 2:length(elements)) lines(x=dvh[,1], y=dvh[,elements[n]], col=el.color, lwd=lwd)
              # plots DVHs to be enhanced
              if (length(enhance)>0)
                for (n in start.enhance:length(enhance)) lines(x=dvh[,1], y=dvh[,enhance[n]+1], col=en.color, lwd=lwd)
            }
            # plot the mean DVH
            if (ncol(dvh)>2) {
              mean.DVH<-apply(X=dvh[,2:ncol(dvh)], MARGIN=1, FUN=mean)
              median.DVH<-apply(X=dvh[,2:ncol(dvh)], MARGIN=1, FUN=median)
            }
            if (mean.median.alone==FALSE) { 
              if ((mean.dvh==TRUE) && (ncol(dvh)>2)) 
                lines(x=dvh[,1], y=mean.DVH, col=mean.color, lwd=lwd)
              if ((median.dvh==TRUE) && (ncol(dvh)>2)) 
                lines(x=dvh[,1], y=median.DVH , col=median.color, lwd=lwd)
            }       
            
            if (mean.median.alone==TRUE) {
              # plot both mean and median
              if ((mean.dvh==TRUE) && (median.dvh==TRUE)) {
                plot(x=dvh[,1], y=mean.DVH, ylim=c(0, ymax), col=mean.color, 
                     lwd=lwd, type="l", xlab="Dose [Gy]", ylab=ylab)
                lines(x=dvh[,1], y=median.DVH , col=median.color, lwd=lwd)
              }
              # plot only mean
              if ((mean.dvh==TRUE) && (median.dvh==FALSE)) {
                plot(x=dvh[,1], y=mean.DVH, ylim=c(0, ymax), col=mean.color, 
                     lwd=lwd, type="l", xlab="Dose [Gy]", ylab=ylab)
              }
              # plot only median
              if ((mean.dvh==FALSE) && (median.dvh==TRUE)) {
                plot(x=dvh[,1], y=median.DVH, ylim=c(0, ymax), col=median.color, 
                     lwd=lwd, type="l", xlab="Dose [Gy]", ylab=ylab)
              }
            }
            # plot of C.I. area, both 95% C.I. of the mean dose and 95% with quantiles are calculated
            if ((C.I.dvh==TRUE) && (ncol(dvh)>2)) {
              message("sampling DVHs for C.I. calculation...")
              DVH.CI<-DVH.baseStat(dvh = x, C.I.width = C.I.dvh.width, n.boot = n.boot)
              if ((mean.dvh==TRUE)&&(median.dvh==FALSE)){            
                lowerCI.mean <- DVH.CI@dvh[,3]
                higherCI.mean<- DVH.CI@dvh[,4]
                xpoly.mean<-c(dvh[,1], rev(dvh[,1]))             # vector of x coordinates of polygon of CI
                ypoly.mean<-c(lowerCI.mean, rev(higherCI.mean))  # vector of y coordinates of polygon of CI of mean
                # draws the polygon of mean CI
                polygon(x=xpoly.mean, y=ypoly.mean, col=adjustcolor(col=mean.color, alpha.f=.25))
              }
              
              if ((mean.dvh==FALSE)&&(median.dvh==TRUE)){
                lowerCI.median <- DVH.CI@dvh[,6]
                higherCI.median<- DVH.CI@dvh[,7]
                xpoly.median<-c(dvh[,1], rev(dvh[,1]))                # vector of x coordinates of polygon of CI
                ypoly.median<-c(lowerCI.median, rev(higherCI.median)) # vector of y coordinates of polygon of CI of median
                # draws the polygon of median CI
                polygon(x=xpoly.median, y=ypoly.median, col=adjustcolor(col=median.color, alpha.f=.25))
              }
              
              if ((median.dvh==TRUE)&&(mean.dvh==TRUE)){
                
                lowerCI.mean <- DVH.CI@dvh[,3]
                higherCI.mean<- DVH.CI@dvh[,4]
                xpoly.mean<-c(dvh[,1], rev(dvh[,1]))             # vector of x coordinates of polygon of CI
                ypoly.mean<-c(lowerCI.mean, rev(higherCI.mean))  # vector of y coordinates of polygon of CI of mean
                lowerCI.median <- DVH.CI@dvh[,6]
                higherCI.median<- DVH.CI@dvh[,7]
                xpoly.median<-c(dvh[,1], rev(dvh[,1]))                # vector of x coordinates of polygon of CI
                ypoly.median<-c(lowerCI.median, rev(higherCI.median)) # vector of y coordinates of polygon of CI of median
                # draws the polygon of mean CI
                polygon(x=xpoly.mean, y=ypoly.mean, col=adjustcolor(col=mean.color, alpha.f=.25))
                # draws the polygon of median CI
                polygon(x=xpoly.median, y=ypoly.median, col=adjustcolor(col=median.color, alpha.f=.25))
              }
            }
            
            # plot the C.I. of data range of DVHs
            if (C.I.dvh.range==TRUE) {
              # quantiles for DVHs data range
              lowerCI <- c()
              higherCI <- c()
              for (n in 1:nrow(dvh)) {
                lowerCI <- c(lowerCI,  quantile(x=dvh[n,2:ncol(dvh)], probs=((1-C.I.dvh.width)/2)))
                higherCI<- c(higherCI, quantile(x=dvh[n,2:ncol(dvh)], probs=((1+C.I.dvh.width)/2)))
              }
              xpoly<-c(dvh[,1], rev(dvh[,1]))      # vector of x coordinates of polygon of CI
              ypoly<-c(lowerCI, rev(higherCI))     # vector of y coordinates of polygon of CI              
              polygon(x=xpoly, y=ypoly, col=adjustcolor(col="grey", alpha.f=.25))  # draws the polygon
            }            
          }
)



#' Function for extracting the V-Dose from cumulative DVH(s)
#' @description Function for calculating the value of V-Dose of a \code{dvhmatrix} object
#' @param dvh A \code{dvhmatrix} class object
#' @param Dose The dose for calculating the V-Dose value(s)
#' @details V-Dose is, in a given Dose Volume Histogram, the value of the volume of the structure receiving
#'          \strong{at least} the chosen level of dose.
#' @return A vector containing the value(s) of V-Dose(s).
#' @export
#' @examples ## create a dvhmatrix class object
#' a<-DVH.generate(dvh.number = 100)
#' DVH.Vdose(dvh = a, Dose = 50)

DVH.Vdose <- function(dvh, Dose) {
  dvh<-DVH.diff.to.cum(dvh = dvh)
  dv<-dvh@dvh[,1]  # vector of doses
  if (ncol(dvh@dvh) > 2) {
    apf <- apply(X = dvh@dvh[,2:ncol(dvh@dvh)], MARGIN = 2, FUN = approxfun, x = dv)  # generate list of approxfun
    return( sapply(X = apf, FUN = function(x) return(x(Dose))) )
  } else { # option with 1 single DVH
    apf <- approxfun(x = dv, y = dvh@dvh[,2])
    return(apf(Dose))
  }
}

#' Function for extracting the D-Volume from cumulative DVH(s)
#' @description Function for calculating the value of D-Volume of a \code{dvhmatrix} object
#' @param dvh A \code{dvhmatrix} class object
#' @param Volume Volume for calculating the D-Volume value(s). Default is 0.001, corresponding to about maximum dose
#'        delivered to the structure.
#' @details D-Volume is, in a given Dose Volume Histogram, the value of the dose received
#'          by a given volume of the structure(s) of interest.
#' @return A vector containing the value(s) of D-Volume(s).
#' @export
#' @examples ## create a dvhmatrix class object
#' a<-DVH.generate(dvh.number = 100)
#' DVH.Dvolume(dvh = a, Volume = 0.5)
DVH.Dvolume <- function(dvh,  Volume=0.001) {
  dvh<-DVH.diff.to.cum(dvh = dvh)
  dv<-dvh@dvh[,1]  # vector of doses
  if (ncol(dvh@dvh) > 2) {
    apf <- apply(X = dvh@dvh[,2:ncol(dvh@dvh)], MARGIN = 2, FUN = approxfun, y = dv)  # generate list of approxfun
    return( sapply(X = apf, FUN = function(x) return(x(Volume))) )
  } else { # option with 1 single DVH
    apf <- approxfun(x = dvh@dvh[,2], y = dv)
    return(apf(Volume))
  }
}

#' Function for calculating the mean and median dvh with related confidence intervals
#' @description This function calculates the mean dvh from a \code{dvhmatrix} class object. The mean dvh is
#'              calculated with its confidence interval that is given by a bootstrapped dvh series from
#'              the dvh given in the \code{dvh} object.
#' @param dvh A \code{dvhmatrix} class object
#' @param C.I.width The width of confidence interval
#' @param n.boot The number of bootstrapped dvhs for computing the quantile in C.I. calculation
#' @return A \code{dvhmatrix} object where 7 columns: Dose, mean DVH, low C.I. of mean DVH, high C.I. of mean DVH, 
#'         median DVH, low C.I. of median DVH, high C.I. of median DVH.
#' @export
#' @useDynLib moddicomV2
#' @examples ## creates a dvhmatrix class object with 100 DVHs
#' a<-DVH.generate(dvh.number = 100, dvh.type = "cumulative", vol.distr = "relative")
#' b<-DVH.baseStat(dvh = a)
#' ## plot the dvhmatrix object "b", showing both mean, 
#' ## median DVH and corresponding confidence intarvals
#' plot(b)
DVH.baseStat<-function(dvh, C.I.width = .95, n.boot = 2000) {
  # calculate mean dvh
  mean.dvh<-apply(X = dvh@dvh[,2:ncol(dvh@dvh)], MARGIN = 1, FUN = mean)
  median.dvh<-apply(X = dvh@dvh[,2:ncol(dvh@dvh)], MARGIN = 1, FUN = median)
  # vectorized DVH
  Vdvh<-as.vector(dvh@dvh[,2:ncol(dvh@dvh)])
  # vector for sampled dvh
  sampledvh<-rep.int(x = 0, times = nrow(x = dvh@dvh) * (ncol(x = dvh@dvh) - 1))
  # create the mean dvh vector
  meanV<-rep.int(x = 0, times = nrow(dvh@dvh) * n.boot)
  medianV<-meanV
  stepMedian<-rep.int(x = 0, times = ncol(dvh@dvh) - 1)
  # number of histograms
  Ndvh <- ncol(dvh@dvh) - 1
  result<-(.C("meanmediandvh", as.double(Vdvh), as.integer(nrow(dvh@dvh)), as.integer(n.boot), 
              as.double(meanV), as.double(sampledvh), as.integer(Ndvh), as.double(medianV), as.double(stepMedian)))
  # create dvh matrix
  meanDvh<-matrix(data = result[[4]], nrow = nrow(dvh@dvh))
  # generate the mean CI DVH
  meanCI<-apply(X = meanDvh, MARGIN = 1, FUN = quantile, probs = c((1-C.I.width)/2, (1+C.I.width)/2))
  # generate the median CI DVH
  medianDvh<-matrix(data = result[[7]], nrow = nrow(dvh@dvh))
  medianCI<-apply(X = medianDvh, MARGIN = 1, FUN = quantile, probs = c((1-C.I.width)/2, (1+C.I.width)/2))
  #  return(new("dvhmatrix", dvh = cbind(dvh@dvh[,1], medianDvh), dvh.type = dvh@dvh.type,
  #             vol.distr = dvh@vol.distr, volume = rep.int(x = mean(dvh@volume), times = n.boot)))
  result.DVH<-new("dvhmatrix", dvh = cbind(dvh@dvh[,1], mean.dvh, meanCI[1,], meanCI[2,], median.dvh, medianCI[1,], medianCI[2,]), dvh.type = dvh@dvh.type,
                  vol.distr = dvh@vol.distr, volume = c(rep.int(x = mean(dvh@volume), times = 3),
                                                        rep.int(x = median(dvh@volume), times = 3)))
  colnames(x = result.DVH@dvh)<-c("Dose [Gy]", "mean", paste((1-C.I.width)/2,"% C.I. of mean"), paste((1+C.I.width)/2,"% C.I. of mean"),
                                  "median", paste((1-C.I.width)/2,"% C.I. of median"), paste((1+C.I.width)/2,"% C.I. of median"))
  return(result.DVH)
}

#' Function that returns the linear quadratic correction of dose bin in a \code{dvhmatrix} object
#' @description This function appies the linear-quadratic correction for the dose in a \code{dvmatrix} object by 
#' converting each dose bins as function of a given \eqn{\alpha\beta} ratio (default value is 3, valid for late outcome).
#' @param dvh A \code{dvhmatrix} object
#' @param ref.frac Reference fractionation, default value 2
#' @param nf Fractions number for the input dvh
#' @param alphabeta \eqn{\alpha\beta} ratio for computing LQ cvonversion, default value 3
#' @return A \code{dvhmatrix} object with converted dose bins
#' @export
DVH.lq.correct<-function(dvh, ref.frac = 2, nf, alphabeta = 3) {
  dvh@dvh[,1] <- dvh@dvh[,1] * (alphabeta + dvh@dvh[,1]/nf) / (alphabeta + ref.frac)
  return(dvh)
}

#' Function for fitting Vdose and Dvolume to a given \code{dvhmatrix} object and outcome vector
#' @description This function can be used to find the best \emph{Vdose} or \emph{Dvolume} values that describe
#' a given clinically observed outcome 
#' @details The \emph{Vdose} is the value of the volume of a structure in a dose volume histogram that receives
#' a given level of dose. It is a widely used indicator and predictor of possible outcome when DVHs
#' have to be clinically evaluated. There are many papers that exploit the correct use of \emph{Vdose} for
#' clinical evaluation of treatment plans. The \code{DR.fit.DoseVolume} function used inside the \pkg{moddicom} package
#' allows to fit clinical data with \code{dvhmatrix} objects in order to detect the dose-volume relationship achievable
#' by dosimetric data. The dose-volume-outcome fitting is performed by iterative search, so, as far as the author know,
#' this is the first example of dose-volume-outcome fitting function that allows to find out \emph{the best fitting value}
#' of \emph{Vdose} or \emph{Dvolume}, rather than using an empyrical search from pre-fetermined values as usually reported in literature.
#' When setting the \code{model} parameter the users have to chose between a \code{\link[stats]{glm}} model and a \code{\link[survival]{coxph}} model
#' by writing the formula of the model \emph{without} the \code{dvhmatrix} object
#' @param dvh A \code{dvhmatrix} object
#' @param outcome A vector of binary values (if not it will coerced to binary values) showing the observed outcome
#' @param model The model that will be used to fit the dose-volume-response data (see examples). Models with available 
#' the \code{\link[stats]{logLik}} method can be used belonging to classes \code{lm}, \code{glm} and \code{coxph}
#' @param type A character value representing the two type of dose-volume fitting as \code{Vdose} or \code{Dvolume}
#' @param CI Logcal Value for calculating the confidence interval of the dose-volume fitting parameter and for the \eqn{\alpha\beta} if modLQ is optioned \code{TRUE}
#' @param CI.width Width of confidence interval
#' @param epsilon Error limit for calculating the Log Likelihood function in determining the confidence interval of dose-volume parameter
#' @examples ## NOT RUN DR.fit.DoseVolume(dvh = D, outcome = Dout, model = glm(formula = Dout ~ 1 , family = binomial(link = 'logit')), type = 'Vdose')
#' @references \emph{Special Considerations Regarding Absorbed-Dose and Dose-Volume Prescribing and Reporting in IMRT}. J ICRU. 2010 Apr;10(1):27-40. doi: 10.1093/jicru/ndq008. PubMed PMID: 24173325.
#' @references Graham MV, Purdy JA, Emami B, Harms W, Bosch W, Lockett MA, Perez CA. \emph{Clinical dose-volume histogram analysis for pneumonitis after 3D treatment for non-small cell lung cancer (NSCLC)}. Int J Radiat Oncol Biol Phys. 1999 Sep 1;45(2):323-9. PubMed PMID: 10487552.
#' @export
DR.fit.DoseVolume<-function(dvh, outcome, model, type = c("Vdose", "Dvolume"), CI = FALSE, CI.width = 0.95, epsilon = 1e-6) {
  type <- match.arg(arg = type)
  if (dvh@dvh.type == "differential") dvh<-DVH.diff.to.cum(dvh = dvh) # transform in cumulative if differential
  
  if (type == "Vdose")   { # create the approxfun and the fitting Vdose function (FUN)
    value<-seq(from = 1, to = max(floor(dvh@dvh[,1])), by = .5) # create the dose vector
    apf <- apply(X = dvh@dvh[, 2:ncol(dvh@dvh)], MARGIN = 2, FUN = approxfun, x = dvh@dvh[, 1]) # approxfun list
    update.model<-function(model, Vdose) {  # function for interactively update the model
      Vdose <<- Vdose                       # .GlobalEnv needed for scoping of variable
      new.model<-update(object = model, formula. = ". ~ . + Vdose")
      return(new.model)    
    }     
    f <- function(val)   return(sapply(X = apf, FUN = function(x) return(x(val))))  # function that calculates the Vdoses
    output.matrix<-c()
    for (n in 1:length(value)) { 
      Vd<-f(value[n])      
      temp.model<-update.model(model = model, Vdose = Vd)
      coeff<-summary(temp.model)$coefficients[,4]      
      if ( length(coeff) == length(temp.model$coefficients) ) # check length to prevent joining models with only one covariate
        output.matrix<-rbind(output.matrix, c(value[n], -logLik(temp.model), coeff)) else 
          warning(paste("Vdose =", value[n], "Gy doesn't fit a model"))
    }     
    
    obj.FUN<-function(x) {   # objective function for optimizing the best model
      Vdose<<-f(x[1])
      new.model<-update(object = model, formula. = ". ~ . + Vdose")
      return(-logLik(new.model))
    }    
    
    opt.model<-nlminb(start = output.matrix[which(output.matrix[,2] == min(output.matrix[,2])), 1], objective = obj.FUN)
    optimized.model<-update.model(model = model, Vdose = f(opt.model$par[1]))  
    coeff<-summary(optimized.model)$coefficients[,4]
    
    output.matrix<-rbind(output.matrix, c(opt.model$par[1], -logLik(optimized.model), coeff))
    output.matrix<-output.matrix[order(output.matrix[,1]),]
    rm(Vdose, envir = .GlobalEnv)
    fit.par<-opt.model$par[1]    
  }
  
  if (type == "Dvolume") { # create the approxfun and the fitting Dvolume function (FUN)
    maxVol<-min(apply(X = dvh@dvh[,2:ncol(dvh@dvh)], FUN = max, MARGIN = 2)) # the minimum of max Volume
    minVol<-max(apply(X = dvh@dvh[,2:ncol(dvh@dvh)], FUN = min, MARGIN = 2)) # the maximum of min Volume
    # browser()
    # dVol<-(maxVol - minVol)/200  # delta Volume
    value <- dvh@dvh[2:(nrow(dvh@dvh) - 1), 1]  # create volume vector removing extreme values
    apf <- apply(X = dvh@dvh[2:(nrow(dvh@dvh) - 1), 2:ncol(dvh@dvh)], MARGIN = 2, FUN = approxfun, x = value)    
    update.model<-function(model, Dvolume) {
      Dvolume <<- Dvolume
      new.model<-update(object = model, formula. = ". ~ . + Dvolume")
      return(new.model)
    }
    browser()
    f <- function(val)   return(sapply(X = apf, FUN = function(x) return(x(val))))  
    output.matrix<-c()
    for (n in 1:length(value)) {
      Dv<-f(value[n])
      temp.model<-update.model(model = model, Dvolume = Dv)
      coeff<-summary(temp.model)$coefficients[,4]      
      if ( length(coeff) == length(temp.model$coefficients) ) # check length to prevent joining models with only one covariate
        output.matrix<-rbind(output.matrix, c(value[n], -logLik(temp.model), coeff)) else 
          warning(paste("Dvolume =", value[n], "doesn't fit a model"))      
    }
    obj.FUN<-function(x) {   # objective function for optimizing the best model
      Dvolume<<-f(x[1])
      new.model<-update(object = model, formula. = ". ~ . + Dvolume")
      return(-logLik(new.model))
    }
    
    opt.model<-nlminb(start = output.matrix[which(output.matrix[,2] == min(output.matrix[,2])), 1], objective = obj.FUN)
    optimized.model<-update.model(model = model, Dvolume = f(opt.model$par[1]))  
    output.matrix<-rbind(output.matrix, c(opt.model$par[1], -logLik(optimized.model), summary(optimized.model)$coefficients[,4]))
    output.matrix<-output.matrix[order(output.matrix[,1]),]
    rm(Dvolume, envir = .GlobalEnv)
    fit.par<-opt.model$par[1]
  }    
  names(fit.par)<-type
  colnames(output.matrix)<-c("Dose", "-LL", names(summary(optimized.model)$coefficients[,4]))
  # calculating confidence interval for Dvolume and Vdose
  CI.val <- NULL
  if (CI == TRUE) { 
    bound <- qchisq(CI.width,1)/2
    ### confidence interval for Vdose
    if (type == "Vdose") {
      #### start caluclation for low bound of CI ####
      start.value<-fit.par  # start with optimal parameter value
      delta.value<-(fit.par - value[1])/2  # delta for steps going down in the CI boundary
      while ( delta.value > epsilon ) {
        start.value <- start.value - delta.value  # create step for the value to be fitted
        while (( start.value >= 0 ) && ((logLik(object = update.model(model = model, Vdose = f(start.value))) -logLik(optimized.model) + bound) > 0 )) 
          start.value <- start.value - delta.value
        start.value <- start.value + delta.value
        delta.value <- delta.value / 2
      }
      low.bound.CI <- start.value - delta.value / 2
      
      #### start calculation for highbound of CI ####
      start.value<-fit.par  # start with optimal parameter value
      delta.value<-(value[length(value)] - fit.par)/2  # delta for steps going up in the CI boundary
      while ( delta.value > epsilon ) {
        start.value <- start.value + delta.value  # create step for the value to be fitted
        while (( start.value <= value[length(value)]) && ((logLik(object = update.model(model = model, Vdose = f(start.value))) -logLik(optimized.model) + bound) > 0 )) 
          start.value <- start.value + delta.value
        start.value <- start.value - delta.value
        delta.value <- delta.value / 2
      }
      high.bound.CI <- start.value + delta.value / 2
    }
    
    ### confidence interval for Dvolume
    if (type == "Dvolume") {
      #### start caluclation for low bound of CI ####
      start.value<-fit.par  # start with optimal parameter value
      delta.value<-(fit.par - value[1])/2  # delta for steps going down in the CI boundary
      while ( delta.value > epsilon ) {
        start.value <- start.value - delta.value  # create step for the value to be fitted
        while (( start.value >= 0 ) && ((logLik(object = update.model(model = model, Dvolume = f(start.value))) -logLik(optimized.model) + bound) > 0))
          start.value <- start.value - delta.value
        start.value <- start.value + delta.value
        delta.value <- delta.value / 2
      }
      low.bound.CI <- start.value - delta.value / 2
      
      #### start calculation for highbound of CI ####
      start.value<-fit.par  # start with optimal parameter value
      delta.value<-(value[length(value)] - fit.par)/2  # delta for steps going up in the CI boundary
      while ( delta.value > epsilon ) {
        start.value <- start.value + delta.value  # create step for the value to be fitted
        while (( start.value <= value[length(value)]) && ((logLik(object = update.model(model = model, Dvolume = f(start.value))) -logLik(optimized.model) + bound) > 0 )) 
          start.value <- start.value + delta.value
        start.value <- start.value - delta.value
        delta.value <- delta.value / 2
      }
      high.bound.CI <- start.value + delta.value / 2
    }
    CI.val <- c(low.bound.CI, high.bound.CI)
    names(CI.val) <- c(paste((1-CI.width)/2, "%"), paste((1+CI.width)/2, "%"))
  }
  
  result <- new('DoseVolumeModel')
  result@output.matrix <- output.matrix
  result@fitted.model <- optimized.model
  result@fitted.value <- fit.par
  result@CI <- CI.val
  result@fitted.parameter <- type
  # return(list(output.matrix = output.matrix, optimized.model = optimized.model, fit.par = fit.par, CI = CI.val))
  return(result)
}
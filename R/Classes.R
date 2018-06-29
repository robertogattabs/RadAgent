#' dvmmatrix S4 class
#' @description This class handles the matrices containing the DVH(s)
#' @slot dvh The numeric matrix containing the DVH(s). The 1st column represents the dose steps, the other(s) represent(s) the volume values.
#' @slot dvh.type A character slot giving the \code{"cumulative"} or \code{"differential"} representation of the volume(s).
#' @slot vol.distr A character slot setting the \code{"absolute"} or \code{"relative"} value of the volume(s).
#' @slot volume A numeric vector giving the value(s) of the total volume(s) of the structure(s) in the DVH(s).
#' @exportClass dvhmatrix

setClass("dvhmatrix", 
         representation(
           dvh="matrix",                # matrix of DVHs, 1st column is dose, other columns are the volume
           dvh.type="character",        # cumulative or differential DVH
           vol.distr="character",       # absolute or relative distribution of the volume
           volume="numeric"             # vector of absolute volumes of DVH
           )
         )

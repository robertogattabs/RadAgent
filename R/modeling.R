#' Function for generating a series of outcomes related to a \code{dvhmatrix} class object
#' @description This function can be used for simulating outcome related to a given vector of doses or
#' a \code{dvhmatrix} class object after setting the outcome related parameters \eqn{TD_{50}}{TD50}, \eqn{\gamma_{50}}{\gamma 50}
#' and \eqn{a}. If \code{doses} isn't a \code{dvhmatrix} class object and it is a vector of nominal doses
#' the \code{a} parameter is ignored.
#' @param doses a \code{dvhmatrix} class object or a vector of nominal doses
#' @param a factor for calculating EUD value according \code{\link{DVH.eud}} function
#' @param TD50 Dose that gives the 50\% probability of outcome
#' @param gamma50 Slope of dose-response curve at TD50 of administered dose
#' @param DR.fun Dose/Response function, a character vector containing the name of one of the function in the package \pkg{moddicom}:
#' \code{Lyman}, \code{Niemierko}, \code{Bentzen}, \code{Goitein}, \code{Munro}, \code{Okunieff}, \code{Warkentin}.
#' @export
#' @return A vector of binary events (1 or 0).
DR.generate<-function(doses, a = 1, TD50 = 45, gamma50 = 1.5, DR.fun = c("Lyman", "Niemierko", "Bentzen", "Goitein",
                                                                         "Munro", "Okunieff", "Warkentin")) {
  DR.fun<-match.arg(arg = DR.fun)
  FUN<-match.fun(FUN = paste("DR.", DR.fun, sep = ""))
  p <-FUN(TD50 = TD50, gamma50 = gamma50, a = a, doses = doses)
  # m is the number of the cases
  if ((class(doses)=="numeric") || (class(doses)=="integer")) m<-length(doses)
  if (class(doses)=="dvhmatrix") m<-ncol(doses@dvh) - 1
  # seed for probability calculation
  prob<-runif(n = m)
  return(as.numeric(FUN(doses = doses, TD50 = TD50, gamma50 = gamma50, a = a)>prob))
}

## Dose/Response according Lyman 1985, Kutcher and Burman 1989, Deasy 2000 (probit) ##
#' Function that calculates NTCP according Lyman/Kutcher/Burman model
#' @description This function calculates the Normal Tissue Complication Probability according the
#' Lyman/Kutcher/Burman model. 
#' @details This model is a reformulation of probit model by using other parameters than mean and standard deviation.
#' In its original formulation the NTCP is given by the following equation: \deqn{NTCP=1/\sqrt{2\pi}\int_{-\infty}^{t}exp(-x^2/2)dx}
#' where: \deqn{v=\frac{V}{V_{ref}}} \deqn{t=(D-TD_{50}(v))/(m*TD_{50}(v))} \deqn{TD(v)=TD(1)*v-n} In previous equations \eqn{V} is the
#' irradiated fraction volume and \eqn{V_{ref}} is the referenced volume for the given outcome, \eqn{m} is another way
#' to describe the slope of dose-response curve (see following lines).
#' In the \pkg{moddicom} implementation the parameters to be set in the model are \eqn{TD_{50}} and \eqn{\gamma_{50}}.
#' The model has been coded by adapting the original formula using these parameters in this way:
#' \deqn{t=(EUD-TD_{50})/(m*TD_{50})} with:
#' \deqn{m=\frac{1}{\gamma_{50}\sqrt{2\pi}}} The relationship between Dose and Volume has been achieved by using the
#' equivalent uniform dose as calculated by \code{\link{DVH.eud}} function.
#' @param doses Either a \code{dvhmatrix} class object or a vector with nominal doses
#' @param TD50 The value of dose that gives the 50\% of probability of outcome
#' @param gamma50 The slope of dose/response curve at 50\% of probability
#' @param a Value for parallel-serial correlation in radiobiological response
#' @export
#' @return A vector with NTCP(s) calculated according LKB model.
#' @references Burman C, Kutcher GJ, Emami B, Goitein M. \emph{Fitting of normal tissue tolerance data to an analytic function}. Int J Radiat Oncol Biol Phys. 1991 May 15;21(1):123-35. PubMed PMID: 2032883.
DR.Lyman <- function (doses, TD50 = 45, gamma50 = 1.5, a = 1) {
  if ((class(doses)=="numeric") || (class(doses)=="integer")) p <- pnorm(q=((doses - TD50)*gamma50*sqrt(2*pi))/TD50)
  if (class(doses)=="dvhmatrix") p <- pnorm(q=((DVH.eud(dvh = doses, a=a) - TD50)*gamma50*sqrt(2*pi))/TD50)
  return(p)
}

## Dose/Response according Goitein 1979 (logprobit) ##
#' Function that calculates NTCP according Goitein (Bentzen) model
#' @description This function calculates the Normal Tissue Complication Probability according the
#' Goitein (Bentzen) model. 
#' @details This model is similar to Lyman model but it is function of \emph{logarithm} of the Dose.
#' Starting from the implementation of Lyman model: \deqn{NTCP=1/\sqrt{2\pi}\int_{-\infty}^{t}exp(-x^2/2)dx}
#' where: \deqn{v=\frac{V}{V_{ref}}} \deqn{t=(log(D)-TD_{50}(v))/(m*TD_{50}(v))} \deqn{TD(v)=TD(1)*v-n} In previous equations \eqn{V} is the
#' irradiated fraction volume and \eqn{V_{ref}} is the referenced volume for the given outcome, \eqn{m} is another way
#' to describe the slope of dose-response curve (see following lines).
#' In the \pkg{moddicom} implementation the parameters to be set in the model are \eqn{TD_{50}} and \eqn{\gamma_{50}}.
#' The model has been coded by adapting the original formula using these parameters in this way:
#' \deqn{t=(log(EUD)-TD_{50})/(m*TD_{50})} with:
#' \deqn{m=\frac{1}{\gamma_{50}\sqrt{2\pi}}} 
#' \eqn{D} can be either the nominal dose or the \eqn{EUD} as calculated by \code{\link{DVH.eud}} function.
#' @param doses Either a \code{dvhmatrix} class object or a vector with nominal doses
#' @param TD50 The value of dose that gives the 50\% of probability of outcome
#' @param gamma50 The slope of dose/response curve at 50\% of probability
#' @param a Value for parallel-serial correlation in radiobiological response
#' @export
#' @return A vector with NTCP(s) calculated according Goitein model.
#' @references Shipley WU, Tepper JE, Prout GR Jr, Verhey LJ, Mendiondo OA, Goitein M, Koehler AM, Suit HD. \emph{Proton radiation as boost therapy for localized prostatic carcinoma}. JAMA. 1979 May 4;241(18):1912-5. PubMed PMID: 107338.
#' @references Bentzen SM, Tucker SL. \emph{Quantifying the position and steepness of radiation dose-response curves}. Int J Radiat Biol. 1997 May;71(5):531-42. Review. PubMed PMID: 9191898.
DR.Goitein <- function (doses, TD50=45, gamma50=1.5, a=1) {
  if ((class(doses)=="numeric") || (class(doses)=="integer")) p <- pnorm(q=(log(doses/TD50)*gamma50*sqrt(2*pi)))
  if (class(doses)=="dvhmatrix") p <- pnorm(q=(log(DVH.eud(dvh = doses, a=a)/TD50)*gamma50*sqrt(2*pi)))
  return(p)
}

## Dose/Response according Niemierko 1991 (loglogit) ##
#' Function that calculates (N)TCP according Niemierko model
#' @description This function calculates the Normal Tissue Complication Probability according the
#' Niemierko model. 
#' @details This model can be used to compute Tumor Control Probability (TCP) too.
#' It is the translation of a \emph{loglogit} generalized linear model as function of \eqn{TD_{50}} and \eqn{\gamma_{50}}.
#' The model equation is: \deqn{(N)TCP=\frac{1}{1+ ({\frac{TD_{50}}{D}})^{4\gamma_{50}}}}
#' \eqn{D} can be either the nominal dose or the \eqn{EUD} as calculated by \code{\link{DVH.eud}} function.
#' @param doses Either a \code{dvhmatrix} class object or a vector with nominal doses
#' @param TD50 The value of dose that gives the 50\% of probability of outcome
#' @param gamma50 The slope of dose/response curve at 50\% of probability
#' @param a Value for parallel-serial correlation in radiobiological response 
#' @export
#' @return A vector with NTCP(s) calculated according Niemierko model.
#' @references Gay HA, Niemierko A. \emph{A free program for calculating EUD-based NTCP and TCP in external beam radiotherapy}. Phys Med. 2007 Dec;23(3-4):115-25. Epub 2007 Sep 7. PubMed PMID: 17825595.
DR.Niemierko <- function (doses, TD50=45, gamma50=1.5, a=1) {
  if ((class(doses)=="numeric") || (class(doses)=="integer")) p <- 1/(1+(TD50/doses)^(4*gamma50))
  if (class(doses)=="dvhmatrix") p <- 1/(1+(TD50/DVH.eud(dvh = doses, a=a))^(4*gamma50))
  return(p)
}

## Dose/Response according Munro, Gilbert, Kallman 1992 (Poisson approximation) ##
#' Function that calculates TCP according Munro/Gilbert/Kallman model
#' @description This function calculates the Tumor Control Probability according the
#' Munro/Gilbert/Kallman model. 
#' @details This model is an empyrical dose/response curve that fits experimental data. In their
#' paper authors assume this curve to be equivalent to a Poisson model. The original model equation is:
#' \deqn{TCP=e^{-EN_{0}e^{\frac{-D}{D_{0}}}}}
#' \eqn{E} is a numerical parameter that is related to tumor radiosensitivity, \eqn{N_{0}} is the total initial number of
#' tumor clonogenic cells, \eqn{D} is the delivered dose and \eqn{D_{0}} is the increment of dose that lowers survival 
#' to 37 per cent. In our implementation Munro/Gilbert/Kallman model has been referenced to \eqn{TD_{50}} and \eqn{\gamma_{50}} as follows:
#' \deqn{TCP=2^{e^{e\gamma_{50}(1-\frac{D}{TD_{50}})}}}
#' In the model equation \eqn{D} can be either the nominal dose or the \eqn{EUD} as calculated by \code{\link{DVH.eud}} function.
#' @param doses Either a \code{dvhmatrix} class object or a vector with nominal doses
#' @param TD50 The value of dose that gives the 50\% of probability of outcome
#' @param gamma50 The slope of dose/response curve at 50\% of probability
#' @param a Value for parallel-serial correlation in radiobiological response 
#' @export
#' @return A vector with TCP calculated according Munro/Gilbert/Kallman model.
#' @references Munro TR, Gilbert CW. \emph{The relation between tumour lethal doses and the radiosensitivity of tumour cells}. Br J Radiol. 1961 Apr;34:246-51. PubMed PMID: 13726846.
DR.Munro <- function (doses, TD50=45, gamma50=1.5, a=1) {
  if ((class(doses)=="numeric") || (class(doses)=="integer")) p <- 2^(-(exp(exp(1)*gamma50*(1-doses/TD50))))
  if (class(doses)=="dvhmatrix") p <- 2^(-(exp(exp(1)*gamma50*(1-DVH.eud(dvh = doses, a=a)/TD50))))
  return(p) 
}

## Dose/Response according Okunieff 1995 (logit) ##
#' Function that calculates TCP according Okunieff model
#' @description This function calculates the Tumor Control Probability according the
#' Okunieff model.
#' @details This model is the equivalent of the \emph{logistic} generalized linear model where the covariates and their coefficients
#' have been reported as function of \eqn{TD_{50}} and \eqn{\gamma_{50}}. The original Okunieff formula is the following:
#' \deqn{TCP=\frac{e^{\frac{D-TD_{50}}{k}}}{1+e^{\frac{D-TD_{50}}{k}}}}
#' where \eqn{k=\gamma_{50}/(4*TD_{50})} and so giving the final model as direct function of \eqn{TD_{50}} and \eqn{\gamma_{50}}:
#' \deqn{TCP=\frac{1}{1+e^{4\gamma_{50}(1-\frac{D}{TD_{50}})}}}
#' In the model equation \eqn{D} can be either the nominal dose or the \eqn{EUD} as calculated by \code{\link{DVH.eud}} function.
#' @param doses Either a \code{dvhmatrix} class object or a vector with nominal doses
#' @param TD50 The value of dose that gives the 50\% of probability of outcome
#' @param gamma50 The slope of dose/response curve at 50\% of probability
#' @param a Value for parallel-serial correlation in radiobiological response 
#' @export
#' @return A vector with TCP calculated according Munro/Gilbert/Kallman model.
#' @references Okunieff P, Morgan D, Niemierko A, Suit HD. \emph{Radiation dose-response of human tumors}. Int J Radiat Oncol Biol Phys. 1995 Jul 15;32(4):1227-37. PubMed PMID: 7607946.
DR.Okunieff <- function (doses, TD50=45, gamma50=1.5, a=1) {
  if ((class(doses)=="numeric") || (class(doses)=="integer")) p <- 1/(1+exp(4*gamma50*(1-(doses/TD50))))
  if (class(doses)=="dvhmatrix") p <- 1/(1+exp(4*gamma50*(1-(DVH.eud(dvh = doses, a=a)/TD50))))
  return(p) 
}

## Dose/Response according Warkentin 2004 (Poisson) ##
#' Function that calculates TCP according Warkentin (Poisson) model
#' @description This function calculates the Tumor Control Probability according the
#' Warkentin model. 
#' @details This model is the equivalent of the \emph{Poisson} model where the covariates and their coefficients
#' have been reported as function of \eqn{TD_{50}} and \eqn{\gamma_{50}}. In Warkentin paper model computation starts with formula:
#' \deqn{TCP=e^{-Np_{s}(D)}}
#' where \eqn{N} is the initial number of clonogens, \eqn{p_{s}(D)} is the survival fraction after the dose \eqn{D}.
#' The previous equation can be rewritten as function of \eqn{TD_{50}} and \eqn{\gamma_{50}}:
#' \deqn{TCP=\left (\frac{1}{2}  \right )^{e^{\left [2\gamma_{50}(1-D/D_{50})/ln2  \right ]}}}
#' Using the assumption of independent subvolumes, for the case of heterogeneous irradiation, 
#' the overall probability of tumor control is the product of the probabilities of killing all 
#' clonogens in each tumor subvolume described by the differential DVH:
#' \deqn{TCP=\prod_{i}TCP(D_{i},v_{i})}
#' Thus, for a given differential DVH \eqn{\left \{ D_{i},v_{i} \right \}}, the TCP can be calculated using the following two-parameter TCP formula:
#' \deqn{TCP=\left (\frac{1}{2}  \right )^{\sum_{i}v_{i}e^{\left [2\gamma_{50}(1-D/D_{50})/ln2  \right ]}}}
#' In the model equation \eqn{D} can be either the nominal dose or the \eqn{EUD} as calculated by \code{\link{DVH.eud}} function.
#' @param doses Either a \code{dvhmatrix} class object or a vector with nominal doses
#' @param TD50 The value of dose that gives the 50\% of probability of outcome
#' @param gamma50 The slope of dose/response curve at 50\% of probability
#' @param a Value for parallel-serial correlation in radiobiological response 
#' @export
#' @return A vector with TCP calculated according Warkentin (Poisson) model.
#' @references Warkentin B, Stavrev P, Stavreva N, Field C, Fallone BG. \emph{A TCP-NTCP estimation module using DVHs and known radiobiological models and parameter sets}. J Appl Clin Med Phys. 2004 Winter;5(1):50-63. Epub 2004 Jan 1. PubMed PMID: 15753933.
DR.Warkentin <- function (doses, TD50=45, gamma50=1.5, a=1) {
  if ((class(doses)=="numeric") || (class(doses)=="integer")) p <- 0.5^(exp(2*gamma50/0.693147181*(1-doses/TD50)))
  if (class(doses)=="dvhmatrix") p <- 0.5^(exp(2*gamma50/0.693147181*(1-DVH.eud(dvh = doses, a=a)/TD50)))
  return(p) 
}

## Dose/Response according Bentzen 1997 (log Poisson) ##
#' Function that calculates TCP according Bentzen (log Poisson) model
#' @description This function calculates the Tumor Control Probability according the
#' Bentzen (log Poisson) model.
#' @details In Bentzen paper dose/response curves were described both for tumor and normal tissue response.
#' In \pkg{moddicom} implementation \emph{log Poisson} model starts from the \code{\link{DR.Warkentin}} equation but 
#' calculating the relationship of response with the \emph{log} Dose as follows:
#' \deqn{TCP=\left (\frac{1}{2}  \right )^{\left (\frac{TD_{50}}{D}  \right )^{\frac{2\gamma_{50}}{log2}}}}
#' In the model equation \eqn{D} can be either the nominal dose or the \eqn{EUD} as calculated by \code{\link{DVH.eud}} function.
#' @param doses Either a \code{dvhmatrix} class object or a vector with nominal doses
#' @param TD50 The value of dose that gives the 50\% of probability of outcome
#' @param gamma50 The slope of dose/response curve at 50\% of probability
#' @param a Value for parallel-serial correlation in radiobiological response 
#' @export
#' @references Bentzen SM, Tucker SL. \emph{Quantifying the position and steepness of radiation dose-response curves}. Int J Radiat Biol. 1997 May;71(5):531-42. Review. PubMed PMID: 9191898.
DR.Bentzen <- function (doses, TD50=45, gamma50=1.5, a=1) {
  if ((class(doses)=="numeric") || (class(doses)=="integer")) p <- 0.5^(TD50/doses)^(2*gamma50/0.693147181)
  if (class(doses)=="dvhmatrix") p <- 0.5^(TD50/DVH.eud(dvh = doses, a=a))^(2*gamma50/0.693147181)
  return(p) 
}

#' Fit a Dose/Response function
#' @description This function fits an object of class \code{dvhatrix} or \code{numeric} to a given series of
#' outcome.
#' @details This function is a wrapper for fitting dose/response functions to a given set of binary outcomes. Currently
#' outcomes must be referred as a numeric vector  representing cases "showing the outcome" with the number 1, 
#' or "not showing the outcome" with the number 0. Model fitting is performed using Maximum Likelihood Estimation method,
#' implemented by wrapping R MLE functions set for specific Log Likelihood functions related to dose/response formulas.
#' If \code{doses} is a vector of nominal doses its length must be equal to the length of \code{outcome}, if \code{doses} is
#' a \code{dvhmatrix} class object the number of histograms must be equal to the length of \code{outcome}.
#' @param doses Either a \code{dvhmatrix} class object or a vector with nominal doses
#' @param outcome A numeric vector of cases showing (1) and not showing (0) the outcome
#' @param DR.fun Dose/Response function, a character vector containing the name of one of the function in the package \pkg{moddicom}:
#' \code{Lyman}, \code{Niemierko}, \code{Bentzen}, \code{Goitein}, \code{Munro}, \code{Okunieff}, \code{Warkentin}.
#' @param type Function type: \code{NTCP}, Normal Tissue Complication Probability, or \code{TCP}, Tumor Control Probability
#' @param CI If \code{TRUE} it returns the value of confidence interval calulated by profile likelihood method
#' @param CI.width The value of width of confidence interval to be returned if \code{CI = TRUE}
#' @param epsilon Lower interval for bisecting the MLE
#' @export
DR.fit <- function (doses, outcome, DR.fun = c("Lyman", "Niemierko", "Bentzen", "Goitein", "Munro", "Okunieff", "Warkentin"),
                    type = c("NTCP", "TCP"), CI = TRUE, C.I.width = .95, epsilon = 1e-6) {
  type<-match.arg(type)
  DR.fun<-match.arg(DR.fun)
  ## fitting two parameters dose/response model
  if ((class(doses)=="numeric") || (class(doses)=="integer")) {
    ## Akaike Information Criterion function
    DR.AIC  <-function(LL) return(4 - 2 * LL)
    ## Akaike Information Criterion corrected by the cases number
    DR.AICc <-function() return(aic + 12/(length(doses) - 3))
    ## define LL functions
    if (DR.fun=="Lyman")      FUN<-function(doses, TD50, gamma50) return(pnorm(q=((doses - TD50)*gamma50*sqrt(2*pi))/TD50))
    if (DR.fun=="Niemierko")  FUN<-function(doses, TD50, gamma50) return(1/(1+(TD50/doses)^(4*gamma50)))
    if (DR.fun=="Bentzen")    FUN<-function(doses, TD50, gamma50) return(0.5^(TD50/doses)^(2*gamma50/0.693147181))
    if (DR.fun=="Goitein")    FUN<-function(doses, TD50, gamma50) return(pnorm(q=(log(doses/TD50)*gamma50*sqrt(2*pi))))
    if (DR.fun=="Munro")      FUN<-function(doses, TD50, gamma50) return(2^(-(exp(exp(1)*gamma50*(1-doses/TD50)))))
    if (DR.fun=="Okunieff")   FUN<-function(doses, TD50, gamma50) return(1/(1+exp(4*gamma50*(1-(doses/TD50)))))
    if (DR.fun=="Warkentin")  FUN<-function(doses, TD50, gamma50) return(0.5^(exp(2*gamma50/0.693147181*(1-doses/TD50))))
    ## define the negative LL function
    nLL<-function(par, doses, outcome){
      TD50<-par[1]
      gamma50<-par[2]
      return(-sum(outcome*log(FUN(doses, TD50, gamma50))+(1 - outcome)*log(1 - FUN(doses, TD50, gamma50))))
    }
    ## LL functions for Profile Likelihood
    nLL.TD50<-function(par, doses, outcome, TD50){
      gamma50<-par[1]
      return(-sum(outcome*log(FUN(doses, TD50, gamma50))+(1 - outcome)*log(1 - FUN(doses, TD50, gamma50))))
    }
    nLL.gamma50<-function(par, doses, outcome, gamma50){
      TD50<-par[1]
      return(-sum(outcome*log(FUN(doses, TD50, gamma50))+(1 - outcome)*log(1 - FUN(doses, TD50, gamma50))))
    }
    # check the function type
    if (type=="NTCP")
      if (any(DR.fun==c("Bentzen", "Munro", "Okunieff", "Warkentin"))) warning(paste("Trying to use", DR.fun, "model as NTCP function!"))
    if (type=="TCP")
      if (any(DR.fun==c("Goitein", "Lyman"))) warning(paste("Trying to use", DR.fun, "model as TCP function!"))
    fit<-nlminb(start = c(45, 1.5), objective = nLL, lower = c(10, .01), upper = c(150, 5), doses = doses, outcome = outcome)
  }
  
  ## fitting three parameters dose/response model
  if (class(doses)=="dvhmatrix") {
    doses<-DVH.relative(dvh = DVH.cum.to.diff(dvh = doses))
    # Akaike Information Criterion function
    DR.AIC <-function(LL) return(6 - 2 * LL)
    ## Akaike Information Criterion corrected by the cases number
    DR.AICc <-function() return(aic + 24/(ncol(doses@dvh) - 5))
    ## define LL functions
    if (DR.fun=="Lyman")      FUN<-function(doses, TD50, gamma50, a) return(pnorm(q=((DVH.eud(dvh = doses, a = a) - TD50)*gamma50*sqrt(2*pi))/TD50))
    if (DR.fun=="Niemierko")  FUN<-function(doses, TD50, gamma50, a) return(1/(1+(TD50/DVH.eud(dvh = doses, a = a))^(4*gamma50)))
    if (DR.fun=="Bentzen")    FUN<-function(doses, TD50, gamma50, a) return(0.5^(TD50/DVH.eud(dvh = doses, a = a))^(2*gamma50/0.693147181))
    if (DR.fun=="Goitein")    FUN<-function(doses, TD50, gamma50, a) return(pnorm(q=(log(DVH.eud(dvh = doses, a = a)/TD50)*gamma50*sqrt(2*pi))))
    if (DR.fun=="Munro")      FUN<-function(doses, TD50, gamma50, a) return(2^(-(exp(exp(1)*gamma50*(1-DVH.eud(dvh = doses, a = a)/TD50)))))
    if (DR.fun=="Okunieff")   FUN<-function(doses, TD50, gamma50, a) return(1/(1+exp(4*gamma50*(1-(DVH.eud(dvh = doses, a = a)/TD50)))))
    if (DR.fun=="Warkentin")  FUN<-function(doses, TD50, gamma50, a) return(0.5^(exp(2*gamma50/0.693147181*(1-DVH.eud(dvh = doses, a = a)/TD50))))
    ## define the negative LL function
    nLL<-function(par, doses, outcome){
      TD50<-par[1]
      gamma50<-par[2]
      a<-par[3]
      return(-sum(outcome*log(FUN(doses, TD50, gamma50, a))+(1 - outcome)*log(1 - FUN(doses, TD50, gamma50, a))))
    }
    ## LL functions for Profile Likelihood
    nLL.TD50<-function(par, doses, outcome, TD50){
      gamma50<-par[1]
      a<-par[2]
      return(-sum(outcome*log(FUN(doses, TD50, gamma50, a))+(1 - outcome)*log(1 - FUN(doses, TD50, gamma50, a))))
    }
    nLL.gamma50<-function(par, doses, outcome, gamma50){
      TD50<-par[1]
      a<-par[2]
      return(-sum(outcome*log(FUN(doses, TD50, gamma50, a))+(1 - outcome)*log(1 - FUN(doses, TD50, gamma50, a))))
    }
    nLL.a<-function(par, doses, outcome, a){
      TD50<-par[1]
      gamma50<-par[2]
      return(-sum(outcome*log(FUN(doses, TD50, gamma50, a))+(1 - outcome)*log(1 - FUN(doses, TD50, gamma50, a))))
    }
    # check the function type
    if (type=="NTCP") {
      Upper<-c(150, 5, 40)
      Lower<-c(10, .01, .01)
      Start<-c(45, 1.5, 2)
      # print warning when fitting a TCP model to NTCP data
      if (any(DR.fun==c("Bentzen", "Munro", "Okunieff", "Warkwntin"))) warning(paste("Trying to use", DR.fun, "model as NTCP function!"))
    }
    if (type=="TCP") {
      Upper<-c(150, 5, -.01)
      Lower<-c(10, .01, -40)
      Start<-c(45, 1.5, 2)
      # print warning when fitting a NTCP model to TCP data
      if (any(DR.fun==c("Goitein", "Lyman"))) warning(paste("Trying to use", DR.fun, "model as TCP function!"))
    }
    fit<-nlminb(start = Start, objective = nLL, lower = Lower, upper = Upper, doses = doses, outcome = outcome)
    #fit<-optimx(par = Start, fn = nLL, lower = Lower, upper = Upper, doses = doses, outcome = outcome)
  }

  MLE<- -fit$objective
  aic <- DR.AIC(LL = MLE)
  aicc <- DR.AICc()
  # set the limit of C.I. with 1/2*chisquare for 1 degree of freedom
  bound <- qchisq(C.I.width,1)/2  
  # set NULL values of CI before calculating
  TD50low<-NULL
  TD50high<-NULL
  gamma50low<-NULL
  gamma50high<-NULL
  alow<-NULL
  ahigh<-NULL
  if (CI == TRUE) { 
    # check the kind of doses and start Profile Likelihood sequence
    if (class(doses)=="dvhmatrix") {
      TD50.CI <- function(start.TD50) {
        temp.fit<-nlminb(start = c(fit$par[2], fit$par[3]), objective = nLL.TD50, lower = Lower[2:3], upper = Upper[2:3], doses = doses, outcome = outcome, TD50 = start.TD50)
        return(list(MLE=-temp.fit$objective, gamma50 = temp.fit$par[1], a = temp.fit$par[2]))
      }
      
      message("Profiling TD50 low C.I. limit...")
      start.TD50<-fit$par[1]
      TD50.step <- 1
      while (TD50.step>epsilon) {
        start.TD50<- start.TD50 - TD50.step
        while ((TD50.CI(start.TD50)$MLE + fit$objective + bound)>0) {       
          #message("TD50=",start.TD50, "   dMLE=" ,TD50.CI(start.TD50)$MLE + fit$objective + bound )
          start.TD50 <- start.TD50 - TD50.step       
        }
        start.TD50<- start.TD50 + TD50.step
        TD50.step<-TD50.step/2
      }
      TD50low<-start.TD50-TD50.step/2
      
      message("Profiling TD50 high C.I. limit...")
      start.TD50<-fit$par[1]
      TD50.step<- 1
      while (TD50.step>epsilon) {
        start.TD50<- start.TD50 + TD50.step
        while ((TD50.CI(start.TD50)$MLE + fit$objective + bound)>0) {
          start.TD50 <- start.TD50 + TD50.step       
        }
        start.TD50<- start.TD50 - TD50.step
        TD50.step<-TD50.step/2
      }
      TD50high<-start.TD50+TD50.step/2
      
      gamma50.CI <- function(start.gamma50) {
        temp.fit<-nlminb(start = c(fit$par[1], fit$par[3]), objective = nLL.gamma50, lower = Lower[c(1,3)], upper = Upper[c(1,3)], doses = doses, outcome = outcome, gamma50 = start.gamma50)
        return(list(MLE=-temp.fit$objective, TD50 = temp.fit$par[1], a = temp.fit$par[2]))
      }
      
      message("Profiling gamma50 low C.I. limit...")
      start.gamma50<-fit$par[2]
      gamma50.step <- .1
      while (gamma50.step>epsilon) {
        start.gamma50<- start.gamma50 - gamma50.step
        while ((gamma50.CI(start.gamma50)$MLE + fit$objective + bound)>0) {       
          start.gamma50 <- start.gamma50 - gamma50.step       
        }
        start.gamma50<- start.gamma50 + gamma50.step
        gamma50.step<-gamma50.step/2
      }
      gamma50low<-start.gamma50-gamma50.step/2
      
      message("Profiling gamma50 high C.I. limit...")
      start.gamma50<-fit$par[2]
      gamma50.step<- .1
      while (gamma50.step>epsilon) {
        start.gamma50<- start.gamma50 + gamma50.step
        while ((gamma50.CI(start.gamma50)$MLE + fit$objective + bound)>0) {
          start.gamma50 <- start.gamma50 + gamma50.step       
        }
        start.gamma50<- start.gamma50 - gamma50.step
        gamma50.step<-gamma50.step/2
      }
      gamma50high<-start.gamma50+gamma50.step/2
      
      a.CI <- function(start.a) {
        temp.fit<-nlminb(start = c(fit$par[1], fit$par[2]), objective = nLL.a, lower = Lower[1:2], upper = Upper[1:2], doses = doses, outcome = outcome, a = start.a)
        return(list(MLE=-temp.fit$objective, TD50 = temp.fit$par[1], gamma50 = temp.fit$par[2]))
      }
      
      message("Profiling a low C.I. limit...")
      start.a<-fit$par[3]
      a.step <- .2
      while (a.step>epsilon) {
        start.a<- start.a - a.step
        while ((a.CI(start.a)$MLE + fit$objective + bound)>0) {   
          start.a <- start.a - a.step       
        }
        start.a<- start.a + a.step
        a.step<-a.step/2
      }
      alow <- start.a - a.step/2
      
      message("Profiling a high C.I. limit...")
      start.a<-fit$par[3]
      a.step<- 1
      while (a.step>epsilon) {
        start.a<- start.a + a.step
        while ((a.CI(start.a)$MLE + fit$objective + bound)>0) {
          start.a <- start.a + a.step       
        }
        start.a <- start.a - a.step
        a.step <- a.step/2
      }
      ahigh <- start.a + a.step/2
    }
    
    # modeling for numeric vector of doses
    if ((class(doses)=="numeric") || (class(doses)=="integer")) {
      TD50.CI <- function(start.TD50) {
        temp.fit<-nlminb(start = fit$par[2], objective = nLL.TD50, lower = .01, upper = 5, doses = doses, outcome = outcome, TD50 = start.TD50)
        return(list(MLE=-temp.fit$objective, gamma50 = temp.fit$par[1]))
      }
      
      message("Profiling TD50 low C.I. limit...")
      start.TD50<-fit$par[1]
      TD50.step <- 1
      while (TD50.step>epsilon) {
        start.TD50<- start.TD50 - TD50.step
        while ((TD50.CI(start.TD50)$MLE + fit$objective + bound)>0) {       
          #message("TD50=",start.TD50, "   dMLE=" ,TD50.CI(start.TD50)$MLE + fit$objective + bound )
          start.TD50 <- start.TD50 - TD50.step       
        }
        start.TD50<- start.TD50 + TD50.step
        TD50.step<-TD50.step/2
      }
      TD50low<-start.TD50-TD50.step/2
      
      message("Profiling TD50 high C.I. limit...")
      start.TD50<-fit$par[1]
      TD50.step<- 1
      while (TD50.step>epsilon) {
        start.TD50<- start.TD50 + TD50.step
        while ((TD50.CI(start.TD50)$MLE + fit$objective + bound)>0) {
          start.TD50 <- start.TD50 + TD50.step       
        }
        start.TD50<- start.TD50 - TD50.step
        TD50.step<-TD50.step/2
      }
      TD50high<-start.TD50+TD50.step/2
      
      gamma50.CI <- function(start.gamma50) {
        temp.fit<-nlminb(start = fit$par[1], objective = nLL.gamma50, lower = 10, upper = 150, doses = doses, outcome = outcome, gamma50 = start.gamma50)
        return(list(MLE=-temp.fit$objective, TD50 = temp.fit$par[1]))
      }
      
      message("Profiling gamma50 low C.I. limit...")
      start.gamma50<-fit$par[2]
      gamma50.step <- .1
      while (gamma50.step>epsilon) {
        start.gamma50<- start.gamma50 - gamma50.step
        while ((gamma50.CI(start.gamma50)$MLE + fit$objective + bound)>0) {       
          start.gamma50 <- start.gamma50 - gamma50.step       
        }
        start.gamma50<- start.gamma50 + gamma50.step
        gamma50.step<-gamma50.step/2
      }
      gamma50low<-start.gamma50-gamma50.step/2
      
      message("Profiling gamma50 high C.I. limit...")
      start.gamma50<-fit$par[2]
      gamma50.step<- .1
      while (gamma50.step>epsilon) {
        start.gamma50<- start.gamma50 + gamma50.step
        while ((gamma50.CI(start.gamma50)$MLE + fit$objective + bound)>0) {
          start.gamma50 <- start.gamma50 + gamma50.step       
        }
        start.gamma50<- start.gamma50 - gamma50.step
        gamma50.step<-gamma50.step/2
      }
      gamma50high<-start.gamma50+gamma50.step/2
    }    
  }
    
  return(list(fit = fit, AIC = aic, AICc = aicc, MLE = MLE, TD50low = TD50low, TD50high = TD50high, 
              gamma50low = gamma50low, gamma50high = gamma50high, alow = alow, ahigh = ahigh, model = DR.fun))
}

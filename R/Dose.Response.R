###########################################################################
## script generating dose/response models over two  over two parameters  ##
## (TD50, gamma50), or three parameters (TD50, gamma50, a)               ##
###########################################################################

## Dose/Response according Lyman 1985, Kutcher and Burman 1989, Deasy 2000 (probit) ##
DR.Lyman <- function (TD50=45, gamma50=1.5, aa=1, diffdvh=NULL, dose=NULL) {
  # check the single choice between dvh matrix or dose series
  if (!is.null(dose) && !is.null(diffdvh)) stop("Select either a DVH or a point dose to calculate NTCP")
  # check the single choice between dvh matrix or dose series 
  if (!is.null(dose)) p <- pnorm(q=((dose - TD50)*gamma50*sqrt(2*pi))/TD50)
  if (!is.null(diffdvh)) p <- pnorm(q=((EUD(dvh.matrix=diffdvh, a=aa) - TD50)*gamma50*sqrt(2*pi))/TD50)
  return(p)
}

## Dose/Response according Goitein 1979 (logprobit) ##
DR.Goitein <- function (TD50=45, gamma50=1.5, aa=1, diffdvh=NULL, dose=NULL) {
  # check the single choice between dvh matrix or dose series
  if (!is.null(dose) && !is.null(diffdvh)) stop("Select either a DVH or a point dose to calculate NTCP")
  # check the single choice between dvh matrix or dose series 
  if (!is.null(dose)) p <- pnorm(q=(log(dose/TD50)*gamma50*sqrt(2*pi)))
  if (!is.null(diffdvh)) p <- pnorm(q=(log(EUD(dvh.matrix=diffdvh, a=aa)/TD50)*gamma50*sqrt(2*pi)))
  return(p)
}

## Dose/Response according Niemierko 1991 (loglogit) ##
DR.Niemierko <- function (TD50=45, gamma50=1.5, aa=1, diffdvh=NULL, dose=NULL) {
  # check the single choice between dvh matrix or dose series
  if (!is.null(dose) && !is.null(diffdvh)) stop("Select either a DVH or a point dose to calculate NTCP")
  # check the single choice between dvh matrix or dose series  
  if (is.vector(dose) && is.null(diffdvh)) p <- 1/(1+(TD50/dose)^(4*gamma50))
  if (is.null(dose) && is.matrix(diffdvh)) p <- 1/(1+(TD50/EUD(dvh.matrix=diffdvh, a=aa))^(4*gamma50))
  return(p)
}

## Dose/Response according Munro, Gilbert, Kallman 1992 (Poisson approximation) ##
DR.Munro <- function (TD50=45, gamma50=1.5, aa=1, diffdvh=NULL, dose=NULL) {
  # check the single choice between dvh matrix or dose series
  if (!is.null(dose) && !is.null(diffdvh)) stop("Select either a DVH or a point dose to calculate NTCP")
  # check the single choice between dvh matrix or dose series  
  if (is.vector(dose) && is.null(diffdvh)) p <- 2^(-(exp(exp(1)*gamma50*(1-dose/TD50))))
  if (is.null(dose) && is.matrix(diffdvh)) p <- 2^(-(exp(exp(1)*gamma50*(1-EUD(dvh.matrix=diffdvh, a=aa)/TD50))))
  return(p) 
}

## Dose/Response according Okunieff 1995 (logit) ##
DR.Okunieff <- function (TD50=45, gamma50=1.5, aa=1, diffdvh=NULL, dose=NULL) {
  # check the single choice between dvh matrix or dose series
  if (!is.null(dose) && !is.null(diffdvh)) stop("Select either a DVH or a point dose to calculate NTCP")
  # check the single choice between dvh matrix or dose series  
  if (is.vector(dose) && is.null(diffdvh)) p <- 1/(1+exp(4*gamma50*(1-(dose/TD50))))
  if (is.null(dose) && is.matrix(diffdvh)) p <- 1/(1+exp(4*gamma50*(1-(EUD(dvh.matrix=diffdvh, a=aa)/TD50))))
  return(p) 
}

## Dose/Response according Warkentin 2004 (Poisson) ##
DR.Warkentin <- function (TD50=45, gamma50=1.5, aa=1, diffdvh=NULL, dose=NULL) {
  # check the single choice between dvh matrix or dose series
  if (!is.null(dose) && !is.null(diffdvh)) stop("Select either a DVH or a point dose to calculate NTCP")
  # check the single choice between dvh matrix or dose series  
  if (is.vector(dose) && is.null(diffdvh)) p <- 0.5^(exp(2*gamma50/log(2)*(1-dose/TD50)))
  if (is.null(dose) && is.matrix(diffdvh)) p <- 0.5^(exp(2*gamma50/log(2)*(1-EUD(dvh.matrix=diffdvh, a=aa)/TD50)))
  return(p) 
}

## Dose/Response according Bentzen 1997 (log Poisson) ##
DR.Bentzen <- function (TD50=45, gamma50=1.5, aa=1, diffdvh=NULL, dose=NULL) {
  # check the single choice between dvh matrix or dose series
  if (!is.null(dose) && !is.null(diffdvh)) stop("Select either a DVH or a point dose to calculate NTCP")
  # check the single choice between dvh matrix or dose series  
  if (is.vector(dose) && is.null(diffdvh)) p <- 0.5^(TD50/dose)^(2*gamma50/log(2))
  if (is.null(dose) && is.matrix(diffdvh)) p <- 0.5^(TD50/EUD(dvh.matrix=diffdvh, a=aa))^(2*gamma50/log(2))
  return(p) 
}

## definition of outcome based by dose
outcome.def <- function (dvh.matrix=NULL, dose=NULL, TD50=45, gamma50=2, a=1, 
                         model=c("Lyman", "Niemierko", "Okunieff", "Munro")) {
  # checks for not empty both dose and dvh.matrix 
  if (is.matrix(dvh.matrix) && is.vector(dose)) 
    stop("You can choice to simulate the outcome EITHER for a dvh based patients series ", "\n", 
         "OR for a single doses vector patients series")
  # check matrix input in dvh.matrix
  if (!is.matrix(dvh.matrix) && is.null(dose))    
    stop("dvh.matrix MUST be a matrix containing relative differential DVHs")
  # check vector input in dose
  if (!is.vector(dose) && is.null(dvh.matrix))
    stop("dose MUST be a vector containing patients doses")
  # generating outcome simulation for a dvh matrix
  if (is.matrix(dvh.matrix)) { 
    if ((dvh.matrix[1,2]==1) && (dvh.matrix[2,2]>0)) 
      warning("dvh.matrix seems to be a cumulative DVH, ", "\n", 
              "  If so, and subsequent fit.NTCP converges,", "\n",
              "  results are not consisntent.", "\n",
              "  If dvh.matrix is a cumulative DVH matrix", "\n",
              "  convert into a differential matrix ", "\n" ,
              "  by using diff.dvh function.")
    dvh.size <- dim(dvh.matrix)                                       # size of NTCP series
    temp.outcome <- seq(1, dvh.size[2] - 1,1)                         # sequence of temporary outcomes
    eqdoses <- EUD(dvh.matrix=dvh.matrix, a=a)                        # sequence of EUD
    model<-match.arg(model)
    if (model=="Lyman") {
      NTCP.Series <- DR.Lyman(TD50=TD50, gamma50=gamma50, aa=a, diffdvh=dvh.matrix)
      for (m in 1:dvh.size[2]) {
        ifelse(runif(1,0,1) < NTCP.Series[m],temp.outcome[m] <- 1, temp.outcome[m] <-0) # outcome definition for Lyman model 
      }
    }
    if (model=="Niemierko") {
      NTCP.Series <- DR.Niemierko(TD50=TD50, gamma50=gamma50, aa=a, diffdvh=dvh.matrix)
      for (m in 1:dvh.size[2]) {
        ifelse(runif(1,0,1) < NTCP.Series[m],temp.outcome[m] <- 1, temp.outcome[m] <-0) # outcome definition for Niemierko 
      }
    }
    outcome <- matrix(nrow=length(temp.outcome), ncol=3)              # create the outcome matrix 
    outcome[,1] <- eqdoses                                            # 1st column: EUDs
    outcome[,2] <- NTCP.Series                                        # 2nd column: NTCP
    outcome[,3] <- temp.outcome                                       # 3rd column: Outcome
    colnames(outcome) <- c("Dose", "NTCP", "Outcome")
    return(outcome)
  }
  if (is.vector(dose)) {
    l <- length(x=dose)
    temp.outcome <- seq(from=1, to=l, 1)
    model<-match.arg(model)
    if (model=="Lyman") {
      NTCP.Series <- DR.Lyman(TD50=TD50, gamma50=gamma50, dose=dose)
      for (m in 1:l) {
        ifelse(runif(1,0,1) < NTCP.Series[m],temp.outcome[m] <- 1, temp.outcome[m] <-0) # outcome definition for Lyman model 
      }
    }
    if (model=="Niemierko") {
      NTCP.Series <- DR.Niemierko(TD50=TD50, gamma50=gamma50, dose=dose)
      for (m in 1:dvh.size[2]) {
        ifelse(runif(1,0,1) < NTCP.Series[m],temp.outcome[m] <- 1, temp.outcome[m] <-0) # outcome definition for Niemierko 
      }
    }
    outcome <- matrix(nrow=length(temp.outcome), ncol=3)              # create the outcome matrix 
    outcome[,1] <- dose                                               # 1st column: vector dose
    outcome[,2] <- NTCP.Series                                        # 2nd column: NTCP
    outcome[,3] <- temp.outcome                                       # 3rd column: Outcome
    colnames(outcome) <- c("Dose", "NTCP", "Outcome")
    return(outcome)
  }
}

## create NTCP curve
NTCP.curve <- function(TD50=45,gamma50=1.5,model=c("Lyman","Niemierko"),dmin=0,dmax=70,add=TRUE) {
  if (model=="Lyman") {
    curve(expr=pnorm(q=((x - TD50)*gamma50*sqrt(2*pi)/TD50)),from=dmin,to=dmax,ylim=c(0,1),
          add=add,xlab="Dose [Gy]",ylab="NTCP [*100%]")
  }
  if (model=="Niemierko") {
    curve(expr=1/(1+(TD50/EUD(dvh.matrix=diffdvh, a=aa))^(4*gamma50)),from=dmin,to=dmax,ylim=c(0,1),
          add=add,xlab="Dose [Gy]",ylab="NTCP [*100%]")
  }  
}



# Log-Likelihood of Lyman function
probit.nLL.Lyman <- function (z,dvh.matrix, outcome) {
  TD50 <- z[1]
  gamma50 <- z[2]
  aa <- z[3]
  p <- pnorm(gamma50*sqrt(2*pi)*(EUD(dvh.matrix, a=aa)/TD50-1))
  return(-sum(outcome*log(p)+(1-outcome)*log(1-p)))
}



# Log-Likelihood of Niemierko function
logit.nLL.Niemierko <- function (z,dvh.matrix, outcome) {
  TD50 <- z[1]
  gamma50 <- z[2]
  aa <- z[3]
  p <- 1/(1+(TD50/EUD(dvh.matrix, a=aa))^(4*gamma50))
  return(-sum(outcome*log(p)+(1-outcome)*log(1-p)))
}

# function for fitting DVH data
fit.NTCP <- function(model=c("Lyman","Niemierko"),dvh.matrix, outcome, 
                     calcCI=TRUE, verbose=FALSE, C.I.width=.95, n.boot=1000,
                     min.TD50=0, max.TD50=150, min.gamma50=0, max.gamma50=15,
                     min.a=-100, max.a=100, C.I.type=c("ProfLik", "Boot")) {
  model<-match.arg(model)
  C.I.type<-match.arg(C.I.type)
  ca1 <<- qchisq(C.I.width,1)/2  # set the bounds for C.I. calculation
  if (verbose==TRUE) message("Fitting ", model, " model")
  if (!is.vector(outcome)) stop("'outcome' MUST be a vector")
  if (verbose==TRUE) message("Outcome series contains ", length(outcome), " observations", "\n")
  if (model=="Lyman") {
    out <- nlminb(c(40,2,2), probit.nLL.Lyman,lower=c(10,0.01,-100),upper=c(180,15,100),
                  dvh.matrix=dvh.matrix,outcome=outcome)
  }
  if (model=="Niemierko") {
    out <- nlminb(c(40,2,2), logit.nLL.Niemierko,lower=c(10,0.01,-100),upper=c(180,15,100),
                  dvh.matrix=dvh.matrix,outcome=outcome)
  }
  # MLE calculation
  MLE <<- -out$objective
  # modelAIC<- 6 -2*out$objective # Akaike information criterion
  modelAIC <- -2*out$objective + 6 +(24 / (ncol(dvh.matrix)-5)) # Akaike information criterion with num cases correction
  # fitted parameters
  D50 <<- out$par[1]
  g50 <<- out$par[2]
  aaa <<- out$par[3]
  iter <- out$iterations
  if (verbose==TRUE) {
    message(model, " model successfully converged")
    message("      TD50 = ", format(D50, digits=3, nsmall=3, justify="right", width=8))
    message("   gamma50 = ", format(g50, digits=3, nsmall=3, justify="right", width=8))
    message("         a = ", format(aaa, digits=3, nsmall=3, justify="right", width=8))
    message("Iterations = ", format(iter, justify="right", width=8), "\n")
  }
  euddose<-EUD(dvh.matrix=dvh.matrix,a=aaa)     ## dose of patients
  # calculation of predicted probabilities
  # outcome MUST be a vector
  predicted<-matrix(nrow=length(outcome),ncol=2) 
  predicted[,2]<-outcome
  
  if (model=="Lyman") {
    predicted[,1]<-DR.Lyman(TD50=D50,gamma50=g50,aa=aaa,diffdvh=dvh.matrix)    
  }
  if (model=="Niemierko") {
    predicted[,1]<-DR.Niemierko(TD50=D50,gamma50=g50,aa=aaa,diffdvh=dvh.matrix)
  }

  # calculation of deviance http://http://stats.stackexchange.com/questions/6581/what-is-deviance-specifically-in-cart-rpart
  res<-outcome - predicted[,1]  # first calculates the residuals
  moddev<-t(res)%*%res          # calculates the deviance 
  # exclusion of values '1' that give NaN in LL function
  predicted[predicted[,1]!=1,]
  # initialize variables for C.I. calculation
  bmatrix<-NULL
  if (calcCI==TRUE) {
    if (C.I.type=="ProfLik") {
      # Confidence interval calculation with old procedure
      # outcome values rounded at 2nd decimal place
      d50L <- CI(MLE=MLE,arg="d50L",dvh.matrix=dvh.matrix,outcome=outcome,model=model,
                 min.TD50=min.TD50, max.TD50=max.TD50, min.gamma50=min.gamma50,
                 max.gamma50=max.gamma50, min.a=min.a, max.a=max.a)
      d50U <- CI(MLE=MLE,arg="d50U",dvh.matrix=dvh.matrix,outcome=outcome,model=model,
                 min.TD50=min.TD50, max.TD50=max.TD50, min.gamma50=min.gamma50,
                 max.gamma50=max.gamma50, min.a=min.a, max.a=max.a)
      g50L <- CI(MLE=MLE,arg="g50L",dvh.matrix=dvh.matrix,outcome=outcome,model=model,
                 min.TD50=min.TD50, max.TD50=max.TD50, min.gamma50=min.gamma50,
                 max.gamma50=max.gamma50, min.a=min.a, max.a=max.a)
      g50U <- CI(MLE=MLE,arg="g50U",dvh.matrix=dvh.matrix,outcome=outcome,model=model,
                 min.TD50=min.TD50, max.TD50=max.TD50, min.gamma50=min.gamma50,
                 max.gamma50=max.gamma50, min.a=min.a, max.a=max.a)
      aL   <- CI(MLE=MLE,arg="aL",dvh.matrix=dvh.matrix,outcome=outcome,model=model,
                 min.TD50=min.TD50, max.TD50=max.TD50, min.gamma50=min.gamma50,
                 max.gamma50=max.gamma50, min.a=min.a, max.a=max.a)
      aU   <- CI(MLE=MLE,arg="aU",dvh.matrix=dvh.matrix,outcome=outcome,model=model,
                 min.TD50=min.TD50, max.TD50=max.TD50, min.gamma50=min.gamma50,
                 max.gamma50=max.gamma50, min.a=min.a, max.a=max.a)
    }
    if (C.I.type=="Boot") {
      # Confidence interval calculation using bootstrapping method
      temp.boot<-bootstrap.CI(model=model, TD50=D50, gamma50=g50, a=aaa, n=n.boot, 
                              dvh.matrix=dvh.matrix, C.I.width=C.I.width)
      d50L<-temp.boot$TD50L
      d50U<-temp.boot$TD50U
      g50L<-temp.boot$gamma50L
      g50U<-temp.boot$gamma50U
      aL  <-temp.boot$aL
      aU  <-temp.boot$aU
      # calculate max interation
      iter<-iter + temp.boot$TD50L[4] + temp.boot$TD50U[4] + 
        temp.boot$gamma50L[4] + temp.boot$gamma50U[4] +
        temp.boot$aL[4] + temp.boot$aU[4]
      # insert the bootstrap matrix
      bmatrix <- temp.boot$bootstrap.matrix
    }
  }
  
  # set NA values in CI if CI are not calculated
  if (calcCI==FALSE) {
    d50L <- c(NA,NA,NA,NA)
    d50U <- c(NA,NA,NA,NA)
    g50L <- c(NA,NA,NA,NA)
    g50U <- c(NA,NA,NA,NA)
    aL <- c(NA,NA,NA,NA)
    aU <- c(NA,NA,NA,NA)
  }
  
  output <- list("TD50"=D50,"TD50L"=d50L,"TD50U"=d50U,
                 "gamma50"=g50,"gamma50L"=g50L,"gamma50U"=g50U,
                 "a"=aaa,"aL"=aL,"aU"=aU,"model"=model,"Iterations"=iter,
                 "MLE"=MLE,"aic"=modelAIC,"outcome"=outcome,"dose"=euddose,
                 "deviance"=moddev,"predicted"=predicted[,1],
                 "DVH.matrix"=dvh.matrix,"bootstrap.matrix"=bmatrix)
  class(output)<-"NTCP"
  # remove global objects form .GlobalEnv
  if (!is.null(highD50)) rm(highD50, lowD50, highg50, lowg50, higha, lowa, envir=.GlobalEnv)
  rm(D50, aaa, g50, MLE, envir=.GlobalEnv)
  return(output)
}

# light version of fit.NTCP it's a function for fitting
# DVH data without C.I., used for bootstrap C.I. calculation
fit.NTCP.CI <- function(model=c("Lyman","Niemierko"),dvh.matrix, outcome, 
                     min.TD50=0, max.TD50=150, min.gamma50=0, max.gamma50=15,
                     min.a=-100, max.a=100) {
  model<-match.arg(model)
  if (model=="Lyman") {
    out <- nlminb(c(40,2,2), probit.nLL.Lyman, lower=c(min.TD50,min.gamma50,min.a),
                  upper=c(max.TD50,max.gamma50,max.a), dvh.matrix=dvh.matrix, outcome=outcome)
  }
  if (model=="Niemierko") {
    out <- nlminb(c(40,2,2), logit.nLL.Niemierko, lower=c(min.TD50,min.gamma50,min.a),
                  upper=c(max.TD50,max.gamma50,max.a), dvh.matrix=dvh.matrix, outcome=outcome)
  }
  # generate output
  output <- list("TD50"=out$par[1], "gamma50"=out$par[2], "a"=out$par[3],
                 "Iterations"=out$iterations, "MLE"=-out$objective)
  return(output)
}

# function for printing basic outcome of fitted functions
print.NTCP<-function(fittedmodel) {
  cat("                           -- NTCP model estimated parameters --","\n\n")
  cat("Model type:",fittedmodel$model,"\n\n")
  cat("Parameter:             Value:                         95% C.I.","\n\n")
  cat("TD50:   ",format(x=c(fittedmodel$TD50,fittedmodel$TD50L[1],fittedmodel$TD50U[1]),
                     justify="right",digits=3,nsmall=3,width=20),"\n")
  cat("g50:    ",format(x=c(fittedmodel$gamma50,fittedmodel$gamma50L[1],fittedmodel$gamma50U[1]),
                     justify="right",digits=3,nsmall=3,width=20),"\n")
  cat("a:      ",format(x=c(fittedmodel$a,fittedmodel$aL[1],fittedmodel$aU[1]),
                     justify="right",digits=3,nsmall=3,width=20),"\n")
  cat("n (1/a):",format(x=c(1/fittedmodel$a,1/fittedmodel$aL[1],1/fittedmodel$aU[1]),
                     justify="right",digits=3,nsmall=3,width=20),"\n")
}

# function that shows model detailed features and alternate parameters
# found during profile likelihood procedures for confidence intervals calculation
summary.NTCP<-function(fittedmodel) {
  cat("                           -- NTCP model estimated parameters --","\n")
  cat("                           -- with alternate parameter values --","\n")
  cat("                           -- corresponding to C.I. limits    --","\n\n\n")
  cat("Model type:",fittedmodel$model,"\n")
  cat("MLE:",fittedmodel$MLE,"\n")
  cat("AIC:",fittedmodel$aic,"\n")
  cat("Deviance:",fittedmodel$deviance,"\n")
  cat("Model iterations number:",fittedmodel$Iterations,"\n\n\n")

  cat("Parameter:             Value:                         95% C.I.","\n\n")
  cat("TD50:   ",format(x=c(fittedmodel$TD50,fittedmodel$TD50L[1],fittedmodel$TD50U[1]),
                        justify="right",digits=3,nsmall=3,width=20),"\n")
  cat("Corresponding g50:           ",format(x=c(fittedmodel$TD50L[2],fittedmodel$TD50U[2]),
                        justify="right",digits=3,nsmall=3,width=20),"\n")
  cat("Corresponding a:             ",format(x=c(fittedmodel$TD50L[3],fittedmodel$TD50U[3]),
                        justify="right",digits=3,nsmall=3,width=20),"\n")
  cat("corresponding n=(1/a):       ",format(x=c(1/fittedmodel$TD50L[3],1/fittedmodel$TD50U[3]),
                        justify="right",digits=3,nsmall=3,width=20),"\n")
  it<-fittedmodel$TD50L[4]+fittedmodel$TD50U[4]
  cat("C.I. Iterations number: ",it,"\n\n")

  cat("g50:    ",format(x=c(fittedmodel$gamma50,fittedmodel$gamma50L[1],fittedmodel$gamma50U[1]),
                        justify="right",digits=3,nsmall=3,width=20),"\n")
  cat("Corresponding TD50:          ",format(x=c(fittedmodel$gamma50L[2],fittedmodel$gamma50U[2]),
                                             justify="right",digits=3,nsmall=3,width=20),"\n")
  cat("Corresponding a:             ",format(x=c(fittedmodel$gamma50L[3],fittedmodel$gamma50U[3]),
                                             justify="right",digits=3,nsmall=3,width=20),"\n")
  cat("corresponding n=(1/a):       ",format(x=c(1/fittedmodel$gamma50L[3],1/fittedmodel$gamma50U[3]),
                                             justify="right",digits=3,nsmall=3,width=20),"\n")
  it<-fittedmodel$gamma50L[4]+fittedmodel$gamma50U[4]
  cat("C.I. Iterations number: ",it,"\n\n")
  
  cat("a:      ",format(x=c(fittedmodel$a,fittedmodel$aL[1],fittedmodel$aU[1]),
                        justify="right",digits=3,nsmall=3,width=20),"\n")
  cat("n=(1/a):",format(x=c(1/fittedmodel$a,1/fittedmodel$aL[1],1/fittedmodel$aU[1]),
                        justify="right",digits=3,nsmall=3,width=20),"\n")
  cat("Corresponding TD50:          ",format(x=c(fittedmodel$aL[2],fittedmodel$aU[2]),
                                             justify="right",digits=3,nsmall=3,width=20),"\n")
  cat("Corresponding gamma50:       ",format(x=c(fittedmodel$aL[3],fittedmodel$aU[3]),
                                             justify="right",digits=3,nsmall=3,width=20),"\n")

  it<-fittedmodel$aL[4]+fittedmodel$aU[4]
  cat("C.I. Iterations number: ",it,"\n\n")  
}

# function for predicting outcome in a given NTCP model
# with a specific level of dose
outcome.NTCP <- function(fittedmodel, dose) {
  x <- dose
  if (fittedmodel$model=="Lyman") {
    # y values of min C.I. curve
    ymin <- min(pnorm(q=((x - fittedmodel$TD50L[1])*fittedmodel$TD50L[2]*sqrt(2*pi))/fittedmodel$TD50L[1]),
                 pnorm(q=((x - fittedmodel$TD50U[1])*fittedmodel$TD50U[2]*sqrt(2*pi))/fittedmodel$TD50U[1]),
                 pnorm(q=((x - fittedmodel$gamma50L[2])*fittedmodel$gamma50L[1]*sqrt(2*pi))/fittedmodel$gamma50L[2]),
                 pnorm(q=((x - fittedmodel$gamma50U[2])*fittedmodel$gamma50U[1]*sqrt(2*pi))/fittedmodel$gamma50U[2]),
                 pnorm(q=((x - fittedmodel$aL[2])*fittedmodel$aL[3]*sqrt(2*pi))/fittedmodel$aL[2]),
                 pnorm(q=((x - fittedmodel$aU[2])*fittedmodel$aU[3]*sqrt(2*pi))/fittedmodel$aU[2]),
                 pnorm(q=((x - fittedmodel$TD50)*fittedmodel$gamma50*sqrt(2*pi))/fittedmodel$TD50))
    # y values of max C.I. curve
    ymax <- max(pnorm(q=((x - fittedmodel$TD50L[1])*fittedmodel$TD50L[2]*sqrt(2*pi))/fittedmodel$TD50L[1]),
                 pnorm(q=((x - fittedmodel$TD50U[1])*fittedmodel$TD50U[2]*sqrt(2*pi))/fittedmodel$TD50U[1]),
                 pnorm(q=((x - fittedmodel$gamma50L[2])*fittedmodel$gamma50L[1]*sqrt(2*pi))/fittedmodel$gamma50L[2]),
                 pnorm(q=((x - fittedmodel$gamma50U[2])*fittedmodel$gamma50U[1]*sqrt(2*pi))/fittedmodel$gamma50U[2]),
                 pnorm(q=((x - fittedmodel$aL[2])*fittedmodel$aL[3]*sqrt(2*pi))/fittedmodel$aL[2]),
                 pnorm(q=((x - fittedmodel$aU[2])*fittedmodel$aU[3]*sqrt(2*pi))/fittedmodel$aU[2]),
                 pnorm(q=((x - fittedmodel$TD50)*fittedmodel$gamma50*sqrt(2*pi))/fittedmodel$TD50))
    # predicted outcome with optimal parameters
    y <- pnorm(q=((x - fittedmodel$TD50)*fittedmodel$gamma50*sqrt(2*pi)/fittedmodel$TD50))
  }
  if (fittedmodel$model=="Niemierko") {
    # y values of min C.I. curve
    ymin <- min(1/(1+(fittedmodel$TD50L[1]/x)^(4*fittedmodel$TD50L[2])),
                 1/(1+(fittedmodel$TD50U[1]/x)^(4*fittedmodel$TD50U[2])),
                 1/(1+(fittedmodel$gamma50L[2]/x)^(4*fittedmodel$gamma50L[1])),
                 1/(1+(fittedmodel$gamma50U[2]/x)^(4*fittedmodel$gamma50U[1])),
                 1/(1+(fittedmodel$aL[2]/x)^(4*fittedmodel$aL[3])),
                 1/(1+(fittedmodel$aU[2]/x)^(4*fittedmodel$aU[3])),
                 1/(1+(fittedmodel$TD50/x)^(4*fittedmodel$gamma50)))
    # y values of max C.I. curve
    ymax <- max(1/(1+(fittedmodel$TD50L[1]/x)^(4*fittedmodel$TD50L[2])),
                 1/(1+(fittedmodel$TD50U[1]/x)^(4*fittedmodel$TD50U[2])),
                 1/(1+(fittedmodel$gamma50L[2]/x)^(4*fittedmodel$gamma50L[1])),
                 1/(1+(fittedmodel$gamma50U[2]/x)^(4*fittedmodel$gamma50U[1])),
                 1/(1+(fittedmodel$aL[2]/x)^(4*fittedmodel$aL[3])),
                 1/(1+(fittedmodel$aU[2]/x)^(4*fittedmodel$aU[3])),
                 1/(1+(fittedmodel$TD50/x)^(4*fittedmodel$gamma50)))
    #predicted outcome with optimal parameters
    y <- 1/(1+(fittedmodel$TD50/x)^(4*fittedmodel$gamma50))
  }  
  cat("Dose:          Probability:                     95% C.I.","\n\n")
  cat(x, format(x=c(y, ymin, ymax),
                        justify="right",digits=3,nsmall=3,width=20),"\n")
}

# function for plotting NTCP curves:
# cases=TRUE: shows all cases on horizontal axes NTCP=0 and NTCP=1
# C.I.=TRUE: shows the confidence interval of the NTCP curve
# quantiles=FALSE: creates quantiles of cases for showing dots along the NTCP curve
# quantiles.no: quantiles number to be shown
# C.I.quantiles: shows the C.I. for NTCP of given quantiles
# C.I.width: width of confidence interval for quatiles
plot.NTCP <- function(fittedmodel,cases=TRUE,C.I.=TRUE,quantiles=FALSE,xmin=0,xmax=0,
                      quantiles.no=4,C.I.quantiles=FALSE,C.I.width=0.95,bounds=FALSE,...) {  
  # plot cases in the graph
  if (cases==TRUE){
    if (xmax==0) {
      plot(x=fittedmodel$dose,y=fittedmodel$outcome,xlim=c(xmin,max(fittedmodel$dose)+5),xlab="Dose [Gy]",ylab="NTCP [*100%]")
    } 
    else {
      plot(x=fittedmodel$dose,y=fittedmodel$outcome,xlim=c(xmin,xmax),xlab="Dose [Gy]",ylab="NTCP [*100%]")  
    }
  }
  # plot Lyman model
  if (fittedmodel$model=="Lyman") {
    if (xmax==0) {
      curve(expr=pnorm(q=((x - fittedmodel$TD50)*fittedmodel$gamma50*sqrt(2*pi)/fittedmodel$TD50)),from=xmin,to=max(fittedmodel$dose)+5,ylim=c(0,1),
            xlab="Dose [Gy]",ylab="NTCP [*100%]",add=cases,lwd=2)
    } else {
      curve(expr=pnorm(q=((x - fittedmodel$TD50)*fittedmodel$gamma50*sqrt(2*pi)/fittedmodel$TD50)),from=xmin,to=xmax,ylim=c(0,1),
            xlab="Dose [Gy]",ylab="NTCP [*100%]",add=cases,lwd=2)  
    }
  }
  # plot Niemierko model
  if (fittedmodel$model=="Niemierko") {
    if (xmax==0) {
      curve(expr=1/(1+(fittedmodel$TD50/x)^(4*fittedmodel$gamma50)),from=xmin,to=max(fittedmodel$dose)+5,ylim=c(0,1),
            xlab="Dose [Gy]",ylab="NTCP [*100%]",add=cases,lwd=2)
    } else {
      curve(expr=1/(1+(fittedmodel$TD50/x)^(4*fittedmodel$gamma50)),from=xmin,to=xmax,ylim=c(0,1),
            xlab="Dose [Gy]",ylab="NTCP [*100%]",add=cases,lwd=2)
    }
  }
  # plot confidence intervals of the models using the extreme values in the parameters matrix
  if (C.I.==TRUE){
    if (xmax==0) {
      x<-seq(0,max(fittedmodel$dose)+5,.25)
    } else {
      x<-seq(0,xmax,.25)
    }
    if (fittedmodel$model=="Lyman") {
      # y values of min C.I. curve
      ymin <- pmin(pnorm(q=((x - fittedmodel$TD50L[1])*fittedmodel$TD50L[2]*sqrt(2*pi))/fittedmodel$TD50L[1]),
                   pnorm(q=((x - fittedmodel$TD50U[1])*fittedmodel$TD50U[2]*sqrt(2*pi))/fittedmodel$TD50U[1]),
                   pnorm(q=((x - fittedmodel$gamma50L[2])*fittedmodel$gamma50L[1]*sqrt(2*pi))/fittedmodel$gamma50L[2]),
                   pnorm(q=((x - fittedmodel$gamma50U[2])*fittedmodel$gamma50U[1]*sqrt(2*pi))/fittedmodel$gamma50U[2]),
                   pnorm(q=((x - fittedmodel$aL[2])*fittedmodel$aL[3]*sqrt(2*pi))/fittedmodel$aL[2]),
                   pnorm(q=((x - fittedmodel$aU[2])*fittedmodel$aU[3]*sqrt(2*pi))/fittedmodel$aU[2]),
                   pnorm(q=((x - fittedmodel$TD50)*fittedmodel$gamma50*sqrt(2*pi))/fittedmodel$TD50))
      # y values of max C.I. curve
      ymax <- pmax(pnorm(q=((x - fittedmodel$TD50L[1])*fittedmodel$TD50L[2]*sqrt(2*pi))/fittedmodel$TD50L[1]),
                   pnorm(q=((x - fittedmodel$TD50U[1])*fittedmodel$TD50U[2]*sqrt(2*pi))/fittedmodel$TD50U[1]),
                   pnorm(q=((x - fittedmodel$gamma50L[2])*fittedmodel$gamma50L[1]*sqrt(2*pi))/fittedmodel$gamma50L[2]),
                   pnorm(q=((x - fittedmodel$gamma50U[2])*fittedmodel$gamma50U[1]*sqrt(2*pi))/fittedmodel$gamma50U[2]),
                   pnorm(q=((x - fittedmodel$aL[2])*fittedmodel$aL[3]*sqrt(2*pi))/fittedmodel$aL[2]),
                   pnorm(q=((x - fittedmodel$aU[2])*fittedmodel$aU[3]*sqrt(2*pi))/fittedmodel$aU[2]),
                   pnorm(q=((x - fittedmodel$TD50)*fittedmodel$gamma50*sqrt(2*pi))/fittedmodel$TD50))
    }
    if (fittedmodel$model=="Niemierko") {
      # y values of min C.I. curve
      ymin <- pmin(1/(1+(fittedmodel$TD50L[1]/x)^(4*fittedmodel$TD50L[2])),
                   1/(1+(fittedmodel$TD50U[1]/x)^(4*fittedmodel$TD50U[2])),
                   1/(1+(fittedmodel$gamma50L[2]/x)^(4*fittedmodel$gamma50L[1])),
                   1/(1+(fittedmodel$gamma50U[2]/x)^(4*fittedmodel$gamma50U[1])),
                   1/(1+(fittedmodel$aL[2]/x)^(4*fittedmodel$aL[3])),
                   1/(1+(fittedmodel$aU[2]/x)^(4*fittedmodel$aU[3])),
                   1/(1+(fittedmodel$TD50/x)^(4*fittedmodel$gamma50)))
      # y values of max C.I. curve
      ymax <- pmax(1/(1+(fittedmodel$TD50L[1]/x)^(4*fittedmodel$TD50L[2])),
                   1/(1+(fittedmodel$TD50U[1]/x)^(4*fittedmodel$TD50U[2])),
                   1/(1+(fittedmodel$gamma50L[2]/x)^(4*fittedmodel$gamma50L[1])),
                   1/(1+(fittedmodel$gamma50U[2]/x)^(4*fittedmodel$gamma50U[1])),
                   1/(1+(fittedmodel$aL[2]/x)^(4*fittedmodel$aL[3])),
                   1/(1+(fittedmodel$aU[2]/x)^(4*fittedmodel$aU[3])),
                   1/(1+(fittedmodel$TD50/x)^(4*fittedmodel$gamma50)))
    }    
    #lines(x,ymin,lty="dotted",lwd=2)
    #lines(x,ymax,lty="dotted",lwd=2)
    lines(x,ymin)
    lines(x,ymax)
  }
  # calculation of quantiles of cases in the plot for plotting quantiles
  if (quantiles==TRUE) {
    fq<-1/quantiles.no                    # fraction of each quantile
    meandoses<-c(1)                       # vector of mean doses
    length(meandoses)<-quantiles.no       # set length of vector of means
    meanoutcome<-c(1)                     # vector of mean outcomes
    length(meanoutcome)<-quantiles.no
    CImatrix<-matrix(nrow=quantiles.no,ncol=2)
    colnames(CImatrix)<-c(
      paste((1-C.I.width)/2,"%",sep=""),
      paste(C.I.width+(1-C.I.width)/2,"%",sep="")
    )
    series<-matrix(nrow=length(fittedmodel$dose), ncol=2)
    series[,1]<-fittedmodel$dose
    series[,2]<-fittedmodel$outcome
    for (m in 1:quantiles.no) {           # loop through quantiles
      doseV<-(1)                          # empty vector of mean doses
      outcomeV<-(1)                       # empty vector of mean outcome
      count<-0                            # counter of elements for each vector of means and results
      for  (n in 1:nrow(series)) {        # loop through all elements in the list        
        if (m==1){                        # 1st quantile
          if (series[n,1][[1]]<=quantile(x=fittedmodel$dose,fq)) {
            count<-count+1
            length(doseV)<-count
            length(outcomeV)<-count
            doseV[count]<-series[n,1][[1]]    # set the element of the dose vector
            outcomeV[count]<-series[n,2][[1]] # set the element of the outcome vector            
          }
        }
        if ((series[n,1][[1]]<=quantile(x=fittedmodel$dose,fq*m))&&(series[n,1][[1]]>quantile(x=fittedmodel$dose,fq*(m-1)))) {
          count<-count+1
          length(doseV)<-count
          length(outcomeV)<-count
          doseV[count]<-series[n,1][[1]]    # set the element of the dose vector
          outcomeV[count]<-series[n,2][[1]] # set the element of the outcome vector
        }          
      }
      meandoses[m]<-mean(doseV)        # set the element of vector of mean doses
      meanoutcome[m]<-mean(outcomeV)   # set the element of vector of mean outcomes
      bstrap <- c(1)                   # vector of bootstrap series for outcome C.I. calculation
      for (i in 1:1000) {
        bsample <- sample(outcomeV,length(outcomeV),replace=T) # sample generation
        bestimate <- mean(bsample)                             # bootstrapped sample mean
        bstrap <- c(bstrap,bestimate)
      }
      CImatrix[m,1]<-quantile(bstrap,(1-C.I.width)/2)          # lower bound of C.I. for mean
      CImatrix[m,2]<-quantile(bstrap,1-(1-C.I.width)/2)        # higher bound of C.I. for mean
    }
    for (p in 1:quantiles.no) {        # plot patients quantiles
      points(x=meandoses,y=meanoutcome,pch=15,cex=1.5)
    }
    if (C.I.quantiles==TRUE) {         # plot outcome C.I. for quantiles
      for (p in 1:quantiles.no) {
        lines(x=c(meandoses[p],meandoses[p]),y=c(CImatrix[p,1],CImatrix[p,2]))      # lines of C.I.
        lines(x=c(meandoses[p]-1,meandoses[p]+1),y=c(CImatrix[p,1],CImatrix[p,1]))  # lower bracket of C.I.
        lines(x=c(meandoses[p]-1,meandoses[p]+1),y=c(CImatrix[p,2],CImatrix[p,2]))  # higher bracket of C.I.
      }
      list("quantiles number"=quantiles.no,"quantiles mean doses"=meandoses,
           "quantiles mean outcomes"=meanoutcome,"quantiles confidence intervals"=CImatrix)
    }
  }  
}

# function for plotting the 3D density kernel of NTCP values of a single differential DVH
plot3D.NTCP <- function(dvh=NULL, model=NULL) {
  label<-c("\n")
  err.state<-0  # no error in getting varibles
  if (sum(dvh[,2])!=1) {
    err.state<-1
    label<-c(label,"dvh doesn't seem to be a relative differential dvh","\n")
  }
  if ((is.null(dvh)) || (is.matrix(dvh)==FALSE)) {
    err.state<-1
    label<-c(label,"You MUST provide a valid differential dvh for calculating 3D plot","\n")
  }
  if (ncol(dvh)>2) {
    err.state<-1
    label<-c(label,"dvh must be a matrix with two columns, 1st as dose bins, 2nd as volume bins", "\n")
  }
  if (is.null(model)) {
    err.state<-1
    label<-c(label, "You MUST provide a valid NTCP model for calculating 3D plot", "\n")
  }
  if (err.state==1) stop(label)
  # calculate the column of EUDs and NTCPs
  doses<-c()
  NTCPs<-c()
  if (model$model=="Lyman")
    for (n in 1:nrow(model$bootstrap.matrix)) {
      doses<-c(doses, EUD(dvh.matrix=dvh, a=model$bootstrap.matrix$a[n]))
      NTCPs<-c(NTCPs, DR.Lyman(TD50=model$bootstrap.matrix$TD50[n], 
                                 gamma50=model$bootstrap.matrix$gamma50[n], 
                                 dose=doses[n]))
    }
  if (model$model=="Niemierko") 
    for (n in 1:nrow(model$bootstrap.matrix)) {
      doses<-c(doses, EUD(dvh.matrix=dvh, a=model$bootstrap.matrix$a[n]))
      NTCPs<-c(NTCPs, DR.Niemierko(TD50=model$bootstrap.matrix$TD50[n], 
                                     gamma50=model$bootstrap.matrix$gamma50[n], 
                                     dose=doses[n]))
    }
  # creates 3D density kernel
  require(MASS)
  # density estimate
  k<-kde2d(x=doses, y=NTCPs, n=100)
  # plot perspective graph
  persp(k, box=TRUE, xlab="EUD (Gy)", ylab="NTCP", r=30, theta=135, phi=45, ticktype="detailed")  
  # calculate final values for EUD and NTCP for optimized parameters
  final.EUD<-EUD(dvh.matrix=dvh, a=model$a)
  if (model$model=="Lyman") final.NTCP<-DR.Lyman(TD50=model$TD50, gamma50=model$gamma50, aa=model$a, diffdvh=dvh)
  if (model$model=="Niemierko") final.NTCP<-DR.Niemierko(TD50=model$TD50, gamma50=model$gamma50, aa=model$a, diffdvh=dvh)
  points(final.EUD, y=final.NTCP, pch=20, col="red")
  plot(doses, NTCPs, pch="+")    
  points(final.EUD, y=final.NTCP, pch=20, col="red")
  abline(v=EUD(dvh.matrix=dvh, a=model$a), col="red", lty=2)
  abline(h=DR.Lyman(TD50=model$TD50, gamma50=model$gamma50, aa=model$a, diffdvh=dvh), col="red", lty=2)
  # calculates the filled contour of the kerneldensity  matrix
  # first calculate the volume of the 3D surface and normalizes to 1
  deltaX<-k$x[2]-k$x[1]
  deltaY<-k$y[2]-k$y[1]
  z.mat<-deltaX*deltaY*k$z       # matrix of volumes of bins of dose/NTCP
  total.abs.volume<-sum(z.mat)   # absolute total volume under the surface (3D integral)
  filled.contour(k)
  return(list("kernel.matrix"=k,"dose"=doses,"NTCP"=NTCPs))
}


# Log-Likelihood for Lyman TD50
probit.nLL.Lyman.TD50 <- function(z,C.I.TD50,dvh.matrix,outcome) { 
  gamma50 <- z[1]
  a <- z[2]
  pp <- pnorm(gamma50*sqrt(2*pi)*(EUD(dvh.matrix, a)/C.I.TD50-1))
  p <- c(1:length(pp))
  m<-length(pp)
  for (n in 1:m) {
    if (pp[n]==1) p[n]<-0 else p[n]<-outcome[n]*log(pp[n])+(1-outcome[n])*log(1-pp[n])
  }
  return(-sum(p))
}

# Log-Likelihood for Niemierko TD50
logit.nLL.Niemierko.TD50 <- function(z,C.I.TD50,dvh.matrix,outcome) { 
  gamma50 <- z[1]
  a <- z[2]
  pp <- 1/(1+(C.I.TD50/EUD(dvh.matrix, a))^(4*gamma50))
  p <- c(1:length(pp))
  m<-length(pp)
  for (n in 1:m) {
    if (pp[n]==1) p[n]<-0 else p[n]<-outcome[n]*log(pp[n])+(1-outcome[n])*log(1-pp[n])
  }
  return(-sum(p))
}

# Log-Likelihood for Lyman gamma50
probit.nLL.Lyman.gamma50 <- function(z,C.I.gamma50,dvh.matrix,outcome) { 
  TD50 <- z[1]
  a <- z[2]
  pp <- pnorm(C.I.gamma50*sqrt(2*pi)*(EUD(dvh.matrix, a)/TD50-1))
  p <- c(1:length(pp))
  m<-length(pp)
  for (n in 1:m) {
    if (pp[n]==1) p[n]<-0 else p[n]<-outcome[n]*log(pp[n])+(1-outcome[n])*log(1-pp[n])
  }
  return(-sum(p))
}

# Log-Likelihood for Niemierko gamma50
logit.nLL.Niemierko.gamma50 <- function(z,C.I.gamma50,dvh.matrix,outcome) { 
  TD50 <- z[1]
  a <- z[2]
  pp <- 1/(1+(TD50/EUD(dvh.matrix, a))^(4*C.I.gamma50))
  p <- c(1:length(pp))
  m<-length(pp)
  for (n in 1:m) {
    if (pp[n]==1) p[n]<-0 else p[n]<-outcome[n]*log(pp[n])+(1-outcome[n])*log(1-pp[n])
  }
  return(-sum(p))
}

# Log-Likelihood for Lyman a
probit.nLL.Lyman.a <- function(z,C.I.a,dvh.matrix,outcome) {
  TD50 <- z[1]
  gamma50 <- z[2]
  pp <- pnorm(gamma50*sqrt(2*pi)*(EUD(dvh.matrix, C.I.a)/TD50 - 1))
  p <- c(1:length(pp))
  m<-length(pp)
  for (n in 1:m) {
    if (pp[n]==1) p[n]<-0 else p[n]<-outcome[n]*log(pp[n])+(1-outcome[n])*log(1-pp[n])
  }
  return(-sum(p))
}

# Log-Likelihood for Niemierko a
logit.nLL.Niemierko.a <- function(z,C.I.a,dvh.matrix,outcome) {
  TD50 <- z[1]
  gamma50 <- z[2]
  pp <- 1/(1+(TD50/EUD(dvh.matrix, C.I.a))^(4*gamma50))
  p <- c(1:length(pp))
  m<-length(pp)
  for (n in 1:m) {
    if (pp[n]==1) p[n]<-0 else p[n]<-outcome[n]*log(pp[n])+(1-outcome[n])*log(1-pp[n])
  }
  return(-sum(p))
}

# building bounds of 95% Confidence Interval: profile likelihood
CI <- function(MLE, arg, dvh.matrix, outcome, model,
               min.TD50=0, max.TD50=150, min.gamma50=0, 
               max.gamma50=10, min.a=-100, max.a=100, verbose=TRUE) {
  if (arg=="d50L") {
    lab <-"lower TD50"
    par0<-"   Parameter Lower TD50 = "
    par1<-"  Corresponding gamma50 = "
    par2<-"        Corresponding a = "
  }
  if (arg=="d50U") {
    lab <-"upper TD50"
    par0<-"   Parameter Upper TD50 = "
    par1<-"  Corresponding gamma50 = "
    par2<-"        Corresponding a = "
  }
  if (arg=="g50L") {
    lab<- "lower gamma50"
    par0<-"Parameter Lower gamma50 = "
    par1<-"     Corresponding TD50 = "
    par2<-"        Corresponding a = "
  }
  if (arg=="g50U") {
    lab<- "upper gamma50"
    par0<-"Parameter Upper gamma50 = "
    par1<-"     Corresponding TD50 = "
    par2<-"        Corresponding a = "
  }
  if (arg=="aL") {
    lab<- "lower a"
    par0<-"      Parameter Lower a = "
    par1<-"     Corresponding TD50 = "
    par2<-"  Corresponding gamma50 = "
    
  }
  if (arg=="aU") {
    lab<- "upper a"
    par0<-"      Parameter Upper a = "
    par1<-"     Corresponding TD50 = "
    par2<-"  Corresponding gamma50 = "
  }
  if (verbose==TRUE) message("Calculating boundaries of Confidence Intervals: ", lab)
  it<-0 # number of iterations performed
  # set the model type
  if (model=="Lyman") {
    nLL.TD50 <- probit.nLL.Lyman.TD50
    nLL.gamma50 <- probit.nLL.Lyman.gamma50
    nLL.a <- probit.nLL.Lyman.a
  }
  if (model=="Niemierko") {
    nLL.TD50 <- logit.nLL.Niemierko.TD50
    nLL.gamma50 <- logit.nLL.Niemierko.gamma50
    nLL.a <- logit.nLL.Niemierko.a
  }
  # lower bound of TD50 C.I.
  if (arg=="d50L") {
    highD50 <- D50                 # starting point for interval search
    lowD50 <- D50 - 10              # ending point for interval search
    lowout <- nlminb(start=c(g50,aaa),objective=nLL.TD50,lower=c(min.gamma50,min.a), upper=c(max.gamma50,max.a),
                     C.I.TD50=lowD50,dvh.matrix=dvh.matrix,outcome=outcome)
    it <- lowout$iterations
    while((-lowout$objective-MLE+ca1)>0) { # function to be repeated until LL-MLE+Chi2/2>0
      highD50 <- lowD50
      lowD50 <- lowD50 - 10
      lowout <- nlminb(start=c(g50,aaa),objective=nLL.TD50, lower=c(min.gamma50,min.a), upper=c(max.gamma50,max.a),
                       C.I.TD50=lowD50,dvh.matrix=dvh.matrix,outcome=outcome)
      it <- it + lowout$iterations
      if (it > 10000) stop("Failed convergence for C.I.")
    }
    # starting point for bisection
    highD50 <<- highD50
    lowD50 <<- lowD50
    # start bisection search
    newD50<-(highD50+lowD50)/2
    newout <- nlminb(start=c(g50,aaa),objective=nLL.TD50, lower=c(min.gamma50,min.a), upper=c(max.gamma50,max.a),
                     C.I.TD50=newD50,dvh.matrix=dvh.matrix,outcome=outcome)
    it <- it + newout$iterations
    while (abs(-newout$objective-MLE+ca1)>1e-5) {
      if (highD50-lowD50 < 1e-8) break
      # grid search of D50 value      
      if ((-newout$objective-MLE+ca1)>0) {   # moving around the zero
        highD50 <- newD50
      } else lowD50 <- newD50
      newD50<-(highD50+lowD50)/2
      newout <- nlminb(start=c(g50,aaa),objective=nLL.TD50, lower=c(min.gamma50,min.a), upper=c(max.gamma50,max.a),
                       C.I.TD50=newD50,dvh.matrix=dvh.matrix,outcome=outcome)
      it<-it+newout$iterations  
    }
    result <- (highD50+lowD50)/2 # final value of lower bound of D50 C.I. and iterations number
  }
  
  # higher bound of TD50 C.I.
  else if (arg=="d50U") {
    lowD50 <- D50                   # starting point for interval search
    highD50 <- D50 + 10             # ending point for interval search    
    highout <- nlminb(start=c(g50,aaa),objective=nLL.TD50, lower=c(min.gamma50,min.a), upper=c(max.gamma50,max.a),
                      C.I.TD50=highD50,dvh.matrix=dvh.matrix,outcome=outcome) 
    it<-highout$iterations
    while((-highout$objective-MLE+ca1)>0) { # routine to be repeated until LL-MLE+Chi2/2>0
      lowD50 <- highD50
      highD50 <- highD50 + 10
      highout <- nlminb(start=c(g50,aaa),objective=nLL.TD50, lower=c(min.gamma50,min.a), upper=c(max.gamma50,max.a),
                        C.I.TD50=highD50,dvh.matrix=dvh.matrix,outcome=outcome)
      it<-it+highout$iterations
      if (it > 10000) stop("Failed convergence for C.I.")
    }
    # starting point for bisection
    lowD50 <<- lowD50
    highD50 <<- highD50
    # start bisection search
    newD50<-(highD50+lowD50)/2
    newout <- nlminb(start=c(g50,aaa),objective=nLL.TD50, lower=c(min.gamma50,min.a), upper=c(max.gamma50,max.a),
                     C.I.TD50=newD50,dvh.matrix=dvh.matrix,outcome=outcome)
    it<-it+newout$iterations
    while (abs(-newout$objective-MLE+ca1)>1e-5) {
      if (highD50-lowD50 < 1e-8) break
      # grid search of D50 value      
      if ((-newout$objective-MLE+ca1)>0) {   # moving around the zero
        lowD50 <- newD50
      } else highD50 <- newD50
      newD50<-(highD50+lowD50)/2
      newout <- nlminb(start=c(g50,aaa),objective=nLL.TD50, lower=c(min.gamma50,min.a), upper=c(max.gamma50,max.a),
                       C.I.TD50=newD50,dvh.matrix=dvh.matrix,outcome=outcome)
      it<-it+newout$iterations      
    }
    result <- (highD50+lowD50)/2 # final value of higher bound of D50 C.I.    
  }
  
  # lower bound of gamma50 C.I.
  else if (arg=="g50L") {
    highg50 <- g50                 # starting point for interval search
    lowg50 <- g50 - 0.5            # ending point for interval search
    lowout <- nlminb(start=c(D50,aaa),objective=nLL.gamma50, lower=c(min.TD50,min.a), upper=c(max.TD50,max.a),
                     C.I.gamma50=lowg50,dvh.matrix=dvh.matrix,outcome=outcome)
    it<-lowout$iterations
    while((-lowout$objective-MLE+ca1)>0) { # routine to be repeated until LL-MLE+Chi2/2>0
      highg50 <- lowg50
      lowg50 <- lowg50 - 0.5
      lowout <- nlminb(start=c(D50,aaa),objective=nLL.gamma50, lower=c(min.TD50,min.a), upper=c(max.TD50,max.a),
                       C.I.gamma50=lowg50,dvh.matrix=dvh.matrix,outcome=outcome)
      it<-it+lowout$iterations
      if (it > 10000) stop("Failed convergence for C.I.")
    }
    # starting point for bisection
    highg50 <<- highg50
    lowg50 <<- lowg50
    # start bisection search
    newg50 <- (highg50+lowg50)/2
    newout <- nlminb(start=c(D50,aaa),objective=nLL.gamma50, lower=c(min.TD50,min.a), upper=c(max.TD50,max.a),
                     C.I.gamma50=newg50,dvh.matrix=dvh.matrix,outcome=outcome)
    it<-it+newout$iterations
    while (abs(-newout$objective-MLE+ca1)>1e-5) {
      # grid search of D50 value      
      if ((-newout$objective-MLE+ca1)>0) {   # moving around the zero
        highg50 <- newg50
      } else lowg50 <- newg50
      newg50<-(highg50+lowg50)/2
      newout <- nlminb(start=c(D50,aaa),objective=nLL.gamma50, lower=c(min.TD50,min.a), upper=c(max.TD50,max.a),
                       C.I.gamma50=newg50,dvh.matrix=dvh.matrix,outcome=outcome)
      it<-it+newout$iterations      
    }
    result <- (highg50+lowg50)/2 # final value of lower bound of gamma50 C.I.
  }
  
  # higher bound of gamma50 C.I.
  else if (arg=="g50U") {
    lowg50 <- g50                   # starting point for interval search
    highg50 <- g50 + 0.5            # ending point for interval search
    highout <- nlminb(start=c(D50,aaa),objective=nLL.gamma50, lower=c(min.TD50,min.a), upper=c(max.TD50,max.a),
                      C.I.gamma50=highg50,dvh.matrix=dvh.matrix,outcome=outcome)
    it<-highout$iterations
    while((-highout$objective-MLE+ca1)>0) { # routine to be repeated until LL-MLE+Chi2/2>0
      lowg50 <- highg50
      highg50 <- highg50 + 0.5
      highout <- nlminb(start=c(D50,aaa),objective=nLL.gamma50, lower=c(min.TD50,min.a), upper=c(max.TD50,max.a),
                        C.I.gamma50=highg50,dvh.matrix=dvh.matrix,outcome=outcome)
      it<-it+highout$iterations
      if (it > 10000) stop("Failed convergence for C.I.")
    }
    # starting point for bisection
    lowg50 <<- lowg50
    highg50 <<- highg50
    # start bisection search
    newg50<-(highg50+lowg50)/2
    newout <- nlminb(start=c(D50,aaa),objective=nLL.gamma50, lower=c(min.TD50,min.a), upper=c(max.TD50,max.a),
                     C.I.gamma50=newg50,dvh.matrix=dvh.matrix,outcome=outcome)
    it<-it+newout$iterations
    while (abs(-newout$objective-MLE+ca1)>1e-5) {
      # grid search of D50 value      
      if ((-newout$objective-MLE+ca1)>0) {   # moving around the zero
        lowg50 <- newg50
      } else highg50 <- newg50
      newg50<-(highg50+lowg50)/2
      newout <- nlminb(start=c(D50,aaa),objective=nLL.gamma50, lower=c(min.TD50,min.a), upper=c(max.TD50,max.a),
                       C.I.gamma50=newg50,dvh.matrix=dvh.matrix,outcome=outcome)
      it<-it+newout$iterations      
    }
    result <- (highg50+lowg50)/2 # final value of higher bound of D50 C.I.
  }
  
  # lower bound of a C.I.
  else if (arg=="aL") {
    if (aaa>5)  inc <- 1 else inc <- 0.25
    higha <- aaa                 # starting point for interval search
    lowa <- aaa - inc            # ending point for interval search
    lowout <- nlminb(start=c(D50,g50),objective=nLL.a, lower=c(min.TD50,min.gamma50), upper=c(max.TD50,max.gamma50),
                     C.I.a=lowa,dvh.matrix=dvh.matrix,outcome=outcome)
    it<-lowout$iterations
    while((-lowout$objective-MLE+ca1)>0) { # routine to be repeated until LL-MLE+Chi2/2>0
      higha <- lowa
      lowa <- lowa - inc
      lowout <- nlminb(start=c(D50,g50),objective=nLL.a, lower=c(min.TD50,min.gamma50), upper=c(max.TD50,max.gamma50),
                       C.I.a=lowa,dvh.matrix=dvh.matrix,outcome=outcome)
      it<-it+lowout$iterations
      if (it > 10000) stop("Failed convergence for C.I.")
    }
    # starting point for bisection
    higha <<- higha
    lowa <<- lowa
    # start bisection search
    newa <- (higha+lowa)/2
    newout <- nlminb(start=c(D50,g50),objective=nLL.a, lower=c(min.TD50,min.gamma50), upper=c(max.TD50,max.gamma50),
                     C.I.a=newa,dvh.matrix=dvh.matrix,outcome=outcome)
    it<-it+newout$iterations
    while (abs(-newout$objective-MLE+ca1)>1e-5) {
      # grid search of D50 value      
      if ((-newout$objective-MLE+ca1)>0) {   # moving around the zero
        higha <- newa
      } else lowa <- newa
      newa<-(higha+lowa)/2
      newout <- nlminb(start=c(D50,g50),objective=nLL.a, lower=c(min.TD50,min.gamma50), upper=c(max.TD50,max.gamma50),
                       C.I.a=newa,dvh.matrix=dvh.matrix,outcome=outcome)
      it<-it+newout$iterations      
    }
    result <- (higha+lowa)/2 # final value of lower bound of a C.I.
  }
  
  # higher bound of a C.I.
  else if (arg=="aU") {
    if (aaa>5)  inc <- 1 else inc <- 0.25
    lowa <- aaa                   # starting point for interval search
    higha <- aaa + inc            # ending point for interval search
    highout <- nlminb(start=c(D50,g50),objective=nLL.a, lower=c(min.TD50,min.gamma50), upper=c(max.TD50,max.gamma50),
                      C.I.a=higha,dvh.matrix=dvh.matrix,outcome=outcome) 
    it<-highout$iterations
    while((-highout$objective-MLE+ca1)>0) { # routine to be repeated until LL-MLE+Chi2/2>0
      lowa <- higha
      higha <- higha + inc
      highout <- nlminb(start=c(D50,g50),objective=nLL.a, lower=c(min.TD50,min.gamma50), upper=c(max.TD50,max.gamma50),
                        C.I.a=higha,dvh.matrix=dvh.matrix,outcome=outcome)
      it<-it+highout$iterations
      if (it > 10000) stop("Failed convergence for C.I.")
    }
    # starting point for bisection
    lowa <<- lowa
    higha <<- higha
    # start bisection search
    newa<-(higha+lowa)/2
    newout <- nlminb(start=c(D50,g50),objective=nLL.a, lower=c(min.TD50,min.gamma50), upper=c(max.TD50,max.gamma50),
                     C.I.a=newa,dvh.matrix=dvh.matrix,outcome=outcome)
    it<-it+newout$iterations
    while (abs(-newout$objective-MLE+ca1)>1e-5) {
      # grid search of D50 value      
      if ((-newout$objective-MLE+ca1)>0) {   # moving around the zero
        lowa <- newa
      } else higha <- newa
      newa<-(higha+lowa)/2
      newout <- nlminb(start=c(D50,g50),objective=nLL.a, lower=c(min.TD50,min.gamma50), upper=c(max.TD50,max.gamma50),
                       C.I.a=newa,dvh.matrix=dvh.matrix,outcome=outcome)
      it<-it+newout$iterations      
    }
    result <- (higha+lowa)/2 # final value of higher bound of D50 C.I.
  }
  # vector of result for one side C.I. and corresponding other parameters and iteraions number
  bound <- c(result, newout$par[1], newout$par[2], it)
  # return one value of C.I. for each time and corresponding other parameters optimized
  # displays messages during C.I. calculation
  if (verbose==TRUE) {
    message(par0, format(result, digits=3, nsmall=3, justify="right", width=8))
    message(par1, format(newout$par[1], digits=3, nsmall=3, justify="right", width=8))
    message(par2, format(newout$par[2], digits=3, nsmall=3, justify="right", width=8))
    message("             Iterations = ", format(it, justify="right", width=8), "\n")
  }
  return(bound) 
}

# function for finding root of not linear equations using bisection method
bisect <- function(f, low, high, delta=1e-8) {
  a<-low
  b<-high
  d<-(a + b)/2
  it<-0
  errorbound <- (abs(b - a))/2
  if (f(a) == 0) return(list("root"=a, "iterations"=it))    
  if (f(b) == 0) return(list("root"=b, "iterations"=it))    
  while (errorbound > delta) {
    it<-it+1    
    if (f(d) == 0) return(list("root"=d, "iterations"=it))    
    if ((sign(f(a))*sign(f(d))) < 0) b<-d else a<-d
    d<-(a+b)/2
    errorbound<-errorbound/2    
  }
  if ((abs(d-low)<delta) ||(abs(b-high)<delta)) return(list("root"=NA, "iterations"=it)) else
    return(list("root"=d, "iterations"=it))
}


# function for calculating profile likelihood of DVH based model
proflik <- function(fitted.model, dvh.matrix, min.TD50=40, max.TD50=50,
                    min.gamma50=0.1, max.gamma50=5, min.a=-100, max.a=100, nbins=10) {
  if (fitted.model$model=="Lyman") {
    nLL.TD50 <- probit.nLL.Lyman.TD50
    nLL.gamma50 <- probit.nLL.Lyman.gamma50
    nLL.a <- probit.nLL.Lyman.a
  }
  if (fitted.model$model=="Niemierko") {
    nLL.TD50 <- logit.nLL.Niemierko.TD50
    nLL.gamma50 <- logit.nLL.Niemierko.gamma50
    nLL.a <- logit.nLL.Niemierko.a
  }
  message("Generating profile likelihood for TD50")  
  # create the matrix for TD50 profile likelihood
  proflik.TD50.matrix <- matrix(nrow=(nbins + 2), ncol=2)
  # create the vector of the TD50 and put into the matrix
  proflik.TD50.matrix[,1] <- sort(c(seq(from=min.TD50, to=max.TD50, by=((max.TD50 - min.TD50)/nbins)), fitted.model$TD50))
  # creates the progress bar
  pb <- txtProgressBar(min = 0, max = nbins + 2, style = 3)
  # creating the values of MLE for fixed parameter
  for (m in 1:nrow(proflik.TD50.matrix)) {
    setTxtProgressBar(pb, m)
    out<-nlminb(start=c(fitted.model$gamma50,fitted.model$a),objective=nLL.TD50, lower=c(0,0.01), upper=c(5,1000),
                C.I.TD50=proflik.TD50.matrix[m, 1],dvh.matrix=dvh.matrix,outcome=fitted.model$outcome)
    if (out$objective==0) proflik.TD50.matrix[m,2]<-NA else proflik.TD50.matrix[m,2]<- -out$objective
  }
  close(pb)  # closes the progress bar
  message("Generating profile likelihood for gamma50")  
  # create the matrix for gamma50 profile likelihood
  proflik.gamma50.matrix <- matrix(nrow=(nbins + 2), ncol=2)
  # create the vector of the TD50 and put into the matrix
  proflik.gamma50.matrix[,1] <- sort(c(seq(from=min.gamma50, to=max.gamma50, by=((max.gamma50 - min.gamma50)/nbins)), fitted.model$gamma50))
  # creates the progress bar
  pb <- txtProgressBar(min = 0, max = nbins + 2, style = 3)
  # creating the values of MLE for fixed parameter
  for (m in 1:nrow(proflik.gamma50.matrix)) {
    setTxtProgressBar(pb, m)
    out<-nlminb(start=c(fitted.model$TD50,fitted.model$a),objective=nLL.gamma50, lower=c(0,0.01), upper=c(100,1000),
                C.I.gamma50=proflik.gamma50.matrix[m, 1],dvh.matrix=dvh.matrix,outcome=fitted.model$outcome)
    if (out$objective==0) proflik.gamma50.matrix[m,2]<-NA else proflik.gamma50.matrix[m,2]<- -out$objective
  }
  close(pb)  # closes the progress bar
  message("Generating profile likelihood for a")
  # create the matrix for a profile likelihood
  proflik.a.matrix <- matrix(nrow=(nbins + 2), ncol=2)
  # create the vector of the TD50 and put into the matrix
  proflik.a.matrix[,1] <- sort(c(seq(from=min.a, to=max.a, by=((max.a - min.a)/nbins)), fitted.model$a))
  # creates the progress bar
  pb <- txtProgressBar(min = 0, max = nbins + 2, style = 3)
  # creating the values of MLE for fixed parameter
  for (m in 1:nrow(proflik.a.matrix)) {
    setTxtProgressBar(pb, m)
    out<-nlminb(start=c(fitted.model$TD50,fitted.model$gamma50),objective=nLL.a, lower=c(0,0), upper=c(100,5),
                C.I.a=proflik.a.matrix[m, 1],dvh.matrix=dvh.matrix,outcome=fitted.model$outcome)
    if (out$objective==0) proflik.a.matrix[m,2]<- NA else proflik.a.matrix[m,2]<- -out$objective
  }
  close(pb)  # closes the progress bar
  # output results
  return(list("TD50"=proflik.TD50.matrix, "gamma50"=proflik.gamma50.matrix, "a"=proflik.a.matrix))
}

# function for calculation of each element simulated in bootstrap during multicore calculation
calc.boot.el <- function(model, dvhnumber, maxdose, dosebin, TD50, gamma50, a) {
  boot.dvh.matrix <- gen.dvh(dvhnumber=dvhnumber, type="random", mindose=0, 
                             maxdose=maxdose, dosebin=dosebin, relative=TRUE)
  boot.outcome <- outcome.def(dvh.matrix=boot.dvh.matrix, TD50=TD50, gamma50=gamma50, 
                              a=a, model=model)
  boot.NTCP <- fit.NTCP.CI(model=model, dvh.matrix=boot.dvh.matrix, outcome=boot.outcome[,3])
  return(list("TD50"=boot.NTCP$TD50, "gamma50"=boot.NTCP$gamma50, "a"=boot.NTCP$a, "Iterations"=boot.NTCP$Iterations))
}

# function for calculation of C.I. of model using bootstrapping method
bootstrap.CI <- function(model, TD50, gamma50, a, n=1000, dvh.matrix, C.I.width=.95) {
  dvh.number<-ncol(dvh.matrix) - 1                          # find the number of DVH in the input DVH matrix 
  result.matrix<-matrix(nrow = n, ncol=4)                   # create the result matrix
  colnames(result.matrix)<-c("TD50", "gamma50", "a", "Iterations") # set the names of result matrix 
  max.dose<-max(dvh.matrix[,1])                             # set the maximum dose for bootstrapped DVHs
  dbin<-dvh.matrix[2,1]-dvh.matrix[1,1]                     # set the dose.bin for bootstrapped DVHs  
  if (Sys.info()["sysname"]=="Linux") {                     # check OS for using parallel computing under Linux
    require(foreach)
    require(doMC)
    require(parallel)
    registerDoMC(cores=detectCores())                       # compute the numbers of cores in CPU and register them
    boot.list <- foreach(m = 1:n) %dopar% calc.boot.el(model=model, dvhnumber=dvh.number, maxdose=max.dose, 
                                                       dosebin=dbin, TD50=TD50, gamma50=gamma50, a=a)
    for (p in 1:n) {
      result.matrix[p,1] <- boot.list[[p]]$TD50
      result.matrix[p,2] <- boot.list[[p]]$gamma50
      result.matrix[p,3] <- boot.list[[p]]$a
      result.matrix[p,4] <- boot.list[[p]]$Iterations
    }  
  } else {
    pb <- txtProgressBar(min = 0, max = n, style = 3)         # creates progressbar
    for (m in 1:n) {
      boot.dvh.matrix <- gen.dvh(dvhnumber=dvh.number, type="random", mindose=0, 
                                 maxdose=max.dose, dosebin=dbin, relative=TRUE)
      boot.outcome <- outcome.def(dvh.matrix=boot.dvh.matrix, TD50=TD50, gamma50=gamma50, 
                                  a=a, model=model)
      boot.NTCP <- fit.NTCP.CI(model=model, dvh.matrix=boot.dvh.matrix, outcome=boot.outcome[,3])
      # fills the value in the result matrix
      result.matrix[m,1] <- boot.NTCP$TD50
      result.matrix[m,2] <- boot.NTCP$gamma50
      result.matrix[m,3] <- boot.NTCP$a
      result.matrix[m,4] <- boot.NTCP$Iterations
      setTxtProgressBar(pb, m)
    }
    close(pb)
  }
  # output of confidence intervals by accessing into the result.matrix
  # find the bootstrapped series that is closest to the given C.I. values  
  TD50L<-which.min(abs(result.matrix[,1] - quantile(x=result.matrix[,1], probs=(1-C.I.width)/2)))
  TD50L<-c(as.numeric(result.matrix[TD50L,1]), as.numeric(result.matrix[TD50L,2]), 
                      as.numeric(result.matrix[TD50L,3]), as.numeric(result.matrix[TD50L,4]))
  TD50U<-which.min(abs(result.matrix[,1] - quantile(x=result.matrix[,1], probs=1-(1-C.I.width)/2)))
  TD50U<-c(as.numeric(result.matrix[TD50U,1]), as.numeric(result.matrix[TD50U,2]),
           as.numeric(result.matrix[TD50U,3]), as.numeric(result.matrix[TD50U,4]))
  gamma50L<-which.min(abs(result.matrix[,2] - quantile(x=result.matrix[,2], probs=(1-C.I.width)/2)))
  gamma50L<-c(as.numeric(result.matrix[gamma50L,2]), as.numeric(result.matrix[gamma50L,1]),
             as.numeric(result.matrix[gamma50L,3]), as.numeric(result.matrix[gamma50L,4]))
  gamma50U<-which.min(abs(result.matrix[,2] - quantile(x=result.matrix[,2], probs=1-(1-C.I.width)/2)))
  gamma50U<-c(as.numeric(result.matrix[gamma50U,2]), as.numeric(result.matrix[gamma50U,1]),
              as.numeric(result.matrix[gamma50U,3]), as.numeric(result.matrix[gamma50U,4]))
  aL<-which.min(abs(result.matrix[,3] - quantile(x=result.matrix[,3], probs=(1-C.I.width)/2)))
  aL<-c(as.numeric(result.matrix[aL,3]), as.numeric(result.matrix[aL,1]), 
                 as.numeric(result.matrix[aL,2]), as.numeric(result.matrix[aL,4]))
  aU<-which.min(abs(result.matrix[,3] - quantile(x=result.matrix[,3], probs=1-(1-C.I.width)/2)))
  aU<-c(as.numeric(result.matrix[aU,3]), as.numeric(result.matrix[aU,1]),
                 as.numeric(result.matrix[aU,2]), as.numeric(result.matrix[aU,4]))
  # generate output list
  output<-list("TD50L"=TD50L, "TD50U"=TD50U, "gamma50L"=gamma50L, "gamma50U"=gamma50U, 
               "aL"=aL, "aU"=aU, "bootstrap.matrix"=as.data.frame(result.matrix))
  return(output)
}


# function for calculating NTCP kernel from bootstrap matrix for a single DVH
NTCP.kernel <- function (model, dvh) {  
  doses <- c()   # empty vector of EUDs
  NTCP  <- c()   # empty vector of NTCP
  for (n in 1:nrow(model$bootstrap.matrix)) {
    doseV <- dvh[,1]^(model$bootstrap.matrix$a[n])  # calculates the doses
    addenda.EUD<-dvh[,2]*doseV
    doses<- c(doses, (sum(addenda.EUD))^(1/model$bootstrap.matrix$a[n])) # (apply(X=addenda.EUD,MARGIN=2,FUN=sum))^(1/a)
  }
  if (model$model == "Lyman") {
    for (n in 1:nrow(model$bootstrap.matrix))
    NTCP <- c(NTCP, DR.Lyman(TD50=model$bootstrap.matrix$TD50[n], 
                               gamma50=model$bootstrap.matrix$gamma50[n], 
                               dose=doses[n])) 
  }
  if (model$model == "Niemierko") {
    for (n in 1:nrow(model$bootstrap.matrix))
      NTCP <- c(NTCP, DR.Niemierko(TD50=model$bootstrap.matrix$TD50[n], 
                                 gamma50=model$bootstrap.matrix$gamma50[n], 
                                 dose=doses[n])) 
  }
  return(as.matrix(cbind(doses, NTCP)))
}

# function for fitting generic dose response function
fit.DR <- function(dvh, dose, outcome, start.TD50=50, start.gamma50=1.5, start.a=2,
                   DR.FUN=c("Lyman","Goitein","Niemierko","Munro","Okunieff","Warkentin","Bentzen"),
                   calcCI=TRUE, verbose=TRUE, C.I.width=.95, n.bins=200) {
  require(bbmle)
  # entry checks
  if (missingArg(dvh) && missingArg(dose)) stop("Please enter EITHER a dvhmatrix OR a dose vector")
  if (missingArg(outcome)) stop("You MUST provide an outcome vector for dose-response model fitting")
  if (!missingArg(dvh)) {
    if (class(dvh)!="dvhmatrix") warning("dvh isn't a dvhmatrix class object. If dvh isn't a relative differential dvh\n", 
                                         "  achieved results are wrong")
    if (class(dvh)=="dvhmatrix") dvh<-dvh@dvh
    if (ncol(dvh)!=(length(outcome)+1)) stop("Number of observation in dvhmatrix differs from outcome vector length")
  }
  fname<-DR.FUN
  DR.FUN=match.arg(arg=DR.FUN)
  DR.FUN=match.fun(FUN=paste("DR.", DR.FUN, sep=""))   # match the dose-response function
  # fits three parameters model
  if (missingArg(dose)) {
    if (verbose==TRUE) message("fitting 3 parameters ", fname, " model...")
    negLogLik<- function(TD50, gamma50, aa)      
      -sum(outcome*log(DR.FUN(diffdvh=dvh, TD50=TD50, gamma50=gamma50, aa=aa)) +
             (1-outcome)*log(1-DR.FUN(diffdvh=dvh, TD50=TD50, gamma50=gamma50, aa=aa)))
    model<-mle2(minuslogl=negLogLik, start=list(TD50=start.TD50, gamma50=start.gamma50, aa=start.a))
    predicted<-DR.FUN(diffdvh=dvh, TD50=model@coef[1], gamma50=model@coef[2], aa=model@coef[3])
    if (verbose==TRUE) message(fname, " model fitting successful")
    # calculate confidence interval of model
    if (calcCI==TRUE) {
      # calculates "fast" C.I. using confint function
      CI.model<-confint(model, level=C.I.width)
      # profiling Likelihood for a wider range of parameters
      # TD50 profiling
      low.TD50<-CI.model[1,1] - (CI.model[1,2]-CI.model[1,1])/20
      high.TD50<-CI.model[1,2] + (CI.model[1,2]-CI.model[1,1])/20
      TD50.seq<-seq(from=low.TD50, to=high.TD50, length.out=n.bins)    
      nLL.TD50.seq<-c()    
      negLogLikTD50<- function(gamma50, aa)      
        -sum(outcome*log(DR.FUN(diffdvh=model$dvh, TD50=TD50.seq[n], gamma50=gamma50, aa=aa)) +
               (1-outcome)*log(1-DR.FUN(diffdvh=model$dvh, TD50=TD50.seq[n], gamma50=gamma50, aa=aa)))
      cat("Profiling TD50...\n")
      pb <- txtProgressBar(min = 0, max = length(TD50.seq), style = 3)
      for (n in 1:length(TD50.seq)) {
        setTxtProgressBar(pb, n)
        nLL.TD50.seq<-c(nLL.TD50.seq, mle2(minuslogl=negLogLikTD50, start=list(gamma50=model$model@coef[2], aa=model$model@coef[3]))@min)
      }
      close(pb)
      # exploring right side of profile for TD50
      n<-which.min(nLL.TD50.seq)
      while (n<length(nLL.TD50.seq)) {
        if (((min(nLL.TD50.seq) + qchisq(p = .95, df = 1)/2 - nLL.TD50.seq[n])>0) &&  
               ((min(nLL.TD50.seq) + qchisq(p = .95, df = 1)/2 - nLL.TD50.seq[n + 1])<0)) {
          low<-n
          high<-n+1
          break
        }
      }
    }
  }
  # fits two parameters model
  if (!missingArg(dose)) {
    if (verbose==TRUE) message("fitting 2 parameters ", fname, " model...")
    negLogLik<- function(TD50, gamma50) 
      -sum(outcome*log(DR.FUN(dose=dose, TD50=TD50, gamma50=gamma50)) +
             (1-outcome)*log(1-DR.FUN(dose=dose, TD50=TD50, gamma50=gamma50)))
    model<-mle2(minuslogl=negLogLik, start=list(TD50=start.TD50, gamma50=start.gamma50))
    predicted<-DR.FUN(dose=dose, TD50=model@coef[1], gamma50=model@coef[2])
    if (verbose==TRUE) message(fname, " model fitting successful")
    # calculate confidence interval of model
    if (calcCI==TRUE) {
      # calculates "fast" C.I. using confint function
      CI.model<-confint(model, level=C.I.width)
      # start manual profiling
    }
  }
  # fills empty arguments
  if (missingArg(dose)) dose<-NULL
  if (missingArg(dvh))  dvh<-NULL
  # creates the output
  output <- list("model"=model, "CI.model"=CI.model, "Predicted"=predicted, 
                 "dvh.matrix"=dvh, "dose"=dose, "outcome"=outcome, "DR.FUN"=DR.FUN)
  # S3 class for DoseRespone
  class(output)<-"DoseResponse"
  return(output)  
}

# function for creating profile likelihood of a DoseResponse object
ProfLik.DR<-function(model, n.bins=200) {
  if (class(model)!="DoseResponse") stop("model must be a DoseResponse object")
  # profile likelihood for 3 parameters model
  DR.FUN<-model$DR.FUN
  outcome<-model$outcome
  if (is.null(model$dose)) {
    # profile likelihood of TD50
    low.TD50<-model$CI.model[1,1] - (model$CI.model[1,2]-model$CI.model[1,1])/20
    high.TD50<-model$CI.model[1,2] + (model$CI.model[1,2]-model$CI.model[1,1])/20
    TD50.seq<-seq(from=low.TD50, to=high.TD50, length.out=n.bins)    
    nLL.TD50.seq<-c()    
    negLogLikTD50<- function(gamma50, aa)      
      -sum(outcome*log(DR.FUN(diffdvh=model$dvh, TD50=TD50.seq[n], gamma50=gamma50, aa=aa)) +
             (1-outcome)*log(1-DR.FUN(diffdvh=model$dvh, TD50=TD50.seq[n], gamma50=gamma50, aa=aa)))
    cat("Profiling TD50...\n")
    pb <- txtProgressBar(min = 0, max = length(TD50.seq), style = 3)
    for (n in 1:length(TD50.seq)) {
      setTxtProgressBar(pb, n)
      nLL.TD50.seq<-c(nLL.TD50.seq, mle2(minuslogl=negLogLikTD50, start=list(gamma50=model$model@coef[2], aa=model$model@coef[3]))@min)
    }
    close(pb)
    # profile likelihood of gamma50
    low.gamma50<-model$CI.model[2,1] - (model$CI.model[2,2]-model$CI.model[2,1])/20
    high.gamma50<-model$CI.model[2,2] + (model$CI.model[2,2]-model$CI.model[2,1])/20
    gamma50.seq<-seq(from=low.gamma50, to=high.gamma50, length.out=n.bins)    
    nLL.gamma50.seq<-c()    
    negLogLikgamma50<- function(TD50, aa)      
      -sum(outcome*log(DR.FUN(diffdvh=model$dvh, gamma50=gamma50.seq[n], TD50=TD50, aa=aa)) +
             (1-outcome)*log(1-DR.FUN(diffdvh=model$dvh, gamma50=gamma50.seq[n], TD50=TD50, aa=aa)))
    cat("Profiling gamma50...\n")
    pb <- txtProgressBar(min = 0, max = length(gamma50.seq), style = 3)
    for (n in 1:length(gamma50.seq)) {
      setTxtProgressBar(pb, n)
      nLL.gamma50.seq<-c(nLL.gamma50.seq, mle2(minuslogl=negLogLikgamma50, start=list(TD50=model$model@coef[1], aa=model$model@coef[3]))@min)
    }
    close(pb)
    # profile likelihood of a
    low.a<-model$CI.model[3,1] - (model$CI.model[3,2]-model$CI.model[3,1])/20
    high.a<-model$CI.model[3,2] + (model$CI.model[3,2]-model$CI.model[3,1])/20
    a.seq<-seq(from=low.a, to=high.a, length.out=n.bins)
    nLL.a.seq<-c()    
    negLogLika<- function(TD50, gamma50)      
      -sum(outcome*log(DR.FUN(diffdvh=model$dvh, aa=a.seq[n], TD50=TD50, gamma50=gamma50)) +
             (1-outcome)*log(1-DR.FUN(diffdvh=model$dvh, aa=a.seq[n], TD50=TD50, gamma50=gamma50)))
    cat("Profiling a...\n")
    pb <- txtProgressBar(min = 0, max = length(a.seq), style = 3)
    for (n in 1:length(a.seq)) {
      setTxtProgressBar(pb, n)
      nLL.a.seq<-c(nLL.a.seq, mle2(minuslogl=negLogLika, start=list(TD50=model$model@coef[1], gamma50=model$model@coef[2]))@min)
    }
    close(pb)
  }
  return(list("TD50.ProfLik"=cbind(TD50.seq, nLL.TD50.seq), "gamma50.ProfLik"=cbind(gamma50.seq, nLL.gamma50.seq),
              "a.ProfLik"=cbind(a.seq, nLL.a.seq)))
}
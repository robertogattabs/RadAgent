# function that calculates a logit model for predicting the outcome
# depending from a Vdose, used by fit.Vdose
glm.Vdose <- function(dvh.matrix, InputDose, formula) {
  VDose<<-Vdose(dvh.matrix=dvh.matrix, Dose=InputDose)
  f <- update.formula(formula, ~ . + VDose)
  if (type.binary=="binomial") out.model <<- glm(formula=f, family=binomial(link = "logit"))
  if (type.binary=="gaussian") out.model <<- glm(formula=f, family=gaussian(link = "identity"))
  if (type.binary=="poisson")  out.model <<- glm(formula=f, family=poisson(link  = "log"))  
  output<-out.model$aic
  return(output)
}

# function that finds Vdose and other covariates affecting a given outcome
fit.Vdose <- function(dvh.matrix, formula, step.dose=.5, type.binary=c("binomial", "gaussian", "poisson"),
                      range.min=NULL, range.max=NULL) {
  require(data.table)
  type.binary<<-match.arg(type.binary)
  formula <- as.formula(formula)
  # define the dose range for Vdose calculation
  if (is.null(range.min) || (range.min < 1)) range.min<-1
  if (is.null(range.max)) range.max<-round(max(dvh.matrix[,1]) - 1)
  # create sequence of doses
  seq.dose <- seq(from=range.min, to=range.max, by=step.dose)
  # create matrix of Vdoses
  seq.Vdose <- c()
  for (n in 1:length(seq.dose)) {
    seq.Vdose <- cbind(seq.Vdose, Vdose(dvh.matrix=dvh.matrix, Dose=seq.dose[n]))
  }
  # create sequence of glm models for getting the minimum aic and starting optimization
  temp.formula <- update.formula(formula, ~ . + V.dose)
  result.matrix <- c()
  for (n in 1:length(seq.dose)) {
    V.dose <<- seq.Vdose[,n]
    if (type.binary=="binomial") family=binomial(link = "logit")
    if (type.binary=="gaussian") family=gaussian(link = "identity")
    if (type.binary=="poisson")  family=poisson(link  = "log")
    temp.model <- glm(temp.formula, family)
    temp.p <- c()
    # create the temp vector of p values
    if (type.binary=="gaussian") for (m in 1:length(summary(temp.model)$coefficients[, "Pr(>|t|)"])) {
      temp.p <- c(temp.p, summary(temp.model)$coefficients[, "Pr(>|t|)"][m])
    } else for (m in 1:length(summary(temp.model)$coefficients[, "Pr(>|z|)"])) {
      temp.p <- c(temp.p, summary(temp.model)$coefficients[, "Pr(>|z|)"][m])
    }
    # add row in the matrix: Vdose, aic, list of P values
    temp.result<-c(seq.dose[n], temp.model$aic, temp.p)
    result.matrix<-rbind(result.matrix, temp.result)
  }
  # typecast result.matrix
  result.matrix<-as.matrix(result.matrix)
  rm(V.dose,envir=.GlobalEnv)
  # set the start Vdose at minimum aic found from series of aic
  start.Vdose <- which.min(x=result.matrix[,2]) * step.dose
  # finds the minimum value of aic for best modeling of Vdose
  out<-nlminb(start=start.Vdose, objective=glm.Vdose, dvh.matrix=dvh.matrix, formula=formula, 
              lower=range.min, upper=range.max)
  rm(VDose, envir=.GlobalEnv)
  # unname the matrix
  result.matrix<-unname(result.matrix)
  cnames<-c("Dose", "AIC")
  temp.p <- c()
  # create the temp vector of p values
  if (type.binary=="gaussian") for (m in 1:length(summary(out.model)$coefficients[, "Pr(>|t|)"])) {
    temp.p <- c(temp.p, summary(out.model)$coefficients[, "Pr(>|t|)"][m])
    # add the column names
    cnames <- c(cnames, paste("P Value ", names(out.model$coefficients[m])))
  } else for (m in 1:length(summary(out.model)$coefficients[, "Pr(>|z|)"])) {
    temp.p <- c(temp.p, summary(out.model)$coefficients[, "Pr(>|z|)"][m])
    # add the column names
    cnames <- c(cnames, paste("P Value ", names(out.model$coefficients[m])))
  }
  colnames(result.matrix) <- cnames
  result.matrix <- rbind(result.matrix, c(out$par, out.model$aic, temp.p))
  result.matrix <- data.table(result.matrix, key="Dose")
  # new typecast after ordering table
  result.matrix <- as.matrix(result.matrix)
  output<-(list("model"=out.model, "Vdose"=out$par, "dvhmatrix"=dvh.matrix, 
                "resultmatrix"=result.matrix))
  class(output)<-"Vdose.fit"
  # clean environment from global variables
  rm(type.binary, out.model, envir=.GlobalEnv)
  return(output)
}


# output of fit.Vdose
print.Vdose.fit <- function (fittedmodel) {
  cat("Vdose = ", fittedmodel$Vdose, " Gy", "\n")
  print(fittedmodel$model)
}

# summary of fit.Vdose
summary.Vdose.fit <- function (fittedmodel) {
  cat("Vdose = ", fittedmodel$Vdose, " Gy", "\n")
#   # add the prediction for given levels of the Volume
#   # minimum volume
#   minV<-min(fittedmodel$dvhmatrix[,2:ncol(fittedmodel$dvhmatrix)])
#   # maximum volume
#   maxV<-max(fittedmodel$dvhmatrix[,2:ncol(fittedmodel$dvhmatrix)])
#   # sequence of volumes for prediction
#   seqV<-seq(from=minV, to=maxV, length.out=11)
#   # sequence of predicted probabilities
#   seqP<-predict(object=fittedmodel$model, newdata=data.frame(VDose=seqV), type="response")
#   seqV<-round(seqV, digits=2)
#   seqP<-round(seqP, digits=2)
#   cat("Volume =      ", seqV, "\n")
#   cat("Probability = ", seqP, "\n")
  summary(fittedmodel$model)
}

# function for plotting the series of P values of covariates
# and AIC related to Vdoses
plot.Vdose.fit <- function (fittedmodel, plot.legend=FALSE) {
  # draw the frame
  par(mar=c(3, 14, 0, 0) + 0.1)
  # plot the AIC curve
  plot(fittedmodel$resultmatrix[,1], fittedmodel$resultmatrix[,2], 
       axes=F, ylim=c(min(fittedmodel$resultmatrix[,2]),max(fittedmodel$resultmatrix[,2])), 
       xlab="", ylab="", type="l",col="black", main="",xlim=c(min(fittedmodel$resultmatrix[,1]),max(fittedmodel$resultmatrix[,1])))
  axis(2, ylim=c(min(fittedmodel$resultmatrix[,2],max(fittedmodel$resultmatrix[,2]))), col="black", lwd=1)
  mtext(2,text=colnames(fittedmodel$resultmatrix)[2] ,line=2)
  abline(v=fittedmodel$Vdose, lty=1, col="red")
  text(x=fittedmodel$Vdose, y=max(fittedmodel$resultmatrix[,2]), 
       labels=paste(round(fittedmodel$Vdose, digits=2), " Gy"), pos=4)
  legend.labels<-"AIC"
  legend.lty <- 1
  # draw the remaining covariates P values
  linepos<-3.5
  for (n in 3:ncol(fittedmodel$resultmatrix)){
    par(new=T)
    plot(fittedmodel$resultmatrix[,1], fittedmodel$resultmatrix[,n], lty=n-1, log="y", 
         axes=F, ylim=c(min(fittedmodel$resultmatrix[,n]),max(fittedmodel$resultmatrix[,n])), 
         xlab="", ylab="", type="l",col="black", main="",xlim=c(min(fittedmodel$resultmatrix[,1]),max(fittedmodel$resultmatrix[,1])))
    axis(2, ylim=c(min(fittedmodel$resultmatrix[,n]),max(fittedmodel$resultmatrix[,n])),lwd=1,line=linepos)
    mtext(2,text=colnames(fittedmodel$resultmatrix)[n] ,line=2+linepos)
    legend.labels<-c(legend.labels, colnames(fittedmodel$resultmatrix)[n])
    legend.lty <- c(legend.lty, n - 1)
    linepos<-linepos + 3.5
  }
  # plot the bottom axis with text
  axis(1, pretty(range(fittedmodel$resultmatrix[,1]),10), lwd=1)
  mtext("V-Dose Value (Gy)", side=1, col="black", line=2)
  # plot the legend if required
  if (plot.legend==TRUE) {
    legend(x=min(fittedmodel$resultmatrix[,1]),y=median(fittedmodel$resultmatrix[,n]),
           legend=legend.labels, lty=legend.lty)
  }
}

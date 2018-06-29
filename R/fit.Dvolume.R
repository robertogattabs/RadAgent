# function that calculates a model for predicting the outcome
# depending from a Dvolume, used by fit.Dvolume
glm.Dvolume <- function(dvh.matrix, InputVolume, formula) {
  DVolume<<-Dvolume(dvh.matrix=dvh.matrix, Volume=InputVolume)
  f <- update.formula(formula, ~ . + DVolume)
  if (type.binary=="binomial") out.model <<- glm(formula=f, family=binomial(link = "logit"))
  if (type.binary=="gaussian") out.model <<- glm(formula=f, family=gaussian(link = "identity"))
  if (type.binary=="poisson")  out.model <<- glm(formula=f, family=poisson(link  = "log"))  
  output<-out.model$aic
  return(output)
}

# function that finds Dvolume and other covariates affecting a given outcome
fit.Dvolume <- function(dvh.matrix, formula, num.steps.volume=100, type.binary=c("binomial", "gaussian", "poisson"),
                        range.min=NULL, range.max=NULL) {
  require(data.table)
  type.binary<<-match.arg(type.binary)
  formula <- as.formula(formula)
  # define the volume range for Dvolume calculation
  # maximum Volume is the minimum volume of the maxima in the dvh matrix
  if (is.null(range.max)) range.max<-min(apply(X=dvh.matrix[,2:ncol(dvh.matrix)], MARGIN=2, FUN=max))
  # fix the number of the steps for calculating the dose range
  if (is.null(range.min)) range.min<-min(dvh.matrix[,2:ncol(dvh.matrix)])
  # create sequence of volumes, find the minimum of the maxima in the volume values
  seq.volume <- seq(from=range.min, to=range.max, length.out=num.steps.volume)
  # delete the zero if present in the sequence and steps forward to the following Dvolume
  if (range.min==0) seq.volume<-seq.volume[2:length(seq.volume)]
  # create matrix of Dvolumes
  seq.Dvolume <- c()
  for (n in 1:length(seq.volume)) {
    seq.Dvolume <- cbind(seq.Dvolume, Dvolume(dvh.matrix=dvh.matrix, Volume=seq.volume[n]))
  }
  # create sequence of glm models for getting the minimum aic and starting optimization
  temp.formula <- update.formula(formula, ~ . + D.volume)
  result.matrix <- c()
  for (n in 1:length(seq.volume)) {
    D.volume <<- seq.Dvolume[,n]
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
    # add row in the matrix: Dvolume, aic, list of P values
    temp.result<-c(seq.volume[n], temp.model$aic, temp.p)
    result.matrix<-rbind(result.matrix, temp.result)
  }
  # typecast result.matrix
  result.matrix<-as.matrix(result.matrix)
  rm(D.volume, envir=.GlobalEnv)
  # set the start Dvolume at minimum aic found from series of aic
  start.Dvolume <- result.matrix [(which.min(x=result.matrix[,2])),1]
  # finds the minimum value of aic for best modeling of Dvolume
  out<-nlminb(start=start.Dvolume, objective=glm.Dvolume, dvh.matrix=dvh.matrix, formula=formula, 
              lower=range.min, upper=range.max)
  # unname the matrix
  result.matrix<-unname(result.matrix)
  cnames<-c("Volume", "AIC")
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
  result.matrix <- data.table(result.matrix, key="Volume")
  # new typecast after ordering table
  result.matrix <- as.matrix(result.matrix)
  output<-(list("model"=out.model, "Dvolume"=out$par, "dvhmatrix"=dvh.matrix, 
                "resultmatrix"=result.matrix))
  class(output)<-"Dvolume.fit"
  # clean environment from global variables
  rm(type.binary, out.model, envir=.GlobalEnv)
  return(output)
}


# output of fit.Dvolume
print.Dvolume.fit <- function (fittedmodel) {
  cat("Dvolume = ", fittedmodel$Dvolume, " cc", "\n")
  print(fittedmodel$model)
}

# summary of fit.Dvolume
summary.Dvolume.fit <- function (fittedmodel) {
  cat("Dvolume = ", fittedmodel$Dvolume, " cc", "\n")
  summary(fittedmodel$model)
}

# function for plotting the series of P values of covariates
# and AIC related to Dvolumes
plot.Dvolume.fit <- function (fittedmodel, plot.legend=FALSE) {
  # draw the frame
  par(mar=c(3, 14, 0, 0) + 0.1)
  # plot the AIC curve
  plot(fittedmodel$resultmatrix[,1], fittedmodel$resultmatrix[,2], 
       axes=F, ylim=c(min(fittedmodel$resultmatrix[,2]),max(fittedmodel$resultmatrix[,2])), 
       xlab="", ylab="", type="l",col="black", main="",xlim=c(min(fittedmodel$resultmatrix[,1]),max(fittedmodel$resultmatrix[,1])))
  axis(2, ylim=c(min(fittedmodel$resultmatrix[,2],max(fittedmodel$resultmatrix[,2]))), col="black", lwd=1)
  mtext(2,text=colnames(fittedmodel$resultmatrix)[2] ,line=2)
  abline(v=fittedmodel$Dvolume, lty=1, col="red")
  text(x=fittedmodel$Dvolume, y=max(fittedmodel$resultmatrix[,2]), 
       labels=paste(round(fittedmodel$Dvolume, digits=2), " cc"), pos=4)
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
  mtext("D-Volume Value (cc)", side=1, col="black", line=2)
  # plot the legend if required
  if (plot.legend==TRUE) {
    legend(x=min(fittedmodel$resultmatrix[,1]),y=median(fittedmodel$resultmatrix[,n]),
           legend=legend.labels, lty=legend.lty)
  }
}
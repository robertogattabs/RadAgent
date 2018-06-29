#' class for glcmTexturalFeatures
#'
#' @description  glcmTexturalFeatures
#' @export
#' @import radiomics data.table

glcmTexturalFeatures <- function(imgObj){

  # compute number of non-NA voxels

  nVoxel <- dim(imgObj)[1]*dim(imgObj)[2]*dim(imgObj)[3] - sum(is.na(imgObj))

  ### compute Gray Levels Cooccurrence Matrices

  G_list <- list()

  #Compute grey level cooccurrence matrices for 4 different directions within each slice
  for (i in 1:dim(imgObj)[3]){

    if (min(imgObj[,,i],na.rm = T) < 0) {imgObj[,,i] <- imgObj[,,i] + abs(min(imgObj[,,i],na.rm = T))}
    #if (max(imgObj[,,i],na.rm = T) <= 10) stop("Controlla virgole")
    #else 
    imgObj[,,i] <- as.integer(imgObj[,,i])

    if(length(imgObj[,,i]) - sum(is.na(imgObj[,,i]))<50) next
    G_list[[(i-1)*4+1]] <- as.matrix(glcm(imgObj[,,i], angle = 0, normalize = F,verbose=F))
    G_list[[(i-1)*4+2]] <- as.matrix(glcm(imgObj[,,i], angle = 45, normalize = F,verbose=F))
    G_list[[(i-1)*4+3]] <- as.matrix(glcm(imgObj[,,i], angle = 90, normalize = F,verbose=F))
    G_list[[(i-1)*4+4]] <- as.matrix(glcm(imgObj[,,i], angle = 135, normalize = F,verbose=F))
  }

if(length(which(sapply(G_list, is.null)))!=0){
  G_list = G_list[-which(sapply(G_list, is.null))]
}

  #Initialise data table for storing GLCM features; I have added a few
  featNames <- c("F_cm.joint.max", "F_cm.joint.avg", "F_cm.joint.var", "F_cm.joint.entr",
                 "F_cm.diff.avg", "F_cm.diff.var", "F_cm.diff.entr",
                 "F_cm.sum.avg", "F_cm.sum.var", "F_cm.sum.entr", "F_cm.energy","F_cm.contrast","F_cm.dissimilarity",
                 "F_cm.inv.diff","F_cm.inv.diff.norm","F_cm.inv.diff.mom","F_cm.inv.diff.mom.norm","F_cm.inv.var",
                 "F_cm.corr","F_cm.auto.corr","F_cm.clust.tend","F_cm.clust.shade","F_cm.clust.prom",
                 "F_cm.info.corr.1","F_cm.info.corr.2")
  F_cm <- data.table(matrix(NA, nrow=length(G_list), ncol=length(featNames)))
  colnames(F_cm) <- featNames

  #Iterate over grey level cooccurrence matrices
  #The idea is that basically every GLCM is the same, i.e. we can just perform the same operations on every glcm.
  for (iter in seq_len(length(G_list))){

    Ng <- ncol(G_list[[iter]])

    #Convert matrix to data table
    df.G    <- data.table(G_list[[iter]])

    #Add row grey level intensity
    df.G$i <- as.numeric(row.names(G_list[[iter]]))

    #Convert from wide to long table. This is the preferred format for data tables and data frames
    df.G   <- melt(df.G, id.vars="i", variable.name="j", value.name="n", variable.factor=FALSE)

    #Convert j from string to numeric
    df.G$j <- as.numeric(df.G$j)

    #Remove combinations with 0 counts
    df.G   <- df.G[n>0,]

    #Convert Grey level coccurrence matrix to joint probability
    df.p_ij <- df.G[,.(p_ij=n/sum(df.G$n)), by=.(i,j)] #joint probability
    df.p_i  <- df.p_ij[,.(p_i=sum(p_ij)), by=i]        #marginal probability over columns
    df.p_j  <- df.p_ij[,.(p_j=sum(p_ij)), by=j]        #marginal probability over rows

    #Diagonal probabilities (p(i-j))
    #First, we create a new column k which contains the absolute value of i-j.
    #Second, we sum the joint probability where k is the same.
    #This can written as one line by chaining the operations.
    df.p_imj <- copy(df.p_ij)
    df.p_imj <- df.p_imj[,"k":=abs(i-j)][,.(p_imj=sum(p_ij)), by=k]

    #Cross-diagonal probabilities (p(i+j))
    #Again, we first create a new column k which contains i+j
    #Second, we sum the probability where k is the same.
    #This is written in one line by chaining the operations.
    df.p_ipj <- copy(df.p_ij)
    df.p_ipj <- df.p_ipj[,"k":=i+j][,.(p_ipj=sum(p_ij)), by=k]

    #Merger of df.p_ij, df.p_i and df.p_j
    df.p_ij <- merge(x=df.p_ij, y=df.p_i, by="i")
    df.p_ij <- merge(x=df.p_ij, y=df.p_j, by="j")

    #Thus we have five probabiF_cmlity matrices
    #Joint probability:          df.p_ij with probabilities p_ij, p_i and p_j, and indices i, j
    #Marginal probability:       df.p_i with probability p_i, and index i
    #Marginal probability:       df.p_j with probability p_j, and index j
    #Diagonal probability:       df.p_imj with probability p_imj and index k
    #Cross-diagonal probability: df.p_ipj with probability p_ipj and index k

    #Calculate features

    #Joint maximum
    F_cm$F_cm.joint.max[iter]         <- max(df.p_ij$p_ij)

    #Joint average
    F_cm$F_cm.joint.avg[iter]         <- sum(df.p_ij$i * df.p_ij$p_ij)

    #Joint variance
    mu <- sum(df.p_ij$i * df.p_ij$p_ij)
    F_cm$F_cm.joint.var[iter]         <- sum((df.p_ij$i-mu)^2 * df.p_ij$p_ij)

    #Joint entropy
    F_cm$F_cm.joint.entr[iter]        <- - sum(df.p_ij$p_ij * log2(df.p_ij$p_ij))

    #Difference average
    F_cm$F_cm.diff.avg[iter]          <- sum((df.p_imj$k) * df.p_imj$p_imj)

    #Difference variance
    mu <- sum((df.p_imj$k) * df.p_imj$p_imj)
    F_cm$F_cm.diff.var[iter]          <- sum((df.p_imj$k-mu)^2 * df.p_imj$p_imj)

    #Difference entropy
    F_cm$F_cm.diff.entr[iter]         <- - sum(df.p_imj$p_imj * log2(df.p_imj$p_imj))

    #Sum average
    F_cm$F_cm.sum.avg[iter]           <- sum(df.p_ipj$k * df.p_ipj$p_ipj)

    #Sum variance
    mu <- sum(df.p_ipj$k * df.p_ipj$p_ipj)
    F_cm$F_cm.sum.var[iter]           <- sum((df.p_ipj$k-mu)^2 * df.p_ipj$p_ipj)

    #Sum entropy
    F_cm$F_cm.sum.entr[iter]          <- - sum(df.p_ipj$p_ipj * log2(df.p_ipj$p_ipj))

    #Angular second moment
    F_cm$F_cm.energy[iter]          <- sum(df.p_ij$p_ij^2)

    #Contrast
    F_cm$F_cm.contrast[iter]          <- sum((df.p_ij$i - df.p_ij$j)^2 * df.p_ij$p_ij)

    #Dissimilarity
    F_cm$F_cm.dissimilarity[iter]          <- sum(abs(df.p_ij$i - df.p_ij$j) * df.p_ij$p_ij)

    #Inverse difference
    F_cm$F_cm.inv.diff[iter]          <- sum( df.p_ij$p_ij/(1+abs(df.p_ij$i - df.p_ij$j)) )

    #Inverse difference normalized
    D <- abs(df.p_ij$i - df.p_ij$j)/Ng
    F_cm$F_cm.inv.diff.norm[iter]          <- sum(df.p_ij$p_ij / (1+D) )

    #Inverse difference moment
    F_cm$F_cm.inv.diff.mom[iter]          <- sum(df.p_ij$p_ij/(1+(df.p_ij$i - df.p_ij$j)^2))

    #Inverse difference moment normalized
    DD <- ((df.p_ij$i - df.p_ij$j)/Ng)^2
    F_cm$F_cm.inv.diff.mom.norm[iter]          <- sum(df.p_ij$p_ij / (1+DD))

    #Inverse variance
    df.jmorethani <- df.p_ij[which(df.p_ij$j>df.p_ij$i),]
    F_cm$F_cm.inv.var[iter]          <- 2*sum(df.jmorethani$p_ij/(df.jmorethani$i-df.jmorethani$j)^2)

    #Correlation
    mu.i <- sum(df.p_i$i * df.p_i$p_i)
    sigma.i <- sqrt(sum((df.p_i$i - mu.i)^2 * df.p_i$p_i))
    F_cm$F_cm.corr[iter] <- (1/sigma.i)^2 * sum((df.p_ij$i - mu.i) * (df.p_ij$j - mu.i) * df.p_ij$p_ij)

    #Autocorrelation
    F_cm$F_cm.auto.corr[iter] <- sum(df.p_ij$i * df.p_ij$j * df.p_ij$p_ij)

    #Cluster tendency
    F_cm$F_cm.clust.tend[iter] <- sum((df.p_ij$i + df.p_ij$j -2*mu.i)^2 * df.p_ij$p_ij)

    #Cluster shade
    F_cm$F_cm.clust.shade[iter] <- sum((df.p_ij$i + df.p_ij$j -2*mu.i)^3 * df.p_ij$p_ij)

    #Cluster prominence
    F_cm$F_cm.clust.prom[iter] <- sum((df.p_ij$i + df.p_ij$j -2*mu.i)^4 * df.p_ij$p_ij)

    #First measure of information correlation
    HXY <- - sum(df.p_ij$p_ij * log2(df.p_ij$p_ij))
    HX <- - sum(df.p_i$p_i * log2(df.p_i$p_i))
    HXY1 <- - sum(df.p_ij$p_ij * log2(df.p_ij$p_i * df.p_ij$p_j))
    F_cm$F_cm.info.corr.1[iter] <- (HXY - HXY1) / HX

    #Second measure of information correlation
    HXY2 <- - sum(df.p_ij$p_i * df.p_ij$p_j * log2(df.p_ij$p_i * df.p_ij$p_j))
    if (HXY > HXY2 ) F_cm$F_cm.info.corr.2[iter] <- 0
    else if (HXY <= HXY2 ) F_cm$F_cm.info.corr.2[iter] <- sqrt(1-exp(-2*(HXY2 - HXY)))

  }

  return(F_cm)
}


######################################
### END FUNCION ######################
######################################

######################################
### NOW RUN IT ######################
######################################

# imgObj <- array(data=0,dim=c(5,4,4))
# imgObj[,,1] <- matrix(c(1,4,4,1,1,   1,4,6,1,1,  4,1,6,4,1,  4,4,6,4,1),nrow=4,ncol=5)
# imgObj[,,2] <- matrix(c(1,4,4,1,1,   1,1,6,1,1,  NA,1,3,1,1, 4,4,6,1,1),nrow=4,ncol=5)
# imgObj[,,3] <- matrix(c(1,4,4,NA,NA, 1,1,1,1,1,  1,1,NA,1,1, 1,1,6,1,1),nrow=4,ncol=5)
# imgObj[,,4] <- matrix(c(1,4,4,NA,NA, 1,1,1,1,1,  1,1,1,1,1,  1,1,6,1,1),nrow=4,ncol=5)
#
# F_cm <- glcmTexturalFeatures(imgObj)
#
# #compute mean on slice
# mean(F_cm$F_cm.diff.avg)
# mean(F_cm$F_cm.sum.avg)

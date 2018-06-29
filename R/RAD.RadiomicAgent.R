#' class for loading and presenting DICOM datajjjj
#' 
#' @description  Iljlkjlk
#' @export
#' @import corrplot pROC
#' @useDynLib moddicomV2 
RAD.RadiomicAgent <- function ( ){
  
  param.cross.correlation.threshold <- NA
  param.wilcoxon.pValue.threshold <- NA
  param.ignorance.threshold<-NA
  
  global.ROIName <- NA
  global.feature.family <- NA
  global.arr.sigma <- NA
  global.extracted.cov <- NA            # Contiene le matrici dei valori delle features

  # --------------------------------------------------
  # get.matrice.piatta.covariate
  # - arr.clinicalOutcome : array degli outcome clinici (array posizionale con nei nomi colonna gli ID dei pazienti)
  # - CalcoloSigma : la lista di matrici, come restituita da 'f.extractor.pluri.LoG.par'
  # --------------------------------------------------
  get.matrice.piatta.covariate<-function( arr.clinicalOutcome , CalcoloSigma ){
    # Inizio sezione  un po' misteriosa che dovrebbe appiattire la lista di matrici ottenuta
    # precedentemente in un unica matriciona, inclusiva dell'outcome clinico 
    Dataset <- data.frame("CS"= as.character( names(arr.clinicalOutcome) ), "outcome"=as.character(arr.clinicalOutcome))
    Dataset$CS <- as.character(Dataset$CS)
    # prendi i nomi delle cartelle (ultima parte del path)
    arr.ID <- lapply( strsplit(CalcoloSigma[[1]][,1],"//") , function(x){  x[2] } )
    # Ora prendi la lista degli outcome, facendo il match con gli ID estratti dai nomi delle cartelle
    arr.outcome.riposizionati <- unlist(
      lapply(arr.ID,
             function(x){  
               tmp.riga <- which( names(arr.clinicalOutcome)==x )
               if(length(tmp.riga)!=1) stop("Nell'array dei pazienti passati, uno non ha il nome corretto")
               return( arr.clinicalOutcome[ tmp.riga ]   )
             }) )
    
    # inizia a popolare la matrice finale
    bigImplodedMatrix <- c()
    bigImplodedMatrix <- cbind(bigImplodedMatrix, CalcoloSigma[[1]][,1])
    bigImplodedMatrix <- cbind(bigImplodedMatrix, arr.ID)
    bigImplodedMatrix <- cbind(bigImplodedMatrix, arr.outcome.riposizionati)
    colnames(bigImplodedMatrix)<-c("path","ID","outcome")
    
    # Appendi le colonne di tutti gli altri fogli contenenti le feature-values ai vari sigma
    # (passa da una lista di matrici ad una matriciona)
    for( foglio in names(CalcoloSigma)) {
      # Estrai la matrice dal foglio di interesse (escludi la prima colonna perche' e' l'ID)
      tmp.subMatrix <- CalcoloSigma[[foglio]][ ,  seq(2,ncol(CalcoloSigma[[foglio]])) ]
      # CAlcola i nuovi nomi (anteponendo il sigma al nome delle colonne)
      arr.nuovi.nomi <- lapply(X = colnames(tmp.subMatrix), FUN = function(x){ str_c(foglio,"_",x) } )
      colnames(tmp.subMatrix) <- arr.nuovi.nomi
      # appendilo alla matrice
      bigImplodedMatrix <- cbind(bigImplodedMatrix,tmp.subMatrix)
    }
    return(bigImplodedMatrix)
  }

  # --------------------------------------------------
  # get.matrice.miglior.predittore.singola.features
  # - arr.clinicalOutcome : array degli outcome clinici (array posizionale con nei nomi colonna gli ID dei pazienti)
  # - CalcoloSigma : la lista di matrici, come restituita da 'f.extractor.pluri.LoG.par'
  # ---------------------------------------------------
  # accorpa le colonne corrispondenti, dei vari fogli, al fine di costruire una lista 
  # di matrici una per ogni covariata (anzichè per ogni sigma!). Già che ci sei, fai sì che le colonne
  # siano elencate per valori di sigma crescenti (auspicabilmente)
  # --------------------------------------------------
  get.matrice.miglior.predittore.singola.features<-function( arr.clinicalOutcome , CalcoloSigma ) {

    # prendi i nomi delle cartelle (ultima parte del path)
    arr.ID <- lapply( strsplit(CalcoloSigma[[1]][,1],"//") , function(x){  x[2] } )
    # Ora prendi la lista degli outcome, facendo il match con gli ID estratti dai nomi delle cartelle
    arr.outcome.riposizionati <- unlist(
      lapply(arr.ID,
             function(x){  
               tmp.riga <- which( names(arr.clinicalOutcome)==x )
               if(length(tmp.riga)!=1) stop("Nell'array dei pazienti passati, uno non ha il nome corretto")
               return( arr.clinicalOutcome[ tmp.riga ]   )
             }) )    
    
    
    # inizia a popolare la matrice che contiene le prime colonne "tipo" di ogni futura sottomatrice
    tmp.matrice.base <- c()
    tmp.matrice.base <- cbind(tmp.matrice.base, CalcoloSigma[[1]][,1])
    tmp.matrice.base <- cbind(tmp.matrice.base, arr.ID)
    tmp.matrice.base <- cbind(tmp.matrice.base, arr.outcome.riposizionati)
    colnames(tmp.matrice.base)<-c("path","ID","outcome")
    
    nomi.colonne <- colnames(CalcoloSigma[[1]]);
    nomi.colonne <- nomi.colonne[2:length(nomi.colonne)]
    lista.matrici <- list()
    lista.correlazioni <- list()
    lista.wilcox.pValue <- list()

    # Fai una prima passata per iniziare a popolare le matrici con le prime colonne (ID, outcome, path)
    for(nome.feature in nomi.colonne) { lista.matrici[[nome.feature]]<-tmp.matrice.base }
    # Ora passa per andare a popolare opportunamente le colonne dei dati,
    # fare la cross-correlazione ed il wilcoxon
    for( nome.feature in nomi.colonne ) {
      colonne.estratte <- lapply( CalcoloSigma , FUN = function(x){ x[, nome.feature] }  )
      colonne.estratte <- matrix(unlist(colonne.estratte),nrow=nrow(CalcoloSigma[[1]]))
      lista.matrici[[nome.feature]]<-cbind(lista.matrici[[nome.feature]] , colonne.estratte)
      colnames(lista.matrici[[nome.feature]]) <- c(colnames(lista.matrici[[nome.feature]])[1:3],names(CalcoloSigma)   )
      # costruisci la matrice di cross-correlazione
      lista.correlazioni[[nome.feature]] <- cor(matrix(as.numeric(lista.matrici[[nome.feature]][,4:ncol(lista.matrici[[nome.feature]])]),nrow =nrow(lista.matrici[[nome.feature]]),byrow = F))
      # Fai il modello predittivo rispetto all'outcome
      m.negativi <- lista.matrici[[nome.feature]][ which(lista.matrici[[nome.feature]][,"outcome"]==1), 4:ncol(lista.matrici[[nome.feature]]) ]
      m.positivi <- lista.matrici[[nome.feature]][ which(lista.matrici[[nome.feature]][,"outcome"]==0), 4:ncol(lista.matrici[[nome.feature]]) ]
      # facciamo i Wilcoxon, uno per ogni colonna (cioe' per ogni sigma)
      lista.wilcox.pValue[[nome.feature]]<-c()
      cat("\n ",nome.feature)
      # if( nome.feature == "F_szm.z.entr") browser()
      for( sigma in colnames(m.negativi))  { 
        # cat("\n \t",sigma)
        # if(sigma == "5" & nome.feature == "F_szm.z.entr") browser()
        tmp.nneg <- as.numeric(m.negativi[,as.character(sigma)])
        tmp.ppos <- as.numeric(m.positivi[,as.character(sigma)])
        # if(nome.feature=="F_stat.old.skewness")  browser();
        if ( sum(tmp.nneg== Inf) >=1 | sum(tmp.ppos == Inf) >= 1 |
             sum(is.na(tmp.nneg)) >=1 | sum(is.na(tmp.ppos)) >=1   ) {
          lista.wilcox.pValue[[nome.feature]] <- c(lista.wilcox.pValue[[nome.feature]], 1 )
        }
        else { 
          lista.wilcox.pValue[[nome.feature]] <- c(lista.wilcox.pValue[[nome.feature]], wilcox.test(x = as.numeric(m.negativi[,as.character(sigma)]),y = as.numeric(m.positivi[,as.character(sigma)]))$p.value )
        }
      }
    }
    return(list(
      "lista.matrici"=lista.matrici,
      "lista.correlazioni"=lista.correlazioni,
      "lista.wilcox.pValue"=lista.wilcox.pValue
    ))
  }
  
  # --------------------------------------------------
  # scoutFeatures
  # --------------------------------------------------
  scoutFeatures <- function( pathName, ROIName, feature.family, 
                             arr.sigma, arr.clinicalOutcome, 
                             px="", py="", forceRecalculus =FALSE,
                             matrice.altre.covariate = c(),
                             cache.fileName = "tmp.f.extractor.pluri.par.RData") {

    # Copia i parametri di lancio cosi' che siano disponibili per l'applicaizone
    # di metodi a seguire
    global.ROIName <<- ROIName
    global.feature.family <<- feature.family
    global.arr.sigma <<- arr.sigma
    
    # Estrai le covariate al variare del sigma, come indicato dai parametri in ingresso
    CalcoloSigma <- f.extractor.pluri.LoG.par(path = pathName, ROIName = ROIName, 
                                              feature.family = feature.family,
                                              forceRecalculus = forceRecalculus, 
                                              sigma.array = arr.sigma, strategy = "none",
                                              fileName = cache.fileName)
    
    mm.1 <- get.matrice.piatta.covariate( arr.clinicalOutcome =arr.clinicalOutcome ,  CalcoloSigma = CalcoloSigma )
    mm.2 <- get.matrice.miglior.predittore.singola.features( arr.clinicalOutcome =arr.clinicalOutcome ,  CalcoloSigma = CalcoloSigma )
    # browser()
    # costruisci il data.frame per la logistica, partendo da una versione "numerica" di 
    # mm.1 cui aggiunge le variabili cliniche
    dataFrameDiTest <- mm.1
    dataFrameDiTest <- as.data.frame(mm.1)
    dataFrameDiTest$outcome <- as.factor(as.integer(dataFrameDiTest$outcome))
    for(i in seq(4,ncol(dataFrameDiTest)))  {	
      dataFrameDiTest[[i]] <- as.numeric(dataFrameDiTest[[i]])
    }
    # Aggiungi eventuali covariate cliniche
    mm.1.data.frame <- cbind(dataFrameDiTest,as.data.frame(matrice.altre.covariate))
    # browser()
    
    # anteponi a tutte le colonne "s_" o in fase di costruzione della formula il tutro si inchioda
    # quando trova una colonna il cui nome comincia per un numero (maremma che palle!)
    colnames(mm.1.data.frame) <- unlist(lapply( colnames(mm.1.data.frame) , function(x) str_c("s_",x)))

    # ----------------------------------------
    # NON CANCELLARE!!!
    # ----------------------------------------
    # ora plotta una cosa figa: la cross correlazione di tutte le covariate
    # (non serve ad una fava ma e' tanto fico!)
    # mm.1 <- m.features$mm.1
    # tmp.mtrx <- matrix(as.numeric(mm.1[,4:ncol(mm.1)]),nrow =nrow(mm.1),byrow = F)
    # colnames(tmp.mtrx)<- colnames(mm.1)[4:ncol(mm.1)]
    # colnames(tmp.mtrx)<- as.character(seq(1,ncol(tmp.mtrx)))
    # ooo <- cor(tmp.mtrx)
    # corrplot(cor(tmp.mtrx),method = "square")
    # corrplot(cor(tmp.mtrx),method = "color")
    # corrplot(cor(tmp.mtrx,use = "complete.obs"),method = "color",order = "AOE")
    # plot(density((ooo),na.rm = T))
    # ----------------------------------------

    global.extracted.cov  <<- list(
      "mm.1" = mm.1,
      "mm.1.data.frame"=mm.1.data.frame,
      "mm.2" = mm.2,
      "mm.other.clinical.cov"= matrice.altre.covariate
    )
    return( global.extracted.cov )
  }
  
  # --------------------------------------------------
  # featureSelection
  # --------------------------------------------------
  featureSelection <- function( number.of.cov = NA ) {
    
    # usa quanto precdentemente calcolato ed archiviato nell'attributo globale 
    # interno all'oggetto "this"
    mm.2 <- global.extracted.cov$mm.2
    mm.1 <- global.extracted.cov$mm.1
    mm.1.data.frame <- global.extracted.cov$mm.1.data.frame
    matrice.altre.covariate <- global.extracted.cov$mm.other.clinical.cov
    
    # carica la soglia di dichiarata ignoranza
    ignorance.threshold <- param.ignorance.threshold
    
    # in questa variabile ci finiranno i vincitori alla prima selezione
    selected.features <- c()
    
    # cicla per tutte le covariate, analizzandone la capacità di separare gli spazi rispetto all'
    # outcome in esame. Questo loop avviente per ogni covariata, nella matrice che ne contiene le 
    # performance al variare dei sigma (per la stessa feature)
    for( feature in names(mm.2$lista.wilcox.pValue) ) {
      
      riga <- c()
      # ordina in funzione del vincitore
      
      # - im 
      mm.2$lista.wilcox.pValue[[feature]][ which(is.na(mm.2$lista.wilcox.pValue[[feature]]))] <- 1
      # - mf
      
      ordinamento <- sort(mm.2$lista.wilcox.pValue[[feature]],index.return = T)
      
      # cat("\n feature: ",feature )
      # if( feature == "F_cm.info.corr.2" ) browser()
      # browser()
      # Vedi se il migliore ha performance sotto la soglia di p.value indicata
      if(ordinamento$x[1]<=param.wilcoxon.pValue.threshold)  {
        # Arruola il primo
        riga <- c( feature, colnames(mm.2$lista.matrici$F_stat.mean)[(3+ordinamento$ix[1])],ordinamento$x[1] )
        # appendilo alla matrice dei vincitori di questa prima selezione
        selected.features <- rbind(selected.features, riga)
        
        # Vedi se anche per il secondo. Analisi future più raffinate potranno anche considerare 
        # più sigma della stessa feature (per questa proof ci fermiamo a due)
        if(ordinamento$x[2]<=param.wilcoxon.pValue.threshold) {
          
          # estrai la p di pearson, per vedere la dipendenza fra primo e secondo
          p.pearson <-abs( mm.2$lista.correlazioni[[ feature]][ ordinamento$ix[1] , ordinamento$ix[2] ])
          if(p.pearson <= param.cross.correlation.threshold) {
            # Arruola il secondo
            riga <- c( feature, colnames(mm.2$lista.matrici$F_stat.mean)[(3+ordinamento$ix[2])],ordinamento$x[1] )
            # appendilo alla matrice dei vincitori di questa prima selezione
            selected.features <- rbind(selected.features, riga)
          }
        }
      }
    }
    
    # Costruisci la colonna con gli stessi nomi della matrice delle covariate complessive
    nomi.completi <- eee <- apply( selected.features, 1, function(x){  str_c(x["sigma"],"_",x["feature"]) } )  
    colnames(selected.features)<-c("feature","sigma","p.value")

    # ed appendila alla matrice appena costruita, nominando le colonne
    nomi.completi <- eee <- apply( selected.features, 1, function(x){  str_c(x["sigma"],"_",x["feature"]) } )  
    selected.features <- cbind(selected.features,"fullName"=nomi.completi)
    
    # OTTIMO! Ora prendi la prima matrice, con le cross correlazioni fra tutte le covariate, e 
    # guarda un po' come sono messe fra loro le veincitrici.
    # Partiamo da un data.frame, che è più facile da gestire in caso di ordinamento
    
    df.selected.features <- data.frame( 
      "feature" = selected.features[, "feature"],
      "sigma" = as.character(selected.features[, "sigma"]), 
      "p.value"=as.numeric(selected.features[, "p.value"]),
      "fullName" = selected.features[, "fullName"]   
    )
    df.selected.features$feature <- as.character(df.selected.features$feature)
    df.selected.features$sigma <- as.character(df.selected.features$sigma)
    df.selected.features$fullName <- as.character(df.selected.features$fullName)
    
    # Ordina il data.frame in funzione delle performance di p.value
    df.selected.features <- df.selected.features[order(df.selected.features$p.value),]

    # il numero di covariate da tener e' pari al numero di pazienti della classe di outcome meno numerosa
    # diviso 15. (in altre parole, per ogni covariata vogliamo che ci siano almeno 15 pazienti nella classe 
    # meno numerosa)
    if(!is.na(number.of.cov)) numero.covariate.da.tenere <- number.of.cov
    else numero.covariate.da.tenere <- as.integer(min(table(unlist(mm.2$lista.matrici[[1]][,"outcome"])))/15)
  
    covariate.reclutate <- c()
    tmp.str.covariate <- ""
    def.str.covariate <- ""
    
    # popola il dataFrame di test, pronti per la costruzione dei modelli
    dataFrameDiTest <- mm.1.data.frame
    
    # loopa fino a che non hai il numero di covariate desiderato
    covariate.incluse <- c();
    arr.acc <- c()
    lst.FP<-list(); lst.FN<-list()
    for( i in seq(1,numero.covariate.da.tenere)) {
      
      # testa ogni possibile covariata, da aggiungere al modello
      # (sia le radiomiche che le cliniche passate dall'utente)
      arr.auc<-c()
      insieme.complessivo.covariate <- c(df.selected.features$fullName, colnames(matrice.altre.covariate))
      lista.modelli <- list()
      for( covariata in insieme.complessivo.covariate ) { 
        
        # Non testare covariate gia' arruolate
        if(! (covariata %in% covariate.incluse) ) { 
          
          # testa il set di covariate in esame
          res <- test.these.cov( c(covariate.incluse,covariata) , ignorance.threshold = ignorance.threshold )        

          # popola l'array dell'accuratezza. In seguito diventera' una lista con anche 
          # altri indicatori di performance 
          arr.acc <- c(arr.acc, res$res.confusion.matrix$accuracy )
          auc.value <- res$auc ; names(auc.value) <- covariata
          arr.auc <- c(arr.auc, auc.value)

          # tieni via la lista dei modelli testati
          local.tmp.str.covariate <- res$local.tmp.str.covariate
          lista.modelli[[ local.tmp.str.covariate ]] <- list()
          lista.modelli[[ local.tmp.str.covariate ]]$modello <- res$modello
          lista.modelli[[ local.tmp.str.covariate ]]$roc <- res$roc.modello
          lista.modelli[[ local.tmp.str.covariate ]]$auc <- res$auc.value
          lista.modelli[[ local.tmp.str.covariate ]]$confusion.matrix <- res$confusion.matrix$conf.matrix
          # if(i == numero.covariate.da.tenere) lst.FP[[ local.tmp.str.covariate ]]<-res$confusion.matrix$conf.matrix["1","0"]
          # if(i == numero.covariate.da.tenere) lst.FN[[ local.tmp.str.covariate ]]<-res$confusion.matrix$conf.matrix["0","1"]

        }
      }
      # aggiusta per prendere quello con AUC migliore
      res.arr.auc <- sort(arr.auc, index.return = T,decreasing = T)
      covariate.incluse <- c(covariate.incluse, names(arr.auc)[res.arr.auc$ix[1]] )
      # aggiorna la stringa che poi diventera' la formula base per la logistica
      tmp.str.covariate <- paste(unlist(lapply( covariate.incluse , function(x) str_c("s_",x))),collapse = '+')
    }
    
    # Ora calcola (avrei potuto farlo anche prima, ma ero un po' pigro)
    # gli array degli altri indicatori di performance (accuracy, specif ....)
    # browser()
    aaa <- calcola.altri.indicatori.diperformance( models = lista.modelli)
    # browser()
    # Bon, ritona il tutto
    return( list(
                  "tmp.str.covariate" = tmp.str.covariate,
                  "covariate.incluse" = covariate.incluse,
                  "lista.modelli" = lista.modelli,
                  "arr.accuracy" = aaa$arr.accuracy,
                  "arr.sensitivity" = aaa$arr.sensitivity,
                  "arr.specificity" = aaa$arr.specificity,
                  "arr.AUC" = aaa$arr.AUC
                  # "lst.FP" = lst.FP,
                  # "lst.FN"  = lst.FN
                )
    )
  }
  calcola.altri.indicatori.diperformance<-function( models  ) {
     
    # Ora ricava le performance di tutti i modelli e buttale in un array
    lst.performance <- list()
    # browser()
    arr.accuracy <- c(); arr.AUC<-c(); arr.sensitivity <-c(); arr.specificity<-c()
    for( signature in names( models ) ) {
      
      cov.sporche <- unlist(str_split(string = signature,pattern = "\\+"))
      cov.pulite <- str_sub(string = cov.sporche, 3,str_length(cov.sporche))
      
      lst.performance[[signature]] <- test.these.cov( cov.pulite ) 
      
      arr.accuracy <- c(arr.accuracy, lst.performance[[signature]]$confusion.matrix$accuracy)
      arr.sensitivity <- c(arr.sensitivity, lst.performance[[signature]]$confusion.matrix$true.positive.rate)
      arr.specificity <- c(arr.specificity, lst.performance[[signature]]$confusion.matrix$true.negative.rate)
      arr.AUC <- c(arr.AUC, lst.performance[[signature]]$auc.value)
      
    }
    names(arr.accuracy) <- names( models )
    names(arr.sensitivity) <- names( models )
    names(arr.specificity) <- names( models )
    names(arr.AUC) <- names( models )    
    
    return(list( 
                  "arr.accuracy"=arr.accuracy,
                  "arr.sensitivity"=arr.sensitivity,
                  "arr.specificity"=arr.specificity,
                  "arr.AUC"=arr.AUC
                )
          )
  }
  # --------------------------------------------------
  # test.this.formula
  # Testa uno specifico set di covariate e ne restituisce la performance
  # --------------------------------------------------
  test.these.cov<-function( arr.cov ,ignorance.threshold=NA, debug.mode = FALSE, n.folder = 10 ) {
    if(debug.mode == TRUE) browser()
    # browser()
    mm.1 <- global.extracted.cov$mm.1
    dataFrameDiTest <- global.extracted.cov$mm.1.data.frame
    if(is.na(ignorance.threshold)) ignorance.threshold <- param.ignorance.threshold
    
    arr.s.cov <- str_c("s_",arr.cov)
    if(length(arr.cov)>1) {
      s.stringa <- paste( arr.s.cov , collapse = "+") 
    } else s.stringa <- arr.s.cov
    
    local.tmp.str.covariate <- s.stringa
    
    # if(tmp.str.covariate=="") local.tmp.str.covariate<-str_c("s_",covariata)
    # else local.tmp.str.covariate <- str_c(tmp.str.covariate,"+", str_c("s_",covariata ) ) 
    tmp.formula <- str_c("s_outcome~",s.stringa)
    
    # tmp.formula <- stringa
    tmp.formula <- as.formula(tmp.formula)
    # browser()
    if(is.factor(dataFrameDiTest$s_cT))
      dataFrameDiTest$s_cT <- as.numeric(levels(dataFrameDiTest$s_cT)[ as.numeric(dataFrameDiTest$s_cT) ])
    if(is.factor(dataFrameDiTest$s_cN))
      dataFrameDiTest$s_cN <- as.numeric(levels(dataFrameDiTest$s_cN)[ as.numeric(dataFrameDiTest$s_cN) ])
    # -im
    # k-folder-validation
    df     <- dataFrameDiTest[sample(nrow(dataFrameDiTest)),]
    folds  <- cut(seq(1,nrow(df)),breaks=n.folder,labels=FALSE)
    result <- list()
    arr.auc.modello <- c()
    # browser()
    for(i in 1:n.folder){
      testIndexes <- which(folds==i,arr.ind=TRUE)
      testData    <- df[testIndexes, ]
      trainData   <- df[-testIndexes, ]
      # fallo solo se hai almeno 1 positivo, nel training!
      if( sum(trainData$s_outcome==1)!=0 & sum(testData$s_outcome==1)!=0 &
          sum(trainData$s_outcome==0)!=0 & sum(testData$s_outcome==0)!=0 ) { 
        modello     <- glm(tmp.formula , family=binomial(link='logit'),data=trainData)
        # cat("\n",sum(trainData$s_outcome==1),"=>",sum(testData$s_outcome==1))
        # cat("\n",sum(trainData$s_outcome==0),"=>",sum(testData$s_outcome==0))
        result[[i]] <- predict(modello, testData) 
        roc.modello <- roc(testData[["s_outcome"]], predict(modello,testData),ci=TRUE )
        arr.auc.modello <- c(arr.auc.modello,roc.modello$auc)
      }
      # auc.modello <- roc.modello$auc
    }
    # -fm 
    
    
    # Ora riaddestra il modello su TUTTO
    # fai la logistica e calcola prendi la ROC
    modello <- glm( tmp.formula , family=binomial(link='logit'), data=dataFrameDiTest )
    roc.modello <- roc(dataFrameDiTest[["s_outcome"]], predict(modello),ci=TRUE )
    auc.value <- roc.modello$auc; # names(auc.value) <- covariata
    # arr.auc <- c(arr.auc, auc.value)
    
    # facciamo anche la cross correlation matrix! (calibrata su uno 0.5)
    esiti <- predict(modello,type = "response")
    # applica la soglia di indecisione
    ammessi <- rep(FALSE, length(esiti))
    if(ignorance.threshold>0) { 
      ammessi[  which( (esiti>=(.5+(ignorance.threshold/2))) | ( .5-(ignorance.threshold/2) >esiti)) ] <- TRUE
    } else { ammessi <- rep(TRUE, length(esiti)) }
    
    # prosegui a calcolare la matrice di confusione
    esiti <- esiti > 0.5
    esiti[ which(esiti==FALSE) ] <- "0"
    esiti[ which(esiti==TRUE) ] <- "1"
    res.confusion.matrix <- calcola.confusion.matrix( unlist(esiti) , unlist(mm.1[,"outcome"]) , ammessi) 
    
    # restituisci le performance
    return(
      list( 
        "confusion.matrix" = res.confusion.matrix,
        "modello" = modello,
        "roc.modello" = roc.modello,
        "auc.value"=auc.value,
        "local.tmp.str.covariate"=local.tmp.str.covariate,
        "arr.auc.modello"=arr.auc.modello
      )
    )
  }  
  old.test.these.cov<-function( arr.cov ,ignorance.threshold=NA, debug.mode = FALSE ) {
    if(debug.mode == TRUE) browser()
    mm.1 <- global.extracted.cov$mm.1
    dataFrameDiTest <- global.extracted.cov$mm.1.data.frame
    if(is.na(ignorance.threshold)) ignorance.threshold <- param.ignorance.threshold
    
    arr.s.cov <- str_c("s_",arr.cov)
    if(length(arr.cov)>1) {
      s.stringa <- paste( arr.s.cov , collapse = "+") 
    } else s.stringa <- arr.s.cov
    
    local.tmp.str.covariate <- s.stringa
    
    # if(tmp.str.covariate=="") local.tmp.str.covariate<-str_c("s_",covariata)
    # else local.tmp.str.covariate <- str_c(tmp.str.covariate,"+", str_c("s_",covariata ) ) 
    tmp.formula <- str_c("s_outcome~",s.stringa)
    
    # tmp.formula <- stringa
    tmp.formula <- as.formula(tmp.formula)

    # fai la logistica e calcola prendi la ROC
    modello <- glm( tmp.formula , family=binomial(link='logit'), data=dataFrameDiTest )
    roc.modello <- roc(dataFrameDiTest[["s_outcome"]], predict(modello),ci=TRUE )
    auc.value <- roc.modello$auc; # names(auc.value) <- covariata
    # arr.auc <- c(arr.auc, auc.value)
    
    # facciamo anche la cross correlation matrix! (calibrata su uno 0.5)
    esiti <- predict(modello,type = "response")
    # applica la soglia di indecisione
    ammessi <- rep(FALSE, length(esiti))
    if(ignorance.threshold>0) { 
      ammessi[  which( (esiti>=(.5+(ignorance.threshold/2))) | ( .5-(ignorance.threshold/2) >esiti)) ] <- TRUE
    } else { ammessi <- rep(TRUE, length(esiti)) }
    
    # prosegui a calcolare la matrice di confusione
    esiti <- esiti > 0.5
    esiti[ which(esiti==FALSE) ] <- "0"
    esiti[ which(esiti==TRUE) ] <- "1"
    res.confusion.matrix <- calcola.confusion.matrix( unlist(esiti) , unlist(mm.1[,"outcome"]) , ammessi) 

    # restituisci le performance
    return(
      list( 
        "confusion.matrix" = res.confusion.matrix,
        "modello" = modello,
        "roc.modello" = roc.modello,
        "auc.value"=auc.value,
        "local.tmp.str.covariate"=local.tmp.str.covariate
      )
    )
  }
  plot.roc <- function( roc, col="back" ){ 
    

  }
  # --------------------------------------------------
  # calcola.confusion.matrix
  # Presi in ingresso due array, calcola la matrice di confusione. Accetta in ingresso anche
  # un arreay di booleani che indica quali elementi devo essere considerati e quali no
  # --------------------------------------------------
  calcola.confusion.matrix<-function( arr.1, arr.2, arr.ammessi ) {
    arr.1 <- arr.1[ which(arr.ammessi==TRUE) ]
    arr.2 <- arr.2[ which(arr.ammessi==TRUE) ]    
    conf.matrix = table( arr.1 , arr.2 )
    
    aaa <- matrix(c(0,0,0,0),ncol=2)
    colnames(aaa) <- c("0","1")
    rownames(aaa) <- c("0","1")
    if( ("0" %in% colnames(conf.matrix))  & ("0" %in% rownames(conf.matrix)) ) aaa["0","0"] <-  conf.matrix["0","0"]
    if( ("0" %in% colnames(conf.matrix))  & ("1" %in% rownames(conf.matrix)) ) aaa["1","0"] <-  conf.matrix["1","0"]
    if( ("1" %in% colnames(conf.matrix))  & ("0" %in% rownames(conf.matrix)) ) aaa["0","1"] <-  conf.matrix["0","1"]
    if( ("1" %in% colnames(conf.matrix))  & ("1" %in% rownames(conf.matrix)) ) aaa["1","1"] <-  conf.matrix["1","1"]
    conf.matrix <- aaa    
    
    accuracy <- sum(diag(conf.matrix)) / sum(conf.matrix)
    true.positive.rate <- conf.matrix["1","1"] / sum(conf.matrix[,"1"])
    true.negative.rate <- conf.matrix["0","0"] / sum(conf.matrix[,"0"])
    
    return( list( "conf.matrix" = conf.matrix, 
                  "accuracy" = accuracy,
                  "true.positive.rate" = true.positive.rate,
                  "true.negative.rate" = true.negative.rate) )
  }
  # --------------------------------------------------
  # posiziona.nuovo.caso
  # posiziona un nuovo caso nello spazio ed esegue la predizione
  # il posizionamento avviene tramite piu' modelli, e le strategie sono 
  # le seguenti: 'AUC', 'accuracy',''sensitivity','specificity'
  # --------------------------------------------------
  posiziona.nuovo.caso <- function( objGeoLet, lst.clinical.cov=list() , models ) {

    # eredita alcuni parametri dai lanci precedenti
    ROIName <- global.ROIName
    feature.family <<- global.feature.family
    arr.sigma <<- global.arr.sigma
    
    # inizializza la lista di matrici che conterranno il risultato
    big.matrice<-list()
    if(file.exists("matricissima.RData")) load(file = "matricissima.RData")
    else {  
      # looppa per ogni sigma
      for( sigma.value in arr.sigma) { 
        # prepara la filtering pipeline
        filterPipeline<- list()
        filterPipeline[[1]]<-list("kernel.type"="LoG", "sigma"=sigma.value)    
        # Estrai le features dal paziente proposto
        tmp.matrice <- computeFeatures.geoLet( obj.geoLet = objGeoLet,
                                ROIName = ROIName , feature.family=feature.family, 
                                filterPipeline=c(), px = "", py ="", error = .1) 
        big.matrice[[as.character(sigma.value)]]<-tmp.matrice      
      }
    }
    

    # Estrai il nome delle covariate dal modello, per andare a riprenderle opportunamente con
    # aaa <- estrai.nome.covariate.e.sigma.da.stringa( models$tmp.str.covariate )

    # browser()
    arr.covariate<-c()
    for( sigma in names(big.matrice) ) { 
      nuovi.nomi <- str_c(  rep( str_c(sigma,"_"),length(names(big.matrice[[ sigma ]])) ), names(big.matrice[[ sigma ]] ))
      da.accodare <- unlist(big.matrice[[ sigma ]])
      names(da.accodare) <- unlist(nuovi.nomi)
      arr.covariate <- c(arr.covariate, da.accodare)
    } 
    arr.covariate <- c(arr.covariate, unlist(lst.clinical.cov))
    
    # Ora manipola il risultato per andare ad estrarre i valori delle covariate indicate in aaa
    # dalla matricona delle features di quello specifico paziente
    # arr.covariate <- c()
    # for( riga in seq(1,nrow(aaa))) {
    #   if( !is.na(aaa[riga,1]) ) {
    #     tmpValore <- big.matrice[[ aaa[riga,1] ]][  aaa[riga,2] ]
    #   }
    #   else { 
    #     tmpValore <- lst.clinical.cov[[ aaa[riga,2] ]]
    #     names(tmpValore)<- aaa[riga,2]
    #   }
    #   arr.covariate <- c(arr.covariate,  tmpValore )
    # }
    
    # browser()
    # Ora ordiniamo tutti i modelli in funzione delle performances!
    win.model.acc <- names(models$arr.accuracy)[  sort(models$arr.accuracy,decreasing = T, index.return = T)$ix ][1]
    win.model.sen <- names(models$arr.sensitivity)[  sort(models$arr.sensitivity,decreasing = T, index.return = T)$ix ][1]
    win.model.spec <- names(models$arr.specificity)[  sort(models$arr.specificity,decreasing = T, index.return = T)$ix ][1]
    win.model.auc <- names(models$arr.AUC)[  sort(models$arr.AUC,decreasing = T, index.return = T)$ix ][1]
    
    # Roba vecchia, legata ad un solo obiettivo
    # Ora applica il modello "vincente" e vediamo quale e' la probabilita'
    names(arr.covariate) <- unlist(str_split(string = win.model.auc,pattern = "\\+"))
    matrice.da.predire <- c()
    matrice.da.predire <- rbind(arr.covariate,arr.covariate)
    nuovoDato <- data.frame(matrice.da.predire)
    probabilita <- predict(models$lista.modelli[[ models$tmp.str.covariate ]]$modello, newdata = nuovoDato,"response")
    probabilita <- probabilita[1]
    # browser()
    if(probabilita<=0.5) classe = 0
    else classe = 1
    if(probabilita>(0.5-param.ignorance.threshold) & probabilita<(0.5+param.ignorance.threshold)) classe = NA
    
    return(list( 
      "p" = probabilita,
      "classe" = classe,
      "modello" = win.model.auc
    ))
  }
  estrai.nome.covariate.e.sigma.da.stringa<-function( stringa ) {
    
    covariate.sporche <- unlist((str_split(string = stringa,"\\+")))
    covariate <- unlist(lapply(covariate.sporche, function(x){ str_sub(x,3,str_length(x)) } ) )    
    aaa <- unlist( lapply(covariate, 
                          function(x){ 
                            posizione <- str_locate(string = x,pattern = "_")[[1]][1]; 
                            tmp.sigma.val <- str_sub(x,0, posizione-1)
                            tmp.nome.var <- str_sub(x,posizione+1,str_length(x))
                            if( is.na(tmp.nome.var)) tmp.nome.var <- x
                            return( list("sigma"=tmp.sigma.val,"nome"=tmp.nome.var ) ) 
                          } ) )
    aaa <- matrix(aaa,ncol=2,byrow = T)
    rownames(aaa) <- covariate.sporche
    return(aaa)
  }
  set.param<-function( cross.correlation.threshold = .4, 
                       wilcoxon.pValue.threshold = .05,
                       ignorance.threshold=.00)  {
    param.cross.correlation.threshold <<- cross.correlation.threshold
    param.wilcoxon.pValue.threshold <<- wilcoxon.pValue.threshold
    param.ignorance.threshold <<- ignorance.threshold    
  }
  constructor<-function(path) {
    set.param()
    
    global.ROIName <<- NA
    global.feature.family <<- c("stat","morph","glcm","rlm","szm")
    global.arr.sigma <<- c(0.49, 0.59, 0.69, 0.79)
    global.extracted.cov <<- list()
    # param.ignorance.threshold <<- 0
  }
  constructor(path = path)
  
  return(
    list(
      "scoutFeatures"=scoutFeatures,
      "featureSelection"=featureSelection,
      "posiziona.nuovo.caso" = posiziona.nuovo.caso,
      "test.these.cov" = test.these.cov,
      "set.param"=set.param,
      "plot.roc"=plot.roc
      )
    )
  
}

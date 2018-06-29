#' class for handling logs/warnings/errorss
#' 
#' @description  It handles messages from script to a chosen output (screen, file, etc.)
#' @export
#' @import stringr XML 
logHandler<-function() {
  sendLog<-function( msg ) {
    if(length(msg)>1) messaggio<-paste(msg,collapse='')
    else messaggio<-msg
    cat("\n",messaggio)
  }
  handle<-function( type , msg) {
    if(length(msg)>1) messaggio<-paste(msg,collapse='')
    else messaggio<-msg    
    cat("\n",type, ":: ", messaggio)
    if(type=="error") stop()
  }
  costructor<-function() {
  }
  costructor();
  return(
    list(
      "sendLog"=sendLog,
      "handle"=handle
    )
  )
}



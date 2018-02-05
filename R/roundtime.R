#' Title
#'
#' @param sec 
#'
#' @return
#' @export
#'
#' @examples
roundtime<-function(sec)
{
  sec<-round(sec)
  if (sec==0)return("<1 seconds")
  min<-0
  h<-0
  d<-0
  if (sec>59){
    min<-floor(sec/60)
    sec<-sec%%60
  }
  if (min>59){
    h<-floor(min/60)
    min<-min%%60
  }
  if (h>23){
    d<-floor(h/24)
    h<-h%%60
  }
  if (sec>1)result <- paste(sec,"seconds")
  if (sec==1)result <- paste(sec,"second")
  if (sec==0)result <- ""
  if (min>1)result <- paste(min,"minutes",result)
  if (min==1)result <- paste(min,"minute",result)
  if (h>1)result <- paste(h,"hours",result)
  if (h==1)result <- paste(h,"hour",result)
  if (d>1)result <- paste(h,"days",result)
  if (d==1)result <- paste(h,"day",result)
  return(trimws(result))
}

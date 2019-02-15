#' Critical Outlier Identification Function
#'
#' This function identifies the critical outliers in the predictor variable. Critical outliers are defined as data points which would maximize the accuracy of the simple linear regression model when handled
#' @param
#' x    Predictor variable
#' @param
#' y    Response Variable
#' @param
#' n 	Number of critical outliers to be identified
#' m 	Residual Error threshold. Takes values greater than 0, ideally between 0 and 1
#' @keywords Critical Outlier Identification
#' @export
#' @examples
#' CriticalOutliers(predictor_variable,response_variable,4,0.1)
CriticalOutliers <- function(x,y,n,m) {
  outinput <- data.frame(x=x, y=y)
  mod <- lm(y ~ x)
  outinput$Diff = abs(mod$residuals)
  dv_threshold = mean(y)*m
  outinput$outlier<-ifelse(outinput$Diff>dv_threshold,1,0)
  s=sum(outinput$outlier)
  prev=0
  for(i in 1:nrow(outinput))
  {
    if(outinput[i,4]!=0)
    {
      outinput[i,4]=outinput[i,4]+prev
      prev=outinput[i,4]
    }
  }
  outinput_new=outinput
  prevnum=0
  for(i in n){
    num= dim(combn(s,i))[2]+prevnum
    prevnum=num

  }
  df=as.data.frame(matrix(nrow=num,ncol=n+4))
  names(df)=c("Correlation","RSQ","Rownum","Number of Removed Outlier")

  count = 1
  count1=0

  for(z in n)
  {
    m= as.data.frame(t(combn(s,z)))
    for(i in 1:nrow(m))
    {
      count1= count1+1
      for(j in 1:ncol(m))
      {
        for (k in 1:nrow(outinput_new))
        {

          if(outinput_new[k,4]==m[i,j])
          {
            temp1 = as.character(outinput_new[k,1])
            outinput_new=outinput_new[-c(k),]
            break
          }
        }
        df[count1,j+4]=temp1
        colnames(df)[j+4] <- paste("Removed values",j,sep=" ")
      }
      corr=cor(outinput_new[,2],outinput_new[,1])
      mod=lm(outinput_new[,2]~outinput_new[,1],data=outinput_new)
      df[count1,"Correlation"]=corr
      df[count1,"RSQ"]=corr^2
      df[count1,"Rownum"]=i
      df[count1,"Number of Removed Outlier"]=z
      outinput_new = outinput
    }
  }
  View(df)
}

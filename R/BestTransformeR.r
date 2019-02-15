#' Best data transformation identification function
#'
#' This function identifies the best data transformation for one variable to maximize the strength of its relationship with another variable. The different transformation techniques compared in this function are MinMax, Zscore, Square, Cube, Square Root, Cube Root,Box-Cox, Exponential & Arcsine and Growth Models: Gompertz, Chapman Richards, Weibull & Morgan-Mercer-Flodin. The best transformation technique is identified based on maximum R Squared Value.
#' @param
#' x    Non-Normalized input variable that needs to be transformed
#' @param
#' y    Response Variable
#' @param
#' type Character string specifying the type of data transformation to use. Value "all" considers all the data transformation techniques. "Non-growth" – considers all the data transformation techniques other than growth models. "growth" - considers all the growth model techniques
#' @keywords Transformation Technique
#' @export
#' @examples
#' BestTranform(input_variable,response_variable,'growth')
#' BestTranform(input_variable,response_variable,'all')
#' BestTranform(input_variable,response_variable,'non-growth')
#' @section Warning:
#' Normalized x variable shouldn’t be passed to the function as normalization happens within the function for relevant transformation techniques; Do not pass the x variable with zero value to this function.

BestTransformeR <- function(x,y,type) {
  i=1
  #Defining Input Dataset
  input1 <- as.data.frame(matrix(nrow=length(x),ncol = 20))
  result <- data.frame()
  if(type=='all')
  {
    #Min Max Transformation
    min=min(x)
    max=max(x)
    input1[,i]=(x-min)/(max-min)
    result[i,1]="MinMax"
    #capturing RSQ
    result[i,2]=summary(lm(y~input1[,i]))$r.squared
    names(result)=c("Transformation","Value")
    i=i+1

    #zscore Transformation
    input1[,i]=(x-mean(x))/sd(x)
    result[i,1]="ZScore"
    #capturing RSQ
    result[i,2]=summary(lm(y~input1[,i]))$r.squared
    names(result)=c("Transformation","Value")
    i=i+1

    #Square root transformation
    input1[,i]=sqrt(x)
    result[i,1]="Square Root"
    #capturing RSQ
    result[i,2]=summary(lm(y~input1[,i]))$r.squared
    i=i+1

    #Cube root transformation
    input1[,i]=x^(1/3)
    result[i,1]="Cube Root"
    #capturing RSQ
    result[i,2]=summary(lm(y~input1[,i]))$r.squared
    i=i+1

    #"x^2" : x square transformation
    input1[,i]=x^2
    result[i,1]="Square"
    #capturing RSQ
    result[i,2]=summary(lm(y~input1[,i]))$r.squared
    i=i+1

    #"x^3" : x^3 square transformation
    input1[,i]=x^3
    result[i,1]="Cube"
    #capturing RSQ
    result[i,2]=summary(lm(y~input1[,i]))$r.squared
    i=i+1

    #sin Transformation
    input1[,i]=sin(x)
    result[i,1]="Sin"
    #capturing RSQ
    result[i,2]=summary(lm(y~input1[,i]))$r.squared
    i=i+1

    #cos Transformation
    input1[,i]=cos(x)
    result[i,1]="Cos"
    #capturing RSQ
    result[i,2]=summary(lm(y~input1[,i]))$r.squared
    i=i+1

    #Exponential Transformation
    #Need to pass normalized numbers as input
    input1[,i]=exp(input1[,1])
    result[i,1]="Exponential"
    #capturing RSQ
    result[i,2]=summary(lm(y~input1[,i]))$r.squared
    i=i+1

    #Arcsine transformation
    #Need to pass normalized numbers as input. Suitable for percentages or proportions
    input1[,i]=asin((input1[,1]))
    result[i,1]="Arcsine"
    #capturing RSQ
    result[i,2]=summary(lm(y~input1[,i]))$r.squared
    i=i+1

    ##Box-Cox transformation
    Box = boxcox(x ~ 1,              # Transform Turbidity as a single vector
                 lambda = seq(-6.0,6.0,0.1)      # Try values -6 to 6 by 0.1
    )
    Cox = data.frame(Box$x, Box$y)            # Create a data frame with the results
    Cox2 = Cox[with(Cox, order(-Cox$Box.y)),] # Order the new data frame by decreasing y
    lambda = Cox2[1, "Box.x"]                 # Extract that lambda
    input1[,i] = (x ^ lambda - 1)/lambda   # Transform the original data
    result[i,1]="Box-Cox"
    #capturing RSQ
    result[i,2]=summary(lm(y~input1[,i]))$r.squared
    i=i+1
    dat <- data.frame(x=input1[,1], y=y)
    ##### Chapman Richards
    pars <- expand.grid(a=seq(1:20),
                        b=seq(0.1,0.9,by=0.1),
                        c=seq(0.1,0.9,by=0.1),
                        m=seq(0.1,0.9, by=0.1))
    fo = y ~ a * (1 - b * exp(-c * x))^(1 / (1 - m))
    res = NULL
    res <- nls2(fo, data=dat,trace=F,
                start=pars, algorithm='brute-force')
    rss <- sum((predict(res) - y) ^ 2)  ## residual sum of squares
    tss <- sum((y - mean(y)) ^ 2)  ## total sum of squares
    rsq <- 1 - rss/tss
    result[i,1]="Chapman Richards"
    result[i,2]=rsq
    mat=coef(res)
    i=i+1

    ### Gompertz
    pars <- expand.grid(a=seq(1,20),b=seq(1, 10),c=seq(1, 10))
    fo = y ~ a * exp(-b * exp(-c * x))
    res = NULL
    res <- nls2(fo, data=dat,trace=F,
                start=pars, algorithm='brute-force')
    rss <- sum((predict(res) - y) ^ 2)  ## residual sum of squares
    tss <- sum((y - mean(y)) ^ 2)  ## total sum of squares
    rsq <- 1 - rss/tss
    result[i,1]="Gompertz"
    result[i,2]=rsq
    mat=coef(res)
    i=i+1

    ### Weibull
    pars <- expand.grid(a=seq(1,20,by=1),
                        b=seq(1, 10),
                        c=seq(1, 10),
                        m=seq(0.1, 0.9, by=0.1))
    fo = y ~ a - b * exp(-c * x^m)
    res = NULL
    res <- nls2(fo, data=dat,trace=F,
                start=pars, algorithm='brute-force')
    rss <- sum((predict(res) - y) ^ 2)  ## residual sum of squares
    tss <- sum((y - mean(y)) ^ 2)  ## total sum of squares
    rsq <- 1 - rss/tss
    result[i,1]="Weibull"
    result[i,2]=rsq
    mat=coef(res)
    i=i+1

    ### Morgan-Mercer-Flodin
    pars <- expand.grid(w0=seq(0.1,0.9,by=0.1),
                        gamma=seq(1, 10, by=1),
                        a=seq(1:20),
                        m=seq(0.1, 1, by=0.1))
    fo = y ~ (w0 * gamma + a * y^m) / (gamma + y^m)
    res = NULL
    res <- nls2(fo, data=dat,trace=F,
                start=pars, algorithm='brute-force')
    rss <- sum((predict(res) - y) ^ 2)  ## residual sum of squares
    tss <- sum((y - mean(y)) ^ 2)  ## total sum of squares
    rsq <- 1 - rss/tss
    result[i,1]="Morgan-Mercer-Flodin"
    result[i,2]=rsq
    mat=coef(res)
    i=i+1

    #Finding Maximum Value
    y_max <- max(result[,2])
    # Corresponding transformation value
    x_val_associated <- result[result$Value == y_max, "Transformation"]
    #Printing Final Result
    if(x_val_associated=="Chapman Richards"||x_val_associated=="Gompertz"||x_val_associated=="Weibull"||x_val_associated=="Morgan-Mercer-Flodin")
    {
      print(x_val_associated)
      print(y_max)
      print(mat)
    }else{
      print(x_val_associated)
      print(y_max)
    }
  }
  if(type=='non-growth')
  {
    #Min Max Transformation
    min=min(x)
    max=max(x)
    input1[,i]=(x-min)/(max-min)
    result[i,1]="MinMax"
    #capturing RSQ
    result[i,2]=summary(lm(y~input1[,i]))$r.squared
    names(result)=c("Transformation","Value")
    i=i+1

    #zscore Transformation
    input1[,i]=(x-mean(x))/sd(x)
    result[i,1]="ZScore"
    #capturing RSQ
    result[i,2]=summary(lm(y~input1[,i]))$r.squared
    names(result)=c("Transformation","Value")
    i=i+1

    #Square root transformation
    input1[,i]=sqrt(x)
    result[i,1]="Square Root"
    #capturing RSQ
    result[i,2]=summary(lm(y~input1[,i]))$r.squared
    i=i+1

    #Cube root transformation
    input1[,i]=x^(1/3)
    result[i,1]="Cube Root"
    #capturing RSQ
    result[i,2]=summary(lm(y~input1[,i]))$r.squared
    i=i+1

    #"x^2" : x square transformation
    input1[,i]=x^2
    result[i,1]="Square"
    #capturing RSQ
    result[i,2]=summary(lm(y~input1[,i]))$r.squared
    i=i+1

    #"x^3" : x^3 square transformation
    input1[,i]=x^3
    result[i,1]="Cube"
    #capturing RSQ
    result[i,2]=summary(lm(y~input1[,i]))$r.squared
    i=i+1

    #sin Transformation
    input1[,i]=sin(x)
    result[i,1]="Sin"
    #capturing RSQ
    result[i,2]=summary(lm(y~input1[,i]))$r.squared
    i=i+1

    #cos Transformation
    input1[,i]=cos(x)
    result[i,1]="Cos"
    #capturing RSQ
    result[i,2]=summary(lm(y~input1[,i]))$r.squared
    i=i+1

    #Exponential Transformation
    #Need to pass normalized numbers as input
    input1[,i]=exp(input1[,1])
    result[i,1]="Exponential"
    #capturing RSQ
    result[i,2]=summary(lm(y~input1[,i]))$r.squared
    i=i+1

    #Arcsine transformation
    #Need to pass normalized numbers as input. Suitable for percentages or proportions
    input1[,i]=asin((input1[,1]))
    result[i,1]="Arcsine"
    #capturing RSQ
    result[i,2]=summary(lm(y~input1[,i]))$r.squared
    i=i+1

    ##Box-Cox transformation
    Box = boxcox(x ~ 1,              # Transform Turbidity as a single vector
                 lambda = seq(-6.0,6.0,0.1)      # Try values -6 to 6 by 0.1
    )
    Cox = data.frame(Box$x, Box$y)            # Create a data frame with the results
    Cox2 = Cox[with(Cox, order(-Cox$Box.y)),] # Order the new data frame by decreasing y
    lambda = Cox2[1, "Box.x"]                 # Extract that lambda
    input1[,i] = (x ^ lambda - 1)/lambda   # Transform the original data
    result[i,1]="Box-Cox"
    #capturing RSQ
    result[i,2]=summary(lm(y~input1[,i]))$r.squared

    #Finding Maximum Value
    y_max <- max(result[,2])
    # Corresponding transformation value
    x_val_associated <- result[result$Value == y_max, "Transformation"]
    #Printing Final Result
    print(x_val_associated)
    print(y_max)
  }
  if(type=='growth')
  {
    #Min Max Transformation
    min=min(x)
    max=max(x)
    input1[,i]=(x-min)/(max-min)
    dat <- data.frame(x=input1[,1], y=y)
    ##### Chapman Richards
    pars <- expand.grid(a=seq(1:20),
                        b=seq(0.1,0.9,by=0.1),
                        c=seq(0.1,0.9,by=0.1),
                        m=seq(0.1,0.9, by=0.1))
    fo = y ~ a * (1 - b * exp(-c * x))^(1 / (1 - m))
    res = NULL
    res <- nls2(fo, data=dat,trace=F,
                start=pars, algorithm='brute-force')
    rss <- sum((predict(res) - y) ^ 2)  ## residual sum of squares
    tss <- sum((y - mean(y)) ^ 2)  ## total sum of squares
    rsq <- 1 - rss/tss
    result[i,1]="Chapman Richards"
    result[i,2]=rsq
    mat=coef(res)
    names(result)=c("Transformation","Value")
    i=i+1

    ### Gompertz
    pars <- expand.grid(a=seq(1,20),b=seq(1, 10),c=seq(1, 10))
    fo = y ~ a * exp(-b * exp(-c * x))
    res = NULL
    res <- nls2(fo, data=dat,trace=F,
                start=pars, algorithm='brute-force')
    rss <- sum((predict(res) - y) ^ 2)  ## residual sum of squares
    tss <- sum((y - mean(y)) ^ 2)  ## total sum of squares
    rsq <- 1 - rss/tss
    result[i,1]="Gompertz"
    result[i,2]=rsq
    mat=coef(res)
    i=i+1

    ### Weibull
    pars <- expand.grid(a=seq(1,20,by=1),
                        b=seq(1, 10),
                        c=seq(1, 10),
                        m=seq(0.1, 0.9, by=0.1))
    fo = y ~ a - b * exp(-c * x^m)
    res = NULL
    res <- nls2(fo, data=dat,trace=F,
                start=pars, algorithm='brute-force')
    rss <- sum((predict(res) - y) ^ 2)  ## residual sum of squares
    tss <- sum((y - mean(y)) ^ 2)  ## total sum of squares
    rsq <- 1 - rss/tss
    result[i,1]="Weibull"
    result[i,2]=rsq
    mat=coef(res)
    i=i+1

    ### Morgan-Mercer-Flodin
    pars <- expand.grid(w0=seq(0.1,0.9,by=0.1),
                        gamma=seq(1, 10, by=1),
                        a=seq(1:20),
                        m=seq(0.1, 1, by=0.1))
    fo = y ~ (w0 * gamma + a * y^m) / (gamma + y^m)
    res = NULL
    res <- nls2(fo, data=dat,trace=F,
                start=pars, algorithm='brute-force')
    rss <- sum((predict(res) - y) ^ 2)  ## residual sum of squares
    tss <- sum((y - mean(y)) ^ 2)  ## total sum of squares
    rsq <- 1 - rss/tss
    result[i,1]="Morgan-Mercer-Flodin"
    result[i,2]=rsq
    mat=coef(res)
    i=i+1

    #Finding Maximum Value
    y_max <- max(result[,2])
    # Corresponding transformation value
    x_val_associated <- result[result$Value == y_max, "Transformation"]
    #Printing Final Result
    print(x_val_associated)
    print(y_max)
    print(mat)
  }
}

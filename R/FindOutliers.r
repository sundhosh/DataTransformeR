#' Outlier identification Method
#'
#' This function identifies the outliers in the predictor variable based on residual error from a simple linear regression model.
#' @param
#' x    Predictor variable
#' @param
#' y    Response Variable
#' @param
#' n 	Residual Error threshold. Takes values greater than 0, ideally between 0 and 1
#' @keywords Outlier Detection
#' @export
#' @examples
#' FindOutliers(predictor_variable,response_variable,0.1)

FindOutliers <- function(x,y,n) {
outinput1 <- data.frame()
mod <- lm(y ~ x)
outinput1 = abs(mod$residuals)
dv_threshold = mean(y)*n
outinput2<-ifelse(outinput1>dv_threshold,1,0)
output = outinput1[ outinput2 == 1]
print(output)
}
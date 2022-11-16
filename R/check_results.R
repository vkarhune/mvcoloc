#' Helper function to check multivariable colocalization results
#'
#' @param res
#' @param type
#' @param thresh1
#' @param thresh2
#'
#' @return
#' @export
#'
#' @examples
check_results <- function(res, type = "cont", thresh1 = 0.5, thresh2 = 0.5){

  x1y <- switch(type, "cont" = "x1y", "bin" = "x1ybin", NULL)
  x2y <- switch(type, "cont" = "x2y", "bin" = "x2ybin", NULL)
  x1y_adjx2 <- switch(type, "cont" = "x1y_adjx2", "bin" = "x1ybin_adjx2", NULL)
  x2y_adjx1 <- switch(type, "cont" = "x2y_adjx1", "bin" = "x2ybin_adjx1", NULL)


  if(res[[x1y]][2] < thresh1 & res[[x2y]][2] < thresh1){
    cat(sprintf("No marginal colocalization\n"))
    out <- "nomargins"
  }

  if(res[[x1y]][2] > thresh1 & res[[x2y]][2] < thresh1){
    if(res[[x2y_adjx1]][2] < thresh2){
      cat(sprintf("X1 colocalizes, X2 adjusted does not colocalize\n"))
      out <- "x1only"
    } else {
      # cat(sprintf("X1 colocalizes, X2 adjusted colocalizes (PP = %.2f)\n", x$x2ybin_adjx1[2]))
      cat(sprintf("X1 colocalizes, X2 adjusted colocalizes\n"))
      out <- "x1_x2cond"
    }
  }

  if(res[[x1y]][2] < thresh1 & res[[x2y]][2] > thresh1){
    if(res[[x1y_adjx2]][2] < thresh2){
      cat(sprintf("X2 colocalizes, X1 adjusted does not colocalize\n"))
      out <- "x2only"
    } else {
      # cat(sprintf("X1 colocalizes, X2 adjusted colocalizes (PP = %.2f)\n", x$x2ybin_adjx1[2]))
      cat(sprintf("X2 colocalizes, X1 adjusted colocalizes\n"))
      out <- "x2_x1cond"
    }
  }

  if(res[[x1y]][2] > thresh1 & res[[x2y]][2] > thresh1){
    cat(sprintf("Both marginals colocalize\n"))
    out <- "bothmargins"
  }

  return(out)
}


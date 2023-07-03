#' UKL colors
#'
#' @description Provides a color palette based on the UKL corporate identity
#'
#' @param n number of colors to return
#'
#' @return list of n colors
#' @export
#'
#' @examples ukl_colors(3)

ukl_colors = function(n){

  # these are the colors from the ukl_corporate identity
  turqoise = rgb(0, 138, 201, maxColorValue = 255)
  darkblue = rgb(36, 51, 87, maxColorValue = 255)
  plum = rgb(163, 0, 80, maxColorValue = 255)
  lila = rgb(128, 100, 162, maxColorValue = 255)
  aquamarin = rgb(75, 172, 198, maxColorValue = 255)
  # maybe plum and lila need to be switched

  ukl_color_function = colorRampPalette(c(aquamarin, turqoise, plum, lila, darkblue))
  ukl_colors = ukl_color_function(n)

  return(ukl_colors)
}

#' check if null
#'
#' @description Checks if input variable is null. If not, returns variable as symbol.
#' This is needed for the evaluation within ggplot2, where NULL can be used as placeholders.
#'
#' @param input input variable. String or NULL
#'
#' @return If input is NULL, returns NULL. If input is string, returns input as symbol
#' (like "get(input)", but better for ggplot)
#' @export
#'
#' @examples
#' checkifnull(NULL)

# helper function to evaluate NULL as placeholders
# if not NULL then return the variable as a symbol (like "get()")
checkifnull = function(input){
  if(is.null(input)){
    return(NULL)
  } else {
    return(sym(input))
  }
}

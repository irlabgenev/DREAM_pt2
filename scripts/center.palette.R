#' Calculates a color palette where zero values are centered at white.
#'
#' @param data A matrix of integers.
#' @param palette_length The size of the color palette. Default = 100
#' @param color1 The color for the low bound. Default = "blue"
#' @param color2 The upper bound color.  Default = "red"
#' @keywords color, heatmap, palette
#' @return A vector containing the new color palette.
#' @examples
#' my_palette <- center.palette(data, palette_length = 200, color1 = "black", color2 = "pink")
#' @export center.palette
#' @author Amy L. Olex \email{alolex@@vcu.edu}
#'
center.palette <- function(data, palette_length = 100, color1 = "blue", color2 = "red"){
  
  my_range <- range(data)
  if((my_range[1] >= 0) | (my_range[2] <= 0)){
    stop("Range does not cross zero. Cannot center color palette.")
  }
  
  diff <- my_range[2] - my_range[1]
  length1 <- floor((abs(my_range[1])/diff)*palette_length)
  length2 <- floor((abs(my_range[2])/diff)*palette_length)
  
  
  color1=colorRampPalette(colors=c(color1, "white"))(length1)
  color2=colorRampPalette(colors=c("white", color2))(length2)
  my_colors <- c(color1, color2)
  
  my_breaks <- seq(from = range(data)[1], to = range(data)[2], length.out = length(my_colors)+1)
  
  return(list(colors = my_colors, breaks = my_breaks))
}

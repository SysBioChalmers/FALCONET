#' This function is used to adjust size of non base metabolite in the map
#'
#'
#' @param positionID the position of the component
#' @param uniqueID the unique ID of the metabolite
#' @param onemap a map in xml format
#'
#' @return
#' @export onemap a refined map
#'
#' @examples
changeFontSize <- function(positionID, uniqueID, onemap) {

  index <- which(str_detect(onemap, "celldesigner:speciesAlias") & str_detect(onemap, positionID) & str_detect(onemap, uniqueID))
  # bounds line
  pos_bounds <- index + 2
  # size line
  pos_size <- index + 3
  # box size line
  pos_box <- index + 7
  # change the size of metabolite component in the graph
  onemap[pos_bounds] <- str_replace_all(onemap[pos_bounds], "w=\"70.0\" h=\"25.0\"/>", "w=\"35.0\" h=\"12.0\"/>")
  onemap[pos_size] <- str_replace_all(onemap[pos_size], "<celldesigner:font size=\"12\"/>", "<celldesigner:font size=\"8\"/>")
  onemap[pos_box] <- str_replace_all(onemap[pos_box], "<celldesigner:boxSize width=\"70.0\" height=\"25.0\"/>", "<celldesigner:boxSize width=\"35.0\" height=\"12.0\"/>")
  return(onemap)

}


#' function to define the color based on the fold change
#'
#' @param omic_fold a vector contains the fold change values of omics data
#' @param up a value define the significantly increased fold changes
#' @param down a value define the significantly decreased fold changes
#'
#' @return
#' @export color a vector contains the color of each reaction based on the flux fold changes
#'
#' @examples
defineFluxFoldColor <- function(omic_fold, up=2, down=0.5) {
  color <- vector()
  for (i in seq_along(omic_fold)) {
    if (omic_fold[i] >= up) {
      color[i] <- "ffff0000" # red
    } else if (omic_fold[i] < down) {
      color[i] <- "ff0eb10e" # green
    } else {
      color[i] <- "ff000000"
    }
  }

  return(color)
}







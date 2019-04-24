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

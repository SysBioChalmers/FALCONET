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



#' function to define the color based on the fold change
#'
#' @param omic_fold a vector contains the fold change values of omics data
#' @param up a value define the significantly increased fold changes
#' @param down a value define the significantly decreased fold changes
#'
#' @return
#' @export color a vector contains the color of each reaction based on the transcriptomics or protein level fold changes
#'
#' @examples
defineGeneFoldColor <- function(omic_fold, up=2, down=0.5) {
  color <- vector()
  for (i in seq_along(omic_fold)) {
    if (omic_fold[i] >= up) {
      color[i] <- "ffff0000" # red
    } else if (omic_fold[i] < down) {
      color[i] <- "ff0eb10e" # green
    } else {
      color[i] <- "ffccff66"
    }
  }

  return(color)
}



# Function to mapping the flux data onto the map
#' Title
#'
#' @param input_map diretory of xml file
#' @param flux_inf a dataframe contains the flux value for each reaction
#' @param output_map diretory of xml file
#'
#' @return
#' @export
#'
#' @examples
fluxMapping <- function(input_map, flux_inf, output_map) {
  yeast_map <- readLines(file(input_map))
  index_rxn <- which(str_detect(yeast_map, "reaction metaid"))
  rxnid <- yeast_map[index_rxn]
  rxnid0 <- str_split(rxnid, " ")
  rxnid1 <- vector()
  for (i in 1:length(rxnid0)) {
    rxnid1[i] <- rxnid0[[i]][3]
  }
  rxnid1 <- str_replace_all(rxnid1, "id=", "") %>%
    str_replace_all(., "\"", "")

  flux_data <- data.frame(rxnid = rxnid1, stringsAsFactors = FALSE)
  flux_data$value <- getSingleReactionFormula(flux_inf$flux_ref, flux_inf$Abbreviation, flux_data$rxnid)
  flux_data$value <- as.numeric(flux_data$value)

  # mapping the fluxs
  min_line_wide <- 1
  max_line_wide <- 7
  min_flux <- abs(min(flux_data$value))
  max_flux <- abs(max(flux_data$value))
  line_wide <- vector()
  print('Calculate the line width')
  for (i in 1:nrow(flux_data)) {
    s0 <- 1 + abs(flux_data$value[i]) / max_flux * max_line_wide
    print(s0)
    flux_data$line_wide[i] <- s0
  }


  # from index_rxn choose the two
  index_line <- which(str_detect(yeast_map, "celldesigner:editPoints"))
  index_line1 <- index_line + 1

  print('Start the data mapping')
  for (i in 1:length(index_line1)) {
    oldline <- "line width=\"1.0\" color=\"ff000000\"/>"
    newline <- paste("line width=\"", flux_data$line_wide[i], "\" color=\"ff000000\"/>", sep = "")
    yeast_map[index_line1[i]] <- str_replace_all(yeast_map[index_line1[i]], oldline, newline)
  }
  writeLines(yeast_map, file(output_map))
}




#' Function to mapping the flux fold changes onto the map
#'
#' @param input_map diretory of xml file
#' @param flux_inf a dataframe contains the flux fold changes value for each reaction under two different conditions
#' @param output_map diretory of xml file
#'
#' @return
#' @export
#'
#' @examples
fluxFoldMapping <- function(input_map, flux_inf, output_map) {
  yeast_map <- readLines(file(input_map))
  index_rxn <- which(str_detect(yeast_map, "reaction metaid"))
  # extract the rxns
  rxnid <- yeast_map[index_rxn]
  rxnid0 <- str_split(rxnid, " ")
  rxnid1 <- vector()
  for (i in 1:length(rxnid0)) {
    rxnid1[i] <- rxnid0[[i]][3]
  }
  rxnid1 <- str_replace_all(rxnid1, "id=", "") %>%
    str_replace_all(., "\"", "")

  flux_data <- data.frame(rxnid = rxnid1, stringsAsFactors = FALSE)
  flux_data$value <- getSingleReactionFormula(flux_inf$flux_fold, flux_inf$Abbreviation, flux_data$rxnid)
  flux_data$value <- as.numeric(flux_data$value)

  # define the line color
  line_color <- defineFluxFoldColor(omic_fold = flux_data$value)

  # from index_rxn choose the two
  index_line <- which(str_detect(yeast_map, "celldesigner:editPoints"))
  index_line1 <- index_line + 1

  for (i in 1:length(index_line1)) {
    oldline <- "line width=\"1.0\" color=\"ff000000\"/>"
    if (line_color[i] != "ff000000") {
      newline <- paste("line width=\"3.0\"", " color=\"", line_color[i], "\"/>", sep = "")
    } else {
      newline <- paste("line width=\"1.0\"", " color=\"", line_color[i], "\"/>", sep = "")
    }
    yeast_map[index_line1[i]] <- str_replace_all(yeast_map[index_line1[i]], oldline, newline)
  }

  writeLines(yeast_map, file(output_map))
}




#' Function to mapping the omics fold changes onto the map
#'
#' @param input_map diretory of xml file
#' @param flux_inf a dataframe contains the omic(gene or protein) fold changes value for each reaction under two different conditions
#' @param output_map diretory of xml file
#'
#' @return
#' @export
#'
#' @examples
omicsFoldMapping <- function(input_map, flux_inf, output_map) {
  yeast_map <- readLines(file(input_map))
  index_rxn <- which(str_detect(yeast_map, "reaction metaid"))
  # extract the rxns
  rxnid <- yeast_map[index_rxn]
  rxnid0 <- str_split(rxnid, " ")
  rxnid1 <- vector()
  for (i in 1:length(rxnid0)){
    rxnid1[i] <- rxnid0[[i]][3]
  }
  rxnid1 <- str_replace_all(rxnid1, "id=","") %>%
    str_replace_all(.,"\"","")

  flux_data <- data.frame(rxnid = rxnid1, stringsAsFactors = FALSE)
  flux_data$gene_id <- str_replace_all(flux_data$rxnid, 'r_','g_')
  flux_data$protein_id <- str_replace_all(flux_data$rxnid, 'r_','p_')
  flux_data$g_value <- getSingleReactionFormula(flux_map$gene_fold, flux_map$Abbreviation,flux_data$rxnid)
  flux_data$p_value <- getSingleReactionFormula(flux_map$protein_fold, flux_map$Abbreviation,flux_data$rxnid)

  color_gene <- defineGeneFoldColor(omic_fold = flux_data$g_value, up=1.5, down = 0.6)
  color_protein <- defineGeneFoldColor(omic_fold = flux_data$p_value)

  #"ff0eb10e\" scheme=\"Color\"/>"   # green
  #ffff0000\" scheme=\"Color\"/>"   # red
  #ffccff66 # original background

  # for gene
  # firstly we get the coordinate of gene and protein on the map based on gene_id and protein_id
  gene0 <- flux_data$gene_id
  gene_index <- vector()
  for(i in seq_along(gene0)){
    gene_index[i] <- which(str_detect(yeast_map, "<species") & str_detect(yeast_map, gene0[i]))
  }

  gene_mapid <- yeast_map[gene_index]
  gene_mapid0 <- str_split(gene_mapid, " ")

  gene_mapid1 <- vector()
  for (i in 1:length(gene_mapid0)){
    gene_mapid1[i] <- gene_mapid0[[i]][3]
  }

  gene_mapid1 <- str_replace_all(gene_mapid1, "id=","") %>%
    str_replace_all(.,"\"","")

  # get the coordinate of gene
  gene_position <- vector()
  for (i in seq_along(gene_mapid1)){
    gene_position[i] <- which(str_detect(yeast_map,"<celldesigner:speciesAlias") & str_detect(yeast_map, gene_mapid1[i]))
  }

  # change the color
  color_position <- gene_position+9 #fill color

  new_color <- paste("<celldesigner:paint color=\"",color_gene,"\" scheme=\"Color\"/>", sep = "")
  yeast_map[color_position] <- new_color

  # for protein
  # firstly we get the coordinate of protein and protein on the map based on protein_id and protein_id
  protein0 <- flux_data$protein_id
  protein_index <- vector()
  for(i in seq_along(protein0)){
    protein_index[i] <- which(str_detect(yeast_map, "<species") & str_detect(yeast_map, protein0[i]))
  }

  protein_mapid <- yeast_map[protein_index]
  protein_mapid0 <- str_split(protein_mapid, " ")

  protein_mapid1 <- vector()
  for (i in 1:length(protein_mapid0)){
    protein_mapid1[i] <- protein_mapid0[[i]][3]
  }

  protein_mapid1 <- str_replace_all(protein_mapid1, "id=","") %>%
    str_replace_all(.,"\"","")

  # get the coordinate of protein
  protein_position <- vector()
  for (i in seq_along(protein_mapid1)){
    protein_position[i] <- which(str_detect(yeast_map,"<celldesigner:speciesAlias") & str_detect(yeast_map, protein_mapid1[i]))
  }

  # change the color
  color_position <- protein_position+9 #fill color
  new_color <- paste("<celldesigner:paint color=\"",color_protein,"\" scheme=\"Color\"/>", sep = "")
  yeast_map[color_position] <- new_color

  # change all other background into the same color
  yeast_map <- str_replace_all(yeast_map, "ffccff66", "ff33ffff")

  writeLines(yeast_map, file(output_map))
}





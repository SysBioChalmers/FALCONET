# An example to show how we can map the fluxes fold changes on to the map
# the xml file obtained based on the automatic layout.
# otherwise the followed code can't be used
# input the flux data


# Function to mapping the flux fold changes onto the map
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


flux_map <- read_excel("data/flux_map.xlsx")
fluxFoldMapping(input_map ="result/model_test_check.xml",
            flux_inf = flux_map,
            output_map = "result/map_with_flux_fold_changes.xml")

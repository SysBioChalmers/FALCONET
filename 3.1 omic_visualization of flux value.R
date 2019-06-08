# An example to show how we can map the fluxes on to the map
# the xml file obtained based on the automatic layout.
# otherwise the followed code can't be used

# Function to mapping the flux data onto the map
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




# input the flux data
flux_map <- read_excel("data/flux_map.xlsx")

fluxMapping(input_map ="result/model_test_check.xml",
              flux_inf = flux_map,
              output_map = "result/map_with_flux_value.xml")
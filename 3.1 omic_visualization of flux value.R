# An example to show how we can map the fluxes on to the map
# the xml file obtained based on the automatic layout.
# otherwise the followed code can't be used

yeast_map <- readLines(file("result/model_test_check.xml"))
index_rxn <- which(str_detect(yeast_map, "reaction metaid"))
# input the flux data
flux_map <- read_excel("data/flux_map.xlsx")


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
flux_data$value <- getSingleReactionFormula(flux_map$flux_ref,flux_map$Abbreviation,flux_data$rxnid)
flux_data$value <- as.numeric(flux_data$value)

# mapping the fluxs
min_line_wide <- 1
max_line_wide <- 7
min_flux <- abs(min(flux_data$value))
max_flux <- abs(max(flux_data$value))
line_wide <- vector()
for (i in 1:nrow(flux_data)){
  s0 <- 1 + abs(flux_data$value[i])/max_flux * max_line_wide
  print(s0)
  flux_data$line_wide[i] <- s0
}


# from index_rxn choose the two
index_line <- which(str_detect(yeast_map, "celldesigner:editPoints"))
index_line1 <- index_line + 1


for (i in 1:length(index_line1)){
oldline <- "line width=\"1.0\" color=\"ff000000\"/>"
newline <- paste("line width=\"", flux_data$line_wide[i], "\" color=\"ff000000\"/>", sep = "")
yeast_map[index_line1[i]] <- str_replace_all(yeast_map[index_line1[i]], oldline, newline)

}
writeLines(yeast_map, file("result/map_with_flux_value.xml"))
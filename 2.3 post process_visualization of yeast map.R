# An example to show how we can map the fluxes on to the map
yeast_map <- readLines(file("result/model_test_check.xml"))
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
flux_data$value <- getSingleReactionFormula(rxn$Flux,rxn$Abbreviation,flux_data$rxnid)
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
writeLines(yeast_map, file("result/model_test_check.xml"))










# test code
model_visualization <- readLines(file("data/model_visualization.xml"))
#gpr_position_in_map <- gpr

# extabolish the position of gene or protein

"<celldesigner:speciesAlias id=\"sa8\" species=\"s8\">"
gene_list <- "<celldesigner:speciesAlias id=\"sa8\" species=\"s8\">"
color_position_g <- which(str_detect(model_visualization, gene_list)==TRUE)+9 #fill color

#color_position_g <- which(str_detect(model_visualization, gene_list)==TRUE)+15 #line color
model_visualization[color_position_g] <- "<celldesigner:paint color=\"ff0eb10e\" scheme=\"Color\"/>"   #green

model_visualization[color_position_g] <- "<celldesigner:paint color=\"ffff0000\" scheme=\"Color\"/>"   #green


# estabolish the postion of reaction

"<reaction metaid=\"r_0200\" id=\"r_0200\" reversible=\"false\">"
size_reaction <- which(str_detect(model_visualization,"<celldesigner:line width=\"1.0\" color=\"ff000000\"/>")==TRUE)
model_visualization[size_reaction[1]] <- "<celldesigner:line width=\"6.0\" color=\"ff000000\"/>"
model_visualization[size_reaction[4]] <- "<celldesigner:line width=\"6.0\" color=\"ff000000\"/>"
model_visualization[size_reaction[7]] <- "<celldesigner:line width=\"6.0\" color=\"ff000000\"/>"
model_visualization[size_reaction[10]] <- "<celldesigner:line width=\"6.0\" color=\"ff000000\"/>"
model_visualization[size_reaction[13]] <- "<celldesigner:line width=\"6.0\" color=\"ff000000\"/>"
model_visualization[size_reaction[16]] <- "<celldesigner:line width=\"6.0\" color=\"ff000000\"/>"
writeLines(model_visualization, file("result/model_visualization_changed_color.xml")) # save




